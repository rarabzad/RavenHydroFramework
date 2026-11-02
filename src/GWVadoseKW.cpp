/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadoseKW.cpp
  :SoilZoneVadose KINEMATIC: gravity-only kinematic wave on the unsaturated column (see GWVadose.h; design note
  soil_zone_kinematic.md). "Richards without capillarity": the column's own Campbell conductivity shifted to vanish at the
  residual (drained) content thr = ths - Sy,
      K(th) = Ks [(th/ths)^e - x0^e] / (1 - x0^e),  e = 2b + 3,  x0 = thr/ths  (0 at thr, Ks at ths),
  the reference profile is thr (gwvz::SetReference(1)). Implicit (backward Euler) first-order upwind finite volumes,
  node by node from the top: (th_i - th_i^0) dz_i/h = q_in,i - K(th_i), one monotone scalar equation per node (Newton,
  started on the side where it converges monotonically, bracketed). A saturated node passes on what it receives up to
  its Ks; the excess backs up into the nodes above, and what the top node cannot take is rejected. Substeps keep every
  content change within the accuracy (SetAccuracy, default 0.02). Mass-exact: the flux into the water table is the
  inflow less the storage change.
  ----------------------------------------------------------------*/
#include "GWVadose.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace gwvz
{
  namespace
  {
    ///< kinematic parameters of one node's material
    struct kw_par
    {
      double ths,thr,e,x0e,A,Ks; ///< saturated and residual contents, exponent, x0^e, normalising factor Ks/(1-x0^e), Ks
      bool   step;               ///< no mobile water (Sy <= 0 or a degenerate factor): K = 0 below ths, Ks at ths
    };
    kw_par Par(const vz_mat &m)
    {
      kw_par p; p.ths=std::max(m.ths,0.0); p.Ks=std::max(m.Ks,0.0); p.e=2.0*m.b+3.0;
      p.thr=std::max(p.ths-std::max(m.Sy,0.0),0.0);
      p.x0e=0.0; p.A=0.0; p.step=true;
      if ((p.ths>0.0) && (p.ths-p.thr>0.0)){                     //(Sy = 0 gives x0 = 1 and A = Ks/0)
        p.x0e=std::pow(p.thr/p.ths,p.e);
        const double d=1.0-p.x0e;
        if ((d>1e-12) && std::isfinite(p.x0e)){p.A=p.Ks/d; p.step=false;}
      }
      return p;
    }
    inline double K(const kw_par &p,const double th,double *dK)
    {
      if (p.step){if (dK){*dK=0.0;} return (th>=p.ths)?p.Ks:0.0;}
      if (th<=p.thr){if (dK){*dK=0.0;} return 0.0;}
      if (th>=p.ths){if (dK){*dK=0.0;} return p.Ks;}
      const double x=th/p.ths,xe1=std::pow(x,p.e-1.0);
      if (dK){*dK=p.A*p.e*xe1/p.ths;}
      return std::max(p.A*(xe1*x-p.x0e),0.0);
    }
    ///< content th of a node after h [d] with inflow qin [m/d] from above: root of g(th) = (th - th0) dz/h + K(th) - qin;
    ///< the node saturates when g(ths) <= 0: th = ths, the excess (qin less what it stores and passes on, Ks) is returned
    ///< in exc [m/d]. Returns the outflow [m/d] (qin - storage rate: >= 0 by the choice of the bracket)
    double SolveNode(const kw_par &p,const double th0,const double dz,const double h,const double qin,double &th,double &exc,int &nit)
    {
      exc=0.0; const double c=dz/h;
      if (p.step){ //no mobile water (Sy <= 0): stores up to ths, then passes on up to Ks (no inverse of K)
        if (qin<=(p.ths-th0)*c){th=th0+qin/c; return 0.0;}
        th=p.ths; double out=qin-(p.ths-th0)*c;
        if (out>p.Ks){exc=out-p.Ks; out=p.Ks;}
        return std::max(out,0.0);
      }
      const double K0=K(p,th0,0),g0=K0-qin;
      if ((g0==0.0) && (th0<=p.ths)){th=th0; return qin;}
      double lo,hi;
      if (g0<0.0){ //wetting
        if ((p.ths-th0)*c+p.Ks-qin<=0.0){ //saturated: stores up to ths, passes on Ks, the rest (exc) backs up
          th=p.ths; double out=qin-(p.ths-th0)*c;
          if (out>p.Ks){exc=out-p.Ks; out=p.Ks;}
          return std::max(out,0.0);
        }
        lo=th0; hi=std::min(th0+(qin-K0)/c,p.ths);           //(the explicit guess has g >= 0, K growing with th: an upper bound)
      }
      else {       //draining: the explicit guess th0 - (K0 - qin) h/dz has g <= 0: a lower bound
        hi=th0; lo=std::max(th0-g0/c,std::min(p.thr,th0));
      }
      //Newton from the upper bound (g convex and increasing below ths: monotone from above), bisection if it leaves the bracket
      const double tol=1e-12*std::max(p.ths,1e-3);                 //(the content to 1e-12: the balance is exact anyway)
      double x=hi;
      for (int it=0;it<100;it++){
        nit++;
        double dk; const double g=(x-th0)*c+K(p,x,&dk)-qin;
        if (g>0.0){hi=x;} else if (g<0.0){lo=x;} else {break;}
        const double dx=g/(c+dk);
        if (std::fabs(dx)<=tol){x-=dx; break;}                 //(converged: the Newton correction itself is below the tolerance)
        double xn=x-dx;
        if (!(xn>lo) || !(xn<hi)){xn=0.5*(lo+hi);}
        if (hi-lo<=tol){x=xn; break;}
        x=xn;
      }
      th=std::min(std::max(x,lo),std::min(hi,p.ths));
      return std::max(qin-(th-th0)*c,0.0);
    }
  }

  double KofKinematic(const vz_mat &m,const double th,double *dK){return K(Par(m),th,dK);}

  double StepKinematic(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info)
  {
    rejected=0.0;
    if (info){info->substeps=0; info->newton=0; info->fails=0;}
    if (!c.active){rejected=qtop*dt; return 0.0;}                               //(no column: rejected, nothing passed on)
    if (dt<=0.0){return 0.0;}
    const double q=std::max(qtop,0.0);
    if (qtop<0.0){rejected=qtop*dt;}                                            //(an upward demand: the column moves water only down)
    const size_t N=c.th.size();
    std::vector<kw_par> M(S.mat.size()); for (size_t j=0;j<M.size();j++){M[j]=Par(S.mat[j]);}
    std::vector<const kw_par*> P(N);
    bool still=(q==0.0);
    for (size_t i=0;i<N;i++){P[i]=&M[c.im[i]]; if (still && (c.th[i]>P[i]->thr)){still=false;}}
    if (still){return 0.0;}                                                     //(no inflow and nothing mobile: nothing moves)
    const double S0=Storage(c),DTH=GetAccuracy();
    std::vector<double> th0,thS=c.th;
    double t=0.0,h=dt,applied=0.0; int guard=0;
    if ((q>0.0) && (c.dz[0]*DTH<q*h)){h=std::max(c.dz[0]*DTH/q,1e-6);}            //(first substep: the top node's change within the accuracy)
    while (t<dt*(1.0-1e-12)){
      h=std::min(h,dt-t);
      th0=c.th;
      int nit=0; double qin=q,rej=0.0,dthmax=0.0;
      for (size_t i=0;i<N;i++){
        double th,exc;
        const double out=SolveNode(*P[i],th0[i],c.dz[i],h,qin,th,exc,nit);
        c.th[i]=th;
        if (exc>0.0){ //a saturated node: the excess backs up into the nodes above (their outflow was that much smaller)
          for (size_t j=i;(j-->0) && (exc>0.0);){
            const double room=std::max(P[j]->ths-c.th[j],0.0)*c.dz[j]/h,take=std::min(exc,room);
            if (take>0.0){c.th[j]=std::min(c.th[j]+take*h/c.dz[j],std::max(P[j]->ths,c.th[j])); exc-=take;}
          }
          if (exc>0.0){rej+=exc;}                                               //(the top node cannot take it: rejected)
        }
        qin=out;
      }
      for (size_t i=0;i<N;i++){dthmax=std::max(dthmax,std::fabs(c.th[i]-th0[i]));}
      if (info){info->newton+=nit;}
      if ((dthmax>DTH) && (h>1e-6) && (++guard<100000)){c.th=th0; h*=0.5; if (info){info->fails++;} continue;}
      rej=std::min(rej,q);
      applied+=(q-rej)*h; rejected+=rej*h;
      t+=h; if (info){info->substeps++;}
      if (dthmax<0.5*DTH){h*=1.5;}
    }
    for (size_t i=0;i<N;i++){if (c.th[i]!=thS[i]){c.psi[i]=(P[i]->ths>0.0)?PsiOf(S.mat[c.im[i]],c.th[i]):0.0;}} //(suction consistent)
    c.satTop=c.dryTop=false;
    //water into the water table: whatever the inflow did not leave in storage (exact by construction; >= 0 up to round-off)
    return std::max((applied-(Storage(c)-S0))/dt,0.0);
  }

  double InitSteadyKinematic(vz_column &c,const vz_strat &S,const double q)
  {
    if (!c.active){return 0.0;}
    const size_t N=c.th.size();
    std::vector<kw_par> P(N);
    double f=std::max(q,0.0); size_t ilim=N;                                   //(the deepest node of the smallest Ks below q)
    for (size_t i=0;i<N;i++){P[i]=Par(S.mat[c.im[i]]); if (P[i].Ks<=f){f=P[i].Ks; ilim=i;}}
    if (!(f<q)){ilim=N;}                                                       //(q at most every Ks: no node backs up)
    for (size_t i=0;i<N;i++){
      const kw_par &p=P[i]; double th;
      if ((ilim<N) && (i<=ilim)){th=p.ths;}                                    //(backed up above the limiting node: saturated)
      else if (p.step){th=(f>0.0)?p.ths:std::max(c.th[i],p.thr);}              //(no mobile water: ths; no inverse)
      else if (f<=0.0){th=p.thr;}
      else if (f>=p.Ks){th=p.ths;}
      else {th=std::min(p.ths*std::pow(p.x0e+f/p.A,1.0/p.e),p.ths);}           //(K(th) = f exactly)
      c.th[i]=th; c.psi[i]=(p.ths>0.0)?PsiOf(S.mat[c.im[i]],th):0.0;
    }
    c.satTop=c.dryTop=false;
    return f;
  }
}
