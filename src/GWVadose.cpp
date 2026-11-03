/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadose.cpp
  One-dimensional unsaturated column (see GWVadose.h and soil_zone_columns_plan.md section 1)
  ----------------------------------------------------------------*/
#include "GWVadose.h"
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>

static double DTHMAX=0.02;  //largest content change in one substep (time accuracy; SetAccuracy)
static int    REFMODE=0;    //reference contents: 0 hydrostatic, 1 residual ths - Sy (SetReference; honoured by Build,
                            // MoveBottom, RelStorage and ExcessIntegral through ThRef)
static const double PSIMIN=1.0e5; //[m] the driest suction considered (about 1000 MPa, beyond air-dry)
namespace gwvz
{
  //-- materials ---------------------------------------------------------------------------------------------
  int vz_strat::Find(const double z) const
  {
    for (size_t i=0;i<mat.size();i++){if (z>=zbot[i]){return (int)i;}}
    return mat.empty()?-1:(int)mat.size()-1;
  }
  // Campbell (1974) retention s = psia S^-b (s suction, S = th/ths), with the Hutson and Cass (1987) parabola for
  // S > Si = 2b/(1+2b): s = s_i sqrt((1-S)/(1-Si)), s_i = psia Si^-b, which meets the Campbell curve with the same slope
  // and reaches s = 0 at saturation; above it th = ths + Ss psi
  static inline double Si_of(const vz_mat &m){return 2.0*m.b/(1.0+2.0*m.b);}
  double Theta(const vz_mat &m,const double psi,double *dth)
  {
    if (psi>0.0){if (dth){*dth=m.Ss;} return m.ths+m.Ss*psi;}
    const double s=-psi,Si=Si_of(m),si=m.psia*std::pow(Si,-m.b);
    if (s>=si){
      const double S=std::pow(s/m.psia,-1.0/m.b);
      if (dth){*dth=m.ths*S/(m.b*s);}                    //d th/d psi = -d th/d s
      return m.ths*S;
    }
    const double x=s/si,S=1.0-(1.0-Si)*x*x;
    if (dth){*dth=m.ths*2.0*(1.0-Si)*x/si;}
    return m.ths*S;
  }
  double PsiOf(const vz_mat &m,const double th)
  {
    const double S=th/m.ths;
    if (S>=1.0){return (m.Ss>0.0)?(th-m.ths)/m.Ss:0.0;}
    const double Si=Si_of(m),si=m.psia*std::pow(Si,-m.b);
    if (S>=Si){return -si*std::sqrt((1.0-S)/(1.0-Si));}
    return std::max(-m.psia*std::pow(std::max(S,1e-12),-m.b),-PSIMIN);            //(no drier than -1e5 m: air-dry)
  }
  double Kof(const vz_mat &m,const double psi,double *dK)
  {
    double dth; const double th=Theta(m,psi,&dth);
    const double S=std::min(th/m.ths,1.0),e=2.0*m.b+3.0;
    const double K=m.Ks*std::pow(S,e);
    if (dK){*dK=(th>=m.ths)?0.0:m.Ks*e*std::pow(S,e-1.0)*dth/m.ths;}
    return K;
  }
  double ThEq(const vz_mat &m,const double height){return Theta(m,-std::max(height,0.0));}
  //reference content of a node height [m] above the water table: the hydrostatic profile (REFMODE 0: RICHARDS, LAG) or
  //the residual content ths - Sy (REFMODE 1: KINEMATIC, whose no-flow state it is; ths for Sy = 0)
  double ThRef(const vz_mat &m,const double height){return (REFMODE==1)?std::max(m.ths-std::max(m.Sy,0.0),0.0):ThEq(m,height);}
  double SyRise(const vz_strat &S,const double z0,const double zt)
  {
    double v=0.0; if (!(zt>z0)){return 0.0;}
    for (size_t i=0;i<S.mat.size();i++){
      const double lo=std::max(z0,S.zbot[i]),hi=std::min(zt,S.ztop[i]);
      if (hi>lo){v+=(hi-lo)*S.mat[i].Sy;}
    }
    return v;
  }
  void   SetAccuracy(const double dthmax){if (dthmax>0.0){DTHMAX=dthmax;}}
  void   SetReference(const int ref){REFMODE=ref;}
  int    GetReference(){return REFMODE;}
  double GetAccuracy(){return DTHMAX;}

  //-- grid --------------------------------------------------------------------------------------------------
  static void MakeGrid(vz_column &c,const vz_strat &S,const int Nmax)
  {
    const double L=c.ztop-c.zbot;
    int N=(int)std::ceil(L/0.05); N=std::max(3,std::min(N,std::max(Nmax,3)));
    std::vector<double> w(N);
    //growth 1.3 from both ends, capped at 1.3^8: from N = 19 on the extra nodes fill the middle at
    //equal thickness, so raising N refines the whole column; N <= 18 never reaches the cap (the grid is unchanged)
    const double r=1.3; int h=std::min((N+1)/2,8);
    for (int k=0;k<N;k++){int j=std::min(k,N-1-k); w[k]=std::pow(r,(double)std::min(j,h));}
    double sw=0.0; for (int k=0;k<N;k++){sw+=w[k];}
    c.zc.assign(N,0.0); c.dz.assign(N,0.0); c.im.assign(N,0);
    double z=c.ztop;
    for (int k=0;k<N;k++){c.dz[k]=L*w[k]/sw; c.zc[k]=z-0.5*c.dz[k]; z-=c.dz[k]; c.im[k]=std::max(S.Find(c.zc[k]),0);}
  }
  double Storage(const vz_column &c){double s=0.0; for (size_t i=0;i<c.th.size();i++){s+=c.th[i]*c.dz[i];} return s;}
  //storage relative to the reference profile (the hydrostatic profile above the water table, or ths - Sy with
  //SetReference(1)): the water table's own fluctuations are MODFLOW's (its specific yield); the column carries what is
  //in excess of (or missing from) the reference: recharge in transit, and deficits left by evaporation and capillary rise
  double RelStorage(const vz_column &c,const vz_strat &S)
  {
    double s=0.0; for (size_t i=0;i<c.th.size();i++){s+=(c.th[i]-ThRef(S.mat[c.im[i]],c.zc[i]-c.zbot))*c.dz[i];} return s;
  }
  double Build(vz_column &c,const vz_strat &S,const double ztop,const double zbot,const int Nmax)
  {
    c.ztop=ztop; c.zbot=zbot; c.active=(ztop>zbot);
    if (!c.active){c.zc.clear(); c.dz.clear(); c.th.clear(); c.psi.clear(); c.im.clear(); return 0.0;}
    MakeGrid(c,S,Nmax);
    const size_t N=c.zc.size(); c.th.assign(N,0.0); c.psi.assign(N,0.0);
    if (REFMODE==1){for (size_t i=0;i<N;i++){const vz_mat &m=S.mat[c.im[i]]; c.th[i]=ThRef(m,0.0); c.psi[i]=(m.ths>0.0)?PsiOf(m,c.th[i]):0.0;}}
    else {for (size_t i=0;i<N;i++){c.psi[i]=zbot-c.zc[i]; c.th[i]=Theta(S.mat[c.im[i]],c.psi[i]);}}
    return RelStorage(c,S);
  }
  // integral of the excess over the reference (piecewise constant per node) over [a,b] (a<b); zero outside the column
  // (below its bottom the ground was saturated: no excess)
  static double ExcessIntegral(const vz_column &c,const vz_strat &S,double a,double b)
  {
    double v=0.0,zt=c.ztop;
    for (size_t i=0;i<c.zc.size();i++){
      const double zb=zt-c.dz[i],lo=std::max(a,zb),hi=std::min(b,zt);
      if (hi>lo){v+=(hi-lo)*(c.th[i]-ThRef(S.mat[c.im[i]],c.zc[i]-c.zbot));}
      zt=zb;
    }
    return v;
  }
  double MoveBottom(vz_column &c,const vz_strat &S,const double zbot,const int Nmax)
  {
    if (!c.active){return 0.0;}
    if (zbot>=c.ztop){ //the column is gone: its excess over equilibrium joins the aquifer
      const double ab=RelStorage(c,S);
      c.active=false; c.satTop=c.dryTop=false; c.zbot=zbot; c.zc.clear(); c.dz.clear(); c.th.clear(); c.psi.clear(); c.im.clear();
      return ab;
    }
    if (std::fabs(zbot-c.zbot)<1e-9){return 0.0;}                               //(not moved: the state stays as it is)
    double absorbed=(zbot>c.zbot)?ExcessIntegral(c,S,c.zbot,zbot):0.0;          //a rising water table takes the excess below it
    vz_column o=c;
    c.zbot=zbot; MakeGrid(c,S,Nmax);
    const size_t N=c.zc.size(); c.th.assign(N,0.0); c.psi.assign(N,0.0);
    double zt=c.ztop;
    for (size_t i=0;i<N;i++){ //the excess keeps its elevation (conservative remap); the reference follows the new water table
      const double zb=zt-c.dz[i]; const vz_mat &m=S.mat[c.im[i]];               //(ths - Sy with SetReference(1): new nodes open at it
      const double eq=ThRef(m,c.zc[i]-zbot),ex=ExcessIntegral(o,S,zb,zt)/c.dz[i]; // and nothing is exchanged, no air-dry clip)
      double th=eq+ex;
      const double thmax=m.ths,thmin=(REFMODE==1)?0.0:Theta(m,-PSIMIN);           //(a node straddling materials may average beyond
      if (th>thmax){absorbed+=(th-thmax)*c.dz[i]; th=thmax;}                      // its porosity: the rest drains to the water table;
      if (th<thmin){absorbed-=(thmin-th)*c.dz[i]; th=thmin;}                      // nor drier than air-dry: taken from it)
      c.th[i]=th; c.psi[i]=PsiOf(m,th);
      zt=zb;
    }
    c.satTop=c.dryTop=false;
    return absorbed;
  }

  //-- one step ----------------------------------------------------------------------------------------------
  static void Thomas(std::vector<double> &a,std::vector<double> &b,std::vector<double> &cc,std::vector<double> &d,std::vector<double> &x)
  {
    const size_t n=b.size();
    for (size_t i=1;i<n;i++){double m=a[i]/b[i-1]; b[i]-=m*cc[i-1]; d[i]-=m*d[i-1];}
    x[n-1]=d[n-1]/b[n-1];
    for (size_t i=n-1;i-->0;){x[i]=(d[i]-cc[i]*x[i+1])/b[i];}
  }
  // one implicit substep. top: 0 flux qtop, 1 psi = 0 (a saturated top that cannot take the inflow), 2 psi = -PSIMIN
  // (an air-dry top that cannot give the upward demand); noUp: no upward flow through the bottom (a top cell MODFLOW
  // has dried: nothing to draw from). Returns false if Newton did not converge; qin is the top inflow [m/d]
  static bool Substep(vz_column &c,const vz_strat &S,const double dt,const double qtop,const int top,const bool noUp,
                      const std::vector<double> &th0,double &qin,int &nit)
  {
    const size_t N=c.zc.size();
    std::vector<double> th(N),C(N),K(N),dK(N),R(N),a(N,0.0),b(N,0.0),cc(N,0.0),x(N);
    std::vector<double> psi=c.psi;
    for (int it=0;it<20;it++){
      nit++;
      for (size_t i=0;i<N;i++){
        const vz_mat &m=S.mat[c.im[i]]; th[i]=Theta(m,psi[i],&C[i]); K[i]=Kof(m,psi[i],&dK[i]);
        if (psi[i]>-m.psia){C[i]=std::max(C[i],std::max(m.Ss,1e-7));}  //(near saturation the Jacobian's capacity is floored at Ss:
      }                                                                     // at psi = 0 the retention slope is 0, and Newton would cycle)
      for (size_t i=0;i<N;i++){R[i]=c.dz[i]*(th[i]-th0[i])/dt; a[i]=0.0; b[i]=c.dz[i]*C[i]/dt; cc[i]=0.0;}
      //interior faces (downward positive)
      for (size_t i=0;i+1<N;i++){
        const double Dz=c.zc[i]-c.zc[i+1],Kf=0.5*(K[i]+K[i+1]),G=(psi[i]-psi[i+1])/Dz+1.0,q=Kf*G;
        const double dqi=0.5*dK[i]*G+Kf/Dz,dqj=0.5*dK[i+1]*G-Kf/Dz;
        R[i]+=q;   b[i]+=dqi;    cc[i]+=dqj;      //leaves node i
        R[i+1]-=q; a[i+1]-=dqi;  b[i+1]-=dqj;     //enters node i+1
      }
      //bottom face: the water table (psi = 0) half a node below
      {
        const size_t i=N-1; const vz_mat &m=S.mat[c.im[i]];
        const double Dz=c.zc[i]-c.zbot,Kf=0.5*(K[i]+m.Ks),G=psi[i]/Dz+1.0,q=Kf*G;
        if (!(noUp && (q<0.0))){R[i]+=q; b[i]+=0.5*dK[i]*G+Kf/Dz;}
      }
      //top face
      if (top>0){
        const vz_mat &m=S.mat[c.im[0]]; const double pt=(top==1)?0.0:-PSIMIN;
        const double Kt=(top==1)?m.Ks:Kof(m,pt),Dz=0.5*c.dz[0],Kf=0.5*(K[0]+Kt),G=(pt-psi[0])/Dz+1.0;
        qin=Kf*G; R[0]-=qin; b[0]-=0.5*dK[0]*G-Kf/Dz;
      }
      else {qin=qtop; R[0]-=qtop;}
      double rmax=0.0; for (size_t i=0;i<N;i++){rmax=std::max(rmax,std::fabs(R[i])*dt);}
      if ((it>0) && (rmax<1e-11)){
        c.psi=psi; for (size_t i=0;i<N;i++){c.th[i]=th[i];}
        return true;
      }
      for (size_t i=0;i<N;i++){R[i]=-R[i];}
      Thomas(a,b,cc,R,x);
      //damping: a node may move towards saturation by most of its suction (dry nodes recover in a few iterations, in
      //log steps), away from it by 30% of it (at least 0.05 m)
      double scale=1.0;
      for (size_t i=0;i<N;i++){
        const bool wetter=(x[i]>0.0) && (psi[i]<0.0);
        const double lim=wetter?std::max(0.5,0.9*std::fabs(psi[i])):std::max(0.05,0.3*std::fabs(psi[i]));
        if (std::fabs(x[i])*scale>lim){scale=lim/std::fabs(x[i]);}
      }
      for (size_t i=0;i<N;i++){psi[i]=std::max(psi[i]+scale*x[i],-PSIMIN); if (!std::isfinite(psi[i])){return false;}}
    }
    return false;
  }
  double Step(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info,const bool noUp)
  {
    rejected=0.0;
    if (info){info->substeps=0; info->newton=0; info->fails=0;}
    if (!c.active){rejected=qtop*dt; return 0.0;}                               //(no column: all of qtop rejected, nothing passed on)
    if (dt<=0.0){return 0.0;}
    const size_t N=c.zc.size();
    //fast path: no inflow and the profile hydrostatic (nothing moves)
    if (qtop==0.0){
      bool eq=true; for (size_t i=0;(i<N) && eq;i++){if (std::fabs(c.psi[i]-(c.zbot-c.zc[i]))>1e-9){eq=false;}}
      if (eq){return 0.0;}
    }
    if (qtop<=0.0){c.satTop=false;}
    if (qtop>=0.0){c.dryTop=false;}
    const double S0=Storage(c);
    double t=0.0,h=dt,inflow=0.0; int guard=0;
    std::vector<double> th0;
    while (t<dt*(1.0-1e-12)){
      if (++guard>100000){rejected+=qtop*(dt-t); t=dt; break;}                 //(the rest of the inflow is rejected, not lost)
      h=std::min(h,dt-t);
      th0=c.th; std::vector<double> psi0=c.psi;
      double qin=0.0; int nit=0; bool ok=false; int top=0;
      //a top that was saturated (inflow) or air-dry (upward demand) is tried with its fixed head first, and keeps it
      //while the column cannot take (or give) more than that
      if ((qtop>0.0) && c.satTop){
        double qd=0.0; ok=Substep(c,S,h,qtop,1,noUp,th0,qd,nit);
        if (ok && (qd<qtop)){qin=qd; top=1;} else {c.th=th0; c.psi=psi0; ok=false;}
      }
      else if ((qtop<0.0) && c.dryTop){
        double qd=0.0; ok=Substep(c,S,h,qtop,2,noUp,th0,qd,nit);
        if (ok && (qd>qtop)){qin=qd; top=2;} else {c.th=th0; c.psi=psi0; ok=false;}
      }
      if (top==0){
        ok=Substep(c,S,h,qtop,0,noUp,th0,qin,nit);
        if (ok && (qtop>0.0) && (c.psi[0]>0.0)){ //the top saturates: it takes only what can flow on (psi = 0 at the top); a
                                                 //negative qd is water the column pushes back out of its top
          std::vector<double> thF=c.th,psiF=c.psi;
          c.th=th0; c.psi=psi0;
          double qd=0.0; bool okD=Substep(c,S,h,qtop,1,noUp,th0,qd,nit);
          if (!okD){ok=false;}                                   //(a failed substep: halved below)
          else if (qd<qtop){qin=qd; top=1;}
          else {c.th=thF; c.psi=psiF;}                           //it can take it all: the flux solution stands
        }
        else if (!ok && (qtop<0.0)){ //an upward demand the top cannot meet: psi = -PSIMIN at the top
          c.th=th0; c.psi=psi0;
          double qd=0.0; bool okD=Substep(c,S,h,qtop,2,noUp,th0,qd,nit);
          if (okD && (qd>qtop)){ok=true; qin=qd; top=2;} else {c.th=th0; c.psi=psi0;}
        }
      }
      double dthmax=0.0; if (ok){for (size_t i=0;i<N;i++){dthmax=std::max(dthmax,std::fabs(c.th[i]-th0[i]));}}
      if (info){info->newton+=nit;}
      if (!ok || ((dthmax>DTHMAX) && (h>1e-6))){ //retry with half the substep
        c.th=th0; c.psi=psi0; h*=0.5; if (info){info->fails++;}
        if (h<1e-10*dt){ //give up refining: keep the state; the inflow of the rest of the step is rejected (an upward
                         //demand that cannot be met: rejected < 0)
          rejected+=qtop*(dt-t); t=dt; break;
        }
        continue;
      }
      c.satTop=(top==1); c.dryTop=(top==2);
      inflow+=qin*h; rejected+=(qtop-qin)*h;
      t+=h; if (info){info->substeps++;}
      if ((nit<=4) && (dthmax<0.5*DTHMAX)){h*=1.5;}
    }
    //water into the water table: whatever the inflow did not leave in storage (exact by construction)
    return (inflow-(Storage(c)-S0))/dt;
  }

  //-- capillary rise ----------------------------------------------------------------------------------------
  double CapRise(const vz_column &c,const double Ks,const double psia,const double b,const double psis,const double dzs)
  {
    if (!c.active || c.zc.empty()){return 0.0;}
    const double dz=dzs+(c.ztop-c.zc[0]);                 //from the column's top node to the soil layer's middle
    return RiseBetween(Ks,psia,b,dz,std::max(-c.psi[0],0.0),psis);
  }
  //---------------------------------------------------------------------------------------------------------------------
  // exact steady upward flux. Dimensionless: u = psi/psia, Q = q/Ks, K/Ks = 1 (u <= 1), u^-m (u > 1), m = 2 + 3/b. Steady
  // upward flux q over a height dz between suctions psi1 (bottom) and psis (top): dpsi/dz = 1 + q/K, so
  //   dz/psia = G(Q,us) - G(Q,u1),  G(Q,u) = min(u,1)/(1+Q) + I(Q,u),  I(Q,u) = int_1^max(u,1) du'/(1+Q u'^m)
  // G(Q,us) - G(Q,u1) decreases from us - u1 (Q = 0) to 0 (Q -> infinity): one root when us - u1 > dz/psia.
  static double RiseInt(const double Q,const double s,const double m)                 //I(Q,s) for s > 1; s <= 0: infinity
  {
    static const double gx[10]={-0.9739065285171717,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,
                                 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    static const double gw[10]={0.0666713443086881,0.1494513491001158,0.2190863625370051,0.2692667193099963,0.2955242247147529,
                                0.2955242247147529,0.2692667193099963,0.2190863625370051,0.1494513491001158,0.0666713443086881};
    const double U=std::pow(1e7/Q,1.0/m);                                             //beyond U: 1/(Q u^m) to 1e-7, integrated exactly
    const bool inf=(s<=0.0);
    const double up=inf?std::max(U,1.0):std::min(s,std::max(U,1.0));
    double I=0.0;
    if (up>1.0){
      const double V=std::log(up); const int NS=12; const double hV=V/NS;
      for (int g=0;g<NS;g++){
        const double c0=(g+0.5)*hV;
        for (int i=0;i<10;i++){const double v=c0+0.5*hV*gx[i],u=std::exp(v); I+=0.5*hV*gw[i]*u/(1.0+Q*std::pow(u,m));}
      }
    }
    if (inf || (s>up)){
      const double a=std::pow(up,1.0-m),bb=inf?0.0:std::pow(s,1.0-m);
      I+=(a-bb)/(Q*(m-1.0));
    }
    return I;
  }
  static double GExact(const double Q,const double u,const double m)                   //u < 0: infinity
  {
    if ((u>=0.0) && (u<=1.0)){return u/(1.0+Q);}
    return 1.0/(1.0+Q)+RiseInt(Q,(u<0.0)?-1.0:u,m);
  }
  //G(Q,u2) - G(Q,u1) directly, for 1 <= u1 < u2 <= e^2 u1 (where the difference of two table values would lose accuracy)
  static double GDiffNear(const double Q,const double u1,const double u2,const double m)
  {
    static const double gx[10]={-0.9739065285171717,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,
                                 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    static const double gw[10]={0.0666713443086881,0.1494513491001158,0.2190863625370051,0.2692667193099963,0.2955242247147529,
                                0.2955242247147529,0.2692667193099963,0.2190863625370051,0.1494513491001158,0.0666713443086881};
    const double a=std::log(u1),h=std::log(u2)-a; double I=0.0;
    for (int i=0;i<10;i++){const double v=a+0.5*h*(1.0+gx[i]),u=std::exp(v); I+=0.5*h*gw[i]*u/(1.0+Q*std::pow(u,m));}
    return I;
  }
  //root in ln Q of f (decreasing) on [a,b] (Brent); f(a) > 0 > f(b) checked by the caller
  template<class F> static double RootLnQ(F f,double a,double b,double fa,double fb,const double rtol=1e-11)
  {
    double c=a,fc=fa,dd=b-a,e=dd;
    for (int it=0;it<100;it++){
      if (fb*fc>0.0){c=a; fc=fa; dd=e=b-a;}
      if (std::fabs(fc)<std::fabs(fb)){a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;}
      const double tol=rtol*std::max(1.0,std::fabs(b)),xm=0.5*(c-b);
      if ((std::fabs(xm)<=tol) || (fb==0.0)){break;}
      if ((std::fabs(e)>=tol) && (std::fabs(fa)>std::fabs(fb))){
        double p,q,r,ss=fb/fa;
        if (a==c){p=2.0*xm*ss; q=1.0-ss;}
        else {q=fa/fc; r=fb/fc; p=ss*(2.0*xm*q*(q-r)-(b-a)*(r-1.0)); q=(q-1.0)*(r-1.0)*(ss-1.0);}
        if (p>0.0){q=-q;} p=std::fabs(p);
        if (2.0*p<std::min(3.0*xm*q-std::fabs(tol*q),std::fabs(e*q))){e=dd; dd=p/q;} else {dd=xm; e=dd;}
      }
      else {dd=xm; e=dd;}
      a=b; fa=fb;
      b+=(std::fabs(dd)>tol)?dd:((xm>0.0)?tol:-tol);
      fb=f(b);
    }
    return b;
  }
  static const double LQ0=-46.0,LQ1=14.0;                                               //range of ln Q
  template<class G> static double SolveRise(G g,const double x,const double u1,const double us,const double m,const double rtol=1e-11)
  {
    if ((us>=0.0) && (us-u1<=x)){return 0.0;}                                          //hydrostatic or downward
    const bool near=(us>0.0) && (us<=std::exp(2.0)*std::max(u1,1.0));               //(the integral directly: no cancellation)
    auto f=[&](double lq){
      const double Q=std::exp(lq);
      if (!near){return g(Q,us)-g(Q,u1)-x;}
      const double sat=(u1<1.0)?(std::min(us,1.0)-u1)/(1.0+Q):0.0;                    //(the part below the air entry)
      return sat+((us>1.0)?GDiffNear(Q,std::max(u1,1.0),us,m):0.0)-x;
    };
    //bracket: from Gardner's estimate for a dry top (Q = (1 + 1.5/(m-1)) x^-m beyond the air entry, 1/x - 1 within it),
    //widened by 3 in ln Q until the sign changes (f falls with Q), within [LQ0, LQ1]
    const double xe=std::max(x+u1,1e-12);
    const double Qg=(xe>1.0)?(1.0+1.5/(m-1.0))*std::pow(xe,-m):std::max(1.0/xe-1.0,1e-3);
    double a=std::min(std::max(std::log(Qg),LQ0),LQ1),fa=f(a),b,fb;
    if (fa>0.0){b=a; fb=fa; do {a=b; fa=fb; b=std::min(b+3.0,LQ1); fb=f(b);} while ((fb>0.0) && (b<LQ1));
                if (fb>=0.0){return std::exp(LQ1);}}
    else       {b=a; fb=fa; do {b=a; fb=fa; a=std::max(a-3.0,LQ0); fa=f(a);} while ((fa<0.0) && (a>LQ0));
                if (fa<=0.0){return std::exp(LQ0);}}
    if (fb==0.0){return std::exp(b);}
    return std::exp(RootLnQ(f,a,b,fa,fb,rtol));
  }
  double RiseBetweenExact(const double Ks,const double psia,const double b,const double dz,const double psi1,const double psis)
  {
    if (!(Ks>0.0) || !(psia>0.0) || !(b>0.0) || !(dz>0.0)){return 0.0;}
    const double m=2.0+3.0/b;
    auto g=[&](double Q,double u){return GExact(Q,u,m);};
    return Ks*SolveRise(g,dz/psia,std::max(psi1,0.0)/psia,(psis>0.0)?psis/psia:-1.0,m);
  }
  double RiseExact(const double Ks,const double psia,const double b,const double d,const double psis)
  {
    return RiseBetweenExact(Ks,psia,b,d,0.0,psis);
  }
  //table per m of ln J(ln Q, t) = ln( I(Q,u)/(u-1) ), u = e^t, t on [0, TT] (J = 1/(1+Q) at t = 0, decreasing in t), and
  //ln I(Q,infinity); bilinear in (ln Q, t), linear in ln Q
  struct rise_table { std::vector<double> lJ,lIinf; double m; };
  static const int    NQ=301,NT=185;
  static const double TT=18.4;
  static const rise_table &RiseTable(const double b)
  {
    static std::map<long,rise_table> cache;
    const long key=(long)std::floor(b*1e4+0.5);
    std::map<long,rise_table>::iterator it=cache.find(key);
    if (it!=cache.end()){return it->second;}
    static const double gx[10]={-0.9739065285171717,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,
                                 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    static const double gw[10]={0.0666713443086881,0.1494513491001158,0.2190863625370051,0.2692667193099963,0.2955242247147529,
                                0.2955242247147529,0.2692667193099963,0.2190863625370051,0.1494513491001158,0.0666713443086881};
    rise_table T; T.m=2.0+3.0/(key*1e-4); T.lJ.resize(NQ*NT); T.lIinf.resize(NQ);
    const double m=T.m,ht=TT/(NT-1);
    for (int i=0;i<NQ;i++){
      const double Q=std::exp(LQ0+(LQ1-LQ0)*i/(NQ-1));
      double I=0.0;
      T.lJ[i*NT]=-std::log(1.0+Q);
      for (int j=1;j<NT;j++){                                                          //I accumulated cell by cell in t
        const double c0=(j-0.5)*ht;
        for (int g=0;g<10;g++){const double v=c0+0.5*ht*gx[g],u=std::exp(v); I+=0.5*ht*gw[g]*u/(1.0+Q*std::pow(u,m));}
        T.lJ[i*NT+j]=std::log(I/(std::exp(j*ht)-1.0));
      }
      T.lIinf[i]=std::log(RiseInt(Q,-1.0,m));
    }
    return cache[key]=T;
  }
  static double GTab(const rise_table &T,const double Q,const double u)            //u < 0: infinity
  {
    if ((u>=0.0) && (u<=1.0)){return u/(1.0+Q);}
    const double fi=std::min(std::max((std::log(Q)-LQ0)/(LQ1-LQ0)*(NQ-1),0.0),(double)(NQ-1));
    const int i=std::min((int)fi,NQ-2); const double wi=fi-i;
    if (u<0.0){return 1.0/(1.0+Q)+std::exp((1.0-wi)*T.lIinf[i]+wi*T.lIinf[i+1]);}
    const double t=std::log(u);
    if (t>=TT){                                                                        //(beyond the table: its end plus the tail 1/(Q u^m))
      const double ue=std::exp(TT),Ie=(ue-1.0)*std::exp((1.0-wi)*T.lJ[i*NT+NT-1]+wi*T.lJ[(i+1)*NT+NT-1]);
      return 1.0/(1.0+Q)+Ie+(std::pow(ue,1.0-T.m)-std::pow(u,1.0-T.m))/(Q*(T.m-1.0));
    }
    const double fj=t/TT*(NT-1); const int j=std::min((int)fj,NT-2); const double wj=fj-j;
    const double lJ=(1.0-wi)*((1.0-wj)*T.lJ[i*NT+j]+wj*T.lJ[i*NT+j+1])+wi*((1.0-wj)*T.lJ[(i+1)*NT+j]+wj*T.lJ[(i+1)*NT+j+1]);
    return 1.0/(1.0+Q)+(u-1.0)*std::exp(lJ);
  }
  double RiseBetween(const double Ks,const double psia,const double b,const double dz,const double psi1,const double psis)
  {
    if (!(Ks>0.0) || !(psia>0.0) || !(b>0.0) || !(dz>0.0)){return 0.0;}
    const rise_table &T=RiseTable(b);
    auto g=[&](double Q,double u){return GTab(T,Q,u);};
    return Ks*SolveRise(g,dz/psia,std::max(psi1,0.0)/psia,(psis>0.0)?psis/psia:-1.0,T.m,1e-6); //(ln Q to 1e-6: the table is within 2 %)
  }
  double Rise(const double Ks,const double psia,const double b,const double d,const double psis)
  {
    return RiseBetween(Ks,psia,b,d,0.0,psis);
  }
  //////////////////////////////////////////////////////////////////
  /// steady RICHARDS profile of column c for the downward flux q [m/d] (:SoilZoneVadoseInitialFlux), in the finite
  ///  volumes of gwvz::Step: every face carries q, K_f (dpsi/dz + 1) = q with K_f the mean of the two nodes' conductivities
  ///  (the water table's face: the bottom node and Ks, psi = 0 at the water table). Solved node by node from the bottom up
  ///  (each face's flux increases with the upper node's pressure head where the gradient is downward): a node above a layer
  ///  whose Ks is less than q builds up positive pressure (a saturated, perched zone), as the steady state requires. When q
  ///  exceeds what the column can pass with psi = 0 at its top (the profile for q then has positive pressure at the top
  ///  node and a top face that cannot carry q: Step would saturate the top and reject part of q), the column starts at
  ///  the backed-up steady state instead: the largest flux f < q whose profile's top face, at psi = 0 (Step's saturated
  ///  top), carries f (bisection on f), with the top saturated; the top then rejects q - f, as InitSteadyKinematic does
  ///. The discrete start is steady to the solver's tolerance in both cases. Returns the column's water
  ///  above the reference [m]; *fw (optional) the flux into the water table [m/d] (q, or f when backed up)
  //
  double InitSteadyRichards(vz_column &c,const vz_strat &S,const double q,double *fw)
  {
    if (fw){*fw=0.0;}
    if (!c.active || !(q>0.0)){return 0.0;}
    const int N=(int)c.zc.size(); const double PMIN=-1e5;
    //root of f(p) = Kf(p) ((p - pl)/Dz + 1) - qq for the node's pressure head p over the one below (pl, conductivity Kl)
    auto solve=[&](const vz_mat &m,const double pl,const double Kl,const double Dz,const double qq)->double{
      auto f=[&](const double p)->double{return 0.5*(Kof(m,p)+Kl)*((p-pl)/Dz+1.0)-qq;};
      double a=std::max(pl-Dz,PMIN),b=std::max(a,0.0)+1.0;                                   //(f(a) <= 0: no downward gradient)
      if (f(a)>=0.0){return a;}
      for (int it=0;(it<200) && (f(b)<0.0);it++){b=a+2.0*(b-a);}
      for (int it=0;it<200;it++){                                                  //(bisection to round-off: f increases on [a,b])
        const double mid=0.5*(a+b); if ((mid<=a) || (mid>=b)){break;}
        if (f(mid)<0.0){a=mid;} else {b=mid;}
      }
      return 0.5*(a+b);
    };
    auto profile=[&](const double qq){
      double pl=0.0,Kl=S.mat[c.im[N-1]].Ks,zl=c.zbot;                              //(the water table: psi = 0, its Ks)
      for (int i=N-1;i>=0;i--){
        const vz_mat &m=S.mat[c.im[i]];
        const double p=solve(m,pl,Kl,c.zc[i]-zl,qq);
        c.psi[i]=p; c.th[i]=Theta(m,p);
        pl=p; Kl=Kof(m,p); zl=c.zc[i];
      }
    };
    auto topFlux=[&]()->double{                                                    //(Step's top face with psi = 0 at the column's top)
      const vz_mat &m=S.mat[c.im[0]];
      return 0.5*(Kof(m,c.psi[0])+m.Ks)*((0.0-c.psi[0])/(0.5*c.dz[0])+1.0);
    };
    profile(q);
    double f=q; bool backed=false;
    if ((c.psi[0]>0.0) && (topFlux()<q)){                                          //(beyond the column's capacity: backed up)
      double lo=0.0,hi=q;                                                          //(topFlux - f: > 0 at lo, < 0 at hi)
      for (int it=0;it<200;it++){
        const double mid=0.5*(lo+hi); if ((mid<=lo) || (mid>=hi)){break;}
        profile(mid);
        if (topFlux()>mid){lo=mid;} else {hi=mid;}
        if (hi-lo<=1e-15*q){break;}
      }
      f=0.5*(lo+hi); profile(f); backed=true;
    }
    c.satTop=backed || (c.psi[0]>=-1e-9); c.dryTop=false;
    if (fw){*fw=f;}
    return RelStorage(c,S);
  }
}
