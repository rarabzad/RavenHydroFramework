/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Soil zone <-> groundwater interface: pure functions (see GWSoilZone.h)
----------------------------------------------------------------*/
#include "GWSoilZone.h"
#include <cmath>
#include <algorithm>

namespace gwsz
{
  //-- level -----------------------------------------------------------------------------------------
  int StackTop(const std::vector<sz_layer> &L,const std::vector<double> &phi)
  {
    int n=(int)L.size();
    for (int j=0;j<n;j++){if (phi[j]<L[j].phimax){return j;}}
    return n-1;
  }
  // lambda = sum_j H_j s_j prod_{i<j} s_i, with s_j = clip((phi_j-phifc_j)/(phimax_j-phifc_j),0,1) the saturated share of
  // layer j. With every layer below full this is the stack definition (level rises through a layer once the layers below
  // are full); free water perched above a partly saturated layer counts in proportion to how saturated the layers below
  // are. This keeps lambda continuous and non-decreasing in every phi_j (a hard "count only once the layer below is full"
  // jumps when the layer below fills, which leaves the local solve without a root).
  double SatThickness(const std::vector<sz_layer> &L,const std::vector<double> &phi)
  {
    double lam=0.0,conn=1.0;
    int n=(int)L.size();
    for (int j=0;(j<n) && (conn>0.0);j++){
      const double fr=L[j].phimax-L[j].phifc;
      double s=(fr>0.0)?(phi[j]-L[j].phifc)/fr:((phi[j]>=L[j].phimax)?1.0:0.0); //no free-water capacity: full or not
      s=std::min(std::max(s,0.0),1.0);
      lam+=L[j].H*s*conn;
      conn*=s;
    }
    return lam;
  }
  double StackFreeWater(const std::vector<sz_layer> &L,const std::vector<double> &phi)
  {
    int top=StackTop(L,phi); double fw=0.0;
    for (int j=top;j>=0;j--){fw+=std::max(phi[j]-L[j].phifc,0.0);}
    return fw;
  }
  double Room(const std::vector<sz_layer> &L,const std::vector<double> &phi)
  {
    double r=0.0;
    for (size_t j=0;j<L.size();j++){r+=std::max(L[j].phimax-phi[j],0.0);}
    return r;
  }

  //-- fill and drain ----------------------------------------------------------------------------------
  double Fill(const std::vector<sz_layer> &L,std::vector<double> &phi,double V)
  {
    for (size_t j=0;(j<L.size()) && (V>0.0);j++){
      double room=std::max(L[j].phimax-phi[j],0.0);
      double a=std::min(room,V);
      phi[j]+=a; V-=a;
      if ((V>0.0) && (phi[j]<L[j].phimax)){phi[j]=L[j].phimax;} //exactly full (round-off only; water above capacity is kept)
    }
    return std::max(V,0.0);
  }
  double Drain(const std::vector<sz_layer> &L,std::vector<double> &phi,double V)
  {
    if (V<=0.0){return 0.0;}
    int top=StackTop(L,phi); double removed=0.0;
    for (int j=top;(j>=0) && (V>0.0);j--){
      double fw=std::max(phi[j]-L[j].phifc,0.0);
      double r=std::min(fw,V);
      phi[j]-=r; V-=r; removed+=r;
    }
    return removed;
  }
  double Apply(const std::vector<sz_layer> &L,std::vector<double> &phi,double dV,double &shortfall)
  {
    shortfall=0.0;
    if (dV>=0.0){return Fill(L,phi,dV);}
    double want=-dV;
    double got=Drain(L,phi,want);
    double rest=want-got;
    for (size_t j=0;(j<L.size()) && (rest>0.0);j++){ //only round-off should reach here: take it below field capacity
      double r=std::min(std::max(phi[j],0.0),rest);
      phi[j]-=r; rest-=r;
    }
    shortfall=std::max(rest,0.0);
    return 0.0;
  }

  //-- local implicit solve ----------------------------------------------------------------------------
  static double LambdaAfter(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const double dV,
                            std::vector<double> &work)
  {
    work=phiStar; double sh;
    Apply(L,work,dV,sh);
    return SatThickness(L,work);
  }
  double LambdaSlope(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const double dV)
  {
    std::vector<double> work;
    const double lo=-StackFreeWater(L,phiStar),hi=Room(L,phiStar);
    const double d=1e-6*std::max(1.0,hi-lo);                 //[mm] well inside one linear piece for any real layer
    double a=dV,b=dV+d;
    if (b>hi){a=dV-d; b=dV;}
    if (a<lo){a=lo;}
    if (b-a<=0.0){return 0.0;}
    return (LambdaAfter(L,phiStar,b,work)-LambdaAfter(L,phiStar,a,work))/(b-a);
  }
  //-- sub-grid soil water level ------------------------------------------------------------------------
  // m(Lv) = sum_s clip(Lv - zbs_s, 0, D) is piecewise linear in Lv, with breakpoints at zbs_s (a sample starts to
  // fill) and zbs_s + D (it is full); its slope is the number of samples in between.
  double SubgridLevel(const sz_cell &s,const double lam,double *dLv)
  {
    const size_t N=s.zbs.size();
    if (N==0){if (dLv){*dLv=1.0;} return lam;}
    if (lam>=s.D){if (dLv){*dLv=0.0;} return s.zbs[N-1]+s.D;}
    if (lam<=0.0){
      if (dLv){size_t k=0; while ((k<N) && (s.zbs[k]<=s.zbs[0])){k++;} *dLv=(double)N/(double)k;}
      return s.zbs[0];
    }
    const double target=lam*(double)N;
    size_t i=0,j=0; int active=0; double Lc=s.zbs[0],acc=0.0;
    while ((i<N) || (j<N)){
      double nbi=(i<N)?s.zbs[i]:1e300,nbj=(j<N)?s.zbs[j]+s.D:1e300;
      double nb=std::min(nbi,nbj);
      double add=active*(nb-Lc);
      if ((active>0) && (acc+add>target)){ //strict: at a breakpoint dLv is the slope above it, as SubgridFlux counts
        if (dLv){*dLv=(double)N/(double)active;}
        return Lc+(target-acc)/(double)active;
      }
      acc+=add; Lc=nb;
      if (nbi<=nbj){active++; i++;} else {active--; j++;}
    }
    if (dLv){*dLv=0.0;}
    return s.zbs[N-1]+s.D;
  }
  double SubgridThickness(const sz_cell &s,const double Lv)
  {
    const size_t N=s.zbs.size(); if (N==0){return 0.0;}
    double m=0.0;
    for (size_t i=0;i<N;i++){m+=std::min(std::max(Lv-s.zbs[i],0.0),s.D);}
    return m/(double)N;
  }
  double SubgridFlux(const sz_cell &s,const double lam,const double h,const double hc,double &F,double &q_h,double &q_hc,
                     double &q_lam)
  {
    const size_t N=s.zbs.size();
    F=q_h=q_hc=q_lam=0.0;
    if (N==0){return 0.0;}
    double dLv; const double Lv=SubgridLevel(s,lam,&dLv);
    double sum=0.0,shc=0.0,sl=0.0,sw=0.0,sh=0.0;
    for (size_t i=0;i<N;i++){
      const double z=s.zbs[i];
      double w,wp;
      if (s.ds<=0.0){w=(hc>z)?1.0:0.0; wp=0.0;}
      else {double x=(hc-z)/s.ds; w=SmoothStep(x); wp=((x>0.0) && (x<1.0))?6.0*x*(1.0-x)/s.ds:0.0;}
      if ((w==0.0) && (wp==0.0)){continue;}
      const double l=std::min(std::max(Lv,z),z+s.D);                  //soil water level at the sample
      //the aquifer side of the sample is the water table, or its base where the water table is below it (only while hc
      //lags above h): such a sample drains just its own saturated water, never the empty soil beneath it
      const double ha=std::max(h,z);
      sum+=w*(l-ha); sw+=w; shc+=wp*(l-ha);
      if (h>z){sh+=w;}
      if ((z<=Lv) && (Lv<z+s.D)){sl+=w*dLv;}
    }
    F=sw/(double)N;
    q_h=-s.C*sh/(double)N; q_hc=s.C*shc/(double)N; q_lam=s.C*sl/(double)N;
    return s.C*sum/(double)N;
  }

  //-- local implicit solve ----------------------------------------------------------------------------
  // g(dV) = dV + c q(lambda(T(phiStar,dV)), h): q [m3/d] times dt [d] over A [m2] is metres, times 1000 mm. Water
  // into the aquifer leaves the column, so g = 0 is the column's backward-Euler balance; g is strictly increasing.
  double LocalSolve(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const sz_cell &s,const double h,
                    const double hc,const double A_m2,const double dt,double &dV,double &overflow,double &dlam_dV)
  {
    std::vector<double> work;
    double F,qh,qc,ql;
    overflow=0.0; dlam_dV=0.0;
    if (!(A_m2>0.0)){dV=0.0; return SatThickness(L,phiStar);}
    const double c=1000.0*dt/A_m2;
    SubgridFlux(s,0.0,h,hc,F,qh,qc,ql);
    if ((F<=0.0) || (s.C<=0.0)){dV=0.0; dlam_dV=LambdaSlope(L,phiStar,0.0); return SatThickness(L,phiStar);} //no contact
    double lo=-StackFreeWater(L,phiStar),hi=Room(L,phiStar);
    double llo=LambdaAfter(L,phiStar,lo,work),lhi=LambdaAfter(L,phiStar,hi,work);
    double glo=lo+c*SubgridFlux(s,llo,h,hc,F,qh,qc,ql);
    double ghi=hi+c*SubgridFlux(s,lhi,h,hc,F,qh,qc,ql);
    if (ghi<0.0){dV=hi; overflow=-ghi; return lhi;}  //root beyond full: the column fills and the rest overflows
    if (glo>0.0){dV=lo; return llo;}                  //cannot drain more than the free water (not reached: q<=0 at lambda=0,
                                                      //as every sample's level is then its base, at or below its aquifer side)
    for (int it=0;it<200;it++){                       //bisection: robust on the piecewise-smooth g; ~1e-13 relative
      double mid=0.5*(lo+hi);
      if ((mid<=lo) || (mid>=hi)){break;}
      double g=mid+c*SubgridFlux(s,LambdaAfter(L,phiStar,mid,work),h,hc,F,qh,qc,ql);
      if (g>0.0){hi=mid;} else {lo=mid;}
      if (hi-lo<=1e-13*std::max(1.0,std::fabs(mid))){break;}
    }
    dV=0.5*(lo+hi);
    dlam_dV=LambdaSlope(L,phiStar,dV);
    return LambdaAfter(L,phiStar,dV,work);
  }

  //-- saturated fraction --------------------------------------------------------------------------------
  double SmoothStep(const double x)
  {
    if (x<=0.0){return 0.0;}
    if (x>=1.0){return 1.0;}
    return x*x*(3.0-2.0*x);
  }
  double SatFraction(const double h,const std::vector<double> &zbs,const double ds)
  {
    if (zbs.empty()){return 0.0;}
    double s=0.0;
    for (size_t i=0;i<zbs.size();i++){
      if (ds<=0.0){s+=(h>zbs[i])?1.0:0.0;}
      else {s+=SmoothStep((h-zbs[i])/ds);}
    }
    return s/(double)zbs.size();
  }
  void Distribute(const std::vector<sz_layer> &L,std::vector<std::vector<double> > &cols,const std::vector<double> &w,
                  const std::vector<double> &target)
  {
    const size_t nc=cols.size();
    for (size_t j=0;j<L.size();j++){
      double mean=0.0,ws=0.0; for (size_t c=0;c<nc;c++){mean+=w[c]*cols[c][j]; ws+=w[c];}
      double d=target[j]-mean;
      if (d>0.0){
        double R=0.0; for (size_t c=0;c<nc;c++){R+=w[c]*std::max(L[j].phimax-cols[c][j],0.0);}
        double f=(R>0.0)?std::min(d/R,1.0):0.0;
        for (size_t c=0;c<nc;c++){cols[c][j]+=f*std::max(L[j].phimax-cols[c][j],0.0);}
        d-=f*R;
      }
      else if (d<0.0){
        //a loss comes first from free water (above field capacity: the water Raven's drainage acts on sits in the
        //saturated columns), then in proportion to the water left
        double S1=0.0; for (size_t c=0;c<nc;c++){S1+=w[c]*std::max(cols[c][j]-L[j].phifc,0.0);}
        double f1=(S1>0.0)?std::min(-d/S1,1.0):0.0;
        for (size_t c=0;c<nc;c++){cols[c][j]-=f1*std::max(cols[c][j]-L[j].phifc,0.0);}
        d+=f1*S1;
        if (d<0.0){
          double S=0.0; for (size_t c=0;c<nc;c++){S+=w[c]*std::max(cols[c][j],0.0);}
          double f=(S>0.0)?std::min(-d/S,1.0):0.0;
          for (size_t c=0;c<nc;c++){cols[c][j]-=f*std::max(cols[c][j],0.0);}
          d+=f*S;
        }
      }
      if ((d!=0.0) && (ws>0.0)){for (size_t c=0;c<nc;c++){cols[c][j]+=d/ws;}} //target outside [0,phimax] (round-off or overfill)
    }
  }
}
