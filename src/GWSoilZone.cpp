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
      double s=(phi[j]-L[j].phifc)/(L[j].phimax-L[j].phifc);
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
      if (V>0.0){phi[j]=L[j].phimax;} //exactly full (no round-off left below capacity)
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
  // g(dV) = dV + c*sum_l Ceff_l*(zb_l + lambda(T(phiStar,dV)) - h_l),  c = 1000*dt/A [mm d / m3 ... per (m2/d * m)]
  // (Ceff*(head difference) is m3/d; times dt is m3; /A is m; *1000 is mm.)  g is continuous and strictly increasing.
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
  double LocalSolve(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const std::vector<sz_link> &links,
                    const double A_m2,const double dt,double &dV,double &overflow,double *dlam_dV)
  {
    std::vector<double> work;
    if (dlam_dV){*dlam_dV=0.0;}
    const double c=1000.0*dt/A_m2;
    double sC=0.0,sR=0.0;                                   //sum Ceff, sum Ceff*(zb-h)
    for (size_t l=0;l<links.size();l++){sC+=links[l].Ceff; sR+=links[l].Ceff*(links[l].zb-links[l].h);}
    overflow=0.0;
    if (sC<=0.0){dV=0.0; if (dlam_dV){*dlam_dV=LambdaSlope(L,phiStar,0.0);} return SatThickness(L,phiStar);}  //no contact: no exchange
    double lo=-StackFreeWater(L,phiStar), hi=Room(L,phiStar);
    double glo=lo+c*(sR+sC*LambdaAfter(L,phiStar,lo,work));
    double ghi=hi+c*(sR+sC*LambdaAfter(L,phiStar,hi,work));
    if (ghi<0.0){ //root beyond full: the soil zone fills and the rest overflows
      dV=hi; overflow=-ghi;
      return LambdaAfter(L,phiStar,hi,work);
    }
    if (glo>0.0){ //cannot drain more than the free water (only possible with contact below a link's base)
      dV=lo;
      return LambdaAfter(L,phiStar,lo,work);
    }
    for (int it=0;it<200;it++){ //bisection: robust on the piecewise-linear g; ~1e-15 relative in 60 halvings
      double mid=0.5*(lo+hi);
      if ((mid<=lo) || (mid>=hi)){break;}
      double g=mid+c*(sR+sC*LambdaAfter(L,phiStar,mid,work));
      if (g>0.0){hi=mid;} else {lo=mid;}
      if (hi-lo<=1e-13*std::max(1.0,std::fabs(mid))){break;}
    }
    dV=0.5*(lo+hi);
    if (dlam_dV){*dlam_dV=LambdaSlope(L,phiStar,dV);}
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
      double mean=0.0; for (size_t c=0;c<nc;c++){mean+=w[c]*cols[c][j];}
      double d=target[j]-mean;
      if (d>0.0){
        double R=0.0; for (size_t c=0;c<nc;c++){R+=w[c]*std::max(L[j].phimax-cols[c][j],0.0);}
        double f=(R>0.0)?std::min(d/R,1.0):0.0;
        for (size_t c=0;c<nc;c++){cols[c][j]+=f*std::max(L[j].phimax-cols[c][j],0.0);}
        d-=f*R;
      }
      else if (d<0.0){
        double S=0.0; for (size_t c=0;c<nc;c++){S+=w[c]*std::max(cols[c][j],0.0);}
        double f=(S>0.0)?std::min(-d/S,1.0):0.0;
        for (size_t c=0;c<nc;c++){cols[c][j]-=f*std::max(cols[c][j],0.0);}
        d+=f*S;
      }
      if (d!=0.0){for (size_t c=0;c<nc;c++){cols[c][j]+=d;}} //target outside [0,phimax] (round-off or overfill)
    }
  }
  double ConnectedBase(const double h,const std::vector<double> &zbs,const double ds)
  {
    if (zbs.empty()){return h;}
    double sw=0.0,swz=0.0,zmin=zbs[0];
    for (size_t i=0;i<zbs.size();i++){
      double w=(ds<=0.0)?((h>zbs[i])?1.0:0.0):SmoothStep((h-zbs[i])/ds);
      sw+=w; swz+=w*zbs[i]; zmin=std::min(zmin,zbs[i]);
    }
    return (sw>0.0)?swz/sw:zmin;
  }
  double SatFractionSlope(const double h,const std::vector<double> &zbs,const double ds)
  {
    if (zbs.empty() || (ds<=0.0)){return 0.0;}
    double s=0.0;
    for (size_t i=0;i<zbs.size();i++){
      double x=(h-zbs[i])/ds;
      if ((x>0.0) && (x<1.0)){s+=6.0*x*(1.0-x)/ds;}
    }
    return s/(double)zbs.size();
  }
}
