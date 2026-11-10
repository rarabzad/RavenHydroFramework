/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadoseLag.cpp
  :SoilZoneVadose LAG {t_m} {N}: the unsaturated column as a linear storage cascade (design note soil_zone_lag.md).
  The column, its grid (MakeGrid), contents, hydrostatic reference profile (SetReference(0), as RICHARDS), MoveBottom,
  RelStorage and hotstart are those of the Richards column; only the step differs. The water a node holds above its
  equilibrium content, S_i = (th_i - th_eq,i) dz_i, is the store of a linear reservoir with residence time
  k_i = t_m dz_i (t_m [d/m] the travel time per metre of the node's material): outflow S_i/k_i (only S_i > 0) goes to
  the next node, the last node's into the water table. Mean travel time through a column of one material t_m L; the
  node count sets the spreading (effective number of reservoirs (sum k)^2 / sum k^2).
  ----------------------------------------------------------------*/
#include "GWVadose.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace gwvz
{
  static const double LAG_QREF  =0.001; //[m/d] reference flux of the derived travel time (1 mm/d)
  static const int    LAG_MAXSUB=50;    //most substeps in one step

  //////////////////////////////////////////////////////////////////
  /// \brief travel time per metre t_m [d/m] of a material: par if > 0, otherwise from the kinematic-wave celerity at the
  ///  reference flux q_r = 1 mm/d with the conductivity K = Ks Se^eps, Se = (th - (ths - Sy))/Sy, eps = 2b + 3:
  ///  c = dK/dth = eps q/(Sy Se_r) at K = q, so t_m = 1/c = Sy Se_r/(eps q), Se_r = (q/Ks)^(1/eps), q = min(q_r, Ks)
  ///  (a material with Ks < q_r cannot carry q_r: then Se_r = 1 and t_m = Sy/(eps Ks))
  //
  double LagTravelTime(const vz_mat &m,const double par)
  {
    if (par>0.0){return par;}
    const double eps=2.0*m.b+3.0;
    if (!(m.Sy>0.0) || !(eps>0.0)){return 0.0;}                                 //(no drainable water: nothing held, passed on at once)
    const double Ks=(m.Ks>0.0)?m.Ks:1e-12;
    const double q=std::min(LAG_QREF,Ks),Se=std::pow(q/Ks,1.0/eps);
    return m.Sy*Se/(eps*q);
  }

  //////////////////////////////////////////////////////////////////
  /// \brief one linear reservoir over a substep h with the constant inflow rate I >= 0 [m/d]: store S [m] (in/out),
  ///  residence time k [d], room Smax >= 0 [m] (no content above ths). Returns the outflow [m] over the substep, exactly
  ///  S0 + I h - S1. Below equilibrium (S < 0) it does not drain: the inflow first fills it up to equilibrium. Beyond its
  ///  room the excess passes on within the substep (a full node passes everything it receives).
  //
  static double Reservoir(double &S,const double I,const double h,const double k,const double Smax)
  {
    const double S0=S,V=I*h;
    double Sr=S0,hr=h;                                                          //store and time left once at or above equilibrium
    if (S0<0.0){
      if (S0+V<=0.0){S=S0+V; return 0.0;}                                       //(still below equilibrium at the end)
      hr=h*(S0+V)/V; Sr=0.0;                                                    //(time left after filling up: (S0 + V)/I)
    }
    double S1;
    if (k<=1e-12*hr){S1=I*k;}                                                   //(a negligible residence time: the store follows I k)
    else {const double em=std::expm1(-hr/k); S1=Sr*(1.0+em)-I*k*em;}            //S1 = Sr e + I k (1 - e), e = exp(-hr/k)
    if (S1>Smax){S1=std::max(Smax,0.0);}                                        //(full: the excess passes on)
    double out=S0+V-S1;
    if (out<0.0){out=0.0; S1=S0+V;}                                             //(round-off only)
    S=S1;
    return out;
  }

  //////////////////////////////////////////////////////////////////
  /// \brief advances the column over dt with the inflow qtop [m/d] at its top (contract in GWVadose.h): storage change =
  ///  (qtop_applied - qbot) dt to round-off, qbot >= 0. The column takes any inflow qtop >= 0 (a full node passes on what
  ///  it receives, the last one into the water table: the capillary fringe, where th_eq is close to ths and the nodes have
  ///  no room, is crossed at once), so rejected = 0 on an active column; an upward demand qtop < 0 is rejected entirely
  ///  (capillary rise of LAG comes from the water table, steady); an inactive column rejects qtop and returns no qbot.
  ///  Exact exponential solution per node and substep with piecewise-constant inflow (the mean outflow of the node above in
  ///  the same substep): substeps no longer than the shortest residence time, at most 50; their number is computed in
  ///  double and clamped before the cast (a tiny residence time, N = 200 or LAG 1e-6, would overflow the int)
  //
  double StepLag(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info,const double par)
  {
    rejected=0.0; if (info){info->substeps=0; info->newton=0; info->fails=0;}
    if (!c.active){rejected=qtop*dt; return 0.0;}                               //(inactive: rejected, nothing to the water table)
    if (!(dt>0.0)){return 0.0;}
    const double qin=std::max(qtop,0.0);
    if (qtop<0.0){rejected=qtop*dt;}                                            //(an upward demand: not from this column)
    const size_t N=c.zc.size();
    std::vector<double> St(N),Sm(N),k(N),eq(N);
    bool drains=false; double kmin=1e300;
    for (size_t i=0;i<N;i++){
      const vz_mat &m=S.mat[c.im[i]];
      eq[i]=ThEq(m,c.zc[i]-c.zbot);
      St[i]=(c.th[i]-eq[i])*c.dz[i]; Sm[i]=std::max(m.ths-eq[i],0.0)*c.dz[i];
      k[i]=LagTravelTime(m,par)*c.dz[i]; kmin=std::min(kmin,k[i]);
      if (St[i]>0.0){drains=true;}
    }
    if ((qin==0.0) && !drains){return 0.0;}                                     //(nothing moves: the state stays as it is)
    const double S0=Storage(c);
    const double ns=std::min(std::max(std::ceil(dt/std::max(kmin,1e-300)-1e-9),1.0),(double)LAG_MAXSUB);
    const int    n=(int)ns;
    const double h=dt/n;
    for (int s=0;s<n;s++){
      double I=qin;                                                             //[m/d] inflow of the node
      for (size_t i=0;i<N;i++){I=Reservoir(St[i],I,h,k[i],Sm[i])/h;}
    }
    for (size_t i=0;i<N;i++){
      const vz_mat &m=S.mat[c.im[i]];
      c.th[i]=eq[i]+St[i]/c.dz[i]; c.psi[i]=PsiOf(m,c.th[i]);
    }
    c.satTop=c.dryTop=false;
    if (info){info->substeps=n;}
    //water into the water table: whatever the inflow did not leave in storage (exact by construction)
    return std::max((qin*dt-(Storage(c)-S0))/dt,0.0);
  }

  //////////////////////////////////////////////////////////////////
  /// \brief steady state of the cascade under the flux q [m/d] (:SoilZoneVadoseInitialFlux): S_i = min(q k_i, Smax_i)
  ///  above the equilibrium profile of the column's current bottom; returns the column's water above the reference [m]
  //
  double InitSteadyLag(vz_column &c,const vz_strat &S,const double q,const double par)
  {
    if (!c.active){return 0.0;}
    const double qq=std::max(q,0.0);
    for (size_t i=0;i<c.zc.size();i++){
      const vz_mat &m=S.mat[c.im[i]];
      const double eq=ThEq(m,c.zc[i]-c.zbot),Sm=std::max(m.ths-eq,0.0)*c.dz[i];
      const double Si=std::min(qq*LagTravelTime(m,par)*c.dz[i],Sm);
      c.th[i]=eq+Si/c.dz[i]; c.psi[i]=PsiOf(m,c.th[i]);
    }
    c.satTop=c.dryTop=false;
    return RelStorage(c,S);
  }
}
