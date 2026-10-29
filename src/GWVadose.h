/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadose.h
  One-dimensional unsaturated (vadose) column between the soil-zone base and the water table, for :SoilZone PHYSICAL
  with :SoilZoneVadose RICHARDS (design note soil_zone_columns_plan.md, section 1). Mixed-form Richards equation, Campbell
  retention with the Hutson-Cass parabola near saturation, implicit Euler, Newton with a tridiagonal solve, adaptive
  sub-steps. Lengths m, times d, contents volumetric, pressure heads m (negative unsaturated), elevations m.
  ----------------------------------------------------------------*/
#ifndef GWVADOSE_H
#define GWVADOSE_H

#include <vector>

namespace gwvz
{
  ///< unsaturated properties of one material (an aquifer class)
  struct vz_mat
  {
    double ths;    ///< saturated water content (porosity) [-]
    double b;      ///< Campbell exponent [-]
    double psia;   ///< air-entry suction (positive) [m]
    double Ks;     ///< saturated (vertical) conductivity [m/d]
    double Sy;     ///< specific yield of the aquifer class [-] (sets the reference content ths-Sy)
    double Ss;     ///< specific storage for positive pressure heads [1/m]
  };
  ///< material by elevation: layer i spans [zbot[i], ztop[i]); listed top down (MODFLOW layers)
  struct vz_strat
  {
    std::vector<double> ztop,zbot;
    std::vector<vz_mat> mat;
    int  Find(const double z) const;   ///< layer of elevation z (the last one below the bottom)
  };
  ///< the column: nodes top down
  struct vz_column
  {
    double ztop;                 ///< soil-zone base [m]
    double zbot;                 ///< water table [m]
    std::vector<double> zc,dz;   ///< node centres and thicknesses [m]
    std::vector<double> th;      ///< water contents [-]
    std::vector<double> psi;     ///< pressure heads [m] (kept consistent with th)
    std::vector<int>    im;      ///< material index of each node (into vz_strat::mat)
    bool active;                 ///< the column exists (length above the minimum)
    bool satTop;                 ///< the top was saturated at the end of the last substep (try psi = 0 there first)
    bool dryTop;                 ///< the top was air-dry under an upward demand (try its driest suction first)
    vz_column():ztop(0),zbot(0),active(false),satTop(false),dryTop(false){}
  };

  // -- retention and conductivity ------------------------------------------------------------------------
  double Theta (const vz_mat &m,const double psi,double *dth_dpsi=0);  ///< content of pressure head
  double PsiOf (const vz_mat &m,const double th);                       ///< pressure head of content (<=0 for th<=ths)
  double Kof   (const vz_mat &m,const double psi,double *dK_dpsi=0);   ///< conductivity of pressure head
  double ThEq  (const vz_mat &m,const double height);                  ///< hydrostatic content at a height above the water table
  void   SetAccuracy(const double dthmax);  ///< largest change of a content in one substep (default 0.02)
  ///< reference contents of the columns (what RelStorage, MoveBottom and Build count from): 0 the hydrostatic profile
  ///< (RICHARDS, LAG; default), 1 the residual content ths - Sy everywhere (KINEMATIC, whose no-flow state it is)
  void   SetReference(const int ref);
  int    GetReference();
  double GetAccuracy();

  // -- column management --------------------------------------------------------------------------------
  ///< builds the grid between ztop and zbot (N nodes at most, geometric from both ends), contents at the hydrostatic
  ///< profile above the water table (psi = zbot - z), returns the column's water relative to the reference [m]
  double Build   (vz_column &c,const vz_strat &S,const double ztop,const double zbot,const int Nmax);
  ///< moves the bottom to the new water table: water below a rising one is returned (absorbed, relative to the
  ///< reference contents) [m]; a falling one opens cells at the reference content ths-Sy (nothing to book).
  ///< The grid is rebuilt and the contents remapped conservatively.
  double MoveBottom(vz_column &c,const vz_strat &S,const double zbot,const int Nmax);
  double Storage  (const vz_column &c);                    ///< water in the column [m]
  double RelStorage(const vz_column &c,const vz_strat &S); ///< water above the hydrostatic (equilibrium) profile [m]

  // -- one time step -------------------------------------------------------------------------------------
  ///< advances the column over dt with the top flux qtop [m/d] (downward positive; negative = upward, taken from the
  ///< top node), the water table fixed at zbot (psi=0 there). Returns the flux into the water table qbot [m/d] such
  ///< that storage change = (qtop_applied - qbot) dt exactly; rejected [m] = qtop dt - applied: inflow the saturated top
  ///< could not take (> 0), or an upward demand (qtop < 0) the top could not supply (< 0).
  struct vz_step_info { int substeps,newton,fails; };
  double Step(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info=0,const bool noUp=false);
  ///< (noUp: no upward flow through the bottom, for a top cell MODFLOW has dried)

  ///< the same contract (storage change = (qtop_applied - qbot) dt exactly, rejected as for Step) for the cheaper methods,
  ///< on the same column (grid, contents, hydrostatic reference profile; MoveBottom and RelStorage apply unchanged). Water
  ///< only moves down (qbot >= 0; an upward demand qtop < 0 is rejected); capillary rise is the water table's (steady).
  ///<   StepKinematic: gravity-only kinematic wave (:SoilZoneVadose KINEMATIC), GWVadoseKW.cpp
  ///<   StepLag:       linear storage cascade (:SoilZoneVadose LAG; par its parameter, 0 = default), GWVadoseLag.cpp
  double StepKinematic(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info=0);
  double StepLag      (vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info=0,const double par=0.0);
  ///< LAG: travel time per metre t_m [d/m] of a material (par if > 0, else derived at q_r = 1 mm/d; GWVadoseLag.cpp), and
  ///< the cascade's steady state under the flux q [m/d] (S_i = min(q k_i, Smax_i); returns the water above the reference [m])
  double LagTravelTime(const vz_mat &m,const double par);
  double InitSteadyLag(vz_column &c,const vz_strat &S,const double q,const double par=0.0);

  ///< capillary rise into the soil zone [m/d] (upward positive): the steady upward flux (RiseBetween) from the column's
  ///< top node to a soil layer at suction psis [m] whose middle is dzs [m] above the column top, through the soil (its
  ///< saturated conductivity Ks [m/d], air entry psia [m] and Clapp-Hornberger b); 0 when the soil is wetter than that
  double CapRise(const vz_column &c,const double Ks,const double psia,const double b,const double psis,const double dzs);
  ///< exact steady upward flux [m/d] over a height dz [m] from suction psi1 >= 0 (bottom) to suction psis (top), Campbell /
  ///< Brooks-Corey conductivity K = Ks (psia/psi)^m above the air entry (m = 2 + 3/b), Ks below it: the root q of
  ///<   dz = integral_psi1^psis dpsi/(1 + q/K(psi))         (0 when psis - psi1 <= dz: hydrostatic or downward)
  ///< psis <= 0 means an infinitely dry top (the largest flux the bottom can sustain over dz)
  double RiseBetweenExact(const double Ks,const double psia,const double b,const double dz,const double psi1,const double psis);
  ///< the same from a table per Clapp-Hornberger b of the integral (built on first use, about 0.02 s each; within 2%, mostly
  ///< within 0.5%; suctions within a factor e^2 of each other are integrated directly)
  double RiseBetween(const double Ks,const double psia,const double b,const double dz,const double psi1,const double psis);
  ///< from a water table (psi1 = 0) d [m] below: RiseBetweenExact / RiseBetween
  double RiseExact(const double Ks,const double psia,const double b,const double d,const double psis);
  double Rise(const double Ks,const double psia,const double b,const double d,const double psis);
}
#endif
