/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Soil zone <-> groundwater interface (:SoilZone PHYSICAL): pure functions of one soil-zone column (the part of
  an HRU's soil zone over one grid cell), independent of Raven and MODFLOW so they can be unit tested.
  Notation as in docs/design/soil_zone_coupling_plan.md of RavenHydroFramework-extra (the design note as built):
    phi [mm] soil water of a layer, phimax = H n (1-SF), phifc = phimax S_fc, H [m] thickness,
    lambda [m] saturated thickness of the column, D [m] soil-zone thickness, F [-] connected (saturated) fraction.
  Layers are always ordered from the bottom of the soil zone (m*) upwards.

  Sub-grid closure. The cell's land surface is sampled at N DEM points; the soil-zone base there is zbs_s (land - D).
  The column's free water is a flat soil water level Lv over these samples (as MODFLOW's head is flat over the cell):
  the saturated thickness at sample s is clip(Lv - zbs_s, 0, D) and lambda is its mean. A sample is in hydraulic
  contact with the aquifer once the water table h is within the depth ds below its base (the capillary fringe), fully at
  its base: w_s = SmoothStep((h - zbs_s)/ds + 1). (Fully at the base, because MODFLOW's top cell, whose top is the
  soil-zone base, has no specific yield above it: an outlet that opened only above the base would leave a confined gap
  in which the head jumps.) The sample is saturated (its recharge rejected) as the water table rises through ds above
  its base: SatFraction, s_s = SmoothStep((h - zbs_s)/ds).
  The exchange into the aquifer is
      q = (C/N) sum_s w_s (clip(Lv, zbs_s, zbs_s + D) - max(h, zbs_s))
  (C the column's soil-aquifer conductance; max(h, zbs_s): a sample in contact while the water table is still below its
  base drains only its own saturated water, and nothing rises into it). With one sample and the water table above the
  base this is q = C (zb + lambda - h). q is continuous in h for ds > 0 (with ds = 0 it jumps where water perches).
----------------------------------------------------------------*/
#ifndef GWSOILZONE_H
#define GWSOILZONE_H

#include <vector>

namespace gwsz
{
  struct sz_layer
  {
    double phimax; ///< [mm] storage capacity H n (1-SF)
    double phifc;  ///< [mm] storage at field capacity phimax*S_fc
    double H;      ///< [m]  thickness
  };
  struct sz_cell   ///< one soil-zone column as seen by its cell
  {
    std::vector<double> zbs; ///< [m] soil-zone base at the DEM samples, ascending
    double ds;               ///< [m] connection smoothing depth (<=0: a step)
    double D;                ///< [m] soil-zone thickness
    double C;                ///< [m2/d] soil-aquifer conductance of the whole column
  };

  /// saturated thickness [m] of the soil zone, lambda = sum_j H_j s_j prod_{i<j} s_i with the saturated share
  /// s_j = clip((phi_j-phifc_j)/(phimax_j-phifc_j),0,1) (layers from the bottom): with every layer below full the level
  /// rises through a layer as its free water grows; continuous and non-decreasing in every phi_j
  double SatThickness(const std::vector<sz_layer> &L,const std::vector<double> &phi);

  /// index of the top layer of the saturated stack (first layer, from the bottom, that is not full; last if all are)
  int    StackTop(const std::vector<sz_layer> &L,const std::vector<double> &phi);

  /// adds V [mm] from the bottom upwards, each layer up to phimax; returns what does not fit (overflow) [mm]
  double Fill(const std::vector<sz_layer> &L,std::vector<double> &phi,double V);

  /// removes up to V [mm] of free water (above field capacity) from the top of the saturated stack downwards;
  /// returns the volume removed [mm] (<= V)
  double Drain(const std::vector<sz_layer> &L,std::vector<double> &phi,double V);

  /// applies a net change dV [mm] (>0 Fill, <0 Drain); returns overflow [mm] (dV>0) and sets short [mm] to any part
  /// of a removal that free water could not supply (that part is then taken below field capacity, bottom up, and
  /// short holds only what even that could not supply, normally 0)
  double Apply(const std::vector<sz_layer> &L,std::vector<double> &phi,double dV,double &shortfall);

  /// free water of the saturated stack [mm] (most that Drain can remove)
  double StackFreeWater(const std::vector<sz_layer> &L,const std::vector<double> &phi);

  /// room left in all layers [mm]
  double Room(const std::vector<sz_layer> &L,const std::vector<double> &phi);

  /// slope dlambda/dV [m/mm] of the saturated thickness after adding dV [mm] to phiStar (one-sided, on the side that
  /// stays inside [-free water, room]); lambda is piecewise linear in dV
  double LambdaSlope(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const double dV);

  /// soil water level Lv [m] of a column with mean saturated thickness lam [m]: (1/N) sum_s clip(Lv-zbs_s,0,D) = lam;
  /// dLv (if given): dLv/dlam, the slope above lam (N/ties at lam = 0, 0 at lam = D)
  double SubgridLevel(const sz_cell &s,const double lam,double *dLv=0);

  /// lam [m] that a soil water level Lv gives (the inverse of SubgridLevel)
  double SubgridThickness(const sz_cell &s,const double Lv);

  /// exchange into the aquifer [m3/d] of a column with saturated thickness lam [m]: h [m] the head that drives the
  /// flow, hc [m] the head that sets the connection weights w_s (the same at a converged interface; the iteration lets
  /// hc lag behind h to damp the switching of the connection); F [-] the contact fraction (mean w_s); partial derivatives
  /// q_h = dq/dh = -(C/N) sum over samples with h above their base of w_s, q_hc = dq/dhc, q_lam = dq/dlam [m2/d]
  double SubgridFlux(const sz_cell &s,const double lam,const double h,const double hc,double &F,double &q_h,double &q_hc,
                     double &q_lam);

  /// backward-Euler local solve of one column for fixed heads h, hc (design note section 5.6):
  ///   dV + c q(lambda(dV), h, hc) = 0,   c = 1000 dt/A [mm per m3]
  /// g is continuous and strictly increasing in dV. Returns lambda [m] at the root; dV [mm] the change of the column;
  /// overflow [mm] when the root lies beyond full; dlam_dV [m/mm] the slope of lambda there (0 when full)
  double LocalSolve(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const sz_cell &s,const double h,
                    const double hc,const double A_m2,const double dt,double &dV,double &overflow,double &dlam_dV);

  /// C1 smooth step: 0 for x<=0, 3x^2-2x^3 on (0,1), 1 for x>=1
  double SmoothStep(const double x);

  /// saturated fraction: mean over the base samples zbs [m] of SmoothStep((h-zbs)/ds) (a step at zbs for ds<=0); the
  /// share of the recharge that is rejected (the contact weights of SubgridFlux reach 1 already at the base)
  double SatFraction(const double h,const std::vector<double> &zbs,const double ds);

  /// brings the columns of one HRU (cols[c][j]: water of layer j in column c [mm], weights w[c], the column areas over
  /// the HRU's) to the HRU totals target[j] = sum_c w[c] cols[c][j]: a loss is taken from each column's free water
  /// (above phifc) first, then in proportion to its water; a gain falls on every column alike, what a full layer cannot
  /// take stays in the layer above in the same column, and what the top layer cannot take overflows. Returns the
  /// overflow [mm over the HRU]: the totals then fall short of the gains by that much, possibly moved between layers.
  /// Optionally layer jR has a known outflow R [mm over the HRU] (the recharge): it is taken by area (byArea, a
  /// constant-rate process) or from free water first, before the layer's other changes are shared, and lossR receives
  /// each column's part of it [mm over the column].
  double Distribute(const std::vector<sz_layer> &L,std::vector<std::vector<double> > &cols,const std::vector<double> &w,
                    const std::vector<double> &target,const int jR=-1,const double R=0.0,const bool byArea=false,
                    std::vector<double> *lossR=0);
}

#endif
