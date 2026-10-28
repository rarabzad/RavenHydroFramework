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
  the saturated thickness at sample s is clip(Lv - zbs_s, 0, D) and lambda is its mean. A sample is connected to the
  aquifer where the water table h is above its base, smoothly over a depth ds: w_s = SmoothStep((h - zbs_s)/ds).
  The exchange into the aquifer is
      q = (C/N) sum_s w_s (clip(Lv, zbs_s, zbs_s + D) - max(h, zbs_s))
  (C the column's soil-aquifer conductance; max(h, zbs_s) differs from h only while the iteration lets the connection
  head lag above h). With one sample this is q = w C (zb + lambda - h). q is continuous in h for ds > 0; a sample
  holding perched water (Lv above its base) adds C/N (Lv - zbs_s) over the depth ds as the water table reaches it, so
  ds must not be 0 where soil water perches (with ds = 0 that is a jump).
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
  /// hc lag behind h to damp the switching of the connection); F [-] the connected fraction; partial derivatives
  /// q_h = dq/dh (-C F at hc = h), q_hc = dq/dhc, q_lam = dq/dlam [m2/d]
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

  /// connected fraction: mean over the base samples zbs [m] of SmoothStep((h-zbs)/ds) (a step for ds<=0)
  double SatFraction(const double h,const std::vector<double> &zbs,const double ds);

  /// brings the columns of one HRU (cols[c][j]: water of layer j in column c [mm], weights w[c] summing to 1) to the
  /// HRU totals target[j] = sum_c w[c] cols[c][j] after it: a gain is shared in proportion to each column's room,
  /// a loss in proportion to each column's free water (above phifc) and, once that is used up, to its water, so no
  /// column leaves [0, phimax] that was inside it; what cannot be
  /// shared that way (only if the target itself lies outside [0, phimax]) is spread uniformly
  void Distribute(const std::vector<sz_layer> &L,std::vector<std::vector<double> > &cols,const std::vector<double> &w,
                  const std::vector<double> &target);
}

#endif
