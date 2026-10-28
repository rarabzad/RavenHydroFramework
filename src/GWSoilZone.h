/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Soil zone <-> groundwater interface (:SoilZone PHYSICAL): pure functions of the soil-zone
  layers of one HRU, independent of Raven and MODFLOW so they can be unit tested.
  Notation as in docs/design/soil_zone_coupling_plan.md (Raven manual conventions):
    phi [mm] soil water of a layer, phimax = H n (1-SF), phifc = phimax S_fc, H [m] thickness,
    lambda [m] saturated thickness of the soil zone above its base, F [-] saturated fraction.
  Layers are always ordered from the bottom of the soil zone (m*) upwards.
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
  struct sz_link   ///< one HRU-cell link as seen by the local solve
  {
    double Ceff;   ///< [m2/d] F*C^sz
    double zb;     ///< [m]    soil-zone base under this link
    double h;      ///< [m]    head of the top active cell
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

  /// backward-Euler local solve of the soil zone of one HRU for fixed heads (plan section 6.1):
  ///   dV + c*sum_l Ceff_l*(zb_l + lambda(dV) - h_l) = 0,   c = 1000*dt/A  [mm per m3/d... see .cpp]
  /// returns lambda [m] at the root; dV [mm] the soil change; overflow [mm] when the root lies beyond full;
  /// dlam_dV (if given) [m/mm] the slope of lambda at the root (0 when the root is pinned at full or at no free water)
  double LocalSolve(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const std::vector<sz_link> &links,
                    const double A_m2,const double dt,double &dV,double &overflow,double *dlam_dV=0);

  /// slope dlambda/dV [m/mm] of the saturated thickness after adding dV [mm] to phiStar (one-sided, on the side that
  /// stays inside [-free water, room]); lambda is piecewise linear in dV
  double LambdaSlope(const std::vector<sz_layer> &L,const std::vector<double> &phiStar,const double dV);

  /// C1 smooth step: 0 for x<=0, 3x^2-2x^3 on (0,1), 1 for x>=1
  double SmoothStep(const double x);

  /// saturated fraction of a link: mean over the base samples zbs [m] of SmoothStep((h-zbs)/ds)
  double SatFraction(const double h,const std::vector<double> &zbs,const double ds);

  /// base [m] of the connected part of a link: sum_s w_s zbs_s / sum_s w_s with w_s = SmoothStep((h-zbs_s)/ds) (a step
  /// for ds<=0), the lowest base while no sample is connected; it lies below h whenever F>0, so the exchange of the
  /// connected part, F C (base + lambda - h), never drives water down across an unsaturated zone
  double ConnectedBase(const double h,const std::vector<double> &zbs,const double ds);

  /// dF/dh [1/m] of SatFraction (0 for ds<=0)
  double SatFractionSlope(const double h,const std::vector<double> &zbs,const double ds);

  /// brings the columns of one HRU (cols[c][j]: water of layer j in column c [mm], weights w[c] summing to 1) to the
  /// HRU totals target[j] = sum_c w[c] cols[c][j] after it: a gain is shared in proportion to each column's room,
  /// a loss in proportion to each column's water, so no column leaves [0, phimax] that was inside it; what cannot be
  /// shared that way (only if the target itself lies outside [0, phimax]) is spread uniformly
  void Distribute(const std::vector<sz_layer> &L,std::vector<std::vector<double> > &cols,const std::vector<double> &w,
                  const std::vector<double> &target);
}

#endif
