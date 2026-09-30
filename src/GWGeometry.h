/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: spatial utilities
  (GeoJSON reading, equal-area projection, ESRI ASCII rasters,
   polygon/cell intersection)
----------------------------------------------------------------*/
#ifndef GWGEOMETRY_H
#define GWGEOMETRY_H

#include <string>
#include <vector>
#include <map>
using namespace std;

struct gw_pt { double x,y; };
typedef vector<gw_pt> gw_ring;   ///< closed or open ring; outer rings CCW (+), holes CW (-)

///////////////////////////////////////////////////////////////////
/// \brief Spherical Lambert azimuthal equal-area projection centred on the model domain
/// \details Equal-area, so overlap areas (and therefore volumes) are preserved
//
class CGWProjection
{
private:
  double _lat0,_lon0,_sinp0,_cosp0;
public:
  CGWProjection();
  void SetOrigin(const double lat0,const double lon0);
  void Forward(const double lon,const double lat,double &x,double &y) const;
  void Inverse(const double x,const double y,double &lon,double &lat) const;
  double GetLat0() const {return _lat0;}
  double GetLon0() const {return _lon0;}
};

///////////////////////////////////////////////////////////////////
/// \brief ESRI ASCII grid in geographic (lon/lat) coordinates, sampled bilinearly
//
class CGWRaster
{
private:
  int _ncols,_nrows;
  double _xll,_yll,_cs,_nodata;
  vector<double> _v;          ///< row-major, north row first
  bool Valid(int r,int c) const;
public:
  CGWRaster();
  bool Read  (const string &filename,string &err);
  bool Sample(const double lon,const double lat,double &val) const;
  bool IsLoaded() const {return (_ncols>0);}
};

///////////////////////////////////////////////////////////////////
/// \brief map projection of an existing MODFLOW model's coordinates (metres): Transverse Mercator (UTM), Albers
/// equal-area conic or Lambert conformal conic on the WGS84 or GRS80 (NAD83) ellipsoid. Specified by an EPSG code
/// (common codes only) or by its parameters; converts model coordinates to longitude/latitude without external libraries
//
class CGWCRS
{
private:
  int    _type;                 ///< 0 none, 1 Transverse Mercator, 2 Albers, 3 Lambert conformal conic
  double _a,_f,_e,_e2;          ///< ellipsoid
  double _lat0,_lon0,_lat1,_lat2,_k0,_fe,_fn; ///< projection parameters [rad] and false easting/northing [m]
  double _n,_A,_b[5],_al[5],_d[5],_xi0; ///< Transverse Mercator (Kruger series, Karney 2011)
  double _nc,_C,_F,_rho0;       ///< conic projections
  double Q  (const double phi) const;
  double M  (const double phi) const;
  double T  (const double phi) const;
public:
  CGWCRS();
  bool   Parse  (const string &spec,string &err); ///< "EPSG:32610" or "TMERC lat0 lon0 k0 FE FN [WGS84|GRS80]", "ALBERS lat1 lat2 lat0 lon0 FE FN [..]", "LCC ..."
  bool   IsSet  () const {return (_type!=0);}
  void   Inverse(const double x,const double y,double &lon,double &lat) const; ///< [m] -> [deg]
  void   Forward(const double lon,const double lat,double &x,double &y) const; ///< [deg] -> [m]
  string Describe() const;
};
namespace GWGeom
{
  bool   ReadGeoJSONPolygons(const string &filename,const string &idField,
                             map<long long,vector<gw_ring> > &shapes,string &err);
  bool   ReadGeoJSONLines   (const string &filename,const string &idField,
                             vector<gw_ring> &lines,vector<long long> &ids,string &err);
  double SignedArea         (const gw_ring &r);
  double ClipRingAreaToRect (const gw_ring &r,const double x0,const double y0,const double x1,const double y1);
}
#endif
