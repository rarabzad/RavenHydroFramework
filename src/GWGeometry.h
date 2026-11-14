/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: spatial utilities
  (GeoJSON reading, equal-area projection, ESRI ASCII and GeoTIFF rasters,
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

class CGWCRS;
///////////////////////////////////////////////////////////////////
/// \brief raster map (ESRI ASCII grid, or GeoTIFF for a .tif/.tiff file) in geographic (lon/lat) coordinates, or in
/// the projected coordinates of an :InputCRS, sampled bilinearly
//
class CGWRaster
{
private:
  int _ncols,_nrows;
  double _xll,_yll,_dx,_dy,_nodata; ///< lower-left corner, cell width and height (equal in an ESRI ASCII grid), NODATA value
  bool _hasNodata;            ///< a NODATA value is given (always for ESRI ASCII grids: -9999 by default)
  vector<double> _v;          ///< row-major, north row first
  const CGWCRS *_crs;         ///< projection of the raster's coordinates (NULL: longitude/latitude)
  string _cmd;                ///< the command and file name as the user wrote them (":LandSurface dem.tif")
  bool _metres;               ///< the user confirmed metres (METRES after the file name): a GeoTIFF band unit Raven does not read as metres is accepted
  string _note;               ///< what the caller should report as a warning (set by Read)
  string _format;             ///< encoding of a GeoTIFF file
  bool Valid(int r,int c) const;
  bool ReadTIFF(const string &filename,string &err);
  bool CheckLonLatRange(const string &filename,string &err);
public:
  CGWRaster();
  void SetCRS(const CGWCRS *crs){_crs=crs;}  ///< (before Read)
  void SetMetres(const bool m){_metres=m;}    ///< (before Read)
  void SetCommand(const string &c){_cmd=c;}   ///< the command and file name as the user wrote them, for messages (before Read)
  bool Read  (const string &filename,string &err);
  bool Sample(const double lon,const double lat,double &val) const;
  bool IsLoaded() const {return (_ncols>0);}
  void Extent(double *b) const;              ///< xmin, ymin, xmax, ymax in the raster's own coordinates
  const string &Note  () const {return _note;}   ///< warnings for the caller, one per line
  const string &Format() const {return _format;} ///< GeoTIFF encoding ("" for ESRI ASCII)
  static const double MAX_PIXELS;            ///< largest raster read [cells]
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
  bool   Same   (const CGWCRS &o) const;        ///< the same projection (parameters within round-off)
  long   EPSG   () const {return _epsg;}        ///< EPSG code it was given by (0: given by its parameters, or ESRI:102001)
  const string &Code() const {return _code;}  ///< the code it was given by ("EPSG:32610", "ESRI:102001"; "": given by its parameters)
private:
  long   _epsg;
  string _code;
};
namespace GWGeom
{
  //(crs: projection of the file's coordinates, :InputCRS; NULL: longitude/latitude. Vertices are returned in lon/lat;
  // rawbox, if given: xmin, ymin, xmax, ymax of the coordinates as written in the file)
  bool   ReadGeoJSONPolygons(const string &filename,const string &idField,
                             map<long long,vector<gw_ring> > &shapes,string &err,const CGWCRS *crs=NULL,double *rawbox=NULL,
                             map<long long,vector<int> > *parts=NULL,  ///< parts: for each ring, the polygon it belongs to (rings of one polygon: its outline and holes)
                             map<long long,string> *noArea=NULL);      ///< noArea: the IDs left out because none of their polygons has an
                                                                       ///< area (rings of fewer than 3 distinct vertices or on one line, empty geometries), with why
  bool   ReadGeoJSONLines   (const string &filename,const string &idField,
                             vector<gw_ring> &lines,vector<long long> &ids,string &err,const CGWCRS *crs=NULL,double *rawbox=NULL);
  double SignedArea         (const gw_ring &r);
  bool   GeographicEPSG     (const long code);  ///< an EPSG code of a geographic (longitude/latitude) coordinate system
  double ClipRingAreaToRect (const gw_ring &r,const double x0,const double y0,const double x1,const double y1);
}
#endif
