/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: GeoTIFF reader for the land-surface and bedrock rasters
  (single band; classic TIFF and BigTIFF, both byte orders, strips and tiles; no compression, LZW, DEFLATE and
   PackBits with Raven's own decoders; predictors 1, 2 and 3; 8-64 bit integers and 32/64 bit floats;
   north-up georeferencing from ModelPixelScale + ModelTiepoint or ModelTransformation; GDAL_NODATA and a GDAL internal
   mask; EPSG code and units of the GeoKeyDirectory). No external library is used.
----------------------------------------------------------------*/
#ifndef GWGEOTIFF_H
#define GWGEOTIFF_H

#include <string>
#include <vector>
using namespace std;

///////////////////////////////////////////////////////////////////
/// \brief a single-band GeoTIFF raster as read by GWGeoTIFF::Read
//
struct gw_geotiff
{
  int    ncols,nrows;
  double xll,yll;      ///< lower-left corner of the lower-left pixel, in the raster's coordinates
  double dx,dy;        ///< pixel width and height (positive)
  bool   hasNodata;    ///< GDAL_NODATA tag given
  double nodata;       ///< its value (a 32-bit float raster: rounded to 32-bit float, as the pixels are)
  vector<double> v;    ///< ncols*nrows values, row-major, north row first (infinite values stored as NaN)
  int    modelType;    ///< GTModelTypeGeoKey: 0 not given, 1 projected, 2 geographic, 3 geocentric, 32767 user-defined
  int    rasterType;   ///< GTRasterTypeGeoKey: 1 PixelIsArea (also when not given), 2 PixelIsPoint
  int    epsg;         ///< EPSG code of the raster's CRS (0: none, or user-defined)
  bool   geographic;   ///< epsg is a geographic (longitude/latitude) CRS
  string format;       ///< description of the encoding, e.g. "BigTIFF, big-endian, tiled 16x16, LZW, predictor 2, 32-bit float"
  string notes;        ///< warnings for the caller, one per line: values that often mark missing data in a raster without
                       ///< GDAL_NODATA, a vertical system of unknown unit, a damaged chain of image directories (no mask used)
  gw_geotiff();
};

namespace GWGeoTIFF
{
  const double             MAX_PIXELS    =5.0e7;                ///< largest raster (pixels; 400 MB of values)
  const unsigned long long MAX_FILE_BYTES=1024ULL*1024ULL*1024ULL; ///< largest file [bytes]

  bool IsTIFFName(const string &filename);    ///< extension .tif, .tiff or .gtiff (any case)
  bool Read      (const string &filename,gw_geotiff &R,string &err);
  bool ReadMemory(const unsigned char *data,const size_t n,const string &name,gw_geotiff &R,string &err);

  //the decoders (exposed for the unit tests): each writes at most `limit` bytes to `out` and returns false with a
  //reason on corrupt input; `full` is set when the output reached `limit` before the end of the encoded data
  bool InflateZlib(const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err);
  bool LZWDecode  (const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err);
  bool PackBits   (const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err);
}
#endif
