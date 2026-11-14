/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: spatial utilities
----------------------------------------------------------------*/
#include "GWGeometry.h"
#include "GWGeoTIFF.h"
#include "GWEPSGTable.h"
#include "GWUtil.h"
#include <cstring>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <limits>

static const double GW_PI      =3.14159265358979323846;
static const double GW_D2R     =GW_PI/180.0;
static const double GW_R_EARTH =6371007.181; //authalic radius [m]

//=================================================================
// Projection
//=================================================================
CGWProjection::CGWProjection(){SetOrigin(0.0,0.0);}
void CGWProjection::SetOrigin(const double lat0,const double lon0)
{
  _lat0=lat0; _lon0=lon0;
  _sinp0=sin(lat0*GW_D2R); _cosp0=cos(lat0*GW_D2R);
}
void CGWProjection::Forward(const double lon,const double lat,double &x,double &y) const
{
  double lam=(lon-_lon0)*GW_D2R, phi=lat*GW_D2R;
  double denom=1.0+_sinp0*sin(phi)+_cosp0*cos(phi)*cos(lam);
  double k=sqrt(2.0/std::max(denom,1e-12));
  x=GW_R_EARTH*k*cos(phi)*sin(lam);
  y=GW_R_EARTH*k*(_cosp0*sin(phi)-_sinp0*cos(phi)*cos(lam));
}
void CGWProjection::Inverse(const double x,const double y,double &lon,double &lat) const
{
  double rho=sqrt(x*x+y*y);
  if (rho<1e-9){lon=_lon0;lat=_lat0;return;}
  double c=2.0*asin(std::min(1.0,rho/(2.0*GW_R_EARTH)));
  double phi=asin(cos(c)*_sinp0+y*sin(c)*_cosp0/rho);
  double lam=atan2(x*sin(c),rho*_cosp0*cos(c)-y*_sinp0*sin(c));
  lat=phi/GW_D2R; lon=_lon0+lam/GW_D2R;
  if (lon>180.0){lon-=360.0;} else if (lon<-180.0){lon+=360.0;} //(an origin near the antimeridian: longitudes within +-180)
}

//=================================================================
// Raster
//=================================================================
const double CGWRaster::MAX_PIXELS=GWGeoTIFF::MAX_PIXELS;
CGWRaster::CGWRaster(){_ncols=_nrows=0;_xll=_yll=0;_dx=_dy=1;_nodata=-9999;_hasNodata=true;_crs=NULL;_metres=false;}
bool CGWRaster::Read(const string &filename,string &err)
{
  _note=""; _format="";
  if (GWGeoTIFF::IsTIFFName(filename)){return ReadTIFF(filename,err);}
  ifstream IN(filename.c_str());
  if (IN.fail()){err="cannot open raster file "+filename;return false;}
  //ESRI ASCII grid: the header keys (any order, any case) up to the first value; GDAL writes dx and dy instead of
  //cellsize for rectangular pixels, and NODATA_value nan for a floating-point raster whose missing pixels are NaN
  if (IN.peek()==0xEF){char b[3]={0,0,0}; IN.read(b,3); if (!((unsigned char)b[1]==0xBB && (unsigned char)b[2]==0xBF)){IN.clear(); IN.seekg(0);}} //(a UTF-8 byte-order mark)
  bool center=false,xc=false,yc=false; string tok; double cs=-1,dx=-1,dy=-1; double ncols=0,nrows=0; int nkey=0;
  bool hasX=false,hasY=false; int seen[10]={0,0,0,0,0,0,0,0,0,0};
  _nodata=-9999; _hasNodata=false; _ncols=_nrows=0; //(no NODATA_value: every number is a value, as GDAL reads the grid)
  while (true)
  {
    streampos pos=IN.tellg();
    if (!(IN>>tok)){err="raster "+filename+" has a header but no values";return false;}
    double v;
    if (gwutil::StrictNumber(tok.c_str(),v) || (tok[0]=='-') || (tok[0]=='+') || (tok[0]=='.') || isdigit((unsigned char)tok[0])){IN.seekg(pos); break;} //(the first value)
    string lk=tok; std::transform(lk.begin(),lk.end(),lk.begin(),::tolower);
    if ((lk=="nan") && (nkey>0)){IN.seekg(pos); break;}                                     //(a first value that is NaN)
    const char *keys[]={"ncols","nrows","xllcorner","yllcorner","xllcenter","yllcenter","cellsize","dx","dy","nodata_value"};
    int kk=-1; for (int i=0;i<10;i++){if (lk==keys[i]){kk=i;}}
    if (kk<0){err="raster "+filename+": unknown header key '"+tok+"' (an ESRI ASCII grid has ncols, nrows, xllcorner or xllcenter, yllcorner or yllcenter, cellsize or dx and dy, and NODATA_value)";return false;}
    if (seen[kk] || ((kk>=2) && (kk<=5) && seen[(kk<4)?kk+2:kk-2])){ //(a key given twice, or a corner and a centre for one axis: which one holds is not clear)
      const string other=seen[kk]?tok:string(keys[(kk<4)?kk+2:kk-2]);
      err="raster "+filename+": the header gives "+(seen[kk]?"the key "+tok+" twice":tok+" and "+other+" (a corner and a centre of the same axis)")+"; give each key once";
      return false;
    }
    seen[kk]=1;
    string sv;
    if (!(IN>>sv)){err="raster "+filename+": the header key "+tok+" has no value";return false;}
    string ls=sv; std::transform(ls.begin(),ls.end(),ls.begin(),::tolower);
    if ((kk==9) && ((ls=="nan") || (ls=="-nan") || (ls=="+nan"))){_nodata=NAN; _hasNodata=false; nkey++; continue;} //(NaN pixels are never valid anyway)
    if (!gwutil::StrictNumber(sv.c_str(),v)){err="raster "+filename+": the value of the header key "+tok+" ('"+sv+"') is not a number";return false;}
    if ((kk<=1) && ((v!=floor(v)) || (v<=0) || (v>2.0e9))){err="raster "+filename+": "+tok+" "+sv+" is not a positive whole number";return false;}
    switch (kk){
      case 0: ncols=v; break;           case 1: nrows=v; break;
      case 2: _xll=v; hasX=true; break; case 3: _yll=v; hasY=true; break;
      case 4: _xll=v; hasX=true; xc=true; break; case 5: _yll=v; hasY=true; yc=true; break;
      case 6: cs=v; break; case 7: dx=v; break; case 8: dy=v; break;
      default: _nodata=v; _hasNodata=true;
    }
    nkey++;
  }
  (void)center;
  if ((ncols<=0) || (nrows<=0) || !hasX || !hasY){err="raster "+filename+": the header lacks "+string((ncols<=0)?"ncols":(nrows<=0)?"nrows":!hasX?"xllcorner (or xllcenter)":"yllcorner (or yllcenter)");return false;}
  if (cs>0){_dx=_dy=cs;}
  else if ((dx>0) && (dy>0)){_dx=dx; _dy=dy;}
  else {err="raster "+filename+": the header needs a positive cellsize, or positive dx and dy";return false;}
  _ncols=(int)ncols; _nrows=(int)nrows; cs=_dx;
  if ((double)_ncols*(double)_nrows>MAX_PIXELS){
    ostringstream s; s<<"raster "<<filename<<" has "<<_ncols<<" x "<<_nrows<<" cells, more than the "<<MAX_PIXELS/1e6<<" million Raven reads: resample it to a "
                      <<"coarser resolution (a few raster cells per model cell are enough; each cell averages 25 samples)";
    err=s.str(); _ncols=0; return false;
  }
  if (xc){_xll-=0.5*_dx;}                              //(each axis by its own key: a centre or a corner)
  if (yc){_yll-=0.5*_dy;}
  if (!CheckLonLatRange(filename,err)){return false;}
  _v.resize((size_t)_ncols*_nrows);
  for (size_t i=0;i<_v.size();i++){
    if (!(IN>>tok)){
      ostringstream s; s<<"raster "<<filename<<" has "<<i<<" values, fewer than ncols*nrows ("<<_ncols<<" x "<<_nrows<<")";
      err=s.str();_ncols=0;return false;
    }
    while (!tok.empty() && ((tok[tok.size()-1]=='\x1a') || (tok[tok.size()-1]=='\0'))){tok.erase(tok.size()-1);} //(a DOS end-of-file mark, as GDAL reads it)
    double v; const char *c=tok.c_str();
    if (!gwutil::StrictNumber(c,v)){
      string l=tok; std::transform(l.begin(),l.end(),l.begin(),::tolower);
      if ((l=="nan") || (l=="-nan") || (l=="+nan")){v=NAN;}
      else {
        ostringstream s; s<<"raster "<<filename<<": the value '"<<tok<<"' (row "<<i/_ncols+1<<", column "<<i%_ncols+1<<") is not a number";
        err=s.str();_ncols=0;return false;
      }
    }
    _v[i]=v;
  }
  while ((IN>>tok) && (tok.find_first_not_of(string("\x1a\0",2))==string::npos)){} //(trailing end-of-file marks)
  if (IN && (tok.find_first_not_of(string("\x1a\0",2))!=string::npos)){
    ostringstream s; s<<"raster "<<filename<<" has more values than ncols*nrows ("<<_ncols<<" x "<<_nrows<<"; the next is '"<<tok<<"'): check ncols and nrows in its header";
    err=s.str();_ncols=0;return false;
  }
  if (!_hasNodata && (_nodata==-9999)){ //(no NODATA_value key: -9999 cells are values, as for GDAL; reported, as they often mark no data)
    size_t c=0; for (size_t i=0;i<_v.size();i++){if (_v[i]==-9999.0){c++;}}
    if (c>0){
      ostringstream s; s<<"raster "<<filename<<" has no NODATA_value but holds "<<c<<((c==1)?" cell":" cells")<<" of -9999, a value that often marks missing data; "
        "they are read as elevations. If they are no data, add the line NODATA_value -9999 to its header";
      _note+=(_note.empty()?"":"\n")+s.str();
    }
  }
  //GDAL's side file <file>.aux.xml, read as for a GeoTIFF (GWGeoTIFF::BandMeta): GDAL's ESRI ASCII driver takes its Scale,
  //Offset and UnitType; it keeps the header's NODATA_value and takes no NoDataValue from the side file
  { string pam,note,pnd; bool has=false,hasPnd=false; double sc=1.0,of=0.0;
    if (!GWGeoTIFF::ReadSideFile(filename,pam,has,err)){_ncols=0; return false;}
    if (!GWGeoTIFF::BandMeta("raster "+filename+": ","",has?&pam:NULL,filename+".aux.xml",_metres,_cmd,sc,of,hasPnd,pnd,note,err)){_ncols=0; return false;}
    if (note!=""){_note+=(_note.empty()?"":"\n")+note;}
    if ((sc!=1.0) || (of!=0.0)){ //(no-data cells become NaN, so that no scaled value is taken for NODATA)
      for (size_t i=0;i<_v.size();i++){
        double &v=_v[i];
        if ((v!=v) || (_hasNodata && (fabs(v-_nodata)<=1e-6))){v=NAN; continue;}
        v=v*sc+of; if (!std::isfinite(v)){v=NAN;}
      }
      _hasNodata=false; _nodata=0.0;
    }
  }
  return true;
}
//without a projection (:InputCRS) the raster's extent must lie in the longitude/latitude range
bool CGWRaster::CheckLonLatRange(const string &filename,string &err)
{
  if ((_crs==NULL) && ((_xll<-360.0) || (_xll+_ncols*_dx>360.0) || (_yll<-90.0-_dy) || (_yll+_nrows*_dy>90.0+_dy))){ //(projected metres read as degrees)
    ostringstream s; s<<"raster "<<filename<<" lies outside the longitude/latitude range (lower-left corner "<<_xll<<", "<<_yll<<", cell size "<<_dx;
    if (_dy!=_dx){s<<" x "<<_dy;}
    s<<"): a raster in projected coordinates needs :InputCRS in the .rvg file";
    err=s.str(); _ncols=0; return false;
  }
  return true;
}
//a GeoTIFF file: read natively (GWGeoTIFF); its EPSG code, if it has one, must agree with :InputCRS (a geographic code
//without :InputCRS; with it, the same code or the same built-in projection); a user-defined CRS or none is treated as an
//ESRI ASCII grid is (the extent decides)
bool CGWRaster::ReadTIFF(const string &filename,string &err)
{
  gw_geotiff T; _ncols=0;
  if (!GWGeoTIFF::Read(filename,T,err,_metres,_cmd)){return false;}
  if (T.epsg>0){
    ostringstream cs; cs<<"EPSG:"<<T.epsg;
    string inp=(_crs==NULL)?"":((_crs->Code()!="")?_crs->Code()+" ("+_crs->Describe()+")":_crs->Describe());
    CGWCRS f; string e2;
    const bool built=!T.geographic && f.Parse(cs.str(),e2);       //(a projection Raven has built in)
    if (T.geographic){
      if (_crs!=NULL){
        err="raster "+filename+" is in longitude/latitude ("+cs.str()+" in its GeoTIFF keys), but :InputCRS is "+inp+
            ": give the raster in the :InputCRS projection (or remove :InputCRS if every input file is in longitude/latitude)";
        return false;
      }
    }
    else if (_crs==NULL){
      err="raster "+filename+" is in "+cs.str()+" (a projection, in its GeoTIFF keys), but there is no :InputCRS, so rasters are expected in longitude/latitude: "+
          (built?"give :InputCRS "+cs.str()+" in the .rvg file":"give :InputCRS with the projection's parameters (TMERC, ALBERS or LCC; "+cs.str()+" is not built in)")+
          ", or reproject the raster to longitude/latitude";
      return false;
    }
    else if (_crs->EPSG()!=T.epsg){
      if (built && !f.Same(*_crs)){
        err="raster "+filename+" is in "+cs.str()+" ("+f.Describe()+", in its GeoTIFF keys), but :InputCRS is "+inp+": all input files must use the projection of :InputCRS";
        return false;
      }
      if (!built && (_crs->Code()!="")){
        err="raster "+filename+" is in "+cs.str()+" (in its GeoTIFF keys), but :InputCRS is "+inp+": all input files must use the projection of :InputCRS";
        return false;
      }
      if (!built){_note="raster "+filename+" is in "+cs.str()+", a projection that is not built in, so it cannot be compared with :InputCRS ("+inp+"); it is read in the :InputCRS coordinates";}
    }
  }
  _ncols=T.ncols; _nrows=T.nrows; _xll=T.xll; _yll=T.yll; _dx=T.dx; _dy=T.dy;
  _hasNodata=T.hasNodata && (T.nodata==T.nodata); _nodata=T.nodata; _format=T.format; //(a NaN NODATA value: NaN pixels, never valid anyway)
  if (T.notes!=""){_note+=(_note.empty()?"":"\n")+T.notes;}
  if (!CheckLonLatRange(filename,err)){return false;}
  _v.swap(T.v);
  return true;
}
void CGWRaster::Extent(double *b) const
{
  b[0]=_xll; b[1]=_yll; b[2]=_xll+_ncols*_dx; b[3]=_yll+_nrows*_dy;
}
bool CGWRaster::Valid(int r,int c) const
{
  if ((r<0) || (c<0) || (r>=_nrows) || (c>=_ncols)){return false;}
  double v=_v[(size_t)r*_ncols+c];
  if (!_hasNodata){return (v==v);}  //(no NODATA value: any number is valid; NaN never is)
  return (fabs(v-_nodata)>1e-6);
}
bool CGWRaster::Sample(const double lon,const double lat,double &val) const
{
  if (_ncols==0){return false;}
  double px=lon,py=lat;
  if (_crs!=NULL){_crs->Forward(lon,lat,px,py);} //(a projected raster: interpolated in its own coordinates)
  double ytop=_yll+_nrows*_dy;
  double fc=(px-_xll)/_dx-0.5;
  if ((_crs==NULL) && ((fc<-1.0) || (fc>_ncols))){ //(a longitude/latitude raster across the antimeridian, written beyond +-180)
    for (int k=-1;k<=1;k+=2){
      const double f=(px+k*360.0-_xll)/_dx-0.5;
      if ((f>=-1.0) && (f<=_ncols)){fc=f; break;}
    }
  }
  double fr=(ytop-py)/_dy-0.5;
  if ((fc<-1.0) || (fr<-1.0) || (fc>_ncols) || (fr>_nrows)){return false;}
  fc=std::max(0.0,std::min((double)(_ncols-1),fc));
  fr=std::max(0.0,std::min((double)(_nrows-1),fr));
  int c0=(int)floor(fc),r0=(int)floor(fr);
  int c1=std::min(c0+1,_ncols-1),r1=std::min(r0+1,_nrows-1);
  double wc=fc-c0,wr=fr-r0;
  int    rr[4]={r0,r0,r1,r1}, cc[4]={c0,c1,c0,c1};
  double ww[4]={(1-wr)*(1-wc),(1-wr)*wc,wr*(1-wc),wr*wc};
  double sum=0,wsum=0;
  for (int i=0;i<4;i++){
    if (Valid(rr[i],cc[i])){sum+=ww[i]*_v[(size_t)rr[i]*_ncols+cc[i]];wsum+=ww[i];}
  }
  if (wsum<1e-9){return false;}
  val=sum/wsum;
  return true;
}

//=================================================================
// Minimal JSON parser (sufficient for GeoJSON)
//=================================================================
namespace
{
  enum jtype{J_NULL,J_BOOL,J_NUM,J_STR,J_ARR,J_OBJ};
  struct JVal
  {
    jtype t; double n; string s;
    vector<JVal> a;
    vector<pair<string,JVal> > o;
    JVal():t(J_NULL),n(0){}
    const JVal *Get(const string &k) const{
      for (size_t i=0;i<o.size();i++){if (o[i].first==k){return &o[i].second;}}
      return NULL;
    }
  };
  class JParser
  {
    const char *p,*e,*b0;
    void ws(){while((p<e) && isspace((unsigned char)*p)){p++;}}
    //the four hexadecimal digits after \u (p at the u; left at the last digit)
    bool Hex4(unsigned long &cp){
      if (e-p<5){return false;}
      cp=0; for (int k=1;k<=4;k++){const char h=p[k]; cp<<=4;
        if ((h>='0') && (h<='9')){cp+=(unsigned long)(h-'0');} else if ((h>='a') && (h<='f')){cp+=(unsigned long)(h-'a'+10);}
        else if ((h>='A') && (h<='F')){cp+=(unsigned long)(h-'A'+10);} else {return false;}}
      p+=4; return true;
    }
    bool str(string &out)
    {
      ws(); if ((p>=e) || (*p!='"')){return false;}
      p++; out.clear();
      while ((p<e) && (*p!='"')){
        if (*p=='\\'){
          p++; if (p>=e){return false;}
          switch(*p){
            case 'n':out+='\n';break; case 't':out+='\t';break; case 'r':out+='\r';break;
            case 'b':out+='\b';break; case 'f':out+='\f';break;
            case 'u':{unsigned long cp; if (!Hex4(cp)){return false;}           //(\uXXXX, a pair of them for a surrogate pair: UTF-8)
              if ((cp>=0xD800) && (cp<=0xDBFF)){
                unsigned long lo; if ((e-p<3) || (p[1]!='\\') || (p[2]!='u')){return false;} p+=2;
                if (!Hex4(lo) || (lo<0xDC00) || (lo>0xDFFF)){return false;}
                cp=0x10000+((cp-0xD800)<<10)+(lo-0xDC00);
              }
              else if ((cp>=0xDC00) && (cp<=0xDFFF)){return false;}
              if (cp<0x80){out+=(char)cp;}
              else if (cp<0x800){out+=(char)(0xC0|(cp>>6)); out+=(char)(0x80|(cp&0x3F));}
              else if (cp<0x10000){out+=(char)(0xE0|(cp>>12)); out+=(char)(0x80|((cp>>6)&0x3F)); out+=(char)(0x80|(cp&0x3F));}
              else {out+=(char)(0xF0|(cp>>18)); out+=(char)(0x80|((cp>>12)&0x3F)); out+=(char)(0x80|((cp>>6)&0x3F)); out+=(char)(0x80|(cp&0x3F));}
              break;}
            default :out+=*p;
          }
          p++;
        }
        else {out+=*p++;}
      }
      if (p>=e){return false;}
      p++; return true;
    }
  public:
    JParser(const char *b,const char *en):p(b),e(en),b0(b){}
    //where parsing stopped (the first token that is not valid JSON): line, column and the text there; for a token that
    //looks like a number JSON does not have, what JSON allows
    string Where() const{
      int line=1,col=1; for (const char *q=b0;(q<p) && (q<e);q++){if (*q=='\n'){line++; col=1;} else {col++;}}
      string t; for (const char *q=p;(q<e) && (t.size()<12) && (*q!='\n') && (*q!='\r');q++){t+=*q;}
      string u=t; if (!u.empty() && ((u[0]=='-') || (u[0]=='+'))){u.erase(0,1);}
      for (size_t i=0;i<u.size();i++){u[i]=(char)tolower((unsigned char)u[i]);}
      string hint;
      const bool after0=(p>b0) && (p<e) && (p[-1]=='0') && ((p-1==b0) || !isdigit((unsigned char)p[-2])) && ((*p=='x') || (*p=='X') || isdigit((unsigned char)*p)); //(0x10, 0122: parsing stops after the 0)
      if (after0){hint=" (hexadecimal numbers and numbers with leading zeros are not JSON)";}
      else if ((u.compare(0,3,"nan")==0) || (u.compare(0,3,"inf")==0)){hint=" (Raven reads non-finite values written as NaN, Infinity and -Infinity only)";}
      return " at line "+to_string((long long)line)+", column "+to_string((long long)col)+((p<e)?" ('"+t+"')":string(" (the end of the file)"))+hint;
    }
    //only blanks after the root value
    bool AtEnd(){ws(); return (p>=e);}
    bool Parse(JVal &v)
    {
      ws(); if (p>=e){return false;}
      char c=*p;
      if (c=='{'){
        v.t=J_OBJ; p++; ws();
        if ((p<e) && (*p=='}')){p++;return true;}
        while (true){
          string k; if (!str(k)){return false;}
          ws(); if ((p>=e) || (*p!=':')){return false;} p++;
          v.o.push_back(make_pair(k,JVal()));
          if (!Parse(v.o.back().second)){return false;}
          ws(); if (p>=e){return false;}
          if (*p==','){p++;continue;}
          if (*p=='}'){p++;return true;}
          return false;
        }
      }
      else if (c=='['){
        v.t=J_ARR; p++; ws();
        if ((p<e) && (*p==']')){p++;return true;}
        while (true){
          v.a.push_back(JVal());
          if (!Parse(v.a.back())){return false;}
          ws(); if (p>=e){return false;}
          if (*p==','){p++;continue;}
          if (*p==']'){p++;return true;}
          return false;
        }
      }
      else if (c=='"'){v.t=J_STR; return str(v.s);}
      else if (c=='t'){if ((e-p<4) || strncmp(p,"true",4)) {return false;} v.t=J_BOOL;v.n=1;p+=4;return true;}
      else if (c=='f'){if ((e-p<5) || strncmp(p,"false",5)){return false;} v.t=J_BOOL;v.n=0;p+=5;return true;}
      else if (c=='n'){if ((e-p<4) || strncmp(p,"null",4)) {return false;} v.t=J_NULL;p+=4;return true;}
      else if ((c=='N') && (e-p>=3) && !strncmp(p,"NaN",3)){v.t=J_NUM; v.n=std::numeric_limits<double>::quiet_NaN(); v.s="NaN"; p+=3; return true;}
      else if ((c=='I') && (e-p>=8) && !strncmp(p,"Infinity",8)){v.t=J_NUM; v.n=HUGE_VAL; v.s="Infinity"; p+=8; return true;}
      else if ((c=='-') && (e-p>=9) && !strncmp(p,"-Infinity",9)){v.t=J_NUM; v.n=-HUGE_VAL; v.s="-Infinity"; p+=9; return true;}
      else { //a number in JSON's grammar: -?(0|[1-9][0-9]*)(.[0-9]+)?([eE][+-]?[0-9]+)? (not hexadecimal or leading zeros);
             //NaN, Infinity and -Infinity (above), which JSON does not have, are read as GDAL and Python's json read them, so
             //that an attribute table with such values is read; a coordinate must be finite (JCoords::Pair)
        const char *q=p;
        if ((q<e) && (*q=='-')){q++;}
        if ((q<e) && (*q=='0')){q++;}
        else if ((q<e) && (*q>='1') && (*q<='9')){while ((q<e) && isdigit((unsigned char)*q)){q++;}}
        else {return false;}
        if ((q<e) && (*q=='.')){q++; if (!((q<e) && isdigit((unsigned char)*q))){return false;} while ((q<e) && isdigit((unsigned char)*q)){q++;}}
        if ((q<e) && ((*q=='e') || (*q=='E'))){
          q++; if ((q<e) && ((*q=='+') || (*q=='-'))){q++;}
          if (!((q<e) && isdigit((unsigned char)*q))){return false;} while ((q<e) && isdigit((unsigned char)*q)){q++;}
        }
        v.t=J_NUM; v.s=string(p,q); v.n=strtod(v.s.c_str(),NULL);
        p=q; return true;
      }
    }
  };
  //coordinates of a file: projected ones (:InputCRS) go to lon/lat as they are read; without a projection the largest
  //|x| and |y| are kept, to recognise projected coordinates read as degrees
  //(the extent of the coordinates as written in the file is kept too, for the longitude/latitude check of :InputCRS)
  struct JCoords
  {
    const CGWCRS *crs; double ax,ay,box[4]; string unknown; //(unknown: a "crs" member naming a code that is not classified)
    string feature,ring,bad; //(the feature and ring being read, for messages; the first coordinate that is not a finite number)
    JCoords(const CGWCRS *c):crs(c),ax(0),ay(0){box[0]=box[1]=1e300; box[2]=box[3]=-1e300;}
    //a coordinate pair [x, y, ...]: both must be finite JSON numbers (a string, null, true or false is no coordinate)
    bool Pair(const JVal &pt,const string &where,double &x,double &y){
      if ((pt.t!=J_ARR) || (pt.a.size()<2)){
        if (bad.empty()){
          bad=feature+", "+where+": it is "+((pt.t==J_ARR)?"an array of "+to_string((unsigned long long)pt.a.size())+" value"+((pt.a.size()==1)?"":"s"):What(pt))+
              ", not a position [x, y] (are the coordinates nested as GeoJSON's geometry type needs?)";
        }
        return false;
      }
      for (int k=0;k<2;k++){
        const JVal &c=pt.a[k];
        if ((c.t==J_NUM) && std::isfinite(c.n)){continue;}
        if (bad.empty()){bad=feature+", "+where+": its "+string(k?"y (latitude)":"x (longitude)")+" is "+What(c)+", not a "+string((c.t==J_NUM)?"finite ":"")+"number";}
        return false;
      }
      x=pt.a[0].n; y=pt.a[1].n; return true;
    }
    //a geometry member that must be an array (of rings, of polygons, of lines)
    bool Array(const JVal &v,const string &where){
      if (v.t==J_ARR){return true;}
      if (bad.empty()){bad=feature+", "+where+": it is "+What(v)+", not an array (are the coordinates nested as GeoJSON's geometry type needs?)";}
      return false;
    }
    static string What(const JVal &c){
      return (c.t==J_STR)?"the string \""+c.s+"\"":(c.t==J_NULL)?string("null"):(c.t==J_BOOL)?string(c.n?"true":"false"):
             (c.t==J_NUM)?(std::isfinite(c.n)?"the number "+c.s:(c.s.size() && (isalpha((unsigned char)c.s[c.s.size()-1])))?c.s:string("a number too large")):
             (c.t==J_ARR)?string("an array"):string("an object");
    }
    gw_pt Point(const double x,const double y){
      gw_pt q; q.x=x; q.y=y;
      box[0]=std::min(box[0],x); box[1]=std::min(box[1],y); box[2]=std::max(box[2],x); box[3]=std::max(box[3],y);
      if (crs!=NULL){crs->Inverse(x,y,q.x,q.y);}
      else {ax=std::max(ax,fabs(x)); ay=std::max(ay,fabs(y));}
      return q;
    }
    void Box(double *b) const {if (b!=NULL){for (int i=0;i<4;i++){b[i]=box[i];}}}
  };
  //a feature's ID: a whole number, or a string that holds one (blanks around it and a decimal point followed by zeros allowed:
  //" 7", "7.0"); null, true, false, a fraction or other text are no ID
  bool FeatureID(const JVal &v,long long &id)
  {
    if (v.t==J_NUM){
      if (!std::isfinite(v.n) || (v.n!=floor(v.n)) || (fabs(v.n)>9.0e15)){return false;}
      id=(long long)v.n; return true;
    }
    if (v.t!=J_STR){return false;}
    const size_t a=v.s.find_first_not_of(" \t"),b=v.s.find_last_not_of(" \t");     //(blanks around the number are left out)
    if (a==string::npos){return false;}
    string s=v.s.substr(a,b-a+1);
    const size_t dp=s.find('.');                                                      //(a decimal point followed by zeros only: "7.0")
    if (dp!=string::npos){if ((dp+1<s.size()) && (s.find_first_not_of('0',dp+1)!=string::npos)){return false;} s.erase(dp);}
    size_t i=((s.size()>0) && ((s[0]=='-') || (s[0]=='+')))?1:0;
    if ((i>=s.size()) || (s.size()>i+18)){return false;}
    for (size_t j=i;j<s.size();j++){if (!isdigit((unsigned char)s[j])){return false;}}
    id=atoll(s.c_str()); return true;
  }
  //the message for an ID that FeatureID refuses (a whole JSON number beyond 9e15, as GDAL writes an int64 field, is not
  //held exactly by the reader's numbers; a string holds up to 18 digits)
  string IDRefused(const JVal &v)
  {
    if ((v.t==J_NUM) && std::isfinite(v.n) && (v.n==floor(v.n)) && (fabs(v.n)>9.0e15)){
      const size_t d=(!v.s.empty() && (v.s[0]=='-'))?1:0;
      const bool digits=(v.s.size()>d) && (v.s.size()<=d+18) && (v.s.find_first_not_of("0123456789",d)==string::npos);
      return "the number "+v.s+", larger than 9e15, the largest ID Raven reads as a number ("+
             (digits?"write it as a string: \""+v.s+"\")":string("write the ID as a string of at most 18 digits)"));
    }
    return JCoords::What(v)+", not a whole number";
  }
  //the message for a missing ID field: GDAL's ID_FIELD creation option moves the field out of "properties" into the
  //feature's "id" member
  string IDMissing(const JVal &F,const string &idField,const string &filename)
  {
    string e="GeoJSON feature missing ID field '"+idField+"' in "+filename;
    if (F.Get("id")!=NULL){e+=" (the feature has an \"id\" member, where GDAL's ID_FIELD option puts the field: write the file without ID_FIELD, so that "+idField+" stays in \"properties\")";}
    return e;
  }
  bool RingFromJSON(const JVal &jr,gw_ring &r,JCoords &C)
  {
    r.clear();
    if (!C.Array(jr,C.ring.substr(0,C.ring.size()-2))){return false;}
    for (size_t i=0;i<jr.a.size();i++){
      const JVal &pt=jr.a[i];
      gw_pt q; if (!C.Pair(pt,C.ring+"vertex "+to_string((unsigned long long)(i+1)),q.x,q.y)){return false;}
      r.push_back(q);
    }
    if ((r.size()>1) && (r.front().x==r.back().x) && (r.front().y==r.back().y)){r.pop_back();}
    //the area is tested in the coordinates as written in the file (before a projected file is taken to longitude/latitude,
    //which would bend a straight ring into a lens):
    set<pair<double,double> > d; for (size_t i=0;(i<r.size()) && (d.size()<3);i++){d.insert(make_pair(r[i].x,r[i].y));}
    if (d.size()<3){return false;} //(a ring of fewer than 3 distinct vertices has no area)
    //(nor has one whose vertices lie on one line, that runs back over itself, or whose loops cancel (a figure eight of
    // equal lobes): its signed area, taken from the first vertex, is at most 1e-9 of the square of its extent, a sliver a
    // billionth of its length wide, well above the rounding of the coordinates)
    double A=0,x0=r[0].x,y0=r[0].y,xmin=x0,xmax=x0,ymin=y0,ymax=y0;
    for (size_t i=0;i<r.size();i++){
      const gw_pt &a=r[i],&b=r[(i+1)%r.size()];
      A+=(a.x-x0)*(b.y-y0)-(b.x-x0)*(a.y-y0);
      xmin=min(xmin,a.x); xmax=max(xmax,a.x); ymin=min(ymin,a.y); ymax=max(ymax,a.y);
    }
    const double ext=max(xmax-xmin,ymax-ymin);
    if (!(fabs(0.5*A)>1e-9*ext*ext)){return false;}
    for (size_t i=0;i<r.size();i++){r[i]=C.Point(r[i].x,r[i].y);}
    return true;
  }
  void AddPolygon(const JVal &jpoly,vector<gw_ring> &rings,JCoords &C,vector<int> *part,const int ipart,const string &where)
  {
    if (!C.Array(jpoly,where)){return;}
    for (size_t i=0;i<jpoly.a.size();i++){
      gw_ring r;
      C.ring="ring "+to_string((unsigned long long)(i+1))+((i==0)?" (outline), ":" (hole), ");
      if (!RingFromJSON(jpoly.a[i],r,C)){
        if (i==0){return;} //(an outline without area: the polygon has none, and its holes are not read)
        continue;
      }
      double A=GWGeom::SignedArea(r);
      bool outer=(i==0);
      if ((outer && (A<0)) || (!outer && (A>0))){std::reverse(r.begin(),r.end());}
      rings.push_back(r);
      if (part!=NULL){part->push_back(ipart);}
    }
  }
  //the GeoJSON "crs" member (GDAL: {"type":"name","properties":{"name":"urn:ogc:def:crs:EPSG::32610"}}) against :InputCRS.
  //The code is the number that ends the name (EPSG:4326, urn:ogc:def:crs:EPSG::4326, urn:ogc:def:crs:EPSG:6.6:4326,
  //http://www.opengis.net/def/crs/EPSG/0/4326); OGC CRS84, CRS83 and CRS27 are longitude/latitude.
  //- a geographic code: longitude/latitude, so it may not have :InputCRS;
  //- a projection that is built in: needs :InputCRS, and the same projection;
  //- any other code cannot be classified: the file is accepted, and its coordinates decide (the longitude/latitude
  //  range test without :InputCRS, the longitude/latitude check of the caller with it)
  bool CheckFileCRS(const JVal &root,const string &filename,const CGWCRS *crs,string &err,string &unknown)
  {
    const JVal *jc=root.Get("crs"); if ((jc==NULL) || (jc->t!=J_OBJ)){return true;}
    const JVal *pr=jc->Get("properties"); if ((pr==NULL) || (pr->t!=J_OBJ)){return true;}
    const JVal *nm=pr->Get("name"); if ((nm==NULL) || (nm->t!=J_STR)){return true;}
    string name=nm->s,u=name; for (size_t i=0;i<u.size();i++){u[i]=(char)toupper((unsigned char)u[i]);}
    while (!u.empty() && isspace((unsigned char)u[u.size()-1])){u.erase(u.size()-1);}
    bool geographic=false; long code=0;
    if ((u.find("CRS84")!=string::npos) || (u.find("CRS83")!=string::npos) || (u.find("CRS27")!=string::npos)){geographic=true;}
    else {
      if (u.find("EPSG")==string::npos){return true;}                      //(another authority: not classified)
      size_t pd=u.size(); while ((pd>0) && isdigit((unsigned char)u[pd-1])){pd--;}
      if (pd==u.size()){return true;}                                        //(no code at the end)
      code=atol(u.c_str()+pd);
      geographic=GWGeom::GeographicEPSG(code);
    }
    ostringstream cs; cs<<"EPSG:"<<code;
    if (geographic && (code>0)){ //(geographic codes whose longitudes are not degrees from Greenwich)
      for (int k=0;k<2;k++){
        const gw_epsg_name *tab=(k==0)?GW_EPSG_GEOG_PRIME_MERIDIAN:GW_EPSG_GEOG_NOT_DEGREES; int nt=(k==0)?GW_EPSG_GEOG_PRIME_MERIDIAN_N:GW_EPSG_GEOG_NOT_DEGREES_N;
        for (int i=0;i<nt;i++){
          if (tab[i].code!=code){continue;}
          err="GeoJSON file "+filename+" is in "+cs.str()+" (its \"crs\" member, "+name+"), whose "+((k==0)?"longitudes count from the prime meridian of "+string(tab[i].name)+" from Greenwich":"angles are in "+string(tab[i].name))+
              "; Raven reads longitude/latitude in degrees from Greenwich: reproject the file (e.g. to EPSG:4326)";
          return false;
        }
      }
    }
    if (geographic){
      if (crs!=NULL){err="GeoJSON file "+filename+" is in longitude/latitude (its \"crs\" member, "+name+"), but :InputCRS says the input files are projected: give the file in the :InputCRS projection, "+crs->Describe()+" (or remove :InputCRS if every input file is in longitude/latitude)";return false;}
      return true;
    }
    for (int i=0;(code>0) && (i<GW_EPSG_PROJECTED_NOT_METRE_N);i++){ //(a projection whose coordinates are not metres)
      if (GW_EPSG_PROJECTED_NOT_METRE[i].code!=code){continue;}
      err="GeoJSON file "+filename+" is in "+cs.str()+" (its \"crs\" member, "+name+"), a projection whose coordinates are in "+string(GW_EPSG_PROJECTED_NOT_METRE[i].name)+
          "; Raven reads projected coordinates in metres: reproject the file to a projection in metres";
      return false;
    }
    CGWCRS f; string e2;
    if (!f.Parse(cs.str(),e2)){unknown=cs.str()+" (its \"crs\" member, "+name+")"; return true;} //(not built in: cannot be compared)
    if (crs==NULL){err="GeoJSON file "+filename+" is in "+cs.str()+" (its \"crs\" member, "+name+"): give :InputCRS "+cs.str()+" in the .rvg file, or reproject the file to longitude/latitude";return false;}
    if (!f.Same(*crs)){err="GeoJSON file "+filename+" is in "+cs.str()+" ("+f.Describe()+"), but :InputCRS is "+crs->Describe()+": all input files must use the projection of :InputCRS";return false;}
    return true;
  }
  bool CheckLonLat(const JCoords &C,const string &filename,string &err)
  {
    if ((C.crs==NULL) && ((C.ax>360.0) || (C.ay>90.0))){
      if (C.unknown!=""){err="GeoJSON file "+filename+" has coordinates outside the longitude/latitude range and is in "+C.unknown+", a projection that is not built in: give :InputCRS with the projection's parameters (TMERC, ALBERS or LCC) in the .rvg file, or reproject the file to longitude/latitude";}
      else {err="GeoJSON file "+filename+" has coordinates outside the longitude/latitude range: a file in projected coordinates needs :InputCRS in the .rvg file";}
      return false;
    }
    return true;
  }
}

//EPSG codes of geographic (longitude/latitude) coordinate systems, 2D and 3D: the table of GWEPSGTable.h (from the EPSG
//registry). Codes of projections and of geocentric systems (metres) are not in it
bool GWGeom::GeographicEPSG(const long code)
{
  return std::binary_search(GW_EPSG_GEOGRAPHIC,GW_EPSG_GEOGRAPHIC+GW_EPSG_GEOGRAPHIC_N,(int)code);
}
double GWGeom::SignedArea(const gw_ring &r)
{
  double A=0; size_t n=r.size();
  for (size_t i=0;i<n;i++){
    const gw_pt &a=r[i],&b=r[(i+1)%n];
    A+=a.x*b.y-b.x*a.y;
  }
  return 0.5*A;
}

bool GWGeom::ReadGeoJSONPolygons(const string &filename,const string &idField,
                                 map<long long,vector<gw_ring> > &shapes,string &err,const CGWCRS *crs,double *rawbox,
                                 map<long long,vector<int> > *parts,map<long long,string> *noArea)
{
  if (noArea!=NULL){noArea->clear();}
  ifstream IN(filename.c_str(),ios::binary);
  if (IN.fail()){err="cannot open GeoJSON file "+filename;return false;}
  stringstream ss; ss<<IN.rdbuf(); string txt=ss.str();
  if (txt.compare(0,3,"\xEF\xBB\xBF")==0){txt.erase(0,3);}   //(a UTF-8 byte-order mark, as some Windows programs write)
  JVal root; JParser P(txt.c_str(),txt.c_str()+txt.size());
  if (!P.Parse(root) || (root.t!=J_OBJ)){err="invalid JSON in "+filename+P.Where();return false;}
  if (!P.AtEnd()){err="invalid JSON in "+filename+": text after the end of the GeoJSON object"+P.Where();return false;}
  const JVal *feats=root.Get("features");
  if ((feats==NULL) || (feats->t!=J_ARR)){err=filename+" is not a GeoJSON FeatureCollection";return false;}
  JCoords C(crs);
  if (!CheckFileCRS(root,filename,crs,err,C.unknown)){return false;}
  int nfound=0,npart=0; map<long long,string> first; //(first: the first feature of an ID read without a ring)
  for (size_t f=0;f<feats->a.size();f++)
  {
    const JVal &F=feats->a[f];
    const JVal *props=F.Get("properties");
    const JVal *geom =F.Get("geometry");
    if ((geom==NULL) || (geom->t!=J_OBJ)){continue;}
    long long id=0;
    C.feature="feature "+to_string((unsigned long long)(f+1));
    const JVal *gt=geom->Get("type"),*co=geom->Get("coordinates");
    if ((gt==NULL) || (co==NULL)){continue;}
    if ((gt->s!="Polygon") && (gt->s!="MultiPolygon")){continue;} //(checked first: features of other types, and their IDs, are not read)
    if (idField!=""){
      if ((props==NULL) || (props->t!=J_OBJ)){err="GeoJSON "+C.feature+" in "+filename+" has no \"properties\" member, so no ID field '"+idField+"'";return false;}
      const JVal *jid=props->Get(idField);
      if (jid==NULL){err=IDMissing(F,idField,filename);return false;}
      if (!FeatureID(*jid,id)){err="GeoJSON file "+filename+", "+C.feature+": its ID field '"+idField+"' is "+IDRefused(*jid);return false;}
    }
    vector<gw_ring> &rings=shapes[id];
    vector<int> *pp=(parts!=NULL)?&(*parts)[id]:NULL;  //(the polygon each ring belongs to: a Polygon, or one part of a MultiPolygon)
    if (idField!=""){C.feature+=" ("+idField+" "+to_string(id)+")";}
    if (gt->s=="Polygon"){AddPolygon(*co,rings,C,pp,npart++,"its coordinates");}
    else if (C.Array(*co,"its coordinates")){for (size_t i=0;i<co->a.size();i++){AddPolygon(co->a[i],rings,C,pp,npart++,"polygon "+to_string((unsigned long long)(i+1)));}}
    if (!C.bad.empty()){err="GeoJSON file "+filename+", "+C.bad;return false;}
    if (rings.empty() && (first.find(id)==first.end())){
      size_t nr=0; //(the rings written: none is an empty geometry)
      if (gt->s=="Polygon"){nr=co->a.size();} else {for (size_t i=0;i<co->a.size();i++){nr+=(co->a[i].t==J_ARR)?co->a[i].a.size():1;}}
      first[id]=C.feature+((nr==0)?": it has no polygon (its coordinates are empty)":": its polygon has no area (its outline has fewer than 3 distinct vertices, lies on one line or encloses no net area)");
    }
    nfound++;
  }
  if (nfound==0){err="no Polygon/MultiPolygon features found in "+filename;return false;}  //(callers may then try lines)
  for (map<long long,string>::const_iterator it=first.begin();it!=first.end();it++){ //(an ID none of whose polygons has an area is left
    if (shapes[it->first].empty()){                                                 // out; the caller decides whether it needs it)
      shapes.erase(it->first); if (parts!=NULL){parts->erase(it->first);}
      if (noArea!=NULL){(*noArea)[it->first]="GeoJSON file "+filename+", "+it->second;}
    }
  }
  C.Box(rawbox);
  return CheckLonLat(C,filename,err);
}

//reads LineString/MultiLineString features; idField may be empty (ids set to -1)
bool GWGeom::ReadGeoJSONLines(const string &filename,const string &idField,
                              vector<gw_ring> &lines,vector<long long> &ids,string &err,const CGWCRS *crs,double *rawbox)
{
  ifstream IN(filename.c_str(),ios::binary);
  if (IN.fail()){err="cannot open GeoJSON file "+filename;return false;}
  stringstream ss; ss<<IN.rdbuf(); string txt=ss.str();
  if (txt.compare(0,3,"\xEF\xBB\xBF")==0){txt.erase(0,3);}   //(a UTF-8 byte-order mark, as some Windows programs write)
  JVal root; JParser P(txt.c_str(),txt.c_str()+txt.size());
  if (!P.Parse(root) || (root.t!=J_OBJ)){err="invalid JSON in "+filename+P.Where();return false;}
  if (!P.AtEnd()){err="invalid JSON in "+filename+": text after the end of the GeoJSON object"+P.Where();return false;}
  const JVal *feats=root.Get("features");
  if ((feats==NULL) || (feats->t!=J_ARR)){err=filename+" is not a GeoJSON FeatureCollection";return false;}
  JCoords C(crs);
  if (!CheckFileCRS(root,filename,crs,err,C.unknown)){return false;}
  for (size_t f=0;f<feats->a.size();f++)
  {
    const JVal &F=feats->a[f];
    const JVal *props=F.Get("properties"),*geom=F.Get("geometry");
    if ((geom==NULL) || (geom->t!=J_OBJ)){continue;}
    long long id=-1;
    C.feature="feature "+to_string((unsigned long long)(f+1));
    const JVal *gt=geom->Get("type"),*co=geom->Get("coordinates");
    if ((gt==NULL) || (co==NULL)){continue;}
    if ((gt->s!="LineString") && (gt->s!="MultiLineString")){continue;} //(features of other types, and their IDs, are not read)
    if ((idField!="") && ((props==NULL) || (props->t!=J_OBJ))){err="GeoJSON "+C.feature+" in "+filename+" has no \"properties\" member, so no ID field '"+idField+"'";return false;}
    if (idField!=""){
      const JVal *jid=props->Get(idField);
      if (jid==NULL){err=IDMissing(F,idField,filename);return false;}
      if (!FeatureID(*jid,id)){err="GeoJSON file "+filename+", "+C.feature+": its ID field '"+idField+"' is "+IDRefused(*jid);return false;}
      C.feature+=" ("+idField+" "+to_string(id)+")";
    }
    vector<const JVal*> parts;
    if      (gt->s=="LineString"     ){parts.push_back(co);}
    else if ((gt->s=="MultiLineString") && C.Array(*co,"its coordinates")){for (size_t i=0;i<co->a.size();i++){parts.push_back(&co->a[i]);}}
    for (size_t i=0;i<parts.size();i++){
      gw_ring L;
      if (!C.Array(*parts[i],(parts.size()>1)?"line "+to_string((unsigned long long)(i+1)):string("its coordinates"))){err="GeoJSON file "+filename+", "+C.bad;return false;}
      for (size_t n=0;n<parts[i]->a.size();n++){
        const JVal &pt=parts[i]->a[n]; double x,y;
        if (C.Pair(pt,((parts.size()>1)?"line "+to_string((unsigned long long)(i+1))+", ":string(""))+"vertex "+to_string((unsigned long long)(n+1)),x,y)){L.push_back(C.Point(x,y));}
        else if (!C.bad.empty()){err="GeoJSON file "+filename+", "+C.bad;return false;}
      }
      if (L.size()>=2){lines.push_back(L); ids.push_back(id);}
    }
  }
  if (lines.size()==0){err="no LineString features found in "+filename;return false;}
  C.Box(rawbox);
  return CheckLonLat(C,filename,err);
}

//Sutherland-Hodgman clipping of an arbitrary ring by an axis-aligned rectangle; returns signed area of result
namespace
{
  inline bool Inside(const gw_pt &p,int edge,double v){
    switch(edge){case 0:return p.x>=v; case 1:return p.x<=v; case 2:return p.y>=v; default:return p.y<=v;}
  }
  inline gw_pt Cross(const gw_pt &a,const gw_pt &b,int edge,double v){
    gw_pt q; double t;
    if (edge<2){t=(v-a.x)/(b.x-a.x); q.x=v; q.y=a.y+t*(b.y-a.y);}
    else       {t=(v-a.y)/(b.y-a.y); q.y=v; q.x=a.x+t*(b.x-a.x);}
    return q;
  }
}
double GWGeom::ClipRingAreaToRect(const gw_ring &r,const double x0,const double y0,const double x1,const double y1)
{
  gw_ring in=r,out;
  double lim[4]={x0,x1,y0,y1};
  for (int edge=0;edge<4;edge++)
  {
    out.clear(); size_t n=in.size();
    if (n==0){break;}
    for (size_t i=0;i<n;i++){
      const gw_pt &cur=in[i],&prev=in[(i+n-1)%n];
      bool ci=Inside(cur,edge,lim[edge]),pi=Inside(prev,edge,lim[edge]);
      if (ci){
        if (!pi){out.push_back(Cross(prev,cur,edge,lim[edge]));}
        out.push_back(cur);
      }
      else if (pi){out.push_back(Cross(prev,cur,edge,lim[edge]));}
    }
    in.swap(out);
  }
  if (in.size()<3){return 0.0;}
  return SignedArea(in);
}

//=====================================================================================================
// CGWCRS: projections of existing models (formulas: Karney 2011 for Transverse Mercator, Snyder 1987 for the conics)
//=====================================================================================================
namespace { const double CRS_PI=3.14159265358979323846,D2R=CRS_PI/180.0; }
CGWCRS::CGWCRS(){_type=0; _epsg=0; _code=""; _a=6378137.0; _f=1.0/298.257223563; _e2=_f*(2-_f); _e=sqrt(_e2); _lat0=_lon0=_lat1=_lat2=0; _k0=1; _fe=_fn=0;
  _n=_A=_xi0=_nc=_C=_F=_rho0=0; for (int i=0;i<5;i++){_b[i]=_al[i]=_d[i]=0;}}
double CGWCRS::M(const double p) const {return cos(p)/sqrt(1-_e2*sin(p)*sin(p));}
double CGWCRS::T(const double p) const {double s=sin(p); return tan(CRS_PI/4-p/2)/pow((1-_e*s)/(1+_e*s),_e/2);}
double CGWCRS::Q(const double p) const {double s=sin(p); return (1-_e2)*(s/(1-_e2*s*s)-1/(2*_e)*log((1-_e*s)/(1+_e*s)));}
bool CGWCRS::Parse(const string &spec0,string &err)
{
  string spec=spec0; for (size_t i=0;i<spec.size();i++){spec[i]=(char)toupper((unsigned char)spec[i]);}
  _epsg=0; _code="";
  istringstream is(spec),is0(spec0); vector<string> w,w0; string t; while (is>>t){w.push_back(t);} while (is0>>t){w0.push_back(t);}
  if (w.empty()){err="empty coordinate system"; return false;}
  string ell="WGS84";
  //every parameter is read strictly (the whole token a finite number) and the count is exact: a typo would otherwise
  //move or shear the model without a word
  const char *pn1[5]={"lat0","lon0","k0","FE","FN"},*pn2[6]={"lat1","lat2","lat0","lon0","FE","FN"};
  double par[6]={0,0,0,0,0,0};
  if ((w[0]=="TMERC") || (w[0]=="ALBERS") || (w[0]=="LCC")){
    const int np=(w[0]=="TMERC")?5:6;
    if ((w.size()<(size_t)np+1) || (w.size()>(size_t)np+2)){
      err=w0[0]+" takes "+((np==5)?"lat0 lon0 k0 FE FN":"lat1 lat2 lat0 lon0 FE FN")+" and an optional ellipsoid (WGS84 or GRS80): "+to_string((long long)w.size()-1)+" values are given";
      return false;
    }
    for (int i=0;i<np;i++){
      if (!gwutil::StrictNumber(w0[i+1].c_str(),par[i])){err=w0[0]+" parameter "+string((np==5)?pn1[i]:pn2[i])+": '"+w0[i+1]+"' is not a number"; return false;}
    }
    if ((np==5) && !(par[2]>0.0)){err="TMERC parameter k0 (the scale factor) must be positive: '"+w0[3]+"'"; return false;}
    for (int i=0;i<np;i++){
      const string nm=(np==5)?pn1[i]:pn2[i];
      if ((nm.compare(0,3,"lat")==0) && (fabs(par[i])>90.0)){err=w0[0]+" parameter "+nm+": "+w0[i+1]+" is not a latitude (-90 to 90 degrees)"; return false;}
      if ((nm=="lon0") && (fabs(par[i])>360.0)){err=w0[0]+" parameter lon0: "+w0[i+1]+" is not a longitude (degrees)"; return false;}
    }
  }
  if ((w[0].compare(0,5,"EPSG:")==0) || (w[0]=="ESRI:102001")){ //(102001 is an ESRI code; EPSG:102001 is accepted as well)
    if (w.size()>1){err="surplus value '"+w0[1]+"' after "+w0[0]+" (an EPSG code takes no parameters)"; return false;}
    const string num=w[0].substr(5);
    if (num.empty() || (num.find_first_not_of("0123456789")!=string::npos) || (num.size()>9)){err="'"+w0[0]+"' is not EPSG: followed by a whole-number code"; return false;}
    long code=atol(num.c_str());
    if ((code>32600) && (code<=32660)){_type=1; _lat0=0; _lon0=-183+6*(code-32600); _k0=0.9996; _fe=500000; _fn=0;}
    else if ((code>32700) && (code<=32760)){_type=1; _lat0=0; _lon0=-183+6*(code-32700); _k0=0.9996; _fe=500000; _fn=10000000;}
    else if ((code>=26901) && (code<=26923)){_type=1; _lat0=0; _lon0=-183+6*(code-26900); _k0=0.9996; _fe=500000; _fn=0; ell="GRS80";}
    else if (code==3005)  {_type=2; _lat1=50; _lat2=58.5; _lat0=45; _lon0=-126; _fe=1000000; _fn=0; ell="GRS80";}           //NAD83 / BC Albers
    else if (code==102001){_type=2; _lat1=50; _lat2=70; _lat0=40; _lon0=-96; _fe=0; _fn=0; ell="GRS80";}                    //Canada Albers (ESRI)
    else if (code==5070)  {_type=2; _lat1=29.5; _lat2=45.5; _lat0=23; _lon0=-96; _fe=0; _fn=0; ell="GRS80";}              //NAD83 / Conus Albers
    else if (code==3978)  {_type=3; _lat1=49; _lat2=77; _lat0=49; _lon0=-95; _fe=0; _fn=0; ell="GRS80";}                    //NAD83 / Canada Atlas Lambert
    else if (code==3347)  {_type=3; _lat1=49; _lat2=77; _lat0=63.390675; _lon0=-91.86666666666666; _fe=6200000; _fn=3000000; ell="GRS80";} //Statistics Canada Lambert
    else {err="EPSG code "+w[0].substr(5)+" is not built in; give the projection by its parameters (TMERC, ALBERS or LCC)"; return false;}
    if ((w[0].compare(0,5,"EPSG:")==0) && (code!=102001)){_epsg=code;}
    _code=(code==102001)?"ESRI:102001":w[0];
  }
  else if (w[0]=="TMERC"){_type=1; _lat0=par[0]; _lon0=par[1]; _k0=par[2]; _fe=par[3]; _fn=par[4]; if (w.size()>6){ell=w[6];}}
  else if ((w[0]=="ALBERS") || (w[0]=="LCC")){
    _type=(w[0]=="ALBERS")?2:3; _lat1=par[0]; _lat2=par[1]; _lat0=par[2]; _lon0=par[3]; _fe=par[4]; _fn=par[5]; if (w.size()>7){ell=w[7];}
  }
  else {err="coordinate system must be EPSG:code, TMERC lat0 lon0 k0 FE FN, ALBERS lat1 lat2 lat0 lon0 FE FN or LCC lat1 lat2 lat0 lon0 FE FN"; return false;}
  if (ell=="GRS80"){_f=1.0/298.257222101;} else if (ell=="WGS84"){_f=1.0/298.257223563;} else {err="ellipsoid must be WGS84 or GRS80, not '"+w0[w.size()-1]+"'"; return false;}
  //a conic projection needs a cone: standard parallels symmetric about the equator (lat1 = -lat2, both 0 included) give
  //none (as PROJ refuses them), and a Lambert conformal parallel at a pole none either
  const string pname=(_type==2)?"ALBERS":"LCC",pars=" (standard parallels lat1 "+gwutil::ToStr(_lat1)+", lat2 "+gwutil::ToStr(_lat2)+")";
  if ((_type>=2) && (fabs(_lat1+_lat2)<1e-10)){
    err=pname+pars+": the standard parallels lie symmetric about the equator, which defines no cone; give two parallels on the same side of the equator";
    return false;
  }
  if ((_type==3) && ((fabs(fabs(_lat1)-90.0)<1e-10) || (fabs(fabs(_lat2)-90.0)<1e-10))){
    err=pname+pars+": a standard parallel at a pole defines no Lambert conformal cone; give parallels between the poles";
    return false;
  }
  _a=6378137.0; _e2=_f*(2-_f); _e=sqrt(_e2);
  _lat0*=D2R; _lon0*=D2R; _lat1*=D2R; _lat2*=D2R;
  if (_type==1){
    double n=_f/(2-_f),n2=n*n,n3=n2*n,n4=n3*n; _n=n;
    _A=_a/(1+n)*(1+n2/4+n4/64);
    _al[1]=n/2-2*n2/3+5*n3/16+41*n4/180; _al[2]=13*n2/48-3*n3/5+557*n4/1440; _al[3]=61*n3/240-103*n4/140; _al[4]=49561*n4/161280;
    _b[1]=n/2-2*n2/3+37*n3/96-n4/360;    _b[2]=n2/48+n3/15-437*n4/1440;      _b[3]=17*n3/480-37*n4/840;    _b[4]=4397*n4/161280;
    _d[1]=2*n-2*n2/3-2*n3+116*n4/45;     _d[2]=7*n2/3-8*n3/5-227*n4/45;      _d[3]=56*n3/15-136*n4/35;     _d[4]=4279*n4/630;
    _xi0=0; if (_lat0!=0){double x0,y0; _type=1; _xi0=0; Forward(_lon0/D2R,_lat0/D2R,x0,y0); _xi0=(y0-_fn)/(_k0*_A);} //meridian arc to lat0
  }
  else if (_type==2){
    double m1=M(_lat1),m2=M(_lat2),q1=Q(_lat1),q2=Q(_lat2),q0=Q(_lat0);
    _nc=(fabs(_lat1-_lat2)>1e-12)?(m1*m1-m2*m2)/(q2-q1):sin(_lat1); _C=m1*m1+_nc*q1; _rho0=_a*sqrt(_C-_nc*q0)/_nc;
  }
  else {
    double m1=M(_lat1),m2=M(_lat2),t1=T(_lat1),t2=T(_lat2),t0=T(_lat0);
    _nc=(fabs(_lat1-_lat2)>1e-12)?(log(m1)-log(m2))/(log(t1)-log(t2)):sin(_lat1); _F=m1/(_nc*pow(t1,_nc));
    _rho0=(fabs(fabs(_lat0)-90.0*D2R)<1e-10)?0.0:_a*_F*pow(t0,_nc); //(an origin at a pole is the cone's apex, as PROJ takes it at either pole)
  }
  if ((_type>=2) && !(std::isfinite(_nc) && (_nc!=0.0) && std::isfinite(_rho0) && std::isfinite((_type==2)?_C:_F))){
    err=pname+pars+" and origin latitude "+gwutil::ToStr(_lat0/D2R)+": these parameters give no usable projection (its cone constants are not finite numbers)";
    return false;
  }
  return true;
}
void CGWCRS::Forward(const double lon,const double lat,double &x,double &y) const
{
  //(the longitude within +-180 degrees first, then its difference from the central meridian, as PROJ takes them (its
  // adjlon, which leaves an angle within pi+1e-12 as it is): the cones are not periodic in it, and a point 180 degrees
  // away keeps the side PROJ gives it)
  struct A{static double Adj(double v){if (fabs(v)<CRS_PI+1e-12){return v;} v+=CRS_PI; v-=2*CRS_PI*floor(v/(2*CRS_PI)); return v-CRS_PI;}};
  const double p=lat*D2R,l=A::Adj(A::Adj(lon*D2R)-_lon0);
  if (_type==1){
    double s=sin(p),c2=2*sqrt(_n)/(1+_n);
    double t=sinh(atanh(s)-c2*atanh(c2*s));
    double xp=atan2(t,cos(l)),ep=atanh(sin(l)/sqrt(1+t*t)),xi=xp,eta=ep;
    for (int j=1;j<=4;j++){xi+=_al[j]*sin(2*j*xp)*cosh(2*j*ep); eta+=_al[j]*cos(2*j*xp)*sinh(2*j*ep);}
    x=_fe+_k0*_A*eta; y=_fn+_k0*_A*(xi-_xi0);
  }
  else if (_type==2){double rho=_a*sqrt(_C-_nc*Q(p))/_nc,th=_nc*l; x=_fe+rho*sin(th); y=_fn+_rho0-rho*cos(th);}
  else if (_type==3){
    double rho=_a*_F*pow(T(p),_nc),th=_nc*l;
    if ((_nc<0) && (fabs(p+CRS_PI/2)<1e-12)){rho=0.0;} //(the apex of a southern cone: T(-pi/2) is no infinity in floating point)
    x=_fe+rho*sin(th); y=_fn+_rho0-rho*cos(th);
  }
  else {x=lon; y=lat;}
}
void CGWCRS::Inverse(const double x,const double y,double &lon,double &lat) const
{
  if (_type==1){
    double xi=(y-_fn)/(_k0*_A)+_xi0,eta=(x-_fe)/(_k0*_A),xp=xi,ep=eta;
    for (int j=1;j<=4;j++){xp-=_b[j]*sin(2*j*xi)*cosh(2*j*eta); ep-=_b[j]*cos(2*j*xi)*sinh(2*j*eta);}
    double chi=asin(sin(xp)/cosh(ep)),p=chi;
    for (int j=1;j<=4;j++){p+=_d[j]*sin(2*j*chi);}
    lat=p/D2R; lon=(_lon0+atan2(sinh(ep),cos(xp)))/D2R;
    return;
  }
  if ((_type==2) || (_type==3)){
    double dx=x-_fe,dy=_rho0-(y-_fn),sg=(_nc<0)?-1.0:1.0;
    double rho=sg*sqrt(dx*dx+dy*dy),th=atan2(sg*dx,sg*dy),p;
    if (_type==2){
      double q=(_C-rho*rho*_nc*_nc/(_a*_a))/_nc; p=asin(max(-1.0,min(1.0,q/2)));
      for (int it=0;it<20;it++){
        double s=sin(p),dp=(1-_e2*s*s)*(1-_e2*s*s)/(2*cos(p))*(q/(1-_e2)-s/(1-_e2*s*s)+1/(2*_e)*log((1-_e*s)/(1+_e*s)));
        p+=dp; if (fabs(dp)<1e-14){break;}
      }
    }
    else {
      double t=pow(rho/(_a*_F),1.0/_nc); p=CRS_PI/2-2*atan(t);
      for (int it=0;it<20;it++){double s=sin(p),pn=CRS_PI/2-2*atan(t*pow((1-_e*s)/(1+_e*s),_e/2)); if (fabs(pn-p)<1e-14){p=pn; break;} p=pn;}
    }
    lat=p/D2R; lon=(_lon0+th/_nc)/D2R;
    return;
  }
  lon=x; lat=y;
}
bool CGWCRS::Same(const CGWCRS &o) const
{
  const double tr=1e-10,tm=1e-4; //[rad], [m]
  if ((_type!=o._type) || (fabs(_f-o._f)>1e-15)){return false;}
  if ((fabs(_lat0-o._lat0)>tr) || (fabs(_lon0-o._lon0)>tr) || (fabs(_fe-o._fe)>tm) || (fabs(_fn-o._fn)>tm)){return false;}
  if (_type==1){return (fabs(_k0-o._k0)<1e-12);}
  return (fabs(_lat1-o._lat1)<tr) && (fabs(_lat2-o._lat2)<tr);
}
string CGWCRS::Describe() const
{
  ostringstream s; s<<setprecision(10);
  const char *nm[4]={"none","Transverse Mercator","Albers equal-area conic","Lambert conformal conic"};
  s<<nm[_type]<<" (lon0 "<<_lon0/D2R<<", lat0 "<<_lat0/D2R;
  if (_type==1){s<<", k0 "<<_k0;} else if (_type>1){s<<", standard parallels "<<_lat1/D2R<<" and "<<_lat2/D2R;}
  s<<", FE "<<_fe<<", FN "<<_fn<<", "<<((fabs(_f-1.0/298.257222101)<1e-15)?"GRS80":"WGS84")<<")";
  return s.str();
}
