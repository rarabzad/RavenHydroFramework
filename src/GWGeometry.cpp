/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: spatial utilities
----------------------------------------------------------------*/
#include "GWGeometry.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <algorithm>

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
}

//=================================================================
// Raster
//=================================================================
CGWRaster::CGWRaster(){_ncols=_nrows=0;_xll=_yll=0;_cs=1;_nodata=-9999;}
bool CGWRaster::Read(const string &filename,string &err)
{
  ifstream IN(filename.c_str());
  if (IN.fail()){err="cannot open raster file "+filename;return false;}
  bool center=false; string key; double val;
  _nodata=-9999;
  for (int i=0;i<6;i++)
  {
    streampos pos=IN.tellg();
    if (!(IN>>key)){err="bad raster header in "+filename;return false;}
    string lk=key; std::transform(lk.begin(),lk.end(),lk.begin(),::tolower);
    if      (lk=="ncols"    ){IN>>_ncols;}
    else if (lk=="nrows"    ){IN>>_nrows;}
    else if (lk=="xllcorner"){IN>>_xll;}
    else if (lk=="yllcorner"){IN>>_yll;}
    else if (lk=="xllcenter"){IN>>_xll;center=true;}
    else if (lk=="yllcenter"){IN>>_yll;center=true;}
    else if (lk=="cellsize" ){IN>>_cs;}
    else if (lk=="nodata_value"){IN>>_nodata;}
    else {IN.seekg(pos);break;} //no NODATA line
  }
  if ((_ncols<=0) || (_nrows<=0) || (_cs<=0)){err="invalid raster header in "+filename;return false;}
  if (center){_xll-=0.5*_cs;_yll-=0.5*_cs;}
  _v.resize((size_t)_ncols*_nrows);
  for (size_t i=0;i<_v.size();i++){
    if (!(IN>>val)){err="raster "+filename+" has fewer values than ncols*nrows";_ncols=0;return false;}
    _v[i]=val;
  }
  return true;
}
bool CGWRaster::Valid(int r,int c) const
{
  if ((r<0) || (c<0) || (r>=_nrows) || (c>=_ncols)){return false;}
  return (fabs(_v[(size_t)r*_ncols+c]-_nodata)>1e-6);
}
bool CGWRaster::Sample(const double lon,const double lat,double &val) const
{
  if (_ncols==0){return false;}
  double ytop=_yll+_nrows*_cs;
  double fc=(lon-_xll)/_cs-0.5;
  double fr=(ytop-lat)/_cs-0.5;
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
    const char *p,*e;
    void ws(){while((p<e) && isspace((unsigned char)*p)){p++;}}
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
            case 'u':if (e-p<5){return false;} out+='?';p+=4;break; //(\uXXXX: non-ASCII names are not needed)
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
    JParser(const char *b,const char *en):p(b),e(en){}
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
      else {
        char *endp; v.t=J_NUM; v.n=strtod(p,&endp);
        if (endp==p){return false;}
        p=endp; return true;
      }
    }
  };
  bool RingFromJSON(const JVal &jr,gw_ring &r)
  {
    r.clear();
    if (jr.t!=J_ARR){return false;}
    for (size_t i=0;i<jr.a.size();i++){
      const JVal &pt=jr.a[i];
      if ((pt.t!=J_ARR) || (pt.a.size()<2)){return false;}
      gw_pt q; q.x=pt.a[0].n; q.y=pt.a[1].n; r.push_back(q);
    }
    if ((r.size()>1) && (r.front().x==r.back().x) && (r.front().y==r.back().y)){r.pop_back();}
    return (r.size()>=3);
  }
  void AddPolygon(const JVal &jpoly,vector<gw_ring> &rings)
  {
    for (size_t i=0;i<jpoly.a.size();i++){
      gw_ring r;
      if (!RingFromJSON(jpoly.a[i],r)){continue;}
      double A=GWGeom::SignedArea(r);
      bool outer=(i==0);
      if ((outer && (A<0)) || (!outer && (A>0))){std::reverse(r.begin(),r.end());}
      rings.push_back(r);
    }
  }
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
                                 map<long long,vector<gw_ring> > &shapes,string &err)
{
  ifstream IN(filename.c_str(),ios::binary);
  if (IN.fail()){err="cannot open GeoJSON file "+filename;return false;}
  stringstream ss; ss<<IN.rdbuf(); string txt=ss.str();
  JVal root; JParser P(txt.c_str(),txt.c_str()+txt.size());
  if (!P.Parse(root) || (root.t!=J_OBJ)){err="invalid JSON in "+filename;return false;}
  const JVal *feats=root.Get("features");
  if ((feats==NULL) || (feats->t!=J_ARR)){err=filename+" is not a GeoJSON FeatureCollection";return false;}
  int nfound=0;
  for (size_t f=0;f<feats->a.size();f++)
  {
    const JVal &F=feats->a[f];
    const JVal *props=F.Get("properties");
    const JVal *geom =F.Get("geometry");
    if ((props==NULL) || (geom==NULL) || (geom->t!=J_OBJ)){continue;}
    long long id=0;
    if (idField!=""){
      const JVal *jid=props->Get(idField);
      if (jid==NULL){err="GeoJSON feature missing ID field '"+idField+"' in "+filename;return false;}
      if (jid->t==J_NUM){id=(long long)llround(jid->n);}
      else if (jid->t==J_STR){id=atoll(jid->s.c_str());}
      else {continue;}
    }
    const JVal *gt=geom->Get("type"),*co=geom->Get("coordinates");
    if ((gt==NULL) || (co==NULL)){continue;}
    if ((gt->s!="Polygon") && (gt->s!="MultiPolygon")){continue;} //(checked first: a line with the same ID adds no empty entry)
    vector<gw_ring> &rings=shapes[id];
    if (gt->s=="Polygon"){AddPolygon(*co,rings);}
    else {for (size_t i=0;i<co->a.size();i++){AddPolygon(co->a[i],rings);}}
    nfound++;
  }
  if (nfound==0){err="no Polygon/MultiPolygon features found in "+filename;return false;}  //(callers may then try lines)
  return true;
}

//reads LineString/MultiLineString features; idField may be empty (ids set to -1)
bool GWGeom::ReadGeoJSONLines(const string &filename,const string &idField,
                              vector<gw_ring> &lines,vector<long long> &ids,string &err)
{
  ifstream IN(filename.c_str(),ios::binary);
  if (IN.fail()){err="cannot open GeoJSON file "+filename;return false;}
  stringstream ss; ss<<IN.rdbuf(); string txt=ss.str();
  JVal root; JParser P(txt.c_str(),txt.c_str()+txt.size());
  if (!P.Parse(root) || (root.t!=J_OBJ)){err="invalid JSON in "+filename;return false;}
  const JVal *feats=root.Get("features");
  if ((feats==NULL) || (feats->t!=J_ARR)){err=filename+" is not a GeoJSON FeatureCollection";return false;}
  for (size_t f=0;f<feats->a.size();f++)
  {
    const JVal &F=feats->a[f];
    const JVal *props=F.Get("properties"),*geom=F.Get("geometry");
    if ((geom==NULL) || (geom->t!=J_OBJ)){continue;}
    long long id=-1;
    if ((idField!="") && (props!=NULL)){
      const JVal *jid=props->Get(idField);
      if (jid==NULL){err="GeoJSON feature missing ID field '"+idField+"' in "+filename;return false;}
      id=(jid->t==J_NUM)?(long long)llround(jid->n):atoll(jid->s.c_str());
    }
    const JVal *gt=geom->Get("type"),*co=geom->Get("coordinates");
    if ((gt==NULL) || (co==NULL)){continue;}
    vector<const JVal*> parts;
    if      (gt->s=="LineString"     ){parts.push_back(co);}
    else if (gt->s=="MultiLineString"){for (size_t i=0;i<co->a.size();i++){parts.push_back(&co->a[i]);}}
    for (size_t i=0;i<parts.size();i++){
      gw_ring L;
      for (size_t n=0;n<parts[i]->a.size();n++){
        const JVal &pt=parts[i]->a[n];
        if ((pt.t==J_ARR) && (pt.a.size()>=2)){gw_pt q; q.x=pt.a[0].n; q.y=pt.a[1].n; L.push_back(q);}
      }
      if (L.size()>=2){lines.push_back(L); ids.push_back(id);}
    }
  }
  if (lines.size()==0){err="no LineString features found in "+filename;return false;}
  return true;
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
CGWCRS::CGWCRS(){_type=0; _a=6378137.0; _f=1.0/298.257223563; _e2=_f*(2-_f); _e=sqrt(_e2); _lat0=_lon0=_lat1=_lat2=0; _k0=1; _fe=_fn=0;
  _n=_A=_xi0=_nc=_C=_F=_rho0=0; for (int i=0;i<5;i++){_b[i]=_al[i]=_d[i]=0;}}
double CGWCRS::M(const double p) const {return cos(p)/sqrt(1-_e2*sin(p)*sin(p));}
double CGWCRS::T(const double p) const {double s=sin(p); return tan(CRS_PI/4-p/2)/pow((1-_e*s)/(1+_e*s),_e/2);}
double CGWCRS::Q(const double p) const {double s=sin(p); return (1-_e2)*(s/(1-_e2*s*s)-1/(2*_e)*log((1-_e*s)/(1+_e*s)));}
bool CGWCRS::Parse(const string &spec0,string &err)
{
  string spec=spec0; for (size_t i=0;i<spec.size();i++){spec[i]=(char)toupper((unsigned char)spec[i]);}
  istringstream is(spec); vector<string> w; string t; while (is>>t){w.push_back(t);}
  if (w.empty()){err="empty coordinate system"; return false;}
  string ell="WGS84";
  if ((w[0].compare(0,5,"EPSG:")==0) || (w[0]=="ESRI:102001")){ //(102001 is an ESRI code; EPSG:102001 is accepted as well)
    long code=atol(w[0].substr(5).c_str());
    if ((code>32600) && (code<=32660)){_type=1; _lat0=0; _lon0=-183+6*(code-32600); _k0=0.9996; _fe=500000; _fn=0;}
    else if ((code>32700) && (code<=32760)){_type=1; _lat0=0; _lon0=-183+6*(code-32700); _k0=0.9996; _fe=500000; _fn=10000000;}
    else if ((code>=26901) && (code<=26923)){_type=1; _lat0=0; _lon0=-183+6*(code-26900); _k0=0.9996; _fe=500000; _fn=0; ell="GRS80";}
    else if (code==3005)  {_type=2; _lat1=50; _lat2=58.5; _lat0=45; _lon0=-126; _fe=1000000; _fn=0; ell="GRS80";}           //NAD83 / BC Albers
    else if (code==102001){_type=2; _lat1=50; _lat2=70; _lat0=40; _lon0=-96; _fe=0; _fn=0; ell="GRS80";}                    //Canada Albers (ESRI)
    else if (code==5070)  {_type=2; _lat1=29.5; _lat2=45.5; _lat0=23; _lon0=-96; _fe=0; _fn=0; ell="GRS80";}              //NAD83 / Conus Albers
    else if (code==3978)  {_type=3; _lat1=49; _lat2=77; _lat0=49; _lon0=-95; _fe=0; _fn=0; ell="GRS80";}                    //NAD83 / Canada Atlas Lambert
    else if (code==3347)  {_type=3; _lat1=49; _lat2=77; _lat0=63.390675; _lon0=-91.86666666666666; _fe=6200000; _fn=3000000; ell="GRS80";} //Statistics Canada Lambert
    else {err="EPSG code "+w[0].substr(5)+" is not built in; give the projection by its parameters (TMERC, ALBERS or LCC)"; return false;}
  }
  else if ((w[0]=="TMERC") && (w.size()>=6)){_type=1; _lat0=atof(w[1].c_str()); _lon0=atof(w[2].c_str()); _k0=atof(w[3].c_str()); _fe=atof(w[4].c_str()); _fn=atof(w[5].c_str()); if (w.size()>6){ell=w[6];}}
  else if (((w[0]=="ALBERS") || (w[0]=="LCC")) && (w.size()>=7)){
    _type=(w[0]=="ALBERS")?2:3; _lat1=atof(w[1].c_str()); _lat2=atof(w[2].c_str()); _lat0=atof(w[3].c_str()); _lon0=atof(w[4].c_str());
    _fe=atof(w[5].c_str()); _fn=atof(w[6].c_str()); if (w.size()>7){ell=w[7];}
  }
  else {err="coordinate system must be EPSG:code, TMERC lat0 lon0 k0 FE FN, ALBERS lat1 lat2 lat0 lon0 FE FN or LCC lat1 lat2 lat0 lon0 FE FN"; return false;}
  if (ell=="GRS80"){_f=1.0/298.257222101;} else if (ell=="WGS84"){_f=1.0/298.257223563;} else {err="ellipsoid must be WGS84 or GRS80"; return false;}
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
    _nc=(fabs(_lat1-_lat2)>1e-12)?(log(m1)-log(m2))/(log(t1)-log(t2)):sin(_lat1); _F=m1/(_nc*pow(t1,_nc)); _rho0=_a*_F*pow(t0,_nc);
  }
  return true;
}
void CGWCRS::Forward(const double lon,const double lat,double &x,double &y) const
{
  double p=lat*D2R,l=lon*D2R-_lon0;
  if (_type==1){
    double s=sin(p),c2=2*sqrt(_n)/(1+_n);
    double t=sinh(atanh(s)-c2*atanh(c2*s));
    double xp=atan2(t,cos(l)),ep=atanh(sin(l)/sqrt(1+t*t)),xi=xp,eta=ep;
    for (int j=1;j<=4;j++){xi+=_al[j]*sin(2*j*xp)*cosh(2*j*ep); eta+=_al[j]*cos(2*j*xp)*sinh(2*j*ep);}
    x=_fe+_k0*_A*eta; y=_fn+_k0*_A*(xi-_xi0);
  }
  else if (_type==2){double rho=_a*sqrt(_C-_nc*Q(p))/_nc,th=_nc*l; x=_fe+rho*sin(th); y=_fn+_rho0-rho*cos(th);}
  else if (_type==3){double rho=_a*_F*pow(T(p),_nc),th=_nc*l; x=_fe+rho*sin(th); y=_fn+_rho0-rho*cos(th);}
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
string CGWCRS::Describe() const
{
  ostringstream s; s<<setprecision(10);
  const char *nm[4]={"none","Transverse Mercator","Albers equal-area conic","Lambert conformal conic"};
  s<<nm[_type]<<" (lon0 "<<_lon0/D2R<<", lat0 "<<_lat0/D2R;
  if (_type==1){s<<", k0 "<<_k0;} else if (_type>1){s<<", standard parallels "<<_lat1/D2R<<" and "<<_lat2/D2R;}
  s<<", FE "<<_fe<<", FN "<<_fn<<", "<<((fabs(_f-1.0/298.257222101)<1e-15)?"GRS80":"WGS84")<<")";
  return s.str();
}
