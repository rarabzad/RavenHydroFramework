/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Class CGroundwaterModel: Raven-native groundwater subsystem on MODFLOW 6
----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include <set>
//optional stage timing for large grids: set RAVEN_GW_PROFILE=1
static double gw_prof_t0=-1;
#define GWPROF(name) { if (getenv("RAVEN_GW_PROFILE")!=NULL){ double tn=double(clock())/CLOCKS_PER_SEC; if (gw_prof_t0<0){gw_prof_t0=tn;} cerr<<"[groundwater] "<<name<<" done at "<<tn<<" s"<<endl; } }
#include "Model.h"
#include "HydroProcessABC.h"
#include "TimeSeries.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#include <ctime>
#include <cstdint>
#ifdef netcdf
  #include <netcdf.h>
#endif
#include "GWUtil.h"
using namespace gwutil;

namespace
{
  const double GW_MIN_LAYER_THICK=0.5;   ///< [m] thinner layers become vertical pass-through cells
  const double GW_PASS_THICK     =0.01;  ///< [m] nominal thickness of pass-through cells

}

//////////////////////////////////////////////////////////////////
/// \brief constructor: creates an empty groundwater model; nothing happens unless :GroundwaterModel is used
//
CGroundwaterModel::CGroundwaterModel(CModel *pModel)
{
  _pModel=pModel;
  _hruIDField="HRU_ID";
  _cellSize=1000.0; _minCoverage=0.25; _headSaveFreq=30; _dvclose=0.00001; _subDaily=false; _lastRelaxed=false;
  _nrow=_ncol=_nlay=_ncpl=0; _x0=_y0=0;
  _gridType=GWGRID_REGULAR; _gridAngle=0.0; _qtMaxLevel=-1; _gridFileID="CELL_ID"; _flowCorr="AUTO"; _gridTypeGiven=false;
  _connectByElevation=false; _disType=0; _meshFidelity=1.0;
  _mname="GWF"; _hf=1.0; _qf=1.0; _rf=1.0;
  _extGridID="CELL_ID"; _streamMap="AUX"; _streamAux="SUBBASIN"; _extSeepage=true; _extTf=_extLf=1.0; _extDvclose=-1.0; _extAreaCheck=0.0;
  _Qother=0.0; _extAux=-1; _extNaux=0; _extOffset=0.0; _extDrop=0; _extNModels=0; _extMoverDropped=0; _extSteady=false; _extNoSto=false; _mf6Version=""; _stoName="STO"; _disName="DIS";
  _iGW=_iReturnSV=DOESNT_EXIST;
  _maxiter=500;
  _X=_rchBound=_rchHcof=_rchRhs=_drnHcof=_drnRhs=_strgss=_strgsy=NULL;
  _rivStage=_rivCond=_rivRbot=_rivHcof=_rivRhs=NULL; _rivNcol=0;
  _evtRate=_evtHcof=_evtRhs=NULL; _evtNcol=0; _Qevt=0.0;
  _welQ=_welHcof=_welRhs=_ghbHead=_ghbCond=_ghbHcof=_ghbRhs=_chdHead=_chdSim=NULL;
  _welNcol=_ghbNcol=_chdNcol=0; _Qwel_spec=_Qwel=_Qghb=_Qchd=0.0;
  _steadyInit=false; _doSteady=false; _ssRecharge=0.0; _initMode="FROM_PROFILE"; _initValue=0.0;
  _hotDims[0]=_hotDims[1]=_hotDims[2]=0; _active=false;
  _rchStoreTot=0.0; _QrchCHD=0.0; _rivBedDepth=-1.0;
  _swxHead=_swxCond=_swxRbot=_swxHcof=_swxRhs=NULL; _swxNcol=0; _Qswx_in=_Qswx_out=_swxAppliedVol=0.0;
  _resExchange=false; _resQ=_resHcof=_resRhs=NULL; _resNcol=0; _Qres_spec=_Qres=_resPendingVol=_resShortVol=0.0;
  _cacheHit=false; _geomSeconds=0.0; _writeNC=false; _ncid=-1; _ncTime=_ncHead=-1; _ncRec=0; _nstp=0; _stepCount=0; _jday0=0; _jyear0=0; _calendar=0;
  _Qriv_in=_Qriv_out=_riverAppliedVol=_riverDeficitVol=0.0; _budgetLineOpen=false;
  _ravenRchVol=_Qrch=_Qseep_tot=_Qsto=_Qerr=_returnVol=0.0;
  _nNonConverged=0;_lastIter=0;_lastConverged=true;
}
CGroundwaterModel::~CGroundwaterModel()
{
  CloseOutputs();
  for (map<long long,CTimeSeries*>::iterator it=_wellTS.begin();it!=_wellTS.end();it++){delete it->second;}
  for (map<string,CTimeSeries*>::iterator it=_headTS.begin();it!=_headTS.end();it++){delete it->second;}
}
//////////////////////////////////////////////////////////////////
/// \brief clears everything built by Initialize(), so the model can be rebuilt (ensemble members)
//
void CGroundwaterModel::ResetState()
{
  CloseOutputs();
  _aqHRUs.clear(); _hruProfile.clear(); _aLinks.clear(); _hruLinks.clear(); _hruLinkedArea.clear(); _hruCells.clear();
  _bndCol.clear(); _evtCol.clear(); _colDomHRU.clear();
  _rivNode.clear(); _rivCol.clear(); _rivSB.clear(); _rivLen.clear(); _rivBed.clear();
  _welNode.clear(); _welWell.clear(); _welFrac.clear(); _ghbNode.clear(); _ghbB.clear(); _chdNode.clear(); _chdB.clear();
  _swxCol.clear(); _swxX.clear(); _swxLinks.clear(); _swxVl.clear(); _swxAppliedVol=0.0;
  _resP.clear(); _resNode.clear(); _resR.clear(); _resW.clear();
  _X=NULL; _active=false; _stepCount=0; _nNonConverged=0; _budgetLineOpen=false;
  _ravenRchVol=_Qrch=_Qseep_tot=_Qsto=_Qerr=_returnVol=_riverAppliedVol=_riverDeficitVol=0.0;
}
void CGroundwaterModel::Reinitialize(const optStruct &Options)
{
  ResetState();
  Initialize(Options);
}
//////////////////////////////////////////////////////////////////
/// \brief date label (yyyy-mm-dd) at model time t, used to stamp outputs at the end of each step
//
string CGroundwaterModel::EndDate(const double &t) const
{
  time_struct te; JulianConvert(t,_jday0,_jyear0,_calendar,te);
  if (!_subDaily){return te.date_string;}
  int minutes=(int)floor((te.julian_day-floor(te.julian_day))*1440.0+0.5); //time of day for sub-daily steps
  if (minutes>=1440){minutes=0;}
  char buf[16]; snprintf(buf,sizeof(buf)," %02d:%02d",minutes/60,minutes%60);
  return te.date_string+string(buf);
}
//////////////////////////////////////////////////////////////////
/// \brief mean of a raster over a grid column (5x5 sample points), rather than the value at its centre
//
double CGroundwaterModel::MeanRaster(const CGWRaster &R,const int ic,bool &ok) const
{
  double sum=0; int n=0; vector<gw_pt> pts; _grid.SamplePoints(ic,5,pts);
  for (size_t k=0;k<pts.size();k++){
    double X,Y,lon,lat,v; _grid.ToWorld(pts[k].x,pts[k].y,X,Y);
    _proj.Inverse(X,Y,lon,lat);
    if (R.Sample(lon,lat,v)){sum+=v;n++;}
  }
  ok=(n>0); return ok?sum/n:0.0;
}
//cell address in the MODFLOW input files: DIS layer row column, DISV layer cell, DISU node
string CGroundwaterModel::CellID(const int l,const int ic) const
{
  ostringstream s;
  if      (_disType==0){s<<l+1<<" "<<ic/_ncol+1<<" "<<ic%_ncol+1;}
  else if (_disType==1){s<<l+1<<" "<<ic+1;}
  else                 {s<<l*_ncpl+ic+1;}
  return s.str();
}
//pieces of a straight segment (projected coordinates) inside grid cells, as parameter ranges [t0,t1] of the segment;
//a stretch lying exactly on a shared cell edge is kept once, in an active cell if there is one
void CGroundwaterModel::SegmentPieces(const double xa,const double ya,const double xb,const double yb,vector<gw_piece> &out) const
{
  out.clear();
  gw_pt A,B; _grid.ToLocal(xa,ya,A.x,A.y); _grid.ToLocal(xb,yb,B.x,B.y);
  vector<int> C; _grid.Candidates(min(A.x,B.x),min(A.y,B.y),max(A.x,B.x),max(A.y,B.y),C);
  vector<gw_piece> all; vector<pair<double,double> > iv;
  for (size_t k=0;k<C.size();k++){
    _grid.SegmentIntervals(A,B,C[k],iv);
    for (size_t j=0;j<iv.size();j++){if (iv[j].second-iv[j].first>1e-12){gw_piece p; p.ic=C[k]; p.t0=iv[j].first; p.t1=iv[j].second; all.push_back(p);}}
  }
  struct byStart { bool operator()(const gw_piece &a,const gw_piece &b) const {return (a.t0<b.t0) || ((a.t0==b.t0) && (a.ic<b.ic));} };
  sort(all.begin(),all.end(),byStart());
  for (size_t k=0;k<all.size();k++){
    if (!out.empty() && (all[k].t0<out.back().t1-1e-9)){ //overlaps the previous piece: segment along a shared edge
      if ((_colProfile[out.back().ic]==DOESNT_EXIST) && (_colProfile[all[k].ic]!=DOESNT_EXIST)){out.back()=all[k];}
      continue;
    }
    out.push_back(all[k]);
  }
}
//active cell at a point of a line: a point lying exactly on a cell edge belongs to either cell, so when the cell
//found there is missing or inactive the point is moved a hair to each side of the line
int CGroundwaterModel::LocateActiveOnLine(const double X,const double Y,const double dx,const double dy) const
{
  int ic=LocateWorld(X,Y);
  if ((ic>=0) && (_colProfile[ic]!=DOESNT_EXIST)){return ic;}
  double L=sqrt(dx*dx+dy*dy); if (L<=0){return -1;}
  double e=1e-6*(_grid.uniform?_cellSize:_grid.MinCellSize()),nx=-dy/L*e,ny=dx/L*e;
  for (int s=-1;s<=1;s+=2){
    int jc=LocateWorld(X+s*nx,Y+s*ny);
    if ((jc>=0) && (_colProfile[jc]!=DOESNT_EXIST)){return jc;}
  }
  return -1;
}
//nearest active cell within one cell size of a point (0 if the point is inside it); -1 if none
int CGroundwaterModel::NearestActiveCell(const double X,const double Y) const
{
  double x,y; _grid.ToLocal(X,Y,x,y);
  double R=_grid.uniform?_cellSize:_grid.MeanCellSize();
  vector<int> C; _grid.Candidates(x-R,y-R,x+R,y+R,C);
  int best=-1; double bd=R;
  for (size_t t=0;t<C.size();t++){
    int ic=C[t]; if (_colProfile[ic]==DOESNT_EXIST){continue;}
    const gw_ring &P=_grid.poly[ic]; double d=1e300;
    if (GWGeom::PointInRing(P,x,y)){d=0.0;}
    else {
      for (size_t k=0;k<P.size();k++){
        const gw_pt &a=P[k],&b=P[(k+1)%P.size()];
        double ex=b.x-a.x,ey=b.y-a.y,L2=ex*ex+ey*ey,t2=(L2>0)?((x-a.x)*ex+(y-a.y)*ey)/L2:0.0;
        t2=max(0.0,min(1.0,t2)); d=min(d,hypot(x-(a.x+t2*ex),y-(a.y+t2*ey)));
      }
    }
    if (d<bd){bd=d; best=ic;}
  }
  return best;
}
int CGroundwaterModel::LocateWorld(const double X,const double Y) const
{
  double x,y; _grid.ToLocal(X,Y,x,y); return _grid.Locate(x,y);
}
void CGroundwaterModel::CellCentreWorld(const int ic,double &X,double &Y) const
{
  if (_grid.uniform && (_grid.angle==0.0)){X=_x0+(ic%_ncol+0.5)*_cellSize; Y=_y0+(_nrow-1-ic/_ncol+0.5)*_cellSize; return;} //original arithmetic
  //(a rotated grid is uniform in its own frame only: its centres must be rotated back to the map)
  _grid.ToWorld(_grid.cx[ic],_grid.cy[ic],X,Y);
}
//grid description for the linkage-cache key
string CGroundwaterModel::GridSpec() const
{
  ostringstream s; s<<setprecision(12)<<"grid:"<<_gridType<<":"<<_cellSize<<":"<<_gridAngle<<":"<<_qtMaxLevel<<":"<<_gridFile<<":"<<_gridFileID;
  for (size_t i=0;i<_refKind.size();i++){s<<":"<<_refKind[i]<<","<<_refFile[i]<<","<<_refSize[i];}
  for (size_t i=0;i<_refKind.size();i++){
    string f=(_refKind[i]=="RIVERS")?_riverFile:_refFile[i];
    if (f!=""){ifstream G(f.c_str(),ios::binary); uint64_t hh=1469598103934665603ULL; char b[65536];
      while (G.read(b,sizeof(b)) || G.gcount()>0){for (streamsize k=0;k<G.gcount();k++){hh^=(unsigned char)b[k]; hh*=1099511628211ULL;}}
      s<<":"<<hex<<hh<<dec;}
    if (_refKind[i]=="WELLS"){
      for (size_t w=0;w<_aWells.size();w++){s<<":"<<_aWells[w].lat<<","<<_aWells[w].lon;}
      for (size_t w=0;w<_aObsWells.size();w++){s<<":o"<<_aObsWells[w].lat<<","<<_aObsWells[w].lon;} //(refinement uses them too)
    }
    if (_refKind[i]=="LAKES"){for (int p=0;p<_pModel->GetNumSubBasins();p++){const CReservoir *pR=_pModel->GetSubBasin(p)->GetReservoir(); if ((pR!=NULL) && (pR->GetHRUIndex()>=0)){s<<":lake"<<_pModel->GetHydroUnit(pR->GetHRUIndex())->GetHRUID();}}}
  }
  if (_extGeomKey!=""){s<<":"<<_extGeomKey;}
  if (_gridFile!=""){ifstream G(_gridFile.c_str(),ios::binary); uint64_t hh=1469598103934665603ULL; char b[65536];
    while (G.read(b,sizeof(b)) || G.gcount()>0){for (streamsize k=0;k<G.gcount();k++){hh^=(unsigned char)b[k]; hh*=1099511628211ULL;}}
    s<<":"<<hex<<hh<<dec;}
  return s.str();
}
//refinement features in the grid's local frame
void CGroundwaterModel::BuildRefinement(vector<gw_refine> &ref) const
{
  ref.clear(); string err;
  for (size_t i=0;i<_refKind.size();i++){
    gw_refine f; f.size=_refSize[i];
    if (_refKind[i]=="WELLS"){
      for (size_t w=0;w<_aWells.size();w++){double X,Y; gw_pt p; _proj.Forward(_aWells[w].lon,_aWells[w].lat,X,Y); _grid.ToLocal(X,Y,p.x,p.y); f.kind=0; f.geom.assign(1,p); ref.push_back(f);}
      for (size_t w=0;w<_aObsWells.size();w++){double X,Y; gw_pt p; _proj.Forward(_aObsWells[w].lon,_aObsWells[w].lat,X,Y); _grid.ToLocal(X,Y,p.x,p.y); f.kind=0; f.geom.assign(1,p); ref.push_back(f);}
    }
    else if ((_refKind[i]=="RIVERS") || (_refKind[i]=="LINES")){
      string file=(_refKind[i]=="RIVERS")?_riverFile:_refFile[i];
      ExitGracefullyIf(file=="","CGroundwaterModel: :GridRefinement RIVERS needs :RiverGeometry",BAD_DATA);
      vector<gw_ring> lines; vector<long long> ids;
      if (!GWGeom::ReadGeoJSONLines(file,(_refKind[i]=="RIVERS")?_riverIDField:"",lines,ids,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
      for (size_t l=0;l<lines.size();l++){
        f.kind=1; f.geom.clear();
        for (size_t k=0;k<lines[l].size();k++){double X,Y; gw_pt p; _proj.Forward(lines[l][k].x,lines[l][k].y,X,Y); _grid.ToLocal(X,Y,p.x,p.y); f.geom.push_back(p);}
        if (f.geom.size()>=2){ref.push_back(f);}
      }
    }
    else if (_refKind[i]=="POLYGONS"){
      map<long long,vector<gw_ring> > polys;
      if (!GWGeom::ReadGeoJSONPolygons(_refFile[i],"",polys,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
      for (map<long long,vector<gw_ring> >::iterator it=polys.begin();it!=polys.end();it++){
        for (size_t r=0;r<it->second.size();r++){
          f.kind=2; f.geom.clear();
          for (size_t k=0;k<it->second[r].size();k++){double X,Y; gw_pt p; _proj.Forward(it->second[r][k].x,it->second[r][k].y,X,Y); _grid.ToLocal(X,Y,p.x,p.y); f.geom.push_back(p);}
          ref.push_back(f);
        }
      }
    }
  }
}
//seeds of the HRU mesh: pairs mirrored across every HRU boundary (so cell faces follow the boundary), seeds on
//refinement lines (rivers run through cell centres), then interior seeds at the base cell size
void CGroundwaterModel::MeshSeeds(const vector<vector<gw_ring> > &rings,const vector<gw_refine> &ref,
                                  const double xa,const double ya,const double xb,const double yb,vector<gw_pt> &seeds) const
{
  seeds.clear();
  double cs=_cellSize,sb=0.5*cs,smin=sb;
  for (size_t i=0;i<ref.size();i++){if (ref[i].kind==1){sb=min(sb,ref[i].size);} smin=min(smin,ref[i].size);}
  double hs=0.25*min(sb,smin); //hash size for spacing checks
  //is a point inside a refinement polygon, or within distance r of its outline?
  auto nearPoly=[&](const gw_ring &P,const double px,const double py,const double r)->bool {
    if (GWGeom::PointInRing(P,px,py)){return true;}
    for (size_t e=0;e<P.size();e++){
      const gw_pt &a=P[e],&b=P[(e+1)%P.size()]; double ex=b.x-a.x,ey=b.y-a.y,L2=ex*ex+ey*ey;
      double t=(L2>0)?max(0.0,min(1.0,((px-a.x)*ex+(py-a.y)*ey)/L2)):0.0;
      if (hypot(px-(a.x+t*ex),py-(a.y+t*ey))<r){return true;}
    }
    return false;
  };
  map<pair<long long,long long>,vector<int> > H;
  struct local { static pair<long long,long long> K(double x,double y,double h){return make_pair((long long)floor(x/h),(long long)floor(y/h));} };
  //nearest distance to accepted seeds within r
  vector<gw_pt> &S=seeds;
  auto nearest=[&](const double px,const double py,const double r)->double { //distance to the nearest accepted seed within r
    double res=1e300; long long ki=(long long)floor(px/hs),kj=(long long)floor(py/hs); int kr=(int)ceil(r/hs);
    for (long long ia=ki-kr;ia<=ki+kr;ia++){for (long long ja=kj-kr;ja<=kj+kr;ja++){
      map<pair<long long,long long>,vector<int> >::iterator it=H.find(make_pair(ia,ja)); if (it==H.end()){continue;}
      for (size_t m=0;m<it->second.size();m++){const gw_pt &q=S[it->second[m]]; res=min(res,hypot(q.x-px,q.y-py));}
    }}
    return res;
  };
  auto addSeed=[&](const double px,const double py){gw_pt q; q.x=px; q.y=py; H[local::K(q.x,q.y,hs)].push_back((int)S.size()); S.push_back(q);};
  //1. mirrored pairs along unique boundary segments
  set<pair<pair<long long,long long>,pair<long long,long long> > > done;
  for (size_t k=0;k<rings.size();k++){
    for (size_t r=0;r<rings[k].size();r++){
      const gw_ring &R=rings[k][r]; size_t m=R.size();
      for (size_t v=0;v<m;v++){ //rings are stored open: the last edge closes back to the first vertex
        const gw_pt &A=R[v],&B=R[(v+1)%m];
        pair<long long,long long> ka((long long)llround(A.x*100),(long long)llround(A.y*100)),kb((long long)llround(B.x*100),(long long)llround(B.y*100));
        if (kb<ka){swap(ka,kb);}
        if (!done.insert(make_pair(ka,kb)).second){continue;}
        double L=hypot(B.x-A.x,B.y-A.y); if (L<1e-6){continue;}
        double sl=sb; //finer spacing where the boundary passes near a refinement polygon (a lake shore)
        for (size_t f=0;f<ref.size();f++){
          if ((ref[f].kind==2) && (ref[f].size<sl) && nearPoly(ref[f].geom,0.5*(A.x+B.x),0.5*(A.y+B.y),ref[f].size)){sl=ref[f].size;}
        }
        double nx=(A.y-B.y)/L,ny=(B.x-A.x)/L,d=0.3*sl;
        int np=max(1,(int)floor(L/sl+0.5));
        for (int j=0;j<np;j++){
          double t=(j+0.5)/np,px=A.x+t*(B.x-A.x),py=A.y+t*(B.y-A.y),dmin;
          dmin=nearest(px,py,0.6*sl); if (dmin<0.6*sl){continue;} //keep pairs whole: skip both seeds
          addSeed(px+d*nx,py+d*ny); addSeed(px-d*nx,py-d*ny);
        }
      }
    }
  }
  //2. seeds on refinement lines, then around refinement points
  for (size_t i=0;i<ref.size();i++){
    const gw_ring &g=ref[i].geom; double s=ref[i].size;
    if (ref[i].kind==1){
      for (size_t v=0;v+1<g.size();v++){
        double L=hypot(g[v+1].x-g[v].x,g[v+1].y-g[v].y); int np=max(1,(int)floor(L/s+0.5));
        for (int j=0;j<np;j++){double t=(j+0.5)/np,px=g[v].x+t*(g[v+1].x-g[v].x),py=g[v].y+t*(g[v+1].y-g[v].y),dmin;
          if ((px<xa) || (px>xb) || (py<ya) || (py>yb)){continue;}
          dmin=nearest(px,py,0.5*s); if (dmin<0.5*s){continue;} addSeed(px,py);}
      }
    }
    else if (ref[i].kind==0){
      for (int ring=0;ring<=2;ring++){int np=(ring==0)?1:6*ring; for (int j=0;j<np;j++){
        double a=2*PI*j/np,px=g[0].x+ring*s*cos(a),py=g[0].y+ring*s*sin(a),dmin;
        if ((px<xa) || (px>xb) || (py<ya) || (py>yb)){continue;} //e.g. a monitoring well outside the model
        dmin=nearest(px,py,0.5*s); if (dmin<0.5*s){continue;} addSeed(px,py);}}
    }
  }
  //2b. refinement polygons (lakes): a staggered lattice at the target spacing inside and within one spacing around them
  for (size_t i=0;i<ref.size();i++){
    if (ref[i].kind!=2){continue;}
    const gw_ring &g=ref[i].geom; double s=ref[i].size,gx0=1e300,gy0=1e300,gx1=-1e300,gy1=-1e300;
    for (size_t v=0;v<g.size();v++){gx0=min(gx0,g[v].x); gx1=max(gx1,g[v].x); gy0=min(gy0,g[v].y); gy1=max(gy1,g[v].y);}
    int row=0;
    for (double y=gy0-s;y<=gy1+s;y+=0.8660254*s,row++){
      for (double x=gx0-s+((row%2)?0.5*s:0.0);x<=gx1+s;x+=s){
        if ((x<xa) || (x>xb) || (y<ya) || (y>yb) || !nearPoly(g,x,y,s)){continue;}
        double dmin=nearest(x,y,0.55*s); if (dmin<0.55*s){continue;} addSeed(x,y);
      }
    }
  }
  //3. interior seeds on a lattice at the base cell size
  for (double y=ya+0.5*cs;y<yb;y+=cs){
    for (double x=xa+0.5*cs;x<xb;x+=cs){
      double dmin; dmin=nearest(x,y,0.55*cs); if (dmin<0.55*cs){continue;} addSeed(x,y);
    }
  }
}

//=================================================================
// parse-time setters
//=================================================================
void CGroundwaterModel::AddAquiferClass(const gw_aquifer_class &c)
{
  ExitGracefullyIf(GetAquiferClassIndex(c.name)!=DOESNT_EXIST,("duplicate aquifer class "+c.name).c_str(),BAD_DATA);
  _aClasses.push_back(c);
}
int CGroundwaterModel::GetAquiferClassIndex(const string &name) const
{
  for (size_t i=0;i<_aClasses.size();i++){if (_aClasses[i].name==name){return (int)i;}}
  return DOESNT_EXIST;
}
void CGroundwaterModel::AddProfile(const gw_profile &p)
{
  ExitGracefullyIf(GetProfileIndex(p.name)!=DOESNT_EXIST,("duplicate aquifer profile "+p.name).c_str(),BAD_DATA);
  ExitGracefullyIf(p.layers.size()==0,("aquifer profile "+p.name+" has no layers").c_str(),BAD_DATA);
  ExitGracefullyIf(p.layers[0].type==GWL_AQUICLUDE,("aquifer profile "+p.name+": top layer cannot be an AQUICLUDE").c_str(),BAD_DATA);
  _aProfiles.push_back(p);
}
int CGroundwaterModel::GetProfileIndex(const string &name) const
{
  for (size_t i=0;i<_aProfiles.size();i++){if (_aProfiles[i].name==name){return (int)i;}}
  return DOESNT_EXIST;
}
gw_profile CGroundwaterModel::DefaultProfile(const string &name)
{
  gw_profile p; p.name=name;
  p.init_head_depth=5.0; p.default_bedrock_depth=50.0;
  p.riverbed_K=0.5; p.riverbed_thick=1.0; p.seep_leakance=1.0;
  p.recharge_delay=0.0; p.soil_zone_depth=0.0; p.extinction_depth=0.0; p.seep_smooth_depth=0.5;
  return p;
}
bool CGroundwaterModel::SetProfileParameter(const string &profile,const string &param,const double &value)
{
  int ip=GetProfileIndex(profile);
  if (ip==DOESNT_EXIST){return false;}
  gw_profile &p=_aProfiles[ip];
  if      (param=="INITIAL_HEAD_DEPTH"   ){p.init_head_depth=value;}
  else if (param=="DEFAULT_BEDROCK_DEPTH"){p.default_bedrock_depth=value;}
  else if (param=="RIVERBED_K"           ){p.riverbed_K=value;}
  else if (param=="RIVERBED_THICKNESS"   ){p.riverbed_thick=value;}
  else if (param=="SEEPAGE_LEAKANCE"     ){p.seep_leakance=value;}
  else if (param=="RECHARGE_DELAY"       ){p.recharge_delay=value;}
  else if (param=="SOIL_ZONE_DEPTH"      ){p.soil_zone_depth=value;}
  else if (param=="EXTINCTION_DEPTH"     ){p.extinction_depth=value;}
  else if (param=="SEEPAGE_SMOOTHING_DEPTH"){p.seep_smooth_depth=value;} //checked at start-up (may not be negative)
  else {
    WriteWarning("CGroundwaterModel: unrecognized aquifer profile parameter "+param+" (ignored)",false);
  }
  return true;
}

//=================================================================
// Initialization
//=================================================================
void CGroundwaterModel::Initialize(const optStruct &Options)
{
  //Raven's ensemble drivers (Monte Carlo, DDS, EnKF) re-run CalculateInitialWaterStorage for every member:
  //each call rebuilds the groundwater model from scratch and restarts MODFLOW 6
  if (_active || (_aqHRUs.size()>0)){ResetState();}
  if (!Options.silent){cout<<"  Initializing groundwater model (MODFLOW 6)..."<<endl;}
  ValidateInputs();
  _jday0=Options.julian_start_day; _jyear0=Options.julian_start_year; _calendar=Options.calendar;
  _nstp=(int)floor(Options.duration/Options.timestep+0.5); _stepCount=0;
  _subDaily=(Options.timestep<1.0-1e-9);
  if (fabs(_nstp*Options.timestep-Options.duration)>1e-6){
    WriteWarning("CGroundwaterModel: :Duration is not a whole number of time steps; MODFLOW 6 runs "+to_string(_nstp)+" steps",Options.noisy);
  }
  ExitGracefullyIf(_aProfiles.size()==0,"CGroundwaterModel: coupled model requires :AquiferProfiles in the .rvp file",BAD_DATA);

  _iGW=_pModel->GetStateVarIndex(GROUNDWATER);
  ExitGracefullyIf(_iGW==DOESNT_EXIST,
    "CGroundwaterModel: no process moves water into GROUNDWATER, so no recharge can reach the aquifer. Route recharge to GROUNDWATER in the .rvi",BAD_DATA);

  //-- which HRUs are coupled -----------------------------------------
  int nHRUs=_pModel->GetNumHRUs();
  _hruProfile.assign(nHRUs,DOESNT_EXIST);
  for (int k=0;k<nHRUs;k++)
  {
    string name=_pModel->GetHydroUnit(k)->GetAquiferProfileName();
    if ((name=="") || (name=="[NONE]") || (name=="NONE")){continue;}
    if (!_pModel->GetHydroUnit(k)->IsEnabled()){
      WriteWarning("CGroundwaterModel: HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+" is disabled; it is not coupled to the groundwater model",false);
      continue;
    }
    int ip=GetProfileIndex(name);
    if (ip==DOESNT_EXIST){
      ExitGracefully(("CGroundwaterModel: HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+
                     " references unknown AQUIFER_PROFILE "+name).c_str(),BAD_DATA);
    }
    _hruProfile[k]=ip;
    _aqHRUs.push_back(k);
  }
  ExitGracefullyIf(_aqHRUs.size()==0,"CGroundwaterModel: no HRU has an AQUIFER_PROFILE in the .rvh file",BAD_DATA);
  _rchStore.assign(nHRUs,0.0); //recharge-delay reservoirs, possibly from hotstart
  for (size_t i=0;i<_aqHRUs.size();i++){
    map<long long,double>::iterator it=_hotRchStore.find(_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID());
    if (it!=_hotRchStore.end()){_rchStore[_aqHRUs[i]]=it->second;}
  }

  //-- where groundwater discharge returns to --------------------------
  if (_seepTarget==""){
    int nsoil=_pModel->GetNumSoilLayers();
    _iReturnSV=(nsoil>0)?_pModel->GetStateVarIndex(SOIL,nsoil-1):_pModel->GetStateVarIndex(SURFACE_WATER);
  }
  else {
    int layer=0;
    sv_type typ=_pModel->GetStateVarInfo()->StringToSVType(_seepTarget,layer,true);
    _iReturnSV=_pModel->GetStateVarIndex(typ,layer);
  }
  ExitGracefullyIf(_iReturnSV==DOESNT_EXIST,"CGroundwaterModel: :SeepageReturnTo storage is not used in this model",BAD_DATA);
  //HRUs without soil (lakes, glaciers, rock) cannot hold seepage in a soil layer: it goes to their surface water
  _hruReturnSV.assign(nHRUs,_iReturnSV);
  int lay_sv=0;
  bool toSoil=(_seepTarget=="") || (_pModel->GetStateVarInfo()->StringToSVType(_seepTarget,lay_sv,false)==SOIL);
  for (size_t i=0;i<_aqHRUs.size();i++){
    int k=_aqHRUs[i]; double th=0;
    for (int m=0;m<_pModel->GetNumSoilLayers();m++){th+=_pModel->GetHydroUnit(k)->GetSoilThickness(m);}
    if (toSoil && (th<=0.0)){_hruReturnSV[k]=_pModel->GetStateVarIndex(SURFACE_WATER);}
  }

  CheckProcesses();
  if (IsExternal()){ExtBuild(Options); GWPROF("ExtBuild");}
  else {
  BuildGeometry(Options); GWPROF("BuildGeometry");
  _disType=(_gridType==GWGRID_REGULAR)?0:(_gridType==GWGRID_QUADTREE)?1:2;
  if (UseXT3D()){
    WriteWarning("CGroundwaterModel: :FlowCorrection XT3D can slow or prevent the Newton convergence where cells dry and rewet; check the 'converged' column of GWBudget.csv",Options.noisy);
  }
  if (_connectByElevation){_disType=2;}
  BuildVertical(Options); GWPROF("BuildVertical");
  BuildLinks(); GWPROF("BuildLinks");
  BuildRivers(); GWPROF("BuildRivers");
  BuildStresses(); GWPROF("BuildStresses");
  } //generated model
  for (map<long long,CTimeSeries*>::iterator it=_wellTS.begin();it!=_wellTS.end();it++){
    it->second->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  for (map<string,CTimeSeries*>::iterator it=_headTS.begin();it!=_headTS.end();it++){
    it->second->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar);
  }
  //MODFLOW ignores recharge on specified-head cells: recharge is spread over each HRU's other cells
  _colCHDTop.assign(_ncpl,0);
  for (size_t e=0;e<_chdNode.size();e++){int ic=_chdNode[e]%_ncpl; bool top=true;
    for (int l=0;l<_chdNode[e]/_ncpl;l++){if (_idomain[UserNode(l,ic)]>0){top=false;}} if (top){_colCHDTop[ic]=1;}}
  _hruRchArea.assign(_pModel->GetNumHRUs(),0.0);
  for (size_t li=0;li<_aLinks.size();li++){if (!_colCHDTop[_aLinks[li].ic]){_hruRchArea[_aLinks[li].k]+=_aLinks[li].area;}}
  //surface-water store <-> aquifer exchanges: one head-dependent cell per column holding group HRUs
  for (size_t x=0;x<_aSWX.size();x++){
    gw_swx &S=_aSWX[x]; int lay=0;
    sv_type typ=_pModel->GetStateVarInfo()->StringToSVType(S.sv,lay,true);
    S.iSV=_pModel->GetStateVarIndex(typ,lay);
    ExitGracefullyIf(S.iSV==DOESNT_EXIST,("CGroundwaterModel: :SurfaceWaterExchange store "+S.sv+" is not used by this model").c_str(),BAD_DATA);
    CHRUGroup *pG=_pModel->GetHRUGroup(S.group);
    ExitGracefullyIf(pG==NULL,("CGroundwaterModel: :SurfaceWaterExchange HRU group "+S.group+" does not exist").c_str(),BAD_DATA);
    S.kk=pG->GetGlobalIndex();
    map<int,vector<int> > bycol;
    for (size_t li=0;li<_aLinks.size();li++){
      if (pG->IsInGroup(_aLinks[li].k) && !_colCHDTop[_aLinks[li].ic]){bycol[_aLinks[li].ic].push_back((int)li);}
    }
    if (bycol.size()==0){WriteWarning("CGroundwaterModel: :SurfaceWaterExchange group "+S.group+" has no coupled HRUs on active cells",false);}
    for (size_t y=0;y<x;y++){ //the same store of the same HRU may not be drawn by two exchanges
      if (_aSWX[y].iSV!=S.iSV){continue;}
      CHRUGroup *pG2=_pModel->GetHRUGroup(_aSWX[y].group);
      for (size_t i=0;i<_aqHRUs.size();i++){
        if (pG->IsInGroup(_aqHRUs[i]) && pG2->IsInGroup(_aqHRUs[i])){
          ExitGracefully(("CGroundwaterModel: HRU "+to_string(_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID())+" is in two :SurfaceWaterExchange groups ("+
                          _aSWX[y].group+", "+S.group+") for the same store "+S.sv).c_str(),BAD_DATA);
        }
      }
    }
    for (map<int,vector<int> >::iterator it=bycol.begin();it!=bycol.end();it++){_swxCol.push_back(it->first); _swxX.push_back((int)x); _swxLinks.push_back(it->second);}
  }
  //reservoirs: seepage of each Raven reservoir enters the top cells under its lake HRU (by overlap)
  if (_resExchange){
    for (int p=0;p<_pModel->GetNumSubBasins();p++){
      CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir();
      if ((pRes==NULL) || (!_pModel->GetSubBasin(p)->IsEnabled())){continue;}
      int k=pRes->GetHRUIndex();
      string nm=to_string(_pModel->GetSubBasin(p)->GetID());
      if ((k==DOESNT_EXIST) || (_hruProfile[k]==DOESNT_EXIST)){WriteWarning("CGroundwaterModel: reservoir in subbasin "+nm+" has no coupled lake HRU (:HRUID with an AQUIFER_PROFILE): not coupled",false); continue;}
      if (pRes->GetSeepageConstant()<=0){WriteWarning("CGroundwaterModel: reservoir in subbasin "+nm+" has no :SeepageParameters: not coupled",false); continue;}
      int r=(int)_resP.size(); double wsum=0; size_t n0=_resNode.size();
      for (size_t j=0;j<_hruLinks[k].size();j++){
        const gw_link &L=_aLinks[_hruLinks[k][j]]; if (_colCHDTop[L.ic]){continue;}
        _resNode.push_back(TopNode(L.ic)); _resR.push_back(r); _resW.push_back(L.area); wsum+=L.area;
      }
      if (wsum<=0){_resNode.resize(n0);_resR.resize(n0);_resW.resize(n0); continue;}
      for (size_t e=n0;e<_resW.size();e++){_resW[e]/=wsum;}
      //the reservoir's stages must be elevations on the aquifer's datum: a stage below the whole aquifer under the
      //lake means relative stages (crest = 0) or another datum, and would drive a huge false inflow
      double zbot=1e300;
      for (size_t e=n0;e<_resNode.size();e++){
        int ic=_resNode[e]%_ncpl; for (int l=_nlay-1;l>=0;l--){int n=UserNode(l,ic); if (_idomain[n]>0){zbot=min(zbot,_botm[n]); break;}}
      }
      ExitGracefullyIf((zbot<1e299) && (pRes->GetResStage()<zbot),("CGroundwaterModel: reservoir in subbasin "+nm+
        " has stage "+to_string(pRes->GetResStage())+" m, below the bottom of the aquifer under its lake ("+to_string(zbot)+
        " m): coupled reservoirs need absolute stages on the aquifer's datum (e.g. :AbsoluteCrestHeight)").c_str(),BAD_DATA);
      _resP.push_back(p);
    }
  }
  _doSteady=_steadyInit;
  if (_doSteady && (_hotHeads.size()>0)){
    WriteAdvisory("CGroundwaterModel: hotstart heads found in .rvc; steady-state initialization skipped",Options.noisy);
    _doSteady=false; //decided per run: the user's setting is kept for later ensemble members
  }
  if (IsExternal()){ExtWritePackages(Options); GWPROF("ExtWritePackages"); if (_owTable.empty()){WriteGridCells(Options); GWPROF("WriteGridCells");}}
  else {
  WriteMF6Files(Options); GWPROF("WriteMF6Files");
  WriteGridCells(Options); GWPROF("WriteGridCells");
  }
  ConnectEngine(); GWPROF("ConnectEngine");
  _active=true;
  if (_doSteady){RunSteadyState();}
  OpenOutputs(Options); GWPROF("OpenOutputs");
  WriteSummary(Options); GWPROF("WriteSummary");
  if (_writeNC){OpenNetCDF(Options); WriteNetCDF(0.0);}
  if (_HEADS.is_open()){ //initial heads (after any steady-state initialization) at t=0
    _HEADS<<"0,"<<EndDate(0.0)<<setprecision(10);
    for (size_t o=0;o<_aObsWells.size();o++){_HEADS<<","<<GetObservationWellHead(_aObsWells[o].id);}
    _HEADS<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief MODFLOW owns groundwater storage in coupled HRUs: no Raven process may draw from GROUNDWATER there
//
//////////////////////////////////////////////////////////////////
/// \brief water-table elevation of a column: head of the highest wet cell, or the bottom of the lowest active cell if
///  all are dry. Under the Newton method a dry cell carries no flow and its head is numerically arbitrary (it can be far
///  below the cell bottom), so raw heads must not be used as physical water levels.
//
double CGroundwaterModel::WaterTable(const int ic) const
{
  double lowest=_top[ic];
  for (int l=0;l<_nlay;l++){
    int n=UserNode(l,ic); if (_idomain[n]<=0){continue;}
    double h=_X[_reduced[n]]*_hf; //[m]
    if (h>_botm[n]){return h;}
    lowest=_botm[n];
  }
  return lowest;
}
void CGroundwaterModel::ClearHotstartData()
{
  _hotHeads.clear(); _hotDims[0]=_hotDims[1]=_hotDims[2]=0;
  _hotRchStore.clear(); _hotResSeep.clear(); _hotRivCarry.clear();
  _initMode="FROM_PROFILE"; _initValue=0.0;
}
void CGroundwaterModel::ValidateInputs() const
{
  //every value that would make the model meaningless, create water, or make MODFLOW fail is rejected with a clear message
  for (size_t c=0;c<_aClasses.size();c++){
    const gw_aquifer_class &A=_aClasses[c]; string n="aquifer class "+A.name+": ";
    ExitGracefullyIf(!(A.Kh>0),(n+"K_HORIZ must be positive").c_str(),BAD_DATA);
    ExitGracefullyIf(!(A.Kv>0),(n+"K_VERT must be positive").c_str(),BAD_DATA);
    ExitGracefullyIf(!(A.Ss>=0),(n+"SPEC_STORAGE may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!((A.Sy>=0) && (A.Sy<=1)),(n+"SPEC_YIELD must lie between 0 and 1").c_str(),BAD_DATA);
    ExitGracefullyIf(!((A.porosity>=0) && (A.porosity<=1)),(n+"POROSITY must lie between 0 and 1").c_str(),BAD_DATA);
  }
  for (size_t p=0;p<_aProfiles.size();p++){
    const gw_profile &P=_aProfiles[p]; string n="aquifer profile "+P.name+": ";
    for (size_t l=0;l<P.layers.size();l++){
      ExitGracefullyIf((!P.layers[l].to_bedrock) && !(P.layers[l].thickness>0),(n+"layer thicknesses must be positive").c_str(),BAD_DATA);
    }
    ExitGracefullyIf(!(P.default_bedrock_depth>0),(n+"DEFAULT_BEDROCK_DEPTH must be positive").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.riverbed_K>=0),(n+"RIVERBED_K may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.riverbed_thick>0),(n+"RIVERBED_THICKNESS must be positive").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.seep_leakance>=0),(n+"SEEPAGE_LEAKANCE may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.recharge_delay>=0),(n+"RECHARGE_DELAY may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.soil_zone_depth>=0),(n+"SOIL_ZONE_DEPTH may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.extinction_depth>=0),(n+"EXTINCTION_DEPTH may not be negative").c_str(),BAD_DATA);
    ExitGracefullyIf(!(P.seep_smooth_depth>=0),(n+"SEEPAGE_SMOOTHING_DEPTH may not be negative").c_str(),BAD_DATA);
  }
  ExitGracefullyIf(!(_cellSize>0),"CGroundwaterModel: :GridCellSize must be positive",BAD_DATA);
  if (IsExternal()){ //existing MODFLOW 6 model: the model supplies grid, hydrogeology, wells and boundaries
    const char *msg="CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove ";
    int nlink=((_extGridFile!="")?1:0)+(_extCRS.IsSet()?1:0)+((_owTable.size()>0)?1:0);
    ExitGracefullyIf(nlink!=1,"CGroundwaterModel: with :MF6Simulation give one of :MF6CRS, :MF6GridFile or :OverlapWeights (to link HRUs to the model's cells)",BAD_DATA);
    ExitGracefullyIf(_extStart.size()<10,"CGroundwaterModel: :MF6Simulation needs :MF6StartDate yyyy-mm-dd",BAD_DATA);
    //(the user's settings are checked, not _gridType: linking an existing model sets it internally, and ensemble
    // members initialise the model again)
    ExitGracefullyIf(_gridTypeGiven || (_gridAngle!=0.0) || (_refKind.size()>0) || (_qtMaxLevel>=0) || (_flowCorr!="AUTO") || _connectByElevation,
      (string(msg)+"the grid options (:GridType, :GridRotation, :GridRefinement, :QuadtreeMaxLevel, :FlowCorrection, :LayerConnection)").c_str(),BAD_DATA);
    ExitGracefullyIf((_aWells.size()>0) || (_aBounds.size()>0),(string(msg)+":Wells and head boundaries (keep the model's packages)").c_str(),BAD_DATA);
    ExitGracefullyIf((_aSWX.size()>0) || _resExchange,(string(msg)+":SurfaceWaterExchange and :ReservoirExchange (not yet available with an existing model)").c_str(),BAD_DATA);
    ExitGracefullyIf(_steadyInit,(string(msg)+":GWInitialization (the model's own stress periods decide)").c_str(),BAD_DATA);
    ExitGracefullyIf(_genOnly!="",(string(msg)+_genOnly+" (the existing model's grid, surfaces and solver settings are used)").c_str(),BAD_DATA);
    ExitGracefullyIf(_writeNC && (_owTable.size()>0),"CGroundwaterModel: :WriteNetCDFHeads needs the cell outlines (:MF6CRS or :MF6GridFile); an :OverlapWeights table has no coordinates",BAD_DATA);
    ExitGracefullyIf(_riverFile!="",(string(msg)+":RiverGeometry (rivers are the model's stream package, :StreamPackage)").c_str(),BAD_DATA);
    ExitGracefullyIf((_aObsWells.size()>0) && (_owTable.size()>0),"CGroundwaterModel: monitoring wells need :MF6CRS or :MF6GridFile (to locate them)",BAD_DATA);
    for (size_t i=0;i<_aProfiles.size();i++){
      ExitGracefullyIf(_aProfiles[i].extinction_depth>0,"CGroundwaterModel: water-table ET (EXTINCTION_DEPTH) is not yet available with an existing model",BAD_DATA);
    }
  }
  else {
    ExitGracefullyIf((_extGridFile!="") || _extCRS.IsSet() || (_owTable.size()>0) || (_takeOver.size()>0) || (_streamPkg!=""),
      "CGroundwaterModel: :MF6CRS, :MF6GridFile, :OverlapWeights, :RavenTakesOver and :StreamPackage need :MF6Simulation",BAD_DATA);
  }
  if (!IsExternal()){ //(an existing model's grid type is set internally when its cells are linked)
    ExitGracefullyIf((_gridType==GWGRID_FILE) && (_gridFile==""),"CGroundwaterModel: :GridType FILE needs :GridFile",BAD_DATA);
    ExitGracefullyIf((_gridType!=GWGRID_FILE) && (_gridFile!=""),"CGroundwaterModel: :GridFile needs :GridType FILE",BAD_DATA);
  }
  ExitGracefullyIf((_gridAngle!=0.0) && (_gridType!=GWGRID_REGULAR) && (_gridType!=GWGRID_QUADTREE),"CGroundwaterModel: :GridRotation applies only to REGULAR and QUADTREE grids",BAD_DATA);
  ExitGracefullyIf(!(fabs(_gridAngle)<=360.0),"CGroundwaterModel: :GridRotation must be between -360 and 360 degrees",BAD_DATA);
  ExitGracefullyIf((_qtMaxLevel>12),"CGroundwaterModel: :QuadtreeMaxLevel may not exceed 12",BAD_DATA);
  ExitGracefullyIf((_qtMaxLevel>=0) && (_gridType!=GWGRID_QUADTREE),"CGroundwaterModel: :QuadtreeMaxLevel needs :GridType QUADTREE",BAD_DATA);
  ExitGracefullyIf((_refKind.size()>0) && (_gridType==GWGRID_FILE),"CGroundwaterModel: :GridRefinement does not apply to imported grids (:GridType FILE)",BAD_DATA);
  for (size_t i=0;i<_refKind.size();i++){
    ExitGracefullyIf(!(_refSize[i]>0),"CGroundwaterModel: :GridRefinement size must be positive",BAD_DATA);
    ExitGracefullyIf(_refSize[i]>=_cellSize,"CGroundwaterModel: :GridRefinement size must be smaller than :GridCellSize",BAD_DATA);
    ExitGracefullyIf((_refKind[i]=="RIVERS") && (_riverFile==""),"CGroundwaterModel: :GridRefinement RIVERS needs :RiverGeometry",BAD_DATA);
    ExitGracefullyIf((_refKind[i]=="WELLS") && (_aWells.size()+_aObsWells.size()==0),"CGroundwaterModel: :GridRefinement WELLS needs :Wells or :ObservationWells",BAD_DATA);
  }
  ExitGracefullyIf(!((_minCoverage>0) && (_minCoverage<=1)),"CGroundwaterModel: :MinCellCoverage must lie in (0,1]",BAD_DATA);
  ExitGracefullyIf(!(_dvclose>0),"CGroundwaterModel: :SolverHeadTolerance must be positive",BAD_DATA);
  ExitGracefullyIf(_maxiter<1,"CGroundwaterModel: :SolverMaxOuterIterations must be at least 1",BAD_DATA);
  ExitGracefullyIf(_headSaveFreq<1,"CGroundwaterModel: :HeadSaveFrequency must be at least 1",BAD_DATA);
  ExitGracefullyIf(_steadyInit && !(_ssRecharge>=0),"CGroundwaterModel: :GWInitialization STEADY_STATE recharge may not be negative",BAD_DATA);
  for (size_t x=0;x<_aSWX.size();x++){ExitGracefullyIf(!(_aSWX[x].leakance>0),("CGroundwaterModel: :SurfaceWaterExchange "+_aSWX[x].group+": LEAKANCE must be positive").c_str(),BAD_DATA);}
}
void CGroundwaterModel::CheckProcesses()
{
  bool warned=false;
  for (int j=0;j<_pModel->GetNumProcesses();j++)
  {
    CHydroProcessABC *pProc=_pModel->GetProcess(j);
    const int *iFrom=pProc->GetFromIndices();
    const int *iTo  =pProc->GetToIndices();
    for (int q=0;q<pProc->GetNumConnections();q++)
    {
      if ((iTo[q]==_iGW) && (iFrom[q]!=_iGW) && !warned){ //water sent to GROUNDWATER in an uncoupled HRU leaves the model
        for (int k=0;k<_pModel->GetNumHRUs();k++){
          if ((_hruProfile[k]==DOESNT_EXIST) && _pModel->GetHydroUnit(k)->IsEnabled() && pProc->ShouldApply(_pModel->GetHydroUnit(k))){
            WriteWarning("CGroundwaterModel: process #"+to_string(j+1)+" sends water to GROUNDWATER in HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+
                         ", which has no AQUIFER_PROFILE: that water leaves the model. Restrict the process with :-->Conditional",false);
            warned=true; break;
          }
        }
      }
      if ((iFrom[q]!=_iGW) || (iTo[q]==_iGW)){continue;}
      for (size_t i=0;i<_aqHRUs.size();i++){
        if (pProc->ShouldApply(_pModel->GetHydroUnit(_aqHRUs[i]))){
          ExitGracefully(("CGroundwaterModel: process #"+to_string(j+1)+" removes water from GROUNDWATER in HRU "+
            to_string(_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID())+
            ", which is coupled to MODFLOW. MODFLOW owns that storage; restrict the process with :-->Conditional").c_str(),BAD_DATA);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief reads HRU polygons, projects them, lays the grid and intersects HRUs with cells
//
void CGroundwaterModel::BuildGeometry(const optStruct &Options)
{
  ExitGracefullyIf(_hruGeomFile=="","CGroundwaterModel: :HRUGeometry is required in the .rvg file",BAD_DATA);
  clock_t tclk=clock();
  map<long long,vector<gw_ring> > shapes;
  string err;
  //-- linkage cache key: FNV-1a hash of the HRU geometry file + everything that shapes the intersection
  uint64_t hsh=1469598103934665603ULL;
  { ifstream G(_hruGeomFile.c_str(),ios::binary); char buf[65536];
    while (G.read(buf,sizeof(buf)) || G.gcount()>0){for (streamsize i=0;i<G.gcount();i++){hsh^=(unsigned char)buf[i]; hsh*=1099511628211ULL;}} }
  ostringstream key; key<<hex<<hsh<<dec<<"|"<<_hruIDField<<"|"<<setprecision(12)<<_cellSize<<"|"<<GridSpec();
  for (size_t i=0;i<_aqHRUs.size();i++){key<<"|"<<_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID();}
  string cachefile=(_cacheFile!="")?_cacheFile:Options.output_dir+"mf6/linkage.cache";
  GW_MKDIR(Options.output_dir.c_str()); GW_MKDIR((Options.output_dir+"mf6").c_str());
  _cacheHit=false;
  { ifstream CF(cachefile.c_str());
    string ckey; if (CF && getline(CF,ckey) && (ckey==key.str())){
      double lat0,lon0; int nH; CF>>lat0>>lon0;
      _proj.SetOrigin(lat0,lon0);
      bool gok=_grid.Read(CF); CF>>nH; if (!gok){CF.setstate(ios::failbit);}
      _ncpl=_grid.ncpl; _nrow=_grid.nrow; _ncol=_grid.ncol; _x0=_grid.x0; _y0=_grid.y0;
      int nHRUs0=_pModel->GetNumHRUs(); _hruPolyArea.assign(nHRUs0,0.0); _hruCells.assign(nHRUs0,map<int,double>());
      bool ok=!CF.fail();
      for (int j=0;(j<nH) && ok;j++){
        long long id; double pa; int nc; CF>>id>>pa>>nc;
        CHydroUnit *pH=_pModel->GetHRUByID(id); ok=(pH!=NULL) && !CF.fail();
        if (!ok){break;}
        int k=pH->GetGlobalIndex(); _hruPolyArea[k]=pa;
        for (int c=0;c<nc;c++){int ic; double a; CF>>ic>>a; _hruCells[k][ic]=a;}
      }
      _cacheHit=ok && !CF.fail();
    }
  }
  if (!_cacheHit){
  if (!GWGeom::ReadGeoJSONPolygons(_hruGeomFile,_hruIDField,shapes,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}

  int nHRUs=_pModel->GetNumHRUs();
  //-- projection origin = mean centroid of coupled HRUs
  double slat=0,slon=0;
  for (size_t i=0;i<_aqHRUs.size();i++){
    location loc=_pModel->GetHydroUnit(_aqHRUs[i])->GetCentroid();
    slat+=loc.latitude; slon+=loc.longitude;
  }
  _proj.SetOrigin(slat/_aqHRUs.size(),slon/_aqHRUs.size());

  double cs;
  //-- project polygons
  vector<vector<gw_ring> > prj(nHRUs);
  _hruPolyArea.assign(nHRUs,0.0);
  double minx=1e30,miny=1e30,maxx=-1e30,maxy=-1e30;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i];
    CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
    map<long long,vector<gw_ring> >::iterator it=shapes.find(pHRU->GetHRUID());
    if (it==shapes.end()){
      ExitGracefully(("CGroundwaterModel: HRU "+to_string(pHRU->GetHRUID())+" has an AQUIFER_PROFILE but no polygon in "+_hruGeomFile).c_str(),BAD_DATA);
      continue;
    }
    for (size_t r=0;r<it->second.size();r++){
      gw_ring ring;
      for (size_t n=0;n<it->second[r].size();n++){
        gw_pt q; _proj.Forward(it->second[r][n].x,it->second[r][n].y,q.x,q.y);
        ring.push_back(q);
        minx=min(minx,q.x);maxx=max(maxx,q.x);miny=min(miny,q.y);maxy=max(maxy,q.y);
      }
      _hruPolyArea[k]+=GWGeom::SignedArea(ring);
      prj[k].push_back(ring);
    }
    double Arvh=pHRU->GetArea()*M2_PER_KM2;
    if (fabs(_hruPolyArea[k]-Arvh)>0.05*Arvh){
      WriteWarning("CGroundwaterModel: polygon area of HRU "+to_string(pHRU->GetHRUID())+" ("+ToStr(_hruPolyArea[k]/M2_PER_KM2)+
                   " km2) differs from .rvh AREA ("+ToStr(pHRU->GetArea())+" km2) by more than 5%. Flux mapping stays mass-conservative.",Options.noisy);
    }
  }

  //-- grid (local frame: identity for unrotated grids, so a uniform regular grid reproduces the original arithmetic)
  cs=_cellSize;
  ExitGracefullyIf(cs<=0,"CGroundwaterModel: :GridCellSize must be positive",BAD_DATA);
  _grid.ox=_grid.oy=0.0; _grid.angle=((_gridType==GWGRID_REGULAR) || (_gridType==GWGRID_QUADTREE))?_gridAngle:0.0;
  if (_grid.angle!=0.0){
    minx=miny=1e30; maxx=maxy=-1e30;
    for (size_t i=0;i<_aqHRUs.size();i++){int k=_aqHRUs[i];
      for (size_t r=0;r<prj[k].size();r++){for (size_t n=0;n<prj[k][r].size();n++){
        double x,y; _grid.ToLocal(prj[k][r][n].x,prj[k][r][n].y,x,y); prj[k][r][n].x=x; prj[k][r][n].y=y;
        minx=min(minx,x);maxx=max(maxx,x);miny=min(miny,y);maxy=max(maxy,y);}}}
  }
  vector<gw_refine> ref; BuildRefinement(ref);
  for (size_t i=0;i<_refKind.size();i++){ //:GridRefinement LAKES - the outline of every reservoir's lake HRU
    if (_refKind[i]!="LAKES"){continue;}
    int nl=0;
    for (int p=0;p<_pModel->GetNumSubBasins();p++){
      const CReservoir *pR=_pModel->GetSubBasin(p)->GetReservoir(); if (pR==NULL){continue;}
      int kh=pR->GetHRUIndex(); if (kh<0){continue;}
      map<long long,vector<gw_ring> >::iterator it=shapes.find(_pModel->GetHydroUnit(kh)->GetHRUID()); if (it==shapes.end()){continue;}
      for (size_t r=0;r<it->second.size();r++){
        gw_refine f; f.kind=2; f.size=_refSize[i];
        for (size_t n=0;n<it->second[r].size();n++){double X,Y; gw_pt q; _proj.Forward(it->second[r][n].x,it->second[r][n].y,X,Y); _grid.ToLocal(X,Y,q.x,q.y); f.geom.push_back(q);}
        if (f.geom.size()>=3){ref.push_back(f); nl++;}
      }
    }
    ExitGracefullyIf(nl==0,"CGroundwaterModel: :GridRefinement LAKES found no reservoir whose lake HRU has a polygon in :HRUGeometry (:Reservoir ... :HRUID)",BAD_DATA);
  }
  if (_gridType==GWGRID_REGULAR){_grid.MakeRegular(minx,miny,maxx,maxy,cs,ref);}
  else if (_gridType==GWGRID_QUADTREE){
    double smin=cs; for (size_t i=0;i<ref.size();i++){smin=min(smin,ref[i].size);}
    int lev=(_qtMaxLevel>=0)?_qtMaxLevel:max(0,(int)ceil(log(cs/smin)/log(2.0)-1e-9));
    _grid.MakeQuadtree(minx,miny,maxx,maxy,cs,ref,lev);
  }
  else if (_gridType==GWGRID_HRUMESH){
    vector<vector<gw_ring> > rr; for (size_t i=0;i<_aqHRUs.size();i++){rr.push_back(prj[_aqHRUs[i]]);}
    vector<gw_pt> seeds; MeshSeeds(rr,ref,minx-cs,miny-cs,maxx+cs,maxy+cs,seeds);
    ExitGracefullyIf(seeds.size()>5000000,"CGroundwaterModel: HRU mesh exceeds 5 million cells per layer; increase :GridCellSize",BAD_DATA);
    _grid.MakeVoronoi(seeds,minx-cs,miny-cs,maxx+cs,maxy+cs);
    for (int i=0;i<_grid.ncpl;i++){ //every seed lies in the box, so every cell must be a proper polygon
      ExitGracefullyIf(_grid.poly[i].size()<3,"CGroundwaterModel: HRU mesh produced a degenerate cell (internal error; please report the input)",RUNTIME_ERR);
    }
  }
  else if (!_extCellsLL.empty()){ //an existing model's cells, from MODFLOW memory and :MF6CRS (lon/lat)
    vector<gw_ring> cells(_extCellsLL.size());
    for (size_t i=0;i<_extCellsLL.size();i++){
      for (size_t n=0;n<_extCellsLL[i].size();n++){gw_pt q; _proj.Forward(_extCellsLL[i][n].x,_extCellsLL[i][n].y,q.x,q.y); cells[i].push_back(q);}
    }
    _grid.MakeFromPolygons(cells);
  }
  else { //imported cells
    map<long long,vector<gw_ring> > cellshapes;
    if (!GWGeom::ReadGeoJSONPolygons(_gridFile,_gridFileID,cellshapes,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
    vector<gw_ring> cells;
    for (map<long long,vector<gw_ring> >::iterator it=cellshapes.begin();it!=cellshapes.end();it++){
      ExitGracefullyIf(it->second.size()!=1,("CGroundwaterModel: grid cell "+to_string(it->first)+" in "+_gridFile+" must be a single polygon without holes").c_str(),BAD_DATA);
      gw_ring R; for (size_t n=0;n<it->second[0].size();n++){gw_pt q; _proj.Forward(it->second[0][n].x,it->second[0][n].y,q.x,q.y); R.push_back(q);}
      cells.push_back(R);
    }
    ExitGracefullyIf(cells.empty(),("CGroundwaterModel: no grid cells in "+_gridFile).c_str(),BAD_DATA);
    _grid.MakeFromPolygons(cells);
  }
  _ncpl=_grid.ncpl; _nrow=_grid.nrow; _ncol=_grid.ncol; _x0=_grid.x0; _y0=_grid.y0;
  ExitGracefullyIf(_ncpl>5000000,"CGroundwaterModel: grid exceeds 5 million cells per layer; increase :GridCellSize",BAD_DATA);

  //-- intersection of every HRU ring with every cell it may touch
  _hruCells.assign(nHRUs,map<int,double>());
  vector<int> cand;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i];
    for (size_t r=0;r<prj[k].size();r++)
    {
      const gw_ring &ring=prj[k][r];
      double rx0=1e30,ry0=1e30,rx1=-1e30,ry1=-1e30;
      for (size_t n=0;n<ring.size();n++){rx0=min(rx0,ring[n].x);rx1=max(rx1,ring[n].x);ry0=min(ry0,ring[n].y);ry1=max(ry1,ring[n].y);}
      _grid.Candidates(rx0,ry0,rx1,ry1,cand);
      for (size_t t=0;t<cand.size();t++){
        double A=_grid.ClipArea(ring,cand[t]);
        if (fabs(A)>1e-6){_hruCells[k][cand[t]]+=A;}
      }
    }
  }

  //-- write the linkage cache
  { ofstream CF(cachefile.c_str());
    CF<<key.str()<<endl<<setprecision(17)<<_proj.GetLat0()<<" "<<_proj.GetLon0()<<endl; _grid.Write(CF);
    CF<<_aqHRUs.size()<<endl;
    for (size_t i=0;i<_aqHRUs.size();i++){
      int k=_aqHRUs[i]; CF<<_pModel->GetHydroUnit(k)->GetHRUID()<<" "<<_hruPolyArea[k]<<" "<<_hruCells[k].size()<<endl;
      for (map<int,double>::iterator it=_hruCells[k].begin();it!=_hruCells[k].end();it++){CF<<it->first<<" "<<it->second<<endl;}
    }
  }
  } //end if !_cacheHit
  _geomSeconds=double(clock()-tclk)/CLOCKS_PER_SEC;
  //-- column cover, dominant profile, fallback land surface from HRU elevations
  vector<double> cover(_ncpl,0.0),elevW(_ncpl,0.0),best(_ncpl,0.0);
  _colProfile.assign(_ncpl,DOESNT_EXIST);
  _colDomHRU.assign(_ncpl,DOESNT_EXIST);
  for (size_t i=0;i<_aqHRUs.size();i++){
    int k=_aqHRUs[i];
    double elev=_pModel->GetHydroUnit(k)->GetElevation();
    for (map<int,double>::iterator it=_hruCells[k].begin();it!=_hruCells[k].end();it++){
      if (it->second<=0){continue;}
      cover[it->first]+=it->second; elevW[it->first]+=it->second*elev;
      if (it->second>best[it->first]){best[it->first]=it->second;_colProfile[it->first]=_hruProfile[k];_colDomHRU[it->first]=k;}
    }
  }
  CGWRaster &dem=_dem;
  if (_demFile!=""){ if (!dem.Read(_demFile,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);} }
  _land.assign(_ncpl,0.0);
  for (int ic=0;ic<_ncpl;ic++)
  {
    if (cover[ic]<_minCoverage*_grid.area[ic]){_colProfile[ic]=DOESNT_EXIST;}
    if (_colProfile[ic]==DOESNT_EXIST){continue;}
    bool ok=false;
    if (dem.IsLoaded()){_land[ic]=MeanRaster(dem,ic,ok);}
    if (!ok){_land[ic]=elevW[ic]/cover[ic];}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief builds model top, layer bottoms, bedrock clipping, pass-through and inactive cells, initial heads
//
void CGroundwaterModel::BuildVertical(const optStruct &Options)
{
  string err;
  CGWRaster bed;
  if (_bedrockFile!=""){ if (!bed.Read(_bedrockFile,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);} }
  _nlay=1;
  for (size_t i=0;i<_aProfiles.size();i++){_nlay=max(_nlay,(int)_aProfiles[i].layers.size());}

  _top.assign(_ncpl,0.0); _bedrock.assign(_ncpl,0.0);
  _botm.assign(_nlay*_ncpl,0.0); _strt.assign(_nlay*_ncpl,0.0); _idomain.assign(_nlay*_ncpl,0);
  int nThin=0;
  for (int ic=0;ic<_ncpl;ic++)
  {
    int ip=_colProfile[ic];
    if (ip!=DOESNT_EXIST)
    {
      const gw_profile &P=_aProfiles[ip];
      double land=_land[ic];
      double top =land-P.soil_zone_depth;
      double bedrock=land-P.default_bedrock_depth;
      bool ok=false;
      if (bed.IsLoaded()){double b=MeanRaster(bed,ic,ok); if (ok){bedrock=b;}}
      if (bedrock>top-GW_MIN_LAYER_THICK){_colProfile[ic]=DOESNT_EXIST; nThin++;} //bedrock at or above model top
      else
      {
        _top[ic]=top; _bedrock[ic]=bedrock;
        double prev=top,lowest=top;
        int nL=(int)P.layers.size();
        for (int l=0;l<_nlay;l++)
        {
          int n=UserNode(l,ic); double bot; int idom;
          if (l<nL)
          {
            const gw_layer &L=P.layers[l];
            bot=L.to_bedrock?bedrock:prev-L.thickness;
            if (bot<bedrock){bot=bedrock;}
            if (l==0){
              if (prev-bot<GW_MIN_LAYER_THICK){bot=prev-GW_MIN_LAYER_THICK;}
              idom=1;
            }
            else if (prev-bot<GW_MIN_LAYER_THICK){bot=prev-GW_PASS_THICK; idom=-1;}
            else {idom=(L.type==GWL_AQUICLUDE)?0:1;}
          }
          else {bot=prev-GW_PASS_THICK; idom=-1;}
          _botm[n]=bot; _idomain[n]=idom;
          if (idom>0){lowest=bot;}
          prev=bot;
        }
        double h0=land-P.init_head_depth;
        if      (_initMode=="DEPTH_BELOW_SURFACE"){h0=land-_initValue;}
        else if (_initMode=="ELEVATION"          ){h0=_initValue;}
        h0=max(h0,lowest+0.1);
        for (int l=0;l<_nlay;l++){_strt[UserNode(l,ic)]=h0;}
        continue;
      }
    }
    //inactive column
    _top[ic]=0.0;
    for (int l=0;l<_nlay;l++){int n=UserNode(l,ic);_botm[n]=-(double)(l+1);_idomain[n]=0;_strt[n]=0.0;}
  }
  if (_hotHeads.size()>0){ //hotstart: cell-by-cell heads from a previous run's solution.rvc
    int hr=(_gridType==GWGRID_REGULAR)?_nrow:1,hc=(_gridType==GWGRID_REGULAR)?_ncol:_ncpl;
    ExitGracefullyIf((_hotDims[0]!=_nlay) || (_hotDims[1]!=hr) || (_hotDims[2]!=hc) || ((int)_hotHeads.size()!=_nlay*_ncpl),
      "CGroundwaterModel: :GWHeads in the .rvc file does not match the groundwater grid (were the .rvg, .rvp or HRU geometry changed?)",BAD_DATA);
    for (int n=0;n<_nlay*_ncpl;n++){if (_idomain[n]!=0){_strt[n]=_hotHeads[n];}}
  }
  if (nThin>0){
    WriteAdvisory("CGroundwaterModel: "+to_string(nThin)+" grid columns deactivated because bedrock is at or above the model top",Options.noisy);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief reduced node numbering and the HRU<->cell linkage over active columns
//
void CGroundwaterModel::BuildLinks()
{
  _reduced.assign(_nlay*_ncpl,-1);
  int n=0;
  for (int i=0;i<_nlay*_ncpl;i++){if (_idomain[i]>0){_reduced[i]=n++;}}

  int nHRUs=_pModel->GetNumHRUs();
  _aLinks.clear(); _hruLinks.assign(nHRUs,vector<int>()); _hruLinkedArea.assign(nHRUs,0.0);
  _colCover.assign(_ncpl,0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i];
    for (map<int,double>::iterator it=_hruCells[k].begin();it!=_hruCells[k].end();it++){
      int ic=it->first;
      if ((it->second<=0) || (_colProfile[ic]==DOESNT_EXIST) || (_idomain[TopNode(ic)]<=0)){continue;}
      gw_link L; L.k=k; L.ic=ic; L.area=it->second;
      _hruLinks[k].push_back((int)_aLinks.size());
      _aLinks.push_back(L);
      _hruLinkedArea[k]+=L.area;
      _colCover[ic]+=L.area;
    }
    if (_hruLinkedArea[k]<=0){
      ExitGracefully(("CGroundwaterModel: HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+
        " does not overlap any active groundwater cell. "+(IsExternal()?
        string("Check that the HRU lies over active cells of the MODFLOW 6 model (projection, :OverlapWeights rows), or remove its AQUIFER_PROFILE."):
        string("Reduce :GridCellSize or :MinCellCoverage, or check bedrock depth."))).c_str(),BAD_DATA);
    }
  }
  _bndCol.clear();
  for (int ic=0;ic<_ncpl;ic++){if ((_colProfile[ic]!=DOESNT_EXIST) && (_colCover[ic]>0)){_bndCol.push_back(ic);}}
  _evtCol.clear();
  for (size_t i=0;i<_bndCol.size();i++){if (_aProfiles[_colProfile[_bndCol[i]]].extinction_depth>0){_evtCol.push_back(_bndCol[i]);}}
  _Qcol.assign(_ncpl,0.0); _Qseep.assign(_ncpl,0.0);
}

//////////////////////////////////////////////////////////////////
/// \brief river cells from river lines: length per cell, channel bed elevation (lowest DEM along the
///  line in the cell), subbasin (ID field, or subbasin of the dominant HRU in the cell), layer holding the bed
//
void CGroundwaterModel::BuildRivers()
{
  _sbGain.assign(_pModel->GetNumSubBasins(),0.0);
  _sbLossRemain.assign(_pModel->GetNumSubBasins(),0.0);
  _sbLossCarry.assign(_pModel->GetNumSubBasins(),0.0); //unsupplied losses carried from a previous run (hotstart)
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    map<long long,double>::iterator it=_hotRivCarry.find(_pModel->GetSubBasin(p)->GetID());
    if (it!=_hotRivCarry.end()){_sbLossCarry[p]=it->second;}
  }
  if (_riverFile==""){return;}
  vector<gw_ring> lines; vector<long long> ids; string err;
  if (!GWGeom::ReadGeoJSONLines(_riverFile,_riverIDField,lines,ids,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
  map<pair<int,int>,pair<double,double> > acc; //(subbasin,column) -> (length, min bed elevation)
  double cs=_cellSize;
  for (size_t i=0;i<lines.size();i++)
  {
    int pfix=DOESNT_EXIST;
    if (ids[i]>=0){pfix=_pModel->GetSubBasinIndex(ids[i]); if (pfix==DOESNT_EXIST){continue;}}
    for (size_t n=0;n+1<lines[i].size();n++)
    {
      const gw_pt &a=lines[i][n],&b=lines[i][n+1];
      double xa,ya,xb,yb; _proj.Forward(a.x,a.y,xa,ya); _proj.Forward(b.x,b.y,xb,yb);
      double len=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));
      if (!_grid.uniform){ //exact pieces of the segment in each cell
        vector<gw_piece> pcs; SegmentPieces(xa,ya,xb,yb,pcs);
        for (size_t q=0;q<pcs.size();q++){
          int ic=pcs[q].ic;
          if ((_colProfile[ic]==DOESNT_EXIST) || (_colCover[ic]<=0)){continue;}
          int p=pfix;
          if (p==DOESNT_EXIST){p=_pModel->GetHydroUnit(_colDomHRU[ic])->GetSubBasinIndex();}
          if (!_pModel->GetSubBasin(p)->IsEnabled()){continue;}
          double plen=(pcs[q].t1-pcs[q].t0)*len,z=_land[ic],val;
          int ns=max(1,(int)ceil(plen/(0.2*sqrt(_grid.area[ic]))));
          for (int s=0;s<ns;s++){double w=pcs[q].t0+(s+0.5)/ns*(pcs[q].t1-pcs[q].t0);
            if (_dem.IsLoaded() && _dem.Sample(a.x+w*(b.x-a.x),a.y+w*(b.y-a.y),val)){z=(s==0)?val:min(z,val);}}
          pair<int,int> key(p,ic);
          if (acc.find(key)==acc.end()){acc[key]=make_pair(0.0,z);}
          acc[key].first+=plen; acc[key].second=min(acc[key].second,z);
        }
        continue;
      }
      int npc=max(1,(int)ceil(len/(cs/20.0)));
      for (int j=0;j<npc;j++){
        double w=(j+0.5)/npc, x=xa+w*(xb-xa), y=ya+w*(yb-ya);
        int ic=LocateActiveOnLine(x,y,xb-xa,yb-ya);
        if (ic<0){continue;}
        if (_colCover[ic]<=0){continue;}
        int p=pfix;
        if (p==DOESNT_EXIST){p=_pModel->GetHydroUnit(_colDomHRU[ic])->GetSubBasinIndex();}
        if (!_pModel->GetSubBasin(p)->IsEnabled()){continue;} //disabled subbasins have no river to exchange with
        double z=_land[ic],val;
        if (_dem.IsLoaded() && _dem.Sample(a.x+w*(b.x-a.x),a.y+w*(b.y-a.y),val)){z=val;}
        pair<int,int> key(p,ic);
        if (acc.find(key)==acc.end()){acc[key]=make_pair(0.0,z);}
        acc[key].first+=len/npc;
        acc[key].second=min(acc[key].second,z);
      }
    }
  }
  for (map<pair<int,int>,pair<double,double> >::iterator it=acc.begin();it!=acc.end();it++)
  {
    int p=it->first.first,ic=it->first.second;
    //the DEM records the water surface/floodplain along the river: channel bed sits one bankfull depth
    //(depth at the Raven reference flow) below it
    CSubBasin *pSB=_pModel->GetSubBasin(p);
    double dbf=(_rivBedDepth>=0)?_rivBedDepth:max(0.0,min(pSB->GetRiverDepthAtFlow(pSB->GetReferenceFlow()),20.0));
    double zbed=min(it->second.second-dbf,_top[ic]);
    int lay=DOESNT_EXIST,lowest=DOESNT_EXIST;
    for (int l=0;l<_nlay;l++){
      int n=UserNode(l,ic); if (_idomain[n]<=0){continue;}
      lowest=l;
      if ((lay==DOESNT_EXIST) && (_botm[n]<zbed)){lay=l;}
    }
    if (lowest==DOESNT_EXIST){continue;}
    if (lay==DOESNT_EXIST){lay=lowest; zbed=_botm[UserNode(lowest,ic)]+0.1;}
    _rivNode.push_back(UserNode(lay,ic)); _rivCol.push_back(ic); _rivSB.push_back(p);
    _rivLen.push_back(it->second.first); _rivBed.push_back(zbed);
  }
  _sbRivLen.assign(_pModel->GetNumSubBasins(),0.0);
  for (size_t e=0;e<_rivNode.size();e++){_sbRivLen[_rivSB[e]]+=_rivLen[e];}
}

//////////////////////////////////////////////////////////////////
/// \brief active columns selected by a GeoJSON file: cell centres inside its polygons, or cells crossed by its lines
//
vector<int> CGroundwaterModel::SelectColumns(const string &file) const
{
  vector<int> cols; vector<char> sel(_ncpl,0); string err;
  map<long long,vector<gw_ring> > polys;
  if (GWGeom::ReadGeoJSONPolygons(file,"",polys,err))
  {
    vector<gw_ring> rings;
    for (map<long long,vector<gw_ring> >::const_iterator it=polys.begin();it!=polys.end();it++){
      for (size_t r=0;r<it->second.size();r++){
        gw_ring R; for (size_t n=0;n<it->second[r].size();n++){gw_pt q; _proj.Forward(it->second[r][n].x,it->second[r][n].y,q.x,q.y); R.push_back(q);}
        rings.push_back(R);
      }
    }
    for (int ic=0;ic<_ncpl;ic++){
      if (_colProfile[ic]==DOESNT_EXIST){continue;}
      double x,y; CellCentreWorld(ic,x,y);
      bool in=false; //even-odd rule over all rings (holes included)
      for (size_t r=0;r<rings.size();r++){
        const gw_ring &R=rings[r]; size_t m=R.size();
        for (size_t i=0,j=m-1;i<m;j=i++){
          if (((R[i].y>y)!=(R[j].y>y)) && (x<(R[j].x-R[i].x)*(y-R[i].y)/(R[j].y-R[i].y)+R[i].x)){in=!in;}
        }
      }
      if (in){sel[ic]=1;}
    }
  }
  else
  {
    vector<gw_ring> lines; vector<long long> ids;
    if (!GWGeom::ReadGeoJSONLines(file,"",lines,ids,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
    for (size_t i=0;i<lines.size();i++){
      for (size_t n=0;n+1<lines[i].size();n++){
        double xa,ya,xb,yb; _proj.Forward(lines[i][n].x,lines[i][n].y,xa,ya); _proj.Forward(lines[i][n+1].x,lines[i][n+1].y,xb,yb);
        if (!_grid.uniform){ //exact pieces; pieces outside the active area attach to the nearest active cell
          vector<gw_piece> pcs; SegmentPieces(xa,ya,xb,yb,pcs);
          for (size_t q=0;q<pcs.size();q++){
            int ic=pcs[q].ic;
            if (_colProfile[ic]==DOESNT_EXIST){double w=0.5*(pcs[q].t0+pcs[q].t1); ic=NearestActiveCell(xa+w*(xb-xa),ya+w*(yb-ya));}
            if (ic>=0){sel[ic]=1;}
          }
          continue;
        }
        double step=_cellSize/20.0;
        int npc=max(1,(int)ceil(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))/step));
        for (int j=0;j<npc;j++){
          double w=(j+0.5)/npc,X=xa+w*(xb-xa),Y=ya+w*(yb-ya); int ic=LocateActiveOnLine(X,Y,xb-xa,yb-ya);
          if (ic<0){ic=NearestActiveCell(X,Y);} //boundary drawn along the edge of the active area
          if (ic>=0){sel[ic]=1;}
        }
      }
    }
  }
  for (int ic=0;ic<_ncpl;ic++){if (sel[ic]){cols.push_back(ic);}}
  return cols;
}

//////////////////////////////////////////////////////////////////
/// \brief wells (screen split over layers by transmissivity), monitoring wells, GHB/CHD boundary cells
//
void CGroundwaterModel::BuildStresses()
{
  for (size_t a=0;a<_aWells.size();a++){for (size_t b=a+1;b<_aWells.size();b++){
    ExitGracefullyIf(_aWells[a].id==_aWells[b].id,("CGroundwaterModel: two wells have ID "+to_string(_aWells[a].id)).c_str(),BAD_DATA);}}
  for (size_t a=0;a<_aObsWells.size();a++){for (size_t b=a+1;b<_aObsWells.size();b++){
    ExitGracefullyIf(_aObsWells[a].id==_aObsWells[b].id,("CGroundwaterModel: two observation wells have ID "+to_string(_aObsWells[a].id)).c_str(),BAD_DATA);}}
  for (size_t a=0;a<_aBounds.size();a++){for (size_t b=a+1;b<_aBounds.size();b++){
    ExitGracefullyIf(_aBounds[a].name==_aBounds[b].name,("CGroundwaterModel: two boundaries are named "+_aBounds[a].name).c_str(),BAD_DATA);}}
  for (size_t w=0;w<_aWells.size();w++)
  {
    gw_well &W=_aWells[w];
    map<long long,CTimeSeries*>::iterator it=_wellTS.find(W.id);
    W.pRate=(it==_wellTS.end())?NULL:it->second;
    if (W.pRate==NULL){WriteWarning("CGroundwaterModel: well "+W.name+" has no :WellRate time series in the .rvt file; rate set to zero",false);}
    double x,y; _proj.Forward(W.lon,W.lat,x,y);
    int ic=LocateWorld(x,y);
    ExitGracefullyIf((ic<0) || (_colProfile[ic]==DOESNT_EXIST),("CGroundwaterModel: well "+W.name+" is outside the active groundwater model").c_str(),BAD_DATA);
    vector<int> nodes; vector<double> T; double Tsum=0; int lowest=-1;
    for (int l=0;l<_nlay;l++){
      int n=UserNode(l,ic); if (_idomain[n]<=0){continue;}
      lowest=n;
      double ov=min(W.screen_top,CellTop(l,ic))-max(W.screen_bot,_botm[n]);
      if (ov>0){
        const gw_profile &P=_aProfiles[_colProfile[ic]];
        double K=(l<(int)P.layers.size())?_aClasses[P.layers[l].iclass].Kh:1.0;
        nodes.push_back(n); T.push_back(K*ov); Tsum+=K*ov;
      }
    }
    if (nodes.size()==0){ //screen outside the aquifer: use the nearest active layer
      nodes.push_back((W.screen_bot>=_top[ic])?TopNode(ic):lowest); T.push_back(1.0); Tsum=1.0;
      WriteWarning("CGroundwaterModel: screen of well "+W.name+" lies outside the aquifer; assigned to the nearest active layer",false);
    }
    for (size_t j=0;j<nodes.size();j++){_welNode.push_back(nodes[j]); _welWell.push_back((int)w); _welFrac.push_back(T[j]/Tsum);}
    for (size_t e=0;e<_aBounds.size();e++){ //MODFLOW ignores wells in fixed-head cells
      if (!_aBounds[e].isCHD){continue;}
      vector<int> cc=SelectColumns(_aBounds[e].file);
      if (std::find(cc.begin(),cc.end(),ic)!=cc.end()){WriteWarning("CGroundwaterModel: well "+W.name+" lies in a cell of fixed-head boundary "+_aBounds[e].name+"; MODFLOW holds the head there and the pumping has no effect",false);}
    }
  }
  for (size_t o=0;o<_aObsWells.size();o++)
  {
    gw_obswell &O=_aObsWells[o]; O.node=-1;
    double x,y; _proj.Forward(O.lon,O.lat,x,y);
    int ic=LocateWorld(x,y);
    if ((ic<0) || (_colProfile[ic]==DOESNT_EXIST)){
      WriteWarning("CGroundwaterModel: observation well "+O.name+" is outside the active groundwater model (no simulated heads)",false); continue;
    }
    int best=-1; double bestd=1e30;
    for (int l=0;l<_nlay;l++){ //layer holding the screen, else the nearest active layer
      int n=UserNode(l,ic); if (_idomain[n]<=0){continue;}
      double top=CellTop(l,ic),bot=_botm[n];
      double d=(O.z>top)?O.z-top:((O.z<bot)?bot-O.z:0.0);
      if (d<bestd){bestd=d;best=n;}
    }
    O.node=best;
  }
  for (size_t b=0;b<_aBounds.size();b++)
  {
    gw_boundary &B=_aBounds[b];
    map<string,CTimeSeries*>::iterator it=_headTS.find(B.name);
    B.pHead=(it==_headTS.end())?NULL:it->second;
    vector<int> cols=SelectColumns(B.file);
    if (cols.size()==0){WriteWarning("CGroundwaterModel: boundary "+B.name+" selects no active cells",false);}
    for (size_t i=0;i<cols.size();i++){
      for (int l=0;l<_nlay;l++){ //MODFLOW requires boundary heads above the cell bottom: skip layers the head lies below
        int n=UserNode(l,cols[i]); if ((_idomain[n]<=0) || (_botm[n]>=B.head)){continue;}
        if (B.isCHD){_chdNode.push_back(n);_chdB.push_back((int)b);} else {_ghbNode.push_back(n);_ghbB.push_back((int)b);}
        if (!B.allLayers){break;}
      }
    }
  }
}


//////////////////////////////////////////////////////////////////
/// \brief XT3D (full conductance tensor) for grids whose connections are not orthogonal
//
bool CGroundwaterModel::UseXT3D() const
{
  if ((_flowCorr=="XT3D") || (_flowCorr=="XT3D_RHS")){return true;}
  if (_flowCorr=="NONE"){return false;}
  return false; //AUTO: no correction - XT3D (even its RHS form) slows Newton badly where cells dry and rewet
}
//////////////////////////////////////////////////////////////////
/// \brief writes every cell of the grid (lon/lat polygons) to mf6/grid_cells.geojson; computes the HRU fidelity of the cells
//
void CGroundwaterModel::WriteGridCells(const optStruct &Options)
{
  double wsum=0,wbest=0;
  for (int ic=0;ic<_ncpl;ic++){
    if ((_colProfile[ic]==DOESNT_EXIST) || (_colCover[ic]<=0)){continue;}
    double best=0; for (size_t i=0;i<_aqHRUs.size();i++){map<int,double>::const_iterator it=_hruCells[_aqHRUs[i]].find(ic); if (it!=_hruCells[_aqHRUs[i]].end()){best=max(best,it->second);}}
    wbest+=best; wsum+=_colCover[ic];
  }
  _meshFidelity=(wsum>0)?wbest/wsum:1.0;
  ofstream G((Options.output_dir+"mf6/grid_cells.geojson").c_str());
  if (!G){return;}
  G<<"{\"type\":\"FeatureCollection\",\"features\":["<<setprecision(10);
  for (int ic=0;ic<_ncpl;ic++){
    G<<((ic>0)?",":"")<<"\n{\"type\":\"Feature\",\"properties\":{\"cell\":"<<ic+1<<",\"active\":"<<((_colProfile[ic]!=DOESNT_EXIST)?1:0);
    if (_colProfile[ic]!=DOESNT_EXIST){G<<",\"profile\":\""<<_aProfiles[_colProfile[ic]].name<<"\",\"hru\":"<<_pModel->GetHydroUnit(_colDomHRU[ic])->GetHRUID();}
    G<<",\"area_m2\":"<<_grid.area[ic]<<"},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[";
    const gw_ring &P=_grid.poly[ic];
    for (size_t k=0;k<=P.size();k++){
      double X,Y,lon,lat; _grid.ToWorld(P[k%P.size()].x,P[k%P.size()].y,X,Y); _proj.Inverse(X,Y,lon,lat);
      G<<((k>0)?",":"")<<"["<<lon<<","<<lat<<"]";
    }
    G<<"]]}}";
  }
  G<<"\n]}\n";
}

//////////////////////////////////////////////////////////////////
/// \brief writes a standard MODFLOW 6 simulation to <output>/mf6/
//
void CGroundwaterModel::WriteMF6Files(const optStruct &Options)
{
  string dir=Options.output_dir+"mf6"; //= main output dir, or the ensemble member's folder
  GW_MKDIR(Options.output_dir.c_str());
  GW_MKDIR(dir.c_str());
  _mf6dir=AbsPath(dir);
  string d=_mf6dir+"/";
  ofstream F;

  int    nstp=(int)floor(Options.duration/Options.timestep+0.5);
  double perlen=nstp*Options.timestep;

  F.open((d+"mfsim.nam").c_str());
  F<<"# MODFLOW 6 simulation generated by Raven (Raven-MODFLOW 6 groundwater coupling)"<<endl;
  F<<"BEGIN OPTIONS\n  CONTINUE\nEND OPTIONS\n\nBEGIN TIMING\n  TDIS6  sim.tdis\nEND TIMING\n\n";
  F<<"BEGIN MODELS\n  GWF6  gwf.nam  GWF\nEND MODELS\n\nBEGIN EXCHANGES\nEND EXCHANGES\n\n";
  F<<"BEGIN SOLUTIONGROUP 1\n  IMS6  sim.ims  GWF\nEND SOLUTIONGROUP\n"; F.close();

  F.open((d+"sim.tdis").c_str());
  F<<"BEGIN OPTIONS\n  TIME_UNITS DAYS\nEND OPTIONS\n\n";
  F<<"BEGIN DIMENSIONS\n  NPER "<<(_doSteady?2:1)<<"\nEND DIMENSIONS\n\nBEGIN PERIODDATA\n";
  if (_doSteady){F<<"  1.0 1 1.0   # steady-state initialization (solved before the Raven run starts)\n";}
  F<<"  "<<setprecision(12)<<perlen<<" "<<nstp<<" 1.0\nEND PERIODDATA\n"; F.close();

  F.open((d+"sim.ims").c_str());
  F<<"BEGIN OPTIONS\n  PRINT_OPTION SUMMARY\n  COMPLEXITY COMPLEX\nEND OPTIONS\n\n";
  //robust Newton settings for cells that dry and rewet: delta-bar-delta under-relaxation, backtracking,
  //deeper ILU preconditioning. Users can raise :SolverMaxOuterIterations in the .rvg file.
  F<<"BEGIN NONLINEAR\n  OUTER_DVCLOSE "<<_dvclose<<"\n  OUTER_MAXIMUM "<<_maxiter<<"\n";
  F<<"  UNDER_RELAXATION DBD\n  UNDER_RELAXATION_THETA 0.7\n  UNDER_RELAXATION_KAPPA 0.07\n  UNDER_RELAXATION_GAMMA 0.1\n  UNDER_RELAXATION_MOMENTUM 0.0\n";
  F<<"  BACKTRACKING_NUMBER 20\n  BACKTRACKING_TOLERANCE 1.05\n  BACKTRACKING_REDUCTION_FACTOR 0.1\n  BACKTRACKING_RESIDUAL_LIMIT 0.002\nEND NONLINEAR\n\n";
  F<<"BEGIN LINEAR\n  INNER_MAXIMUM 500\n  INNER_DVCLOSE "<<0.1*_dvclose<<"\n  INNER_RCLOSE 1.0\n  LINEAR_ACCELERATION BICGSTAB\n"; //inner tolerance tighter than the outer one
  F<<"  PRECONDITIONER_LEVELS 5\n  PRECONDITIONER_DROP_TOLERANCE 0.0001\n  NUMBER_ORTHOGONALIZATIONS 2\nEND LINEAR\n"; F.close();

  F.open((d+"gwf.nam").c_str());
  //Newton without its UNDER_RELAXATION option: that option resets heads that fall below a cell bottom, and in MODFLOW 6.6
  //the Newton change of such a cell then repeats unchanged every iteration, so the step never converges (a drying aquitard of
  //the Liard model failed on most days with 6.6.3); without it 6.6.3 and 6.8.1 converge in the same few iterations
  F<<"BEGIN OPTIONS\n  NEWTON\nEND OPTIONS\n\nBEGIN PACKAGES\n";
  F<<"  "<<((_disType==0)?"DIS6   gwf.dis ":(_disType==1)?"DISV6  gwf.disv":"DISU6  gwf.disu")<<"  DIS\n  NPF6  gwf.npf  NPF\n  STO6  gwf.sto  STO\n  IC6   gwf.ic   IC\n";
  F<<"  RCH6  gwf.rch  RCH_RAVEN\n  DRN6  gwf.drn  DRN_RAVEN\n";
  if (_rivNode.size()>0){F<<"  RIV6  gwf.riv  RIV_RAVEN\n";}
  if (_evtCol.size()>0) {F<<"  EVT6  gwf.evt  EVT_RAVEN\n";}
  if (_welNode.size()>0){F<<"  WEL6  gwf.wel  WEL_RAVEN\n";}
  if (_ghbNode.size()>0){F<<"  GHB6  gwf.ghb  GHB_RAVEN\n";}
  if (_chdNode.size()>0){F<<"  CHD6  gwf.chd  CHD_RAVEN\n";}
  if (_swxCol.size()>0) {F<<"  RIV6  gwf.swx  SWX_RAVEN\n";}
  if (_resNode.size()>0){F<<"  WEL6  gwf.res  RES_RAVEN\n";}
  F<<"  OC6   gwf.oc   OC\nEND PACKAGES\n"; F.close();

  //-- discretization: DIS (regular grids), DISV (quadtree), DISU (meshes, imported cells, connections by elevation)
  const char *disname[3]={"gwf.dis","gwf.disv","gwf.disu"};
  F.open((d+disname[_disType]).c_str());
  F<<"# Lambert azimuthal equal-area projection, origin lat="<<setprecision(10)<<_proj.GetLat0()<<" lon="<<_proj.GetLon0()<<" (sphere R=6371007.181 m)"<<endl;
  if (_disType==0)
  {
    double X0,Y0; _grid.ToWorld(_grid.colX[0],_grid.rowY[_nrow-1],X0,Y0); //lower-left corner of the grid
    bool orig=_grid.uniform && (_grid.angle==0.0); //unrotated uniform grid: the original arithmetic
    F<<"BEGIN OPTIONS\n  LENGTH_UNITS METERS\n  XORIGIN "<<(orig?_x0:X0)<<"\n  YORIGIN "<<(orig?_y0:Y0)<<"\n";
    if (_grid.angle!=0.0){F<<"  ANGROT "<<_grid.angle<<"\n";}
    F<<"END OPTIONS\n\n";
    F<<"BEGIN DIMENSIONS\n  NLAY "<<_nlay<<"\n  NROW "<<_nrow<<"\n  NCOL "<<_ncol<<"\nEND DIMENSIONS\n\n";
    if (_grid.uniform){F<<"BEGIN GRIDDATA\n  DELR\n    CONSTANT "<<_cellSize<<"\n  DELC\n    CONSTANT "<<_cellSize<<endl;}
    else {F<<"BEGIN GRIDDATA\n"; WriteArray(F,"DELR",_grid.colW,10); WriteArray(F,"DELC",_grid.rowH,10);}
    WriteArray(F,"TOP",_top,_ncol); WriteArray(F,"BOTM",_botm,_ncol,_nlay); WriteArray(F,"IDOMAIN",_idomain,_ncol,_nlay,"1");
    F<<"END GRIDDATA\n"; F.close();
  }
  else if (_disType==1)
  {
    F<<"BEGIN OPTIONS\n  LENGTH_UNITS METERS\n  XORIGIN 0.0\n  YORIGIN 0.0\n";
    if (_grid.angle!=0.0){F<<"  ANGROT "<<_grid.angle<<"\n";}
    F<<"END OPTIONS\n\nBEGIN DIMENSIONS\n  NLAY "<<_nlay<<"\n  NCPL "<<_ncpl<<"\n  NVERT "<<_grid.verts.size()<<"\nEND DIMENSIONS\n\n";
    F<<"BEGIN GRIDDATA\n"; WriteArray(F,"TOP",_top,10); WriteArray(F,"BOTM",_botm,10,_nlay); WriteArray(F,"IDOMAIN",_idomain,10,_nlay,"1");
    F<<"END GRIDDATA\n\nBEGIN VERTICES\n"<<setprecision(15);
    for (size_t v=0;v<_grid.verts.size();v++){F<<"  "<<v+1<<" "<<_grid.verts[v].x<<" "<<_grid.verts[v].y<<"\n";}
    F<<"END VERTICES\n\nBEGIN CELL2D\n";
    for (int ic=0;ic<_ncpl;ic++){
      F<<"  "<<ic+1<<" "<<_grid.cx[ic]<<" "<<_grid.cy[ic]<<" "<<_grid.cverts[ic].size();
      for (size_t k=0;k<_grid.cverts[ic].size();k++){F<<" "<<_grid.cverts[ic][k]+1;}
      F<<"\n";
    }
    F<<"END CELL2D\n"; F.close();
  }
  else
  {
    //nodes keep the layered numbering (layer*ncpl+cell); pass-through cells (IDOMAIN -1) are bridged vertically,
    //aquicludes (IDOMAIN 0) are not; horizontal connections join cells of the same layer, or (:LayerConnection ELEVATION)
    //any cells of neighbouring columns whose vertical intervals overlap (staggered connections, IHC 2)
    int nodes=_nlay*_ncpl;
    vector<double> topN(nodes),botN(nodes),areaN(nodes); vector<int> idN(nodes);
    for (int l=0;l<_nlay;l++){for (int ic=0;ic<_ncpl;ic++){int n=UserNode(l,ic);
      topN[n]=(l==0)?_top[ic]:_botm[n-_ncpl]; botN[n]=_botm[n]; areaN[n]=_grid.area[ic]; idN[n]=(_idomain[n]>0)?1:0;}}
    struct cn {int m,ihc; double cl1,cl2,hw,ang;};
    vector<vector<cn> > C(nodes);
    for (int ic=0;ic<_ncpl;ic++){
      for (int l=0;l<_nlay;l++){
        int n=UserNode(l,ic); if (!idN[n]){continue;}
        for (int l2=l+1;l2<_nlay;l2++){ //vertical: next active cell below in the same column, bridging pass-through
          int m=UserNode(l2,ic);                  //cells (IDOMAIN -1) only; an aquiclude (IDOMAIN 0) blocks the column,
          if (_idomain[m]==0){break;}             //as MODFLOW does for DIS and DISV grids
          if (!idN[m]){continue;}
          cn a; a.m=m; a.ihc=0; a.cl1=0.5*(topN[n]-botN[n]); a.cl2=0.5*(topN[m]-botN[m]); a.hw=areaN[n]; a.ang=0; C[n].push_back(a);
          cn b=a; b.m=n; swap(b.cl1,b.cl2); C[m].push_back(b); break;
        }
      }
      for (size_t k=0;k<_grid.nbr[ic].size();k++){
        const gw_nbr &g=_grid.nbr[ic][k]; if (g.j<ic){continue;}
        for (int l=0;l<_nlay;l++){
          int n=UserNode(l,ic); if (!idN[n]){continue;}
          for (int l2=0;l2<_nlay;l2++){
            int m=UserNode(l2,g.j); if (!idN[m]){continue;}
            bool link=(l2==l);
            if (_connectByElevation){link=(min(topN[n],topN[m])-max(botN[n],botN[m])>1e-6);}
            if (!link){continue;}
            //ANGLDEGX is measured in the frame of the vertices, which are written in the grid's own (rotated) frame;
            //ANGROT only places that frame on the map
            cn a; a.m=m; a.ihc=((l2==l) && !_connectByElevation)?1:2; a.cl1=g.di; a.cl2=g.dj; a.hw=g.len; a.ang=fmod(g.ang*180.0/PI+720.0,360.0);
            C[n].push_back(a);
            cn b=a; b.m=n; swap(b.cl1,b.cl2); b.ang=fmod(a.ang+180.0,360.0); C[m].push_back(b);
          }
        }
      }
    }
    struct byNode { bool operator()(const cn &a,const cn &b) const {return a.m<b.m;} };
    for (int n=0;n<nodes;n++){sort(C[n].begin(),C[n].end(),byNode());} //MODFLOW pairs connections n-m and m-n by position
    long long nja=0; for (int n=0;n<nodes;n++){nja+=1+(long long)C[n].size();}
    F<<"BEGIN OPTIONS\n  LENGTH_UNITS METERS\n";
    if (_grid.angle!=0.0){F<<"  ANGROT "<<_grid.angle<<"\n";}
    F<<"END OPTIONS\n\nBEGIN DIMENSIONS\n  NODES "<<nodes<<"\n  NJA "<<nja<<"\n  NVERT "<<_grid.verts.size()<<"\nEND DIMENSIONS\n\n";
    F<<"BEGIN GRIDDATA\n"; WriteArray(F,"TOP",topN,10); WriteArray(F,"BOT",botN,10); WriteArray(F,"AREA",areaN,10); WriteArray(F,"IDOMAIN",idN,10,1,"1");
    F<<"END GRIDDATA\n\nBEGIN CONNECTIONDATA\n  IAC\n    INTERNAL\n"<<setprecision(15);
    for (int n=0;n<nodes;n++){F<<((n%10==0)?"      ":" ")<<1+C[n].size(); if ((n%10==9) || (n==nodes-1)){F<<"\n";}}
    const char *nm[5]={"JA","IHC","CL12","HWVA","ANGLDEGX"};
    for (int q=0;q<5;q++){
      F<<"  "<<nm[q]<<"\n    INTERNAL\n";
      for (int n=0;n<nodes;n++){
        F<<"     ";
        if      (q==0){F<<" "<<n+1;} else if (q==1){F<<" 1";} else {F<<" 0";}
        for (size_t k=0;k<C[n].size();k++){
          const cn &a=C[n][k];
          if      (q==0){F<<" "<<a.m+1;}
          else if (q==1){F<<" "<<a.ihc;}
          else if (q==2){F<<" "<<a.cl1;}
          else if (q==3){F<<" "<<a.hw;}
          else          {F<<" "<<((a.ihc==0)?0.0:fmod(a.ang+720.0,360.0));}
        }
        F<<"\n";
      }
    }
    F<<"END CONNECTIONDATA\n\nBEGIN VERTICES\n"; //cell outlines: needed for face directions (XT3D, specific discharge)
    for (size_t v=0;v<_grid.verts.size();v++){F<<"  "<<v+1<<" "<<_grid.verts[v].x<<" "<<_grid.verts[v].y<<"\n";}
    F<<"END VERTICES\n\nBEGIN CELL2D\n";
    for (int n=0;n<nodes;n++){
      int ic=n%_ncpl; F<<"  "<<n+1<<" "<<_grid.cx[ic]<<" "<<_grid.cy[ic]<<" "<<_grid.cverts[ic].size();
      for (size_t k=0;k<_grid.cverts[ic].size();k++){F<<" "<<_grid.cverts[ic][k]+1;}
      F<<"\n";
    }
    F<<"END CELL2D\n"; F.close();
  }

  int aC=(_disType==0)?_ncol:10,aL=(_disType==2)?1:_nlay; //array layout
  //-- cell properties from each column's profile layer
  vector<int> icelltype(_nlay*_ncpl,1),iconvert(_nlay*_ncpl,1);
  vector<double> K(_nlay*_ncpl,1.0),K33(_nlay*_ncpl,1.0),Ss(_nlay*_ncpl,1e-5),Sy(_nlay*_ncpl,0.1);
  for (int ic=0;ic<_ncpl;ic++){
    int ip=_colProfile[ic]; if (ip==DOESNT_EXIST){continue;}
    const gw_profile &P=_aProfiles[ip];
    for (int l=0;l<(int)P.layers.size();l++){
      int n=UserNode(l,ic);
      const gw_aquifer_class &C=_aClasses[P.layers[l].iclass];
      K[n]=C.Kh; K33[n]=C.Kv; Ss[n]=C.Ss; Sy[n]=C.Sy;
      bool conf=(P.layers[l].type==GWL_CONFINED_AQUIFER);
      icelltype[n]=conf?0:1; iconvert[n]=conf?0:1;
    }
  }
  F.open((d+"gwf.npf").c_str());
  F<<"BEGIN OPTIONS\n  SAVE_FLOWS\n"<<(UseXT3D()?((_flowCorr=="XT3D")?"  XT3D\n":"  XT3D RHS\n"):"")<<((_disType==2)?"  SAVE_SPECIFIC_DISCHARGE\n":"")<<"END OPTIONS\n\nBEGIN GRIDDATA\n";
  WriteArray(F,"ICELLTYPE",icelltype,aC,aL,"1"); WriteArray(F,"K",K,aC,aL); WriteArray(F,"K33",K33,aC,aL);
  F<<"END GRIDDATA\n"; F.close();

  F.open((d+"gwf.sto").c_str());
  F<<"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN GRIDDATA\n";
  WriteArray(F,"ICONVERT",iconvert,aC,aL,"1"); WriteArray(F,"SS",Ss,aC,aL); WriteArray(F,"SY",Sy,aC,aL);
  if (_doSteady){F<<"END GRIDDATA\n\nBEGIN PERIOD 1\n  STEADY-STATE\nEND PERIOD\n\nBEGIN PERIOD 2\n  TRANSIENT\nEND PERIOD\n";}
  else            {F<<"END GRIDDATA\n\nBEGIN PERIOD 1\n  TRANSIENT\nEND PERIOD\n";}
  F.close();

  F.open((d+"gwf.ic").c_str());
  F<<"BEGIN GRIDDATA\n"; WriteArray(F,"STRT",_strt,aC,aL); F<<"END GRIDDATA\n"; F.close();

  //-- recharge and seepage-face drains on the top cell of every linked column
  F.open((d+"gwf.rch").c_str());
  F<<"# rates are overwritten by Raven every time step\nBEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_bndCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
  for (size_t i=0;i<_bndCol.size();i++){int ic=_bndCol[i]; F<<"  "<<CellID(0,ic)<<" 0.0"<<endl;}
  F<<"END PERIOD\n"; F.close();

  F.open((d+"gwf.drn").c_str());
  F<<"# seepage faces at the model top; discharge returns to Raven. DDRN<0: discharge starts DDRN below ELEV and\n";
  F<<"# reaches full conductance at ELEV (smooth onset helps Newton convergence)\n";
  F<<"BEGIN OPTIONS\n  AUXILIARY DDRN\n  AUXDEPTHNAME DDRN\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_bndCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
  F<<setprecision(10);
  for (size_t i=0;i<_bndCol.size();i++){
    int ic=_bndCol[i];
    double cond=_aProfiles[_colProfile[ic]].seep_leakance*_grid.area[ic];
    //MODFLOW refuses a drain whose smoothing interval reaches below the cell bottom: limit it to the top layer
    double dsm=min(_aProfiles[_colProfile[ic]].seep_smooth_depth,_top[ic]-_botm[UserNode(0,ic)]);
    F<<"  "<<CellID(0,ic)<<" "<<_top[ic]<<" "<<cond<<" "<<-dsm<<endl;
  }
  F<<"END PERIOD\n"; F.close();

  if (_rivNode.size()>0){
    F.open((d+"gwf.riv").c_str());
    F<<"# river cells linked to Raven reaches; stage and conductance are set by Raven every time step\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_rivNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_rivNode.size();e++){
      int n=_rivNode[e],l=n/_ncpl,ic=n%_ncpl;
      F<<"  "<<CellID(l,ic)<<" "<<_rivBed[e]<<" 0.0 "<<_rivBed[e]<<endl;
    }
    F<<"END PERIOD\n"; F.close();
  }
  if (_evtCol.size()>0){
    F.open((d+"gwf.evt").c_str());
    F<<"# water-table evapotranspiration; rate = PET left unmet by Raven, set every time step\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_evtCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t i=0;i<_evtCol.size();i++){
      int ic=_evtCol[i];
      F<<"  "<<CellID(0,ic)<<" "<<_land[ic]<<" 0.0 "<<_aProfiles[_colProfile[ic]].extinction_depth<<endl;
    }
    F<<"END PERIOD\n"; F.close();
  }
  if (_welNode.size()>0){
    F.open((d+"gwf.wel").c_str());
    F<<"# wells from .rvg :Wells; rates set by Raven every step. AUTO_FLOW_REDUCE cuts pumping as a cell dries\n";
    F<<"BEGIN OPTIONS\n  AUTO_FLOW_REDUCE 0.1\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_welNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
    for (size_t e=0;e<_welNode.size();e++){int n=_welNode[e],l=n/_ncpl,ic=n%_ncpl; F<<"  "<<CellID(l,ic)<<" 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_ghbNode.size()>0){
    F.open((d+"gwf.ghb").c_str());
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_ghbNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_ghbNode.size();e++){int n=_ghbNode[e],l=n/_ncpl,ic=n%_ncpl; const gw_boundary &B=_aBounds[_ghbB[e]];
      F<<"  "<<CellID(l,ic)<<" "<<B.head<<" "<<B.cond<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_chdNode.size()>0){
    F.open((d+"gwf.chd").c_str());
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_chdNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_chdNode.size();e++){int n=_chdNode[e],l=n/_ncpl,ic=n%_ncpl; F<<"  "<<CellID(l,ic)<<" "<<_aBounds[_chdB[e]].head<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_swxCol.size()>0){
    F.open((d+"gwf.swx").c_str());
    F<<"# surface-water stores (wetlands, lakes) exchanging with the aquifer: river-type boundary with its bed at the\n# land surface; stage = land surface + store depth and conductance are set every step\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_swxCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_swxCol.size();e++){int ic=_swxCol[e]; F<<"  "<<CellID(0,ic)<<" "<<_top[ic]<<" 0.0 "<<_top[ic]<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_resNode.size()>0){
    F.open((d+"gwf.res").c_str());
    F<<"# seepage of Raven reservoirs into the aquifer under their lake HRUs (negative: groundwater feeding the reservoir)\n";
    F<<"BEGIN OPTIONS\n  AUTO_FLOW_REDUCE 0.1\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_resNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
    for (size_t e=0;e<_resNode.size();e++){int n=_resNode[e],ic=n%_ncpl; F<<"  "<<CellID(0,ic)<<" 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  F.open((d+"gwf.oc").c_str());
  F<<"BEGIN OPTIONS\n  HEAD FILEOUT gwf.hds\nEND OPTIONS\n\nBEGIN PERIOD 1\n  SAVE HEAD FREQUENCY "<<max(1,_headSaveFreq)<<"\n  PRINT BUDGET LAST\nEND PERIOD\n"; F.close();
}

//////////////////////////////////////////////////////////////////
/// \brief MODFLOW 6 library: :MF6Library, else the RAVEN_MF6_LIB environment variable, else the first library found next
///  to the Raven executable or in the lib/mf6 folder of the source checkout (where tools/get_mf6.py puts it), else the
///  platform's default name (found through the system's library search path)
//
string CGroundwaterModel::MF6LibraryPath() const
{
  string lib=_libPath;
  if (lib==""){const char *env=getenv("RAVEN_MF6_LIB"); if ((env!=NULL) && (env[0]!='\0')){lib=env;}}
  if (lib==""){
    vector<string> c=CMF6Engine::LibraryCandidates();
    for (size_t i=0;i<c.size();i++){if (CMF6Engine::FileExists(c[i])){lib=c[i]; break;}}
  }
  if (lib==""){lib=CMF6Engine::DefaultLibraryName();}
  return lib;
}

//////////////////////////////////////////////////////////////////
/// \brief message when the library cannot be loaded: where Raven looked and how to get it
//
string CGroundwaterModel::MF6LibraryHelp() const
{
  string s=" MODFLOW 6 library not found or not loadable. Raven looks at :MF6Library in the .rvg file, then the RAVEN_MF6_LIB "
           "environment variable, then";
  vector<string> c=CMF6Engine::LibraryCandidates();
  for (size_t i=0;i<c.size();i++){s+=" "+c[i]+((i+1<c.size())?",":";");}
  s+=" then "+CMF6Engine::DefaultLibraryName()+" on the system library path. To download it, run  python tools/get_mf6.py  "
     "in the Raven source folder (or build with CMake, which fetches it), or get it from "
     "https://github.com/MODFLOW-ORG/modflow6/releases";
  return s;
}

//////////////////////////////////////////////////////////////////
/// \brief loads libmf6, initializes the simulation and binds memory pointers
//
void CGroundwaterModel::ConnectEngine()
{
  string lib=MF6LibraryPath();
  if (!_mf6.Load(lib))        {ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()+"."+MF6LibraryHelp()).c_str(),RUNTIME_ERR);}
  _mf6Version=_mf6.GetVersion();
  if (!_mf6.Initialize(_mf6dir)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()+" (see "+_mf6dir+"/mfsim.lst)").c_str(),RUNTIME_ERR);}
  _X=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"","X"));
  ExitGracefullyIf(_X==NULL,"CGroundwaterModel: cannot access MODFLOW 6 heads",RUNTIME_ERR);
  int *mxit=_mf6.GetIntPtr("SLN_1/MXITER");
  if (mxit!=NULL){_maxiter=*mxit;}
  if (IsExternal()){ExtConnect();}
  string a1=_mf6.GetVarAddress(_mname,_stoName,"STRGSS"),a2=_mf6.GetVarAddress(_mname,_stoName,"STRGSY");
  if (_mf6.HasVar(a1) && _mf6.HasVar(a2)){_strgss=_mf6.GetDoublePtr(a1);_strgsy=_mf6.GetDoublePtr(a2);}
  else {WriteWarning("CGroundwaterModel: MODFLOW 6 storage flows not accessible; GWBudget.csv storage column will be inferred",false);}
}

//bound arrays can be re-associated when stress-period data are read, so pointers are refreshed every step
void CGroundwaterModel::RefreshPointers()
{
  //MODFLOW 6.5+ packages hold user input in named arrays (RECHARGE) that are copied into BOUND each step;
  //older versions only have BOUND
  string rch=_mf6.GetVarAddress(_mname,"RCH_RAVEN","RECHARGE");
  if (!_mf6.HasVar(rch)){rch=_mf6.GetVarAddress(_mname,"RCH_RAVEN","BOUND");}
  _rchBound=_mf6.GetDoublePtr(rch);
  _rchHcof =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RCH_RAVEN","HCOF"));
  _rchRhs  =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RCH_RAVEN","RHS"));
  bool drn=!IsExternal() || _extSeepage; //an existing model may run without Raven's seepage package
  _drnHcof =drn?_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"DRN_RAVEN","HCOF")):NULL;
  _drnRhs  =drn?_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"DRN_RAVEN","RHS")):NULL;
  if (_rivNode.size()>0){
    string st=_mf6.GetVarAddress(_mname,"RIV_RAVEN","STAGE");
    if (_mf6.HasVar(st)){
      _rivNcol=0;
      _rivStage=_mf6.GetDoublePtr(st);
      _rivCond =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RIV_RAVEN","COND"));
      _rivRbot =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RIV_RAVEN","RBOT"));
    }
    else { //older MODFLOW 6: BOUND(stage,cond,rbot)
      _rivNcol=3;
      _rivStage=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RIV_RAVEN","BOUND"));
      _rivCond=_rivStage+1; _rivRbot=_rivStage+2;
    }
    _rivHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RIV_RAVEN","HCOF"));
    _rivRhs =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RIV_RAVEN","RHS"));
    ExitGracefullyIf((_rivStage==NULL) || (_rivCond==NULL) || (_rivRbot==NULL) || (_rivHcof==NULL) || (_rivRhs==NULL),
                     ("CGroundwaterModel: cannot access MODFLOW 6 river arrays: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);
  }
  if (_evtCol.size()>0){
    string ra=_mf6.GetVarAddress(_mname,"EVT_RAVEN","RATE");
    if (_mf6.HasVar(ra)){_evtNcol=0; _evtRate=_mf6.GetDoublePtr(ra);}
    else {_evtNcol=3; double *b=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"EVT_RAVEN","BOUND")); _evtRate=(b==NULL)?NULL:b+1;} //BOUND(surface,rate,depth)
    _evtHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"EVT_RAVEN","HCOF"));
    _evtRhs =_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"EVT_RAVEN","RHS"));
    ExitGracefullyIf((_evtRate==NULL) || (_evtHcof==NULL) || (_evtRhs==NULL),
                     ("CGroundwaterModel: cannot access MODFLOW 6 EVT arrays: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);
  }
  if (_welNode.size()>0){
    string a=_mf6.GetVarAddress(_mname,"WEL_RAVEN","Q");
    if (_mf6.HasVar(a)){_welNcol=0; _welQ=_mf6.GetDoublePtr(a);} else {_welNcol=1; _welQ=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"WEL_RAVEN","BOUND"));}
    _welHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"WEL_RAVEN","HCOF")); _welRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"WEL_RAVEN","RHS"));
    ExitGracefullyIf((_welQ==NULL) || (_welHcof==NULL) || (_welRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 WEL arrays",RUNTIME_ERR);
  }
  if (_ghbNode.size()>0){
    string a=_mf6.GetVarAddress(_mname,"GHB_RAVEN","BHEAD");
    if (_mf6.HasVar(a)){_ghbNcol=0; _ghbHead=_mf6.GetDoublePtr(a); _ghbCond=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"GHB_RAVEN","COND"));}
    else {_ghbNcol=2; _ghbHead=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"GHB_RAVEN","BOUND")); _ghbCond=(_ghbHead==NULL)?NULL:_ghbHead+1;}
    _ghbHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"GHB_RAVEN","HCOF")); _ghbRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"GHB_RAVEN","RHS"));
    ExitGracefullyIf((_ghbHead==NULL) || (_ghbCond==NULL) || (_ghbHcof==NULL) || (_ghbRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 GHB arrays",RUNTIME_ERR);
  }
  if (_chdNode.size()>0){
    string a=_mf6.GetVarAddress(_mname,"CHD_RAVEN","HEAD");
    if (_mf6.HasVar(a)){_chdNcol=0; _chdHead=_mf6.GetDoublePtr(a);} else {_chdNcol=1; _chdHead=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"CHD_RAVEN","BOUND"));}
    _chdSim=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"CHD_RAVEN","SIMVALS"));
    ExitGracefullyIf((_chdHead==NULL) || (_chdSim==NULL),"CGroundwaterModel: cannot access MODFLOW 6 CHD arrays",RUNTIME_ERR);
  }
  if (_swxCol.size()>0){
    string a=_mf6.GetVarAddress(_mname,"SWX_RAVEN","STAGE");
    if (_mf6.HasVar(a)){_swxNcol=0; _swxHead=_mf6.GetDoublePtr(a); _swxCond=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","COND"));
                        _swxRbot=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","RBOT"));}
    else {_swxNcol=3; _swxHead=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","BOUND")); _swxCond=(_swxHead==NULL)?NULL:_swxHead+1; _swxRbot=(_swxHead==NULL)?NULL:_swxHead+2;}
    _swxHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","HCOF")); _swxRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","RHS"));
    ExitGracefullyIf((_swxHead==NULL) || (_swxCond==NULL) || (_swxRbot==NULL) || (_swxHcof==NULL) || (_swxRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 surface-water exchange arrays",RUNTIME_ERR);
  }
  if (_resNode.size()>0){
    string a=_mf6.GetVarAddress(_mname,"RES_RAVEN","Q");
    if (_mf6.HasVar(a)){_resNcol=0; _resQ=_mf6.GetDoublePtr(a);} else {_resNcol=1; _resQ=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RES_RAVEN","BOUND"));}
    _resHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RES_RAVEN","HCOF")); _resRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"RES_RAVEN","RHS"));
    ExitGracefullyIf((_resQ==NULL) || (_resHcof==NULL) || (_resRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 reservoir-exchange arrays",RUNTIME_ERR);
  }
  ExitGracefullyIf((_rchBound==NULL) || (_rchHcof==NULL) || (_rchRhs==NULL) || (drn && ((_drnHcof==NULL) || (_drnRhs==NULL))),
                   ("CGroundwaterModel: cannot access MODFLOW 6 boundary arrays: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);
}

//=================================================================
// Per-step exchange
//=================================================================
//////////////////////////////////////////////////////////////////
/// \brief one coupled step: Raven recharge -> MODFLOW; solve; seepage -> Raven stores
/// \param aPhinew [in/out] Raven state variables at end of step (GROUNDWATER holds this step's recharge [mm])
//
void CGroundwaterModel::Exchange(double **aPhinew,const double &tstep,const optStruct &Options,const time_struct &tt)
{
  FinishBudgetLine(); //river terms of the previous step are final once routing is done
  _stepCount++;
  ExitGracefullyIf(_stepCount>_nstp,"CGroundwaterModel: Raven is running more time steps than the MODFLOW 6 simulation holds (check :Duration and :TimeStep)",RUNTIME_ERR);
  if (IsExternal()){ //an existing model's specified-head cells are known only once MODFLOW has read the step's period data
    if (!_mf6.PrepareTimeStep(tstep)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
    RefreshPointers();
    ExtUpdateCHDTop();
  }
  //-- 1. Raven recharge -> column flows [m3/d], mass-conservative per HRU
  std::fill(_Qcol.begin(),_Qcol.end(),0.0);
  _ravenRchVol=0.0; _QrchCHD=0.0;
  if (_work.size()!=(size_t)_ncpl){_work.assign(_ncpl,0.0);}
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i];
    double V=aPhinew[k][_iGW]/MM_PER_METER*_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2; //[m3]
    _ravenRchVol+=V;
    double tau=_aProfiles[_hruProfile[k]].recharge_delay;
    if (tau>0){ //linear-reservoir lag through a deep vadose zone
      _rchStore[k]+=V;
      V=_rchStore[k]*(1.0-exp(-tstep/tau));
      _rchStore[k]-=V;
    }
    if (_hruRchArea[k]<=0){_QrchCHD+=V/tstep; continue;} //HRU entirely on specified-head cells
    for (size_t j=0;j<_hruLinks[k].size();j++){
      const gw_link &L=_aLinks[_hruLinks[k][j]];
      if (_colCHDTop[L.ic]){continue;}
      _Qcol[L.ic]+=V/tstep*L.area/_hruRchArea[k];
    }
  }

  _rchStoreTot=0.0; for (size_t i=0;i<_aqHRUs.size();i++){_rchStoreTot+=_rchStore[_aqHRUs[i]];}
  //-- 2. MODFLOW 6 time step
  if (!IsExternal()){
    if (!_mf6.PrepareTimeStep(tstep)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
    RefreshPointers();
  }
  for (size_t i=0;i<_bndCol.size();i++){_rchBound[i]=_Qcol[_bndCol[i]]/MF6CellArea(_bndCol[i])*_rf;} //rate in model units, over MODFLOW's own cell area
  //water-table ET rate = PET not already used by Raven this step, area-weighted over the column's HRUs [m/d]
  if (_evtCol.size()>0){
    int iAET=_pModel->GetStateVarIndex(AET);
    std::fill(_work.begin(),_work.end(),0.0); //sum(area*unmet PET)
    for (size_t li=0;li<_aLinks.size();li++){
      const gw_link &L=_aLinks[li];
      if (_aProfiles[_hruProfile[L.k]].extinction_depth<=0){continue;} //only HRUs whose profile has water-table ET
      double pet  =_pModel->GetHydroUnit(L.k)->GetForcingFunctions()->PET*tstep;      //[mm]
      double used =(iAET!=DOESNT_EXIST)?aPhinew[L.k][iAET]:0.0;                         //[mm]
      _work[L.ic]+=L.area*max(pet-used,0.0)/tstep/MM_PER_METER;                        //[m3/d]
    }
    int es=max(_evtNcol,1);
    //MODFLOW multiplies the rate by the whole cell area: spread the unmet-PET volume of the coupled HRUs over the cell,
    //so the maximum ET MODFLOW can remove equals that volume (not more where HRUs cover only part of the cell)
    for (size_t i=0;i<_evtCol.size();i++){int ic=_evtCol[i]; _evtRate[i*es]=_work[ic]/_grid.area[ic];}
  }
  SetBoundaryArrays(tt.model_time,tstep);
  //reservoirs: seepage Raven computed during the previous step's routing enters the aquifer now (one-step lag)
  _Qres_spec=0.0;
  for (size_t e=0;e<_resNode.size();e++){
    CReservoir *pRes=_pModel->GetSubBasin(_resP[_resR[e]])->GetReservoir();
    double V=pRes->GetReservoirGWLosses(tstep);
    if (_stepCount==1){ //first step: nothing is in transit, unless a hotstart file says so (the reservoir object may
                        //still hold the last seepage of a previous ensemble member)
      V=0.0;
      map<long long,double>::iterator it=_hotResSeep.find(_pModel->GetSubBasin(_resP[_resR[e]])->GetID());
      if (it!=_hotResSeep.end()){V=it->second;}
    }
    double q=V/tstep*_resW[e]; //[m3/d] (+: into aquifer)
    _resQ[e*max(_resNcol,1)]=q; _Qres_spec+=q;
  }
  //surface-water stores: head = land surface + store depth; conductance limited so the largest possible
  //leakage in the step cannot exceed the water the stores hold in that cell
  for (size_t e=0;e<_swxCol.size();e++){
    const gw_swx &S=_aSWX[_swxX[e]]; int ic=_swxCol[e];
    double A=0,dA=0,V=0;
    if (_swxVl.size()!=_swxCol.size()){_swxVl.assign(_swxCol.size(),vector<double>());}
    _swxVl[e].assign(_swxLinks[e].size(),0.0);
    for (size_t j=0;j<_swxLinks[e].size();j++){
      const gw_link &L=_aLinks[_swxLinks[e][j]]; double stor=max(aPhinew[L.k][S.iSV],0.0)/MM_PER_METER; //[m]
      A+=L.area; dA+=L.area*stor;
      //water available in this cell = the HRU's real store volume (depth x .rvh area) times the cell's share
      //of the HRU's linked area, so the shares of all cells add up to exactly what Raven holds
      _swxVl[e][j]=stor*_pModel->GetHydroUnit(L.k)->GetArea()*M2_PER_KM2*L.area/_hruLinkedArea[L.k];
      V+=_swxVl[e][j];
    }
    double sbar=(A>0)?dA/A:0.0;                       //mean store depth [m]
    double bhead=_top[ic]+sbar;                         //stage; bed (RBOT) = model top
    double C=S.leakance*A, dhmax=sbar;                  //leakage is at most C*sbar (water table below the bed)
    //supply limit: the largest possible leakage, C*(bhead - cell bottom)*dt, may not exceed the water the stores
    //hold, so a store can never be over-drawn whatever the head does during the solve. (Groundwater reaching the
    //surface under an empty store is still returned to Raven by the seepage drains at the model top.)
    if (dhmax>1e-6){C=min(C,max(V,0.0)/(tstep*dhmax));}
    if ((V<=0.0) || (sbar<=0.0)){C=0.0;}
    int st=max(_swxNcol,1); _swxHead[e*st]=bhead; _swxCond[e*st]=C; _swxRbot[e*st]=_top[ic]; //bed set in memory: exactly equal to stage when the store is empty
  }
  bool conv; int niter;
  //target: _dvclose; a step not converged after 100 iterations continues with MODFLOW's tolerance at 10*_dvclose
  double *dvc=_mf6.GetDoublePtr("SLN_1/DVCLOSE"); //MODFLOW's own tolerance: raised for the rest of a hard step only
  //hard steps: after 100 iterations (or half the iteration limit, if that is smaller) MODFLOW's own tolerance is relaxed
  //tenfold for the rest of the step (in the model's
  //units: with an existing model its own DVCLOSE is the target)
  double relaxTol=(IsExternal() && (_extDvclose>0))?10.0*_extDvclose/_hf:10.0*_dvclose;
  if (!_mf6.SolveTimeStep(_maxiter,conv,niter,dvc,relaxTol,min(100,max(1,_maxiter/2)),&_lastRelaxed)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  _lastConverged=conv; _lastIter=niter;
  if (!conv){
    _nNonConverged++;
    if (_nNonConverged<=10){WriteWarning("CGroundwaterModel: MODFLOW 6 did not converge on "+tt.date_string+" (continuing)",Options.noisy);}
  }
  if (!_mf6.FinalizeTimeStep()){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}

  //-- 3. boundary flows from final heads (q = hcof*h - rhs, positive into aquifer)
  _Qrch=0.0; _Qseep_tot=0.0;
  std::fill(_Qseep.begin(),_Qseep.end(),0.0);
  for (size_t i=0;i<_bndCol.size();i++){
    int ic=_bndCol[i];
    int node=_reduced[TopNode(ic)];
    double h=_X[node];
    _Qrch+=(_rchHcof[i]*h-_rchRhs[i])*_qf; //[m3/d]
    double qd=(_drnHcof==NULL)?0.0:-(_drnHcof[i]*h-_drnRhs[i])*_qf; //discharge out of aquifer [m3/d]
    _Qseep[ic]=qd; _Qseep_tot+=qd;
  }
  _Qriv_in=_Qriv_out=0.0;
  std::fill(_sbGain.begin(),_sbGain.end(),0.0);
  for (size_t e=0;e<_rivNode.size();e++){
    double q=(_rivHcof[e]*_X[_reduced[_rivNode[e]]]-_rivRhs[e])*_qf; //into aquifer [m3/d]
    if (q>0){_Qriv_in+=q;} else {_Qriv_out-=q;}
    _sbGain[_rivSB[e]]-=q;
  }
  if (IsExternal()){ExtStreamFlows();}
  ComputeStressFlows();
  _Qres=0.0; _resShortVol=0.0;
  for (size_t e=0;e<_resNode.size();e++){
    double act=(_resHcof[e]*_X[_reduced[_resNode[e]]]-_resRhs[e])*_qf,spec=_resQ[e*max(_resNcol,1)]; //[m3/d], + into aquifer
    _Qres+=act;
    //a gain Raven already gave the lake (spec<0) that MODFLOW could not fully supply because the cell was drying
    //(AUTO_FLOW_REDUCE): the difference is taken back from the water reaching that subbasin in this step's routing,
    //like a river loss, so no water is created (reported as "reservoir gain not supplied by MODFLOW")
    if ((spec<0.0) && (act>spec) && (_resP[_resR[e]]<(int)_sbLossCarry.size())){
      double S=(act-spec)*tstep; _sbLossCarry[_resP[_resR[e]]]+=S; _resShortVol+=S;
    }
  }
  for (size_t r=0;r<_resP.size();r++){ //local groundwater head seen by each reservoir during this step's routing
    //a lake above a water table deeper than its bed is disconnected: the seepage then depends on the bed, not on the
    //water-table depth (as for MODFLOW's river boundary). The bed is the model top, unless the lake stands below it
    //(a DEM that recorded a higher water surface): then the bed is taken at the lake level, so a deep water table
    //gives no exchange rather than a false gain
    CReservoir *pRes=_pModel->GetSubBasin(_resP[r])->GetReservoir();
    double stage=pRes->GetResStage();
    double h=0; for (size_t e=0;e<_resNode.size();e++){if (_resR[e]==(int)r){int ic=_resNode[e]%_ncpl; h+=_resW[e]*max(WaterTable(ic),min(_top[ic],stage));}}
    pRes->SetGWParameters(pRes->GetSeepageConstant(),h);
  }
  _Qswx_in=_Qswx_out=0.0; _swxAppliedVol=0.0;
  _swxQ.assign(_swxCol.size(),0.0);
  for (size_t e=0;e<_swxCol.size();e++){
    double q=(_swxHcof[e]*_X[_reduced[TopNode(_swxCol[e])]]-_swxRhs[e])*_qf; //[m3/d] into aquifer (MODFLOW's value)
    _swxQ[e]=q;
    if (q>0){_Qswx_in+=q;} else {_Qswx_out-=q;}
  }
  _Qevt=0.0;
  for (size_t i=0;i<_evtCol.size();i++){
    int ic=_evtCol[i];
    _Qevt-=(_evtHcof[i]*_X[_reduced[TopNode(ic)]]-_evtRhs[i])*_qf; //EVT flow is out of the aquifer
  }
  _Qsto=0.0;
  if (_strgss!=NULL){
    int nn=0; for (int n=0;n<_nlay*_ncpl;n++){if (_idomain[n]>0){nn++;}}
    for (int n=0;n<nn;n++){_Qsto+=(_strgss[n]+_strgsy[n])*_qf;}
    _Qerr=_Qrch+_Qsto-_Qseep_tot+_Qriv_in-_Qriv_out-_Qevt+_Qwel+_Qghb+_Qchd+_Qswx_in-_Qswx_out+_Qres; //(_QrchCHD never enters MODFLOW)
    if (IsExternal()){_Qerr+=_Qother;}
  }
  else if (IsExternal() && _extNoSto){ //an existing model without a storage package is steady state: no storage flow, the residual is MODFLOW's error
    _Qsto=0.0; _Qerr=_Qrch-_Qseep_tot+_Qriv_in-_Qriv_out-_Qevt+_Qwel+_Qghb+_Qchd+_Qswx_in-_Qswx_out+_Qres+_Qother;
  }
  else {_Qsto=_Qseep_tot-_Qrch-_Qriv_in+_Qriv_out+_Qevt-_Qwel-_Qghb-_Qchd-_Qswx_in+_Qswx_out-_Qres-_Qother; _Qerr=0.0;}

  //-- guard: an unconverged solution can produce spurious flows. Seepage returned to Raven is then limited to
  //   the water that entered the aquifer this step (recharge + river leakage), scaled down uniformly
  if (!_lastConverged){ //all discharge returned to Raven (seepage + river gains) limited to this step's inflow
    double cap=max(_Qrch+_Qriv_in+_Qswx_in+max(_Qres,0.0)+max(_Qwel,0.0)+max(_Qghb,0.0)+max(_Qchd,0.0)
                   +(IsExternal()?max(_Qother,0.0):0.0),0.0); //all inflows this step (with an existing model, its own packages too)
    double gains=0.0; for (size_t p=0;p<_sbGain.size();p++){gains+=max(_sbGain[p],0.0);}
    for (size_t e=0;e<_swxQ.size();e++){gains+=max(-_swxQ[e],0.0);} //discharge into surface stores
    if (_Qseep_tot+gains>cap){
      double f=(_Qseep_tot+gains>0)?cap/(_Qseep_tot+gains):0.0;
      for (int ic=0;ic<_ncpl;ic++){_Qseep[ic]*=f;}
      for (size_t p=0;p<_sbGain.size();p++){if (_sbGain[p]>0){_sbGain[p]*=f;}}
      for (size_t e=0;e<_swxQ.size();e++){if (_swxQ[e]<0){_swxQ[e]*=f;}}
      WriteWarning("CGroundwaterModel: unconverged step ending "+EndDate(tt.model_time+tstep)+": discharge returned to Raven limited from "+
                   ToStr(_Qseep_tot+gains)+" to "+ToStr(cap)+" m3/d",false);
      _Qseep_tot*=f;
    }
  }
  //-- surface-store exchange applied to Raven stores (after the guard)
  for (size_t e=0;e<_swxCol.size();e++){
    const gw_swx &S=_aSWX[_swxX[e]];
    double vol=_swxQ[e]*tstep, V=0, A=0;
    const vector<double> &Vl=_swxVl[e]; //available water frozen before the solve (the cap was set from it)
    for (size_t j=0;j<_swxLinks[e].size();j++){V+=Vl[j]; A+=_aLinks[_swxLinks[e][j]].area;}
    if ((vol>0) && (vol>V)){ //bounded by the supply-limited conductance; only round-off can exceed it
      if (vol-V>1e-6*V+1e-3){WriteWarning("CGroundwaterModel: surface-store leakage exceeded available water by "+ToStr(vol-V)+" m3 (clipped)",false);}
      vol=V;
    }
    for (size_t j=0;j<_swxLinks[e].size();j++){
      const gw_link &L=_aLinks[_swxLinks[e][j]];
      double share=(vol>0)?((V>0)?Vl[j]/V:0.0):L.area/A; //losses by available water, gains by area
      double dV=-vol*share; //[m3] change of the Raven store
      aPhinew[L.k][S.iSV]+=dV/(_pModel->GetHydroUnit(L.k)->GetArea()*M2_PER_KM2)*MM_PER_METER;
      _swxAppliedVol+=dV;
    }
  }
  //-- 4. seepage -> Raven store of overlapping HRUs (share of column by overlap area)
  _returnVol=0.0;
  for (size_t li=0;li<_aLinks.size();li++){
    const gw_link &L=_aLinks[li];
    if (_Qseep[L.ic]==0.0){continue;}
    double V=_Qseep[L.ic]*tstep*L.area/_colCover[L.ic]; //[m3]
    aPhinew[L.k][_hruReturnSV[L.k]]+=V/(_pModel->GetHydroUnit(L.k)->GetArea()*M2_PER_KM2)*MM_PER_METER;
    _returnVol+=V;
  }
  WriteStepOutputs(tstep,tt);
}

//////////////////////////////////////////////////////////////////
/// \brief sets river stage/conductance, well rates and boundary heads for the step starting at t
//
void CGroundwaterModel::SetBoundaryArrays(const double &t,const double &tstep)
{
  //river stage from Raven reach depth, conductance from current top width (start of step)
  int stride=max(_rivNcol,1);
  for (size_t e=0;e<_rivNode.size();e++){
    CSubBasin *pSB=_pModel->GetSubBasin(_rivSB[e]);
    const gw_profile &P=_aProfiles[_colProfile[_rivCol[e]]];
    double depth=max(pSB->GetRiverDepth(),0.0), width=max(pSB->GetTopWidth(),0.0);
    if (!(depth<1e6)){depth=0.0;} if (!(width<1e6)){width=0.0;}
    _rivStage[e*stride]=min(_rivBed[e]+depth,max(_land[_rivCol[e]],_rivBed[e])); //above bankfull the water is on the floodplain, not in the channel
    _rivRbot [e*stride]=_rivBed[e];
    double C=P.riverbed_K*width*_rivLen[e]/max(P.riverbed_thick,0.01);
    //supply limit: the largest possible loss, C*(stage-rbot), may not exceed this cell's share (by river length)
    //of the flow carried by the reach at the start of the step, so MODFLOW never takes water the river does not have
    //the river cells lie on the channel reach, upstream of any reservoir at the subbasin outlet: use the channel flow
    //(GetOutflowRate would return the reservoir release, which can be zero while the channel flows)
    double qavail=max(pSB->GetChannelOutflowRate(),0.0)*SEC_PER_DAY*_rivLen[e]/max(_sbRivLen[_rivSB[e]],1e-9); //[m3/d]
    double dh=_rivStage[e*stride]-_rivRbot[e*stride];
    if (dh>1e-6){C=min(C,qavail/dh);}
    _rivCond [e*stride]=C;
  }
  _Qwel_spec=0.0;
  for (size_t e=0;e<_welNode.size();e++){
    const gw_well &W=_aWells[_welWell[e]];
    double rate=(W.pRate==NULL)?0.0:W.pRate->GetAvgValue(t,tstep);
    if ((rate==RAV_BLANK_DATA) || !(fabs(rate)<1e20)){rate=0.0;}
    _welQ[e*max(_welNcol,1)]=rate*_welFrac[e];
    _Qwel_spec+=rate*_welFrac[e];
  }
  for (size_t e=0;e<_ghbNode.size();e++){
    const gw_boundary &B=_aBounds[_ghbB[e]];
    double h=B.head; if (B.pHead!=NULL){double v=B.pHead->GetAvgValue(t,tstep); if (v!=RAV_BLANK_DATA){h=v;}}
    _ghbHead[e*max(_ghbNcol,1)]=max(h,_botm[_ghbNode[e]]+0.01); _ghbCond[e*max(_ghbNcol,1)]=B.cond;
  }
  for (size_t e=0;e<_chdNode.size();e++){
    const gw_boundary &B=_aBounds[_chdB[e]];
    double h=B.head; if (B.pHead!=NULL){double v=B.pHead->GetAvgValue(t,tstep); if (v!=RAV_BLANK_DATA){h=v;}}
    _chdHead[e*max(_chdNcol,1)]=max(h,_botm[_chdNode[e]]+0.01);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief flows of wells and head boundaries after a solve (positive into the aquifer) [m3/d]
//
void CGroundwaterModel::ComputeStressFlows()
{
  _Qwel=_Qghb=_Qchd=0.0;
  for (size_t e=0;e<_welNode.size();e++){_Qwel+=(_welHcof[e]*_X[_reduced[_welNode[e]]]-_welRhs[e])*_qf;}
  for (size_t e=0;e<_ghbNode.size();e++){_Qghb+=(_ghbHcof[e]*_X[_reduced[_ghbNode[e]]]-_ghbRhs[e])*_qf;}
  for (size_t e=0;e<_chdNode.size();e++){_Qchd+=_chdSim[e];}
}
//////////////////////////////////////////////////////////////////
/// \brief steady-state first stress period with uniform recharge, solved before Raven's first step
//
void CGroundwaterModel::RunSteadyState()
{
  if (!_mf6.PrepareTimeStep(1.0)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  RefreshPointers();
  for (size_t i=0;i<_bndCol.size();i++){_rchBound[i]=_ssRecharge/MM_PER_METER*_colCover[_bndCol[i]]/_grid.area[_bndCol[i]];}
  SetBoundaryArrays(0.0,1.0);
  //keep the starting heads: a steady state that does not converge is a worse start than they are
  int nx=0; for (size_t i=0;i<_reduced.size();i++){if (_reduced[i]>=0){nx++;}}
  vector<double> x0; if (_X!=NULL){x0.assign(_X,_X+nx);}
  bool conv; int niter;
  if (!_mf6.SolveTimeStep(_maxiter,conv,niter)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  if (!_mf6.FinalizeTimeStep()){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  cout<<"  Groundwater steady-state initialization "<<(conv?"converged":"DID NOT CONVERGE")<<" in "<<niter<<" outer iterations"<<endl;
  if (!conv){
    if (!x0.empty() && (_X!=NULL)){for (int i=0;i<nx;i++){_X[i]=x0[i];}} //MODFLOW takes these as the previous heads of the first step
    WriteWarning("CGroundwaterModel: steady-state initialization did not converge (a steady state may not exist, e.g. wells pumping more "
                 "than the aquifer can supply); the run starts from the initial heads (:GWInitialization) instead",false);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief NetCDF head field output (CF conventions): head(time,layer,y,x) with 2-D lat/lon; requires -Dnetcdf
//
void CGroundwaterModel::OpenNetCDF(const optStruct &Options)
{
#ifdef netcdf
  string fn=FilenamePrepare("GWHeads.nc",Options);
  int dT,dL,dY,dX,vX,vY,vLat,vLon,vL;
  if (nc_create(fn.c_str(),NC_CLOBBER|NC_NETCDF4,&_ncid)!=NC_NOERR){WriteWarning("CGroundwaterModel: cannot create "+fn,false); _writeNC=false; return;}
  bool st=(_gridType==GWGRID_REGULAR); int ny=st?_nrow:1,nx=st?_ncol:_ncpl;
  nc_def_dim(_ncid,"time",NC_UNLIMITED,&dT); nc_def_dim(_ncid,"layer",_nlay,&dL);
  nc_def_dim(_ncid,st?"y":"one",ny,&dY); nc_def_dim(_ncid,st?"x":"cell",nx,&dX);
  nc_def_var(_ncid,"time",NC_DOUBLE,1,&dT,&_ncTime);
  time_struct t0; JulianConvert(0.0,_jday0,_jyear0,_calendar,t0);
  string units="days since "+t0.date_string+" 00:00:00";
  nc_put_att_text(_ncid,_ncTime,"units",units.size(),units.c_str());
  nc_def_var(_ncid,"layer",NC_INT,1,&dL,&vL);
  if (st){nc_def_var(_ncid,"x",NC_DOUBLE,1,&dX,&vX); nc_def_var(_ncid,"y",NC_DOUBLE,1,&dY,&vY);}
  else   {nc_def_var(_ncid,"x",NC_DOUBLE,1,&dX,&vX); nc_def_var(_ncid,"y",NC_DOUBLE,1,&dX,&vY);}
  nc_put_att_text(_ncid,vX,"units",1,"m"); nc_put_att_text(_ncid,vY,"units",1,"m");
  int d2[2]={dY,dX}; nc_def_var(_ncid,"lat",NC_DOUBLE,2,d2,&vLat); nc_def_var(_ncid,"lon",NC_DOUBLE,2,d2,&vLon);
  nc_put_att_text(_ncid,vLat,"units",13,"degrees_north"); nc_put_att_text(_ncid,vLon,"units",12,"degrees_east");
  int d4[4]={dT,dL,dY,dX}; nc_def_var(_ncid,"head",NC_FLOAT,4,d4,&_ncHead);
  float fill=-9999.0f; nc_put_att_float(_ncid,_ncHead,"_FillValue",NC_FLOAT,1,&fill);
  nc_put_att_text(_ncid,_ncHead,"units",1,"m"); nc_put_att_text(_ncid,_ncHead,"coordinates",7,"lat lon");
  nc_put_att_text(_ncid,_ncHead,"long_name",20,"groundwater head (m)");
  if (!st){string g="unstructured grid: cell polygons in mf6/grid_cells.geojson (cell = index+1)"; nc_put_att_text(_ncid,NC_GLOBAL,"grid",g.size(),g.c_str());}
  ostringstream pr; pr<<"Lambert azimuthal equal-area, sphere R=6371007.181 m, lat0="<<setprecision(10)<<_proj.GetLat0()<<" lon0="<<_proj.GetLon0();
  if (_grid.angle!=0.0){pr<<"; grid x/y rotated "<<_grid.angle<<" deg counterclockwise";}
  string ps=pr.str(); nc_put_att_text(_ncid,NC_GLOBAL,"projection",ps.size(),ps.c_str());
  nc_put_att_text(_ncid,NC_GLOBAL,"source",34,"Raven-MODFLOW 6 groundwater coupling");
  nc_enddef(_ncid);
  vector<int> L(_nlay); for (int l=0;l<_nlay;l++){L[l]=l+1;} nc_put_var_int(_ncid,vL,&L[0]);
  vector<double> xs(nx),ys(st?ny:nx),la(_ncpl),lo(_ncpl);
  if (st){
    for (int c=0;c<_ncol;c++){xs[c]=_grid.uniform?_x0+(c+0.5)*_cellSize:_grid.colX[c]+0.5*_grid.colW[c];}
    for (int r=0;r<_nrow;r++){ys[r]=_grid.uniform?_y0+(_nrow-1-r+0.5)*_cellSize:_grid.rowY[r]+0.5*_grid.rowH[r];}
  }
  else {for (int ic=0;ic<_ncpl;ic++){xs[ic]=_grid.cx[ic]; ys[ic]=_grid.cy[ic];}}
  for (int ic=0;ic<_ncpl;ic++){double X,Y; CellCentreWorld(ic,X,Y); _proj.Inverse(X,Y,lo[ic],la[ic]);}
  nc_put_var_double(_ncid,vX,&xs[0]); nc_put_var_double(_ncid,vY,&ys[0]); nc_put_var_double(_ncid,vLat,&la[0]); nc_put_var_double(_ncid,vLon,&lo[0]);
  _ncRec=0;
#else
  (void)Options;
  WriteWarning("CGroundwaterModel: :WriteNetCDFHeads ignored - Raven was compiled without NetCDF (-Dnetcdf)",false);
  _writeNC=false;
#endif
}
void CGroundwaterModel::WriteNetCDF(const double &t)
{
#ifdef netcdf
  if (_ncid<0){return;}
  vector<float> h(_nlay*_ncpl,-9999.0f);
  for (int n=0;n<_nlay*_ncpl;n++){if (_reduced[n]>=0){h[n]=(float)(_X[_reduced[n]]*_hf);}}
  size_t s1[1]={(size_t)_ncRec},c1[1]={1}; nc_put_vara_double(_ncid,_ncTime,s1,c1,&t);
  bool st=(_gridType==GWGRID_REGULAR); size_t s4[4]={(size_t)_ncRec,0,0,0},c4[4]={1,(size_t)_nlay,(size_t)(st?_nrow:1),(size_t)(st?_ncol:_ncpl)}; nc_put_vara_float(_ncid,_ncHead,s4,c4,&h[0]);
  _ncRec++;
#else
  (void)t;
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief simulated head at a monitoring well [m]
//
double CGroundwaterModel::GetObservationWellHead(const long long id) const
{
  if ((!_active) || (_X==NULL)){return RAV_BLANK_DATA;}
  for (size_t o=0;o<_aObsWells.size();o++){
    if ((_aObsWells[o].id==id) && (_aObsWells[o].node>=0)){
      double h=_X[_reduced[_aObsWells[o].node]]*_hf;
      if (h<_botm[_aObsWells[o].node]){return RAV_BLANK_DATA;} //screened cell is dry: no water level to observe
      return h;
    }
  }
  return RAV_BLANK_DATA;
}
//////////////////////////////////////////////////////////////////
/// \brief writes cell-by-cell heads to the solution (.rvc) file for hotstart
//
void CGroundwaterModel::WriteHotstart(ofstream &RVC) const
{
  if ((!_active) || (_X==NULL)){return;}
  if ((_gridType==GWGRID_REGULAR) || (IsExternal() && (_disType==0))){RVC<<":GWHeads "<<_nlay<<" "<<_nrow<<" "<<_ncol<<endl<<setprecision(12);}
  else {RVC<<":GWHeads "<<_nlay<<" 1 "<<_ncpl<<endl<<setprecision(12);} //unstructured: one row of cells
  for (int n=0;n<_nlay*_ncpl;n++){
    double h=(_reduced[n]>=0)?_X[_reduced[n]]:_strt[n];
    RVC<<((n%10==0)?"  ":" ")<<h; if ((n%10==9) || (n==_nlay*_ncpl-1)){RVC<<endl;}
  }
  RVC<<":EndGWHeads"<<endl;
  if (_resP.size()>0){ //reservoir seepage computed in the last step, not yet delivered to the aquifer
    RVC<<":GWReservoirSeepage   # subbasin ID, seepage in transit [m3]"<<endl;
    for (size_t r=0;r<_resP.size();r++){RVC<<"  "<<_pModel->GetSubBasin(_resP[r])->GetID()<<", "<<setprecision(15)<<_pModel->GetSubBasin(_resP[r])->GetReservoir()->GetReservoirGWLosses(1.0)<<endl;}
    RVC<<":EndGWReservoirSeepage"<<endl;
  }
  bool anyCarry=false; for (size_t p=0;p<_sbLossRemain.size();p++){if (_sbLossRemain[p]>0){anyCarry=true;}}
  if (anyCarry){ //river loss not yet supplied at the end of the run
    RVC<<":GWRiverLossCarry   # subbasin ID, river loss to be taken from the river in the next step [m3]"<<endl;
    for (size_t p=0;p<_sbLossRemain.size();p++){if (_sbLossRemain[p]>0){RVC<<"  "<<_pModel->GetSubBasin((int)p)->GetID()<<", "<<setprecision(15)<<_sbLossRemain[p]<<endl;}}
    RVC<<":EndGWRiverLossCarry"<<endl;
  }
  bool anyDelay=false; for (size_t p=0;p<_aProfiles.size();p++){if (_aProfiles[p].recharge_delay>0){anyDelay=true;}}
  if (anyDelay){
    RVC<<":GWRechargeStore   # HRU ID, water in recharge-delay reservoir [m3]"<<endl;
    for (size_t i=0;i<_aqHRUs.size();i++){RVC<<"  "<<_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID()<<", "<<_rchStore[_aqHRUs[i]]<<endl;}
    RVC<<":EndGWRechargeStore"<<endl;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief adds aquifer discharge to (or removes river leakage from) each subbasin's in-catchment inflow
/// \details losses are limited to the water Raven routes to the reach this step; any shortfall is
///  reported as a deficit in GWBudget.csv rather than silently created
//
void CGroundwaterModel::ApplyRiverExchange(double *aRouted,const double &tstep)
{
  _riverAppliedVol=0.0; _riverDeficitVol=0.0;
  for (size_t p=0;p<_sbGain.size();p++){
    _sbLossRemain[p]=0.0;
    if (((_sbGain[p]==0.0) && (_sbLossCarry[p]==0.0)) || (!_pModel->GetSubBasin((int)p)->IsEnabled())){continue;}
    double old=aRouted[p];
    aRouted[p]+=_sbGain[p]*tstep-_sbLossCarry[p]; //loss left unsupplied last step is taken now
    _sbLossCarry[p]=0.0;
    if (aRouted[p]<0.0){_sbLossRemain[p]=-aRouted[p]; _riverDeficitVol+=-aRouted[p]; aRouted[p]=0.0;} //provisional; see TakeRiverLossFromInflow
    _riverAppliedVol+=aRouted[p]-old;
  }
}
//////////////////////////////////////////////////////////////////
/// \brief called in the routing loop: river loss not covered by local runoff is taken from the flow entering the reach
//
double CGroundwaterModel::TakeRiverLossFromInflow(const int p,const double &Qin,const double &tstep)
{
  if ((p>=(int)_sbLossRemain.size()) || (_sbLossRemain[p]<=0.0) || (Qin<=0.0)){return Qin;}
  double avail=Qin*tstep*SEC_PER_DAY;
  double take =min(_sbLossRemain[p],avail);
  _sbLossRemain[p]-=take; _riverDeficitVol-=take; _riverAppliedVol-=take;
  return Qin-take/(tstep*SEC_PER_DAY);
}
void CGroundwaterModel::FinishBudgetLine()
{
  //river loss still unsupplied after routing is carried to the next step (reported as "river loss carried")
  //(added, not assigned: on the first step this keeps a carry read from a hotstart file)
  for (size_t p=0;p<_sbLossRemain.size();p++){_sbLossCarry[p]+=_sbLossRemain[p]; _sbLossRemain[p]=0.0;}
  if (_budgetLineOpen && _BUDGET.is_open()){
    _BUDGET<<","<<_Qriv_in<<","<<_Qriv_out<<","<<_riverAppliedVol<<","<<_riverDeficitVol<<","<<_Qevt<<","<<_Qwel_spec<<","<<_Qwel<<","<<_Qghb<<","<<_Qchd<<","<<_Qswx_in<<","<<_Qswx_out<<","<<_swxAppliedVol;
    _resPendingVol=0.0;
    for (size_t r=0;r<_resP.size();r++){_resPendingVol+=_pModel->GetSubBasin(_resP[r])->GetReservoir()->GetReservoirGWLosses(1.0);}
    _BUDGET<<","<<_Qres_spec<<","<<_Qres<<","<<_resPendingVol<<","<<_resShortVol;
    if (IsExternal()){_BUDGET<<","<<_Qother;}
    _BUDGET<<endl;
  }
  _budgetLineOpen=false;
}

//////////////////////////////////////////////////////////////////
/// \brief diagnostic: wet (saturated) fraction of the top layer per profile and number of connected wet clusters
//
void CGroundwaterModel::WriteConnectivity(const double &t,const string &dstr)
{
  vector<int> wet(_ncpl,0),lab(_ncpl,-1);
  vector<int> nprof(_aProfiles.size(),0),nwet(_aProfiles.size(),0);
  for (size_t i=0;i<_bndCol.size();i++){
    int ic=_bndCol[i]; int n=TopNode(ic);
    wet[ic]=(_X[_reduced[n]]*_hf>_botm[n]+1e-3)?1:0;
    nprof[_colProfile[ic]]++; nwet[_colProfile[ic]]+=wet[ic];
  }
  int nclust=0; vector<int> stack;
  for (int ic=0;ic<_ncpl;ic++){
    if (!wet[ic] || (lab[ic]>=0)){continue;}
    nclust++; lab[ic]=nclust; stack.push_back(ic);
    while (!stack.empty()){
      int c=stack.back(); stack.pop_back();
      for (size_t k=0;k<_grid.nbr[c].size();k++){int nb=_grid.nbr[c][k].j; if (wet[nb] && (lab[nb]<0)){lab[nb]=nclust;stack.push_back(nb);}}
    }
  }
  int tw=0,tn=0; for (size_t p=0;p<nprof.size();p++){tw+=nwet[p];tn+=nprof[p];}
  _CONN<<setprecision(8)<<t<<","<<dstr<<","<<nclust<<","<<(double)tw/max(tn,1);
  for (size_t p=0;p<nprof.size();p++){_CONN<<","<<(double)nwet[p]/max(nprof[p],1);}
  _CONN<<endl;
}

//=================================================================
// Outputs
//=================================================================
void CGroundwaterModel::OpenOutputs(const optStruct &Options)
{
  _BUDGET.open(FilenamePrepare("GWBudget.csv",Options).c_str());
  ExitGracefullyIf(_BUDGET.fail(),"CGroundwaterModel: cannot open GWBudget.csv",FILE_OPEN_ERR);
  _BUDGET<<"time [d],date,Raven recharge [m3/d],MF6 recharge [m3/d],seepage to Raven [m3/d],storage release [m3/d],"
         <<"MF6 balance error [m3/d],MF6 balance error [%],returned to Raven [m3],outer iterations,converged,recharge delay storage [m3],recharge onto CHD cells [m3/d],"
         <<"river leakage to aquifer [m3/d],aquifer discharge to rivers [m3/d],river exchange applied to Raven [m3],river loss carried to next step [m3],water-table ET [m3/d],wells specified [m3/d],wells actual [m3/d],GHB [m3/d],CHD [m3/d],"
         <<"surface-store leakage to aquifer [m3/d],aquifer discharge to surface stores [m3/d],surface-store exchange applied to Raven [m3],"
         <<"reservoir seepage sent [m3/d],reservoir seepage delivered [m3/d],reservoir seepage awaiting delivery [m3],reservoir gain not supplied by MODFLOW [m3]"<<(IsExternal()?",other model packages [m3/d]":"")<<endl;
  if (_aObsWells.size()>0){
    _HEADS.open(FilenamePrepare("GWHeads.csv",Options).c_str());
    _HEADS<<"time [d],date";
    for (size_t o=0;o<_aObsWells.size();o++){_HEADS<<","<<_aObsWells[o].name<<" ("<<_aObsWells[o].id<<") [m]";}
    _HEADS<<endl;
  }
  _CONN.open(FilenamePrepare("GWConnectivity.csv",Options).c_str());
  _CONN<<"time [d],date,connected wet clusters,wet fraction (all)";
  for (size_t p=0;p<_aProfiles.size();p++){_CONN<<",wet fraction "<<_aProfiles[p].name;}
  _CONN<<endl;
  _HRUSTATE.open(FilenamePrepare("GWHRUState.csv",Options).c_str());
  ExitGracefullyIf(_HRUSTATE.fail(),"CGroundwaterModel: cannot open GWHRUState.csv",FILE_OPEN_ERR);
  _HRUSTATE<<"time [d],date";
  for (size_t i=0;i<_aqHRUs.size();i++){_HRUSTATE<<",head_"<<_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID()<<" [m]";}
  for (size_t i=0;i<_aqHRUs.size();i++){_HRUSTATE<<",wtdepth_"<<_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID()<<" [m]";}
  _HRUSTATE<<endl;
}

void CGroundwaterModel::WriteStepOutputs(const double &tstep,const time_struct &tt)
{
  double t=tt.model_time+tstep;
  double scale=max(fabs(_Qrch)+fabs(_Qseep_tot)+fabs(_Qsto)+_Qriv_in+_Qriv_out+fabs(_Qevt)+fabs(_Qwel)+fabs(_Qghb)+fabs(_Qchd)+_Qswx_in+_Qswx_out+fabs(_Qres)+fabs(_Qother),1.0); //gross flow through the aquifer [m3/d]
  string dstr=EndDate(t);
  _BUDGET<<setprecision(10)<<t<<","<<dstr<<","<<_ravenRchVol/tstep<<","<<_Qrch<<","<<_Qseep_tot<<","<<_Qsto<<","
         <<_Qerr<<","<<100.0*_Qerr/scale<<","<<_returnVol<<","<<_lastIter<<","<<(_lastConverged?(_lastRelaxed?2:1):0)<<","<<_rchStoreTot<<","<<_QrchCHD; //line completed in FinishBudgetLine
  _budgetLineOpen=true;
  if (_writeNC && ((_stepCount%max(_headSaveFreq,1))==0)){WriteNetCDF(t);}
  if (_HEADS.is_open()){
    _HEADS<<setprecision(10)<<t<<","<<dstr;
    for (size_t o=0;o<_aObsWells.size();o++){_HEADS<<","<<GetObservationWellHead(_aObsWells[o].id);}
    _HEADS<<endl;
  }
  WriteConnectivity(t,dstr);

  _HRUSTATE<<setprecision(8)<<t<<","<<dstr;
  vector<double> head(_aqHRUs.size(),0.0),depth(_aqHRUs.size(),0.0);
  for (size_t i=0;i<_aqHRUs.size();i++){
    int k=_aqHRUs[i]; double sa=0,sh=0,sd=0;
    for (size_t j=0;j<_hruLinks[k].size();j++){
      const gw_link &L=_aLinks[_hruLinks[k][j]];
      double h=WaterTable(L.ic); //physical water table (not the numerically arbitrary head of a dry cell)
      sa+=L.area; sh+=L.area*h; sd+=L.area*(_land[L.ic]-h);
    }
    head[i]=sh/sa; depth[i]=sd/sa;
  }
  for (size_t i=0;i<head.size(); i++){_HRUSTATE<<","<<head[i];}
  for (size_t i=0;i<depth.size();i++){_HRUSTATE<<","<<depth[i];}
  _HRUSTATE<<endl;
}

void CGroundwaterModel::WriteSummary(const optStruct &Options)
{
  ofstream S(FilenamePrepare("GWModelSummary.txt",Options).c_str());
  int nact=0,npass=0,ninact=0,ncols=0;
  for (int n=0;n<_nlay*_ncpl;n++){if (_idomain[n]>0){nact++;} else if (_idomain[n]<0){npass++;}}
  for (int ic=0;ic<_ncpl;ic++){if (_colProfile[ic]!=DOESNT_EXIST){ncols++;} else {ninact++;}}
  S<<"Raven-MODFLOW 6 groundwater model summary"<<endl;
  S<<"  MODFLOW 6 files       : "<<_mf6dir<<endl;
  S<<"  MODFLOW 6 library     : "<<MF6LibraryPath()<<((_mf6Version=="")?string(""):" (version "+_mf6Version+")")<<endl;
  const char *gname[4]={"regular","quadtree","HRU mesh","imported"}; const char *dname[3]={"DIS","DISV","DISU"};
  if (IsExternal()){ //existing MODFLOW 6 model
    const char *tu=(_extTf==86400.0)?"seconds":(_extTf==1440.0)?"minutes":(_extTf==24.0)?"hours":(_extTf==1.0)?"days":"years";
    const char *lu=(fabs(_extLf-1.0)<1e-12)?"metres":(fabs(_extLf-100.0)<1e-9)?"centimetres":"feet";
    S<<"  existing model        : "<<_mname<<" in "<<_extSim<<" ("<<dname[_disType]<<", "<<_nlay<<" layers x "<<_ncpl<<" cells per layer; "<<lu<<", "<<tu<<")"<<endl;
    S<<"  cells linked by       : "<<(_extCRS.IsSet()?"cell outlines from MODFLOW, "+_extCRS.Describe():(_extGridFile!="")?"cell polygons in "+_extGridFile+" (layout matches MODFLOW's; area ratio "+ToStr(_extAreaCheck)+")":string(":OverlapWeights table"))<<endl;
    S<<"  packages taken over   :"; for (size_t i=0;i<_takeOver.size();i++){S<<" "<<_takeOver[i];} S<<((_takeOver.empty())?" none":"")<<endl;
    S<<"  packages in balance   :"; for (size_t i=0;i<_extOther.size();i++){S<<" "<<_extOther[i];} S<<endl;
    if (_extSteady){S<<"  time                  : steady-state model with one stress period, stretched over Raven's run: every Raven step is a steady-state solution (no storage)"<<endl;}
    if (_extMoverDropped>0){S<<"  water mover           : "<<_extMoverDropped<<" movers from or to packages taken over by Raven removed"<<endl;}
    S<<"  stream package        : "<<((_streamPkg=="")?string("none"):_streamPkg+((_streamMap=="AUX")?" (subbasin from auxiliary "+_streamAux+")":" (subbasin of the dominant HRU)"))<<endl;
    double Ahru=0,Alink=0; int nlow=0;
    bool poly=_owTable.empty(); //with polygons, coverage is measured against each HRU's polygon; otherwise against its .rvh area
    for (size_t i=0;i<_aqHRUs.size();i++){int k=_aqHRUs[i]; double A=poly?_hruPolyArea[k]:_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2; Ahru+=A; Alink+=_hruLinkedArea[k]; if (_hruLinkedArea[k]<0.5*A){nlow++;}}
    S<<"  HRU coverage          : "<<setprecision(4)<<100.0*Alink/max(Ahru,1.0)<<"% of the coupled HRUs' "<<(poly?"polygon":".rvh")<<" area lies on active model cells; "<<nlow<<" HRUs below 50%"<<setprecision(6)<<endl;
    if (_extDvclose>0.01){S<<"  note                  : the model's solver tolerance ("<<_extDvclose<<" m) is loose for exact water accounting; 0.001 m or less is recommended"<<endl;}
  }
  else {
  if (_gridType==GWGRID_REGULAR){S<<"  grid                  : regular, "<<_nlay<<" layers x "<<_nrow<<" rows x "<<_ncol<<" cols, cell size "<<(_grid.uniform?_cellSize:_grid.MinCellSize())<<(_grid.uniform?"":"-"+ToStr(_cellSize))<<" m";}
  else {S<<"  grid                  : "<<gname[_gridType]<<", "<<_nlay<<" layers x "<<_ncpl<<" cells, cell size "<<_grid.MinCellSize()<<"-"<<_cellSize<<" m";}
  S<<((_grid.angle!=0.0)?", rotated "+ToStr(_grid.angle)+" deg":"")<<"; MODFLOW "<<dname[_disType]<<(UseXT3D()?((_flowCorr=="XT3D")?" with XT3D":" with XT3D (RHS)"):"")<<(_connectByElevation?", layers connected by elevation":"")<<endl;
  S<<"  HRU fidelity of cells : "<<setprecision(4)<<100.0*_meshFidelity<<"% of active-cell HRU area lies in each cell's dominant HRU"<<setprecision(6)<<endl;
  }
  S<<"  active columns        : "<<ncols<<" (inactive "<<ninact<<")"<<endl;
  S<<"  active cells          : "<<nact<<"   vertical pass-through cells: "<<npass<<endl;
  S<<"  coupled HRUs          : "<<_aqHRUs.size()<<endl;
  S<<"  recharge/seepage cells: "<<_bndCol.size()<<endl;
  if (IsExternal()){ //rivers, wells and boundaries are the model's own packages (listed above); Raven adds none
    S<<"  river/well/boundary cells: those of the model's own packages listed above (none added by Raven)"<<endl;
  }
  else {
  S<<"  river cells           : "<<_rivNode.size()<<endl;
  S<<"  water-table ET cells  : "<<_evtCol.size()<<endl;
  S<<"  wells / WEL cells     : "<<_aWells.size()<<" / "<<_welNode.size()<<endl;
  S<<"  GHB cells / CHD cells : "<<_ghbNode.size()<<" / "<<_chdNode.size()<<endl;
  }
  S<<"  observation wells     : "<<_aObsWells.size()<<endl;
  if (!IsExternal()){
  S<<"  surface-water exchange cells: "<<_swxCol.size()<<endl;
  S<<"  coupled reservoirs / cells  : "<<_resP.size()<<" / "<<_resNode.size()<<endl;
  }
  S<<"  HRU/grid linkage      : "<<(_cacheHit?"read from cache":"computed")<<" in "<<_geomSeconds<<" s"<<endl;
  for (size_t o=0;o<_aObsWells.size();o++){
    int n=_aObsWells[o].node;
    S<<"    "<<_aObsWells[o].name<<": "<<((n<0)?string("outside model"):"layer "+to_string(n/_ncpl+1)+", MODFLOW cell ("+CellID(n/_ncpl,n%_ncpl)+")")<<endl;
  }
  if (_doSteady){S<<"  initialization        : steady state, recharge "<<_ssRecharge<<" mm/d"<<endl;}
  else if (_hotHeads.size()>0){S<<"  initialization        : hotstart heads from .rvc"<<endl;}
  for (size_t p=0;p<_aProfiles.size();p++){
    int n=0; for (int ic=0;ic<_ncpl;ic++){if (_colProfile[ic]==(int)p){n++;}}
    S<<"  profile "<<_aProfiles[p].name<<": "<<_aProfiles[p].layers.size()<<" layers, "<<n<<" columns"<<endl;
  }
  S<<endl<<"HRU_ID,profile,rvh_area_km2,polygon_area_km2,linked_area_km2,linked_fraction_of_polygon,n_cells"<<endl;
  for (size_t i=0;i<_aqHRUs.size();i++){
    int k=_aqHRUs[i]; CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
    S<<pHRU->GetHRUID()<<","<<_aProfiles[_hruProfile[k]].name<<","<<pHRU->GetArea()<<","<<_hruPolyArea[k]/M2_PER_KM2<<","
     <<_hruLinkedArea[k]/M2_PER_KM2<<","<<_hruLinkedArea[k]/max(_hruPolyArea[k],1e-9)<<","<<_hruLinks[k].size()<<endl;
  }
}

void CGroundwaterModel::CloseOutputs()
{
  FinishBudgetLine();
  if (_BUDGET.is_open()){_BUDGET.close();}
  if (_HRUSTATE.is_open()){_HRUSTATE.close();}
  if (_CONN.is_open()){_CONN.close();}
  if (_HEADS.is_open()){_HEADS.close();}
#ifdef netcdf
  if (_ncid>=0){nc_close(_ncid); _ncid=-1;}
#endif
  if (_nNonConverged>0){
    cout<<"WARNING: MODFLOW 6 failed to converge on "<<_nNonConverged<<" time steps (see GWBudget.csv)"<<endl;
    _nNonConverged=0;
  }
  _mf6.Finalize();
  _X=NULL; _active=false; //MODFLOW memory is released: nothing may read heads after this point
}
