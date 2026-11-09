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
static const double SZ_FAR_RESID=0.01; ///< largest MODFLOW residual [-] of the gross aquifer flow with which an unconverged step holding a head far
                                       ///< outside the model is carried on
static double gw_prof_t0=-1;
#define GWPROF(name) { if (getenv("RAVEN_GW_PROFILE")!=NULL){ double tn=double(clock())/CLOCKS_PER_SEC; if (gw_prof_t0<0){gw_prof_t0=tn;} cerr<<"[groundwater] "<<name<<" done at "<<tn<<" s"<<endl; } }
#include "Model.h"
#include "HydroProcessABC.h"
#include "SoilWaterMovers.h"
#include "LateralExchangeABC.h"
#include "ProcessGroup.h"
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
  _Qother=0.0; _extAux=-1; _extNaux=0; _extOffset=0.0; _extDrop=0; _extNModels=0; _extMoverDropped=0; _extTSShifted=0; _extSteady=false; _extNoSto=false; _mf6Version=""; _stoName="STO"; _disName="DIS";
  _iGW=_iReturnSV=DOESNT_EXIST;
  _maxiter=500;
  _X=_rchBound=_rchHcof=_rchRhs=_drnHcof=_drnRhs=_strgss=_strgsy=NULL;
  _rivStage=_rivCond=_rivRbot=_rivHcof=_rivRhs=NULL; _rivNcol=0;
  _evtRate=_evtHcof=_evtRhs=NULL; _evtNcol=0; _Qevt=0.0;
  _welQ=_welHcof=_welRhs=_ghbHead=_ghbCond=_ghbHcof=_ghbRhs=_chdHead=_chdSim=NULL;
  _welNcol=_ghbNcol=_chdNcol=0; _Qwel_spec=_Qwel=_Qghb=_Qchd=0.0;
  _chdNb=_chdNl=_chdIb=NULL; _chdIbc=0; _chdPacked=false; _bndOffNow=0; _bndOffN=0;
  _steadyInit=false; _doSteady=false; _ssRecharge=0.0; _initMode="FROM_PROFILE"; _initValue=0.0;
  _hotDims[0]=_hotDims[1]=_hotDims[2]=0; _active=false;
  _rchStoreTot=0.0; _QrchCHD=0.0; _rivBedDepth=-1.0;
  _swxHead=_swxCond=_swxRbot=_swxHcof=_swxRhs=NULL; _swxNcol=0; _Qswx_in=_Qswx_out=_swxAppliedVol=0.0;
  _resExchange=false; _resQ=_resHcof=_resRhs=NULL; _resNcol=0; _Qres_spec=_Qres=_resPendingVol=_resShortVol=0.0;
  _cacheHit=false; _geomSeconds=0.0; _writeNC=false; _writeCBC=false; _ncid=-1; _ncTime=_ncHead=-1; _ncRec=0; _nstp=0; _stepCount=0; _jday0=0; _jyear0=0; _calendar=0; _tsDur=_tsDt=0.0;
  _Qriv_in=_Qriv_out=_riverAppliedVol=_riverDeficitVol=0.0; _budgetLineOpen=false;
  _ravenRchVol=_Qrch=_Qseep_tot=_Qsto=_Qerr=_returnVol=0.0;
  _nNonConverged=0;_lastIter=0;_lastConverged=true;
  _szMode=0; _szNLayers=0; _szMaxIter=100; _szTolH=1e-4; _szTolQ=1e-6; _szSubgrid=5;
  _vzMode=0; _vzPar=0.0; _vzQ0=0.0; _vzClsDev=0.0; _vzBackV=_vzBackWV=0.0; _vzClasses=1; _vzAcc=0.02; _vzN=16; _szCapWT=0.0; _szwQ=NULL; _szwNcol=0; _vzUpQ=0.0; _vzSlowed=0; _ibPtr=NULL; _mf6Newton=true; _farWarned=false; _szCapRise=-1; _szTiming=1; _szClasses=1; _szClsWarned=false; _szDebtVol=_szRepayQ=0.0; _QszCHD_in=_QszCHD_out=_szRchCHD=_QrchCHD0=0.0; _vzInVol=_vzOutVol=_vzCapVol=_vzStoreVol=_vzStore0=0.0; _vzSub=_vzNewton=0;
  _szxHead=_szxCond=_szxHcof=_szxRhs=NULL; _szxNcol=0; _szqQ=_szqHcof=_szqRhs=NULL; _szqNcol=0; _szdElev=_szdCond=_szdAux=_szdHcof=_szdRhs=NULL; _szdNcol=0;
  _Qsz_in=_Qsz_out=_szRejVol=_szOvfVol=_szAppliedVol=_szShortVol=0.0; _szIter=0; _szResid=0.0; _szConv=true; _szNotConv=0; _szShortN=0; _szRchByArea.clear(); _szCapUpd=false; _szOneOK=true; _szMovedV=0.0; _szWithheldV=0.0; _szWithheldG=0.0; _szWithheldS=0.0; _szGuardN=0; _ncWithheldVol=0.0; _ncWithheldN=0;
  _szTstep=1.0; _szPrevQ=0.0; _szRatioMin=_szRatioMax=0.0; _szColInit=false; _szInitVol=0.0; _szInitHot=0; _szResE=-1; _szDQ=_szDQtol=0.0; _szExtra=0;
  _szCapN=0; _szCapExcess=0.0; _szDropN=_szDropPlain=0; _szDropSteps=0; _szDropThis=false; _szxDropped.clear(); _szEsum=0.0; _szWhy=0; _szWhyE=-1; _szWhyV=_szWhyFa=_szWhyFb=0.0;
  _szNegCutQ=_szNegUnheldV=0.0; _negCutVol=_negUnheldVol=0.0; _negUnheldN=0; _negWarned=false; _Qgross=1.0; _QabsOther=0.0; _rchLostQ=_rchLostVol=0.0; _rchLostN=0;
  _resUndelVol=_resUndelTot=0.0; _resUndelN=0; _resUnsent.clear(); _resUse.clear(); _resWs.clear(); _resOff.clear(); _resShortTot=0.0; _resShortN=0;
  _vzChdV=0.0; _vzChdN=0; _vzHotUnusedV=0.0; _vzHotUnusedN=0;
  _wetKeptIt=_wetKeptSteps=0; _farSlowN=0; _farWarned=false; _xPrevIt.clear(); _xStart.clear(); _keptCol.clear(); _stepEnd=""; _lastConverged=true;
  _summaryFile="";
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
  _chdNb=_chdNl=_chdIb=NULL; _chdIbc=0; _chdPacked=false; _chdAct.clear(); _bndOffNow=0; _bndOffN=0;
  _swxCol.clear(); _swxX.clear(); _swxLinks.clear(); _swxVl.clear(); _swxAppliedVol=0.0;
  for (size_t r=0;(r<_resP.size()) && (r<_resK0.size());r++){ //(a reservoir's seepage constant is scaled by the share of its lake that can
    CReservoir *pRes=_pModel->GetSubBasin(_resP[r])->GetReservoir();      // take seepage: its own value for the next ensemble member)
    if (pRes!=NULL){pRes->SetGWParameters(_resK0[r],pRes->GetLocalGWHead());}
  }
  _resP.clear(); _resZbot.clear(); _resNode.clear(); _resR.clear(); _resW.clear(); _resK0.clear(); _resA.clear(); _resWk.clear();
  _X=NULL; _active=false; _stepCount=0; _nNonConverged=0; _budgetLineOpen=false;
  _ravenRchVol=_Qrch=_Qseep_tot=_Qsto=_Qerr=_returnVol=_riverAppliedVol=_riverDeficitVol=0.0;
  _szOn.clear(); _szSV.clear(); _szLay.clear(); _szD.clear(); _szLam.clear(); _szRchV.clear(); _szRej.clear(); _szOvfV.clear(); _szNetV.clear(); _szW.clear();
  _linkF.clear(); _cellDem.clear(); _szxLink.clear(); _szxC.clear(); _szxZb.clear(); _szxQ.clear(); _colCoverNS.clear(); _szxOfHRU.clear();
  _szxCell.clear(); _szxQc.clear(); _szxCg.clear(); _szxHc.clear(); _szxHp.clear(); _szxHr.clear(); _szxHq.clear(); _szxQq.clear(); _szxFw.clear(); _szxFc.clear(); _szxOne.clear(); _szSyncOvf.clear(); _szxBa.clear(); _szxBra.clear(); _szxBb.clear(); _szxBrb.clear(); _szxBk.clear(); _szxLam.clear(); _szxRej.clear(); _szxCol.clear(); _linkSzx.clear(); _szColOff.clear(); _szColA.clear(); _szColPhi.clear(); _szColInit=false;
  _szPhiProc.clear(); _szPhiOff.clear(); _evtBase.clear();
  _szxCapQ.clear(); _szxCapL.clear(); _szCapWT=0.0; _szxUp.clear(); _szxUpF.clear(); _szxUpCHD.clear(); _ibPtr=NULL; _vzUp.clear(); _vzUpA.clear(); _vzUpQ=0.0; _szwStart.clear(); _szwNode.clear(); _szwQ=NULL; _szwNcol=0;
  _szCls.clear(); _szClsW.clear(); _szClsPhi.clear(); _szClsRes.clear(); _szClsGW.clear(); _szClsDone.clear();
  _szcLink.clear(); _szcCol.clear(); _szcCell.clear(); _szcRej.clear(); _szcFw.clear(); _szcQ.clear(); _szcCapL.clear(); _szcOfHRU.clear(); _linkSzc.clear(); _QszCHD_in=_QszCHD_out=_szRchCHD=_QrchCHD0=0.0;
  _szxQprev.clear(); _szxLagL.clear(); _szxLagD.clear(); _szLagNet.clear(); _szxDebt.clear(); _szxRepay.clear(); _szDebtVol=_szRepayQ=0.0;
  _vzcCol.clear(); _vzOff.clear(); _vzcA.clear(); _vzcZt.clear(); _vzcZbs.clear(); _vzcUse.clear(); _vzcF0.clear(); _vzcOut.clear(); _vzcRq.clear(); _vzAbsN.clear(); _vzQN.clear(); _vzClsDev=0.0; _vzBackW.clear(); _vzPass.clear(); _vzRl.clear(); _vzBackV=_vzBackWV=0.0;
  _vzStrat.clear(); _vzUse.clear(); _vzQ.clear(); _vzAbsQ.clear(); _vzRejV.clear(); _vzFp.clear(); _vzInV.clear(); _vzCapV.clear(); _vzCapL.clear();
  _vzInVol=_vzOutVol=_vzCapVol=_vzStoreVol=_vzStore0=0.0; _vzSub=_vzNewton=0; _vzSlowed=0;
  _szxHead=_szxCond=_szxHcof=_szxRhs=NULL; _szxNcol=0; _szqQ=_szqHcof=_szqRhs=NULL; _szqNcol=0; _szdElev=_szdCond=_szdAux=_szdHcof=_szdRhs=NULL; _szdNcol=0;
  _Qsz_in=_Qsz_out=_szRejVol=_szOvfVol=_szAppliedVol=_szShortVol=0.0; _szIter=0; _szResid=0.0; _szConv=true; _szNotConv=0; _szShortN=0; _szRchByArea.clear(); _szCapUpd=false; _szOneOK=true; _szMovedV=0.0; _szWithheldV=0.0; _szWithheldG=0.0; _szWithheldS=0.0; _szGuardN=0; _ncWithheldVol=0.0; _ncWithheldN=0;
  //(every soil-zone, column and stop-rule member and run total: ensemble members start as the first did)
  _szCapN=0; _szCapExcess=0.0; _szDropN=_szDropPlain=0; _szDropSteps=0; _szDropThis=false; _szxDropped.clear(); _szEsum=0.0; _szWhy=0; _szWhyE=-1; _szWhyV=_szWhyFa=_szWhyFb=0.0;
  _szNegCutQ=_szNegUnheldV=0.0; _negCutVol=_negUnheldVol=0.0; _negUnheldN=0; _negWarned=false; _Qgross=1.0; _QabsOther=0.0; _rchLostQ=_rchLostVol=0.0; _rchLostN=0;
  _resUndelVol=_resUndelTot=0.0; _resUndelN=0; _resUnsent.clear(); _resUse.clear(); _resWs.clear(); _resOff.clear(); _resShortTot=0.0; _resShortN=0;
  _vzChdV=0.0; _vzChdN=0; _vzHotUnusedV=0.0; _vzHotUnusedN=0;
  _wetKeptIt=_wetKeptSteps=0; _farSlowN=0; _farWarned=false; _xPrevIt.clear(); _xStart.clear(); _keptCol.clear(); _stepEnd=""; _lastConverged=true;
  _szxRs.clear(); _szxRq.clear(); _szxRqP.clear(); _szxRl.clear(); _szxRqc.clear(); _szxF0.clear(); _szxDrn.clear(); _szxKink.clear(); _szxNegCap.clear(); _vzUsed.clear(); _cellConv.clear(); _szResE=-1; _szExtra=0;
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
//////////////////////////////////////////////////////////////////
/// \brief reads a :LandSurface or :BedrockSurface raster (ESRI ASCII, or GeoTIFF for .tif/.tiff) in the coordinates of
/// :InputCRS (else longitude/latitude); stops on a read error or a raster left in longitude/latitude under :InputCRS
//
void CGroundwaterModel::ReadRaster(CGWRaster &R,const string &file,const optStruct &Options) const
{
  string err;
  R.SetCRS(InCRS());
  if (!R.Read(file,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
  istringstream notes(R.Note()); string note;
  while (getline(notes,note)){if (note!=""){WriteWarning("CGroundwaterModel: "+note,Options.noisy);}}
  double rb[4]; R.Extent(rb); CheckInputLonLat(file,rb,"raster");
}
//////////////////////////////////////////////////////////////////
/// \brief a raster must give values for the active columns: none (a raster of another area, a wrong corner or
/// projection, or NODATA over the whole model) stops the run with both extents; some missing gives a warning
//
void CGroundwaterModel::CheckRasterCover(const CGWRaster &R,const string &file,const string &cmd,const vector<int> &cols,const int nCov,
                                         const string &fallback,const optStruct &Options) const
{
  const int nAct=(int)cols.size();
  if ((nAct==0) || (nCov==nAct)){return;}
  if (nCov>0){
    WriteWarning("CGroundwaterModel: "+cmd+" raster "+file+" gives no value for "+to_string(nAct-nCov)+" of the "+to_string(nAct)+
                 " active columns (outside the raster, or NODATA there); these take "+fallback,Options.noisy);
    return;
  }
  double rb[4],mb[4]={1e300,1e300,-1e300,-1e300}; R.Extent(rb);
  const CGWCRS *crs=InCRS();
  for (size_t i=0;i<cols.size();i++){ //(the points the columns' means sample, in the raster's coordinates)
    int ic=cols[i]; vector<gw_pt> pts; _grid.SamplePoints(ic,5,pts);
    for (size_t k=0;k<pts.size();k++){
      double X,Y,lon,lat,x,y; _grid.ToWorld(pts[k].x,pts[k].y,X,Y); _proj.Inverse(X,Y,lon,lat);
      if (crs!=NULL){crs->Forward(lon,lat,x,y);} else {x=lon; y=lat;}
      mb[0]=min(mb[0],x); mb[1]=min(mb[1],y); mb[2]=max(mb[2],x); mb[3]=max(mb[3],y);
    }
  }
  ostringstream S; S<<setprecision(9)<<"CGroundwaterModel: "<<cmd<<" raster "<<file<<" gives no value for any of the "<<nAct<<" active columns: the raster covers x "
    <<rb[0]<<" to "<<rb[2]<<", y "<<rb[1]<<" to "<<rb[3]<<", and the model's active columns lie at x "<<mb[0]<<" to "<<mb[2]<<", y "<<mb[1]<<" to "<<mb[3]
    <<" ("<<((crs!=NULL)?"coordinates of :InputCRS":"longitude/latitude")<<"). Check that the file covers the model area, its corner coordinates and :InputCRS"
    <<" (a raster that does overlap the model has only NODATA there)";
  ExitGracefully(S.str().c_str(),BAD_DATA);
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
  if (_inCRS.IsSet()){s<<":incrs:"<<_inCRS.Describe();}
  return s.str();
}
//////////////////////////////////////////////////////////////////
/// \brief start-up diagnostic: a layer of a generated model whose active convertible cells all start at or below their
///  bottom (initial heads, cold start without a steady-state period, or hotstart heads) is dry as a whole; MODFLOW then has
///  to rewet it in the first steps, which can take many iterations or leave steps unconverged
//
void CGroundwaterModel::CheckDryStart(const optStruct &Options)
{
  _dryStart="";
  if (IsExternal() || _doSteady){return;}
  for (int l=0;l<_nlay;l++){
    int nact=0,ndry=0;
    for (int ic=0;ic<_ncpl;ic++){
      const int ip=_colProfile[ic],n=UserNode(l,ic);
      if ((ip==DOESNT_EXIST) || (_idomain[n]<=0) || (l>=(int)_aProfiles[ip].layers.size())){continue;}
      if (_aProfiles[ip].layers[l].type==GWL_CONFINED_AQUIFER){continue;} //(confined cells never dry)
      nact++; if (_strt[n]<=_botm[n]){ndry++;}
    }
    if ((nact==0) || (ndry<nact)){continue;}
    _dryStart+=(_dryStart==""?"":", ")+to_string(l+1)+" ("+to_string(nact)+" cells)";
    WriteWarning("CGroundwaterModel: layer "+to_string(l+1)+" starts dry: the initial head is at or below the bottom of all its "+to_string(nact)+
                 " active convertible cells"+string((_hotHeads.size()>0)?" (hotstart heads)":" (INITIAL_HEAD_DEPTH or :InitialGWHeads deeper than the layers above it)")+
                 "; MODFLOW has to rewet the whole layer, which can take many iterations or leave the first steps unconverged. Start from "
                 ":GWInitialization STEADY_STATE, a smaller INITIAL_HEAD_DEPTH or :InitialGWHeads",Options.noisy);
  }
}
//wells and monitoring wells given as X,Y (:InputCRS) to longitude/latitude
void CGroundwaterModel::ResolveInputCRS()
{
  for (size_t w=0;w<_aWells.size()+_aObsWells.size();w++){
    bool ow=(w>=_aWells.size()); double &lat=ow?_aObsWells[w-_aWells.size()].lat:_aWells[w].lat,&lon=ow?_aObsWells[w-_aWells.size()].lon:_aWells[w].lon;
    bool &xy=ow?_aObsWells[w-_aWells.size()].xy:_aWells[w].xy;
    if (!xy){continue;}
    ExitGracefullyIf(!_inCRS.IsSet(),("CGroundwaterModel: well "+(ow?_aObsWells[w-_aWells.size()].name:_aWells[w].name)+" is given by X and Y: give :InputCRS (their projection) in the .rvg file, or LATITUDE and LONGITUDE").c_str(),BAD_DATA);
    double x=lon,y=lat; _inCRS.Inverse(x,y,lon,lat); xy=false;
  }
}
//:InputCRS: a GeoJSON or raster input whose coordinates (box: xmin, ymin, xmax, ymax as written in the file) lie in the
//longitude/latitude range, over the HRUs' LATITUDE/LONGITUDE of the .rvh file, and not near them in the :InputCRS
//coordinates, is in longitude/latitude: stop (read as projected metres it would lie by the projection's false origin)
void CGroundwaterModel::CheckInputLonLat(const string &file,const double *box,const string &what) const
{
  if (!_inCRS.IsSet()){return;}
  if ((box[0]>box[2]) || (box[0]<-180.0) || (box[2]>180.0) || (box[1]<-90.0) || (box[3]>90.0)){return;}
  double L[4]={1e300,1e300,-1e300,-1e300},P[4]={1e300,1e300,-1e300,-1e300}; int n=0;
  int nH=(int)_aqHRUs.size(); if (nH==0){nH=_pModel->GetNumHRUs();}
  for (int i=0;i<nH;i++){
    location loc=_pModel->GetHydroUnit(_aqHRUs.empty()?i:_aqHRUs[i])->GetCentroid();
    double lo=loc.longitude,la=loc.latitude,x,y;
    if (!std::isfinite(lo) || !std::isfinite(la) || (fabs(la)>90.0) || (fabs(lo)>360.0)){continue;}
    _inCRS.Forward(lo,la,x,y); if (!std::isfinite(x) || !std::isfinite(y)){continue;}
    L[0]=min(L[0],lo); L[1]=min(L[1],la); L[2]=max(L[2],lo); L[3]=max(L[3],la);
    P[0]=min(P[0],x);  P[1]=min(P[1],y);  P[2]=max(P[2],x);  P[3]=max(P[3],y); n++;
  }
  if (n==0){return;}
  double mL=0.5+0.5*max(L[2]-L[0],L[3]-L[1]),mP=50000.0+0.5*max(P[2]-P[0],P[3]-P[1]); //(HRUs extend beyond their centroids)
  bool overLL=(box[0]<=L[2]+mL) && (box[2]>=L[0]-mL) && (box[1]<=L[3]+mL) && (box[3]>=L[1]-mL);
  bool overP =(box[0]<=P[2]+mP) && (box[2]>=P[0]-mP) && (box[1]<=P[3]+mP) && (box[3]>=P[1]-mP);
  if (overLL && !overP){
    ostringstream S; S<<setprecision(6)<<"CGroundwaterModel: "<<what<<" "<<file<<" looks like longitude/latitude: its coordinates (x "<<box[0]<<" to "<<box[2]
      <<", y "<<box[1]<<" to "<<box[3]<<") lie over the HRUs' LATITUDE/LONGITUDE, while in the coordinates of :InputCRS ("<<_inCRS.Describe()
      <<") the HRUs lie at x "<<setprecision(9)<<P[0]<<" to "<<P[2]<<", y "<<P[1]<<" to "<<P[3]<<". With :InputCRS every GeoJSON and raster input must be "
      <<"in that projection: reproject the file (or remove :InputCRS if every input file is in longitude/latitude)";
    ExitGracefully(S.str().c_str(),BAD_DATA);
  }
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
      double rb[4];
      if (!GWGeom::ReadGeoJSONLines(file,(_refKind[i]=="RIVERS")?_riverIDField:"",lines,ids,err,InCRS(),rb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
      CheckInputLonLat(file,rb,"GeoJSON file");
      for (size_t l=0;l<lines.size();l++){
        f.kind=1; f.geom.clear();
        for (size_t k=0;k<lines[l].size();k++){double X,Y; gw_pt p; _proj.Forward(lines[l][k].x,lines[l][k].y,X,Y); _grid.ToLocal(X,Y,p.x,p.y); f.geom.push_back(p);}
        if (f.geom.size()>=2){ref.push_back(f);}
      }
    }
    else if (_refKind[i]=="POLYGONS"){
      map<long long,vector<gw_ring> > polys;
      double rb[4];
      if (!GWGeom::ReadGeoJSONPolygons(_refFile[i],"",polys,err,InCRS(),rb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
      CheckInputLonLat(_refFile[i],rb,"GeoJSON file");
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
  else { //(as for :AquiferClasses: a misspelt parameter would silently keep its default)
    ExitGracefully(("CGroundwaterModel: unknown aquifer profile parameter "+param+" in :AquiferProfileParameters (INITIAL_HEAD_DEPTH, DEFAULT_BEDROCK_DEPTH, "
      "RIVERBED_K, RIVERBED_THICKNESS, SEEPAGE_LEAKANCE, RECHARGE_DELAY, SOIL_ZONE_DEPTH, EXTINCTION_DEPTH, SEEPAGE_SMOOTHING_DEPTH)").c_str(),BAD_DATA);
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
  _tsDur=Options.duration; _tsDt=Options.timestep;
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
  if (SZPhysical()){SZSetup();}

  CheckProcesses();
  for (map<string,CTimeSeries*>::iterator it=_headTS.begin();it!=_headTS.end();it++){ //(before BuildStresses: the boundaries are placed with the
    it->second->Initialize(Options.julian_start_day,Options.julian_start_year,Options.duration,Options.timestep,false,Options.calendar); // series' heads)
  }
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
  //MODFLOW ignores recharge on specified-head cells: recharge is spread over each HRU's other cells
  _colCHDTop.assign(_ncpl,0);
  for (size_t e=0;e<_chdNode.size();e++){int ic=_chdNode[e]%_ncpl; bool top=true;
    for (int l=0;l<_chdNode[e]/_ncpl;l++){if (_idomain[UserNode(l,ic)]>0){top=false;}} if (top){_colCHDTop[ic]=1;}}
  _hruRchArea.assign(_pModel->GetNumHRUs(),0.0);
  for (size_t li=0;li<_aLinks.size();li++){if (!_colCHDTop[_aLinks[li].ic]){_hruRchArea[_aLinks[li].k]+=_aLinks[li].area;}}
  if (SZPhysical()){SZBuildEntries();}
  //surface-water store <-> aquifer exchanges: one head-dependent cell per column holding group HRUs
  for (size_t x=0;x<_aSWX.size();x++){
    gw_swx &S=_aSWX[x]; int lay=0;
    sv_type typ=_pModel->GetStateVarInfo()->StringToSVType(S.sv,lay,true);
    S.iSV=_pModel->GetStateVarIndex(typ,lay);
    ExitGracefullyIf(S.iSV==DOESNT_EXIST,("CGroundwaterModel: :SurfaceWaterExchange store "+S.sv+" is not used by this model").c_str(),BAD_DATA);
    CHRUGroup *pG=_pModel->GetHRUGroup(S.group);
    ExitGracefullyIf(pG==NULL,("CGroundwaterModel: :SurfaceWaterExchange HRU group "+S.group+" does not exist").c_str(),BAD_DATA);
    if (SZPhysical() && (typ==SOIL)){ //the soil zone's layers exchange with the aquifer through SZX_RAVEN already
      for (size_t i=0;i<_aqHRUs.size();i++){
        int k=_aqHRUs[i]; bool inZone=false;
        for (size_t j=0;j<_szSV[k].size();j++){if (_szSV[k][j]==S.iSV){inZone=true;}}
        ExitGracefullyIf(_szOn[k] && inZone && pG->IsInGroup(k),("CGroundwaterModel: :SurfaceWaterExchange on "+S.sv+" of HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+
          ": with :SoilZone PHYSICAL that layer exchanges with the aquifer through the soil zone already").c_str(),BAD_DATA);
      }
    }
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
  //a reservoir's lake HRU must lie in its own subbasin: Raven sends the HRU's surface water to the reservoir of the HRU's
  //subbasin, and there is none (Raven crashes in its first step, with or without groundwater)
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    const CReservoir *pRes=_pModel->GetSubBasin(p)->GetReservoir(); if (pRes==NULL){continue;}
    const int k=pRes->GetHRUIndex(); if ((k==DOESNT_EXIST) || !_pModel->GetHydroUnit(k)->IsEnabled()){continue;}
    const int pk=_pModel->GetHydroUnit(k)->GetSubBasinIndex();
    ExitGracefullyIf((pk!=p) && (_pModel->GetSubBasin(pk)->GetReservoir()==NULL),("CGroundwaterModel: the lake HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+
      " of the reservoir in subbasin "+to_string(_pModel->GetSubBasin(p)->GetID())+" lies in subbasin "+to_string(_pModel->GetSubBasin(pk)->GetID())+
      ", which has no reservoir: give the reservoir an :HRUID of its own subbasin").c_str(),BAD_DATA);
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
      int r=(int)_resP.size(); double wsum=0,wall=0; size_t n0=_resNode.size();
      for (size_t j=0;j<_hruLinks[k].size();j++){
        const gw_link &L=_aLinks[_hruLinks[k][j]]; wall+=L.area; if (_colCHDTop[L.ic]){continue;}
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
      //(Raven sets the reservoirs' initial stages after the groundwater model is built, or in its first routing: here
      // their dry stage, which is on the same datum, and again with the stage after the first step)
      ExitGracefullyIf((zbot<1e299) && (max(pRes->GetResStage(),pRes->GetDryStage())<zbot),("CGroundwaterModel: reservoir in subbasin "+nm+
        " has stages from "+to_string(pRes->GetDryStage())+" m, below the bottom of the aquifer under its lake ("+to_string(zbot)+
        " m): coupled reservoirs need absolute stages on the aquifer's datum (e.g. :AbsoluteCrestHeight)").c_str(),BAD_DATA);
      if ((zbot<1e299) && (pRes->GetResStage()<zbot)){
        WriteWarning("CGroundwaterModel: reservoir in subbasin "+nm+" starts at stage "+to_string(pRes->GetResStage())+" m, below the aquifer under its lake "
                     "(no :InitialReservoirStage or :InitialReservoirFlow in the .rvc file, or a stage below it): Raven's first step fills it, a jump in Raven's water balance",false);
      }
      _resP.push_back(p); _resZbot.push_back(zbot);
      //the lake's bed over specified-head columns takes no seepage: its seepage constant is scaled by the share left (each step
      //also by the share of its cells that can take seepage), so the cells left get the seepage of their own share of the bed
      _resK0.push_back(pRes->GetSeepageConstant()); _resA.push_back((wall>0.0)?wsum/wall:1.0);
      map<long long,double>::const_iterator hf=_hotResFrac.find(_pModel->GetSubBasin(p)->GetID());
      _resWk.push_back((hf!=_hotResFrac.end())?min(max(hf->second,0.0),1.0):1.0);
    }
  }
  _doSteady=_steadyInit;
  if (_doSteady && (_hotHeads.size()>0)){
    WriteAdvisory("CGroundwaterModel: hotstart heads found in .rvc; steady-state initialization skipped",Options.noisy);
    _doSteady=false; //decided per run: the user's setting is kept for later ensemble members
  }
  CheckDryStart(Options);
  if (IsExternal()){ExtWritePackages(Options); GWPROF("ExtWritePackages"); if (_owTable.empty()){WriteGridCells(Options); GWPROF("WriteGridCells");}}
  else {
  WriteMF6Files(Options); GWPROF("WriteMF6Files");
  WriteGridCells(Options); GWPROF("WriteGridCells");
  }
  ConnectEngine(); GWPROF("ConnectEngine");
  _active=true;
  if (_doSteady){RunSteadyState();}
  if (SZPhysical()){SZInitColumns();} //before Raven's initial storage is taken (CModel::CalculateInitialWaterStorage)
  if ((!_hotVz.empty() || !_hotVzc.empty()) && !(SZPhysical() && (_vzMode>0))){ //unsaturated-column water of the hotstart file a run without columns
    _vzHotUnusedN=0; _vzHotUnusedV=0.0; bool vol=!_hotVzt.empty();                //  cannot hold: dropped, warned about and reported
    for (map<long long,map<int,vector<double> > >::const_iterator it=_hotVz.begin();it!=_hotVz.end();it++){_vzHotUnusedN+=(int)it->second.size();}
    for (map<long long,map<int,map<int,vector<double> > > >::const_iterator it=_hotVzc.begin();it!=_hotVzc.end();it++){
      for (map<int,map<int,vector<double> > >::const_iterator jt=it->second.begin();jt!=it->second.end();jt++){_vzHotUnusedN+=(int)jt->second.size();}
    }
    for (map<long long,map<int,map<int,vector<double> > > >::const_iterator it=_hotVzt.begin();it!=_hotVzt.end();it++){
      for (map<int,map<int,vector<double> > >::const_iterator jt=it->second.begin();jt!=it->second.end();jt++){
        for (map<int,vector<double> >::const_iterator kt=jt->second.begin();kt!=jt->second.end();kt++){if (kt->second.size()>=2){_vzHotUnusedV+=kt->second[1];}}
      }
    }
    WriteWarning("CGroundwaterModel: the hotstart file holds "+to_string(_vzHotUnusedN)+" unsaturated columns (:GWVadoseColumns or :GWVadoseClassColumns), but this run has none ("+
                 string(SZPhysical()?":SoilZoneVadose NONE":"no :SoilZone PHYSICAL")+"): their water"+(vol?(", "+ToStr(_vzHotUnusedV)+" m3 above their reference profiles,"):string(""))+
                 " is dropped (reported in the run totals of GWModelSummary.txt). Restart with the columns of the run that wrote the file to keep it",false);
  }
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
  _hotRchStore.clear(); _hotResSeep.clear(); _hotResFrac.clear(); _hotRivCarry.clear(); _hotSzCols.clear(); _hotVz.clear(); _hotVzc.clear(); _hotVzt.clear(); _hotLag.clear(); _hotDebt.clear();
  _initMode="FROM_PROFILE"; _initValue=0.0;
}
//////////////////////////////////////////////////////////////////
/// \brief well-rate and boundary-head series of the .rvt file. An EnKF member's :ExtraRVTFilename is read after the model
///  was built (CalculateInitialWaterStorage): the series it gives replace the ones linked, are initialised and linked here
//
void CGroundwaterModel::SetWellFlowRateSeries(const long long id,CTimeSeries *pTS)
{
  map<long long,CTimeSeries*>::iterator it=_wellTS.find(id);
  if ((it!=_wellTS.end()) && (it->second!=pTS)){delete it->second;}
  _wellTS[id]=pTS;
  if (!_active){return;}
  pTS->Initialize(_jday0,_jyear0,_tsDur,_tsDt,false,_calendar);
  bool used=false; string ids;
  for (size_t w=0;w<_aWells.size();w++){if (_aWells[w].id==id){_aWells[w].pRate=pTS; used=true;} ids+=((w>0)?", ":"")+to_string(_aWells[w].id);}
  ExitGracefullyIf(!used,("CGroundwaterModel: :WellFlowRate "+to_string(id)+" in the .rvt file names no well of :Wells in the .rvg file ("+
    ((_aWells.size()>0)?("well IDs: "+ids):string("no wells"))+")").c_str(),BAD_DATA);
}
void CGroundwaterModel::SetBoundaryHeadSeries(const string &name,CTimeSeries *pTS)
{
  map<string,CTimeSeries*>::iterator it=_headTS.find(name);
  if ((it!=_headTS.end()) && (it->second!=pTS)){delete it->second;}
  _headTS[name]=pTS;
  if (!_active){return;}
  pTS->Initialize(_jday0,_jyear0,_tsDur,_tsDt,false,_calendar);
  bool used=false; string names;
  for (size_t b=0;b<_aBounds.size();b++){if (_aBounds[b].name==name){_aBounds[b].pHead=pTS; used=true;} names+=((b>0)?", ":"")+_aBounds[b].name;}
  ExitGracefullyIf(!used,("CGroundwaterModel: :BoundaryHead "+name+" in the .rvt file names no :GeneralHeadBoundary or :SpecifiedHeadBoundary of the .rvg file ("+
    ((_aBounds.size()>0)?("boundaries: "+names):string("no boundaries"))+")").c_str(),BAD_DATA);
}
//////////////////////////////////////////////////////////////////
/// \brief the state variable that receives HRU k's groundwater seepage (as Initialize sets _hruReturnSV), from the
///  inputs alone (the .rvc is read before Initialize); DOESNT_EXIST for an HRU that is not coupled
//
int CGroundwaterModel::SeepageTargetSV(const int k) const
{
  const CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
  string name=pHRU->GetAquiferProfileName();
  if ((name=="") || (name=="[NONE]") || (name=="NONE") || !pHRU->IsEnabled()){return DOESNT_EXIST;}
  int iSV=DOESNT_EXIST; bool toSoil=true;
  if (_seepTarget==""){
    int nsoil=_pModel->GetNumSoilLayers();
    iSV=(nsoil>0)?_pModel->GetStateVarIndex(SOIL,nsoil-1):_pModel->GetStateVarIndex(SURFACE_WATER);
  }
  else {
    int layer=0;
    sv_type typ=_pModel->GetStateVarInfo()->StringToSVType(_seepTarget,layer,false);
    if (typ==UNRECOGNIZED_SVTYPE){return DOESNT_EXIST;}
    iSV=_pModel->GetStateVarIndex(typ,layer); toSoil=(typ==SOIL);
  }
  if (toSoil){
    double th=0; for (int m=0;m<_pModel->GetNumSoilLayers();m++){th+=pHRU->GetSoilThickness(m);}
    if (th<=0.0){iSV=_pModel->GetStateVarIndex(SURFACE_WATER);}
  }
  return iSV;
}
//////////////////////////////////////////////////////////////////
/// \brief true when state variable i of HRU k is the store of one of the :SurfaceWaterExchange commands (its HRU group holds k,
///  and k is coupled); groundwater discharge fills that store with no capacity limit, as it does the seepage target
//
bool CGroundwaterModel::SWXStoreSV(const int k,const int i) const
{
  if (_aSWX.empty() || (SeepageTargetSV(k)==DOESNT_EXIST)){return false;}
  for (size_t x=0;x<_aSWX.size();x++){
    int lay=0; sv_type typ=_pModel->GetStateVarInfo()->StringToSVType(_aSWX[x].sv,lay,false);
    if ((typ==UNRECOGNIZED_SVTYPE) || (_pModel->GetStateVarIndex(typ,lay)!=i)){continue;}
    const CHRUGroup *pG=_pModel->GetHRUGroup(_aSWX[x].group);
    if ((pG!=NULL) && pG->IsInGroup(k)){return true;}
  }
  return false;
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
    ExitGracefullyIf(!(A.clapp_b>0) || !(A.psia>0),(n+"CLAPP_B and AIR_ENTRY_PRESSURE must be positive").c_str(),BAD_DATA);
    ExitGracefullyIf((_vzMode>0) && (A.porosity<A.Sy),(n+"with :SoilZoneVadose "+VZName()+", POROSITY (the saturated content of the unsaturated column) may not be smaller than SPEC_YIELD").c_str(),BAD_DATA);
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
    ExitGracefullyIf(_steadyInit,(string(msg)+":GWInitialization (the model's own stress periods decide)").c_str(),BAD_DATA);
    string gen=_genOnly; size_t pl;
    if (SZPhysical()){while ((pl=gen.find(" :LandSurface"))!=string::npos){gen.erase(pl,13);}} //(with the soil zone: its relief within the cells)
    ExitGracefullyIf(gen!="",(string(msg)+gen+" (these apply to models Raven builds: an existing model keeps its own grid, surfaces and solver settings)").c_str(),BAD_DATA);
    ExitGracefullyIf(SZPhysical() && (_szExtClass!="") && (GetAquiferClassIndex(_szExtClass)==DOESNT_EXIST),("CGroundwaterModel: :SoilZoneAquiferClass "+_szExtClass+" is not in :AquiferClasses").c_str(),BAD_DATA);
    if (SZPhysical() && (_vzMode>0) && (_szExtClass=="")){WriteWarning("CGroundwaterModel: :SoilZoneVadose "+VZName()+" with an existing model and no :SoilZoneAquiferClass: the unsaturated columns use CLAPP_B 4, AIR_ENTRY_PRESSURE 0.3 m and a porosity of SY + 0.1",false);}
    ExitGracefullyIf(_writeNC && (_owTable.size()>0),"CGroundwaterModel: :WriteNetCDFHeads needs the cell outlines (:MF6CRS or :MF6GridFile); an :OverlapWeights table has no coordinates",BAD_DATA);
    ExitGracefullyIf(_riverFile!="",(string(msg)+":RiverGeometry (rivers are the model's stream package, :StreamPackage)").c_str(),BAD_DATA);
    ExitGracefullyIf((_aObsWells.size()>0) && (_owTable.size()>0),"CGroundwaterModel: monitoring wells need :MF6CRS or :MF6GridFile (to locate them)",BAD_DATA);
    bool layersGiven=false;
    for (size_t i=0;i<_aProfiles.size();i++){
      if (_aProfiles[i].layers.size()>0){layersGiven=true;}
    }
    if (layersGiven){WriteWarning("CGroundwaterModel: with :MF6Simulation the layers of the aquifer profiles are ignored (the model's own layers are used); a profile can be given as: name, 0",false);}
  }
  else {
    for (size_t i=0;i<_aProfiles.size();i++){ //Raven builds the model from the profiles' layers
      ExitGracefullyIf(_aProfiles[i].layers.size()==0,("CGroundwaterModel: aquifer profile "+_aProfiles[i].name+" has no layers (layers may be left out only with :MF6Simulation)").c_str(),BAD_DATA);
      ExitGracefullyIf(_aProfiles[i].layers[0].type==GWL_AQUICLUDE,("CGroundwaterModel: aquifer profile "+_aProfiles[i].name+": top layer cannot be an AQUICLUDE").c_str(),BAD_DATA);
    }
    ExitGracefullyIf((_extGridFile!="") || _extCRS.IsSet() || (_owTable.size()>0) || (_takeOver.size()>0) || (_streamPkg!=""),
      "CGroundwaterModel: :MF6CRS, :MF6GridFile, :OverlapWeights, :RavenTakesOver and :StreamPackage need :MF6Simulation",BAD_DATA);
    ExitGracefullyIf(_extOnly!="",("CGroundwaterModel:"+_extOnly+" apply only to an existing model (:MF6Simulation); Raven builds this model, "
      "with its seepage drain and without MODFLOW packages of its own: remove them").c_str(),BAD_DATA);
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
  //a time series of the .rvt file that names no well or boundary would be dropped silently (a typo would change the forcing)
  for (map<long long,CTimeSeries*>::const_iterator it=_wellTS.begin();it!=_wellTS.end();it++){
    bool used=false; string ids;
    for (size_t w=0;w<_aWells.size();w++){if (_aWells[w].id==it->first){used=true;} ids+=((w>0)?", ":"")+to_string(_aWells[w].id);}
    ExitGracefullyIf(!used,("CGroundwaterModel: :WellFlowRate "+to_string(it->first)+" in the .rvt file names no well of :Wells in the .rvg file ("+
      ((_aWells.size()>0)?("well IDs: "+ids):string("no wells"))+")").c_str(),BAD_DATA);
  }
  for (map<string,CTimeSeries*>::const_iterator it=_headTS.begin();it!=_headTS.end();it++){
    bool used=false; string names;
    for (size_t b=0;b<_aBounds.size();b++){if (_aBounds[b].name==it->first){used=true;} names+=((b>0)?", ":"")+_aBounds[b].name;}
    ExitGracefullyIf(!used,("CGroundwaterModel: :BoundaryHead "+it->first+" in the .rvt file names no :GeneralHeadBoundary or :SpecifiedHeadBoundary of the .rvg file ("+
      ((_aBounds.size()>0)?("boundaries: "+names):string("no boundaries"))+")").c_str(),BAD_DATA);
  }
}
//////////////////////////////////////////////////////////////////
/// \brief warns once, after the .rvg file is read, for soil-zone options given where they do nothing (given: the options
///  written, as " :Name ..."); :SoilZone PHYSICAL is on
//
void CGroundwaterModel::WarnIgnoredSoilZoneOptions(const string &given,const bool noisy) const
{
  string vz="";
  if ((_vzMode==0) && (given.find(" :SoilZoneVadoseClasses")!=string::npos)){vz+=" :SoilZoneVadoseClasses";}
  if ((_vzMode==0) && (given.find(" :SoilZoneVadoseInitialFlux")!=string::npos)){vz+=" :SoilZoneVadoseInitialFlux";}
  if (vz!=""){WriteWarning("ParseGWFile:"+vz+" ignored: these apply only to unsaturated columns (:SoilZoneVadose RICHARDS, KINEMATIC or LAG); :SoilZoneVadose is NONE",noisy);}
  if (!IsExternal() && (given.find(" :SoilZoneAquiferClass")!=string::npos)){
    WriteWarning("ParseGWFile: :SoilZoneAquiferClass ignored: it applies only to an existing model (:MF6Simulation); the unsaturated columns of a model Raven builds use the aquifer classes of its profiles",noisy);
  }
}
void CGroundwaterModel::CheckProcesses()
{
  bool warned=false;
  if (SZPhysical()){_szRchByArea.assign(_pModel->GetNumHRUs(),0);}
  for (int j=0;j<_pModel->GetNumProcesses();j++)
  {
    CHydroProcessABC *pProc=_pModel->GetProcess(j);
    vector<CHydroProcessABC*> procs(1,pProc);                    //(expanded: process groups hold processes of their own)
    for (size_t pi=0;pi<procs.size();pi++){const CProcessGroup *pG=dynamic_cast<const CProcessGroup*>(procs[pi]);
      if (pG!=NULL){for (int s=0;s<pG->GetGroupSize();s++){procs.push_back(pG->GetSubProcess(s));}}}
    for (size_t pi=0;pi<procs.size();pi++){
      if (!SZPhysical()){continue;}
      //constant-rate percolation SOIL[M-1] -> GROUNDWATER takes the recharge from the columns by area (HRUs it applies to;
      //inside a process group the group's conditions: Raven does not apply a member's own conditions)
      const CmvPercolation *pPerc=dynamic_cast<const CmvPercolation*>(procs[pi]);
      if ((pPerc!=NULL) && (pPerc->GetPercolationType()==PERC_CONSTANT)){
        const int iBot=_pModel->GetStateVarIndex(SOIL,_szNLayers-1);
        const int *pf=procs[pi]->GetFromIndices(),*pt=procs[pi]->GetToIndices();
        for (int q=0;q<procs[pi]->GetNumConnections();q++){
          if ((pf[q]!=iBot) || (pt[q]!=_iGW)){continue;}
          for (size_t i=0;i<_aqHRUs.size();i++){const CHydroUnit *pH=_pModel->GetHydroUnit(_aqHRUs[i]);
            if (pProc->ShouldApply(pH)){_szRchByArea[_aqHRUs[i]]=1;}}
        }
      }
      if (procs[pi]->GetNumLatConnections()<=0){continue;} //lateral processes may not reach the pass-through GROUNDWATER either
      const CLateralExchangeProcessABC *pLat=dynamic_cast<const CLateralExchangeProcessABC*>(procs[pi]);
      if (pLat!=NULL){
        const int *lf=pLat->GetLateralFromIndices(),*lt=pLat->GetLateralToIndices(),*kf=pLat->GetFromHRUIndices(),*kt=pLat->GetToHRUIndices();
        for (int q=0;q<procs[pi]->GetNumLatConnections();q++){
          bool bad=((lt[q]==_iGW) && (_hruProfile[kt[q]]!=DOESNT_EXIST)) || ((lf[q]==_iGW) && (_hruProfile[kf[q]]!=DOESNT_EXIST));
          ExitGracefullyIf(bad,("CGroundwaterModel: with :SoilZone PHYSICAL, lateral process #"+to_string(j+1)+" moves water into or out of GROUNDWATER of a coupled HRU; "
            "recharge must come from the bottom soil-zone layer by a vertical process").c_str(),BAD_DATA);
        }
      }
    }
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
      if (SZPhysical() && (iTo[q]==_iGW) && (iFrom[q]!=_iGW)){ //recharge must leave the soil zone from its bottom layer
        int iBot=_pModel->GetStateVarIndex(SOIL,_szNLayers-1);
        if (iFrom[q]!=iBot){
          for (size_t i=0;i<_aqHRUs.size();i++){
            if (_szOn[_aqHRUs[i]] && pProc->ShouldApply(_pModel->GetHydroUnit(_aqHRUs[i]))){
              ExitGracefully(("CGroundwaterModel: with :SoilZone PHYSICAL, recharge must come from the bottom soil-zone layer SOIL["+to_string(_szNLayers-1)+
                "]; process #"+to_string(j+1)+" sends water to GROUNDWATER from another store in HRU "+to_string(_pModel->GetHydroUnit(_aqHRUs[i])->GetHRUID())).c_str(),BAD_DATA);
            }
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
  double rb[4];
  if (!GWGeom::ReadGeoJSONPolygons(_hruGeomFile,_hruIDField,shapes,err,InCRS(),rb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
  CheckInputLonLat(_hruGeomFile,rb,"GeoJSON file");

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
  if (_gridType!=GWGRID_FILE){ //a cell size far too small for the HRUs' extent (given in km, or split by a thousands separator) would exhaust the memory
    const double maxCells=1.0e7;
    double nc=ceil((maxx-minx)/cs)+1.0,nr=ceil((maxy-miny)/cs)+1.0;
    if (nc*nr>maxCells){
      ostringstream S; S<<"CGroundwaterModel: :GridCellSize "<<cs<<" m gives about "<<setprecision(3)<<nc*nr<<" cells per layer over the HRUs' extent ("
        <<setprecision(4)<<(maxx-minx)/1000.0<<" km by "<<(maxy-miny)/1000.0<<" km), more than "<<maxCells<<": give the cell size in metres as a plain number "
        <<"(for example 2000, not 2,000 or 2km)";
      ExitGracefully(S.str().c_str(),BAD_DATA);
    }
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
  if ((_gridType==GWGRID_REGULAR) && !ref.empty()){ //columns x rows of the refined regular grid (before its 2:1 smoothing), before it is built
    double ncr[2]={0,0};
    for (int a=0;a<2;a++){
      double a0=floor(((a==0)?minx:miny)/cs)*cs; int n=max(1,(int)ceil((((a==0)?maxx:maxy)-a0)/cs));
      for (int i=0;i<n;i++){
        double lo=a0+i*cs,hi=lo+cs,target=cs;
        for (size_t f=0;f<ref.size();f++){
          double mn=1e300,mx=-1e300;
          for (size_t k=0;k<ref[f].geom.size();k++){double v=(a==0)?ref[f].geom[k].x:ref[f].geom[k].y; mn=min(mn,v); mx=max(mx,v);}
          if ((mx+ref[f].size>=lo) && (mn-ref[f].size<=hi)){target=min(target,ref[f].size);}
        }
        ncr[a]+=max(1.0,ceil(cs/target-1e-9));
      }
    }
    if (ncr[0]*ncr[1]>1.0e7){
      ostringstream S; S<<"CGroundwaterModel: the :GridRefinement sizes give a regular grid of about "<<setprecision(3)<<ncr[0]*ncr[1]<<" cells per layer ("
        <<setprecision(9)<<ncr[0]<<" columns by "<<ncr[1]<<" rows), more than 1e+07: give larger refinement sizes in metres as plain numbers, "
        <<"or use :GridType QUADTREE or HRU_MESH, which refine only around the features";
      ExitGracefully(S.str().c_str(),BAD_DATA);
    }
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
    double gb[4];
    if (!GWGeom::ReadGeoJSONPolygons(_gridFile,_gridFileID,cellshapes,err,InCRS(),gb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
    CheckInputLonLat(_gridFile,gb,"GeoJSON file");
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
  dem.SetCRS(InCRS());
  if (_demFile!=""){ReadRaster(dem,_demFile,Options);}
  _land.assign(_ncpl,0.0);
  vector<int> act; int nCov=0; //(active columns, and how many of them the DEM gives a value for)
  for (int ic=0;ic<_ncpl;ic++)
  {
    if (cover[ic]<_minCoverage*_grid.area[ic]){_colProfile[ic]=DOESNT_EXIST;}
    if (_colProfile[ic]==DOESNT_EXIST){continue;}
    bool ok=false;
    if (dem.IsLoaded()){_land[ic]=MeanRaster(dem,ic,ok); act.push_back(ic); if (ok){nCov++;}}
    if (!ok){_land[ic]=elevW[ic]/cover[ic];}
  }
  if (dem.IsLoaded()){CheckRasterCover(dem,_demFile,":LandSurface",act,nCov,"the area-weighted ELEVATION of their HRUs",Options);}
}

//////////////////////////////////////////////////////////////////
/// \brief builds model top, layer bottoms, bedrock clipping, pass-through and inactive cells, initial heads
//
void CGroundwaterModel::BuildVertical(const optStruct &Options)
{
  CGWRaster bed;
  bed.SetCRS(InCRS());
  if (_bedrockFile!=""){ReadRaster(bed,_bedrockFile,Options);}
  vector<int> act; int nCov=0; //(active columns, and how many of them the bedrock raster gives a value for)
  _nlay=1;
  for (size_t i=0;i<_aProfiles.size();i++){_nlay=max(_nlay,(int)_aProfiles[i].layers.size());}

  _top.assign(_ncpl,0.0); _bedrock.assign(_ncpl,0.0);
  _botm.assign(_nlay*_ncpl,0.0); _strt.assign(_nlay*_ncpl,0.0); _idomain.assign(_nlay*_ncpl,0);
  int nThin=0;
  vector<double> colD; //:SoilZone PHYSICAL: overlap-weighted soil-zone thickness of the column (0 for HRUs without soil)
  if (SZPhysical()){
    colD.assign(_ncpl,0.0); vector<double> colA(_ncpl,0.0);
    for (size_t i=0;i<_aqHRUs.size();i++){int k=_aqHRUs[i];
      for (map<int,double>::iterator it=_hruCells[k].begin();it!=_hruCells[k].end();it++){colA[it->first]+=it->second; colD[it->first]+=it->second*_szD[k];}}
    for (int ic=0;ic<_ncpl;ic++){colD[ic]=(colA[ic]>0)?colD[ic]/colA[ic]:0.0;}
  }
  for (int ic=0;ic<_ncpl;ic++)
  {
    int ip=_colProfile[ic];
    if (ip!=DOESNT_EXIST)
    {
      const gw_profile &P=_aProfiles[ip];
      double land=_land[ic];
      double top =SZPhysical()?land-colD[ic]:land-P.soil_zone_depth;
      double bedrock=land-P.default_bedrock_depth;
      bool ok=false;
      if (bed.IsLoaded()){double b=MeanRaster(bed,ic,ok); act.push_back(ic); if (ok){bedrock=b; nCov++;}}
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
  if (bed.IsLoaded()){CheckRasterCover(bed,_bedrockFile,":BedrockSurface",act,nCov,"a bedrock at their profile's DEFAULT_BEDROCK_DEPTH below the land surface",Options);}
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
  double rb[4];
  if (!GWGeom::ReadGeoJSONLines(_riverFile,_riverIDField,lines,ids,err,InCRS(),rb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
  CheckInputLonLat(_riverFile,rb,"GeoJSON file");
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
  map<long long,vector<gw_ring> > polys; double rb[4];
  if (GWGeom::ReadGeoJSONPolygons(file,"",polys,err,InCRS(),rb))
  {
    CheckInputLonLat(file,rb,"GeoJSON file");
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
    ExitGracefullyIf(err.find("no Polygon/MultiPolygon")==string::npos,("CGroundwaterModel: "+err).c_str(),BAD_DATA); //(a polygon file refused by its coordinates)
    vector<gw_ring> lines; vector<long long> ids; err="";
    if (!GWGeom::ReadGeoJSONLines(file,"",lines,ids,err,InCRS(),rb)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
    CheckInputLonLat(file,rb,"GeoJSON file");
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
    if (W.pRate==NULL){WriteWarning("CGroundwaterModel: well "+W.name+" has no :WellFlowRate time series in the .rvt file; rate set to zero",false);}
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
    //the highest and lowest heads the boundary can have: HEAD and every value of its :BoundaryHead series (the whole series,
    //not only the run's steps, so a restarted run gets the same cells; HEAD applies where the series has no value). Each
    //cell whose bottom lies below the highest head is placed; on a step whose head is at or below a cell's bottom that
    //cell is switched off (SetBoundaryArrays). LAYERS TOP: the top layer the head can reach and the layers below it down
    //to the first one the lowest head reaches; each step holds the highest of them whose bottom lies below the step's head
    double hp=B.head,hm=B.head;
    if (B.pHead!=NULL){
      for (int s=0;s<B.pHead->GetNumValues();s++){double v=B.pHead->GetValue(s); if ((v==RAV_BLANK_DATA) || !std::isfinite(v)){continue;} hp=max(hp,v); hm=min(hm,v);}
    }
    size_t nb=_chdNode.size()+_ghbNode.size(); double zmax=-1e300;
    for (size_t i=0;i<cols.size();i++){
      for (int l=0;l<_nlay;l++){ //MODFLOW requires boundary heads above the cell bottom: skip layers the head never rises above
        int n=UserNode(l,cols[i]); if (_idomain[n]<=0){continue;}
        zmax=max(zmax,_botm[n]); if (_botm[n]>=hp){continue;}
        if (B.isCHD){_chdNode.push_back(n);_chdB.push_back((int)b);} else {_ghbNode.push_back(n);_ghbB.push_back((int)b);}
        if (!B.allLayers && (_botm[n]<hm)){break;}
      }
    }
    if ((cols.size()>0) && (_chdNode.size()+_ghbNode.size()==nb)){
      WriteWarning("CGroundwaterModel: boundary "+B.name+" has no cells: its head ("+ToStr(hp)+" m"+((B.pHead!=NULL)?", the highest of HEAD and its :BoundaryHead series":", HEAD")+
                   ") lies below the bottom of every cell it selects (highest "+ToStr(zmax)+" m); MODFLOW needs a boundary head above the cell bottom",false);
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
  F<<"BEGIN OPTIONS\n  NEWTON\n"<<(_writeCBC?"  SAVE_FLOWS   # :WriteGWBudgetFile: every package's flows to gwf.cbc\n":"")<<"END OPTIONS\n\nBEGIN PACKAGES\n";
  F<<"  "<<((_disType==0)?"DIS6   gwf.dis ":(_disType==1)?"DISV6  gwf.disv":"DISU6  gwf.disu")<<"  DIS\n  NPF6  gwf.npf  NPF\n  STO6  gwf.sto  STO\n  IC6   gwf.ic   IC\n";
  F<<"  RCH6  gwf.rch  RCH_RAVEN\n  DRN6  gwf.drn  DRN_RAVEN\n";
  if (_rivNode.size()>0){F<<"  RIV6  gwf.riv  RIV_RAVEN\n";}
  if (_evtCol.size()>0) {F<<"  EVT6  gwf.evt  EVT_RAVEN\n";}
  if (_welNode.size()>0){F<<"  WEL6  gwf.wel  WEL_RAVEN\n";}
  if (_ghbNode.size()>0){F<<"  GHB6  gwf.ghb  GHB_RAVEN\n";}
  if (_chdNode.size()>0){F<<"  CHD6  gwf.chd  CHD_RAVEN\n";}
  if (_swxCol.size()>0) {F<<"  RIV6  gwf.swx  SWX_RAVEN\n";}
  if (_resNode.size()>0){F<<"  WEL6  gwf.res  RES_RAVEN\n";}
  if (_szxLink.size()>0){F<<"  GHB6  gwf.szx  SZX_RAVEN\n  WEL6  gwf.szq  SZQ_RAVEN\n  DRN6  gwf.szd  SZD_RAVEN\n";}
  if (_szwNode.size()>0){F<<"  WEL6  gwf.szw  SZW_RAVEN\n";}
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
  vector<double> K(_nlay*_ncpl,1.0),K33(_nlay*_ncpl,1.0),Ss(_nlay*_ncpl,1e-5),Sy(_nlay*_ncpl,0.1),Por(_nlay*_ncpl,0.3);
  for (int ic=0;ic<_ncpl;ic++){
    int ip=_colProfile[ic]; if (ip==DOESNT_EXIST){continue;}
    const gw_profile &P=_aProfiles[ip];
    for (int l=0;l<(int)P.layers.size();l++){
      int n=UserNode(l,ic);
      const gw_aquifer_class &C=_aClasses[P.layers[l].iclass];
      K[n]=C.Kh; K33[n]=C.Kv; Ss[n]=C.Ss; Sy[n]=C.Sy; Por[n]=C.porosity;
      bool conf=(P.layers[l].type==GWL_CONFINED_AQUIFER);
      icelltype[n]=conf?0:1; iconvert[n]=conf?0:1;
    }
  }
  F.open((d+"gwf.npf").c_str());
  F<<"BEGIN OPTIONS\n  SAVE_FLOWS\n"<<(UseXT3D()?((_flowCorr=="XT3D")?"  XT3D\n":"  XT3D RHS\n"):"")<<(((_disType==2) || _writeCBC)?"  SAVE_SPECIFIC_DISCHARGE\n":"")
   <<(_writeCBC?"  SAVE_SATURATION\n":"")<<"END OPTIONS\n\nBEGIN GRIDDATA\n"; //(:WriteGWBudgetFile: specific discharge and saturation for particle tracking)
  WriteArray(F,"ICELLTYPE",icelltype,aC,aL,"1"); WriteArray(F,"K",K,aC,aL); WriteArray(F,"K33",K33,aC,aL);
  F<<"END GRIDDATA\n"; F.close();

  F.open((d+"gwf.sto").c_str());
  F<<"BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN GRIDDATA\n";
  WriteArray(F,"ICONVERT",iconvert,aC,aL,"1"); WriteArray(F,"SS",Ss,aC,aL); WriteArray(F,"SY",Sy,aC,aL);
  if (_doSteady){F<<"END GRIDDATA\n\nBEGIN PERIOD 1\n  STEADY-STATE\nEND PERIOD\n\nBEGIN PERIOD 2\n  TRANSIENT\nEND PERIOD\n";}
  else            {F<<"END GRIDDATA\n\nBEGIN PERIOD 1\n  TRANSIENT\nEND PERIOD\n";}
  F.close();

  if (_writeCBC){ //:WriteGWBudgetFile: cell porosity for offline particle tracking (the GRIDDATA of a PRT MIP6 file)
    F.open((d+"gwf.porosity").c_str());
    F<<"# porosity of every cell (POROSITY of the aquifer classes) for particle tracking: usable as the MIP6 file of a PRT model\nBEGIN GRIDDATA\n";
    WriteArray(F,"POROSITY",Por,aC,aL); F<<"END GRIDDATA\n"; F.close();
  }
  F.open((d+"gwf.ic").c_str());
  F<<"BEGIN GRIDDATA\n"; WriteArray(F,"STRT",_strt,aC,aL,"1.0",SZPhysical()?17:10); F<<"END GRIDDATA\n"; F.close(); //(full precision with the soil zone, as its hotstart)

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
    //:SoilZone PHYSICAL: where soil-zone HRUs cover the column the exchange is SZX_RAVEN; the drain serves the rest
    if (SZPhysical() && (_colCover[ic]>0)){cond*=_colCoverNS[ic]/_colCover[ic];}
    //MODFLOW refuses a drain whose smoothing interval reaches below the cell bottom: limit it to the top layer (less a
    //margin: written to 10 digits, top minus the full thickness could round below the bottom)
    double dsm=min(_aProfiles[_colProfile[ic]].seep_smooth_depth,max(_top[ic]-_botm[UserNode(0,ic)]-1e-8*max(1.0,fabs(_top[ic])),0.0));
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
      //surface: the land (also with :SoilZone PHYSICAL: the extinction depth counts from the ground, above the soil zone)
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
      F<<"  "<<CellID(l,ic)<<" "<<((B.head>_botm[n])?B.head:_botm[n]+0.01)<<" "<<B.cond<<endl;}  //(a cell placed by a series' higher heads: Raven sets each step's head)
    F<<"END PERIOD\n"; F.close();
  }
  if (_chdNode.size()>0){
    F.open((d+"gwf.chd").c_str());
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_chdNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_chdNode.size();e++){int n=_chdNode[e],l=n/_ncpl,ic=n%_ncpl; const double h=_aBounds[_chdB[e]].head; F<<"  "<<CellID(l,ic)<<" "<<((h>_botm[n])?h:_botm[n]+0.01)<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_swxCol.size()>0){
    F.open((d+"gwf.swx").c_str());
    F<<"# surface-water stores (wetlands, lakes) exchanging with the aquifer: river-type boundary with its bed at the\n# land surface; stage = land surface + store depth and conductance are set every step\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_swxCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_swxCol.size();e++){int ic=_swxCol[e]; F<<"  "<<CellID(0,ic)<<" "<<SWXDatum(e)<<" 0.0 "<<SWXDatum(e)<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_resNode.size()>0){
    F.open((d+"gwf.res").c_str());
    F<<"# seepage of Raven reservoirs into the aquifer under their lake HRUs (negative: groundwater feeding the reservoir)\n";
    F<<"BEGIN OPTIONS\n  AUTO_FLOW_REDUCE 0.1\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_resNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
    for (size_t e=0;e<_resNode.size();e++){int n=_resNode[e],ic=n%_ncpl; F<<"  "<<CellID(0,ic)<<" 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
  }
  if (_szxLink.size()>0){
    F.open((d+"gwf.szx").c_str());
    F<<"# soil zone <-> aquifer (:SoilZone PHYSICAL): one entry per HRU-cell link, the exchange of the link's soil-zone\n";
    F<<"# column with the cell as a linear (Robin) boundary; head and conductance set by Raven every iteration\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_szxLink.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_szxLink.size();e++){int ic=_aLinks[_szxLink[e]].ic; F<<"  "<<CellID(0,ic)<<" "<<max(_szxZb[e],_botm[UserNode(0,ic)]+1e-3)<<" 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
    F.open((d+"gwf.szq").c_str());
    F<<"# soil zone <-> aquifer, second part: where the boundary head of SZX_RAVEN would lie below the cell bottom (which\n";
    F<<"# MODFLOW does not accept) that head is raised to the bottom and this rate carries the difference; set by Raven\n";
    F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_szxLink.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
    for (size_t e=0;e<_szxLink.size();e++){int ic=_aLinks[_szxLink[e]].ic; F<<"  "<<CellID(0,ic)<<" 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
    F.open((d+"gwf.szd").c_str());
    F<<"# soil zone <-> aquifer, third part: a boundary MODFLOW keeps (explicit mode, frozen interface) on links where water\n";
    F<<"# rises into the soil, as a drain: it stops instead of reversing when the head falls (onset smoothed over DDRN)\n";
    F<<"BEGIN OPTIONS\n  AUXILIARY DDRN\n  AUXDEPTHNAME DDRN\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_szxLink.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(10);
    for (size_t e=0;e<_szxLink.size();e++){int ic=_aLinks[_szxLink[e]].ic; F<<"  "<<CellID(0,ic)<<" "<<_top[ic]<<" 0.0 0.0"<<endl;}
    F<<"END PERIOD\n"; F.close();
    if (_szwNode.size()>0){
      F.open((d+"gwf.szw").c_str());
      F<<"# soil zone <-> aquifer, fourth part: upward flows (capillary rise, the unsaturated column's upflow) out of the cell\n";
      F<<"# holding the water table when the top cell has dried, which MODFLOW could not supply; rates set by Raven\n";
      F<<"BEGIN OPTIONS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_szwNode.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
      for (size_t j=0;j<_szwNode.size();j++){int n=_szwNode[j]; F<<"  "<<CellID(n/_ncpl,n%_ncpl)<<" 0.0"<<endl;}
      F<<"END PERIOD\n"; F.close();
    }
  }
  F.open((d+"gwf.oc").c_str());
  if (_writeCBC){ //:WriteGWBudgetFile: cell-by-cell flows at the same steps as the heads
    F<<"BEGIN OPTIONS\n  BUDGET FILEOUT gwf.cbc\n  HEAD FILEOUT gwf.hds\nEND OPTIONS\n\nBEGIN PERIOD 1\n  SAVE HEAD FREQUENCY "<<max(1,_headSaveFreq)<<"\n  SAVE BUDGET FREQUENCY "<<max(1,_headSaveFreq)<<"\n  PRINT BUDGET LAST\nEND PERIOD\n";
  }
  else {F<<"BEGIN OPTIONS\n  HEAD FILEOUT gwf.hds\nEND OPTIONS\n\nBEGIN PERIOD 1\n  SAVE HEAD FREQUENCY "<<max(1,_headSaveFreq)<<"\n  PRINT BUDGET LAST\nEND PERIOD\n";}
  F.close();
}

//////////////////////////////////////////////////////////////////
/// \brief MODFLOW 6 library: :MF6Library, else the RAVEN_MF6_LIB environment variable, else the first library found next
///  to the Raven executable or in the lib/mf6 folder of the source checkout, else the
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
  s+=" then "+CMF6Engine::DefaultLibraryName()+" on the system library path. Download it with get_mf6.py of the RavenHydroFramework-extra "
     "repository (https://github.com/rarabzad/RavenHydroFramework-extra), or from https://github.com/MODFLOW-ORG/modflow6/releases "
     "(the library is in the bin folder of the release zip), and put it next to the Raven executable";
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
  else if (!(IsExternal() && _extNoSto)){WriteWarning("CGroundwaterModel: MODFLOW 6 storage flows not accessible; GWBudget.csv storage column will be inferred",false);}
  //(an existing model without STO6 is steady state: its storage is 0, nothing to infer)
  //convertible cells (KeepColumnsWet): a generated model's from the layer types it was written with (an existing model's
  //from MODFLOW's memory, ExtConnect)
  if (!IsExternal()){
    _cellConv.assign(_nlay*_ncpl,1);
    for (int ic=0;ic<_ncpl;ic++){
      if (_colProfile[ic]==DOESNT_EXIST){continue;}
      const gw_profile &P=_aProfiles[_colProfile[ic]];
      for (int l=0;(l<_nlay) && (l<(int)P.layers.size());l++){if (P.layers[l].type==GWL_CONFINED_AQUIFER){_cellConv[UserNode(l,ic)]=0;}}
    }
  }
  //an outer-iteration limit below 20 leaves many steps unconverged with the soil-zone interface, which is updated within
  //MODFLOW's outer iterations; written once MODFLOW's limit is known
  if (SZPhysical() && (_szxLink.size()>0) && (_maxiter<20)){
    WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL with "+(IsExternal()?"OUTER_MAXIMUM "+to_string(_maxiter)+" in the model's IMS file (MODFLOW's MXITER)":
                 ":SolverMaxOuterIterations "+to_string(_maxiter))+": the soil-zone interface is updated within MODFLOW's outer iterations; a limit below 20 leaves many steps unconverged",false);
  }
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
    _chdNb=_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,"CHD_RAVEN","NBOUND")); _chdNl=_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,"CHD_RAVEN","NODELIST"));
    {string a2=_mf6.GetVarAddress(_mname,"","IBOUND"); _chdIb=_mf6.HasVar(a2)?_mf6.GetIntPtr(a2):NULL;}
  }
  if (_swxCol.size()>0){
    string a=_mf6.GetVarAddress(_mname,"SWX_RAVEN","STAGE");
    if (_mf6.HasVar(a)){_swxNcol=0; _swxHead=_mf6.GetDoublePtr(a); _swxCond=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","COND"));
                        _swxRbot=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","RBOT"));}
    else {_swxNcol=3; _swxHead=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","BOUND")); _swxCond=(_swxHead==NULL)?NULL:_swxHead+1; _swxRbot=(_swxHead==NULL)?NULL:_swxHead+2;}
    _swxHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","HCOF")); _swxRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SWX_RAVEN","RHS"));
    ExitGracefullyIf((_swxHead==NULL) || (_swxCond==NULL) || (_swxRbot==NULL) || (_swxHcof==NULL) || (_swxRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 surface-water exchange arrays",RUNTIME_ERR);
  }
  if (_szxLink.size()>0){ //MODFLOW >= 6.5 re-reads BHEAD/COND every outer iteration (tests/soilzone/spike_xmi.py)
    string a=_mf6.GetVarAddress(_mname,"SZX_RAVEN","BHEAD");
    if (_mf6.HasVar(a)){_szxNcol=0; _szxHead=_mf6.GetDoublePtr(a); _szxCond=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZX_RAVEN","COND"));}
    else {_szxNcol=2; _szxHead=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZX_RAVEN","BOUND")); _szxCond=(_szxHead==NULL)?NULL:_szxHead+1;}
    _szxHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZX_RAVEN","HCOF")); _szxRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZX_RAVEN","RHS"));
    ExitGracefullyIf((_szxHead==NULL) || (_szxCond==NULL) || (_szxHcof==NULL) || (_szxRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 soil-zone exchange arrays",RUNTIME_ERR);
    a=_mf6.GetVarAddress(_mname,"SZQ_RAVEN","Q");
    if (_mf6.HasVar(a)){_szqNcol=0; _szqQ=_mf6.GetDoublePtr(a);} else {_szqNcol=1; _szqQ=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZQ_RAVEN","BOUND"));}
    _szqHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZQ_RAVEN","HCOF")); _szqRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZQ_RAVEN","RHS"));
    ExitGracefullyIf((_szqQ==NULL) || (_szqHcof==NULL) || (_szqRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 soil-zone exchange arrays (SZQ_RAVEN)",RUNTIME_ERR);
    if (_szwNode.size()>0){
      a=_mf6.GetVarAddress(_mname,"SZW_RAVEN","Q");
      if (_mf6.HasVar(a)){_szwNcol=0; _szwQ=_mf6.GetDoublePtr(a);} else {_szwNcol=1; _szwQ=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZW_RAVEN","BOUND"));}
      ExitGracefullyIf(_szwQ==NULL,"CGroundwaterModel: cannot access MODFLOW 6 soil-zone exchange arrays (SZW_RAVEN)",RUNTIME_ERR);
    }
    a=_mf6.GetVarAddress(_mname,"SZD_RAVEN","ELEV");
    if (_mf6.HasVar(a)){_szdNcol=0; _szdElev=_mf6.GetDoublePtr(a); _szdCond=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZD_RAVEN","COND"));}
    else {_szdNcol=2; _szdElev=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZD_RAVEN","BOUND")); _szdCond=(_szdElev==NULL)?NULL:_szdElev+1;}
    //(one auxiliary variable, DDRN; MODFLOW >= 6.4.3 reads it from the input context, AUXVAR_IDM, not from AUXVAR)
    a=_mf6.GetVarAddress(_mname,"SZD_RAVEN","AUXVAR_IDM");
    _szdAux=_mf6.GetDoublePtr(_mf6.HasVar(a)?a:_mf6.GetVarAddress(_mname,"SZD_RAVEN","AUXVAR"));
    _szdHcof=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZD_RAVEN","HCOF")); _szdRhs=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,"SZD_RAVEN","RHS"));
    ExitGracefullyIf((_szdElev==NULL) || (_szdCond==NULL) || (_szdAux==NULL) || (_szdHcof==NULL) || (_szdRhs==NULL),"CGroundwaterModel: cannot access MODFLOW 6 soil-zone exchange arrays (SZD_RAVEN)",RUNTIME_ERR);
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
  _stepCount++; _stepEnd=EndDate(tt.model_time+tstep); //(per-step warnings name the end of the step, as GWBudget.csv)
  ExitGracefullyIf(_stepCount>_nstp,"CGroundwaterModel: Raven is running more time steps than the MODFLOW 6 simulation holds (check :Duration and :TimeStep)",RUNTIME_ERR);
  if (IsExternal()){ //an existing model's specified-head cells are known only once MODFLOW has read the step's period data
    if (!_mf6.PrepareTimeStep(tstep)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
    RefreshPointers();
    ExtUpdateCHDTop();
  }
  if (SZPhysical()){ //soil-zone water after Raven's processes this step: start of every local solve
    //Raven's changes are shared among the columns (gwsz::Distribute); the recharge (end-of-step GROUNDWATER, which only
    //SOIL[M-1] fills) leaves the columns first, as the process takes it (by area for PERC_CONSTANT, else from free water
    //first), and each column's part is its share of the recharge
    int iSW=_pModel->GetStateVarIndex(SURFACE_WATER);
    vector<double> lossR;   //(_szMovedV is reset when the budget line is written: a hotstart sync counts in the first step)
    for (size_t i=0;i<_aqHRUs.size();i++){int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
      for (size_t j=0;j<_szSV[k].size();j++){_szPhiProc[_szPhiOff[k]+j]=aPhinew[k][_szSV[k][j]];}
      double ovf=((SZClassCount(k)>1) && (k<(int)_szClsDone.size()) && _szClsDone[k])?
                 SZSyncClasses(k,&_szPhiProc[_szPhiOff[k]],max(aPhinew[k][_iGW],0.0),lossR):  //Raven's processes ran per column class
                 SZSyncColumns(k,&_szPhiProc[_szPhiOff[k]],max(aPhinew[k][_iGW],0.0),&lossR);
      if (_szTiming==1){SZLagReverse(k);}                                          //LAGGED: the source back out of the columns
      for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){int e=_szxOfHRU[k][jj],c=_szxCol[e]; _szxFw[e]=lossR[c-_szColOff[k]]/MM_PER_METER*_szColA[c];}
      for (size_t jj=0;jj<_szcOfHRU[k].size();jj++){int e=_szcOfHRU[k][jj],c=_szcCol[e]; _szcFw[e]=lossR[c-_szColOff[k]]/MM_PER_METER*_szColA[c];}
      //Raven's soil-zone layers become the columns' mean; gains full columns could not take overflow to SURFACE_WATER
      int nc=SZNCols(k),o=_szColOff[k]; double sw=0.0; for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
      for (size_t j=0;j<_szSV[k].size();j++){
        double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*_szNLayers+j];}
        if (sw>0.0){aPhinew[k][_szSV[k][j]]=_szPhiProc[_szPhiOff[k]+j]=v/sw;}
      }
      if (ovf>0.0){aPhinew[k][iSW]+=ovf;}
      _szSyncOvf[k]=ovf/MM_PER_METER*_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
    }
    _szColInit=true;
    _szTstep=tstep; _szResid=0.0; _szIter=0; _szConv=true; _szExtra=0; _szWhy=0; _szWhyE=-1;
    _szxDropped.assign(_szxLink.size(),0); _szDropThis=false;
  }
  if (_xStart.empty() && (_X!=NULL)){const int nx=_mf6.GetVarSize(_mf6.GetVarAddress(_mname,"","X")); if (nx>0){_xStart.assign(_X,_X+nx);}} //(heads at the start of the run)
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
    if (SZPhysical()){_szRchV[k]=V;} //distributed with rejection by SZRecharge
    if ((V<0.0) && !_negWarned && !(SZPhysical() && _szOn[k])){ //(a negative recharge without a soil zone to limit it: plain recharge)
      _negWarned=true;
      bool conf=false; for (size_t j=0;j<_hruLinks[k].size();j++){const int n=TopNode(_aLinks[_hruLinks[k][j]].ic); if ((n<(int)_cellConv.size()) && !_cellConv[n]){conf=true;}}
      WriteWarning("CGroundwaterModel: HRU "+to_string(_pModel->GetHydroUnit(k)->GetHRUID())+" takes water out of GROUNDWATER (a negative recharge, e.g. PERC_GR4JEXCH) in the step ending "+_stepEnd+
        "; without a soil zone (:SoilZone PHYSICAL on an HRU with soil) it is not limited by what the aquifer cells can give: "+
        (_mf6Newton?string("MODFLOW takes it out of the top cell whatever the cell holds (a cell drained by it can make the solve diverge)"):
                    string("without NEWTON MODFLOW takes it until the cells dry and then applies none of it (booked as recharge not applied)")+
                    (conf?"; its confined cells (ICELLTYPE 0) never dry and take all of it":"")),false);
    }
    if (_hruRchArea[k]<=0){ //HRU entirely on specified-head cells (with a soil zone: its columns there, SZRecharge)
      if (!(SZPhysical() && _szOn[k])){_QrchCHD+=V/tstep;}
      continue;
    }
    for (size_t j=0;j<_hruLinks[k].size();j++){
      const gw_link &L=_aLinks[_hruLinks[k][j]];
      if (_colCHDTop[L.ic]){continue;}
      _Qcol[L.ic]+=V/tstep*L.area/_hruRchArea[k];
    }
  }

  _QrchCHD0=_QrchCHD; _szRchCHD=0.0;
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
    for (size_t i=0;i<_evtCol.size();i++){int ic=_evtCol[i]; _evtRate[i*es]=_work[ic]/MF6CellArea(ic)*_rf;} //(model units: an existing model's own cell area)
    if (SZPhysical()){_evtBase.assign(_evtCol.size(),0.0); for (size_t i=0;i<_evtCol.size();i++){_evtBase[i]=_evtRate[i*es];}}
  }
  SetBoundaryArrays(tt.model_time,tstep); _bndOffN+=_bndOffNow;
  //reservoirs: seepage Raven computed during the previous step's routing enters the aquifer now (one-step lag)
  if (_stepCount==2){ //(the stage Raven's first routing set: with no initial flow in the .rvc file the stage is 0 until
                      // then, and nothing before that routing depends on it)
    for (size_t r=0;r<_resP.size();r++){
      CReservoir *pRes=_pModel->GetSubBasin(_resP[r])->GetReservoir(); const double zbot=(r<_resZbot.size())?_resZbot[r]:1e300;
      ExitGracefullyIf((zbot<1e299) && (pRes->GetResStage()<zbot),("CGroundwaterModel: reservoir in subbasin "+to_string(_pModel->GetSubBasin(_resP[r])->GetID())+
        " has stage "+to_string(pRes->GetResStage())+" m, below the bottom of the aquifer under its lake ("+to_string(zbot)+
        " m): coupled reservoirs need absolute stages on the aquifer's datum (e.g. :AbsoluteCrestHeight)").c_str(),BAD_DATA);
    }
  }
  _Qres_spec=0.0;
  //cells that cannot take the seepage this step (an existing model's specified heads, and inactive cells of a model without
  //NEWTON: MODFLOW's WEL books no flow there) get none; the reservoir's other cells take it, by their overlap weights
  vector<double> &resUse=_resUse,&resWs=_resWs; vector<int> &resOff=_resOff;
  resUse.assign(_resNode.size(),1.0); resWs.assign(_resP.size(),0.0); resOff.assign(_resP.size(),0);
  int *resIb=NULL; if (IsExternal() && !_mf6Newton && (_resNode.size()>0)){string a=_mf6.GetVarAddress(_mname,"","IBOUND"); if (_mf6.HasVar(a)){resIb=_mf6.GetIntPtr(a);}}
  for (size_t e=0;e<_resNode.size();e++){
    const int r=_reduced[_resNode[e]];
    if (_colCHDTop[_resNode[e]%_ncpl] || (r<0) || ((resIb!=NULL) && (resIb[r]<=0))){resUse[e]=0.0; resOff[_resR[e]]=1;}
    resWs[_resR[e]]+=_resW[e]*resUse[e];
  }
  _resUnsent.assign(_resNode.size(),0.0);
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
    _Qres_spec+=q;
    //the seepage was computed for the share wk of the lake's cells (the seepage constant scaled by it, at the read-back): each
    //cell that can take seepage now gets its own part of it, the seepage of wk spread over the share wn of the cells able to
    //take it now, at most its rate when computed (D = max(wn, wk)); the rest is booked at the read-back (returned or not supplied)
    const int r=_resR[e]; const double wn=resOff[r]?resWs[r]:1.0,wk=(r<(int)_resWk.size())?_resWk[r]:1.0,D=max(wn,wk);
    if (D<=0.0){_resUnsent[e]=q; q=0.0;}                                         //(no cell can take it: booked at the read-back)
    else if (resOff[r] || (wk<1.0)){_resUnsent[e]=q*(1.0-wn/D); q=q*resUse[e]/D;}
    _resQ[e*max(_resNcol,1)]=q/_qf;                                              //(model units)
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
    const double zt=SWXDatum(e);                        //bed: the model top in model units (an existing model's own value); with
                                                        //:SoilZone PHYSICAL the ground of the group's HRUs, above the soil zone
    double bhead=zt+sbar/_hf;                           //stage; bed (RBOT) = that datum [model units]
    double C=S.leakance*A, dhmax=sbar;                  //leakage is at most C*sbar (water table below the bed)
    //supply limit: the largest possible leakage, C*(bhead - cell bottom)*dt, may not exceed the water the stores
    //hold, so a store can never be over-drawn whatever the head does during the solve. (Groundwater reaching the
    //surface under an empty store is still returned to Raven by the seepage drains at the model top.)
    if (dhmax>1e-6){C=min(C,max(V,0.0)/(tstep*dhmax));}
    if ((V<=0.0) || (sbar<=0.0)){C=0.0;}
    int st=max(_swxNcol,1); _swxHead[e*st]=bhead; _swxCond[e*st]=C*_hf/_qf; _swxRbot[e*st]=zt; //bed set in memory: exactly equal to stage when the store is empty
  }
  bool conv; int niter;
  std::function<bool(int,bool)> szcb; //:SoilZone PHYSICAL: interface update after every MODFLOW outer iteration
  if (_szxLink.size()>0){
    _szxHq.assign(_szxHq.size(),-HUGE_VAL);                                       //no local solve yet this step
    _szOneOK=(_maxiter>=10);
    SZLinkFractions(false); SZUpNodes();
    if (_vzMode>0){VZStep(tstep);}                                              //unsaturated columns: once per step, heads of the last step
    if (SZCapOn()){SZCapRiseWT(tstep);}                                          //capillary rise from a water table in or near the soil zone
    SZRepayRates(tstep);                                                         //LAGGED: repayment of the lagged source's debt
    SZNegCaps(tstep);                                                            //what the top cells can still give to a negative recharge
    SZRecharge(tstep,false); SZLocalSolves(tstep,false); SZSetArrays(false);   //predictor from the heads of the last step
    _szPrevQ=HUGE_VAL;
    if (SZFreezeAt()==1){ //(one update per step, an outer-iteration limit of 2 or 3: its flow test compares with the predictor's
      _szPrevQ=0.0; for (size_t e=0;e<_szxLink.size();e++){_szPrevQ+=_szxQc[e];} // total exchange)
    }
    szcb=[this,aPhinew](int it,bool c){return SZIterate(aPhinew,it,c,false,_szPrevQ,it>=_maxiter);};
  }
  else if (SZPhysical() && (_szcLink.size()>0)){ //every soil-zone link on specified heads: no interface, but their recharge,
    SZLinkFractions(false); SZUpNodes();          //its rejection and the capillary rise still apply
    if (SZCapOn()){SZCapRiseWT(tstep);}
    SZNegCaps(tstep);
    SZRecharge(tstep,true);
  }
  //target: _dvclose; a step not converged after 100 iterations continues with MODFLOW's tolerance at 10*_dvclose
  double *dvc=_mf6.GetDoublePtr("SLN_1/DVCLOSE"); //MODFLOW's own tolerance: raised for the rest of a hard step only
  //hard steps: after 100 iterations MODFLOW has not converged on (or half the iteration limit, if that is smaller; soil-zone
  //interface updates after converged iterations do not count) MODFLOW's own tolerance is relaxed tenfold for the rest of
  //the step (in the model's
  //units: with an existing model its own DVCLOSE is the target)
  double relaxTol=(IsExternal() && (_extDvclose>0))?10.0*_extDvclose/_hf:10.0*_dvclose;
  //drained columns (KeepColumnsWet) before each solve, with NEWTON and soil-zone links only
  bool wetStep=false; std::function<void(int)> wetcb;
  _keptCol.assign(_ncpl,0);
  if ((_szxLink.size()>0) && _mf6Newton){wetcb=[this,&wetStep](int it){if (KeepColumnsWet(it)>0){_wetKeptIt++; wetStep=true;}};}
  //the heads before the solve of the last outer iteration the limit allows (an unconverged step always ends with it): how far
  //that iteration still moved a head far outside the model decides whether the step diverged (below)
  _xPrevIt.clear();
  std::function<void(int)> solvecb=[this,&wetcb](int it){
    if (wetcb){wetcb(it);}
    if (it==_maxiter){const int nx=_mf6.GetVarSize(_mf6.GetVarAddress(_mname,"","X")); if (nx>0){_xPrevIt.assign(_X,_X+nx);}}
  };
  if (!_mf6.SolveTimeStep(_maxiter,conv,niter,dvc,relaxTol,min(100,max(1,_maxiter/2)),&_lastRelaxed,szcb,solvecb)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  if (wetStep){_wetKeptSteps++;}
  _lastConverged=conv; _lastIter=niter;
  if ((_szxLink.size()>0) && !_szConv){
    _szNotConv++;
    if (_szNotConv<=10){
      string where=""; if (_szResE>=0){const gw_link &L=_aLinks[_szxLink[_szResE]]; where=", largest at HRU "+to_string(_pModel->GetHydroUnit(L.k)->GetHRUID())+" in cell ("+CellID(TopNode(L.ic)/_ncpl,L.ic)+")";}
      //(the soil levels and the total exchange within the tolerances, a later test failed: that test and its link)
      string at=""; if ((_szWhy>0) && (_szWhyE>=0)){const gw_link &L=_aLinks[_szxLink[_szWhyE]]; at=" at HRU "+to_string(_pModel->GetHydroUnit(L.k)->GetHRUID())+" in cell ("+CellID(TopNode(L.ic)/_ncpl,L.ic)+")";}
      const string within="; the last change of the soil level ("+ToStr(_szResid)+" m) and of the total exchange ("+ToStr(_szDQ)+" m3/d) are within the tolerances";
      //(no update in the step: the values of an earlier step are not this step's; one update: no earlier exchange to compare with)
      const string dq=(_szIter==0)?"no interface update in this step (MODFLOW solved with the predictor's boundary)":
                      (_szWhy==1)?"soil-zone shortfall "+ToStr(_szWhyV)+" m3: at MODFLOW's heads the boundary it solved with takes more out of the columns than they hold, most"+at+within:
                      (_szWhy==2)?"connection head "+ToStr(_szWhyV)+" m from MODFLOW's head"+at+" (saturated fraction "+ToStr(_szWhyFa)+" against "+ToStr(_szWhyFb)+" at MODFLOW's head)"+within:
                      "last change of the soil level "+ToStr(_szResid)+" m"+where+(std::isfinite(_szDQ)?"; change of the total exchange "+ToStr(_szDQ)+" m3/d, tolerance "+ToStr(_szDQtol)+" m3/d":
                      "; a single interface update, no change of the total exchange to compare with (tolerance "+ToStr(_szDQtol)+" m3/d)");
      WriteWarning("CGroundwaterModel: soil-zone interface did not converge in the step ending "+_stepEnd+" ("+dq+"; continuing with MODFLOW's flows)",Options.noisy);
    }
  }
  if (!conv){
    _nNonConverged++;
    if (_nNonConverged<=10){WriteWarning("CGroundwaterModel: MODFLOW 6 did not converge in the step ending "+_stepEnd+" (continuing)",Options.noisy);}
  }
  if (!_mf6.FinalizeTimeStep()){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  for (size_t n=0;n<_reduced.size();n++){ //a diverged solve: stop at once (the far-head test, after the flows, below)
    const int r=_reduced[n];
    if ((r>=0) && !std::isfinite(_X[r])){
      ExitGracefully(("CGroundwaterModel: MODFLOW 6 heads are not finite in the step ending "+_stepEnd+" (the solve diverged): shorter time steps or a more robust "
        "solver setting may help; check the model (over-pumping: a well without AUTO_FLOW_REDUCE; negative recharge, Raven's processes taking water out of "
        "GROUNDWATER, from a cell that cannot supply it; a confined cell taking recharge with no outlet, e.g. under a top cell that dried without NEWTON)").c_str(),RUNTIME_ERR);
    }
  }

  //-- 3. boundary flows from final heads (q = hcof*h - rhs, positive into aquifer)
  _Qrch=0.0; _Qseep_tot=0.0;
  std::fill(_Qseep.begin(),_Qseep.end(),0.0);
  double gross=0.0;                                                              //[m3/d] every entry's inflow and outflow apart (below)
  _rchLostQ=0.0; int nDry=0,nChd=0;
  const bool lostOn=IsExternal() && !_mf6Newton;                                 //(recharge MODFLOW does not apply: a model without NEWTON)
  int *ib=NULL; if (lostOn){string a=_mf6.GetVarAddress(_mname,"","IBOUND"); if (_mf6.HasVar(a)){ib=_mf6.GetIntPtr(a);}}
  for (size_t i=0;i<_bndCol.size();i++){
    int ic=_bndCol[i];
    int node=_reduced[TopNode(ic)];
    double h=_X[node];
    const double qr=(_rchHcof[i]*h-_rchRhs[i])*_qf;                              //[m3/d]
    _Qrch+=qr; gross+=fabs(qr);
    double qd=(_drnHcof==NULL)?0.0:-(_drnHcof[i]*h-_drnRhs[i])*_qf; //discharge out of aquifer [m3/d]
    _Qseep[ic]=qd; _Qseep_tot+=qd; gross+=fabs(qd);
    if (lostOn){ //without NEWTON MODFLOW applies no recharge to a column whose cells are all dry (inactive), nor to a specified head
      const double d=_Qcol[ic]-qr;                                               //(signed: a negative recharge there was taken out of GROUNDWATER)
      if (fabs(d)>1e-9*max(fabs(_Qcol[ic]),1e-6)){
        _rchLostQ+=d;
        if (ib!=NULL){
          int top=-1; for (int l=0;l<_nlay;l++){const int r=_reduced[UserNode(l,ic)]; if ((r>=0) && (ib[r]!=0)){top=r; break;}}
          if (top<0){nDry++;} else if (ib[top]<0){nChd++;}
        }
      }
    }
  }
  if (lostOn){
    _rchLostVol+=_rchLostQ*tstep;
    if ((fabs(_rchLostQ)*tstep>1e-6) && (++_rchLostN<=10)){
      string why="";
      if (nDry>0){why+=to_string(nDry)+" column(s) whose cells are all inactive (dry)";}
      if (nChd>0){why+=string((nDry>0)?" and ":"")+to_string(nChd)+" column(s) whose highest active cell is a specified head";}
      if (why==""){why="columns MODFLOW took out of its solve";}
      WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+": "+ToStr(_rchLostQ*tstep)+" m3 of Raven's recharge was not applied by MODFLOW (a model without NEWTON) on "+why+
                   ((_rchLostQ<0.0)?"; negative: taken out of GROUNDWATER without MODFLOW supplying it":"; it leaves the coupled model")+
                   "; booked as 'recharge not applied by MODFLOW' (GWBudget.csv); any water those cells held when they dried is in neither water balance",false);
    }
  }
  _Qrch-=_vzUpQ;                                                                 //(the unsaturated columns' upflow out of lower cells, SZW_RAVEN)
  if (_szwQ!=NULL){const int sw=max(_szwNcol,1); for (size_t j=0;j<_szwNode.size();j++){gross+=fabs(_szwQ[j*sw])*_qf;}}
  _Qriv_in=_Qriv_out=0.0;
  std::fill(_sbGain.begin(),_sbGain.end(),0.0);
  for (size_t e=0;e<_rivNode.size();e++){
    double q=(_rivHcof[e]*_X[_reduced[_rivNode[e]]]-_rivRhs[e])*_qf; //into aquifer [m3/d]
    if (q>0){_Qriv_in+=q;} else {_Qriv_out-=q;}
    _sbGain[_rivSB[e]]-=q;
  }
  for (size_t e=0;e<_welNode.size();e++){gross+=fabs((_welHcof[e]*_X[_reduced[_welNode[e]]]-_welRhs[e])*_qf);}
  for (size_t e=0;e<_ghbNode.size();e++){gross+=fabs((_ghbHcof[e]*_X[_reduced[_ghbNode[e]]]-_ghbRhs[e])*_qf);}
  for (size_t i=0;i<_chdAct.size();i++){gross+=fabs(_chdSim[i]);}
  if (IsExternal()){ExtStreamFlows();}
  ComputeStressFlows();
  _Qres=0.0; _resShortVol=0.0; _resUndelVol=0.0;
  for (size_t e=0;e<_resNode.size();e++){
    double act=(_resHcof[e]*_X[_reduced[_resNode[e]]]-_resRhs[e])*_qf,spec=_resQ[e*max(_resNcol,1)]*_qf; //[m3/d], + into aquifer
    _Qres+=act; gross+=fabs(act);
    const double uns=(e<_resUnsent.size())?_resUnsent[e]:0.0;                   //(seepage no cell could take)
    if (_resP[_resR[e]]>=(int)_sbLossCarry.size()){continue;}
    //a gain Raven already gave the lake (spec<0) that MODFLOW could not fully supply because the cell was drying
    //(AUTO_FLOW_REDUCE): the difference is taken back from the water reaching that subbasin in this step's routing,
    //like a river loss, so no water is created (reported as "reservoir gain not supplied by MODFLOW")
    if ((spec<0.0) && (act>spec)){double S=(act-spec)*tstep; _sbLossCarry[_resP[_resR[e]]]+=S; _resShortVol+=S;}
    if (uns<0.0){double S=-uns*tstep; _sbLossCarry[_resP[_resR[e]]]+=S; _resShortVol+=S;}
    //a loss Raven already took out of the lake (spec>0) that MODFLOW did not take (a cell that became a specified head or
    //inactive during the solve; or no cell could take it): returned to the water reaching that subbasin in this step's
    //routing, which flows on to the reservoir (reported as "reservoir seepage not delivered by MODFLOW")
    double U=0.0;
    if ((spec>0.0) && (act<spec)){U+=(spec-act)*tstep;}
    if (uns>0.0){U+=uns*tstep;}
    if (U>0.0){_sbLossCarry[_resP[_resR[e]]]-=U; _resUndelVol+=U;}
  }
  _resUndelTot+=_resUndelVol; _resShortTot+=_resShortVol;
  if ((_resShortVol>1e-6) && (++_resShortN<=10)){
    WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+": "+ToStr(_resShortVol)+" m3 of reservoir gain was not supplied by MODFLOW (cells under the lake that "
                 "were drying, specified heads or inactive); taken back from the water reaching the reservoir's subbasin and booked as 'reservoir gain not supplied by MODFLOW' (GWBudget.csv)"+
                 ((_resShortN==10)?" (further such warnings suppressed)":""),false);
  }
  if ((_resUndelVol>1e-6) && (++_resUndelN<=10)){
    WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+": "+ToStr(_resUndelVol)+" m3 of reservoir seepage was not delivered by MODFLOW (cells under the "
                 "lake that are specified heads or inactive); returned to the reservoir's subbasin and booked as 'reservoir seepage not delivered by MODFLOW' (GWBudget.csv)"+
                 ((_resUndelN==10)?" (further such warnings suppressed)":""),false);
  }
  for (size_t r=0;r<_resP.size();r++){ //local groundwater head seen by each reservoir during this step's routing
    //a lake above a water table deeper than its bed is disconnected: the seepage then depends on the bed, not on the
    //water-table depth (as for MODFLOW's river boundary). The bed is the ground: the model top, or with :SoilZone PHYSICAL
    //the lake HRU's ground above its soil zone (as for surface stores, SWXDatum); unless the lake stands below it (a DEM
    //that recorded a higher water surface): then the bed is taken at the lake level, so a deep water table gives no
    //exchange rather than a false gain
    CReservoir *pRes=_pModel->GetSubBasin(_resP[r])->GetReservoir();
    double stage=pRes->GetResStage();
    const int kl=pRes->GetHRUIndex(); const bool szk=SZPhysical() && (kl>=0) && (kl<(int)_szOn.size()) && _szOn[kl];
    double h=0; for (size_t e=0;e<_resNode.size();e++){if (_resR[e]==(int)r){int ic=_resNode[e]%_ncpl; const double zg=szk?SZBase(kl,ic)+_szD[kl]:_top[ic]; h+=_resW[e]*max(WaterTable(ic),min(zg,stage));}}
    double F=(r<_resA.size())?_resA[r]:1.0,wsr=1.0;                             //(the share of the lake's bed that can take seepage)
    if ((r<_resOff.size()) && _resOff[r]){ //(cells that could not take seepage this step are left out: the head of the others;
      h=0; for (size_t e=0;e<_resNode.size();e++){if ((_resR[e]==(int)r) && (_resUse[e]>0.0)){int ic=_resNode[e]%_ncpl; const double zg=szk?SZBase(kl,ic)+_szD[kl]:_top[ic]; h+=_resW[e]*max(WaterTable(ic),min(zg,stage));}}
      h=(_resWs[r]>0.0)?h/_resWs[r]:stage;                                       // none: the lake is disconnected, no seepage)
      wsr=max(_resWs[r],0.0); F*=wsr;
    }
    if (r<_resWk.size()){_resWk[r]=wsr;}
    //(never 0: Raven's reservoir keeps its last seepage when the constant is 0; 1e-12 of it with the level at the stage is none)
    pRes->SetGWParameters(((r<_resK0.size())?_resK0[r]:pRes->GetSeepageConstant())*max(F,1e-12),h);
  }
  _Qswx_in=_Qswx_out=0.0; _swxAppliedVol=0.0;
  _swxQ.assign(_swxCol.size(),0.0);
  for (size_t e=0;e<_swxCol.size();e++){
    double q=(_swxHcof[e]*_X[_reduced[TopNode(_swxCol[e])]]-_swxRhs[e])*_qf; //[m3/d] into aquifer (MODFLOW's value)
    _swxQ[e]=q;
    if (q>0){_Qswx_in+=q;} else {_Qswx_out-=q;}
  }
  _Qevt=0.0;
  //MODFLOW's EVT moves an entry whose cell is inactive (a top cell that dried, without NEWTON) to the highest active cell
  //of the column, and computes the ET there: the flow is read back with the head of the cell MODFLOW used (NODELIST)
  int *evtNl=(_evtCol.size()>0)?_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,"EVT_RAVEN","NODELIST")):NULL;
  for (size_t i=0;i<_evtCol.size();i++){
    int ic=_evtCol[i];
    const int r=((evtNl!=NULL) && (evtNl[i]>0))?evtNl[i]-1:_reduced[TopNode(ic)];
    const double qe=(_evtHcof[i]*_X[r]-_evtRhs[i])*_qf;
    _Qevt-=qe; gross+=fabs(qe);                                                  //EVT flow is out of the aquifer
  }
  _Qsz_in=_Qsz_out=0.0;
  if (SZPhysical()){SZCHDSolves(_szTstep);}                                      //columns over specified heads: exact, after the solve
  for (size_t e=0;e<_szxLink.size();e++){ //soil zone <-> aquifer, MODFLOW's booked flow (+ into aquifer)
    if (_colCHDTop[_aLinks[_szxLink[e]].ic]){ //a specified head this step (existing models): the local solve's exchange with it
      const gw_link &L=_aLinks[_szxLink[e]];
      const double q=SZChdExchange(L.k,_szxCol[e],_szxCell[e],_szxRej[e],_X[_reduced[TopNode(L.ic)]]*_hf,_szTstep,(e<_szxLagD.size())?_szxLagD[e]:0.0); _szxQ[e]=q;
      if (q>0){_QszCHD_in+=q;} else {_QszCHD_out-=q;}
      for (int j=0;j<_szNLayers;j++){_QszCHD_out+=_szxCapL[e*_szNLayers+j];}
      continue;
    }
    double h=_X[_reduced[TopNode(_aLinks[_szxLink[e]].ic)]];
    double q=(_szxHcof[e]*h-_szxRhs[e])*_qf+(_szqHcof[e]*h-_szqRhs[e])*_qf+(_szdHcof[e]*h-_szdRhs[e])*_qf; //SZX_RAVEN + SZQ_RAVEN + SZD_RAVEN
    if ((e<_szxUpCHD.size()) && _szxUpCHD[e]){for (int j=0;j<_szNLayers;j++){_QszCHD_out+=_szxCapL[e*_szNLayers+j];}} //(rise out of a specified head below a dried top cell)
    const bool inQ=(e>=_szxUp.size()) || (_szxUp[e]==TopNode(_aLinks[_szxLink[e]].ic)); //(else in SZW_RAVEN, out of a lower cell)
    if (inQ && (e<_szxCapQ.size())){q+=_szxCapQ[e];}                             //(the capillary rise from the water table, a fixed rate in SZQ_RAVEN, is booked apart)
    if (inQ && (e<_szxRepay.size())){q+=_szxRepay[e];}                           //(as is the repayment of the lagged source's debt)
    _szxQ[e]=q; if (q>0){_Qsz_in+=q;} else {_Qsz_out-=q;}
    gross+=fabs(q);                                                              //(each link net)
  }
  _Qsto=0.0;
  for (size_t e=0;e<_swxQ.size();e++){gross+=fabs(_swxQ[e]);}
  gross+=_Qriv_in+_Qriv_out;                                       //(rivers and the stream package: each entry's sign apart)
  if (IsExternal()){gross+=_QabsOther;}
  //MODFLOW's residual: inflows + storage release - outflows, one sum for every branch that knows the storage flow (the steady
  //existing model without STO6 has storage 0)
  auto resid=[this](){
    double r=_Qrch+_Qsto-_Qseep_tot+_Qriv_in-_Qriv_out-_Qevt+_Qwel+_Qghb+_Qchd+_Qswx_in-_Qswx_out+_Qres+_Qsz_in-_Qsz_out-_szCapWT-_szRepayQ; //(_QrchCHD never enters MODFLOW)
    if (IsExternal()){r+=_Qother;}
    return r;
  };
  if (_strgss!=NULL){
    int nn=0; for (int n=0;n<_nlay*_ncpl;n++){if (_idomain[n]>0){nn++;}}
    for (int n=0;n<nn;n++){_Qsto+=(_strgss[n]+_strgsy[n])*_qf; gross+=(fabs(_strgss[n])+fabs(_strgsy[n]))*_qf;}
    _Qerr=resid();
  }
  else if (IsExternal() && _extNoSto){_Qsto=0.0; _Qerr=resid();} //an existing model without a storage package is steady state: no storage flow
  else {_Qsto=0.0; _Qsto=-resid(); _Qerr=0.0; gross+=fabs(_Qsto);}             //(inferred storage)
  //gross flow through the aquifer: every entry's inflow and outflow apart (MODFLOW's TOTAL IN + TOTAL OUT), at least 1 m3/d;
  //computed once, before the guard: the far-head test and GWBudget's "MF6 balance error [%]" (signed
  //package totals cancelled: a well taking 1e4 m3/d resupplied by a specified head gave 1 m3/d)
  _Qgross=max(gross,1.0);

  { //far heads: an unconverged step that holds a head far outside the model's elevations (more than 10 times its relief, at
    //least 1 km, beyond its highest top or lowest bottom) is carried on only if every far head that moved in the step (not
    //within the margin of where it was at the start of the step, the start taken only from converged steps) was moved by the
    //last outer iteration by at most the margin, and MODFLOW's residual is at most 1 % of the gross flow; otherwise the run
    //stops. The same heads on a converged step are what MODFLOW computes from the input
    //(e.g. recharge into a confined cell under a dried top cell without NEWTON): a warning once per run. HDRY is skipped
    double zmin,zmax; const double m=FarMargin(zmin,zmax);
    int nf=-1,n3=-1; double dmax=0.0,mv=0.0; bool moved=false;
    for (size_t n=0;n<_reduced.size();n++){
      const int r=_reduced[n]; if (r<0){continue;}
      if (!(fabs(_X[r])<1e29)){continue;}                                        //(HDRY, in the model's units)
      const double h=_X[r]*_hf,d=FarDistance(h,zmin,zmax);
      if (!(d>m)){continue;}
      if (d>dmax){dmax=d; nf=(int)n;}
      if (!conv){
        bool mvd=true;
        if (r<(int)_xStart.size()){
          const double x0=_xStart[r]*_hf;
          const int sh=(h>zmax)?1:-1,s0=(x0>zmax)?1:((x0<zmin)?-1:0);
          mvd=(fabs(h-x0)>m) && ((sh!=s0) || (d>FarDistance(x0,zmin,zmax)+m));
        }
        if (mvd){
          moved=true;
          const double dm=(r<(int)_xPrevIt.size())?fabs(_X[r]-_xPrevIt[r])*_hf:HUGE_VAL;
          if ((n3<0) || (dm>mv)){mv=dm; n3=(int)n;}
        }
      }
    }
    if (nf>=0){
      const int nn=(n3>=0)?n3:nf; const double hn=_X[_reduced[nn]]*_hf;
      const string what="MODFLOW 6 head "+to_string(hn)+" m in cell ("+CellID(nn/_ncpl,nn%_ncpl)+") in the step ending "+_stepEnd+" lies far outside the model's elevations ("+
                        to_string(zmin)+" to "+to_string(zmax)+" m)";
      const string lim=IsExternal()?"OUTER_MAXIMUM in the model's IMS file":":SolverMaxOuterIterations";
      const string causes="; check the model (over-pumping: a well without AUTO_FLOW_REDUCE; negative recharge, Raven's processes taking water out of GROUNDWATER, "
                          "from a cell that cannot supply it; a confined cell taking recharge with no outlet, e.g. under a top cell that dried without NEWTON)";
      if (conv){
        if (!_farWarned){_farWarned=true; WriteWarning("CGroundwaterModel: "+what+" in a converged step"+causes,Options.noisy);}
      }
      else {
        const double pc=100.0*fabs(_Qerr)/_Qgross;
        const string res="MODFLOW's residual "+ToStr(_Qerr)+" m3/d, "+ToStr(100.0*_Qerr/_Qgross)+" % of the gross flow "+ToStr(_Qgross)+" m3/d";
        string why="";
        if (moved && (mv>m)){why=" and the step did not converge: the solve diverged"+(std::isfinite(mv)?" (its last outer iteration moved it "+to_string(mv)+" m)":string(""))+" ("+res+")";}
        else if (pc>100.0*SZ_FAR_RESID){why=" and the step did not converge with "+res+" (more than "+ToStr(100.0*SZ_FAR_RESID)+" %): its heads are wrong"+
                                         (moved?" (the last outer iteration moved this head "+to_string(mv)+" m)":string(""));}
        if (why!=""){ExitGracefully(("CGroundwaterModel: "+what+why+"; shorter time steps or a larger "+lim+" may help"+causes).c_str(),RUNTIME_ERR);}
        if (++_farSlowN<=10){WriteWarning("CGroundwaterModel: "+what+" in a step MODFLOW did not converge"+(moved?", but its last outer iteration moved the head only "+to_string(mv)+" m (at most the margin, "+to_string(m)+" m)":string(" (as at the start of the step)"))+
                                          "; "+res+" (at most "+ToStr(100.0*SZ_FAR_RESID)+" %) is booked as its balance error; a larger "+lim+" may let such steps converge",Options.noisy);}
      }
    }
    if (conv){const int nx=_mf6.GetVarSize(_mf6.GetVarAddress(_mname,"","X")); if (nx>0){_xStart.assign(_X,_X+nx);}} //(the start of the next step: converged steps only)
  }

  //-- :SoilZone PHYSICAL guard: water that MODFLOW's books do not back is not handed to Raven. Two sources: on an
  //   unconverged step, MODFLOW's residual (more outflow booked than its inflows and storage release supply; the soil
  //   zone otherwise received flows far beyond its local solve, which overflowed: outflow 17 times the
  //   iterated run's with an outer-iteration limit of 12); on any step, the soil-zone shortfall (water MODFLOW took from
  //   columns that did not hold it, which it then discharged or returned as upflow with balanced books). All discharge to Raven (seepage, gains of enabled subbasins, surface-store gains, the exchange into the soil
  //   zone) is reduced in proportion; what it cannot cover (little discharge this step) stays in the aquifer
  //   (reported: MF6 balance error, soil-zone shortfall; the withheld volume in GWBudget.csv)
  _szWithheldV=_szWithheldG=_szWithheldS=0.0;
  if (SZPhysical()){
    double R=SZTrialShortfall(tstep)/tstep;                                       //[m3/d]
    if (!_lastConverged && (_Qerr<0.0)){R-=_Qerr;}
    if (R>0.0){
      double gR=0.0,gS=0.0;
      for (size_t p=0;p<_sbGain.size();p++){if (_pModel->GetSubBasin((int)p)->IsEnabled()){gR+=max(_sbGain[p],0.0);}}
      for (size_t e=0;e<_swxQ.size();e++){gS+=max(-_swxQ[e],0.0);}
      double gains=gR+gS;
      double up=0.0; for (size_t e=0;e<_szxQ.size();e++){if (!_colCHDTop[_aLinks[_szxLink[e]].ic]){up+=max(-_szxQ[e],0.0);}} //(not from specified heads)
      double D=_Qseep_tot+gains+up,seep0=_Qseep_tot;
      if (D>0.0){
        double f=max(1.0-R/D,0.0);
        for (int ic=0;ic<_ncpl;ic++){_Qseep[ic]*=f;} _Qseep_tot*=f;
        for (size_t p=0;p<_sbGain.size();p++){if (_sbGain[p]>0){_sbGain[p]*=f;}}
        for (size_t e=0;e<_swxQ.size();e++){if (_swxQ[e]<0){_swxQ[e]*=f;}}
        _Qsz_out=0.0; for (size_t e=0;e<_szxQ.size();e++){if (_colCHDTop[_aLinks[_szxLink[e]].ic]){continue;} if (_szxQ[e]<0){_szxQ[e]*=f;} if (_szxQ[e]<0){_Qsz_out-=_szxQ[e];}}
        _szWithheldV=(1.0-f)*(seep0+up)*tstep; _szWithheldG=(1.0-f)*gR*tstep; _szWithheldS=(1.0-f)*gS*tstep; //(split: MODFLOW's own budget
                                                    //books the gains apart, and the river and store exchanges are checked apart)
        if ((_szWithheldV+_szWithheldG+_szWithheldS>1e-6) && (++_szGuardN<=10)){
          WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+": "+ToStr(_szWithheldV+_szWithheldG+_szWithheldS)+" m3 of discharge to Raven withheld for "
                       "MODFLOW's residual and the soil-zone shortfall (MODFLOW unconverged or an interface boundary it kept)"+(_szGuardN==10?" (further such warnings suppressed)":""),false);
        }
      }
    }
  }
  if (!_lastConverged && !SZPhysical() && (_Qerr<0.0)){
    //-- without :SoilZone PHYSICAL, on a step MODFLOW did not converge (CONTINUE-like): every flow MODFLOW
    //   booked is given to Raven, except its negative residual (more outflow booked than its inflows and storage release
    //   supply), which is withheld from all discharge to Raven (seepage, gains of enabled subbasins, surface-store gains) in
    //   proportion, as the PHYSICAL guard above does; what that discharge cannot cover stays in the aquifer (MODFLOW's balance
    //   error). The withheld volume is a run total in GWModelSummary.txt; per step it is the GWBudget row's residual, recharge
    //   + storage release + other inflows - discharge - MF6 balance error: the row's discharge columns (seepage to Raven,
    //   aquifer discharge to rivers and to surface stores) are the discharge handed to Raven, after the withholding.
    const double R=-_Qerr;                                                       //[m3/d]
    double gR=0.0,gS=0.0;
    for (size_t p=0;p<_sbGain.size();p++){if (_pModel->GetSubBasin((int)p)->IsEnabled()){gR+=max(_sbGain[p],0.0);}}
    for (size_t e=0;e<_swxQ.size();e++){gS+=max(-_swxQ[e],0.0);}               //(discharge into surface stores)
    const double D=_Qseep_tot+gR+gS;
    if (D>0.0){
      const double f=max(1.0-R/D,0.0);
      for (int ic=0;ic<_ncpl;ic++){_Qseep[ic]*=f;} _Qseep_tot*=f;
      for (size_t p=0;p<_sbGain.size();p++){if (_sbGain[p]>0){_sbGain[p]*=f;}}
      for (size_t e=0;e<_swxQ.size();e++){if (_swxQ[e]<0){_swxQ[e]*=f;}}
      _Qriv_out-=(1.0-f)*gR; _Qswx_out-=(1.0-f)*gS;                              //(GWBudget: the discharge handed to Raven)
      const double w=(1.0-f)*D*tstep;                                            //[m3]
      _ncWithheldVol+=w; _ncWithheldN++;
      if ((w>1e-6) && (++_szGuardN<=10)){
        WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+" (MODFLOW unconverged): "+ToStr(w)+" m3 of discharge to Raven withheld for "
                     "MODFLOW's negative residual ("+ToStr(_Qerr)+" m3/d, booked as its balance error)"+(_szGuardN==10?" (further such warnings suppressed)":""),false);
      }
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
  //-- soil zone: MODFLOW's exchange, rejected recharge and overflow to Raven (:SoilZone PHYSICAL)
  if (SZPhysical()){SZApply(aPhinew,tstep);}
  //-- 4. seepage -> Raven store of overlapping HRUs (share of column by overlap area)
  _returnVol=0.0;
  for (size_t li=0;li<_aLinks.size();li++){
    const gw_link &L=_aLinks[li];
    if (_Qseep[L.ic]==0.0){continue;}
    double cover=_colCover[L.ic];
    if (SZPhysical()){ //the drain serves only the HRUs without a soil zone in the column
      if (_szOn[L.k]){continue;}
      cover=_colCoverNS[L.ic]; if (cover<=0.0){continue;}
    }
    double V=_Qseep[L.ic]*tstep*L.area/cover; //[m3]
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
  //head boundaries: a cell whose bottom lies at or above the step's head is switched off for the step (it neither holds
  //the water table nor feeds water): a GHB gets conductance 0, a CHD entry is left out of CHD_RAVEN's list (its cell is an
  //ordinary active cell for the step). LAYERS TOP: only the highest placed layer whose bottom lies below the head
  _bndOffNow=0;
  int pb=-1,pic=-1; bool taken=false;
  for (size_t e=0;e<_ghbNode.size();e++){
    const gw_boundary &B=_aBounds[_ghbB[e]]; const int n=_ghbNode[e];
    double h=B.head; if (B.pHead!=NULL){double v=B.pHead->GetAvgValue(t,tstep); if (v!=RAV_BLANK_DATA){h=v;}}
    if ((_ghbB[e]!=pb) || (n%_ncpl!=pic)){pb=_ghbB[e]; pic=n%_ncpl; taken=false;}
    const bool on=(h>_botm[n]) && (B.allLayers || !taken); if (on){taken=true;}
    if (!(h>_botm[n])){_bndOffNow++;}
    _ghbHead[e*max(_ghbNcol,1)]=max(h,_botm[n]+0.01); _ghbCond[e*max(_ghbNcol,1)]=on?B.cond:0.0;
  }
  if (_chdNode.size()>0){
    vector<int> act; vector<double> hs(_chdNode.size(),0.0); vector<char> on(_chdNode.size(),0);
    pb=-1; pic=-1; taken=false;
    for (size_t e=0;e<_chdNode.size();e++){
      const gw_boundary &B=_aBounds[_chdB[e]]; const int n=_chdNode[e];
      double h=B.head; if (B.pHead!=NULL){double v=B.pHead->GetAvgValue(t,tstep); if (v!=RAV_BLANK_DATA){h=v;}}
      if ((_chdB[e]!=pb) || (n%_ncpl!=pic)){pb=_chdB[e]; pic=n%_ncpl; taken=false;}
      on[e]=(h>_botm[n]) && (B.allLayers || !taken); if (on[e]){taken=true; act.push_back((int)e);}
      if (!(h>_botm[n])){_bndOffNow++;}
      hs[e]=max(h,_botm[n]+0.01);
    }
    if ((act.size()==_chdNode.size()) && !_chdPacked){ //every entry held, MODFLOW's list as read from gwf.chd
      for (size_t e=0;e<_chdNode.size();e++){_chdHead[e*max(_chdNcol,1)]=hs[e];}
    }
    else { //the list of the entries held this step: NODELIST (reduced node numbers, from 1), heads and NBOUND; IBOUND marks the
      //held cells (MODFLOW's CHD sets it only when it reads a period) and makes the others ordinary active cells
      ExitGracefullyIf((_chdNb==NULL) || (_chdNl==NULL) || (_chdIb==NULL),"CGroundwaterModel: cannot switch specified-head cells off: MODFLOW 6's CHD list or IBOUND is not accessible",RUNTIME_ERR);
      if (_chdIbc==0){const int r0=_reduced[_chdNode[0]]; _chdIbc=((r0>=0) && (_chdIb[r0]<0))?_chdIb[r0]:-1;} //(MODFLOW's value of a held cell, read while its list is whole)
      for (size_t i=0;i<act.size();i++){const int e=act[i]; _chdNl[i]=_reduced[_chdNode[e]]+1; _chdHead[i*max(_chdNcol,1)]=hs[e];}
      *_chdNb=(int)act.size();
      for (size_t e=0;e<_chdNode.size();e++){const int r=_reduced[_chdNode[e]]; if (r>=0){_chdIb[r]=on[e]?_chdIbc:1;}}
      _chdPacked=(act.size()!=_chdNode.size());
    }
    _chdAct=act;
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
  for (size_t i=0;i<_chdAct.size();i++){_Qchd+=_chdSim[i];}
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
  if (_szxLink.size()>0){ //soil zone not yet active: SZX acts as a drain at each link's base (no soil water, no supply)
    SZLinkFractions(false); SZLocalSolves(1.0,true); SZSetArrays(true); _szPrevQ=HUGE_VAL; _szExtra=0; _szIter=0;
    std::function<bool(int,bool)> cb=[this](int it,bool c){return SZIterate(NULL,it,c,true,_szPrevQ,it>=_maxiter);};
    if (!_mf6.SolveTimeStep(_maxiter,conv,niter,NULL,0.0,0,NULL,cb)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
  }
  else if (!_mf6.SolveTimeStep(_maxiter,conv,niter)){ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()).c_str(),RUNTIME_ERR);}
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
  //(12 significant digits; 17 with :SoilZone PHYSICAL, whose restart then starts from MODFLOW's heads as they were, as the
  //starting heads Raven writes for MODFLOW, raven.ic and the generated model's STRT)
  const int hp=SZPhysical()?17:12;
  if ((_gridType==GWGRID_REGULAR) || (IsExternal() && (_disType==0))){RVC<<":GWHeads "<<_nlay<<" "<<_nrow<<" "<<_ncol<<endl<<setprecision(hp);}
  else {RVC<<":GWHeads "<<_nlay<<" 1 "<<_ncpl<<endl<<setprecision(hp);} //unstructured: one row of cells
  for (int n=0;n<_nlay*_ncpl;n++){
    double h=(_reduced[n]>=0)?_X[_reduced[n]]:_strt[n];
    RVC<<((n%10==0)?"  ":" ")<<h; if ((n%10==9) || (n==_nlay*_ncpl-1)){RVC<<endl;}
  }
  RVC<<":EndGWHeads"<<endl<<setprecision(12); //
  if (_resP.size()>0){ //reservoir seepage computed in the last step, not yet delivered to the aquifer
    RVC<<":GWReservoirSeepage   # subbasin ID, seepage in transit [m3], share of the lake's cells able to take seepage when it was computed"<<endl;
    for (size_t r=0;r<_resP.size();r++){RVC<<"  "<<_pModel->GetSubBasin(_resP[r])->GetID()<<", "<<setprecision(15)<<_pModel->GetSubBasin(_resP[r])->GetReservoir()->GetReservoirGWLosses(1.0)<<", "<<((r<_resWk.size())?_resWk[r]:1.0)<<endl;}
    RVC<<":EndGWReservoirSeepage"<<endl;
  }
  bool anyCarry=false; for (size_t p=0;p<_sbLossRemain.size();p++){if (_sbLossRemain[p]>0){anyCarry=true;}}
  if (anyCarry){ //river loss not yet supplied at the end of the run
    RVC<<":GWRiverLossCarry   # subbasin ID, river loss to be taken from the river in the next step [m3]"<<endl;
    for (size_t p=0;p<_sbLossRemain.size();p++){if (_sbLossRemain[p]>0){RVC<<"  "<<_pModel->GetSubBasin((int)p)->GetID()<<", "<<setprecision(15)<<_sbLossRemain[p]<<endl;}}
    RVC<<":EndGWRiverLossCarry"<<endl;
  }
  if (SZPhysical() && _szColInit){ //water of the soil-zone columns (Raven's SOIL[m] holds only their mean)
    RVC<<":GWSoilZoneColumns   # HRU ID, cell (grid column number; 0: rest of the HRU), water of SOIL["<<_szNLayers-1<<"] up to SOIL[0] in the column [mm]"<<endl;
    for (size_t i=0;i<_aqHRUs.size();i++){
      int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
      const vector<int> cell=SZColumnCells(k);
      for (int c=0;c<SZNCols(k);c++){
        RVC<<"  "<<_pModel->GetHydroUnit(k)->GetHRUID()<<", "<<cell[c]<<setprecision(15);
        for (int j=0;j<_szNLayers;j++){RVC<<", "<<_szColPhi[(_szColOff[k]+c)*_szNLayers+j];}
        RVC<<endl;
      }
    }
    RVC<<":EndGWSoilZoneColumns"<<endl;
    bool anyVz=false; for (size_t i=0;i<_vzcCol.size();i++){if (_vzcCol[i].active){anyVz=true;}}
    if (anyVz){ //unsaturated columns (:SoilZoneVadose); with depth classes (:SoilZoneVadoseClasses K > 1) one line per class column
      const bool cls=(_vzClasses>1);
      if (!cls){RVC<<":GWVadoseColumns   # HRU ID, cell (grid column number), water table [m], water contents of the nodes from the top [-]"<<endl;}
      else     {RVC<<":GWVadoseClassColumns   # HRU ID, cell (grid column number), depth class, classes of the cell, water table [m], water contents of the nodes from the top [-]"<<endl;}
      for (size_t e=0;e<_szxLink.size();e++){
        const gw_link &L=_aLinks[_szxLink[e]];
        for (int i=_vzOff[e];i<_vzOff[e+1];i++){
          const gwvz::vz_column &c=_vzcCol[i]; if (!c.active){continue;}
          RVC<<"  "<<_pModel->GetHydroUnit(L.k)->GetHRUID()<<", "<<L.ic+1;
          if (cls){RVC<<", "<<i-_vzOff[e]+1<<", "<<_vzOff[e+1]-_vzOff[e];}
          RVC<<setprecision(15)<<", "<<c.zbot;
          for (size_t q=0;q<c.th.size();q++){RVC<<", "<<c.th[q];}
          RVC<<endl;
        }
      }
      RVC<<(cls?":EndGWVadoseClassColumns":":EndGWVadoseColumns")<<endl;
      //each column's top (matched on a restart: another :SoilZoneSubgrid or ground gives other tops) and its water
      RVC<<":GWVadoseColumnTops   # HRU ID, cell (grid column number), depth class (1 without depth classes), column top [m], water above the column's reference profile [m3]"<<endl;
      for (size_t e=0;e<_szxLink.size();e++){
        const gw_link &L=_aLinks[_szxLink[e]]; const double a=_szColA[_szxCol[e]];
        for (int i=_vzOff[e];i<_vzOff[e+1];i++){
          const gwvz::vz_column &c=_vzcCol[i]; if (!c.active){continue;}
          RVC<<"  "<<_pModel->GetHydroUnit(L.k)->GetHRUID()<<", "<<L.ic+1<<", "<<i-_vzOff[e]+1<<setprecision(15)<<", "<<c.ztop<<", "<<gwvz::RelStorage(c,_vzStrat[e])*(a*_vzcA[i])<<endl;
        }
      }
      RVC<<":EndGWVadoseColumnTops"<<endl;
    }
    if ((_szTiming==1) && (_szxQprev.size()==_szxLink.size())){ //LAGGED: the last step's exchange is the next step's source
      RVC<<":GWSoilZoneLag   # HRU ID, cell (grid column number), exchange of the last step into the aquifer [m3/d], debt the aquifer still repays [m3]"<<endl;
      for (size_t e=0;e<_szxLink.size();e++){
        const gw_link &L=_aLinks[_szxLink[e]];
        //over a specified head in the last step the exchange was exact, not a lagged source: a continued run does not apply it,
        //and a restarted existing model knows its specified heads only after the first step's lagged source
        const double q=_colCHDTop[L.ic]?0.0:_szxQprev[e];
        RVC<<"  "<<_pModel->GetHydroUnit(L.k)->GetHRUID()<<", "<<L.ic+1<<setprecision(15)<<", "<<q<<", "<<((e<_szxDebt.size())?_szxDebt[e]:0.0)<<endl;
      }
      RVC<<":EndGWSoilZoneLag"<<endl;
    }
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
    if (SZPhysical()){_BUDGET<<","<<_Qsz_in<<","<<_Qsz_out<<","<<_szRejVol<<","<<_szOvfVol<<","<<_szAppliedVol<<","<<_szShortVol<<","<<_szIter<<","<<(_szConv?1:0)<<","<<_szResid<<","<<_szMovedV<<","<<_szWithheldV<<","<<_szWithheldG<<","<<_szWithheldS<<","<<_szNegUnheldV;}
    if (SZPhysical() && (_vzMode>0)){_BUDGET<<","<<_vzInVol-_vzBackV-_vzBackWV<<","<<_vzOutVol-_vzBackV<<","<<_vzCapVol<<","<<_vzStoreVol;} //(water handed back: rule (i)
                                                                                 // reduces the columns' inflow and outflow, rule (ii) their inflow)
    if (SZPhysical() && SZCapOn()){_BUDGET<<","<<_szCapWT*_szTstep;}
    if (SZPhysical() && (_szTiming==1)){_BUDGET<<","<<_szDebtVol<<","<<_szRepayQ*_szTstep;}
    if (SZPhysical() && SZChdColumns()){_BUDGET<<","<<_QszCHD_in<<","<<_QszCHD_out;}
    if (IsExternal() && !_mf6Newton){_BUDGET<<","<<_rchLostQ;}
    if (IsExternal() && (_resNode.size()>0)){_BUDGET<<","<<_resUndelVol;}
    _BUDGET<<endl;
  }
  if (_budgetLineOpen){_szMovedV=0.0;} //(called before the first step too: a hotstart sync then counts in step 1)
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
         <<"reservoir seepage sent [m3/d],reservoir seepage delivered [m3/d],reservoir seepage awaiting delivery [m3],reservoir gain not supplied by MODFLOW [m3]"<<(IsExternal()?",other model packages [m3/d]":"")
         <<(SZPhysical()?",soil zone to aquifer [m3/d],aquifer to soil zone [m3/d],rejected recharge [m3],soil-zone overflow [m3],soil-zone exchange applied to Raven [m3],"
                         "soil-zone shortfall [m3],interface updates,interface converged,interface change [m],soil-zone layer redistribution [m3],seepage and soil inflow withheld [m3],river gains withheld [m3],surface-store gains withheld [m3],negative recharge the soil zone did not hold [m3]":"")
         <<((SZPhysical() && (_vzMode>0))?",unsaturated zone inflow [m3],unsaturated zone to aquifer [m3],capillary rise to soil zone [m3],unsaturated zone storage [m3]":"")
         <<((SZPhysical() && SZCapOn())?",capillary rise from the water table [m3]":"")
         <<((SZPhysical() && (_szTiming==1))?",lagged-exchange debt [m3],lagged-exchange repayment [m3]":"")
         <<((SZPhysical() && SZChdColumns())?",soil zone to specified heads [m3/d],specified heads to soil zone [m3/d]":"")
         <<((IsExternal() && !_mf6Newton)?",recharge not applied by MODFLOW [m3/d]":"")
         <<((IsExternal() && (_resNode.size()>0))?",reservoir seepage not delivered by MODFLOW [m3]":"")<<endl;
  if (SZPhysical()){
    _SZOUT.open(FilenamePrepare("GWSoilExchange.csv",Options).c_str());
    ExitGracefullyIf(_SZOUT.fail(),"CGroundwaterModel: cannot open GWSoilExchange.csv",FILE_OPEN_ERR);
    _SZOUT<<"time [d],date,HRU_ID,soil-zone base [m],soil water level [m],water table [m],saturated fraction [-],"
          <<"soil zone to aquifer [m3/d],aquifer to soil zone [m3/d],rejected recharge [mm],overflow [mm],soil-zone water [mm],soil-zone capacity [mm],interface updates,interface change [m]"
          <<(SZChdColumns()?",soil zone to specified heads [m3/d],specified heads to soil zone [m3/d]":"")<<endl;
  }
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
  const double scale=_Qgross;                                                    //gross flow through the aquifer [m3/d] (Exchange, before the guard)
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
  if (SZPhysical()){SZWriteStep(t,dstr);}

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
  _summaryFile=FilenamePrepare("GWModelSummary.txt",Options);
  ofstream S(_summaryFile.c_str());
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
    if (_extTSShifted>0){S<<"  time series           : "<<_extTSShifted<<" TS6 files shifted to Raven's start ("<<_extOffset<<" d after the model's time zero)"<<endl;}
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
    if (_evtCol.size()>0){S<<"  water-table ET cells  : "<<_evtCol.size()<<" (EVT_RAVEN, raven.evt)"<<endl;}
    if (_swxCol.size()>0){S<<"  surface-water exchange cells: "<<_swxCol.size()<<" (SWX_RAVEN, raven.swx)"<<endl;}
    if (_resExchange){S<<"  coupled reservoirs / cells  : "<<_resP.size()<<" / "<<_resNode.size()<<" (RES_RAVEN, raven.res)"<<endl;}
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
  if (_inCRS.IsSet()){S<<"  input coordinates     : "<<_inCRS.Describe()<<" (:InputCRS)"<<endl;}
  if (_writeCBC){S<<"  budget file           : gwf.cbc every "<<max(1,_headSaveFreq)<<" steps (with specific discharge and saturation), porosity in gwf.porosity (:WriteGWBudgetFile)"<<endl;}
  S<<"  HRU/grid linkage      : "<<(_cacheHit?"read from cache":"computed")<<" in "<<_geomSeconds<<" s"<<endl;
  for (size_t o=0;o<_aObsWells.size();o++){
    int n=_aObsWells[o].node;
    S<<"    "<<_aObsWells[o].name<<": "<<((n<0)?string("outside model"):"layer "+to_string(n/_ncpl+1)+", MODFLOW cell ("+CellID(n/_ncpl,n%_ncpl)+")")<<endl;
  }
  if (_doSteady){S<<"  initialization        : steady state, recharge "<<_ssRecharge<<" mm/d"<<endl;}
  else if (_hotHeads.size()>0){S<<"  initialization        : hotstart heads from .rvc"<<endl;}
  if (_dryStart!=""){S<<"  dry at the start      : layer "<<_dryStart<<"; every active convertible cell at or below its bottom"<<endl;}
  for (size_t p=0;p<_aProfiles.size();p++){
    int n=0; for (int ic=0;ic<_ncpl;ic++){if (_colProfile[ic]==(int)p){n++;}}
    S<<"  profile "<<_aProfiles[p].name<<": ";
    if (IsExternal()){S<<"layers from the model, ";} else {S<<_aProfiles[p].layers.size()<<" layers, ";}
    S<<n<<" columns"<<endl;
  }
  if (SZPhysical()){
    int nOn=0; double dmin=1e300,dmax=0,cmin=1e300,cmax=0;
    for (size_t i=0;i<_aqHRUs.size();i++){int k=_aqHRUs[i]; if (_szOn[k]){nOn++; dmin=min(dmin,_szD[k]); dmax=max(dmax,_szD[k]);}}
    for (size_t e=0;e<_szxC.size();e++){double a=_szColA[_szxCol[e]]; if (a>0){cmin=min(cmin,_szxC[e]/a); cmax=max(cmax,_szxC[e]/a);}}
    S<<"  soil zone             : PHYSICAL, SOIL[0.."<<_szNLayers-1<<"]; "<<nOn<<" HRUs with a soil zone, "<<_aqHRUs.size()-nOn<<" without (seepage drain)"<<endl;
    if (nOn>0){
      S<<"  soil-zone thickness   : "<<dmin<<" to "<<dmax<<" m; model top = land surface - overlap-weighted thickness"<<endl;
      if (_szxC.size()>0){S<<"  soil-zone links       : "<<_szxC.size()<<" (SZX_RAVEN), conductance per area "<<cmin<<" to "<<cmax<<" 1/d; Sy/drainable porosity "<<_szRatioMin<<" to "<<_szRatioMax<<endl;}
      if ((_szxC.size()>0) && (_vzMode>0)){S<<"  unsaturated columns   : "<<VZName()<<", "<<_vzN<<" nodes at most, ";
        if (_vzMode==3){ //LAG: the travel times (eq. lagtm) instead of an accuracy it does not use
          double tmin=1e300,tmax=0.0,cmin=1e300,cmax=0.0;
          for (size_t e=0;e<_vzStrat.size();e++){
            for (size_t m=0;m<_vzStrat[e].mat.size();m++){const double t=gwvz::LagTravelTime(_vzStrat[e].mat[m],_vzPar); tmin=min(tmin,t); tmax=max(tmax,t);}
            for (int ii=_vzOff[e];ii<_vzOff[e+1];ii++){const gwvz::vz_column &c=_vzcCol[ii]; if (!c.active){continue;}
              double tc=0.0; for (size_t i=0;i<c.dz.size();i++){tc+=gwvz::LagTravelTime(_vzStrat[e].mat[c.im[i]],_vzPar)*c.dz[i];}
              cmin=min(cmin,tc); cmax=max(cmax,tc);}
          }
          if (tmax>0.0){S<<"travel time "<<tmin<<" to "<<tmax<<" d/m"<<((_vzPar>0.0)?" (given)":" (derived at 1 mm/d)");}
          if (cmax>0.0){S<<", columns "<<cmin<<" to "<<cmax<<" d at the start";}
          S<<" (no Newton iterations: the substeps of the run are in the run totals)"; //(the counts are 0 at the start)
        }
        else {S<<"accuracy "<<_vzAcc<<"; "<<_vzSub<<" substeps and "<<_vzNewton<<" Newton iterations in all";}
        S
                                              <<endl;
        if (_vzQ0>0.0){S<<"  initial flux          : "<<_vzQ0<<" mm/d (steady start of the columns on a cold start)"<<endl;}
        if (_vzClasses>1){S<<"  depth classes         : up to "<<_vzClasses<<" per soil-zone link (least-spread grouping of its DEM samples), "<<_vzcCol.size()<<" class columns on "<<_szxC.size()
                            <<" links; class saturated fractions sum to the link's within "<<_vzClsDev<<endl;}
        S<<"  unsaturated storage   : "<<setprecision(15)<<_vzStore0<<setprecision(6)<<" m3 at the start (above the columns' reference profile; GWBudget.csv 'unsaturated zone storage')"<<endl;
      }
      S<<"  interface iteration   : up to "<<_szMaxIter<<" updates after MODFLOW has converged, tolerances "<<_szTolH<<" m and "<<_szTolQ<<"; saturated fraction and contact from "<<_szSubgrid<<"x"<<_szSubgrid<<" DEM samples per cell"<<(_dem.IsLoaded()?"":" (no DEM: cell mean)")<<endl;
      if ((_szxC.size()>0) && _mf6Newton){
        int nk=0; for (size_t i=0;i<_cellConv.size();i++){if (!_cellConv[i] && (_idomain[i]>0)){nk++;}}
        S<<"  drained columns       : kept 1 mm above their bottoms in MODFLOW's iterates (NEWTON; columns whose active cells are all convertible; "<<nk<<" confined cells)"<<endl;
      }
      S<<"  soil-zone columns     : "<<_szColA.size()<<" (one per link, plus the rest of each HRU); at the start "<<_szInitHot<<" HRUs from the hotstart file, "
       <<_szInitVol<<" m3 added to saturate the soil zone below the initial water table"<<endl;
    }
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
  AppendRunTotals();
  if (_BUDGET.is_open()){_BUDGET.close();}
  if (_HRUSTATE.is_open()){_HRUSTATE.close();}
  if (_CONN.is_open()){_CONN.close();}
  if (_HEADS.is_open()){_HEADS.close();}
  if (_SZOUT.is_open()){_SZOUT.close();}
  if (_szNotConv>0){
    cout<<"WARNING: the soil-zone interface did not converge on "<<_szNotConv<<" time steps (see GWBudget.csv, 'interface converged')"<<endl;
    _szNotConv=0;
  }
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

//////////////////////////////////////////////////////////////////
/// \brief run totals appended to GWModelSummary.txt at the end of the run (once)
//
void CGroundwaterModel::AppendRunTotals()
{
  if (_summaryFile.empty() || (_stepCount<=0)){return;}
  ofstream S(_summaryFile.c_str(),ios::app);
  if (!S.good()){return;}
  S<<endl<<"Run totals (end of run)"<<endl;
  S<<"  time steps            : "<<_stepCount<<"; MODFLOW 6 unconverged on "<<_nNonConverged<<endl;
  if (_szxLink.size()>0){S<<"  soil-zone interface   : unconverged on "<<_szNotConv<<" steps"<<endl;}
  if ((_szxLink.size()>0) && (_vzMode>0)){S<<"  unsaturated columns   : "<<_vzSub<<" substeps"<<((_vzMode==3)?"":(" and "+to_string(_vzNewton)+" Newton iterations"))<<" in all"<<endl;}
  if (_szxLink.size()>0){S<<"  drained columns       : kept 1 mm above their bottoms in MODFLOW's iterates in "<<_wetKeptIt<<" outer iterations of "<<_wetKeptSteps<<" steps"<<endl;}
  if ((_szxLink.size()>0) && (_szMaxIter>0) && (SZFreezeAt()>=1)){S<<"  supply-limited updates: "<<_szCapN<<" link updates; largest take at the lowest head beyond the column's free water "<<_szCapExcess<<" m3/d (exchange line and rejection term)"<<endl;}
  if ((_szxLink.size()>0) && (_szMaxIter>0)){S<<"  connection-test drops : "<<_szDropN<<" link updates on "<<_szDropSteps<<" steps (connection head further than d_s from MODFLOW's head: bracket and last residual dropped); next update plain on "<<_szDropPlain<<", the rest on the step's last update"<<endl;}
  if (_ncWithheldN>0){S<<"  discharge withheld    : "<<_ncWithheldVol<<" m3 on "<<_ncWithheldN<<" unconverged steps (MODFLOW's negative residual, withheld from the discharge to Raven in proportion: seepage, river and surface-store gains; it stays in the aquifer, booked as MF6 balance error in GWBudget.csv)"<<endl;}
  if (_negCutVol>0.0){S<<"  negative recharge cut : "<<_negCutVol<<" m3 the top cells could not supply, taken back out of the soil zone (negative rejected recharge in GWBudget.csv)"<<endl;}
  if (_negUnheldVol>0.0){S<<"  negative recharge not held: "<<_negUnheldVol<<" m3 of that cut the soil zone did not hold when it was taken back (the HRU's soil zone held less, at the end of the step, than the cut: a process after the import moved it on before the coupled exchange): no store supplied it ('negative recharge the soil zone did not hold' in GWBudget.csv)"<<endl;}
  if (_vzHotUnusedN>0){S<<"  hotstart columns unused: "<<_vzHotUnusedN<<" unsaturated columns of the hotstart file, "<<_vzHotUnusedV<<" m3 above their reference profiles, dropped (this run has no unsaturated columns)"<<endl;}
  if (IsExternal() && (_vzMode>0) && SZPhysical()){S<<"  columns over specified heads: "<<_vzChdN<<" unsaturated columns ended on cells that were specified heads, their water ("<<_vzChdV<<" m3 above their reference profiles) to the specified heads (recharge onto CHD cells in GWBudget.csv)"<<endl;}
  bool series=false; for (size_t b=0;b<_aBounds.size();b++){if (_aBounds[b].pHead!=NULL){series=true;}}
  if (series){S<<"  head boundaries off   : "<<_bndOffN<<" cell-steps (a :BoundaryHead series at or below the cell bottom: GHB conductance 0, CHD cell not held)"<<endl;}
  if (_resNode.size()>0){S<<"  reservoir gain not supplied: "<<_resShortTot<<" m3 (lake gains MODFLOW did not supply: cells under the lakes that were drying, specified heads or inactive; taken back from the reservoirs' subbasins; GWBudget.csv)"<<endl;}
  if (IsExternal() && (_resNode.size()>0)){S<<"  reservoir seepage not delivered: "<<_resUndelTot<<" m3 (cells under the lakes that were specified heads or inactive: returned to the reservoirs' subbasins; GWBudget.csv)"<<endl;}
  if (IsExternal() && !_mf6Newton){S<<"  recharge not applied  : "<<_rchLostVol<<" m3 (columns whose cells were all dry or whose highest active cell is a specified head, a model without NEWTON; GWBudget.csv; net: a negative recharge there, taken out of GROUNDWATER without MODFLOW supplying it, counts negative)"<<endl;}
  _summaryFile="";
}

//////////////////////////////////////////////////////////////////
/// \brief divergence margin [m]: 10 times the model's relief, at least 1 km; zmin, zmax the lowest bottom and highest top
//
double CGroundwaterModel::FarMargin(double &zmin,double &zmax) const
{
  zmax=-1e300; zmin=1e300;
  for (int ic=0;ic<_ncpl;ic++){for (int l=0;l<_nlay;l++){const int nd=UserNode(l,ic); if (_idomain[nd]>0){zmax=max(zmax,_top[ic]); zmin=min(zmin,_botm[nd]);}}}
  if (zmax<zmin){zmax=zmin=0.0;}
  return max(10.0*(zmax-zmin),1000.0);
}
double CGroundwaterModel::FarDistance(const double x,const double zmin,const double zmax)
{
  if (x>zmax){return x-zmax;}
  if (x<zmin){return zmin-x;}
  return 0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief drained columns (NEWTON, before each solve): a column is drained when every active cell is convertible (ICELLTYPE
///  != 0) and its head is at or below the cell's bottom; such a column's lowest cell is set to its bottom + 1 mm in MODFLOW's
///  iterate (only the iterate is changed: a column with no wet cell leaves MODFLOW's Newton matrix singular). Once per step,
///  and (rule b) once more on any later solve of that step that leaves the lowest cell more than 1 m below its bottom, and
///  after a singular solve (heads not finite); never in the last two outer iterations except for heads that are not finite.
///  A column with an active confined cell is never drained (MODFLOW keeps its full transmissivity at any head; its dry convertible cells are left as they are). Only
///  non-finite heads or the lowest cell are tested. Returns the number of columns set
//
int CGroundwaterModel::KeepColumnsWet(const int it)
{
  if (_X==NULL){return 0;}
  if ((int)_keptCol.size()!=_ncpl){_keptCol.assign(_ncpl,0);}
  const bool late=(it>_maxiter-2);
  int nset=0;
  for (int ic=0;ic<_ncpl;ic++){
    if (_colProfile.size() && !IsExternal() && (_colProfile[ic]==DOESNT_EXIST)){continue;}
    int low=-1; bool conf=false,nonfin=false;
    for (int l=0;l<_nlay;l++){
      const int n=UserNode(l,ic),r=_reduced[n]; if ((r<0) || (_idomain[n]<=0)){continue;}
      if ((n<(int)_cellConv.size()) && !_cellConv[n]){conf=true; break;}
      if (!std::isfinite(_X[r])){nonfin=true;}
      low=n;
    }
    if (conf || (low<0)){continue;}
    const int r=_reduced[low]; const double x=_X[r]*_hf;
    if (!nonfin){
      if (!(x<=_botm[low])){continue;}                                           //(the lowest cell holds water: not drained)
      if (late){continue;}
      if ((_keptCol[ic]>=2) || ((_keptCol[ic]==1) && !(x<_botm[low]-1.0))){continue;}
    }
    if (nonfin){for (int l=0;l<_nlay;l++){const int n=UserNode(l,ic),rr=_reduced[n]; if ((rr>=0) && !std::isfinite(_X[rr])){_X[rr]=(_botm[n]+1e-3)/_hf;}}}
    _X[r]=(_botm[low]+1e-3)/_hf; _keptCol[ic]++; nset++;
  }
  return nset;
}
