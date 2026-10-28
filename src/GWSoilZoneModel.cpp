/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Soil zone <-> groundwater interface of CGroundwaterModel (:SoilZone PHYSICAL).
  Design, notation and equations: docs/design/soil_zone_coupling_plan.md (RavenHydroFramework-extra).

  Raven owns the soil zone SOIL[0..M-1] of every coupled HRU with soil; MODFLOW owns everything below its base. The
  soil zone of an HRU is split into columns, one per HRU-cell link (plus the rest of the HRU); each column keeps its
  own water and Raven's SOIL[m] is their area-weighted mean (as GSFLOW keeps a gravity reservoir per HRU-cell
  intersection). Each link has one GHB entry (package SZX_RAVEN) on the top active cell; with one DEM sample and the
  water table above the base it exchanges
     q = C^sz(k,i) (z_base(k,i) + lambda(k,i) - h)   into the aquifer
  (lambda: saturated thickness of the link's column; GWSoilZone.h has the sub-grid form, with contact through ds below
  the base). Recharge on saturated ground is rejected (returned to the column); column overflow goes to
  SURFACE_WATER the same step. lambda is the backward-Euler local solve of the column for the current head
  (GWSoilZone.cpp), updated inside MODFLOW's outer iterations, and MODFLOW sees the exchange as a Robin condition
  with the column's true slope, until heads and soil levels agree.
----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include <cstdlib>
#include "Model.h"
#include "HydroProcessABC.h"
#include "GWUtil.h"
#include <cmath>
#include <iomanip>
#include <algorithm>

//////////////////////////////////////////////////////////////////
/// \brief soil-zone layers of every coupled HRU, thickness, and the checks of the plan (section 3)
//
void CGroundwaterModel::SZSetup()
{
  int nHRUs=_pModel->GetNumHRUs(),nsoil=_pModel->GetNumSoilLayers();
  ExitGracefullyIf(_szNLayers>nsoil,("CGroundwaterModel: :SoilZone PHYSICAL lists "+to_string(_szNLayers)+" soil layers but the model has "+to_string(nsoil)).c_str(),BAD_DATA);
  _szOn.assign(nHRUs,0); _szSV.assign(nHRUs,vector<int>()); _szLay.assign(nHRUs,vector<gwsz::sz_layer>());
  _szD.assign(nHRUs,0.0);
  int nOn=0,nThick=0; double Dmax=0;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; const CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
    string id=to_string(pHRU->GetHRUID());
    double th=0; for (int m=0;m<_szNLayers;m++){th+=pHRU->GetSoilThickness(m);}           //[mm]
    HRU_type ty=pHRU->GetHRUType();
    if ((th<1e-3) || (ty==HRU_LAKE) || (ty==HRU_WATER) || (ty==HRU_GLACIER) || (ty==HRU_ROCK)){continue;} //no soil (Raven gives
                                                       //profiles without horizons a vanishing thickness): keeps the seepage drain
    for (int j=_szNLayers-1;j>=0;j--) //bottom (m*) first
    {
      double H=pHRU->GetSoilThickness(j)/MM_PER_METER;                 //[m]
      ExitGracefullyIf(!(H>0.0),("CGroundwaterModel: HRU "+id+": SOIL["+to_string(j)+"] has zero thickness; with :SoilZone PHYSICAL every coupled HRU with soil needs all soil-zone layers").c_str(),BAD_DATA);
      gwsz::sz_layer L;
      L.H=H;
      L.phimax=pHRU->GetSoilCapacity(j);                                //[mm] H n (1-SF)
      L.phifc =L.phimax*pHRU->GetSoilProps(j)->field_capacity;          //[mm]
      ExitGracefullyIf(!(L.phimax-L.phifc>1e-9),("CGroundwaterModel: HRU "+id+": SOIL["+to_string(j)+"] has no drainable porosity (POROSITY 0 or FIELD_CAPACITY 1); the water table cannot rise into it").c_str(),BAD_DATA);
      _szLay[k].push_back(L);
      _szSV[k].push_back(_pModel->GetStateVarIndex(SOIL,j));
      _szD[k]+=H;
    }
    map<const void*,string>::const_iterator nk=_szNoKsat.find((const void*)pHRU->GetSoilProps(_szNLayers-1));
    if (nk!=_szNoKsat.end()){
      ExitGracefully(("CGroundwaterModel: HRU "+id+": soil class "+nk->second+" of SOIL["+to_string(_szNLayers-1)+
        "] (the bottom of the soil zone) needs HYDRAUL_COND in the .rvp file with :SoilZone PHYSICAL (it sets the soil-aquifer conductance)").c_str(),BAD_DATA);
    }
    ExitGracefullyIf(_aProfiles[_hruProfile[k]].recharge_delay>0.0,("CGroundwaterModel: HRU "+id+": RECHARGE_DELAY of aquifer profile "+_aProfiles[_hruProfile[k]].name+
      " cannot be used with :SoilZone PHYSICAL (the delay store would hold recharge that the saturated soil zone rejects); set it to 0").c_str(),BAD_DATA);
    double ks=pHRU->GetSoilProps(_szNLayers-1)->hydraul_cond;
    ExitGracefullyIf(!(ks>0.0),("CGroundwaterModel: HRU "+id+": HYDRAUL_COND of SOIL["+to_string(_szNLayers-1)+"] (the bottom of the soil zone) must be positive with :SoilZone PHYSICAL").c_str(),BAD_DATA);
    ExitGracefullyIf(_szD[k]>20.0,("CGroundwaterModel: HRU "+id+": the soil zone SOIL[0.."+to_string(_szNLayers-1)+"] is "+to_string(_szD[k])+
      " m thick; :SoilZone PHYSICAL is meant for physical soil layers (conceptual reservoirs of tens of metres are not a soil zone)").c_str(),BAD_DATA);
    if (_szD[k]>5.0){nThick++;}
    Dmax=max(Dmax,_szD[k]);
    _szOn[k]=1; nOn++;
  }
  if (nThick>0){WriteWarning("CGroundwaterModel: "+to_string(nThick)+" HRUs have a soil zone thicker than 5 m (up to "+to_string(Dmax)+" m); check that SOIL[0.."+to_string(_szNLayers-1)+"] are physical soil layers",false);}
  if ((nOn>0) && (_maxiter<20)){WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL with :SolverMaxOuterIterations "+to_string(_maxiter)+": the soil-zone interface is updated within MODFLOW's outer iterations; a limit below 20 leaves many steps unconverged",false);}
  if (nOn==0){WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: no coupled HRU has soil layers; the setting has no effect",false);}
  for (size_t p=0;p<_aProfiles.size();p++){
    if (_aProfiles[p].soil_zone_depth>0){WriteWarning("CGroundwaterModel: SOIL_ZONE_DEPTH of aquifer profile "+_aProfiles[p].name+" is ignored with :SoilZone PHYSICAL (the soil-zone base comes from the soil layers)",false);}
  }
  if (_seepTarget!=""){WriteWarning("CGroundwaterModel: with :SoilZone PHYSICAL, :SeepageReturnTo applies only to coupled HRUs without soil (lakes, rock, glaciers)",false);}
  if (_szOutGroup!=""){ExitGracefullyIf(_pModel->GetHRUGroup(_szOutGroup)==NULL,("CGroundwaterModel: :SoilZoneOutput HRU group "+_szOutGroup+" does not exist").c_str(),BAD_DATA);}
  //work arrays
  _szPhiOff.assign(nHRUs,-1); int off=0;
  for (int k=0;k<nHRUs;k++){if (_szOn[k]){_szPhiOff[k]=off; off+=(int)_szSV[k].size();}}
  _szPhiProc.assign(off,0.0);
  _szW.assign(nHRUs,0.0); _szLam.assign(nHRUs,0.0); _szRchV.assign(nHRUs,0.0); _szRej.assign(nHRUs,0.0); _szOvfV.assign(nHRUs,0.0); _szNetV.assign(nHRUs,0.0); _szSyncOvf.assign(nHRUs,0.0);
}

double CGroundwaterModel::SZBase(const int k,const int ic) const
{
  if (IsExternal()){return _top[ic];}                                            //an existing model: the soil zone lies on its top
  return _land[ic]-_szD[k];
}
double CGroundwaterModel::SZCellKv(const int ic) const
{
  if (IsExternal()){return _extKv[TopNode(ic)];}
  return _aClasses[_aProfiles[_colProfile[ic]].layers[0].iclass].Kv;
}
double CGroundwaterModel::SZCellSy(const int ic) const
{
  if (IsExternal()){return _extSy[TopNode(ic)];}
  return _aClasses[_aProfiles[_colProfile[ic]].layers[0].iclass].Sy;
}

double CGroundwaterModel::SZNodeSy(const int n) const
{
  if (IsExternal()){return _extSy[n];}
  const int ic=n%_ncpl,l=n/_ncpl; const gw_profile &P=_aProfiles[_colProfile[ic]];
  if (l<(int)P.layers.size()){return _aClasses[P.layers[l].iclass].Sy;}
  return SZCellSy(ic);
}

//////////////////////////////////////////////////////////////////
/// \brief the cell each SZX entry's upward flows come from this step (heads of the last step): its top cell while that
///  holds the water table more than max(5 % of its thickness, 5 cm) above its bottom; else the highest cell below that
///  does (through SZW_RAVEN; MODFLOW cannot supply a rate out of a cell it has dried); when none does, the highest cell
///  with any water, the rates scaled down smoothly to 0 as it dries (_szxUpF); -1 when no cell has water
//
void CGroundwaterModel::SZUpNodes()
{
  for (size_t e=0;e<_szxLink.size();e++){
    const int ic=_aLinks[_szxLink[e]].ic,top=TopNode(ic);
    int up=-1,wet=-1; double fw=0.0;                                             //(wet: the highest cell with any water, its scale)
    for (int j=(int)_szwStart[e]-1;(j<_szwStart[e+1]) && (up<0);j++){             //(j = start-1: the top cell; then the cells below it)
      const bool isTop=(j<_szwStart[e]);
      if (!isTop && (_szwQ==NULL)){break;}
      const int n=isTop?top:_szwNode[j],nu=isTop?-1:n-_ncpl;                    //(its top: the bottom of the cell directly above, active or not)
      const double h=_X[_reduced[n]]*_hf,b=_botm[n],t=isTop?_top[ic]:_botm[nu],mg=max(0.05*(t-b),0.05);
      if (h>=b+mg){up=n; break;}
      if ((h>b) && (wet<0)){wet=n; const double x=(h-b)/mg; fw=x*x*(3.0-2.0*x);} //(smooth: 0 at the bottom, 1 at the margin)
    }
    if (up>=0){_szxUp[e]=up; _szxUpF[e]=1.0;}
    else if ((wet>=0) && ((wet==top) || (_szwQ!=NULL))){_szxUp[e]=wet; _szxUpF[e]=fw;} //(a cell near drying: the rates fade out with
    else {_szxUp[e]=-1; _szxUpF[e]=0.0;}                                         // its water, no switch)
    _szxUpCHD[e]=(_szxUp[e]>=0) && (_szxUp[e]!=top) && NodeCHD(_szxUp[e]);       //(a specified head below: MODFLOW ignores a well there)
    if (_szxUpCHD[e]){_szxUpF[e]=1.0;}
  }
}
double CGroundwaterModel::SZUpCap(const int e,const double &tstep)
{ //half the water the source cell holds above its bottom under the link, faded out as it nears drying; a specified head: no limit
  const int up=_szxUp[e]; if (up<0){return 0.0;}
  if (_szxUpCHD[e]){return 1e300;}
  const double hu=_X[_reduced[up]]*_hf,a=_szColA[_szxCol[e]];
  return _szxUpF[e]*0.5*max(SZNodeSy(up),0.01)*a*max(hu-_botm[up],0.0)/tstep;
}
bool CGroundwaterModel::NodeCHD(const int n)
{
  if (_ibPtr==NULL){string a=_mf6.GetVarAddress(_mname,"","IBOUND"); if (_mf6.HasVar(a)){_ibPtr=_mf6.GetIntPtr(a);}}
  if (_ibPtr!=NULL){const int r=_reduced[n]; return (r>=0) && (_ibPtr[r]<0);}
  return std::find(_chdNode.begin(),_chdNode.end(),n)!=_chdNode.end();
}

//////////////////////////////////////////////////////////////////
/// \brief one SZX entry per HRU-cell link of a soil-zone HRU (top active cell), conductances (half the bottom soil
///  layer and half the top aquifer cell in series), DEM samples per cell; links on specified-head cells get a CHD entry
///  instead (their column exchanges with the fixed head, no MODFLOW boundary)
//
void CGroundwaterModel::SZBuildEntries()
{
  int nHRUs=_pModel->GetNumHRUs();
  _linkF.assign(_aLinks.size(),0.0);
  _colCoverNS.assign(_ncpl,0.0);
  _szxLink.clear(); _szxC.clear(); _szxZb.clear(); _szxQ.clear(); _szxCell.clear(); _szxOfHRU.assign(nHRUs,vector<int>());
  _szcLink.clear(); _szcCell.clear(); _szcOfHRU.assign(nHRUs,vector<int>());
  _cellDem.assign(_ncpl,vector<double>());
  int nRatio=0,nMix=0; double rmin=1e300,rmax=0,dmix=0.0;
  for (size_t li=0;li<_aLinks.size();li++)
  {
    const gw_link &L=_aLinks[li];
    if (!_szOn[L.k]){_colCoverNS[L.ic]+=L.area; continue;}
    const bool chd=_colCHDTop[L.ic];
    int ic=L.ic;
    const gw_profile &P=_aProfiles[_colProfile[ic]];
    const double Kv1=SZCellKv(ic),Sy1=SZCellSy(ic);
    ExitGracefullyIf(!(Kv1>0.0),("CGroundwaterModel: :SoilZone PHYSICAL: the top active cell of grid column "+to_string(ic+1)+" has no vertical conductivity").c_str(),BAD_DATA);
    double b1=_top[ic]-_botm[TopNode(ic)];                                         //[m] top cell
    const gwsz::sz_layer &Lb=_szLay[L.k][0];                                      //bottom soil-zone layer
    double ks=_pModel->GetHydroUnit(L.k)->GetSoilProps(_szNLayers-1)->hydraul_cond/MM_PER_METER; //[m/d]
    double res=Lb.H/(2.0*ks)+b1/(2.0*Kv1);                                       //[d]
    double C=L.area/res;                                                          //[m2/d]
    if (!chd){
      _szxOfHRU[L.k].push_back((int)_szxLink.size());
      _szxLink.push_back((int)li); _szxC.push_back(C); _szxZb.push_back(SZBase(L.k,ic)); _szxQ.push_back(0.0);
    }
    //storage consistency: specific yield of the top aquifer layer vs drainable porosity of the bottom soil layer
    double Sd=(Lb.phimax-Lb.phifc)/(Lb.H*MM_PER_METER);
    double r=Sy1/Sd; rmin=min(rmin,r); rmax=max(rmax,r);
    if ((r<0.5) || (r>2.0)){nRatio++;}
    //the model top is the overlap-weighted soil-zone base of the cell: an HRU whose soil zone is thicker or thinner
    //than the cell's mean (cells shared with lakes or other soil profiles) has its base below or above the top
    double dm=fabs(_szD[L.k]-(_land[ic]-_top[ic]));
    if (!IsExternal() && (dm>max(0.05,0.1*_szD[L.k]))){nMix++; dmix=max(dmix,dm);}
    //DEM sample points of the cell (sub-grid saturated fraction)
    if (_cellDem[ic].empty()){
      if (_dem.IsLoaded() && (_szSubgrid>1) && !(IsExternal() && !_owTable.empty())){ //(an :OverlapWeights table has no cell outlines)
        vector<gw_pt> pts; _grid.SamplePoints(ic,_szSubgrid,pts);
        for (size_t s=0;s<pts.size();s++){
          double X,Y,lon,lat,v; _grid.ToWorld(pts[s].x,pts[s].y,X,Y); _proj.Inverse(X,Y,lon,lat);
          if (_dem.Sample(lon,lat,v)){_cellDem[ic].push_back(v);}
        }
      }
      if (_cellDem[ic].empty()){_cellDem[ic].push_back(_land[ic]);}
    }
    gwsz::sz_cell cell;                                                           //sub-grid closure of the column (GWSoilZone.h)
    if (IsExternal()){ //an existing model: the base is the model top, its relief within the cell from the DEM (about the cell's mean)
      double mz=0.0; for (size_t q=0;q<_cellDem[ic].size();q++){mz+=_cellDem[ic][q];} mz/=max((double)_cellDem[ic].size(),1.0);
      for (size_t q=0;q<_cellDem[ic].size();q++){cell.zbs.push_back(_top[ic]+(_cellDem[ic][q]-mz));}
    }
    else {for (size_t q=0;q<_cellDem[ic].size();q++){cell.zbs.push_back(_cellDem[ic][q]-_szD[L.k]);}}
    std::sort(cell.zbs.begin(),cell.zbs.end());
    cell.ds=P.seep_smooth_depth; cell.D=_szD[L.k]; cell.C=C;
    if (chd){_szcOfHRU[L.k].push_back((int)_szcLink.size()); _szcLink.push_back((int)li); _szcCell.push_back(cell);}
    else    {_szxCell.push_back(cell);}
  }
  if (nMix>0){
    WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: on "+to_string(nMix)+" links the soil-zone base of the HRU lies up to "+to_string(dmix)+
                 " m from the model top (the cell is shared with HRUs without soil or with a different soil-zone thickness); between the base "
                 "and the top the soil zone and the aquifer both store water (base below the top) or neither does (base above it)",false);
  }
  if (nRatio>0){
    WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: on "+to_string(nRatio)+" links the specific yield of the top aquifer layer differs from the drainable porosity "
                 "n(1-SF)(1-FIELD_CAPACITY) of the bottom soil-zone layer by more than a factor 2 (ratio "+to_string(rmin)+" to "+to_string(rmax)+
                 "); the water table then rises at a different rate in the soil than below it",false);
  }
  _szRatioMin=(rmin<1e299)?rmin:0.0; _szRatioMax=rmax;
  //soil-zone columns (design note section 4): the soil zone of an HRU is split into one column per SZX entry and the rest
  //of the HRU (cells with specified heads, parts outside the grid). A column's area is its share of the HRU polygon
  //times the HRU's AREA (the polygon and the .rvh AREA rarely agree exactly). Each column keeps its own water, as
  //GSFLOW keeps a gravity reservoir per HRU-cell intersection; Raven's SOIL[m] is their area-weighted mean.
  _linkSzx.assign(_aLinks.size(),-1);
  for (size_t e=0;e<_szxLink.size();e++){_linkSzx[_szxLink[e]]=(int)e;}
  _linkSzc.assign(_aLinks.size(),-1);
  for (size_t e=0;e<_szcLink.size();e++){_linkSzc[_szcLink[e]]=(int)e;}
  _szColOff.assign(nHRUs,-1); _szColA.clear(); _szxCol.assign(_szxLink.size(),-1); _szcCol.assign(_szcLink.size(),-1);
  for (int k=0;k<nHRUs;k++){
    if (!_szOn[k]){continue;}
    double A=_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2,S=0.0,sum=0.0;    //[m2]
    for (size_t j=0;j<_hruLinks[k].size();j++){S+=_aLinks[_hruLinks[k][j]].area;}
    double P=max(((int)_hruPolyArea.size()>k)?_hruPolyArea[k]:0.0,S);         //[m2] polygon area (at least the linked part)
    double sc=(P>0.0)?A/P:1.0;
    _szColOff[k]=(int)_szColA.size();
    for (size_t j=0;j<_szxOfHRU[k].size();j++){
      int e=_szxOfHRU[k][j]; double a=_aLinks[_szxLink[e]].area*sc;
      _szxCol[e]=(int)_szColA.size(); _szColA.push_back(a); sum+=a;
    }
    for (size_t j=0;j<_szcOfHRU[k].size();j++){                                   //columns over specified heads
      int e=_szcOfHRU[k][j]; double a=_aLinks[_szcLink[e]].area*sc;
      _szcCol[e]=(int)_szColA.size(); _szColA.push_back(a); sum+=a;
    }
    _szColA.push_back(max(A-sum,0.0));                                            //rest of the HRU (outside the grid)
  }
  //the conductance of a link acts on its column: on the HRU's AREA share, as the column's water does (with the
  //overlap area instead, an HRU whose polygon is larger than its AREA would exchange its soil water that much faster)
  for (size_t e=0;e<_szxLink.size();e++){
    double la=_aLinks[_szxLink[e]].area,sc=(la>0.0)?_szColA[_szxCol[e]]/la:1.0;
    _szxC[e]*=sc; _szxCell[e].C*=sc;
  }
  for (size_t e=0;e<_szcLink.size();e++){
    double la=_aLinks[_szcLink[e]].area,sc=(la>0.0)?_szColA[_szcCol[e]]/la:1.0;
    _szcCell[e].C*=sc;
  }
  _szcRej.assign(_szcLink.size(),0.0); _szcFw.assign(_szcLink.size(),0.0); _szcQ.assign(_szcLink.size(),0.0); _szcCapL.assign(_szcLink.size()*_szNLayers,0.0);
  _szxCapQ.assign(_szxLink.size(),0.0); _szxCapL.assign(_szxLink.size()*_szNLayers,0.0); //(capillary rise, zero unless it is on)
  //SZW_RAVEN: upward flows (capillary rise, the unsaturated column's upflow) come from the cell holding the water table
  //when the top cell has (nearly) dried: one rate per active cell below each entry's top cell
  _szwStart.assign(_szxLink.size()+1,0); _szwNode.clear();
  for (size_t e=0;e<_szxLink.size();e++){
    _szwStart[e]=(int)_szwNode.size();
    const int ic=_aLinks[_szxLink[e]].ic,top=TopNode(ic);
    for (int l=top/_ncpl+1;l<_nlay;l++){const int n=UserNode(l,ic); if (_idomain[n]>0){_szwNode.push_back(n);}}
  }
  _szwStart[_szxLink.size()]=(int)_szwNode.size();
  _szxUp.assign(_szxLink.size(),-1); _szxUpF.assign(_szxLink.size(),0.0); _szxUpCHD.assign(_szxLink.size(),0); _vzUp.assign(_szxLink.size(),0.0); _vzUpA.assign(_szxLink.size(),0.0); _vzUpQ=0.0;
  { double dsmin=1e300; for (size_t e=0;e<_szxCell.size();e++){dsmin=min(dsmin,_szxCell[e].ds);}
    if ((_szxCell.size()>0) && (dsmin<0.1)){
      WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: SEEPAGE_SMOOTHING_DEPTH of the cells' profiles is as small as "+to_string(dsmin)+
                   " m; the soil-zone contact then switches almost abruptly and the interface may not converge (0.1 to 0.5 m is typical)",false);}
  }
  _szColPhi.assign(_szColA.size()*_szNLayers,0.0); _szColInit=false;
  _szxLam.assign(_szxLink.size(),0.0); _szxRej.assign(_szxLink.size(),0.0); _szxFw.assign(_szxLink.size(),0.0); _szxFc.assign(_szxLink.size(),0.0); _szxOne.assign(_szxLink.size(),0);
  _szxQc.assign(_szxLink.size(),0.0); _szxCg.assign(_szxLink.size(),0.0);
  size_t ne=_szxLink.size();
  _szxHc.assign(ne,0.0); _szxHp.assign(ne,0.0); _szxHr.assign(ne,-HUGE_VAL); _szxHq.assign(ne,-HUGE_VAL); _szxQq.assign(ne,0.0);
  _szxBa.assign(ne,0.0); _szxBra.assign(ne,0.0); _szxBb.assign(ne,0.0); _szxBrb.assign(ne,0.0); _szxBk.assign(ne,0);
  if (_vzMode>0){VZSetup();}
}

//////////////////////////////////////////////////////////////////
/// \brief soil-zone columns at the start of the run (after the initial heads are final, before Raven takes its
///  initial storage): from the hotstart file where it has them; otherwise from Raven's initial soil water, each column
///  then filled (never drained) until its soil water level reaches the initial water table of its cell (lambda =
///  gwsz::SubgridThickness(h0)), the state the exchange would otherwise reach in its first steps by pulling water out
///  of the aquifer. Raven's SOIL[m] is set to the columns' mean.
//
void CGroundwaterModel::SZInitColumns()
{
  const int n=_szNLayers;
  _szInitVol=0.0; _szInitHot=0; int nWarn=0;
  vector<double> tg(n),work;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
    double A=pHRU->GetArea()*M2_PER_KM2;
    int nc=SZNCols(k),o=_szColOff[k];
    for (int j=0;j<n;j++){tg[j]=pHRU->GetStateVarValue(_szSV[k][j]);}
    vector<vector<double> > hc;
    int hot=SZHotColumns(k,hc);
    if ((hot<0) && (++nWarn<=10)){
      WriteWarning("CGroundwaterModel: :GWSoilZoneColumns of HRU "+to_string(pHRU->GetHRUID())+" do not match its soil-zone columns ("+to_string(nc)+
                   " columns of "+to_string(n)+" layers; the grid or :SoilZone changed?): the columns start from the water table instead",false);
    }
    if (hot>0){
      for (int c=0;c<nc;c++){for (int j=0;j<n;j++){_szColPhi[(o+c)*n+j]=hc[c][j];}}
      _szInitHot++;
      _szColInit=true; double ovf=SZSyncColumns(k,&tg[0]); //hotstart columns and the .rvc soil water normally agree; any difference is shared
      double sw=0.0; for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
      for (int j=0;j<n;j++){
        double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*n+j];}
        pHRU->SetStateVarValue(_szSV[k][j],(sw>0.0)?v/sw:tg[j]);
      }
      if (ovf>0.0){int iSW=_pModel->GetStateVarIndex(SURFACE_WATER); pHRU->SetStateVarValue(iSW,pHRU->GetStateVarValue(iSW)+ovf);}
      if (ovf>1e-6){WriteWarning("CGroundwaterModel: HRU "+to_string(pHRU->GetHRUID())+": the .rvc soil water exceeds what its :GWSoilZoneColumns hold by "+gwutil::ToStr(ovf)+
                                 " mm where the columns are full; that water starts in SURFACE_WATER",false);}
      continue;
    }
    for (int c=0;c<nc;c++){for (int j=0;j<n;j++){_szColPhi[(o+c)*n+j]=tg[j];}}
    double added=0.0;                                                           //[m3]
    const int nx=(int)_szxOfHRU[k].size(),nch=(int)_szcOfHRU[k].size();
    for (int jj=0;jj<nx+nch;jj++){                                              //(exchange columns, then those over specified heads)
      const bool ch=(jj>=nx); const int e=ch?_szcOfHRU[k][jj-nx]:_szxOfHRU[k][jj],c=ch?_szcCol[e]:_szxCol[e];
      const gw_link &L=_aLinks[ch?_szcLink[e]:_szxLink[e]];
      double h0=_X[_reduced[TopNode(L.ic)]]*_hf;                                //[m]
      const gwsz::sz_cell &cell=ch?_szcCell[e]:_szxCell[e];
      if (!ch){_szxHc[e]=h0;}
      if (gwsz::SatFraction(h0,cell.zbs,cell.ds)<=0.0){continue;}               //not connected: the soil keeps Raven's water
      double lamT=gwsz::SubgridThickness(cell,h0);                              //[m] saturated thickness at equilibrium
      vector<double> phi(_szColPhi.begin()+c*n,_szColPhi.begin()+(c+1)*n);
      if (gwsz::SatThickness(_szLay[k],phi)>=lamT){continue;}
      double lo=0.0,hi=gwsz::Room(_szLay[k],phi);                               //[mm] water to add: bisection on lambda(Fill)
      for (int itb=0;itb<200 && (hi-lo>1e-10*max(1.0,hi));itb++){
        double mid=0.5*(lo+hi); work=phi; gwsz::Fill(_szLay[k],work,mid);
        if (gwsz::SatThickness(_szLay[k],work)<lamT){lo=mid;} else {hi=mid;}
      }
      gwsz::Fill(_szLay[k],phi,hi);
      for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
      added+=hi/MM_PER_METER*_szColA[c];
    }
    //Raven's soil-zone layers: area-weighted mean of the columns
    double sw=0.0; for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
    for (int j=0;j<n;j++){
      double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*n+j];}
      pHRU->SetStateVarValue(_szSV[k][j],(sw>0.0)?v/sw:tg[j]);
    }
    _szInitVol+=added;
    (void)A;
  }
  _szColInit=true;
  _szxQprev.assign(_szxLink.size(),0.0);                                         //LAGGED source of the first step: from the hotstart file
  _szxDebt.assign(_szxLink.size(),0.0); _szxRepay.assign(_szxLink.size(),0.0);
  if (!_hotLag.empty()){
    for (size_t e=0;e<_szxLink.size();e++){
      const gw_link &L=_aLinks[_szxLink[e]];
      map<long long,map<int,double> >::const_iterator it=_hotLag.find(_pModel->GetHydroUnit(L.k)->GetHRUID());
      if (it==_hotLag.end()){continue;}
      map<int,double>::const_iterator jt=it->second.find(L.ic+1);
      if (jt!=it->second.end()){_szxQprev[e]=jt->second; _szxDebt[e]=max(_hotDebt[it->first][L.ic+1],0.0);}
    }
  }
  if (_vzMode>0){VZInit();}
}

//////////////////////////////////////////////////////////////////
/// \brief cell of column c (local) of HRU k: the grid column number (1-based) of the link that owns it, 0 for the rest
///  column (written to and matched in the hotstart file, so columns follow their cells if the links are reordered)
//
vector<int> CGroundwaterModel::SZColumnCells(const int k) const
{
  vector<int> cell(SZNCols(k),0);
  for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
    int e=_szxOfHRU[k][jj];
    cell[_szxCol[e]-_szColOff[k]]=_aLinks[_szxLink[e]].ic+1;
  }
  if (k<(int)_szcOfHRU.size()){
    for (size_t jj=0;jj<_szcOfHRU[k].size();jj++){
      int e=_szcOfHRU[k][jj];
      cell[_szcCol[e]-_szColOff[k]]=_aLinks[_szcLink[e]].ic+1;
    }
  }
  return cell;
}
int CGroundwaterModel::SZHotColumns(const int k,vector<vector<double> > &cols) const
{
  map<long long,map<int,vector<double> > >::const_iterator it=_hotSzCols.find(_pModel->GetHydroUnit(k)->GetHRUID());
  if (it==_hotSzCols.end()){return 0;}
  const int nc=SZNCols(k),n=_szNLayers;
  if ((int)it->second.size()!=nc){return -1;}
  cols.assign(nc,vector<double>());
  const vector<int> cell=SZColumnCells(k);
  for (int c=0;c<nc;c++){
    map<int,vector<double> >::const_iterator jt=it->second.find(cell[c]);
    if ((jt==it->second.end()) || ((int)jt->second.size()!=n)){return -1;}
    cols[c]=jt->second;
  }
  return 1;
}

//////////////////////////////////////////////////////////////////
/// \brief brings the columns of HRU k to its soil-zone water target[j] (layers bottom up) [mm]: Raven's processes
///  changed the HRU totals since the columns were last set (gwsz::Distribute: recharge as the process takes it,
///  percolation within each column, evaporation and rain alike on every column, overflow where full).
///  On the first call the columns start from the hotstart file or, without one, uniform.
//
double CGroundwaterModel::SZSyncColumns(const int k,const double *target,const double R,vector<double> *lossR)
{
  int nc=SZNCols(k),n=_szNLayers,o=_szColOff[k];
  vector<vector<double> > cols(nc,vector<double>(n));
  vector<double> w(nc),tg(target,target+n);
  double sw=0.0; for (int c=0;c<nc;c++){w[c]=_szColA[o+c]; sw+=w[c];}
  for (int c=0;c<nc;c++){w[c]=(sw>0.0)?w[c]/sw:1.0/nc;}
  for (int c=0;c<nc;c++){for (int j=0;j<n;j++){cols[c][j]=_szColPhi[(o+c)*n+j];}}
  if (!_szColInit){
    vector<vector<double> > hc;
    if (SZHotColumns(k,hc)>0){cols=hc;}
    else {for (int c=0;c<nc;c++){cols[c]=tg;}}
  }
  double mv=0.0;
  double ovf=gwsz::Distribute(_szLay[k],cols,w,tg,0,R,(k<(int)_szRchByArea.size()) && _szRchByArea[k],lossR,&mv);
  _szMovedV+=mv/MM_PER_METER*_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
  for (int c=0;c<nc;c++){for (int j=0;j<n;j++){_szColPhi[(o+c)*n+j]=cols[c][j];}}
  return ovf;
}

//////////////////////////////////////////////////////////////////
/// \brief connected (saturated) fraction of every soil-zone link (design note section 5.4): the share of the link's recharge
///  that is rejected, and of its cell's water-table ET that is left to Raven.
///  The connection follows a head hc, the root of r(hc) = h(hc) - hc, h the head MODFLOW returns for that connection.
///  The predictor and the updates after outer iterations on which MODFLOW has not converged take hc = h; after converged
///  ones, each link takes hc = h (a plain update, extrapolated by the secant while r keeps its sign and shrinks) until
///  r changes sign, and
///  from then on regula falsi with the Illinois correction inside the bracket (Dowell & Jarratt 1971): near a soil base
///  the connection can switch between two states that each send the head to the other side of the base, and a plain
///  update then cycles; the bracketed update settles between them. At a converged interface hc = h.
//
void CGroundwaterModel::SZLinkFractions(const bool relax,const bool mf6conv)
{
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    double h=_X[_reduced[TopNode(L.ic)]]*_hf;                                    //[m]
    if (!relax || !mf6conv){_szxHc[e]=h; _szxHr[e]=-HUGE_VAL; _szxBk[e]=0;} //heads MODFLOW has not converged on are no residuals of hc
    else {
      double hc=_szxHc[e],r=h-hc,hn=h;
      if ((_szxHr[e]!=-HUGE_VAL) && (r!=0.0)){
        double hp=_szxHp[e],rp=_szxHr[e];
        if (_szxBk[e]==0){
          if (rp*r<0.0){_szxBa[e]=hp; _szxBra[e]=rp; _szxBb[e]=hc; _szxBrb[e]=r; _szxBk[e]=2;}
          else { //no sign change yet: a shrinking residual (0 < rho < 1, a slow drift) is extrapolated by the secant
                 //through the last two points (Aitken on the plain update); a growing one is left to the plain update,
                 //which soon changes sign (turning back along the secant destabilizes neighbouring links)
            double rho=r/rp;
            if ((rho>0.0) && (rho<0.95)){hn=hc+r/(1.0-rho);}            //|step| <= 20 |r|
          }
        }
        else { //replace the end whose residual has the same sign; Illinois: halve the other end's residual when the
               //same end is replaced twice in a row (_szxBk: the end replaced last, 1 = a, 2 = b)
          if (r*_szxBrb[e]>0.0){_szxBb[e]=hc; _szxBrb[e]=r; if (_szxBk[e]==2){_szxBra[e]*=0.5;} _szxBk[e]=2;}
          else                 {_szxBa[e]=hc; _szxBra[e]=r; if (_szxBk[e]==1){_szxBrb[e]*=0.5;} _szxBk[e]=1;}
        }
        if (_szxBk[e]!=0){
          double a=_szxBa[e],ra=_szxBra[e],b=_szxBb[e],rb=_szxBrb[e];
          hn=(rb!=ra)?b-rb*(b-a)/(rb-ra):0.5*(a+b);
          double lo=min(a,b),hi=max(a,b);
          if (!(hn>lo) || !(hn<hi)){hn=0.5*(a+b);}
          if (hi-lo<1e-9*max(1.0,fabs(hc))){hn=h; _szxBk[e]=0;} //collapsed on a root that moved (other links changed): start again
        }
      }
      _szxHp[e]=hc; _szxHr[e]=r; _szxHc[e]=hn;
    }
    _linkF[_szxLink[e]]=gwsz::SatFraction(_szxHc[e],_szxCell[e].zbs,_szxCell[e].ds);
  }
  for (size_t e=0;e<_szcLink.size();e++){ //over specified heads: the fixed head
    const gw_link &L=_aLinks[_szcLink[e]];
    _linkF[_szcLink[e]]=gwsz::SatFraction(_X[_reduced[TopNode(L.ic)]]*_hf,_szcCell[e].zbs,_szcCell[e].ds);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief recharge to MODFLOW with rejection on saturated ground (design note section 6.1): _Qcol, RCH rates, _szxRej, _szRej
//
void CGroundwaterModel::SZRecharge(const double &tstep,const bool write)
{
  std::fill(_Qcol.begin(),_Qcol.end(),0.0);
  _szRchCHD=0.0;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; double V=_szRchV[k],rej=0.0;
    double Aall=_hruRchArea[k];                                                  //(with a soil zone: the columns over specified heads too)
    if (_szOn[k]){ //(all its links' areas: also those on cells that are specified heads this step)
      Aall=0.0; for (size_t j=0;j<_hruLinks[k].size();j++){Aall+=_aLinks[_hruLinks[k][j]].area;}
    }
    if (Aall<=0){if (_szOn[k]){_szRej[k]=0.0;} continue;} //all on specified-head cells (counted in _QrchCHD)
    //an HRU with a soil zone: recharge leaves the columns that hold free water in the bottom layer (Raven's drainage acts
    //on the HRU mean, and its loss is taken from the columns' free water first, gwsz::Distribute), not by area: drainage
    //from saturated columns spread over unsaturated ones would recharge the aquifer upslope of where it came from
    double W=0.0;
    if (_szOn[k]){
      for (size_t j=0;j<_hruLinks[k].size();j++){int li=_hruLinks[k][j]; if (_linkSzx[li]>=0){W+=_szxFw[_linkSzx[li]];}}
      for (size_t j=0;j<_szcOfHRU[k].size();j++){W+=_szcFw[_szcOfHRU[k][j]];}
    }
    for (size_t j=0;j<_hruLinks[k].size();j++){
      int li=_hruLinks[k][j]; const gw_link &L=_aLinks[li];
      if (_szOn[k] && (_linkSzc[li]>=0)){ //a column over a specified head: its recharge goes to the specified head, rejected where saturated
        const int ec=_linkSzc[li];
        const double share=(W>0.0)?_szcFw[ec]/W:L.area/Aall,F=_linkF[li];
        _szcRej[ec]=V*share*F; rej+=_szcRej[ec]; _szRchCHD+=V/tstep*share*(1.0-F);
        continue;
      }
      if (_szOn[k] && _colCHDTop[L.ic] && (_linkSzx[li]>=0)){ //the same for an exchange column whose cell is a specified head this
        const int e2=_linkSzx[li];                                                //step (existing models)
        const double share=(W>0.0)?_szxFw[e2]/W:L.area/Aall,F=gwsz::SatFraction(_X[_reduced[TopNode(L.ic)]]*_hf,_szxCell[e2].zbs,_szxCell[e2].ds);
        _szxRej[e2]=V*share*F; rej+=_szxRej[e2]; _szRchCHD+=V/tstep*share*(1.0-F);
        continue;
      }
      if (_colCHDTop[L.ic]){continue;}
      double share=(W>0.0)?((_linkSzx[li]>=0)?_szxFw[_linkSzx[li]]/W:0.0):L.area/(_szOn[k]?Aall:_hruRchArea[k]),F=_linkF[li];
      const int e=_linkSzx[li];
      if ((_vzMode>0) && (e>=0) && (e<(int)_vzUse.size())){
        if ((_vzAbsQ[e]<0.0) && (e<(int)_szxUpCHD.size()) && _szxUpCHD[e]){_szRchCHD+=_vzAbsQ[e];} //water a risen water table took from the column; the
        else if (!(_vzUpA[e]>0.0)){_Qcol[L.ic]+=_vzAbsQ[e];}                     // new nodes of a fallen one (out of a specified head below; out of a
                                                                                 // lower cell: SZW_RAVEN)
        if (_vzUse[e]){ //the unsaturated column carried the unsaturated share this step: its bottom flux recharges the
                        //aquifer, and the saturated share is rejected with the saturated fraction its inflow used
          double rj=V*share*_vzFp[e]+_vzRejV[e];
          if ((_vzQ[e]<0.0) && (e<(int)_szxUpCHD.size()) && _szxUpCHD[e]){_szRchCHD+=_vzQ[e];} //(upflow out of a specified head below: supplied by it,
          else {_Qcol[L.ic]+=(_vzUp[e]>0.0)?0.0:_vzQ[e];}                          // MODFLOW does not see it; from a lower cell: SZW_RAVEN)
          rej+=rj; _szxRej[e]=rj;
          continue;
        }
      }
      _Qcol[L.ic]+=V/tstep*share*(1.0-F);
      rej+=V*share*F;
      if (e>=0){_szxRej[e]=V*share*F;}
    }
    if (_szOn[k]){_szRej[k]=rej;}
  }
  _QrchCHD=_QrchCHD0+_szRchCHD;
  if (write){for (size_t i=0;i<_bndCol.size();i++){_rchBound[i]=_Qcol[_bndCol[i]]/MF6CellArea(_bndCol[i])*_rf;}}
}

//////////////////////////////////////////////////////////////////
/// \brief local implicit solve of every exchange column for the current heads (design note sections 5.6, 7.1), and the linear
///  boundary MODFLOW gets for it: _szxLam, _szxQc (exchange at the current head), _szxCg (its slope)
///  (steady: the soil zone is not active yet, lambda = 0: each connected sample drains to its own base)
//
void CGroundwaterModel::SZLocalSolves(const double &tstep,const bool steady)
{
  vector<double> phi;
  double maxd=0.0; const int n=_szNLayers; _szResE=-1;
  for (size_t e=0;e<_szxLink.size();e++)
  {
    const gw_link &L=_aLinks[_szxLink[e]]; int k=L.k,c=_szxCol[e];
    const gwsz::sz_cell &cell=_szxCell[e];
    double h=_X[_reduced[TopNode(L.ic)]]*_hf;                                    //[m] (the soil zone works in metres and m3/d)
    double a=_szColA[c];                                                         //[m2] column area
    double lam=0.0,lamp=0.0;
    if (!steady){
      phi.assign(_szColPhi.begin()+c*n,_szColPhi.begin()+(c+1)*n);
      gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);                       //rejected recharge back to the column
      if ((e<_szxLagD.size()) && (_szxLagD[e]!=0.0)){                            //LAGGED: the part of the source not taken back
        if (_szxLagD[e]>0.0){gwsz::Fill(_szLay[k],phi,_szxLagD[e]);} else {double sh=0.0; gwsz::Apply(_szLay[k],phi,_szxLagD[e],sh);}
      }
      double dV,ovf;
      if (_colCHDTop[L.ic]){ //a specified head this step (existing models): no MODFLOW boundary; the exchange with the head is
                             //found after MODFLOW's solve (SZChdExchange, in the booking)
        _szxLam[e]=gwsz::SatThickness(_szLay[k],phi); _szxQc[e]=0.0; _szxCg[e]=0.0; _szxOne[e]=0; _szxFc[e]=gwsz::SatFraction(h,cell.zbs,cell.ds); _szxHc[e]=h;
        continue;
      }
      lam=gwsz::LocalSolve(_szLay[k],phi,cell,h,_szxHc[e],a,tstep,dV,ovf,lamp);
      if (fabs(lam-_szxLam[e])>maxd){maxd=fabs(lam-_szxLam[e]); _szResE=(int)e;}
    }
    _szxLam[e]=lam;
    //MODFLOW gets q = Cg (Hg - h): q the exchange at this head, Cg = -dq/dh with the column's own response included
    //(implicit function theorem on the local balance dV + c q(lambda(dV),h) = 0, c = 1000 dt/a, the connection hc
    //held; its change is found by the update of hc):
    //   dq/dh = q_h / (1 + c q_lam lambda') = -F C / (1 + c q_lam lambda'),   lambda' = dlambda/dV [m/mm]
    //Cg sets only how fast the interface converges. For water rising into the soil it is kept at least the secant to
    //the level where the contact ends (the lowest base minus ds), below which nothing can rise: the linear boundary then
    //never extracts water from a cell that has fallen below every base, as a too small slope would (a fixed-rate well
    //drying its cell). (A floor taken at the base itself is far too stiff just above it, and the iteration creeps.)
    double F,qh,qc,ql;
    double q=gwsz::SubgridFlux(cell,lam,h,_szxHc[e],F,qh,qc,ql);
    if (_szxHc[e]!=h){ //the connection still lags the head: the exchange it misses, in head units (q/C)
      double F1,a1,a2,a3; double d=fabs(gwsz::SubgridFlux(cell,lam,h,h,F1,a1,a2,a3)-q)/cell.C;
      if (d>maxd){maxd=d; _szResE=(int)e;}
    }
    double x=steady?0.0:MM_PER_METER*tstep/a*ql*lamp;
    double Cg=-qh/(1.0+x),sec=0.0; (void)qc;
    //never below the secant through the last local solve of this step: the exchange as a function of the head steepens
    //where the column fills (lambda' falls to 0, the slope jumps from -q_h/(1 + c q_lam lambda') to -q_h) or a sample
    //connects; a slope taken on the flat side lets MODFLOW overshoot across that kink and cycle (a slope too large only
    //slows the convergence, one too small can diverge)
    if (!steady){
      if ((_szxHq[e]!=-HUGE_VAL) && (fabs(h-_szxHq[e])>1e-6)){
        double s=-(q-_szxQq[e])/(h-_szxHq[e]);
        if (s>Cg){Cg=min(s,max(-qh,Cg));}
      }
      _szxHq[e]=h; _szxQq[e]=q;
    }
    if (q<0.0){sec=-q/max(h-(cell.zbs[0]-max(cell.ds,0.0)),1e-3);} //no contact below the lowest base minus ds
    Cg=max(max(Cg,sec),1e-9*cell.C);
    //explicit mode (no interface iteration) and the last update possible before the interface freezes (iterations or
    //update budget): the boundary MODFLOW keeps may not take more from the column than it holds, whatever the head does.
    //Where the soil drains to the aquifer (q >= 0) the slope is limited so that the flow never exceeds the column's free
    //water over the step down to the cell bottom, as the river and store exchanges are limited. Where water rises into
    //the soil (q < 0) a linear boundary would reverse when the head falls below Hg and drain columns without free water
    //(water created), and a slope small enough to prevent that pumps water up from below the contact
    //: the line is given to MODFLOW as a drain instead (SZSetArrays), which stops at Hg, and the secant
    //floor keeps Hg at or above the end of the contact.
    //(Iterated, the converged flow is the local solve's, which respects both limits.)
    _szxOne[e]=0;
    if (!steady && ((_szMaxIter==0) || (SZFreezeAt()<1) || _szCapUpd)){
      double avail=gwsz::StackFreeWater(_szLay[k],phi)/MM_PER_METER*a/tstep;           //[m3/d]
      double bot=_botm[TopNode(L.ic)];
      double cap=max(avail-q,0.0)/max(h-bot,1e-3);
      if (q>=0.0){Cg=max(min(Cg,cap),1e-9*cell.C);}
      else if (Cg>cap){ //only where the line could take more than the column holds: one-sided, if at least 10 MODFLOW
                        //solves follow; with fewer MODFLOW did not settle on the changed boundary, and the line stays,
                        //as limited as the secant floor allows (a reversal can then still overdraw)
        //(not in a cell dry at the head of linearization: a drain there is inactive at that head, MODFLOW's first Newton
        //iterate sees no boundary and its head runs away)
        if (_szOneOK && (h>bot+1e-3)){_szxOne[e]=1;}
        else {Cg=max(cap,sec);}
      }
    }
    _szxFc[e]=F; _szxQc[e]=q; _szxCg[e]=Cg;
  }
  _szResid=maxd;
}

//////////////////////////////////////////////////////////////////
/// \brief SZX heads and conductances, RCH rates (with rejection) and EVT rates (outside saturated ground)
//
void CGroundwaterModel::SZSetArrays(const bool steady)
{
  //MODFLOW refuses a GHB head below the cell bottom: there the head is raised to the bottom and SZQ_RAVEN (a rate)
  //carries Cg (Hg - bottom), which leaves the boundary q = Cg (Hg - h) unchanged for every h. Hg - h = q/Cg, so at a
  //converged interface the flow MODFLOW books is exactly the flow of the local solve (a Robin interface condition).
  //A boundary MODFLOW keeps where water rises into the soil (_szxOne) is the same line as a drain (SZD_RAVEN): it stops
  //at Hg instead of turning into flow from the soil. Its onset is smoothed (MODFLOW's DDRN) over half the distance to
  //Hg, at most the link's smoothing depth (the flow at the head of linearization is then the line's), and at least a
  //quarter of it (MODFLOW's Newton iteration converges badly on a sharper kink).
  int st=max(_szxNcol,1),sq=max(_szqNcol,1),sd=max(_szdNcol,1);
  vector<char> pureRate(_szxLink.size(),0);
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    double h=_X[_reduced[TopNode(L.ic)]]*_hf;                                    //[m]
    double bot=_botm[TopNode(L.ic)]+1e-3;                                         //[m] margin: MODFLOW reads the bottom rounded
    double Cg=_szxCg[e];
    const bool inQ=(e>=_szxUp.size()) || (_szxUp[e]==TopNode(L.ic));            //(else SZW_RAVEN, out of a lower cell: below)
    const double cap=inQ?(((e<_szxCapQ.size())?_szxCapQ[e]:0.0)+((e<_szxRepay.size())?_szxRepay[e]:0.0))/_qf:0.0; //capillary rise from the water table (SZCapRiseWT) and
                                                                                 //the repayment of the lagged source's debt (SZRepayRates), out of the aquifer
    //(MODFLOW's arrays in the model's units: heads /_hf, conductances *_hf/_qf, rates /_qf; 1 for generated models)
    const double cf=_hf/_qf;
    _szdElev[e*sd]=_top[L.ic]/_hf; _szdCond[e*sd]=0.0; _szdAux[e]=0.0;
    if ((Cg<=1e-6*_szxC[e]) || _colCHDTop[L.ic]){ //exchange (nearly) independent of the head: a rate (no GHB head far outside the
                                                   //cell); on a specified-head cell (existing models) no boundary at all
      _szxHead[e*st]=bot/_hf; _szxCond[e*st]=0.0; _szqQ[e*sq]=_colCHDTop[L.ic]?0.0:_szxQc[e]/_qf-cap;
      if (!_colCHDTop[L.ic]){ //an outlet for a water table that rises past the soil-zone base during MODFLOW's solve: a rate has no
                              //slope, and the excess would go into the elastic storage above the cell top (heads far above ground);
                              //the drain starts at the lowest base sample, where the soil begins to saturate (zero flow below it,
                              //so a solve that stays below the base is unchanged); its conductance below, after the recharge
        const double zl=_szxCell[e].zbs.empty()?_szxZb[e]:_szxCell[e].zbs.front();
        _szdElev[e*sd]=max(zl,bot)/_hf; _szdCond[e*sd]=_szxC[e]*cf; _szdAux[e]=max(_szxCell[e].ds,1e-3)/_hf; pureRate[e]=1;
      }
      continue;
    }
    double Hg=h+_szxQc[e]/Cg;
    double H=max(Hg,bot);
    if (_szxOne[e] && (_szxQc[e]<0.0)){ //(a drain level below the cell bottom: the line is steepened to stop at the bottom,
                                        //not continued by a fixed rate that would pump the dried cell; at most to the link
                                        //conductance, then the flow at h is less than the line's)
      if (Hg<bot){Cg=min(-_szxQc[e]/max(h-bot,1e-6),max(Cg,_szxC[e])); H=bot;}
      _szxHead[e*st]=bot/_hf; _szxCond[e*st]=0.0; _szqQ[e*sq]=-cap;
      _szdElev[e*sd]=H/_hf; _szdCond[e*sd]=Cg*cf; _szdAux[e]=max(max(min(_szxCell[e].ds,0.5*(h-H)),0.25*_szxCell[e].ds),1e-3)/_hf;
      continue;
    }
    _szxHead[e*st]=H/_hf; _szxCond[e*st]=Cg*cf;
    _szqQ[e*sq]=Cg*(Hg-H)/_qf-cap;
  }
  if (_szwQ!=NULL){ //upward flows out of the cell holding the water table below a dried top cell (fixed rates for the step)
    const int sw=max(_szwNcol,1);
    for (size_t j=0;j<_szwNode.size();j++){_szwQ[j*sw]=0.0;}
    for (size_t e=0;(e<_szxLink.size()) && (e+1<_szwStart.size());e++){
      const int up=(e<_szxUp.size())?_szxUp[e]:-1; if ((up<0) || (up==TopNode(_aLinks[_szxLink[e]].ic)) || _szxUpCHD[e]){continue;}
      const double q=((e<_szxCapQ.size())?_szxCapQ[e]:0.0)+((e<_szxRepay.size())?_szxRepay[e]:0.0)+((e<_vzUp.size())?_vzUp[e]+_vzUpA[e]:0.0);
      for (int j=_szwStart[e];j<_szwStart[e+1];j++){if (_szwNode[j]==up){_szwQ[j*sw]-=q/_qf; break;}}
    }
  }
  if (steady){return;}
  SZRecharge(_szTstep,true);
  for (size_t e=0;e<_szxLink.size();e++){ //the outlet above: the link's conductance, at least what carries the link's inflow (its
    if (!pureRate[e]){continue;}            //rate and its share of the cell's recharge) out over one soil-zone thickness above the
    const gw_link &L=_aLinks[_szxLink[e]];  //base, so that a solve without interface updates keeps the head within the soil zone
    const double Qin=max(_szxQc[e],0.0)+max(_Qcol[L.ic],0.0)*min(L.area/max(MF6CellArea(L.ic),1e-12),1.0);
    const double Cd=max(_szxC[e],Qin/max(_szD[L.k],0.05));
    _szdCond[e*sd]=Cd*_hf/_qf;
  }
  if (_evtCol.size()>0){ //water-table ET only on the unsaturated share of the cell: saturated soil already evaporates in Raven
    vector<double> Fa(_ncpl,0.0);
    for (size_t e=0;e<_szxLink.size();e++){const gw_link &L=_aLinks[_szxLink[e]]; Fa[L.ic]+=L.area*_linkF[_szxLink[e]];}
    int es=max(_evtNcol,1);
    for (size_t i=0;i<_evtCol.size();i++){
      int ic=_evtCol[i]; double f=(_colCover[ic]>0)?Fa[ic]/_colCover[ic]:0.0;
      _evtRate[i*es]=_evtBase[i]*(1.0-min(max(f,0.0),1.0));
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief one interface update after a MODFLOW outer iteration; returns whether the interface has converged
//
bool CGroundwaterModel::SZIterate(double **aPhinew,const int iter,const bool mf6conv,const bool steady,double &prevQ,const bool last)
{
  (void)aPhinew;
  //Schedule: in the first outer iterations the interface is updated after every one (on easy steps both converge
  //together); after that only when MODFLOW has converged on the interface it was given, so that the changing
  //boundary does not keep MODFLOW's own iteration (its under-relaxation and backtracking) from converging. The update
  //budget counts the updates made after MODFLOW had converged; when it is used up the values are frozen.
  const int nEarly=min(5,max(1,_maxiter/5));
  if (_szMaxIter==0){return true;}                                            //explicit: the predictor only
  //the last outer iterations (a tenth, at least 2) are left to MODFLOW on the boundary it has: an update just before
  //the limit would leave MODFLOW's final solve unconverged, and its residual would be water created or lost. The
  //update after iteration SZFreezeAt() and the last one the update budget allows make their boundary supply-limited,
  //as in explicit mode. Once frozen nothing changes any more (a supply-limited update after the freeze left too few
  //solves for MODFLOW to absorb it when MODFLOW could not converge). When MODFLOW does not converge in
  //the iterations after the first ones, the last update can be an earlier, unlimited one; with a low outer-iteration
  //limit its boundary can then take more from a column than it holds (reported as soil-zone shortfall, like MODFLOW's
  //own residual on such steps)
  const int freezeAt=SZFreezeAt();
  bool frozen=(iter>freezeAt);
  if (!frozen && !(!mf6conv && (iter>nEarly))){if (mf6conv){_szExtra++;} frozen=(_szExtra>_szMaxIter);}
  if (frozen){
    if (_szIter==0){_szConv=false;}                                           //no update at all this step: not iterated
    return true;
  }
  if (!mf6conv && (iter>nEarly)){return false;}                               //let MODFLOW converge first
  //state MODFLOW just solved with: kept if this update shows the interface has converged, so that MODFLOW's budget
  //(computed from the arrays it solved with) and the flows Raven books stay exactly consistent
  vector<double> F0(_linkF),hc0(_szxHc),lam0(_szxLam),q0(_szxQc),cg0(_szxCg),rej0(_szxRej),srej0(_szRej),fc0(_szxFc); vector<char> one0(_szxOne);
  SZLinkFractions(true,mf6conv);
  if (!steady){SZRecharge(_szTstep,false);}
  _szCapUpd=(iter==freezeAt) || (_szExtra>=_szMaxIter); //the last update possible (iterations or budget): its boundary stays,
                                                        //supply-limited as in explicit mode
  _szOneOK=(_maxiter-iter>=10);
  SZLocalSolves(_szTstep,steady);
  _szCapUpd=false;
  //total exchange of the updated interface at the current heads (into the aquifer) [m3/d] (the boundary MODFLOW gets
  //reproduces it at these heads, whatever Cg)
  double Q=0.0,G=0.0,sC=0.0;
  for (size_t e=0;e<_szxLink.size();e++){Q+=_szxQc[e]; G+=fabs(_szxQc[e]); sC+=_szxFc[e]*_szxC[e];}
  //flow test: relative change, or a change no larger than the head tolerances allow on every link
  //(sum of conductances x the larger of the interface tolerance and MODFLOW's DVCLOSE)
  double dQ=fabs(Q-prevQ);
  _szDQ=dQ; _szDQtol=max(_szTolQ*max(G,1e-6),sC*max(_szTolH,(IsExternal() && (_extDvclose>0))?_extDvclose:_dvclose));
  bool conv=(_szResid<_szTolH) && (dQ<=_szDQtol);
  prevQ=Q; _szIter++; _szConv=conv; (void)iter;
  if (conv || last){ //keep what MODFLOW solved with (the update changed it by less than the tolerances, or no solve follows:
                     //the step's flows must be those of the arrays of MODFLOW's last solve)
    _linkF=F0; _szxHc=hc0; _szxLam=lam0; _szxQc=q0; _szxCg=cg0; _szxRej=rej0; _szRej=srej0; _szxFc=fc0; _szxOne=one0;
    if (!steady){SZRecharge(_szTstep,false);} //_Qcol of those fractions
  }
  else {SZSetArrays(steady);}
  return conv;
}

//////////////////////////////////////////////////////////////////
/// \brief applies MODFLOW's booked exchange (_szxQ, after the non-convergence guard), the rejected recharge and the
///  overflow to the soil-zone columns, Raven's soil-zone layers (their area-weighted mean) and SURFACE_WATER
///  (design note sections 7.6, 9)
//
void CGroundwaterModel::SZApply(double **aPhinew,const double &tstep)
{
  int iSW=_pModel->GetStateVarIndex(SURFACE_WATER);
  const int n=_szNLayers;
  _szRejVol=_szOvfVol=_szAppliedVol=_szShortVol=0.0; _szDebtVol=0.0;
  if (_szxDebt.size()!=_szxLink.size()){_szxDebt.assign(_szxLink.size(),0.0);}
  vector<double> phi(n),tg(n);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    double A=_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
    for (int j=0;j<n;j++){tg[j]=aPhinew[k][_szSV[k][j]];}
    double ovfV=SZSyncColumns(k,&tg[0])/MM_PER_METER*A;                         //no-op unless SOIL changed since the exchange began
    double netV=0.0,shV=0.0,rejV=0.0;                                            //[m3]
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      int e=_szxOfHRU[k][jj],c=_szxCol[e]; double a=_szColA[c];
      for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
      double ovf=gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);           //[mm over the column]
      double dV=-_szxQ[e]*tstep/a*MM_PER_METER;                                  //[mm] change of the column
      double sh=0.0;
      ovf+=gwsz::Apply(_szLay[k],phi,dV,sh);
      if ((e<(int)_szxLagD.size()) && (_szxLagD[e]!=0.0)){                       //LAGGED: the part of the source not taken back; what the
        const double D=_szxLagD[e]; double shD=0.0;                              //column cannot give back (water Raven's processes used that
        if (D>0.0){ovf+=gwsz::Fill(_szLay[k],phi,D);}                            //the step's exchange did not supply) is a debt the aquifer
        else {gwsz::Apply(_szLay[k],phi,D,shD);}                                 //repays in the next steps
        dV+=D+shD;
        if (shD>0.0){const double v=shD/MM_PER_METER*a; _szxDebt[e]+=v; _szDebtVol+=v;}
      }
      if ((_vzMode>0) && (e<(int)_vzCapV.size()) && (_vzCapV[e]>0.0)){             //capillary rise from the unsaturated column into its layers
        for (int j=0;j<n;j++){
          const double v=_vzCapV[e]*_vzCapL[e*n+j]/a*MM_PER_METER; if (!(v>0.0)){continue;}
          const double add=min(v,max(_szLay[k][j].phimax-phi[j],0.0));
          phi[j]+=add; ovf+=gwsz::Fill(_szLay[k],phi,v-add); dV+=v;
        }
      }
      if ((e+1)*n<=(int)_szxCapL.size()){                                          //capillary rise from the water table into its layers
        for (int j=0;j<n;j++){
          const double v=_szxCapL[e*n+j]*tstep/a*MM_PER_METER; if (!(v>0.0)){continue;}
          const double add=min(v,max(_szLay[k][j].phimax-phi[j],0.0));
          phi[j]+=add; ovf+=gwsz::Fill(_szLay[k],phi,v-add); dV+=v;
        }
      }
      for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
      ovfV+=ovf/MM_PER_METER*a; netV+=(dV+sh)/MM_PER_METER*a; shV+=sh/MM_PER_METER*a; rejV+=_szxRej[e];
    }
    for (size_t jj=0;jj<_szcOfHRU[k].size();jj++){ //columns over specified heads: the local solve's exchange with the fixed head
      int e=_szcOfHRU[k][jj],c=_szcCol[e]; double a=_szColA[c];
      if (!(a>0.0)){continue;}
      for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
      double ovf=gwsz::Fill(_szLay[k],phi,_szcRej[e]/a*MM_PER_METER);
      double dV=-_szcQ[e]*tstep/a*MM_PER_METER,sh=0.0;
      ovf+=gwsz::Apply(_szLay[k],phi,dV,sh);
      for (int j=0;j<n;j++){
        const double v=_szcCapL[e*n+j]*tstep/a*MM_PER_METER; if (!(v>0.0)){continue;}
        const double add=min(v,max(_szLay[k][j].phimax-phi[j],0.0));
        phi[j]+=add; ovf+=gwsz::Fill(_szLay[k],phi,v-add); dV+=v;
      }
      for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
      ovfV+=ovf/MM_PER_METER*a; netV+=(dV+sh)/MM_PER_METER*a; shV+=sh/MM_PER_METER*a; rejV+=_szcRej[e];
    }
    //Raven's soil-zone layers: area-weighted mean of the columns (the rest column did not change)
    int nc=SZNCols(k),o=_szColOff[k]; double sw=0.0,W=0.0;
    for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
    for (int j=0;j<n;j++){
      double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*n+j];}
      v=(sw>0.0)?v/sw:tg[j];
      aPhinew[k][_szSV[k][j]]=v; W+=v;
    }
    _szW[k]=W;
    if (ovfV>0.0){aPhinew[k][iSW]+=ovfV/A*MM_PER_METER;}
    if (k<(int)_szLagNet.size()){netV+=_szLagNet[k];}                           //LAGGED: source and its reversal (given to Raven's soil zone at the sync)
    _szAppliedVol+=rejV+netV;                                                    //net water added to Raven (overflow stays in Raven)
    _szRejVol+=rejV; _szShortVol+=shV;
    ovfV+=_szSyncOvf[k];                                                         //(sent to SURFACE_WATER when the exchange began)
    _szOvfVol+=ovfV; _szOvfV[k]=ovfV; _szNetV[k]=netV;
  }
  _szxQprev=_szxQ;                                                               //(the LAGGED source of the next step)
  for (size_t e=0;e<_szxRepay.size() && e<_szxDebt.size();e++){_szxDebt[e]=max(_szxDebt[e]-_szxRepay[e]*tstep,0.0);}
  if ((_szShortVol>1e-6) && (++_szShortN<=10)){
    WriteWarning("CGroundwaterModel: soil-zone exchange: "+to_string(_szShortVol)+" m3 that MODFLOW took from the soil zone was not in it (step "+to_string(_stepCount)+")",false);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief shortfall [m3] the booked exchange (_szxQ) and rejected recharge would leave in the columns: SZApply's
///  booking without changing anything (for the guard in Exchange, which withholds it from the discharge to Raven)
//
double CGroundwaterModel::SZTrialShortfall(const double &tstep) const
{
  const int n=_szNLayers; double S=0.0;
  vector<double> phi(n);
  for (size_t e=0;e<_szxLink.size();e++){
    int c=_szxCol[e],k=_aLinks[_szxLink[e]].k; double a=_szColA[c];
    if (_colCHDTop[_aLinks[_szxLink[e]].ic]){continue;}                          //(exchange with a specified head: not MODFLOW's water)
    double dV=-_szxQ[e]*tstep/a*MM_PER_METER;
    if (!(dV<0.0)){continue;}                                                    //only flow out of the column can fall short (LAGGED: the
                                                                                 //source not taken back is settled after it, its shortfall a debt)
    for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
    gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);
    double sh=0.0; gwsz::Apply(_szLay[k],phi,dV,sh);
    S+=sh/MM_PER_METER*a;
  }
  return S;
}

//////////////////////////////////////////////////////////////////
/// \brief GWSoilExchange.csv rows of this step (HRU values: area-weighted over its exchange columns)
//
void CGroundwaterModel::SZWriteStep(const double &t,const string &dstr)
{
  if (!_SZOUT.is_open()){return;}
  CHRUGroup *pG=(_szOutGroup=="")?NULL:_pModel->GetHRUGroup(_szOutGroup);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    if ((pG!=NULL) && !pG->IsInGroup(k)){continue;}
    double A=_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
    double sa=0,zb=0,lam=0,hc=0,Fk=0,qa=0,qs=0;
    for (size_t j=0;j<_szxOfHRU[k].size();j++){
      int e=_szxOfHRU[k][j]; const gw_link &L=_aLinks[_szxLink[e]];
      vector<double> phi(_szColPhi.begin()+_szxCol[e]*_szNLayers,_szColPhi.begin()+(_szxCol[e]+1)*_szNLayers); //column after the step
      sa+=L.area; zb+=L.area*_szxZb[e]; lam+=L.area*gwsz::SubgridLevel(_szxCell[e],gwsz::SatThickness(_szLay[L.k],phi)); hc+=L.area*WaterTable(L.ic); Fk+=L.area*_linkF[_szxLink[e]];
      if (_szxQ[e]>0){qa+=_szxQ[e];} else {qs-=_szxQ[e];}
    }
    double lv=0.0;
    if (sa>0){zb/=sa; lv=lam/sa; hc/=sa; Fk/=sa;}
    _szLam[k]=lv-zb;
    double Wmax=0; for (size_t j=0;j<_szSV[k].size();j++){Wmax+=_szLay[k][j].phimax;}
    double W=_szW[k];
    _SZOUT<<setprecision(10)<<t<<","<<dstr<<","<<_pModel->GetHydroUnit(k)->GetHRUID()<<","<<zb<<","<<lv<<","<<hc<<","<<Fk<<","
          <<qa<<","<<qs<<","<<_szRej[k]/A*MM_PER_METER<<","<<_szOvfV[k]/A*MM_PER_METER<<","<<W<<","<<Wmax<<","<<_szIter<<","<<_szResid<<endl;
  }
}
