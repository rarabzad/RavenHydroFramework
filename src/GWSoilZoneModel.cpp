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

static const double SZ_WSREL=0.1; ///< wrong-sign exemption of the per-link boundary test: at most this share of the link's flow [-]
static const double SZ_WSH  =1.0; ///< ... or at least this many head tolerances times the link's K

//////////////////////////////////////////////////////////////////
/// \brief soil-zone layers of every coupled HRU, thickness, and the checks of the soil-zone set-up
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
  //(the warning about an outer-iteration limit below 20 is written once MODFLOW's limit is known, after ConnectEngine: an
  // existing model's is its IMS OUTER_MAXIMUM)
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
//////////////////////////////////////////////////////////////////
/// \brief bed and stage datum of surface-store exchange entry e (SWX_RAVEN) [model units]: the model top; with :SoilZone PHYSICAL
///  the ground of the group's HRUs in the column (the soil zone lies between the model top and the ground), area-weighted:
///  an HRU with a soil zone at its soil-zone base plus its thickness, one without at the model top (where its seepage leaves)
//
double CGroundwaterModel::SWXDatum(const size_t e) const
{
  const int ic=_swxCol[e];
  const double zt=IsExternal()?_extTopRaw[ic]:_top[ic];
  if (!SZPhysical() || (e>=_swxLinks.size())){return zt;}
  double A=0.0,zA=0.0;
  for (size_t j=0;j<_swxLinks[e].size();j++){
    const gw_link &L=_aLinks[_swxLinks[e][j]];
    const double g=((L.k<(int)_szOn.size()) && _szOn[L.k])?SZBase(L.k,ic)+_szD[L.k]:_top[ic];
    A+=L.area; zA+=L.area*g;
  }
  if (A<=0.0){return zt;}
  return zt+(zA/A-_top[ic])/_hf;
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
  int nRatio=0,nMix=0; double rmin=1e300,rmax=0,dmix=0.0,fmin=1e300,fmax=0.0; //(fmin, fmax: the ratios of the links outside a factor 2)
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
    if ((r<0.5) || (r>2.0)){nRatio++; fmin=min(fmin,r); fmax=max(fmax,r);}
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
    //no base below the top cell's bottom + 1 mm, the lowest level MODFLOW's drains of the link can have (SZSetArrays): in
    //the 1 mm band below it the local solve's upflow and the drain MODFLOW gets would differ (in a rough-terrain case the
    //booked flows would be 68 % short)
    {const double zf=_botm[TopNode(ic)]+1e-3; for (size_t q=0;q<cell.zbs.size();q++){cell.zbs[q]=max(cell.zbs[q],zf);}}
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
                 "n(1-SF)(1-FIELD_CAPACITY) of the bottom soil-zone layer by more than a factor 2 ("+((fmin==fmax)?"ratio "+to_string(fmin):"ratio "+to_string(fmin)+" to "+to_string(fmax))+
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
  _szxRs.assign(ne,0.0); _szxRq.assign(ne,0.0); _szxRqP.assign(ne,0.0); _szxRl.assign(ne,0.0); _szxRqc.assign(ne,0.0); _szxF0.assign(ne,0.0);
  _szxDrn.assign(ne,0); _szxNegCap.assign(ne,1e300); _vzUsed.assign(ne,0.0); _szxDropped.assign(ne,0);
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
      _szColInit=true;
      bool same=true; //(a file written at full precision: the .rvc soil water is the columns' mean, and the columns continue as they were)
      {double sw=0.0; for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
       for (int j=0;j<n;j++){double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*n+j];} if ((sw>0.0) && (fabs(v/sw-tg[j])>1e-9)){same=false;}}}
      if (same){continue;}
      double ovf=SZSyncColumns(k,&tg[0]); //hotstart columns and the .rvc soil water normally agree; any difference is shared
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
      if (jt!=it->second.end()){_szxQprev[e]=jt->second; if (_szTiming==1){_szxDebt[e]=max(_hotDebt[it->first][L.ic+1],0.0);}} //(only LAGGED repays a debt)
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
    if ((e<_szxDropped.size()) && _szxDropped[e]){                               //(the update after a failed connection test: plain?
      if (!relax || !mf6conv || (_szxHr[e]==-HUGE_VAL)){_szDropPlain++;}        // GWModelSummary.txt, run totals)
      _szxDropped[e]=0;
    }
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
/// \brief what the top cell of every exchange link can still give to a negative recharge this step [m3/d] (_szxNegCap):
///  the cell's limit for upward flows (SZUpCap: half its water above its bottom under the link, faded out as it nears
///  drying) less what the unsaturated columns (_vzUsed, VZStep), the capillary rise from the water table (SZCapRiseWT)
///  and the repayment of the lagged source's debt (SZRepayRates) took out of it before, in that order, so that the total
///  out of a cell never exceeds its limit. A top cell that has dried (upward flows come from a lower cell, SZW_RAVEN) can
///  give nothing; a specified head everything
//
void CGroundwaterModel::SZNegCaps(const double &tstep)
{
  const size_t ne=_szxLink.size();
  _szxNegCap.assign(ne,0.0);
  for (size_t e=0;e<ne;e++){
    const int ic=_aLinks[_szxLink[e]].ic;
    if (_colCHDTop[ic]){_szxNegCap[e]=1e300; continue;}
    if ((e>=_szxUp.size()) || (_szxUp[e]!=TopNode(ic))){continue;}
    double used=(e<_vzUsed.size())?_vzUsed[e]:0.0;
    if (e<_szxCapQ.size()){used+=_szxCapQ[e];}
    if (e<_szxRepay.size()){used+=_szxRepay[e];}
    _szxNegCap[e]=max(SZUpCap((int)e,tstep)-used,0.0);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief the rejected recharge of SZX entry e that follows the head h [m3/d] (and its slope dR/dh [m2/d]): the link's
///  recharge where no unsaturated column carries it, R_l F(h), and the water its column gives back to the soil zone on
///  ground that saturated since the start of the step (rule (i)), R_q clip(F(h) - F0, 0, 1 - F0) with R_q the column's
///  outputs over 1 - F0 (one depth class: F of the link)
//
double CGroundwaterModel::SZRejRate(const int e,const double &h,double *slope) const
{
  if (slope!=NULL){*slope=0.0;}
  if ((e<0) || (e>=(int)_szxRl.size())){return 0.0;}
  const double Rl=_szxRl[e],Rq=(e<(int)_szxRqc.size())?_szxRqc[e]:0.0;
  if ((Rl<=0.0) && (Rq<=0.0)){return 0.0;}
  if ((_vzMode>0) && (e<(int)_vzUse.size()) && _vzUse[e] && (e+1<(int)_vzOff.size())){
    //unsaturated columns carried the link's share this step, class by class (depth classes, :SoilZoneVadoseClasses): rule (i)
    //per class column (VZBackRate: R_q,c clip(F_c - F0_c, 0, 1 - F0_c)), and classes without a column reject their share
    //with their own fraction, R_l a_c F_c(h)
    double R=VZBackRate(e,h);
    const bool expl=(_szMaxIter==0) || (SZFreezeAt()<1);
    for (int i=_vzOff[e];i<_vzOff[e+1];i++){
      const double Fc=VZClassF(e,i,h);
      const double dF=(slope==NULL)?0.0:gwsz::SatFractionSlope(h,(_vzOff[e+1]-_vzOff[e]<=1)?_szxCell[e].zbs:_vzcZbs[i],_szxCell[e].ds);
      if (!_vzcUse[i]){R+=Rl*_vzcA[i]*Fc; if (slope!=NULL){*slope+=Rl*_vzcA[i]*dF;} continue;}
      const double F0=_vzcF0[i];
      if ((slope!=NULL) && !expl && (_vzcRq[i]>0.0) && (F0<1.0-1e-12) && (Fc>F0) && (Fc<1.0)){*slope+=_vzcRq[i]*dF;}
    }
    return R;
  }
  const gwsz::sz_cell &cell=_szxCell[e];
  const double F=gwsz::SatFraction(h,cell.zbs,cell.ds);
  double R=Rl*F,dF=0.0;
  if (slope!=NULL){dF=gwsz::SatFractionSlope(h,cell.zbs,cell.ds); *slope=Rl*dF;}
  if (Rq>0.0){
    const double F0=_szxF0[e],x=F-F0;
    if (x>0.0){R+=Rq*min(x,1.0-F0); if ((slope!=NULL) && (x<1.0-F0)){*slope+=Rq*dF;}}
  }
  return R;
}
bool CGroundwaterModel::SZRejFollows(const int e) const
{
  if ((e<0) || (e>=(int)_szxRl.size())){return false;}
  if (_colCHDTop[_aLinks[_szxLink[e]].ic]){return false;}
  return (_szxRl[e]>0.0) || ((e<(int)_szxRqc.size()) && (_szxRqc[e]>0.0));
}

//////////////////////////////////////////////////////////////////
/// \brief recharge to MODFLOW with rejection on saturated ground (design note section 6.1): _Qcol, RCH rates, _szxRej, _szRej;
///  the head-dependent part of the rejection (_szxRl, _szxRqc, _szxF0: SZRejRate) and its value and slope at the connection
///  head (_szxRq, _szxRqP), from which SZLocalSolves makes the rejection's linear term
//
void CGroundwaterModel::SZRecharge(const double &tstep,const bool write)
{
  std::fill(_Qcol.begin(),_Qcol.end(),0.0);
  _szRchCHD=0.0; _vzBackV=_vzBackWV=0.0;
  const size_t ne=_szxLink.size();
  _szxRl.assign(ne,0.0); _szxRqc.assign(ne,0.0); _szxRq.assign(ne,0.0); _szxRqP.assign(ne,0.0);
  if (_szxF0.size()!=ne){_szxF0.assign(ne,0.0);}
  if (_szxNegCap.size()!=ne){_szxNegCap.assign(ne,1e300);}
  _szNegCutQ=0.0;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; double V=_szRchV[k],rej=0.0;
    double Aall=_hruRchArea[k];                                                  //(with a soil zone: the columns over specified heads too)
    if (_szOn[k]){ //(all its links' areas: also those on cells that are specified heads this step)
      Aall=0.0; for (size_t j=0;j<_hruLinks[k].size();j++){Aall+=_aLinks[_hruLinks[k][j]].area;}
    }
    if (Aall<=0){if (_szOn[k]){_szRej[k]=0.0;} continue;} //all on specified-head cells (counted in _QrchCHD)
    //a negative recharge (Raven's processes take water out of GROUNDWATER, e.g. PERC_GR4JEXCH with GR4J_X2 > 0) is not split
    //by the saturated fraction: all of it is taken from MODFLOW (the top cell's recharge; on a column link the column's share,
    //VZStep), nothing is rejected for saturation. MODFLOW's recharge would take it out of the top cell whatever the cell holds
    //(a cell drained by it can fall to heads of -1e4 m and stop the run), so on an exchange link it is at
    //most what the top cell can still give this step (_szxNegCap, SZNegCaps); the rest is taken back out of the soil zone,
    //as a negative rejected recharge (SZApply). HRUs without a soil zone keep the plain recharge
    const bool neg=(V<0.0);
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
        const double share=(W>0.0)?_szcFw[ec]/W:L.area/Aall,F=neg?0.0:_linkF[li];
        _szcRej[ec]=V*share*F; rej+=_szcRej[ec]; _szRchCHD+=V/tstep*share*(1.0-F);
        continue;
      }
      if (_szOn[k] && _colCHDTop[L.ic] && (_linkSzx[li]>=0)){ //the same for an exchange column whose cell is a specified head this
        const int e2=_linkSzx[li];                                                //step (existing models)
        const double share=(W>0.0)?_szxFw[e2]/W:L.area/Aall,F=neg?0.0:gwsz::SatFraction(_X[_reduced[TopNode(L.ic)]]*_hf,_szxCell[e2].zbs,_szxCell[e2].ds);
        _szxRej[e2]=V*share*F; rej+=_szxRej[e2]; _szRchCHD+=V/tstep*share*(1.0-F);
        if ((_vzMode>0) && (e2<(int)_vzAbsQ.size())){_szRchCHD+=_vzAbsQ[e2];}     //(the water of its unsaturated columns, which ended: VZStep)
        continue;
      }
      if (_colCHDTop[L.ic]){continue;}
      //(an HRU without a soil zone partly on specified-head cells: shares of its whole linked area, the rest went to the
      // specified heads in Exchange)
      const double Ano=(_hruRchArea[k]<_hruLinkedArea[k])?_hruLinkedArea[k]:_hruRchArea[k];
      double share=(W>0.0)?((_linkSzx[li]>=0)?_szxFw[_linkSzx[li]]/W:0.0):L.area/(_szOn[k]?Aall:Ano),F=_linkF[li];
      const int e=_linkSzx[li];
      if ((_vzMode>0) && (e>=0) && (e<(int)_vzRl.size())){
        _vzRl[e]=V/tstep*share;                                                  //(for SZRejRate)
        if (_vzBackW[e]>0.0){rej+=_vzBackW[e]; _vzBackWV+=_vzBackW[e];}          //rule (ii): water of vanished columns beyond Sy dh, back to the soil zone
        _Qcol[L.ic]+=_vzPass[e];                                                 //(a negative recharge share passed on to the cell)
      }
      if ((_vzMode>0) && (e>=0) && (e<(int)_vzUse.size())){
        //water a risen water table took from the columns (> 0, the top cell's recharge); the new nodes of a fallen one (< 0:
        //out of a specified head below, which supplies it; out of a lower cell: SZW_RAVEN; else the top cell). With depth
        //classes the link's columns may do both in one step: the takes (_vzAbsN) are booked by their source, the rest is the
        //top cell's
        const bool chdUp=(e<(int)_szxUpCHD.size()) && _szxUpCHD[e];
        if (_vzOff[e+1]-_vzOff[e]<=1){
          if ((_vzAbsQ[e]<0.0) && chdUp){_szRchCHD+=_vzAbsQ[e];}
          else if (!(_vzUpA[e]>0.0)){_Qcol[L.ic]+=_vzAbsQ[e];}
        }
        else {
          const double take=_vzAbsN[e],give=_vzAbsQ[e]+take;                      //[m3/d] taken and released
          if (chdUp){_szRchCHD-=take; _Qcol[L.ic]+=give;}
          else if (_vzUpA[e]>0.0){_Qcol[L.ic]+=give;}
          else {_Qcol[L.ic]+=_vzAbsQ[e];}
        }
        if (_vzUse[e]){ //the unsaturated columns carried the unsaturated share this step: their bottom flux recharges the
                        //aquifer, and the saturated share is rejected with the saturated fraction their inflow used
          //(a negative share is not split by F and went to the cell, VZStep: what the cell could not give is in _vzRejV)
          double rj=(neg?0.0:V*share*_vzFp[e])+_vzRejV[e];
          if (neg && (_vzRejV[e]<0.0)){_szNegCutQ-=_vzRejV[e]/tstep;}
          if (_vzOff[e+1]-_vzOff[e]<=1){
            if ((_vzQ[e]<0.0) && chdUp){_szRchCHD+=_vzQ[e];}                        //(upflow out of a specified head below: supplied by it,
            else {_Qcol[L.ic]+=(_vzUp[e]>0.0)?0.0:_vzQ[e];}                         // MODFLOW does not see it; from a lower cell: SZW_RAVEN)
          }
          else {
            const double up=_vzQN[e],down=_vzQ[e]+up;                             //[m3/d] upflow and downflow of the link's columns
            if (chdUp){_szRchCHD-=up; _Qcol[L.ic]+=down;}
            else if (_vzUp[e]>0.0){_Qcol[L.ic]+=down;}
            else {_Qcol[L.ic]+=_vzQ[e];}
            //depth classes without a column this step: their unsaturated share reaches the water table in the step and
            //their saturated share is rejected, with their own fraction at the connection head (as a link without columns)
            double au=0.0,fu=0.0;
            for (int i=_vzOff[e];i<_vzOff[e+1];i++){if (!_vzcUse[i]){au+=_vzcA[i]; fu+=_vzcA[i]*VZClassF(e,i,_szxHc[e]);}}
            if (au>0.0){
              if (neg){ //(a negative recharge: not split by F; at most what the top cell can still give after the classes with a
                        // column, the capillary rise and the repayment took (_szxNegCap): the rest back to the soil zone, as on a link
                        // without columns)
                const double q=V/tstep*share*au,take=max(q,-_szxNegCap[e]);
                _Qcol[L.ic]+=take; rj+=(q-take)*tstep; _szNegCutQ+=take-q;
              }
              else {_Qcol[L.ic]+=V/tstep*share*(au-fu); rj+=V*share*fu;}
            }
          }
          const double bk=VZBackRate(e,_szxHc[e]);                               //rule (i): the columns' outputs on ground saturated since the
          if (bk>0.0){_Qcol[L.ic]-=bk; rj+=bk*tstep; _vzBackV+=bk*tstep;}        // start of the step go back to the soil zone
          rej+=rj; _szxRej[e]=rj+((_vzBackW[e]>0.0)?_vzBackW[e]:0.0);
          //the head-dependent part (SZRejRate, class by class): the share of classes without a column, and rule (i)
          bool unused=false; for (int i=_vzOff[e];i<_vzOff[e+1];i++){if (!_vzcUse[i]){unused=true;}}
          if (!neg && unused){_szxRl[e]=V/tstep*share;}                          //(only classes without a column follow with R_l)
          double rq=0.0; for (int i=_vzOff[e];i<_vzOff[e+1];i++){if (_vzcUse[i]){rq+=_vzcRq[i]; _szxF0[e]=_vzcF0[i];}}
          if (!((_szMaxIter==0) || (SZFreezeAt()<1))){_szxRqc[e]=rq;}            //(explicit coupling: nothing follows the head)
          if (SZRejFollows(e)){double sl=0.0; _szxRq[e]=SZRejRate(e,_szxHc[e],&sl); _szxRqP[e]=sl;}
          continue;
        }
      }
      if (neg){ //all of the negative recharge to MODFLOW, at most what the top cell can still give; the cut back to the soil zone
        const double q=V/tstep*share,cap=(e>=0)?_szxNegCap[e]:1e300;           //[m3/d]
        const double take=max(q,-cap);
        _Qcol[L.ic]+=take;
        if (e>=0){
          _szxRej[e]=(q-take)*tstep; rej+=_szxRej[e]; _szNegCutQ+=take-q;
          if ((_vzMode>0) && (e<(int)_vzBackW.size()) && (_vzBackW[e]>0.0)){_szxRej[e]+=_vzBackW[e];} //rule (ii): the water of a column that vanished this
        }                                                                        // step goes back to the soil zone (in rej already, above)
        else if (q<take){_Qcol[L.ic]+=q-take;}                                 //(no soil zone: unlimited)
        continue;
      }
      _Qcol[L.ic]+=V/tstep*share*(1.0-F);
      rej+=V*share*F;
      if (e>=0){
        _szxRej[e]=V*share*F+(((_vzMode>0) && (e<(int)_vzBackW.size()) && (_vzBackW[e]>0.0))?_vzBackW[e]:0.0);
        _szxRl[e]=V/tstep*share;                                                 //(the rejection follows the connection head: SZRejRate)
        double sl=0.0; _szxRq[e]=SZRejRate(e,_szxHc[e],&sl); _szxRqP[e]=sl;
      }
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
      {double shR=0.0; gwsz::Apply(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER,shR);} //rejected recharge back to the column (negative: the cut
                                                                                 // taken back out, as SZApply does)
      if ((e<_szxLagD.size()) && (_szxLagD[e]!=0.0)){                            //LAGGED: the part of the source not taken back
        if (_szxLagD[e]>0.0){gwsz::Fill(_szLay[k],phi,_szxLagD[e]);} else {double sh=0.0; gwsz::Apply(_szLay[k],phi,_szxLagD[e],sh);}
      }
      double dV,ovf;
      if (_colCHDTop[L.ic]){ //a specified head this step (existing models): no MODFLOW boundary; the exchange with the head is
                             //found after MODFLOW's solve (SZChdExchange, in the booking)
        _szxLam[e]=gwsz::SatThickness(_szLay[k],phi); _szxQc[e]=0.0; _szxCg[e]=0.0; _szxOne[e]=0; _szxFc[e]=gwsz::SatFraction(h,cell.zbs,cell.ds); _szxHc[e]=h;
        _szxRs[e]=0.0;
        continue;
      }
      if (_szxKink.size()!=_szxLink.size()){_szxKink.assign(_szxLink.size(),0);}
      _szxKink[e]=(gwsz::StackFreeWater(_szLay[k],phi)<=1e-9)?1:0;              //(no free water: below the kink nothing flows)
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
    //(water created), and a slope small enough to prevent that pumps water up from below the contact: the line is
    //given to MODFLOW as a drain instead (SZSetArrays), which stops at Hg, and the secant
    //floor keeps Hg at or above the end of the contact.
    //(Iterated, the converged flow is the local solve's, which respects both limits.)
    //the rejected recharge's linear term s (hc - h) (SZSetArrays): the slope of the head-dependent rejection at the connection
    //head, at least the chord to the head MODFLOW returned last (a tangent taken on a flat end of F would let MODFLOW cross
    //the whole ramp), never negative
    double s=0.0;
    if (!steady && SZRejFollows((int)e)){
      s=max(_szxRqP[e],0.0);
      if (fabs(h-_szxHc[e])>1e-9){s=max(s,(SZRejRate((int)e,h)-SZRejRate((int)e,_szxHc[e]))/(h-_szxHc[e]));}
    }
    _szxOne[e]=0;
    if (!steady && ((_szMaxIter==0) || (SZFreezeAt()<1) || _szCapUpd)){
      //down to the lowest head the cell can reach: the bottom of the column's lowest active cell (under NEWTON the head of a
      //dried top cell falls below its bottom, and a column limited down to the top cell's bottom only can lose thousands
      //of m3/d more than it holds)
      double avail=gwsz::StackFreeWater(_szLay[k],phi)/MM_PER_METER*a/tstep;           //[m3/d]
      double zlow=_botm[TopNode(L.ic)];
      for (int l=_nlay-1;l>=0;l--){const int nd=UserNode(l,L.ic); if (_idomain[nd]>0){zlow=_botm[nd]; break;}}
      const double bot=_botm[TopNode(L.ic)];
      //the rejection's term first: at zlow it may take at most the head-dependent rejection at hc (what it stands for), and
      //no more than the free water the line leaves; then the line, with what the term
      //leaves: line and term together never take more than the column's free water
      if (s>0.0){
        const double r0=max(_szxRq[e],0.0),dz=max(_szxHc[e]-zlow,1e-3);
        s=min(s,r0/dz); s=min(s,max(avail-max(q,0.0),0.0)/dz);
      }
      const double ts=s*max(_szxHc[e]-zlow,0.0);                                  //[m3/d] the term's take at zlow
      double cap=max(avail-ts-q,0.0)/max(h-zlow,1e-3);
      if (q>=0.0){Cg=max(min(Cg,cap),1e-9*cell.C);}
      else if (Cg>cap){ //only where the line could take more than the column holds: one-sided, if at least 10 MODFLOW
                        //solves follow; with fewer MODFLOW does not settle on the changed boundary, and the line stays,
                        //as limited as the secant floor allows (a reversal can then still overdraw)
        //(not in a cell dry at the head of linearization: a drain there is inactive at that head, MODFLOW's first Newton
        //iterate sees no boundary and its head runs away)
        if (_szOneOK && (h>bot+1e-3)){_szxOne[e]=1;}
        else {Cg=max(cap,sec);}
      }
      if (_szCapUpd){ //(reported: the largest take at zlow beyond the free water, line and term together; 0 up to round-off)
        const double take=(q>=0.0)?q+((Cg>1e-6*cell.C)?Cg*max(h-zlow,0.0):0.0)+ts:ts; //(a slope at its floor is given as a rate, SZSetArrays)
        _szCapN++; _szCapExcess=max(_szCapExcess,take-avail);
      }
    }
    _szxRs[e]=s;
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
  if (_szxDrn.size()!=_szxLink.size()){_szxDrn.assign(_szxLink.size(),0);}
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
    //the rejection's linear term s (hc - h) (SZLocalSolves): a GHB of head hc and conductance s, combined with the exchange's
    const double s=((e<_szxRs.size()) && !steady && !_colCHDTop[L.ic])?_szxRs[e]:0.0,hr=_szxHc[e];
    _szxDrn[e]=0;
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
        if (s>0.0){const double H=max(hr,bot); _szxHead[e*st]=H/_hf; _szxCond[e*st]=s*cf; _szqQ[e*sq]+=s*(hr-H)/_qf;}
      }
      continue;
    }
    double Hg=h+_szxQc[e]/Cg;
    const bool kink=_szxOne[e] || ((e<_szxKink.size()) && _szxKink[e]);
    if ((_szxQc[e]<0.0) && (h>bot) && kink){ //water rising into a soil without free water: the kink drain. Nothing rises below
      //the level where the contact ends, and a column without free water gives nothing either, so the exchange has a kink
      //there, across which a line overshoots (it reverses into flow out of the column, which holds none) and the iteration
      //cycles (test S2); a column with free water drains below the line's zero, and keeps the line (a drain there would turn
      //the exchange's slope: the Liard link of exemption (i) would creep for 100 updates). MODFLOW gets the exchange as a drain
      //ending at the line's zero, at least
      //the cell bottom + 1 mm (no sample base lies below that level, SZBuildEntries), steepened there to carry the exchange
      //at the head of linearization (at most to the larger of Cg and the link conductance). Its onset is smoothed (DDRN)
      //over at most half the distance to the head of linearization (the flow there is the line's), at least 1 mm (a quarter
      //of d_s for the one-sided drain of a supply-limited update, MODFLOW's Newton iteration converges badly on a sharper
      //kink). (Not in a cell dry at the head of linearization: there the drain would be inactive and MODFLOW's first Newton
      //iterate would see no boundary)
      double H=max(Hg,bot),Cd=Cg;
      if (Hg<bot){Cd=min(-_szxQc[e]/max(h-bot,1e-6),max(Cg,_szxC[e]));}
      const double ds=_szxCell[e].ds;
      const double dep=_szxOne[e]?max(max(min(ds,0.5*(h-H)),0.25*ds),1e-3):max(min(0.25*ds,0.5*(h-H)),1e-3);
      _szdElev[e*sd]=H/_hf; _szdCond[e*sd]=Cd*cf; _szdAux[e]=dep/_hf; _szxDrn[e]=1;
      if (s>0.0){const double Hs=max(hr,bot); _szxHead[e*st]=Hs/_hf; _szxCond[e*st]=s*cf; _szqQ[e*sq]=s*(hr-Hs)/_qf-cap;}
      else {_szxHead[e*st]=bot/_hf; _szxCond[e*st]=0.0; _szqQ[e*sq]=-cap;}
      continue;
    }
    //the line q = Cg (Hg - h) with the term s (hc - h): (Cg + s) (Ht - h), Ht = (Cg Hg + s hc)/(Cg + s)
    const double Ct=Cg+s,Ht=(Cg*Hg+s*hr)/Ct;
    const double H=max(Ht,bot);
    _szxHead[e*st]=H/_hf; _szxCond[e*st]=Ct*cf;
    _szqQ[e*sq]=Ct*(Ht-H)/_qf-cap;
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
  //Schedule: in the first outer iterations (nEarly) the interface is updated after every one (on easy steps both converge
  //together); after that only when MODFLOW has converged on the interface it was given, so that the changing
  //boundary does not keep MODFLOW's own iteration (its under-relaxation and backtracking) from converging. The update
  //budget counts the updates made after MODFLOW had converged; when it is used up the values are frozen.
  const int nEarly=min(5,max(1,_maxiter/5));
  if ((_szMaxIter==0) || (SZFreezeAt()<1)){return true;}                      //explicit: the predictor only (also with an outer-iteration
                                                                              //limit of 1, which leaves no update: the predictor's boundary is
                                                                              //supply-limited as in explicit mode, SZLocalSolves)
  //the last outer iterations (a tenth, at least 2) are left to MODFLOW on the boundary it has: an update just before
  //the limit would leave MODFLOW's final solve unconverged, and its residual would be water created or lost. The
  //update after iteration SZFreezeAt() and the last one the update budget allows make their boundary supply-limited,
  //as in explicit mode. Once frozen nothing changes (a supply-limited update after the freeze would leave too few
  //solves for MODFLOW to absorb it when MODFLOW cannot converge). When MODFLOW does not converge in
  //the iterations after the first ones, the last update can be an earlier, unlimited one; with a low outer-iteration
  //limit its boundary can then take more from a column than it holds (reported as soil-zone shortfall, like MODFLOW's
  //own residual on such steps)
  const int freezeAt=SZFreezeAt();
  bool frozen=(iter>freezeAt);
  if (!frozen && !(!mf6conv && (iter>nEarly))){if (mf6conv){_szExtra++;} frozen=(_szExtra>_szMaxIter) || (!mf6conv && (_szExtra>=_szMaxIter));}
  //(the second: the budget's supply-limited update was made in an earlier one of the first outer iterations, which update
  //without MODFLOW's convergence and do not count: it stays the last one)
  if (frozen){
    if (_szIter==0){_szConv=false;}                                           //no update at all this step: not iterated
    return true;
  }
  if (!mf6conv && (iter>nEarly)){return false;}                               //let MODFLOW converge first
  if (!mf6conv && SZDivergedIterate()){return false;} //no update from a diverged iterate: its local solves give absurd boundaries
                                                      //(GHB heads of 1e45 m, exchanges of -5e49 m3/d) from which MODFLOW does not
                                                      //recover, whereas on the boundary it has MODFLOW alone does)
  //state MODFLOW just solved with: kept if this update shows the interface has converged, so that MODFLOW's budget
  //(computed from the arrays it solved with) and the flows Raven books stay exactly consistent
  _szEsum=0.0;
  vector<double> F0(_linkF),hc0(_szxHc),lam0(_szxLam),q0(_szxQc),cg0(_szxCg),rej0(_szxRej),srej0(_szRej),fc0(_szxFc),rs0(_szxRs); vector<char> one0(_szxOne);
  const vector<double> hq0(_szxHq),qq0(_szxQq);                               //(head and exchange of the last local solve: the secant below)
  SZLinkFractions(true,mf6conv);
  if (!steady){SZRecharge(_szTstep,false);}
  _szCapUpd=(iter==freezeAt) || (_szExtra>=_szMaxIter); //the last update possible (iterations or budget): its boundary stays,
                                                        //supply-limited as in explicit mode
  _szOneOK=(_maxiter-iter>=10);
  SZLocalSolves(_szTstep,steady);
  _szCapUpd=false;
  //the rejected recharge follows the connection too: the recharge by which the saturated fraction MODFLOW solved with
  //differs from the fraction at the current head must be within the head tolerance as well (at a converged interface
  //hc = h), in head units: over the larger of the link's conductance (as the exchange it misses, SZLocalSolves) and the
  //column's storage over the step (Sy a / dt), either of which takes up such a flow. The exchange does not see the
  //connection where the soil water stands at the water table (every sample's level at its aquifer side), and without
  //this test the interface could pass as converged with the fraction of an earlier iterate (in a rough-terrain case 91 %
  //of a link's recharge rejected with the water table below every base point of its cell)
  if (!steady){
    for (size_t e=0;e<_szxLink.size();e++){
      if (!SZRejFollows((int)e)){continue;}
      const int li=_szxLink[e]; const gw_link &L=_aLinks[li];
      if (_colCHDTop[L.ic]){continue;}
      const double K=max(_szxCell[e].C,max(SZCellSy(L.ic),0.01)*_szColA[_szxCol[e]]/_szTstep);          //[m2/d]
      if (!(K>0.0)){continue;}
      const double d=fabs(SZRejRate((int)e,_X[_reduced[TopNode(L.ic)]]*_hf)-SZRejRate((int)e,hc0[e]))/K;
      if (d>_szResid){_szResid=d; _szResE=(int)e;}
    }
    //nor may the boundary MODFLOW solved with (the exchange's line, drain or rate, and the rejection's linear term, from
    //MODFLOW's HCOF and RHS as the booking reads them) differ at the current head from what it stands for: the local
    //solve's exchange there plus the change of the rejected recharge from the connection head it was built at. The tests
    //above compare successive local solves only; a line linearized where the exchange is steep (just above the soil-zone
    //base) and continued below it, where the exchange does not depend on the head, would pass as converged although at the
    //head MODFLOW returns it carries several times the local solve's flow (existing-model round trips of test_vz_lag L5
    //and test_vz_kinematic K6: generated and existing model stop on different sides of that kink); and the rejection's
    //term s (hc - h), solved with a chord s from an earlier iterate, can carry 50 m3/d into a cell whose head has fallen
    //1.3 m below hc, where the rejection is 0 at both heads (rough suite x43 LAG). Per link, in head units:
    //the difference E over the larger of the link's conductance and the column's storage over the step (as the tests
    //above) must be below the head tolerance. Two parts that no update removes are left out, and only they:
    // - a line link (not a drain) whose local solves' exchange rises with the head (upflow falling as the head rises: the secant through
    //   the last two local solves has the wrong sign for a boundary line, whose conductance cannot be negative), the
    //   head having moved by more than the head tolerance between them with the connection held (below): the line's
    //   error from that slope (the exchange part of E; the rejection's part is still tested), up to B = max(10 % of the
    //   link's flow, the head tolerance times K); beyond B it is tested. A per-link test without the exemption stalls the
    //   Liard run on such a link (MODFLOW's head creeps 3 mm per update for 100 updates, 7 % of the link's flow). Its exchange part
    //   stays in the summed test (_szEsum, held to the flow test's tolerance), and the connection test below allows it
    //   the same B (the connection's lag is that drift);
    // - on a link MODFLOW has as a drain (the kink drain: _szxDrn) whose head is within 1 mm of the level of the drain
    //   MODFLOW solved with (and within that drain's smoothing depth): MODFLOW's smoothing of the drain's onset (DDRN)
    //   takes the difference between the unsmoothed drain and the smoothed one there, which no update removes (rough
    //   suite e106: the interface would stall on 3 to 4 steps). Only that difference, computed from MODFLOW's own drain
    //   factor, is left out: the booked flow then differs from the local solve's by at most the drain's conductance
    //   (at most C) times 1 mm (up to all of a small flow).
    //Drains are held to the column's supply by SZSolvedShortfall
    _szEsum=0.0;
    if ((_szxHcof!=NULL) && (_szqHcof!=NULL) && (_szdHcof!=NULL)){
      const int sd=max(_szdNcol,1);
      for (size_t e=0;e<_szxLink.size();e++){
        const gw_link &L=_aLinks[_szxLink[e]];
        if (_colCHDTop[L.ic]){continue;}
        const bool drn=(e<_szxDrn.size()) && _szxDrn[e];
        const double hx=_X[_reduced[TopNode(L.ic)]];                                                     //(model units)
        const double hm=hx*_hf;                                                                          //[m]
        double qs=(_szxHcof[e]*hx-_szxRhs[e])*_qf+(_szqHcof[e]*hx-_szqRhs[e])*_qf+(_szdHcof[e]*hx-_szdRhs[e])*_qf;
        const bool inQ=(e>=_szxUp.size()) || (_szxUp[e]==TopNode(L.ic));
        if (inQ && (e<_szxCapQ.size())){qs+=_szxCapQ[e];}
        if (inQ && (e<_szxRepay.size())){qs+=_szxRepay[e];}
        double E=qs-_szxQc[e],Er=0.0;                                                                    //[m3/d] (Er: the rejection's part)
        if (SZRejFollows((int)e)){Er=-(SZRejRate((int)e,hc0[e])-SZRejRate((int)e,hm));}
        if (e<rs0.size()){Er+=rs0[e]*(hc0[e]-hm);}
        E-=SZRejFollows((int)e)?SZRejRate((int)e,hc0[e])-SZRejRate((int)e,hm):0.0;
        if (!drn){_szEsum+=fabs(E);}
        if (drn){ //MODFLOW's smoothing of the drain's onset at this head: its flow less that of the unsmoothed drain, within
                  //1 mm of the level of the drain MODFLOW solved with (and within that drain's smoothing depth). (A one-sided
                  //drain, smoothed over at least d_s/4, never reaches this test: it is made only by the supply-limited
                  //last update, after which nothing is updated or tested)
          const double el=_szdElev[e*sd]*_hf,dep=_szdAux[e]*_hf,cd=_szdCond[e*sd]*_qf/_hf;               //[m], [m], [m2/d]
          if ((hm>el) && (hm-el<min(dep,1e-3)+1e-9)){E-=(_szdHcof[e]*hx-_szdRhs[e])*_qf+cd*(hm-el);} //(smoothed - unsmoothed flow into the aquifer)
          //(a drain raised to the cell bottom + 1 mm with its conductance at its limit max(C_g, C) carries all the upflow
          //the local solve gives: no sample base lies below that level (SZBuildEntries), and a sample gives upflow only where the
          //head is above its base, at most C (h - base) in all; so no upflow beyond cd (h - bottom) needs an exemption)
        }
        const double K=max(_szxCell[e].C,max(SZCellSy(L.ic),0.01)*_szColA[_szxCol[e]]/_szTstep);       //[m2/d]
        if (!(K>0.0)){continue;}
        //wrong-sign secant (exemption (i)): only through a real head change (more than the head tolerance) with the
        //connection held at the previous local solve (hc0 = hq0: the last update set the connection to its head, not to an
        //Aitken or regula-falsi estimate), so that its point is one of the converged exchange q(h; h); the current local
        //solve uses the connection this update just set, which may be extrapolated, so the secant is the slope of q(h; h)
        //up to the current solve's connection lag (on an accepted step within the head tolerance times C, the lag test of
        //SZLocalSolves); with
        //MODFLOW's head barely moving, the sign would come from the change of the connection alone (on a Liard step: the head
        //moving 1.6e-5 m, hc0 3.2 cm above it, secant -4.5e7 m2/d, 47 % of the link's flow left out). The
        //part left out is bounded per link: up to SZ_WSREL (10 %) of the link's flow (at least the head tolerance times
        //K), beyond which it is tested as on any other link; the connection test below allows the same B
        double Ex=0.0,sec=0.0;                                                                           //[m3/d] (Ex: exempt)
        const bool hasSec=(e<hq0.size()) && (hq0[e]!=-HUGE_VAL) && (fabs(hm-hq0[e])>_szTolH) && (fabs(hc0[e]-hq0[e])<=1e-9);
        if (hasSec){sec=-(_szxQc[e]-qq0[e])/(hm-hq0[e]);}                                                //[m2/d] local solves' secant
        const bool ws=hasSec && (sec<0.0) && !drn;
        const double B=ws?max(SZ_WSREL*fabs(_szxQc[e]),SZ_WSH*_szTolH*K):0.0;                            //[m3/d] exempt at most
        if (ws){const double x=E-Er; Ex=max(min(x,B),-B);}
        const double Et=E-Ex;
        //the connection MODFLOW solved with (hc0, not the one this update moved on to, which a plain update sets to the
        //head) must give the exchange of a converged connection (hc = h) at the current head within the head tolerance,
        //as in SZLocalSolves (otherwise hc0 can be 1.3 m from MODFLOW's head on a converged step); on every link (on a
        //wrong-sign link the lag is that drift: there beyond the same bound B)
        if ((fabs(hc0[e]-hm)>1e-12) && (_szxCell[e].C>0.0)){
          double F1,a1,a2,a3;
          const double dq=fabs(gwsz::SubgridFlux(_szxCell[e],_szxLam[e],hm,hc0[e],F1,a1,a2,a3)-gwsz::SubgridFlux(_szxCell[e],_szxLam[e],hm,hm,F1,a1,a2,a3));
          const double dl=max(dq-B,0.0)/_szxCell[e].C;
          if (dl>_szResid){_szResid=dl; _szResE=(int)e;}
        }
        const double d=fabs(Et)/K;
        if (d>_szResid){_szResid=d; _szResE=(int)e;}
      }
    }
  }
  //total exchange of the updated interface at the current heads (into the aquifer) [m3/d] (the boundary MODFLOW gets
  //reproduces it at these heads, whatever Cg)
  double Q=0.0,G=0.0,sC=0.0;
  for (size_t e=0;e<_szxLink.size();e++){Q+=_szxQc[e]; G+=fabs(_szxQc[e]); sC+=_szxFc[e]*_szxC[e];}
  //flow test: relative change, or a change no larger than the head tolerances allow on every link
  //(sum of conductances x the larger of the interface tolerance and MODFLOW's DVCLOSE)
  double dQ=fabs(Q-prevQ);
  _szDQ=dQ; _szDQtol=max(_szTolQ*max(G,1e-6),sC*max(_szTolH,(IsExternal() && (_extDvclose>0))?_extDvclose:_dvclose));
  bool conv=(_szResid<_szTolH) && (dQ<=_szDQtol) && (steady || (_szEsum<=_szDQtol));
  _szWhy=0; _szWhyE=-1;
  //nor may the boundary MODFLOW solved with take, at its current heads, more out of a column than the column holds (a
  //line steeper than the column can feed, continued below its zero level, reverses into flow out of a column without free
  //water when MODFLOW's head ends a little below where it was linearized: within the head tolerance, but with a
  //conductance of 1e4 m2/d that is water the soil zone does not have; rough suite: up to 10 m3 on a step judged converged otherwise). The
  //next update linearizes at the current head, where the line carries the local solve's flow
  if (conv && !steady){ //(with the rejected recharge booked at convergence)
    int es=-1; const double sf=SZSolvedShortfall(_szTstep,rej0,&es);
    if (sf>1e-6){conv=false; _szWhy=1; _szWhyE=es; _szWhyV=sf;}
  }
  //nor may the connection MODFLOW solved with lie further from the returned head than the smoothing depth d_s (at least 10
  //times the head tolerance) where the saturated fraction differs between them: at a converged interface hc = h, and the
  //fraction of hc is what GWSoilExchange.csv reports and what splits the water-table ET. The tests above are of flows, and
  //where none depends on the connection (an empty soil zone, no recharge to reject, columns giving nothing back) they pass
  //whatever hc is: a regula falsi bracket formed after an early iterate raised a cell's head 7.5 m can stay while the
  //head falls back (Illinois halves the stale end's residual until the bracket collapses onto it), and the interface
  //would converge with hc 2.5 m above the head, a saturated fraction of 0.52 reported over a water table below every base
  //sample (rough suite R20, LAG with :SoilZoneVadoseClasses 3). Such a link's bracket is dropped, so that the next update
  //is the plain one. (Otherwise converged connections lag the head by far less than d_s: by 1 to 2 cm on the instrumented
  //R20 runs; the quick rough suite, all methods, one and three classes, has no other run this test changes). The link's
  //last residual goes too (_szxHr), so that the next update is the plain one, not a secant or a new bracket with that
  //stale point
  if (conv && !steady){
    for (size_t e=0;e<_szxLink.size();e++){
      const gw_link &L=_aLinks[_szxLink[e]]; if (_colCHDTop[L.ic]){continue;}
      const double hx=_X[_reduced[TopNode(L.ic)]]; if (!(std::fabs(hx)<1e29)){continue;}   //(HDRY: no head to follow)
      const double hm=hx*_hf;
      if (!(fabs(hc0[e]-hm)>max(_szxCell[e].ds,10.0*_szTolH))){continue;}
      const double dF=fabs(gwsz::SatFraction(hc0[e],_szxCell[e].zbs,_szxCell[e].ds)-gwsz::SatFraction(hm,_szxCell[e].zbs,_szxCell[e].ds));
      if (dF>1e-6){
        conv=false; _szxBk[e]=0; _szxHr[e]=-HUGE_VAL;
        if (_szxDropped.size()!=_szxLink.size()){_szxDropped.assign(_szxLink.size(),0);}
        _szxDropped[e]=1; _szDropN++; if (!_szDropThis){_szDropSteps++; _szDropThis=true;}
        if ((_szWhy!=2) || (fabs(hc0[e]-hm)>fabs(_szWhyV))){
          _szWhy=2; _szWhyE=(int)e; _szWhyV=hc0[e]-hm;
          _szWhyFa=gwsz::SatFraction(hc0[e],_szxCell[e].zbs,_szxCell[e].ds); _szWhyFb=gwsz::SatFraction(hm,_szxCell[e].zbs,_szxCell[e].ds);
        }
      }
    }
  }
  prevQ=Q; _szIter++; _szConv=conv; (void)iter;
  if (conv || last){ //keep what MODFLOW solved with (the update changed it by less than the tolerances, or no solve follows:
                     //the step's flows must be those of the arrays of MODFLOW's last solve)
    _linkF=F0; _szxHc=hc0; _szxLam=lam0; _szxQc=q0; _szxCg=cg0; _szxRej=rej0; _szRej=srej0; _szxFc=fc0; _szxOne=one0; _szxRs=rs0;
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
  _szRejVol=_szOvfVol=_szAppliedVol=_szShortVol=0.0; _szDebtVol=0.0; _szNegUnheldV=0.0;
  if (_szxDebt.size()!=_szxLink.size()){_szxDebt.assign(_szxLink.size(),0.0);}
  vector<double> phi(n),tg(n);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    double A=_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
    for (int j=0;j<n;j++){tg[j]=aPhinew[k][_szSV[k][j]];}
    double ovfV=SZSyncColumns(k,&tg[0])/MM_PER_METER*A;                         //no-op unless SOIL changed since the exchange began
    double netV=0.0,shV=0.0,rejV=0.0,nuV=0.0;                                    //[m3]
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      int e=_szxOfHRU[k][jj],c=_szxCol[e]; double a=_szColA[c];
      for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
      //rejected recharge back to the column; a negative one (the part of a negative recharge the top cell could not give)
      //taken out of it, first: what the column does not hold is the "not held" part (below)
      double ovf=0.0;
      if (_szxRej[e]>=0.0){ovf=gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);}           //[mm over the column]
      else {double shR=0.0; gwsz::Apply(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER,shR); nuV+=shR/MM_PER_METER*a;}
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
      double ovf=0.0;
      if (_szcRej[e]>=0.0){ovf=gwsz::Fill(_szLay[k],phi,_szcRej[e]/a*MM_PER_METER);}
      else {double shR=0.0; gwsz::Apply(_szLay[k],phi,_szcRej[e]/a*MM_PER_METER,shR); nuV+=shR/MM_PER_METER*a;}
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
    int nc=SZNCols(k),o=_szColOff[k];
    //negative recharge not held: the cut a column did not hold is taken out of the HRU's other columns (the same Raven
    //soil store; the rest column and those over specified heads included), in proportion to their water; only what the
    //whole HRU's soil zone cannot supply is booked, and nothing when it holds more than that (no round-off booked).
    //The soil zone then holds less, at the end of the step, than the cut: a process after the import (e.g. after
    //PERC_GR4JEXCH) has moved the water on before the coupled exchange
    if (nuV>0.0){
      double Wt=0.0; for (int c=0;c<nc;c++){double w=0.0; for (int j=0;j<n;j++){w+=max(_szColPhi[(o+c)*n+j],0.0);} Wt+=w/MM_PER_METER*_szColA[o+c];}
      const double f=(Wt>0.0)?nuV/Wt:1.0;
      double taken=0.0;
      for (int c=0;c<nc;c++){
        const double a=_szColA[o+c]; if (!(a>0.0)){continue;}
        double w=0.0; for (int j=0;j<n;j++){phi[j]=_szColPhi[(o+c)*n+j]; w+=max(phi[j],0.0);}
        if (!(w>0.0)){continue;}
        if (f>=1.0){for (int j=0;j<n;j++){taken+=max(phi[j],0.0)/MM_PER_METER*a; phi[j]=min(phi[j],0.0);}}
        else {for (int j=0;j<n;j++){const double t=f*max(phi[j],0.0); phi[j]-=t; taken+=t/MM_PER_METER*a;}}
        for (int j=0;j<n;j++){_szColPhi[(o+c)*n+j]=phi[j];}
      }
      netV+=nuV-taken;                                                           //(the soil lost the cut less what it did not hold, and what the
                                                                                 // other columns gave)
      nuV=(f>=1.0)?max(nuV-taken,0.0):0.0;                                       //(booked: what the whole HRU could not supply)
    }
    //Raven's soil-zone layers: area-weighted mean of the columns (the rest column did not change)
    double sw=0.0,W=0.0;
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
    _szRejVol+=rejV; _szShortVol+=shV; _szNegUnheldV+=nuV;
    ovfV+=_szSyncOvf[k];                                                         //(sent to SURFACE_WATER when the exchange began)
    _szOvfVol+=ovfV; _szOvfV[k]=ovfV; _szNetV[k]=netV;
  }
  _szxQprev=_szxQ;                                                               //(the LAGGED source of the next step; none from a link over
  for (size_t e=0;e<_szxQprev.size();e++){if (_colCHDTop[_aLinks[_szxLink[e]].ic]){_szxQprev[e]=0.0;}} // a specified head this step, whose exchange is exact, as the hotstart writes)
  for (size_t e=0;e<_szxRepay.size() && e<_szxDebt.size();e++){_szxDebt[e]=max(_szxDebt[e]-_szxRepay[e]*tstep,0.0);}
  _negCutVol+=_szNegCutQ*tstep; _negUnheldVol+=_szNegUnheldV;
  if ((_szNegUnheldV>1e-6) && (++_negUnheldN<=10)){
    WriteWarning("CGroundwaterModel: in the step ending "+_stepEnd+": "+to_string(_szNegUnheldV)+" m3 of the negative recharge the top cells could not "
                 "supply was not held by the soil zone when it was taken back: the soil zone holds less, at the end of the step, than the cut "
                 "(a process after the import, e.g. after PERC_GR4JEXCH, moved it on before the coupled exchange); booked as 'negative recharge "
                 "the soil zone did not hold' in GWBudget.csv",false);
  }
  if ((_szShortVol>1e-6) && (++_szShortN<=10)){
    WriteWarning("CGroundwaterModel: soil-zone exchange: "+to_string(_szShortVol)+" m3 that MODFLOW took from the soil zone was not in it in the step ending "+_stepEnd,false);
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
    double shR=0.0; gwsz::Apply(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER,shR);   //(negative: taken out first, SZApply)
    double sh=0.0; gwsz::Apply(_szLay[k],phi,dV,sh);
    S+=sh/MM_PER_METER*a;                                                        //(the exchange's: shR is fixed for the step, SZApply books it)
  }
  return S;
}

//////////////////////////////////////////////////////////////////
/// \brief shortfall [m3] that the boundary MODFLOW last solved with (SZX + SZQ + SZD, its HCOF and RHS) would leave in the
///  columns at MODFLOW's current heads: SZTrialShortfall with that flow instead of the booked one
//
double CGroundwaterModel::SZSolvedShortfall(const double &tstep,const vector<double> &rej,int *emax) const
{
  if ((_szxHcof==NULL) || (_szqHcof==NULL) || (_szdHcof==NULL)){return 0.0;}
  const int n=_szNLayers; double S=0.0,big=0.0;
  vector<double> phi(n);
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    if (_colCHDTop[L.ic]){continue;}
    const double h=_X[_reduced[TopNode(L.ic)]];
    double q=(_szxHcof[e]*h-_szxRhs[e])*_qf+(_szqHcof[e]*h-_szqRhs[e])*_qf+(_szdHcof[e]*h-_szdRhs[e])*_qf; //(as booked after the solve)
    const bool inQ=(e>=_szxUp.size()) || (_szxUp[e]==TopNode(L.ic));
    if (inQ && (e<_szxCapQ.size())){q+=_szxCapQ[e];}
    if (inQ && (e<_szxRepay.size())){q+=_szxRepay[e];}
    const int c=_szxCol[e]; const double a=_szColA[c];
    const double dV=-q*tstep/a*MM_PER_METER;
    if (!(dV<0.0)){continue;}
    for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
    double shR=0.0; gwsz::Apply(_szLay[L.k],phi,rej[e]/a*MM_PER_METER,shR);     //(negative: taken out first, SZApply)
    double sh=0.0; gwsz::Apply(_szLay[L.k],phi,dV,sh);
    S+=sh/MM_PER_METER*a;                                                        //(the exchange's: shR is fixed for the step, SZApply books it)
    if ((emax!=NULL) && (sh>0.0) && (sh/MM_PER_METER*a>big)){big=sh/MM_PER_METER*a; *emax=(int)e;}
  }
  return S;
}

//////////////////////////////////////////////////////////////////
/// \brief an iterate no interface update may be made from (MODFLOW has not converged on it): a head that is not finite, or
///  a link whose top-cell head (the only heads the local solves read) lies far outside the model's elevations (the divergence
///  margin of the stop after the solve: 10 times the relief, at least 1 km) and is not where the link's connection head of
///  the last update or the last step was: more than the margin from it, and on the other side of the model or more than the
///  margin further outside (position-based). The local solves at such heads give MODFLOW absurd boundaries (GHB heads of
///  1e45 m) from which it does not recover; a far head a converged state legitimately holds (an over-pumped confined
///  cell) does not skip the updates
//
bool CGroundwaterModel::SZDivergedIterate() const
{
  for (size_t n=0;n<_reduced.size();n++){const int r=_reduced[n]; if ((r>=0) && !std::isfinite(_X[r])){return true;}}
  double zmin,zmax; const double m=FarMargin(zmin,zmax);
  for (size_t e=0;e<_szxLink.size();e++){
    const double x=_X[_reduced[TopNode(_aLinks[_szxLink[e]].ic)]];
    if (!(fabs(x)<1e29)){continue;}                                              //(HDRY)
    const double h=x*_hf;
    if (!(FarDistance(h,zmin,zmax)>m)){continue;}
    const double ref=_szxHc[e];
    if (fabs(h-ref)<=m){continue;}
    const int sh=(h>zmax)?1:((h<zmin)?-1:0),sr=(ref>zmax)?1:((ref<zmin)?-1:0);
    if ((sh!=sr) || (FarDistance(h,zmin,zmax)>FarDistance(ref,zmin,zmax)+m)){return true;}
  }
  return false;
}

//////////////////////////////////////////////////////////////////
/// \brief GWSoilExchange.csv rows of this step (HRU values area-weighted over its exchange columns; the saturated fraction
///  over all its columns on active cells, those over specified heads included, as the rejected recharge it goes with)
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
    //base, level and water table: area means over all the HRU's columns, those over specified heads included (their water
    //table is the specified head); the exchange split as GWBudget.csv books it: with the aquifer (its active cells), and
    //with specified heads (columns over them, and exchange columns whose cell is a specified head this step, with the
    //capillary rise a specified head supplies)
    const int n=_szNLayers;
    double sa=0,zb=0,lam=0,hc=0,Fk=0,qa=0,qs=0,ca=0,cs=0,xa=0,xzb=0,xlam=0;
    for (size_t j=0;j<_szxOfHRU[k].size();j++){
      int e=_szxOfHRU[k][j]; const gw_link &L=_aLinks[_szxLink[e]];
      vector<double> phi(_szColPhi.begin()+_szxCol[e]*n,_szColPhi.begin()+(_szxCol[e]+1)*n); //column after the step
      const double lev=gwsz::SubgridLevel(_szxCell[e],gwsz::SatThickness(_szLay[L.k],phi));
      sa+=L.area; zb+=L.area*_szxZb[e]; lam+=L.area*lev; hc+=L.area*WaterTable(L.ic); Fk+=L.area*_linkF[_szxLink[e]];
      xa+=L.area; xzb+=L.area*_szxZb[e]; xlam+=L.area*lev;
      double rise=0.0; if (_colCHDTop[L.ic] || (((size_t)e<_szxUpCHD.size()) && _szxUpCHD[e])){for (int jj=0;jj<n;jj++){rise+=_szxCapL[e*n+jj];}}
      if (_colCHDTop[L.ic]){if (_szxQ[e]>0){ca+=_szxQ[e];} else {cs-=_szxQ[e];}}
      else                 {if (_szxQ[e]>0){qa+=_szxQ[e];} else {qs-=_szxQ[e];}}
      cs+=rise;
    }
    if (k<(int)_szcOfHRU.size()){
      for (size_t j=0;j<_szcOfHRU[k].size();j++){
        const int e=_szcOfHRU[k][j],li=_szcLink[e]; const gw_link &L=_aLinks[li];
        vector<double> phi(_szColPhi.begin()+_szcCol[e]*n,_szColPhi.begin()+(_szcCol[e]+1)*n);
        const double z0=SZBase(L.k,L.ic);
        sa+=L.area; zb+=L.area*z0; lam+=L.area*gwsz::SubgridLevel(_szcCell[e],gwsz::SatThickness(_szLay[L.k],phi)); hc+=L.area*WaterTable(L.ic); Fk+=L.area*_linkF[li];
        if (_szcQ[e]>0){ca+=_szcQ[e];} else {cs-=_szcQ[e];}
        for (int jj=0;jj<n;jj++){cs+=_szcCapL[e*n+jj];}
      }
    }
    if (xa>0){_szLam[k]=xlam/xa-xzb/xa;} else {_szLam[k]=0.0;}
    double Wmax=0; for (size_t j=0;j<_szSV[k].size();j++){Wmax+=_szLay[k][j].phimax;}
    double W=_szW[k];
    _SZOUT<<setprecision(10)<<t<<","<<dstr<<","<<_pModel->GetHydroUnit(k)->GetHRUID()<<",";
    if (sa>0){_SZOUT<<zb/sa<<","<<lam/sa<<","<<hc/sa<<","<<Fk/sa<<",";} else {_SZOUT<<",,,,";} //(empty: an HRU with no soil-zone column)
    _SZOUT<<qa<<","<<qs<<","<<_szRej[k]/A*MM_PER_METER<<","<<_szOvfV[k]/A*MM_PER_METER<<","<<W<<","<<Wmax<<","<<_szIter<<","<<_szResid;
    if (SZChdColumns()){_SZOUT<<","<<ca<<","<<cs;}
    _SZOUT<<endl;
  }
}
