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
  return _land[ic]-_szD[k];
}

//////////////////////////////////////////////////////////////////
/// \brief one SZX entry per HRU-cell link of a soil-zone HRU (top active cell, not on specified-head cells),
///  conductances (half the bottom soil layer and half the top aquifer cell in series), DEM samples per cell
//
void CGroundwaterModel::SZBuildEntries()
{
  int nHRUs=_pModel->GetNumHRUs();
  _linkF.assign(_aLinks.size(),0.0);
  _colCoverNS.assign(_ncpl,0.0);
  _szxLink.clear(); _szxC.clear(); _szxZb.clear(); _szxQ.clear(); _szxCell.clear(); _szxOfHRU.assign(nHRUs,vector<int>());
  _cellDem.assign(_ncpl,vector<double>());
  int nRatio=0,nMix=0; double rmin=1e300,rmax=0,dmix=0.0;
  for (size_t li=0;li<_aLinks.size();li++)
  {
    const gw_link &L=_aLinks[li];
    if (!_szOn[L.k]){_colCoverNS[L.ic]+=L.area; continue;}
    if (_colCHDTop[L.ic]){continue;}
    int ic=L.ic;
    const gw_profile &P=_aProfiles[_colProfile[ic]];
    const gw_aquifer_class &C1=_aClasses[P.layers[0].iclass];
    double b1=_top[ic]-_botm[UserNode(0,ic)];                                    //[m] top cell
    const gwsz::sz_layer &Lb=_szLay[L.k][0];                                      //bottom soil-zone layer
    double ks=_pModel->GetHydroUnit(L.k)->GetSoilProps(_szNLayers-1)->hydraul_cond/MM_PER_METER; //[m/d]
    double res=Lb.H/(2.0*ks)+b1/(2.0*C1.Kv);                                     //[d]
    double C=L.area/res;                                                          //[m2/d]
    _szxOfHRU[L.k].push_back((int)_szxLink.size());
    _szxLink.push_back((int)li); _szxC.push_back(C); _szxZb.push_back(SZBase(L.k,ic)); _szxQ.push_back(0.0);
    //storage consistency: specific yield of the top aquifer layer vs drainable porosity of the bottom soil layer
    double Sd=(Lb.phimax-Lb.phifc)/(Lb.H*MM_PER_METER);
    double r=C1.Sy/Sd; rmin=min(rmin,r); rmax=max(rmax,r);
    if ((r<0.5) || (r>2.0)){nRatio++;}
    //the model top is the overlap-weighted soil-zone base of the cell: an HRU whose soil zone is thicker or thinner
    //than the cell's mean (cells shared with lakes or other soil profiles) has its base below or above the top
    double dm=fabs(_szD[L.k]-(_land[ic]-_top[ic]));
    if (dm>max(0.05,0.1*_szD[L.k])){nMix++; dmix=max(dmix,dm);}
    //DEM sample points of the cell (sub-grid saturated fraction)
    if (_cellDem[ic].empty()){
      if (_dem.IsLoaded() && (_szSubgrid>1)){
        vector<gw_pt> pts; _grid.SamplePoints(ic,_szSubgrid,pts);
        for (size_t s=0;s<pts.size();s++){
          double X,Y,lon,lat,v; _grid.ToWorld(pts[s].x,pts[s].y,X,Y); _proj.Inverse(X,Y,lon,lat);
          if (_dem.Sample(lon,lat,v)){_cellDem[ic].push_back(v);}
        }
      }
      if (_cellDem[ic].empty()){_cellDem[ic].push_back(_land[ic]);}
    }
    gwsz::sz_cell cell;                                                           //sub-grid closure of the column (GWSoilZone.h)
    for (size_t q=0;q<_cellDem[ic].size();q++){cell.zbs.push_back(_cellDem[ic][q]-_szD[L.k]);}
    std::sort(cell.zbs.begin(),cell.zbs.end());
    cell.ds=P.seep_smooth_depth; cell.D=_szD[L.k]; cell.C=C;
    _szxCell.push_back(cell);
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
  _szColOff.assign(nHRUs,-1); _szColA.clear(); _szxCol.assign(_szxLink.size(),-1);
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
    _szColA.push_back(max(A-sum,0.0));                                            //rest of the HRU
  }
  //the conductance of a link acts on its column: on the HRU's AREA share, as the column's water does (with the
  //overlap area instead, an HRU whose polygon is larger than its AREA would exchange its soil water that much faster)
  for (size_t e=0;e<_szxLink.size();e++){
    double la=_aLinks[_szxLink[e]].area,sc=(la>0.0)?_szColA[_szxCol[e]]/la:1.0;
    _szxC[e]*=sc; _szxCell[e].C*=sc;
  }
  { double dsmin=1e300; for (size_t e=0;e<_szxCell.size();e++){dsmin=min(dsmin,_szxCell[e].ds);}
    if ((_szxCell.size()>0) && (dsmin<0.1)){
      WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: SEEPAGE_SMOOTHING_DEPTH of the cells' profiles is as small as "+to_string(dsmin)+
                   " m; the soil-zone contact then switches almost abruptly and the interface may not converge (0.1 to 0.5 m is typical)",false);}
  }
  _szColPhi.assign(_szColA.size()*_szNLayers,0.0); _szColInit=false;
  _szxLam.assign(_szxLink.size(),0.0); _szxRej.assign(_szxLink.size(),0.0); _szxFw.assign(_szxLink.size(),0.0); _szxFc.assign(_szxLink.size(),0.0);
  _szxQc.assign(_szxLink.size(),0.0); _szxCg.assign(_szxLink.size(),0.0);
  size_t ne=_szxLink.size();
  _szxHc.assign(ne,0.0); _szxHp.assign(ne,0.0); _szxHr.assign(ne,-HUGE_VAL); _szxHq.assign(ne,-HUGE_VAL); _szxQq.assign(ne,0.0);
  _szxBa.assign(ne,0.0); _szxBra.assign(ne,0.0); _szxBb.assign(ne,0.0); _szxBrb.assign(ne,0.0); _szxBk.assign(ne,0);
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
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      int e=_szxOfHRU[k][jj],c=_szxCol[e]; const gw_link &L=_aLinks[_szxLink[e]];
      double h0=_X[_reduced[TopNode(L.ic)]];
      const gwsz::sz_cell &cell=_szxCell[e];
      _szxHc[e]=h0;
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
    double h=_X[_reduced[TopNode(L.ic)]];
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
}

//////////////////////////////////////////////////////////////////
/// \brief recharge to MODFLOW with rejection on saturated ground (design note section 6.1): _Qcol, RCH rates, _szxRej, _szRej
//
void CGroundwaterModel::SZRecharge(const double &tstep,const bool write)
{
  std::fill(_Qcol.begin(),_Qcol.end(),0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; double V=_szRchV[k],rej=0.0;
    if (_hruRchArea[k]<=0){if (_szOn[k]){_szRej[k]=0.0;} continue;} //all on specified-head cells (counted in _QrchCHD)
    //an HRU with a soil zone: recharge leaves the columns that hold free water in the bottom layer (Raven's drainage acts
    //on the HRU mean, and its loss is taken from the columns' free water first, gwsz::Distribute), not by area: drainage
    //from saturated columns spread over unsaturated ones would recharge the aquifer upslope of where it came from
    double W=0.0;
    if (_szOn[k]){for (size_t j=0;j<_hruLinks[k].size();j++){int li=_hruLinks[k][j]; if (!_colCHDTop[_aLinks[li].ic] && (_linkSzx[li]>=0)){W+=_szxFw[_linkSzx[li]];}}}
    for (size_t j=0;j<_hruLinks[k].size();j++){
      int li=_hruLinks[k][j]; const gw_link &L=_aLinks[li];
      if (_colCHDTop[L.ic]){continue;}
      double share=(W>0.0)?((_linkSzx[li]>=0)?_szxFw[_linkSzx[li]]/W:0.0):L.area/_hruRchArea[k],F=_linkF[li];
      _Qcol[L.ic]+=V/tstep*share*(1.0-F);
      rej+=V*share*F;
      if (_linkSzx[li]>=0){_szxRej[_linkSzx[li]]=V*share*F;}
    }
    if (_szOn[k]){_szRej[k]=rej;}
  }
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
    double h=_X[_reduced[TopNode(L.ic)]];
    double a=_szColA[c];                                                         //[m2] column area
    double lam=0.0,lamp=0.0;
    if (!steady){
      phi.assign(_szColPhi.begin()+c*n,_szColPhi.begin()+(c+1)*n);
      gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);                       //rejected recharge back to the column
      double dV,ovf;
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
    //update budget): the boundary MODFLOW keeps may not drain more than the column's free water over the step whatever
    //the head does (down to the cell bottom), as the river and store exchanges are limited. The limit wins over the
    //secant floor for water rising into the soil: with q < 0 a head falling below Hg reverses the linear boundary, and
    //the floor let it drain columns that held no free water (water created). The cost: such a boundary
    //can let more water rise than the contact allows (it fills the column and overflows; no water is created), and
    //keeps extracting when the head falls below every base. MODFLOW needs iterations to settle on that: with fewer
    //than 5 solves left (a very low outer-iteration limit) the floor is kept, as before.
    //(Iterated, the converged flow is the local solve's, which respects both; the limit would only make the slope
    //too small. Limiting a boundary MODFLOW already iterates on, after the fact, left MODFLOW unconverged.)
    if (!steady && ((_szMaxIter==0) || (SZFreezeAt()<1) || _szCapUpd)){
      double avail=gwsz::StackFreeWater(_szLay[k],phi)/MM_PER_METER*a/tstep;           //[m3/d]
      double bot=_botm[UserNode(0,L.ic)];
      Cg=min(Cg,max(avail-q,0.0)/max(h-bot,1e-3));
      if ((q<0.0) && !_szCapWins){Cg=max(Cg,sec);}
      Cg=max(Cg,1e-9*cell.C);
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
  int st=max(_szxNcol,1),sq=max(_szqNcol,1);
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    double h=_X[_reduced[TopNode(L.ic)]];
    double bot=_botm[UserNode(0,L.ic)]+1e-3;                                    //[m] margin: MODFLOW reads the bottom rounded
    double Cg=_szxCg[e];
    if (Cg<=1e-6*_szxC[e]){ //exchange (nearly) independent of the head: a rate (no GHB head far outside the cell)
      _szxHead[e*st]=bot; _szxCond[e*st]=0.0; _szqQ[e*sq]=_szxQc[e]; continue;
    }
    double Hg=h+_szxQc[e]/Cg;
    double H=max(Hg,bot);
    _szxHead[e*st]=H; _szxCond[e*st]=Cg; _szqQ[e*sq]=Cg*(Hg-H);
  }
  if (steady){return;}
  SZRecharge(_szTstep,true);
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
  vector<double> F0(_linkF),hc0(_szxHc),lam0(_szxLam),q0(_szxQc),cg0(_szxCg),rej0(_szxRej),srej0(_szRej),fc0(_szxFc);
  SZLinkFractions(true,mf6conv);
  if (!steady){SZRecharge(_szTstep,false);}
  _szCapUpd=(iter==freezeAt) || (_szExtra>=_szMaxIter); //the last update possible (iterations or budget): its boundary stays,
                                                        //supply-limited as in explicit mode
  _szCapWins=(_maxiter-iter>=5);
  SZLocalSolves(_szTstep,steady);
  _szCapUpd=false;
  //total exchange of the updated interface at the current heads (into the aquifer) [m3/d] (the boundary MODFLOW gets
  //reproduces it at these heads, whatever Cg)
  double Q=0.0,G=0.0,sC=0.0;
  for (size_t e=0;e<_szxLink.size();e++){Q+=_szxQc[e]; G+=fabs(_szxQc[e]); sC+=_szxFc[e]*_szxC[e];}
  //flow test: relative change, or a change no larger than the head tolerances allow on every link
  //(sum of conductances x the larger of the interface tolerance and MODFLOW's DVCLOSE)
  double dQ=fabs(Q-prevQ);
  _szDQ=dQ; _szDQtol=max(_szTolQ*max(G,1e-6),sC*max(_szTolH,_dvclose));
  bool conv=(_szResid<_szTolH) && (dQ<=_szDQtol);
  prevQ=Q; _szIter++; _szConv=conv; (void)iter;
  if (conv || last){ //keep what MODFLOW solved with (the update changed it by less than the tolerances, or no solve follows:
                     //the step's flows must be those of the arrays of MODFLOW's last solve)
    _linkF=F0; _szxHc=hc0; _szxLam=lam0; _szxQc=q0; _szxCg=cg0; _szxRej=rej0; _szRej=srej0; _szxFc=fc0;
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
  _szRejVol=_szOvfVol=_szAppliedVol=_szShortVol=0.0;
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
      for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
      ovfV+=ovf/MM_PER_METER*a; netV+=(dV+sh)/MM_PER_METER*a; shV+=sh/MM_PER_METER*a; rejV+=_szxRej[e];
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
    _szAppliedVol+=rejV+netV;                                                    //net water added to Raven (overflow stays in Raven)
    _szRejVol+=rejV; _szShortVol+=shV;
    ovfV+=_szSyncOvf[k];                                                         //(sent to SURFACE_WATER when the exchange began)
    _szOvfVol+=ovfV; _szOvfV[k]=ovfV; _szNetV[k]=netV;
  }
  if ((_szShortVol>1e-6) && (++_szShortN<=10)){
    WriteWarning("CGroundwaterModel: soil-zone exchange: "+to_string(_szShortVol)+" m3 that MODFLOW took from the soil zone was not in it (step "+to_string(_stepCount)+")",false);
  }
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
