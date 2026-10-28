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
  intersection). Each link has one GHB entry (package SZX_RAVEN) on the top active cell, exchanging
     q = F(k,i) C^sz(k,i) (z_base(k,i) + lambda(k,i) - h)   into the aquifer
  (lambda: saturated thickness of the link's column; F: saturated fraction of the link, 0 while the water table is
  below the base). Recharge on saturated ground is rejected (returned to the column); column overflow goes to
  SURFACE_WATER the same step. lambda is the backward-Euler local solve of the column for the current head
  (GWSoilZone.cpp), updated inside MODFLOW's outer iterations, and MODFLOW sees the exchange as a Robin condition
  with the column's true slope, until heads and soil levels agree.
----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include "Model.h"
#include "HydroProcessABC.h"
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
    double ks=pHRU->GetSoilProps(_szNLayers-1)->hydraul_cond;
    ExitGracefullyIf(!(ks>0.0),("CGroundwaterModel: HRU "+id+": HYDRAUL_COND of SOIL["+to_string(_szNLayers-1)+"] (the bottom of the soil zone) must be positive with :SoilZone PHYSICAL").c_str(),BAD_DATA);
    ExitGracefullyIf(_szD[k]>20.0,("CGroundwaterModel: HRU "+id+": the soil zone SOIL[0.."+to_string(_szNLayers-1)+"] is "+to_string(_szD[k])+
      " m thick; :SoilZone PHYSICAL is meant for physical soil layers (conceptual reservoirs of tens of metres are not a soil zone)").c_str(),BAD_DATA);
    if (_szD[k]>5.0){nThick++;}
    Dmax=max(Dmax,_szD[k]);
    _szOn[k]=1; nOn++;
  }
  if (nThick>0){WriteWarning("CGroundwaterModel: "+to_string(nThick)+" HRUs have a soil zone thicker than 5 m (up to "+to_string(Dmax)+" m); check that SOIL[0.."+to_string(_szNLayers-1)+"] are physical soil layers",false);}
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
  _szW.assign(nHRUs,0.0); _szLam.assign(nHRUs,0.0); _szRchV.assign(nHRUs,0.0); _szRej.assign(nHRUs,0.0); _szOvfV.assign(nHRUs,0.0); _szNetV.assign(nHRUs,0.0);
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
  _linkF.assign(_aLinks.size(),0.0); _linkFp.assign(_aLinks.size(),0.0); _linkZb.assign(_aLinks.size(),0.0); _linkW.assign(_aLinks.size(),_szOmega); _linkDF.assign(_aLinks.size(),0.0);
  _colCoverNS.assign(_ncpl,0.0);
  _szxLink.clear(); _szxC.clear(); _szxZb.clear(); _szxQ.clear(); _szxDl.clear(); _szxOfHRU.assign(nHRUs,vector<int>());
  _cellDem.assign(_ncpl,vector<double>());
  int nRatio=0; double rmin=1e300,rmax=0;
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
    _szxLink.push_back((int)li); _szxC.push_back(C); _szxZb.push_back(SZBase(L.k,ic)); _szxQ.push_back(0.0); _szxDl.push_back(0.0);
    //storage consistency: specific yield of the top aquifer layer vs drainable porosity of the bottom soil layer
    double Sd=(Lb.phimax-Lb.phifc)/(Lb.H*MM_PER_METER);
    double r=C1.Sy/Sd; rmin=min(rmin,r); rmax=max(rmax,r);
    if ((r<0.5) || (r>2.0)){nRatio++;}
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
  }
  if (nRatio>0){
    WriteWarning("CGroundwaterModel: :SoilZone PHYSICAL: on "+to_string(nRatio)+" links the specific yield of the top aquifer layer differs from the drainable porosity "
                 "n(1-SF)(1-FIELD_CAPACITY) of the bottom soil-zone layer by more than a factor 2 (ratio "+to_string(rmin)+" to "+to_string(rmax)+
                 "); the water table then rises at a different rate in the soil than below it",false);
  }
  _szRatioMin=(rmin<1e299)?rmin:0.0; _szRatioMax=rmax;
  //soil-zone columns (plan section 5.3): the soil zone of an HRU is split into one column per SZX entry and the rest
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
  _szColPhi.assign(_szColA.size()*_szNLayers,0.0); _szColInit=false;
  _szxLam.assign(_szxLink.size(),0.0); _szxRej.assign(_szxLink.size(),0.0);
}

//////////////////////////////////////////////////////////////////
/// \brief brings the columns of HRU k to its soil-zone water target[j] (layers bottom up) [mm]: Raven's processes
///  changed the HRU totals since the columns were last set (gains shared by room, losses by water; gwsz::Distribute).
///  On the first call the columns start from the hotstart file or, without one, uniform.
//
void CGroundwaterModel::SZSyncColumns(const int k,const double *target)
{
  int nc=SZNCols(k),n=_szNLayers,o=_szColOff[k];
  vector<vector<double> > cols(nc,vector<double>(n));
  vector<double> w(nc),tg(target,target+n);
  double sw=0.0; for (int c=0;c<nc;c++){w[c]=_szColA[o+c]; sw+=w[c];}
  for (int c=0;c<nc;c++){w[c]=(sw>0.0)?w[c]/sw:1.0/nc;}
  for (int c=0;c<nc;c++){for (int j=0;j<n;j++){cols[c][j]=_szColPhi[(o+c)*n+j];}}
  if (!_szColInit){
    bool hot=false;
    map<long long,vector<vector<double> > >::const_iterator it=_hotSzCols.find(_pModel->GetHydroUnit(k)->GetHRUID());
    if ((it!=_hotSzCols.end()) && ((int)it->second.size()==nc)){
      hot=true; for (int c=0;c<nc;c++){if ((int)it->second[c].size()!=n){hot=false;}}
      if (hot){for (int c=0;c<nc;c++){cols[c]=it->second[c];}}
    }
    if (!hot){for (int c=0;c<nc;c++){cols[c]=tg;}}
  }
  gwsz::Distribute(_szLay[k],cols,w,tg);
  for (int c=0;c<nc;c++){for (int j=0;j<n;j++){_szColPhi[(o+c)*n+j]=cols[c][j];}}
}

//////////////////////////////////////////////////////////////////
/// \brief saturated fraction of every soil-zone link from the current heads (plan section 5.4)
//
void CGroundwaterModel::SZLinkFractions(const bool relax)
{
  vector<double> zbs;
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    double h=_X[_reduced[TopNode(L.ic)]];
    const vector<double> &dem=_cellDem[L.ic];
    zbs.resize(dem.size());
    for (size_t s=0;s<dem.size();s++){zbs[s]=dem[s]-_szD[L.k];}
    double ds=_aProfiles[_colProfile[L.ic]].seep_smooth_depth;
    double F=gwsz::SatFraction(h,zbs,ds);
    double Fp=gwsz::SatFractionSlope(h,zbs,ds);
    _linkZb[_szxLink[e]]=min(gwsz::ConnectedBase(h,zbs,ds),h); //h itself while nothing is connected (only a relaxed F>0 then)
    //F depends steeply on the head inside the smoothing band: inside the iteration it is under-relaxed per link, the
    //relaxation halved (down to 1/16) whenever the update of that link reverses direction (an oscillation)
    int li=_szxLink[e];
    if (relax){
      double dF=F-_linkF[li];
      if (dF*_linkDF[li]<0.0){_linkW[li]=max(0.5*_linkW[li],1.0/16.0);}
      _linkDF[li]=dF;
      F=_linkF[li]+_linkW[li]*dF; Fp*=_linkW[li];
    }
    _linkF[li]=F; _linkFp[li]=Fp;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief recharge to MODFLOW with rejection on saturated ground (plan section 5.6): _Qcol, RCH rates, _szxRej, _szRej
//
void CGroundwaterModel::SZRecharge(const double &tstep,const bool write)
{
  std::fill(_Qcol.begin(),_Qcol.end(),0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; double V=_szRchV[k],rej=0.0;
    if (_hruRchArea[k]<=0){if (_szOn[k]){_szRej[k]=0.0;} continue;} //all on specified-head cells (counted in _QrchCHD)
    for (size_t j=0;j<_hruLinks[k].size();j++){
      int li=_hruLinks[k][j]; const gw_link &L=_aLinks[li];
      if (_colCHDTop[L.ic]){continue;}
      double share=L.area/_hruRchArea[k],F=_linkF[li];
      _Qcol[L.ic]+=V/tstep*share*(1.0-F);
      rej+=V*share*F;
      if (_linkSzx[li]>=0){_szxRej[_linkSzx[li]]=V*share*F;}
    }
    if (_szOn[k]){_szRej[k]=rej;}
  }
  if (write){for (size_t i=0;i<_bndCol.size();i++){_rchBound[i]=_Qcol[_bndCol[i]]/MF6CellArea(_bndCol[i])*_rf;}}
}

//////////////////////////////////////////////////////////////////
/// \brief local implicit solve of every exchange column for the current heads and fractions (plan section 6.1);
///  a column exchanges with its own cell only, so the slope of lambda against that head is exact
//
void CGroundwaterModel::SZLocalSolves(const double &tstep)
{
  vector<gwsz::sz_link> links(1); vector<double> phi;
  double maxd=0.0; const int n=_szNLayers;
  for (size_t e=0;e<_szxLink.size();e++)
  {
    const gw_link &L=_aLinks[_szxLink[e]]; int k=L.k,c=_szxCol[e];
    double a=_szColA[c];                                                         //[m2] column area
    phi.assign(_szColPhi.begin()+c*n,_szColPhi.begin()+(c+1)*n);
    gwsz::Fill(_szLay[k],phi,_szxRej[e]/a*MM_PER_METER);                         //rejected recharge back to the column
    links[0].Ceff=_linkF[_szxLink[e]]*_szxC[e]; links[0].zb=_linkZb[_szxLink[e]]; links[0].h=_X[_reduced[TopNode(L.ic)]];
    double dV,ovf,lamp;
    double lam=gwsz::LocalSolve(_szLay[k],phi,links,a,tstep,dV,ovf,&lamp);
    maxd=max(maxd,fabs(lam-_szxLam[e]));
    _szxLam[e]=lam;
    //response of lambda to the head (implicit function theorem on g(dV;h)=0): dlambda/dh = c Ceff lambda'/(1 + c Ceff lambda'),
    //c = 1000 dt/a [mm per m3], lambda' = dlambda/dV [m/mm]
    double x=MM_PER_METER*tstep/a*links[0].Ceff*lamp;
    _szxDl[e]=x/(1.0+x);
  }
  _szResid=maxd;
}

//////////////////////////////////////////////////////////////////
/// \brief SZX heads and conductances, RCH rates (with rejection) and EVT rates (outside saturated ground)
//
void CGroundwaterModel::SZSetArrays(const bool steady)
{
  //Each link exchanges q = F C^sz (z_b + lambda - h) with the aquifer (into the aquifer), lambda from its column's
  //local solve. MODFLOW is given this flow as a linear boundary q = Cg (Hg - h) that reproduces it at the current head
  //and has the slope the soil zone really has: dq/dh = -F C^sz (1 - dlambda/dh) + C^sz dF/dh (z_b + lambda - h). The
  //first term is the column's storage seen from the aquifer (its level follows the head), the second the growth of the
  //saturated fraction, kept only where it adds conductance (water moving up into the soil). Hg - h = q/Cg, so at a
  //converged interface the flow MODFLOW books is exactly the flow of the local solve; Cg sets only how fast the
  //interface converges (a Robin interface condition; with a fixed soil level it converges slowly on stiff contacts).
  //MODFLOW refuses a GHB head below the cell bottom: there the head is raised to the bottom and SZQ_RAVEN (a rate)
  //carries Cg (Hg - bottom), which leaves the boundary q = Cg (Hg - h) unchanged for every h.
  int st=max(_szxNcol,1),sq=max(_szqNcol,1);
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]]; int li=_szxLink[e];
    double hb=_linkZb[li]+(steady?0.0:_szxLam[e]);
    double Ce=_linkF[li]*_szxC[e],Hg=hb,Cg=Ce;
    if (!steady){
      double h=_X[_reduced[TopNode(L.ic)]];
      Cg=Ce*(1.0-_szxDl[e])+max(0.0,_linkFp[li]*_szxC[e]*(h-hb));
      if (Cg>1e-12*_szxC[e]){Hg=h+Ce*(hb-h)/Cg;} else {Cg=Ce;}
    }
    double bot=_botm[UserNode(0,L.ic)];
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
bool CGroundwaterModel::SZIterate(double **aPhinew,const int iter,const bool steady,double &prevQ)
{
  (void)aPhinew;
  if (iter>=_szMaxIter){return true;} //update budget used: values frozen, MODFLOW finishes the step
  //state MODFLOW just solved with: kept if this update shows the interface has converged, so that MODFLOW's budget
  //(computed from the arrays it solved with) and the flows Raven books stay exactly consistent
  vector<double> F0(_linkF),Fp0(_linkFp),Zb0(_linkZb),lam0(_szxLam),Dl0(_szxDl),rej0(_szxRej),srej0(_szRej);
  SZLinkFractions(true);
  if (!steady){
    SZRecharge(_szTstep,false);
    SZLocalSolves(_szTstep);
  }
  //total exchange of the updated interface at the current heads (into the aquifer) [m3/d] (the boundary MODFLOW gets
  //reproduces it at these heads, whatever Cg)
  double Q=0.0,G=0.0,sC=0.0;
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]];
    double Ce=_linkF[_szxLink[e]]*_szxC[e];
    double q=Ce*(_linkZb[_szxLink[e]]+(steady?0.0:_szxLam[e])-_X[_reduced[TopNode(L.ic)]]);
    Q+=q; G+=fabs(q); sC+=Ce;
  }
  //flow test: relative change, or a change no larger than the head tolerances allow on every link
  //(sum of conductances x the larger of the interface tolerance and MODFLOW's DVCLOSE)
  double dQ=fabs(Q-prevQ);
  if (steady){_szResid=0.0;}
  bool conv=(steady || (_szResid<_szTolH)) && ((dQ<=_szTolQ*max(G,1e-6)) || (dQ<=sC*max(_szTolH,_dvclose)));
  prevQ=Q; _szIter=iter; _szConv=conv;
  if (conv){ //keep what MODFLOW solved with (the update changed it by less than the tolerances)
    _linkF=F0; _linkFp=Fp0; _linkZb=Zb0; _szxLam=lam0; _szxDl=Dl0; _szxRej=rej0; _szRej=srej0;
    if (!steady){SZRecharge(_szTstep,false);} //_Qcol of those fractions
  }
  else {SZSetArrays(steady);}
  return conv;
}

//////////////////////////////////////////////////////////////////
/// \brief applies MODFLOW's booked exchange (_szxQ, after the non-convergence guard), the rejected recharge and the
///  overflow to the soil-zone columns, Raven's soil-zone layers (their area-weighted mean) and SURFACE_WATER
///  (plan sections 5.2, 7)
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
    SZSyncColumns(k,&tg[0]);                                                    //no-op unless SOIL changed since the snapshot
    double ovfV=0.0,netV=0.0,shV=0.0,rejV=0.0;                                  //[m3]
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
    _szRejVol+=rejV; _szOvfVol+=ovfV; _szShortVol+=shV;
    _szOvfV[k]=ovfV; _szNetV[k]=netV;
  }
  if (_szShortVol>1e-6){
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
      sa+=L.area; zb+=L.area*_szxZb[e]; lam+=L.area*_szxLam[e]; hc+=L.area*WaterTable(L.ic); Fk+=L.area*_linkF[_szxLink[e]];
      if (_szxQ[e]>0){qa+=_szxQ[e];} else {qs-=_szxQ[e];}
    }
    if (sa>0){zb/=sa; lam/=sa; hc/=sa; Fk/=sa;}
    _szLam[k]=lam;
    double Wmax=0; for (size_t j=0;j<_szSV[k].size();j++){Wmax+=_szLay[k][j].phimax;}
    double W=_szW[k];
    _SZOUT<<setprecision(10)<<t<<","<<dstr<<","<<_pModel->GetHydroUnit(k)->GetHRUID()<<","<<zb<<","<<zb+lam<<","<<hc<<","<<Fk<<","
          <<qa<<","<<qs<<","<<_szRej[k]/A*MM_PER_METER<<","<<_szOvfV[k]/A*MM_PER_METER<<","<<W<<","<<Wmax<<","<<_szIter<<","<<_szResid<<endl;
  }
}
