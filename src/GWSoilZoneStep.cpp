/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWSoilZoneStep.cpp
  :SoilZone PHYSICAL, before Raven's processes of a step (soil_zone_columns_plan.md sections 3 and 4):
  - :SoilZoneExchangeTiming LAGGED: the columns' exchange of the last step is given to them (and to Raven's soil-zone
    layers) before Raven's processes, so that evaporation and percolation see the water that rises in; it is taken
    back out of the columns when the exchange begins, and the step's exchange is booked as before (MODFLOW's flow).
  - :SoilZoneClasses K: the columns of each HRU are grouped into up to K classes by wetness; Raven's processes run once
    per class with the class's soil water (Solvers.cpp), and each class's result goes to its own columns.
  ----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include "Model.h"

//////////////////////////////////////////////////////////////////
/// \brief before Raven's processes of the step (called from MassEnergyBalance)
//
void CGroundwaterModel::SZStepBegin(double **aPhinew,double **aPhi,const double &tstep,const optStruct &Options)
{
  const int nH=_pModel->GetNumHRUs();
  _szCls.assign(nH,vector<int>()); _szClsDone.assign(nH,0);
  _szxLagL.assign(_szxLink.size(),0.0); _szxLagD.assign(_szxLink.size(),0.0); _szLagNet.assign(nH,0.0);
  if (!SZPhysical() || !_szColInit){return;}
  if (_szTiming==1){SZLagApply(aPhinew,aPhi,tstep);}
  if (_szClasses>1){
    if (Options.sol_method==ORDERED_SERIES){SZClassify();}
    else if (!_szClsWarned){
      _szClsWarned=true;
      WriteWarning("CGroundwaterModel: :SoilZoneClasses needs :Method ORDERED_SERIES; Raven's processes act on the HRU's mean soil water",false);
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief column classes of every soil-zone HRU: the columns sorted by saturated thickness, then total soil water,
///  are split into up to K groups of about equal area; classes whose soil water is the same are merged
//
void CGroundwaterModel::SZClassify()
{
  const int n=_szNLayers,nH=_pModel->GetNumHRUs();
  _szClsW.assign(nH,vector<double>()); _szClsPhi.assign(nH,vector<double>());
  _szClsRes.assign(nH,vector<double>()); _szClsGW.assign(nH,vector<double>());
  vector<double> phi(n);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    const int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const int nc=SZNCols(k),o=_szColOff[k],K=min(_szClasses,nc);
    if (K<=1){continue;}
    vector<pair<pair<double,double>,int> > key(nc);
    double A=0.0;
    for (int c=0;c<nc;c++){
      double W=0.0; for (int j=0;j<n;j++){phi[j]=_szColPhi[(o+c)*n+j]; W+=phi[j];}
      key[c]=make_pair(make_pair(gwsz::SatThickness(_szLay[k],phi),W),c);
      A+=_szColA[o+c];
    }
    if (!(A>0.0)){continue;}
    std::sort(key.begin(),key.end());
    vector<int> cls(nc,0); double cum=0.0;
    for (int r=0;r<nc;r++){
      const int c=key[r].second; const double a=_szColA[o+c];
      cls[c]=min((int)((cum+0.5*a)/A*K),K-1); cum+=a;
    }
    //class means; drop empty classes and merge classes with the same soil water (the processes would give the same result)
    vector<double> w(K,0.0),m(K*n,0.0);
    for (int c=0;c<nc;c++){const double a=_szColA[o+c]; w[cls[c]]+=a; for (int j=0;j<n;j++){m[cls[c]*n+j]+=a*_szColPhi[(o+c)*n+j];}}
    vector<int> remap(K,-1); vector<double> W2,M2;
    for (int q=0;q<K;q++){
      if (!(w[q]>0.0)){continue;}
      for (int j=0;j<n;j++){m[q*n+j]/=w[q];}
      int same=-1;
      for (size_t p=0;p<W2.size();p++){
        bool eq=true; for (int j=0;j<n;j++){if (fabs(M2[p*n+j]-m[q*n+j])>1e-9*max(1.0,fabs(m[q*n+j]))){eq=false; break;}}
        if (eq){same=(int)p; break;}
      }
      if (same>=0){ //merge into an earlier class (the means are equal: the merged mean is the same)
        remap[q]=same; W2[same]+=w[q];
      }
      else {remap[q]=(int)W2.size(); W2.push_back(w[q]); for (int j=0;j<n;j++){M2.push_back(m[q*n+j]);}}
    }
    if (W2.size()<=1){continue;}
    for (int c=0;c<nc;c++){cls[c]=remap[cls[c]];}
    for (size_t p=0;p<W2.size();p++){W2[p]/=A;}
    _szCls[k]=cls; _szClsW[k]=W2; _szClsPhi[k]=M2;
    _szClsRes[k].assign(M2.size(),0.0); _szClsGW[k].assign(W2.size(),0.0);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief class c's soil water into the state vector aPhi (a copy of HRU k's)
//
void CGroundwaterModel::SZClassSoil(const int k,const int c,double *aPhi) const
{
  const int n=_szNLayers;
  for (int j=0;j<n;j++){aPhi[_szSV[k][j]]=_szClsPhi[k][c*n+j];}
}

//////////////////////////////////////////////////////////////////
/// \brief class c's state after Raven's processes: its soil water and its recharge (GROUNDWATER pass-through)
//
void CGroundwaterModel::SZClassResult(const int k,const int c,const double *aPhi)
{
  const int n=_szNLayers;
  for (int j=0;j<n;j++){_szClsRes[k][c*n+j]=aPhi[_szSV[k][j]];}
  _szClsGW[k][c]=(_iGW>=0)?aPhi[_iGW]:0.0;
  if (c==(int)_szClsW[k].size()-1){_szClsDone[k]=1;}
}

//////////////////////////////////////////////////////////////////
/// \brief brings each class's columns to the class's result (Exchange, in place of SZSyncColumns): Raven's changes
///  to HRU k's soil zone after the class runs (lateral processes) are shared by area among all columns; returns the
///  overflow [mm over the HRU], the recharge loss of every column in lossR [mm over the column]. Falls back to the
///  HRU-wide sync when a class target would be negative.
//
double CGroundwaterModel::SZSyncClasses(const int k,const double *target,const double R,vector<double> &lossR)
{
  const int n=_szNLayers,nc=SZNCols(k),o=_szColOff[k],K=(int)_szClsW[k].size();
  vector<double> d(n,0.0); double dR=R;
  for (int j=0;j<n;j++){d[j]=target[j]; for (int c=0;c<K;c++){d[j]-=_szClsW[k][c]*_szClsRes[k][c*n+j];}}
  for (int c=0;c<K;c++){dR-=_szClsW[k][c]*max(_szClsGW[k][c],0.0);}
  //(round-off of the class runs and of the lateral share, about 1e-15 mm, must not decide whether the classes apply:
  // targets a hair below zero are zero)
  const double tol=1e-9;                                                          //[mm]
  for (int c=0;c<K;c++){
    for (int j=0;j<n;j++){if (_szClsRes[k][c*n+j]+d[j]<-tol){return SZSyncColumns(k,target,R,&lossR);}}
    if (max(_szClsGW[k][c],0.0)+dR<-tol){return SZSyncColumns(k,target,R,&lossR);}
  }
  lossR.assign(nc,0.0);
  double ovf=0.0,mvT=0.0;
  const bool byArea=(k<(int)_szRchByArea.size()) && _szRchByArea[k];
  for (int c=0;c<K;c++){
    vector<int> idx; double sw=0.0;
    for (int q=0;q<nc;q++){if (_szCls[k][q]==c){idx.push_back(q); sw+=_szColA[o+q];}}
    if (idx.empty()){continue;}
    vector<vector<double> > cols(idx.size(),vector<double>(n)); vector<double> w(idx.size()),tg(n),lr;
    for (size_t i=0;i<idx.size();i++){
      w[i]=(sw>0.0)?_szColA[o+idx[i]]/sw:1.0/idx.size();
      for (int j=0;j<n;j++){cols[i][j]=_szColPhi[(o+idx[i])*n+j];}
    }
    for (int j=0;j<n;j++){tg[j]=max(_szClsRes[k][c*n+j]+d[j],0.0);}
    double mv=0.0;
    double ov=gwsz::Distribute(_szLay[k],cols,w,tg,0,max(max(_szClsGW[k][c],0.0)+dR,0.0),byArea,&lr,&mv);
    ovf+=ov*_szClsW[k][c]; mvT+=mv*_szClsW[k][c];
    for (size_t i=0;i<idx.size();i++){
      for (int j=0;j<n;j++){_szColPhi[(o+idx[i])*n+j]=cols[i][j];}
      lossR[idx[i]]=(i<lr.size())?lr[i]:0.0;
    }
  }
  _szMovedV+=mvT/MM_PER_METER*_pModel->GetHydroUnit(k)->GetArea()*M2_PER_KM2;
  return ovf;
}

//////////////////////////////////////////////////////////////////
/// \brief LAGGED: each column receives its last step's exchange (as much as it can hold or give) before Raven's
///  processes; Raven's soil-zone layers become the columns' mean
//
void CGroundwaterModel::SZLagApply(double **aPhinew,double **aPhi,const double &tstep)
{
  const int n=_szNLayers;
  if (_szxQprev.size()!=_szxLink.size()){return;}
  vector<double> phi(n);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    const int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    bool any=false;
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      const int e=_szxOfHRU[k][jj],c=_szxCol[e]; const double a=_szColA[c];
      if (!(a>0.0) || (_szxQprev[e]==0.0) || _colCHDTop[_aLinks[_szxLink[e]].ic]){continue;} //(over a specified head the exchange is exact:
      for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}                                       // no lag, as for the CHD entries)
      double dV=-_szxQprev[e]*tstep/a*MM_PER_METER;                             //[mm] into the column
      if (dV>0.0){
        double room=0.0; for (int j=0;j<n;j++){room+=max(_szLay[k][j].phimax-phi[j],0.0);}
        dV=min(dV,room); gwsz::Fill(_szLay[k],phi,dV);
      }
      else {
        double avail=0.0; for (int j=0;j<n;j++){avail+=max(phi[j],0.0);}
        dV=max(dV,-avail); double sh=0.0; gwsz::Apply(_szLay[k],phi,dV,sh); dV+=sh;
      }
      for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
      _szxLagL[e]=dV; _szLagNet[k]+=dV/MM_PER_METER*a; any=true;
    }
    if (!any){continue;}
    const int nc=SZNCols(k),o=_szColOff[k]; double sw=0.0; for (int c=0;c<nc;c++){sw+=_szColA[o+c];}
    if (!(sw>0.0)){continue;}
    for (int j=0;j<n;j++){
      double v=0.0; for (int c=0;c<nc;c++){v+=_szColA[o+c]*_szColPhi[(o+c)*n+j];}
      aPhinew[k][_szSV[k][j]]=v/sw;
      if (aPhi!=NULL){aPhi[k][_szSV[k][j]]=v/sw;}                               //(the start of the step for the EULER and HEUN solvers)
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief LAGGED: after the sync of HRU k's columns to Raven's processes, the lagged source is taken back out of each
///  column (the step's exchange is then booked in full, as without the lag); what cannot be taken back (water the
///  processes used) stays and is settled in SZApply (_szxLagD)
//
void CGroundwaterModel::SZLagReverse(const int k)
{
  const int n=_szNLayers;
  if (_szxLagL.size()!=_szxLink.size()){return;}
  vector<double> phi(n);
  for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
    const int e=_szxOfHRU[k][jj],c=_szxCol[e]; const double a=_szColA[c],L=_szxLagL[e];
    if ((L==0.0) || !(a>0.0)){continue;}
    for (int j=0;j<n;j++){phi[j]=_szColPhi[c*n+j];}
    double r=0.0,D=0.0;
    if (L>0.0){double sh=0.0; gwsz::Apply(_szLay[k],phi,-L,sh); r=-L+sh; D=-sh;}
    else      {double ov=gwsz::Fill(_szLay[k],phi,-L); r=-L-ov; D=ov;}
    for (int j=0;j<n;j++){_szColPhi[c*n+j]=phi[j];}
    _szxLagD[e]=D; _szLagNet[k]+=r/MM_PER_METER*a;
  }
}

//////////////////////////////////////////////////////////////////
/// \brief LAGGED: the aquifer repays each link's debt (lagged source Raven's processes used beyond the exchange) over
///  this step, as a fixed rate through SZQ_RAVEN, unless the cell is near its bottom (MODFLOW would reduce the rate)
//
void CGroundwaterModel::SZRepayRates(const double &tstep)
{
  const size_t ne=_szxLink.size();
  _szxRepay.assign(ne,0.0); _szRepayQ=0.0;
  if ((_szTiming!=1) || (_szxDebt.size()!=ne)){return;}
  for (size_t e=0;e<ne;e++){
    if (!(_szxDebt[e]>0.0)){continue;}
    const int ic=_aLinks[_szxLink[e]].ic;
    if (_colCHDTop[ic]){continue;}                                               //(a specified head this step: MODFLOW would ignore the rate)
    if ((e>=_szxUp.size()) || (_szxUp[e]<0) || (_szxUpF[e]<1.0) || _szxUpCHD[e]){continue;}                         //(no cell holds the water table clear of its bottom: later;
    double used=(e<_szxCapQ.size())?_szxCapQ[e]:0.0;                           //(what the source cell has given this step: the rise, and the
    if (_szxUp[e]!=TopNode(ic)){used+=_vzUp[e]+_vzUpA[e];}                      // column's upflow and uptake, through SZW_RAVEN or, out of the
    else if ((_vzMode>0) && (e<_vzUse.size())){used+=max(-_vzAbsQ[e],0.0)+(_vzUse[e]?max(-_vzQ[e],0.0):0.0);} // top cell, its recharge)
    const double left=SZUpCap(e,tstep)-used;                                     //(what it can still give)
    _szxRepay[e]=max(min(_szxDebt[e]/tstep,left),0.0); _szRepayQ+=_szxRepay[e];                    // from a lower cell when the top one has dried, SZW_RAVEN)
  }
}

//////////////////////////////////////////////////////////////////
/// \brief exchange of column c (HRU k) with a specified head h [m] over the step, its rejected recharge rejV [m3]
///  returned to it first: the local solve with the fixed head, exact (the head does not respond) [m3/d into the head]
//
double CGroundwaterModel::SZChdExchange(const int k,const int c,const gwsz::sz_cell &cell,const double rejV,const double h,const double &tstep,const double D) const
{
  const int n=_szNLayers; const double a=_szColA[c];
  if (!(a>0.0)){return 0.0;}
  vector<double> phi(_szColPhi.begin()+c*n,_szColPhi.begin()+(c+1)*n);
  gwsz::Fill(_szLay[k],phi,rejV/a*MM_PER_METER);
  if (D>0.0){gwsz::Fill(_szLay[k],phi,D);} else if (D<0.0){double sh=0.0; gwsz::Apply(_szLay[k],phi,D,sh);} //(a lagged source not taken back, mm)
  double dV=0.0,ovf=0.0,lamp=0.0;
  gwsz::LocalSolve(_szLay[k],phi,cell,h,h,a,tstep,dV,ovf,lamp);
  return -(dV+ovf)/MM_PER_METER*a/tstep;                                       //(what overflows came from the head)
}

//////////////////////////////////////////////////////////////////
/// \brief columns over specified-head cells (soil_zone_columns_plan.md section 5): the local solve with the fixed head is the
///  exact exchange of the step (no MODFLOW boundary: the specified head takes or gives the water, as it takes the
///  recharge onto it); booked as soil zone to / from specified heads, with their capillary rise (SZCapRiseWT). Called
///  after MODFLOW's solve (the specified head is in MODFLOW's heads only from its first iteration on)
//
void CGroundwaterModel::SZCHDSolves(const double &tstep)
{
  const int n=_szNLayers;
  _QszCHD_in=_QszCHD_out=0.0;
  for (size_t e=0;e<_szcLink.size();e++){
    const gw_link &L=_aLinks[_szcLink[e]];
    const double q=SZChdExchange(L.k,_szcCol[e],_szcCell[e],_szcRej[e],_X[_reduced[TopNode(L.ic)]]*_hf,tstep);
    _szcQ[e]=q;
    if (q>0.0){_QszCHD_in+=q;} else {_QszCHD_out-=q;}
  }
  for (size_t e=0;e<_szcLink.size();e++){for (int j=0;j<n;j++){_QszCHD_out+=_szcCapL[e*n+j];}}
}
