/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadoseModel.cpp
  :SoilZone PHYSICAL with :SoilZoneVadose RICHARDS: the unsaturated column of every exchange column between the
  soil-zone base and the water table (soil_zone_columns_plan.md, sections 1 and 2): materials, initial state, one solve per
  Raven step before MODFLOW's (with the heads of the last step), capillary rise into the soil zone, recharge to MODFLOW.
  ----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include "Model.h"

static const double VZ_LMIN=0.05;   //[m] shortest unsaturated column (one node of the finest grid)

//////////////////////////////////////////////////////////////////
/// \brief materials of every exchange column, by MODFLOW layer (after SZBuildEntries)
//
void CGroundwaterModel::VZSetup()
{
  const size_t ne=_szxLink.size();
  _vzCol.assign(ne,gwvz::vz_column()); _vzStrat.assign(ne,gwvz::vz_strat());
  _vzUse.assign(ne,0); _vzQ.assign(ne,0.0); _vzAbsQ.assign(ne,0.0); _vzRejV.assign(ne,0.0); _vzFp.assign(ne,0.0);
  _vzInV.assign(ne,0.0); _vzCapV.assign(ne,0.0);
  gwvz::SetAccuracy(_vzAcc);
  vector<char> warned(_aClasses.size(),0); string noUnsat="";
  for (size_t e=0;e<ne;e++){
    const gw_link &L=_aLinks[_szxLink[e]]; const int ic=L.ic;
    const gw_profile &P=_aProfiles[_colProfile[ic]];
    gwvz::vz_strat &S=_vzStrat[e];
    for (int l=0;l<_nlay;l++){
      int n=UserNode(l,ic);
      if ((n>=(int)_idomain.size()) || (_idomain[n]<=0)){continue;}            //(an inactive layer: the one below continues)
      gwvz::vz_mat m;
      if (IsExternal()){ //an existing model: conductivity and specific yield of its cells, retention of :SoilZoneAquiferClass
        const int ix=GetAquiferClassIndex(_szExtClass);
        m.Ks=_extKv[n]; m.Sy=min(max(_extSy[n],0.01),0.6); m.Ss=1e-5;
        if (ix!=DOESNT_EXIST){const gw_aquifer_class &C=_aClasses[ix]; m.b=C.clapp_b; m.psia=C.psia; m.ths=max(C.porosity,m.Sy+0.01);}
        else {m.b=4.0; m.psia=0.3; m.ths=m.Sy+0.1;}
        if (!(m.Ks>0.0)){continue;}
      }
      else {
        if (l>=(int)P.layers.size()){break;}
        const gw_aquifer_class &C=_aClasses[P.layers[l].iclass];
        m.ths=C.porosity; m.b=C.clapp_b; m.psia=C.psia; m.Ks=C.Kv; m.Sy=C.Sy; m.Ss=max(C.Ss,1e-7);
        if (!C.unsatGiven && !warned[P.layers[l].iclass]){warned[P.layers[l].iclass]=1; noUnsat+=(noUnsat==""?"":", ")+C.name;}
      }
      double zt=(S.mat.empty())?1e9:S.zbot.back();                             //(the first material extends up without limit)
      S.ztop.push_back(zt); S.zbot.push_back(_botm[n]); S.mat.push_back(m);
    }
    if (!S.zbot.empty()){S.zbot.back()=-1e9;}                                  //(and the last one down)
  }
  if (noUnsat!=""){
    WriteWarning("CGroundwaterModel: :SoilZoneVadose RICHARDS: aquifer classes "+noUnsat+" have no CLAPP_B and AIR_ENTRY_PRESSURE; the unsaturated column "
                 "uses 4 and 0.3 m (give both in :AquiferClasses)",false);
  }
}

//////////////////////////////////////////////////////////////////
/// \brief link li's share of soil-zone HRU k's recharge, by the rules of SZRecharge: the columns' losses to recharge
///  (_szxFw, _szcFw; W their sum over the HRU), else by area over all the HRU's links
//
double CGroundwaterModel::VZRchW(const int k,double &Aall) const
{
  double W=0.0; Aall=0.0;
  for (size_t j=0;j<_hruLinks[k].size();j++){int l2=_hruLinks[k][j]; Aall+=_aLinks[l2].area; if (_linkSzx[l2]>=0){W+=_szxFw[_linkSzx[l2]];}}
  for (size_t j=0;j<_szcOfHRU[k].size();j++){W+=_szcFw[_szcOfHRU[k][j]];}
  return W;
}
double CGroundwaterModel::VZLinkShare(const int k,const int li) const
{
  double Aall; const double W=VZRchW(k,Aall);
  if (W>0.0){return (_linkSzx[li]>=0)?_szxFw[_linkSzx[li]]/W:((_linkSzc[li]>=0)?_szcFw[_linkSzc[li]]/W:0.0);}
  return (Aall>0.0)?_aLinks[li].area/Aall:0.0;
}

//////////////////////////////////////////////////////////////////
/// \brief water above the reference contents (ths - Sy, what a falling water table leaves) in all columns [m3]
//
double CGroundwaterModel::VZStorage() const
{
  double v=0.0;
  for (size_t e=0;e<_vzCol.size();e++){if (_vzCol[e].active){v+=gwvz::RelStorage(_vzCol[e],_vzStrat[e])*_szColA[_szxCol[e]];}}
  return v;
}

//////////////////////////////////////////////////////////////////
/// \brief columns at the start of the run: from the hotstart file (:GWVadoseColumns) where it has them, otherwise at the
///  hydrostatic profile above the initial water table (no flow until recharge arrives)
//
void CGroundwaterModel::VZInit()
{
  int nHot=0,nBad=0;
  for (size_t e=0;e<_vzCol.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]]; const int ic=L.ic;
    const double h=_X[_reduced[TopNode(ic)]]*_hf,zb=_szxZb[e],Lmin=VZ_LMIN;
    gwvz::vz_column &c=_vzCol[e]; const gwvz::vz_strat &S=_vzStrat[e];
    c=gwvz::vz_column(); c.ztop=zb; c.zbot=h;
    map<long long,map<int,vector<double> > >::const_iterator it=_hotVz.find(_pModel->GetHydroUnit(L.k)->GetHRUID());
    if (it!=_hotVz.end()){
      map<int,vector<double> >::const_iterator jt=it->second.find(ic+1);
      if ((jt!=it->second.end()) && (jt->second.size()>=2)){
        const vector<double> &v=jt->second;                                    //zbot, then the contents top down
        gwvz::Build(c,S,zb,v[0],_vzN);
        if (c.active && (c.th.size()==v.size()-1)){
          for (size_t i=0;i<c.th.size();i++){c.psi[i]=gwvz::PsiOf(S.mat[c.im[i]],v[i+1]); c.th[i]=gwvz::Theta(S.mat[c.im[i]],c.psi[i]);} //(content and suction consistent)
          nHot++;
          continue;
        }
        if (zb-h>Lmin){nBad++;}                                                 //(another node count: the column starts at equilibrium;
        c=gwvz::vz_column(); c.ztop=zb; c.zbot=h;                                //  a water table now at the base: no column)
      }
    }
    if (zb-h>Lmin){gwvz::Build(c,S,zb,h,_vzN);}
  }
  _vzStore0=_vzStoreVol=VZStorage();
  if (nHot>0){WriteAdvisory("CGroundwaterModel: "+to_string(nHot)+" unsaturated columns read from the hotstart file",false);}
  if (nBad>0){WriteWarning("CGroundwaterModel: "+to_string(nBad)+" unsaturated columns of the hotstart file do not match the columns (another :SoilZoneVadose node count?); they start at equilibrium",false);}
}

//////////////////////////////////////////////////////////////////
/// \brief the unsaturated columns over this step, with the water table at the start of the step (soil_zone_columns_plan.md
///  section 1.3): the column follows the water table (water below a risen one goes to the aquifer this step), takes the
///  unsaturated share of the column's recharge at its top (less the capillary rise into the soil zone), and gives its
///  bottom flux to MODFLOW's recharge (SZRecharge). Solved once per step, before MODFLOW's solve.
//
void CGroundwaterModel::VZStep(const double &tstep)
{
  const int n=_szNLayers;
  const bool cap=(_szCapRise==1) || ((_szCapRise<0) && (_vzMode>0));
  _vzInVol=_vzOutVol=_vzCapVol=0.0;
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const soil_struct *sp=_pModel->GetHydroUnit(k)->GetSoilProps(_szNLayers-1);   //the bottom soil-zone layer's soil
    const gwsz::sz_layer &Lb=_szLay[k][0];
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      const int e=_szxOfHRU[k][jj],li=_szxLink[e],c=_szxCol[e]; const gw_link &L=_aLinks[li];
      const double a=_szColA[c];
      _vzUse[e]=0; _vzQ[e]=0.0; _vzAbsQ[e]=0.0; _vzRejV[e]=0.0; _vzInV[e]=0.0; _vzCapV[e]=0.0;
      if (!(a>0.0) || _colCHDTop[L.ic]){continue;}                              //(a specified head this step: its recharge goes there)
      gwvz::vz_column &col=_vzCol[e]; const gwvz::vz_strat &S=_vzStrat[e];
      const double h=_X[_reduced[TopNode(L.ic)]]*_hf,zb=_szxZb[e],Lmin=VZ_LMIN;
      //follow the water table: a rising one takes the excess over equilibrium below it, a falling one opens cells at
      //equilibrium; a column shorter than Lmin (5 cm) is the soil zone's contact exchange and capillary rise instead
      if (col.active){
        double ab=(zb-h<=Lmin)?gwvz::MoveBottom(col,S,zb,_vzN):gwvz::MoveBottom(col,S,h,_vzN);
        _vzAbsQ[e]=ab*a/tstep;
      }
      else if (zb-h>Lmin){
        gwvz::Build(col,S,zb,h,_vzN);                                            //(at equilibrium: no excess, nothing exchanged)
      }
      if (!col.active){continue;}
      //inflow: the unsaturated share of the column's recharge (the saturated share is rejected as before, with this
      //saturated fraction for the whole step)
      const double share=VZLinkShare(k,li),F=_linkF[li];
      const double P=_szRchV[k]*share*(1.0-F);                                  //[m3] over the step
      _vzFp[e]=F; _vzInV[e]=P;
      //capillary rise into the bottom soil-zone layer: the steady upward flux from the column's top node to the layer's
      //middle at its Campbell suction (Raven's AIR_ENTRY_PRESSURE, CLAPP_B, HYDRAUL_COND of that soil)
      double R=0.0;
      if (cap && (sp!=NULL)){
        const double phi=_szColPhi[c*n+0];
        const double Ssat=min(max(phi/max(Lb.phimax,1e-12),0.01),1.0),bs=max(sp->clapp_b,0.1),psia=max(sp->air_entry_pressure,0.0)/MM_PER_METER;
        const double psis=psia*pow(Ssat,-bs);                                    //[m]
        if (phi<Lb.phifc){R=gwvz::CapRise(col,max(sp->hydraul_cond,0.0)/MM_PER_METER,psia,bs,psis,0.5*Lb.H);}
        if (R>0.0){
          const gwvz::vz_mat &m0=S.mat[col.im[0]];
          double room=max(Lb.phifc-phi,0.0)/MM_PER_METER;                       //[m] the layer's deficit to field capacity
          double avail=max(col.th[0]-gwvz::Theta(m0,-min(psis,1e5)),0.0)*col.dz[0]; //[m] the top node down to the soil's suction
          R=min(R,min(room,avail)/tstep);
          _vzCapV[e]=R*tstep*a; _vzCapVol+=_vzCapV[e];                           //(added to the layer in SZApply, after the columns' sync)
        }
      }
      double rej=0.0; gwvz::vz_step_info inf;
      const double bot=_botm[TopNode(L.ic)],ctop=_top[L.ic];
      const bool noUp=(h<bot+max(0.05*(ctop-bot),0.05));                         //(a top cell MODFLOW has (nearly) dried: no upward flow
      double qb=gwvz::Step(col,S,tstep,P/(a*tstep)-R,rej,&inf,noUp);             // out of it, which MODFLOW could not supply)
      if (noUp && (qb<0.0)){ //(the Newton tolerance leaves a trace of upward flow: it is taken back out of the bottom node,
        const size_t ib=col.th.size()-1;                                         // not drawn from a dry cell)
        col.th[ib]=max(col.th[ib]+qb*tstep/col.dz[ib],0.0); col.psi[ib]=gwvz::PsiOf(S.mat[col.im[ib]],col.th[ib]); qb=0.0;
      }
      _vzSub+=inf.substeps; _vzNewton+=inf.newton;
      if (rej<0.0){_vzCapV[e]=max(_vzCapV[e]+rej*a,0.0); _vzCapVol+=rej*a; rej=0.0;} //(an upward demand the column could not supply)
      _vzRejV[e]=rej*a; _vzQ[e]=qb*a; _vzUse[e]=1;
      _vzInVol+=P-_vzRejV[e]; _vzOutVol+=(_vzQ[e]+_vzAbsQ[e])*tstep;           //(inflow the column took; its rejection is in the rejected recharge)
    }
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){int e=_szxOfHRU[k][jj]; if (!_vzUse[e]){_vzOutVol+=_vzAbsQ[e]*tstep;}}
  }
  _vzStoreVol=VZStorage();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief capillary rise from the water table into the soil zone where there is no unsaturated column (no
///  :SoilZoneVadose, or the water table within the smoothing depth of the soil-zone base or above it), soil_zone_columns_plan.md
///  section 2. Into every soil-zone layer below field capacity whose middle is above the water table: the exact steady
///  upward flux (gwvz::Rise) from the water table d below the layer's middle to the layer's suction psi_s (its CLAPP_B,
///  AIR_ENTRY_PRESSURE, HYDRAUL_COND; none when psi_s <= d, hydrostatic or wetter), on the unsaturated share of the
///  column, at most the layer's deficit to field capacity; all layers together at most the largest steady rise to the
///  lowest of them (psi_s infinite: the water reaching the upper layers passes through it). With
///  the heads of the last step, a fixed rate for the step: taken from the aquifer through SZQ_RAVEN (SZSetArrays) and
///  added to the layers in SZApply.
//
void CGroundwaterModel::SZCapRiseWT(const double &tstep)
{
  const int n=_szNLayers; const size_t ne=_szxLink.size();
  _szxCapQ.assign(ne,0.0); _szxCapL.assign(ne*n,0.0); _szCapWT=0.0;
  _szcCapL.assign(_szcLink.size()*n,0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
    const int nx=(int)_szxOfHRU[k].size(),nch=(int)_szcOfHRU[k].size();
    for (int jj=0;jj<nx+nch;jj++){                                              //(exchange columns, then those over specified heads,
      const bool sc=(jj>=nx);                                                    // whose rise the specified head supplies)
      const int e=sc?_szcOfHRU[k][jj-nx]:_szxOfHRU[k][jj],li=sc?_szcLink[e]:_szxLink[e],c=sc?_szcCol[e]:_szxCol[e]; const gw_link &L=_aLinks[li];
      const bool ch=sc || _colCHDTop[L.ic];                                      //(also an exchange column on a cell that is a specified head this step)
      const double a=_szColA[c];
      if (!(a>0.0)){continue;}
      if (!ch && (_vzMode>0) && (e<(int)_vzCol.size()) && _vzCol[e].active){continue;}  //(the column's own rise, VZStep)
      const double h=_X[_reduced[TopNode(L.ic)]]*_hf;
      const double bot=_botm[TopNode(L.ic)],top=_top[L.ic];
      if (!ch && (h<bot+max(0.05*(top-bot),0.05))){continue;}                   //(MODFLOW reduces a rate taken near the cell bottom)
      vector<double> &CL=sc?_szcCapL:_szxCapL;
      const double F=min(max(_linkF[li],0.0),1.0); if (F>=1.0){continue;}
      double z=sc?SZBase(L.k,L.ic):_szxZb[e],Qmax=-1.0,Qs=0.0;
      for (int j=0;j<n;j++){
        const gwsz::sz_layer &Lj=_szLay[k][j]; const double phi=_szColPhi[c*n+j];
        const double zm=z+0.5*Lj.H; z+=Lj.H;
        if ((zm<=h) || (phi>=Lj.phifc)){continue;}                               //(below the water table, or wet enough)
        const soil_struct *sp=pHRU->GetSoilProps(n-1-j);
        if (sp==NULL){continue;}
        const double bs=max(sp->clapp_b,0.1),psia=max(sp->air_entry_pressure,0.0)/MM_PER_METER,Ks=max(sp->hydraul_cond,0.0)/MM_PER_METER;
        const double d=zm-h,Ssat=min(max(phi/max(Lj.phimax,1e-12),0.01),1.0);
        const double psis=psia*pow(Ssat,-bs);                                    //[m] the layer's suction
        if (Qmax<0.0){Qmax=gwvz::Rise(Ks,psia,bs,d,-1.0)*a*(1.0-F);}             //[m3/d] the largest rise to the lowest layer that takes water
        if (!(psis>d)){continue;}                                                //(at or wetter than hydrostatic: no rise)
        double Q=min(gwvz::Rise(Ks,psia,bs,d,psis)*a*(1.0-F),(Lj.phifc-phi)/MM_PER_METER*a/tstep);
        if (Q>0.0){CL[e*n+j]=Q; Qs+=Q;}
      }
      if (!(Qs>0.0)){continue;}
      if (Qs>Qmax){const double r=Qmax/Qs; for (int j=0;j<n;j++){CL[e*n+j]*=r;} Qs=Qmax;}
      if (!ch){_szxCapQ[e]=Qs; _szCapWT+=Qs;}
    }
  }
}
