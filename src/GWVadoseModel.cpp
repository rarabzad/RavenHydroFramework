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
  _vzInV.assign(ne,0.0); _vzCapV.assign(ne,0.0); _vzCapL.assign(ne*_szNLayers,0.0);
  gwvz::SetAccuracy(_vzAcc);
  vector<char> warned(_aClasses.size(),0),warnSy(_aClasses.size(),0); string noUnsat="",badSy=""; bool extBadSy=false;
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
      { //Sy against the drainable porosity of the retention curve 1 m above the water table (soil_zone_columns_plan.md 1.5): a
        //rising water table releases the column's water above equilibrium to MODFLOW, which stores it with Sy; with Sy
        //much smaller each step's rise releases more than the last (a feedback that stops only at the soil-zone base)
        const double dp=m.ths-gwvz::Theta(m,-1.0);
        if ((m.Sy>2.0*dp) || (m.Sy<0.5*dp)){
          if (IsExternal()){extBadSy=true;}
          else {const int ic2=P.layers[l].iclass; if (!warnSy[ic2]){warnSy[ic2]=1; badSy+=(badSy==""?"":", ")+_aClasses[ic2].name+" (Sy "+to_string(m.Sy)+", drainable porosity "+to_string(dp)+")";}}
        }
      }
      double zt=(S.mat.empty())?1e9:S.zbot.back();                             //(the first material extends up without limit)
      S.ztop.push_back(zt); S.zbot.push_back(_botm[n]); S.mat.push_back(m);
    }
    if (!S.zbot.empty()){S.zbot.back()=-1e9;}                                  //(and the last one down)
  }
  if ((badSy!="") || extBadSy){
    WriteWarning("CGroundwaterModel: :SoilZoneVadose RICHARDS: the specific yield differs by more than a factor 2 from the drainable porosity "
                 "theta_s - theta(psi = -1 m) of POROSITY, CLAPP_B and AIR_ENTRY_PRESSURE"+(badSy!=""?(" for "+badSy):string(" in cells of the model and :SoilZoneAquiferClass"))+
                 "; a rising water table then releases more (or less) water from the unsaturated column than MODFLOW stores with Sy, and with Sy too small "
                 "the water table rises step by step to the soil-zone base",false);
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
    double h=_X[_reduced[TopNode(ic)]]*_hf; const double zb=_szxZb[e],Lmin=VZ_LMIN;
    if (!(std::fabs(h)<1e29)){h=WaterTable(ic);}                                 //(a dry top cell without NEWTON: the water table below it)
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
          if (!c.th.empty()){c.dryTop=(c.psi[0]<=-0.999999*1e5); c.satTop=(c.psi[0]>=-1e-9);} //(the top's memory: a top held at air-dry or
                                                                                 // saturation by the last step, as in a continuous run)
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
  _vzInVol=_vzOutVol=_vzCapVol=0.0; _vzUpQ=0.0;
  vector<double> QL(max(n,1),0.0),QW(max(n,1),0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const soil_struct *sp=_pModel->GetHydroUnit(k)->GetSoilProps(_szNLayers-1);   //the bottom soil-zone layer's soil
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      const int e=_szxOfHRU[k][jj],li=_szxLink[e],c=_szxCol[e]; const gw_link &L=_aLinks[li];
      const double a=_szColA[c];
      _vzUse[e]=0; _vzQ[e]=0.0; _vzAbsQ[e]=0.0; _vzRejV[e]=0.0; _vzInV[e]=0.0; _vzCapV[e]=0.0; _vzUp[e]=0.0; _vzUpA[e]=0.0;
      for (int j=0;j<n;j++){_vzCapL[e*n+j]=0.0;}
      if (!(a>0.0) || _colCHDTop[L.ic]){continue;}                              //(a specified head this step: its recharge goes there)
      gwvz::vz_column &col=_vzCol[e]; const gwvz::vz_strat &S=_vzStrat[e];
      const double h=_X[_reduced[TopNode(L.ic)]]*_hf,zb=_szxZb[e],Lmin=VZ_LMIN;
      if (!(std::fabs(h)<1e29)){continue;}                                       //(a top cell out of MODFLOW's solve, HDRY without NEWTON: the
                                                                                 // column is left as it is until the cell rewets)
      //follow the water table: a rising one takes the excess over equilibrium below it, a falling one opens cells at
      //equilibrium; a column shorter than Lmin (5 cm) is the soil zone's contact exchange and capillary rise instead
      if (col.active){
        //a rising water table takes the column's water above equilibrium below it, at most what MODFLOW stores for that
        //rise (Sy dh): with Sy smaller than the drainable porosity each step's rise would release more than the last (a
        //feedback with gain eps/Sy > 1, design note 1.5); the column's bottom then follows more slowly, releasing the rest
        //in the next steps (with consistent parameters eps <= Sy and this does not apply)
        const double zt=(zb-h<=Lmin)?zb:h;
        if (zt>col.zbot+1e-9){ //a rising water table (to the base: the column goes)
          const double z0=col.zbot,Sy=S.mat[S.Find(z0+1e-6)].Sy,capV=Sy*(zt-z0);   //[m]
          gwvz::vz_column t=col; const double abT=gwvz::MoveBottom(t,S,zt,_vzN);
          const double capU=_szxUpCHD[e]?1e300:SZUpCap(e,tstep)*tstep/a;         //[m] a deficit below it (a column dried by evaporation or
          if ((zt<zb) && (abT>capV*(1.0+1e-9))){                                 // capillary rise) is filled out of the source cell, at most what
            double lo=z0,hi=min(zt,col.ztop);                                    // it can give this step (SZUpCap)
            for (int it=0;it<60;it++){ //(the excess below a level grows with the level)
              const double zm=0.5*(lo+hi); t=col;
              if (gwvz::MoveBottom(t,S,zm,_vzN)>capV){hi=zm;} else {lo=zm;}
              if (hi-lo<1e-7){break;}
            }
            const double zn=(lo>col.ztop-0.01)?col.ztop:lo;                      //(a remainder thinner than 1 cm: the column goes)
            _vzAbsQ[e]=gwvz::MoveBottom(col,S,zn,_vzN)*a/tstep;
            if (_vzSlowed++==0){WriteWarning("CGroundwaterModel: :SoilZoneVadose RICHARDS: a rising water table met more water above equilibrium in an unsaturated "
                                             "column than MODFLOW stores with Sy; the column releases at most Sy dh per step (Sy below the drainable porosity?)",false);}
          }
          else if ((abT<0.0) && (-abT>capU)){                                    //(the bottom then follows more slowly; the column stays)
            double lo=z0,hi=min(zt,col.ztop-1e-3);
            for (int it=0;it<60;it++){
              const double zm=0.5*(lo+hi); t=col;
              if (-gwvz::MoveBottom(t,S,zm,_vzN)>capU){hi=zm;} else {lo=zm;}
              if (hi-lo<1e-7){break;}
            }
            _vzAbsQ[e]=(lo>z0+1e-9)?gwvz::MoveBottom(col,S,lo,_vzN)*a/tstep:0.0;
          }
          else {col=t; _vzAbsQ[e]=abT*a/tstep;}
          if ((_vzAbsQ[e]<0.0) && (_szxUp[e]>=0) && (_szxUp[e]!=TopNode(L.ic)) && !_szxUpCHD[e]){_vzUpA[e]=-_vzAbsQ[e]; _vzUpQ+=_vzUpA[e];}
        }
        else if (zt<col.zbot-1e-9){ //a falling water table: the new nodes start at equilibrium, the water they take comes out of the
          //cell holding the water table, at most what it can give this step (SZUpCap: none without such a cell, the bottom then
          //stays); the column's bottom follows more slowly otherwise
          const double capV=(_szxUpCHD[e]?1e300:SZUpCap(e,tstep)*tstep/a);        //[m]
          gwvz::vz_column t=col; const double abT=gwvz::MoveBottom(t,S,zt,_vzN);
          if (-abT>capV){
            double lo=zt,hi=col.zbot;                                            //(lo: too far down; hi: not moved)
            for (int it=0;it<60;it++){
              const double zm=0.5*(lo+hi); t=col;
              if (-gwvz::MoveBottom(t,S,zm,_vzN)>capV){lo=zm;} else {hi=zm;}
              if (hi-lo<1e-7){break;}
            }
            _vzAbsQ[e]=(hi<col.zbot-1e-9)?gwvz::MoveBottom(col,S,hi,_vzN)*a/tstep:0.0;
          }
          else {col=t; _vzAbsQ[e]=abT*a/tstep;}
          if ((_vzAbsQ[e]<0.0) && (_szxUp[e]>=0) && (_szxUp[e]!=TopNode(L.ic)) && !_szxUpCHD[e]){_vzUpA[e]=-_vzAbsQ[e]; _vzUpQ+=_vzUpA[e];}
        }
        //(unmoved: nothing)
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
      //capillary rise into the soil zone (SZRiseLayers, as from the water table where there is no column): the steady
      //upward flux from the column's top node (its suction) into every layer below field capacity, on the unsaturated
      //share. Per step at most the larger of the top node's water down to the suction of the lowest receiving layer and
      //what the water table itself could supply with its drawdown: as the column shortens the top node empties and the
      //rise becomes the water table's (the two meet at the shortest column, no jump), and what the column cannot pass
      //on within the step comes back as an unmet upward demand below
      double R=0.0;
      if (cap && (sp!=NULL)){
        double psL=0.0;
        const double Ae=a*(1.0-min(max(F,0.0),1.0));
        const double Qc=SZRiseLayers(k,c,zb,col.zc[0],max(-col.psi[0],0.0),a,Ae,0.0,tstep,&QL[0],&psL);
        if (Qc>0.0){
          const gwvz::vz_mat &m0=S.mat[col.im[0]];
          const double avail=max(col.th[0]-gwvz::Theta(m0,-min(psL,1e5)),0.0)*col.dz[0]*a; //[m3] the top node down to the lowest layer's suction
          const int up=_szxUp[e];                                                //(as SZCapRiseWT: from the cell holding the water table; only
          double Qw=0.0;                                                         // needed when the top node cannot give the rise)
          if ((up>=0) && (Qc>avail/tstep)){
            const double hu=_X[_reduced[up]]*_hf,Syu=_szxUpCHD[e]?0.0:max(SZNodeSy(up),0.01);
            Qw=SZRiseLayers(k,c,zb,hu,0.0,a,Ae,Syu,tstep,&QW[0],NULL,VZGapMat(e));
            if (!_szxUpCHD[e]){Qw=min(Qw*_szxUpF[e],SZUpCap(e,tstep));}            //(as SZCapRiseWT)
          }
          const double Q=min(Qc,max(avail/tstep,Qw));                            //[m3/d]
          if (Q>0.0){
            R=Q/a;                                                               //[m/d] out of the column's top
            _vzCapV[e]=Q*tstep; _vzCapVol+=_vzCapV[e];                           //(added to the layers in SZApply, after the columns' sync)
            for (int j=0;j<n;j++){_vzCapL[e*n+j]=QL[j]/Qc;}
          }
        }
      }
      double rej=0.0; gwvz::vz_step_info inf;
      const bool noUp=(_szxUp[e]<0);                                            //(no cell holds any water table: no upward
      double qb=gwvz::Step(col,S,tstep,P/(a*tstep)-R,rej,&inf,noUp);             // flow, which MODFLOW could not supply)
      if (qb<0.0){ //upflow through the bottom: out of the cell holding the water table, scaled down as it nears drying and at
        //most half the water it holds above its bottom under the column per step (none without such a cell); the part
        //MODFLOW does not supply is taken back out of the column, bottom node first, down to air-dry
        const int up=_szxUp[e];
        double allow=0.0;                                                        //[m/d]
        if (!noUp && _szxUpCHD[e]){allow=-qb;}                                   //(a specified head supplies all of it)
        else if (!noUp){                                                         //(the source cell's limit, less what the new nodes took)
          allow=min(_szxUpF[e]*(-qb),max(SZUpCap(e,tstep)+min(_vzAbsQ[e],0.0),0.0)/a);
        }
        double need=(-qb-allow)*tstep;                                           //[m]
        for (size_t ib=col.th.size();(ib-->0) && (need>0.0);){
          const gwvz::vz_mat &mb=S.mat[col.im[ib]];
          const double take=min(need,max(col.th[ib]-gwvz::Theta(mb,-1e5),0.0)*col.dz[ib]);
          if (take>0.0){col.th[ib]-=take/col.dz[ib]; col.psi[ib]=gwvz::PsiOf(mb,col.th[ib]); need-=take; qb+=take/tstep;}
        }
      }
      _vzSub+=inf.substeps; _vzNewton+=inf.newton;
      if (rej<0.0){_vzCapV[e]=max(_vzCapV[e]+rej*a,0.0); _vzCapVol+=rej*a; rej=0.0;} //(an upward demand the column could not supply)
      _vzRejV[e]=rej*a; _vzQ[e]=qb*a; _vzUse[e]=1;
      if ((qb<0.0) && (_szxUp[e]>=0) && (_szxUp[e]!=TopNode(L.ic)) && !_szxUpCHD[e]){_vzUp[e]=-qb*a; _vzUpQ+=_vzUp[e];} //(upflow from a lower
                                                                                 // cell: SZW_RAVEN; otherwise, also a remainder without a source cell, the top cell's recharge)
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
///  lowest of them (psi_s infinite: the water reaching the upper layers passes through it). Implicit in the drawdown it
///  causes (SZRiseLayers): the rate at the end-of-step water table, lowered by the risen volume over Sy and the share's
///  area. With the heads of the last step otherwise, a fixed rate for the step: taken from the aquifer through
///  SZQ_RAVEN (SZSetArrays) and added to the layers in SZApply.
//
void CGroundwaterModel::SZCapRiseWT(const double &tstep)
{
  const int n=_szNLayers; const size_t ne=_szxLink.size();
  _szxCapQ.assign(ne,0.0); _szxCapL.assign(ne*n,0.0); _szCapWT=0.0;
  _szcCapL.assign(_szcLink.size()*n,0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const int nx=(int)_szxOfHRU[k].size(),nch=(int)_szcOfHRU[k].size();
    for (int jj=0;jj<nx+nch;jj++){                                              //(exchange columns, then those over specified heads,
      const bool sc=(jj>=nx);                                                    // whose rise the specified head supplies)
      const int e=sc?_szcOfHRU[k][jj-nx]:_szxOfHRU[k][jj],li=sc?_szcLink[e]:_szxLink[e],c=sc?_szcCol[e]:_szxCol[e]; const gw_link &L=_aLinks[li];
      const bool ch=sc || _colCHDTop[L.ic];                                      //(also an exchange column on a cell that is a specified head this step)
      const double a=_szColA[c];
      if (!(a>0.0)){continue;}
      if (!ch && (_vzMode>0) && (e<(int)_vzCol.size()) && _vzCol[e].active){continue;}  //(the column's own rise, VZStep)
      const int up=ch?TopNode(L.ic):_szxUp[e];                                  //(the cell holding the water table: the top cell, or one below
      if (up<0){continue;}                                                       // a dried top cell through SZW_RAVEN; none: no rise)
      const bool chu=ch || _szxUpCHD[e];                                         //(a specified head supplies it, top cell or below)
      const double h=_X[_reduced[up]]*_hf;
      vector<double> &CL=sc?_szcCapL:_szxCapL;
      const double F=min(max(_linkF[li],0.0),1.0); if (F>=1.0){continue;}
      const double z0=sc?SZBase(L.k,L.ic):_szxZb[e];
      const double Sy=chu?0.0:max(SZNodeSy(up),0.01);                           //(a specified head does not draw down)
      double Qs=SZRiseLayers(k,c,z0,h,0.0,a,a*(1.0-F),Sy,tstep,&CL[e*n],NULL,sc?NULL:VZGapMat(e));
      if (!chu){ //out of an aquifer cell: fades out as it nears drying, and at most half the water it holds above its bottom
        const double cap=SZUpCap(e,tstep);                                       // under the link per step (a thin cell would be pumped dry)
        const double sc2=(Qs*_szxUpF[e]>cap)?cap/Qs:_szxUpF[e];
        if (sc2<1.0){Qs*=sc2; for (int j=0;j<n;j++){CL[e*n+j]*=sc2;}}
      }
      if (!(Qs>0.0)){continue;}
      if (!chu){_szxCapQ[e]=Qs; _szCapWT+=Qs;}
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief capillary rise into the soil-zone layers of column c (HRU k, soil-zone base zb) from a source at elevation zs
///  with suction s1 (the water table: s1 = 0; a column's top node: its suction): into every layer below field capacity
///  whose middle is above zs, the exact steady upward flux over the distance from zs to the layer's middle up to the
///  layer's suction (gwvz::RiseBetween; none when hydrostatic or wetter), on the area Ae, at most the layer's deficit to
///  field capacity over the column area a; all layers together at most the rise to the lowest of them with an infinitely
///  dry top. With Sy > 0 the source is the water table and the rate is implicit in its drawdown: Q = T(Q dt/(Sy Ae)),
///  the rate with the water table lowered by the risen volume over the share's area (no lateral resupply: the bound of
///  a local column, stable for every step length; the explicit rate overshoots, e.g. 0.1 m of rise and 0.7 m of
///  drawdown in one day). With gap (the aquifer material just below the base, :SoilZoneVadose) the water from the water
///  table first crosses that material: at most its steady rise to an infinitely dry top over the gap (the column above
///  the shortest column length sees the same material). QL gets each layer's rate [m3/d]; psisLow the suction of the
///  lowest receiving layer [m].
//
double CGroundwaterModel::SZRiseLayers(const int k,const int c,const double zb,const double zs,const double s1,const double a,
                                       const double Ae,const double Sy,const double tstep,double *QL,double *psisLow,const gwvz::vz_mat *gap) const
{
  const int n=_szNLayers; const CHydroUnit *pHRU=_pModel->GetHydroUnit(k);
  //jl: the lowest layer whose middle lies above the source, through which all rising water passes (also when it is at
  //field capacity itself: the cap below and the column's limit use it, so a layer reaching field capacity switches nothing)
  vector<double> d0(n,-1.0),Ks(n,0.0),pa(n,0.0),bs(n,1.0),ps(n,0.0),room(n,0.0); int jl=-1; bool any=false;
  double z=zb;
  for (int j=0;j<n;j++){
    QL[j]=0.0; d0[j]=-1.0;
    const gwsz::sz_layer &Lj=_szLay[k][j]; const double phi=_szColPhi[c*n+j];
    const double zm=z+0.5*Lj.H; z+=Lj.H;
    if (zm<=zs){continue;}                                                       //(below the source)
    const soil_struct *sp=pHRU->GetSoilProps(n-1-j);
    if (sp==NULL){continue;}
    bs[j]=max(sp->clapp_b,0.1); pa[j]=max(sp->air_entry_pressure,0.0)/MM_PER_METER; Ks[j]=max(sp->hydraul_cond,0.0)/MM_PER_METER;
    if (!(Ks[j]>0.0) || !(pa[j]>0.0)){continue;}
    const double Ssat=min(max(phi/max(Lj.phimax,1e-12),0.01),1.0);
    ps[j]=pa[j]*pow(Ssat,-bs[j]);
    if (jl<0){jl=j; d0[j]=-(zm-zs);}                                             //(its distance kept for the cap: -d)
    if (phi>=Lj.phifc){continue;}                                                //(wet enough: takes none)
    d0[j]=zm-zs; room[j]=(Lj.phifc-phi)/MM_PER_METER*a/tstep; any=true;
  }
  if (psisLow!=NULL){*psisLow=(jl>=0)?ps[jl]:0.0;}
  if ((jl<0) || !any || !(Ae>0.0)){return 0.0;}
  const double dl=fabs(d0[jl]);                                                  //(distance to the lowest layer's middle)
  vector<double> buf(n,0.0),bufLo(n,0.0);
  auto T=[&](const double D,double *fill)->double{                               //all layers' rate with the source lowered by D
    double Qs=0.0;
    for (int j=jl;j<n;j++){
      if (d0[j]<0.0){if (fill){fill[j]=0.0;} continue;}                          //(below the source, or at field capacity)
      const double Q=min(gwvz::RiseBetween(Ks[j],pa[j],bs[j],d0[j]+D,s1,ps[j])*Ae,room[j]);
      if (fill){fill[j]=Q;} Qs+=Q;
    }
    if (!(Qs>0.0)){return 0.0;}
    double Qmax=gwvz::RiseBetween(Ks[jl],pa[jl],bs[jl],dl+D,s1,-1.0)*Ae;
    if ((gap!=NULL) && (zb>zs-D)){                                               //(through the aquifer between the water table and the base,
      Qmax=min(Qmax,gwvz::RiseBetween(gap->Ks,gap->psia,gap->b,zb-zs+D,0.0,-1.0)*Ae); // at most its own rise to an infinitely dry top)
    }
    if (Qs>Qmax){if (fill){for (int j=jl;j<n;j++){fill[j]*=Qmax/Qs;}} Qs=Qmax;}
    return Qs;
  };
  auto Put=[&](const double Q,const double Tq,const double *f)->double{        //the layers' shares of T at Q's drawdown, scaled to Q
    for (int j=jl;j<n;j++){QL[j]=(Tq>0.0)?f[j]*Q/Tq:0.0;}
    return Q;
  };
  const double T0=T(0.0,&buf[0]);
  if (!(T0>0.0)){return 0.0;}
  if (!(Sy>0.0)){return Put(T0,T0,&buf[0]);}
  //implicit in the drawdown: g(Q) = Q - T(Q w), w = dt/(Sy Ae), increases; T falls with D, so the root lies between
  //lo = T(T0 w) and T0, with g(T0) = T0 - lo known. Illinois on that bracket (a geometric bisection step after one that
  //does not halve it) to 0.2 % of the root (the rise table itself is within about 2 %): 2 to 5 evaluations of T; none
  //beyond the first when the drawdown changes the rate by less than that. The result is the bracket's low side (never
  //more than the root)
  const double RTOL=2e-3,w=tstep/(Sy*Ae);
  double lo=T(T0*w,&buf[0]),hi=T0;
  if (hi-lo<=RTOL*lo){return Put(lo,lo,&buf[0]);}                              //(the shares at T0's drawdown: within the tolerance)
  double Tlo=T(lo*w,&bufLo[0]);
  double glo=lo-Tlo,ghi=hi-lo;
  //the best lower bound Lb of the root: lo, and T(Q w) for every Q above the root (T falls); its layer shares
  double Lb=lo,LbT=Tlo; vector<double> bufLb(bufLo);
  int side=0; bool bis=false; double wprev=hi-lo;
  for (int it=0;(it<60) && (hi-Lb>RTOL*Lb) && (glo<0.0);it++){
    const double Q=bis?((lo>0.0)?sqrt(lo*hi):0.5*(lo+hi)):(lo*ghi-hi*glo)/(ghi-glo);
    const double Tq=T(Q*w,&buf[0]),g=Q-Tq;
    if (g>0.0){hi=Q; ghi=g; if (side==1){glo*=0.5;} side=1; if (Tq>Lb){Lb=Tq; LbT=Tq; bufLb=buf;}}
    else      {lo=Q; glo=g; if (side==-1){ghi*=0.5;} side=-1; if (Q>Lb){Lb=Q; LbT=Tq; bufLb=buf;}}
    bis=(hi-lo>0.5*wprev); wprev=hi-lo;                                         //(no halving, e.g. at a layer's room: bisect next)
  }
  return Put(Lb,LbT,&bufLb[0]);
}

//////////////////////////////////////////////////////////////////
/// \brief aquifer material just below the soil-zone base of SZX entry e (with :SoilZoneVadose), else NULL
//
const gwvz::vz_mat *CGroundwaterModel::VZGapMat(const int e) const
{
  if ((_vzMode<=0) || (e<0) || (e>=(int)_vzStrat.size()) || _vzStrat[e].mat.empty()){return NULL;}
  const int i=_vzStrat[e].Find(_szxZb[e]-1e-6);
  return (i>=0)?&_vzStrat[e].mat[i]:NULL;
}
