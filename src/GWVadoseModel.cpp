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
#include <set>
#include "Model.h"

static const double VZ_LMIN=0.05;   //[m] shortest unsaturated column (one node of the finest grid)
#include "GWUtil.h"
using gwutil::ToStr;

//////////////////////////////////////////////////////////////////
/// \brief materials of every exchange column, by MODFLOW layer (after SZBuildEntries)
//
void CGroundwaterModel::VZSetup()
{
  const size_t ne=_szxLink.size();
  _vzStrat.assign(ne,gwvz::vz_strat());
  _vzUse.assign(ne,0); _vzQ.assign(ne,0.0); _vzAbsQ.assign(ne,0.0); _vzRejV.assign(ne,0.0); _vzFp.assign(ne,0.0);
  _vzInV.assign(ne,0.0); _vzCapV.assign(ne,0.0); _vzCapL.assign(ne*_szNLayers,0.0);
  gwvz::SetAccuracy(_vzAcc); gwvz::SetReference((_vzMode==2)?1:0);
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
        ExitGracefullyIf(!(m.ths>0.0),("CGroundwaterModel: :SoilZoneVadose "+VZName()+": aquifer class "+C.name+" has POROSITY 0 (the saturated content "
                         "of the unsaturated columns); give every aquifer class the columns pass through a positive POROSITY").c_str(),BAD_DATA);
        if (!C.unsatGiven && !warned[P.layers[l].iclass]){warned[P.layers[l].iclass]=1; noUnsat+=(noUnsat==""?"":", ")+C.name;}
      }
      if (_vzMode!=2){ //Sy against the drainable porosity of the retention curve 1 m above the water table (soil_zone_columns_plan.md
        //1.5): a rising water table releases the column's water above equilibrium to MODFLOW, which stores it with Sy; with Sy
        //much smaller each step's rise releases more than the last (a feedback that stops only at the soil-zone base).
        //(KINEMATIC: its reference is ths - Sy, a rise releases at most Sy dh; the warning is meaningless)
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
    WriteWarning("CGroundwaterModel: :SoilZoneVadose "+VZName()+": the specific yield differs by more than a factor 2 from the drainable porosity "
                 "theta_s - theta(psi = -1 m) of POROSITY, CLAPP_B and AIR_ENTRY_PRESSURE"+(badSy!=""?(" for "+badSy):string(" in cells of the model and :SoilZoneAquiferClass"))+
                 "; a rising water table then releases more (or less) water from the unsaturated column than MODFLOW stores with Sy, and with Sy too small "
                 "the water table rises step by step to the soil-zone base",false);
  }
  if (noUnsat!=""){
    WriteWarning("CGroundwaterModel: :SoilZoneVadose "+VZName()+": aquifer classes "+noUnsat+" have no CLAPP_B and AIR_ENTRY_PRESSURE; the unsaturated column "
                 "uses 4 and 0.3 m (give both in :AquiferClasses)",false);
  }
  VZClassSetup();                                                                //the columns: one per depth class of every entry
}

//////////////////////////////////////////////////////////////////
/// \brief advances one unsaturated column over dt by the method of :SoilZoneVadose (gwvz::Step, StepKinematic, StepLag);
///  all share the column (grid, contents, equilibrium reference, MoveBottom) and return the flux into the water table
//
double CGroundwaterModel::VZColStep(gwvz::vz_column &c,const gwvz::vz_strat &S,const double dt,const double qtop,double &rej,gwvz::vz_step_info *inf,const bool noUp) const
{
  if (_vzMode==2){return gwvz::StepKinematic(c,S,dt,qtop,rej,inf);}
  if (_vzMode==3){return gwvz::StepLag(c,S,dt,qtop,rej,inf,_vzPar);}
  return gwvz::Step(c,S,dt,qtop,rej,inf,noUp);
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
/// \brief depth classes of a link (:SoilZoneVadoseClasses K): its DEM samples' soil-zone bases z (ascending) grouped into at
///  most K contiguous classes with the least spread (Fisher's optimal partition: the smallest sum over the classes of the
///  squared deviations from the class mean, by dynamic programming over the sorted bases). Equal bases (within 1e-9 m) are
///  never split, so K is also capped by the number of distinct bases. lo gets the class boundaries (lo[0] = 0 ... lo[K] = N).
//
static void VZGroups(const vector<double> &z,const int Kmax,vector<int> &lo)
{
  const int N=(int)z.size();
  lo.assign(1,0);
  if (N<=1){lo.push_back(N); return;}
  vector<char> cut(N+1,0); int nd=1;                                            //cut[j]: a class may start at sample j
  for (int j=1;j<N;j++){if (z[j]-z[j-1]>1e-9){cut[j]=1; nd++;}}
  const int K=max(1,min(Kmax,nd));
  if (K==1){lo.push_back(N); return;}
  vector<double> s1(N+1,0.0),s2(N+1,0.0);                                       //(prefix sums, about the first base: no cancellation)
  for (int j=0;j<N;j++){const double d=z[j]-z[0]; s1[j+1]=s1[j]+d; s2[j+1]=s2[j]+d*d;}
  auto sse=[&](const int a,const int b)->double{const double n=b-a,m=s1[b]-s1[a]; return max(s2[b]-s2[a]-m*m/n,0.0);}; //samples a..b-1
  const double INF=1e300;
  vector<vector<double> > D(K+1,vector<double>(N+1,INF)); vector<vector<int> > B(K+1,vector<int>(N+1,-1));
  D[0][0]=0.0;
  for (int k=1;k<=K;k++){
    for (int j=1;j<=N;j++){
      if ((j<N) && !cut[j]){continue;}                                           //(a class ends only where the next may start)
      for (int i=k-1;i<j;i++){
        if ((i>0) && !cut[i]){continue;}
        if (D[k-1][i]>=INF){continue;}
        const double v=D[k-1][i]+sse(i,j);
        if (v<D[k][j]){D[k][j]=v; B[k][j]=i;}
      }
    }
  }
  vector<int> b(K+1,N);
  for (int k=K;k>0;k--){b[k-1]=B[k][b[k]];}
  lo=b;
}

//////////////////////////////////////////////////////////////////
/// \brief depth classes of every exchange column (after VZSetup's materials): class columns, their area shares, tops and
///  samples. With K = 1 (or one distinct base) one column per entry, at the entry's mean soil-zone base, as without classes.
///  A class column starts at the mean base of its class, shifted with the entry's mean base (z_c = zb + mean_c - mean):
///  the area-weighted mean of the class tops is the entry's mean base, so the mean depth to the water table is kept and
///  K = 1 is the single column.
//
void CGroundwaterModel::VZClassSetup()
{
  const size_t ne=_szxLink.size();
  _vzOff.assign(ne+1,0); _vzcA.clear(); _vzcZt.clear(); _vzcZbs.clear(); _vzClsDev=0.0;
  vector<int> lo;
  for (size_t e=0;e<ne;e++){
    _vzOff[e]=(int)_vzcA.size();
    const vector<double> &z=_szxCell[e].zbs; const int N=(int)z.size();
    if ((_vzClasses<=1) || (N<=1)){lo.assign(1,0); lo.push_back(N);} else {VZGroups(z,_vzClasses,lo);}
    const int K=(int)lo.size()-1;
    if (K<=1){_vzcA.push_back(1.0); _vzcZt.push_back(_szxZb[e]); _vzcZbs.push_back(vector<double>()); continue;}
    double m=0.0; for (int q=0;q<N;q++){m+=z[q];} m/=N;
    for (int c=0;c<K;c++){
      double mc=0.0; for (int q=lo[c];q<lo[c+1];q++){mc+=z[q];} mc/=(lo[c+1]-lo[c]);
      _vzcA.push_back((double)(lo[c+1]-lo[c])/(double)N);
      _vzcZt.push_back(_szxZb[e]+(mc-m));
      _vzcZbs.push_back(vector<double>(z.begin()+lo[c],z.begin()+lo[c+1]));
    }
    //check: the class fractions, area-weighted, are the entry's (the same smooth step on the class's samples): heads from
    //below the lowest base's contact to above the highest base
    const double ds=_szxCell[e].ds;
    for (int t=0;t<=60;t++){
      const double h=z.front()-1.5*ds+(z.back()-z.front()+3.0*ds)*t/60.0;
      double Fs=0.0; for (int c=0;c<K;c++){Fs+=_vzcA[_vzOff[e]+c]*VZClassF((int)e,_vzOff[e]+c,h);}
      _vzClsDev=max(_vzClsDev,fabs(Fs-gwsz::SatFraction(h,z,ds)));
    }
  }
  _vzOff[ne]=(int)_vzcA.size();
  const size_t nc=_vzcA.size();
  _vzcCol.assign(nc,gwvz::vz_column()); _vzcUse.assign(nc,0); _vzcF0.assign(nc,0.0); _vzcOut.assign(nc,0.0); _vzcRq.assign(nc,0.0);
  _vzAbsN.assign(ne,0.0); _vzQN.assign(ne,0.0); _vzBackW.assign(ne,0.0); _vzPass.assign(ne,0.0); _vzRl.assign(ne,0.0);
}

//////////////////////////////////////////////////////////////////
/// \brief saturated fraction of class column i of SZX entry e at head h: gwsz::SatFraction on the class's own samples (same
///  smoothing depth), so that sum_c a_c F_c = F of the entry (up to round-off); one class: the entry's F
//
double CGroundwaterModel::VZClassF(const int e,const int i,const double h) const
{
  if (_vzOff[e+1]-_vzOff[e]<=1){return gwsz::SatFraction(h,_szxCell[e].zbs,_szxCell[e].ds);}
  return gwsz::SatFraction(h,_vzcZbs[i],_szxCell[e].ds);
}

//////////////////////////////////////////////////////////////////
/// \brief rule (i) of the water-table following: where the water table
///  reaches ground that was unsaturated at the start of the step, the outputs of the columns on that ground go back to the
///  soil zone (rejected recharge) instead of the aquifer: class by class R_q,c clip(F_c(h) - F0_c, 0, 1 - F0_c), R_q,c the
///  class's outputs per unit unsaturated share; with one class F_c = F at the connection head (_linkF). [m3/d]
//
double CGroundwaterModel::VZBackRate(const int e,const double h) const
{
  if ((_vzMode<=0) || (e<0) || (e+1>=(int)_vzOff.size())){return 0.0;}
  double r=0.0;
  const bool expl=(_szMaxIter==0) || (SZFreezeAt()<1);
  for (int i=_vzOff[e];i<_vzOff[e+1];i++){
    if (!_vzcUse[i]){continue;}
    const double F0=_vzcF0[i];
    if (F0>=1.0-1e-12){r+=max(_vzcOut[i],0.0)/_szTstep; continue;}             //(all saturated at the start: all of it, also in explicit coupling)
    if (expl || !(_vzcRq[i]>0.0)){continue;}                                     //(explicit coupling: no head to follow)
    const double F=VZClassF(e,i,h);                                              //(one class: the link's F at h; at the connection head = _linkF)
    r+=_vzcRq[i]*min(max(F-F0,0.0),1.0-F0);
  }
  return r;
}

//////////////////////////////////////////////////////////////////
/// \brief water above the reference contents (ths - Sy, what a falling water table leaves) in all columns [m3]
//
double CGroundwaterModel::VZStorage() const
{
  double v=0.0;
  for (size_t e=0;e+1<_vzOff.size();e++){
    const double a=_szColA[_szxCol[e]];
    for (int i=_vzOff[e];i<_vzOff[e+1];i++){if (_vzcCol[i].active){v+=gwvz::RelStorage(_vzcCol[i],_vzStrat[e])*(a*_vzcA[i]);}}
  }
  return v;
}


//////////////////////////////////////////////////////////////////
/// \brief columns at the start of the run: from the hotstart file (:GWVadoseColumns, or :GWVadoseClassColumns with depth
///  classes) where it has them; otherwise at their reference profile above the initial water table (no flow until recharge
///  arrives), or, with :SoilZoneVadoseInitialFlux on a cold start, at their steady profile for that flux. On a groundwater
///  hotstart (heads, soil-zone or unsaturated columns read from the .rvc) the initial flux is ignored: a column the hotstart
///  does not hold formed at rest in the run that wrote it (or did not exist), so it starts at rest here too.
//
void CGroundwaterModel::VZInit()
{
  int nHot=0,nSteady=0,nRest=0,nBack=0;
  const bool cls=(_vzClasses>1);
  const bool gwHot=!_hotHeads.empty() || !_hotSzCols.empty() || !_hotVz.empty() || !_hotVzc.empty();
  const double q0=_vzQ0/MM_PER_METER;                                            //[m/d]
  //the methods count the columns' water from different reference profiles (KINEMATIC from ths - Sy, RICHARDS and LAG from
  //the hydrostatic profile): contents another reference wrote would be other water here, so such lines are not used
  const bool refOther=!_hotVzMethod.empty() && ((_hotVzMethod=="KINEMATIC")!=(_vzMode==2));
  set<pair<long long,pair<int,int> > > used;                                    //(HRU ID, cell, class) of the lines used
  for (size_t e=0;e<_szxLink.size();e++){
    const gw_link &L=_aLinks[_szxLink[e]]; const int ic=L.ic;
    if (!IsExternal() && _colCHDTop[ic]){ //a generated model's cell held as a specified head at the start (a :BoundaryHead series): no
      for (int i=_vzOff[e];i<_vzOff[e+1];i++){                                  //  column (one forms, at rest, when the cell is released)
        _vzcCol[i]=gwvz::vz_column(); _vzcCol[i].ztop=_vzcZt[i]; _vzcCol[i].zbot=_vzcZt[i];
      }
      continue;
    }
    double h=_X[_reduced[TopNode(ic)]]*_hf; const double Lmin=VZ_LMIN;
    if (!(std::fabs(_X[_reduced[TopNode(ic)]])<1e29)){h=WaterTable(ic);}       //(HDRY is in the model's units: a dry top cell without NEWTON, the water table below it)
    const double hF=h;                                                           //(the initial head: the saturated fractions of the first step)
    const gwvz::vz_strat &S=_vzStrat[e];
    double zlow=-1e300; for (int l=_nlay-1;l>=0;l--){const int nd=UserNode(l,ic); if (_idomain[nd]>0){zlow=_botm[nd]; break;}}
    h=max(h,zlow);                                                               //(columns not below the aquifer's bottom)
    const long long id=_pModel->GetHydroUnit(L.k)->GetHRUID();
    const int nK=_vzOff[e+1]-_vzOff[e];
    //the hotstart lines of this link: of its own kind (:GWVadoseColumns with one class, :GWVadoseClassColumns with K > 1);
    //lines of the other kind were written with another number of classes and are not used (nor lines written by a method
    //of another reference profile). With depth classes a class column inactive at the end of the run that wrote the file
    //has no line (the writer skips it; the water table fell below its base in the last step, and the column forms here, at
    //rest, as it would have in the continuous run's next step). Lines not used are reported after the loop
    const vector<double> *v1=NULL; const map<int,vector<double> > *vc=NULL;
    if (!refOther && !cls){
      map<long long,map<int,vector<double> > >::const_iterator it=_hotVz.find(id);
      if (it!=_hotVz.end()){map<int,vector<double> >::const_iterator jt=it->second.find(ic+1); if (jt!=it->second.end()){v1=&jt->second;}}
    }
    else if (!refOther){
      map<long long,map<int,map<int,vector<double> > > >::const_iterator kt=_hotVzc.find(id);
      if (kt!=_hotVzc.end()){map<int,map<int,vector<double> > >::const_iterator jt=kt->second.find(ic+1); if (jt!=kt->second.end()){vc=&jt->second;}}
    }
    //the column tops written with the lines (:GWVadoseColumnTops; files without them are not checked): a top that differs by
    //more than 1 mm was written for other columns (another :SoilZoneSubgrid, ground or soil-zone thickness)
    const map<int,vector<double> > *vt=NULL;
    {map<long long,map<int,map<int,vector<double> > > >::const_iterator tt=_hotVzt.find(id); if (tt!=_hotVzt.end()){map<int,map<int,vector<double> > >::const_iterator ut=tt->second.find(ic+1); if (ut!=tt->second.end()){vt=&ut->second;}}}
    for (int i=_vzOff[e];i<_vzOff[e+1];i++){
      gwvz::vz_column &c=_vzcCol[i]; const double zb=_vzcZt[i];
      c=gwvz::vz_column(); c.ztop=zb; c.zbot=h;
      vector<double> v; bool have=false;
      if (v1!=NULL){if (v1->size()>=2){v=*v1; have=true;}}
      else if (vc!=NULL){
        map<int,vector<double> >::const_iterator ct=vc->find(i-_vzOff[e]+1);
        if ((ct!=vc->end()) && (ct->second.size()>=3) && ((int)std::floor(ct->second[0]+0.5)==nK)){v.assign(ct->second.begin()+1,ct->second.end()); have=true;}
      }
      if (have && (vt!=NULL)){
        map<int,vector<double> >::const_iterator ct=vt->find(i-_vzOff[e]+1);
        if ((ct!=vt->end()) && (ct->second.size()>=1) && (std::fabs(ct->second[0]-zb)>1e-3)){have=false;}
      }
      if (have){
        gwvz::Build(c,S,zb,v[0],_vzN);                                            //zbot, then the contents top down
        if (c.active && (c.th.size()==v.size()-1)){
          for (size_t q=0;q<c.th.size();q++){ //the content as written; floored only at what the method can hold: air-dry
            const gwvz::vz_mat &M=S.mat[c.im[q]];                                 //(RICHARDS, LAG), 0 (KINEMATIC, whose contents
            double th=max(v[q+1],(gwvz::GetReference()==1)?0.0:gwvz::Theta(M,-1e5)); // go down to ths - Sy, below air-dry when Sy
            if (!(M.Ss>0.0)){th=min(th,M.ths);}                                   // is near the porosity), and at most saturation
            c.th[q]=th; c.psi[q]=gwvz::PsiOf(M,th);                               // without specific storage
          }
          if (!c.th.empty()){c.dryTop=(c.psi[0]<=-0.999999*1e5); c.satTop=(c.psi[0]>=-1e-9);} //(the top's memory: a top held at air-dry or
                                                                                 // saturation by the last step, as in a continuous run)
          nHot++; used.insert(make_pair(id,make_pair(ic+1,cls?(i-_vzOff[e]+1):1)));
          continue;
        }
        c=gwvz::vz_column(); c.ztop=zb; c.zbot=h;                                //(another node count: the column starts at its reference;
                                                                                 //  a water table now at the base: no column)
      }
      if (zb-h>Lmin){
        gwvz::Build(c,S,zb,h,_vzN);
        if (c.active && (q0>0.0)){
          if (gwHot){nRest++;}
          else {
            //the flux the column receives: the unsaturated share of q, (1 - F) q, F the saturated fraction of its class (one
            //class: the link's) at the initial head, as VZStep gives it in the first step
            const double qc=q0*(1.0-min(max(VZClassF((int)e,i,hF),0.0),1.0));   //[m/d]
            if (qc>0.0){
              double fw=qc;
              if (_vzMode==2){fw=gwvz::InitSteadyKinematic(c,S,qc);}          //(steady starts of the column library)
              else if (_vzMode==3){gwvz::InitSteadyLag(c,S,qc,_vzPar);}
              else {gwvz::InitSteadyRichards(c,S,qc,&fw);}
              if (fw<qc*(1.0-1e-9)){nBack++;}                                    //(more than the column can carry: backed up, its top rejects the rest)
              nSteady++;
            }
          }
        }
      }
    }
  }
  _vzStore0=_vzStoreVol=VZStorage();
  if (nHot>0){WriteAdvisory("CGroundwaterModel: "+to_string(nHot)+" unsaturated columns read from the hotstart file",false);}
  //the hotstart lines no column used: written for other columns (another method's reference profile, node count, number of
  //classes or column top: those columns start at rest), or on cells that are no longer exchange links of a soil-zone HRU
  //(now specified heads, another HRU or grid): their water is dropped, warned about and reported in the run totals
  int nUn=0,nBad=0; double vUn=0.0; bool vol=!_hotVzt.empty();
  map<long long,set<int> > linkCells; for (size_t e=0;e<_szxLink.size();e++){const gw_link &L=_aLinks[_szxLink[e]]; linkCells[_pModel->GetHydroUnit(L.k)->GetHRUID()].insert(L.ic+1);}
  vector<pair<long long,pair<int,int> > > lines;
  for (map<long long,map<int,vector<double> > >::const_iterator it=_hotVz.begin();it!=_hotVz.end();it++){
    for (map<int,vector<double> >::const_iterator jt=it->second.begin();jt!=it->second.end();jt++){lines.push_back(make_pair(it->first,make_pair(jt->first,1)));}
  }
  for (map<long long,map<int,map<int,vector<double> > > >::const_iterator it=_hotVzc.begin();it!=_hotVzc.end();it++){
    for (map<int,map<int,vector<double> > >::const_iterator jt=it->second.begin();jt!=it->second.end();jt++){
      for (map<int,vector<double> >::const_iterator kt=jt->second.begin();kt!=jt->second.end();kt++){lines.push_back(make_pair(it->first,make_pair(jt->first,kt->first)));}
    }
  }
  for (size_t j=0;j<lines.size();j++){
    if (used.count(lines[j])){continue;}
    nUn++;
    map<long long,set<int> >::const_iterator lt=linkCells.find(lines[j].first);
    const bool atLink=(lt!=linkCells.end()) && lt->second.count(lines[j].second.first) && (IsExternal() || !_colCHDTop[lines[j].second.first-1]);
    if (atLink){nBad++;}
    map<long long,map<int,map<int,vector<double> > > >::const_iterator tt=_hotVzt.find(lines[j].first);
    if (tt==_hotVzt.end()){continue;}
    map<int,map<int,vector<double> > >::const_iterator ut=tt->second.find(lines[j].second.first); if (ut==tt->second.end()){continue;}
    map<int,vector<double> >::const_iterator vt=ut->second.find(lines[j].second.second); if ((vt!=ut->second.end()) && (vt->second.size()>=2)){vUn+=vt->second[1];}
  }
  _vzHotUnusedN=nUn; _vzHotUnusedV=vUn; _vzHotBadN=nBad;
  if (nUn>0){
    string why;
    if (nBad>0){why+=to_string(nBad)+" written for other columns ("+(refOther?(":SoilZoneVadose "+_hotVzMethod+", whose reference profile differs from this run's: they start at rest"):
                     string("another :SoilZoneVadose node count, another :SoilZoneVadoseClasses, or another column top: :SoilZoneSubgrid, ground or soil-zone thickness?; they start at rest"))+")";}
    if (nUn>nBad){why+=string((nBad>0)?", ":"")+to_string(nUn-nBad)+" on cells that are no longer exchange links of a soil-zone HRU (specified heads, another HRU or grid?)";}
    WriteWarning("CGroundwaterModel: "+to_string(nUn)+" unsaturated columns of the hotstart file are not used: "+why+"; their water"+(vol?(", "+ToStr(vUn)+" m3 above their reference profiles,"):string(""))+
                 " is dropped (reported in the run totals of GWModelSummary.txt)",false);
  }
  if (nSteady>0){WriteAdvisory("CGroundwaterModel: "+to_string(nSteady)+" unsaturated columns start at their steady profile for :SoilZoneVadoseInitialFlux "+ToStr(_vzQ0)+" mm/d "
                               "(the unsaturated share (1 - F) of it each column receives)",false);}
  if (nBack>0){WriteWarning("CGroundwaterModel: "+to_string(nBack)+" unsaturated columns cannot carry the flux of :SoilZoneVadoseInitialFlux "+ToStr(_vzQ0)+" mm/d (their layers "
                            "pass less with a saturated top); they start at their backed-up steady profile for the flux they can pass, and their tops reject the rest",false);}
  if (nRest>0){WriteWarning("CGroundwaterModel: :SoilZoneVadoseInitialFlux ignored: the run starts from a groundwater hotstart; "+to_string(nRest)+
                            " unsaturated columns the hotstart file does not hold start at rest (their reference profile)",false);}
}

//////////////////////////////////////////////////////////////////
/// \brief the unsaturated columns over this step, with the water table at the start of the step (soil_zone_columns_plan.md
///  section 1.3): every column follows the water table (water below a risen one goes to the aquifer this step), takes the
///  unsaturated share of its recharge at its top (less the capillary rise into the soil zone), and gives its bottom flux to
///  MODFLOW's recharge (SZRecharge). Solved once per step, before MODFLOW's solve. With depth classes
///  (:SoilZoneVadoseClasses) an exchange column has one unsaturated column per class (area a a_c, top at its class's base,
///  inflow the class's unsaturated share a_c (1 - F_c)); the classes of a link are processed in turn and SHARE the limits of
///  the link's source cell (what earlier classes took is subtracted) and the room of its soil-zone layers; the per-link
///  arrays (_vzQ, _vzAbsQ, _vzRejV, _vzInV, _vzCapV, _vzUp, _vzUpA, _vzCapL, _vzUse) are sums over its classes.
//
void CGroundwaterModel::VZStep(const double &tstep)
{
  const int n=_szNLayers;
  const bool cap=(_szCapRise==1) || ((_szCapRise<0) && (_vzMode>0));
  _vzInVol=_vzOutVol=_vzCapVol=0.0; _vzUpQ=0.0;
  vector<double> QL(max(n,1),0.0),QW(max(n,1),0.0),capLV(max(n,1),0.0),capQL(max(n,1),0.0);
  for (size_t i=0;i<_aqHRUs.size();i++)
  {
    int k=_aqHRUs[i]; if (!_szOn[k]){continue;}
    const soil_struct *sp=_pModel->GetHydroUnit(k)->GetSoilProps(_szNLayers-1);   //the bottom soil-zone layer's soil
    vector<double> pend;                                                         //[m3] releases of columns that went (or carried nothing)
    for (size_t jj=0;jj<_szxOfHRU[k].size();jj++){
      const int e=_szxOfHRU[k][jj],li=_szxLink[e],c=_szxCol[e]; const gw_link &L=_aLinks[li];
      const double a=_szColA[c];
      const int nK=_vzOff[e+1]-_vzOff[e];
      _vzUse[e]=0; _vzQ[e]=0.0; _vzAbsQ[e]=0.0; _vzRejV[e]=0.0; _vzInV[e]=0.0; _vzCapV[e]=0.0; _vzUp[e]=0.0; _vzUpA[e]=0.0; _vzFp[e]=0.0;
      _vzAbsN[e]=0.0; _vzQN[e]=0.0; _vzBackW[e]=0.0; _vzPass[e]=0.0;
      if (_vzUsed.size()!=_szxLink.size()){_vzUsed.assign(_szxLink.size(),0.0);}
      _vzUsed[e]=0.0;
      for (int j=0;j<n;j++){_vzCapL[e*n+j]=0.0; capLV[j]=0.0; capQL[j]=0.0;}
      for (int ii=_vzOff[e];ii<_vzOff[e+1];ii++){_vzcUse[ii]=0; _vzcF0[ii]=0.0; _vzcOut[ii]=0.0; _vzcRq[ii]=0.0;}
      if (!(a>0.0)){continue;}
      if (_colCHDTop[L.ic]){ //a specified head this step (an existing model's own, known once MODFLOW has read the period): its recharge
        //goes there (SZRecharge), and the link's unsaturated columns end, their water above their reference profiles to the
        //specified head (as recharge onto it, _vzAbsQ); a column forms again, at rest, when the cell is no longer one
        for (int ii=_vzOff[e];ii<_vzOff[e+1];ii++){
          gwvz::vz_column &col=_vzcCol[ii]; if (!col.active){continue;}
          const double R=gwvz::RelStorage(col,_vzStrat[e])*(a*_vzcA[ii]);       //[m3]
          const double zt0=col.ztop; col=gwvz::vz_column(); col.ztop=zt0; col.zbot=zt0;
          _vzAbsQ[e]+=R/tstep; _vzOutVol+=R; _vzChdV+=R; if (R!=0.0){_vzChdN++;} //(a column at rest holds no water: not counted)
        }
        continue;
      }
      const gwvz::vz_strat &S=_vzStrat[e];
      const double h=_X[_reduced[TopNode(L.ic)]]*_hf,Lmin=VZ_LMIN;
      if (!(std::fabs(_X[_reduced[TopNode(L.ic)]])<1e29)){continue;}           //(HDRY is in the model's units: a top cell out of MODFLOW's
                                                                                 // solve, HDRY without NEWTON: the columns are left as they are until the cell rewets)
      double zlow=-1e300; for (int l=_nlay-1;l>=0;l--){const int nd=UserNode(l,L.ic); if (_idomain[nd]>0){zlow=_botm[nd]; break;}}
      const double share=VZLinkShare(k,li);
      double used=0.0;                                                           //[m3/d] taken from the source cell by this link's columns so far
      double capV=0.0;                                                           //[m3] the link's capillary rise so far
      for (int ii=_vzOff[e];ii<_vzOff[e+1];ii++){
        gwvz::vz_column &col=_vzcCol[ii];
        const double zb=_vzcZt[ii],ac=_vzcA[ii],Ac=a*ac;                         //(one class: the link's mean base and area)
        double absQ=0.0;                                                         //[m3/d] this column's release (> 0) or uptake (< 0)
        //follow the water table: a rising one takes the excess over equilibrium below it, a falling one opens cells at
        //equilibrium; a column shorter than Lmin (5 cm) is the soil zone's contact exchange and capillary rise instead
        if (col.active){
          //a rising water table takes the column's water above equilibrium below it, at most what MODFLOW stores for that
          //rise (Sy dh): with Sy smaller than the drainable porosity each step's rise would release more than the last (a
          //feedback with gain eps/Sy > 1, design note 1.5); the column's bottom then follows more slowly, releasing the rest
          //in the next steps (with consistent parameters eps <= Sy and this does not apply)
          const double zt=(zb-h<=Lmin)?zb:max(h,zlow);                           //(not below the aquifer's bottom: a dry column's head can be)
          if (zt>col.zbot+1e-9){ //a rising water table (to the base: the column goes)
            const double z0=col.zbot;
            const int tn=TopNode(L.ic); const double ztc=CellTop(tn/_ncpl,L.ic);
            const double capS=min(gwvz::SyRise(S,z0,zt),gwvz::SyRise(S,h,min(zb,ztc))); //[m] Sy dh across the layers the rise crosses, and the
                                                                                 // room the top cell has left above the water table
            gwvz::vz_column t=col; const double abT=gwvz::MoveBottom(t,S,zt,_vzN);
            const double capU=_szxUpCHD[e]?1e300:(SZUpCap(e,tstep)-used)*tstep/Ac; //[m] a deficit below it (a column dried by evaporation or
            if ((zt<zb) && (abT>capS*(1.0+1e-9))){                               // capillary rise) is filled out of the source cell, at most what
              double lo=z0,hi=min(zt,col.ztop);                                  // it can still give this step (SZUpCap less what the link's other
              for (int it=0;it<60;it++){ //(the excess below a level grows with the level)       classes took)
                const double zm=0.5*(lo+hi); t=col;
                if (gwvz::MoveBottom(t,S,zm,_vzN)>capS){hi=zm;} else {lo=zm;}
                if (hi-lo<1e-7){break;}
              }
              const double zn=(lo>col.ztop-0.01)?col.ztop:lo;                    //(a remainder thinner than 1 cm: the column goes)
              absQ=gwvz::MoveBottom(col,S,zn,_vzN)*Ac/tstep;
              if ((_vzSlowed++==0) && (_vzMode!=2)){WriteWarning("CGroundwaterModel: :SoilZoneVadose "+VZName()+": a rising water table met more water above equilibrium in an unsaturated "
                                               "column than MODFLOW stores with Sy; the column releases at most Sy dh per step (Sy below the drainable porosity?)",false);}
            }
            else if ((abT<0.0) && (-abT>capU)){                                  //(the bottom then follows more slowly; the column stays)
              double lo=z0,hi=min(zt,col.ztop-1e-3);
              for (int it=0;it<60;it++){
                const double zm=0.5*(lo+hi); t=col;
                if (-gwvz::MoveBottom(t,S,zm,_vzN)>capU){hi=zm;} else {lo=zm;}
                if (hi-lo<1e-7){break;}
              }
              absQ=(lo>z0+1e-9)?gwvz::MoveBottom(col,S,lo,_vzN)*Ac/tstep:0.0;
            }
            else if ((zt>=zb) && (abT>capS*(1.0+1e-9))){                       //(rule (ii): a column that vanishes releases at most Sy dh to
              col=t; absQ=capS*Ac/tstep; _vzBackW[e]+=(abT-capS)*Ac;             // the aquifer; the rest of its water goes back to the soil zone,
            }                                                                    // as rejected recharge)
            else {col=t; absQ=abT*Ac/tstep;}
          }
          else if (zt<col.zbot-1e-9){ //a falling water table: the new nodes start at equilibrium, the water they take comes out of the
            //cell holding the water table, at most what it can still give this step (SZUpCap: none without such a cell, the bottom
            //then stays); the column's bottom follows more slowly otherwise
            const double capF=(_szxUpCHD[e]?1e300:(SZUpCap(e,tstep)-used)*tstep/Ac); //[m]
            gwvz::vz_column t=col; const double abT=gwvz::MoveBottom(t,S,zt,_vzN);
            if (-abT>capF){
              double lo=zt,hi=col.zbot;                                          //(lo: too far down; hi: not moved)
              for (int it=0;it<60;it++){
                const double zm=0.5*(lo+hi); t=col;
                if (-gwvz::MoveBottom(t,S,zm,_vzN)>capF){lo=zm;} else {hi=zm;}
                if (hi-lo<1e-7){break;}
              }
              absQ=(hi<col.zbot-1e-9)?gwvz::MoveBottom(col,S,hi,_vzN)*Ac/tstep:0.0;
            }
            else {col=t; absQ=abT*Ac/tstep;}
          }
          //(unmoved: nothing)
          if (absQ<0.0){
            used-=absQ; _vzAbsN[e]-=absQ;
            if ((_szxUp[e]>=0) && (_szxUp[e]!=TopNode(L.ic)) && !_szxUpCHD[e]){_vzUpA[e]-=absQ; _vzUpQ-=absQ;}
          }
          _vzAbsQ[e]+=absQ;
        }
        else if (zb-max(h,zlow)>Lmin){
          gwvz::Build(col,S,zb,max(h,zlow),_vzN);                                //(at equilibrium: no excess, nothing exchanged; not below the
        }                                                                        // aquifer's bottom)
        if (!col.active){pend.push_back(absQ*tstep); _vzcOut[ii]=absQ*tstep; continue;} //(booked after the HRU's columns, as without classes)
        //inflow: the unsaturated share of the class's recharge (the saturated share is rejected, as on a link without columns, with this
        //saturated fraction for the whole step)
        const double F=(nK==1)?_linkF[li]:VZClassF(e,ii,_szxHc[e]);
        const bool negP=(_szRchV[k]<0.0);                                        //(a negative recharge is not split by F)
        double P=negP?_szRchV[k]*share*ac:_szRchV[k]*share*ac*(1.0-F);           //[m3] over the step
        if (negP){ //a negative recharge (a process taking water out of GROUNDWATER) never goes through the column (it would take
          //water out of the column's top, booked as negative capillary rise): passed on to the cell, out of its source's water,
          //at most what the top cell can still give this step (counted in used/_vzUsed, so the capillary rise and the repayment
          //do not draw the same water again); the rest is taken back out of the soil zone as a negative rejected recharge
          double capN=0.0;
          if (_szxUpCHD[e] && (_szxUp[e]==TopNode(L.ic))){capN=1e300;}
          else if (_szxUp[e]==TopNode(L.ic)){capN=max(SZUpCap(e,tstep)-used,0.0);}
          const double passV=max(P,-capN*tstep);                                 //[m3]
          _vzPass[e]+=passV/tstep; used-=passV/tstep; _vzRejV[e]+=P-passV; P=0.0;
        }
        _vzcF0[ii]=F; _vzcUse[ii]=1; _vzFp[e]+=ac*F; _vzInV[e]+=P;
        //capillary rise into the soil zone (SZRiseLayers, as from the water table where there is no column): the steady
        //upward flux from the column's top node (its suction) into every layer below field capacity, on the unsaturated
        //share. Per step at most the larger of the top node's water down to the suction of the lowest receiving layer and
        //what the water table itself could supply with its drawdown: as the column shortens the top node empties and the
        //rise becomes the water table's (the two meet at the shortest column, no jump), and what the column cannot pass
        //on within the step comes back as an unmet upward demand below. The classes of a link share its layers' room
        double R=0.0,capC=0.0; bool rose=false;
        if (cap && (sp!=NULL) && VZOwnRise()){                                    //(other methods: from the water table, SZCapRiseWT)
          double psL=0.0;
          const double Ae=Ac*(1.0-min(max(F,0.0),1.0));
          const double Qc=SZRiseLayers(k,c,zb,col.zc[0],max(-col.psi[0],0.0),a,Ae,0.0,tstep,&QL[0],&psL,NULL,(nK>1)?&capQL[0]:NULL);
          if (Qc>0.0){
            const gwvz::vz_mat &m0=S.mat[col.im[0]];
            const double avail=max(col.th[0]-gwvz::Theta(m0,-min(psL,1e5)),0.0)*col.dz[0]*Ac; //[m3] the top node down to the lowest layer's suction
            const int up=_szxUp[e];                                              //(as SZCapRiseWT: from the cell holding the water table; only
            double Qw=0.0;                                                       // needed when the top node cannot give the rise)
            if ((up>=0) && (Qc>avail/tstep)){
              const double hu=_X[_reduced[up]]*_hf,Syu=_szxUpCHD[e]?0.0:max(SZNodeSy(up),0.01);
              Qw=SZRiseLayers(k,c,zb,hu,0.0,a,Ae,Syu,tstep,&QW[0],NULL,VZGapMat(e,(nK>1)?zb:-HUGE_VAL),(nK>1)?&capQL[0]:NULL);
              if (!_szxUpCHD[e]){Qw=min(Qw*_szxUpF[e],SZUpCap(e,tstep));}          //(as SZCapRiseWT)
            }
            const double Q=min(Qc,max(avail/tstep,Qw));                          //[m3/d]
            if (Q>0.0){
              R=Q/Ac;                                                            //[m/d] out of the column's top
              capC=Q*tstep; rose=true;
              if (nK==1){_vzCapV[e]=capC; _vzCapVol+=capC; for (int j=0;j<n;j++){_vzCapL[e*n+j]=QL[j]/Qc;}} //(added to the layers in SZApply,
              else {_vzCapVol+=capC; for (int j=0;j<n;j++){capQL[j]+=QL[j]*Q/Qc;}}                           // after the columns' sync)
            }
          }
          if (rose && (nK>1)){for (int j=0;j<n;j++){QL[j]=QL[j]/Qc;}}             //(layer fractions of this class's rise)
          (void)Qc;
        }
        double rej=0.0,rejX=0.0; gwvz::vz_step_info inf;                         //(rejX [m3]: the negative rejected recharge below)
        const bool noUp=(_szxUp[e]<0);                                          //(no cell holds any water table: no upward
        double qb=VZColStep(col,S,tstep,P/(Ac*tstep)-R,rej,&inf,noUp);          // flow, which MODFLOW could not supply)
        if (qb<0.0){ //upflow through the bottom: out of the cell holding the water table, scaled down as it nears drying and at
          //most half the water it holds above its bottom under the column per step (none without such a cell), less what the
          //link's columns already took; the part MODFLOW does not supply is taken back out of the column, bottom node first,
          //down to air-dry
          double allow=0.0;                                                      //[m/d]
          if (!noUp && _szxUpCHD[e]){allow=-qb;}                                 //(a specified head supplies all of it)
          else if (!noUp){                                                       //(the source cell's limit, less what was taken)
            allow=min(_szxUpF[e]*(-qb),max(SZUpCap(e,tstep)-used,0.0)/Ac);
          }
          double need=(-qb-allow)*tstep;                                         //[m]
          for (size_t ib=col.th.size();(ib-->0) && (need>0.0);){
            const gwvz::vz_mat &mb=S.mat[col.im[ib]];
            const double take=min(need,max(col.th[ib]-gwvz::Theta(mb,-1e5),0.0)*col.dz[ib]);
            if (take>0.0){col.th[ib]-=take/col.dz[ib]; col.psi[ib]=gwvz::PsiOf(mb,col.th[ib]); need-=take; qb+=take/tstep;}
          }
          if ((need>0.0) && rose){ //what the bottom node cannot give either is taken off what the column delivered at its top: its
            //capillary rise (the column's storage change is unchanged)
            double &cv=(nK==1)?_vzCapV[e]:capC;
            const double t2=min(need,cv/Ac);
            if (t2>0.0){cv-=t2*Ac; _vzCapVol-=t2*Ac; need-=t2; qb+=t2/tstep;}
          }
          if ((need>0.0) && !noUp && !_szxUpCHD[e] && (_szxUp[e]==TopNode(L.ic))){ //then the top cell as far as it can still give
            //(its limit less what the link's columns took, this upflow's allowance included)
            const double topLeft=max(SZUpCap(e,tstep)-used-allow*Ac,0.0)*tstep/Ac; //[m]
            need-=min(need,topLeft);                                             //(stays in qb: upflow out of the top cell)
          }
          if (need>0.0){ //the rest handed back to the soil zone as a negative rejected recharge (taken out by SZApply, possibly
            //booked as not held), so that the cell never gives more than it can
            rejX=need*Ac; _vzRejV[e]-=rejX; qb+=need/tstep; need=0.0;
          }
        }
        _vzSub+=inf.substeps; _vzNewton+=inf.newton;
        if (nK==1){
          if (rej<0.0){_vzCapV[e]=max(_vzCapV[e]+rej*a,0.0); _vzCapVol+=rej*a; rej=0.0;} //(an upward demand the column could not supply)
        }
        else {
          if (rej<0.0){capC=max(capC+rej*Ac,0.0); _vzCapVol+=rej*Ac; rej=0.0;}
          if (rose){capV+=capC; for (int j=0;j<n;j++){capLV[j]+=capC*QL[j];}}
        }
        const double qV=qb*Ac;                                                   //[m3/d]
        _vzRejV[e]+=rej*Ac; _vzQ[e]+=qV; _vzUse[e]=1;
        if (qb<0.0){
          used-=qV; _vzQN[e]-=qV;                                                //(the source cell's water: the link's later classes have less)
          if ((_szxUp[e]>=0) && (_szxUp[e]!=TopNode(L.ic)) && !_szxUpCHD[e]){_vzUp[e]-=qV; _vzUpQ-=qV;} //(upflow from a lower cell: SZW_RAVEN;
        }                                                                        // otherwise, also a remainder without a source cell, the top cell's recharge)
        _vzcOut[ii]=(qV+absQ)*tstep;
        _vzcRq[ii]=(F<1.0)?max(qV+absQ,0.0)/(1.0-F):0.0;                         //[m3/d] the outputs per unit unsaturated share (rule (i))
        _vzInVol+=P-rej*Ac+rejX; _vzOutVol+=(qV+absQ)*tstep;                         //(inflow the column took; its rejection is in the rejected recharge)
      }
      if ((nK>1) && (capV>0.0)){_vzCapV[e]=capV; for (int j=0;j<n;j++){_vzCapL[e*n+j]=capLV[j]/capV;}}
      _vzUsed[e]=used;                                                           //(what the link's columns took from the source cell: SZNegCaps)
    }
    for (size_t q=0;q<pend.size();q++){_vzOutVol+=pend[q];}
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
      const int nK=(!ch && (_vzMode>0) && (e+1<(int)_vzOff.size()))?_vzOff[e+1]-_vzOff[e]:1;
      if (!ch && VZOwnRise() && (nK==1) && (e+1<(int)_vzOff.size()) && _vzcCol[_vzOff[e]].active){continue;}  //(the column's own rise, VZStep)
      const int up=ch?TopNode(L.ic):_szxUp[e];                                  //(the cell holding the water table: the top cell, or one below
      if (up<0){continue;}                                                       // a dried top cell through SZW_RAVEN; none: no rise)
      const bool chu=ch || _szxUpCHD[e];                                         //(a specified head supplies it, top cell or below)
      const double h=_X[_reduced[up]]*_hf;
      vector<double> &CL=sc?_szcCapL:_szxCapL;
      const double F=min(max(_linkF[li],0.0),1.0); if (F>=1.0){continue;}
      const double z0=sc?SZBase(L.k,L.ic):_szxZb[e];
      const double Sy=chu?0.0:max(SZNodeSy(up),0.01);                           //(a specified head does not draw down)
      double Qs=0.0;
      if (nK==1){Qs=SZRiseLayers(k,c,z0,h,0.0,a,a*(1.0-F),Sy,tstep,&CL[e*n],NULL,sc?NULL:VZGapMat(e));}
      else { //depth classes: from the water table to each class's own base, on its unsaturated share a a_c (1 - F_c), where the
        //class has no column (RICHARDS) or for every class (the other methods); the classes share the layers' room, of which
        //the columns' own rise (VZStep) already took its part
        vector<double> tk(n,0.0),QLc(n,0.0);
        if (e<(int)_vzCapV.size()){for (int j=0;j<n;j++){tk[j]=_vzCapV[e]*_vzCapL[e*n+j]/tstep;}}
        for (int j=0;j<n;j++){CL[e*n+j]=0.0;}
        for (int ii=_vzOff[e];ii<_vzOff[e+1];ii++){
          if (VZOwnRise() && _vzcCol[ii].active){continue;}
          const double Fc=min(max(VZClassF(e,ii,_szxHc[e]),0.0),1.0),Ae=a*_vzcA[ii]*(1.0-Fc);
          if (!(Ae>0.0)){continue;}
          const double Qc=SZRiseLayers(k,c,_vzcZt[ii],h,0.0,a,Ae,Sy,tstep,&QLc[0],NULL,VZGapMat(e,_vzcZt[ii]),&tk[0]);
          if (!(Qc>0.0)){continue;}
          Qs+=Qc; for (int j=0;j<n;j++){CL[e*n+j]+=QLc[j]; tk[j]+=QLc[j];}
        }
      }
      if (!chu){ //out of an aquifer cell: fades out as it nears drying, and at most half the water it holds above its bottom
        double cap=SZUpCap(e,tstep);                                             // under the link per step (a thin cell would be pumped dry), less
        if ((_vzMode>0) && (e<(int)_vzUsed.size())){                             // what the link's columns took out of it in this step
          cap=max(cap-_vzUsed[e],0.0);                                           // (uptakes of moving bottoms, passed negative recharge and
        }                                                                        //  the columns' bottom upflow, as SZNegCaps counts them)
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
                                       const double Ae,const double Sy,const double tstep,double *QL,double *psisLow,const gwvz::vz_mat *gap,const double *taken) const
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
    room[j]=(Lj.phifc-phi)/MM_PER_METER*a/tstep;
    if (taken!=NULL){room[j]-=taken[j]; if (!(room[j]>0.0)){room[j]=0.0; continue;}} //(other depth classes filled it this step)
    d0[j]=zm-zs; any=true;
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
const gwvz::vz_mat *CGroundwaterModel::VZGapMat(const int e,const double zb) const
{
  if ((_vzMode<=0) || (e<0) || (e>=(int)_vzStrat.size()) || _vzStrat[e].mat.empty()){return NULL;}
  const int i=_vzStrat[e].Find(((zb==-HUGE_VAL)?_szxZb[e]:zb)-1e-6);
  return (i>=0)?&_vzStrat[e].mat[i]:NULL;
}
