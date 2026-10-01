/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: an existing MODFLOW 6 model coupled to Raven.
  The user's simulation is copied into the output folder and kept intact except for the time discretization (split
  into Raven time steps), the packages Raven takes over (removed) and two packages added (RCH_RAVEN, DRN_RAVEN).
  HRUs are linked to the cells of the model's top active layer through overlap areas - from cell polygons, or from an
  :OverlapWeights table as in the MODFLOW-USG coupling - normalised per HRU so that each HRU's volume is conserved.
----------------------------------------------------------------*/
#include "GroundwaterModel.h"
#include "Model.h"
#include "GWUtil.h"
using namespace gwutil;
#include <algorithm>
#include <set>
#include <cstdint>
#include <climits>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#ifdef _WIN32
  #include <windows.h>
  #include <direct.h>
#else
  #include <dirent.h>
#endif

namespace
{
  string Upper(string s){for (size_t i=0;i<s.size();i++){s[i]=(char)toupper((unsigned char)s[i]);} return s;}
  vector<string> Words(const string &line)
  {
    string l=line; size_t h=l.find('#'); if (h!=string::npos){l=l.substr(0,h);} h=l.find('!'); if (h!=string::npos){l=l.substr(0,h);}
    istringstream is(l); vector<string> w; string t; while (is>>t){w.push_back(t);} return w;
  }
  string DirOf(const string &f){size_t p=f.find_last_of("/\\"); return (p==string::npos)?string("./"):f.substr(0,p+1);}
  bool IsDir(const string &p){struct stat st; return (stat(p.c_str(),&st)==0) && ((st.st_mode&S_IFMT)==S_IFDIR);}
  void MakeDir(const string &p)
  {
#ifdef _WIN32
    _mkdir(p.c_str());
#else
    mkdir(p.c_str(),0755);
#endif
  }
  bool CopyFileBytes(const string &a,const string &b)
  {
    ifstream I(a.c_str(),ios::binary); if (!I){return false;}
    ofstream O(b.c_str(),ios::binary); if (!O){return false;}
    if (I.peek()!=EOF){O<<I.rdbuf();} //(an empty file stays empty: streaming nothing would set the fail bit)
    O.close(); return !O.fail();
  }
  //copies a simulation folder (recursively), leaving out old model results; returns the number of files copied
  int CopyTree(const string &src,const string &dst,const string &skipDir)
  {
    MakeDir(dst); int n=0;
    vector<string> names;
#ifdef _WIN32
    WIN32_FIND_DATAA fd; HANDLE h=FindFirstFileA((src+"*").c_str(),&fd);
    if (h!=INVALID_HANDLE_VALUE){do {names.push_back(fd.cFileName);} while (FindNextFileA(h,&fd)); FindClose(h);}
#else
    DIR *d=opendir(src.c_str()); if (d!=NULL){struct dirent *e; while ((e=readdir(d))!=NULL){names.push_back(e->d_name);} closedir(d);}
#endif
    for (size_t i=0;i<names.size();i++){
      const string &nm=names[i]; if ((nm==".") || (nm=="..")){continue;}
      string s=src+nm,t=dst+nm;
      if (IsDir(s)){ //never copy the working copy into itself (the output folder may lie inside the model folder)
        if (gwutil::InsideDir(skipDir,s)){continue;} //s is the working copy, or holds it
        n+=CopyTree(s+"/",t+"/",skipDir); continue;
      }
      string U=Upper(nm);
      if ((U.size()>4) && ((U.substr(U.size()-4)==".HDS") || (U.substr(U.size()-4)==".CBC") || (U.substr(U.size()-4)==".LST") || (U.substr(U.size()-4)==".GRB") ||
                          (U.substr(U.size()-4)==".BUD") || (U.substr(U.size()-4)==".CBB"))){continue;} //old results
      if (!CopyFileBytes(s,t)){ExitGracefully(("CGroundwaterModel: could not copy "+s+" to "+t).c_str(),BAD_DATA);}
      n++;
    }
    return n;
  }
  bool IsTakeOverType(const string &ftype){string f=Upper(ftype); return (f.find("RCH")==0) || (f.find("EVT")==0) || (f.find("UZF")==0);}
  bool IsStreamType(const string &ftype){string f=Upper(ftype); return (f=="RIV6") || (f=="DRN6") || (f=="GHB6");}
}

//////////////////////////////////////////////////////////////////
/// \brief reads mfsim.nam and the model name file; checks the declared packages (double-counting guard)
//
void CGroundwaterModel::ExtReadSimulation()
{
  ifstream S(_extSim.c_str());
  ExitGracefullyIf(!S,("CGroundwaterModel: cannot open MODFLOW 6 simulation file "+_extSim).c_str(),BAD_DATA);
  string line,block,tdis; vector<string> mnam,mname; int nsol=0; _extNModels=0;
  while (getline(S,line)){
    vector<string> w=Words(line); if (w.empty()){continue;}
    string W0=Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); continue;}
    if (W0=="END"){block=""; continue;}
    if ((block=="TIMING") && (W0=="TDIS6") && (w.size()>1)){tdis=w[1];}
    if ((block=="MODELS") && (w.size()>2)){_extNModels++;}
    if ((block=="MODELS") && (W0=="GWF6") && (w.size()>2)){mnam.push_back(w[1]); mname.push_back(w[2]);}
    if ((block=="SOLUTIONGROUP") && (W0!="MXITER") && (W0.size()>1) && (W0[W0.size()-1]=='6')){nsol++;} //IMS6, EMS6, ...
  }
  ExitGracefullyIf(tdis=="","CGroundwaterModel: no TDIS6 entry in the MODFLOW 6 simulation file",BAD_DATA);
  ExitGracefullyIf(mnam.empty(),"CGroundwaterModel: no GWF6 model in the MODFLOW 6 simulation file",BAD_DATA);
  //Raven drives one solution (SLN_1): a second solution (e.g. a transport model) would not be solved
  ExitGracefullyIf(nsol>1,"CGroundwaterModel: the MODFLOW 6 simulation has more than one solution (e.g. a GWT or GWE model with its own IMS); "
    "only simulations with a single groundwater-flow solution can be coupled at present",BAD_DATA);
  int im=-1;
  if (_extModel==""){
    ExitGracefullyIf(mnam.size()>1,"CGroundwaterModel: the simulation holds several GWF models; name the one to couple with :MF6Model",BAD_DATA);
    im=0;
  }
  else {for (size_t i=0;i<mname.size();i++){if (Upper(mname[i])==Upper(_extModel)){im=(int)i;}}}
  ExitGracefullyIf(im<0,("CGroundwaterModel: :MF6Model "+_extModel+" is not a GWF model of the simulation").c_str(),BAD_DATA);
  _extTdisFile=tdis; _extNamFile=mnam[im]; _mname=Upper(mname[im]);

  string dir=DirOf(_extSim);
  ifstream N((dir+_extNamFile).c_str());
  ExitGracefullyIf(!N,("CGroundwaterModel: cannot open the model name file "+dir+_extNamFile).c_str(),BAD_DATA);
  _extPkgs.clear(); block=""; map<string,int> count;
  while (getline(N,line)){
    vector<string> w=Words(line); if (w.empty()){continue;}
    string W0=Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); continue;}
    if (W0=="END"){block=""; continue;}
    if ((block=="PACKAGES") && (w.size()>=2)){
      gw_ext_pkg P; P.ftype=W0; P.fname=w[1];
      string base=W0.substr(0,W0.size()-1); //without the trailing 6
      if (w.size()>=3){P.pname=Upper(w[2]);}
      else {int c=++count[base]; P.pname=((base=="DIS") || (base=="DISV") || (base=="DISU") || (base=="NPF") || (base=="STO") || (base=="IC") || (base=="OC"))?base:base+"-"+to_string(c);}
      _extPkgs.push_back(P);
    }
  }
  _extDisFtype=""; _stoName="STO"; _disName="DIS";
  for (size_t i=0;i<_extPkgs.size();i++){
    const string &f=_extPkgs[i].ftype;
    if ((f=="DIS6") || (f=="DISV6") || (f=="DISU6")){_extDisFtype=f; _disName=_extPkgs[i].pname;}
    if (f=="STO6"){_stoName=_extPkgs[i].pname;}
  }
  ExitGracefullyIf(_extDisFtype=="","CGroundwaterModel: the MODFLOW 6 model has no discretization package",BAD_DATA);
  ExitGracefullyIf(_extDisFtype=="DISU6",
    "CGroundwaterModel: models with DISU discretization are not yet supported in existing-model mode (DIS and DISV are)",BAD_DATA);

  //every name the user gave must exist; every recharge/ET/UZF package must be declared taken over or kept
  vector<string> all=_takeOver; all.insert(all.end(),_keepPkgs.begin(),_keepPkgs.end()); if (_streamPkg!=""){all.push_back(_streamPkg);}
  for (size_t j=0;j<all.size();j++){
    bool found=false; for (size_t i=0;i<_extPkgs.size();i++){if (_extPkgs[i].pname==Upper(all[j])){found=true;}}
    ExitGracefullyIf(!found,("CGroundwaterModel: package "+all[j]+" is not in the model name file "+_extNamFile).c_str(),BAD_DATA);
  }
  for (size_t i=0;i<_extPkgs.size();i++){
    const gw_ext_pkg &P=_extPkgs[i];
    bool tk=(find(_takeOver.begin(),_takeOver.end(),P.pname)!=_takeOver.end()),kp=(find(_keepPkgs.begin(),_keepPkgs.end(),P.pname)!=_keepPkgs.end());
    ExitGracefullyIf(tk && kp,("CGroundwaterModel: package "+P.pname+" is listed both in :RavenTakesOver and :KeepPackages").c_str(),BAD_DATA);
    if (IsTakeOverType(P.ftype) && !tk && !kp){
      ExitGracefully(("CGroundwaterModel: the model's "+P.ftype+" package "+P.pname+" would add recharge or evapotranspiration beside Raven's. "
                      "List it in :RavenTakesOver (Raven replaces it) or :KeepPackages (it stays)").c_str(),BAD_DATA);
    }
    if (tk && !IsTakeOverType(P.ftype)){
      ExitGracefully(("CGroundwaterModel: :RavenTakesOver applies to recharge, evapotranspiration and UZF packages; "+P.pname+" is "+P.ftype).c_str(),BAD_DATA);
    }
    if ((P.pname==Upper(_streamPkg)) && !IsStreamType(P.ftype)){
      ExitGracefully(("CGroundwaterModel: :StreamPackage must be a RIV, DRN or GHB package; "+P.pname+" is "+P.ftype).c_str(),BAD_DATA);
    }
    if ((P.pname==Upper(_streamPkg)) && tk){ExitGracefully("CGroundwaterModel: the stream package cannot be taken over",BAD_DATA);}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief rewrites the copied TDIS: each stress period keeps its length and is split into Raven time steps;
///  the simulation ends with Raven's run
//
void CGroundwaterModel::ExtRewriteTDIS(const optStruct &Options)
{
  string f=_mf6dir+"/"+_extTdisFile;
  ifstream I(f.c_str()); ExitGracefullyIf(!I,("CGroundwaterModel: cannot open "+f).c_str(),BAD_DATA);
  string line,block,units="DAYS"; vector<double> perlen; vector<string> opts;
  while (getline(I,line)){
    vector<string> w=Words(line); if (w.empty()){continue;}
    string W0=Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); continue;}
    if (W0=="END"){block=""; continue;}
    if (block=="OPTIONS"){
      ExitGracefullyIf(W0=="ATS6","CGroundwaterModel: the MODFLOW 6 model uses adaptive time stepping (ATS6); Raven needs fixed time steps. Remove ATS6 from the TDIS options",BAD_DATA);
      if ((W0=="START_DATE_TIME") && (_extOffset>0)){opts.push_back("  START_DATE_TIME "+_extRavenStart+"  # Raven's start (restart)"); continue;}
      opts.push_back(line); if ((W0=="TIME_UNITS") && (w.size()>1)){units=Upper(w[1]);}
    }
    if (block=="PERIODDATA"){
      ExitGracefullyIf((W0=="OPEN/CLOSE") || (W0=="INTERNAL"),"CGroundwaterModel: TDIS period data must be listed in the TDIS file itself",BAD_DATA);
      perlen.push_back(atof(w[0].c_str()));
    }
  }
  I.close();
  double Tf=1.0; //model time units per day
  if      (units=="SECONDS"){Tf=86400.0;}
  else if (units=="MINUTES"){Tf=1440.0;}
  else if (units=="HOURS")  {Tf=24.0;}
  else if (units=="DAYS")   {Tf=1.0;}
  else if (units=="YEARS")  {Tf=1.0/365.25;}
  else {WriteWarning("CGroundwaterModel: MODFLOW 6 TIME_UNITS undefined; days assumed",Options.noisy);}
  _extTf=Tf;
  ExitGracefullyIf(perlen.empty(),"CGroundwaterModel: no stress periods in the TDIS file",BAD_DATA);
  //a steady-state model (no storage package, or storage steady in every period it names) with a single stress period has no
  //physical time: its period is stretched over Raven's run, so that every Raven step is a steady-state solution
  _extSteady=false;
  _extNoSto=true; for (size_t i=0;i<_extPkgs.size();i++){if (Upper(_extPkgs[i].ftype)=="STO6"){_extNoSto=false;}}
  if (perlen.size()==1){
    bool hasSto=false,trans=false,named=false;
    for (size_t i=0;i<_extPkgs.size();i++){
      if (Upper(_extPkgs[i].ftype)!="STO6"){continue;}
      hasSto=true;
      ifstream S((_mf6dir+"/"+_extPkgs[i].fname).c_str()); string l,b;
      while (getline(S,l)){
        vector<string> w=Words(l); if (w.empty()){continue;}
        string W0=Upper(w[0]);
        if ((W0=="BEGIN") && (w.size()>1)){b=Upper(w[1]); continue;}
        if (W0=="END"){b=""; continue;}
        if (b=="PERIOD"){if (W0=="TRANSIENT"){trans=true;} if (W0=="STEADY-STATE"){named=true;}}
      }
    }
    if (!hasSto || (named && !trans)){
      _extSteady=true;
      perlen[0]=(Options.duration+_extOffset)*Tf;
    }
  }
  double dt=Options.timestep,remain=Options.duration,total=0,skip=_extOffset;
  vector<double> pl; vector<int> ns; _extDrop=0;
  for (size_t p=0;p<perlen.size();p++){total+=perlen[p]/Tf;}
  ExitGracefullyIf(skip>=total-1e-9,"CGroundwaterModel: Raven starts after the end of the MODFLOW 6 model",BAD_DATA);
  total=0;
  for (size_t p=0;(p<perlen.size()) && (remain>1e-9);p++){
    double days=perlen[p]/Tf; total+=days;
    if (skip>0){ //elapsed before Raven's start (restart): drop whole periods, shorten the one holding the start
      if (skip>=days-1e-9){skip-=days; _extDrop++; continue;}
      double left=days-skip,k=left/dt;
      ExitGracefullyIf(fabs(k-floor(k+0.5))>1e-6,"CGroundwaterModel: Raven's start does not fall on a Raven time step within its stress period",BAD_DATA);
      days=left; skip=0;
    }
    int n=(int)floor(days/dt+0.5);
    ExitGracefullyIf((n<1) || (fabs(n*dt-days)>1e-6*max(1.0,days)),
      ("CGroundwaterModel: stress period "+to_string(p+1)+" ("+ToStr(days)+" d) is not a whole number of Raven time steps ("+ToStr(dt)+" d)").c_str(),BAD_DATA);
    if (days>remain+1e-9){
      n=(int)floor(remain/dt+0.5); days=n*dt;
      ExitGracefullyIf(n<1,"CGroundwaterModel: Raven's :Duration is not a whole number of time steps, so the last MODFLOW 6 stress period would be empty; make :Duration a multiple of :TimeStep",BAD_DATA);
    }
    pl.push_back(days*Tf); ns.push_back(n); remain-=days;
  }
  ExitGracefullyIf(remain>1e-6,("CGroundwaterModel: the MODFLOW 6 model ends before Raven's run ("+ToStr(total-_extOffset)+" days after Raven's start, :Duration "+ToStr(Options.duration)+")").c_str(),BAD_DATA);
  ofstream O(f.c_str());
  ExitGracefullyIf(!O,("CGroundwaterModel: cannot write "+f).c_str(),BAD_DATA);
  O<<"# rewritten by Raven: stress periods kept, split into Raven time steps, ending with Raven's run\nBEGIN OPTIONS\n";
  for (size_t i=0;i<opts.size();i++){O<<opts[i]<<"\n";}
  O<<"END OPTIONS\n\nBEGIN DIMENSIONS\n  NPER "<<pl.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIODDATA\n"<<setprecision(15);
  _nstp=0;
  for (size_t p=0;p<pl.size();p++){O<<"  "<<pl[p]<<" "<<ns[p]<<" 1.0\n"; _nstp+=ns[p];}
  O<<"END PERIODDATA\n";
}

//////////////////////////////////////////////////////////////////
/// \brief days since 1970-01-01 of a yyyy-mm-dd date (proleptic Gregorian); LLONG_MIN if malformed
//
long long CGroundwaterModel::CivilDays(const string &date)
{
  int y,m,d; char c1,c2; istringstream is(date.substr(0,10));
  if (!(is>>y>>c1>>m>>c2>>d) || (c1!='-') || (c2!='-') || (m<1) || (m>12) || (d<1) || (d>31)){return LLONG_MIN;}
  y-=(m<=2)?1:0; long long era=(y>=0?y:y-399)/400; long long yoe=y-era*400;
  long long doy=(153*(m+(m>2?-3:9))+2)/5+d-1,doe=yoe*365+yoe/4-yoe/100+doy;
  return era*146097+doe-719468;
}
//////////////////////////////////////////////////////////////////
/// \brief restart inside the model period: the first _extDrop stress periods are gone, so every kept package's
///  PERIOD blocks are renumbered; the new first period takes the block in effect at Raven's start
//
void CGroundwaterModel::ExtRenumberPeriods()
{
  //list-type packages: the last period block stays in effect until the next one
  const char *latest[12]={"WEL6","RIV6","GHB6","CHD6","DRN6","RCH6","EVT6","STO6","OC6","RCHA6","EVTA6","MVR6"};
  //advanced packages: a period block changes only the settings it names, the others carry on (like arrays)
  const char *accum[4]={"SFR6","LAK6","MAW6","UZF6"};
  //only the coupled model's packages are renumbered: other models' period data would stay on the old numbers
  ExitGracefullyIf((_extDrop>0) && (_extNModels>1),"CGroundwaterModel: Raven starts after the model's time zero, which renumbers the stress periods, "
    "but the simulation holds several models; start Raven at :MF6StartDate",BAD_DATA);
  for (size_t i=0;i<_extPkgs.size();i++){
    const gw_ext_pkg &P=_extPkgs[i];
    if (find(_takeOver.begin(),_takeOver.end(),P.pname)!=_takeOver.end()){continue;}
    string f=_mf6dir+"/"+P.fname;
    ifstream I(f.c_str()); if (!I){continue;}
    vector<string> pre; map<int,vector<string> > blk; string line; int k=-1; bool arrays=false,ts=false,tv=false;
    while (getline(I,line)){
      vector<string> w=Words(line); string W0=w.empty()?"":Upper(w[0]);
      if ((W0=="BEGIN") && (w.size()>2) && (Upper(w[1])=="PERIOD")){k=atoi(w[2].c_str()); blk[k].clear(); continue;}
      if ((k>=0) && (W0=="END")){k=-1; continue;}
      if (k>=0){blk[k].push_back(line); continue;}
      if (W0=="READASARRAYS"){arrays=true;}
      if ((W0=="TS6") || (W0=="TAS6") || (W0=="TIMESERIESFILE")){ts=true;}
      if ((W0=="TVK6") || (W0=="TVS6")){tv=true;} //time-varying K or storage: their own period blocks accumulate
      pre.push_back(line);
    }
    I.close();
    //time series count from the model's time zero: any later start shifts them, even if no period is dropped
    ExitGracefullyIf(ts,("CGroundwaterModel: "+P.pname+" uses time series, whose times count from the model's time zero; start Raven at :MF6StartDate").c_str(),BAD_DATA);
    if (_extDrop==0){continue;} //start inside the first period: period numbers are unchanged
    ExitGracefullyIf(tv,("CGroundwaterModel: "+P.pname+" uses time-varying properties (TVK6/TVS6), whose period data cannot be renumbered safely; start Raven at :MF6StartDate").c_str(),BAD_DATA);
    if (blk.empty()){continue;}
    bool ok=false; for (int j=0;j<12;j++){if (P.ftype==latest[j]){ok=true;}}
    if ((P.ftype=="RCHA6") || (P.ftype=="EVTA6")){arrays=true;} //array-based by type
    for (int j=0;j<4;j++){if (P.ftype==accum[j]){ok=true; arrays=true;}} //settings accumulate: as for arrays
    if (arrays){ //array blocks accumulate (arrays not given keep earlier values): safe only when one block precedes the start
      int nbefore=0; for (map<int,vector<string> >::iterator it=blk.begin();it!=blk.end();it++){if (it->first<=_extDrop+1){nbefore++;}}
      ok=ok && (nbefore<=1); arrays=!ok;
    }
    ExitGracefullyIf(!ok || arrays,("CGroundwaterModel: Raven starts after the model's time zero, which needs the stress periods renumbered; "
      "the period data of "+P.pname+" ("+P.ftype+(arrays?"; its period blocks carry earlier settings forward":"")+") cannot be renumbered safely. Start Raven at :MF6StartDate").c_str(),BAD_DATA);
    map<int,vector<string> > nb;
    int first=-1; for (map<int,vector<string> >::iterator it=blk.begin();it!=blk.end();it++){if (it->first<=_extDrop+1){first=it->first;}} //in effect at the start
    if (first>0){nb[1]=blk[first];}
    for (map<int,vector<string> >::iterator it=blk.begin();it!=blk.end();it++){if (it->first>_extDrop+1){nb[it->first-_extDrop]=it->second;}}
    ofstream O(f.c_str());
    ExitGracefullyIf(!O,("CGroundwaterModel: cannot write "+f).c_str(),BAD_DATA);
    for (size_t j=0;j<pre.size();j++){O<<pre[j]<<"\n";}
    for (map<int,vector<string> >::iterator it=nb.begin();it!=nb.end();it++){
      O<<"BEGIN PERIOD "<<it->first<<(((it->first==1) && (first!=_extDrop+1))?"   # in effect at Raven's start (was period "+to_string(first)+")":string(""))<<"\n";
      for (size_t j=0;j<it->second.size();j++){O<<it->second[j]<<"\n";}
      O<<"END PERIOD\n";
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief the model's own specified-head cells (IBOUND<0) at the top of a column take no recharge in MODFLOW: Raven's
///  recharge is spread over each HRU's other cells, as for generated models (updated when the stress period changes them)
//
void CGroundwaterModel::ExtUpdateCHDTop()
{
  string a=_mf6.GetVarAddress(_mname,"","IBOUND"); if (!_mf6.HasVar(a)){return;}
  int *ib=_mf6.GetIntPtr(a); if (ib==NULL){return;}
  bool changed=false;
  for (int ic=0;ic<_ncpl;ic++){
    int n=_topNode[ic]; if (n<0){continue;}
    int r=_reduced[n]; int c=((r>=0) && (ib[r]<0))?1:0;
    if (c!=_colCHDTop[ic]){_colCHDTop[ic]=c; changed=true;}
  }
  if (!changed){return;}
  _hruRchArea.assign(_pModel->GetNumHRUs(),0.0);
  for (size_t li=0;li<_aLinks.size();li++){if (!_colCHDTop[_aLinks[li].ic]){_hruRchArea[_aLinks[li].k]+=_aLinks[li].area;}}
}

//////////////////////////////////////////////////////////////////
/// \brief writes the working copy's simulation file as mfsim.nam (the name MODFLOW's library opens), with CONTINUE
///  set so that a step that does not converge is flagged by the coupling instead of stopping the whole run
//
void CGroundwaterModel::ExtWriteSimFile()
{
  string base=_extSim.substr(_extSim.find_last_of("/\\")==string::npos?0:_extSim.find_last_of("/\\")+1);
  string src=_mf6dir+"/"+base;
  ifstream I(src.c_str()); ExitGracefullyIf(!I,("CGroundwaterModel: cannot read the copied simulation file "+src).c_str(),BAD_DATA);
  vector<string> L; string line; bool inOpt=false,hasOpt=false,hasCont=false;
  while (getline(I,line)){
    if (!line.empty() && (line[line.size()-1]=='\r')){line.erase(line.size()-1);}
    vector<string> w=Words(line); string W0=w.empty()?"":Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1) && (Upper(w[1])=="OPTIONS")){inOpt=true; hasOpt=true;}
    else if (inOpt && (W0=="END")){inOpt=false;}
    else if (inOpt && (W0=="CONTINUE")){hasCont=true;}
    L.push_back(line);
  }
  I.close();
  ofstream O((_mf6dir+"/mfsim.nam").c_str());
  ExitGracefullyIf(!O,("CGroundwaterModel: cannot write "+_mf6dir+"/mfsim.nam").c_str(),BAD_DATA);
  if (!hasOpt){O<<"# CONTINUE added by Raven: a step that does not converge is flagged in GWBudget.csv\nBEGIN OPTIONS\n  CONTINUE\nEND OPTIONS\n\n";}
  for (size_t i=0;i<L.size();i++){
    O<<L[i]<<"\n";
    vector<string> w=Words(L[i]);
    if (hasOpt && !hasCont && (w.size()>1) && (Upper(w[0])=="BEGIN") && (Upper(w[1])=="OPTIONS")){O<<"  CONTINUE  # added by Raven\n";}
  }
  ExitGracefullyIf(!O.good(),("CGroundwaterModel: failed writing "+_mf6dir+"/mfsim.nam").c_str(),BAD_DATA);
}

//////////////////////////////////////////////////////////////////
/// \brief the model's water mover (MVR): movers from or to a package Raven takes over are removed from the copied MVR file
///  (MODFLOW stops if the MVR names a package that is no longer in the model); a stream package that sends its water
///  through the mover is refused, because Raven would take the same water from the aquifer a second time
//
void CGroundwaterModel::ExtEditMover()
{
  string fn;
  for (size_t i=0;i<_extPkgs.size();i++){if (Upper(_extPkgs[i].ftype).find("MVR")==0){fn=_extPkgs[i].fname;}}
  if (fn==""){return;}
  string f=_mf6dir+"/"+fn;
  ifstream I(f.c_str()); ExitGracefullyIf(!I,("CGroundwaterModel: cannot read the copied water mover file "+f).c_str(),BAD_DATA);
  vector<string> L; string line,block;
  while (getline(I,line)){L.push_back(line);}
  I.close();
  set<string> tk(_takeOver.begin(),_takeOver.end());
  vector<string> keepPk; bool anyTk=false;
  for (size_t i=0;i<L.size();i++){ //the packages the mover connects
    vector<string> w=Words(L[i]); if (w.empty() || (w[0][0]=='#')){continue;}
    string W0=Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); continue;}
    if (W0=="END"){block=""; continue;}
    if (block=="PACKAGES"){
      string pk=Upper(w.back()); //(a model-level MVR lists package names; a simulation-level one would add the model name)
      if ((_streamPkg!="") && (pk==Upper(_streamPkg))){
        ExitGracefully(("CGroundwaterModel: the stream package "+_streamPkg+" sends water through the model's water mover ("+fn+
          "); Raven would take that water from the aquifer a second time. Choose a stream package that is not a mover provider").c_str(),BAD_DATA);
      }
      if (tk.count(pk)){anyTk=true;} else {keepPk.push_back(pk);}
    }
  }
  if (!anyTk){return;}
  ExitGracefullyIf(keepPk.size()<2,("CGroundwaterModel: with "+string("the packages Raven takes over removed, the water mover ")+fn+
    " connects fewer than two packages; remove the MVR package from the model name file").c_str(),BAD_DATA);
  ofstream O(f.c_str()); ExitGracefullyIf(!O,("CGroundwaterModel: cannot write "+f).c_str(),BAD_DATA);
  int ndrop=0; block="";
  for (size_t i=0;i<L.size();i++){
    vector<string> w=Words(L[i]);
    if (w.empty() || (w[0][0]=='#')){O<<L[i]<<"\n"; continue;}
    string W0=Upper(w[0]);
    if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); O<<L[i]<<"\n"; continue;}
    if (W0=="END"){block=""; O<<L[i]<<"\n"; continue;}
    if ((block=="DIMENSIONS") && (W0=="MAXPACKAGES")){O<<"  MAXPACKAGES  "<<keepPk.size()<<"   # "<<L[i]<<" (packages taken over by Raven removed)\n"; continue;}
    if ((block=="PACKAGES") && tk.count(Upper(w.back()))){O<<"# "<<L[i]<<"   # taken over by Raven\n"; continue;}
    if (block=="PERIOD"){
      ExitGracefullyIf((W0=="OPEN/CLOSE") || (W0=="INTERNAL"),("CGroundwaterModel: the water mover "+fn+
        " reads its period data from another file; Raven cannot remove the movers of the packages it takes over").c_str(),BAD_DATA);
      //provider id receiver id type value (model-level): a mover from or to a taken-over package goes
      if ((w.size()>=4) && (tk.count(Upper(w[0])) || tk.count(Upper(w[2])))){ndrop++; continue;}
    }
    O<<L[i]<<"\n";
  }
  ExitGracefullyIf(!O.good(),("CGroundwaterModel: failed writing "+f).c_str(),BAD_DATA);
  _extMoverDropped=ndrop;
}

//////////////////////////////////////////////////////////////////
/// \brief edits the copied model name file: taken-over packages removed, Raven's packages (and initial heads) added
//
void CGroundwaterModel::ExtEditNameFile(const bool addRaven)
{
  string f=_mf6dir+"/"+_extNamFile;
  ifstream I(f.c_str()); vector<string> L; string line,block;
  ExitGracefullyIf(!I,("CGroundwaterModel: cannot read the copied model name file "+f).c_str(),BAD_DATA);
  while (getline(I,line)){L.push_back(line);}
  I.close();
  ExitGracefullyIf(L.empty(),("CGroundwaterModel: the copied model name file "+f+" is empty").c_str(),BAD_DATA);
  ofstream O(f.c_str());
  ExitGracefullyIf(!O,("CGroundwaterModel: cannot write "+f).c_str(),BAD_DATA);
  for (size_t i=0;i<L.size();i++){
    vector<string> w=Words(L[i]);
    if (!w.empty() && (Upper(w[0])=="BEGIN") && (w.size()>1)){block=Upper(w[1]);}
    if (!w.empty() && (Upper(w[0])=="END") && (block=="PACKAGES") && addRaven){
      O<<"  RCH6  raven.rch  RCH_RAVEN   # added by Raven\n";
      if (_extSeepage){O<<"  DRN6  raven.drn  DRN_RAVEN   # added by Raven\n";}
    }
    if (!w.empty() && (Upper(w[0])=="END")){block="";}
    if ((block=="PACKAGES") && (w.size()>=2) && (Upper(w[0])!="BEGIN")){
      bool tk=false; for (size_t j=0;j<_extPkgs.size();j++){if ((_extPkgs[j].fname==w[1]) && (find(_takeOver.begin(),_takeOver.end(),_extPkgs[j].pname)!=_takeOver.end())){tk=true;}}
      if (tk){O<<"# "<<L[i]<<"   # taken over by Raven\n"; continue;}
      if (addRaven && (_hotHeads.size()>0) && (Upper(w[0])=="IC6")){O<<"  IC6  raven.ic  "<<((w.size()>=3)?w[2]:string("IC"))<<"   # hotstart heads from Raven\n"; continue;}
    }
    O<<L[i]<<"\n";
  }
}

//////////////////////////////////////////////////////////////////
/// \brief reads the model's geometry, active cells and units from MODFLOW 6 memory (probe initialisation)
//
void CGroundwaterModel::ExtProbe()
{
  string lib=MF6LibraryPath();
  if (!_mf6.Load(lib))        {ExitGracefully(("CGroundwaterModel: "+_mf6.GetLastError()+"."+MF6LibraryHelp()).c_str(),RUNTIME_ERR);}
  _mf6Version=_mf6.GetVersion();
  { //MODFLOW 6.6 (and older) stalls when a model uses NEWTON with its UNDER_RELAXATION keyword and a cell dries: the same Newton
    //change repeats every iteration and the step never converges (6.8.1 converges); warn, as the model is the user's
    int ma=0,mi=0; bool uo=false;
    if (sscanf(_mf6Version.c_str(),"%d.%d",&ma,&mi)==2){uo=(ma<6) || ((ma==6) && (mi<8));}
    ifstream N((_mf6dir+"/"+_extNamFile).c_str()); string l,b; bool nu=false;
    while (getline(N,l)){
      vector<string> w=Words(l); if (w.empty()){continue;}
      string W0=Upper(w[0]);
      if ((W0=="BEGIN") && (w.size()>1)){b=Upper(w[1]); continue;}
      if (W0=="END"){b=""; continue;}
      if ((b=="OPTIONS") && (W0=="NEWTON")){for (size_t j=1;j<w.size();j++){if (Upper(w[j])=="UNDER_RELAXATION"){nu=true;}}}
    }
    if (uo && nu){
      WriteWarning("CGroundwaterModel: the model uses NEWTON UNDER_RELAXATION with MODFLOW "+_mf6Version+"; with versions before 6.8 "
        "a drying cell can make every later step fail to converge. Use MODFLOW 6.8 or newer, or remove UNDER_RELAXATION from the NEWTON option",false);
    }
  }
  if (!_mf6.Initialize(_mf6dir)){ExitGracefully(("CGroundwaterModel: the MODFLOW 6 model does not initialise: "+_mf6.GetLastError()+" (see "+_mf6dir+"/mfsim.lst)").c_str(),RUNTIME_ERR);}
  string M=_mname,D=_disName;
  int *p=_mf6.GetIntPtr(_mf6.GetVarAddress(M,D,"NLAY")); ExitGracefullyIf(p==NULL,"CGroundwaterModel: cannot read NLAY of the MODFLOW 6 model",RUNTIME_ERR);
  if (p==NULL){return;}
  _nlay=*p;
  if (_extDisFtype=="DIS6"){
    int *pr=_mf6.GetIntPtr(_mf6.GetVarAddress(M,D,"NROW")),*pc=_mf6.GetIntPtr(_mf6.GetVarAddress(M,D,"NCOL"));
    ExitGracefullyIf((pr==NULL) || (pc==NULL),"CGroundwaterModel: cannot read NROW/NCOL of the MODFLOW 6 model",RUNTIME_ERR);
    _nrow=*pr; _ncol=*pc; _ncpl=_nrow*_ncol; _disType=0;
  }
  else {
    int *pn=_mf6.GetIntPtr(_mf6.GetVarAddress(M,D,"NCPL"));
    ExitGracefullyIf(pn==NULL,"CGroundwaterModel: cannot read NCPL of the MODFLOW 6 model",RUNTIME_ERR);
    _ncpl=*pn; _nrow=1; _ncol=_ncpl; _disType=1;
  }
  int nuser=_nlay*_ncpl;
  string ar=_mf6.GetVarAddress(M,D,"NODEREDUCED");
  int *nr=(_mf6.HasVar(ar) && (_mf6.GetVarSize(ar)==nuser))?_mf6.GetIntPtr(ar):NULL;
  _reduced.assign(nuser,-1);
  int nodes=0;
  for (int n=0;n<nuser;n++){_reduced[n]=(nr==NULL)?n:((nr[n]>0)?nr[n]-1:-1); if (_reduced[n]>=0){nodes++;}}
  double *top=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"TOP")),*bot=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"BOT")),*area=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"AREA"));
  ExitGracefullyIf((top==NULL) || (bot==NULL) || (area==NULL),"CGroundwaterModel: cannot read TOP, BOT or AREA of the MODFLOW 6 model",RUNTIME_ERR);
  int *lu=_mf6.GetIntPtr(_mf6.GetVarAddress(M,D,"LENUNI"));
  int lenuni=(lu==NULL)?0:*lu;
  double Lf=1.0; //model length units per metre
  if      (lenuni==1){Lf=1.0/0.3048;}
  else if (lenuni==3){Lf=100.0;}
  else if (lenuni==0){WriteWarning("CGroundwaterModel: MODFLOW 6 LENGTH_UNITS undefined; metres assumed",false);}
  _extLf=Lf;
  _hf=1.0/Lf; _qf=_extTf/(Lf*Lf*Lf); _rf=Lf/_extTf;
  //per column: top active cell, its top and area; cell bottoms (m)
  _topNode.assign(_ncpl,-1); _top.assign(_ncpl,0.0); _extArea.assign(_ncpl,0.0); _extTopRaw.assign(_ncpl,0.0);
  _botm.assign(nuser,0.0); _idomain.assign(nuser,0); _strt.assign(nuser,0.0); _nodeTop.assign(nuser,0.0);
  for (int ic=0;ic<_ncpl;ic++){
    for (int l=0;l<_nlay;l++){
      int n=UserNode(l,ic),r=_reduced[n]; if (r<0){continue;}
      _idomain[n]=1; _botm[n]=bot[r]*_hf; _nodeTop[n]=top[r]*_hf;
      if (_topNode[ic]<0){_topNode[ic]=n; _top[ic]=top[r]*_hf; _extTopRaw[ic]=top[r]; _extArea[ic]=area[r]*_hf*_hf;}
    }
  }
  _userOf.assign(nodes,-1); for (int n=0;n<nuser;n++){if (_reduced[n]>=0){_userOf[_reduced[n]]=n;}}
  string dv="SLN_1/DVCLOSE"; if (_mf6.HasVar(dv)){double *d=_mf6.GetDoublePtr(dv); if (d!=NULL){_extDvclose=*d*_hf;}}
  if (_extCRS.IsSet()){ExtCellsFromMemory();}
  _extConn.clear(); _extLocalC.clear();
  if (_extGridFile!=""){ //cell centres in the model's frame (metres), to check the layout of the grid file later
    double X0,Y0,ca,sa; vector<gw_ring> local; ExtLocalCells(local,X0,Y0,ca,sa);
    _extLocalC.assign(_ncpl,gw_pt());
    for (int ic=0;ic<_ncpl;ic++){
      double sx=0,sy=0; for (size_t k=0;k<local[ic].size();k++){sx+=local[ic][k].x; sy+=local[ic][k].y;}
      double n=max((double)local[ic].size(),1.0),x=sx/n,y=sy/n; //vertex mean: enough for the check
      _extLocalC[ic].x=(X0+x*ca-y*sa)*_hf; _extLocalC[ic].y=(Y0+x*sa+y*ca)*_hf;
    }
  }
  if (!_owTable.empty() || (_extGridFile!="")){ //neighbours of the top active cells from MODFLOW's own connections
    int *IA=_mf6.GetIntPtr(_mf6.GetVarAddress(M,"CON","IA")),*JA=_mf6.GetIntPtr(_mf6.GetVarAddress(M,"CON","JA"));
    int *IHC=_mf6.GetIntPtr(_mf6.GetVarAddress(M,"CON","IHC")),*JAS=_mf6.GetIntPtr(_mf6.GetVarAddress(M,"CON","JAS"));
    if (!_owTable.empty()){_grid.nbr.assign(_ncpl,vector<gw_nbr>());}
    if ((IA!=NULL) && (JA!=NULL) && (IHC!=NULL) && (JAS!=NULL)){
      for (int ic=0;ic<_ncpl;ic++){
        if (_topNode[ic]<0){continue;}
        int r=_reduced[_topNode[ic]];
        for (int q=IA[r]-1;q<IA[r+1]-1;q++){
          int m=JA[q]-1; if ((m==r) || (IHC[JAS[q]-1]==0)){continue;}
          int jc=_userOf[m]%_ncpl; if (jc==ic){continue;}
          if (jc>ic){_extConn.push_back(make_pair(ic,jc));}
          if (_owTable.empty()){continue;}
          gw_nbr e; e.j=jc; e.len=e.di=e.dj=e.ang=0.0;
          bool dup=false; for (size_t k=0;k<_grid.nbr[ic].size();k++){if (_grid.nbr[ic][k].j==e.j){dup=true;}}
          if (!dup){_grid.nbr[ic].push_back(e);}
        }
      }
    }
  }
  _mf6.Finalize();
}

//////////////////////////////////////////////////////////////////
/// \brief top-layer cell outlines from MODFLOW memory (DIS: DELR/DELC; DISV: VERTICES/IAVERT/JAVERT), placed with
///  XORIGIN/YORIGIN/ANGROT, converted to metres and to longitude/latitude with :MF6CRS
//
void CGroundwaterModel::ExtCellsFromMemory()
{
  double X0,Y0,ca,sa; vector<gw_ring> local; ExtLocalCells(local,X0,Y0,ca,sa);
  //model coordinates -> projected metres -> lon/lat; fingerprint for the linkage cache
  _extCellsLL.assign(_ncpl,gw_ring()); uint64_t hsh=1469598103934665603ULL;
  for (int ic=0;ic<_ncpl;ic++){
    for (size_t k=0;k<local[ic].size();k++){
      double xm=(X0+local[ic][k].x*ca-local[ic][k].y*sa)*_hf,ym=(Y0+local[ic][k].x*sa+local[ic][k].y*ca)*_hf; gw_pt q;
      _extCRS.Inverse(xm,ym,q.x,q.y); _extCellsLL[ic].push_back(q);
      unsigned char *b=(unsigned char*)&xm; for (int j=0;j<8;j++){hsh^=b[j]; hsh*=1099511628211ULL;}
      b=(unsigned char*)&ym;                for (int j=0;j<8;j++){hsh^=b[j]; hsh*=1099511628211ULL;}
    }
  }
  ostringstream k; k<<"crs:"<<_extCRS.Describe()<<":"<<hex<<hsh; _extGeomKey=k.str();
}
//////////////////////////////////////////////////////////////////
/// \brief checks that the polygons of :MF6GridFile are the model's cells, in the right order.
/// \details A map projection scales areas and turns directions smoothly, so the checks allow one overall scale and one
///  overall rotation but nothing else: (1) the ratio of each cell's area to MODFLOW's is the same for all cells (a wrong
///  length unit or projection moves the typical ratio far from 1); (2) for every pair of cells MODFLOW connects, the
///  polygons are as far apart, and in the same direction, as MODFLOW's cells (a shuffled numbering changes distances;
///  a mirrored one reverses directions). Some bending is allowed: 20% in distance, 15 degrees in direction.
//
void CGroundwaterModel::ExtCheckGridFile()
{
  const string f=_extGridFile;
  vector<double> ratio;
  for (int ic=0;ic<_ncpl;ic++){if ((_topNode[ic]>=0) && (_extArea[ic]>0)){ratio.push_back(_grid.area[ic]/_extArea[ic]);}}
  ExitGracefullyIf(ratio.empty(),"CGroundwaterModel: the MODFLOW 6 model has no active top cells",BAD_DATA);
  vector<double> sr=ratio; sort(sr.begin(),sr.end()); double med=sr[sr.size()/2];
  ExitGracefullyIf((med<0.8) || (med>1.25),("CGroundwaterModel: the cells of "+f+" are typically "+ToStr(med)+" times MODFLOW's cell area. "
    "Check the length units of the model and the projection in which the polygons were made").c_str(),BAD_DATA);
  int nbadA=0,wc=-1; double worst=0;
  for (int ic=0;ic<_ncpl;ic++){
    if ((_topNode[ic]<0) || !(_extArea[ic]>0)){continue;}
    double d=fabs(_grid.area[ic]/_extArea[ic]/med-1.0);
    if (d>0.05){nbadA++;} if (d>worst){worst=d; wc=ic;}
  }
  //connections: distance ratio and direction difference against MODFLOW's own cell centres
  int nbadC=0,nC=0,wcC=-1; double sc=sqrt(med),mc=0,ms=0;
  vector<double> dang(_extConn.size(),0.0),drat(_extConn.size(),1.0);
  for (size_t e=0;e<_extConn.size();e++){
    int i=_extConn[e].first,j=_extConn[e].second;
    double fx=_grid.cx[j]-_grid.cx[i],fy=_grid.cy[j]-_grid.cy[i],mx=_extLocalC[j].x-_extLocalC[i].x,my=_extLocalC[j].y-_extLocalC[i].y;
    double dm=hypot(mx,my); if (!(dm>0)){continue;}
    drat[e]=hypot(fx,fy)/(dm*sc);
    double a=atan2(fy,fx)-atan2(my,mx); dang[e]=a; mc+=cos(a); ms+=sin(a); nC++;
  }
  double amean=atan2(ms,mc);
  for (size_t e=0;e<_extConn.size();e++){
    double a=dang[e]-amean; a=atan2(sin(a),cos(a));
    if ((drat[e]<0.8) || (drat[e]>1.25) || (fabs(a)>15.0*3.14159265358979323846/180.0)){nbadC++; if (wcC<0){wcC=_extConn[e].first;}}
  }
  ExitGracefullyIf(nbadC>0,("CGroundwaterModel: in "+f+", "+to_string(nbadC)+" of "+to_string(nC)+" pairs of neighbouring cells do not lie "
    "like MODFLOW's cells (first at cell "+to_string(wcC+1)+"). The cell numbering (or its orientation: MODFLOW's row 1 is the "
    "northern row) or the cell shapes differ from the model's").c_str(),BAD_DATA);
  ExitGracefullyIf(nbadA>0,("CGroundwaterModel: "+to_string(nbadA)+" cells of "+f+" differ in area from MODFLOW's by more than 5% beyond "
    "the common scale (worst: cell "+to_string(wc+1)+", "+ToStr(100*worst)+"%). Check the cell numbering and shapes").c_str(),BAD_DATA);
  _extAreaCheck=med;
}

//////////////////////////////////////////////////////////////////
/// \brief the top-layer cell outlines in the model's own coordinates (model length units, before XORIGIN/ANGROT),
///  read from MODFLOW's memory: DIS from DELR/DELC (row 1 north), DISV from its vertices
//
void CGroundwaterModel::ExtLocalCells(vector<gw_ring> &local,double &X0,double &Y0,double &ca,double &sa)
{
  string M=_mname,D=_disName;
  double *xo=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"XORIGIN")),*yo=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"YORIGIN")),*an=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"ANGROT"));
  X0=(xo!=NULL)?*xo:0.0; Y0=(yo!=NULL)?*yo:0.0; double ang=((an!=NULL)?*an:0.0)*3.14159265358979323846/180.0; ca=cos(ang); sa=sin(ang);
  local.assign(_ncpl,gw_ring());
  if (_disType==0){
    double *dr=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"DELR")),*dc=_mf6.GetDoublePtr(_mf6.GetVarAddress(M,D,"DELC"));
    ExitGracefullyIf((dr==NULL) || (dc==NULL),"CGroundwaterModel: cannot read DELR/DELC of the MODFLOW 6 model",RUNTIME_ERR);
    vector<double> xe(_ncol+1,0.0),ye(_nrow+1,0.0); //column edges from the west; row edges from the south (row 1 is north)
    for (int c=0;c<_ncol;c++){xe[c+1]=xe[c]+dr[c];}
    for (int r=_nrow-1;r>=0;r--){ye[_nrow-r]=ye[_nrow-1-r]+dc[r];}
    for (int r=0;r<_nrow;r++){for (int c=0;c<_ncol;c++){
      double ya=ye[_nrow-1-r],yb=ye[_nrow-r]; gw_ring &R=local[r*_ncol+c]; gw_pt p;
      p.x=xe[c];   p.y=ya; R.push_back(p); p.x=xe[c+1]; p.y=ya; R.push_back(p);
      p.x=xe[c+1]; p.y=yb; R.push_back(p); p.x=xe[c];   p.y=yb; R.push_back(p);
    }}
  }
  else {
    string av=_mf6.GetVarAddress(M,D,"VERTICES"),ai=_mf6.GetVarAddress(M,D,"IAVERT"),aj=_mf6.GetVarAddress(M,D,"JAVERT");
    double *V=_mf6.GetDoublePtr(av); int *IA=_mf6.GetIntPtr(ai),*JA=_mf6.GetIntPtr(aj);
    ExitGracefullyIf((V==NULL) || (IA==NULL) || (JA==NULL),"CGroundwaterModel: cannot read the cell vertices of the MODFLOW 6 model",RUNTIME_ERR);
    for (int ic=0;ic<_ncpl;ic++){
      for (int k=IA[ic]-1;k<IA[ic+1]-1;k++){gw_pt p; int v=JA[k]-1; p.x=V[2*v]; p.y=V[2*v+1]; local[ic].push_back(p);}
      if ((local[ic].size()>1) && (local[ic].front().x==local[ic].back().x) && (local[ic].front().y==local[ic].back().y)){local[ic].pop_back();}
    }
  }
}

//////////////////////////////////////////////////////////////////
/// \brief links HRUs to the top-layer cells: overlap areas from cell polygons (checked against the model's own cell
///  areas) or from an :OverlapWeights table; columns without an active cell are left out
//
void CGroundwaterModel::ExtLinkHRUs(const optStruct &Options)
{
  int nHRUs=_pModel->GetNumHRUs();
  _hruCells.assign(nHRUs,map<int,double>()); _hruPolyArea.assign(nHRUs,0.0);
  vector<double> gridArea;
  if (_extCRS.IsSet()){ //cells from MODFLOW memory: geometry is the model's own, so no area check is needed
    int ncpl=_ncpl,nlay=_nlay,nrow=_nrow,ncol=_ncol;
    _gridType=GWGRID_FILE; _gridFile=""; BuildGeometry(Options);
    _ncpl=ncpl; _nlay=nlay; _nrow=nrow; _ncol=ncol;
  }
  else if (_extGridFile!=""){
    map<long long,vector<gw_ring> > shapes; string err;
    if (!GWGeom::ReadGeoJSONPolygons(_extGridFile,_extGridID,shapes,err)){ExitGracefully(("CGroundwaterModel: "+err).c_str(),BAD_DATA);}
    ExitGracefullyIf((int)shapes.size()!=_ncpl,("CGroundwaterModel: "+_extGridFile+" holds "+to_string(shapes.size())+" cells; the model layer has "+to_string(_ncpl)).c_str(),BAD_DATA);
    long long expect=1;
    for (map<long long,vector<gw_ring> >::iterator it=shapes.begin();it!=shapes.end();it++,expect++){
      ExitGracefullyIf(it->first!=expect,("CGroundwaterModel: cell IDs in "+_extGridFile+" must be 1.."+to_string(_ncpl)+" (cell number within a layer)").c_str(),BAD_DATA);
    }
    int ncpl=_ncpl,nlay=_nlay,nrow=_nrow,ncol=_ncol; //BuildGeometry sets the grid of the polygons
    _gridType=GWGRID_FILE; _gridFile=_extGridFile; _gridFileID=_extGridID;
    BuildGeometry(Options);
    _ncpl=ncpl; _nlay=nlay; _nrow=nrow; _ncol=ncol;
    ExtCheckGridFile();
  }
  else {
    //:OverlapWeights: w = A_overlap / A_cell for the model's top layer (MODFLOW-USG coupling format)
    _grid.ncpl=_ncpl; _grid.area=_extArea;
    for (size_t e=0;e<_owTable.size();e++){
      CHydroUnit *pH=_pModel->GetHRUByID(_owTable[e].hru);
      ExitGracefullyIf(pH==NULL,("CGroundwaterModel: :OverlapWeights refers to unknown HRU "+to_string(_owTable[e].hru)).c_str(),BAD_DATA);
      int ic=_owTable[e].cell-1;
      ExitGracefullyIf((ic<0) || (ic>=_ncpl),("CGroundwaterModel: :OverlapWeights cell "+to_string(_owTable[e].cell)+" is outside the model layer (1.."+to_string(_ncpl)+")").c_str(),BAD_DATA);
      int k=pH->GetGlobalIndex();
      if (_hruProfile[k]==DOESNT_EXIST){continue;}
      double A=_owTable[e].w*((_topNode[ic]>=0)?_extArea[ic]:0.0);
      if (A>0){_hruCells[k][ic]+=A; _hruPolyArea[k]+=A;}
    }
    for (size_t i=0;i<_aqHRUs.size();i++){_hruPolyArea[_aqHRUs[i]]=max(_hruPolyArea[_aqHRUs[i]],1e-300);}
  }
  //columns: active in the model and covered enough by coupled HRUs; dominant HRU and its profile
  vector<double> cover(_ncpl,0.0),best(_ncpl,0.0);
  _colProfile.assign(_ncpl,DOESNT_EXIST); _colDomHRU.assign(_ncpl,DOESNT_EXIST); _land.assign(_ncpl,0.0);
  for (size_t i=0;i<_aqHRUs.size();i++){
    int k=_aqHRUs[i];
    for (map<int,double>::iterator it=_hruCells[k].begin();it!=_hruCells[k].end();it++){
      int ic=it->first; if ((ic<0) || (ic>=_ncpl) || (it->second<=0)){continue;}
      cover[ic]+=it->second;
      if (it->second>best[ic]){best[ic]=it->second; _colProfile[ic]=_hruProfile[k]; _colDomHRU[ic]=k;}
    }
  }
  for (int ic=0;ic<_ncpl;ic++){
    double Acell=(_owTable.empty())?_grid.area[ic]:_extArea[ic]; //same area measure as the overlaps
    if ((_topNode[ic]<0) || (cover[ic]<_minCoverage*Acell)){_colProfile[ic]=DOESNT_EXIST;}
    _land[ic]=_top[ic];
  }
}

//////////////////////////////////////////////////////////////////
/// \brief builds the coupling to an existing model (replaces grid generation and file writing)
//
void CGroundwaterModel::ExtBuild(const optStruct &Options)
{
  ExtReadSimulation();
  GW_MKDIR(Options.output_dir.c_str()); GW_MKDIR((Options.output_dir+"mf6").c_str());
  _mf6dir=AbsPath(Options.output_dir+"mf6");
  ExitGracefullyIf(gwutil::InsideDir(DirOf(_extSim),_mf6dir),("CGroundwaterModel: the MODFLOW 6 simulation "+_extSim+
    " lies in Raven's working folder "+_mf6dir+", which Raven overwrites; keep the model in a folder of its own").c_str(),BAD_DATA);
  int nf=CopyTree(DirOf(_extSim),_mf6dir+"/",_mf6dir+"/");
  ExitGracefullyIf(nf==0,"CGroundwaterModel: nothing was copied from the MODFLOW 6 simulation folder",BAD_DATA);
  if (nf>200){WriteAdvisory("CGroundwaterModel: "+to_string(nf)+" files were copied from the folder of "+_extSim+
                            "; keep the MODFLOW 6 model in a folder of its own so that only its files are copied",Options.noisy);}
  //start date: the model's time zero must be Raven's start
  time_struct t0; JulianConvert(0.0,Options.julian_start_day,Options.julian_start_year,Options.calendar,t0);
  ExitGracefullyIf((Options.calendar!=CALENDAR_PROLEPTIC_GREGORIAN) && (Options.calendar!=CALENDAR_GREGORIAN),
    "CGroundwaterModel: an existing MODFLOW 6 model needs a Gregorian :Calendar (dates are compared day by day)",BAD_DATA);
  long long dm=CivilDays(_extStart),dr=CivilDays(t0.date_string);
  _extRavenStart=t0.date_string.substr(0,10);
  ExitGracefullyIf((dm==LLONG_MIN) || (dr==LLONG_MIN),"CGroundwaterModel: :MF6StartDate must be yyyy-mm-dd",BAD_DATA);
  ExitGracefullyIf(dr<dm,("CGroundwaterModel: Raven's :StartDate "+t0.date_string+" is before the model's time zero (:MF6StartDate "+_extStart+")").c_str(),BAD_DATA);
  //[d] Raven starts this long after the model's time zero (a restart), including a time of day in :StartDate
  _extOffset=(double)(dr-dm)+(Options.julian_start_day-floor(Options.julian_start_day));
  ExitGracefullyIf((_extOffset>0) && (Options.timestep<1.0-1e-9),"CGroundwaterModel: runs starting after the model's time zero need daily or longer time steps",BAD_DATA);
  ExtWriteSimFile();
  ExtRewriteTDIS(Options);
  ExtEditMover();
  if (_extOffset>0){ExtRenumberPeriods();} //(also checks what a later start makes unsafe when no period is dropped)
  ExtEditNameFile(false);
  ExtProbe();
  if ((_hotHeads.size()>0) && ((_hotDims[0]!=_nlay) || (_hotDims[1]!=_nrow) || (_hotDims[2]!=_ncol) || ((int)_hotHeads.size()!=_nlay*_ncpl))){
    ExitGracefully("CGroundwaterModel: :GWHeads in the .rvc file does not match the MODFLOW 6 model",BAD_DATA);
  }
  ExtLinkHRUs(Options);
  BuildLinks();
  //stream package: auxiliary index of the subbasin ID
  _extAux=-1; _extNaux=0;
  if ((_streamPkg!="") && (_streamMap=="AUX")){
    string fn; for (size_t i=0;i<_extPkgs.size();i++){if (_extPkgs[i].pname==Upper(_streamPkg)){fn=_extPkgs[i].fname;}}
    ifstream P((_mf6dir+"/"+fn).c_str()); string line,block;
    while (getline(P,line)){
      vector<string> w=Words(line); if (w.empty()){continue;}
      string W0=Upper(w[0]);
      if ((W0=="BEGIN") && (w.size()>1)){block=Upper(w[1]); continue;}
      if (W0=="END"){block=""; continue;}
      if ((block=="OPTIONS") && ((W0=="AUXILIARY") || (W0=="AUX"))){
        _extNaux=(int)w.size()-1;
        for (size_t j=1;j<w.size();j++){if (Upper(w[j])==Upper(_streamAux)){_extAux=(int)j-1;}}
      }
    }
    ExitGracefullyIf(_extAux<0,("CGroundwaterModel: stream package "+_streamPkg+" has no auxiliary variable "+_streamAux).c_str(),BAD_DATA);
  }
  _sbGain.assign(_pModel->GetNumSubBasins(),0.0);
  _sbLossRemain.assign(_pModel->GetNumSubBasins(),0.0);
  _sbLossCarry.assign(_pModel->GetNumSubBasins(),0.0);
  for (int p=0;p<_pModel->GetNumSubBasins();p++){
    map<long long,double>::iterator it=_hotRivCarry.find(_pModel->GetSubBasin(p)->GetID());
    if (it!=_hotRivCarry.end()){_sbLossCarry[p]=it->second;}
  }
  _sbRivLen.assign(_pModel->GetNumSubBasins(),0.0);
  if (_owTable.empty()){BuildStresses();} //monitoring wells (located through the cell polygons)
}

//////////////////////////////////////////////////////////////////
/// \brief writes Raven's packages (in the model's units) and, on hotstart, the initial heads
//
void CGroundwaterModel::ExtWritePackages(const optStruct &Options)
{
  (void)Options;
  ofstream F((_mf6dir+"/raven.rch").c_str());
  F<<"# recharge from Raven: rates set every time step\nBEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  MAXBOUND "<<_bndCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n";
  for (size_t i=0;i<_bndCol.size();i++){int ic=_bndCol[i],n=TopNode(ic); F<<"  "<<CellID(n/_ncpl,ic)<<" 0.0\n";}
  F<<"END PERIOD\n"; F.close();
  if (_extSeepage){
    F.open((_mf6dir+"/raven.drn").c_str());
    F<<"# seepage to Raven at the top of the top active cell (model units)\nBEGIN OPTIONS\n  AUXILIARY DDRN\n  AUXDEPTHNAME DDRN\n  SAVE_FLOWS\nEND OPTIONS\n\n";
    F<<"BEGIN DIMENSIONS\n  MAXBOUND "<<_bndCol.size()<<"\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n"<<setprecision(12);
    for (size_t i=0;i<_bndCol.size();i++){
      int ic=_bndCol[i],n=TopNode(ic); const gw_profile &P=_aProfiles[_colProfile[ic]];
      double cond=P.seep_leakance*_extArea[ic]*_extLf*_extLf/_extTf; //[L2/T]
      double d=min(P.seep_smooth_depth,_top[ic]-_botm[n]); //MODFLOW refuses a smoothing interval below the cell bottom
      F<<"  "<<CellID(n/_ncpl,ic)<<" "<<_extTopRaw[ic]<<" "<<cond<<" "<<-d*_extLf<<"\n";
    }
    F<<"END PERIOD\n"; F.close();
  }
  if (_hotHeads.size()>0){
    F.open((_mf6dir+"/raven.ic").c_str());
    F<<"# initial heads from Raven's hotstart (:GWHeads)\nBEGIN GRIDDATA\n";
    WriteArray(F,"STRT",_hotHeads,(_disType==0)?_ncol:10,_nlay);
    F<<"END GRIDDATA\n"; F.close();
  }
  ExtEditNameFile(true);
}

//////////////////////////////////////////////////////////////////
/// \brief flows of the stream package to Raven's reaches, and of all other kept packages (for the balance)
//
void CGroundwaterModel::ExtStreamFlows()
{
  _Qriv_in=_Qriv_out=0.0; _Qother=0.0;
  if (_streamPkg!=""){
    string P=Upper(_streamPkg);
    int *nb=_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,P,"NBOUND"));
    int *nl=_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,P,"NODELIST"));
    double *hc=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,P,"HCOF")),*rh=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,P,"RHS"));
    double *ax=NULL;
    if (_extAux>=0){ //MODFLOW 6.5+ keeps the values read from the file in its input context; older versions in the package
      string a=_mf6.GetVarAddress(_mname,P,"AUXVAR"),ia="__INPUT__/"+a,in="__INPUT__/"+_mf6.GetVarAddress(_mname,P,"NBOUND");
      int *inb=_mf6.HasVar(in)?_mf6.GetIntPtr(in):NULL;
      if (_mf6.HasVar(ia) && (inb!=NULL) && (nb!=NULL) && (*inb==*nb)){ax=_mf6.GetDoublePtr(ia);}
      else {ax=_mf6.GetDoublePtr(a);}
    }
    ExitGracefullyIf((nb==NULL) || (nl==NULL) || (hc==NULL) || (rh==NULL) || ((_extAux>=0) && (ax==NULL)),
      "CGroundwaterModel: cannot read the stream package from MODFLOW 6 memory",RUNTIME_ERR);
    for (int i=0;i<*nb;i++){
      int r=nl[i]-1; if (r<0){continue;}
      double q=(hc[i]*_X[r]-rh[i])*_qf; //into the aquifer [m3/d]
      int p=DOESNT_EXIST;
      if (_extAux>=0){
        long long id=llround(ax[(size_t)i*_extNaux+_extAux]);
        map<long long,int>::iterator it=_extSBIndex.find(id);
        if (it==_extSBIndex.end()){p=_pModel->GetSubBasinIndex(id); _extSBIndex[id]=p;} else {p=it->second;}
        ExitGracefullyIf(p<0,("CGroundwaterModel: stream package entry with subbasin ID "+to_string(id)+" - no such subbasin").c_str(),BAD_DATA);
      }
      else {
        int ic=_userOf[r]%_ncpl,k=_colDomHRU[ic];
        ExitGracefullyIf(k==DOESNT_EXIST,"CGroundwaterModel: a stream cell lies on a column without coupled HRUs; use AUX subbasin IDs",BAD_DATA);
        p=_pModel->GetHydroUnit(k)->GetSubBasinIndex();
      }
      ExitGracefullyIf(!_pModel->GetSubBasin(p)->IsEnabled(),"CGroundwaterModel: a stream cell exchanges with a disabled subbasin",BAD_DATA);
      _sbGain[p]-=q;
      if (q>0){_Qriv_in+=q;} else {_Qriv_out-=q;}
    }
  }
  for (size_t j=0;j<_extOther.size();j++){
    int *nb=_mf6.GetIntPtr(_mf6.GetVarAddress(_mname,_extOther[j],"NBOUND"));
    double *sv=_mf6.GetDoublePtr(_mf6.GetVarAddress(_mname,_extOther[j],"SIMVALS"));
    if ((nb==NULL) || (sv==NULL)){continue;}
    for (int i=0;i<*nb;i++){_Qother+=sv[i]*_qf;}
    //a package that sends water to the model's mover (MVR) reports that part separately (the list file's "-TO-MVR" term):
    //it also leaves the aquifer
    //(only when the package's mover is on: otherwise the array is not sized to the package)
    string im=_mf6.GetVarAddress(_mname,_extOther[j],"IMOVER"),mv=_mf6.GetVarAddress(_mname,_extOther[j],"SIMTOMVR");
    int *pim=_mf6.HasVar(im)?_mf6.GetIntPtr(im):NULL;
    if ((pim!=NULL) && (*pim!=0) && _mf6.HasVar(mv) && (_mf6.GetVarSize(mv)>=*nb)){double *sm=_mf6.GetDoublePtr(mv); if (sm!=NULL){for (int i=0;i<*nb;i++){_Qother+=sm[i]*_qf;}}}
  }
}

//////////////////////////////////////////////////////////////////
/// \brief after the engine is connected: packages whose simulated flows enter the balance check
//
void CGroundwaterModel::ExtConnect()
{
  _extOther.clear();
  for (size_t i=0;i<_extPkgs.size();i++){
    const gw_ext_pkg &P=_extPkgs[i];
    if ((find(_takeOver.begin(),_takeOver.end(),P.pname)!=_takeOver.end()) || (P.pname==Upper(_streamPkg))){continue;}
    string a=_mf6.GetVarAddress(_mname,P.pname,"SIMVALS");
    if (_mf6.HasVar(a)){_extOther.push_back(P.pname);}
  }
}
