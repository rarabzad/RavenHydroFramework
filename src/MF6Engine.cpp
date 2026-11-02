#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: libmf6 wrapper
----------------------------------------------------------------*/
#include "MF6Engine.h"
#include "RavenInclude.h"   // WriteWarning, ExitGracefully (version guard)
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <climits>
#ifdef _WIN32
  #include <windows.h>
  #include <direct.h>
  #define GW_CHDIR _chdir
  #define GW_GETCWD _getcwd
#else
  #include <dlfcn.h>
  #include <unistd.h>
  #include <sys/stat.h>
  #ifdef __APPLE__
    #include <mach-o/dyld.h>
  #endif
  #define GW_CHDIR chdir
  #define GW_GETCWD getcwd
#endif

namespace
{
  /// Scoped change of working directory: MODFLOW 6 resolves all of its file names
  /// relative to the directory holding mfsim.nam, while Raven uses its own run directory
  //-- exit handler: MODFLOW stops the whole program (a Fortran STOP, exit code 2) on an input error it finds while Raven calls
  //   it; the handler, registered once, reports MODFLOW's error report as Raven's error (only while a MODFLOW call is in
  //   progress), from the run folder, and lets the exit go on
  const char *g_call=NULL;          ///< the MODFLOW call in progress (NULL: none)
  char        g_runDir[4096]={0};   ///< Raven's folder before the call
  string      g_mf6Dir;             ///< MODFLOW's folder (mfsim.lst)
  bool        g_registered=false;
  void MF6ExitHandler()
  {
    if (g_call==NULL){return;}
    const string call(g_call); g_call=NULL;
    string rep="",line; bool on=false;
    ifstream L((g_mf6Dir+"/mfsim.lst").c_str());
    while (getline(L,line)){
      if (line.find("ERROR REPORT")!=string::npos){on=true; continue;}
      if (!on){continue;}
      size_t a=line.find_first_not_of(" \t\r"); if (a==string::npos){continue;}
      size_t b=line.find_last_not_of(" \t\r"); line=line.substr(a,b-a+1);
      if ((line.find("UNIT NUMBER")!=string::npos) || (line.find("STOP")==0)){break;}
      rep+=(rep.empty()?"":" ")+line;
      if (rep.size()>1500){break;}
    }
    if (rep.empty()){rep="see "+g_mf6Dir+"/mfsim.lst";}
    if (g_runDir[0]!='\0'){if (GW_CHDIR(g_runDir)!=0){}}
    const string msg="CGroundwaterModel: MODFLOW 6 stopped the run during "+call+": "+rep;
    ExitGracefully(msg.c_str(),BAD_DATA_WARN);                                 //(writes the error; no exit from within exit)
    cerr<<endl<<"============== Exiting Gracefully ==========================="<<endl<<"Exiting Gracefully: "<<msg<<endl
        <<"Error Type: Runtime Error"<<endl<<"============================================================="<<endl;
  }
  struct CallGuard //marks a MODFLOW call in progress
  {
    CallGuard(const char *name,const string &dir){
      if (!g_registered){g_registered=true; atexit(MF6ExitHandler);}
      g_call=name; g_mf6Dir=dir;
    }
    ~CallGuard(){g_call=NULL;}
  };
  class CDGuard
  {
    char _old[4096]; bool _ok;
  public:
    CDGuard(const string &dir){
      _ok=(GW_GETCWD(_old,sizeof(_old))!=NULL);
      if (_ok){strncpy(g_runDir,_old,sizeof(g_runDir)-1);}
      if (!dir.empty()){ if (GW_CHDIR(dir.c_str())!=0){_ok=false;} }
    }
    ~CDGuard(){ if (_ok){ if (GW_CHDIR(_old)!=0){} } }
    bool ok() const {return _ok;}
  };
}

CMF6Engine::CMF6Engine():_lib(NULL),_initialized(false),
  _initialize(NULL),_finalize(NULL),_prepare_time_step(NULL),_finalize_time_step(NULL),
  _prepare_solve(NULL),_solve(NULL),_finalize_solve(NULL),_get_value_ptr(NULL),_get_var_address(NULL),
  _get_var_rank(NULL),_get_var_shape(NULL),_get_last_bmi_error(NULL),_get_current_time(NULL),_get_version(NULL){}
CMF6Engine::~CMF6Engine(){Finalize();}

void *CMF6Engine::Sym(const char *name)
{
#ifdef _WIN32
  void *f=(void*)GetProcAddress((HMODULE)_lib,name);
#else
  void *f=dlsym(_lib,name);
#endif
  if (f==NULL){_lasterr+=string(" missing symbol ")+name+";";}
  return f;
}

//////////////////////////////////////////////////////////////////
/// \brief file name of the MODFLOW 6 library on this platform
//
string CMF6Engine::DefaultLibraryName()
{
#if defined(_WIN32)
  return "libmf6.dll";
#elif defined(__APPLE__)
  return "libmf6.dylib";
#else
  return "libmf6.so";
#endif
}
//////////////////////////////////////////////////////////////////
/// \brief folder holding the running Raven executable ("" if it cannot be determined)
//
string CMF6Engine::ExecutableDir()
{
  string p;
#if defined(_WIN32)
  char buf[4096]; DWORD n=GetModuleFileNameA(NULL,buf,sizeof(buf));
  if ((n>0) && (n<sizeof(buf))){p=string(buf,n);}
#elif defined(__APPLE__)
  char buf[4096]; uint32_t sz=sizeof(buf);
  if (_NSGetExecutablePath(buf,&sz)==0){char r[PATH_MAX]; p=(realpath(buf,r)!=NULL)?string(r):string(buf);}
#else
  char buf[4096]; ssize_t n=readlink("/proc/self/exe",buf,sizeof(buf)-1);
  if (n>0){p=string(buf,(size_t)n);}
#endif
  size_t k=p.find_last_of("/\\");
  return (k==string::npos)?string(""):p.substr(0,k);
}
//////////////////////////////////////////////////////////////////
/// \brief places searched for the library when neither :MF6Library nor RAVEN_MF6_LIB is set:
///  next to the executable, then in lib/mf6 beside it or one to three folders up (src/Raven.exe, build/Raven,
///  build/Release/Raven.exe and Visual Studio's src/x64/Release/Raven.exe of a source checkout all reach the checkout's
///  lib/mf6)
//
vector<string> CMF6Engine::LibraryCandidates()
{
  vector<string> c; string d=ExecutableDir(),n=DefaultLibraryName();
  if (d==""){return c;}
  const char *sub[5]={"","lib/mf6/","../lib/mf6/","../../lib/mf6/","../../../lib/mf6/"};
  for (int i=0;i<5;i++){c.push_back(d+"/"+sub[i]+n);}
  return c;
}
bool CMF6Engine::FileExists(const string &path)
{
#if defined(_WIN32)
  DWORD a=GetFileAttributesA(path.c_str()); return (a!=INVALID_FILE_ATTRIBUTES) && !(a & FILE_ATTRIBUTE_DIRECTORY);
#else
  struct stat s; return (stat(path.c_str(),&s)==0) && S_ISREG(s.st_mode);
#endif
}

bool CMF6Engine::Load(const string &libpath)
{
  _lasterr="";
#ifdef _WIN32
  _lib=(void*)LoadLibraryA(libpath.c_str());
  if (_lib==NULL){_lasterr="cannot load MODFLOW 6 library "+libpath;return false;}
#else
  _lib=dlopen(libpath.c_str(),RTLD_NOW|RTLD_LOCAL);
  if (_lib==NULL){_lasterr="cannot load MODFLOW 6 library "+libpath+": "+string(dlerror());return false;}
#endif
  _initialize        =(f_str)   Sym("initialize");
  _finalize          =(f_void)  Sym("finalize");
  _prepare_time_step =(f_dbl)   Sym("prepare_time_step");
  _finalize_time_step=(f_void)  Sym("finalize_time_step");
  _prepare_solve     =(f_int)   Sym("prepare_solve");
  _solve             =(f_int2)  Sym("solve");
  _finalize_solve    =(f_int)   Sym("finalize_solve");
  _get_value_ptr     =(f_ptr)   Sym("get_value_ptr");
  _get_var_address   =(f_addr)  Sym("get_var_address");
  _get_var_rank      =(f_rank)  Sym("get_var_rank");
  _get_var_shape     =(f_rank)  Sym("get_var_shape");
  _get_last_bmi_error=(f_errmsg)Sym("get_last_bmi_error");
  _get_current_time  =(f_dbl)   Sym("get_current_time");
  bool ok=(_lasterr=="");
#ifdef _WIN32
  _get_version=(f_errmsg)GetProcAddress((HMODULE)_lib,"get_version");
#else
  _get_version=(f_errmsg)dlsym(_lib,"get_version");
#endif
  if (ok){CheckVersion(libpath);}
  return ok;
}

//////////////////////////////////////////////////////////////////
/// \brief numeric parts of a MODFLOW version string: "6.8.1" -> {6,8,1}. Leading blanks and a sign are read as strtol
/// reads them; reading stops at the first character that does not continue the number, so " 6.8.1", "6.8.1-dev",
/// "6.8.1+abc.dirty" and "6.08.01" read as 6.8.1, "6.5.0rc1" as 6.5.0, "6..9" as 6 and "-1.0" as -1.0
/// \return false if the string does not start with a number ("", "abc", "v6.8.1")
//
bool CMF6Engine::ParseVersion(const string &v,vector<long> &c)
{
  c.clear();
  const char *p=v.c_str(); char *e=NULL;
  long x=strtol(p,&e,10);
  if (e==p){return false;}
  c.push_back(x); p=e;
  while ((*p=='.') && isdigit((unsigned char)(p[1]))){
    x=strtol(p+1,&e,10); c.push_back(x); p=e;
  }
  return true;
}
/// \brief compares versions part by part as numbers (6.10.0 > 6.8.1; missing parts count as 0, so 6.5 == 6.5.0)
int CMF6Engine::CompareVersions(const vector<long> &a,const vector<long> &b)
{
  for (size_t i=0;i<max(a.size(),b.size());i++){
    long x=(i<a.size())?a[i]:0, y=(i<b.size())?b[i]:0;
    if (x<y){return -1;}
    if (x>y){return  1;}
  }
  return 0;
}

//////////////////////////////////////////////////////////////////
/// \brief version guard, right after the library is loaded: a library older than MF6_MIN_VERSION stops the run; any
/// version other than MF6_TESTED_VERSION, no version, or a version that cannot be read gives one warning per run
//
void CMF6Engine::CheckVersion(const string &libpath)
{
  static bool warned=false;
  string v=GetVersion();
  vector<long> have,lo,tested;
  ParseVersion(MF6_MIN_VERSION,lo); ParseVersion(MF6_TESTED_VERSION,tested);
  string need=string("the coupling needs MODFLOW 6 ")+MF6_MIN_VERSION+" or newer and was tested with "+MF6_TESTED_VERSION;
  string warn="";
  if (_get_version==NULL){
    warn="the MODFLOW 6 library "+libpath+" does not report its version (no get_version); "+need;
  }
  else if (v==""){
    warn="the MODFLOW 6 library "+libpath+" does not report its version; "+need;
  }
  else if (!ParseVersion(v,have)){
    warn="the MODFLOW 6 library "+libpath+" reports a version that cannot be read (\""+v+"\"); "+need;
  }
  else if (CompareVersions(have,lo)<0){
    ExitGracefully(("CGroundwaterModel: the MODFLOW 6 library "+libpath+" is version "+v+"; the coupling needs MODFLOW 6 "+
                    MF6_MIN_VERSION+" or newer (tested with "+MF6_TESTED_VERSION+"). Get the tested version with "
                    "tools/get_mf6.py of RavenHydroFramework-extra, or build Raven with -DRAVEN_MF6=PREBUILT or SOURCE, "
                    "which puts it next to the Raven executable, and remove :MF6Library or RAVEN_MF6_LIB if they name "
                    "the old library").c_str(),BAD_DATA);
    return;
  }
  else if (CompareVersions(have,tested)!=0){
    warn="the MODFLOW 6 library "+libpath+" is version "+v+"; the coupling was tested with MODFLOW 6 "+MF6_TESTED_VERSION+
         " (another version can give slightly different results; tools/get_mf6.py of RavenHydroFramework-extra gets "
         "the tested one)";
  }
  if ((warn!="") && (!warned)){
    WriteWarning("CGroundwaterModel: "+warn,true);
    warned=true;
  }
}

string CMF6Engine::GetVersion()
{
  if (_get_version==NULL){return "";}
  char buf[512]; memset(buf,0,sizeof(buf));
  if (_get_version(buf)!=0){return "";}
  string v(buf); size_t e=v.find_last_not_of(" \t\r\n"); return (e==string::npos)?"":v.substr(0,e+1);
}

bool CMF6Engine::Check(int rc,const string &what)
{
  if (rc==0){return true;}
  char buf[4096]; buf[0]='\0';
  if (_get_last_bmi_error!=NULL){_get_last_bmi_error(buf);}
  _lasterr="MODFLOW 6 "+what+" failed: "+string(buf);
  return false;
}

bool CMF6Engine::Initialize(const string &workdir)
{
  _workdir=workdir;
  CDGuard g(_workdir);
  if (!g.ok()){_lasterr="cannot enter MODFLOW 6 folder "+_workdir; return false;}
  CallGuard cg("initialize",_workdir);
  if (!Check(_initialize("mfsim.nam"),"initialize")){return false;}
  _initialized=true;
  return true;
}
bool CMF6Engine::PrepareTimeStep(double dt)
{
  CDGuard g(_workdir);
  if (!g.ok()){_lasterr="cannot enter MODFLOW 6 folder "+_workdir; return false;}
  CallGuard cg("prepare_time_step",_workdir);
  return Check(_prepare_time_step(&dt),"prepare_time_step");
}
/// Newton loop: the step ends when MODFLOW has converged (OUTER_DVCLOSE) and, with afterIter, the coupling iteration
/// too. Optionally (dvclose = MODFLOW's tolerance variable), a step on which MODFLOW has not converged in relaxedAfter
/// outer iterations (counting only iterations MODFLOW did not converge on, not those spent on the coupling) continues
/// with the tolerance raised to relaxedTol for the rest of that step only (*relaxed = true if it then converges); the
/// tolerance is restored afterwards.
bool CMF6Engine::SolveTimeStep(int maxiter,bool &converged,int &niter,
                               double *dvclose,double relaxedTol,int relaxedAfter,bool *relaxed,
                               const std::function<bool(int,bool)> &afterIter,
                               const std::function<void(int)> &beforeSolve)
{
  CDGuard g(_workdir);
  CallGuard cg("solve",_workdir);
  if (!g.ok()){_lasterr="cannot enter MODFLOW 6 folder "+_workdir; return false;}
  int sol=1,conv=0,nhard=0; bool isrelaxed=false;
  converged=false; niter=0; if (relaxed!=NULL){*relaxed=false;}
  bool tiered=(dvclose!=NULL) && (relaxedTol>0) && (relaxedAfter>0);
  double target=tiered?*dvclose:0.0;
  if (!Check(_prepare_solve(&sol),"prepare_solve")){return false;}
  for (niter=1;niter<=maxiter;niter++){
    if (tiered && !isrelaxed && (nhard==relaxedAfter)){*dvclose=relaxedTol; isrelaxed=true;}
    if (beforeSolve){beforeSolve(niter);}                                       //(may change MODFLOW's iterate: KeepColumnsWet)
    if (!Check(_solve(&sol,&conv),"solve")){if (tiered){*dvclose=target;} return false;}
    if (conv!=1){nhard++;}
    bool cpl=true;
    if (afterIter){cpl=afterIter(niter,conv==1);}
    if ((conv==1) && cpl){converged=true; if (isrelaxed && (relaxed!=NULL)){*relaxed=true;} break;}
    if (niter==maxiter){converged=(conv==1);} //coupling not converged at the limit: report MODFLOW's own state
  }
  if (tiered){*dvclose=target;}
  if (niter>maxiter){niter=maxiter;}
  return Check(_finalize_solve(&sol),"finalize_solve");
}
bool CMF6Engine::FinalizeTimeStep()
{
  CDGuard g(_workdir);
  if (!g.ok()){_lasterr="cannot enter MODFLOW 6 folder "+_workdir; return false;}
  CallGuard cg("finalize_time_step",_workdir);
  return Check(_finalize_time_step(),"finalize_time_step");
}
void CMF6Engine::Finalize()
{
  if (_initialized && (_finalize!=NULL)){
    CDGuard g(_workdir);
    if (g.ok()){_finalize();}
  }
  _initialized=false;
}
string CMF6Engine::GetVarAddress(const string &comp,const string &sub,const string &var)
{
  char buf[1024]; buf[0]='\0';
  if (_get_var_address(comp.c_str(),sub.c_str(),var.c_str(),buf)!=0){return "";}
  return string(buf);
}
bool CMF6Engine::HasVar(const string &address)
{
  int rank=-1;
  if (address.empty()){return false;}
  return (_get_var_rank(address.c_str(),&rank)==0);
}
int CMF6Engine::GetVarSize(const string &address)
{
  int rank=0;
  if (_get_var_rank(address.c_str(),&rank)!=0){return -1;}
  if (rank==0){return 1;}
  int shape[8]={0,0,0,0,0,0,0,0};
  if (_get_var_shape(address.c_str(),shape)!=0){return -1;}
  int n=1; for (int i=0;i<rank;i++){n*=shape[i];}
  return n;
}
double *CMF6Engine::GetDoublePtr(const string &address)
{
  void *p=NULL;
  if (!Check(_get_value_ptr(address.c_str(),&p),"get_value_ptr("+address+")")){return NULL;}
  return (double*)p;
}
int *CMF6Engine::GetIntPtr(const string &address)
{
  void *p=NULL;
  if (!Check(_get_value_ptr(address.c_str(),&p),"get_value_ptr("+address+")")){return NULL;}
  return (int*)p;
}
double CMF6Engine::GetModelTime()
{
  double t=0; if (_get_current_time!=NULL){_get_current_time(&t);} return t;
}
