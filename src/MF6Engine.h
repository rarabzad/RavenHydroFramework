/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: thin wrapper around the
  MODFLOW 6 shared library (libmf6) XMI/BMI interface, loaded at runtime
----------------------------------------------------------------*/
#ifndef MF6ENGINE_H
#define MF6ENGINE_H
#include <string>
#include <vector>
using namespace std;

class CMF6Engine
{
private:
  void  *_lib;
  string _workdir;      ///< absolute path of MODFLOW 6 simulation directory
  string _lasterr;
  bool   _initialized;

  typedef int (*f_str)   (const char*);
  typedef int (*f_void)  ();
  typedef int (*f_dbl)   (double*);
  typedef int (*f_int)   (int*);
  typedef int (*f_int2)  (int*,int*);
  typedef int (*f_ptr)   (const char*,void**);
  typedef int (*f_addr)  (const char*,const char*,const char*,char*);
  typedef int (*f_rank)  (const char*,int*);
  typedef int (*f_errmsg)(char*);

  f_str  _initialize; f_void _finalize;
  f_dbl  _prepare_time_step; f_void _finalize_time_step;
  f_int  _prepare_solve; f_int2 _solve; f_int _finalize_solve;
  f_ptr  _get_value_ptr; f_addr _get_var_address;
  f_rank _get_var_rank;  f_rank _get_var_shape;
  f_errmsg _get_last_bmi_error;
  f_dbl  _get_current_time;
  f_errmsg _get_version;   ///< optional (MODFLOW 6 exports get_version)

  void *Sym(const char *name);
  bool  Check(int rc,const string &what);

public:
  CMF6Engine();
  ~CMF6Engine();
  bool   Load            (const string &libpath);
  bool   Initialize      (const string &workdir);
  bool   PrepareTimeStep (double dt);
  bool   SolveTimeStep   (int maxiter,bool &converged,int &niter,
                          double *dvclose=NULL,double relaxedTol=0.0,int relaxedAfter=0,bool *relaxed=NULL);
  bool   FinalizeTimeStep();
  void   Finalize        ();

  string GetVarAddress   (const string &comp,const string &sub,const string &var);
  bool   HasVar          (const string &address);
  int    GetVarSize      (const string &address);
  double*GetDoublePtr    (const string &address);
  int   *GetIntPtr       (const string &address);
  double GetModelTime  ();
  string GetLastError    () const {return _lasterr;}
  string GetVersion      ();
  static string         DefaultLibraryName(); ///< libmf6.dll / libmf6.dylib / libmf6.so
  static string         ExecutableDir     (); ///< folder of the running executable
  static vector<string> LibraryCandidates (); ///< automatic search places for the library
  static bool           FileExists        (const string &path);       ///< the library's MODFLOW 6 version ("" if the library does not say)
  bool   IsInitialized   () const {return _initialized;}
};
#endif
