/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: small helpers shared by the coupling's source files
----------------------------------------------------------------*/
#ifndef GWUTIL_H
#define GWUTIL_H
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#ifdef _WIN32
  #include <direct.h>
  #define GW_MKDIR(d) _mkdir(d)
#else
  #define GW_MKDIR(d) mkdir(d,0755)
#endif
namespace gwutil
{
  inline std::string AbsPath(const std::string &p)
  {
#ifdef _WIN32
    char buf[4096]; if (_fullpath(buf,p.c_str(),sizeof(buf))!=NULL){return std::string(buf);} return p;
#else
    char buf[4096]; if (realpath(p.c_str(),buf)!=NULL){return std::string(buf);} return p;
#endif
  }
  /// absolute path in one comparable form: forward slashes, no trailing slash (Windows: lower case)
  inline std::string NormPath(const std::string &p)
  {
    std::string a=AbsPath(p);
    for (size_t i=0;i<a.size();i++){if (a[i]=='\\'){a[i]='/';}}
    while ((a.size()>1) && (a[a.size()-1]=='/')){a.erase(a.size()-1);}
#ifdef _WIN32
    for (size_t i=0;i<a.size();i++){a[i]=(char)tolower((unsigned char)a[i]);}
#endif
    return a;
  }
  /// true if folder a is folder b or lies inside it
  inline bool InsideDir(const std::string &a,const std::string &b)
  {
    std::string x=NormPath(a),y=NormPath(b);
    return (x==y) || ((x.size()>y.size()) && (x.compare(0,y.size(),y)==0) && (x[y.size()]=='/'));
  }
  /// writes a MODFLOW 6 array; 3-D arrays (nlay>1) in LAYERED form. MODFLOW reads one grid row per
  /// record, so every row must start on a new line (wrapped at 10 values within the row)
  template <class T> void WriteArray(std::ofstream &OUT,const std::string &name,const std::vector<T> &v,const int ncol,const int nlay=1,const char *fac="1.0")
  {
    size_t n=v.size()/nlay;
    OUT<<"  "<<name<<((nlay>1)?"  LAYERED":"")<<std::endl<<std::setprecision(10);
    for (int l=0;l<nlay;l++){
      OUT<<"    INTERNAL FACTOR "<<fac<<std::endl;
      for (size_t i=0;i<n;i++){
        size_t j=i%ncol;
        OUT<<((j%10==0)?"      ":" ")<<v[l*n+i];
        if ((j%10==9) || (j==(size_t)ncol-1) || (i==n-1)){OUT<<std::endl;} //every row, and the last value, end a line
      }
    }
  }
  inline std::string ToStr(double v){std::ostringstream s; s<<v; return s.str();}
}
#endif
