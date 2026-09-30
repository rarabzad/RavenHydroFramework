#include "GWGeometry.h"
#include <cstdio>
#include <cstdlib>
int main(int argc,char**argv){ CGWCRS c; std::string err; if(!c.Parse(argv[1],err)){printf("ERR %s\n",err.c_str());return 1;}
  double x,y,lon,lat; while(scanf("%lf %lf",&x,&y)==2){c.Inverse(x,y,lon,lat); printf("%.12f %.12f\n",lon,lat);} }
