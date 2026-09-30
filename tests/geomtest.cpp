// Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
// SPDX-License-Identifier: Artistic-2.0
// Unit test harness for GWGeometry: clipping areas, projection, raster sampling
#include "GWGeometry.h"   // compile: g++ -std=c++11 -O2 -I../src geomtest.cpp ../src/GWGeometry.cpp -o geomtest; then python3 geomtest_check.py
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(){
  srand(12345);
  FILE *f=fopen("clip.txt","w");
  for (int t=0;t<400;t++){               // random star-shaped (often concave) polygons, some with a hole
    int n=5+rand()%12; double cx=rand()%1000/100.0, cy=rand()%1000/100.0;
    gw_ring r;
    for (int i=0;i<n;i++){double a=2*M_PI*i/n, rad=0.5+3.0*(rand()%1000)/1000.0; gw_pt p; p.x=cx+rad*cos(a); p.y=cy+rad*sin(a); r.push_back(p);}
    double x0=rand()%1000/100.0-1, y0=rand()%1000/100.0-1, w=0.3+rand()%500/100.0;
    double A=GWGeom::ClipRingAreaToRect(r,x0,y0,x0+w,y0+w);
    fprintf(f,"%d %.17g %.17g %.17g %.17g %.17g",n,x0,y0,x0+w,y0+w,A);
    for (int i=0;i<n;i++){fprintf(f," %.17g %.17g",r[i].x,r[i].y);} fprintf(f,"\n");
  }
  fclose(f);
  CGWProjection P; P.SetOrigin(59.0,-125.3);
  f=fopen("proj.txt","w");
  for (int t=0;t<300;t++){double lon=-129+8.0*(rand()%10000)/10000.0, lat=56+7.0*(rand()%10000)/10000.0,x,y,lo,la;
    P.Forward(lon,lat,x,y); P.Inverse(x,y,lo,la); fprintf(f,"%.10f %.10f %.6f %.6f %.12f %.12f\n",lon,lat,x,y,lo,la);}
  fclose(f);
  CGWRaster R; std::string err; if (!R.Read("test.asc",err)){printf("raster read failed: %s\n",err.c_str()); return 1;}
  f=fopen("raster.txt","w");
  for (int t=0;t<200;t++){double lon=-2+4.0*(rand()%10000)/10000.0, lat=10+3.0*(rand()%10000)/10000.0,v=-1; bool ok=R.Sample(lon,lat,v); fprintf(f,"%.8f %.8f %d %.10f\n",lon,lat,ok?1:0,v);}
  fclose(f);
  printf("done\n"); return 0;
}
