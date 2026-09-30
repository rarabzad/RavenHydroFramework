#include "GWGrid.h"   // compile: g++ -std=c++11 -I../src gridtest.cpp ../src/GWGrid.cpp ../src/GWGeometry.cpp
#include <cstdio>
#include <cmath>
#include <random>
#include <sstream>
using namespace std;
static int fails=0;
void check(const char *name,bool ok,double v){printf("%-58s %s  (%.3g)\n",name,ok?"ok":"FAIL",v); if(!ok)fails++;}
void invariants(const char *tag,CGWGrid &G,double domA){
  double A=0; for(int i=0;i<G.ncpl;i++)A+=G.area[i];
  char b[200]; sprintf(b,"%s: cells tile the domain (area sum / domain - 1)",tag); check(b,fabs(A/domA-1)<1e-9,A/domA-1);
  double asym=0; long nn=0;
  for(int i=0;i<G.ncpl;i++)for(auto &e:G.nbr[i]){nn++; bool f=false; for(auto &g:G.nbr[e.j]) if(g.j==i){f=true; asym=max(asym,fabs(g.len-e.len)/e.len+fabs(g.di-e.dj)/max(e.dj,1e-9));} if(!f)asym=1e9;}
  sprintf(b,"%s: neighbours symmetric (%ld links)",tag,nn); check(b,asym<1e-6,asym);
  double per=0,sh=0; for(int i=0;i<G.ncpl;i++){auto &r=G.poly[i]; for(size_t k=0;k<r.size();k++)per+=hypot(r[(k+1)%r.size()].x-r[k].x,r[(k+1)%r.size()].y-r[k].y); for(auto&e:G.nbr[i])sh+=e.len;}
  double outer=0; for(int i=0;i<G.ncpl;i++){ } 
  int miss=0; for(int i=0;i<G.ncpl;i++){int L=G.Locate(G.cx[i],G.cy[i]); if(L!=i && G.convex[i])miss++;}
  sprintf(b,"%s: centroid locates its own cell (misses)",tag); check(b,miss==0,miss);
  (void)per;(void)sh;(void)outer;
}
int main(){
  vector<gw_refine> none;
  CGWGrid R; R.MakeRegular(0.3,0.7,10000.2,7000.9,1000,none);
  invariants("regular uniform",R,R.ncol*1000.0*R.nrow*1000.0);
  check("regular uniform: bounds exactly x0+c*cs",R.bx0[5]==R.x0+5*R.cs && R.bx1[5]==R.bx0[5]+R.cs,0);
  vector<gw_refine> ref(1); ref[0].kind=1; ref[0].size=250; {gw_pt a{4200,3300},b{5800,3900}; ref[0].geom={a,b};}
  CGWGrid V; V.MakeRegular(0,0,10000,7000,1000,ref);
  double dom=10000.0*7000.0; invariants("regular variable",V,dom);
  double wmin=1e9,wmax=0,rat=0; for(int c=0;c<V.ncol;c++){wmin=min(wmin,V.colW[c]);wmax=max(wmax,V.colW[c]); if(c) rat=max(rat,max(V.colW[c]/V.colW[c-1],V.colW[c-1]/V.colW[c]));}
  check("regular variable: widths 250..1000, neighbour ratio <= 2",fabs(wmin-250)<1e-9&&fabs(wmax-1000)<1e-9&&rat<=2.0001,rat);
  CGWGrid Q; Q.MakeQuadtree(0,0,10000,7000,1000,ref,3);
  invariants("quadtree",Q,dom);
  double lv=0; for(int i=0;i<Q.ncpl;i++)for(auto&e:Q.nbr[i]){double si=Q.bx1[i]-Q.bx0[i],sj=Q.bx1[e.j]-Q.bx0[e.j]; lv=max(lv,max(si/sj,sj/si));}
  check("quadtree: 2:1 balance (max size ratio of neighbours)",lv<=2.0001,lv);
  int hang=0; for(int i=0;i<Q.ncpl;i++) if(Q.poly[i].size()>4) hang++;
  printf("  quadtree cells %d (smallest %.0f m), cells with hanging vertices %d\n",Q.ncpl,Q.MinCellSize(),hang);
  CGWGrid F; F.MakeFromPolygons(Q.poly); invariants("imported (quadtree polygons)",F,dom);
  double dn=0; for(int i=0;i<Q.ncpl;i++){ if(F.nbr[i].size()!=Q.nbr[i].size()) dn=1e9; double a=0,b=0; for(auto&e:F.nbr[i])a+=e.len; for(auto&e:Q.nbr[i])b+=e.len; dn=max(dn,fabs(a-b));}
  check("imported: same neighbours and face lengths as quadtree",dn<1e-6,dn);
  mt19937 rng(7); uniform_real_distribution<double> ux(0,10000),uy(0,7000); vector<gw_pt> s; for(int i=0;i<600;i++) s.push_back({ux(rng),uy(rng)});
  CGWGrid M; M.MakeVoronoi(s,0,0,10000,7000); invariants("voronoi",M,dom);
  double orth=0; for(int i=0;i<M.ncpl;i++)for(auto&e:M.nbr[i]) orth=max(orth,fabs(e.di+e.dj-hypot(s[e.j].x-s[i].x,s[e.j].y-s[i].y)));
  check("voronoi: centres equidistant from faces (orthogonal)",orth<1e-9,orth);
  int nc=0; for(int i=0;i<M.ncpl;i++) nc+=M.convex[i]; check("voronoi: all cells convex",nc==M.ncpl,nc);
  { std::vector<gw_pt> L; for(double y=125;y<7000;y+=250) for(double x=125;x<10000;x+=250) L.push_back({x,y});
    CGWGrid G; G.MakeVoronoi(L,0,0,10000,7000); invariants("voronoi (lattice seeds, cocircular)",G,dom);
    gw_ring cov={{-5000,-5000},{20000,-5000},{20000,15000},{-5000,15000}}; double e=0; int zero=0;
    for(int i=0;i<G.ncpl;i++){double a=G.ClipArea(cov,i); e=std::max(e,fabs(a-G.area[i])/G.area[i]); if(a<1e-6) zero++;}
    check("voronoi (lattice): every cell clips to its full area",e<1e-9&&zero==0,e);
    int nv=0; for(int i=0;i<G.ncpl;i++) nv=std::max(nv,(int)G.poly[i].size()); check("voronoi (lattice): squares have 4 vertices",nv==4,nv); }
  // clipping: a polygon covering the whole domain overlaps each cell by its area
  gw_ring big={{-5000,-5000},{20000,-5000},{20000,15000},{-5000,15000}}; double e1=0;
  for(CGWGrid*G:{&R,&V,&Q,&M,&F}) for(int i=0;i<G->ncpl;i+=7) e1=max(e1,fabs(G->ClipArea(big,i)-G->area[i])/G->area[i]);
  check("clipping: covering polygon returns each cell area",e1<1e-9,e1);
  // concave cell: L-shape clipped by a square
  gw_ring L={{0,0},{2,0},{2,1},{1,1},{1,2},{0,2}}; CGWGrid C1; C1.MakeFromPolygons({L});
  gw_ring sq={{0.5,0.5},{1.5,0.5},{1.5,1.5},{0.5,1.5}}; double a=C1.ClipArea(sq,0);
  check("concave cell: L-shape by square = 0.75",fabs(a-0.75)<1e-12,a);
  stringstream ss; Q.Write(ss); CGWGrid Q2; bool ok=Q2.Read(ss); double d=0; for(int i=0;i<Q.ncpl;i++) d+=fabs(Q2.area[i]-Q.area[i])+fabs((double)Q2.nbr[i].size()-Q.nbr[i].size());
  check("cache round trip (quadtree)",ok&&d==0,d);
  stringstream s2; R.Write(s2); CGWGrid R2; ok=R2.Read(s2); d=0; for(int i=0;i<R.ncpl;i++) d+=fabs(R2.bx0[i]-R.bx0[i])+fabs(R2.by1[i]-R.by1[i]);
  check("cache round trip (regular uniform, bit-exact)",ok&&d==0&&R2.uniform,d);
  printf("%s: %d failures\n",fails?"FAILED":"ALL PASSED",fails); return fails;
}
