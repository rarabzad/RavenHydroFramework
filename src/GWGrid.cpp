/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: model grids
----------------------------------------------------------------*/
#include "GWGrid.h"
#include <cmath>
#include <algorithm>
#include <map>
#include <iomanip>

namespace
{
  const double GW_PI=3.14159265358979323846;
  inline double Cross(const gw_pt &o,const gw_pt &a,const gw_pt &b){return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);}
  void DropClosing(gw_ring &r){if ((r.size()>1) && (r.front().x==r.back().x) && (r.front().y==r.back().y)){r.pop_back();}}
  gw_ring RectRing(double xa,double ya,double xb,double yb){gw_ring r(4); r[0].x=xa;r[0].y=ya; r[1].x=xb;r[1].y=ya; r[2].x=xb;r[2].y=yb; r[3].x=xa;r[3].y=yb; return r;}

  //subdivides consecutive intervals of base width cs where refinement features ask for smaller cells; 2:1 smoothing
  void RefineAxis(const double a0,const int n,const double cs,const vector<gw_refine> &ref,const bool xAxis,
                  vector<double> &start,vector<double> &width)
  {
    vector<double> w;
    for (int i=0;i<n;i++){
      double lo=a0+i*cs,hi=lo+cs,target=cs;
      for (size_t f=0;f<ref.size();f++){
        double mn=1e300,mx=-1e300;
        for (size_t k=0;k<ref[f].geom.size();k++){double v=xAxis?ref[f].geom[k].x:ref[f].geom[k].y; mn=min(mn,v); mx=max(mx,v);}
        if ((mx+ref[f].size>=lo) && (mn-ref[f].size<=hi)){target=min(target,ref[f].size);}
      }
      int k=max(1,(int)ceil(cs/target-1e-9));
      for (int j=0;j<k;j++){w.push_back(cs/k);}
    }
    bool changed=true;
    while (changed){
      changed=false;
      for (size_t i=0;i+1<w.size();i++){
        if (w[i]>2.0001*w[i+1]){double h=w[i]/2; w[i]=h; w.insert(w.begin()+i,h); changed=true; break;}
        if (w[i+1]>2.0001*w[i]){double h=w[i+1]/2; w[i+1]=h; w.insert(w.begin()+i+1,h); changed=true; break;}
      }
    }
    start.resize(w.size()); width=w;
    double s=a0; for (size_t i=0;i<w.size();i++){start[i]=s; s+=w[i];}
  }

  //does a feature touch the box? (points inside, polylines crossing, polygons overlapping)
  bool FeatureTouches(const gw_refine &f,double xa,double ya,double xb,double yb)
  {
    const gw_ring &g=f.geom;
    for (size_t k=0;k<g.size();k++){if ((g[k].x>=xa) && (g[k].x<=xb) && (g[k].y>=ya) && (g[k].y<=yb)){return true;}}
    if (f.kind>=1){
      size_t ns=(f.kind==2)?g.size():g.size()-1;
      for (size_t k=0;k<ns;k++){if (GWGeom::SegmentRectLength(g[k],g[(k+1)%g.size()],xa,ya,xb,yb)>0){return true;}}
    }
    if (f.kind==2){if (GWGeom::PointInRing(g,0.5*(xa+xb),0.5*(ya+yb))){return true;}}
    return false;
  }
}

//=====================================================================================================
namespace GWGeom
{
  //Sutherland-Hodgman clipping of any ring by a convex counterclockwise polygon; signed area (holes stay negative)
  double ClipRingAreaToConvex(const gw_ring &r,const gw_ring &clip)
  {
    gw_ring out=r; DropClosing(out);
    size_t m=clip.size();
    double scale=0; for (size_t e=0;e<m;e++){scale=max(scale,max(fabs(clip[e].x-clip[0].x),fabs(clip[e].y-clip[0].y)));}
    for (size_t e=0;(e<m) && (out.size()>0);e++){
      const gw_pt &A=clip[e],&B=clip[(e+1)%m];
      if (hypot(B.x-A.x,B.y-A.y)<=1e-9*scale){continue;} //degenerate edge: its inside is undefined
      gw_ring in=out; out.clear();
      for (size_t k=0;k<in.size();k++){
        const gw_pt &P=in[k],&Q=in[(k+1)%in.size()];
        double cp=Cross(A,B,P),cq=Cross(A,B,Q);
        bool pin=(cp>=0),qin=(cq>=0);
        if (pin){out.push_back(P);}
        if (pin!=qin){double t=cp/(cp-cq); gw_pt I; I.x=P.x+t*(Q.x-P.x); I.y=P.y+t*(Q.y-P.y); out.push_back(I);}
      }
    }
    return (out.size()<3)?0.0:SignedArea(out);
  }
  bool PointInRing(const gw_ring &r,const double x,const double y)
  {
    bool in=false; size_t n=r.size();
    for (size_t i=0,j=n-1;i<n;j=i++){
      if (((r[i].y>y)!=(r[j].y>y)) && (x<(r[j].x-r[i].x)*(y-r[i].y)/(r[j].y-r[i].y)+r[i].x)){in=!in;}
    }
    return in;
  }
  bool IsConvex(const gw_ring &r)
  {
    size_t n=r.size(); if (n<3){return false;}
    for (size_t i=0;i<n;i++){if (Cross(r[i],r[(i+1)%n],r[(i+2)%n])<-1e-9){return false;}}
    return true;
  }
  //ear clipping of a simple counterclockwise polygon
  void Triangulate(const gw_ring &r0,vector<gw_ring> &tris)
  {
    tris.clear(); gw_ring r=r0; DropClosing(r);
    vector<int> idx(r.size()); for (size_t i=0;i<r.size();i++){idx[i]=(int)i;}
    int guard=0;
    while ((idx.size()>3) && (guard<100000)){
      guard++; bool cut=false;
      for (size_t i=0;i<idx.size();i++){
        int a=idx[(i+idx.size()-1)%idx.size()],b=idx[i],c=idx[(i+1)%idx.size()];
        if (Cross(r[a],r[b],r[c])<=1e-12){continue;}
        bool empty=true;
        for (size_t k=0;k<idx.size();k++){
          int p=idx[k]; if ((p==a) || (p==b) || (p==c)){continue;}
          if ((Cross(r[a],r[b],r[p])>=0) && (Cross(r[b],r[c],r[p])>=0) && (Cross(r[c],r[a],r[p])>=0)){empty=false;break;}
        }
        if (!empty){continue;}
        gw_ring t(3); t[0]=r[a]; t[1]=r[b]; t[2]=r[c]; tris.push_back(t);
        idx.erase(idx.begin()+i); cut=true; break;
      }
      if (!cut){break;}
    }
    if (idx.size()==3){gw_ring t(3); t[0]=r[idx[0]]; t[1]=r[idx[1]]; t[2]=r[idx[2]]; tris.push_back(t);}
  }
  //length of segment ab inside an axis-aligned rectangle (Liang-Barsky)
  double SegmentRectLength(const gw_pt &a,const gw_pt &b,const double x0,const double y0,const double x1,const double y1)
  {
    double t0=0,t1=1,dx=b.x-a.x,dy=b.y-a.y;
    double p[4]={-dx,dx,-dy,dy},q[4]={a.x-x0,x1-a.x,a.y-y0,y1-a.y};
    for (int i=0;i<4;i++){
      if (p[i]==0){if (q[i]<0){return 0.0;} continue;}
      double t=q[i]/p[i];
      if (p[i]<0){t0=max(t0,t);} else {t1=min(t1,t);}
      if (t0>t1){return 0.0;}
    }
    return (t1-t0)*sqrt(dx*dx+dy*dy);
  }
}

//=====================================================================================================
CGWGrid::CGWGrid()
{
  type=GWGRID_REGULAR; ox=oy=angle=0.0; nrow=ncol=ncpl=0; uniform=false; x0=y0=cs=0.0;
  _bs=1.0; _bxa=_bya=0.0; _bnx=_bny=0;
}
void CGWGrid::ToLocal(const double X,const double Y,double &x,double &y) const
{
  if (angle==0.0){x=X-ox; y=Y-oy; return;}
  double a=-angle*GW_PI/180.0,dx=X-ox,dy=Y-oy;
  x=dx*cos(a)-dy*sin(a); y=dx*sin(a)+dy*cos(a);
}
void CGWGrid::ToWorld(const double x,const double y,double &X,double &Y) const
{
  if (angle==0.0){X=x+ox; Y=y+oy; return;}
  double a=angle*GW_PI/180.0;
  X=ox+x*cos(a)-y*sin(a); Y=oy+x*sin(a)+y*cos(a);
}

//-- regular grid: one cell size, or variable spacing refined around features ---------------------------
void CGWGrid::MakeRegular(const double xa,const double ya,const double xb,const double yb,const double cellsize,
                          const vector<gw_refine> &ref)
{
  type=GWGRID_REGULAR; cs=cellsize;
  x0=floor(xa/cs)*cs; y0=floor(ya/cs)*cs;
  int nc=max(1,(int)ceil((xb-x0)/cs)),nr=max(1,(int)ceil((yb-y0)/cs));
  uniform=ref.empty();
  if (uniform){
    ncol=nc; nrow=nr;
    colX.resize(ncol); colW.assign(ncol,cs); rowY.resize(nrow); rowH.assign(nrow,cs);
    for (int c=0;c<ncol;c++){colX[c]=x0+c*cs;}
    for (int r=0;r<nrow;r++){rowY[r]=y0+(nrow-1-r)*cs;}
  }
  else {
    vector<double> ys,hs;
    RefineAxis(x0,nc,cs,ref,true ,colX,colW);
    RefineAxis(y0,nr,cs,ref,false,ys,hs);
    ncol=(int)colX.size(); nrow=(int)ys.size();
    rowY.resize(nrow); rowH.resize(nrow);
    for (int r=0;r<nrow;r++){rowY[r]=ys[nrow-1-r]; rowH[r]=hs[nrow-1-r];}
  }
  ncpl=nrow*ncol;
  poly.resize(ncpl); nbr.assign(ncpl,vector<gw_nbr>());
  for (int r=0;r<nrow;r++){
    for (int c=0;c<ncol;c++){
      int ic=r*ncol+c;
      poly[ic]=RectRing(colX[c],rowY[r],colX[c]+colW[c],rowY[r]+rowH[r]);
    }
  }
  Finish(false);
  //exact rectangle bounds (the uniform case reproduces the original arithmetic)
  for (int r=0;r<nrow;r++){for (int c=0;c<ncol;c++){int ic=r*ncol+c; bx0[ic]=colX[c]; by0[ic]=rowY[r]; bx1[ic]=colX[c]+colW[c]; by1[ic]=rowY[r]+rowH[r];}}
}

//-- quadtree: squares split where features ask for smaller cells, balanced so neighbours differ by one level at most
namespace
{
  struct qnode { double x,y,s; int level; int child[4]; };
  int QFind(const vector<qnode> &Q,int n,double x,double y)
  {
    while (Q[n].child[0]>=0){
      double h=Q[n].s/2; int k=((x>=Q[n].x+h)?1:0)+((y>=Q[n].y+h)?2:0);
      n=Q[n].child[k];
    }
    return n;
  }
  void QSplit(vector<qnode> &Q,int n)
  {
    double h=Q[n].s/2;
    for (int k=0;k<4;k++){
      qnode c; c.x=Q[n].x+((k&1)?h:0); c.y=Q[n].y+((k&2)?h:0); c.s=h; c.level=Q[n].level+1;
      for (int i=0;i<4;i++){c.child[i]=-1;}
      Q.push_back(c); Q[n].child[k]=(int)Q.size()-1;
    }
  }
}
void CGWGrid::MakeQuadtree(const double xa,const double ya,const double xb,const double yb,const double cellsize,
                           const vector<gw_refine> &ref,const int maxLevel)
{
  type=GWGRID_QUADTREE; cs=cellsize; uniform=false;
  x0=floor(xa/cs)*cs; y0=floor(ya/cs)*cs;
  int nc=max(1,(int)ceil((xb-x0)/cs)),nr=max(1,(int)ceil((yb-y0)/cs));
  ncol=nc; nrow=nr;
  vector<qnode> Q; Q.reserve(nc*nr*4);
  for (int r=0;r<nr;r++){for (int c=0;c<nc;c++){
    qnode q; q.x=x0+c*cs; q.y=y0+(nr-1-r)*cs; q.s=cs; q.level=0; for (int i=0;i<4;i++){q.child[i]=-1;} Q.push_back(q);}}
  int nroot=nc*nr;
  //refine to feature sizes. For speed, polylines are tested segment by segment, and a square tests only the features
  //that touched its parent (a child's padded box lies inside its parent's); the squares split are exactly the same
  vector<gw_refine> F; vector<double> fx0,fy0,fx1,fy1;
  for (size_t f=0;f<ref.size();f++){
    if ((ref[f].kind==1) && (ref[f].geom.size()>2)){
      for (size_t k=0;k+1<ref[f].geom.size();k++){gw_refine s=ref[f]; s.geom.assign(ref[f].geom.begin()+k,ref[f].geom.begin()+k+2); F.push_back(s);}
    }
    else {F.push_back(ref[f]);}
  }
  for (size_t f=0;f<F.size();f++){
    double a=1e300,b=1e300,c=-1e300,d=-1e300,pad=0.5*F[f].size;
    for (size_t k=0;k<F[f].geom.size();k++){a=min(a,F[f].geom[k].x); b=min(b,F[f].geom[k].y); c=max(c,F[f].geom[k].x); d=max(d,F[f].geom[k].y);}
    fx0.push_back(a-pad); fy0.push_back(b-pad); fx1.push_back(c+pad); fy1.push_back(d+pad);
  }
  vector<vector<int> > cand(Q.size());
  vector<int> all(F.size()); for (size_t f=0;f<F.size();f++){all[f]=(int)f;}
  vector<int> stack; for (int i=0;i<nroot;i++){stack.push_back(i); cand[i]=all;}
  while (!stack.empty()){
    int n=stack.back(); stack.pop_back();
    vector<int> C; C.swap(cand[n]);
    if (Q[n].level>=maxLevel){continue;}
    double target=1e300; vector<int> touch;
    double qa=Q[n].x,qb=Q[n].y,qc=Q[n].x+Q[n].s,qd=Q[n].y+Q[n].s; //the square
    for (size_t t=0;t<C.size();t++){
      int f=C[t]; if ((fx1[f]<qa) || (fx0[f]>qc) || (fy1[f]<qb) || (fy0[f]>qd)){continue;} //bounding boxes apart
      double pad=0.5*F[f].size;
      if (FeatureTouches(F[f],qa-pad,qb-pad,qc+pad,qd+pad)){target=min(target,F[f].size); touch.push_back(f);}
    }
    if (target<Q[n].s*0.999){
      QSplit(Q,n); if (cand.size()<Q.size()){cand.resize(Q.size());}
      for (int k=0;k<4;k++){stack.push_back(Q[n].child[k]); cand[Q[n].child[k]]=touch;}
    }
  }
  //2:1 balance
  double dom_x1=x0+nc*cs,dom_y1=y0+nr*cs;
  bool changed=true;
  while (changed){
    changed=false;
    size_t nQ=Q.size();
    for (size_t n=0;n<nQ;n++){
      if (Q[n].child[0]>=0){continue;}
      double s=Q[n].s,e=s*1e-3;
      double px[8]={Q[n].x+0.25*s,Q[n].x+0.75*s,Q[n].x+0.25*s,Q[n].x+0.75*s,Q[n].x-e,Q[n].x-e,Q[n].x+s+e,Q[n].x+s+e};
      double py[8]={Q[n].y-e,Q[n].y-e,Q[n].y+s+e,Q[n].y+s+e,Q[n].y+0.25*s,Q[n].y+0.75*s,Q[n].y+0.25*s,Q[n].y+0.75*s};
      for (int k=0;k<8;k++){
        if ((px[k]<x0) || (py[k]<y0) || (px[k]>=dom_x1) || (py[k]>=dom_y1)){continue;}
        int rc=(int)floor((px[k]-x0)/cs),rr=nr-1-(int)floor((py[k]-y0)/cs);
        int m=QFind(Q,rr*nc+rc,px[k],py[k]);
        if (Q[m].level<Q[n].level-1){QSplit(Q,m); changed=true;}
      }
    }
  }
  //leaves in a deterministic order: roots row by row (north first), children south-west, south-east, north-west, north-east
  poly.clear();
  vector<int> leaves;
  for (int i=0;i<nroot;i++){
    vector<int> st(1,i);
    while (!st.empty()){int n=st.back(); st.pop_back();
      if (Q[n].child[0]<0){leaves.push_back(n); continue;}
      for (int k=3;k>=0;k--){st.push_back(Q[n].child[k]);}
    }
  }
  ncpl=(int)leaves.size();
  for (size_t L=0;L<leaves.size();L++){
    const qnode &q=Q[leaves[L]]; double s=q.s,e=s*1e-3;
    gw_ring r; gw_pt p;
    //each side: add the midpoint when the neighbour across it is smaller (two different leaves at its quarter points)
    double sx[4]={q.x,q.x+s,q.x+s,q.x},sy[4]={q.y,q.y,q.y+s,q.y+s};
    double ox4[4]={0,e,0,-e},oy4[4]={-e,0,e,0};
    for (int k=0;k<4;k++){
      p.x=sx[k]; p.y=sy[k]; r.push_back(p);
      double ax=sx[k],ay=sy[k],bx=sx[(k+1)%4],by=sy[(k+1)%4];
      double q1x=ax+0.25*(bx-ax)+ox4[k],q1y=ay+0.25*(by-ay)+oy4[k],q3x=ax+0.75*(bx-ax)+ox4[k],q3y=ay+0.75*(by-ay)+oy4[k];
      if ((q1x<x0) || (q1y<y0) || (q1x>=dom_x1) || (q1y>=dom_y1)){continue;}
      int r1=QFind(Q,(nr-1-(int)floor((q1y-y0)/cs))*nc+(int)floor((q1x-x0)/cs),q1x,q1y);
      int r3=QFind(Q,(nr-1-(int)floor((q3y-y0)/cs))*nc+(int)floor((q3x-x0)/cs),q3x,q3y);
      if (r1!=r3){p.x=0.5*(ax+bx); p.y=0.5*(ay+by); r.push_back(p);}
    }
    poly.push_back(r);
  }
  nbr.assign(ncpl,vector<gw_nbr>());
  Finish(false);
  for (int ic=0;ic<ncpl;ic++){const qnode &q=Q[leaves[ic]]; bx0[ic]=q.x; by0[ic]=q.y; bx1[ic]=q.x+q.s; by1[ic]=q.y+q.s; rect[ic]=1;}
}

//-- Voronoi mesh by half-plane clipping; every edge remembers the neighbouring seed that produced it -------------
void CGWGrid::MakeVoronoi(const vector<gw_pt> &seeds,const double xa,const double ya,const double xb,const double yb)
{
  type=GWGRID_HRUMESH; uniform=false; nrow=ncol=0;
  int n=(int)seeds.size(); ncpl=n;
  double bs=sqrt((xb-xa)*(yb-ya)/max(n,1));
  int bnx=max(1,(int)ceil((xb-xa)/bs)),bny=max(1,(int)ceil((yb-ya)/bs));
  vector<vector<int> > B(bnx*bny);
  for (int i=0;i<n;i++){
    int bxi=min(bnx-1,max(0,(int)floor((seeds[i].x-xa)/bs))),byi=min(bny-1,max(0,(int)floor((seeds[i].y-ya)/bs)));
    B[byi*bnx+bxi].push_back(i);
  }
  poly.assign(n,gw_ring()); nbr.assign(n,vector<gw_nbr>());
  vector<vector<int> > labs(n);
  for (int i=0;i<n;i++){
    const gw_pt &P=seeds[i];
    gw_ring c=RectRing(xa,ya,xb,yb); vector<int> lab(4,-1);
    int bxi=min(bnx-1,max(0,(int)floor((P.x-xa)/bs))),byi=min(bny-1,max(0,(int)floor((P.y-ya)/bs)));
    for (int k=0;k<=bnx+bny;k++){
      for (int yy=byi-k;yy<=byi+k;yy++){
        for (int xx=bxi-k;xx<=bxi+k;xx++){
          if ((max(abs(xx-bxi),abs(yy-byi))!=k) || (xx<0) || (yy<0) || (xx>=bnx) || (yy>=bny)){continue;}
          const vector<int> &L=B[yy*bnx+xx];
          for (size_t t=0;t<L.size();t++){
            int j=L[t]; if (j==i){continue;}
            const gw_pt &Qj=seeds[j];
            double dx=Qj.x-P.x,dy=Qj.y-P.y,mx=0.5*(P.x+Qj.x),my=0.5*(P.y+Qj.y);
            gw_ring nc; vector<int> nl;
            size_t m=c.size();
            for (size_t v=0;v<m;v++){
              const gw_pt &A=c[v],&Bp=c[(v+1)%m];
              double fa=(A.x-mx)*dx+(A.y-my)*dy,fb=(Bp.x-mx)*dx+(Bp.y-my)*dy;
              bool ain=(fa<=0),bin=(fb<=0);
              if (ain){nc.push_back(A); nl.push_back(lab[v]);}
              if (ain!=bin){
                double s=fa/(fa-fb); gw_pt I; I.x=A.x+s*(Bp.x-A.x); I.y=A.y+s*(Bp.y-A.y);
                nc.push_back(I); nl.push_back(ain?j:lab[v]);
              }
            }
            c=nc; lab=nl;
            if (c.size()<3){break;}
          }
        }
      }
      double rmax=0; for (size_t v=0;v<c.size();v++){rmax=max(rmax,hypot(c[v].x-P.x,c[v].y-P.y));}
      if (k*bs>2.0*rmax){break;}
    }
    //remove near-duplicate vertices (cocircular seeds leave zero-length edges); keep the label of the edge that survives
    gw_ring cc; vector<int> ll;
    for (size_t v=0;v<c.size();v++){
      const gw_pt &A=c[v],&Bn=c[(v+1)%c.size()];
      if ((c.size()>3) && (hypot(Bn.x-A.x,Bn.y-A.y)<1e-9*bs)){continue;}
      cc.push_back(A); ll.push_back(lab[v]);
    }
    poly[i]=cc; labs[i]=ll;
  }
  Finish(true);
  //neighbours from edge labels: exact shared faces, centres equidistant from the face (orthogonal connections)
  map<pair<int,int>,double> shared;
  for (int i=0;i<n;i++){
    size_t m=poly[i].size();
    for (size_t v=0;v<m;v++){
      int j=labs[i][v]; if (j<0){continue;}
      double len=hypot(poly[i][(v+1)%m].x-poly[i][v].x,poly[i][(v+1)%m].y-poly[i][v].y);
      if (len<1e-9*bs){continue;}
      shared[make_pair(min(i,j),max(i,j))]+=0.5*len;
    }
  }
  for (map<pair<int,int>,double>::iterator it=shared.begin();it!=shared.end();it++){
    int i=it->first.first,j=it->first.second;
    double d=0.5*hypot(seeds[j].x-seeds[i].x,seeds[j].y-seeds[i].y),a=atan2(seeds[j].y-seeds[i].y,seeds[j].x-seeds[i].x);
    gw_nbr e; e.len=it->second; e.di=d; e.dj=d;
    e.j=j; e.ang=a; nbr[i].push_back(e);
    e.j=i; e.ang=a+GW_PI; nbr[j].push_back(e);
  }
}

//-- imported cells ---------------------------------------------------------------------------------------------
void CGWGrid::MakeFromPolygons(const vector<gw_ring> &cells)
{
  type=GWGRID_FILE; uniform=false; nrow=ncol=0;
  ncpl=(int)cells.size(); poly=cells;
  for (int i=0;i<ncpl;i++){
    DropClosing(poly[i]);
    //GIS exports often repeat a vertex; repeated vertices break the convexity test and the triangulation used for
    //clipping, so consecutive (near-)duplicates are removed first
    double tol=1e-7*sqrt(fabs(GWGeom::SignedArea(poly[i])))+1e-9;
    gw_ring q;
    for (size_t v=0;v<poly[i].size();v++){
      if (q.empty() || (hypot(poly[i][v].x-q.back().x,poly[i][v].y-q.back().y)>tol)){q.push_back(poly[i][v]);}
    }
    while ((q.size()>1) && (hypot(q.front().x-q.back().x,q.front().y-q.back().y)<=tol)){q.pop_back();}
    poly[i]=q;
    if (GWGeom::SignedArea(poly[i])<0){reverse(poly[i].begin(),poly[i].end());}
  }
  nbr.assign(ncpl,vector<gw_nbr>());
  Finish(false);
}

//-- geometry, bins, vertices, neighbours ----------------------------------------------------------------------------
void CGWGrid::Finish(const bool neighboursKnown)
{
  ComputeGeometry();
  BuildBins();
  WeldVertices(1e-6*MeanCellSize());
  if (!neighboursKnown && (type!=GWGRID_REGULAR)){NeighboursFromEdges(1e-6*MeanCellSize());}
  if (type==GWGRID_REGULAR){
    nbr.assign(ncpl,vector<gw_nbr>());
    for (int r=0;r<nrow;r++){for (int c=0;c<ncol;c++){
      int ic=r*ncol+c;
      int jj[4]={(c<ncol-1)?ic+1:-1,(r>0)?ic-ncol:-1,(c>0)?ic-1:-1,(r<nrow-1)?ic+ncol:-1};
      for (int k=0;k<4;k++){
        if (jj[k]<0){continue;}
        int rj=jj[k]/ncol,cj=jj[k]%ncol;
        gw_nbr e; e.j=jj[k]; e.ang=k*GW_PI/2;
        if ((k==0) || (k==2)){e.len=rowH[r]; e.di=0.5*colW[c]; e.dj=0.5*colW[cj];}
        else                 {e.len=colW[c]; e.di=0.5*rowH[r]; e.dj=0.5*rowH[rj];}
        nbr[ic].push_back(e);
      }
    }}
  }
}
void CGWGrid::ComputeGeometry()
{
  area.assign(ncpl,0.0); cx.assign(ncpl,0.0); cy.assign(ncpl,0.0);
  bx0.assign(ncpl,0.0); by0.assign(ncpl,0.0); bx1.assign(ncpl,0.0); by1.assign(ncpl,0.0);
  rect.assign(ncpl,0); convex.assign(ncpl,0);
  for (int i=0;i<ncpl;i++){
    const gw_ring &r=poly[i]; size_t m=r.size();
    double A=0,X=0,Y=0,xa=1e300,ya=1e300,xb=-1e300,yb=-1e300;
    for (size_t k=0;k<m;k++){
      const gw_pt &p=r[k],&q=r[(k+1)%m];
      double c=p.x*q.y-q.x*p.y; A+=c; X+=(p.x+q.x)*c; Y+=(p.y+q.y)*c;
      xa=min(xa,p.x); xb=max(xb,p.x); ya=min(ya,p.y); yb=max(yb,p.y);
    }
    A*=0.5; area[i]=A; cx[i]=(A!=0)?X/(6*A):xa; cy[i]=(A!=0)?Y/(6*A):ya;
    bx0[i]=xa; by0[i]=ya; bx1[i]=xb; by1[i]=yb;
    convex[i]=GWGeom::IsConvex(r)?1:0;
    rect[i]=((type==GWGRID_REGULAR) || (type==GWGRID_QUADTREE) || ((m==4) && (fabs(A-(xb-xa)*(yb-ya))<1e-9*fabs(A))))?1:0;
  }
  if (type==GWGRID_REGULAR){for (int i=0;i<ncpl;i++){area[i]=colW[i%ncol]*rowH[i/ncol];}}
}
double CGWGrid::MeanCellSize() const
{
  if (uniform){return cs;}
  double A=0; for (int i=0;i<ncpl;i++){A+=area[i];}
  return (ncpl>0)?sqrt(A/ncpl):1.0;
}
double CGWGrid::MinCellSize() const
{
  if (uniform){return cs;}
  double s=1e300; for (int i=0;i<ncpl;i++){s=min(s,min(bx1[i]-bx0[i],by1[i]-by0[i]));}
  return s;
}
void CGWGrid::BuildBins()
{
  double xa=1e300,ya=1e300,xb=-1e300,yb=-1e300;
  for (int i=0;i<ncpl;i++){xa=min(xa,bx0[i]); ya=min(ya,by0[i]); xb=max(xb,bx1[i]); yb=max(yb,by1[i]);}
  _bs=max(MeanCellSize(),1e-9); _bxa=xa; _bya=ya;
  _bnx=max(1,(int)ceil((xb-xa)/_bs)); _bny=max(1,(int)ceil((yb-ya)/_bs));
  while ((double)_bnx*_bny>4.0*max(ncpl,1)+1000.0){_bs*=1.5; _bnx=max(1,(int)ceil((xb-xa)/_bs)); _bny=max(1,(int)ceil((yb-ya)/_bs));}
  _bins.assign(_bnx*_bny,vector<int>());
  for (int i=0;i<ncpl;i++){
    int ia=max(0,(int)floor((bx0[i]-_bxa)/_bs)),ib=min(_bnx-1,(int)floor((bx1[i]-_bxa)/_bs));
    int ja=max(0,(int)floor((by0[i]-_bya)/_bs)),jb=min(_bny-1,(int)floor((by1[i]-_bya)/_bs));
    for (int j=ja;j<=jb;j++){for (int k=ia;k<=ib;k++){_bins[j*_bnx+k].push_back(i);}}
  }
}
void CGWGrid::Candidates(const double xa,const double ya,const double xb,const double yb,vector<int> &out) const
{
  out.clear();
  if (uniform){ //regular uniform grid: direct index range
    int c0=max(0,(int)floor((xa-x0)/cs)),c1=min(ncol-1,(int)floor((xb-x0)/cs));
    int rlo=max(0,nrow-1-(int)floor((yb-y0)/cs)),rhi=min(nrow-1,nrow-1-(int)floor((ya-y0)/cs));
    for (int r=rlo;r<=rhi;r++){for (int c=c0;c<=c1;c++){out.push_back(r*ncol+c);}}
    return;
  }
  int ia=max(0,(int)floor((xa-_bxa)/_bs)),ib=min(_bnx-1,(int)floor((xb-_bxa)/_bs));
  int ja=max(0,(int)floor((ya-_bya)/_bs)),jb=min(_bny-1,(int)floor((yb-_bya)/_bs));
  for (int j=ja;j<=jb;j++){for (int k=ia;k<=ib;k++){
    const vector<int> &L=_bins[j*_bnx+k];
    for (size_t t=0;t<L.size();t++){int i=L[t]; if ((bx1[i]>=xa) && (bx0[i]<=xb) && (by1[i]>=ya) && (by0[i]<=yb)){out.push_back(i);}}
  }}
  sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
}
int CGWGrid::Locate(const double x,const double y) const
{
  if (uniform){
    int col=(int)floor((x-x0)/cs),row=nrow-1-(int)floor((y-y0)/cs);
    return ((col>=0) && (col<ncol) && (row>=0) && (row<nrow))?row*ncol+col:-1;
  }
  vector<int> C; Candidates(x,y,x,y,C);
  for (size_t t=0;t<C.size();t++){
    int i=C[t];
    if (rect[i]){if ((x>=bx0[i]) && (x<bx1[i]) && (y>=by0[i]) && (y<by1[i])){return i;} continue;}
    if (GWGeom::PointInRing(poly[i],x,y)){return i;}
  }
  return -1;
}
double CGWGrid::ClipArea(const gw_ring &r,const int ic) const
{
  if (rect[ic]){return GWGeom::ClipRingAreaToRect(r,bx0[ic],by0[ic],bx1[ic],by1[ic]);}
  if (convex[ic]){return GWGeom::ClipRingAreaToConvex(r,poly[ic]);}
  vector<gw_ring> tris; GWGeom::Triangulate(poly[ic],tris);
  double A=0; for (size_t t=0;t<tris.size();t++){A+=GWGeom::ClipRingAreaToConvex(r,tris[t]);}
  return A;
}
void CGWGrid::SamplePoints(const int ic,const int n,vector<gw_pt> &pts) const
{
  pts.clear(); gw_pt p;
  if (uniform){ //identical to the original raster sampling
    for (int i=0;i<n;i++){for (int j=0;j<n;j++){
      p.x=x0+(ic%ncol+(i+0.5)/n)*cs; p.y=y0+(nrow-1-ic/ncol+(j+0.5)/n)*cs; pts.push_back(p);}}
    return;
  }
  for (int i=0;i<n;i++){for (int j=0;j<n;j++){
    p.x=bx0[ic]+(i+0.5)/n*(bx1[ic]-bx0[ic]); p.y=by0[ic]+(j+0.5)/n*(by1[ic]-by0[ic]);
    if (rect[ic] || GWGeom::PointInRing(poly[ic],p.x,p.y)){pts.push_back(p);}
  }}
  if (pts.empty()){p.x=cx[ic]; p.y=cy[ic]; pts.push_back(p);}
}
namespace
{
  //Cyrus-Beck: parameter range of segment AB inside a convex counterclockwise polygon; false if outside
  bool ClipSegmentConvex(const gw_pt &A,const gw_pt &B,const gw_ring &P,double &t0,double &t1)
  {
    t0=0.0; t1=1.0; double dx=B.x-A.x,dy=B.y-A.y; size_t m=P.size();
    double scale=0; for (size_t e=0;e<m;e++){scale=max(scale,max(fabs(P[e].x-P[0].x),fabs(P[e].y-P[0].y)));}
    for (size_t e=0;e<m;e++){
      const gw_pt &p=P[e],&q=P[(e+1)%m];
      double ex=q.x-p.x,ey=q.y-p.y; if (hypot(ex,ey)<=1e-9*scale){continue;}
      double num=ex*(A.y-p.y)-ey*(A.x-p.x); //>=0 inside (left of the edge)
      double den=ex*dy-ey*dx;
      if (den==0){if (num<0){return false;} continue;}
      double t=-num/den;
      if (den>0){t0=max(t0,t);} else {t1=min(t1,t);}
      if (t0>t1){return false;}
    }
    return true;
  }
}
void CGWGrid::SegmentIntervals(const gw_pt &A,const gw_pt &B,const int ic,vector<pair<double,double> > &iv) const
{
  iv.clear(); double t0,t1;
  if (rect[ic]){
    double dx=B.x-A.x,dy=B.y-A.y; t0=0; t1=1;
    double p[4]={-dx,dx,-dy,dy},q[4]={A.x-bx0[ic],bx1[ic]-A.x,A.y-by0[ic],by1[ic]-A.y};
    for (int i=0;i<4;i++){
      if (p[i]==0){if (q[i]<0){return;} continue;}
      double t=q[i]/p[i]; if (p[i]<0){t0=max(t0,t);} else {t1=min(t1,t);}
      if (t0>t1){return;}
    }
    iv.push_back(make_pair(t0,t1)); return;
  }
  if (convex[ic]){if (ClipSegmentConvex(A,B,poly[ic],t0,t1)){iv.push_back(make_pair(t0,t1));} return;}
  vector<gw_ring> tris; GWGeom::Triangulate(poly[ic],tris);
  for (size_t k=0;k<tris.size();k++){if (ClipSegmentConvex(A,B,tris[k],t0,t1) && (t1>t0)){iv.push_back(make_pair(t0,t1));}}
}
void CGWGrid::WeldVertices(const double tol)
{
  verts.clear(); cverts.assign(ncpl,vector<int>());
  map<pair<long long,long long>,int> key;
  for (int i=0;i<ncpl;i++){
    for (int k=(int)poly[i].size()-1;k>=0;k--){ //clockwise
      const gw_pt &p=poly[i][k];
      pair<long long,long long> kk((long long)llround(p.x/tol),(long long)llround(p.y/tol));
      map<pair<long long,long long>,int>::iterator it=key.find(kk);
      int v; if (it==key.end()){v=(int)verts.size(); verts.push_back(p); key[kk]=v;} else {v=it->second;}
      if (cverts[i].empty() || (cverts[i].back()!=v)){cverts[i].push_back(v);}
    }
    if ((cverts[i].size()>1) && (cverts[i].front()==cverts[i].back())){cverts[i].pop_back();}
  }
}
//shared faces between cells whose edges lie on the same line with opposite direction (handles T-junctions)
void CGWGrid::NeighboursFromEdges(const double tol)
{
  nbr.assign(ncpl,vector<gw_nbr>());
  vector<int> C;
  for (int i=0;i<ncpl;i++){
    Candidates(bx0[i]-tol,by0[i]-tol,bx1[i]+tol,by1[i]+tol,C);
    for (size_t t=0;t<C.size();t++){
      int j=C[t]; if (j<=i){continue;}
      //a shared face may be made of several collinear or slightly bent segments (T-junctions, cells re-projected
      //from another map projection): each overlapping pair of opposite edges adds its length, and the face's position
      //and normal are the length-weighted averages over those pieces
      double len=0,sx=0,sy=0,px=0,py=0;
      size_t mi=poly[i].size(),mj=poly[j].size();
      for (size_t a=0;a<mi;a++){
        const gw_pt &A=poly[i][a],&B=poly[i][(a+1)%mi];
        double ex=B.x-A.x,ey=B.y-A.y,L=hypot(ex,ey); if (L<tol){continue;}
        ex/=L; ey/=L;
        double ctol=max(tol,1e-3*L); //off-line distance allowed: 0.1% of the edge length (re-projected straight edges bend)
        for (size_t b=0;b<mj;b++){
          const gw_pt &P=poly[j][b],&Q=poly[j][(b+1)%mj];
          if ((fabs((P.x-A.x)*ey-(P.y-A.y)*ex)>ctol) || (fabs((Q.x-A.x)*ey-(Q.y-A.y)*ex)>ctol)){continue;}
          if (((Q.x-P.x)*ex+(Q.y-P.y)*ey)>=0){continue;}
          double s1=(P.x-A.x)*ex+(P.y-A.y)*ey,s2=(Q.x-A.x)*ex+(Q.y-A.y)*ey;
          double lo=max(0.0,min(s1,s2)),hi=min(L,max(s1,s2));
          if (hi-lo>tol){
            double w=hi-lo,sm=0.5*(lo+hi);
            len+=w; sx+=w*ey; sy+=w*(-ex); px+=w*(A.x+sm*ex); py+=w*(A.y+sm*ey);
          }
        }
      }
      if (len<=tol){continue;}
      double sn=hypot(sx,sy); if (sn<=0){continue;}
      double nx=sx/sn,ny=sy/sn,lx=px/len,ly=py/len; //mean outward normal of cell i and a point on the mean face line
      //distances from the two centres to the shared face line
      double di=fabs((cx[i]-lx)*nx+(cy[i]-ly)*ny),dj=fabs((cx[j]-lx)*nx+(cy[j]-ly)*ny);
      gw_nbr e; e.len=len; e.di=di; e.dj=dj;
      e.j=j; e.ang=atan2(ny,nx);       nbr[i].push_back(e);
      swap(e.di,e.dj); e.j=i; e.ang=atan2(-ny,-nx); nbr[j].push_back(e);
    }
  }
}

//-- serialization for the linkage cache --------------------------------------------------------------------------
void CGWGrid::Write(ostream &O) const
{
  O<<setprecision(17)<<(int)type<<" "<<ox<<" "<<oy<<" "<<angle<<" "<<(uniform?1:0)<<" "<<x0<<" "<<y0<<" "<<cs<<" "<<nrow<<" "<<ncol<<" "<<ncpl<<"\n";
  if (type==GWGRID_REGULAR){
    for (int c=0;c<ncol;c++){O<<colX[c]<<" "<<colW[c]<<"\n";}
    for (int r=0;r<nrow;r++){O<<rowY[r]<<" "<<rowH[r]<<"\n";}
    return;
  }
  for (int i=0;i<ncpl;i++){
    O<<poly[i].size(); for (size_t k=0;k<poly[i].size();k++){O<<" "<<poly[i][k].x<<" "<<poly[i][k].y;}
    O<<" "<<bx0[i]<<" "<<by0[i]<<" "<<bx1[i]<<" "<<by1[i]<<" "<<(int)rect[i]<<" "<<nbr[i].size();
    for (size_t k=0;k<nbr[i].size();k++){O<<" "<<nbr[i][k].j<<" "<<nbr[i][k].len<<" "<<nbr[i][k].di<<" "<<nbr[i][k].dj<<" "<<nbr[i][k].ang;}
    O<<"\n";
  }
}
bool CGWGrid::Read(istream &I)
{
  int t,u; I>>t>>ox>>oy>>angle>>u>>x0>>y0>>cs>>nrow>>ncol>>ncpl; if (I.fail()){return false;}
  type=(gw_grid_type)t; uniform=(u==1);
  if (type==GWGRID_REGULAR){
    colX.resize(ncol); colW.resize(ncol); rowY.resize(nrow); rowH.resize(nrow);
    for (int c=0;c<ncol;c++){I>>colX[c]>>colW[c];}
    for (int r=0;r<nrow;r++){I>>rowY[r]>>rowH[r];}
    if (I.fail()){return false;}
    poly.resize(ncpl); nbr.assign(ncpl,vector<gw_nbr>());
    for (int r=0;r<nrow;r++){for (int c=0;c<ncol;c++){poly[r*ncol+c]=RectRing(colX[c],rowY[r],colX[c]+colW[c],rowY[r]+rowH[r]);}}
    Finish(false);
    for (int r=0;r<nrow;r++){for (int c=0;c<ncol;c++){int ic=r*ncol+c; bx0[ic]=colX[c]; by0[ic]=rowY[r]; bx1[ic]=colX[c]+colW[c]; by1[ic]=rowY[r]+rowH[r];}}
    return true;
  }
  poly.assign(ncpl,gw_ring()); nbr.assign(ncpl,vector<gw_nbr>());
  vector<double> a0(ncpl),b0(ncpl),a1(ncpl),b1(ncpl); vector<char> rc(ncpl);
  for (int i=0;i<ncpl;i++){
    size_t m; I>>m; if (I.fail() || (m>100000)){return false;}
    poly[i].resize(m); for (size_t k=0;k<m;k++){I>>poly[i][k].x>>poly[i][k].y;}
    int r; size_t nn; I>>a0[i]>>b0[i]>>a1[i]>>b1[i]>>r>>nn; rc[i]=(char)r; if (I.fail() || (nn>100000)){return false;}
    nbr[i].resize(nn);
    for (size_t k=0;k<nn;k++){I>>nbr[i][k].j>>nbr[i][k].len>>nbr[i][k].di>>nbr[i][k].dj>>nbr[i][k].ang;}
  }
  if (I.fail()){return false;}
  vector<vector<gw_nbr> > keep=nbr;
  Finish(true); nbr=keep;
  for (int i=0;i<ncpl;i++){bx0[i]=a0[i]; by0[i]=b0[i]; bx1[i]=a1[i]; by1[i]=b1[i]; rect[i]=rc[i];}
  return true;
}
