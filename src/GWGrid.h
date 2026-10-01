/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: model grids
  (regular with rotation and variable spacing, quadtree, Voronoi mesh, imported cells)
----------------------------------------------------------------*/
#ifndef GWGRID_H
#define GWGRID_H
#include "GWGeometry.h"
#include <iostream>

enum gw_grid_type {GWGRID_REGULAR=0,GWGRID_QUADTREE=1,GWGRID_HRUMESH=2,GWGRID_FILE=3};

///< neighbour of a cell: shared face length, distances from both cell centres to the face, direction of the face normal
struct gw_nbr { int j; double len,di,dj,ang; };  ///< [m],[m],[m],[rad, local frame, pointing from this cell to j]

///< grid refinement feature, local frame: kind 0 point, 1 polyline, 2 polygon; target cell size [m]
struct gw_refine { int kind; double size; gw_ring geom; };

///////////////////////////////////////////////////////////////////
/// \brief a two-dimensional model grid (one layer of cells) in a local, possibly rotated, frame
/// \details local = rotation(-angle) of (projected world - origin). Every cell is a counterclockwise polygon.
/// Regular and quadtree cells are axis-aligned rectangles in the local frame, so overlaps use exact rectangle
/// clipping; Voronoi cells are convex; imported cells may be non-convex (clipped by triangles).
//
class CGWGrid
{
public:
  gw_grid_type type;
  double ox,oy,angle;               ///< local frame: world origin [m], rotation [deg, counterclockwise]
  //-- regular grids
  int    nrow,ncol;                 ///< row 0 is the northernmost row
  bool   uniform;                   ///< regular grid with one cell size (reproduces the original arithmetic exactly)
  double x0,y0,cs;                  ///< lower-left corner and cell size of a uniform grid (local frame)
  vector<double> colX,colW,rowY,rowH; ///< column west edge/width (west->east); row south edge/height (row 0 north)
  //-- every grid
  int    ncpl;                      ///< cells per layer
  vector<gw_ring>  poly;            ///< cell polygons, counterclockwise, local frame
  vector<double>   area,cx,cy,bx0,by0,bx1,by1;
  vector<char>     rect,convex;
  vector<vector<gw_nbr> > nbr;
  vector<gw_pt>    verts;           ///< welded vertices (DISV)
  vector<vector<int> > cverts;      ///< per cell: vertex indices, clockwise (DISV order)

  CGWGrid();
  void   ToLocal(const double X,const double Y,double &x,double &y) const;
  void   ToWorld(const double x,const double y,double &X,double &Y) const;

  void   MakeRegular (const double xa,const double ya,const double xb,const double yb,const double cellsize,
                      const vector<gw_refine> &ref);
  void   MakeQuadtree(const double xa,const double ya,const double xb,const double yb,const double cellsize,
                      const vector<gw_refine> &ref,const int maxLevel);
  void   MakeVoronoi (const vector<gw_pt> &seeds,const double xa,const double ya,const double xb,const double yb);
  void   MakeFromPolygons(const vector<gw_ring> &cells);

  int    Locate      (const double x,const double y) const;
  void   Candidates  (const double xa,const double ya,const double xb,const double yb,vector<int> &out) const;
  double ClipArea    (const gw_ring &r,const int ic) const;
  void   SamplePoints(const int ic,const int n,vector<gw_pt> &pts) const;
  void   SegmentIntervals(const gw_pt &A,const gw_pt &B,const int ic,vector<pair<double,double> > &iv) const; ///< parameter ranges of segment AB inside a cell
  double MinCellSize () const;
  double MeanCellSize() const;
  bool   IsNeighbourCVFD() const {return (type==GWGRID_HRUMESH);}

  void   Write(ostream &O) const;
  bool   Read (istream &I);

private:
  double _bs,_bxa,_bya; int _bnx,_bny;
  vector<vector<int> > _bins;
  void   Finish(const bool neighboursKnown);
  void   ComputeGeometry();
  void   BuildBins();
  void   WeldVertices(const double tol);
  void   NeighboursFromEdges(const double tol);
};

namespace GWGeom
{
  double ClipRingAreaToConvex(const gw_ring &r,const gw_ring &clip);
  bool   PointInRing        (const gw_ring &r,const double x,const double y);
  bool   IsConvex           (const gw_ring &r);
  void   Triangulate        (const gw_ring &r,vector<gw_ring> &tris);
  double SegmentRectLength  (const gw_pt &a,const gw_pt &b,const double x0,const double y0,const double x1,const double y1);
}
#endif
