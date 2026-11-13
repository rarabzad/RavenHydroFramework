/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Class CGroundwaterModel: Raven-native groundwater subsystem
  driven by MODFLOW 6 (libmf6) through its XMI interface.

  Raven owns the definition (aquifer classes/profiles in .rvp,
  AQUIFER_PROFILE in .rvh, geometry and stresses in .rvg); this class
  builds a standard MODFLOW 6 model from it, runs it in lock-step
  with Raven, and exchanges fluxes every time step.
----------------------------------------------------------------*/
#ifndef GROUNDWATER_H
#define GROUNDWATER_H

#include "RavenInclude.h"
#include "GWGeometry.h"
#include "GWGrid.h"
#include "MF6Engine.h"
#include "GWSoilZone.h"
#include "GWVadose.h"
#include <fstream>

class CModel;
class CTimeSeries;

struct gw_well       ///< pumping/injection well
{
  long long id; string name;
  double lat,lon,screen_top,screen_bot;  ///< [deg],[m]
  bool   xy;                             ///< given as X,Y in :InputCRS (lat,lon hold Y,X until ResolveInputCRS)
  CTimeSeries *pRate;                    ///< [m3/d], negative = withdrawal (from .rvt :WellFlowRate)
};
struct gw_obswell    ///< monitoring well
{
  long long id; string name;
  double lat,lon,z;                      ///< screen elevation [m]
  bool   xy;                             ///< given as X,Y in :InputCRS (lat,lon hold Y,X until ResolveInputCRS)
  int    node;                           ///< user cell, -1 if outside the active model
};
struct gw_swx        ///< surface-water store <-> aquifer exchange (wetlands, lakes)
{
  string sv; string group; double leakance;   ///< Raven store, HRU group, leakance [1/d]
  int iSV; int kk;                            ///< resolved state-variable and HRU-group indices
};
struct gw_boundary   ///< general-head (GHB) or specified-head (CHD) boundary
{
  string name,file; bool isCHD; bool allLayers;
  double head,cond;                      ///< [m], per-cell conductance [m2/d] (GHB)
  CTimeSeries *pHead;                    ///< optional head time series (.rvt :BoundaryHead)
};

enum gw_layer_type
{
  GWL_AQUIFER,          ///< convertible (unconfined when head below top)
  GWL_CONFINED_AQUIFER, ///< confined, constant transmissivity
  GWL_AQUITARD,         ///< explicit low-K leaky layer
  GWL_AQUICLUDE         ///< no-flow barrier (IDOMAIN=0)
};

struct gw_aquifer_class
{
  string name;
  double Kh;        ///< horizontal hydraulic conductivity [m/d]
  double Kv;        ///< vertical hydraulic conductivity [m/d]
  double Ss;        ///< specific storage [1/m]
  double Sy;        ///< specific yield [-]
  double porosity;  ///< total porosity [-] (transport, reporting; saturated content of the unsaturated column)
  double clapp_b;   ///< Campbell exponent of the unsaturated column [-] (:SoilZoneVadose RICHARDS)
  double psia;      ///< air-entry suction of the unsaturated column [m]
  bool   unsatGiven;///< CLAPP_B and AIR_ENTRY_PRESSURE were given
  bool   syGiven;   ///< SPEC_YIELD was given
};

struct gw_layer
{
  int           iclass;     ///< index of aquifer class
  double        thickness;  ///< [m]; ignored if to_bedrock
  bool          to_bedrock; ///< layer extends down to bedrock surface
  gw_layer_type type;
};

struct gw_profile
{
  string           name;
  vector<gw_layer> layers;
  double init_head_depth;       ///< default initial water table depth below land surface [m]
  double default_bedrock_depth; ///< depth to bedrock below land surface when no bedrock raster [m]
  double riverbed_K;            ///< riverbed hydraulic conductivity [m/d] (RIV)
  double riverbed_thick;        ///< riverbed thickness [m] (RIV)
  double seep_leakance;         ///< seepage-face leakance where water table reaches model top [1/d]
  double recharge_delay;        ///< linear-reservoir recharge lag [d] (reserved)
  double soil_zone_depth;       ///< depth of model top below land surface [m]
  double extinction_depth;      ///< water-table ET extinction depth below land surface [m]; 0 = no water-table ET
  double seep_smooth_depth;     ///< seepage starts this far below the model top and reaches full conductance at it [m]; 0 = sharp
};

struct gw_link
{
  int    k;     ///< HRU index
  int    ic;    ///< grid column (cell in plan view) index
  double area;  ///< overlap area [m2]
};

class CGroundwaterModel
{
private:/*----------------------------------------------------*/
  CModel                  *_pModel;

  //-- user definition -----------------------------------------
  vector<gw_aquifer_class> _aClasses;
  vector<gw_profile>       _aProfiles;
  string _hruGeomFile,_hruIDField;
  string _demFile,_bedrockFile,_libPath,_seepTarget,_riverFile,_riverIDField;
  double _cellSize;          ///< [m]
  double _minCoverage;       ///< fraction of cell covered by aquifer HRUs needed to activate it
  int    _headSaveFreq;      ///< head file save frequency [steps]
  double _dvclose;           ///< Newton head-change tolerance [m]
  bool   _lastRelaxed;       ///< last step converged (MODFLOW's own test) at the relaxed tolerance, 10x the target, after 100 iterations
  bool   _subDaily;          ///< time step shorter than a day: output labels include the time of day

  //-- grid ------------------------------------------------------
  CGWProjection _proj;
  CGWCRS        _inCRS;                 ///< :InputCRS projection of the geometry and raster input files (not set: lon/lat)
  CGWGrid       _grid;                  ///< cells of one layer (regular, quadtree, HRU mesh or imported)
  int    _gridType;                    ///< gw_grid_type
  double _gridAngle;                   ///< rotation of regular and quadtree grids [deg, counterclockwise]
  int    _qtMaxLevel;                  ///< quadtree: maximum number of halvings of the base cell size
  string _gridFile,_gridFileID;        ///< imported cell polygons (GeoJSON) and their ID field
  vector<string> _refKind; vector<string> _refFile; vector<double> _refSize; ///< :GridRefinement requests
  string _flowCorr;                    ///< AUTO, XT3D or NONE
  bool   _connectByElevation;          ///< :LayerConnection ELEVATION (unstructured connections between overlapping layers)
  string _mname;                       ///< MODFLOW model name in memory addresses (GWF for generated models)
  double _hf,_qf,_rf;                  ///< unit factors: m per model length, m3/d per model flow, model rate per m/d (1 for generated models)
  vector<int> _topNode;                ///< per column: user node of the top active cell (empty: layer 1)
  int    _disType;
  //-- existing MODFLOW 6 model (see GWExternal.cpp)
  struct gw_ext_pkg { string ftype,fname,pname; };
  vector<double> _extKv,_extSy;   ///< existing model: vertical conductivity [m/d] and specific yield of every user node (soil zone)
  struct gw_ow      { long long hru; int cell; double w; };  ///< :OverlapWeights row: HRU ID, cell (1..NCPL), A_overlap/A_cell
  string _extSim,_extModel,_extGridFile,_extGridID,_extStart;   ///< .rvg inputs
  string _streamPkg,_streamMap,_streamAux;                      ///< stream package, mapping (AUX or DOMINANT_HRU), aux name
  vector<string> _takeOver,_keepPkgs;                           ///< packages Raven replaces / declared kept
  vector<gw_ow>  _owTable;
  CGWCRS         _extCRS;                                       ///< :MF6CRS projection of the model's coordinates
  vector<gw_ring> _extCellsLL;                                  ///< top-layer cells from MODFLOW memory, lon/lat (with :MF6CRS)
  string         _extGeomKey;                                   ///< fingerprint of those cells (linkage cache)
  bool   _extSeepage;                                           ///< add DRN_RAVEN (seepage to Raven)
  string _extTdisFile,_extNamFile,_extDisFtype,_stoName,_disName;
  vector<gw_ext_pkg> _extPkgs;
  vector<string> _extOther;                                     ///< kept packages whose simulated flows enter the balance
  double _extTf,_extLf,_extDvclose,_extAreaCheck,_Qother;       ///< model time units per day, length units per m; solver tolerance [m]
  vector<double> _extArea,_extTopRaw;                           ///< per column: cell area [m2]; top of the top active cell [model units]
  vector<int>    _userOf;                                       ///< reduced node -> user node
  int    _extAux,_extNaux;
  double _extOffset;                                           ///< [d] Raven's start after the model's time zero (restart)
  int    _extNModels;
  string _mf6Version;                                           ///< version reported by the MODFLOW 6 library
  bool   _extNoSto;                                             ///< the existing model has no storage package (steady state)
  bool   _extSteady;                                            ///< steady-state model: its one period stretched over Raven's run
  int    _extMoverDropped;                                      ///< movers of taken-over packages removed from the MVR copy
  string _extRavenStart;                                        ///< Raven's start date (yyyy-mm-dd)                                           ///< models listed in the simulation (all types)
  int    _extDrop;                                             ///< stress periods elapsed before Raven's start
  map<long long,int> _extSBIndex;                     ///< 0 DIS, 1 DISV, 2 DISU
  double _meshFidelity;                ///< HRU mesh: share of active-cell area inside the dominant HRU
  int    _nrow,_ncol,_nlay,_ncpl;
  double _x0,_y0;                        ///< lower-left corner, projected [m]
  vector<double> _land,_top,_bedrock;    ///< per column [m]
  vector<double> _botm,_strt;            ///< per cell (layer-major) [m]
  vector<int>    _idomain;               ///< per cell
  vector<int>    _reduced;               ///< user cell -> reduced node (0-based), -1 if removed
  vector<int>    _colProfile;            ///< per column: profile index, -1 inactive
  vector<double> _colCover;              ///< per column: aquifer-HRU overlap area on active cells [m2]
  vector<map<int,double> > _hruCells;    ///< per HRU: column -> overlap area [m2]
  vector<double> _hruPolyArea;           ///< per HRU: polygon area [m2]
  vector<gw_link>        _aLinks;
  vector<vector<int> >   _hruLinks;      ///< per HRU: indices into _aLinks
  vector<double>         _hruLinkedArea; ///< per HRU: overlap with active columns [m2]
  vector<int>            _hruProfile;    ///< per HRU: profile index or -1
  vector<int>            _aqHRUs;        ///< indices of coupled HRUs
  vector<int>            _bndCol;        ///< column of each RCH/DRN entry
  vector<double>         _Qcol,_Qseep;   ///< work arrays per column [m3/d]
  vector<int>            _colDomHRU;     ///< per column: HRU with largest overlap
  CGWRaster              _dem;
  //-- river cells (RIV) --
  vector<int>            _rivNode;       ///< user cell (layer-major) of each RIV entry
  vector<int>            _rivCol,_rivSB; ///< column and subbasin index of each RIV entry
  vector<double>         _rivLen,_rivBed;///< river length in cell [m], channel bed elevation [m]
  vector<double>         _sbGain;        ///< per subbasin: net aquifer discharge to river [m3/d]
  vector<int>            _evtCol;        ///< column of each EVT entry
  //-- stresses, observations, initial conditions --
  vector<gw_well>        _aWells;
  vector<gw_obswell>     _aObsWells;
  vector<gw_boundary>    _aBounds;
  map<long long,CTimeSeries*> _wellTS;   ///< well rates parsed from .rvt, by well ID
  map<string,CTimeSeries*>    _headTS;   ///< boundary heads parsed from .rvt, by boundary name
  vector<int>            _welNode,_welWell; vector<double> _welFrac;
  vector<int>            _ghbNode,_ghbB,_chdNode,_chdB;
  double                *_welQ,*_welHcof,*_welRhs,*_ghbHead,*_ghbCond,*_ghbHcof,*_ghbRhs,*_chdHead,*_chdSim;
  int                    _welNcol,_ghbNcol,_chdNcol;
  int                   *_chdNb,*_chdNl,*_chdIb; int _chdIbc; bool _chdPacked; ///< CHD_RAVEN's NBOUND and NODELIST, the model's IBOUND; IBOUND of a held cell; list without its switched-off entries
  vector<int>            _chdAct;    ///< entries of _chdNode in CHD_RAVEN's list this step, in list order (all, unless some are switched off)
  int                    _bndOffNow; long long _bndOffN; ///< head-boundary cells switched off (the step's head at or below the cell bottom): this step, run total of cell-steps
  double                 _Qwel_spec,_Qwel,_Qghb,_Qchd; ///< [m3/d]
  bool                   _steadyInit; double _ssRecharge; ///< steady-state first period requested (.rvg), recharge [mm/d]
  bool                   _doSteady;  ///< steady-state period in this run (off on hotstart; decided per run/member)
  string                 _dryStart;  ///< layers whose active convertible cells all start dry (summary; empty: none)
  string                 _initMode;  double _initValue;   ///< .rvc :InitialGWHeads
  vector<double>         _hotHeads;  int _hotDims[3];     ///< .rvc :GWHeads (hotstart)
  bool                   _active;
  vector<double>         _rchStore;  ///< per HRU: water held in the recharge-delay reservoir [m3]
  map<long long,double>  _hotRchStore; ///< from .rvc :GWRechargeStore, by HRU ID [m3]
  double                 _rchStoreTot; ///< [m3] recharge-delay storage after the step
  vector<int>            _hruReturnSV; ///< per HRU: state variable receiving seepage (soil-less HRUs: SURFACE_WATER)
  map<long long,double>  _hotResSeep;  ///< from .rvc :GWReservoirSeepage, by subbasin ID [m3] (hotstart of seepage in transit)
  double                 _rivBedDepth; ///< fixed channel depth below DEM [m]; <0: depth at Raven reference flow
  int                    _nstp,_stepCount;
  double                 _jday0; int _jyear0,_calendar; ///< for end-of-step date labels
  double                 _tsDur,_tsDt;   ///< run duration and time step [d] (time series read after Initialize)
  vector<double>         _work;      ///< work array per column
  //-- reservoir <-> aquifer (Raven reservoirs with :SeepageParameters on coupled lake HRUs)
  bool                   _resExchange;
  vector<int>            _resP;      ///< per coupled reservoir: subbasin index
  vector<double>         _resZbot;   ///< per coupled reservoir: lowest aquifer bottom under its lake [m] (stage check in the first step)
  vector<int>            _resNode,_resR; vector<double> _resW; ///< per RES entry: cell, reservoir, weight
  double                *_resQ,*_resHcof,*_resRhs; int _resNcol;
  double                 _Qres_spec,_Qres,_resPendingVol; ///< [m3/d] seepage sent, delivered; [m3] seepage awaiting delivery
  double                 _resShortVol;                     ///< [m3] reservoir gain MODFLOW could not supply (drying cell), taken from the reach
  vector<double>         _resUnsent;  ///< per RES entry: seepage of a reservoir with no cell able to take it this step (specified heads, inactive cells) [m3/d]
  vector<double>         _resUse,_resWs; vector<int> _resOff; ///< per RES entry: 1 if its cell can take seepage this step, else 0; per reservoir: weight of those cells, 1 if any cell cannot
  double                 _resUndelVol,_resUndelTot; ///< [m3] reservoir seepage MODFLOW did not deliver, returned to the reservoir's subbasin: this step, run total
  int                    _resUndelN;  ///< warnings about it
  vector<double>         _resK0;      ///< per coupled reservoir: its own seepage constant (:SeepageParameters) [m3/s/m]
  vector<double>         _resA;       ///< per coupled reservoir: share of its lake's linked area its entries cover (generated models: not on specified-head columns)
  vector<double>         _resWk;      ///< per coupled reservoir: share (of its entries' weight) of the cells that could take seepage when Raven last computed it
  double                 _resShortTot; int _resShortN; ///< [m3] run total of reservoir gains MODFLOW did not supply; warnings about them
  map<long long,double>  _hotResFrac; ///< from .rvc :GWReservoirSeepage, by subbasin ID: _resWk of the seepage in transit
  vector<gw_swx>         _aSWX;      ///< surface-water exchanges
  vector<int>            _swxCol,_swxX; ///< per SWX entry: column, exchange index
  vector<vector<int> >   _swxLinks;  ///< per SWX entry: links (HRU,cell) of group HRUs in the column
  vector<vector<double> > _swxVl;    ///< per SWX entry and link: water available at the start of the step [m3]
  double                *_swxHead,*_swxCond,*_swxRbot,*_swxHcof,*_swxRhs; int _swxNcol;
  double                 _Qswx_in,_Qswx_out; ///< [m3/d] store leakage into aquifer; aquifer discharge into stores
  double                 _swxAppliedVol;     ///< [m3] net water added to Raven stores by the exchange this step
  string                 _cacheFile; bool _cacheHit; double _geomSeconds;
  bool                   _writeNC; int _ncid,_ncTime,_ncHead,_ncRec;
  bool                   _writeCBC;      ///< :WriteGWBudgetFile: MODFLOW budget file gwf.cbc, specific discharge, saturation and gwf.porosity
  void                   OpenNetCDF(const optStruct &Options);
  void                   WriteNetCDF(const double &t);
  vector<char>           _colCHDTop; ///< per column: top active cell is a specified-head (CHD) cell
  vector<char>           _colCHDSw;  ///< per column (generated models): its specified-head cells all follow a :BoundaryHead series; _colCHDTop holds on the steps one of them is held
  bool                   _chdSwAny,_resSw; ///< some column is such a column; some lake lies over one (its seepage on their held steps is returned, as in an existing model)
  vector<double>         _hruRchArea;///< per HRU: linked area excluding CHD-topped columns [m2]
  double                 _QrchCHD;   ///< [m3/d] recharge of HRUs lying entirely on CHD cells (passes to the CHD boundary)
  string                 EndDate(const double &t) const;
  void                   ResetState();
  string CellID       (const int l,const int ic) const;
  int    LocateWorld  (const double X,const double Y) const;
  int    LocateActiveOnLine(const double X,const double Y,const double dx,const double dy) const;
  int    NearestActiveCell(const double X,const double Y) const;
  struct gw_piece { int ic; double t0,t1; };
  void   SegmentPieces(const double xa,const double ya,const double xb,const double yb,vector<gw_piece> &out) const;
  void   CellCentreWorld(const int ic,double &X,double &Y) const;
  void   BuildRefinement(vector<gw_refine> &ref) const;
  void   MeshSeeds    (const vector<vector<gw_ring> > &rings,const vector<gw_refine> &ref,const double xa,const double ya,const double xb,const double yb,vector<gw_pt> &seeds) const;
  string GridSpec     () const;
  bool   UseXT3D      () const;
  string MF6LibraryPath() const;
  string MF6LibraryHelp() const;
  void   ExtReadSimulation();
  void   ExtRewriteTDIS   (const optStruct &Options);
  void   ExtEditNameFile  (const bool addRaven);
  void   ExtEditMover     ();
  vector<string> ExtReadOpenClose(const string &name) const;
  void   ExtShiftTimeSeries(const string &name);
  int    _extTSShifted;    ///< time-series files shifted to Raven's later start
  void   ExtRenumberPeriods();
  void   ExtWriteSimFile   ();
  void   ExtUpdateCHDTop   ();
  void   ExtCheckGridFile  ();
  void   ExtLocalCells     (vector<gw_ring> &local,double &X0,double &Y0,double &ca,double &sa);
  vector<pair<int,int> > _extConn; ///< top-layer cell pairs MODFLOW connects (grid-file check)
  vector<gw_pt>  _extLocalC;       ///< MODFLOW's top-layer cell centres in the model frame [m] (grid-file check)
  static long long CivilDays(const string &date);
  void   ExtProbe         ();
  void   ExtCellsFromMemory();
  void   ExtLinkHRUs      (const optStruct &Options);
  void   ExtBuild         (const optStruct &Options);
  void   ExtWritePackages (const optStruct &Options);
  void   ExtStreamFlows   ();
  void   ExtConnect       ();
  void   WriteGridCells(const optStruct &Options);
  double                 MeanRaster(const CGWRaster &R,const int ic,bool &ok) const;
  void                   ReadRaster(CGWRaster &R,const string &file,const optStruct &Options) const;
  void                   CheckRasterCover(const CGWRaster &R,const string &file,const string &cmd,const vector<int> &cols,const int nCov,const string &fallback,const optStruct &Options) const;
  ofstream               _HEADS;
  vector<double>         _sbRivLen;      ///< per subbasin: total river length on active cells [m]
  vector<double>         _sbLossCarry;   ///< per subbasin: river loss not yet supplied, taken in the next step [m3]
  map<long long,double>  _hotRivCarry;   ///< from .rvc :GWRiverLossCarry, by subbasin ID [m3]
  vector<double>         _swxQ;          ///< per SWX entry: flow into the aquifer this step [m3/d]
  vector<double>         _sbLossRemain;  ///< per subbasin: river loss not yet taken from Raven this step [m3]
  vector<double>         _sbInTake,_sbOutTake; ///< per subbasin: river loss taken from the reach's inflow and from its outflow this step [m3]
  vector<double>         _sbOutTakePrev; ///< per subbasin: river loss taken from the reach's outflow in the last step [m3]
  vector<double>         _sbQoutDeliv;   ///< per subbasin: the reach's outflow delivered at the start of the step [m3/s], while its
                                         ///< routing reads the outflow before the last step's take (BeginReachRouting)
  vector<double>         _sbQoutRec;     ///< per subbasin: the reach's stored outflow right after Raven's UpdateOutflows [m3/s] (NaN: not known)
  vector<char>           _sbAddBack;     ///< per subbasin: the last take was added back for this step's routing
  int                    _sbOverN;       ///< steps on which a reach's outflow was set from outside the routing after a take (warnings)
  double                 _rivTake,_rivTakePrev; ///< river loss taken from reach inflows and outflows this step and in the last [m3] (Raven's
                                                ///< channel storage and outflow book such a take half in its step and half in the next)
  map<long long,pair<double,double> > _hotRivTake; ///< from .rvc :GWRiverTake, by subbasin ID: inflow and outflow take of the last step [m3]
  map<long long,double>  _hotQoutRec;    ///< from .rvc :GWRiverTake: the reach's stored outflow after Raven's routing in the last step [m3/s]
  double                 _rivCarryEnd;   ///< river loss still carried at the end of the run [m3]
  int                    _rivCarryN;     ///< steps on which river loss was carried (warnings)
  void                   InitRiverTakes();
  bool                   _budgetLineOpen;
  int    _iGW,_iReturnSV;

  //-- soil zone <-> groundwater (:SoilZone PHYSICAL; docs/design/soil_zone_coupling_plan.md in RavenHydroFramework-extra) --
  int                    _szMode;        ///< 0 CONCEPTUAL (default, unchanged behaviour), 1 PHYSICAL
  int                    _szNLayers;     ///< soil-zone layers SOIL[0..M-1]; the interface layer m* = SOIL[M-1]
  int                    _szMaxIter;     ///< interface updates per step once MODFLOW itself has converged (0 = explicit predictor only)
  int                    _szExtra;       ///< interface updates this step made after MODFLOW had converged
  double                 _szTolH,_szTolQ;///< [m], [-] interface tolerances
  int                    _szSubgrid;     ///< DEM samples per cell side for the saturated fraction (1 = cell mean)
  string                 _szOutGroup;    ///< HRU group written to GWSoilExchange.csv ("" = all)
  vector<char>           _szOn;          ///< per HRU: soil-zone coupled (coupled HRU with soil layers)
  vector<vector<int> >   _szSV;          ///< per HRU: state-variable indices of the soil-zone layers, bottom -> top
  vector<vector<gwsz::sz_layer> > _szLay; ///< per HRU: capacities of those layers
  vector<double>         _szD;           ///< per HRU: soil-zone thickness [m]
  vector<double>         _szLam;         ///< per HRU: area-weighted saturated thickness of its exchange columns [m] (output)
  vector<double>         _szRchV;        ///< per HRU: recharge released to the aquifer this step, before rejection [m3]
  vector<double>         _szRej;         ///< per HRU: rejected recharge of the current iterate [m3]
  vector<double>         _szOvfV,_szNetV;///< per HRU: overflow [m3] and net exchange added to the soil [m3] this step
  vector<double>         _linkF;         ///< per link (_aLinks): saturated fraction of the current iterate (0: not soil-zone)
  vector<vector<double> > _cellDem;      ///< per column: DEM sample elevations for the saturated fraction [m]
  vector<int>            _szxLink;       ///< per SZX entry: link index
  vector<double>         _szxC,_szxZb;   ///< per SZX entry: conductance C^sz [m2/d], soil-zone base [m]
  vector<double>         _szxQ;          ///< per SZX entry: flow into the aquifer of the final solution [m3/d]
  vector<gwsz::sz_cell>  _szxCell;       ///< per SZX entry: sub-grid description of its column (bases at the DEM samples, ds, D, C)
  vector<double>         _szxQc,_szxCg;  ///< per SZX entry: exchange at the current head [m3/d] and its slope -dq/dh [m2/d]
  vector<double>         _szxHc;         ///< per SZX entry: head that sets the connection weights (lags the head inside the iteration) [m]
  vector<double>         _szxHp,_szxHr;  ///< per SZX entry: previous hc and its residual h - hc [m] (-HUGE_VAL: none yet this step)
  vector<double>         _szxFc;         ///< per SZX entry: contact fraction (mean contact weight) of the last local solve
  vector<double>         _szSyncOvf;     ///< per HRU: overflow when Raven's gains were shared among full columns this step [m3]
  bool                   _szCapUpd;      ///< the current update is the last before the interface freezes (boundary supply-limited)
  double                 _szMovedV;      ///< water moved up between soil-zone layers when Raven's changes were shared this step [m3]
  double                 _szWithheldV;   ///< discharge to Raven withheld this step for MODFLOW's residual and the soil-zone shortfall: from seepage and the soil-zone inflow [m3]
  double                 _szWithheldG;   ///< the same, from river gains [m3]
  double                 _szWithheldS;   ///< the same, from surface-store gains [m3]
  int                    _szGuardN;      ///< warnings written by that guard
  double                 _ncWithheldVol; ///< without :SoilZone PHYSICAL: discharge to Raven withheld on unconverged steps for MODFLOW's negative residual, run total [m3]
  int                    _ncWithheldN;   ///< steps on which it was withheld
  vector<char>           _szRchByArea;   ///< per HRU: recharge leaves the bottom layer at a constant rate (PERC_CONSTANT): columns give it up by area
  vector<double>         _szxFw;         ///< per SZX entry: the recharge its column gave up this step [m3] (recharge shares)
  vector<double>         _szxHq,_szxQq;  ///< per SZX entry: head [m] and exchange [m3/d] of the last local solve this step (_szxHq -HUGE_VAL: none)
  vector<double>         _szxBa,_szxBra,_szxBb,_szxBrb; ///< per SZX entry: bracket of the root of h - hc (ends, residuals) [m]
  vector<int>            _szxBk;         ///< per SZX entry: bracket state (0 none, else the end replaced last: 1 a, 2 b)
  vector<double>         _szxLam;        ///< per SZX entry: saturated thickness of its column, current iterate [m]
  vector<double>         _szxRej;        ///< per SZX entry: rejected recharge of the current iterate [m3]
  vector<int>            _szxCol;        ///< per SZX entry: its soil-zone column
  vector<int>            _linkSzx;       ///< per link: SZX entry (-1: none)
  vector<int>            _szColOff;      ///< per HRU: first soil-zone column (exchange columns in _szxOfHRU order, then the rest)
  vector<double>         _szColA;        ///< per soil-zone column: area [m2]
  vector<double>         _szColPhi;      ///< per soil-zone column and layer (bottom up): water [mm]
  bool                   _szColInit;     ///< columns hold a state (set at the first exchange)
  map<const void*,string> _szNoKsat;   ///< soil classes (their soil_struct) whose HYDRAUL_COND was not given (.rvp)
  map<long long,map<int,vector<double> > > _hotSzCols; ///< from .rvc :GWSoilZoneColumns, by HRU ID and cell (column number of the grid, 0 = rest of the HRU)
  vector<double>         _colCoverNS;    ///< per column: overlap area of coupled HRUs without soil zone [m2]
  double                *_szxHead,*_szxCond,*_szxHcof,*_szxRhs; int _szxNcol;
  double                *_szqQ,*_szqHcof,*_szqRhs; int _szqNcol; ///< SZQ_RAVEN (WEL): the part of an SZX flow whose head would lie below the cell bottom
  double                *_szdElev,*_szdCond,*_szdAux,*_szdHcof,*_szdRhs; int _szdNcol; ///< SZD_RAVEN (DRN, smoothed): a boundary MODFLOW keeps (explicit, frozen) where water rises into the soil: one-sided
  vector<char>           _szxOne;        ///< per SZX entry: the boundary is the one-sided SZD_RAVEN
  bool                   _szOneOK;       ///< the current update may give one-sided boundaries (at least 10 MODFLOW solves follow)
  double                 _Qsz_in,_Qsz_out,_szRejVol,_szOvfVol,_szAppliedVol,_szShortVol; ///< step totals [m3/d],[m3]
  int                    _szIter;        ///< interface updates in the last step
  int                    _szShortN;      ///< steps with a soil-zone shortfall (warnings)
  double                 _szResid;       ///< last interface change of the saturated thickness [m]
  int                    _szResE;        ///< SZX entry of that largest change (-1: none)
  double                 _szDQ,_szDQtol; ///< last change of the total exchange and the change the tolerances allow [m3/d]
  bool                   _szConv;        ///< interface converged in the last step
  int                    _szNotConv;     ///< steps whose interface did not converge
  double                 _szTstep;       ///< step length of the current interface iteration [d]
  vector<double>         _szPhiProc;     ///< work: soil water after Raven's processes (flattened per HRU) [mm]
  vector<int>            _szPhiOff;      ///< per HRU: offset into _szPhiProc
  vector<vector<int> >   _szxOfHRU;      ///< per HRU: its SZX entries
  vector<double>         _evtBase;       ///< per EVT entry: rate before the saturated-fraction reduction [m/d]
  double                 _szPrevQ;       ///< total interface flow of the previous iterate [m3/d]
  vector<double>         _szW;           ///< per HRU: soil-zone water after the step [mm]
  double                 _szRatioMin,_szRatioMax; ///< range of S_y(top layer)/drainable porosity(bottom soil layer)
  ofstream               _SZOUT;
  //-- soil zone: unsaturated columns, exchange timing, column classes (docs/design/soil_zone_columns_plan.md) --
  int                    _vzMode;        ///< unsaturated column between the soil-zone base and the water table: 0 NONE, 1 RICHARDS,
                                         ///<  2 KINEMATIC (gravity-only kinematic wave), 3 LAG (linear storage cascade)
  double                 _vzPar;         ///< parameter of the method (LAG: see SetSoilZoneVadose; 0 = default)
  int                    _vzClasses;     ///< unsaturated columns per exchange column, by depth class of its DEM samples (1: one, at the mean base)
  double                 _vzAcc;         ///< largest content change of a column in one substep [-]
  int                    _vzN;           ///< nodes of a column at most
  int                    _szCapRise;     ///< capillary rise into the soil zone: -1 default (on with RICHARDS), 0 off, 1 on
  int                    _szTiming;      ///< 0 AFTER (the exchange after Raven's processes), 1 LAGGED (default: last step's exchange as a source)
  int                    _szClasses;     ///< column classes per HRU for Raven's processes (1: the HRU mean)
  vector<vector<int> >   _szCls;         ///< per HRU: class of each of its columns this step (empty: one class)
  vector<vector<double> > _szClsW;       ///< per HRU: area weight of each class (sum 1)
  vector<vector<double> > _szClsPhi;     ///< per HRU: class mean soil water at the start of the step [mm], c*n+j (j bottom up)
  vector<vector<double> > _szClsRes;     ///< per HRU: class soil water after Raven's processes [mm], c*n+j
  vector<vector<double> > _szClsGW;      ///< per HRU: class GROUNDWATER pass-through (recharge) after the processes [mm over the class]
  vector<char>           _szClsDone;     ///< per HRU: class results recorded this step
  vector<double>         _szxQprev;      ///< per SZX entry: last step's booked exchange [m3/d] (the LAGGED source; hotstart)
  vector<double>         _szxLagL;       ///< per SZX entry: lagged source given to the column before Raven's processes [mm over the column]
  vector<double>         _szxLagD;       ///< per SZX entry: part of the lagged source not taken back at the sync, applied in SZApply [mm]
  vector<double>         _szLagNet;      ///< per HRU: water the lagged source and its reversal left in Raven's soil zone [m3]
  vector<double>         _szxDebt;       ///< per SZX entry: lagged source Raven's processes used that the step's exchange did not supply [m3], repaid by the aquifer
  vector<double>         _szxRepay;      ///< per SZX entry: repayment of _szxDebt this step [m3/d] (a fixed rate out of the aquifer through SZQ_RAVEN)
  double                 _szDebtVol;     ///< new debt this step [m3]
  double                 _szRepayQ;      ///< repayment this step, all links [m3/d]
  map<long long,map<int,double> > _hotLag; ///< hotstart :GWSoilZoneLag: HRU ID -> cell (1-based) -> last step's exchange [m3/d]
  map<long long,map<int,double> > _hotDebt; ///< hotstart :GWSoilZoneLag: HRU ID -> cell (1-based) -> outstanding debt [m3]
  double                 _szDebtHotUnused; ///< debt of the hotstart file this run does not repay (no LAGGED timing, or no such link) [m3]
  void   SZClassify     ();                                     ///< column classes of every soil-zone HRU for this step
  void   SZLagApply     (double **aPhinew,double **aPhi,const double &tstep);  ///< LAGGED: last step's exchange to the columns and Raven's soil zone
  void   SZLagReverse   (const int k);                          ///< LAGGED: takes the source back out of the columns after the sync
  void   SZRepayRates   (const double &tstep);                  ///< LAGGED: the aquifer's repayment of the debt this step
  //columns over specified-head cells (soil_zone_columns_plan.md section 5): the exchange with the fixed head is the local
  //solve's, booked against the specified head (no MODFLOW boundary)
  vector<int>            _szcLink;       ///< per CHD entry: link
  vector<int>            _szcCol;        ///< per CHD entry: its column
  vector<gwsz::sz_cell>  _szcCell;       ///< per CHD entry: sub-grid closure
  vector<double>         _szcRej;        ///< per CHD entry: rejected recharge this step [m3]
  vector<double>         _szcFw;         ///< per CHD entry: recharge share (the column's loss to recharge) [m3]
  vector<double>         _szcQ;          ///< per CHD entry: exchange into the specified head this step [m3/d]
  vector<double>         _szcCapL;       ///< per CHD entry and layer: capillary rise from the specified head [m3/d], e*n+j
  vector<vector<int> >   _szcOfHRU;      ///< per HRU: its CHD entries
  vector<int>            _linkSzc;       ///< per link: CHD entry or -1
  double SZChdExchange  (const int k,const int c,const gwsz::sz_cell &cell,const double rejV,const double h,const double &tstep,const double D=0.0) const;
  double                 _QszCHD_in;     ///< soil zone to specified heads this step [m3/d]
  double                 _QszCHD_out;    ///< specified heads to soil zone this step [m3/d] (contact and capillary rise)
  double                 _szRchCHD;      ///< recharge of soil-zone columns onto specified heads this step [m3/d]
  double                 _QrchCHD0;      ///< recharge onto specified heads of HRUs without a soil zone [m3/d]
  void   SZCHDSolves    (const double &tstep);                  ///< local solves of the columns over specified heads
  double VZRchW         (const int k,double &Aall) const;       ///< sum of the recharge shares of HRU k's columns; Aall its links' area
  double SZCellKv       (const int ic) const;                   ///< vertical conductivity of a column's top active cell [m/d]
  double SZCellSy       (const int ic) const;                   ///< specific yield of a column's top active cell [-]
  double SZNodeSy       (const int n) const;                    ///< specific yield of a cell (user node) [-]
  void   SZUpNodes      ();                                     ///< _szxUp from the heads of the last step
  const gwvz::vz_mat *VZGapMat(const int e,const double zb=-HUGE_VAL) const; ///< aquifer material just below the soil-zone base of SZX entry e (or below zb), else NULL
  double SZRiseLayers   (const int k,const int c,const double zb,const double zs,const double s1,const double a,const double Ae,
                         const double Sy,const double tstep,double *QL,double *psisLow,const gwvz::vz_mat *gap=NULL,const double *taken=NULL) const; ///< capillary rise into the soil-zone layers [m3/d]
                                                                ///<  (taken: what other depth classes give each layer this step [m3/d])
  double SZSyncClasses  (const int k,const double *target,const double R,vector<double> &lossR); ///< sync class by class
  bool                   _szClsWarned;   ///< warned that classes need ORDERED_SERIES
  string                 _szExtClass;    ///< :SoilZoneAquiferClass (existing models: retention of the cells for the unsaturated columns)
  vector<gwvz::vz_column> _vzcCol;       ///< unsaturated columns, one per depth class of every SZX entry e: _vzOff[e].._vzOff[e+1]-1
                                         ///<  (:SoilZoneVadoseClasses K; with K = 1 one column per entry, at its mean base)
  vector<int>            _vzOff;         ///< per SZX entry (and one past the last): its first class column in _vzcCol
  vector<double>         _vzcA;          ///< per class column: its share of the entry's area (n_c/N of the DEM samples) [-]
  vector<double>         _vzcZt;         ///< per class column: its top (the soil-zone base of its class) [m]
  vector<vector<double> > _vzcZbs;       ///< per class column: soil-zone bases of its DEM samples, ascending [m] (K > 1)
  vector<char>           _vzcUse;        ///< per class column: it carried its share of this step's recharge
  vector<double>         _vzcF0;         ///< per class column: saturated fraction its inflow used this step [-]
  vector<double>         _vzcOut;        ///< per class column: its outputs to the aquifer this step [m3] (bottom flux and release)
  vector<double>         _vzcRq;         ///< per class column: rate of its outputs on the unsaturated share, _vzcOut/(1-F0_c)/dt [m3/d]
  vector<double>         _vzAbsN;        ///< per SZX entry: water the moving bottoms of its columns took from the source cell this step [m3/d] (>= 0)
  vector<double>         _vzQN;          ///< per SZX entry: upflow through the bottoms of its columns this step [m3/d] (>= 0)
  vector<double>         _vzBackW;       ///< per SZX entry: water of its vanished columns beyond Sy dh, given back to the soil zone this step [m3] (rule (ii))
  vector<double>         _vzPass;        ///< per SZX entry: negative recharge share of its columns passed on to the cell this step [m3/d] (never into a column)
  vector<double>         _vzRl;          ///< per SZX entry: its recharge share this step, V share / dt [m3/d] (SZRecharge; for SZRejRate)
  double                 _vzBackV,_vzBackWV; ///< this step, all entries: column outputs handed back on newly saturated ground (rule (i)) and _vzBackW [m3]
  double VZBackRate     (const int e,const double h) const;    ///< rule (i): outputs of entry e's columns handed back to the soil zone at head h [m3/d]
  double                 _vzClsDev;      ///< largest |sum_c a_c F_c - F| over the entries and a range of heads (check of the class fractions)
  double                 _vzQ0;          ///< :SoilZoneVadoseInitialFlux [mm/d]: steady start of the columns on a cold start (0: at rest)
  void   VZClassSetup   ();                                     ///< depth classes and class columns of every SZX entry (VZSetup)
  double VZClassF       (const int e,const int i,const double h) const; ///< saturated fraction of class column i of entry e at head h [-]
  vector<gwvz::vz_strat>  _vzStrat;      ///< per SZX entry: materials by elevation
  vector<char>           _vzUse;         ///< per SZX entry: the column carried this step's recharge
  vector<double>         _vzQ;           ///< per SZX entry: recharge to the aquifer from the column this step [m3/d] (bottom flux)
  vector<double>         _vzAbsQ;        ///< per SZX entry: water a rising water table took from the column, given to the aquifer this step [m3/d]
  vector<double>         _vzRejV;        ///< per SZX entry: inflow the saturated column top could not take [m3]
  vector<double>         _vzFp;          ///< per SZX entry: saturated fraction the column's inflow used this step
  vector<double>         _vzInV,_vzCapV; ///< per SZX entry: inflow [m3] and capillary rise into the soil zone [m3] this step
  long                   _vzSlowed;       ///< steps in which a column's bottom followed a rising water table more slowly (released at most Sy dh)
  vector<double>         _vzCapL;        ///< per SZX entry and soil-zone layer (bottom up, e*n+j): its fraction of _vzCapV
  vector<double>         _szxCapQ;       ///< per SZX entry: capillary rise from the water table into the soil zone this step [m3/d] (a fixed rate taken through SZQ_RAVEN)
  vector<int>            _szwStart;      ///< SZW_RAVEN: entries of SZX entry e are _szwStart[e].._szwStart[e+1]-1 (the active cells below its top cell)
  vector<int>            _szwNode;       ///< SZW_RAVEN: user node of each entry
  double                *_szwQ;          ///< SZW_RAVEN rates (MODFLOW's array; NULL without the package)
  int                    _szwNcol;       ///< stride of _szwQ (0: Q array, 1: BOUND)
  vector<int>            _szxUp;         ///< per SZX entry this step: the cell upward flows come from (its top cell; a lower wet cell when the top cell is dry, through SZW_RAVEN; -1 none)
  vector<char>           _szxUpCHD;      ///< per SZX entry this step: _szxUp is a specified-head cell (the head supplies the rise; no well there)
  mutable int           *_ibPtr;         ///< MODFLOW's IBOUND (reduced nodes), fetched on first use
  bool   NodeCHD        (const int n);                          ///< user node n is a specified head now (IBOUND < 0)
  bool                   _farWarned;     ///< the warning about heads far outside the model was written
  bool                   _mf6Newton;     ///< the MODFLOW model uses NEWTON (generated models always; existing ones as their name file says)
  vector<double>         _szxUpF;        ///< per SZX entry this step: scale of the upward flows out of _szxUp (1; less when that cell is within max(5 %, 5 cm) of drying)
  vector<double>         _vzUpA;         ///< per SZX entry: water a falling water table's new column nodes took from a lower cell this step [m3/d] (through SZW_RAVEN)
  double SZUpCap        (const int e,const double &tstep);      ///< most an upward flow may take out of SZX entry e's source cell this step [m3/d]
  vector<double>         _vzUp;          ///< per SZX entry: upward flow from a lower wet cell into the unsaturated column this step [m3/d] (through SZW_RAVEN)
  double                 _vzUpQ;         ///< all of _vzUp [m3/d] (MODFLOW recharge, negative)
  vector<double>         _szxCapL;       ///< per SZX entry and soil-zone layer (bottom up, e*n+j): its share of _szxCapQ [m3/d]
  double                 _szCapWT;       ///< capillary rise from the water table, all links [m3/d]
  void   SZCapRiseWT    (const double &tstep);                ///< capillary rise from the water table where there is no unsaturated column
  bool   SZCapOn        () const {return (_szCapRise==1) || ((_szCapRise<0) && (_vzMode>0));}
  bool   VZOwnRise      () const {return (_vzMode==1);}  ///< capillary rise comes out of the unsaturated column's top (RICHARDS);
                                                         ///<  otherwise from the water table through it (SZCapRiseWT), steady
  string VZName         () const {return (_vzMode==1)?"RICHARDS":((_vzMode==2)?"KINEMATIC":((_vzMode==3)?"LAG":"NONE"));}
  double VZColStep      (gwvz::vz_column &c,const gwvz::vz_strat &S,const double dt,const double qtop,double &rej,gwvz::vz_step_info *inf,const bool noUp) const;
  double                 _vzInVol,_vzOutVol,_vzCapVol,_vzStoreVol,_vzStore0; ///< step totals [m3] and storage above the reference [m3]
  long long              _vzSub,_vzNewton; ///< substeps and Newton iterations of the columns (summary)
  map<long long,map<int,vector<double> > > _hotVz; ///< from .rvc :GWVadoseColumns: by HRU ID and cell: zbot, contents
  map<long long,map<int,map<int,vector<double> > > > _hotVzc; ///< from .rvc :GWVadoseClassColumns: by HRU ID, cell and class: classes of the cell, zbot, contents
  map<long long,map<int,map<int,vector<double> > > > _hotVzt; ///< from .rvc :GWVadoseColumnTops: by HRU ID, cell and class (1 without classes): column top [m], water above its reference profile [m3]
  string                 _hotVzMethod;   ///< from .rvc: the :SoilZoneVadose method that wrote the unsaturated columns ('' if not given)
  int                    _vzHotBadN;     ///< hotstart column lines written for other columns (method's reference, node count, classes, top)
  double                 _vzChdV; int _vzChdN; ///< [m3] run total of the water of unsaturated columns ended on cells that are specified heads (existing models); columns ended
  double                 _vzHotUnusedV; int _vzHotUnusedN; ///< [m3] water of hotstart column lines a run without columns does not use; lines
  void   VZSetup        ();                                    ///< materials of every exchange column (after SZBuildEntries)
  void   VZInit         ();                                    ///< columns from the initial heads or the hotstart
  void   VZStep         (const double &tstep);                 ///< the columns over the step (predictor heads), capillary rise
  double VZLinkShare    (const int k,const int li) const;      ///< link li's share of HRU k's recharge (as SZRecharge)
  double VZStorage      () const;                              ///< water above the reference contents in all columns [m3]
  bool   SZPhysical() const {return _szMode==1;}
  void   SZSetup        ();                                    ///< after HRUs are known: layers, thickness, checks
  void   SZBuildEntries ();                                    ///< after links: SZX entries, conductances, bases
  double SZBase         (const int k,const int ic) const;      ///< soil-zone base of HRU k in column ic [m]
  double SWXDatum       (const size_t e) const;                 ///< bed and stage datum of surface-store exchange entry e [model units]
  void   SZLinkFractions(const bool relax,const bool mf6conv=true); ///< _szxHc and _linkF from the current heads (hc iterated inside the step)
  void   SZRecharge     (const double &tstep,const bool write);///< _Qcol, _szRej, _szxRej from _szRchV and _linkF (write: RCH array too)
  void   SZLocalSolves  (const double &tstep,const bool steady);///< _szxLam, _szxQc, _szxCg per exchange column
  int    SZNCols        (const int k) const {return (int)_szxOfHRU[k].size()+((k<(int)_szcOfHRU.size())?(int)_szcOfHRU[k].size():0)+1;}
  int    SZFreezeAt     () const {return (_maxiter>=2)?max(1,_maxiter-max(2,_maxiter/10)):0;} ///< last outer iteration after which the interface may be updated (<1: predictor only)
  double SZSyncColumns  (const int k,const double *target,const double R=0.0,vector<double> *lossR=NULL); ///< columns of HRU k to its soil-zone water target[j] [mm] (R: recharge out of the bottom layer [mm], lossR: each column's part); returns the overflow [mm over the HRU]
  vector<int> SZColumnCells(const int k) const;               ///< cell of each column of HRU k (column number of the grid, 1-based; 0 = rest)
  int    SZHotColumns   (const int k,vector<vector<double> > &cols) const; ///< hotstart columns of HRU k: 1 found, 0 none, -1 not matching
  void   SZInitColumns  ();                                    ///< columns at the start: hotstart, else saturated below the initial water table
  double                 _szInitVol;     ///< water added to the soil zone at the start to saturate it below the initial water table [m3]
  int                    _szInitHot;     ///< HRUs whose columns came from the hotstart file
  void   SZSetArrays    (const bool steady);                   ///< SZX head/conductance, RCH and EVT rates
  bool   SZIterate      (double **aPhinew,const int iter,const bool mf6conv,const bool steady,double &prevQ,const bool last);
  void   SZApply        (double **aPhinew,const double &tstep);///< final exchange, rejection and overflow to Raven
  double SZTrialShortfall(const double &tstep) const;          ///< shortfall [m3] the booked exchange would leave (SZApply without changing anything)
  void   SZWriteStep    (const double &t,const string &dstr);
  //-- soil zone interface: the head-dependent rejection term, the kink drain, the
  //   per-link convergence tests, negative recharge, drained columns, divergence and the MODFLOW side
  vector<double>         _szxRs;         ///< per SZX entry: slope s of the rejection's linear term s (hc - h) MODFLOW solves with [m2/d]
  vector<double>         _szxRq;         ///< per SZX entry: head-dependent rejected recharge of the link at its connection head [m3/d] (SZRejRate(hc))
  vector<double>         _szxRqP;        ///< per SZX entry: its slope dR/dh at the connection head [m2/d]
  vector<double>         _szxRl;         ///< per SZX entry: recharge rate of the link whose rejection follows F(h) (no column: R_l) [m3/d]
  vector<double>         _szxRqc;        ///< per SZX entry: rate of the columns' outputs that follow the head (rule (i), R_q = out/(1-F0)) [m3/d]
  vector<double>         _szxF0;         ///< per SZX entry: saturated fraction the columns' inflow used this step (F0)
  vector<char>           _szxDrn;        ///< per SZX entry: MODFLOW has the exchange as the kink drain (SZD_RAVEN, water rising into the soil)
  vector<char>           _szxKink;       ///< per SZX entry: the column has no free water at the start of the local solve (its exchange cannot reverse: kink drain)
  vector<double>         _szxNegCap;     ///< per SZX entry: what the top cell can still give to a negative recharge this step [m3/d]
  vector<double>         _vzUsed;        ///< per SZX entry: water the unsaturated columns took out of the source cell this step [m3/d]
  vector<char>           _szxDropped;    ///< per SZX entry: the connection test failed on the last update (its next update is counted)
  vector<char>           _cellConv;      ///< per user node: convertible cell (ICELLTYPE != 0; generated: layer type; existing: NPF memory)
  vector<double>         _xPrevIt;       ///< MODFLOW's iterate before the solve of the last outer iteration the limit allows (model units)
  vector<double>         _xStart;        ///< heads at the start of the step, taken only from converged steps (far-head test; model units)
  vector<double>         _hotXStart;     ///< from .rvc :GWStartHeads: _xStart of the run that wrote it, when its last step did not converge
  vector<int>            _keptCol;       ///< per grid column: times KeepColumnsWet set the column in this step
  long                   _szCapN;        ///< supply-limited link updates (last update before the freeze) in the run
  double                 _szCapExcess;   ///< largest take of such a boundary beyond the column's free water [m3/d] (reported, must be 0)
  long                   _szDropN,_szDropPlain; ///< link updates the connection test failed (bracket dropped) in the run, and of them those whose next update in the step was the plain one
  int                    _szDropSteps;   ///< steps with such a drop
  bool                   _szDropThis;    ///< a drop in the current step
  double                 _szEsum;        ///< summed |E| of the non-drain links of the last per-link boundary test [m3/d]
  int                    _szWhy;         ///< test that left the last update unconverged with soil levels and total exchange within the tolerances: 0 none, 1 soil-zone shortfall, 2 connection head
  int                    _szWhyE;        ///< SZX entry of that test's largest failure (-1: none)
  double                 _szWhyV,_szWhyFa,_szWhyFb; ///< shortfall [m3], or connection head less MODFLOW's head [m] and the saturated fractions at both
  double                 _szNegCutQ;     ///< negative recharge cut by the top cells' limit this step [m3/d]
  double                 _szNegUnheldV;  ///< of that cut, what the HRUs' soil zones did not hold this step [m3]
  double                 _negCutVol,_negUnheldVol; ///< run totals of the two [m3]
  int                    _negUnheldN;    ///< warnings about negative recharge not held
  bool                   _negWarned;     ///< warned about unlimited negative recharge (HRUs without a soil zone)
  double                 _Qgross;        ///< gross flow through the aquifer this step: every entry's inflow and outflow apart [m3/d] (floor 1)
  double                 _QabsOther;     ///< the existing model's other packages: every entry's |q| [m3/d]
  double                 _rchLostQ,_rchLostVol; ///< recharge MODFLOW did not apply (existing model without NEWTON) this step [m3/d] and in the run [m3] (signed)
  int                    _rchLostN;      ///< warnings about it
  long                   _wetKeptIt,_wetKeptSteps; ///< outer iterations and steps in which KeepColumnsWet set a drained column
  int                    _farSlowN;      ///< warnings about far heads on unconverged steps
  string                 _stepEnd;       ///< end of the current step (date, and time for sub-daily steps), for the per-step warnings
  string                 _summaryFile;   ///< GWModelSummary.txt (run totals appended at the end; "" once written)
  void   SZNegCaps      (const double &tstep);                 ///< _szxNegCap: what each exchange link's top cell can still give this step
  double SZRejRate      (const int e,const double &h,double *slope=NULL) const; ///< head-dependent rejected recharge of SZX entry e at head h [m3/d]
  bool   SZRejFollows   (const int e) const;                   ///< SZX entry e has a head-dependent rejection
  bool   SZDivergedIterate() const;                            ///< an iterate no interface update may be made from (not finite, or far outside at a link)
  double SZSolvedShortfall(const double &tstep,const vector<double> &rej,int *emax=NULL) const; ///< (emax: the SZX entry of the largest one) shortfall [m3] the boundary MODFLOW just solved with (and the rejected recharge rej [m3] it goes with) would leave in the columns at its current heads
  int    KeepColumnsWet (const int it);                        ///< before a solve (NEWTON): drained columns set 1 mm above their bottoms; returns the number set
  double FarMargin      (double &zmin,double &zmax) const;     ///< divergence margin max(10 x relief, 1 km) [m] and the model's elevation range
  static double FarDistance(const double x,const double zmin,const double zmax); ///< how far elevation x lies outside [zmin, zmax] [m]
  void   AppendRunTotals();                                    ///< run totals at the end of GWModelSummary.txt

  //-- engine -----------------------------------------------------
  CMF6Engine _mf6;
  string     _mf6dir;
  int        _maxiter;
  double    *_X,*_rchBound,*_rchHcof,*_rchRhs,*_drnHcof,*_drnRhs,*_strgss,*_strgsy;
  double    *_rivStage,*_rivCond,*_rivRbot,*_rivHcof,*_rivRhs;
  double    *_evtRate,*_evtHcof,*_evtRhs; int _evtNcol;
  int        _rivNcol;   ///< 0 if named arrays (STAGE, COND, RBOT) exist, else 3 (BOUND columns)

  //-- step budget -------------------------------------------------
  double _ravenRchVol;  ///< [m3] recharge Raven sent this step
  double _Qrch,_Qseep_tot,_Qsto,_Qerr; ///< [m3/d]
  double _returnVol;    ///< [m3] water returned to Raven stores this step
  double _Qevt;                       ///< [m3/d] water-table evapotranspiration
  double _Qriv_in,_Qriv_out;          ///< [m3/d] river leakage into aquifer, aquifer discharge to rivers
  double _riverAppliedVol,_riverDeficitVol; ///< [m3] net river exchange applied to Raven reaches; loss Raven could not supply
  ofstream _CONN;
  int    _nNonConverged,_lastIter;
  bool   _lastConverged;
  ofstream _BUDGET,_HRUSTATE;

  int  UserNode(const int lay,const int ic) const {return lay*_ncpl+ic;}
  int  TopNode (const int ic) const {return _topNode.empty()?ic:_topNode[ic];} ///< user node of the top active cell of a column
  void BuildGeometry  (const optStruct &Options);
  void BuildVertical  (const optStruct &Options);
  void BuildLinks     ();
  void BuildRivers    ();
  void BuildStresses  ();
  double _gridEst; string _gridWhat; ///< cells per layer of the grid being built (estimated before it is built) and the command that sets it
  void SetBoundaryArrays(const double &t,const double &tstep);
  int  CHDEntriesOn(const double &t,const double &tstep,vector<char> &on,vector<double> *hs) const;
  bool ResUndelivered() const {return (_resNode.size()>0) && (IsExternal() || _resSw);} ///< seepage the cells under a lake cannot take is returned and booked
  void UpdateSwitchedCHDTop(const double &t,const double &tstep);
  void RunSteadyState ();
  void ComputeStressFlows();
  vector<int> SelectColumns(const string &file) const;
  bool   _gridTypeGiven;
  string _extOnly;        ///< commands given that apply only to existing models (:MF6Model, :MF6StartDate, :KeepPackages, :SeepageToRaven; checked with a generated model)
  string _genOnly;        ///< commands given that apply only to generated models (checked with an existing model)  ///< :GridType or :GridFile was written in the .rvg (existing models set the grid type internally)
  vector<double> _nodeTop; ///< [m] top of every active cell, from MODFLOW (existing models, where cells above may be inactive)
  double CellTop(const int l,const int ic) const {
    if (!_nodeTop.empty()){return _nodeTop[l*_ncpl+ic];}
    return (l==0)?_top[ic]:_botm[(l-1)*_ncpl+ic];
  }
  void WriteConnectivity(const double &t,const string &dstr);
  void WriteMF6Files  (const optStruct &Options);
  void ConnectEngine  ();
  void RefreshPointers();
  void CheckProcesses ();
  double WaterTable   (const int ic) const; ///< physical water-table elevation of a column [m] (highest wet cell; column bottom if dry)
  void ValidateInputs () const;
  void OpenOutputs    (const optStruct &Options);
  void WriteSummary   (const optStruct &Options);
  void WriteStepOutputs(const double &tstep,const time_struct &tt);

public:/*-------------------------------------------------------*/
  CGroundwaterModel(CModel *pModel);
  ~CGroundwaterModel();

  //-- parse-time setters (.rvp, .rvg) ----------------------------
  void   AddAquiferClass      (const gw_aquifer_class &c);
  int    GetAquiferClassIndex (const string &name) const;
  void   AddProfile           (const gw_profile &p);
  int    GetProfileIndex      (const string &name) const;
  int    GetNumProfiles       () const {return (int)_aProfiles.size();}
  bool   SetProfileParameter  (const string &profile,const string &param,const double &value);
  static gw_profile DefaultProfile(const string &name);

  void   SetHRUGeometry       (const string &file,const string &idfield){_hruGeomFile=file;_hruIDField=idfield;}
  void   SetLandSurfaceFile   (const string &file){_demFile=file; _genOnly+=" :LandSurface";}
  void   SetBedrockFile       (const string &file){_bedrockFile=file; _genOnly+=" :BedrockSurface";}
  void   SetCellSize          (const double &cs){_cellSize=cs; _genOnly+=" :GridCellSize";}
  void   SetMinCoverage       (const double &f){_minCoverage=f;}
  void   SetLibraryPath       (const string &path){_libPath=path;}
  void   SetSeepageTarget     (const string &sv){_seepTarget=sv;}
  void   SetRiverGeometry     (const string &file,const string &idfield){_riverFile=file;_riverIDField=idfield;}
  void   SetHeadSaveFrequency (const int n){_headSaveFreq=n;}
  void   SetMaxOuterIterations(const int n){_maxiter=n; _genOnly+=" :SolverMaxOuterIterations";}
  void   SetHeadTolerance     (const double &d){_dvclose=d; _genOnly+=" :SolverHeadTolerance";}
  void   SetRiverBedDepth     (const double &d){_rivBedDepth=d;}
  void   AddSurfaceWaterExchange(const string &sv,const string &group,const double &leak){gw_swx s; s.sv=sv; s.group=group; s.leakance=leak; s.iSV=s.kk=-1; _aSWX.push_back(s);}
  void   SetLinkageCache      (const string &f){_cacheFile=f;}
  void   SetWriteNetCDFHeads  (const bool b){_writeNC=b;}
  void   SetWriteBudgetFile   (){_writeCBC=true; _genOnly+=" :WriteGWBudgetFile";}
  void   SetReservoirExchange (const bool b){_resExchange=b;}
  void   SetHotstartRechargeStore(const long long hru,const double &v){_hotRchStore[hru]=v;}
  void   SetHotstartReservoirSeepage(const long long sbid,const double &v){_hotResSeep[sbid]=v;}
  void   SetHotstartReservoirFraction(const long long sbid,const double &v){_hotResFrac[sbid]=v;}
  void   SetHotstartRiverLossCarry(const long long sbid,const double &v){_hotRivCarry[sbid]=v;}
  void   SetHotstartStartHeads(const vector<double> &x){_hotXStart=x;}
  void   SetHotstartRiverTake (const long long sbid,const double &vin,const double &vout,const double &qrec){_hotRivTake[sbid]=make_pair(vin,vout); if (qrec==qrec){_hotQoutRec[sbid]=qrec;}}
  void   SetHotstartSoilZoneColumn(const long long hru,const int cell,const vector<double> &phi){if (cell<0){return;} _hotSzCols[hru][cell]=phi;}
  void   Reinitialize         (const optStruct &Options); ///< fresh start for a new ensemble member
  void   AddWell              (const gw_well &w){_aWells.push_back(w);}
  void   AddObservationWell   (const gw_obswell &o){_aObsWells.push_back(o);}
  void   AddBoundary          (const gw_boundary &b){_aBounds.push_back(b);}
  void   SetWellFlowRateSeries    (const long long id,CTimeSeries *pTS); ///< (also after Initialize: an ensemble member's extra .rvt)
  void   SetBoundaryHeadSeries(const string &name,CTimeSeries *pTS);
  int    SeepageTargetSV      (const int k) const; ///< state variable receiving HRU k's seepage (DOESNT_EXIST: HRU not coupled); usable before Initialize
  bool   SWXStoreSV           (const int k,const int i) const; ///< state variable i of HRU k is the store of a :SurfaceWaterExchange; usable before Initialize
  void   SetSteadyStateInit   (const double &rch_mm_d){_steadyInit=true;_ssRecharge=rch_mm_d;}
  void   SetInitialHeads      (const string &mode,const double &value){_initMode=mode;_initValue=value;}
  void   ClearHotstartData    (); ///< before an ensemble member's .rvc is read: forget the previous member's values
  void   SetGridType          (const int t)            {_gridType=t; _gridTypeGiven=true;}
  void   SetMF6Simulation     (const string &f)        {_extSim=f;}
  void   SetMF6Model          (const string &m)        {_extModel=m; _extOnly+=" :MF6Model";}
  void   SetMF6GridFile       (const string &f,const string &id) {_extGridFile=f; _extGridID=id;}
  void   SetMF6StartDate      (const string &d)        {_extStart=d; _extOnly+=" :MF6StartDate";}
  bool   SetMF6CRS            (const string &spec,string &err) {return _extCRS.Parse(spec,err);}
  bool   SetInputCRS          (const string &spec,string &err) {return _inCRS.Parse(spec,err);}
  const CGWCRS *InCRS         () const {return _inCRS.IsSet()?&_inCRS:NULL;}
  void   ResolveInputCRS      ();   ///< wells given as X,Y: to lon/lat (after the .rvg file is read)
  void   CheckInputLonLat     (const string &file,const double *box,const string &what) const; ///< :InputCRS: stop on a file in longitude/latitude
  void   CheckDryStart        (const optStruct &Options);   ///< warns of layers that start completely dry
  void   AddTakeOver          (const string &p)        {string u=p; for (size_t i=0;i<u.size();i++){u[i]=(char)toupper((unsigned char)u[i]);} _takeOver.push_back(u);}
  void   AddKeepPackage       (const string &p)        {string u=p; for (size_t i=0;i<u.size();i++){u[i]=(char)toupper((unsigned char)u[i]);} _keepPkgs.push_back(u); if (_extOnly.find(" :KeepPackages")==string::npos){_extOnly+=" :KeepPackages";}}
  void   SetStreamPackage     (const string &p,const string &map,const string &aux) {_streamPkg=p; _streamMap=map; _streamAux=aux;}
  void   SetExtSeepage        (const bool on)          {_extSeepage=on; _extOnly+=" :SeepageToRaven";}
  void   AddOverlapWeight     (const long long hru,const int cell,const double w) {gw_ow o; o.hru=hru; o.cell=cell; o.w=w; _owTable.push_back(o);}
  bool   IsExternal           () const                 {return _extSim!="";}
  bool   SZChdColumns         () const {return (_szcLink.size()>0) || IsExternal() || !_chdNode.empty();} ///< soil-zone exchange with specified heads is possible (its output columns)
  double MF6CellArea          (const int ic) const     {return IsExternal()?_extArea[ic]:_grid.area[ic];} ///< [m2] as MODFLOW measures it
  void   SetGridRotation      (const double a)         {_gridAngle=a;}
  void   SetQuadtreeMaxLevel  (const int n)            {_qtMaxLevel=n;}
  void   SetGridFile          (const string &f,const string &id) {_gridFile=f; _gridFileID=id; _gridTypeGiven=true;}
  void   AddGridRefinement    (const string &kind,const string &file,const double size){_refKind.push_back(kind); _refFile.push_back(file); _refSize.push_back(size);}
  void   SetFlowCorrection    (const string &s)        {_flowCorr=s;}
  void   SetLayerConnection   (const bool byElevation) {_connectByElevation=byElevation;}
  void   SetSoilZone          (const bool physical,const int nlayers){_szMode=physical?1:0; _szNLayers=nlayers;   //the last :SoilZone counts
                                  size_t p; while ((p=_genOnly.find(" :SoilZone"))!=string::npos){_genOnly.erase(p,10);}}  //(available with existing models too)
  void   SetSoilZoneIteration (const int maxit,const double tolh,const double tolq){_szMaxIter=maxit;_szTolH=tolh;_szTolQ=tolq;}
  void   SetSoilZoneSubgrid   (const int n){_szSubgrid=n;}
  void   SetSoilZoneOutput    (const string &group){_szOutGroup=group;}
  void   SetSoilZoneVadose    (const int mode,const double acc,const int n,const double par=0.0){_vzMode=mode; if (acc>0){_vzAcc=acc;} if (n>0){_vzN=n;} _vzPar=par;}
  void   SetSoilZoneVadoseClasses(const int K){_vzClasses=K;}
  void   SetSoilZoneVadoseInitialFlux(const double q){_vzQ0=q;}
  void   SetSoilZoneCapillaryRise(const int on){_szCapRise=on;}
  void   SetSoilZoneTiming    (const int t){_szTiming=t;}
  void   SetSoilZoneClasses   (const int k){_szClasses=k;}
  void   SetSoilZoneAquiferClass(const string &c){_szExtClass=c;}
  void   SetHotstartSoilZoneLag(const long long hru,const int cell,const double q,const double debt){_hotLag[hru][cell]=q; _hotDebt[hru][cell]=debt;}
  void   SetHotstartVadose    (const long long hru,const int cell,const vector<double> &v){if (cell<=0){return;} _hotVz[hru][cell]=v;}
  void   SetHotstartVadoseClass(const long long hru,const int cell,const int cls,const vector<double> &v){if ((cell<=0) || (cls<=0)){return;} _hotVzc[hru][cell][cls]=v;}
  void   SetHotstartVadoseTop (const long long hru,const int cell,const int cls,const vector<double> &v){if ((cell<=0) || (cls<=0)){return;} _hotVzt[hru][cell][cls]=v;}
  void   SetHotstartVadoseMethod(const string &m){_hotVzMethod=m;}
  bool   UsesPhysicalSoilZone () const {return _szMode==1;}
  void   WarnIgnoredSoilZoneOptions(const string &given,const bool noisy) const;
  void   NoteSoilClassWithoutKsat(const void *soil,const string &tag){_szNoKsat[soil]=tag;}
  void   SetHotstartHeads     (const int nlay,const int nrow,const int ncol,const vector<double> &h)
                              {_hotDims[0]=nlay;_hotDims[1]=nrow;_hotDims[2]=ncol;_hotHeads=h;}

  //-- run-time ------------------------------------------------------
  void   Initialize           (const optStruct &Options);
  void   Exchange             (double **aPhinew,const double &tstep,const optStruct &Options,const time_struct &tt);
  void   SZStepBegin          (double **aPhinew,double **aPhi,const double &tstep,const optStruct &Options); ///< before Raven's processes (:SoilZone PHYSICAL)
  int    SZClassCount         (const int k) const {return ((k<(int)_szCls.size()) && !_szCls[k].empty())?(int)_szClsW[k].size():1;}
  double SZClassWeight        (const int k,const int c) const {return _szClsW[k][c];}
  void   SZClassSoil          (const int k,const int c,double *aPhi) const;   ///< class c's soil water into a state vector
  void   SZClassResult        (const int k,const int c,const double *aPhi);   ///< class c's state after Raven's processes
  void   ApplyRiverExchange   (double *aRouted,const double &tstep);
  double TakeRiverLossFromInflow(const int p,const double &Qin,const double &tstep); ///< returns reduced upstream inflow [m3/s]
  void   BeginReachRouting    (const int p,const double &tstep);                   ///< the routing state before the last step's outflow take
  void   TakeRiverLossFromOutflow(const int p,double &Qout,const double &tstep);    ///< reduces the reach's routed outflow [m3/s]
  void   EndReachRouting      (const int p);                                       ///< records the reach's stored outflow after Raven's UpdateOutflows
  void   FinishBudgetLine     ();
  ///  [m3] net water added to Raven this step, as Raven's balance books it: a river loss taken from a reach's inflow or
  ///  outflow enters Raven's channel storage and outflow (trapezoidal in time) half in its step and half in the next
  double GetStepReturnVolume  () const {return _returnVol+_riverAppliedVol+_swxAppliedVol+_szAppliedVol+0.5*(_rivTake-_rivTakePrev);}
  void   CloseOutputs         ();
  bool   IsActive             () const {return _active;}
  double GetObservationWellHead(const long long id) const; ///< [m]; RAV_BLANK_DATA if unknown
  void   WriteHotstart        (ofstream &RVC) const;
};

#endif
