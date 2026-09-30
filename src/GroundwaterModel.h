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
#include <fstream>

class CModel;
class CTimeSeries;

struct gw_well       ///< pumping/injection well
{
  long long id; string name;
  double lat,lon,screen_top,screen_bot;  ///< [deg],[m]
  CTimeSeries *pRate;                    ///< [m3/d], negative = withdrawal (from .rvt :WellRate)
};
struct gw_obswell    ///< monitoring well
{
  long long id; string name;
  double lat,lon,z;                      ///< screen elevation [m]
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
  GWL_AQUITARD,         ///< explicit low-K leaky layer (also CONFINING_LAYER)
  GWL_AQUICLUDE         ///< no-flow barrier (IDOMAIN=0)
};

struct gw_aquifer_class
{
  string name;
  double Kh;        ///< horizontal hydraulic conductivity [m/d]
  double Kv;        ///< vertical hydraulic conductivity [m/d]
  double Ss;        ///< specific storage [1/m]
  double Sy;        ///< specific yield [-]
  double porosity;  ///< total porosity [-] (transport, reporting)
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
  double riverbed_K;            ///< riverbed hydraulic conductivity [m/d] (Phase 2 RIV)
  double riverbed_thick;        ///< riverbed thickness [m] (Phase 2 RIV)
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
  //-- Phase 3 stresses, observations, initial conditions --
  vector<gw_well>        _aWells;
  vector<gw_obswell>     _aObsWells;
  vector<gw_boundary>    _aBounds;
  map<long long,CTimeSeries*> _wellTS;   ///< well rates parsed from .rvt, by well ID
  map<string,CTimeSeries*>    _headTS;   ///< boundary heads parsed from .rvt, by boundary name
  vector<int>            _welNode,_welWell; vector<double> _welFrac;
  vector<int>            _ghbNode,_ghbB,_chdNode,_chdB;
  double                *_welQ,*_welHcof,*_welRhs,*_ghbHead,*_ghbCond,*_ghbHcof,*_ghbRhs,*_chdHead,*_chdSim;
  int                    _welNcol,_ghbNcol,_chdNcol;
  double                 _Qwel_spec,_Qwel,_Qghb,_Qchd; ///< [m3/d]
  bool                   _steadyInit; double _ssRecharge; ///< steady-state first period requested (.rvg), recharge [mm/d]
  bool                   _doSteady;  ///< steady-state period in this run (off on hotstart; decided per run/member)
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
  vector<double>         _work;      ///< work array per column
  //-- reservoir <-> aquifer (Raven reservoirs with :SeepageParameters on coupled lake HRUs)
  bool                   _resExchange;
  vector<int>            _resP;      ///< per coupled reservoir: subbasin index
  vector<int>            _resNode,_resR; vector<double> _resW; ///< per RES entry: cell, reservoir, weight
  double                *_resQ,*_resHcof,*_resRhs; int _resNcol;
  double                 _Qres_spec,_Qres,_resPendingVol; ///< [m3/d] seepage sent, delivered; [m3] seepage awaiting delivery
  double                 _resShortVol;                     ///< [m3] reservoir gain MODFLOW could not supply (drying cell), taken from the reach
  vector<gw_swx>         _aSWX;      ///< surface-water exchanges
  vector<int>            _swxCol,_swxX; ///< per SWX entry: column, exchange index
  vector<vector<int> >   _swxLinks;  ///< per SWX entry: links (HRU,cell) of group HRUs in the column
  vector<vector<double> > _swxVl;    ///< per SWX entry and link: water available at the start of the step [m3]
  double                *_swxHead,*_swxCond,*_swxRbot,*_swxHcof,*_swxRhs; int _swxNcol;
  double                 _Qswx_in,_Qswx_out; ///< [m3/d] store leakage into aquifer; aquifer discharge into stores
  double                 _swxAppliedVol;     ///< [m3] net water added to Raven stores by the exchange this step
  string                 _cacheFile; bool _cacheHit; double _geomSeconds;
  bool                   _writeNC; int _ncid,_ncTime,_ncHead,_ncRec;
  void                   OpenNetCDF(const optStruct &Options);
  void                   WriteNetCDF(const double &t);
  vector<char>           _colCHDTop; ///< per column: top active cell is a specified-head (CHD) cell
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
  ofstream               _HEADS;
  vector<double>         _sbRivLen;      ///< per subbasin: total river length on active cells [m]
  vector<double>         _sbLossCarry;   ///< per subbasin: river loss not yet supplied, taken in the next step [m3]
  map<long long,double>  _hotRivCarry;   ///< from .rvc :GWRiverLossCarry, by subbasin ID [m3]
  vector<double>         _swxQ;          ///< per SWX entry: flow into the aquifer this step [m3/d]
  vector<double>         _sbLossRemain;  ///< per subbasin: river loss not yet taken from Raven this step [m3]
  bool                   _budgetLineOpen;
  int    _iGW,_iReturnSV;

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
  void SetBoundaryArrays(const double &t,const double &tstep);
  void RunSteadyState ();
  void ComputeStressFlows();
  vector<int> SelectColumns(const string &file) const;
  bool   _gridTypeGiven;
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
  void   SetReservoirExchange (const bool b){_resExchange=b;}
  void   SetHotstartRechargeStore(const long long hru,const double &v){_hotRchStore[hru]=v;}
  void   SetHotstartReservoirSeepage(const long long sbid,const double &v){_hotResSeep[sbid]=v;}
  void   SetHotstartRiverLossCarry(const long long sbid,const double &v){_hotRivCarry[sbid]=v;}
  void   Reinitialize         (const optStruct &Options); ///< fresh start for a new ensemble member
  void   AddWell              (const gw_well &w){_aWells.push_back(w);}
  void   AddObservationWell   (const gw_obswell &o){_aObsWells.push_back(o);}
  void   AddBoundary          (const gw_boundary &b){_aBounds.push_back(b);}
  void   SetWellRateSeries    (const long long id,CTimeSeries *pTS){_wellTS[id]=pTS;}
  void   SetBoundaryHeadSeries(const string &name,CTimeSeries *pTS){_headTS[name]=pTS;}
  void   SetSteadyStateInit   (const double &rch_mm_d){_steadyInit=true;_ssRecharge=rch_mm_d;}
  void   SetInitialHeads      (const string &mode,const double &value){_initMode=mode;_initValue=value;}
  void   ClearHotstartData    (); ///< before an ensemble member's .rvc is read: forget the previous member's values
  void   SetGridType          (const int t)            {_gridType=t; _gridTypeGiven=true;}
  void   SetMF6Simulation     (const string &f)        {_extSim=f;}
  void   SetMF6Model          (const string &m)        {_extModel=m;}
  void   SetMF6GridFile       (const string &f,const string &id) {_extGridFile=f; _extGridID=id;}
  void   SetMF6StartDate      (const string &d)        {_extStart=d;}
  bool   SetMF6CRS            (const string &spec,string &err) {return _extCRS.Parse(spec,err);}
  void   AddTakeOver          (const string &p)        {string u=p; for (size_t i=0;i<u.size();i++){u[i]=(char)toupper((unsigned char)u[i]);} _takeOver.push_back(u);}
  void   AddKeepPackage       (const string &p)        {string u=p; for (size_t i=0;i<u.size();i++){u[i]=(char)toupper((unsigned char)u[i]);} _keepPkgs.push_back(u);}
  void   SetStreamPackage     (const string &p,const string &map,const string &aux) {_streamPkg=p; _streamMap=map; _streamAux=aux;}
  void   SetExtSeepage        (const bool on)          {_extSeepage=on;}
  void   AddOverlapWeight     (const long long hru,const int cell,const double w) {gw_ow o; o.hru=hru; o.cell=cell; o.w=w; _owTable.push_back(o);}
  bool   IsExternal           () const                 {return _extSim!="";}
  double MF6CellArea          (const int ic) const     {return IsExternal()?_extArea[ic]:_grid.area[ic];} ///< [m2] as MODFLOW measures it
  void   SetGridRotation      (const double a)         {_gridAngle=a;}
  void   SetQuadtreeMaxLevel  (const int n)            {_qtMaxLevel=n;}
  void   SetGridFile          (const string &f,const string &id) {_gridFile=f; _gridFileID=id; _gridTypeGiven=true;}
  void   AddGridRefinement    (const string &kind,const string &file,const double size){_refKind.push_back(kind); _refFile.push_back(file); _refSize.push_back(size);}
  void   SetFlowCorrection    (const string &s)        {_flowCorr=s;}
  void   SetLayerConnection   (const bool byElevation) {_connectByElevation=byElevation;}
  void   SetHotstartHeads     (const int nlay,const int nrow,const int ncol,const vector<double> &h)
                              {_hotDims[0]=nlay;_hotDims[1]=nrow;_hotDims[2]=ncol;_hotHeads=h;}

  //-- run-time ------------------------------------------------------
  void   Initialize           (const optStruct &Options);
  void   Exchange             (double **aPhinew,const double &tstep,const optStruct &Options,const time_struct &tt);
  void   ApplyRiverExchange   (double *aRouted,const double &tstep);
  double TakeRiverLossFromInflow(const int p,const double &Qin,const double &tstep); ///< returns reduced upstream inflow [m3/s]
  void   FinishBudgetLine     ();
  double GetStepReturnVolume  () const {return _returnVol+_riverAppliedVol+_swxAppliedVol;}  ///< [m3] net water added to Raven this step
  void   CloseOutputs         ();
  bool   IsActive             () const {return _active;}
  double GetObservationWellHead(const long long id) const; ///< [m]; RAV_BLANK_DATA if unknown
  void   WriteHotstart        (ofstream &RVC) const;
};

#endif
