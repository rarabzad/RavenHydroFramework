# Raven–MODFLOW 6 groundwater coupling: user guide

This chapter describes how to add a groundwater system to any Raven model. Raven builds a standard
MODFLOW 6 groundwater-flow model from information in the usual Raven input files, runs it in lock-step
with Raven through the MODFLOW 6 shared library (`libmf6`), and exchanges water every time step.
You never edit MODFLOW files; Raven writes them to `<output>/mf6/` so you can inspect them with FloPy or ModelMuse.

## 1. How it works

- HRUs with a groundwater system are marked by `AQUIFER_PROFILE` in the `.rvh` file.
- HRU polygons (GeoJSON) are intersected with a grid that Raven lays over those HRUs (regular by default; rotated, refined, quadtree, HRU-following mesh or imported cells, section 4.3a).
- Aquifer profiles in the `.rvp` file define the layer stack and hydraulic properties; land surface and
  bedrock come from rasters. One MODFLOW 6 model holds all profiles, so lateral flow, leakage and flow
  through aquitards between different groundwater systems are computed by MODFLOW itself.
- Every time step, after Raven's vertical and lateral processes and before routing:
  1. whatever Raven processes put into `GROUNDWATER` in coupled HRUs becomes MODFLOW recharge;
  2. river stages come from Raven reach depths; well rates and boundary heads from `.rvt` time series;
  3. MODFLOW solves (Newton–Raphson, so cells dry and rewet automatically);
  4. seepage where the water table reaches the model top returns to a Raven store, and aquifer–river
     exchange is added to or taken from the reaches.

Raven owns the soil and vadose zone; MODFLOW owns the saturated zone. In coupled HRUs `GROUNDWATER` is a
pass-through store: its end-of-step value is the recharge sent to MODFLOW, and Raven's water balance
counts it as leaving Raven. Groundwater storage itself is reported by the files in section 5.

## 2. Requirements and building

- The MODFLOW 6 shared library: `libmf6.so` (Linux), `libmf6.dll` (Windows) or `libmf6.dylib` (macOS); 6.5 or later,
  tested with 6.8.1 and 6.6.3, 6.8 or newer recommended. Get it with one command in the source folder:

      python tools/get_mf6.py            # MODFLOW 6.8.1 for this computer, into lib/mf6/

  (`--latest` or `--version 6.x.y` for another release, `--zip file` for a release zip downloaded by hand). Building with
  CMake does the same automatically (option `RAVEN_GET_MF6`, on by default) and copies the library next to Raven.
  The script downloads the release zip from https://github.com/MODFLOW-ORG/modflow6/releases and keeps only the library.
  MODFLOW 6.8 is released for Linux x86_64, Windows 64-bit and Apple Silicon Macs; Intel Macs need an earlier release.
- Build Raven as usual: CMake (`mkdir build; cd build; cmake ..; make`) or `make` in `src/` (both link `-ldl`); on Windows
  `src/Raven.sln` (the project includes the new files).
  Raven builds and runs without MODFLOW; the library is loaded only when `:GroundwaterModel MODFLOW6` is used.
- Raven finds the library by itself: next to its executable, or in `lib/mf6/` of the source folder (one or two folders
  above the executable, as for `src/Raven.exe`, `build/Raven` or `build/Release/Raven.exe`). `:MF6Library` in the `.rvg`
  file or the `RAVEN_MF6_LIB` environment variable override this; the system library path is searched last.
  `GWModelSummary.txt` names the library used and its version.

## 3. Quick start: four edits

1. `.rvh`: fill the `AQUIFER_PROFILE` column (`[NONE]` = no groundwater under this HRU).
2. `.rvp`: add `:AquiferClasses`, `:AquiferProfiles` and (optionally) `:AquiferProfileParameters`.
3. `.rvg`: at minimum `:HRUGeometry`; usually also `:LandSurface`, `:BedrockSurface`, `:RiverGeometry`, `:GridCellSize`.
4. `.rvi`: add `:GroundwaterModel MODFLOW6` and route recharge into `GROUNDWATER` for coupled HRUs, e.g.

```
:GroundwaterModel MODFLOW6
:DefineHRUGroups  AquiferHRUs
:HydrologicProcesses
  ...
  :Percolation PERC_CONSTANT FAST_RESERVOIR SLOW_RESERVOIR
    :-->Conditional HRU_GROUP IS_NOT AquiferHRUs
  :Percolation PERC_CONSTANT FAST_RESERVOIR GROUNDWATER
    :-->Conditional HRU_GROUP IS AquiferHRUs
:EndHydrologicProcesses
```

Raven stops with an error if a process removes water from `GROUNDWATER` in a coupled HRU (MODFLOW owns
that storage) or if no process sends water to `GROUNDWATER`.

## 4. Input reference

The groundwater file is `<model>.rvg` by default; `:rvg_Filename file` in the `.rvi` names another file (relative to the `.rvi`).

### 4.1 `.rvh`

`AQUIFER_PROFILE` (existing HRU column): profile name from `:AquiferProfiles`, or `[NONE]`.

### 4.2 `.rvp`

```
:AquiferClasses
  :Attributes, K_HORIZ, K_VERT, SPEC_STORAGE, SPEC_YIELD, POROSITY
  :Units,      m/d,     m/d,    1/m,          none,       none
  GRAVEL,      50.0,    5.0,    1.0E-5,       0.25,       0.30
:EndAquiferClasses
```

| Attribute | Default | Notes |
| --- | --- | --- |
| K_HORIZ | required | horizontal hydraulic conductivity [m/d] |
| K_VERT | K_HORIZ/10 | vertical hydraulic conductivity [m/d] |
| SPEC_STORAGE | 1e-5 | specific storage [1/m] |
| SPEC_YIELD | 0.1 | specific yield [-] |
| POROSITY | 0.3 | not used by flow; reserved for transport and particle tracking |

```
:AquiferProfiles
  # name, nLayers, {class, thickness [m] | TO_BEDROCK, layer type} x nLayers
  VALLEY_ALLUVIUM, 3, GRAVEL, 15.0, AQUIFER, CLAY_TILL, 5.0, AQUITARD, SAND, TO_BEDROCK, CONFINED_AQUIFER
:EndAquiferProfiles
```

| Layer type | MODFLOW 6 | Behaviour |
| --- | --- | --- |
| AQUIFER (UNCONFINED_AQUIFER) | ICELLTYPE 1, ICONVERT 1 | convertible: Sy and saturated-thickness transmissivity when unconfined |
| CONFINED_AQUIFER | ICELLTYPE 0, ICONVERT 0 | constant transmissivity, Ss only |
| AQUITARD (CONFINING_LAYER) | convertible, low K from class | explicit leaky layer |
| AQUICLUDE | IDOMAIN 0 | no-flow barrier; not allowed as the top layer |

Layers thinner than 0.5 m after clipping at bedrock become vertical pass-through cells. Layer k of one
profile connects laterally to layer k of a neighbouring profile, so list layers in the same stratigraphic
order where profiles meet.

```
:AquiferProfileParameters
  :Attributes,     INITIAL_HEAD_DEPTH, SEEPAGE_LEAKANCE, RECHARGE_DELAY
  :Units,          m,                  1/d,              d
  VALLEY_ALLUVIUM, 3.0,                1.0,              0.0
:EndAquiferProfileParameters
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| INITIAL_HEAD_DEPTH | 5.0 m | initial water-table depth below land surface |
| DEFAULT_BEDROCK_DEPTH | 50.0 m | bedrock depth where no bedrock raster is given |
| SOIL_ZONE_DEPTH | 0.0 m | model top below land surface (0 = land surface) |
| SEEPAGE_LEAKANCE | 1.0 1/d | seepage-face conductance per unit area |
| SEEPAGE_SMOOTHING_DEPTH | 0.5 m | seepage starts this far below the model top and reaches full conductance at it |
| RIVERBED_K | 0.5 m/d | riverbed hydraulic conductivity |
| RIVERBED_THICKNESS | 1.0 m | riverbed thickness |
| EXTINCTION_DEPTH | 0.0 m | water-table ET extinction depth (0 = no water-table ET) |
| RECHARGE_DELAY | 0.0 d | linear-reservoir lag between Raven and the water table (0 = direct) |

### 4.3 `.rvg`

| Command | Default | Meaning |
| --- | --- | --- |
| `:HRUGeometry file {id field}` | required; field `HRU_ID` | HRU polygons (GeoJSON, lon/lat) |
| `:RiverGeometry file {subbasin field}` | none | river lines; without a field, each piece belongs to the subbasin of the dominant HRU in its cell |
| `:LandSurface file` | HRU ELEVATION | DEM, ESRI ASCII grid in lon/lat; averaged over each cell |
| `:BedrockSurface file` | DEFAULT_BEDROCK_DEPTH | bedrock elevation raster, same format |
| `:GridCellSize m` | 1000 | base cell size [m] of the grid |
| `:MinCellCoverage f` | 0.25 | fraction of a cell covered by coupled HRUs needed to activate it |
| `:MF6Library path` | RAVEN_MF6_LIB | MODFLOW 6 shared library |
| `:SeepageReturnTo SV` | lowest soil layer | Raven store receiving groundwater seepage (HRUs without soil layers: SURFACE_WATER) |
| `:RiverBedDepth m` or `REFERENCE_FLOW` | REFERENCE_FLOW | channel bed depth below the DEM along rivers |
| `:GWInitialization STEADY_STATE r` | off | steady-state first stress period with uniform recharge r [mm/d] If the steady state does not converge (e.g. wells pumping more than the aquifer supplies), the run starts from the initial heads, with a warning. |
| `:SurfaceWaterExchange SV group LEAKANCE L` | none | Raven store (e.g. DEPRESSION, LAKE_STORAGE) of an HRU group exchanging with the aquifer; leakance L [1/d]; never draws more water than the store holds |
| `:ReservoirExchange` | off | couple Raven reservoirs that have :SeepageParameters and a lake HRU (:HRUID) with an AQUIFER_PROFILE: the water table under the lake (never below the lake bed: a lake above a deep water table is disconnected) sets the reservoir's local groundwater head, and its seepage enters the aquifer the next step; stages must be absolute (:AbsoluteCrestHeight) |
| `:WriteNetCDFHeads` | off | full head field to GWHeads.nc every :HeadSaveFrequency steps (Raven built with -Dnetcdf) |
| `:LinkageCache file` | <output>/mf6/linkage.cache | cached HRU/grid intersection reused by later runs |
| `:SolverMaxOuterIterations n` | 500 | Newton outer iterations per step |
| `:SolverHeadTolerance m` | 0.00001 | Newton head-change tolerance (target); after 100 iterations MODFLOW continues with 10x this tolerance for the rest of that step (GWBudget `converged`=2) |
| `:HeadSaveFrequency n` | 30 | steps between saves of the full head field (`mf6/gwf.hds`) |
| `:Wells` table | none | NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT [m]; rates from `.rvt` |
| `:ObservationWells` table | none | NAME, LATITUDE, LONGITUDE, SCREEN_ELEV [m] |
| `:GeneralHeadBoundary name file HEAD h CONDUCTANCE c {LAYERS ALL\|TOP}` | none | head-dependent boundary; c per cell [m2/d] |
| `:SpecifiedHeadBoundary name file HEAD h {LAYERS ALL\|TOP}` | none | fixed-head boundary |

Boundary files are GeoJSON: polygons select cells whose centre lies inside; lines select cells they cross.
A boundary is placed only in layers whose bottom lies below its head. Recharge falling on specified-head
cells is spread over the HRU's other cells.

```
:Wells
  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT
  1, PW_1, 59.2731, -124.7346, 1199.7, 1188.7
:EndWells
:ObservationWells
  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV
  101, OW_NEAR, 59.2731, -124.7306, 1192.7
:EndObservationWells
```

Well pumping is split over the screened layers by transmissivity; MODFLOW's AUTO_FLOW_REDUCE lowers
pumping as a cell dries (GWBudget.csv reports specified and actual rates).

### 4.3a Grid options (`.rvg`)

| Command | Meaning (default) |
|---|---|
| `:GridType REGULAR\|QUADTREE\|HRU_MESH\|FILE` | kind of grid (`REGULAR`) |
| `:GridRotation deg` | rotation of regular and quadtree grids, counterclockwise (0) |
| `:GridRefinement RIVERS s` / `WELLS s` | cells of size `s` [m] along rivers / around wells and monitoring wells |
| `:GridRefinement LAKES 250` | finer cells in and around every reservoir's lake HRU (outline from `:HRUGeometry`) |
| `:GridRefinement LINES file s` / `POLYGONS file s` | cells of size `s` along lines / in and around polygons (GeoJSON, lon/lat) |
| `:GridFile cells.geojson ID` | imported cells (one single-ring polygon per cell), with `:GridType FILE` |
| `:QuadtreeMaxLevel n` | maximum halvings of `:GridCellSize` (from the refinement sizes) |
| `:FlowCorrection AUTO\|XT3D\|XT3D_RHS\|NONE` | XT3D conductance tensor; `AUTO` = off (XT3D can slow convergence where cells dry and rewet) |
| `:LayerConnection INDEX\|ELEVATION` | connect neighbouring layers by number or by overlapping elevation (DISU) |

MODFLOW formats: regular grids are written as DIS, quadtrees as DISV, HRU meshes, imported cells and
`:LayerConnection ELEVATION` as DISU. Every run writes `mf6/grid_cells.geojson` (all cells, lon/lat) for mapping.
Rivers and boundary lines are clipped exactly by every cell; a boundary drawn along the edge of the active area is
attached to the nearest active cell. `RAVEN_GW_PROFILE=1` prints the time of each model-building stage.

### 4.3b Coupling an existing MODFLOW 6 model (`.rvg`)

Raven links its coupled HRUs (those with an `AQUIFER_PROFILE`) to the top active layer of the user's model and leaves the
model intact (it works on a copy in `<output>/mf6/`).

| Command | Meaning |
|---|---|
| `:MF6Simulation model/mfsim.nam` | the existing simulation (DIS or DISV) |
| `:MF6Model name` | the GWF model to couple, if the simulation holds several |
| `:MF6StartDate yyyy-mm-dd` | calendar date of the model's time zero (Raven may start later: restarts) |
| `:MF6CRS EPSG:32610` | projection of the model's coordinates; Raven reads the cell outlines from the model. Built in: UTM (EPSG 326xx/327xx, NAD83 269xx), 3005, ESRI:102001, 5070, 3978, 3347, or `TMERC`/`ALBERS`/`LCC` parameters |
| `:MF6GridFile cells.geojson CELL_ID` | alternative: top-layer cell polygons (lon/lat), checked against MODFLOW's cell areas |
| `:OverlapWeights ... :EndOverlapWeights` | alternative: `[HRU ID] [cell] [A_overlap/A_cell]` rows (MODFLOW-USG coupling format) |
| `:RavenTakesOver p ...` / `:KeepPackages p ...` | every recharge, ET or UZF package must be declared one way or the other |
| `:StreamPackage p {AUX name / DOMINANT_HRU}` | RIV/DRN/GHB package whose flows go to Raven's reaches (subbasin from an auxiliary variable, default SUBBASIN) |
| `:SeepageToRaven ON/OFF` | DRN_RAVEN at the model top (default ON) |

Each HRU's recharge is split by its relative overlaps, so its volume reaches MODFLOW exactly whatever the delineation
quality. Units (feet, seconds, ...) are converted automatically. The model's own solver settings are used; allow at
least 500 outer iterations (OUTER_MAXIMUM), as the first coupled step (often a steady-state period) is the hardest. `GWBudget.csv` gains a column with the net flow of the
model's own packages (including the water they pass to the model's water mover). See the manual, chapter "Coupling an
existing MODFLOW 6 model".

- **Water mover (MVR):** movers from or to a package Raven takes over (for example a UZF package whose rejected
  infiltration goes to SFR) are removed from the working copy and counted in `GWModelSummary.txt`; a stream package that
  sends water through the mover is refused (its water would be counted twice).
- **Steady-state models** (no storage package, one stress period of any length): the period is stretched over Raven's
  run, so every Raven step is a steady-state solution; the storage column stays zero.
- Tested on USGS example models (Sagehen Creek with UZF, SFR and mover; a steady-state stream-capture model) with MODFLOW
  6.8.1 and 6.6.3: `tests/existing/usgs_sagehen.py`, `usgs_capture.py`.

### 4.4 `.rvt`

```
:WellRate 1 m3/d                 # well ID; negative = withdrawal
  1985-10-01 00:00:00 1.0 1826
  0.0
  ...
:EndWellRate
:BoundaryHead GHB1 m             # overrides HEAD of the named boundary
  ...
:EndBoundaryHead
:ObservationData GW_HEAD 101 m   # observation-well ID; scored in Diagnostics.csv
  ...
:EndObservationData
```

Observed heads at date D are compared with the simulated head at D (end of the previous step).

### 4.5 `.rvc`

```
:InitialGWHeads DEPTH_BELOW_SURFACE 5.0   # or FROM_PROFILE (default), or ELEVATION [m]
```

Raven's `solution.rvc` also contains `:GWHeads nlay nrow ncol` (the head of every cell), and when needed
`:GWRechargeStore` (recharge-delay stores), `:GWReservoirSeepage` (reservoir seepage in transit) and
`:GWRiverLossCarry` (river loss still to be taken), each closed by its `:End...` line. Using `solution.rvc` as the `.rvc`
of the next run continues exactly where the previous run stopped.

## 5. Outputs

| File | Content |
| --- | --- |
| `GWBudget.csv` | per step: Raven recharge, MODFLOW recharge, seepage, storage release, balance error, river leakage and discharge, river exchange applied to Raven, river loss carried to next step, water-table ET, wells (specified, actual), GHB, CHD, recharge-delay storage, recharge onto CHD cells |
| `GWHRUState.csv` | per coupled HRU: mean water-table head and depth below land surface |
| `GWHeads.csv` | head at each observation well, starting at t = 0 (a monitoring well whose screened cell is dry reports the missing value -1.2345) |
| `GWConnectivity.csv` | saturated fraction of the top layer per profile; number of connected saturated areas |
| `GWHeads.nc` | head(time, layer, y, x) with latitude/longitude (optional, :WriteNetCDFHeads) |
| `GWModelSummary.txt` | grid, active cells, profiles, river/well/boundary cells, HRU-to-grid linkage per HRU |
| `mf6/` | the generated MODFLOW 6 model, list file (`mfsim.lst`) and head file (`gwf.hds`) |

All rows are stamped with the date at the end of the step (with the time of day, e.g. `1985-10-01 12:00`, for sub-daily time steps). Groundwater discharge returned to Raven is
counted as input to Raven's water balance (`WatershedStorage.csv`), recharge as output.

## 6. Calibration with Ostrich

Aquifer class and profile parameters live in the `.rvp` file, so Ostrich templates work as for any Raven
parameter. Example template lines (`Liard.rvp.tpl`):

```
  GRAVEL,      par_K_gravel, 5.0, 1.0E-5, par_Sy_gravel, 0.30
  UPLAND_TILL, par_K_till,   0.02, 1.0E-4, 0.08,          0.30
```

Head observations are scored like hydrographs (`GW_HEAD` rows in `Diagnostics.csv`), so they can enter
an Ostrich objective alongside streamflow. An example `ostIn.txt` is provided with the calibration test.
Calibrate hydraulic conductivity in log space.

## 7. Ensembles

Raven's Monte Carlo, DDS and EnKF modes re-initialize the model for every member; the groundwater model is
rebuilt and MODFLOW 6 restarted each time, with its files in each member's output folder.

## 8. Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| "cannot load MODFLOW 6 library" | run `python tools/get_mf6.py` in the source folder, or set `:MF6Library` or `RAVEN_MF6_LIB` to the full path of libmf6; the message lists the places searched |
| Raven stops with no Raven message during initialization | MODFLOW rejected its input: read `<output>/mf6/mfsim.lst` |
| "HRU ... has an AQUIFER_PROFILE but no polygon" | HRU IDs in the GeoJSON do not match the `.rvh`; check the ID field name |
| "does not overlap any active groundwater cell" | reduce `:GridCellSize` or `:MinCellCoverage`, or check bedrock depth |
| Warnings "polygon area ... differs from .rvh AREA" | informational; fluxes are mapped mass-conservatively regardless |
| "MODFLOW 6 did not converge" | start from `:GWInitialization STEADY_STATE` or realistic `:InitialGWHeads`; raise `:SolverMaxOuterIterations`; increase `SEEPAGE_SMOOTHING_DEPTH` |
| "process ... sends water to GROUNDWATER in HRU ... which has no AQUIFER_PROFILE" | that water leaves the model; restrict the process with :-->Conditional |
| well rate is zero ("has no :WellRate time series") | the `:WellRate` block is missing or its well ID does not match an ID in `:Wells` |
| "running more time steps than the MODFLOW 6 simulation holds" | `:Duration` or `:TimeStep` changed after the model files were built; rerun |

## 9. Limitations

Refinement of a regular grid spans whole rows and columns; XT3D (optional, off by default) costs more per iteration and can slow convergence;
rasters in latitude/longitude (ESRI ASCII);
river stage from the start of the step; reservoirs exchange with a one-step lag (:ReservoirExchange); groundwater transport, particle tracking, demand-driven wells and frozen-ground recharge control
are planned extensions. (Coupling an existing MODFLOW 6 model is available: section 4.3b and `examples/Liard_existing_model/`.)
