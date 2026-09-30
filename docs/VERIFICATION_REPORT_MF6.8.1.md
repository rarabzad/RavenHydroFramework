# Verification report (Raven 4.15 with the MODFLOW 6 coupling)

- Raven: `src/Raven.exe` (branch mf6-groundwater-coupling, on Raven 4.15 v593; GNU make build)
- MODFLOW 6 library: `libmf6.so (MODFLOW 6.8.1, Linux)`
- scratch folder: `<scratch>`
- date: 2026-09-29
- mode: full

The suite ran in one pass; the last step (`ext_fuzz`) was interrupted by a restart of the machine and was repeated,
together with `edge`, after a rebuild that only changed one error message and restored `Raven_BMI.cpp` to 4.15.
With MODFLOW 6.6.3, `geomtest` first stopped while deleting a temporary file that the parallel 6.8.1 run had already removed
(all its checks had passed); the cleanup now tolerates that, and the step was repeated.

| Step | Result | Time [s] | Summary |
|---|---|---|---|
| gridtest | PASS | 7 | ALL PASSED: 0 failures |
| geomtest | PASS | 3 | ALL PASSED |
| crs | PASS | 3 | 1 worst position difference [m], largest 3.0e-06 (limit 1e-03) |
| mf6files | PASS | 119 | 16 cases, 127 checks, 0 failed |
| edge_nolib | PASS | 21 | 53 of 53 |
| edge | PASS | 254 | 53 of 53 |
| thiem | PASS | 29 | 16 grids fitted (1 with no cells in the window), 0 runs with failed steps |
| restart | PASS | 52 | 5 restart deviations, largest 9.9e-07 (limit 1e-04) |
| lake | PASS | 14 | no FAIL |
| liard_grids | PASS | 741 | no FAILED results in <scratch>/liard_grids/results.json |
| behaviour | PASS | 175 | behaviour checks: ALL PASSED |
| budget_xcheck | PASS | 1 | 1 worst differences, largest 1.1e-07 (limit 1e-06) |
| fuzz_grids | PASS | 683 | 72 of 72 set-ups pass, 20 with a restart (results_11.json, results_21.json, results_22.json) |
| fuzz_lakes | PASS | 179 | 24 of 24 set-ups pass, 6 with a restart (results_12_lakes.json, results_13_lakes.json) |
| ext_guards | PASS | 104 | 14 of 14 |
| ext_roundtrip | PASS | 54 | 5 round-trip differences, largest 1.7e-06 (limit 1e-04) |
| ext_restart | PASS | 310 | 2 restart deviations, largest 1.4e-07 (limit 1e-04) |
| usgs_sagehen | PASS | 970 | Sagehen checks: ALL PASSED |
| usgs_capture | PASS | 20 | capture checks: ALL PASSED |
| ext_fuzz | PASS | 1001 | 20 of 20 set-ups pass, 6 with a restart (fz_results_5.json) |

All 20 steps pass.


## gridtest (PASS)

```
regular variable: widths 250..1000, neighbour ratio <= 2   ok  (2)
quadtree: cells tile the domain (area sum / domain - 1)    ok  (0)
quadtree: neighbours symmetric (478 links)                 ok  (0)
quadtree: centroid locates its own cell (misses)           ok  (0)
quadtree: 2:1 balance (max size ratio of neighbours)       ok  (2)
  quadtree cells 121 (smallest 250 m), cells with hanging vertices 23
imported (quadtree polygons): cells tile the domain (area sum / domain - 1) ok  (0)
imported (quadtree polygons): neighbours symmetric (478 links) ok  (0)
imported (quadtree polygons): centroid locates its own cell (misses) ok  (0)
imported: same neighbours and face lengths as quadtree     ok  (0)
voronoi: cells tile the domain (area sum / domain - 1)     ok  (-2.22e-16)
voronoi: neighbours symmetric (3406 links)                 ok  (0)
voronoi: centroid locates its own cell (misses)            ok  (0)
voronoi: centres equidistant from faces (orthogonal)       ok  (0)
voronoi: all cells convex                                  ok  (600)
voronoi (lattice seeds, cocircular): cells tile the domain (area sum / domain - 1) ok  (0)
voronoi (lattice seeds, cocircular): neighbours symmetric (4344 links) ok  (0)
voronoi (lattice seeds, cocircular): centroid locates its own cell (misses) ok  (0)
voronoi (lattice): every cell clips to its full area       ok  (0)
voronoi (lattice): squares have 4 vertices                 ok  (4)
clipping: covering polygon returns each cell area          ok  (1.62e-13)
concave cell: L-shape by square = 0.75                     ok  (0.75)
cache round trip (quadtree)                                ok  (0)
cache round trip (regular uniform, bit-exact)              ok  (0)
ALL PASSED: 0 failures
```

## geomtest (PASS)

```
clipping: 400 polygons against shapely (worst relative error) ok  (4.6e-12)
projection: forward against PROJ (worst position error, m)   ok  (6.74e-07)
projection: inverse round trip (worst error, degrees)        ok  (0)
raster: bilinear sampling of a plane (worst error)           ok  (3.55e-15)
ALL PASSED
```

## crs (PASS)

```
EPSG:32610   largest position difference 3.05e-06 m
EPSG:32609   largest position difference 2.95e-06 m
EPSG:32718   largest position difference 2.97e-06 m
EPSG:26910   largest position difference 2.72e-06 m
EPSG:3005    largest position difference 7.44e-07 m
EPSG:5070    largest position difference 7.52e-07 m
EPSG:3978    largest position difference 7.40e-07 m
EPSG:3347    largest position difference 6.91e-07 m
ESRI:102001  largest position difference 7.79e-07 m
worst 3.0e-06 m
```

## mf6files (PASS)

```
regular (default)                            PASS 
rotated 20, rivers 500                       PASS 
quadtree, rivers 500                         PASS 
HRU mesh, rivers 1000                        PASS 
layers by elevation                          PASS 
rotated 20 uniform, boundaries               PASS 
regular, boundaries                          PASS 
HRU mesh, aquiclude                          PASS 
elevation, aquiclude                         PASS 
rotated 30, elevation                        PASS 
thin top layer, smoothing 1 m                PASS 
HRU mesh, wells, monitoring well outside     PASS 
HRU mesh, lakes 500                          PASS 
quadtree, lakes 500                          PASS 
rotated 20, lakes 750                        PASS 
imported cells (FILE)                        PASS
16 cases, 127 checks, 0 failed; results in <scratch>/mf6files/results.json
```

## edge_nolib (PASS)

```
N30_reservoir_absolute_stages            ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
N31_mesh_obswell_outside_refined         ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
N32_reservoir_hru_of_other_subbasin      error     PASS  CModel::Initialize: the :HRUID of reservoir TestLake is an HRU of another subbasin; it must be an HRU of subba
P01_cell1000                             ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
P02_cell2000                             ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
P03_cell4000                             ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
P04_halfday                              ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
E01_two_link_methods                     error     PASS  CGroundwaterModel: with :MF6Simulation give one of :MF6CRS, :MF6GridFile or :OverlapWeights (to link HRUs to t
E02_generated_setting                    error     PASS  CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove  :GridCellS
E03_grid_type                            error     PASS  CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove the grid op
E04_netcdf_with_weights                  error     PASS  CGroundwaterModel: :WriteNetCDFHeads needs the cell outlines (:MF6CRS or :MF6GridFile); an :OverlapWeights tab
E05_start_before_model                   error     PASS  CGroundwaterModel: Raven's :StartDate 1985-09-01 is before the model's time zero (:MF6StartDate 1985-10-01)
E06_two_solutions                        error     PASS  CGroundwaterModel: the MODFLOW 6 simulation has more than one solution (e.g. a GWT or GWE model with its own I
E07_adaptive_time_steps                  error     PASS  CGroundwaterModel: the MODFLOW 6 model uses adaptive time stepping (ATS6); Raven needs fixed time steps. Remov
E08_noleap_calendar                      error     PASS  CGroundwaterModel: an existing MODFLOW 6 model needs a Gregorian :Calendar (dates are compared day by day)
E09_model_in_working_folder              error     PASS  CGroundwaterModel: the MODFLOW 6 simulation out/mf6/mfsim.nam lies in Raven's working folder <scratch>
E10_restart_with_time_series             error     PASS  CGroundwaterModel: WEL_PUMPS uses time series, whose times count from the model's time zero; start Raven at :M
E11_restart_with_tvk                     error     PASS  CGroundwaterModel: NPF uses time-varying properties (TVK6/TVS6), whose period data cannot be renumbered safely
E12_usg_command                          error     PASS  ParseGWFile: :GWRiverConnection belongs to the MODFLOW-USG coupling of Raven 4.15 and earlier, which this vers
E13_period_not_whole_steps               error     PASS  CGroundwaterModel: stress period 1 (1.5 d) is not a whole number of Raven time steps (1 d)
E20_valid                                ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
E21_valid_restart                        ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
E22_valid_sim_file_other_name            ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
E23_valid_output_inside_model_folder     ok        PASS  CGroundwaterModel: cannot load MODFLOW 6 library /nonexistent/libmf6.so: /nonexistent/libmf6.so: cannot open s
53 of 53 cases pass (without the MODFLOW library); results in <scratch>/edge/results.json
```

## edge (PASS)

```
N30_reservoir_absolute_stages            ok        PASS  Successful Simulation
N31_mesh_obswell_outside_refined         ok        PASS  Successful Simulation
N32_reservoir_hru_of_other_subbasin      error     PASS  CModel::Initialize: the :HRUID of reservoir TestLake is an HRU of another subbasin; it must be an HRU of subba
P01_cell1000                             ok        PASS  Successful Simulation
P02_cell2000                             ok        PASS  Successful Simulation
P03_cell4000                             ok        PASS  Successful Simulation
P04_halfday                              ok        PASS  Successful Simulation
E01_two_link_methods                     error     PASS  CGroundwaterModel: with :MF6Simulation give one of :MF6CRS, :MF6GridFile or :OverlapWeights (to link HRUs to t
E02_generated_setting                    error     PASS  CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove  :GridCellS
E03_grid_type                            error     PASS  CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove the grid op
E04_netcdf_with_weights                  error     PASS  CGroundwaterModel: :WriteNetCDFHeads needs the cell outlines (:MF6CRS or :MF6GridFile); an :OverlapWeights tab
E05_start_before_model                   error     PASS  CGroundwaterModel: Raven's :StartDate 1985-09-01 is before the model's time zero (:MF6StartDate 1985-10-01)
E06_two_solutions                        error     PASS  CGroundwaterModel: the MODFLOW 6 simulation has more than one solution (e.g. a GWT or GWE model with its own I
E07_adaptive_time_steps                  error     PASS  CGroundwaterModel: the MODFLOW 6 model uses adaptive time stepping (ATS6); Raven needs fixed time steps. Remov
E08_noleap_calendar                      error     PASS  CGroundwaterModel: an existing MODFLOW 6 model needs a Gregorian :Calendar (dates are compared day by day)
E09_model_in_working_folder              error     PASS  CGroundwaterModel: the MODFLOW 6 simulation out/mf6/mfsim.nam lies in Raven's working folder <scratch>
E10_restart_with_time_series             error     PASS  CGroundwaterModel: WEL_PUMPS uses time series, whose times count from the model's time zero; start Raven at :M
E11_restart_with_tvk                     error     PASS  CGroundwaterModel: NPF uses time-varying properties (TVK6/TVS6), whose period data cannot be renumbered safely
E12_usg_command                          error     PASS  ParseGWFile: :GWRiverConnection belongs to the MODFLOW-USG coupling of Raven 4.15 and earlier, which this vers
E13_period_not_whole_steps               error     PASS  CGroundwaterModel: stress period 1 (1.5 d) is not a whole number of Raven time steps (1 d)
E20_valid                                ok        PASS  Successful Simulation
E21_valid_restart                        ok        PASS  Successful Simulation
E22_valid_sim_file_other_name            ok        PASS  Successful Simulation
E23_valid_output_inside_model_folder     ok        PASS  Successful Simulation
53 of 53 cases pass; results in <scratch>/edge/results.json
```

## thiem (PASS)

```
HRU mesh 500 m                             dh(600-1200) err   +7.80%   dh(1200-2400) err   +0.78%   residual -6.4e-10  failed 0
HRU mesh 250 m                             dh(600-1200) err  +11.81%   dh(1200-2400) err   +0.28%   residual 5.8e-10  failed 0
HRU mesh 250 m without XT3D                dh(600-1200) err  +11.81%   dh(1200-2400) err   +0.28%   residual 5.8e-10  failed 0
HRU mesh 250 m with XT3D                   dh(600-1200) err   +9.90%   dh(1200-2400) err   -0.15%   residual -5.7e-09  failed 0
HRU mesh 250 m lattice only                dh(600-1200) err   +3.05%   dh(1200-2400) err  +17.93%   residual -4.9e-12  failed 0
regular 250 m via DISU                     dh(600-1200) err   +3.07%   dh(1200-2400) err  +17.97%   residual 3.0e-11  failed 0
imported (quadtree cells)                  dh(600-1200) err   +5.94%   dh(1200-2400) err   +0.54%   residual -3.8e-09  failed 0
analytic head difference per doubling of distance: 1.1032 m
hru_mesh_1000_m_125_m_at_wells           slope error  -7.80%   fit RMS 0.0901 m   cells used 39
hru_mesh_250_m                           slope error  +0.29%   fit RMS 0.0071 m   cells used 437
hru_mesh_250_m_lattice_only              slope error  +0.00%   fit RMS 0.0040 m   cells used 393
hru_mesh_250_m_with_xt3d                 slope error  -0.06%   fit RMS 0.0016 m   cells used 437
hru_mesh_250_m_without_xt3d              slope error  +0.29%   fit RMS 0.0071 m   cells used 437
hru_mesh_500_m                           slope error  +0.57%   fit RMS 0.0318 m   cells used 103
imported_quadtree_cells                  slope error  +0.25%   fit RMS 0.0082 m   cells used 484
quadtree_1000_m_to_62_5_m                slope error  +0.80%   fit RMS 0.0121 m   cells used 120
quadtree_250_m_to_62_5_m                 slope error  +0.26%   fit RMS 0.0082 m   cells used 485
quadtree_500_m_to_62_5_m                 slope error  +0.59%   fit RMS 0.0144 m   cells used 174
regular_1000_m                           slope error -100.00%   fit RMS nan m   cells used 0
regular_1000_m_variable_spacing_to_125_m slope error  +1.01%   fit RMS 0.0048 m   cells used 774
regular_125_m                            slope error  +0.01%   fit RMS 0.0021 m   cells used 1745
regular_250_m                            slope error  +0.04%   fit RMS 0.0043 m   cells used 393
regular_250_m_rotated_17_deg             slope error  +0.03%   fit RMS 0.0042 m   cells used 393
regular_250_m_via_disu                   slope error  +0.04%   fit RMS 0.0043 m   cells used 393
regular_500_m                            slope error  +0.60%   fit RMS 0.0098 m   cells used 65
```

## restart (PASS)

```
regular (default)            restart deviation 9.9e-07 (60 days)  residual -1.6e-11  recharge 0e+00  river 9e-12  failed 0
rotated + variable spacing   restart deviation 1.9e-07 (60 days)  residual -1.4e-11  recharge 0e+00  river 7e-12  failed 0
quadtree                     restart deviation 4.4e-07 (60 days)  residual -1.6e-11  recharge 0e+00  river 2e-12  failed 0
HRU mesh                     restart deviation 1.9e-07 (60 days)  residual 4.1e-12  recharge 1e-11  river 4e-11  failed 0
layers by elevation          restart deviation 8.3e-08 (60 days)  residual 1.6e-11  recharge 0e+00  river 2e-11  failed 0
```

## lake (PASS)

```
lake HRU 3795 area 126.9 km2 (rvh), elevation 1880 m
regular, rotated, LAKES 750  cells near the lake:  31, median  750 m | elsewhere:  544, median 1500 m | recharge 0e+00 residual -2.9e-13 river 2e-11 failed 0
      reservoir checks: {'reservoir: sent(t) - awaiting(t-1) (rel)': 1.5282051174084972e-10, 'reservoir: delivered / sent': 1.0, 'reservoir: (delivered + not supplied) / sent': 1.0}
quadtree, LAKES 375          cells near the lake:  81, median  375 m | elsewhere:  280, median 3000 m | recharge 0e+00 residual -1.2e-11 river 1e-11 failed 0
      reservoir checks: {'reservoir: sent(t) - awaiting(t-1) (rel)': 1.5529382274134303e-10, 'reservoir: delivered / sent': 1.0, 'reservoir: (delivered + not supplied) / sent': 1.0}
HRU mesh, LAKES 500          cells near the lake:  45, median  491 m | elsewhere:  547, median 1717 m | recharge 0e+00 residual 3.1e-13 river 3e-12 failed 0
      reservoir checks: {'reservoir: sent(t) - awaiting(t-1) (rel)': 1.5623214009672461e-10, 'reservoir: delivered / sent': 1.0, 'reservoir: (delivered + not supplied) / sent': 1.0}
HRU mesh, POLYGONS 500       cells near the lake:  45, median  491 m | elsewhere:  547, median 1717 m | recharge 0e+00 residual 3.1e-13 river 3e-12 failed 0
      reservoir checks: {'reservoir: sent(t) - awaiting(t-1) (rel)': 1.5623214009672461e-10, 'reservoir: delivered / sent': 1.0, 'reservoir: (delivered + not supplied) / sent': 1.0}
HRU mesh, no refinement      cells near the lake:  14, median 1210 m | elsewhere:  546, median 1717 m | recharge 0e+00 residual -2.4e-13 river 1e-11 failed 0
      reservoir checks: {'reservoir: sent(t) - awaiting(t-1) (rel)': 1.5574037499193725e-10, 'reservoir: delivered / sent': 1.0, 'reservoir: (delivered + not supplied) / sent': 1.0}
```

## liard_grids (PASS)

```
regular 2 km (default) ok 58 s
regular 2 km, rotated 20 deg, 1 km along rivers ok 323 s
quadtree 2 km to 500 m along rivers ok 151 s
HRU mesh 2 km, 1 km along rivers ok 147 s
regular 2 km, layers by elevation ok 60 s
grid                                                 cells   time  NSE 10BE004/005/010      seepage  rivers | recharge  residual failed
regular 2 km (default)                                8840     58  0.728 / 0.861 / 0.635       9.48    0.82 |    0e+00  -2.4e-11 0
regular 2 km, rotated 20 deg, 1 km along rivers      35058    323  0.726 / 0.861 / 0.638       9.92    0.66 |    0e+00  -2.0e-11 0
quadtree 2 km to 500 m along rivers                  16276    151  0.719 / 0.858 / 0.611       9.49    0.01 |    0e+00  -1.6e-11 0
HRU mesh 2 km, 1 km along rivers                     14610    147  0.717 / 0.859 / 0.627       9.94    0.35 |    8e-12  -1.3e-11 0
regular 2 km, layers by elevation                     8840     60  0.666 / 0.859 / 0.606       4.11    0.49 |    0e+00  -1.7e-11 0
results in <scratch>/liard_grids/results.json
```

## behaviour (PASS)

```
reservoir: lake on HRU 3795, 50 km2, 731 days: Successful Simulation
  recharge identity                                  0.00e+00  ok
  MODFLOW cumulative residual                        1.68e-11  ok
  river bookkeeping                                  9.48e-12  ok
  reservoir: sent(t) - awaiting(t-1)                 0.00e+00  ok
  reservoir: (delivered + not supplied) / sent - 1   0.00e+00  ok
  unconverged steps                                  0.00e+00  ok
  Raven balance error [mm]                           5.15e-01  ok
  seepage sent: mean 716777 m3/d, total 5.24e+08 m3; delivered/sent 1.000000000000; wells delivered/specified 1.000
ensemble: 2 Monte Carlo members: Successful Simulation
  largest relative difference of a member from the single run: 0.0e+00  ok
behaviour checks: ALL PASSED
```

## budget_xcheck (PASS)

```
<scratch>/behaviour/reservoir/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.6489e+07/-1.6489e+07 (1.8e-11) | RIV 587887/587887 (1.7e-10) | WEL -20000/-20000 (0.0e+00) | RES 719427/719427 (0.0e+00)
<scratch>/lake/lk_hru_mesh_lakes_500/out RCH 0/0 (0.0e+00) | DRN -1.23035e+06/-1.23035e+06 (1.6e-10) | RIV -613345/-613345 (0.0e+00) | WEL -5000/-5000 (0.0e+00) | RES 636841/636841 (0.0e+00)
<scratch>/lake/lk_hru_mesh_no_refinement/out RCH 0/0 (0.0e+00) | DRN -1.21568e+06/-1.21568e+06 (0.0e+00) | RIV -609819/-609819 (0.0e+00) | WEL -5000/-5000 (0.0e+00) | RES 638837/638837 (0.0e+00)
<scratch>/lake/lk_hru_mesh_polygons_500/out RCH 0/0 (0.0e+00) | DRN -1.23035e+06/-1.23035e+06 (1.6e-10) | RIV -613345/-613345 (0.0e+00) | WEL -5000/-5000 (0.0e+00) | RES 636841/636841 (0.0e+00)
<scratch>/lake/lk_quadtree_lakes_375/out RCH 0/0 (0.0e+00) | DRN -1.12697e+06/-1.12697e+06 (3.5e-10) | RIV -883987/-883987 (4.5e-11) | WEL -353.858/-353.858 (1.1e-07) | RES 640389/640389 (0.0e+00)
<scratch>/lake/lk_regular_rotated_lakes_750/out RCH 0/0 (0.0e+00) | DRN -1.14834e+06/-1.14834e+06 (0.0e+00) | RIV -718586/-718586 (0.0e+00) | WEL -5000/-5000 (0.0e+00) | RES 652184/652184 (0.0e+00)
<scratch>/liard_grids/hru_mesh_2_km_1_km_along_rivers/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.56241e+07/-1.56241e+07 (1.2e-10) | RIV 1.14089e+06/1.14089e+06 (1.8e-10)
<scratch>/liard_grids/quadtree_2_km_to_500_m_along_rivers/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.57703e+07/-1.57703e+07 (2.4e-10) | RIV 1.1825e+06/1.1825e+06 (1.7e-10)
<scratch>/liard_grids/regular_2_km_default/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.57888e+07/-1.57888e+07 (1.1e-10) | RIV 651322/651322 (3.1e-10)
<scratch>/liard_grids/regular_2_km_layers_by_elevation/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.16675e+07/-1.16675e+07 (3.5e-10) | RIV 812275/812275 (6.2e-10)
<scratch>/liard_grids/regular_2_km_rotated_20_deg_1_km_along_rivers/out RCH 1.2974e+07/1.2974e+07 (1.2e-10) | DRN -1.57845e+07/-1.57845e+07 (2.3e-10) | RIV 896537/896537 (4.5e-10)
45 package budgets compared; WORST relative difference coupler vs MODFLOW list-file budget: 1.08e-07
```

## fuzz_grids (PASS)

```
results in <scratch>/fuzz/results_11.json
results in <scratch>/fuzz/results_21.json
results in <scratch>/fuzz/results_22.json
```

## fuzz_lakes (PASS)

```
results in <scratch>/fuzz/results_12_lakes.json
results in <scratch>/fuzz/results_13_lakes.json
```

## ext_guards (PASS)

```
model written: 117 x 98 x 4 layers; active columns 5835 ; river cells 239 ; periods 25 ; units meters_days
model written: 117 x 98 x 4 layers; active columns 5835 ; river cells 239 ; periods 25 ; units feet_seconds
7714 weights written to <scratch>/extmodel/md/weights.rvg
ok   valid_aux                          Successful Simulation
ok   valid_dominant_hru                 Successful Simulation
ok   undeclared_evt                     CGroundwaterModel: the model's EVT6 package EVT_MODEL would add recharge or evapotranspiration beside Raven's. List it in :RavenTakesOver (Raven repla
ok   takeover_of_wells                  CGroundwaterModel: :RavenTakesOver applies to recharge, evapotranspiration and UZF packages; WEL_PUMPS is WEL6
ok   stream_not_stream_type             CGroundwaterModel: :StreamPackage must be a RIV, DRN or GHB package; WEL_PUMPS is WEL6
ok   unknown_package                    CGroundwaterModel: package RIV_NOPE is not in the model name file gwf_liard.nam
ok   start_date_mismatch                CGroundwaterModel: Raven's :StartDate 1985-10-01 is before the model's time zero (:MF6StartDate 1985-10-02)
ok   longer_than_model                  CGroundwaterModel: the MODFLOW 6 model ends before Raven's run (731 days after Raven's start, :Duration 800)
ok   period_not_whole_steps             CGroundwaterModel: stress period 1 (1.5 d) is not a whole number of Raven time steps (1 d)
ok   two_models_unnamed                 CGroundwaterModel: the simulation holds several GWF models; name the one to couple with :MF6Model
ok   grid_in_wrong_units                CGroundwaterModel: the cells of ext_cells.geojson are typically 0.0924295 times MODFLOW's cell area. Check the length units of the model and the proje
ok   generated_option_left              CGroundwaterModel: with :MF6Simulation the model's own grid, wells and boundaries are used; remove the grid options (:GridType, :GridRotation, :GridRe
ok   both_link_methods                  CGroundwaterModel: with :MF6Simulation give one of :MF6CRS, :MF6GridFile or :OverlapWeights (to link HRUs to the model's cells)
ok   no_link_method                     CGroundwaterModel: with :MF6Simulation give one of :MF6CRS, :MF6GridFile or :OverlapWeights (to link HRUs to the model's cells)
14 of 14 behave as required
```

## ext_roundtrip (PASS)

```
existing-model run: Successful Simulation
  Raven recharge [m3/d]      max relative difference 0.0e+00
  MF6 recharge [m3/d]        max relative difference 0.0e+00
  seepage to Raven [m3/d]    max relative difference 5.4e-07
  storage release [m3/d]     max relative difference 8.0e-07
  outlet flow SUB_43              max relative difference 1.7e-06
  audit: recharge 0e+00 residual -3.0e-11 failed 0
```

## ext_restart (PASS)

```
continuous: Successful Simulation
split at day  30 (1985-10-31): 170 days compared, max deviation 1.4e-07; wells/boundaries flow change seen in restarted run: yes; audit recharge 0e+00 residual 6.9e-12 failed 0
   renumbered first period block note: []
split at day 150 (1986-02-28): 50 days compared, max deviation 9.0e-08; wells/boundaries flow change seen in restarted run: yes; audit recharge 0e+00 residual 9.9e-12 failed 0
   renumbered first period block note: []
```

## usgs_sagehen (PASS)

```
      Raven recharge 2.12e+04 m3/d, model packages -5070 m3/d, storage -1.613e+04 m3/d, seepage to Raven 0 m3/d (means)
takeover      expected ok   : Successful Simulation
      recharge identity              0.00e+00 ok
      MODFLOW cumulative residual    2.82e-10 ok
      river bookkeeping              0.00e+00 ok
      unconverged steps              0.00e+00 ok
      seepage returned               0.00e+00 ok
      Raven recharge 2.12e+04 m3/d, model packages -5069 m3/d, storage -1.613e+04 m3/d, seepage to Raven 0 m3/d (means)
seepage       expected ok   : Successful Simulation
      recharge identity              0.00e+00 ok
      MODFLOW cumulative residual    1.27e-11 ok
      river bookkeeping              0.00e+00 ok
      unconverged steps              0.00e+00 ok
      seepage returned               0.00e+00 ok
      Raven recharge 2.12e+04 m3/d, model packages -4599 m3/d, storage -1.612e+04 m3/d, seepage to Raven 481.8 m3/d (means)
stream_mover  expected error: CGroundwaterModel: the stream package DRN-1 sends water through the model's water mover (ex-gwf-sagehen.mvr); Raven would take that water from the aqu
keep_uzf      expected ok   : Successful Simulation
      recharge identity              0.00e+00 ok
      MODFLOW cumulative residual    5.31e-09 ok
      river bookkeeping              0.00e+00 ok
      unconverged steps              1 (flagged; accounts still exact)
      seepage returned               0.00e+00 ok
      Raven recharge 2.12e+04 m3/d, model packages -1.274e+04 m3/d, storage -8460 m3/d, seepage to Raven 0 m3/d (means)
restart       day 200 -> 199 days compared, largest deviation 3.3e-08 ok
Sagehen checks: ALL PASSED
```

## usgs_capture (PASS)

```
published  Successful Simulation
      steady-state model recognised: True
      recharge identity              0.00e+00 ok
      MODFLOW residual               1.71e-07 on the converged days ok
      river bookkeeping              6.45e-11 ok
      unconverged steps              4 (the model's own tolerance of 1e-8 m; flagged, as designed)
      seepage returned               0.00e+00 ok
      storage release (must be 0 in a steady model): 0.0e+00 m3/d; river to Raven -365 m3/d, Raven recharge 1631 m3/d (means)
solver     Successful Simulation
      steady-state model recognised: True
      recharge identity              0.00e+00 ok
      MODFLOW cumulative residual    2.78e-09 ok
      river bookkeeping              1.10e-11 ok
      unconverged steps              0.00e+00 ok
      seepage returned               0.00e+00 ok
      storage release (must be 0 in a steady model): 0.0e+00 m3/d; river to Raven -365 m3/d, Raven recharge 1631 m3/d (means)
restart    day 30 -> 30 days compared, largest deviation 1.0e-07 ok
capture checks: ALL PASSED
```

## ext_fuzz (PASS)

```
f00 PASS {'units': 'fs', 'link': 'file', 'stream': 'dominant', 'seep': True, 'evt': 'keep', 'cov': 0.5, 'days': 45}  
f01 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.1, 'days': 90}  
f02 PASS {'units': 'fs', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 45}  
f03 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'keep', 'cov': 0.25, 'days': 45} (5, 6.906379219003669e-08) 
f04 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': False, 'evt': 'takeover', 'cov': 0.25, 'days': 45}  
f05 PASS {'units': 'fs', 'link': 'crs', 'stream': 'aux', 'seep': False, 'evt': 'keep', 'cov': 0.1, 'days': 90} (51, 5.523285604327903e-08) 
f06 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 45}  
f07 PASS {'units': 'fs', 'link': 'file', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 90} (65, 9.561063236612013e-08) 
f08 PASS {'units': 'fs', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.5, 'days': 90}  
f09 PASS {'units': 'md', 'link': 'weights', 'stream': 'dominant', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 45}  
f10 PASS {'units': 'fs', 'link': 'crs', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 90}  
f11 PASS {'units': 'md', 'link': 'weights', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.25, 'days': 45}  
f12 PASS {'units': 'md', 'link': 'crs', 'stream': 'dominant', 'seep': True, 'evt': 'takeover', 'cov': 0.1, 'days': 90} (74, 2.3515050308229927e-07) 
f13 PASS {'units': 'md', 'link': 'file', 'stream': 'dominant', 'seep': True, 'evt': 'takeover', 'cov': 0.5, 'days': 45}  
f14 PASS {'units': 'md', 'link': 'weights', 'stream': 'aux', 'seep': True, 'evt': 'takeover', 'cov': 0.5, 'days': 90} (81, 2.687272397099456e-07) 
f15 PASS {'units': 'md', 'link': 'file', 'stream': 'none', 'seep': True, 'evt': 'takeover', 'cov': 0.1, 'days': 45}  
f16 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': False, 'evt': 'keep', 'cov': 0.1, 'days': 90} (60, 3.6357161203359714e-08) 
f17 PASS {'units': 'md', 'link': 'crs', 'stream': 'aux', 'seep': False, 'evt': 'takeover', 'cov': 0.1, 'days': 90}  
f18 PASS {'units': 'fs', 'link': 'crs', 'stream': 'none', 'seep': False, 'evt': 'keep', 'cov': 0.5, 'days': 90}  
f19 PASS {'units': 'md', 'link': 'weights', 'stream': 'none', 'seep': False, 'evt': 'takeover', 'cov': 0.25, 'days': 45}
```
