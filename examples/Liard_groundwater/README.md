# Liard River with groundwater: Raven + MODFLOW 6 example

A hypothetical aquifer is added under three nested Central Liard subbasins (52 -> 48 -> 43; 85 HRUs).
The other 62 subbasins run exactly as in the original Raven model. Adding groundwater took four edits to
ordinary Raven files plus one new `.rvg` file; this README walks through them so you can repeat them on
your own model.

## 0. Prerequisites
- Raven built from this branch (CMake, or `make` in `src/`; add `-Dnetcdf` for NetCDF head output).
- The MODFLOW 6 library: `python tools/get_mf6.py` in the source folder (MODFLOW 6.8.1, the tested version, into
  `lib/mf6/`, where Raven finds it), or a CMake build, which fetches it. To use a library elsewhere, set
  `RAVEN_MF6_LIB=/path/to/libmf6.so` (Windows: `set RAVEN_MF6_LIB=C:\mf6\libmf6.dll`).

## 1. Run it
```
cd examples/Liard_groundwater
../../Raven.exe Liard -o output/    # 20 years, daily; about 10 minutes on one core
python3 tools_audit.py output       # optional: automated water-balance and consistency audit
```
Groundwater results: `output/GWBudget.csv`, `GWHRUState.csv`, `GWConnectivity.csv`, `GWModelSummary.txt`,
and the generated MODFLOW 6 model in `output/mf6/`. Streamflow and diagnostics are in the usual Raven files.

## 2. What was changed, step by step

**Edit 1 - `Liard.rvh`: say which HRUs have groundwater.** The existing `AQUIFER_PROFILE` column holds a profile
name (`VALLEY_ALLUVIUM` for forest/shrub/wetland HRUs, `UPLAND_TILL` for barren HRUs) or `[NONE]`. A group
`AquiferHRUs` lists the same HRUs, for use in `.rvi` conditionals.

**Edit 2 - `Liard.rvp`: describe the hydrogeology** (end of the file):
```
:AquiferClasses          # materials: K_HORIZ, K_VERT [m/d], SPEC_STORAGE [1/m], SPEC_YIELD, POROSITY
:AquiferProfiles         # layer stacks: class, thickness or TO_BEDROCK, layer type, per layer
:AquiferProfileParameters# initial water-table depth, bedrock fallback, seepage leakance, ...
```

**Edit 3 - `Liard.rvi`: switch groundwater on and route recharge into it:**
```
:GroundwaterModel MODFLOW6
:DefineHRUGroups AquiferHRUs
...
:Percolation PERC_CONSTANT FAST_RESERVOIR SLOW_RESERVOIR
  :-->Conditional HRU_GROUP IS_NOT AquiferHRUs      # (existing process, now restricted)
:Percolation PERC_CONSTANT FAST_RESERVOIR GROUNDWATER
  :-->Conditional HRU_GROUP IS AquiferHRUs          # new: this water becomes MODFLOW recharge
```
Groundwater seepage returns to the lowest soil layer (`SLOW_RESERVOIR`) and reaches the streams through the
model's existing baseflow process; aquifer-river exchange enters the reaches directly.

**Edit 4 - `Liard.rvg` (new): geometry and settings.** HRU polygons, river lines, DEM, bedrock, grid size.
Optional features (steady-state start, wetland exchange, wells, monitoring wells, boundaries, NetCDF heads)
are listed in the file, commented out, ready to switch on.

No MODFLOW file is written by hand: Raven builds the grid (3 layers x 85 x 54 cells at 2 km, 3,252 active
columns, 178 river cells), writes `output/mf6/`, and runs MODFLOW 6 alongside Raven every day.

## 3. About the input data in this example
The Liard package has no HRU polygons, DEM or bedrock data, so stand-ins were generated for testing:
- `Liard_HRUs_standin.geojson` - Voronoi cells of the `.rvh` HRU centroids clipped to each subbasin
- `Liard_dem_standin.asc`, `Liard_bedrock_standin.asc` - DEM interpolated from HRU elevations; bedrock
  60 m (valley) / 15 m (upland) / 5 m (glacier) below it
- `Liard_rivers.geojson` - the package's river network; `Liard_subbasins.geojson` - its subbasin polygons

`tools_build_standin_inputs.py` made these files and the four edits from the original package (run it only
on a fresh copy of the original files). For a real study, replace the stand-ins with your HRU polygons (from
BasinMaker or your GIS), a DEM and a bedrock or aquifer-base surface; keep the HRU IDs in the GeoJSON equal
to those in the `.rvh`.

Because the stand-in polygons are Voronoi cells, their areas do not match the `.rvh` `AREA` values: Raven warns
"polygon area ... differs from .rvh AREA" for every coupled HRU, and some differ by a factor of 50-100. Water volumes
are still exchanged exactly (each HRU's recharge is spread over the cells its polygon covers), but recharge is then
concentrated on fewer or more cells than the real HRU would cover, which distorts local heads. With real HRU polygons
these warnings disappear; treat them in your own model as a sign that the polygons and the `.rvh` do not match.

The original package was built for Raven 2.9.1; for Raven 4.x, Windows `\` paths were changed to `/` and
`DEP_THRESHHOLD` was renamed `DEP_THRESHOLD`.

## 4. Results you should see (final audit)
Recharge transferred exactly; MODFLOW cumulative balance residual ~3e-6 of gross flow; worst day 0.04%;
every step converged; no river asked for more water than it carried; Raven balance error 0.52 mm.
NSE at gauges 10BE004 / 10BE005 / 10BE010: 0.850 / 0.821 / 0.634 (uncoupled 0.867 / 0.835 / 0.681;
not recalibrated, stand-in geometry).

## 5. Adapting this to your own model
1. Build HRU polygons with the same HRU IDs as your `.rvh` (GeoJSON, lon/lat).
2. Choose aquifer profiles; set `AQUIFER_PROFILE` per HRU; define classes and profiles in the `.rvp`.
3. Route the water your model sends to deep groundwater into `GROUNDWATER` for the coupled HRUs, and make
   sure no process draws from `GROUNDWATER` there (Raven checks this).
4. Write a `.rvg` (start from this one), pick a cell size (a few hundred to a few thousand metres).
5. Run, read `GWModelSummary.txt` and `GWBudget.csv`, run `tools_audit.py`, then calibrate with Ostrich
   (see `docs/ostrich_example`).
