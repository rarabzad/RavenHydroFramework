# Release audit of the Raven–MODFLOW 6 coupling

This report sums up the audit of the complete package (code, manual, examples and tests) made before release: five phases
without the MODFLOW library, a sixth that repeated every coupled test with MODFLOW 6.8.1, a seventh with two USGS example
models and MODFLOW 6.6.3, and an eighth that carried the coupling over to Raven 4.15. The manual has the same
material in its verification chapter (section "Release audit" and the findings table).

## How the audit could work without MODFLOW

The MODFLOW 6 library could not be downloaded where the audit ran: the network blocked GitHub and USGS. The audit
depends on one property of the coupling: **Raven writes every MODFLOW input file, and `grid_cells.geojson`, before it
loads the library.** In the existing-model mode it also copies the model and rewrites the simulation, time and name
files first. If `RAVEN_MF6_LIB` points at a file that does not exist, the run stops with a message just after that point.
So everything except the solve itself could be checked. When the library became available (MODFLOW 6.8.1), every coupled
test was repeated on the final code (phase 6).

## Phase 1: code review

Three independent line-by-line reviews covered the coupling and every changed Raven file. The defects they found were all
corrected. The full list is in the manual's findings table under *Release audit*. The main ones:

| Area | Defect | Correction |
|---|---|---|
| Raven core | a reservoir drying under a minimum/override flow booked no evaporation (114 million m³ unbalanced in a 20-year test; 158 with 4.1.1) | books the evaporation the outflow formula removed, and no seepage |
| Raven core | reservoir constraint read before the first step without a value (4.1.1; found by UBSan) | initialised; outputs unchanged |
| Raven core | `GW_HEAD` observations bypassed Raven's observation matching; ensemble snapshot lost the last lateral inflow | corrected |
| Generated grids | rotated grids without refinement: cell centres and DIS origin in the unrotated frame (a boundary polygon then selected no cell, so no CHD package was written) | rotated frame |
| Generated grids | DISU: aquicludes bridged by vertical connections (2 794 in one test); face angles had the rotation added twice (30° off) | corrected |
| Generated grids | HRU mesh: a well outside the mesh crashed the run; the closing edge of each HRU outline got no seeds | corrected |
| Generated grids | imported cells: repeated vertices; faces at T-junctions of re-projected grids lost | de-duplicated; length-relative tolerance |
| Generated grids | drain smoothing below the top cell (2 794 of 3 252 drains with a thin top layer) | limited to the top layer |
| Reservoir exchange | bed at the model top when the lake stood below it (false gains); unsupplied gains created water; lake below the aquifer bottom accepted | bed at the lake level; gains taken back from the reach (new budget column); error |
| Existing models | `mfsim.nam` without `CONTINUE`; recharge on the model's own specified-head cells; restarts with time series, ATS, TVK/TVS or several models accepted; START_DATE_TIME/time of day on restart; non-Gregorian calendars; `:MF6GridFile` check by area only; second initialisation in ensembles failed | corrected or refused with a message |
| Windows | build failed (`GetCurrentTime` clashes with a `windows.h` macro) | renamed |

## Phase 2: builds and static analysis

- GCC and Clang, `-O0` and `-O3`: no warnings in the coupling's files. The warnings that remain are in untouched 4.1.1 code.
- cppcheck and clang-tidy: false positives only.
- Windows: MinGW build run under Wine. Results match Linux apart from line endings and the last digits of two
  near-zero storage columns (values around 10⁻¹⁷).

## Phase 3: MODFLOW input files and numbers

- `tests/mf6files_check.py` runs 16 generated set-ups covering every grid type, rotation, refinement, layers by
  elevation, aquicludes, boundaries, wells, lakes and imported cells. It checks that:
  - FloPy reads every model;
  - cells are valid and do not overlap;
  - the DIS origin and rotation match the cell file (0.5 m);
  - DISV/DISU vertices are right, DISU connections are symmetric and face angles match the vertices (1°);
  - no vertical connection crosses an aquiclude;
  - drains sit inside the top cell;
  - boundary polygons select the right cells.

  **127 of 127 checks pass.** The build before the audit fails five of the set-ups.
- Geometry against shapely and PROJ: clipping within 1.3·10⁻¹⁰, projections within 3 µm (9 coordinate systems).
- Scoping of the manual's numbers: the MODFLOW files of the regular, rotated and refined, quadtree and layers-by-elevation
  Liard runs, and of every Thiem grid except the HRU mesh, are identical before and after the audit. The rotated uniform
  grid differs only in its origin. Their results stand. Results that depend on changed code were repeated in phase 6.

## Phase 4: documentation

- All 70 `.rvg`/`.rvi` commands shown in the manual, user guide and examples are recognised by Raven.
- The three complete `.rvg` listings in the manual (including the full example of "Inputs and outputs" and the
  one in "Applying the coupling") build their models up to the library load.
- One snippet note added: a `:WellRate` series must cover the whole run.
- Visual review of every page: diagrams redrawn in a sans-serif font, overflowing labels and tables fixed, breakable
  appendix tables.
- Spell check clean.
- 106 pages, no overfull boxes, no undefined references.

## Phase 5: regression, sanitizers, fresh-user walkthrough

- The 20-year Liard model without groundwater gives outputs identical byte for byte to Raven 4.1.1. `Raven_errors.txt`
  differs only in the build date.
- AddressSanitizer and UndefinedBehaviorSanitizer were run on the 16 file-check set-ups, the 52 input set-ups, seven
  reservoir and ensemble cases and three years of the Liard model. They found one uninitialised value, in 4.1.1's
  reservoir code, now corrected. Otherwise there was no error.
- `tests/edge.py --nolib` runs 35 generated-model and 17 existing-model set-ups (36 and 17 after phase 6). **All behave as required**: each
  faulty one stops with a message naming the problem, and each valid one prepares its files. The valid existing-model
  preparations were checked in detail: stress periods split at the start date, restart renumbering, and the original
  model left untouched.
- A reviewer who had not seen the package followed `START_HERE.md` from the zip:
  - the build took 49 s, with one warning in 4.1.1 code;
  - both examples stop with a clear library message after writing their MODFLOW files;
  - `run_all.py --nolib` passes every step;
  - `make_repo.sh` builds the two-branch repository, and the branch compiles;
  - the Overleaf folder compiles to the same PDF as the one shipped.

  The reviewer's findings were all corrected:
  - the example READMEs called `Raven.exe` without its path;
  - `START_HERE` said every audit line should be near 0, but ratios should be near 1;
  - two stale lines in the user guide;
  - a compiled Python file in the example;
  - the existing-model library message lacked the `RAVEN_MF6_LIB` hint;
  - `tools_audit.py` gave a traceback on an unfinished run;
  - the example README did not explain the area warnings caused by its stand-in HRU polygons.

## Phase 6: coupled verification with MODFLOW 6.8.1

`python3 tests/run_all.py` with the MODFLOW 6.8.1 library (Linux build), final code, 45 minutes on one core. Every step
passes (full report: `docs/VERIFICATION_REPORT_MF6.8.1.md`):

| Step | Result |
|---|---|
| unit tests (grids, geometry, projections) | all pass; projections within 3 µm of PROJ |
| MODFLOW input files (`mf6files_check.py`) | 127 of 127 checks |
| faulty and unusual input (`edge.py`, without and with MODFLOW) | 53 of 53 |
| Thiem well test | every grid converges to the analytical slope; HRU mesh +0.57 / +0.29 / −0.06 % |
| restarts, five grid types | largest deviation 9.9·10⁻⁷ |
| lake refinement | reservoir identity exact (sent = previous step's volume to 1.6·10⁻¹⁰, delivered = sent) |
| Liard on every grid (`liard_grids.py`) | regular, rotated, quadtree and elevation rows reproduce the manual; HRU mesh: 14 610 cells, NSE 0.717 / 0.859 / 0.627 |
| reservoir and ensemble behaviour (`behaviour_check.py`, new) | 50 km² lake delivered exactly one step later; ensemble members identical to the single run |
| MODFLOW's list-file budget (`budget_xcheck.py`) | 45 package budgets, worst difference 1.1·10⁻⁷ (printing precision) |
| randomized grid rounds | 72 of 72 (20 restarts, largest deviation 3.1·10⁻⁶) and 24 of 24 with reservoir and lake refinement |
| existing-model tests | safeguards 14/14, round trip 1.7·10⁻⁶, restarts 1.4·10⁻⁷, randomized 20/20 (6 restarts, within 2.7·10⁻⁷) |

The 20-year Liard run and the two-year existing-model example give outputs identical to the runs made before the audit, so
the data figures and the behaviour-test numbers drawn from them stand. All numbers the manual had marked "before the
release audit" were replaced with the final-code results.

Phase 6 found and fixed:
- **a Raven 4.1.1 crash:** a reservoir whose `:HRUID` is an HRU of another subbasin crashed the run (segmentation fault).
  It now stops with a message, and `edge.py` has a case for it;
- **the randomized harness ran on the wrong model:** it used the full three-subbasin model instead of the subbasin-52
  model, and set the crest from `.rvh` elevations that lie below the stand-in DEM. Both are corrected: the crest is now the
  highest land inside the lake polygon. The restart criterion now scales by each flow's largest value, as the
  existing-model harness does. One set-up differed by 0.16 m³/d of storage release on the first restarted day, which the
  old criterion reported as 1.5·10⁻⁴;
- **the budget cross-check compared nothing:** newer FloPy names list-file budget lines by package type. The check now
  reads MODFLOW's list file by package name, and it fails if nothing is compared;
- **the existing-model harness:** its pass flag was stored as text, and a restart comparison counted as run when only one
  of its two runs had finished;
- **enhancements:**
  - `run_all.py` merges partial reruns (`--only`) into the report;
  - `liard_grids.py --table` prints the manual's table;
  - `tools_audit.py` audits ensemble members;
  - the existing-model summary no longer shows zero river/well/boundary counts for the model's own packages.

## Phase 7: USGS example models and a second MODFLOW version

Two models of the USGS MODFLOW 6 example collection were coupled in the existing-model mode with `:OverlapWeights`
(`tests/existing/usgs_sagehen.py`, `usgs_capture.py`), and the whole test suite was run with MODFLOW 6.8.1 and 6.6.3.

- **Sagehen Creek** (real watershed; Newton, UZF, SFR, drains, water mover, 399 daily periods). Raven takes over the UZF
  package; the recharge identity, river bookkeeping and seepage returned are exact in every case. MODFLOW's cumulative
  residual:
  - 8.5·10⁻⁷ with the model's own loose solver;
  - 2.8·10⁻¹⁰ with a tight solver;
  - 1.3·10⁻¹¹ with Raven's seepage drains.

  A restart inside the daily periods reproduces the continuous run to 3.3·10⁻⁸. A stream package feeding the mover is
  refused.
- **Stream capture** (steady state, time in seconds, river, wells). The one-second steady period is stretched over
  Raven's run. The accounts are exact (residual 2.8·10⁻⁹) and a restart matches within 1.0·10⁻⁷. With the model's own
  tolerance of 10⁻⁸ m, four days do not converge; Raven flags them.
- **MODFLOW 6.6.3.** Every step of `run_all.py` passes (all 20 steps, both USGS models included; full reports:
  `docs/VERIFICATION_REPORT_MF6.8.1.md` and `docs/VERIFICATION_REPORT_MF6.6.3.md`). The Liard grid table agrees with 6.8.1
  to the last digit shown.

Found and fixed:
- **Water sent to the mover was missing from the balance.** The water the model's packages pass to the mover (MVR)
  was not counted: the Sagehen balance was off by 0.1 %, and by up to 30 % with the model's UZF kept. It is now counted.
- **A missing package stopped MODFLOW.** A mover that named a package Raven takes over stopped MODFLOW, and Raven with
  it, without a message. Those movers are now removed from the working copy, and the count is reported.
- **Double counting through the mover.** A stream package that feeds the mover would have been counted twice. It is now
  refused.
- **Steady-state models could not be coupled.** Models with one short stress period (often one second) did not split
  into Raven days. The period is now stretched over Raven's run. Without a storage package, MODFLOW's error is reported
  as error, not as storage.
- **Restarts refused with a mover or advanced packages.** The mover is now renumbered like a list package; SFR, LAK, MAW
  and UZF are renumbered like arrays.
- **MODFLOW 6.6 stalled on drying cells.** The `UNDER_RELAXATION` keyword of Raven's `NEWTON` option made steps with a
  drying cell stall, and the default Liard grid failed on 295 of 344 days. `NEWTON` is now written alone; the results
  with 6.8.1 are identical. Raven reports the library's version, and warns when an existing model uses the keyword with
  a version before 6.8.

## Phase 8: carried over to Raven 4.15

The coupling was developed on Raven 4.1.1. For the pull request it was carried over to Raven 4.15 (v593, the `main` of
the fork), where the source files live in `src/`.

- **Merge.** A three-way merge left 26 overlapping edits in 22 files, each resolved by hand (listed in the manual,
  section "Carrying the coupling over to Raven 4.15"). Raven 4.15 had already fixed the BMI index in the same way, so
  its version was kept. Its new `RebootTimeVariables` (ensemble reset of the balance arrays) is kept beside the
  coupling's restoring of routing memory and reservoir states.
- **Models without groundwater, against unmodified 4.15 (both built with NetCDF):**
  - the Liard model: identical byte for byte;
  - Raven's benchmark models: the 11 that run are identical, except that Lake of the Woods, whose reservoir dries
    out on 45 days, books less evaporation on those days (hydrographs and stages identical; balance error 5.44 → 4.97
    mm);
  - the other 14 benchmark set-ups stop with the same input error in unmodified 4.15.
- **Tests.** Raven 4.15 prints its closing message to the error stream, so the test scripts now read both streams.
  The whole suite passes on the ported code with MODFLOW 6.8.1 and 6.6.3, with the same results as on 4.1.1 (the
  Liard grid table agrees to the last digit shown).
- **Builds.** GNU make (`src/Makefile`) and CMake (NetCDF found, libdl linked) both build and pass.
- **Getting MODFLOW 6.** `tools/get_mf6.py` downloads the tested MODFLOW 6.8.1 library for the computer from the MODFLOW 6
  releases into `lib/mf6/`, and a CMake build does the same. Raven finds the library there, or next to its executable,
  without any setting. The downloaded Linux library is identical to the one all tests ran with. The Windows and Apple
  Silicon zips were checked to hold the library in the same place; 6.8.1 has no Intel-Mac build. A new test step
  (`mf6lib`) checks the offline install, both lookup places and the message when no library is found.

## Still open

- **Native Windows (Visual Studio) and macOS builds** of the 4.15 port.
- **Linking by map projection on a real model.** `:MF6CRS` and `:MF6GridFile` have not been tried on a real,
  map-referenced USGS model. The example models are not georeferenced, so they were linked with `:OverlapWeights`.
- **Models with several groundwater models or exchanges** (for example LGR). These are refused for restarts and were not
  coupled.

## Open item for the Raven developers

In a test where a reservoir is drawn empty by a minimum flow, 1 072 of the 1 086 dry days balance exactly. The remaining
13.7 million m³ falls almost entirely (13.69 million) on the day after the reservoir empties. On that day the reported
outflow still averages in the previous day's end-of-step release rate, although the lake is empty. This is in Raven
4.1.1's reservoir routine, not in the coupling. For comparison, the same test gives 158 million m³ with 4.1.1.

## Reproducing the audit

```
(cd src && make)                            # or CMake: mkdir build; cd build; cmake ..; make
python3 tests/run_all.py --nolib            # unit tests, file checks, input checks (minutes; no MODFLOW needed)
export RAVEN_MF6_LIB=/path/to/libmf6.so
export RAVEN_USGS_EXAMPLES=/path/to/mf6examples.zip   # optional: the USGS example-model steps
python3 tests/run_all.py                    # everything, one report in tests/_work/run_all_report.md
```
