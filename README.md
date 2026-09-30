# RavenHydroFramework

The code repository for the Raven Hydrological Modelling Framework developed at the University of Waterloo.

Release versions, tutorials, documentation, and more distributed at [the Raven website](https://raven.uwaterloo.ca/Main.html).

Intended for use with Visual Studio Community Edition 2022, but also provided with Windows/linux/unix/MacOS g++ makefile and CMake configuration file.

Note unconventional two-space tabbing conventions.
If you would like to work with the active development of Raven's core, please do so on a branch and coordinate commits to the trunk with the Raven development team.

Please contact us while you're at it - we love to have people helping out.

## Building Raven

Raven is built using CMake. A suggested sequence of commands to build Raven executable file is:

```bash
mkdir build
cd build
cmake [OPTIONS] ../
make
```

The `cmake` command can be configured with the following optional arguments (```[OPTIONS]``` above):

* `-DCOMPILE_LIB` to build Raven as a library (default: `OFF`)
* `-DCOMPILE_EXE` to build Raven as an executable (default: `ON`)

So that the ```cmake``` command to build Raven as solely a dynamic library becomes:

```bash
cmake -DCOMPILE_LIB=ON -DCOMPILE_EXE=OFF ../
```
Raven can alternately be bullt in unix/MacOS using the makefile provided with the source code (g++ must be installed on the machine). Lastly, it may be compiled within Visual Studio Community Edition 2022.

## MODFLOW 6 groundwater coupling

Raven can simulate groundwater with MODFLOW 6: it either builds a MODFLOW 6 model from its HRUs
(`:GroundwaterModel MODFLOW6`) or couples an existing one (`:MF6Simulation`), and exchanges recharge, seepage, river
leakage and other fluxes every time step. MODFLOW 6 is loaded at run time; `python tools/get_mf6.py` (or a CMake build)
downloads the library into `lib/mf6/`, where Raven finds it. Raven builds and runs unchanged without it.

* User guide: `docs/Raven_MF6_Groundwater_UserGuide.md`; manual: `docs/manual/` (LaTeX, with `main.pdf`)
* Examples: `examples/Liard_groundwater` (Raven builds the model), `examples/Liard_existing_model` (existing model)
* Tests: `python3 tests/run_all.py` (see `tests/README.md`); audit: `docs/AUDIT_REPORT.md`
