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
* `-DRAVEN_MF6=SOURCE|PREBUILT|PATH|OFF`: the MODFLOW 6 library for the groundwater coupling, put next to Raven
  (`cmake/RavenMF6.cmake`). Default: `SOURCE` (built from MODFLOW 6.8.1 with its meson build) when a Fortran compiler,
  meson >= 1.7, ninja and git are found, else `PREBUILT` (official release library, checked against pinned SHA-256)
  where MODFLOW releases one, else `OFF` with a warning. `PATH` uses `-DRAVEN_MF6_PATH=/path/to/libmf6.so`; offline:
  `-DRAVEN_MF6_ZIP=mf6.8.1_linux.zip` or `-DRAVEN_MF6_SOURCE_DIR=...`. With `OFF` the build is unchanged. `cmake --install`
  and `cpack` give one folder with Raven, the library and its Fortran runtime. The coupling's manual, user guide,
  examples and tests are in https://github.com/rarabzad/RavenHydroFramework-extra (branch `main`; clone it beside this
  repository).

So that the ```cmake``` command to build Raven as solely a dynamic library becomes:

```bash
cmake -DCOMPILE_LIB=ON -DCOMPILE_EXE=OFF ../
```
Raven can alternately be bullt in unix/MacOS using the makefile provided with the source code (g++ must be installed on the machine). Lastly, it may be compiled within Visual Studio Community Edition 2022.
