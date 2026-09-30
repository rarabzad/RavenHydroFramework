# MODFLOW 6 library

The Raven–MODFLOW 6 groundwater coupling loads the MODFLOW 6 shared library at run time. It is not stored in the
repository; get it with

    python tools/get_mf6.py            # MODFLOW 6.8.1 (the tested version) for this computer

or build Raven with CMake, which fetches it (option `RAVEN_GET_MF6`, on by default). Either way it lands here as
`libmf6.so` (Linux), `libmf6.dll` (Windows) or `libmf6.dylib` (macOS), with `libmf6_version.txt` naming its version.

Raven finds it here, or next to the Raven executable, without any setting. `:MF6Library` in the `.rvg` file or the
`RAVEN_MF6_LIB` environment variable override this. MODFLOW 6 is public-domain software of the U.S. Geological
Survey: https://github.com/MODFLOW-ORG/modflow6
