# RavenMF6.cmake: the MODFLOW 6 shared library (libmf6) for Raven's groundwater coupling.
#
# Raven loads libmf6 at run time (dlopen/LoadLibrary) and finds it next to its executable. This module puts it there, in
# the build tree and in the install tree, in one of four ways chosen by the cache option RAVEN_MF6:
#
#   SOURCE    build libmf6 from the MODFLOW 6 source with MODFLOW's own meson build (needs a Fortran compiler, meson >= 1.7,
#             ninja, and git unless RAVEN_MF6_SOURCE_DIR or RAVEN_MF6_SOURCE_URL is given). Source: git clone of release
#             tag RAVEN_MF6_VERSION (6.8.1 checked against its pinned commit), or RAVEN_MF6_SOURCE_DIR (a source folder),
#             or RAVEN_MF6_SOURCE_URL (an archive; a pinned release zip works, its src/ equals the tag).
#   PREBUILT  the official USGS release zip for this platform, downloaded at configure time and checked against the pinned
#             SHA256 (6.8.1: linux, macarm, win64), or a local zip given by RAVEN_MF6_ZIP (offline).
#   PATH      a library you already have: RAVEN_MF6_PATH=/path/to/libmf6.so (copied and installed as libmf6.so).
#   OFF       Raven alone (the coupling then needs a library at run time: :MF6Library, RAVEN_MF6_LIB, or next to Raven).
#
# Default (RAVEN_MF6 empty, re-chosen at every configure; the cache keeps only a value you give): an offline option names
# its mode (RAVEN_MF6_PATH -> PATH, RAVEN_MF6_ZIP -> PREBUILT, RAVEN_MF6_SOURCE_DIR/SOURCE_URL -> SOURCE); otherwise SOURCE
# when the source tools are found, else PREBUILT when MODFLOW releases a library for this platform, else OFF with a
# warning (-DRAVEN_MF6=OFF silences it). A default that cannot be completed (download failed) falls back to OFF with a
# warning; an explicit (or inferred) mode that cannot be completed stops.
#
# Files: RavenMF6Verify.cmake (git clone + pinned-commit check, build time), RavenMF6Copy.cmake (copies the library and
# writes libmf6_version.txt next to Raven, build time), RavenMF6Bundle.cmake (install time: the Fortran runtime the
# library needs, honouring DESTDIR).
#
# Results for CMakeLists.txt (normal variables): RAVEN_MF6_MODE (SOURCE/PREBUILT/PATH/OFF), RAVEN_MF6_LIB_FILE (the
# library before renaming), RAVEN_MF6_LIB_VERSION (MODFLOW version, "" when unknown), RAVEN_MF6_FORTRAN_USED (source
# builds); internal cache RAVEN_MF6_LIBNAME (libmf6.so/.dylib/.dll); function raven_mf6_rpath(target).

include_guard(GLOBAL)

set(RAVEN_MF6_TESTED_VERSION 6.8.1)
set(RAVEN_MF6_PINNED_COMMIT_6.8.1 5b2ebd1a95720d5c3732dc482e05fbf559e6ef41)        # tag 6.8.1 of MODFLOW-ORG/modflow6
set(RAVEN_MF6_PIN_mf6.8.1_linux.zip  439a29f3e60f8e88b0ac083db10c92e3474b203d18ca6583863d8464b811ba90)
set(RAVEN_MF6_PIN_mf6.8.1_macarm.zip 2167d553d8f0c57c6c356156cb158b14e60684199785a82031b5c6cccf0b1e4c)
set(RAVEN_MF6_PIN_mf6.8.1_win64.zip  de54d0d78f2c3e5aa2104ec688c458e85ac4cbc770feff24015df4ba7d027837)
set(_RMF6_DIR "${CMAKE_CURRENT_LIST_DIR}")

# RPATH of the Raven executable: $ORIGIN (macOS @loader_path) so that a moved install finds what lies beside it, used in
# the build tree too (no relink at install), plus the folders of linked libraries outside the system folders (netCDF from
# conda or spack), without which a build-tree Raven would no longer start.
function(raven_mf6_rpath target)
  if(APPLE)
    set(_origin "@loader_path")
  else()
    set(_origin "$ORIGIN")
  endif()
  set_target_properties(${target} PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${_origin}"
                        INSTALL_RPATH_USE_LINK_PATH TRUE)
endfunction()
if(RAVEN_MF6_FUNCTIONS_ONLY)       # tests/cmake_checks.py includes this file for raven_mf6_rpath() alone
  return()
endif()

# ---------------------------------------------------------------------------------------------------------------------
# 1. the mode
# ---------------------------------------------------------------------------------------------------------------------
set(_rmf6_given "")
if(DEFINED RAVEN_MF6 AND NOT "${RAVEN_MF6}" STREQUAL "")
  string(TOUPPER "${RAVEN_MF6}" _rmf6_given)
  if(NOT _rmf6_given MATCHES "^(SOURCE|PREBUILT|PATH|OFF)$")
    message(FATAL_ERROR "RAVEN_MF6=${RAVEN_MF6} is not a valid choice: use SOURCE, PREBUILT, PATH or OFF, or leave it "
                        "empty for the default (-DRAVEN_MF6=)")
  endif()
endif()

set(_rmf6_offline RAVEN_MF6_PATH RAVEN_MF6_ZIP RAVEN_MF6_ZIP_SHA256 RAVEN_MF6_SOURCE_DIR RAVEN_MF6_SOURCE_URL)
macro(_rmf6_report_ignored mode)
  set(_rmf6_used_PATH RAVEN_MF6_PATH)
  set(_rmf6_used_PREBUILT RAVEN_MF6_ZIP RAVEN_MF6_ZIP_SHA256)
  set(_rmf6_used_SOURCE RAVEN_MF6_SOURCE_DIR RAVEN_MF6_SOURCE_URL)
  set(_rmf6_used_OFF "")
  foreach(_v IN LISTS _rmf6_offline)
    if(NOT "${${_v}}" STREQUAL "" AND NOT _v IN_LIST _rmf6_used_${mode})
      message(STATUS "RAVEN_MF6: ${_v} is ignored with RAVEN_MF6=${mode}")
    endif()
  endforeach()
endmacro()

# A plain -DRAVEN_MF6=OFF leaves the build exactly as without this module: no cache entries, targets, RPATH or install
# rules are added (the compile and link commands are those of Raven before the MODFLOW 6 build integration).
if(_rmf6_given STREQUAL "OFF")
  _rmf6_report_ignored(OFF)
  message(STATUS "RAVEN_MF6=OFF: Raven is built without the MODFLOW 6 library (a coupled model then needs one at run time:"
                 " :MF6Library, RAVEN_MF6_LIB, or libmf6 next to Raven)")
  set(RAVEN_MF6_MODE OFF)
  return()
endif()

if(NOT DEFINED CACHE{RAVEN_MF6})
  set(RAVEN_MF6 "" CACHE STRING "MODFLOW 6 library: SOURCE, PREBUILT, PATH or OFF; empty = default, chosen at every configure")
else()
  set_property(CACHE RAVEN_MF6 PROPERTY TYPE STRING)
endif()
set_property(CACHE RAVEN_MF6 PROPERTY STRINGS ";SOURCE;PREBUILT;PATH;OFF")
set(RAVEN_MF6_VERSION "${RAVEN_MF6_TESTED_VERSION}" CACHE STRING "MODFLOW 6 release for SOURCE (git tag) and PREBUILT")
set(RAVEN_MF6_PATH "" CACHE FILEPATH "PATH mode: the MODFLOW 6 library to use")
set(RAVEN_MF6_ZIP "" CACHE FILEPATH "PREBUILT offline: a MODFLOW 6 release zip on disk")
set(RAVEN_MF6_ZIP_SHA256 "" CACHE STRING "expected SHA256 of RAVEN_MF6_ZIP or of the downloaded release zip")
set(RAVEN_MF6_SOURCE_DIR "" CACHE PATH "SOURCE offline: a MODFLOW 6 source folder (holding meson.build)")
set(RAVEN_MF6_SOURCE_URL "" CACHE STRING "SOURCE: URL of a MODFLOW 6 source archive (tar.gz/zip; file:// works)")
set(RAVEN_MESON "" CACHE FILEPATH "meson for SOURCE (empty: searched on PATH)")
set(RAVEN_MF6_FORTRAN "" CACHE FILEPATH "Fortran compiler for SOURCE (empty: FC, CMake's check, gfortran/ifx on PATH)")
set(RAVEN_MF6_JOBS "" CACHE STRING "parallel jobs of the MODFLOW 6 source build (empty: meson's default)")
set(RAVEN_MF6_RELEASE_URL "https://github.com/MODFLOW-ORG/modflow6/releases/download" CACHE STRING
    "base URL of the MODFLOW 6 release zips (<url>/<version>/mf<version>_<platform>.zip)")
set(RAVEN_MF6_GIT_URL "https://github.com/MODFLOW-ORG/modflow6.git" CACHE STRING "MODFLOW 6 git repository for SOURCE")
mark_as_advanced(RAVEN_MESON RAVEN_MF6_FORTRAN RAVEN_MF6_JOBS RAVEN_MF6_RELEASE_URL RAVEN_MF6_GIT_URL)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(RAVEN_MF6_LIBNAME libmf6.dll CACHE INTERNAL "")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(RAVEN_MF6_LIBNAME libmf6.dylib CACHE INTERNAL "")
else()
  set(RAVEN_MF6_LIBNAME libmf6.so CACHE INTERNAL "")
endif()
set(_rmf6_work "${CMAKE_BINARY_DIR}/mf6")
file(MAKE_DIRECTORY "${_rmf6_work}")

# explicit, inferred from an offline option, or default
set(_rmf6_explicit FALSE)
set(RAVEN_MF6_MODE "")
if(_rmf6_given)
  set(_rmf6_explicit TRUE)
  set(RAVEN_MF6_MODE ${_rmf6_given})
  set(_rmf6_how "RAVEN_MF6=${_rmf6_given}")
elseif(NOT "${RAVEN_MF6_PATH}" STREQUAL "")
  set(_rmf6_explicit TRUE)
  set(RAVEN_MF6_MODE PATH)
  set(_rmf6_how "PATH (from RAVEN_MF6_PATH)")
elseif(NOT "${RAVEN_MF6_ZIP}" STREQUAL "")
  set(_rmf6_explicit TRUE)
  set(RAVEN_MF6_MODE PREBUILT)
  set(_rmf6_how "PREBUILT (from RAVEN_MF6_ZIP)")
elseif(NOT "${RAVEN_MF6_SOURCE_DIR}" STREQUAL "" OR NOT "${RAVEN_MF6_SOURCE_URL}" STREQUAL "")
  set(_rmf6_explicit TRUE)
  set(RAVEN_MF6_MODE SOURCE)
  set(_rmf6_how "SOURCE (from RAVEN_MF6_SOURCE_DIR/RAVEN_MF6_SOURCE_URL)")
endif()

# ---------------------------------------------------------------------------------------------------------------------
# 2. helpers
# ---------------------------------------------------------------------------------------------------------------------
# stop (explicit mode) or fall back to OFF with a warning (default)
macro(_rmf6_giveup msg)
  if(_rmf6_explicit)
    message(FATAL_ERROR "RAVEN_MF6: ${msg}\nChoose another way to get the library (-DRAVEN_MF6=SOURCE, PREBUILT or PATH "
                        "with -DRAVEN_MF6_PATH=/path/to/${RAVEN_MF6_LIBNAME}), or build Raven without it: -DRAVEN_MF6=OFF")
  endif()
  message(WARNING "RAVEN_MF6: ${msg}\nRaven is built without the MODFLOW 6 library (a coupled model then needs one at run "
                  "time: :MF6Library, RAVEN_MF6_LIB, or ${RAVEN_MF6_LIBNAME} next to Raven). To build it from source, "
                  "install a Fortran compiler (gfortran or ifx), meson >= 1.7, ninja and git (-DRAVEN_MF6=SOURCE); to use "
                  "a library you have: -DRAVEN_MF6_PATH=/path/to/${RAVEN_MF6_LIBNAME}; to silence this warning: "
                  "-DRAVEN_MF6=OFF")
  set(RAVEN_MF6_MODE OFF)
  return()
endmacro()

# MODFLOW version in a meson.build (project(... version: 'X.Y.Z' ...)), "" if none
function(_rmf6_meson_version file out)
  set(_v "")
  if(EXISTS "${file}")
    file(READ "${file}" _t LIMIT 4000)
    if(_t MATCHES "project\\([^)]*version[ \t\r\n]*:[ \t\r\n]*'([0-9][0-9A-Za-z.+-]*)'")
      set(_v "${CMAKE_MATCH_1}")
    endif()
  endif()
  set(${out} "${_v}" PARENT_SCOPE)
endfunction()

# target architecture: on Windows the generator platform (-A ARM64), else the compiler's, else the processor (); elsewhere a single CMAKE_OSX_ARCHITECTURES entry, else CMAKE_SYSTEM_PROCESSOR
function(_rmf6_target_arch out)
  set(_a "")
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    if(CMAKE_GENERATOR_PLATFORM)
      set(_a "${CMAKE_GENERATOR_PLATFORM}")
    elseif(CMAKE_CXX_COMPILER_ARCHITECTURE_ID)
      set(_a "${CMAKE_CXX_COMPILER_ARCHITECTURE_ID}")
    else()
      set(_a "${CMAKE_SYSTEM_PROCESSOR}")
    endif()
  elseif(APPLE OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    list(LENGTH CMAKE_OSX_ARCHITECTURES _n)
    if(_n EQUAL 1)
      set(_a "${CMAKE_OSX_ARCHITECTURES}")
    else()
      set(_a "${CMAKE_SYSTEM_PROCESSOR}")
    endif()
  else()
    set(_a "${CMAKE_SYSTEM_PROCESSOR}")
  endif()
  set(${out} "${_a}" PARENT_SCOPE)
endfunction()

# release zip suffix for this platform ("" + reason when MODFLOW releases no library for it)
function(_rmf6_release_suffix version out why)
  _rmf6_target_arch(_arch)
  string(TOLOWER "${_arch}" _la)
  set(_s "")
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(_sys Windows)
    if(_la MATCHES "^(x64|amd64|x86_64|em64t)$")
      set(_s win64)
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(_sys macOS)
    if(_la MATCHES "^(arm64|aarch64)$")
      set(_s macarm)
    elseif(_la MATCHES "^(x86_64|amd64)$" AND version VERSION_LESS 6.8)
      set(_s mac)            # Intel Macs: release zips up to 6.7 (none from 6.8 on; 404)
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(_sys Linux)
    if(_la MATCHES "^(x86_64|amd64)$")
      set(_s linux)
    endif()
  else()
    set(_sys "${CMAKE_SYSTEM_NAME}")
  endif()
  set(${out} "${_s}" PARENT_SCOPE)
  set(${why} "MODFLOW 6 ${version} releases no library for ${_sys} ${_arch}" PARENT_SCOPE)
endfunction()

# the Fortran compiler meson will use: RAVEN_MF6_FORTRAN, then FC (meson's own choice; read here because CMake's check
# ignores FC with the Visual Studio generators), then CMake's check_language(Fortran), then gfortran/ifx on
# PATH
function(_rmf6_find_fortran out)
  set(_fc "")
  if(NOT "${RAVEN_MF6_FORTRAN}" STREQUAL "")
    find_program(_RMF6_FC_PROBE NAMES "${RAVEN_MF6_FORTRAN}")
    set(_fc "${_RMF6_FC_PROBE}")
    unset(_RMF6_FC_PROBE CACHE)
  elseif(NOT "$ENV{FC}" STREQUAL "")
    separate_arguments(_fcl NATIVE_COMMAND "$ENV{FC}")
    list(GET _fcl 0 _fc0)
    find_program(_RMF6_FC_PROBE NAMES "${_fc0}")
    set(_fc "${_RMF6_FC_PROBE}")
    unset(_RMF6_FC_PROBE CACHE)
  endif()
  if(NOT _fc)
    if("${CMAKE_Fortran_COMPILER}" MATCHES "NOTFOUND$")
      unset(CMAKE_Fortran_COMPILER CACHE)          # look again at every configure (a compiler installed since)
    endif()
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER AND EXISTS "${CMAKE_Fortran_COMPILER}")
      set(_fc "${CMAKE_Fortran_COMPILER}")
    endif()
  endif()
  if(NOT _fc)
    find_program(_RMF6_FC_PROBE NAMES gfortran ifx ifort)
    if(_RMF6_FC_PROBE)
      set(_fc "${_RMF6_FC_PROBE}")
    endif()
    unset(_RMF6_FC_PROBE CACHE)
  endif()
  set(${out} "${_fc}" PARENT_SCOPE)
endfunction()

# source-build tools; sets _rmf6_fc, _rmf6_meson, _rmf6_ninja, _rmf6_git (parent) and returns the missing ones in why
function(_rmf6_source_tools ok why)
  set(_miss "")
  _rmf6_find_fortran(_fc)
  if(NOT _fc)
    list(APPEND _miss "no Fortran compiler found (gfortran or ifx; FC or -DRAVEN_MF6_FORTRAN=...)")
  endif()
  set(_meson "")
  if(NOT "${RAVEN_MESON}" STREQUAL "")
    if(EXISTS "${RAVEN_MESON}" AND NOT IS_DIRECTORY "${RAVEN_MESON}")
      set(_meson "${RAVEN_MESON}")
    else()
      list(APPEND _miss "meson not found (RAVEN_MESON=${RAVEN_MESON} does not exist)")
    endif()
  else()
    find_program(_RMF6_MESON_PROBE NAMES meson meson.py)
    if(_RMF6_MESON_PROBE)
      set(_meson "${_RMF6_MESON_PROBE}")
    else()
      list(APPEND _miss "meson not found (pip install meson, or -DRAVEN_MESON=...)")
    endif()
    unset(_RMF6_MESON_PROBE CACHE)
  endif()
  if(_meson)
    execute_process(COMMAND "${_meson}" --version OUTPUT_VARIABLE _mv ERROR_VARIABLE _me RESULT_VARIABLE _mr
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE TIMEOUT 60)
    if(NOT _mr EQUAL 0 OR NOT _mv MATCHES "^[ \t]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)")
      string(SUBSTRING "${_mv}${_me}" 0 80 _mshow)
      list(APPEND _miss "meson ${_meson} gives no readable version ('${_mshow}')")
      set(_meson "")
    elseif(CMAKE_MATCH_1 VERSION_LESS 1.7)
      list(APPEND _miss "meson ${_meson} is version ${CMAKE_MATCH_1}; MODFLOW 6 needs meson 1.7 or newer")
      set(_meson "")
    endif()
  endif()
  set(_ninja "")
  if(CMAKE_GENERATOR MATCHES "Ninja" AND CMAKE_MAKE_PROGRAM)
    set(_ninja "${CMAKE_MAKE_PROGRAM}")
  else()
    find_program(_RMF6_NINJA_PROBE NAMES ninja ninja-build)
    if(_RMF6_NINJA_PROBE)
      set(_ninja "${_RMF6_NINJA_PROBE}")
    endif()
    unset(_RMF6_NINJA_PROBE CACHE)
  endif()
  if(NOT _ninja)
    list(APPEND _miss "ninja not found (pip install ninja)")
  endif()
  set(_git "")
  if("${RAVEN_MF6_SOURCE_DIR}" STREQUAL "" AND "${RAVEN_MF6_SOURCE_URL}" STREQUAL "")
    find_package(Git QUIET)                        # the source is cloned
    if(GIT_FOUND AND GIT_EXECUTABLE)
      set(_git "${GIT_EXECUTABLE}")
    else()
      list(APPEND _miss "git not found (needed to clone the MODFLOW 6 source; or give RAVEN_MF6_SOURCE_DIR/SOURCE_URL)")
    endif()
  endif()
  set(_rmf6_fc "${_fc}" PARENT_SCOPE)
  set(_rmf6_meson "${_meson}" PARENT_SCOPE)
  set(_rmf6_ninja "${_ninja}" PARENT_SCOPE)
  set(_rmf6_git "${_git}" PARENT_SCOPE)
  list(JOIN _miss "; " _w)
  if(_miss)
    set(${ok} FALSE PARENT_SCOPE)
  else()
    set(${ok} TRUE PARENT_SCOPE)
  endif()
  set(${why} "${_w}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------------------------------------------------
# 3. the default
# ---------------------------------------------------------------------------------------------------------------------
_rmf6_release_suffix("${RAVEN_MF6_VERSION}" _rmf6_suffix _rmf6_nozip_why)
if(NOT RAVEN_MF6_MODE)
  _rmf6_source_tools(_rmf6_tools_ok _rmf6_tools_why)
  if(_rmf6_tools_ok)
    set(RAVEN_MF6_MODE SOURCE)
    set(_rmf6_how "SOURCE (default: Fortran compiler, meson, ninja and git found)")
  elseif(_rmf6_suffix)
    set(RAVEN_MF6_MODE PREBUILT)
    set(_rmf6_how "PREBUILT (default; not SOURCE: ${_rmf6_tools_why})")
  else()
    _rmf6_giveup("no way to get the MODFLOW 6 library by default: ${_rmf6_nozip_why}, and it cannot be built from source "
                 "here (${_rmf6_tools_why}).")
  endif()
endif()
_rmf6_report_ignored(${RAVEN_MF6_MODE})
message(STATUS "RAVEN_MF6: ${_rmf6_how}")

set(RAVEN_MF6_LIB_FILE "")
set(RAVEN_MF6_LIB_VERSION "")
set(RAVEN_MF6_FORTRAN_USED "")
set(_rmf6_desc "")
set(_rmf6_meson_build "")

# ---------------------------------------------------------------------------------------------------------------------
# 4a. PATH
# ---------------------------------------------------------------------------------------------------------------------
if(RAVEN_MF6_MODE STREQUAL "PATH")
  if("${RAVEN_MF6_PATH}" STREQUAL "")
    _rmf6_giveup("RAVEN_MF6=PATH needs the library: -DRAVEN_MF6_PATH=/path/to/${RAVEN_MF6_LIBNAME}")
  elseif(IS_DIRECTORY "${RAVEN_MF6_PATH}")
    _rmf6_giveup("RAVEN_MF6_PATH=${RAVEN_MF6_PATH} is a folder; give the library file itself (${RAVEN_MF6_LIBNAME})")
  elseif(NOT EXISTS "${RAVEN_MF6_PATH}")
    _rmf6_giveup("RAVEN_MF6_PATH=${RAVEN_MF6_PATH} does not exist")
  endif()
  get_filename_component(RAVEN_MF6_LIB_FILE "${RAVEN_MF6_PATH}" REALPATH)    # libmf6.so -> libmf6.so.6.8.1 etc.
  set(_rmf6_desc "library given by RAVEN_MF6_PATH (${RAVEN_MF6_PATH}")
  if(NOT RAVEN_MF6_LIB_FILE STREQUAL RAVEN_MF6_PATH)
    string(APPEND _rmf6_desc " -> ${RAVEN_MF6_LIB_FILE}")
  endif()
  string(APPEND _rmf6_desc ")")
  message(STATUS "RAVEN_MF6: PATH, ${_rmf6_desc}; installed as ${RAVEN_MF6_LIBNAME}; MODFLOW version unknown (given library)")
endif()

# ---------------------------------------------------------------------------------------------------------------------
# 4b. PREBUILT
# ---------------------------------------------------------------------------------------------------------------------
if(RAVEN_MF6_MODE STREQUAL "PREBUILT")
  set(_pre "${_rmf6_work}/prebuilt")
  set(_stamp "${_rmf6_work}/prebuilt.stamp")
  string(TOLOWER "${RAVEN_MF6_ZIP_SHA256}" _usersha)
  string(STRIP "${_usersha}" _usersha)
  if(_usersha AND NOT _usersha MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR "RAVEN_MF6_ZIP_SHA256=${RAVEN_MF6_ZIP_SHA256} is not a SHA256 (64 hexadecimal digits)")
  endif()
  set(_zip "")
  if(NOT "${RAVEN_MF6_ZIP}" STREQUAL "")
    # a zip on disk: stamped by its own SHA256 (a replaced zip, or another with the same name, is unpacked again) and compared with the expected SHA256 at EVERY configure
    if(NOT EXISTS "${RAVEN_MF6_ZIP}" OR IS_DIRECTORY "${RAVEN_MF6_ZIP}")
      _rmf6_giveup("RAVEN_MF6_ZIP=${RAVEN_MF6_ZIP} does not exist (or is not a file)")
    endif()
    set(_zip "${RAVEN_MF6_ZIP}")
    get_filename_component(_zname "${_zip}" NAME)
    file(SHA256 "${_zip}" _actual)
    set(_key "${_actual}")
    if(_usersha)
      set(_expected "${_usersha}")
      set(_whose "RAVEN_MF6_ZIP_SHA256")
      set(_rmf6_desc "zip ${_zname} given by RAVEN_MF6_ZIP (bin/${RAVEN_MF6_LIBNAME}; SHA256 of RAVEN_MF6_ZIP_SHA256)")
    elseif(DEFINED RAVEN_MF6_PIN_${_zname})
      set(_expected "${RAVEN_MF6_PIN_${_zname}}")
      set(_whose "the official MODFLOW release ${_zname}")
      set(_rmf6_desc "official release ${_zname} given by RAVEN_MF6_ZIP (SHA256 checked against the pinned checksum)")
    else()
      set(_expected "")
      message(WARNING "RAVEN_MF6_ZIP: ${_zname} is not a MODFLOW release zip Raven knows the checksum of, and "
                      "RAVEN_MF6_ZIP_SHA256 is not set: it is used unchecked. Give -DRAVEN_MF6_ZIP_SHA256=<its SHA256> "
                      "to check it.")
      set(_rmf6_desc "zip ${_zname} given by RAVEN_MF6_ZIP (bin/${RAVEN_MF6_LIBNAME}; not checked)")
    endif()
    if(_expected AND NOT _actual STREQUAL _expected)
      message(FATAL_ERROR "RAVEN_MF6_ZIP: ${_zip} has SHA256 ${_actual}, but ${_whose} has ${_expected}: a damaged or "
                          "different file (truncated download?). Download it again, or correct RAVEN_MF6_ZIP_SHA256.")
    endif()
    set(_src_desc "${_zip}")
  else()
    if(NOT _rmf6_suffix)
      _rmf6_giveup("${_rmf6_nozip_why} (no release zip); build it from source (-DRAVEN_MF6=SOURCE) or give one "
                   "(-DRAVEN_MF6=PATH -DRAVEN_MF6_PATH=...).")
    endif()
    set(_zname "mf${RAVEN_MF6_VERSION}_${_rmf6_suffix}.zip")
    set(_url "${RAVEN_MF6_RELEASE_URL}/${RAVEN_MF6_VERSION}/${_zname}")
    if(_usersha)
      set(_expected "${_usersha}")
      set(_whose "RAVEN_MF6_ZIP_SHA256")
      set(_rmf6_desc "release zip ${_url} (downloaded, SHA256 of RAVEN_MF6_ZIP_SHA256)")
    elseif(DEFINED RAVEN_MF6_PIN_${_zname})
      set(_expected "${RAVEN_MF6_PIN_${_zname}}")
      set(_whose "the official MODFLOW release ${_zname}")
      set(_rmf6_desc "official release ${_url} (downloaded, SHA256 checked against the pinned checksum)")
    else()
      set(_expected "")
      message(WARNING "RAVEN_MF6: Raven has no pinned checksum for ${_zname} (MODFLOW ${RAVEN_MF6_VERSION}); the download is "
                      "not checked. Give -DRAVEN_MF6_ZIP_SHA256=<its SHA256> to check it.")
      set(_rmf6_desc "release zip ${_url} (downloaded, not checked)")
    endif()
    # stamp named by the expected SHA256: a new expected SHA256 means a new download and check
    if(_expected)
      set(_key "${_expected}")
    else()
      set(_key "url:${_url}")
    endif()
    set(_src_desc "${_url}")
  endif()

  set(_old_key "")
  if(EXISTS "${_stamp}")
    file(STRINGS "${_stamp}" _sl)
    foreach(_l IN LISTS _sl)
      if(_l MATCHES "^([a-z]+)=(.*)$")
        set(_st_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
      endif()
    endforeach()
    set(_old_key "${_st_key}")
  endif()
  if(_old_key STREQUAL _key AND EXISTS "${_pre}/${_st_lib}" AND NOT "${_st_lib}" STREQUAL "")
    set(RAVEN_MF6_LIB_FILE "${_pre}/${_st_lib}")
    set(RAVEN_MF6_LIB_VERSION "${_st_version}")
    set(_member "${_st_lib}")
  else()
    if(NOT _zip)                                   # download (configure time)
      set(_zip "${_rmf6_work}/download/${_zname}")
      file(REMOVE "${_zip}")
      message(STATUS "RAVEN_MF6: downloading ${_url}")
      file(DOWNLOAD "${_url}" "${_zip}.part" STATUS _dst TLS_VERIFY ON INACTIVITY_TIMEOUT 120)
      list(GET _dst 0 _dcode)
      if(NOT _dcode EQUAL 0)
        list(GET _dst 1 _dmsg)
        file(REMOVE "${_zip}.part")
        _rmf6_giveup("could not download ${_url} (${_dmsg}). Offline: download the zip elsewhere and give "
                     "-DRAVEN_MF6_ZIP=/path/to/${_zname}.")
      endif()
      file(SHA256 "${_zip}.part" _actual)
      if(_expected AND NOT _actual STREQUAL _expected)
        file(REMOVE "${_zip}.part")
        message(FATAL_ERROR "RAVEN_MF6: the download ${_url} has SHA256 ${_actual}, but ${_whose} has ${_expected}: "
                            "damaged or truncated download, or a different file. Try again, or correct RAVEN_MF6_ZIP_SHA256.")
      endif()
      file(RENAME "${_zip}.part" "${_zip}")
    endif()
    # unpack only the library (and the meson.build files that give the version) into an emptied folder
    file(REMOVE_RECURSE "${_pre}")
    file(MAKE_DIRECTORY "${_pre}")
    file(ARCHIVE_EXTRACT INPUT "${_zip}" DESTINATION "${_pre}" PATTERNS "*${RAVEN_MF6_LIBNAME}" "*meson.build")
    file(GLOB_RECURSE _libs RELATIVE "${_pre}" "${_pre}/*${RAVEN_MF6_LIBNAME}")
    set(_member "")
    foreach(_l IN LISTS _libs)
      get_filename_component(_ln "${_l}" NAME)
      if(_ln STREQUAL RAVEN_MF6_LIBNAME)
        if(NOT _member OR (_l MATCHES "/bin/" AND NOT _member MATCHES "/bin/"))
          set(_member "${_l}")
        endif()
      endif()
    endforeach()
    if(NOT _member)
      file(REMOVE_RECURSE "${_pre}")
      _rmf6_giveup("${_src_desc} holds no ${RAVEN_MF6_LIBNAME}: is it a MODFLOW 6 release zip for ${CMAKE_SYSTEM_NAME}?")
    endif()
    # version: the zip's top-level meson.build, else the zip name (mfX.Y.Z_...), else unknown
    file(GLOB_RECURSE _mbs RELATIVE "${_pre}" "${_pre}/*meson.build")
    set(_mb "")
    foreach(_m IN LISTS _mbs)
      string(LENGTH "${_m}" _mlen)
      if(NOT _mb OR _mlen LESS _mblen)
        set(_mb "${_m}")
        set(_mblen ${_mlen})
      endif()
    endforeach()
    set(RAVEN_MF6_LIB_VERSION "")
    if(_mb)
      _rmf6_meson_version("${_pre}/${_mb}" RAVEN_MF6_LIB_VERSION)
    endif()
    if(NOT RAVEN_MF6_LIB_VERSION AND _zname MATCHES "^mf([0-9]+\\.[0-9]+\\.[0-9]+)_")
      set(RAVEN_MF6_LIB_VERSION "${CMAKE_MATCH_1}")
    endif()
    foreach(_m IN LISTS _mbs)
      file(REMOVE "${_pre}/${_m}")
    endforeach()
    set(RAVEN_MF6_LIB_FILE "${_pre}/${_member}")
    file(WRITE "${_stamp}" "key=${_key}\nlib=${_member}\nversion=${RAVEN_MF6_LIB_VERSION}\n")
    message(STATUS "RAVEN_MF6: unpacked ${_member} from ${_zname}")
  endif()
  string(APPEND _rmf6_desc " (${_member})")
  if(RAVEN_MF6_LIB_VERSION)
    set(_vshow "MODFLOW ${RAVEN_MF6_LIB_VERSION}")
  else()
    set(_vshow "MODFLOW version unknown")
  endif()
  message(STATUS "RAVEN_MF6: PREBUILT ${_vshow}, ${_rmf6_desc}")
endif()

# ---------------------------------------------------------------------------------------------------------------------
# 4c. SOURCE
# ---------------------------------------------------------------------------------------------------------------------
if(RAVEN_MF6_MODE STREQUAL "SOURCE")
  if(NOT DEFINED _rmf6_tools_ok)
    _rmf6_source_tools(_rmf6_tools_ok _rmf6_tools_why)
  endif()
  if(NOT _rmf6_tools_ok)
    _rmf6_giveup("cannot build MODFLOW 6 from source: ${_rmf6_tools_why}.")
  endif()
  set(RAVEN_MF6_FORTRAN_USED "${_rmf6_fc}")
  execute_process(COMMAND "${_rmf6_fc}" --version OUTPUT_VARIABLE _fcv ERROR_QUIET TIMEOUT 60)
  string(REGEX MATCH "^[^\r\n]*" _fcv "${_fcv}")
  set(_ep_prefix "${_rmf6_work}/source")
  set(_ep_bin "${_ep_prefix}/build")
  set(_ep_args "")
  if(NOT "${RAVEN_MF6_SOURCE_DIR}" STREQUAL "")
    get_filename_component(_ep_src "${RAVEN_MF6_SOURCE_DIR}" ABSOLUTE)
    if(NOT EXISTS "${_ep_src}/meson.build")
      _rmf6_giveup("RAVEN_MF6_SOURCE_DIR=${RAVEN_MF6_SOURCE_DIR} holds no meson.build: not a MODFLOW 6 source folder")
    endif()
    _rmf6_meson_version("${_ep_src}/meson.build" RAVEN_MF6_LIB_VERSION)
    list(APPEND _ep_args DOWNLOAD_COMMAND "" BUILD_ALWAYS TRUE)     # your folder: ninja rebuilds what changed
    set(_rmf6_desc "built from source ${_ep_src}")
  elseif(NOT "${RAVEN_MF6_SOURCE_URL}" STREQUAL "")
    set(_ep_src "${_ep_prefix}/src/mf6_source")
    get_filename_component(_uname "${RAVEN_MF6_SOURCE_URL}" NAME)
    list(APPEND _ep_args URL "${RAVEN_MF6_SOURCE_URL}" DOWNLOAD_NO_PROGRESS TRUE)
    if(DEFINED RAVEN_MF6_PIN_${_uname})
      list(APPEND _ep_args URL_HASH "SHA256=${RAVEN_MF6_PIN_${_uname}}")
      set(_rmf6_desc "built from the official release ${RAVEN_MF6_SOURCE_URL} (SHA256 checked against the pinned checksum)")
    else()
      message(WARNING "RAVEN_MF6_SOURCE_URL: ${_uname} has no pinned checksum; the download is not checked")
      set(_rmf6_desc "built from source ${RAVEN_MF6_SOURCE_URL} (not checked)")
    endif()
    # version: the source already downloaded by an earlier build, else the archive name (a release zip mfX.Y.Z_*.zip,
    # or a tag archive (modflow6-)vX.Y.Z.tar.gz/.tgz/.zip), else unknown until built
    _rmf6_meson_version("${_ep_src}/meson.build" RAVEN_MF6_LIB_VERSION)
    if(NOT RAVEN_MF6_LIB_VERSION)
      if(_uname MATCHES "^mf([0-9]+\\.[0-9]+\\.[0-9]+)_[A-Za-z0-9]+\\.zip$")
        set(RAVEN_MF6_LIB_VERSION "${CMAKE_MATCH_1}")
      elseif(_uname MATCHES "^(modflow6-)?v?([0-9]+\\.[0-9]+\\.[0-9]+)\\.(tar\\.gz|tgz|zip)$")
        set(RAVEN_MF6_LIB_VERSION "${CMAKE_MATCH_2}")
      endif()
    endif()
  else()
    set(_ep_src "${_ep_prefix}/src/mf6_source")
    set(_pin "")
    if(DEFINED RAVEN_MF6_PINNED_COMMIT_${RAVEN_MF6_VERSION})
      set(_pin "${RAVEN_MF6_PINNED_COMMIT_${RAVEN_MF6_VERSION}}")
      set(_rmf6_desc "built from source, git tag ${RAVEN_MF6_VERSION} of ${RAVEN_MF6_GIT_URL} (commit ${_pin}, checked)")
    else()
      message(WARNING "RAVEN_MF6: no pinned commit for MODFLOW ${RAVEN_MF6_VERSION}; the clone of tag ${RAVEN_MF6_VERSION} "
                      "is not checked")
      set(_rmf6_desc "built from source, git tag ${RAVEN_MF6_VERSION} of ${RAVEN_MF6_GIT_URL} (not checked)")
    endif()
    list(APPEND _ep_args DOWNLOAD_COMMAND "${CMAKE_COMMAND}" "-DGIT=${_rmf6_git}" "-DURL=${RAVEN_MF6_GIT_URL}"
         "-DTAG=${RAVEN_MF6_VERSION}" "-DCOMMIT=${_pin}" "-DDIR=${_ep_src}" -P "${_RMF6_DIR}/RavenMF6Verify.cmake")
    set(RAVEN_MF6_LIB_VERSION "${RAVEN_MF6_VERSION}")
  endif()
  string(APPEND _rmf6_desc " with ${_rmf6_fc}")
  if(_fcv)
    string(APPEND _rmf6_desc " (${_fcv})")
  endif()
  set(_rmf6_meson_build "${_ep_src}/meson.build")

  # meson: release build of the shared library only; the Fortran compiler passed as FC (recorded in the configure
  # command, mf6_source-cfgcmd.txt); the library's own run path is its folder, so that the Fortran runtime installed
  # beside it is found
  set(_env "${CMAKE_COMMAND}" -E env "FC=${_rmf6_fc}" "NINJA=${_rmf6_ninja}")
  set(_setup_args setup --reconfigure "${_ep_bin}" "${_ep_src}" --buildtype=release --wrap-mode=nodownload
      -Ddebug=false)
  if(APPLE)
    list(APPEND _setup_args "-Dfortran_link_args=-Wl,-rpath,@loader_path")
  elseif(NOT WIN32)
    list(APPEND _setup_args "-Dfortran_link_args=-Wl,-rpath,$ORIGIN")
  endif()
  set(_jobs "")
  if(NOT "${RAVEN_MF6_JOBS}" STREQUAL "")
    set(_jobs -j ${RAVEN_MF6_JOBS})
  endif()
  include(ExternalProject)
  ExternalProject_Add(mf6_source
    PREFIX "${_ep_prefix}"
    SOURCE_DIR "${_ep_src}"
    BINARY_DIR "${_ep_bin}"
    ${_ep_args}
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${_env} "${_rmf6_meson}" ${_setup_args}
    BUILD_COMMAND ${_env} "${_rmf6_meson}" compile -C "${_ep_bin}" ${_jobs} mf6:shared_library
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "${_ep_bin}/srcbmi/${RAVEN_MF6_LIBNAME}"
    USES_TERMINAL_BUILD TRUE)
  set(RAVEN_MF6_LIB_FILE "${_ep_bin}/srcbmi/${RAVEN_MF6_LIBNAME}")
  if(RAVEN_MF6_LIB_VERSION)
    set(_vshow "MODFLOW ${RAVEN_MF6_LIB_VERSION}")
  else()
    set(_vshow "MODFLOW version from the source's meson.build when built")
  endif()
  message(STATUS "RAVEN_MF6: SOURCE ${_vshow}, ${_rmf6_desc}; meson ${_rmf6_meson}, ninja ${_rmf6_ninja}")
endif()

# ---------------------------------------------------------------------------------------------------------------------
# 5. next to Raven at build time, with libmf6_version.txt (first line "MODFLOW X.Y.Z", as tools/get_mf6.py writes;
#    "MODFLOW unknown (given library)" for PATH)
# ---------------------------------------------------------------------------------------------------------------------
if(RAVEN_MF6_MODE STREQUAL "PATH")
  set(_vtext "unknown (given library)")
elseif(RAVEN_MF6_LIB_VERSION)
  set(_vtext "${RAVEN_MF6_LIB_VERSION}")
else()
  set(_vtext "")                                   # source build: read from meson.build when built
endif()
set(MF6_DEST "$<TARGET_FILE_DIR:Raven>")
add_custom_target(raven_mf6_lib ALL
  COMMAND "${CMAKE_COMMAND}" "-DSRC=${RAVEN_MF6_LIB_FILE}" "-DDEST=${MF6_DEST}" "-DNAME=${RAVEN_MF6_LIBNAME}"
          "-DVERSION=${_vtext}" "-DMESON_BUILD=${_rmf6_meson_build}" "-DDESC=${_rmf6_desc}"
          "-DSTAGE=${_rmf6_work}/libmf6_version.txt" -P "${_RMF6_DIR}/RavenMF6Copy.cmake"
  COMMENT "MODFLOW 6 library -> ${MF6_DEST}/${RAVEN_MF6_LIBNAME}"
  VERBATIM)
if(TARGET mf6_source)
  add_dependencies(raven_mf6_lib mf6_source)
endif()
