cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing
# modified for GitHub and C++11 by Trevor James Smith of Ouranos in Montreal
# modified for lpsolve support by Maxim Krassovski of DFO

# Set the C++ standard to C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# optional cmake command line arguments (e.g, "cmake -D COMPILE_LIB=ON" .)
option(COMPILE_LIB "If ON, will create a dynamic lib file (default: OFF)" OFF)
option(COMPILE_EXE "If ON, will create a executable file (default: ON)" ON)
option(PYTHON, "If ON, will create a share library for python (default: OFF)" OFF)
option(LPSOLVE, "If ON, will link to lp_solve optimization library (default: OFF)" OFF)

# Setup Project
PROJECT(Raven CXX)

# Single-configuration generators (Unix Makefiles, Ninja) ignore "cmake --build --config Release": without a build type
# Raven would be compiled without optimisation (about 3 times slower). Default: Release; -DCMAKE_BUILD_TYPE=... overrides
# it. Multi-configuration generators (Visual Studio, Xcode, Ninja Multi-Config) choose with --config as before.
get_property(_raven_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _raven_multi_config AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release (default), Debug, RelWithDebInfo or MinSizeRel" FORCE)
  message(STATUS "CMAKE_BUILD_TYPE not given: Release")
endif()

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find NetCDF
find_package(NetCDF) #may also be 'netCDF'
if (NOT NetCDF_FOUND)
  message(STATUS "NetCDF not found, trying netCDF")
  find_package(netCDF)
endif()

# Find HDF5
find_package(HDF5)

# find header & source
file(GLOB HEADER "src/*.h")
file(GLOB SOURCE "src/*.cpp")

# Create library with Python bindings
# To work, install pybind11, making sure it comes with cmake files (conda install -c conda-forge pybind11)
if(PYTHON)
  SET(PYBIND11_NEWPYTHON ON)
  find_package(Python COMPONENTS REQUIRED Interpreter Development)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(libraven MODULE src/py/libraven.cpp)
  include_directories("src")
  target_compile_features(libraven PUBLIC cxx_std_11)
endif()

# creates a shared library - file extension is OS dependent (Linux: .so, Windows: .dll)
if(COMPILE_LIB)
  add_library(ravenbmi SHARED ${SOURCE})
  target_compile_definitions(ravenbmi PUBLIC BMI_LIBRARY)
  target_link_libraries(ravenbmi ${CMAKE_DL_LIBS})
endif()

# creates an executable - file extension is OS dependent (Linux: none, Windows: .exe)
if(COMPILE_EXE)
  add_executable(Raven
    ${SOURCE}
    ${HEADER}
  )
  set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)

  # MODFLOW 6 groundwater coupling: libmf6 is loaded at run time (dlopen on Linux/macOS)
  target_link_libraries(Raven ${CMAKE_DL_LIBS})

  # Remove deprecation warnings for GCC
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(Raven PRIVATE "-Wno-deprecated")
    message(STATUS "Modified compile flags with '-Wno-deprecated'")
  ENDIF()

  if(LPSOLVE)
    target_link_directories(Raven PRIVATE lib/lp_solve)  # where liblpsolve55.so is located
    target_link_libraries(Raven lpsolve55)
    add_definitions(-D_LPSOLVE_)
  endif()
endif()

if(NETCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven NetCDF::NetCDF)
elseif(netCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven netcdf)
else()
  message(STATUS "NetCDF not found")
endif()

source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

# MODFLOW 6 library for the groundwater coupling, placed next to Raven (cmake/RavenMF6.cmake):
#   -DRAVEN_MF6=SOURCE|PREBUILT|PATH|OFF (default: SOURCE when a Fortran compiler, meson >= 1.7, ninja and git are found
#   and the MODFLOW 6 git repository answers at configure time, else PREBUILT when MODFLOW releases a library for this
#   platform and it can be downloaded, else OFF with a warning). With OFF the compile and link commands are the same as
#   without this block, and cmake --install installs Raven alone.
if(TARGET Raven)
  include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/RavenMF6.cmake")
  # cmake --install / cpack: one bin folder with Raven, and unless RAVEN_MF6 is OFF libmf6 (always under its standard
  # name), libmf6_version.txt and the Fortran runtime the library needs
  install(TARGETS Raven RUNTIME DESTINATION bin COMPONENT Raven)
  if(RAVEN_MF6_MODE STREQUAL "OFF")
    message(STATUS "RAVEN_MF6: OFF, so cmake --install installs Raven alone (bin/Raven, no MODFLOW 6 library); a coupled "
                   "model then needs one at run time (:MF6Library, RAVEN_MF6_LIB, or libmf6 next to Raven)")
    install(CODE "message(STATUS \"RAVEN_MF6 was OFF: Raven installed without the MODFLOW 6 library (a coupled model needs \"
                 \"one at run time: :MF6Library, RAVEN_MF6_LIB, or libmf6 next to Raven)\")" COMPONENT Raven)
  else()
    raven_mf6_rpath(Raven)

    install(PROGRAMS "${RAVEN_MF6_LIB_FILE}" DESTINATION bin RENAME "${RAVEN_MF6_LIBNAME}" COMPONENT mf6)
    install(FILES "${CMAKE_BINARY_DIR}/mf6/libmf6_version.txt" DESTINATION bin COMPONENT mf6)
    install(CODE "set(RAVEN_MF6_BUNDLE_LIB \"${RAVEN_MF6_LIB_FILE}\")
      set(RAVEN_MF6_BUNDLE_FC \"${RAVEN_MF6_FORTRAN_USED}\")
      set(RAVEN_MF6_BUNDLE_SYSTEM \"${CMAKE_SYSTEM_NAME}\")
      set(RAVEN_MF6_BUNDLE_NAME \"${RAVEN_MF6_LIBNAME}\")
      include(\"${CMAKE_CURRENT_SOURCE_DIR}/cmake/RavenMF6Bundle.cmake\")" COMPONENT mf6)

    # package Raven-<Raven version>-mf6-<MODFLOW version>-<platform> (no MODFLOW version when it is unknown)
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/src/RavenInclude.h" _raven_ver REGEX "__RAVEN_VERSION__ *=")
    string(REGEX REPLACE ".*\"([^\"]*)\".*" "\\1" _raven_ver "${_raven_ver}")
    _rmf6_target_arch(_raven_arch)
    string(TOLOWER "${_raven_arch}" _raven_arch)
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      set(_raven_sys windows)
      if(_raven_arch MATCHES "^(x64|amd64|x86_64)$")
        set(_raven_arch x64)
      endif()
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      set(_raven_sys macos)
    else()
      string(TOLOWER "${CMAKE_SYSTEM_NAME}" _raven_sys)
    endif()
    if(RAVEN_MF6_LIB_VERSION)
      set(_raven_mf6tag "-mf6-${RAVEN_MF6_LIB_VERSION}")
    else()
      set(_raven_mf6tag "")
    endif()
    set(CPACK_PACKAGE_NAME Raven)
    set(CPACK_PACKAGE_VENDOR "University of Waterloo")
    set(CPACK_PACKAGE_VERSION "${_raven_ver}")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Raven hydrological modelling framework with the MODFLOW 6 library")
    set(CPACK_PACKAGE_FILE_NAME "Raven-${_raven_ver}${_raven_mf6tag}-${_raven_sys}-${_raven_arch}")
    if(WIN32)
      set(CPACK_GENERATOR ZIP)
    else()
      set(CPACK_GENERATOR TGZ)
    endif()
    message(STATUS "RAVEN_MF6: package ${CPACK_PACKAGE_FILE_NAME}")
    include(CPack)
  endif()
endif()

# unset cmake variables to avoid polluting the cache
unset(COMPILE_LIB CACHE)
unset(COMPILE_EXE CACHE)
