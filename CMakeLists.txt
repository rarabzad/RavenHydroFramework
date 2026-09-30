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
option(RAVEN_GET_MF6 "If ON, downloads the MODFLOW 6 library (libmf6) for the groundwater coupling into lib/mf6 (default: ON)" ON)
set(RAVEN_MF6_VERSION "6.8.1" CACHE STRING "MODFLOW 6 release fetched when RAVEN_GET_MF6 is ON (6.8.1 is the tested one)")

# Setup Project
PROJECT(Raven CXX)

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

# MODFLOW 6 library for the groundwater coupling (loaded at run time; Raven builds and runs without it). Fetched once
# from the MODFLOW 6 releases into lib/mf6 (as tools/get_mf6.py does) and copied next to the executable, where Raven
# finds it without any setting. A failed download only warns: the library can be fetched later with tools/get_mf6.py.
if(WIN32)
  set(MF6_LIBNAME "libmf6.dll")
  set(MF6_ASSET "win64")
elseif(APPLE)
  set(MF6_LIBNAME "libmf6.dylib")
  if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    set(MF6_ASSET "macarm")
  else()
    set(MF6_ASSET "mac")
  endif()
else()
  set(MF6_LIBNAME "libmf6.so")
  set(MF6_ASSET "linux")
endif()
set(MF6_DIR "${CMAKE_SOURCE_DIR}/lib/mf6")
if(RAVEN_GET_MF6 AND NOT EXISTS "${MF6_DIR}/${MF6_LIBNAME}")
  set(MF6_URL "https://github.com/MODFLOW-ORG/modflow6/releases/download/${RAVEN_MF6_VERSION}/mf${RAVEN_MF6_VERSION}_${MF6_ASSET}.zip")
  set(MF6_ZIP "${CMAKE_BINARY_DIR}/_mf6/mf${RAVEN_MF6_VERSION}_${MF6_ASSET}.zip")
  message(STATUS "Downloading the MODFLOW 6 library: ${MF6_URL}")
  file(DOWNLOAD "${MF6_URL}" "${MF6_ZIP}" STATUS MF6_STATUS TLS_VERIFY ON)
  list(GET MF6_STATUS 0 MF6_CODE)
  if(MF6_CODE EQUAL 0)
    file(ARCHIVE_EXTRACT INPUT "${MF6_ZIP}" DESTINATION "${CMAKE_BINARY_DIR}/_mf6/x")
    file(GLOB_RECURSE MF6_FOUND "${CMAKE_BINARY_DIR}/_mf6/x/*/${MF6_LIBNAME}")
    if(MF6_FOUND)
      list(GET MF6_FOUND 0 MF6_SRC)
      file(MAKE_DIRECTORY "${MF6_DIR}")
      file(COPY "${MF6_SRC}" DESTINATION "${MF6_DIR}")
      file(WRITE "${MF6_DIR}/libmf6_version.txt" "MODFLOW ${RAVEN_MF6_VERSION}\nsource: ${MF6_URL}\n")
      message(STATUS "MODFLOW ${RAVEN_MF6_VERSION} library: ${MF6_DIR}/${MF6_LIBNAME}")
    else()
      message(WARNING "${MF6_URL} holds no ${MF6_LIBNAME}")
    endif()
  else()
    message(WARNING "Could not download the MODFLOW 6 library (${MF6_STATUS}). Raven builds anyway; for the groundwater "
                    "coupling run  python tools/get_mf6.py  or set RAVEN_MF6_LIB.")
  endif()
  file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_mf6")
endif()
if(COMPILE_EXE AND EXISTS "${MF6_DIR}/${MF6_LIBNAME}")
  add_custom_command(TARGET Raven POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MF6_DIR}/${MF6_LIBNAME}" "$<TARGET_FILE_DIR:Raven>")
endif()

# unset cmake variables to avoid polluting the cache
unset(COMPILE_LIB CACHE)
unset(COMPILE_EXE CACHE)
