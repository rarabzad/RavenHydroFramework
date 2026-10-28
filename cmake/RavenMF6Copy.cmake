# RavenMF6Copy.cmake: puts the MODFLOW 6 library next to the Raven executable (build time; target raven_mf6_lib).
#
#   cmake -DSRC=<library> -DDEST=<folder of Raven> -DNAME=<libmf6.so|.dylib|.dll> -DVERSION=<X.Y.Z, text or empty>
#         -DMESON_BUILD=<meson.build of a source build or empty> -DDESC=<where the library came from>
#         -DSTAGE=<build-tree copy of libmf6_version.txt, installed with the library> -P RavenMF6Copy.cmake
#
# The library is copied under its standard name (a library given as libmf6_663.so or libmf6.so.6.8.1 becomes libmf6.so),
# and libmf6_version.txt is written beside it in the format of tools/get_mf6.py:
#   MODFLOW 6.8.1
#   source: <where it came from>
#   sha256: <SHA256 of the library>
cmake_minimum_required(VERSION 3.20)
if(NOT EXISTS "${SRC}")
  message(FATAL_ERROR "MODFLOW 6 library ${SRC} not found")
endif()
file(MAKE_DIRECTORY "${DEST}")
set(_dst "${DEST}/${NAME}")
file(SHA256 "${SRC}" _sha)
set(_same FALSE)
if(EXISTS "${_dst}")
  file(SHA256 "${_dst}" _dsha)
  if(_dsha STREQUAL _sha)
    set(_same TRUE)
  endif()
endif()
if(NOT _same)
  file(REMOVE "${_dst}")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E copy "${SRC}" "${_dst}" RESULT_VARIABLE _r)   # (file(COPY_FILE) needs 3.21)
  if(NOT _r EQUAL 0)
    message(FATAL_ERROR "cannot copy ${SRC} to ${_dst}")
  endif()
  message(STATUS "MODFLOW 6 library: ${_dst}")
endif()

set(_v "${VERSION}")
if("${_v}" STREQUAL "" AND NOT "${MESON_BUILD}" STREQUAL "" AND EXISTS "${MESON_BUILD}")
  file(READ "${MESON_BUILD}" _t LIMIT 4000)
  if(_t MATCHES "project\\([^)]*version[ \t\r\n]*:[ \t\r\n]*'([0-9][0-9A-Za-z.+-]*)'")
    set(_v "${CMAKE_MATCH_1}")
  endif()
endif()
if("${_v}" STREQUAL "")
  set(_v "unknown")
endif()
set(_text "MODFLOW ${_v}\nsource: ${DESC}\nsha256: ${_sha}\n")
foreach(_f "${DEST}/libmf6_version.txt" "${STAGE}")
  if(NOT "${_f}" STREQUAL "")
    set(_old "")
    if(EXISTS "${_f}")
      file(READ "${_f}" _old)
    endif()
    if(NOT _old STREQUAL _text)
      file(WRITE "${_f}" "${_text}")
    endif()
  endif()
endforeach()
