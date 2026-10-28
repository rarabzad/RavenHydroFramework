# RavenMF6Verify.cmake: download step of the MODFLOW 6 source build (RAVEN_MF6=SOURCE without RAVEN_MF6_SOURCE_DIR/URL).
#
#   cmake -DGIT=<git> -DURL=<repository> -DTAG=<release tag> -DCOMMIT=<pinned commit or empty> -DDIR=<folder>
#         -P RavenMF6Verify.cmake
#
# Shallow clone of the release tag (the tag archive could not be pinned: the GitHub tag tarball is generated on demand),
# then the checked-out commit is compared with the pinned one; a clone that does not match is removed and the build stops.
cmake_minimum_required(VERSION 3.20)
foreach(_v GIT URL TAG DIR)
  if("${${_v}}" STREQUAL "")
    message(FATAL_ERROR "RavenMF6Verify.cmake: ${_v} not given")
  endif()
endforeach()

function(_head out)
  execute_process(COMMAND "${GIT}" -C "${DIR}" rev-parse HEAD OUTPUT_VARIABLE _h RESULT_VARIABLE _r
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT _r EQUAL 0)
    set(_h "")
  endif()
  set(${out} "${_h}" PARENT_SCOPE)
endfunction()

if(EXISTS "${DIR}/.git")
  _head(_have)
  if(COMMIT AND _have STREQUAL COMMIT)
    message(STATUS "MODFLOW 6 source: ${DIR} already at tag ${TAG} (commit ${COMMIT})")
    return()
  endif()
endif()
file(REMOVE_RECURSE "${DIR}")
get_filename_component(_parent "${DIR}" DIRECTORY)
file(MAKE_DIRECTORY "${_parent}")
message(STATUS "MODFLOW 6 source: git clone --depth 1 --branch ${TAG} ${URL}")
execute_process(COMMAND "${GIT}" -c advice.detachedHead=false clone --depth 1 --branch "${TAG}" "${URL}" "${DIR}"
                RESULT_VARIABLE _r)
if(NOT _r EQUAL 0)
  file(REMOVE_RECURSE "${DIR}")
  message(FATAL_ERROR "MODFLOW 6 source: git clone of tag ${TAG} from ${URL} failed. Offline: give a source folder "
                      "(-DRAVEN_MF6_SOURCE_DIR=...) or a release zip as source (-DRAVEN_MF6_SOURCE_URL=file:///.../mf${TAG}_linux.zip), "
                      "or use -DRAVEN_MF6=PREBUILT with -DRAVEN_MF6_ZIP=...")
endif()
_head(_have)
if(COMMIT)
  if(NOT _have STREQUAL COMMIT)
    file(REMOVE_RECURSE "${DIR}")
    message(FATAL_ERROR "MODFLOW 6 source: tag ${TAG} of ${URL} is commit ${_have}, but Raven pins commit ${COMMIT}: "
                        "not the official release. Check RAVEN_MF6_GIT_URL.")
  endif()
  message(STATUS "MODFLOW 6 source: tag ${TAG} is the pinned commit ${COMMIT}")
else()
  message(STATUS "MODFLOW 6 source: tag ${TAG} is commit ${_have} (no pinned commit; not checked)")
endif()
