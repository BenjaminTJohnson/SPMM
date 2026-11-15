cmake_minimum_required(VERSION 3.20)
project(SPMM LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SPMM_NATIVE "Tune for the build host CPU (binary may not run on other nodes)" OFF)
option(SPMM_BUILD_LEGACY "Also build the deprecated program spmm_legacy (src/domain_melt5.f90)" OFF)

set(id ${CMAKE_Fortran_COMPILER_ID})
if(id STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -Wextra -fcheck=all -fbacktrace")
  set(native -march=native)
elseif(id STREQUAL "Intel" OR id STREQUAL "IntelLLVM")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -warn all -check all -traceback")
  set(native -xHost)
endif()

add_subdirectory(src)
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER} (${id}), build type ${CMAKE_BUILD_TYPE}")
