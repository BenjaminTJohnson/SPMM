# CMakeLists files in this project can
cmake_minimum_required (VERSION 3.0.0)
project (SPMM LANGUAGES Fortran )

# Display the Fortran compiler informatio
set(CMAKE_Fortran_COMPILER "ifort")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ipo -xHost -qopt-zmm-usage=high -fp-model fast=2 -align array64byte -unroll")

message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")



add_subdirectory (src)




