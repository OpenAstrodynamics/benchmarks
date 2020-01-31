cmake_minimum_required(VERSION 3.3)

project(dopri Fortran)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -Wall")

set(SOURCE_FILES
    dopri.f90
    dop853.f
    dopri5.f)

add_library(dopri SHARED ${SOURCE_FILES})
add_library(dopri_static STATIC ${SOURCE_FILES})
