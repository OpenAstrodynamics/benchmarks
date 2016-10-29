cmake_minimum_required(VERSION 3.3)
project(ICATT C)

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")

include_directories(${PROJECT_SOURCE_DIR}/../../deps/dopri)

add_subdirectory(../../deps/dopri ${CMAKE_BINARY_DIR}/dopri)

matlab_add_mex(NAME propagator
    SRC propagator.c
    LINK_TO dopri)
