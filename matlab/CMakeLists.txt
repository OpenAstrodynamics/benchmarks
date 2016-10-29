cmake_minimum_required(VERSION 3.3)
project(ICATT Fortran)

enable_language(CXX)
enable_testing()

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MAIN_PROGRAM)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -g")

add_subdirectory(../lib/dopri ${CMAKE_BINARY_DIR}/dopri)

matlab_add_mex(NAME propagator
    SRC propagator.f90 propagator_module.f90
    LINK_TO dopri_static)
matlab_add_unit_test(NAME proptest
    ADDITIONAL_PATH ${CMAKE_SOURCE_DIR}
    UNITTEST_FILE ${CMAKE_SOURCE_DIR}/proptest.m
    NO_UNITTEST_FRAMEWORK)

install(TARGETS propagator DESTINATION ${CMAKE_SOURCE_DIR})
