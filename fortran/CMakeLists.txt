cmake_minimum_required(VERSION 3.3)

project(icatt Fortran)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

add_subdirectory(../../deps/dopri ${PROJECT_BINARY_DIR}/dopri)

set(SRCS src/main.f90
    src/elements.f90
    src/fdopri.f90
    src/kepler.f90
    src/newton.f90
    src/lambert.f90)

add_executable(icatt ${SRCS})
target_link_libraries(icatt dopri)
