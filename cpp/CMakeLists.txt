cmake_minimum_required(VERSION 3.3)
project(ICATT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)

add_subdirectory(../lib/dopri ${CMAKE_BINARY_DIR}/dopri)

include_directories(
        ${EIGEN3_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/src
        ../../deps/dopri
        ../../deps/spice/src/cspice/include)
link_directories(${PROJECT_SOURCE_DIR})

set(SOURCE_FILES src/main.cpp src/elements.cpp src/kepler.cpp src/cppdopri.cpp src/lambert.cpp)
add_executable(icatt ${SOURCE_FILES})
target_link_libraries(icatt dopri cspice)
