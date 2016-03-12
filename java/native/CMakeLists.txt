cmake_minimum_required(VERSION 3.3)

add_subdirectory(../../../deps/dopri ${CMAKE_BINARY_DIR}/dopri)

if(DEFINED ENV{JAVA_HOME})
    set(JAVA_HOME ENV{JAVA_HOME})
else()
    set(JAVA_HOME /Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home/)
endif()

include_directories(${CMAKE_SOURCE_DIR}
    ${JAVA_HOME}/include
    ${JAVA_HOME}/include/darwin
    ${CMAKE_SOURCE_DIR}/../../lib/dopri)

set(SRCS jdopri.c)

add_library(jnilib SHARED ${SRCS})
target_link_libraries(jnilib dopri)
