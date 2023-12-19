cmake_minimum_required(VERSION 3.5)
project(STACCATO)
set(CMAKE_BUILD_TYPE Debug)
set (SRCS
    src/DSDDecompose.c
    src/DSDInterface.c
    src/DSDManager.c
    src/DSDManipulations.c
    src/DSDNewDecompose.c
    src/DSDOrDecompose.c
    src/DSDPrimeDecompose.c
    src/DSDUtilities.c
    src/DSDXorDecompose.c
    src/fixheap.c
)

add_library(STACCATO SHARED ${SRCS})

target_include_directories(STACCATO
    PUBLIC
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd/cudd
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd/st
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd/mtr
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd/epd
    $ENV{UTOPIA_HOME}/lib/bill/lib/cudd/util
)
    
if(DEFINED ENV{CUDD_DIR})
    set(CUDD_LIBRARY_PATH "$ENV{CUDD_DIR}/libcudd.so")
else()
    set(CUDD_LIBRARY_PATH "/usr/local/lib/libcudd.so")
endif()

target_link_libraries(STACCATO PUBLIC ${CUDD_LIBRARY_PATH})
