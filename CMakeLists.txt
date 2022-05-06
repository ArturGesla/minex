cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(CMAKE_BUILD_TYPE Debug)

project(minex VERSION 1.0)

file(GLOB_RECURSE SRC_FILES .)

add_executable(minex main.cpp ${SRC_FILES})

set(INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/library/inc")
#list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/lean_vtk")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../eigen")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/UMFPACK/Include")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/AMD/Include")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/SuiteSparse_config")

#message (STATUS "${INC_DIRS}")


target_include_directories(minex PUBLIC ${INC_DIRS})

#umfpack search
find_library(UMFPACK_DIR NAMES libumfpack.so libumfpack.so.5)
message (STATUS "umfpack found at: ${UMFPACK_DIR}")
target_link_libraries(minex ${UMFPACK_DIR})



