cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(CMAKE_BUILD_TYPE Debug)

project(minex VERSION 1.0)

file(GLOB_RECURSE SRC_FILES .)

add_executable(minex main.cpp ${SRC_FILES})

set(INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
# list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/lean_vtk")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../eigen")
#  list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/UMFPACK/Include")
#  list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/AMD/Include")
#  list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/SuiteSparse_config")
 list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../superlu/SRC")
#  list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../arpackpp/include")
#  list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../arpackpp/examples/matrices/nonsym")
 list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/inc/arpackpp-for-eigen/include")
 list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/inc/arpackpp-for-eigen/examples/matrices/nonsym")




#message (STATUS "${INC_DIRS}")


target_include_directories(minex PUBLIC ${INC_DIRS})

# lapack search
find_library(LAPACK_DIR NAMES openblas) #lapack is built in
message (STATUS "superlu found at: ${LAPACK_DIR}")
target_link_libraries(minex ${LAPACK_DIR})

# superlu search
find_library(SUPERLU_DIR NAMES libsuperlu.so libsuperlu.so.5)
message (STATUS "superlu found at: ${SUPERLU_DIR}")
target_link_libraries(minex ${SUPERLU_DIR})

# arpack search (probably in some local dir)
set(LOC_LIB_DIRS "/people/gesla/Documents/git/others/my-libs/usr/lib/x86_64-linux-gnu")
find_library(ARPACK_DIR NAMES arpack PATHS ${LOC_LIB_DIRS})
message (STATUS "superlu found at: ${ARPACK_DIR}")
target_link_libraries(minex ${ARPACK_DIR})



