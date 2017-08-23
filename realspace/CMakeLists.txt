# Specify the minimum version for CMake

cmake_minimum_required(VERSION 2.4.0)

# Project's name
PROJECT(LEO2D)
SET(CMAKE_BUILD_TYPE Debug)

# Define name and version of package
SET(PACKAGE_NAME "LEO2D" CACHE STRING "package_name" FORCE)
SET(PACKAGE_VERSION "1.0" CACHE STRING "package_version" FORCE)

# set code language
ENABLE_LANGUAGE(CXX)
SET(CMAKE_CXX_COMPILER mpicxx)

# setup packages

# (WITH MKL)
# SET(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-std=c++11 -mkl -DUSE_MKL -Wno-unused-variable -Wno-unused-parameter")

# (WITHOUT MKL)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-variable -Wno-unused-parameter")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/module/")

#MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )

find_package(FFTW REQUIRED)
find_package(MPI)
find_package(Eigen3 REQUIRED)

# include required directories
include_directories(${FFTW_INCLUDES} ${MPI_INCLUDE_PATH} ${EIGEN3_INCLUDE_DIR})
include(CMakeBackwardCompatibilityCXX)
include(CheckLibraryExists)

# link necessary subdirectory library paths

add_subdirectory(src)
add_subdirectory(app)
