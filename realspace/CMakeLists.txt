# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.5.1)

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

SET(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-std=c++1z -Wno-unused-variable -Wno-unused-parameter")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/module/")

find_package(FFTW)
find_package(MPI)
find_package(Eigen3 3.3 NO_MODULE)

# include required directories
include_directories(${FFTW_INCLUDES} ${MPI_INCLUDE_PATH} ${EIGEN3_INCLUDES})
include(CMakeBackwardCompatibilityCXX)
include(CheckLibraryExists)

# link necessary subdirectory library paths

add_subdirectory(src)
add_subdirectory(app)
