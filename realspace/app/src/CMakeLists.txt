cmake_minimum_required(VERSION 2.6) # Pour éviter un warning au lancement de CMake, on doit préciser la version minimale
######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)
if(MKL_GNU)
	INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
endif()

######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
if(MKL_GNU)
	LINK_DIRECTORIES(${MKL_LIBRARIES})
endif()

######################################################################################
### 3) Define your target: just an executable here
######################################################################################

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/app)

ADD_EXECUTABLE(LEO2D main.cpp)
ADD_DEPENDENCIES(LEO2D leo2d)


######################################################################################
### 4) Link the librairies for your executable
######################################################################################

if(MKL_GNU)
	TARGET_LINK_LIBRARIES(LEO2D leo2d mkl_intel_lp64 mkl_sequential mkl_core)
else()
	TARGET_LINK_LIBRARIES(LEO2D leo2d)
endif()
######################################################################################
