cmake_minimum_required(VERSION 2.6) 
######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)
if(MKL_OFF)
else()
	INCLUDE_DIRECTORIES(${FFTW_INCLUDES} ${MKL_INCLUDE_DIRS})
endif()

######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
if(MKL_OFF)
else()
	LINK_DIRECTORIES(${FFTW_LIBRARIES} ${MKL_LIBRARIES})
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
	#TARGET_LINK_LIBRARIES(LEO2D leo2d fftw3 mkl_intel_lp64 mkl_sequential mkl_core)
	TARGET_LINK_LIBRARIES(LEO2D leo2d fftw3 ${MKL_LIBRARIES})
elseif(MKL_INTEL)
	TARGET_LINK_LIBRARIES(LEO2D leo2d fftw3 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
else()
	TARGET_LINK_LIBRARIES(LEO2D leo2d)
endif()
######################################################################################
