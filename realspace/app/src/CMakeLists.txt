cmake_minimum_required(VERSION 2.6)
######################################################################################
### 1) Include the sources
######################################################################################


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/)



if (LINK_FFTW)
	INCLUDE_DIRECTORIES(${FFTW_INCLUDES})
endif()

if(MKL_OFF)
else()
	INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
endif()

######################################################################################
### 2) Specify where CMake can find the libraries
######################################################################################

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

if(MKL_OFF)
else()
	LINK_DIRECTORIES(${MKL_LIBRARIES})
endif()


######################################################################################
### 3) Define your target: just an executable here
######################################################################################

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/app)

ADD_EXECUTABLE(LEO2D main.cpp)
#ADD_EXECUTABLE(LEO2D_serial main_serial.cpp)

ADD_DEPENDENCIES(LEO2D leo2d)
#ADD_DEPENDENCIES(LEO2D_serial leo2d)


######################################################################################
### 4) Link the librairies for your executable
######################################################################################

if (LINK_FFTW)
	if(MKL_GNU)
		#TARGET_LINK_LIBRARIES(LEO2D leo2d fftw3 mkl_intel_lp64 mkl_sequential mkl_core)
		TARGET_LINK_LIBRARIES(LEO2D leo2d ${FFTW_LIBRARIES} ${MKL_LIBRARIES})
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d fftw3 ${MKL_LIBRARIES})
	elseif(MKL_INTEL)
		TARGET_LINK_LIBRARIES(LEO2D leo2d ${FFTW_LIBRARIES} mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d fftw3 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
	else()
		TARGET_LINK_LIBRARIES(LEO2D ${FFTW_LIBRARIES}  leo2d)
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d)
	endif()
else()
	if(MKL_GNU)
		#TARGET_LINK_LIBRARIES(LEO2D leo2d fftw3 mkl_intel_lp64 mkl_sequential mkl_core)
		TARGET_LINK_LIBRARIES(LEO2D leo2d ${MKL_LIBRARIES})
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d fftw3 ${MKL_LIBRARIES})
	elseif(MKL_INTEL)
		TARGET_LINK_LIBRARIES(LEO2D leo2d mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d fftw3 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
	else()
		TARGET_LINK_LIBRARIES(LEO2D leo2d)
		#TARGET_LINK_LIBRARIES(LEO2D_serial leo2d)
	endif()
endif()
######################################################################################
