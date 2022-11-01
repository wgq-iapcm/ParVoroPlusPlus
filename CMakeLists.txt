####################################################################
#
# This a CMake project file that controls which source files must
# be built for the executable and which third-party libraries it 
# needs to be linked with.
#
# Use the cmake utility to create a Makefile from this project file.
# See http://www.cmake.org/ for details.
#
####################################################################

PROJECT(Parvoro++)
CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -std=c++11")
MESSAGE(CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS})

###################################################################

SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)
# CMAKE_BUILD_TYPE defaults to Debug
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
ENDIF()

##################################################################

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# MPI
IF(NOT DEFINED MPI_ROOT)
    SET(MPI_ROOT /home/wgq/local/mpich) #notebook
    #SET(MPI_ROOT /usr/local/mpi3_gnu) 
    #SET(MPI_ROOT /usr/local/mpi3-gcc830-shared)
    MESSAGE(STATUS "MPI_ROOT not set, use ${MPI_ROOT} as default path.  You can change the value with cmake -DMPI_ROOT=<dir>")
ENDIF()
FIND_PACKAGE(MPI REQUIRED)
IF(MPI_FOUND)
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_DIRS})
    SET(libraries ${libraries} ${MPI_CXX_LIBRARIES})
    ADD_DEFINITIONS(-DPARALLEL)
ELSE()
    MESSAGE(FATAL_ERROR "MPI library not found")
ENDIF()

#Zlib
IF(NOT DEFINED zlib_PATH)
	SET(zlib_PATH /usr/) #notebook 
	MESSAGE(WARNING "zlib_PATH not set, use ${zlib_PATH} as default path.  You can change the value with cmake -Dzlib_PATH=<dir>")
ENDIF()
FIND_PACKAGE(zlib)
IF(zlib_FOUND)
	INCLUDE_DIRECTORIES(${zlib_INCLUDES})
	SET(libraries ${libraries} ${zlib_LIBRARIES})
	ADD_DEFINITIONS(-DHAVE_ZLIB)
ELSE()
	MESSAGE(FATAL_ERROR "zlib library not found. Please install this library on your system or specify the location of the library manually.")
ENDIF()

#HDF5
IF(NOT DEFINED HDF5_ROOT)
	SET(HDF5_ROOT /home/wgq/local/hdf5) #notebook
	#SET(HDF5_ROOT /vol8/home/jsg045819y/local/hdf5) 
	MESSAGE(STATUS "HDF5_ROOT not set, use ${HDF5_ROOT} as default path.  You can change the value with cmake -DHDF5_ROOT=<dir>")
ENDIF()
FIND_PACKAGE(HDF5 REQUIRED)
IF(NOT HDF5_FOUND)
	MESSAGE("HDF5 library not found.")
	return()
ELSE()
	MESSAGE(STATUS "HDF5_VERSION: ${HDF5_VERSION}")
	INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
	SET(libraries ${libraries} ${HDF5_LIBRARIES} -lz)
	ADD_DeFINITIONS(-DHAVE_HDF5)
	ADD_DeFINITIONS(-DH5_HAVE_PARALLEL)
ENDIF()

#Silo
IF(NOT DEFINED silo_PATH)
	SET(silo_PATH /home/wgq/local/silo) #notebook
	#SET(silo_PATH /vol8/home/jsg045819y/local/silo) 
	#SET(silo_PATH /public/home/ch160458/local/silo) 
	MESSAGE(WARNING "silo_PATH not set, use ${silo_PATH} as default path.  You can change the value with cmake -Dsilo_PATH=<dir>")
ENDIF()
FIND_PACKAGE(silo)
IF(silo_FOUND)
	INCLUDE_DIRECTORIES(${silo_INCLUDES})
	SET(libraries ${libraries} ${silo_LIBRARIES})
	ADD_DEFINITIONS(-DHAVE_SILO)
ELSE()
	MESSAGE(FATAL_ERROR "silo library not found. Please install this library on your system or specify the location of the library manually.")
ENDIF()

####################################################################

# Define install
MESSAGE(STATUS)
MESSAGE(STATUS "Install Prefix (CMAKE_INSTALL_PREFIX):     ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS "To override install options, add -D{OPTION_NAME}=... to the cmake command")
MESSAGE(STATUS "To build and install MDanalysis, enter \"make\" and \"make install\"")
MESSAGE(STATUS)

# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Include install
INSTALL (DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

###################################################################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Subdirectories
ADD_SUBDIRECTORY	    (3rdparty)
SET(libraries ${libraries} diy)
SET(libraries ${libraries} VoroPlusPlus)

ADD_SUBDIRECTORY	    (bin)

###################################################################





