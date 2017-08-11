#
# Project settings
#

# Minimum required version of cmake
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

# Project name and programming languages used
PROJECT( MY_COOL_PROJECT CXX )

# CMake module(s) path
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR} )

#
# Prerequisites
#
FIND_PACKAGE( ACADO REQUIRED )

#
# Include directories
#
INCLUDE_DIRECTORIES( . ${ACADO_INCLUDE_DIRS} )

# Use C++11
# ADD_DEFINITIONS(-std=C++11)

#
# Build an executable
#
ADD_EXECUTABLE(        my_exe my_code )
# Linux and OS X
TARGET_LINK_LIBRARIES( my_exe ${ACADO_SHARED_LIBRARIES} )
# Windows
#TARGET_LINK_LIBRARIES( my_exe ${ACADO_STATIC_LIBRARIES} )
SET_TARGET_PROPERTIES( my_exe
                       PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )