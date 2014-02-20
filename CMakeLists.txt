#
#  Copyright 2013 CNRS
#  Author: Nirmal Giftsun
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# laas part
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/lapack.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

# set up project
SET(PROJECT_NAME sot_solver)
project(${PROJECT_NAME})
SET(PROJECT_VERSION 1.0)
SET(PROJECT_DESCRIPTION "the best solver ever written")
SET(PROJECT_URL "")
SET(${PROJECT_NAME}_HEADERS
    include/sot_solver/sot_solver.hh
)
SET(CXX_DISABLE_WERROR True)

SETUP_PROJECT()
SEARCH_FOR_EIGEN()
#set
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/src)
set(CMAKE_INSTALL_RPATH "${LIBRARY_OUTPUT_PATH}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_PATH}")

#  Dependencies
ADD_REQUIRED_DEPENDENCY("dynamic-graph")
ADD_REQUIRED_DEPENDENCY("dynamic-graph-python")



ADD_SUBDIRECTORY(src)

PKG_CONFIG_APPEND_CFLAGS(${_Eigen_CFLAGS})
SETUP_PROJECT_FINALIZE()



