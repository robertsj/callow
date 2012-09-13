#------------------------------------------------------------------------------#
# callow - a naive sparse linear algebra library
#------------------------------------------------------------------------------#

project(callow)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#------------------------------------------------------------------------------#
# OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: Debug Release."
    FORCE)
endif()

# Options for Python bindings. On or Off.
SET(CALLOW_ENABLE_PYTHON NO CACHE BOOL
  "Enable Python bindings via SWIG."
)

# Options for OpenMP directives. On or Off.
SET(CALLOW_ENABLE_OPENMC NO CACHE BOOL
  "Enable OpenMP for multithreaded solves."
)

# Write config file, install it, and include that directory for all source
set(CALLOW_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/config/callow_config.hh.in")
set(CALLOW_CONFIG_OUT "${PROJECT_BINARY_DIR}/config/callow_config.hh")
configure_file(${CALLOW_CONFIG_IN} ${CALLOW_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/config")
install(FILES ${PROJECT_BINARY_DIR}/config/callow_config.hh
        DESTINATION include)

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(config)
include_directories(utils)
include_directories(vector)
include_directories(matrix)
include_directories(solver)
add_subdirectory(utils)
add_subdirectory(vector)
add_subdirectory(matrix)
add_subdirectory(solver)
if(CALLOW_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

#------------------------------------------------------------------------------#
# TESTING
#------------------------------------------------------------------------------#

OPTION (USE_TESTING "Build testing" ON)
if (CALLOW_ENABLE_TEST)
    enable_testing()
    include(CTest)
    include_directories(test)
    add_subdirectory(test)
endif()

#------------------------------------------------------------------------------#
# FINAL LIBRARY
#------------------------------------------------------------------------------#

ADD_LIBRARY(callow
            SHARED
            ${UTILS_SRC}
            ${VECTOR_SRC}
            ${MATRIX_SRC}
            ${SOLVER_SRC}
)
#ADD_LIBRARY(callow_static
#            STATIC
#            ${SRC} 
#)

#------------------------------------------------------------------------------#
# INSTALL
#------------------------------------------------------------------------------#

if( NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
      "Installation prefix; default is build directory."
      FORCE)
endif()

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("---------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION")
message("---------------------------------------------------------------")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")

message("---------------------------------------------------------------")
if(CALLOW_ENABLE_PYTHON)
message("++++ Python:      enabled")
else()
message("++++ Python:      disabled")
endif()
if(CALLOW_ENABLE_OPENMP)
message("++++ OpenMP:      enabled")
else()
message("++++ OpenMP:      disabled")
endif()
if(CALLOW_ENABLE_TEST )
message("++++ Testing:     enabled")
else()
message("++++ Testing:     disabled")
endif()
message("---------------------------------------------------------------")


