##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013-2015, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Builds the seissol_kernels unit tests.

cmake_minimum_required( VERSION 2.8 )
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#option( GENERATED_KERNELS "Generated kernels to use: swsm, dwsm, ssnb, dsnb, sknc, dknc, shsw, dhsw" )
#option( CONVERGENCE_ORDER "Convergence order of the ADER-DG scheme" )
#option( GENERATED_CODE    "Generated code to use" )
option( SEISSOL_PROXY     "Location to the the seissol proxy" )

macro (test_cache_enum var)
	get_property (value CACHE ${var} PROPERTY VALUE)
	get_property (enum CACHE ${var} PROPERTY STRINGS)
	list (FIND enum "${value}" val_in_enum)
	if (val_in_enum LESS 0)
		message (FATAL_ERROR "Set ${var} to one of the following: ${enum}")
	endif ()
endmacro ()

set (GENERATED_KERNELS "" CACHE STRING "Generated kernels to use")
set_property (CACHE GENERATED_KERNELS PROPERTY STRINGS swsm dwsm ssnb dsnb sknc dknc shsw dhsw)
test_cache_enum (GENERATED_KERNELS)

set (CONVERGENCE_ORDER "" CACHE STRING "Convergence order")
set_property (CACHE CONVERGENCE_ORDER PROPERTY STRINGS 2 3 4 5 6 7 8)
test_cache_enum (CONVERGENCE_ORDER)

set (GENERATED_CODE "" CACHE PATH "Path to generated code")
if (NOT EXISTS "${GENERATED_CODE}/matrix_kernels/")
	message("Could not find ${GENERATED_CODE}/matrix_kernels/. Did you set GENERATED_CODE?")
endif ()


# set convergence order
add_definitions( -DCONVERGENCE_ORDER=${CONVERGENCE_ORDER} )

# set number of quantities
add_definitions( -DNUMBER_OF_QUANTITIES=9 )

# set include to generated code
if( GENERATED_CODE )
  include_directories( ${GENERATED_CODE} )
endif()

# split by architecture and precision
string( SUBSTRING ${GENERATED_KERNELS} 0 1  PRECISION_PLAIN )
string( SUBSTRING ${GENERATED_KERNELS} 1 -1 ARCHITECTURE_PLAIN )

# build the architecture dependent matrix kernels as library
if( ${ARCHITECTURE_PLAIN} STREQUAL wsm )
  add_definitions( -DALIGNMENT=16
                   -msse3 )
elseif( ${ARCHITECTURE_PLAIN} STREQUAL snb )
  add_definitions( -DALIGNMENT=32
                   -mavx )
elseif( ${ARCHITECTURE_PLAIN} STREQUAL hsw )
  add_definitions( -DALIGNMENT=32
                   -xCORE-AVX2
                   -fma )
elseif( ${ARCHITECTURE_PLAIN} STREQUAL knc )
  add_definitions( -DALIGNMENT=64
                   -mmic
                   -fma )
  set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mmic" )
else()
  message( FATAL_ERROR "GENERATED_KERNELS ${GENERATED_KERNELS} not supported" )
endif()

# add matrix kernel library
add_library( matrix_kernels
             submodules/seissol/Monitoring/FlopCounter.cpp
             ${GENERATED_CODE}/matrix_kernels/sparse_${GENERATED_KERNELS}.cpp
             ${GENERATED_CODE}/matrix_kernels/${PRECISION_PLAIN}gemm_${ARCHITECTURE_PLAIN}.cpp )


# convert to upper case for definition
string( TOUPPER ${GENERATED_KERNELS} GENERATED_KERNELS )

# define the architecture
add_definitions( -D${GENERATED_KERNELS} )
# debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# add submodules
include_directories( src
                     unit_tests
                     submodules
                     submodules/seissol
                   )

# add pugixml library
add_library( pugiXml
             submodules/pugixml/src/pugixml.cpp
           )

# add integration kernels
add_library( time_kernel
             src/Time.cpp )

add_library( volume_kernel
             src/Volume.cpp )

add_library( flux_kernel
             src/Boundary.cpp )

if(NOT SEISSOL_PROXY STREQUAL "OFF" )
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()

  add_definitions( -DNDEBUG )

  add_executable( seissol_proxy ${SEISSOL_PROXY}/proxy_seissol.cpp )

  # enable ipo
  set_property(TARGET matrix_kernels PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
  set_property(TARGET flux_kernel PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
  set_property(TARGET volume_kernel PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
  set_property(TARGET time_kernel PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
  set_property(TARGET seissol_proxy PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

  target_link_libraries( seissol_proxy time_kernel volume_kernel flux_kernel matrix_kernels )
endif()

# add unit tests
find_package( CxxTest )
if(CXXTEST_FOUND)
  set( CXXTEST_USE_PYTHON TRUE )
  include_directories( ${CXXTEST_INCLUDE_DIR} )
  enable_testing()

  # unit tests for XML
  CXXTEST_ADD_TEST( setupUnitTests
                    setup.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/XmlParserTestSuite.t.h
                  )
  target_link_libraries( setupUnitTests pugiXml )

  # unit tests for the matrix kernels
#  CXXTEST_ADD_TEST( matrixKernelsUnitTests
#                    matrixKernels.cpp
#                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/DenseMatrixKernelTestSuite.t.h
#                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/BlockedMatrixKernelsTestSuite.t.h
#                    #${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/MatrixKernelsTestSuite.t.h
#                  )
#  target_link_libraries( matrixKernelsUnitTests pugiXml )

  # unit tests for the time kernel
  CXXTEST_ADD_TEST( time_unit_tests
                    timeUnitTests.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/TimeKernelTestSuite.t.h  
                  )

  # unit tests for the volume kernel
  CXXTEST_ADD_TEST( volume_unit_tests
                    volumeUnitTests.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/VolumeKernelTestSuite.t.h
                  )

  # unit tests for the flux kernel
  CXXTEST_ADD_TEST( flux_unit_tests
                    fluxUnitTests.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/BoundaryKernelTestSuite.t.h
                  )

  target_link_libraries( time_unit_tests   time_kernel   pugiXml matrix_kernels )
  target_link_libraries( volume_unit_tests volume_kernel pugiXml matrix_kernels )  
  target_link_libraries( flux_unit_tests   flux_kernel   pugiXml matrix_kernels )  
endif()
