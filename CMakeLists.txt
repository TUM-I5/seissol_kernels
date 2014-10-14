##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013-2014, SeisSol Group
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

cmake_minimum_required( VERSION 2.6 )

option( VECTOR_INSTRUCTION_SET "Vectorization instruction set: avx, avx2 or mic" )
option( CONVERGENCE_ORDER      "Convergence order of the ADER-DG scheme" )

add_definitions( -DCONVERGENCE_ORDER=${CONVERGENCE_ORDER} )

# set number of quantities
add_definitions( -DNUMBER_OF_QUANTITIES=9 )

# build the architecture dependent matrix kernels as library
if( ${VECTOR_INSTRUCTION_SET} STREQUAL avx )
  message( "-- using 32-bit aligned matrix kernels" )

  add_library( matrix_kernels
               submodules/seissol/Monitoring/FlopCounter.cpp
               preprocessing/generated_code/matrix_kernels/dgemm_32.cpp )
  add_definitions( -DALIGNMENT=32
                   -m${VECTOR_INSTRUCTION_SET} )
elseif( ${VECTOR_INSTRUCTION_SET} STREQUAL avx2 OR ${VECTOR_INSTRUCTION_SET} STREQUAL mic )
  message( "-- using 64-bit aligned matrix kernels" )

  add_library( matrix_kernels
               submodules/seissol/Monitoring/FlopCounter.cpp
               preprocessing/generated_code/matrix_kernels/dgemm_64.cpp )
#TODO: add proper avx2 command
  add_definitions( -DALIGNMENT=64
                   -m${VECTOR_INSTRUCTION_SET} )
else()
  message( WARNING "vector instruction set ${VECTOR_INSTRUCTION_SET} not supported" )
  message( "-- using 32-bit alignment and compiler generated kernels" )
  add_definitions( -DALIGNMENT=32 )
  add_library( matrix_kernels
               submodules/seissol/Monitoring/FlopCounter.cpp
               preprocessing/generated_code/matrix_kernels/dgemm_32.cpp )
endif()

# enable C++11
ADD_DEFINITIONS( -std=c++11 )

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

# add kernels
add_library( time_kernel
             src/Time.cpp )

add_library( volume_kernel
             src/Volume.cpp )

add_library( flux_kernel
             src/Boundary.cpp )

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
