##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
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

# set number of basis functions
add_definitions( -DNUMBEROFBASISFUNCTIONS=20
                 -DNUMBEROFVARIABLES=9
                 -DGENERATEDKERNELS )

# use intel compilers
include(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")
CMAKE_FORCE_Fortran_COMPILER(ifort "Intel Fortran Compiler")

# enable C++11
ADD_DEFINITIONS( -std=c++11 )

# add submodules
include_directories( submodules
                     submodules/seissol
                   )

# add pugixml library
add_library( pugiXml
             submodules/pugixml/src/pugixml.cpp
           )

# add memory manager
add_library( memoryManager
             submodules/seissol/Initializer/MemoryAllocator.cpp
             submodules/seissol/Initializer/MemoryManager.cpp )

# add integrators
add_library( integrators
             submodules/seissol/Monitoring/FlopCounter.cpp
             src/TimeIntegrator.cpp
             src/VolumeIntegrator.cpp
             src/BoundaryIntegrator.cpp )

# add unit tests
find_package( CxxTest )
if(CXXTEST_FOUND)
  set( CXXTEST_USE_PYTHON TRUE )
  include_directories( ${CXXTEST_INCLUDE_DIR} )
  enable_testing()

  # unit tests for XML and memory
  CXXTEST_ADD_TEST( setupUnitTests
                    setup.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/XmlParserTestSuite.t.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/MemoryAllocatorTestSuite.t.h
                  )
  target_link_libraries( setupUnitTests memoryManager pugiXml )

  # unit tests for the matrix kernels
  CXXTEST_ADD_TEST( matrixKernelsUnitTests
                    matrixKernels.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/DenseMatrixKernelTestSuite.t.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/BlockedMatrixKernelsTestSuite.t.h
                    #${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/MatrixKernelsTestSuite.t.h
                  )
  target_link_libraries( matrixKernelsUnitTests memoryManager pugiXml )

  # unit tests for the integrators
  CXXTEST_ADD_TEST( integratorsUnitTests
                    integrators.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/TimeIntegratorTestSuite.t.h 
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/VolumeIntegratorTestSuite.t.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/BoundaryIntegratorTestSuite.t.h     
                  )
  target_link_libraries( integratorsUnitTests integrators memoryManager pugiXml )  
endif()
