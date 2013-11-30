# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
#
# According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
#
# The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
#
# Copyright (c) 2013
# Technische Universitaet Muenchen
# Department of Informatics
# Chair of Scientific Computing
# http://www5.in.tum.de/
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
# Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Builds the seissol_kernels unit tests.

cmake_minimum_required( VERSION 2.6 )

# set number of basis functions
add_definitions( -DNUMBEROFBASISFUNCTIONS=10
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
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/VolumeIntegratorTestSuite.t.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/BoundaryIntegratorTestSuite.t.h     
                  )
  target_link_libraries( integratorsUnitTests integrators memoryManager pugiXml )  
endif()
