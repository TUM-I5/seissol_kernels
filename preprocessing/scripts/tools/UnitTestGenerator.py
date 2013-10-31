#!/bin/env python
# @file
# This file is part of SeisSol.
# 
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# 
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
#
# Generates unit tests for matrix kernels.
#

import tools.Logger as l_logger
import tools.SeisSolGen as l_seisSolGen

class UnitTestGenerator():
  ###
  ### Class Variables
  ###
  m_pathToMatrices = ''

  # Constructor
  #
  # @param i_pathToMatrices path, where the matrices are stored.
  def __init__( self, i_pathToMatrices ):
    self.m_pathToMatrices = i_pathToMatrices
  
  # Write matrix unit tests source code
  def writeMatrixUnitTests( self,
                            i_matrices,
                            i_pathToOutputFile,
                            i_kernelType ):
    # assert a supported matrix kernel type
    assert( i_kernelType == 'sparse_volume_boundary' or i_kernelType == 'sparse_time' or i_kernelType == 'dense')

    # holds the generated source code
    l_sourceCode = ''

    # remove old file contents
    open(i_pathToOutputFile, 'w').close()

    # write license
    l_logger.writeFileHeader(i_pathToOutputFile, '// ')

    l_firstIteration = True

    for l_matrix in i_matrices:
      # generate if-clause
      if( l_firstIteration == True ):
        l_firstIteration = False
        l_sourceCode = l_sourceCode + 'if'
      else:
        l_sourceCode = l_sourceCode + 'else if'

      if( i_kernelType == 'sparse_volume_boundary' ):
        l_sourceCode = l_sourceCode + '( numberOfRows == ' + str(l_matrix['numberOfDenseRows']) + ' && '\
            'std::string::npos != pathToMatrixMarketFile.find("'+ str(l_matrix['pathToMatrixMarketFile']) +'") ){ \n'
      elif( i_kernelType == 'sparse_time' ):
        l_sourceCode = l_sourceCode + '(  l_numberOfBasisFunctions == ' + str(l_matrix['numberOfDenseRows']) + ' && '\
            'l_matrixName=="'+ str(l_matrix['name']) +'" ){ \n'
      elif( i_kernelType == 'dense' ):
        l_sourceCode = l_sourceCode + '( i_m == ' + str(l_matrix['m']) + ' && '+\
                                        'i_n == ' + str(l_matrix['n']) + ' && '+\
                                        'i_k == ' + str(l_matrix['k']) + ' ){ \n'

      # generate matrix call
      if( i_kernelType == 'sparse_volume_boundary' ):
        l_sourceCode = l_sourceCode + '  '+ l_matrix['routineNameOfGeneratedKernel'] + '( flatColumnMajor,l_testMatrix, l_generatedResultMatrix ); \n'
      elif( i_kernelType == 'sparse_time' ):
        l_sourceCode = l_sourceCode + '  '+ l_matrix['routineNameOfGeneratedKernel'] + '( l_a, l_b, l_c1, l_nonZeroBlockSize ); \n'
      elif( i_kernelType == 'dense' ):
        l_sourceCode = l_sourceCode + '  '+ l_matrix['routineNameOfGeneratedKernel'] + '( i_a, i_b, io_c ); \n'
      # close conditional statement
      l_sourceCode = l_sourceCode + '} \n'

    # assert that we found our matrix call
    l_sourceCode = l_sourceCode + 'else \n'
    
    if( i_kernelType == 'sparse_volume_boundary' ):
      l_sourceCode = l_sourceCode + '  TS_FAIL("invalid number of basis function or matrix file name"); \n'
    elif( i_kernelType == 'sparse_time' ):
      l_sourceCode = l_sourceCode + '  TS_FAIL("invalid number of basis functions or matrix name"); \n'
    elif( i_kernelType == 'dense' ):
      l_sourceCode = l_sourceCode + '   TS_FAIL("invalid matrix dimensions (m, n, k)"); \n'

    # write code to file
    l_file = open(i_pathToOutputFile ,'a')
    l_file.write(l_sourceCode)
    l_file.close()

  # Generates Unit tests for matrix kernels.
  #
  # @param i_pathToOutputDirectory directory to which the code will be written
  def generateMatrixUnitTests( self, i_pathToOutputDirectory ):
    l_logger.log('generating unit tests for sparse matrix kernels', 2)

    # where to save the sparse unit tests for volume and boundary integration
    l_pathToOutputFile = i_pathToOutputDirectory+'/volume_boundary_sparse_matrix_kernels.hpp_include'

    # get sparse matrices
    l_sparseMatrices = l_seisSolGen.getSparseMatrices( self.m_pathToMatrices )

    # filter matrices of the volume and boundary integration which don't have special handling of zero block
    l_sparseVolumeBoundary = [l_sparseMatrix for l_sparseMatrix in l_sparseMatrices if l_sparseMatrix['aderZeroBlocks'] == False]

    # write unit tests for sparse matrices, without zero block appearing in the ADER time integration
    self.writeMatrixUnitTests( i_matrices = l_sparseVolumeBoundary,
                               i_pathToOutputFile = l_pathToOutputFile,
                               i_kernelType = 'sparse_volume_boundary' )


    # where to save the sparse unit tests for time integration
    l_pathToOutputFile = i_pathToOutputDirectory+'/time_sparse_matrix_kernels.hpp_include'

    # filter matrices of the time integration which have special handling of zero block
    l_sparseTime = [l_sparseMatrix for l_sparseMatrix in l_sparseMatrices if l_sparseMatrix['aderZeroBlocks'] == True]

    # write unit tests for sparse matrices, without zero block appearing in the ADER time integration
    self.writeMatrixUnitTests( i_matrices = l_sparseTime,
                               i_pathToOutputFile = l_pathToOutputFile,
                               i_kernelType = 'sparse_time' )


    # where to save the dense unit tests
    l_pathToOutputFile = i_pathToOutputDirectory+'/dense_matrix_kernels.hpp_include'

    # get dense matrices
    l_denseMatrices = l_seisSolGen.getDenseMatrices()

    # write unit tests for dense matrices
    self.writeMatrixUnitTests( i_matrices = l_denseMatrices,
                               i_pathToOutputFile = l_pathToOutputFile,
                               i_kernelType = 'dense' )
