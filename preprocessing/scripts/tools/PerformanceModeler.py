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
# Copyright (c) 2012
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
# Generates theoretical performance models for SeisSol.
#

import os

import tools.Logger as l_logger
from scipy.io.mmio import *
import numpy as numpy
import csv
from elementtree.ElementTree import parse
import csv
import ntpath

class PerformanceModeler():
  # mapping: #basis functions -> polynomial degree
  m_mapping = { 4:  1,
                10: 2,
                20: 3,
                35: 4,
                56: 5
               }

  # Constructor
  # @param i_pathToOutputDirectory directory, where the performance models will be stored.
  def __init__( self, i_pathToOutputDirectory ):
    # set output directory
    self.pathToOutputDirectory = i_pathToOutputDirectory

  # Deprecated: This function is part of the CK-only performance model and not used anymore.
  #
  # Generates a theoretical performance model for a sparse matrix kernel (spase*dense or dense*sparse).
  #
  # \param i_pathToSparseMatrix path to the sparse matrix (format: matrix market).
  # \param i_numberOfDenseRows #rows of dense matrix.
  # \param i_numberOfDenseColumns #columns of dense matrix
  # \param i_kernelType type of kernel: 'SparseDense' or 'DenseSparse'
  # \return #non zeros of sparse matrx, #scalar multiplications for a single call
  def generateSparseMatrixPerformanceModel( self,
                                            i_pathToSparseMatrix,
                                            i_numberOfDenseRows,
                                            i_numberOfDenseColumns,
                                            i_kernelType ):
    l_logger.log('Generating sparse matrix performance model for '+ i_pathToSparseMatrix, 3)
    
    # assert correct kernel type
    assert( i_kernelType == 'SparseDense' or i_kernelType == 'DenseSparse' )
    
    # read the full matrix
    l_sparseMatrix = mmread(i_pathToSparseMatrix)
    
    # assert correct dimensions
    if( i_kernelType == 'SparseDense' ):
      assert( len(l_sparseMatrix[0]) == i_numberOfDenseRows )
    elif( i_kernelType == 'DenseSparse' ):
      assert( i_numberOfDenseColumns == len(l_sparseMatrix) )
    else:
      assert ( False )
    
    # get number of nonzeros in sparse matrix
    l_numberOfNonZeros = numpy.count_nonzero(l_sparseMatrix)
    
    # get number of scalar multiplications
    if( i_kernelType == 'SparseDense' ):
      l_numberOfScalarMultiplications = l_numberOfNonZeros * i_numberOfDenseColumns
    elif( i_kernelType == 'DenseSparse' ):
      l_numberOfScalarMultiplications = l_numberOfNonZeros * i_numberOfDenseRows
    else:
      assert( False )
    
    # print information  
    l_logger.log('Number of nonzeros: '+str(l_numberOfNonZeros), 4)
    l_logger.log('Number of scalar multiplications: '+str(l_numberOfScalarMultiplications), 4)
    
    return l_numberOfNonZeros, l_numberOfScalarMultiplications
  
  # Compute roofline model.
  # @param io_performanceModel performance model containing loads, multiplications and stores. will be extended by computational intensity.
  def computeRooflineModel( self,
                            io_performanceModel ):
    l_logger.log('computing roofline model', 3)
    
    # compute computational intensity (flops per byte)
    io_performanceModel['computationalIntensity'] = io_performanceModel['numberOfScalarMultiplications'] / \
                                                    ( (io_performanceModel['numberOfMainMemoryLoads'] + io_performanceModel['numberOfMainMemoryStores'])
                                                    * io_performanceModel['precision'] )
    
    # compute peak performance for each CPU (in terms of multiplications)
    io_performanceModel['theoreticalPeakPerformance'] = []
    for l_i in xrange( len(self.architectures['clock_speed_ghz']) ):
      io_performanceModel['theoreticalPeakPerformance'] = io_performanceModel['theoreticalPeakPerformance'] + [ self.architectures['clock_speed_ghz'][l_i] * self.architectures['single_precision_sse3_pipelined_multiply_adds_per_core_and_cycle'][l_i] ]
      
      if (io_performanceModel['precision'] == 8):
        io_performanceModel['theoreticalPeakPerformance'][l_i] = io_performanceModel['theoreticalPeakPerformance'][l_i] / 2.0
      else:
        assert( io_performanceModel['precision'] == 4)
      
      # mininum of memory bandwidth and fpu performance
      io_performanceModel['theoreticalPeakPerformance'][l_i] = min( io_performanceModel['theoreticalPeakPerformance'][l_i],
                                                                    self.architectures['main_memory_bandwidth_gb_per_second'][l_i] / self.architectures['number_of_cores'][l_i] * io_performanceModel['computationalIntensity'] )
 
  # Deprecated: This function is part of the CK-only performance model and not used anymore.
  #
  # Generated a theoretical performance model for the Cauchy Kovalewski kernel.
  #
  # \param i_pathToMatrices location of stiffness and star matrices
  # \param i_numberOfBasisFunctions #basis functions
  # \param i_numberOfVariables #variables
  def generateCauchyKovalewskiPerformanceModel( self,
                                                i_pathToMatrices,
                                                i_numberOfBasisFunctions,
                                                i_numberOfVariables = 9 ):
    l_logger.log('Generating Cauchy Kovalewski performance model for #basis functions='+\
                 str(i_numberOfBasisFunctions)+', #variables='+str(i_numberOfVariables)+'.', 2)
    
    # performace model, default: double precision (8 byte)
    l_performanceModel = { 'numberOfMainMemoryLoads': 0,
                           'numberOfMainMemoryStores': 0,
                           'numberOfScalarMultiplications': 0,
                           'precision': 8.0
                         }
    
    # mapping polynomial degree of <-> #basis functions
    l_mapping = {'polynomialDegree': [], 'numberOfBasisFunctions': []}
    l_mapping['polynomialDegree'] = [1, 2, 3, 4, 5]
    l_mapping['numberOfBasisFunctions'] = [4, 10, 20, 35, 56]
    
    # assert a valid number of basis functions
    assert( i_numberOfBasisFunctions in l_mapping['numberOfBasisFunctions'] )
    
    # get polynomial degree
    l_mappingIndex = l_mapping['numberOfBasisFunctions'].index(i_numberOfBasisFunctions)
    l_polynomialDegree = l_mapping['polynomialDegree'][l_mappingIndex]
    
    # matrix properties
    l_matrices = { 'baseName': [],\
                   'fileName': [],\
                   'kernelType': [],\
                   'numberOfNonZeros': [],\
                   'numberOfScalarMultiplications': []\
                 }
    
    # add star matrix, representative for (A*, B*, C*)
    l_matrices['baseName'] = ['starMatrix']
    l_matrices['fileName'] = l_matrices['fileName'] + ['starMatrix_3D_maple.mtx']
    l_matrices['kernelType'] = l_matrices['kernelType'] + ['DenseSparse']
    l_matrices['numberOfNonZeros'] = l_matrices['numberOfNonZeros'] + [-1]
    l_matrices['numberOfScalarMultiplications'] = l_matrices['numberOfScalarMultiplications'] + [-1]
    
    for l_baseName in 'kXi', 'kEta', 'kZeta':
      l_matrices['baseName'] = l_matrices['baseName'] + [l_baseName]
      l_matrices['fileName'] = l_matrices['fileName'] + [ l_baseName+\
                                                          'DivMT_3D_' +\
                                                          str(l_polynomialDegree) +\
                                                          '_maple.mtx']
      l_matrices['kernelType'] = l_matrices['kernelType'] + ['SparseDense']
      # initilize currently unknown values
      l_matrices['numberOfNonZeros'] = l_matrices['numberOfNonZeros'] + [-1]
      l_matrices['numberOfScalarMultiplications'] = l_matrices['numberOfScalarMultiplications'] + [-1]
    
    # get files in matrix folder
    l_matrixFiles = os.listdir('matrices')
    l_matrixFiles.sort()
    
    # iterate of over files
    for l_file in l_matrixFiles:
      if( l_file in l_matrices['fileName'] ):
        # get index in matrix list
        l_matrixIndex = l_matrices['fileName'].index(l_file)
        
        # get matrix local performance data
        l_matrices['numberOfNonZeros'][l_matrixIndex], \
        l_matrices['numberOfScalarMultiplications'][l_matrixIndex]   = \
          self.generateSparseMatrixPerformanceModel( i_pathToSparseMatrix = 'matrices/'+l_file,
                                                     i_numberOfDenseRows = i_numberOfBasisFunctions,
                                                     i_numberOfDenseColumns = i_numberOfVariables,
                                                     i_kernelType = l_matrices['kernelType'][l_matrixIndex])
    #                                     
    # generate performance model
    #
    
    l_numberOfUnknowns = i_numberOfBasisFunctions * i_numberOfVariables
    
    # load unknowns and time integrated unknowns
    l_performanceModel['numberOfMainMemoryLoads'] = 2*( l_performanceModel['numberOfMainMemoryLoads'] + \
                                                        l_numberOfUnknowns )
    
    # load three star matrices into cache (A*, B*, C*)
    l_matrixIndex = l_matrices['baseName'].index('starMatrix')
    l_performanceModel['numberOfMainMemoryLoads'] = l_performanceModel['numberOfMainMemoryLoads'] + \
                                                    3*l_matrices['numberOfNonZeros'][l_matrixIndex]
    
    # compute first taylor step
    l_performanceModel['numberOfScalarMultiplications'] = l_performanceModel['numberOfScalarMultiplications'] +\
                                                          l_numberOfUnknowns
    
    # iterate over talyor expansion
    for l_taylorSeriesOrder in range(1, l_polynomialDegree+1):
      # compute stiffness matrix multiplications
      for l_baseName in 'kXi', 'kEta', 'kZeta': 
        l_matrixIndex = l_matrices['baseName'].index(l_baseName)
        l_performanceModel['numberOfScalarMultiplications'] = l_performanceModel['numberOfScalarMultiplications'] +\
                                                              l_matrices['numberOfScalarMultiplications'][l_matrixIndex]
      l_matrixIndex = l_matrices['baseName'].index('starMatrix')
      # compute star matrix multiplications
      l_performanceModel['numberOfScalarMultiplications'] = l_performanceModel['numberOfScalarMultiplications'] +\
                                                            3*(l_matrices['numberOfScalarMultiplications'][l_matrixIndex])
      
      # update time integrated unknowns
      l_performanceModel['numberOfScalarMultiplications'] = l_performanceModel['numberOfScalarMultiplications'] +\
                                                            l_numberOfUnknowns
    # write time integrated unknowns to main memory
    l_performanceModel['numberOfMainMemoryStores'] = l_performanceModel['numberOfMainMemoryStores'] + \
                                                     l_numberOfUnknowns
    
    # derive roofline model
    self.computeRooflineModel( io_performanceModel = l_performanceModel )
    
    l_logger.log(str(l_performanceModel), 3)
    
    # compute time per iteration
    l_performanceModel['timePerCauchyKovalewskiIteration'] = []
    l_logger.log('Minimal single core time per Cauchy Kovalewski iteration (roofline model for main memory bandwidth)', 3)
    for l_i in xrange( len(l_performanceModel['theoreticalPeakPerformance']) ):
      l_performanceModel['timePerCauchyKovalewskiIteration'] = l_performanceModel['timePerCauchyKovalewskiIteration'] +\
                                                               [ l_performanceModel['numberOfScalarMultiplications'] / l_performanceModel['theoreticalPeakPerformance'][l_i] ]
      l_performanceModel['timePerCauchyKovalewskiIteration'][l_i] = l_performanceModel['timePerCauchyKovalewskiIteration'][l_i] / (1E+9)
      
      l_logger.log( self.architectures['cluster'][l_i]+ ' - ' +
                    self.architectures['full_name'][l_i]+
                    ' (' + self.architectures['codename'][l_i] +') '+
                    '@' + str(self.architectures['clock_speed_ghz'][l_i]) + 'GHz, '+ str(self.architectures['main_memory_bandwidth_gb_per_second'][l_i]) + 'GB/s: '+ 
                    str(l_performanceModel['timePerCauchyKovalewskiIteration'][l_i]), 4)
    
    # write results to disk
    with open(self.pathToCauchyKovalewskiPerformanceModel, 'a') as l_csvFile:
      for l_i in xrange( len(l_performanceModel['theoreticalPeakPerformance']) ):
        l_csvFile.write('ck_time_integration' + ',' )
        
        # write number of basis functions and variables
        l_csvFile.write( str(i_numberOfBasisFunctions) + ',' )
        l_csvFile.write( str(i_numberOfVariables) + ',' )
      
        # write architecture information
        l_csvFile.write( str(self.architectures['cluster'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['full_name'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['codename'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['number_of_cores'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['clock_speed_ghz'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['single_precision_sse3_pipelined_multiply_adds_per_core_and_cycle'][l_i]) + ',' )
        l_csvFile.write( str(self.architectures['main_memory_bandwidth_gb_per_second'][l_i]) + ',' )
        
        # write derived data
        l_csvFile.write( str(l_performanceModel['numberOfMainMemoryLoads']) + ',' )
        l_csvFile.write( str(l_performanceModel['numberOfMainMemoryStores']) + ',' )
        l_csvFile.write( str(l_performanceModel['precision']) + ',' )
        l_csvFile.write( str(l_performanceModel['numberOfScalarMultiplications']) + ',' )
        l_csvFile.write( str(l_performanceModel['computationalIntensity']) + ',' )
        l_csvFile.write( str(l_performanceModel['theoreticalPeakPerformance'][l_i]) + ',' )
        l_csvFile.write( str(l_performanceModel['timePerCauchyKovalewskiIteration'][l_i]) + '\n' )

  # Read matrix information from a XML-file.
  #
  # i_pathToMatricesFile path to matrices xml-file.
  def readMatricesFile( self, i_pathToMatricesFile ):
    l_logger.log( 'reading xml-file ' + i_pathToMatricesFile, 2 )
    # parse XML-file
    l_xmlRoot = parse( i_pathToMatricesFile ).getroot()

    l_matrices = {}

    # add global matrices
    for l_globalMatrix in  l_xmlRoot.find( 'global' ):
      # check for valid numbers
      assert( int(l_globalMatrix.attrib['rows']) == int(l_globalMatrix.attrib['columns']) > 0 )
      assert( len( list(l_globalMatrix) ) <=  int(l_globalMatrix.attrib['rows']) *  int(l_globalMatrix.attrib['columns']) )

      l_matrices[ l_globalMatrix.attrib['name'] ] = {}
      l_matrices[ l_globalMatrix.attrib['name'] ]['#non_zeros'] = len( list(l_globalMatrix) )
      l_matrices[ l_globalMatrix.attrib['name'] ]['#rows'] = int(l_globalMatrix.attrib['rows'])
      l_matrices[ l_globalMatrix.attrib['name'] ]['#columns'] = int(l_globalMatrix.attrib['columns'])
      l_matrices[ l_globalMatrix.attrib['name'] ]['sparse'] = l_globalMatrix.attrib['sparse']
      l_matrices[ l_globalMatrix.attrib['name'] ]['type'] = l_globalMatrix.tag

      # build recursive structure of the ADER time integration
      if l_globalMatrix.attrib['name'] in ['kXiDivMT', 'kEtaDivMT', 'kZetaDivMT']:
        # initialize #nnz
        l_matrices[ l_globalMatrix.attrib['name'] ]['#non_zeros'] = {}

        # order of the taylor series expansion for this #(basis functions)
        l_orderOfTaylorSeriesExpansion = self.m_mapping[ l_matrices[ l_globalMatrix.attrib['name'] ]['#rows'] ] + 1

        # iterate over all sub-matrices of the orders
        for l_order in range(1, l_orderOfTaylorSeriesExpansion ):
          # compute size of the non-zero stiffness-block/submatrices
          l_nonZeroBlockSizeStiffness = l_orderOfTaylorSeriesExpansion - l_order + 1
          l_nonZeroBlockSizeStiffness = l_nonZeroBlockSizeStiffness * (l_nonZeroBlockSizeStiffness + 1) * (l_nonZeroBlockSizeStiffness + 2) / 6

          # set #nnz to zero
          l_matrices[ l_globalMatrix.attrib['name'] ]['#non_zeros'][l_nonZeroBlockSizeStiffness] = 0

          # count the #nnz
          for l_element in l_globalMatrix:
            # check if the nnz is part of non-zero block, if yes: add
            if int(l_element.attrib['row']) <= l_nonZeroBlockSizeStiffness and int(l_element.attrib['column']) <= l_nonZeroBlockSizeStiffness:
              l_matrices[ l_globalMatrix.attrib['name'] ]['#non_zeros'][l_nonZeroBlockSizeStiffness] = l_matrices[ l_globalMatrix.attrib['name'] ]['#non_zeros'][l_nonZeroBlockSizeStiffness] + 1

    # add local matrices
    # TODO: Hardcoded check if more general information is necessary
    l_matrices[ 'flux_solver' ] = { '#non_zeros': 9*9,
                                    '#rows':     9,
                                    '#columns':  9,
                                    'sparse':    'false',
                                    'type':      'jacobian'
                                 }
    l_matrices[ 'star_matrix' ] = { '#non_zeros': 24,
                                    '#rows':     9,
                                    '#columns':  9,
                                    'sparse':    'true',
                                    'type':      'jacobian'
                                  }

    return l_matrices  

  # Generates a performance model for the time, volume and boundary kernel.
  #
  # i_pathToMatricesFile path to matrices xml-file.
  # i_pathToOutputDirectory path to the output directory.
  def generatePerformanceModel( self,
                                i_pathToMatricesFile,
                                i_pathToOutputDirectory='performance_models' ):
    l_logger.log( 'generating volume performance model for ' + i_pathToMatricesFile )
    l_matrices = self.readMatricesFile( i_pathToMatricesFile = i_pathToMatricesFile )

    # get number of basis functions from the stiffness matrix and variables from the star matrix
    l_numberOfBasisFunctions = int(l_matrices['kXiDivM']['#rows'])
    l_numberOfVariables = int(l_matrices['star_matrix']['#rows'])

    # determine corr. polynomial degree
    l_polynomialDegree = self.m_mapping[ l_numberOfBasisFunctions ]

    l_logger.log( '#basis functions: ' + str(l_numberOfBasisFunctions) +\
                  ', #variables: ' + str(l_numberOfVariables)+\
                  ', polynomial degree: '+ str(l_polynomialDegree), 2 )

    # initialize performance model
    l_performanceModel = {}

    # generate names for flux matrices
    l_fluxMatrices = ['fM1', 'fM2', 'fM3', 'fM4']
    for l_localSide in ['1','2','3','4']:
      for l_neighboringSide in ['1', '2', '3', '4']:
        for l_vertexCombination in ['1', '2', '3']:
          l_fluxMatrices = l_fluxMatrices+['fP'+l_localSide+l_neighboringSide+l_vertexCombination]

    # collect #scalar floating point multiplications for volume and boundary matrices
    for l_matrixName in ['kXiDivM',  'kEtaDivM',  'kZetaDivM',
                         'kXiDivMT', 'kEtaDivMT', 'kZetaDivMT',
                         'flux_solver', 'star_matrix'] + l_fluxMatrices:
      # TODO: take care of this guy for dense matrices in ADER time integration
      if( l_matrices[l_matrixName]['sparse'] == 'true' ):
        l_numberOfNonZeros = l_matrices[l_matrixName]['#non_zeros']
      else:
        assert(  l_matrices[l_matrixName]['sparse'] == 'false' )
        l_numberOfNonZeros = int(l_matrices[l_matrixName]['#rows']) * int(l_matrices[l_matrixName]['#columns'])

      # determine type of multiplication (left * unknowns, unknowns * right)
      if( l_matrices[l_matrixName]['type'] in ['flux', 'stiffness'] ):
        l_repeatsPerEntry = l_numberOfVariables
      else:
        assert( l_matrices[l_matrixName]['type'] == 'jacobian' )
        l_repeatsPerEntry = l_numberOfBasisFunctions

      # default case: volume or boundary integration
      if( l_matrixName not in [ 'kXiDivMT', 'kEtaDivMT', 'kZetaDivMT' ] ):
        # create new performance model
        l_performanceModel[l_matrixName] = {}

        # each element is multiplied by a complete row of B (stiffness, flux) / column of A (jacobian)
        l_performanceModel[l_matrixName]['FLOPS'] = ( l_numberOfNonZeros * l_repeatsPerEntry ) * 2      
      # time integration, each appearing submatrix has to be taken care of
      else:
        # iterate over submatrices
        for l_key in l_numberOfNonZeros:
          # create new performance model
          l_performanceModel[l_matrixName+'_'+str(l_key)] = {}

          # add #(floating point operations)
          l_performanceModel[l_matrixName+'_'+str(l_key)]['FLOPS'] = ( l_numberOfNonZeros[l_key] * l_repeatsPerEntry ) * 2

    # order of the taylor series expansion
    l_orderOfTaylorSeriesExpansion = l_polynomialDegree+1

    # set up ader star matrices
    for l_order in xrange(1, l_orderOfTaylorSeriesExpansion):
      # compute size of the non-zero star-dense-blocks
      l_nonZeroBlockSizeStar = l_orderOfTaylorSeriesExpansion - l_order
      l_nonZeroBlockSizeStar = l_nonZeroBlockSizeStar * (l_nonZeroBlockSizeStar + 1) * (l_nonZeroBlockSizeStar + 2) / 6

      # create new performance model
      l_performanceModel['star_matrix_ader_'+str(l_nonZeroBlockSizeStar)] = {}

      # set #(floating point operations)
      l_performanceModel['star_matrix_ader_'+str(l_nonZeroBlockSizeStar)]['FLOPS'] = ( l_matrices['star_matrix']['#non_zeros'] * l_nonZeroBlockSizeStar ) * 2

    #initialize time, volume and element boundary kernel
    l_performanceModel['time_kernel'] = {}
    l_performanceModel['volume_kernel']= {}
    l_performanceModel['boundary_kernel_element_local'] = {}

    #
    # Cauchy Kovaleski kernel
    #    
    l_numberOfUnknowns = l_numberOfBasisFunctions * l_numberOfVariables
    
    # compute first taylor step
    l_performanceModel['time_kernel']['FLOPS'] = l_numberOfUnknowns

    # iterate over talyor expansion
    for l_order in range(1, l_orderOfTaylorSeriesExpansion):
      # compute size of the non-zero blocks/submatrices
      l_nonZeroBlockSizeStiffness = l_orderOfTaylorSeriesExpansion - l_order + 1
      l_nonZeroBlockSizeStiffness = l_nonZeroBlockSizeStiffness * (l_nonZeroBlockSizeStiffness + 1) * (l_nonZeroBlockSizeStiffness + 2) / 6

      l_nonZeroBlockSizeStar = l_orderOfTaylorSeriesExpansion - l_order
      l_nonZeroBlockSizeStar = l_nonZeroBlockSizeStar * (l_nonZeroBlockSizeStar + 1) * (l_nonZeroBlockSizeStar + 2) / 6

      # compute stiffness matrix multiplications
      l_performanceModel['time_kernel']['FLOPS'] = l_performanceModel['time_kernel']['FLOPS'] +\
                                                   l_performanceModel['kXiDivMT_'   + str(l_nonZeroBlockSizeStiffness)]['FLOPS'] +\
                                                   l_performanceModel['kEtaDivMT_'  + str(l_nonZeroBlockSizeStiffness)]['FLOPS'] +\
                                                   l_performanceModel['kZetaDivMT_' + str(l_nonZeroBlockSizeStiffness)]['FLOPS']
      # compute star matrix multiplications
      l_performanceModel['time_kernel']['FLOPS'] =   l_performanceModel['time_kernel']['FLOPS'] +\
                                                     l_performanceModel['star_matrix_ader_'+str(l_nonZeroBlockSizeStar)]['FLOPS'] * 3
      
      # update time integrated unknowns
      l_performanceModel['time_kernel']['FLOPS'] =   l_performanceModel['time_kernel']['FLOPS'] +\
                                                     ( l_numberOfUnknowns ) * 2

    #
    # volume kernel
    #
    l_performanceModel['volume_kernel']['FLOPS'] =   l_performanceModel['kXiDivM']['FLOPS']            +\
                                                     l_performanceModel['kEtaDivM']['FLOPS']           +\
                                                     l_performanceModel['kZetaDivM']['FLOPS']          +\
                                                     l_performanceModel['star_matrix']['FLOPS'] * 3
    #
    # element local boundary integral
    # Remark: We can't determin the part of the boundary kernel, which accounts for the contribution of the
    #         neighboring element, because the choice of flux matrices at the faces is mesh dependent.
    l_performanceModel['boundary_kernel_element_local']['FLOPS'] = l_performanceModel['fM1']['FLOPS'] +\
                                                                   l_performanceModel['fM2']['FLOPS'] +\
                                                                   l_performanceModel['fM3']['FLOPS'] +\
                                                                   l_performanceModel['fM4']['FLOPS'] +\
                                                                   l_performanceModel['flux_solver']['FLOPS'] * 4

    # get base name from input file
    l_baseName = ntpath.basename( i_pathToMatricesFile )

    # open csv writer
    l_writer = csv.writer(open(i_pathToOutputDirectory+'/'+l_baseName.replace('xml', 'csv'), 'w'))

    # add header
    l_writer.writerow( ['kernel', 'FLOPS'] )

    # write performance model
    for l_key in sorted(l_performanceModel):
      l_writer.writerow( [l_key, l_performanceModel[l_key]['FLOPS'] ] )

