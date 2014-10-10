#!/bin/env python
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
#
# Assembles the matrix operations in SeisSol.
#

###
### Libraries
###
import argparse
import os

import re

import tools.Logger as l_logger

l_availableModules = dict()

try:
  import tools.MatrixConverter as MatrixConverter
  import numpy
  l_availableModules['matrix_converter'] = True
except ImportError:
  l_availableModules['matrix_converter'] = False

try:
  import tools.Vectorizer as l_vectorizer
  l_availableModules['vectorizer'] = True
except ImportError:
  l_availableModules['vectorizer'] = False

try:
  import tools.SeisSolGen as l_seisSolGen
  l_availableModules['seissol_gen'] = True
except ImportError:
  l_availableModules['seissol_gen'] = False

try:
  import tools.PerformanceModeler as PerformanceModeler
  l_availableModules['performance_modeler'] = True
except ImportError:
  l_availableModules['performance_modeler'] = False

try:
  import tools.UnitTestGenerator as UnitTestGenerator
  l_availableModules['unit_test_generator'] = True
except ImportError:
  l_availableModules['unit_test_generator'] = False

###
### Command line arguments
###
l_commandLineParser = argparse.ArgumentParser( description=
  'This script is the frontend for generated matrix kernels in SeisSol; \
   Available modules in your python configuration:\n' + str(l_availableModules) +'.'
)

l_commandLineParser.add_argument('--convertFullToSparse',
                                 action='store_true',
                                 help='convert full matrices (generated by Maple) to sparse "CSC-/CSR-Matrix Market" matrices.')

l_commandLineParser.add_argument('--convertToXml',
                                 action='store_true',
                                 help='convert full matrices (generated by Maple) to XML with sparse index storage.')

l_commandLineParser.add_argument('--plotSparsityPatterns',
                                 action='store_true',
                                 help='plots sparsity patterns of the matrices.')

l_commandLineParser.add_argument('--runBenchmarks',
                                 action='store_true',
                                 help='run units tests.')

l_commandLineParser.add_argument('--runVectorization',
                                 action='store_true',
                                 help='vectorizes the columns of the matrices.')

l_commandLineParser.add_argument('--generatePreprocessorCode',
                                 action='store_true',
                                 help='generates preprocessor code which contains the dimensions of the matrices.')

l_commandLineParser.add_argument('--generateStarMatrixInitializationCode',
                                 action='store_true',
                                 help='generates Fortran Code, which initializes a flat star matrix given the full matrix.')

l_commandLineParser.add_argument('--generateMatrixKernelsInitializationCode',
                                 action='store_true',
                                 help='generates C Code, which initializes the function pointers to the kernels of the matrices appearing in the boundary integration.')
                                 
l_commandLineParser.add_argument('--generateMatrixKernels',
                                 action='store_true',
                                 help='generates the C++ matrix kernels.')

l_commandLineParser.add_argument('--generateMatrixUnitTests',
                                 action='store_true',
                                 help='generates unit tests for the C++ matrix kernels.')

l_commandLineParser.add_argument('--generatePerformanceModel',
                                 action='store_true',
                                 help='generates a theoretical performance model.')

l_commandLineArguments = l_commandLineParser.parse_args()

###
### Main
###

l_logger.printWelcomeMessage()

# check for valid settings
if l_commandLineArguments.runBenchmarks:
  l_logger.log('Benchmark implementation outdated. Aborting')
  exit()

# construct classes
if l_availableModules['matrix_converter']:
  l_matrixConverter = MatrixConverter.MatrixConverter()

if l_availableModules['performance_modeler']:
  l_performanceModeler = PerformanceModeler.PerformanceModeler( i_pathToOutputDirectory = 'performance_models' )

if l_availableModules['unit_test_generator']:
  l_unitTestGenerator = UnitTestGenerator.UnitTestGenerator( i_pathToMatrices = 'matrices' )

#if l_commandLineArguments.convertFullToSparse or l_commandLineArguments.plotSparsityPatterns or l_commandLineArguments.runVectorization or l_commandLineArguments.generatePreprocessorCode or #l_commandLineArguments.generateStarMatrixInitializationCode:
# get the matrix files and sort them
l_matrixFiles = os.listdir('matrices')
l_matrixFiles.sort()

for l_file in l_matrixFiles:
  if l_file.endswith('_maple.mtx'):
    l_baseName = l_file.replace('_maple.mtx','')
    if l_commandLineArguments.convertFullToSparse:
      l_matrixConverter.convertFullToSparse( 'matrices/'+l_file, l_baseName, 'matrices' )
    if l_commandLineArguments.plotSparsityPatterns:
      l_matrixConverter.plotSparsityPattern( i_fullMatrix = 'matrices/'+l_file,
                                             i_baseName = l_baseName,
                                             i_pathToOutputDirectory = 'matrices',
                                             i_readMatrix = True )
    if l_commandLineArguments.runVectorization:
      l_vectorizer.vectorizeMatrix( 'matrices/'+l_file, l_baseName, 'vectorization', True)
    if l_commandLineArguments.generatePreprocessorCode:
      l_matrixConverter.addMatrixToPreProcessorCode( i_pathToFullMatrix = 'matrices/'+l_file,
                                                     i_baseName = l_baseName)
    if l_baseName == 'starMatrix_3D':
      if l_commandLineArguments.generateStarMatrixInitializationCode:
        l_matrixConverter.generateMatrixInitializationCode('matrices/'+l_file, l_baseName, 'initializeFlatStarMatrixColumnMajor', 'csc' ,'generated_code/initialization')

if l_commandLineArguments.convertToXml:
  l_sparseDenseSwitches = [1.0] #numpy.arange(0, 1.01, 0.01).tolist()

  for l_sparseDenseSwitch in l_sparseDenseSwitches:
    l_matrixConverter.convertToXml( i_pathToMatrices='matrices',
                                    i_pathToOutputDirectory='matrices',
                                    i_sparseDenseSwitch=l_sparseDenseSwitch  )

if l_commandLineArguments.generatePerformanceModel:
  for l_matrixFile in ['matrices_4.xml', 'matrices_10.xml', 'matrices_20.xml', 'matrices_35.xml', 'matrices_56.xml']:
    l_pathToMatricesFile = 'matrices/'+l_matrixFile
    l_performanceModeler.generatePerformanceModel( i_pathToMatricesFile = l_pathToMatricesFile )
  #l_matrixFiles = [ 'matrices_0.18_20.xml', 'matrices_0.25_35.xml', 'matrices_0.15_56.xml' ]

  # pure sparse (execpt for 'flux solver') performance model
  #l_matrixFiles = l_matrixFiles +\
  #                  ['matrices_1.01_4.xml', 'matrices_1.01_10.xml', 'matrices_1.01_20.xml',
  #                   'matrices_1.01_35.xml', 'matrices_1.01_56.xml']

  #for l_matrix in l_matrixFiles:
  #l_performanceModeler.generatePerformanceModel( i_pathToMatricesFile = 'matrices_56.xml' )
  
  if l_commandLineArguments.generatePreprocessorCode:
    l_matrixConverter.writePreProcessorCode('generated_code/pre_processor')

if l_commandLineArguments.runBenchmarks:
  l_logger.log('running benchmarks')
  l_seisSolGen.runBenchmarks('MatrixGen', 'matrices', 'build_log.txt')

if l_commandLineArguments.generateMatrixKernelsInitializationCode:
  l_seisSolGen.generateMatrixKernelsInitializationCode( i_pathToMatrices='matrices',
                                                        i_pathToOutputFile='generated_code/initialization/bind_matrix_kernels.hpp_include' )
  
if l_commandLineArguments.generateMatrixKernels:
  l_seisSolGen.generateMatrixKernels( i_pathToSeisSolGen = 'SeisSolGen/generator.exe',
                                      i_pathToMatrices = 'matrices',
                                      i_pathToOutputDirectory = 'generated_code/matrix_kernels' )

if l_commandLineArguments.generateMatrixUnitTests:
  l_unitTestGenerator.generateMatrixUnitTests( i_pathToOutputDirectory = 'generated_code/unit_tests' )

l_logger.printFinishMessage()
