#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# @author Alexander Heinecke (alexander.heinecke AT mytum.de)
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
#
# Interface to the code generator SeisSolGen.
#

import re
import subprocess
import tempfile
import tools.Logger as l_logger
import MatrixConverter
import xmltodict
import logging
import pprint
import os
import scipy.io
import numpy
import sys

import Configuration

class SeisSolGen:
  m_configuration = 0
  m_matrixSetup = 0

  # Constuctor
  #
  # @param i_configuration Configuration
  def __init__( self,
                i_matrixSetup ):
    self.m_configuration = i_matrixSetup.m_configuration
    self.m_matrixSetup = i_matrixSetup

  # Executes SeisSolGen with the given command line parameter
  #
  # @param i_pathToSeisSolGen location of SeisSolGen.
  # @param i_commandLineParameters command line parameters used in the execution of SeisSolGen.
  def executeSeisSolgen( self,
                         i_pathToSeisSolGen,
                         i_commandLineParameters ):
      l_bashCommand = i_pathToSeisSolGen + ' ' + i_commandLineParameters
      l_logger.log('generating matrix kernels using '+l_bashCommand, 2)
      l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd='.', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      l_logger.log( l_bashProcess.communicate()[0], 3)

  # Generates matrix kernels.
  #
  # @param i_pathToSeisSolGen path to executable of the code generator.
  # @param i_pathToMatrices path to directory, which contains the matrices in MatrixMarket-format.
  # @param i_pathToOutputDirectory path to directory where the generated code will be written to.
  def generateMatrixKernels( self ):

    l_outputDir = self.m_configuration.m_pathToGeneratedCodeDir+'/matrix_kernels'

    l_logger.log('generating matrix kernels with'+
                 '\nSeisSolGen='+self.m_configuration.m_matricesDir+
                 '\n matrix directory='+self.m_configuration.m_pathToGemmCodeGenerator+
                 '\n output directory='+l_outputDir)

    # create directory if not existent
    if not os.path.exists(l_outputDir):
      os.makedirs(l_outputDir)

    for l_precision in ['s', 'd']:
      l_fileNames = []

      for l_architecture in self.m_configuration.m_architectures:
        l_baseName = l_precision + 'gemm_' + l_architecture

        l_fileNames = l_fileNames + [l_baseName+".h"]
        l_fileNames = l_fileNames + [l_baseName+".cpp"]

      for l_architecture in self.m_configuration.m_architectures:
        l_baseName = 'sparse_' + l_precision + l_architecture

        l_fileNames = l_fileNames + [l_baseName+".h"]
        l_fileNames = l_fileNames + [l_baseName+".cpp"]


      # create new files and write header information
      for l_file in l_fileNames:
        # output file
        l_pathToOutputFile = l_outputDir+'/'+l_file
      
        # remove old file contents
        open(l_pathToOutputFile, 'w').close()

        # write license
        l_logger.writeFileHeader(l_pathToOutputFile, '// ')
        
        # add include guards, add include for intrinsics
        l_includeGuardName = re.sub("[^a-zA-Z]","",  l_file).upper() # remove non letters, uppercase
        l_includeCommand = '#ifndef ' + l_includeGuardName + '\n' \
                           '#define ' + l_includeGuardName + '\n\n' \
                           '#if defined( __SSE3__) || defined(__MIC__)\n#include <immintrin.h>\n#endif\n\n'\
                           '#include <cstddef>\n'\
                           '#ifndef NDEBUG\n'\
                           'extern long long libxsmm_num_total_flops;\n'\
                           '#endif\n\n'
        l_outputFile = open(l_pathToOutputFile ,'a')

        l_outputFile.write(l_includeCommand)
        l_outputFile.close()

      # get dense matrices
      l_denseMatrices = self.m_matrixSetup.getDenseMatrices( i_architectures = self.m_configuration.m_architectures, i_precision = l_precision )

      # get sparse matrices
      l_sparseMatrices = self.m_matrixSetup.getSparseMatrices( i_architectures = self.m_configuration.m_architectures, i_precision = l_precision )

      # generate code for all matrices
      for l_matrix in l_denseMatrices + l_sparseMatrices:
        # write header
        l_header = open(l_outputDir+'/'+l_matrix['fileNameOfGeneratedKernel']+'.h', 'a')
        l_header.write( 'void ' + l_matrix['routine_name']+'(const real *i_A, const real *i_B, real *io_C,'\
                                                            'const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );\n' )
        l_header.close()

        # output file
        l_pathToOutputFile = l_outputDir+'/'+l_matrix['fileNameOfGeneratedKernel']+'.cpp'

        if( l_matrix['type'] == 'dense' ):
          l_commandLineParameters =     l_matrix['type']          +\
                                    ' '+l_pathToOutputFile        +\
                                    ' '+l_matrix['routine_name']  +\
                                    ' '+str(l_matrix['m'])        +\
                                    ' '+str(l_matrix['n'])        +\
                                    ' '+str(l_matrix['k'])        +\
                                    ' '+str(l_matrix['ld_a'])     +\
                                    ' '+str(l_matrix['ld_b'])     +\
                                    ' '+str(l_matrix['ld_c'])     +\
                                    ' 1 1'                        +\
                                    ' '+str(int(l_matrix['add']))
          # 1 1: we always rely on alignment based on LDA and LDC
        else:
          # convert to temporary csc file
          l_cscFile = tempfile.NamedTemporaryFile(); l_csrFile = tempfile.NamedTemporaryFile();
          MatrixConverter.MatrixConverter.convertFullToSparse( i_pathToFullMatrix = l_matrix['matrix_market'],
                                                               i_cscFile          = l_cscFile,
                                                               i_csrFile          = l_csrFile )


          l_commandLineParameters =     l_matrix['type']          +\
                                    ' '+l_pathToOutputFile        +\
                                    ' '+l_matrix['routine_name']  +\
                                    ' '+l_cscFile.name            +\
                                    ' '+str(l_matrix['m'])        +\
                                    ' '+str(l_matrix['n'])        +\
                                    ' '+str(l_matrix['k'])        +\
                                    ' '+str(l_matrix['ld_a'])     +\
                                    ' '+str(l_matrix['ld_b'])     +\
                                    ' '+str(l_matrix['ld_c'])     +\
                                    ' 1 1'                        +\
                                    ' '+str(int(l_matrix['add']))
          # 1 1: we always rely on alignment based on LDA and LDC

        l_commandLineParameters += ' '+str(l_matrix['arch']) + ' ' + str(l_matrix['prefetch']) + ' ' + l_precision.upper() + 'P'

        # generate code C++-code for the current matrix
        self.executeSeisSolgen( self.m_configuration.m_pathToGemmCodeGenerator,\
                                l_commandLineParameters )

      # close include guards
      for l_file in l_fileNames:
        l_header = open(l_outputDir+'/'+l_file, 'a' )
        l_header.write('#endif\n')
        l_header.close()

    # write sparse global incluce
    l_globalInclude = l_outputDir+'/dense.h'
    open( l_globalInclude, 'w').close()
    l_logger.writeFileHeader(l_globalInclude, '// ')
    l_globalInclude = open( l_globalInclude, 'a')
    for l_architecture in self.m_configuration.m_architectures:
      for l_precision in ['s', 'd']:
        l_command =   '#ifdef ' + (str(l_precision+l_architecture)).upper()     + '\n' \
                    + '#include <matrix_kernels/' + l_precision + 'gemm_' + l_architecture + '.h>'+ '\n' \
                    + '#endif\n'
        l_globalInclude.write( l_command )
    l_globalInclude.close()

    # write sparse global incluce
    l_globalInclude = l_outputDir+'/sparse.h'
    open( l_globalInclude, 'w').close()
    l_logger.writeFileHeader(l_globalInclude, '// ')
    l_globalInclude = open( l_globalInclude, 'a')
    for l_architecture in self.m_configuration.m_architectures:
      for l_precision in ['s', 'd']:
        l_command =   '#ifdef ' + (str(l_precision+l_architecture)).upper()     + '\n' \
                    + '#include <matrix_kernels/sparse_' + l_precision +  l_architecture+ '.h>'+ '\n' \
                    + '#endif\n'
        l_globalInclude.write( l_command )
    l_globalInclude.close()

  # Generates code, which initializes the matrix function pointers for in the time, volume and flux kernel.
  def generateMatrixKernelsInitializationCode( self ):
    l_logger.log('generating dense matrix initialization code' )

    l_outputDir = self.m_configuration.m_pathToGeneratedCodeDir+'/initialization'

    # create directory if not existent
    if not os.path.exists(l_outputDir):
      os.makedirs(l_outputDir)

    l_globalInclude = ''

    for l_precision in [ 's', 'd' ]:
      for l_architecture in self.m_configuration.m_architectures:
        logging.info( 'generating bindings for '+l_precision+l_architecture )

        l_outFile = l_outputDir + '/' + 'bind_' + l_precision + l_architecture + '.h'

        l_globalInclude = l_globalInclude + '#ifdef ' + (str(l_precision+l_architecture)).upper() + '\n'                              \
                                          + '#include <initialization/' + 'bind_' + l_precision + l_architecture + '.h' + '>\n' \
                                          + '#endif\n'

        # holds the generated source code
        l_sourceCode = ''

        # remove old file contents
        open(l_outFile, 'w').close()

        # write license
        l_logger.writeFileHeader(l_outFile, '// ')

        l_alignment = self.m_configuration.m_alignments[l_architecture]

        l_sourceCode = l_sourceCode + '#if ALIGNMENT!=' + str(l_alignment) + '\n' \
                                    + '#error alignment-architecture mismatch\n'   \
                                    + '#endif\n\n'

        l_pathToSparseDenseSwitch = self.m_configuration.m_pathToSparseDenseConfigs + "/" + l_precision + l_architecture + ".xml"

        # iterate over convergence orders
        for l_order in range(2, self.m_configuration.m_maximumDegree + 2):
          l_sourceCode = l_sourceCode + '#if CONVERGENCE_ORDER==' + str(l_order) + '\n'

          # all sparse matrices
          l_allSparse = []

          # iterate over kernels
          for l_kernel in ['time', 'volume', 'boundary']:
            if( l_kernel == 'time' ):
              l_numberOfBinds = 4 * (l_order-1)

              # get the order and name of the matrices
              l_matrixNames = {}
              for l_bind in range(0, l_numberOfBinds):
                l_matrixNames[l_bind] = { 'order': (l_numberOfBinds - l_bind - 1)/4 + 2,
                                          'name' : self.m_configuration.m_matrixBinds['time'][l_bind%4] }

                # decrease the order for star matrices which operate on reduced time derivatives
                if (l_bind+1)%4 == 0:
                  l_matrixNames[l_bind]['order'] = l_matrixNames[l_bind]['order'] - 1

              l_sparse      = self.m_matrixSetup.getSparseTimeMatrices(            i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_pathToSparseDenseSwitch = l_pathToSparseDenseSwitch,
                                                                                   i_precision               = l_precision )

              l_globalDgemm = self.m_matrixSetup.getDenseStiffTimeMatrices(        i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )
              l_localDgemm  = self.m_matrixSetup.getDenseStarSolverMatrices(       i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = reversed(range(l_order-1)),
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )
              assert( len(l_globalDgemm) == len(l_localDgemm) )

            if( l_kernel == 'volume' ):
              l_numberOfBinds = 4

              # get the order and name of the matrices
              l_matrixNames = {}
              for l_bind in range(0, l_numberOfBinds):
                l_matrixNames[l_bind] = { 'order': l_order,
                                          'name' : self.m_configuration.m_matrixBinds['volume'][l_bind] }

              l_sparse =  self.m_matrixSetup.getSparseVolumeAndFluxMatrices(       i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_pathToSparseDenseSwitch = l_pathToSparseDenseSwitch,
                                                                                   i_integrationKernels      = ['volume'],
                                                                                   i_precision               = l_precision )

              l_globalDgemm = self.m_matrixSetup.getDenseStiffVolumeMatrices(      i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )

              l_localDgemm  = self.m_matrixSetup.getDenseStarSolverMatrices(       i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )

            if( l_kernel == 'boundary' ):
              l_numberOfBinds = 54

              # get the order and name of the matrices
              l_matrixNames = {}
              for l_bind in range(0, l_numberOfBinds):
                l_matrixNames[l_bind] = { 'order': l_order,
                                          'name' : self.m_configuration.m_matrixBinds['boundary'][l_bind] }

              l_sparse      =  self.m_matrixSetup.getSparseVolumeAndFluxMatrices(  i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_pathToSparseDenseSwitch = l_pathToSparseDenseSwitch,
                                                                                   i_integrationKernels      = ['boundary'],
                                                                                   i_precision               = l_precision )

              l_globalDgemm = self.m_matrixSetup.getDenseFluxMatrices(             i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )

              if l_architecture in ['wsm', 'snb', 'hsw']:
                l_fluxMatrix_prefetch = 'BL2viaC'
              else:
                l_fluxMatrix_prefetch = 'pfsigonly'

              l_globalDgemm = l_globalDgemm + self.m_matrixSetup.getDenseFluxMatrices(             i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision,
                                                                                   i_prefetch                = l_fluxMatrix_prefetch )

              l_localDgemm  = self.m_matrixSetup.getDenseStarSolverMatrices(       i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision )

              if l_architecture in ['wsm', 'snb', 'hsw']:
                l_starSolver_prefetch = 'pfsigonly'
              else:
                l_starSolver_prefetch = 'pfsigonly'

              l_localDgemm  = l_localDgemm + self.m_matrixSetup.getDenseStarSolverMatrices(       i_alignment               = l_alignment,
                                                                                   i_degreesOfBasisFunctions = [l_order-1],
                                                                                   i_numberOfQuantities      = 9,
                                                                                   i_precision               = l_precision,
                                                                                   i_prefetch                = l_starSolver_prefetch )

            l_sourceCode = l_sourceCode + '\n#ifdef ' + l_kernel.upper() + '_KERNEL\n'

            # store current sparse matrices
            l_allSparse = l_allSparse + l_sparse

            # GEMM matrices
            l_gemmId = 0

            for l_bind in range( 0, l_numberOfBinds ):
              # check if we have a matching sparse matrix
              l_match = [ l_sparseMatrix for l_sparseMatrix in l_sparse if l_sparseMatrix['bind'] == l_bind ]
              assert( len( l_match ) < 2 )

              # get matrix bind of the local matrix
              if(   l_kernel == 'time' or l_kernel == 'volume' ):
                l_localBind = self.m_configuration.m_matrixBinds[l_kernel]['starMatrix']
              elif( l_kernel == 'boundary' ):
                l_localBind = self.m_configuration.m_matrixBinds['boundary']['fluxSolver']
              else:
                assert( False )           

              # get nonzeros
              l_matrixOrder = l_matrixNames[l_bind]['order']
              l_nonZeros = self.m_configuration.m_nonZeros[l_matrixOrder][l_matrixNames[l_bind]['name']]

              # this a local matrix
              if( (l_bind+1) % (l_localBind+1) == 0 ):
                l_gemmMatrix = l_localDgemm[l_gemmId]

                # increase gemm id for time kernel
                l_gemmId = l_gemmId + 1
                if (l_kernel != 'boundary' ):
                  l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*self.m_matrixSetup.getNumberOfBasisFunctions(l_matrixOrder)*2)   + ';\n'
              else:
                l_gemmMatrix = l_globalDgemm[l_gemmId]
                if (l_kernel != 'boundary' ):
                  l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*l_gemmMatrix['n']*2)   + ';\n'

              # sparse setup
              if( len(l_match) == 1 ):
                if (l_kernel == 'boundary' ):
                  l_gemmMatrix = l_globalDgemm[0]
                  l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*l_gemmMatrix['n']*2)   + ';\n'

                l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_match[0]['flops'])   + ';\n'
                l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_match[0]['routine_name'] + ';\n'
              # dense setup
              else:
                if (l_kernel == 'boundary' ):
                  if (l_bind < 4):
                    l_gemmMatrix = l_globalDgemm[0]
                    l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*l_gemmMatrix['n']*2)   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_gemmMatrix['flops'])   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                  elif (l_bind == 52):
                    l_gemmMatrix = l_localDgemm[0]
                    l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*self.m_matrixSetup.getNumberOfBasisFunctions(l_matrixOrder)*2)   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_gemmMatrix['flops'])   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                  elif (l_bind > 3 and l_bind < 52):
                    l_gemmMatrix = l_globalDgemm[1]
                    l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*l_gemmMatrix['n']*2)   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_gemmMatrix['flops'])   + ';\n'
                    l_sourceCode = l_sourceCode + '#ifdef ENABLE_MATRIX_PREFETCH\n'
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                    l_sourceCode = l_sourceCode + '#else\n'
                    l_gemmMatrix = l_globalDgemm[0]
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                    l_sourceCode = l_sourceCode + '#endif\n'
                  elif (l_bind == 53):
                    l_gemmMatrix = l_localDgemm[1]
                    l_sourceCode = l_sourceCode + 'm_nonZeroFlops[' + str(l_bind) + '] = ' + str(l_nonZeros*self.m_matrixSetup.getNumberOfBasisFunctions(l_matrixOrder)*2)   + ';\n'
                    l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_gemmMatrix['flops'])   + ';\n'
                    l_sourceCode = l_sourceCode + '#ifdef ENABLE_MATRIX_PREFETCH\n'
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                    l_sourceCode = l_sourceCode + '#else\n'
                    l_gemmMatrix = l_localDgemm[0]
                    l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                    l_sourceCode = l_sourceCode + '#endif\n'
                  else:
                    sys.exit("generation error in binding for boundary kernel")
                else:
                  l_sourceCode = l_sourceCode + 'm_matrixKernels[' + str(l_bind) + '] = ' + l_gemmMatrix['routine_name'] + ';\n'
                  l_sourceCode = l_sourceCode + 'm_hardwareFlops[' + str(l_bind) + '] = ' + str(l_gemmMatrix['flops'])   + ';\n'


            l_sourceCode = l_sourceCode + '#endif\n\n'

          #
          # setup sparse-dense decision matrix at compile time
          #
          l_sourceCode = l_sourceCode + '#ifdef SPARSE_SWITCH\n'
          for l_matrix in range(0, 60):
            # check if we have a matching sparse matrix
            l_match = [ l_sparseMatrix for l_sparseMatrix in l_allSparse if l_sparseMatrix['name'] == self.m_configuration.m_matrixNames[l_matrix] ]

            # sparse setup
            if( len(l_match) > 0 ):
               l_nnz = len(numpy.nonzero(scipy.io.mmread( l_match[0]['matrix_market'] ))[0])
               l_sourceCode = l_sourceCode + 'm_sparseSwitch[' + str(l_matrix) + '] = ' + str(l_nnz) + '; \n'
            # dense setup
            else:
               l_sourceCode = l_sourceCode + 'm_sparseSwitch[' + str(l_matrix) + '] = -1; \n'
          l_sourceCode = l_sourceCode + '#endif\n'

          # define a compile switch if the star matrix is sparse
          if( len( [ l_sparseMatrix for l_sparseMatrix in l_allSparse if 'starMatrix' in l_sparseMatrix['name'] ] ) > 0 ):
            l_sourceCode = l_sourceCode + '\n#define STAR_NNZ 24\n'
          else:
            l_sourceCode = l_sourceCode + '\n#define STAR_NNZ 81\n'

          l_sourceCode = l_sourceCode + '\n#endif\n\n'

        # write code to file
        l_outFile = open(l_outFile ,'a')
        l_outFile.write(l_sourceCode)
        l_outFile.close()

    # write global header
    l_outFile = l_outputDir + '/' + 'bind.h'
    open(l_outFile ,'w').close()
    l_logger.writeFileHeader(l_outFile, '// ')
    l_outFile = open(l_outFile ,'a')
    l_outFile.write(l_globalInclude)
    l_outFile.close()

    # setup precision switch
    l_sourceCode = ''
    for l_precision in ['s', 'd']:
       for l_architecture in self.m_configuration.m_architectures:
         l_sourceCode = l_sourceCode + '\n#ifdef ' + str( l_precision + l_architecture ).upper() + "\n"
         if l_precision == 's':
           l_sourceCode = l_sourceCode + "#define SINGLE_PRECISION\n"
         else:
           l_sourceCode = l_sourceCode + "#define DOUBLE_PRECISION\n"
         l_sourceCode = l_sourceCode + "#endif\n"

    # precision switch
    l_outFilePath = l_outputDir + '/' + 'precision.h'
    open(l_outFilePath ,'w').close()
    l_outFile = open(l_outFilePath ,'a')
    l_outFile.write( '#if 0\n' )
    l_outFile.close()
    l_logger.writeFileHeader(l_outFilePath, '// ')
    l_outFile = open(l_outFilePath ,'a')
    l_outFile.write( '#endif\n' )
    l_outFile.write(l_sourceCode)
    l_outFile.close()

