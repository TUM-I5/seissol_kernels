#!/bin/env python

# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# @author Alexander Heinecke (alexander.heinecke AT mytum.de)
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
# Interface to the code generator SeisSolGen.
#

import os
import re
import subprocess
import tools.Logger as l_logger

# Run a single SeisSolGen benchmark
# TODO: This routine is outdated and references to an old version of the code generator.
#
# \param i_pathToSeisSolGen path to the code generator.
# \param i_pathToBuildLogFile path to the log file of the build processes.
# \param i_commandLineParameters command line parameters given to the code generator.
# \param i_archiveFileName name of the file, where the generated code is stored.
def runSingleBenchmark( i_pathToSeisSolGen,
                        i_buildLogFile,
                        i_commandLineParameters,
                        i_archiveFileName ):
  # initialize build log
  l_buildLog = ''

  # generate code C++-code for the current matrix
  l_bashCommand = './generator.exe ' + i_commandLineParameters + ' 1'
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  l_buildLog += l_bashProcess.communicate()[0]

  # archive the generated code
  l_bashCommand = 'cp generated_matmul.imp generated_code/'+i_archiveFileName
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  l_buildLog += l_bashProcess.communicate()[0]
        
  # clean the current build
  l_bashCommand = 'make benchclean'
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,)
  l_buildLog += l_bashProcess.communicate()[0]

  # build the benchmark for the current matrix
  l_bashCommand = 'make bench'
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,)
  l_buildLog += l_bashProcess.communicate()[0]

  i_buildLogFile.write(l_buildLog)

  # run the benchmark
  l_bashCommand = './benchmark.exe ' + i_commandLineParameters
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  print l_bashProcess.communicate()[0]

# Run SeisSolGen benchmarks
# TODO: This routine is outdated and references to an old version of the code generator.
#
# \param i_pathToSeisSolGen path to the code generator.
# \param i_pathToMatrices path to the matrices in Matrix Market format
def runBenchmarks( i_pathToSeisSolGen,
                   i_pathToMatrices,
                   i_pathToBuildLogFile ):
  
  # open the build log file
  l_buildLogFile = open(i_pathToBuildLogFile, 'w')

  log('** bulding the code generator')

  # log for the build commands
  l_buildLog = '*********** generator ***********\n'

  l_bashCommand = 'make clean'
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  l_buildLog += l_bashProcess.communicate()[0]

  l_bashCommand = 'make'
  l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd=i_pathToSeisSolGen, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  l_buildLog += l_bashProcess.communicate()[0]

  # write code generator log
  l_buildLogFile.write(l_buildLog)

  #old single right routine
  '''
  # star matrix
  l_file = 'starMatrix_3D_csc.mtx'
  l_pathToFile = i_pathToMatrices+'/'+l_file
  for l_degree in range(1,6):
    l_numberOfBasisFunctions = str((l_degree+1)*(l_degree+2)*(l_degree+3)/6)

    log('*** generating code and running benchmarks for '+l_file+', numerical order: ' + str(l_degree) )

    l_buildLog = '*********** ' + l_file + ' - '+ str(l_degree) + ' ***********\n'
    l_buildLogFile.write(l_buildLog)

    # command line parameter (filename, isCSR(0/1), nDenseRows, nDenseCols, nDenseLD, isLeft(0/1))
    l_commandLineParameters = ' ../'+l_pathToFile+' 0 '+l_numberOfBasisFunctions+' 9 '+l_numberOfBasisFunctions+' 0'

    # run benchmark for this matrix
    runSingleBenchmark( i_pathToSeisSolGen, l_buildLogFile, l_commandLineParameters, l_file.replace('_csc.mtx','')+'_'+str(l_degree)+'.imp' )
   '''

  # stiffness matrices
  # get and sort the matrix files
  l_matrixFiles = os.listdir(i_pathToMatrices)
  l_matrixFiles.sort()

  for l_degree in range(1,8):
    l_numberOfQuantities = '9'
    l_numberOfBasisFunctions = str((l_degree+1)*(l_degree+2)*(l_degree+3)/6)

    log( '*** generating code and running benchmarks for: '+str(l_degree)+' (order of basis), '+l_numberOfBasisFunctions+' (#basis functions), '+l_numberOfQuantities+' (#quantities)')

    l_buildLog = '*********** order - ' + str(l_degree) + ' ***********\n'
    l_buildLogFile.write(l_buildLog)

    # name of the stiffness matrices
    l_kXi =   'kXiDivMT_3D_'+str(l_degree)+'_csc.mtx'
    l_kEta =  'kEtaDivMT_3D_'+str(l_degree)+'_csc.mtx'
    l_kZeta = 'kZetaDivMT_3D_'+str(l_degree)+'_csc.mtx'

    # name of the star matrices
    l_aStar = 'starMatrix_3D_csc.mtx'
    l_bStar = 'starMatrix_3D_csc.mtx'
    l_cStar = 'starMatrix_3D_csc.mtx'
    
    # assert existance
    assert( l_kXi in l_matrixFiles )
    assert( l_kEta in l_matrixFiles )
    assert( l_kZeta in l_matrixFiles )

    assert( l_aStar in l_matrixFiles )
    assert( l_bStar in l_matrixFiles )
    assert( l_cStar in l_matrixFiles )

    # build complete paths
    l_kXi =   i_pathToMatrices+'/'+l_kXi
    l_kEta =  i_pathToMatrices+'/'+l_kEta
    l_kZeta = i_pathToMatrices+'/'+l_kZeta

    l_aStar = i_pathToMatrices+'/'+l_aStar
    l_bStar = i_pathToMatrices+'/'+l_bStar
    l_cStar = i_pathToMatrices+'/'+l_cStar

    # assert correct dimensions within the files
    l_matrixDimension = (open(l_kXi, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEta, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEta, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )

    l_matrixDimension = (open(l_aStar, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfQuantities )
    assert( l_matrixDimension[1] == l_numberOfQuantities )
    l_matrixDimension = (open(l_bStar, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfQuantities )
    assert( l_matrixDimension[1] == l_numberOfQuantities )
    l_matrixDimension = (open(l_cStar, 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfQuantities )
    assert( l_matrixDimension[1] == l_numberOfQuantities )

    # command line parameters (    filename-k-one, filename-star-one,
    #                              filename-k-two, filename-star-two
    #                              filename-k-three, filename-star-three
    #                              nDenseRows, nDenseCols, nDenseLD )
    l_commandLineParameters = '../' + l_kXi   + ' ' + '../' + l_aStar + ' ' +\
                              '../' + l_kEta  + ' ' + '../' + l_bStar + ' ' +\
                              '../' + l_kZeta + ' ' + '../' + l_cStar + ' ' +\
                              l_numberOfBasisFunctions + ' ' + l_numberOfQuantities + ' ' + l_numberOfBasisFunctions

    # run the actual benchmark
    runSingleBenchmark( i_pathToSeisSolGen, l_buildLogFile, l_commandLineParameters, 'gen_'+str(l_degree)+'.imp' )

  # old single left routine
  '''
  for l_file in l_matrixFiles:
    # get matrices in csc format
    if l_file.endswith('_csc.mtx'):
      # get stiffness matrices
      if re.match('k[A-Za-z0-9_]*DivM_3D[A-Za-z0-9_]*', l_file):
        log('*** generating code and running benchmarks for '+l_file)

        l_buildLog = '*********** ' + l_file + ' ***********\n'
        l_buildLogFile.write(l_buildLog)

        # generate full path
        l_pathToFile = i_pathToMatrices+'/'+l_file

        # read the dimension og the current matrix
        l_matrixDimension = open(l_pathToFile, 'r')
        l_matrixDimension = l_matrixDimension.readlines()[2]
        l_matrixDimension = l_matrixDimension.split()

        # command line parameters (    filename-k-one, filename-star-one,
        #                              filename-k-two, filename-star-two
        #                              filename-k-three, filename-star-three
        #                              nDenseRows, nDenseCols, nDenseLD )
        l_commandLineParameters = ' ../'+l_pathToFile+' 0 '+l_matrixDimension[0]+' 9 '+l_matrixDimension[1]+' 1'

        # run benchmark for this matrix
        # runSingleBenchmark( i_pathToSeisSolGen, l_buildLogFile, l_commandLineParameters, l_file.replace('_csc.mtx','')+'.imp' )
     '''

# Executes SeisSolGen with the given command line parameter
#
# \param i_pathToSeisSolGen location of SeisSolGen.
# \param i_commandLineParameters command line parameters used in the execution of SeisSolGen.
def executeSeisSolgen( i_pathToSeisSolGen,
                       i_commandLineParameters ):
    l_bashCommand = i_pathToSeisSolGen + ' ' + i_commandLineParameters
    l_logger.log('generating matrix kernels using '+l_bashCommand, 2)
    l_bashProcess = subprocess.Popen(l_bashCommand.split(), cwd='.', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    l_logger.log( l_bashProcess.communicate()[0], 3)

# Returns a list of matrix directories contataining star, stiffness and flux matrices.
#   Each matrix is specified by:
#     fileNameOfGeneratedKernel    - base name of generated matrix kernel
#     routineNameOfGeneratedKernel - name of the routine inside the generated kernel
#     pathToMatrixMarketFile       - location of the matrix in CSC-MatrixMarkeitFormat
#     multiplicationSide           - side on which the matrix appears within the kernel call
#                                    left: sparse * dense, right: dense * matrix
#     numberOfDenseRows            - #(dense rows) of the matrix
#     numberOfDenseColumns         - #(dense columns) of the matrix
#     sizeOfDenseLeadingDimension  - size of the leading dimension
#     aderZeroBlocks               - recursively appearing blocks of zeros appearing in the ADER time integration
#                                    True: matrix is part of the ADER time integration and zero blocks occur recursively
#                                    False: otherwise
#     add                          - True: C += A.B, False: C = A.B
#
# \param i_pathToMatrices path to the directory, which contains the matrices.
# \param i_numberOfQuantities number of quantities (elastics = 9, attenuation > 9)
# \param i_maximumDegreeOfBasisFunctions maximum order of the involved basis functions
def getSparseMatrices( i_pathToMatrices,
                       i_numberOfQuantities = 9,
                       i_maximumDegreeOfBasisFunctions = 6 ): #TODO: 7
  l_logger.log('getting matrices', 2)

  # list which holds the different matrix structures
  l_sparseMatrices = []

  # get and sort the matrix files
  l_matrixFiles = os.listdir(i_pathToMatrices)
  l_matrixFiles.sort()

  # convert input parameters to str
  l_numberOfQuantities = str(i_numberOfQuantities)

  ###
  # star matrices
  ###  

  # name of the star matrices
  l_starMatrix = 'starMatrix_3D_csc.mtx'

  # assert existance
  assert( l_starMatrix in l_matrixFiles )

  # build complete path
  l_starMatrix = i_pathToMatrices+'/'+l_starMatrix
  
  # assert correct dimensions within the file
  l_matrixDimension = (open(l_starMatrix, 'r').readlines()[2]).split()
  assert( l_matrixDimension[0] == l_numberOfQuantities )
  assert( l_matrixDimension[1] == l_numberOfQuantities )

  ###
  # stiffness and flux matrices
  ##
  for l_degree in range(1,i_maximumDegreeOfBasisFunctions):
    #each matric is a dictionary containing information about it
    l_matrix = dict()

    l_numberOfBasisFunctions = str((l_degree+1)*(l_degree+2)*(l_degree+3)/6)

    l_logger.log( 'adding star, stiffness and flux matrices to dictionaries for: '+str(l_degree)+' (order of basis), '+l_numberOfBasisFunctions+' (#basis functions), '+l_numberOfQuantities+' (#quantities)', 3)

    # name and filename of the stiffness matrices, TODO: avoid copy and paste code..
    l_kXi =   dict( matrixName = 'kXiDivM',   matrixMarketFileName='kXiDivM_3D_'   + str(l_degree) + '_csc.mtx' )
    l_kEta =  dict( matrixName = 'kEtaDivM',  matrixMarketFileName='kEtaDivM_3D_'  + str(l_degree) + '_csc.mtx' )
    l_kZeta = dict( matrixName = 'kZetaDivM', matrixMarketFileName='kZetaDivM_3D_' + str(l_degree) + '_csc.mtx' )
    
    l_kXiTransposed =   dict( matrixName = 'kXiDivMT',   matrixMarketFileName='kXiDivMT_3D_'   + str(l_degree) + '_csc.mtx' )
    l_kEtaTransposed =  dict( matrixName = 'kEtaDivMT',  matrixMarketFileName='kEtaDivMT_3D_'  + str(l_degree) + '_csc.mtx' )
    l_kZetaTransposed = dict( matrixName = 'kZetaDivMT', matrixMarketFileName='kZetaDivMT_3D_' + str(l_degree) + '_csc.mtx' )

    # assert existance
    assert( l_kXi['matrixMarketFileName']   in l_matrixFiles )
    assert( l_kEta['matrixMarketFileName']  in l_matrixFiles )
    assert( l_kZeta['matrixMarketFileName'] in l_matrixFiles )
    
    assert( l_kXiTransposed['matrixMarketFileName']   in l_matrixFiles )
    assert( l_kEtaTransposed['matrixMarketFileName']  in l_matrixFiles )
    assert( l_kZetaTransposed['matrixMarketFileName'] in l_matrixFiles )  

    # build complete paths
    l_kXi['pathToMatrixMarketFile']   = i_pathToMatrices+'/'+l_kXi['matrixMarketFileName']
    l_kEta['pathToMatrixMarketFile']  = i_pathToMatrices+'/'+l_kEta['matrixMarketFileName']
    l_kZeta['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_kZeta['matrixMarketFileName']
    
    l_kXiTransposed['pathToMatrixMarketFile']   = i_pathToMatrices+'/'+l_kXiTransposed['matrixMarketFileName']
    l_kEtaTransposed['pathToMatrixMarketFile']  = i_pathToMatrices+'/'+l_kEtaTransposed['matrixMarketFileName']
    l_kZetaTransposed['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_kZetaTransposed['matrixMarketFileName']

    # assert correct dimensions within the files
    l_matrixDimension = (open(l_kXi['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEta['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEta['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    
    l_matrixDimension = (open(l_kXiTransposed['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEtaTransposed['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
    l_matrixDimension = (open(l_kEtaTransposed['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
    assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
    assert( l_matrixDimension[1] == l_numberOfBasisFunctions )

    ###
    # Flux matrices
    ##
   
    # initialize empty list for flux matrices accounting for the elements contribution $F^{-,i}, \quad i \in \{1..4\}$
    l_fluxPlus = []   
    
    # initialize empty list for flux matrices accounting for the contribution of the neighboring elemetns $F^{+,i,j,h}, \quad i,j \in \{1..4\}, \; h \in \{1..3\}$
    l_fluxMinus = []

    # generate dictionaries, assert existance, build complete path and assert correct dimensions
    # for the flux matrices (i <-> local face, j <-> neighboring face, h <-> vertex combination)
    for l_localFace in range(0,4):
      ## element local flux matrix
      # generate dictionary
      #   matrixIds:
      #     0:  \f$ F^{-, 1} \f$
      #     1:  \f$ F^{-, 2} \f$
      #     2:  \f$ F^{-, 3} \f$
      #     3:  \f$ F^{-, 4} \f$
      l_fluxMinus = l_fluxMinus + [ dict( matrixName = 'fM'+str(l_localFace+1),
                                          matrixId = l_localFace,
                                          matrixMarketFileName='fM'+str(l_localFace+1)+'DivM_3D_'   + str(l_degree) + '_csc.mtx' ) ]
      
      # assert existance
      assert( l_fluxMinus[-1]['matrixMarketFileName']   in l_matrixFiles )

      # generate complete path
      l_fluxMinus[-1]['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxMinus[l_localFace]['matrixMarketFileName']

      # assert correct dimensions
      l_matrixDimension = (open(l_fluxMinus[-1]['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )

      for l_neighboringFace in range(0,4):
        for l_vertexCombination in range (0,3):
          ## neigboring element flux matrix
          # matrix ids
          #   4:  \f$ F^+{+, 1, 1, 1} \f$
          #   5:  \f$ F^+{+, 1, 1, 2} \f$
          #   6:  \f$ F^+{+, 1, 1, 3} \f$
          #   7:  \f$ F^+{+, 1, 1, 1} \f$
          #   [...]
          #  51:  \f$ F^+{+, 4, 4, 3} \f$
          # jump over elements contribute \f$ f^{-,i} \f$
          l_matrixId = 4
          # local face contribution to the id
          l_matrixId += 12*l_localFace
          # neighboring face contribution to the id
          l_matrixId += 3*l_neighboringFace
          # contribution of vertex orientation to the id
          l_matrixId += l_vertexCombination

          # generate dictionary
          l_multiIndex = str(l_localFace+1)+str(l_neighboringFace+1)+str(l_vertexCombination+1)
          l_fluxPlus = l_fluxPlus + [ dict( matrixName = 'fP'+l_multiIndex,
                                            matrixId   = l_matrixId,
                                            matrixMarketFileName='fP'+l_multiIndex+'DivM_3D_'   + str(l_degree) + '_csc.mtx' ) ]

          # assert existance
          assert( l_fluxPlus[-1]['matrixMarketFileName']   in l_matrixFiles )

          # generate complete path
          l_fluxPlus[-1]['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxPlus[-1]['matrixMarketFileName']

          # assert correct dimensions
          l_matrixDimension = (open(l_fluxPlus[-1]['pathToMatrixMarketFile'], 'r').readlines()[2]).split()
          assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
          assert( l_matrixDimension[1] == l_numberOfBasisFunctions )

    ###
    # Combine all matrices
    ##

    # add star matrices for the volume integration to list of matrices (TODO: this loop needs to be extended for attenuation)
    l_sparseMatrices.append( dict( fileNameOfGeneratedKernel    = 'star_matrices_3d.hpp_include',
                             name='volumeStarMatrix',
                             id = 59,
                             routineNameOfGeneratedKernel = 'generatedMatrixMultiplication_volumeStarMatrix_3D_'+l_numberOfQuantities+'_'+l_numberOfBasisFunctions,
                             pathToMatrixMarketFile       = l_starMatrix,
                             multiplicationSide           = '0',
                             numberOfDenseRows            = l_numberOfBasisFunctions,
                             numberOfDenseColumns         = l_numberOfQuantities,
                             sizeOfDenseLeadingDimension  = l_numberOfBasisFunctions,
                             aderZeroBlocks               = False,
                             add                          = True
                           )
                     )
    # star matrix multiplications in the time integration kernel
    l_sparseMatrices.append( dict(l_sparseMatrices[-1]) )
    l_sparseMatrices[-1]['name']                         = 'aderStarMatrix'
    l_sparseMatrices[-1]['routineNameOfGeneratedKernel'] = 'generatedMatrixMultiplication_aderStarMatrix_3D_'+l_numberOfQuantities+'_'+l_numberOfBasisFunctions
    l_sparseMatrices[-1]['aderZeroBlocks']               = True
    
  
    # add stiffness matrices to list of matrices
    l_stiffnessMatrixId = 53 # matrix id of \f$ K^\xi \f$, rest follows ascending
    for l_stiffnessMatrix in [l_kXi, l_kEta, l_kZeta, l_kXiTransposed, l_kEtaTransposed, l_kZetaTransposed]:
      l_sparseMatrices.append( dict( fileNameOfGeneratedKernel    = 'stiffness_matrices_3d.hpp_include',
                                     name                         = l_stiffnessMatrix['matrixName'],
                                     id                           = l_stiffnessMatrixId,
                                     routineNameOfGeneratedKernel = 'generatedMatrixMultiplication_'+l_stiffnessMatrix['matrixName']+'_'+l_numberOfQuantities+'_'+l_numberOfBasisFunctions+'',
                                     pathToMatrixMarketFile       = l_stiffnessMatrix['pathToMatrixMarketFile'],
                                     multiplicationSide           = '1',
                                     numberOfDenseRows            = l_numberOfBasisFunctions,
                                     numberOfDenseColumns         = l_numberOfQuantities,
                                     sizeOfDenseLeadingDimension  = l_numberOfBasisFunctions,
                                     aderZeroBlocks               = False,
                                     add                          = False
                                   )
                             )
      # increase matrix id
      l_stiffnessMatrixId = l_stiffnessMatrixId + 1

      # transposed stiffness matrices are used recursively in the ADER time integration
      if l_stiffnessMatrix in [l_kXiTransposed, l_kEtaTransposed, l_kZetaTransposed]:
        l_sparseMatrices[-1]['aderZeroBlocks'] = True

    # add flux matrices
    for l_fluxMatrix in l_fluxMinus+l_fluxPlus:
      l_sparseMatrices.append( dict( fileNameOfGeneratedKernel    = 'flux_matrices_3d.hpp_include',
                               name=l_fluxMatrix['matrixName'],
                               id = l_fluxMatrix['matrixId'],
                               routineNameOfGeneratedKernel = 'generatedMatrixMultiplication_'+l_fluxMatrix['matrixName']+'_'+l_numberOfQuantities+'_'+l_numberOfBasisFunctions+'',
                               pathToMatrixMarketFile       = l_fluxMatrix['pathToMatrixMarketFile'],
                               multiplicationSide           = '1',
                               numberOfDenseRows            = l_numberOfBasisFunctions,
                               numberOfDenseColumns         = l_numberOfQuantities,
                               sizeOfDenseLeadingDimension  = l_numberOfBasisFunctions,
                               aderZeroBlocks               = False,
                               add                          = False
                             )
                       )
  # done, return the matrices
  return l_sparseMatrices

# Returns a list of dense*dense matrix kernels for near dense matrices in SeisSol.
#   Each matrix is specified by:
#     fileNameOfGeneratedKernel    - base name of generated matrix kernel
#     routineNameOfGeneratedKernel - name of the routine inside the generated kernel
#     m                            - #(dense rows) of the left matrix
#     n                            - #(dense columns) of the right matrix
#     k                            - #(dense columns) of the left matrix
#     ldA                          - size of the leading dimension of the left matrix
#     ldB                          - size of the leading dimension of the right matrix
#     ldC                          - size of the leading dimnesion of the result matrix
#     add                          - True: C += A.B, False: C = A.B
#
#   Sketch (C = A.B)
#                    B
#            ****************
#            *     k x n    *
#            ****************
#      A             C
#   *******  ****************
#   *     *  *              *
#   *  m  *  *              *
#   *  x  *  *     m x n    *
#   *  k  *  *              *
#   *     *  *              *
#   *******  ****************
#
# \param i_numberOfQuantities number of quantities (elastics = 9, attenuation > 9)
# \param i_maximumDegreeOfBasisFunctions maximum order of the involved basis functions
# \return dictionary conatinig the dense matrices described abover.
def getDenseMatrices( i_numberOfQuantities = 9,
                      i_maximumDegreeOfBasisFunctions = 6 ):
  
  # list which holds the different matrix structures
  l_denseMatrices = []

  # file where the generated code is stored
  l_fileNameOfGeneratedKernel = 'dense_matrices.hpp_include'

  # iterate of basis functions degree
  for l_degree in range(1,i_maximumDegreeOfBasisFunctions):
    # compute number of basis functions
    l_numberOfBasisFunctions = str((l_degree+1)*(l_degree+2)*(l_degree+3)/6)

    l_logger.log( 'adding dense matrices: '+str(l_degree)+' (degree of basis), '+l_numberOfBasisFunctions+' (#basis functions), '+str(i_numberOfQuantities)+' (#quantities)', 3)

    # generate kernels for flux matrices
    #   m = #(basis functions)
    #   n = #(quantities)
    #   k = #(basis functions)
    l_m = l_numberOfBasisFunctions
    l_n = i_numberOfQuantities
    l_k = l_numberOfBasisFunctions

    # where to store the dense matrix kernels
    l_fileNameOfGeneratedKernel    = 'dense_matrices.hpp_include'

    # generate a dense kernel for boundary and volume (dense) and one for the time integration (denseCK)
    for l_denseType in ['dense', 'denseCK']:
      # generate function name
      l_routineNameOfGeneratedKernel = 'generatedMatrixMultiplication_'+ 'dense'
      if l_denseType == 'denseCK':
         l_routineNameOfGeneratedKernel =  l_routineNameOfGeneratedKernel + 'Ader'

      l_routineNameOfGeneratedKernel = l_routineNameOfGeneratedKernel + '_'+str(l_m)+'_'+str(l_n)+'_'+str(l_k)

      # Add matrix to dictionary
      #   TODO: We use leading dimensions equal to standard column wise storage
      l_denseMatrices = l_denseMatrices + [ dict(\
                                              fileNameOfGeneratedKernel    = l_fileNameOfGeneratedKernel, \
                                              routineNameOfGeneratedKernel = l_routineNameOfGeneratedKernel, \
                                              m          = l_m, \
                                              n          = l_n, \
                                              k          = l_k, \
                                              ldA        = l_m, \
                                              ldB        = l_k, \
                                              ldC        = l_m, \
                                              denseType  = l_denseType, \
                                              add        = False
                                          )]

    # generate kernels for the flux solver
    # m = #(basis functions)
    # n = #(quantities)
    # k = #(quantities)   
    l_m = l_numberOfBasisFunctions
    l_n = i_numberOfQuantities
    l_k = i_numberOfQuantities

    # generate function name
    l_routineNameOfGeneratedKernel = 'generatedMatrixMultiplication_'+ 'dense_'+str(l_m)+'_'+str(l_n)+'_'+str(l_k)

    # Add matrix to dictionary
    #   TODO: We use leading dimensions equal to standard column wise storage
    l_denseMatrices = l_denseMatrices + [ dict(\
                                            fileNameOfGeneratedKernel    = l_fileNameOfGeneratedKernel, \
                                            routineNameOfGeneratedKernel = l_routineNameOfGeneratedKernel, \
                                            m          = l_m, \
                                            n          = l_n, \
                                            k          = l_k, \
                                            ldA        = l_m, \
                                            ldB        = l_k, \
                                            ldC        = l_m, \
                                            denseType  = 'dense', \
                                            add        = True
                                        )]

  return l_denseMatrices

# Generates matrix kernels.
#
# @param i_pathToSeisSolGen path to executable of the code generator.
# @param i_pathToMatrices path to directory, which contains the matrices in CSC-MatrixMarket-format.
# @param i_pathToOutputDirectory path to directory where the generated code will be written to.
def generateMatrixKernels( i_pathToSeisSolGen,
                           i_pathToMatrices,
                           i_pathToOutputDirectory ):
  l_logger.log('generating matrix kernels with'+
               '\nSeisSolGen='+i_pathToSeisSolGen+
               '\n matrix directory='+i_pathToMatrices+
               '\n output directory='+i_pathToOutputDirectory)

  # get dense matrices
  l_denseMatrices = getDenseMatrices()

  # get sparse matrices
  l_sparseMatrices = getSparseMatrices(i_pathToMatrices)

  l_fileNames = [ 'dense_matrices.hpp_include',
                  'star_matrices_3d.hpp_include',
                  'flux_matrices_3d.hpp_include',
                  'stiffness_matrices_3d.hpp_include' ]

  # create new files and write header information
  for l_file in l_fileNames:
    # output file
    l_pathToOutputFile = i_pathToOutputDirectory+'/'+l_file
  
    # remove old file contents
    open(l_pathToOutputFile, 'w').close()

    # write license
    l_logger.writeFileHeader(l_pathToOutputFile, '// ')
    
    # add include guards, add include for intrinsics
    l_includeGuardName = re.sub("[^a-zA-Z]","",  l_file).upper() # remove non letters, uppercase
    l_includeCommand = '#ifndef ' + l_includeGuardName + '\n' \
                       '#define ' + l_includeGuardName + '\n\n' \
                       '#if defined( __SSE3__) || defined(__MIC__)\n#include <immintrin.h>\n#endif\n\n'
    l_outputFile = open(l_pathToOutputFile ,'a')

    # TODO: Remove this code once we have different compile units
    if( l_file == 'flux_matrices_3d.hpp_include' ):
      l_outputFile.write('#ifndef __MIC__\n')

    l_outputFile.write(l_includeCommand)
    l_outputFile.close()

  # generate code for all dense matrices
  for l_matrix in l_denseMatrices:
    # output file
    l_pathToOutputFile = i_pathToOutputDirectory+'/'+l_matrix['fileNameOfGeneratedKernel']

    l_type = l_matrix['denseType']
    l_commandLineParameters = l_type + ' ' + l_pathToOutputFile+\
                              ' '+l_matrix['routineNameOfGeneratedKernel']+\
                              ' '+str(l_matrix['m'])+\
                              ' '+str(l_matrix['n'])+\
                              ' '+str(l_matrix['k'])+\
                              ' '+str(l_matrix['ldA'])+\
                              ' '+str(l_matrix['ldB'])+\
                              ' '+str(l_matrix['ldC'])+\
                              ' '+str(int(l_matrix['add']))
    # generate code C++-code for the current matrix
    executeSeisSolgen( i_pathToSeisSolGen,\
                       l_commandLineParameters )

  # generate code for all sparse matrices
  for l_matrix in l_sparseMatrices:
    # output file
    l_pathToOutputFile = i_pathToOutputDirectory+'/'+l_matrix['fileNameOfGeneratedKernel']

    # check if the matrix qualifies for the reduced ADER time integration, which avoids multiplications with zero blocks for polynomial degrees above zero
    if l_matrix['aderZeroBlocks'] == False:    
      l_type = 'sparse'
    else:
      l_type = 'sparseCK'

    l_commandLineParameters = l_type + ' ' + l_pathToOutputFile+\
                              ' '+l_matrix['routineNameOfGeneratedKernel']+\
                              ' '+l_matrix['pathToMatrixMarketFile']+\
                              ' '+l_matrix['multiplicationSide']+\
                              ' '+l_matrix['numberOfDenseRows']+\
                              ' '+l_matrix['numberOfDenseColumns']+\
                              ' '+l_matrix['sizeOfDenseLeadingDimension']+\
                              ' '+str(int(l_matrix['add']))
    # generate code C++-code for the current matrix
    executeSeisSolgen( i_pathToSeisSolGen,\
                       l_commandLineParameters )

  # TODO: Remove this code once we have different compile units
  l_pathToOutputFile = i_pathToOutputDirectory+'/flux_matrices_3d.hpp_include'
  l_outputFile = open(l_pathToOutputFile ,'a')
  l_outputFile.write('#else \n#warning sparse flux matrices are not defined\n#endif\n')
  l_outputFile.close()

  # close include guards
  for l_file in l_fileNames:
    # output file
    l_pathToOutputFile = i_pathToOutputDirectory+'/'+l_file
    l_outputFile = open(l_pathToOutputFile ,'a')
    l_outputFile.write('\n#endif')
    l_outputFile.close()

# Generates code, which initializes the flux matrix function pointers for the boundary integrator.
#
# @param i_pathToMatrices path to directory containing the matrices.
# @param i_pathToOutputFile path to file, which will hold the initialization code.
def generateBoundaryMatrixKernelsInitializationCode( i_pathToMatrices,
                                          i_pathToOutputFile ):
  l_logger.log('generating flux matrix initialization code' )

  # holds the generated source code
  l_sourceCode = ''

  # remove old file contents
  open(i_pathToOutputFile, 'w').close()

  # write license
  l_logger.writeFileHeader(i_pathToOutputFile, '// ')
  
  # get sparse matrices of degree 2 only, because we are interested in the base names only
  l_matrices = getSparseMatrices( i_pathToMatrices = i_pathToMatrices,
                                  i_numberOfQuantities = 9,
                                  i_maximumDegreeOfBasisFunctions = 2 )

  # use dense kernel if not sparse
  l_sourceCode = l_sourceCode + '// preprocessor helping function to convert integer constants to strings\n'\
                              +  '#define CONCAT_HELPER_4(a,b,c,d) a ## b ## c ## d\n'\
                              +  '#define CONCAT_4(a,b,c,d) CONCAT_HELPER_4(a,b,c,d)\n'\
                              +  '#define CONCAT_HELPER_6(a,b,c,d,e,f) a ## b ## c ## d ## e ## f\n'\
                              +  '#define CONCAT_6(a,b,c,d,e,f) CONCAT_HELPER_6(a,b,c,d, e, f)\n'\
                              +  '#ifndef __MIC__\n'\
                              +  '// dense flux kernel\n'\
                              +  'if( i_sparse == false && i_id < 52 ){\n'\
                              +  '  m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );\n'\
                              +  '}\n'\
  # iterate of flux matrice
  for l_matrix in l_matrices:
    if 'fP' in l_matrix['name'] or 'fM' in l_matrix['name']:
      l_sourceCode = l_sourceCode + '// sparse flux kernel\n'\
                                  + 'else if( i_id == ' + str(l_matrix['id']) + '){\n'\
                                  + '  m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_'+ l_matrix['name'] +'_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );\n'\
                                  + '}\n'

  # add code for dense jacobian
  l_sourceCode = l_sourceCode + '// dense jacobian kernel\n'\
                              + 'else if( i_id == 52 ){\n'\
                              + '  m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFVARIABLES );\n'\
                              + '}\n'

  # we didn't find our matrix
  l_sourceCode = l_sourceCode + '// something went wrong: we dont have a kernel for this matrix\n'\
                              + 'else{\n'\
                              + '  assert(false);\n'\
                              + '}\n'
  
  # add dense only for Intel MIC
  l_sourceCode = l_sourceCode + '#else\n'\
                              + '// dense only for Intel MIC\n'\
                              + 'if( i_id == 52 ){ // jacobian\n'\
                              + '  m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFVARIABLES );\n'\
                              + '}\n'\
                              + 'else{ // flux matrix\n'\
                              + '  // assert only sparse flux kernels are requested\n'\
                              + '  assert(i_sparse == false);\n'\
                              + '  m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );\n'\
                              + '}\n'\
                              + '#endif\n'

  # undefine helper functions
  l_sourceCode = l_sourceCode +  '#undef CONCAT_HELPER_4\n'\
                              +  '#undef CONCAT_4\n'\
                              +  '#undef CONCAT_HELPER_6\n'\
                              +  '#undef CONCAT_6\n'

  # write code to file
  l_file = open(i_pathToOutputFile ,'a')
  l_file.write(l_sourceCode)
  l_file.close()
