#!/bin/env python
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
#
# Converts different matrices to CSC and CSR.
#

from scipy.io.mmio import *
from scipy.sparse import *

from matplotlib import pyplot, colors
from pylab import savefig
from numpy import arange, nonzero, sort

import tools.Logger as l_logger

#from elementtree.SimpleXMLWriter import XMLWriter
import tools.SimpleXMLWriter as XMLWriter

import os

class MatrixConverter():
###
### Class Variables
###
  # stores preprocessor commands
  preProcessorCommands2D = {'c++':'', 'fortran':''}
  preProcessorCommands3D = {'c++':'', 'fortran':''}

###
### Local functions
###

  # Generates a single preprocessor command for a given matrix,
  # which checks again a specified dimension and contains its number of nonzeros.
  #
  # \param i_pathToFullMatrix path to the full matrix (format: matrix market).
  # \param i_baseName base name of the matrix.
  # \param i_dimensionSize size of the dimesion to check again (i.e. 9 for #variables in elastics).
  # \param i_dimensionName name of the dimension to check again (i.e. NUMBEROFBASISFUNCTIONS for stiffness matrices).
  def generatePreProcessorCode( self,
                                i_pathToFullMatrix,
                                i_baseName,
                                i_dimensionSize,
                                i_dimensionName ):
    l_logger.log('generating preprocessor code for matrix: '+i_pathToFullMatrix, 2)

    # read the full matrix
    l_fullMatrix = mmread(i_pathToFullMatrix)

    # convert the full matrix to a list of lists
    l_fullMatrix = lil_matrix(l_fullMatrix)

    # create preprocessor command
    l_preProcessorCommand = '#if ' + i_dimensionName  + ' == ' + str(i_dimensionSize) + '\n' 
    l_preProcessorCommand = l_preProcessorCommand + '#define ' + i_baseName.upper() + '_NUMBEROFNONZEROS ' + str(len(nonzero(l_fullMatrix)[0])) + '\n'
    l_preProcessorCommand = l_preProcessorCommand + '#endif \n\n'

    return l_preProcessorCommand

  ###
  ### Global functions
  ###

  # Returns a list of matrix directories contataining star, stiffness and flux matrices.
  #   Each matrix is specified by:
  #     fileNameOfGeneratedKernel    - base name of generated matrix kernel
  #     routineNameOfGeneratedKernel - name of the routine inside the generated kernel
  #     pathToMatrixMarketFile       - location of the matrix in MatrixMarkeitFormat
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
  def getSparseMatrices( self,
                         i_pathToMatrices,
                         i_numberOfQuantities = 9,
                         i_maximumDegreeOfBasisFunctions = 7 ):
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
    l_starMatrix = 'starMatrix_3D_maple.mtx'

    # assert existance
    assert( l_starMatrix in l_matrixFiles )

    # build complete path
    l_starMatrix = i_pathToMatrices+'/'+l_starMatrix
    
    # assert correct dimensions within the file
    l_matrixDimension = (open(l_starMatrix, 'r').readlines()[1]).split()
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
      l_kXi =   dict( matrixName = 'kXiDivM',   matrixMarketFileName='kXiDivM_3D_'   + str(l_degree) + '_maple.mtx' )
      l_kEta =  dict( matrixName = 'kEtaDivM',  matrixMarketFileName='kEtaDivM_3D_'  + str(l_degree) + '_maple.mtx' )
      l_kZeta = dict( matrixName = 'kZetaDivM', matrixMarketFileName='kZetaDivM_3D_' + str(l_degree) + '_maple.mtx' )
      
      l_kXiTransposed =   dict( matrixName = 'kXiDivMT',   matrixMarketFileName='kXiDivMT_3D_'   + str(l_degree) + '_maple.mtx' )
      l_kEtaTransposed =  dict( matrixName = 'kEtaDivMT',  matrixMarketFileName='kEtaDivMT_3D_'  + str(l_degree) + '_maple.mtx' )
      l_kZetaTransposed = dict( matrixName = 'kZetaDivMT', matrixMarketFileName='kZetaDivMT_3D_' + str(l_degree) + '_maple.mtx' )

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
      l_matrixDimension = (open(l_kXi['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      l_matrixDimension = (open(l_kEta['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      l_matrixDimension = (open(l_kEta['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      
      l_matrixDimension = (open(l_kXiTransposed['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      l_matrixDimension = (open(l_kEtaTransposed['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      l_matrixDimension = (open(l_kEtaTransposed['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
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
                                            matrixMarketFileName='fM'+str(l_localFace+1)+'DivM_3D_'   + str(l_degree) + '_maple.mtx' ) ]
        
        # assert existance
        assert( l_fluxMinus[-1]['matrixMarketFileName']   in l_matrixFiles )

        # generate complete path
        l_fluxMinus[-1]['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxMinus[l_localFace]['matrixMarketFileName']

        # assert correct dimensions
        l_matrixDimension = (open(l_fluxMinus[-1]['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
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
                                              matrixMarketFileName='fP'+l_multiIndex+'DivM_3D_'   + str(l_degree) + '_maple.mtx' ) ]

            # assert existance
            assert( l_fluxPlus[-1]['matrixMarketFileName']   in l_matrixFiles )

            # generate complete path
            l_fluxPlus[-1]['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxPlus[-1]['matrixMarketFileName']

            # assert correct dimensions
            l_matrixDimension = (open(l_fluxPlus[-1]['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
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

  # Converts a given full matrix to compressed sparse row and compressed sparse column format
  #
  # \param i_pathToFullMatrix path to the full matrix (format: matrix market).
  # \param i_cscFile open file object to csc file.
  # \param i_csrFile open file object to csr file.
  @staticmethod
  def convertFullToSparse(  i_pathToFullMatrix,
                            i_cscFile,
                            i_csrFile ):
    l_logger.log('converting to CSR and CSC: '+i_pathToFullMatrix, 2)

    # read the full matrix
    l_fullMatrix = mmread(i_pathToFullMatrix)

    # convert the full matrix to a list of lists
    l_fullMatrix = lil_matrix(l_fullMatrix)
    
    # generate a list containing the matrix as tuples: (row, column, value)
    l_sparseMatrix = l_fullMatrix.todok().keys()
    for l_i in range(len(l_sparseMatrix)):
      l_sparseMatrix[l_i] = l_sparseMatrix[l_i] + (l_fullMatrix.todok().values()[l_i],)
      
    # sort csc and csr in a dictionary
    l_sortedMatrices = {}
    l_sortedMatrices['csc'] = sorted( l_sparseMatrix, key = lambda e: (e[1], e[0] ) )
    l_sortedMatrices['csr'] = sorted( l_sparseMatrix, key = lambda e: (e[0], e[1] ) )
    
    # convert to coordinate format
    for l_format in ['csc', 'csr']:
      # data structure containing a list for rows, columns and data
      l_row = []
      l_col = []
      l_data = []
      
      for l_element in l_sortedMatrices[l_format]:
        l_row = l_row +   [l_element[0]]
        l_col = l_col +   [l_element[1]]
        l_data = l_data + [l_element[2]]
      
      # create coordinate format matrix
      l_sortedMatrices[l_format] = coo_matrix( (l_data, (l_row,l_col)), shape=l_fullMatrix.shape)

    # write csr and csc
    mmwrite(i_csrFile, l_sortedMatrices['csr'])
    mmwrite(i_cscFile, l_sortedMatrices['csc'])

  # Plot the sparsity pattern of a given full matrix.
  #
  # \param i_fullMatrix full matrix or path to the full matrix (format: matrix market).
  # \param i_baseName base name of the plot.
  # \param i_pathToOutputDirectory path to the output directory.
  # \param i_readMatrix false or true, whether the matrix needs be read in or not.
  def plotSparsityPattern( self,
                           i_fullMatrix,
                           i_baseName,
                           i_pathToOutputDirectory,
                           i_readMatrix = False):
    l_logger.log('plotting sparsity pattern of: '+i_fullMatrix, 2)
    
    if( i_readMatrix == True ):
      l_fullMatrix = mmread( i_fullMatrix )
    else:
      l_fullMatrix = i_fullMatrix

    # create black/white colormap
    l_colorMap = colors.ListedColormap(['white', 'black'])

    # set up ticks
    pyplot.xticks( arange(len(l_fullMatrix[0,:])) )
    pyplot.yticks( arange(len(l_fullMatrix[:,0])) )
    pyplot.tick_params(axis='both', which='major', labelsize=5)

    # draw image
    pyplot.imshow((l_fullMatrix!=0), interpolation='nearest', cmap=l_colorMap)

    # save figure to disk
    savefig(i_pathToOutputDirectory+'/'+i_baseName+'_sparsity_pattern.pdf')

  # Adds preprocessor code for the given matrix.
  # The decision about dimensions, base name and dimension is taken by considering the file name.
  # The code is stored in one of the class variables.
  #
  # @param i_pathToFullMatrix path to the full matrix (format: matrix market).
  # @param i_baseName of file (generated by Maple).
  def addMatrixToPreProcessorCode( self,
                                   i_pathToFullMatrix,
                                   i_baseName ):

    l_baseName = i_baseName.split('_')

    # create different statements for 2D and 3D
    if l_baseName[0] == 'starMatrix':
      # TODO: fixed to elastics in 3d
      l_dimensionSize = '9'
      l_localString = self.preProcessorCommands3D
      l_dimensionName = 'NUMBEROFVARIABLES'
    else:
      if l_baseName[1] == '2D':
        l_dimensionSize = (int(l_baseName[2]) + 1) * (int(l_baseName[2]) + 2) / 2
        l_localString = self.preProcessorCommands2D
      elif l_baseName[1] == '3D':
        l_dimensionSize = (int(l_baseName[2]) + 1) * (int(l_baseName[2]) + 2) * (int(l_baseName[2]) + 3) / 6
        l_localString = self.preProcessorCommands3D
      else:
        assert(False)
      l_dimensionName = 'NUMBEROFBASISFUNCTIONS'

    # override base name by true name
    l_baseName = l_baseName[0]

    for l_language in ['c++', 'fortran']:
      l_localString[l_language] =l_localString[l_language] + self.generatePreProcessorCode( i_pathToFullMatrix,
                                                                                            l_baseName,
                                                                                            i_dimensionSize = l_dimensionSize,
                                                                                            i_dimensionName = l_dimensionName )

  # Writes the stored preprocessor code to disk.
  #
  # \param i_pathToOutputDirectory path to the output directory.
  def writePreProcessorCode( self,
                            i_pathToOutputDirectory ):
    l_logger.log('writing preprocessor code')

    l_pathTo2DFortranFile = i_pathToOutputDirectory + '/matrixSizes2D.fpp'
    l_pathTo3DFortranFile = i_pathToOutputDirectory + '/matrixSizes3D.fpp'
    l_pathTo2DCppFile     = i_pathToOutputDirectory + '/matrixSizes2D.h'
    l_pathTo3DCppFile     = i_pathToOutputDirectory + '/matrixSizes3D.h'

    # create empty files
    open(l_pathTo2DFortranFile, 'w').close()
    open(l_pathTo3DFortranFile, 'w').close()
    open(l_pathTo2DCppFile, 'w').close()
    open(l_pathTo3DCppFile, 'w').close()

    # write license
    l_logger.writeFileHeader(l_pathTo2DFortranFile, '! ')
    l_logger.writeFileHeader(l_pathTo3DFortranFile, '! ')
    l_logger.writeFileHeader(l_pathTo2DCppFile, '// ')
    l_logger.writeFileHeader(l_pathTo3DCppFile, '// ')

    # open files
    l_matrixSizes2DFortran = open(l_pathTo2DFortranFile ,'a')
    l_matrixSizes3DFortran = open(l_pathTo3DFortranFile ,'a')
    l_matrixSizes2DCpp     = open(l_pathTo2DCppFile ,'a')
    l_matrixSizes3DCpp     = open(l_pathTo3DCppFile ,'a')

    l_descriptionFortran ='''!
! @section DESCRIPTION
! 
! Defines matrix sizes in the preprocessing phase.
!   Remark: This file was generated. It should not be modified by hand.
! 
'''
    l_descriptionCpp ='''//
// @section DESCRIPTION
// 
// Defines matrix sizes in the preprocessing phase.
//   Remark: This file was generated. It should not be modified by hand.
// 
'''

    # write files
    l_matrixSizes2DFortran.write( l_descriptionFortran+self.preProcessorCommands2D['fortran'] )
    l_matrixSizes3DFortran.write( l_descriptionFortran+self.preProcessorCommands3D['fortran'] )
    l_matrixSizes2DCpp.write(     l_descriptionCpp+self.preProcessorCommands2D['c++']     )
    l_matrixSizes3DCpp.write(     l_descriptionCpp+self.preProcessorCommands3D['c++']     )

    # close files
    l_matrixSizes2DFortran.close()
    l_matrixSizes3DFortran.close()
    l_matrixSizes2DCpp.close()
    l_matrixSizes3DCpp.close()

  def generateMatrixInitializationCode( self,
                                        i_pathToFullMatrix,
                                        i_mapleBaseName,
                                        i_functionName,
                                        i_sparseFormat,
                                        i_pathToOutputDirectory ):
    
    l_logger.log('generating matrix initialization code for '+i_pathToFullMatrix);

    # TODO: only csc is supported
    assert(i_sparseFormat == 'csc')    

    # read the full matrix
    l_fullMatrix = mmread(i_pathToFullMatrix)

    # convert the full matrix to a list of lists
    l_fullMatrix = lil_matrix(l_fullMatrix)

    # convert to csc
    l_fullMatrix = l_fullMatrix.tocsc()

    # generate function header
    l_sourceCode = '''! Initialization of flat matrices.
!   Remark: This file was generated. It should not be modified by hand.
  
!Copies non-zeros (regarding the original Maple pattern) of the full matrix to the flat matrix, Maple base name: ''' + i_mapleBaseName + '''
!   Remark: Zeros will be written to the flat matrix, if present at the given position.
!           Error checking is disabled with the preprocessor variable NDEBUG.
!           This file was generated. It should not be modified by hand.
!
!   The output order is column-major:
!<pre>
!          _                         _
!         | 0  0  0  0  0  0 16  0  0 |
!         | 0  0  0  0  0  0  0 19  0 |
!         | 0  0  0  0  0  0  0  0 22 |
!         | 0  0  0  0  0  0 17 20  0 |
!         | 0  0  0  0  0  0  0 21 23 | -> [1 2 3 4 5 .. 20 21 23 24]
!         | 0  0  0  0  0  0 18  0 24 |
!         | 1  4  7 10  0 14  0  0  0 |
!         | 2  5  8 11 12  0  0  0  0 |
!         |_3  6  9  0 13 15  0  0  0_|
!</pre>
! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!
subroutine '''
    l_sourceCode = l_sourceCode + i_functionName + '''( i_fullMatrix, o_flatMatrix )
implicit none

! function parameters
real, intent(in),  dimension( :,: ) :: i_fullMatrix !< full matrix stored in a two-dimensional array.
real, intent(out), dimension( : )   :: o_flatMatrix !< flat matrix will be set to: non-zeros of the full matrix (col-major).

! local variables
#ifndef NDEBUG
integer :: l_row, l_col   !< loop variables.
#endif

! function
'''
    # assert a proper size of the flat matrix
    l_sourceCode = l_sourceCode +'''#ifndef NDEBUG
! assert a proper size of the flat matrix
if( ('''+str(l_fullMatrix.tocsc().getnnz())+''' .ne. ubound(o_flatMatrix,1)) .or. (1 .ne. lbound(o_flatMatrix,1)) ) then
  write(STDERR, *) '''+'\''+i_functionName+''', o_flatMatrix bounds are not matching.'
  stop
endif
#endif

'''

    # assert that all nonzeros are covered
    l_sourceCode = l_sourceCode + '''#ifndef NDEBUG
! assert that all nonzeros are covered
do l_col = lbound(i_fullMatrix, 2), ubound(i_fullMatrix, 2)
  do l_row = lbound(i_fullMatrix, 1), ubound(i_fullMatrix, 1)
    if ( (abs(i_fullMatrix(l_row, l_col)) .ge. ZEROTOLERANCE ) .and. .not. ( '''
    # get non-zero indices sorted columnwise
    l_cscIndizes =  sorted( l_fullMatrix.todok().keys(), key = lambda e: (e[1], e[0]) )

    # iterate over non zeros and generate an if-statement
    for l_nonZeroEntry in l_cscIndizes:
      l_sourceCode = l_sourceCode + ' ( l_row .eq. ' + str(l_nonZeroEntry[0]+1)\
                                  + ' .and. l_col .eq. ' + str(l_nonZeroEntry[1]+1) + ' )'
      if l_nonZeroEntry != l_cscIndizes[-1]:
        l_sourceCode = l_sourceCode + ' .or. '
    l_sourceCode = l_sourceCode+''') ) then
       write(STDERR, *) '''+'\''+i_functionName+''', o_flatMatrix does not cover all nonzeros.', l_row, l_col, i_fullMatrix(l_row, l_col)
       stop
    endif
  end do
end do
#endif

'''

    # initialization of the flat matrix
    l_flatMatrixIndex = 1
    for l_nonZeroEntry in l_cscIndizes:
      l_sourceCode = l_sourceCode + 'o_flatMatrix(' + str(l_flatMatrixIndex) +\
                                    ') = i_fullMatrix(' + str(l_nonZeroEntry[0]+1) + ','\
                                    + str(l_nonZeroEntry[1]+1) + ')\n'
      l_flatMatrixIndex = l_flatMatrixIndex + 1


    # close function
    l_sourceCode = l_sourceCode + '\nend subroutine '+i_functionName

    # write code to disk
    l_pathToOutputFile = i_pathToOutputDirectory+'/'+i_functionName+'.f90_include'
    
    # create empty files
    open(l_pathToOutputFile , 'w').close()

    # write file header
    l_logger.writeFileHeader(l_pathToOutputFile, '! ')
    

    l_file = open(l_pathToOutputFile ,'a')
    l_file.write(l_sourceCode)
    l_file.close()

  # reads a matrix in market format.
  #
  # \input i_pathToMatrix path to the matrix.
  # \return dictionary with the matrix structure in coordinate format
  def readMatrixMarket( self,
                        i_pathToMatrix,
                      ):
    l_matrixFile = open(i_pathToMatrix)

    # jump over header
    if(    l_matrixFile.readline() != "%%MatrixMarket matrix array real general\n" ):
      print( 'aborting, wrong matrix market format' )
      exit(1)

    l_dimensions = l_matrixFile.readline().split()
    l_numberOfRows    = int( l_dimensions[0] )
    l_numberOfColumns = int( l_dimensions[1] )

    l_row = 1; l_column = 1;

    l_values = []
    l_rows   = []
    l_columns = []

    while( True ):
      l_value = l_matrixFile.readline().rstrip("\n");

      if not l_value:
        assert( l_row    == 1                   );
        assert( l_column == l_numberOfColumns+1 );
        break
      
      if( l_value != "0." ):
        l_values = l_values + [l_value]
        l_rows = l_rows       + [l_row]
        l_columns = l_columns + [l_column]

      l_row = l_row + 1
      if( l_row > l_numberOfRows ):
        l_row = 1
        l_column = l_column + 1

    # return everythin
    return { '#rows':    l_numberOfRows,
             '#columns': l_numberOfColumns,
             '#nnz':     len( l_values ),
             'rows':     l_rows,
             'columns':  l_columns,
             'values':   l_values }

  # Converts a given set of dense matrices to a single xml-file with sparse index storage.
  #
  # \param i_pathToMatrices path to the dense matrix (format: matrix market).
  def convertToXml(  self,
                     i_pathToMatrices,
                     i_pathToOutputDirectory ):
    l_logger.log('converting matrices in folder \''+i_pathToMatrices )

    # #(basis functions) we write XML files for
    l_numberOfBasisFunctionsList = [4, 10, 20, 35, 56, 84];

    # get sparse matrices
    l_matrices = self.getSparseMatrices( i_pathToMatrices                = i_pathToMatrices,
                                         i_numberOfQuantities            = 9,
                                         i_maximumDegreeOfBasisFunctions = 8 )

    # iterate over basis functions
    #   Remark: Each matrix is a subset of the matrix for the next degree.
    #           By generating different XML-file we are able to decide per degree if
    #           the same matrix should be handled sparse or dense.
    for l_numberOfBasisFunctions in l_numberOfBasisFunctionsList:
      # generate file name
      l_pathToOutputFile = i_pathToOutputDirectory+'/matrices_'+str(l_numberOfBasisFunctions)+'.xml'

      l_logger.log('writing: '+ l_pathToOutputFile, 2);

      # open xml file
      l_xmlFile = XMLWriter(l_pathToOutputFile)

      # root element
      l_xmlFile.start('matrices')

      #
      # add global matrices for this degree
      #
      l_xmlFile.start("global")

      for l_matrix in l_matrices:
        # read the matrix structure
        l_matrixStructure = self.readMatrixMarket( i_pathToMatrix = l_matrix['pathToMatrixMarketFile'] )

        # consider only global matrices (stiffness/flux)
        if l_matrixStructure['#rows'] != l_numberOfBasisFunctions:
          continue;

        # write xml for this matrix
        if 'fP' in l_matrix['name'] or 'fM' in l_matrix['name']:
          l_matrixType = "flux"
        elif 'kXi' in l_matrix['name'] or 'kEta' in l_matrix['name'] or 'kZeta' in l_matrix['name']:
          l_matrixType = "stiffness"
        else:
          assert(False)

        # assert dimensions match
        assert( len( l_matrixStructure['columns'] ) == len( l_matrixStructure['rows']   ) )
        assert( len( l_matrixStructure['columns'] ) == len( l_matrixStructure['values'] ) )

        l_xmlFile.start( l_matrixType,
          # add matrix meta information
                         name    = l_matrix['name'],
                         id      = str(l_matrix['id']),
                         rows    = str(l_matrixStructure['#rows']),
                         columns = str(l_matrixStructure['#columns']) )
      
        for l_entry in xrange( len( l_matrixStructure['rows'] ) ):
          # add this matrix entry
          l_xmlFile.element( "entry",
                             row    = str(l_matrixStructure['rows'][l_entry]),
                             column = str(l_matrixStructure['columns'][l_entry]),
                             value  = l_matrixStructure['values'][l_entry] )

        #print l_matrixEntries
        l_xmlFile.end(l_matrixType)
    
      l_xmlFile.end('global')

      #
      # add local matrices for this degree
      # TODO: Hardcoded #variables for elastics
      #
      l_xmlFile.start('local')
      
      # Add flux solver
      # \f$  N_{k,i} A_k^+ N_{k,i}^{-1} \f$ and \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      l_xmlFile.element( 'fluxSolver',
                         id='52',
                         rows='9',
                         columns='9' )

      # Add star matrices
      for l_matrix in l_matrices:
        if( l_matrix['name'] in ['volumeStarMatrix'] ):
          # read the matrix structure
          l_matrixStructure = self.readMatrixMarket( i_pathToMatrix = l_matrix['pathToMatrixMarketFile'] )
          
          # add matrix meta information
          l_xmlFile.start( 'starMatrix',
                           id      = str(l_matrix['id']),
                           rows    = str(l_matrixStructure['#rows']),
                           columns = str(l_matrixStructure['#columns']) )

      
          for l_entry in xrange( len(l_matrixStructure['rows']) ):
            # add this matrix entry
            l_xmlFile.element( "entry",
                               row    = str(l_matrixStructure['rows'][l_entry]),
                               column = str(l_matrixStructure['columns'][l_entry])
                             )
          l_xmlFile.end('starMatrix')
          
          # add only a single star matrix
          break
      l_xmlFile.end('local')

      l_xmlFile.end('matrices')
