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
# Converts different matrices to CSC and CSR.
#

from scipy.io.mmio import *
from scipy.sparse import *

from matplotlib import pyplot, colors
from pylab import savefig
from numpy import arange, nonzero, sort

import tools.Logger as l_logger
import tools.SeisSolGen as l_seisSolGen

from elementtree.SimpleXMLWriter import XMLWriter

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

  # Converts a given full matrix to compressed sparse row and compressed sparse column format
  #   An "_csr" and "_csc" will be appended to the base name.
  #
  # \param i_pathToFullMatrix path to the full matrix (format: matrix market).
  # \param i_baseName base name of the sparse output.
  # \param i_pathToOutputDirectory path to the output directory.
  def convertFullToSparse(  self,
                            i_pathToFullMatrix,
                            i_baseName,
                            i_pathToOutputDirectory ):
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
    mmwrite(i_pathToOutputDirectory+'/'+i_baseName+'_csr', l_sortedMatrices['csr'])
    mmwrite(i_pathToOutputDirectory+'/'+i_baseName+'_csc', l_sortedMatrices['csc'])

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
    # read the full matrix
    l_matrixEntries = mmread( i_pathToMatrix )

    # convert matrix to coordinate format
    l_matrixEntries = coo_matrix(l_matrixEntries)

    # get #rows and #columns
    l_numberOfRows = l_matrixEntries.shape[0]
    l_numberOfColumns = l_matrixEntries.shape[1]

    # get #nnz
    l_numberOfNonZeros = l_matrixEntries.nnz

    # get row, columns and values
    l_rows    = l_matrixEntries.row
    l_columns = l_matrixEntries.col
    l_values  = l_matrixEntries.data

    # return everythin
    return { '#rows':    l_numberOfRows,
             '#columns': l_numberOfColumns,
             '#nnz':     l_numberOfNonZeros,
             'rows':     l_rows,
             'columns':  l_columns,
             'values':   l_values }

  # Converts a given set of dense matrices to a single xml-file with sparse index storage.
  #
  # \param i_pathToMatrices path to the dense matrix (format: matrix market).
  def convertToXml(  self,
                     i_pathToMatrices,
                     i_pathToOutputDirectory,
                     i_sparseDenseSwitch = 0.2 ):
    l_logger.log('converting matrices in folder \''+i_pathToMatrices )

    # #(basis functions) we write XML files for
    l_numberOfBasisFunctionsList = [4, 10, 20, 35, 56];

    # get sparse matrices
    l_matrices = l_seisSolGen.getSparseMatrices( i_pathToMatrices = i_pathToMatrices,
                                                 i_numberOfQuantities = 9,
                                                 i_maximumDegreeOfBasisFunctions = 7 )

    # iterate over basis functions
    #   Remark: Each matrix is a subset of the matrix for the next degree.
    #           By generating different XML-file we are able to decide per degree if
    #           the same matrix should be handled sparse or dense.
    for l_numberOfBasisFunctions in l_numberOfBasisFunctionsList:
      # generate file name
      l_pathToOutputFile = i_pathToOutputDirectory+'/matrices_'+("%0.2f" % i_sparseDenseSwitch)+'_'+str(l_numberOfBasisFunctions)+'.xml'

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
        assert( l_matrixStructure['rows'].size      == l_matrixStructure['rows'].size )
        assert( l_matrixStructure['columns'].size   == l_matrixStructure['values'].size  )

        # TODO: add proper sparse/dense switch here
        if( l_matrixStructure['#nnz'] / float( l_matrixStructure['#rows'] * l_matrixStructure['#columns'] ) > i_sparseDenseSwitch ):
          l_sparse = 'false'
        else:
          l_sparse = 'true'

        l_xmlFile.start( l_matrixType,
          # add matrix meta information
                         name    = l_matrix['name'],
                         id      = str(l_matrix['id']),
                         sparse  = l_sparse,
                         rows    = str(l_matrixStructure['#rows']),
                         columns = str(l_matrixStructure['#columns']) )

      
        for l_entry in xrange( l_matrixStructure['rows'].size ):
          # add this matrix entry
          l_xmlFile.element( "entry",
                             row    = str(l_matrixStructure['rows'][l_entry]+1),
                             column = str(l_matrixStructure['columns'][l_entry]+1),
                             value  = repr(l_matrixStructure['values'][l_entry]) )

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
                         sparse='false',
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
                           sparse  = 'true', # TODO: Always sparse, dense switch?
                           rows    = str(l_matrixStructure['#rows']),
                           columns = str(l_matrixStructure['#columns']) )

      
          for l_entry in xrange( l_matrixStructure['rows'].size ):
            # add this matrix entry
            l_xmlFile.element( "entry",
                               row    = str(l_matrixStructure['rows'][l_entry]+1),
                               column = str(l_matrixStructure['columns'][l_entry]+1)
                             )
          l_xmlFile.end('starMatrix')
          
          # add only a single star matrix
          break
      l_xmlFile.end('local')

      l_xmlFile.end('matrices')
