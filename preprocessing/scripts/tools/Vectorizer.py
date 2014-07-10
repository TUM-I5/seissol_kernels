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
# Performs a vectorization of matrix rows or columns.
#   Please note: A simple brute force is used in the scripts.
#                There's much space left for optimization.
#

from scipy.io.mmio import *
from scipy.sparse import *
from numpy import *

from datetime import datetime

from matplotlib import pyplot, colors
from pylab import *

import tools.MatrixConverter as l_matrixConverter
import tools.Logger as l_logger

###
### Local functions
###

# Place a single interval to all the possible places.
#   This function works recursive, therefore at each level new intervals will be placed until the threshold is exceeded.
#
# \param i_boolVector vector, which should be covered.
# \param i_intervalVector vector, which is currently covered by the intervals.
# \param i_minimumStartIndex minimum position where the end of the interval at the current level should be placed.
# \param i_numberOfIntervals number of intervals which still need a position.
# \param i_intervalSize size of each interval.
# \param o_result will be set to the result, once a solution is found.
#
def placeInterval( i_boolVector, i_intervalVector , i_minimumStartIndex, i_numberOfIntervals, i_intervalSize, o_result):

  # compute remaining intervals which need to be added
  l_numberOfIntervals = i_numberOfIntervals -1

  l_currentStartIndex = i_minimumStartIndex

  l_foundSolution = False

  while( l_currentStartIndex < len(i_boolVector) ):
    # create a temporary copy of the interval vector
    l_tempVector = copy(i_intervalVector)

    l_tempVector[max(0,l_currentStartIndex):l_currentStartIndex+i_intervalSize] = True

    # add new intervals if necessary or check for optimality
    if( l_numberOfIntervals != 0 ):
      # the sub intervals found an optimal solution
      if( placeInterval( i_boolVector, l_tempVector, l_currentStartIndex+1, l_numberOfIntervals, i_intervalSize, o_result) ):
        l_foundSolution = True
    # all intervals added check for optimality
    else:
      if( (l_tempVector + i_boolVector).all() == True ):
       o_result[max(0,l_currentStartIndex):l_currentStartIndex+i_intervalSize] = i_numberOfIntervals
       l_foundSolution = True
    
    if(l_foundSolution == True):
       o_result[max(0,l_currentStartIndex):l_currentStartIndex+i_intervalSize] = i_numberOfIntervals
       return True

    l_currentStartIndex += 1
  
  # nothing found (shift parent interval or try with more intervals)
  return False

# Compute minimum number of intervals with the given size which cover all non-zero elements of a vector.
#
# \param i_vector vector which should be covered.
# \param size of each interval.
def computeCoveringIntervals( i_vector, i_intervalSize):
  if (count_nonzero(i_vector) == 0):
    return zeros(len(i_vector))

  # convert the given numeric vector to a vector of booleans whether there is a zero entry or not
  l_boolVector = i_vector == 0

  # start with a minimum number of intervals
  l_numberOfIntervals = ceil( count_nonzero(i_vector) / float(i_intervalSize) )

  # where to store the results
  l_result = zeros(len(i_vector))
  
  # iterate of the number of intervals
  while( True ):
    l_logger.log( 'number of intervals: '+str(l_numberOfIntervals), 4 )

    # exit if we found the solution
    if ( placeInterval( l_boolVector, zeros(len(l_boolVector)).astype(bool), 0, l_numberOfIntervals, i_intervalSize, l_result) ):
      return l_result
    l_numberOfIntervals+=1

# Write the vectorized matrices to sparse format (CSC and CSR). Additional zero entries are added where vector operations are used.
# \param i_denseMatrix original dense matrix (contains all the initial zero entries).
# \param i_vectorizedMatrix vectorized matrix, where the original values are replaced by the correspond number of the vector operation. Therefore original zeros might be non-zeros due to the vectorization procedure
# \param i_baseName base name of the sparse output.
# \param i_pathToOutputDirectory path to the output directory.
def writeSparseVectorizedMatrices( i_denseMatrix,
                                   i_vectorizedMatrix,
                                   i_baseName,
                                   i_pathToOutputDirectory ):
  l_denseMatrix = copy(i_denseMatrix)

  # assert that we have matrices of the same shape
  assert(len(l_denseMatrix[:,0]) == len(i_vectorizedMatrix[:,0]))
  assert(len(l_denseMatrix[0,:]) == len(i_vectorizedMatrix[0,:]))

  # set overlapping (dense/vector) elements to nan
  for l_i in range(len(i_denseMatrix[:,0])):
    for l_j in range(len(i_denseMatrix[0,:])):
      if( l_denseMatrix[l_i][l_j] == 0. and i_vectorizedMatrix[l_i][l_j] != 0. ):
        l_denseMatrix[l_i][l_j] = np.nan

  # convert the dense matrix to a list of lists
  l_denseMatrix = lil_matrix(l_denseMatrix)

  # define 'base'-paths
  l_pathToCsrFile = i_pathToOutputDirectory+'/'+i_baseName+'_csr'
  l_pathToCscFile = i_pathToOutputDirectory+'/'+i_baseName+'_csc'

  # write csr and csc
  mmwrite(l_pathToCsrFile, l_denseMatrix.tocsr())
  mmwrite(l_pathToCscFile, l_denseMatrix.tocsc())

  # make the paths complete
  l_pathToCsrFile += '.mtx'
  l_pathToCscFile += '.mtx'

  # replace the nan's by zeros again
  for l_sparseFile in [l_pathToCsrFile, l_pathToCscFile]:
    # read the current lines
    l_lines = open(l_sparseFile, 'r').readlines()
    
    # replace the file with a new empty one
    l_file = open(l_sparseFile, 'w')
    for l_line in l_lines:
      # print all the modified strings
      l_file.write( l_line.replace('nan', '0.0') )
    l_file.close()

# Plot the vectorization of a given vectorized matrix.
#
# \param i_vectorMatriox vectorized matrix.
# \param i_baseName base name of the plot.
# \param i_pathToOutputDirectory path to the output directory.
def plotVectorization( i_vectorMatrix,
                       i_baseName,
                       i_pathToOutputDirectory ):
  # convert list to matrix
  l_vectorMatrix = matrix(i_vectorMatrix)

  # replace zeros by nan's and set corr. bad color
  for l_i in xrange(len(l_vectorMatrix[:])):
    for l_j in xrange(len(l_vectorMatrix[:])):
      if(l_vectorMatrix[l_i, l_j] == 0.):
        l_vectorMatrix[l_i, l_j] = np.nan
  l_colorMap = matplotlib.cm.get_cmap('spectral')
  l_colorMap.set_bad('w')

  # set up the axis
  #pyplot.xticks( arange(len(l_denseMatrix[0,:])) )
  #pyplot.yticks( arange(len(l_denseMatrix[:,0])) )
  pyplot.tick_params(axis='both', which='major', labelsize=5)
  
  # draw the vectorized result
  pyplot.imshow(transpose(matrix(l_vectorMatrix)), interpolation='nearest', cmap=l_colorMap)

  # save figure to disk
  savefig(i_pathToOutputDirectory+'/'+i_baseName+'_vectorized.pdf')

# Performs a column-wise vectorization of the given matrix.
#
# \param i_pathToDenseMatrix path to the matrix, which will be vectorized (format: matrix market).
# \param i_baseName base name of the vectorized output.
# \param i_pathToOutputDirectory path to the output directory.
# \param i_drawPlots If true: draw plots of the input matrices and vectorized results.

def vectorizeMatrix( i_pathToDenseMatrix,
                     i_baseName,
                     i_pathToOutputDirectory,
                     i_drawPlots ):
  l_logger.log('vectorizing: '+ i_pathToDenseMatrix, 2)

  # read matrix
  l_denseMatrix = mmread( i_pathToDenseMatrix )

  # create list for the vectorized result.
  l_vectorMatrix = list()

  # do the vectorization
  for l_i in xrange(len(l_denseMatrix[:,0])):
    l_logger.log('current column: ' + str(l_i), 3)
    l_vectorMatrix.append( computeCoveringIntervals(l_denseMatrix[:,l_i], 4) )

  # write vectorized matrix to disk
  writeSparseVectorizedMatrices( l_denseMatrix, transpose(l_vectorMatrix), i_baseName+'_vec', i_pathToOutputDirectory)

  # draw plots
  if i_drawPlots:
    l_matrixConverter.plotSparsityPattern( l_denseMatrix,
                                           i_baseName,
                                           i_pathToOutputDirectory )

    plotVectorization( l_vectorMatrix,
                       i_baseName,
                       i_pathToOutputDirectory )
