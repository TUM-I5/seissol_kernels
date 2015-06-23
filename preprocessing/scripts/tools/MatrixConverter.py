#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
from numpy import arange, nonzero, sort, float64
from lxml import etree

import Logger as l_logger

import os

class MatrixConverter():  
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
                         i_numberOfViscoelasticQuantities = 6,
                         i_maximumDegreeOfBasisFunctions = 7 ):
    l_logger.log('getting matrices', 2)

    # list which holds the different matrix structures
    l_sparseMatrices = []

    # get and sort the matrix files
    l_matrixFiles = os.listdir(i_pathToMatrices)
    l_matrixFiles.sort()

    # convert input parameters to str
    l_numberOfQuantities = str(i_numberOfQuantities)
    l_numberOfViscoelasticQuantities = str(i_numberOfViscoelasticQuantities)

    ###
    # star matrices
    ###  

    # name of the star matrices
    l_starMatrix = dict( name = 'starMatrix', id = 59, matrixMarketFileName='starMatrix_3D_maple.mtx' )
    assert( l_starMatrix['matrixMarketFileName'] in l_matrixFiles )
    l_starMatrix['pathToMatrixMarketFile'] = i_pathToMatrices + '/' + l_starMatrix['matrixMarketFileName']
    l_matrixDimension = (open(l_starMatrix['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
    assert( l_matrixDimension[0] == l_numberOfQuantities )
    assert( l_matrixDimension[1] == l_numberOfQuantities )
    l_sparseMatrices.append(l_starMatrix)

    # viscoelastic star matrix
    l_viscoelasticStarMatrix = dict( name = 'viscoelasticStarMatrix', id = 62, matrixMarketFileName='viscoelasticStarMatrix_3D_maple.mtx' )
    assert( l_viscoelasticStarMatrix['matrixMarketFileName'] in l_matrixFiles )
    l_viscoelasticStarMatrix['pathToMatrixMarketFile'] = i_pathToMatrices + '/' + l_viscoelasticStarMatrix['matrixMarketFileName']
    l_matrixDimension = (open(l_viscoelasticStarMatrix['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
    assert( l_matrixDimension[0] == l_numberOfQuantities)
    assert( l_matrixDimension[1] == l_numberOfViscoelasticQuantities )
    l_sparseMatrices.append(l_viscoelasticStarMatrix)
    
    # viscoelastic source matrix
    l_et = dict( name = 'viscoelasticSourceMatrix', id = 60, matrixMarketFileName='eT_3D_maple.mtx' )
    assert( l_et['matrixMarketFileName'] in l_matrixFiles )
    l_et['pathToMatrixMarketFile'] = i_pathToMatrices + '/' + l_et['matrixMarketFileName']
    l_matrixDimension = (open(l_et['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
    assert( l_matrixDimension[0] == l_numberOfViscoelasticQuantities )
    assert( l_matrixDimension[1] == l_numberOfQuantities )
    l_sparseMatrices.append(l_et)

    ###
    # stiffness, flux, and mass matrices
    ##
    for l_degree in range(1,i_maximumDegreeOfBasisFunctions):
      # each matrix is a dictionary containing information about it
      l_matrix = dict()

      l_numberOfBasisFunctions = str((l_degree+1)*(l_degree+2)*(l_degree+3)/6)

      l_logger.log( 'adding star, stiffness, flux, and mass matrices to dictionaries for: '+str(l_degree+1)+' (order of basis), '+l_numberOfBasisFunctions+' (#basis functions), '+l_numberOfQuantities+' (#quantities)', 3)

      # name and filename of the stiffness matrices, TODO: avoid copy and paste code..
      l_kXi =   dict( name = 'kXiDivM',   id = 53, matrixMarketFileName='kXiDivM_3D_'   + str(l_degree) + '_maple.mtx' )
      l_kEta =  dict( name = 'kEtaDivM',  id = 54, matrixMarketFileName='kEtaDivM_3D_'  + str(l_degree) + '_maple.mtx' )
      l_kZeta = dict( name = 'kZetaDivM', id = 55, matrixMarketFileName='kZetaDivM_3D_' + str(l_degree) + '_maple.mtx' )
      
      l_kXiTransposed =   dict( name = 'kXiDivMT',   id = 56, matrixMarketFileName='kXiDivMT_3D_'   + str(l_degree) + '_maple.mtx' )
      l_kEtaTransposed =  dict( name = 'kEtaDivMT',  id = 57, matrixMarketFileName='kEtaDivMT_3D_'  + str(l_degree) + '_maple.mtx' )
      l_kZetaTransposed = dict( name = 'kZetaDivMT', id = 58, matrixMarketFileName='kZetaDivMT_3D_' + str(l_degree) + '_maple.mtx' )

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
      
      l_sparseMatrices.append(l_kXi)
      l_sparseMatrices.append(l_kEta)
      l_sparseMatrices.append(l_kZeta)
      l_sparseMatrices.append(l_kXiTransposed)
      l_sparseMatrices.append(l_kEtaTransposed)
      l_sparseMatrices.append(l_kZetaTransposed)

      ###
      # Flux matrices
      ##

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
        l_fluxMinus = dict( name = 'fM'+str(l_localFace+1),
                            id = l_localFace,
                            matrixMarketFileName='fM'+str(l_localFace+1)+'DivM_3D_'   + str(l_degree) + '_maple.mtx' )
        
        # assert existance
        assert( l_fluxMinus['matrixMarketFileName']   in l_matrixFiles )

        # generate complete path
        l_fluxMinus['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxMinus['matrixMarketFileName']

        # assert correct dimensions
        l_matrixDimension = (open(l_fluxMinus['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
        assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
        assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
        
        l_sparseMatrices.append(l_fluxMinus)

      for l_localFace in range(0,4):
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
            l_fluxPlus = dict( name = 'fP'+l_multiIndex,
                               id   = l_matrixId,
                               matrixMarketFileName='fP'+l_multiIndex+'DivM_3D_'   + str(l_degree) + '_maple.mtx' )
            # assert existance
            assert( l_fluxPlus['matrixMarketFileName']   in l_matrixFiles )

            # generate complete path
            l_fluxPlus['pathToMatrixMarketFile'] = i_pathToMatrices+'/'+l_fluxPlus['matrixMarketFileName']

            # assert correct dimensions
            l_matrixDimension = (open(l_fluxPlus['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
            assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
            assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
            
            l_sparseMatrices.append(l_fluxPlus)
      ###
      # Mass matrices
      ##
      
      l_m = dict( name = 'm', id = 63, matrixMarketFileName='m_3D_' + str(l_degree) + '_maple.mtx')
      assert( l_m['matrixMarketFileName']   in l_matrixFiles )
      l_m['pathToMatrixMarketFile'] = i_pathToMatrices + '/' + l_m['matrixMarketFileName']
      l_matrixDimension = (open(l_m['pathToMatrixMarketFile'], 'r').readlines()[1]).split()
      assert( l_matrixDimension[0] == l_numberOfBasisFunctions )
      assert( l_matrixDimension[1] == l_numberOfBasisFunctions )
      l_sparseMatrices.append(l_m)
      

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
    
  # Plot the sparsity pattern of all matrices in a directory.
  #
  # \param i_pathToMatrices path that contains the matrices.
  # \param i_pathToOutputDirectory path to the output directory.
  def plotSparsityPatterns( self,
                            i_pathToMatrices,
                            i_pathToOutputDirectory ):
    l_matrixFiles = os.listdir(i_pathToMatrices)
    for l_file in l_matrixFiles:
      if l_file.endswith('_maple.mtx'):
        l_baseName = l_file.replace('_maple.mtx','')
        self.plotSparsityPattern( i_fullMatrix = i_pathToMatrices + '/' + l_file,
                                               i_baseName = l_baseName,
                                               i_pathToOutputDirectory = i_pathToOutputDirectory,
                                               i_readMatrix = True )

  def writeSparsityPattern(self, i_matrix, o_node):
    l_matrixStructure = self.readMatrixMarket( i_pathToMatrix = i_matrix['pathToMatrixMarketFile'] )
    
    # add matrix meta information
    l_star_attributes = {"id" : str(i_matrix['id']), "rows" : str(l_matrixStructure['#rows']), "columns" : str(l_matrixStructure['#columns'])}
    l_star = etree.SubElement(o_node, i_matrix['name'], l_star_attributes)

    for l_entry in l_matrixStructure['matrix']:
      # element attributes
      l_star_entry_attributes = {"row" : str(l_entry[0]), "column" : str(l_entry[1]) }
      # add node to XML
      l_star_entry = etree.SubElement(l_star, "entry", l_star_entry_attributes)


  # Reads a matrix in market format. We do not use mmread here, as we do
  # not want any alteration of the values due to rounding issues.
  #
  # \input i_pathToMatrix path to the matrix.
  # \return dictionary with the matrix structure in coordinate format
  def readMatrixMarket( self,
                        i_pathToMatrix ):
    matrixFile = open(i_pathToMatrix)
    
    if (not matrixFile.readline().startswith('%%MatrixMarket matrix array real general')):
      print('Wrong matrix market format.')
      exit(1)
      
    dimensions = matrixFile.readline().split()
    numberOfRows = int(dimensions[0])
    numberOfColumns = int(dimensions[1])
    
    matrix = numberOfRows * numberOfColumns * [(0, 0, '')]  # pre-allocate matrix
    entry = 0
    for line in matrixFile:
      # format: row, column, value
      matrix[entry] = (entry % numberOfRows + 1, entry / numberOfRows + 1, line.strip());
      entry = entry + 1
      
    sparseMatrix = filter(lambda x: float64(x[2]) != 0.0, matrix)
    
    return { '#rows':     numberOfRows,
             '#columns':  numberOfColumns,
             'matrix':    sparseMatrix
    }

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

      # root element
      l_root = etree.Element("matrices")

      #
      # add global matrices for this degree
      #
      l_global = etree.SubElement(l_root, "global")

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
        elif 'm' in l_matrix['name']:
          l_matrixType = "inverseMass"
          l_matrixStructure['matrix'] = [(x[0], x[1], '%.20f' % (1.0 / float64(x[2]))) for x in l_matrixStructure['matrix']]
        else:
          assert(False)

        # matrix attributes
        l_global_matrix_attributes = {"name": l_matrix['name'], "id": str(l_matrix['id']), "rows" : str(l_matrixStructure['#rows']), "columns" : str(l_matrixStructure['#columns'])}
        # add node to XML
        l_global_matrix = etree.SubElement(l_global, l_matrixType, l_global_matrix_attributes) 
      
        for l_entry in l_matrixStructure['matrix']:
          # element attributes
          l_global_entry_attributes = {"row": str(l_entry[0]), "column": str(l_entry[1]), "value" : l_entry[2]}
          # add node to XML
          l_global_entry = etree.SubElement(l_global_matrix, "entry", l_global_entry_attributes)
         
      #
      # add local matrices for this degree
      # TODO: Hardcoded #variables for elastics
      #
      l_local = etree.SubElement(l_root, "local")
      
      # Add flux solvers
      # \f$  N_{k,i} A_k^+ N_{k,i}^{-1} \f$ and \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      l_flux_attributes = {"id" : "52", "rows" : "9", "columns" : "9"}
      l_flux = etree.SubElement(l_local, "fluxSolver", l_flux_attributes)
      # viscoelastic
      l_flux_attributes = {"id" : "61", "rows" : "9", "columns" : "6"}
      l_flux = etree.SubElement(l_local, "viscoelasticFluxSolver", l_flux_attributes)

      # Add star matrices
      for l_matrix in l_matrices:
        if( l_matrix['name'] in ['starMatrix', 'viscoelasticStarMatrix', 'viscoelasticSourceMatrix'] ):
          self.writeSparsityPattern(l_matrix, l_local)

      #write XML file
      l_xml_tree = etree.ElementTree(l_root)
      l_xml_tree.write(l_pathToOutputFile, pretty_print=True, encoding='utf-8')

