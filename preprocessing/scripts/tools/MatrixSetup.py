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
# Sets up the matrix configurations.
#
import logging
import xmltodict
import pprint

import Configuration

class MatrixSetup:

  m_configuration = 0

  # Constuctor
  #
  # @param i_configuration Configuration
  def __init__( self,
                i_configuration ):
    logging.debug( "Constructed a new MatrixSetup()" )

    self.m_configuration = i_configuration

  # Get the number of basis function for the given convergence order
  #
  # @param i_order convegence order.
  def getNumberOfBasisFunctions( self,
                                 i_order ):
    return i_order*(i_order+1)*(i_order+2)/6

  # Get the number of basis function for the given convergence order and floating point precision
  # aligned to the given alignment
  #
  # @param i_order convergence order.
  # @param i_alignment target alignment.
  # @param i_precision floating point precision. 
  def getAlignedNumberOfBasisFunctions( self,
                                        i_order,
                                        i_alignment,
                                        i_precision ):

    # compute the #remaing non-zeros for the recursive scheme
    l_numberOfNonZeroBasisFunctions = self.getNumberOfBasisFunctions( i_order )

    l_nonZeroBits = l_numberOfNonZeroBasisFunctions * i_precision

    # align the per-derivative non-zeros to memory
    l_alignedBits =  l_nonZeroBits +\
                    (i_alignment-(l_nonZeroBits % i_alignment))%i_alignment

    l_alignedBasisFunctions = l_alignedBits / i_precision

    return l_alignedBasisFunctions

  # Gets the leading dimension for the recursive computation of the time derivatives.
  # Each step of the recursion has an associated order.
  #
  # Example for O = 4, 64 bit alignment and double precision
  #
  # || 20 basis functions | 4 entries alignment | 10 basis functions | 6 entries alignment || 4 basis functions | 4 entries alignment || 1 basis function | 7 entries alignment ||  
  #
  # This leads in total to a leading dimension of 20+4+10+6+4+4+1+7=56 for the associated order 4, 10+6+4+4+1+7=32 for the associated order 3 and so on. 
  # Note that no entry is stored for the 3rd derivative,
  #
  # @param i_associatedOrder associated order of the time derivative
  # @param i_alignment memory alignment in bits.
  # @param i_precision floating point precision.
  #
  def getLeadingDimensionOfTimeDerivatives( self,
                                            i_associatedOrder,
                                            i_alignment,
                                            i_precision = 8 ):
    l_leadingDimension = 0
    # iterate over all remaning time derivatives
    for l_timeDerivative in range(1, i_associatedOrder+1):
      l_leadingDimension = l_leadingDimension + getAlignedNumberOfBasisFunctions( i_order     = l_timeDerivative,
                                                                                  i_alignment = i_alignment,
                                                                                  i_precision = i_precision )

    return l_leadingDimension

  # Gets the dense matrices for the stiffness matrix computation in the time kernel.
  #
  # @param i_alignment assumed memory alignment of the stiffness matrices and time derivatives of the DOFs.
  # @param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # @param i_numberOfQuantities #(quantities).
  # @param i_precision precision either 's' for single precision or 'd' for double precision
  def getDenseStiffTimeMatrices( self,
                                 i_alignment,
                                 i_degreesOfBasisFunctions,
                                 i_numberOfQuantities,
                                 i_precision ):
    l_denseMatrices = []

    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      # generate dense matrix kernels for stiffness and star matrices in the time kernel
      # We have to compute O-1 derivatives; zeroth derivative is equal to the DOFs
      #
      # Remark: We can't exploit the hierachical character of the basis functions
      #         because the leading dimension of the stiffness matrices changes with the decreasing associated order.
      for l_computedDerivatives in range(0,l_order-1):

        # associated order of the derivatives
        l_associatedOrder = l_order - l_computedDerivatives

        # remaining non zero basis function for the associated order
        l_nonZeroBasisFunctions = (l_associatedOrder)*(l_associatedOrder+1)*(l_associatedOrder+2)/6

        # DGEMM specification
        # m     = #(aligned basis functions of the associated order-1)
        # n     = #(quantities)
        # k     = #(basis functions)
        # ld(A) = #(aligned basis function of the order-1)
        # ld(B) = #(aligned basis function of the associated order)
        # ld(C) = #(aligned basis functions of the associated order-1)
        l_m   = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder-1,
                                                       i_alignment       = i_alignment,
                                                       i_precision       = self.m_configuration.m_bytesPerReal[i_precision] )
        l_n   = i_numberOfQuantities
        l_k   = l_nonZeroBasisFunctions
        l_ldA = self.getAlignedNumberOfBasisFunctions( i_order           = l_order-1,
                                                       i_alignment       = i_alignment,
                                                       i_precision       = self.m_configuration.m_bytesPerReal[i_precision] )
        l_ldB = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder,
                                                       i_alignment       = i_alignment,
                                                       i_precision       = self.m_configuration.m_bytesPerReal[i_precision] )
        l_ldC = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder-1,
                                                       i_alignment       = i_alignment,
                                                       i_precision       = self.m_configuration.m_bytesPerReal[i_precision] )

        l_routineNameOfGeneratedKernel = i_precision + "gemm_m"   + str(l_m)   + "_n"   + str(l_n)   + "_k"   + str(l_k)\
                                                     + "_ldA" + str(l_ldA) + "_ldB" + str(l_ldB) + "_ldC" + str(l_ldC)\
                                                     + "_beta0_pfsigonly";

        l_flops = l_m * l_k * l_n * 2;

        # Add matrix to dictionary
        l_denseMatrices = l_denseMatrices + [ dict(\
                                                routine_name = l_routineNameOfGeneratedKernel, \
                                                m            = l_m,                            \
                                                n            = l_n,                            \
                                                k            = l_k,                            \
                                                ld_a         = l_ldA,                          \
                                                ld_b         = l_ldB,                          \
                                                ld_c         = l_ldC,                          \
                                                flops        = l_flops,                        \
                                                add          = False,                          \
                                                bind         = -1,                             \
                                                prefetch     = 'pfsigonly'
                                              )]

    return l_denseMatrices

  # Gets the dense matrices for the stiffness matrix computation in the volume kernel.
  #
  # @param i_alignment assumed memory alignment of the stiffness matrices and time integrated DOFs.
  # @param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # @param i_numberOfQuantities #(quantities).
  def getDenseStiffVolumeMatrices( self,
                                   i_alignment,
                                   i_degreesOfBasisFunctions,
                                   i_numberOfQuantities,
                                   i_precision ):
    l_denseMatrices = []
    
    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      # generate stiffness matrix kernels for the volume kernel
      #
      # DGEMM specification
      # m     = #(aligned basis functions)
      # n     = #(quantities)
      # k     = #(basis functions of the order-1)
      # ld(A) = #(aligned basis functions)
      # ld(B) = #(aligned basis functions)
      # ld(C) = #(aligned basus functions) 
      l_m   = self.getAlignedNumberOfBasisFunctions( i_order = l_order    , i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_n   = i_numberOfQuantities
      l_k   = self.getNumberOfBasisFunctions(        i_order = l_order - 1                                                                                            )
      l_ldA = self.getAlignedNumberOfBasisFunctions( i_order = l_order    , i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_ldB = self.getAlignedNumberOfBasisFunctions( i_order = l_order    , i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_ldC = self.getAlignedNumberOfBasisFunctions( i_order = l_order    , i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )

      l_routineNameOfGeneratedKernel = i_precision + "gemm_m"   + str(l_m)   + "_n"   + str(l_n)   + "_k"   + str(l_k)\
                                                   + "_ldA" + str(l_ldA) + "_ldB" + str(l_ldB) + "_ldC" + str(l_ldC)\
                                                   + "_beta0_pfsigonly";

      l_flops = l_m * l_k * l_n * 2;

      # Add matrix to dictionary
      l_denseMatrices = l_denseMatrices + [ dict(\
                                              routine_name = l_routineNameOfGeneratedKernel, \
                                              m            = l_m,                            \
                                              n            = l_n,                            \
                                              k            = l_k,                            \
                                              ld_a         = l_ldA,                          \
                                              ld_b         = l_ldB,                          \
                                              ld_c         = l_ldC,                          \
                                              flops        = l_flops,                        \
                                              add          = False,                          \
                                              bind         = -1,                             \
                                              prefetch     = 'pfsigonly'
                                           )]

    return l_denseMatrices

  # Gets the dense matrices for flux matrix computation in the flux kernel.
  #
  # \param i_alignment assumed memory alignment of the flux matrices and time integrated DOFs.
  # \param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # \param i_numberOfQuantities #(quantities).
  def getDenseFluxMatrices( self,
                            i_alignment,
                            i_degreesOfBasisFunctions,
                            i_numberOfQuantities,
                            i_precision, 
                            i_prefetch = 'pfsigonly' ):
    l_denseMatrices = []
    
    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      # generate stiffness matrix kernels for the flux kernel
      #
      # DGEMM specification
      # m     = #(aligned basis functions)
      # n     = #(quantities)
      # k     = #(basis functions of the order)
      # ld(A) = #(aligned basis functions)
      # ld(B) = #(aligned basis functions)
      # ld(C) = #(aligned basus functions) 
      l_m   = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_n   = i_numberOfQuantities
      l_k   = self.getNumberOfBasisFunctions(        i_order = l_order                                                                                            )
      l_ldA = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_ldB = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
      l_ldC = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )

      l_routineNameOfGeneratedKernel = i_precision + "gemm_m"   + str(l_m)   + "_n"   + str(l_n)   + "_k"   + str(l_k)\
                                                   + "_ldA" + str(l_ldA) + "_ldB" + str(l_ldB) + "_ldC" + str(l_ldC)\
                                                   + "_beta0_" + i_prefetch;

      l_flops = l_m * l_k * l_n * 2;

      # Add matrix to dictionary
      l_denseMatrices = l_denseMatrices + [ dict(\
                                              routine_name = l_routineNameOfGeneratedKernel, \
                                              m            = l_m,                            \
                                              n            = l_n,                            \
                                              k            = l_k,                            \
                                              ld_a         = l_ldA,                          \
                                              ld_b         = l_ldB,                          \
                                              ld_c         = l_ldC,                          \
                                              flops        = l_flops,                        \
                                              add          = False,                          \
                                              bind         = -1,                             \
                                              prefetch     = i_prefetch
                                          )]
    return l_denseMatrices           

  # Gets the dense matrices for star matrix and flux solver comptation.
  #
  # @param i_alignment assumed memory alignment of (time integrated, derivatives of) DOFs.
  # @param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # @param i_numberOfQuantities #(quantities).
  def getDenseStarSolverMatrices( self,
                                  i_alignment,
                                  i_degreesOfBasisFunctions,
                                  i_numberOfQuantities,
                                  i_precision, 
                                  i_prefetch = 'pfsigonly'):
    l_denseMatrices = []
    
    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      # generate star matrix and flux matrix kernels for the volume and flux kernels
      #
      # DGEMM specification
      # m     = #(aligned basis functions of the order)
      # n     = #(quantities)
      # k     = #(quantities)
      # ld(A) = #(aligned basis functions of the order)
      # ld(B) = #(quantities)
      # ld(C) = #(aligned basis functions of the order)
      l_m   = self.getAlignedNumberOfBasisFunctions( i_order     = l_order,
                                                     i_alignment = i_alignment,
                                                     i_precision = self.m_configuration.m_bytesPerReal[i_precision] ); 
      l_n   = i_numberOfQuantities
      l_k   = i_numberOfQuantities
      l_ldA = self.getAlignedNumberOfBasisFunctions( i_order     = l_order,
                                                     i_alignment = i_alignment,
                                                     i_precision = self.m_configuration.m_bytesPerReal[i_precision] ); 
      l_ldB = i_numberOfQuantities
      l_ldC = self.getAlignedNumberOfBasisFunctions( i_order     = l_order,
                                                     i_alignment = i_alignment,
                                                     i_precision = self.m_configuration.m_bytesPerReal[i_precision] ); 

      l_routineNameOfGeneratedKernel = i_precision + "gemm_m"   + str(l_m)   + "_n"   + str(l_n)   + "_k"   + str(l_k)\
                                                   + "_ldA" + str(l_ldA) + "_ldB" + str(l_ldB) + "_ldC" + str(l_ldC)\
                                                   + "_beta1_" + i_prefetch;

      l_flops = l_m * l_k * l_n * 2;

      # Add matrix to dictionary
      l_denseMatrices = l_denseMatrices + [ dict(\
                                              routine_name = l_routineNameOfGeneratedKernel, \
                                              m            = l_m,                            \
                                              n            = l_n,                            \
                                              k            = l_k,                            \
                                              ld_a         = l_ldA,                          \
                                              ld_b         = l_ldB,                          \
                                              ld_c         = l_ldC,                          \
                                              flops        = l_flops,                        \
                                              add          = True,                           \
                                              bind         = -1,                             \
                                              prefetch     = i_prefetch
                                          )]

    return l_denseMatrices

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
  # @param i_architectures architectures to generate kernels for.
  # @param i_precision machine precisions code should be generated for
  # @param i_numberOfQuantities number of quantities (elastics = 9, attenuation > 9)
  # @param i_maximumDegreeOfBasisFunctions maximum order of the involved basis functions
  # @return dictionary conatinig the dense matrices described abover.
  def getDenseMatrices( self,
                        i_architectures,
                        i_precision,
                        i_numberOfQuantities,
                        i_maximumDegreeOfBasisFunctions = 8 ):
    
    # list which holds the different matrix structures
    l_denseMatrices = []

    # iterate over different architectures:
    for l_architecture in i_architectures:
      l_alignment = self.m_configuration.m_alignments[l_architecture]
      l_alignedGemm = []

      for l_precision in i_precision:
        l_alignedGemm = l_alignedGemm + self.getDenseStiffTimeMatrices(   i_alignment               = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision )

        l_alignedGemm = l_alignedGemm + self.getDenseStiffVolumeMatrices( i_alignment               = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision )

        l_alignedGemm = l_alignedGemm + self.getDenseFluxMatrices(        i_alignment = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision )
        
        if l_architecture in ['wsm', 'snb', 'hsw', 'skx']:
          l_fluxMatrix_prefetch = 'BL2viaC'
        elif l_architecture in ['knl']:
          l_fluxMatrix_prefetch = 'curAL2_BL2viaC'
        else:
          l_fluxMatrix_prefetch = 'pfsigonly'

        l_alignedGemm = l_alignedGemm + self.getDenseFluxMatrices(        i_alignment = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision,
                                                                          i_prefetch                = l_fluxMatrix_prefetch )

        l_alignedGemm = l_alignedGemm + self.getDenseStarSolverMatrices(  i_alignment               = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision )

        if l_architecture in ['wsm', 'snb', 'hsw', 'skx']:
          l_starSolver_prefetch = 'pfsigonly'
        elif l_architecture in ['knl']:
          l_starSolver_prefetch = 'AL2jpst_BL2viaC'
        else:
          l_starSolver_prefetch = 'pfsigonly'

        l_alignedGemm = l_alignedGemm + self.getDenseStarSolverMatrices(  i_alignment               = l_alignment,
                                                                          i_degreesOfBasisFunctions = range(0,i_maximumDegreeOfBasisFunctions),
                                                                          i_numberOfQuantities      = i_numberOfQuantities,
                                                                          i_precision               = l_precision,
                                                                          i_prefetch                = l_starSolver_prefetch )

        # remove all duplicates which might have been generated (recursive time integration)
        l_alignedGemm =  {l_value['routine_name']:l_value for l_value in l_alignedGemm}.values()

      # file where the generated code is stored
      for l_matrix in range(len(l_alignedGemm)):
        l_alignedGemm[l_matrix]['type'] = "dense"
        l_alignedGemm[l_matrix]['arch'] = l_architecture
        if 'dgemm' in l_alignedGemm[l_matrix]['routine_name']:
          l_alignedGemm[l_matrix]['fileNameOfGeneratedKernel'] =  l_fileNameOfGeneratedKernel = 'dgemm_' + str(l_architecture)
        else:
          l_alignedGemm[l_matrix]['fileNameOfGeneratedKernel'] =  l_fileNameOfGeneratedKernel = 'sgemm_' + str(l_architecture)

      # add the alignment to all matrices
      l_denseMatrices = l_denseMatrices + l_alignedGemm 

    return l_denseMatrices

  # Gets the sparse time matrices for flux computation.
  #
  # @param i_alignment assumed memory alignment of the flux matrices and time integrated DOFs.
  # @param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # @param i_numberOfQuantities #(quantities).
  # @param i_pathToSparseDenseSwitch
  def getSparseTimeMatrices( self,
                             i_alignment,
                             i_degreesOfBasisFunctions,
                             i_numberOfQuantities,
                             i_pathToSparseDenseSwitch,
                             i_precision ):
    # open sparse-dense-switch configuration
    l_sparseDense = open(i_pathToSparseDenseSwitch)

    # no configuration present continue
    if( not l_sparseDense ):
      # command line interface
      logging.info( "time matrices: no sparse-dense configuration, using dense only" )
      return []

    # parse xml file
    l_sparseDense = xmltodict.parse( l_sparseDense )['sparse_matrices']

    # generate a single list for the kernels
    l_sparseMatrices = []
    
    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      l_setup = "O"+str(l_order)

      if l_setup in l_sparseDense:
        # generate sparse matrix kernels for stiffness and star matrices in the time kernel
        # We have to compute O-1 derivatives; zeroth derivative is equal to the DOFs
        #
        # Remark: We can't exploit the hierachical character of the basis functions
        #         because the leading dimension of the stiffness matrices changes with the decreasing associated order.
        for l_computedDerivatives in range(0,l_order-1):
          # associated order of the derivatives
          l_associatedOrder = l_order - l_computedDerivatives

          # iterate over all sparse time matrices and generate kernel specifications
          for l_sparseMatrix in ( l_sparseDense[l_setup]['time'].keys() if l_sparseDense[l_setup]['time'] != None else [] ):
            # compute bind id
            l_bindId = l_computedDerivatives * 4 + self.m_configuration.m_matrixBinds['time'][l_sparseMatrix]

            # DGEMM specification
            # m     = #(basis functions of the associated order-1) if the star matrix is sparse
            #                                        or
            #         #(aligned basis function of the associated order-1) if the star matrix is dense -> zeros are imposed for the oncoming computation
            # n     = #(quantities)
            # k     = sparse, not available: -order
            # ld(A) = sparse, not available: -1
            # ld(B) = #(aligned basis function of the associated order)
            # ld(C) = #(aligned basis functions of the associated order-1)

            if( "starMatrix" in ( l_sparseDense[l_setup]['local'].keys() if l_sparseDense[l_setup]['local'] != None else [] ) ):
              l_m   = self.getNumberOfBasisFunctions(        i_order           = l_associatedOrder-1                        )
            else:
              l_m   = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder-1,
                                                             i_alignment       = i_alignment,
                                                             i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
            l_n   = i_numberOfQuantities
            l_k   = self.getNumberOfBasisFunctions(          i_order           = l_associatedOrder                          )
            l_ldA = -(l_degree+1)
            l_ldB = self.getAlignedNumberOfBasisFunctions(   i_order           = l_associatedOrder,
                                                             i_alignment       = i_alignment,
                                                             i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
            l_ldC = self.getAlignedNumberOfBasisFunctions(   i_order           = l_associatedOrder-1,
                                                             i_alignment       = i_alignment,
                                                             i_precision = self.m_configuration.m_bytesPerReal[i_precision] )

            l_routineName = i_precision               +\
                            "sparse"                  +\
                            "_"      + l_sparseMatrix +\
                            "_m"     + str(l_m)       +\
                            "_n"     + str(l_n)       +\
                            "_k"     + str(l_k)       +\
                            "_ldAna" + str(-l_ldA)    +\
                            "_ldB"   + str(l_ldB)     +\
                            "_ldC"   + str(l_ldC)     +\
                            "_beta0_pfsigonly"

            l_flops = self.m_configuration.m_nonZeros[l_associatedOrder][l_sparseMatrix] * l_n *2

            # Add matrix to dictionary
            l_sparseMatrices = l_sparseMatrices + [ dict(
                                                      name          = l_sparseMatrix,
                                                      matrix_market = self.m_configuration.m_matrixMarketFiles[l_order][l_sparseMatrix],
                                                      routine_name  = l_routineName,
                                                      m             = l_m,
                                                      n             = l_n,
                                                      k             = l_k,
                                                      ld_a          = l_ldA,
                                                      ld_b          = l_ldB,
                                                      ld_c          = l_ldC,
                                                      flops         = l_flops,
                                                      add           = False,
                                                      bind          = l_bindId,
                                                      prefetch      = 'pfsigonly'
                                                   )
                                                  ]
          # generate sparse star matrix if requested
          if( "starMatrix" in ( l_sparseDense[l_setup]['local'].keys() if l_sparseDense[l_setup]['local'] != None else [] ) ):
            # compute bind id
            l_bindId = l_computedDerivatives * 4 + self.m_configuration.m_matrixBinds['time']['starMatrix']

            # DGEMM specification
            # m     = #(basis functions of the associated order-1)
            # n     = #(quantities)
            # k     = #(quantities)
            # ld(A) = #(aligned basis function of the order-1)
            # ld(B) = not available (sparse): -order
            # ld(C) = #(aligned basis functions of the associated order-1)
            l_m   = self.getNumberOfBasisFunctions(        i_order           = l_associatedOrder-1 )
            l_n   = i_numberOfQuantities
            l_k   = i_numberOfQuantities
            l_ldA = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder-1,
                                                           i_alignment       = i_alignment,
                                                           i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
            l_ldB = -(l_degree+1)
            l_ldC = self.getAlignedNumberOfBasisFunctions( i_order           = l_associatedOrder-1,
                                                           i_alignment       = i_alignment,
                                                           i_precision = self.m_configuration.m_bytesPerReal[i_precision] )

            l_routineName = i_precision             +\
                            "sparse"                +\
                            "_"      + "starMatrix" +\
                            "_m"     + str(l_m)     +\
                            "_n"     + str(l_n)     +\
                            "_k"     + str(l_k)     +\
                            "_ldA"   + str(l_ldA)   +\
                            "_ldBna" + str(-l_ldB)  +\
                            "_ldC"   + str(l_ldC)   +\
                            "_beta1_pfsigonly"

            l_flops = self.m_configuration.m_nonZeros[l_associatedOrder]['starMatrix'] * l_m * 2

            # Add matrix to dictionary
            l_sparseMatrices = l_sparseMatrices + [ dict(
                                                      name          = 'starMatrix',
                                                      matrix_market = self.m_configuration.m_matrixMarketFiles['starMatrix'],
                                                      routine_name  = l_routineName,
                                                      m             = l_m,
                                                      n             = l_n,
                                                      k             = l_k,
                                                      ld_a          = l_ldA,
                                                      ld_b          = l_ldB,
                                                      ld_c          = l_ldC,
                                                      flops         = l_flops,
                                                      add           = True,
                                                      bind          = l_bindId,
                                                      prefetch      = 'pfsigonly'
                                                   )
                                                  ]
    return l_sparseMatrices

  # Gets the sparse volume and flux matrices for flux computation.
  #
  # @param i_alignment assumed memory alignment of the flux matrices and time integrated DOFs.
  # @param i_degreesOfBasisFunctions list of degrees of the basis functions for which kernels are generated.
  # @param i_numberOfQuantities #(quantities).
  # @param i_pathToSparseDenseSwitch
  # @param i_integrationKernels get either the matrices for both kernels ['volume', 'boundary']  or only one: ['volume'] or ['boundary']
  def getSparseVolumeAndFluxMatrices( self,
                                      i_alignment,
                                      i_degreesOfBasisFunctions,
                                      i_numberOfQuantities,
                                      i_pathToSparseDenseSwitch,
                                      i_precision,
                                      i_integrationKernels = ['volume', 'boundary'] ):

    # open sparse-dense-switch configuration
    l_sparseDense = open(i_pathToSparseDenseSwitch)

    # no configuration present continue
    if( not l_sparseDense ):
      # command line interface
      logging.info( "volume and flux matrices: no sparse-dense configuration, using dense only" )
      return []

    # parse xml file
    l_sparseDense = xmltodict.parse( l_sparseDense )['sparse_matrices']

    # generate a single list for the kernels
    l_sparseMatrices = []
    
    for l_degree in i_degreesOfBasisFunctions:
      l_order = l_degree + 1

      l_setup = "O"+str(l_order)

      if l_setup in l_sparseDense:
        l_sparseNames = []
        for l_kernel in i_integrationKernels:
          if( l_sparseDense[l_setup][l_kernel] ):
            l_sparseNames = list( l_sparseDense[l_setup][l_kernel].keys() )

          for l_sparseMatrix in l_sparseNames:
            # continue if this matrix does not belong to the current kernel
            if( l_sparseMatrix not in self.m_configuration.m_matrixBinds[l_kernel].keys() ):
              continue

            # set up bind id
            l_bindId = self.m_configuration.m_matrixBinds[l_kernel][l_sparseMatrix]

            # generate sparse matrix kernels for the flux kernel
            #
            # DGEMM specification
            # m     = #(aligned basis functions)
            # n     = #(quantities)
            # k     = not available: -order
            # ld(A) = not available: -1
            # ld(B) = #(aligned basis functions)
            # ld(C) = #(aligned basus functions) 
            l_m   = self.getNumberOfBasisFunctions( i_order = l_order )
            l_n   = i_numberOfQuantities
            l_k   = self.getNumberOfBasisFunctions( i_order = l_order )
            l_ldA = -(l_degree+1)
            l_ldB = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
            l_ldC = self.getAlignedNumberOfBasisFunctions( i_order = l_order, i_alignment = i_alignment, i_precision = self.m_configuration.m_bytesPerReal[i_precision] )

            l_routineName = i_precision               +\
                            "sparse"                  +\
                            "_"      + l_sparseMatrix +\
                            "_m"     + str(l_m)       +\
                            "_n"     + str(l_n)       +\
                            "_k"     + str(l_k)       +\
                            "_ldAna" + str(-l_ldA)    +\
                            "_ldB"   + str(l_ldB)     +\
                            "_ldC"   + str(l_ldC)     +\
                            "_beta0_pfsigonly"

            l_flops = self.m_configuration.m_nonZeros[l_order][l_sparseMatrix] * l_n * 2

            # Add matrix to dictionary
            l_sparseMatrices = l_sparseMatrices + [ dict(
                                                      name          = l_sparseMatrix,
                                                      matrix_market = self.m_configuration.m_matrixMarketFiles[l_order][l_sparseMatrix],
                                                      routine_name  = l_routineName,
                                                      m             = l_m,
                                                      n             = l_n,
                                                      k             = l_k,
                                                      ld_a          = l_ldA,
                                                      ld_b          = l_ldB,
                                                      ld_c          = l_ldC,
                                                      flops         = l_flops,
                                                      add           = False,
                                                      bind          = l_bindId,
                                                      prefetch      = 'pfsigonly'
                                                   )
                                                  ]
          # generate sparse star matrix if requested
          if( l_kernel == 'volume' and l_sparseDense[l_setup]['local'] != None and "starMatrix" in l_sparseDense[l_setup]['local'].keys() ):
            # set up bind id
            l_bindId = self.m_configuration.m_matrixBinds[l_kernel]['starMatrix']

            # DGEMM specification
            # m     = #(basis functions of the order)
            # n     = #(quantities)
            # k     = #(quantities)
            # ld(A) = #(aligned basis function of the order)
            # ld(B) = not available (sparse): -order
            # ld(C) = #(aligned basis functions of the order)
            l_m   = self.getNumberOfBasisFunctions(        i_order           = l_order                                    )
            l_n   = i_numberOfQuantities
            l_k   = i_numberOfQuantities
            l_ldA = self.getAlignedNumberOfBasisFunctions( i_order           = l_order,
                                                           i_alignment       = i_alignment,
                                                           i_precision = self.m_configuration.m_bytesPerReal[i_precision] )
            l_ldB = -(l_degree+1)
            l_ldC = self.getAlignedNumberOfBasisFunctions( i_order           = l_order,
                                                           i_alignment       = i_alignment,
                                                           i_precision = self.m_configuration.m_bytesPerReal[i_precision]   )

            l_routineName = i_precision             +\
                            "sparse"                +\
                            "_"      + "starMatrix" +\
                            "_m"     + str(l_m)     +\
                            "_n"     + str(l_n)     +\
                            "_k"     + str(l_k)     +\
                            "_ldA"   + str(l_ldA)   +\
                            "_ldBna" + str(-l_ldB)  +\
                            "_ldC"   + str(l_ldC)   +\
                            "_beta1_pfsigonly"

            l_flops = self.m_configuration.m_nonZeros[l_order]['starMatrix'] * l_m * 2

            # Add matrix to dictionary
            l_sparseMatrices = l_sparseMatrices + [ dict(
                                                      name          = "starMatrix",
                                                      matrix_market = self.m_configuration.m_matrixMarketFiles['starMatrix'],
                                                      routine_name  = l_routineName,
                                                      m             = l_m,
                                                      n             = l_n,
                                                      k             = l_k,
                                                      ld_a          = l_ldA,
                                                      ld_b          = l_ldB,
                                                      ld_c          = l_ldC,
                                                      flops         = l_flops,
                                                      add           = True,
                                                      bind          = l_bindId,
                                                      prefetch      = 'pfsigonly'
                                                   )
                                                  ]
    return l_sparseMatrices

  # Returns a list of sparse matrix kernels for SeisSol.
  # @param i_architectures architectures to get sparse matrices for.
  # @param i_precision machine precisions code should be generated for
  # @param i_numberOfQuantities number of quantities (elastics = 9, attenuation > 9)
  # @param i_maximumDegreeOfBasisFunctions maximum order of the involved basis functions
  # @return dictionary containing the sparse matrices described.
  def getSparseMatrices( self,
                         i_architectures,
                         i_precision,
                         i_numberOfQuantities,
                         i_maximumDegreeOfBasisFunctions = 8 ):

    # list which holds the different matrix structures
    l_sparseMatrices = []

    l_minimumGlobalDegree   = 0
    l_minimumLocalDegree    = 0

    for l_architecture in i_architectures:
      l_alignment = self.m_configuration.m_alignments[l_architecture]

      # path to the sparse dense configuration
      l_pathToSparseDenseSwitch = self.m_configuration.m_pathToSparseDenseConfigs + '/' + i_precision + l_architecture + ".xml"

      for l_precision in i_precision:
        l_alignedSparse = []

        l_alignedSparse = l_alignedSparse + self.getSparseTimeMatrices(          i_alignment               = l_alignment,
                                                                                 i_degreesOfBasisFunctions = range(l_minimumGlobalDegree,i_maximumDegreeOfBasisFunctions),
                                                                                 i_pathToSparseDenseSwitch = l_pathToSparseDenseSwitch,
                                                                                 i_numberOfQuantities      = i_numberOfQuantities,
                                                                                 i_precision               = l_precision )

        l_alignedSparse = l_alignedSparse + self.getSparseVolumeAndFluxMatrices( i_alignment               = l_alignment,
                                                                                 i_degreesOfBasisFunctions = range(l_minimumLocalDegree,i_maximumDegreeOfBasisFunctions),
                                                                                 i_pathToSparseDenseSwitch = l_pathToSparseDenseSwitch,
                                                                                 i_numberOfQuantities      = i_numberOfQuantities,
                                                                                 i_precision               = l_precision )

        # file where the generated code is stored
        l_fileNameOfGeneratedKernel = 'sparse_' + l_precision + str(l_architecture)
        for l_matrix in range(len(l_alignedSparse)):
          l_alignedSparse[l_matrix]['type'] = "sparse"
          l_alignedSparse[l_matrix]['arch'] = l_architecture
          l_alignedSparse[l_matrix]['fileNameOfGeneratedKernel'] = l_fileNameOfGeneratedKernel

        # add the alignment to all matrices
        l_sparseMatrices = l_sparseMatrices + l_alignedSparse 

    return l_sparseMatrices
