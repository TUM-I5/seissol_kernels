/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2013
 * Technische Universitaet Muenchen
 * Department of Informatics
 * Chair of Scientific Computing
 * http://www5.in.tum.de/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
 * Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Volume integration in SeisSol.
 **/
#include "VolumeIntegrator.h"

void seissol::kernels::VolumeIntegrator::setUpMatrixKernel( unsigned int i_id,
                                                            bool i_sparse ) {
  // assert we are not out of bounds
  assert( i_id < 4 );

  // assert no dense kernel for the star matrices is requested
  assert( i_sparse == true || i_id < 3 );

  /*
   * set up the the function pointers
   */
#ifndef __MIC__
  if( i_sparse == false && i_id < 3 ){
    m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse stiffness matrix \f$ K^\xi \f$
  else if( i_id == 0 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kXiDivM_,   NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse stiffness matrix \f$ K^\eta \f$ 
  else if( i_id == 1 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kEtaDivM_,  NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse stiffness matrix \f$ K^\zeta \f$ 
  else if( i_id == 2 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kZetaDivM_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
#else
  // dense stiffness matrices only for Intel MIC
  if( i_id < 3 ){
    m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_dense_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
#endif
  // sparse star matrix
  else{
    assert( i_sparse == true && i_id == 3 );
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_volumeStarMatrix_3D_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
}

seissol::kernels::VolumeIntegrator::VolumeIntegrator( const seissol::XmlParser                   &i_matrixReader,
                                                      const seissol::initializers::MemoryManager &i_memoryManager ) {
  /*
   * Initialize stiffness matrix pointers.
   */
  m_stiffnessMatrixPointers = i_memoryManager.getStiffnessMatrixPointers();

  /*
   * Initialize function pointers to the matrix kernels
   */
  //! vectors, which hold information about our matrices
  std::vector< unsigned int > l_matrixIds;
  std::vector< std::string  > l_matrixNames;
  std::vector< unsigned int > l_matrixNumberOfRows;
  std::vector< unsigned int > l_matrixNumberOfColumns;
  std::vector< bool         > l_matrixSparsities;

  // element information in coordinate format
  std::vector< std::vector<unsigned int> > l_matrixRows;
  std::vector< std::vector<unsigned int> > l_matrixColumns;
  std::vector< std::vector<double>       > l_matrixValues;

  // read the flux matrices
  i_matrixReader.readGlobalMatrices( "stiffness",
                                     l_matrixIds,  l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );

  // assert we have all matrices: 3 stiffness matrices given transposed and non-transposed.
  assert( l_matrixIds.size() == 6 );

  // set up the stiffness matrices
  // Remark: The first three entries correspond to the non-transposed stiffness matrices used in the volume integration.
  for( int l_i = 0; l_i < 3; l_i++) {
    // assert we have the right matrices
    assert( l_matrixIds[l_i] > 52 && l_matrixIds[l_i] < 56 );

    // set up matrix kernels of the stiffness matrices
    setUpMatrixKernel( l_i,
                       l_matrixSparsities[l_i] );
  }

  // set up the star matrix kernel
  setUpMatrixKernel( 3, true );
}

void seissol::kernels::VolumeIntegrator::computeVolumeIntegral( double i_timeIntegratedUnknowns[NUMBEROFUNKNOWNS],
                                                                double *i_stiffnessMatrices[3],
                                                                double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                double i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                double io_unknowns[NUMBEROFUNKNOWNS] ) {
  /*
   * Assert alignments, which are assumed in the matrix kernels.
   */
#if NUMBEROFBASISFUNCTIONS == 1
  // no alignment
#elif NUMBEROFBASISFUNCTIONS == 4
  // 32 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknowns) % 32 == 0 );
  assert( ((uintptr_t)io_unknowns)              % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 10
  // 16 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknowns) % 16 == 0 );
  assert( ((uintptr_t)io_unknowns)              % 16 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 20
  // 32 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknowns) % 32 == 0 );
  assert( ((uintptr_t)io_unknowns)              % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 35
  // 8 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknowns) %  8 == 0 );
  assert( ((uintptr_t)io_unknowns)              %  8 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 56
  // 64 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknowns) % 64 == 0 );
  assert( ((uintptr_t)io_unknowns)              % 64 == 0 );
#else
#error Preprocessor flag NUMBEROFBASISFUNCTIONS is not in {1, 4, 10, 20, 35, 56}.
#endif

#ifndef NDEBUG
  // 64 byte alignment of stiffness matrices
  for( int l_matrix = 0; l_matrix < 3; l_matrix++ ) {
    assert( ((uintptr_t)i_stiffnessMatrices[l_matrix]) % 64 == 0 );
  }
#endif

  /*
   * Computation
   */
  // temporary matrix for two-step multiplications (64 byte aligned)
  double l_partialProduct[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

  // TODO: Switch to a loop when the memory manager is able to take care about the star matrices.

  // reset temporary matrix
  memset( l_partialProduct, 0, NUMBEROFUNKNOWNS*sizeof(*l_partialProduct) );

  // calculate $K_\xi.I(Q_k, t^n, t^{n+1}$ and $(K_\xi.I(Q_k, t^n, t^{n+1}).A^*$
  m_matrixKernels[0] ( i_stiffnessMatrices[0], i_timeIntegratedUnknowns, l_partialProduct,
                       l_partialProduct,       i_aStar,                  io_unknowns       ); // prefetches
  m_matrixKernels[3] ( l_partialProduct,       i_aStar,                  io_unknowns,
                       i_stiffnessMatrices[1], i_timeIntegratedUnknowns, l_partialProduct  ); // prefetches

  // reset temporary matrix
  memset( l_partialProduct, 0, NUMBEROFUNKNOWNS*sizeof(*l_partialProduct) );

  // calculate $K_\eta.I(Q_k, t^n, t^{n+1}$ and $(K_\eta.I(Q_k, t^n, t^{n+1}).B^*$
  m_matrixKernels[1] ( i_stiffnessMatrices[1], i_timeIntegratedUnknowns, l_partialProduct,
                       l_partialProduct,       i_bStar,                  io_unknowns       ); // prefetches
  m_matrixKernels[3] ( l_partialProduct,       i_bStar,                  io_unknowns,
                       i_stiffnessMatrices[2], i_timeIntegratedUnknowns, l_partialProduct  ); // prefetches

  // reset temporary matrix
  memset( l_partialProduct, 0, NUMBEROFUNKNOWNS*sizeof(*l_partialProduct) );

  // calculate $K_\zeta.I(Q_k, t^n, t^{n+1}$ and $(K_\zeta.I(Q_k, t^n, t^{n+1}).C^*$
  m_matrixKernels[2] ( i_stiffnessMatrices[2], i_timeIntegratedUnknowns, l_partialProduct,
                       l_partialProduct,       i_cStar,                  io_unknowns       ); // prefetches 
  m_matrixKernels[3] ( l_partialProduct,       i_cStar,                  io_unknowns,
                       i_stiffnessMatrices[6], i_timeIntegratedUnknowns, NULL              ); // inter-kernel prefetches for the bnd. int.

#ifndef NDEBUG
  // update flop counter
  addVolumeFlops();
#endif
}
