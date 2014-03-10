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
 * Time integration in SeisSol.
 **/

#include <Monitoring/FlopCounter.hpp>
#include <utils/logger.h>

#ifndef NDEBUG
#pragma message "compiling time integrator with assertions"
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <cstring>
#include <assert.h>

#include <generated_code/matrix_kernels/dense_matrices.hpp_include>
#include <generated_code/matrix_kernels/star_matrices_3d.hpp_include>
#ifndef __MIC__
#include <generated_code/matrix_kernels/stiffness_matrices_3d.hpp_include>
#endif
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include "TimeIntegrator.h"


void seissol::kernels::TimeIntegrator::setUpMatrixKernel( unsigned int i_id,
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
    m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_denseAder_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse transposed stiffness matrix \f$ M^{-1} K^\xi \f$
  else if( i_id == 0 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kXiDivMT_,   NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse transposed stiffness matrix \f$ M^{-1} K^\eta \f$ 
  else if( i_id == 1 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kEtaDivMT_,  NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
  // sparse transposed stiffness matrix \f$ M^{-1} K^\zeta \f$ 
  else if( i_id == 2 ){
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_kZetaDivMT_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
#else
  // dense stiffness matrices only for Intel MIC
  if( i_id < 3 ){
    m_matrixKernels[i_id] = CONCAT_6( generatedMatrixMultiplication_denseAder_, NUMBEROFBASISFUNCTIONS, _, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
#endif
  // sparse star matrix
  else{
    assert( i_sparse == true && i_id == 3 );
    m_matrixKernels[i_id] = CONCAT_4( generatedMatrixMultiplication_aderStarMatrix_3D_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS );
  }
}

#ifndef DIRTY_EXCLUDE_ON_MIC
seissol::kernels::TimeIntegrator::TimeIntegrator( const seissol::XmlParser                   &i_matrixReader,
                                                  const seissol::initializers::MemoryManager &i_memoryManager ) {
  /*
   * Initialize transposed stiffness matrix pointers (jump over the first three non-transposed stiffness matrices).
   */
  m_stiffnessMatrixPointers = i_memoryManager.getStiffnessMatrixPointers()+3;

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
  // Remark: The last three entries correspond to the transposed stiffness matrices used in the time integration.
  for( int l_i = 0; l_i < 3; l_i++) {
    // assert we have the right matrices
    assert( l_matrixIds[l_i+3] > 55 && l_matrixIds[l_i+3] < 59 );

    // set up matrix kernels of the stiffness matrices
    setUpMatrixKernel( l_i,
                       l_matrixSparsities[l_i+3] );
  }

  // set up the star matrix kernel
  setUpMatrixKernel( 3, true );
}
#endif

#ifdef __INTEL_OFFLOAD
seissol::kernels::TimeIntegrator::TimeIntegrator() {
  // set up the stiffness matrices as dense matrices
  for( int l_i = 0; l_i < 3; l_i++) {
    setUpMatrixKernel( l_i, false );
  }

  // set up the star matrix kernel
  setUpMatrixKernel( 3, true );
}
#endif

void seissol::kernels::TimeIntegrator::computeTimeDerivatives( const double  i_unknowns[NUMBEROFUNKNOWNS],
                                                                     double *i_stiffnessMatrices[3],
                                                                     double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                     double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                     double  i_cStar[STARMATRIX_NUMBEROFNONZEROS], 
                                                                     double  o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] ) {
  /*
   * Assert alignments, which are assumed in the matrix kernels.
   */
  //TODO: What do we assume in the matrix kernels??

  /*
   * Computation
   */
  // non-zero block sizes in stiffness and star matrix multiplications
  int l_nonZeroBlockSizeStiffness;
  int l_nonZeroBlockSizeStar = NUMBEROFBASISFUNCTIONS;

  // temporary matrix for summation (64 byte aligned)
  double l_partialSum[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

  // copy unknowns to zeroth derivative 
  for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
    o_timeDerivatives[0][l_entry] = i_unknowns[l_entry];
  }

  // iterate over order in time
  for( int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++) {
    // set the this derivative to zero
    // TODO: Replace by C=A.B matrix operators
    for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
      o_timeDerivatives[l_order][l_entry] = 0;
    }    

    // save non-zero block size of the siffness multiplcations
    l_nonZeroBlockSizeStiffness = l_nonZeroBlockSizeStar;

    /* Compute non-zero block size for the star matrix multiplications in this step and the
     * stiffness matrix multiplications in the next step.
     *
     * A detailed description can be found in the global time stepping functionality of the time integrator (see below).
     */
    l_nonZeroBlockSizeStar  = ORDEROFTAYLORSERIESEXPANSION - l_order;
    l_nonZeroBlockSizeStar  = l_nonZeroBlockSizeStar * (l_nonZeroBlockSizeStar + 1) * (l_nonZeroBlockSizeStar + 2);
    l_nonZeroBlockSizeStar /= 6;

    // calculate $K_\xi.Q_k$ and $(K_\xi.Q_k).A*$
    m_matrixKernels[0] ( i_stiffnessMatrices[0], o_timeDerivatives[l_order-1], l_partialSum,               l_nonZeroBlockSizeStiffness,
                         NULL,                   NULL,                         NULL                                                     ); // TODO: prefetches
    m_matrixKernels[3] ( l_partialSum,           i_aStar,                      o_timeDerivatives[l_order], l_nonZeroBlockSizeStar,
                         NULL,                   NULL,                         NULL                                                     ); // TODO: prefetches

    // calculate $K_\eta.Q_k$ and $(K_\eta.Q_k).B*$
    m_matrixKernels[1] ( i_stiffnessMatrices[1], o_timeDerivatives[l_order-1], l_partialSum,               l_nonZeroBlockSizeStiffness,
                         NULL,                   NULL,                         NULL                                                     ); //TODO: prefetches
    m_matrixKernels[3] ( l_partialSum,           i_bStar,                      o_timeDerivatives[l_order], l_nonZeroBlockSizeStar,
                         NULL,                   NULL,                         NULL                                                     ); //TODO: prefetches

    // calculate $K_\zeta.Q_k$ and calculate $(K_\zeta.Q_k).C*$
    m_matrixKernels[2] ( i_stiffnessMatrices[2], o_timeDerivatives[l_order-1], l_partialSum,               l_nonZeroBlockSizeStiffness,
                         NULL,                   NULL,                         NULL                                                     ); //TODO: prefetches
    m_matrixKernels[3] ( l_partialSum,           i_cStar,                      o_timeDerivatives[l_order], l_nonZeroBlockSizeStar,
                         NULL,                   NULL,                         NULL                                                     ); //TODO: prefetches
  }
}

void seissol::kernels::TimeIntegrator::computeTimeEvaluation( const double  i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                                              const double &i_expansionPoint,
                                                              const double &i_evaluationPoint,
                                                                    double  o_unknowns[NUMBEROFUNKNOWNS] ) {
  // assert that non-backward evaluation in time
  assert( i_expansionPoint >= i_evaluationPoint );

  // initialization of factors in the Taylor series expansion (1st term)
  double l_taylorSeriesDelta  = i_evaluationPoint - i_expansionPoint;
  double l_taylorSeriesFactor = 1.0; 

  // set unknowns to unkowns at time of expansion point (0th derivative)
  for(int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++) {
    o_unknowns[l_entry] = i_timeDerivatives[0][l_entry];
  }
 
 
  // evaluate taylor series expansion at evvaluation point
  for(int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) { 
    //   Remark: the negative sign of the derivative is include here
    l_taylorSeriesFactor = -( l_taylorSeriesFactor*l_taylorSeriesDelta ) / double(l_order);

    // update the unknowns by the contribution of this derivative
    for(int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++) {
      o_unknowns[l_entry] += l_taylorSeriesFactor * i_timeDerivatives[l_order][l_entry];
    }
  }
}

void seissol::kernels::TimeIntegrator::computeTimeIntegral( const double  i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                                            const double &i_deltaTLower,
                                                            const double &i_deltaTUpper,
                                                                  double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
  // assert that this is a forwared integration in time
  assert( i_deltaTUpper > i_deltaTLower );

  // initialization of factors in the Taylor series expansion (0th term)
  double l_taylorSeriesDelta = i_deltaTUpper - i_deltaTLower;
  double l_taylorSeriesFactor = l_taylorSeriesDelta; 

  // update the time integrated unknowns by the contribution of the 0th derivative
  for(int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++) {
    o_timeIntegratedUnknowns[l_entry] = l_taylorSeriesFactor * i_timeDerivatives[0][l_entry];
  }
 
 
  // iterate over order in time
  for(int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) { 
    // compute $\frac{(\Delta t^\text{up})^{j+1} - (\Delta t^\text{lo}^{j+1})}{(j+1)!}$ for the current term (l_order == j in the code)
    //   Remark: the negative sign of the derivative is include here
    l_taylorSeriesFactor = -( l_taylorSeriesFactor*l_taylorSeriesDelta ) / double(l_order + 1);

    // update the time integrated unknowns by the contribution of this derivative
    for(int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++) {
      o_timeIntegratedUnknowns[l_entry] += l_taylorSeriesFactor * i_timeDerivatives[l_order][l_entry];
    }
  }
}

void seissol::kernels::TimeIntegrator::computeTimeIntegral( const double   i_unknowns[NUMBEROFUNKNOWNS],
                                                                  double  *i_stiffnessMatrices[3],
                                                                  double   i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                  double   i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                                                  double   i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                                            const double  &i_deltaT,
                                                                  double   o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
  /*
   * Assert alignments, which are assumed in the matrix kernels.
   */
#if NUMBEROFBASISFUNCTIONS == 1
  // no alignment
#elif NUMBEROFBASISFUNCTIONS == 4
  // 32 byte alignment
  assert( ((uintptr_t)i_unknowns)               % 32 == 0 );
  assert( ((uintptr_t)o_timeIntegratedUnknowns) % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 10
  // 16 byte alignment
  assert( ((uintptr_t)i_unknowns)               % 16 == 0 );
  assert( ((uintptr_t)o_timeIntegratedUnknowns) % 16 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 20
  // 32 byte alignment
  assert( ((uintptr_t)i_unknowns)               % 32 == 0 );
  assert( ((uintptr_t)o_timeIntegratedUnknowns) % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 35
  // 8 byte alignment
  assert( ((uintptr_t)i_unknowns)               %  8 == 0 );
  assert( ((uintptr_t)o_timeIntegratedUnknowns) %  8 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 56
  // 64 byte alignment
  assert( ((uintptr_t)i_unknowns)               % 64 == 0 );
  assert( ((uintptr_t)o_timeIntegratedUnknowns) % 64 == 0 );
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
  // initialization of the factors in the Taylor series expansion (0th term)
  double l_taylorSeriesFactor = i_deltaT;

  // non-zero block sizes in stiffness and star matrix multiplications
  int l_nonZeroBlockSizeStiffness;
  int l_nonZeroBlockSizeStar = NUMBEROFBASISFUNCTIONS;

  // temporary matrices for summation (64 byte aligned)
  double l_firstPartialSum[NUMBEROFUNKNOWNS]        __attribute__((aligned(64)));
  double l_secondPartialSum[NUMBEROFUNKNOWNS]       __attribute__((aligned(64)));
  double l_differentiatedUnknowns[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));


  // initialize differentiated and time integrated unknowns
  for( int l_i = 0; l_i < NUMBEROFUNKNOWNS; l_i++ ) {
    o_timeIntegratedUnknowns[l_i] = i_unknowns[l_i] * l_taylorSeriesFactor;
    l_differentiatedUnknowns[l_i] = i_unknowns[l_i];
  }
#ifndef NDEBUG
  // add the previous multiplications to the FLOP counter
  num_flops += NUMBEROFUNKNOWNS;
#endif

  for( int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++) {
    // save non-zero block size of the siffness multiplcations
    l_nonZeroBlockSizeStiffness = l_nonZeroBlockSizeStar;

    /* Compute non-zero block size for the star matrix multiplications in this step and the
     * stiffness matrix multiplications in the next step.
     *
     * Remark: The size of the non-zero time differentiated unknowns reduces in every step.
     *   Sketch for polynomial degree 2:
     * <pre>
     *   First step:
     *       \f$(k^\Xi)^T\f$     * \f$\frac{\partial^0}{\partial t^0} Q_k\f$ * \f$A^*_k\f$ + [...]
     *    _                   _     _                 _     _                 _             _                 _
     *   | - * - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * - - |           | * * * * * * * * * | <----/ remaining non-zero block
     *   | - - - - * - - - - - |   | * * * * * * * * * |   | - - - - - - - * - |           | * * * * * * * * * | <---/  size:
     *   | - - - - - * - - - - |   | * * * * * * * * * |   | - - - - - - - - * |           | * * * * * * * * * | <--/     #basis functions of polynomial degree 1
     *   | - - - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * * - |           | * * * * * * * * * | <-/
     *   | - - - - - - - - - - | * | * * * * * * * * * | * | - - - - - - - * * | + [...] = | - - - - - - - - - |    Remark: The star matrix multiplication is one
     *   | - - - - - - - - - - |   | * * * * * * * * * |   | - - - - - - * - * |           | - - - - - - - - - |            step ahead because of associativity
     *   | - - - - - - - - - - |   | * * * * * * * * * |   | * * * * - * - - - |           | - - - - - - - - - |            of the matrix multiplication.
     *   | - - - - - - - - - - |   | * * * * * * * * * |   | * * * * * - - - - |           | - - - - - - - - - |
     *   | - - - - - - - - - - |   | * * * * * * * * * |   |_* * * - * * - - -_|           |_- - - - - - - - -_|
     *   |_- - - - - - - - - -_|   |_* * * * * * * * *_|
     *
     *   Second step:
     *        \f$(k^\Xi)^T\f$     * \f$\frac{\partial^1}{\partial t^1} Q_k\f$ * \f$A^*_k\f$ + [...] 
     *    _                   _     _                 _     _                 _
     *   | - * - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * - - |
     *   | - - - - * - - - - - |   | * * * * * * * * * |   | - - - - - - - * - |
     *   | - - - - - * - - - - |   | * * * * * * * * * |   | - - - - - - - - * |
     *   | - - - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * * - |
     *   | - - - - - - - - - - | * | - - - - - - - - - | * | - - - - - - - * * | + [...] = [...]
     *   | - - - - - - - - - - |   | - - - - - - - - - |   | - - - - - - * - * |
     *   | - - - - - - - - - - |   | - - - - - - - - - |   | * * * * - * - - - |
     *   | - - - - - - - - - - |   | - - - - - - - - - |   | * * * * * - - - - |
     *   | - - - - - - - - - - |   | - - - - - - - - - |   |_* * * - * * - - -_|
     *   |_- - - - - - - - - -_|   |_- - - - - - - - -_|
     *             ___________
     *            /     |     \
     *                  |
     *         columns, which hit zeros in the second step
     *      -> the product has only one non-zero row left
     * </pre>
     */
    l_nonZeroBlockSizeStar  = ORDEROFTAYLORSERIESEXPANSION - l_order;
    l_nonZeroBlockSizeStar  = l_nonZeroBlockSizeStar * (l_nonZeroBlockSizeStar + 1) * (l_nonZeroBlockSizeStar + 2);
    l_nonZeroBlockSizeStar /= 6;

    // store differentiated unknowns (each iteration calculates a new derivative)
    if (l_order > 1)
      memcpy( l_differentiatedUnknowns, l_secondPartialSum, NUMBEROFUNKNOWNS*sizeof(*l_secondPartialSum) );

    // reset second partial sum
    memset( l_secondPartialSum, 0, NUMBEROFUNKNOWNS*sizeof(*l_secondPartialSum) );

    // calculate $K_\xi.Q_k$ and $(K_\xi.Q_k).A*$
    m_matrixKernels[0] ( i_stiffnessMatrices[0], l_differentiatedUnknowns, l_firstPartialSum,  l_nonZeroBlockSizeStiffness,
                         l_firstPartialSum,      i_aStar,                  l_secondPartialSum                               ); // prefetches
    m_matrixKernels[3] ( l_firstPartialSum,      i_aStar,                  l_secondPartialSum, l_nonZeroBlockSizeStar,
                         i_stiffnessMatrices[1], l_differentiatedUnknowns, l_firstPartialSum                                ); // prefetches

    // calculate $K_\eta.Q_k$ and $(K_\eta.Q_k).B*$
    m_matrixKernels[1] ( i_stiffnessMatrices[1], l_differentiatedUnknowns, l_firstPartialSum,  l_nonZeroBlockSizeStiffness,
                         l_firstPartialSum,      i_bStar,                  l_secondPartialSum                               ); // prefetches
    m_matrixKernels[3] ( l_firstPartialSum,      i_bStar,                  l_secondPartialSum, l_nonZeroBlockSizeStar,
                         i_stiffnessMatrices[2], l_differentiatedUnknowns, l_firstPartialSum                                ); // prefetches

    // calculate $K_\zeta.Q_k$ and calculate $(K_\zeta.Q_k).C*$
    m_matrixKernels[2] ( i_stiffnessMatrices[2], l_differentiatedUnknowns, l_firstPartialSum,  l_nonZeroBlockSizeStiffness,
                         l_firstPartialSum,      i_cStar,                  l_secondPartialSum                               ); // prefetches
    m_matrixKernels[3] ( l_firstPartialSum,      i_cStar,                  l_secondPartialSum, l_nonZeroBlockSizeStar,
                         i_stiffnessMatrices[0], l_differentiatedUnknowns, l_firstPartialSum                                ); // prefetches; TODO: Upcoming time int.

    // compute $\frac{(t^{n+1} - t^n)^{j+1}}{(j+1)!}$ for the current term (l_order == j in the code)
    //   Remark: the negative sign of the derivative is include here
    l_taylorSeriesFactor = -l_taylorSeriesFactor * i_deltaT / double(l_order+1);

    // update the time integrated unknowns
    for( int l_i = 0; l_i < NUMBEROFUNKNOWNS; l_i++ )
        o_timeIntegratedUnknowns[l_i] += l_secondPartialSum[l_i] * l_taylorSeriesFactor;
#ifndef NDEBUG
    // add the previous multiply-adds to the FLOP counter
    num_flops += NUMBEROFUNKNOWNS * 2;
#endif
  }

#ifndef NDEBUG
  // update flop counter
  addTimeFlops();
#endif
}
