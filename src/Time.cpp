/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Time kernel of SeisSol.
 **/

#include "Time.h"

#if ALIGNMENT==32
#include <generated_code/matrix_kernels/dgemm_32.h>
#elif ALIGNMENT==64
#include <generated_code/matrix_kernels/dgemm_64.h>
#else
#error ALIGNMENT not supported
#endif

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
#endif

#include <cstring>
#include <cassert>
#include <cstdint>

seissol::kernels::Time::Time() {
  // compute the aligned number of basis functions and offsets of the derivatives
  m_derivativesOffsets[0] = 0;
  for( int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
    m_numberOfAlignedBasisFunctions[l_order] = getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order, ALIGNMENT );

    if( l_order > 0 ) {
      m_derivativesOffsets[l_order]  =  m_numberOfAlignedBasisFunctions[l_order-1] * NUMBER_OF_QUANTITIES;
      m_derivativesOffsets[l_order] +=  m_derivativesOffsets[l_order-1];
    }
  }

  // intialize the function pointers to the matrix kernels
#define TIME_KERNEL
#include <generated_code/initialization/bind_matrix_kernels.hpp_include>
#undef TIME_KERNEL
}

void seissol::kernels::Time::computeDerivatives(       real** i_stiffnessMatrices,
                                                 const real*  i_degreesOfFreedom,
                                                       real** i_starMatrices,
                                                       real*  o_timeDerivatives ) {
  // assert alignments
  assert( ((uintptr_t)i_degreesOfFreedom)     % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[0]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[1]) % ALIGNMENT == 0 );
  assert( ((uintptr_t)i_stiffnessMatrices[2]) % ALIGNMENT == 0 );

  // reset time derivatives
  // TODO: inefficient
  memset( o_timeDerivatives, 0, getAlignedTimeDerivativesSize()*sizeof(real) );

  // temporary result
  real l_temporaryResult[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(64)));

  // copy DOFs to 0th derivative
  memcpy( o_timeDerivatives, i_degreesOfFreedom, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real) );

  // compute all derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    // iterate over dimensions 
    for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
      // compute $K_{\xi_c}.Q_k$ and $(K_{\xi_c}.Q_k).A*$
      m_matrixKernels[ (l_derivative-1)*2     ] ( i_stiffnessMatrices[l_c], o_timeDerivatives+m_derivativesOffsets[l_derivative-1],  l_temporaryResult,
                                                  NULL,                     NULL,                                                    NULL                                                  ); // prefetches

      m_matrixKernels[ (l_derivative-1)*2 + 1 ] ( l_temporaryResult,        i_starMatrices[l_c],                                     o_timeDerivatives+m_derivativesOffsets[l_derivative],
                                                  NULL,                     NULL,                                                    NULL                                                  ); // prefetches
    }
  }

}

void seissol::kernels::Time::computeExtrapolation(       real   i_expansionPoint,
                                                         real   i_evaluationPoint,
                                                   const real*  i_timeDerivatives,
                                                         real*  o_timeEvaluated ) {
  // assert alignments
  assert( ((uintptr_t)i_timeDerivatives) % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeEvaluated)   % ALIGNMENT == 0 );

  // assert that non-backward evaluation in time
  assert( i_evaluationPoint >= i_expansionPoint );

  // copy DOFs (scalar==1)
  memcpy( o_timeEvaluated, i_timeDerivatives, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real) );

  // initialize scalars in the taylor series expansion (1st derivative)
  real l_deltaT = i_evaluationPoint - i_expansionPoint;
  real l_scalar = 1.0;
 
  // evaluate taylor series expansion at evvaluation point
  for( unsigned int l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_scalar *= -l_deltaT;
    l_scalar /= (real) l_derivative;

    // update the time evaluated taylor series by the contribution of this derivative
    for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[l_derivative]; l_basisFunction++ ) {
        // compute the different indices of the degrees of freedom due to the compressed storage scheme
        unsigned int l_timeEvaluatedIndex  =                                      l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS             + l_basisFunction;
        unsigned int l_timeDerivativeIndex = m_derivativesOffsets[l_derivative] + l_quantity*m_numberOfAlignedBasisFunctions[l_derivative] + l_basisFunction;

        o_timeEvaluated[l_timeEvaluatedIndex] += l_scalar * i_timeDerivatives[l_timeDerivativeIndex];
      }
    }
  }

}

void seissol::kernels::Time::computeIntegral(       real   i_expansionPoint,
                                                    real   i_integrationStart,
                                                    real   i_integrationEnd,
                                              const real*  i_timeDerivatives,
                                                    real*  o_timeIntegrated ) {
  // assert alignments
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart > i_expansionPoint   );
  assert( i_integrationEnd   > i_integrationStart );

  // reset time integrated degrees of freedom
  memset( o_timeIntegrated, 0, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real) );

  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real)  1;
  real l_secondTerm = (real)  1;
  real l_factorial  = (real) -1;
  real l_scalar;
 
  // iterate over time derivatives
  for(int l_derivative = 0; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    l_firstTerm  *=  l_deltaTUpper;
    l_secondTerm *=  l_deltaTLower;
    l_factorial  *= -(real)(l_derivative+1);

    l_scalar  = l_firstTerm - l_secondTerm;
    l_scalar /= l_factorial;

    // update the time integrated degrees of freedom by the contribution of this derivative
    for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
      for( unsigned int l_basisFunction = 0; l_basisFunction < m_numberOfAlignedBasisFunctions[l_derivative]; l_basisFunction++ ) {
        // compute the different indices of the degrees of freedom due to the compressed storage scheme
        unsigned int l_timeIntegratedIndex =                                      l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS             + l_basisFunction;
        unsigned int l_timeDerivativeIndex = m_derivativesOffsets[l_derivative] + l_quantity*m_numberOfAlignedBasisFunctions[l_derivative] + l_basisFunction;

        o_timeIntegrated[l_timeIntegratedIndex] += l_scalar * i_timeDerivatives[l_timeDerivativeIndex];
      }
    }
  }

}
