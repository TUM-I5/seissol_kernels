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

#ifndef TIME_H_
#define TIME_H_

#include <typedefs.hpp>
#include <cassert>
#include <common.hpp>

namespace seissol {
  namespace kernels {
    class Time;
  }
}

/**
 * Time kernel, which computes the time derivatives, intgrals and extrapolation for an element using the ADER scheme.
 *
 * All functions operate on a compressed memory format, which stores only the (memory aligned) non-zeros of the time derivatives.
 *
 * Numerical motivation: The size of the non-zero time differentiated unknowns reduces in every step.
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
 **/
class seissol::kernels::Time {
  // explicit private for unit tests
  private:
    //! aligned number of basis functions in decreasing order.
    unsigned int m_numberOfAlignedBasisFunctions[CONVERGENCE_ORDER];

    /*
     *! Offsets of the derivatives.
     *
     * * Offset counting starts at the zeroth derivative with o_derivativesOffset[0]=0; increasing derivatives follow:
     *   1st derivative: o_derivativesOffset[1]
     *   2nd derivative: o_derivativesOffset[2]
     *   ...
     * * Offset are always counted from positition zero; for example the sixth derivative will include all jumps over prior derivatives 0 to 5.
     */
    unsigned int m_derivativesOffsets[CONVERGENCE_ORDER];

    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global transposed stiffness matrix (case a) or B a star matrix (case b).
     * Each dense kernel (TODO: sparse implementation) has hardcoded BLAS-specifiers (M, N, K, ld(A), ld(B), ld(C), beta) exploiting the recursive structure.
     * The kernels are ordered as follows:
     *    0:         1st derivative \f$ M^{-1} ( K^\xi )^T \vee M^{-1} ( K^\eta )^T \vee M^{-1} ( K^\zeta )^T \f$
     *    1:         1st derivative \f$ A^* \vee B^* \vee C^* \f
     *    2:         2nd derivative \f$ M^{-1} ( K^\xi )^T \vee M^{-1} ( K^\eta )^T \vee M^{-1} ( K^\zeta )^T \f$
     *    3:         2nd derivative \f$ A^* \vee B^* \vee C^* \f
     *    ...
     *    2*(O-2):   O-1th derivative
     *    2*(O-2)+1: O-1th derivative
     *
     * Remark: The mass matrix \f$ M \f$ is diagonal.
     * 
     * The matrix kernels might prefetch (TODO: not implemented!) matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * In the case of the time integrator this means prefetching the transposed stiffness matrices (multiplied by the inverse mass matrices)
     * or star matrices of the upcoming operations in the recursive computation of the time derivatives.
     * The last operation,
     * \f[
     *   \left( \frac{\partial^{j_{\max-1}}}{\partial t^{j_{\max-1}} Q_k \right) C_k^*
     * \f],
     * already prefetches the stiffness, unknowns and time integrated unknowns matrix of the upcoming time integration of the next element.
     *
     * @param i_A left/transposed stiffness matrix (case a) or derivatives matrix (case b).
     * @param i_B right/derivatives matrix (case a) or star matrix (case b).
     * @param io_C resulting matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[(CONVERGENCE_ORDER-1)*2])( real *i_A,         real *i_B,         real *io_C,
                                                      real *i_APrefetch, real *i_BPrefetch, real *i_CPrefetch );

  public:
    /**
     * Gets the total size of the time derivative for a specific convergence order and alignment.
     *
     * @param i_convergenceOrder convergence order of the ADER-DG method.
     * @param i_numberOfQuantities number of quantities.
     * @param i_alignment memory alignment in bits.
     *
     **/
    static unsigned int getAlignedTimeDerivativesSize( unsigned int i_convergenceOrder   = CONVERGENCE_ORDER,
                                                       unsigned int i_numberOfQuantities = NUMBER_OF_QUANTITIES,
                                                       unsigned int i_alignment          = ALIGNMENT ) {
      assert( i_convergenceOrder >= 1 );
      unsigned int l_alignedDerivativesSize = 0;

      for( unsigned int l_associatedOrder = i_convergenceOrder; l_associatedOrder > 0; l_associatedOrder-- ) {
        l_alignedDerivativesSize += getNumberOfAlignedBasisFunctions( l_associatedOrder, i_alignment ) * i_numberOfQuantities;
      }

      return l_alignedDerivativesSize;
    }

    /**
     * Constructor, which initializes the time kernel.
     **/
    Time();

    /**
     * Computes the time derivatives.
     *   Storage format is compressed (storing only non-zeros).
     *   Use computeTimeIntegral to compute time integrated degrees of freedom from the derivatives.
     *   Use computeTimeEvaluation to evaluate the time prediction of the degrees of freedom based on the derivatives.
     *
     * @param i_stiffnessMatrices stiffness matrices, 0: \f$ K^\xi \f$, 1: \f$ K^\eta\f$, 2: \f K^\zeta \f$.
     * @param i_degreesOfFreedom of the current time step \f$ t^\text{cell} \f$ for which the time derivatives \f$ \frac{\partial^j}{\partial t^j} \f$ will be computed.
     * @param i_starMatrices star matrices, 0: \f$ A^*_k \f$, 1: \f$ B^*_k \f$, 2: \f$ C^*_k \f$.
     * @param o_timeDerivatives time derivatives of the degrees of freedom in compressed format.
     **/
    void computeDerivatives(       real** i_stiffnessMatrices,
                             const real*  i_degreesOfFreedom,
                                   real** i_starMatrices,
                                   real*  o_timeDerivatives );

    /**
     * Evaluates the taylor series expansion at the given evaluation point.
     *
     * @param i_expansionPoint expansion point of the taylor series (point in time where the time derivatives were computed from).
     * @param i_evaluationPoint point in time when the taylor series is evaluated.
     * @param i_timeDerivatives time derivatives.
     * @param o_degreesOfFreedom degrees of freedom at the time of the evaluation point.
     **/
    void computeExtrapolation(       real   i_expansionPoint,
                                     real   i_evaluationPoint,
                               const real*  i_timeDerivatives,
                                     real*  o_degreesOfFreedom );

    /**
     * Computes the time integrated degrees of freedom from previously computed time derivatives.
     *
     * @param i_expansionPoint expansion point (in time) of the Taylor series.
     * @param i_integrationStart start of the integration interval.
     * @param i_integrationEnd end of the integration interval.
     * @param i_timeDerivatives time derivatives.
     * @param o_timeIntegratedDofs time integrated DOFs over the interval: \f$ [ t^\text{start},  t^\text{end} ] \f$ 
     */
    void computeIntegral(       real   i_expansionPoint,
                                real   i_integrationStart,
                                real   i_integrationEnd,
                          const real*  i_timeDerivatives,
                                real*  o_timeIntegratedDofs );

};

#endif

