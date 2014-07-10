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
 * Time integration in SeisSol.
 **/

#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_

#include <Initializer/MemoryManager.h>
#include <Initializer/XmlParser.hpp>
#include <Initializer/preProcessorMacros.fpp>

namespace seissol {
  namespace kernels {
    class TimeIntegrator;
  }
}

/**
 * Time integrator, which computes the time integral for an element using the ADER scheme.
 **/
class seissol::kernels::TimeIntegrator {
  // explicit private for unit tests
  private:
    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global transposed stiffness matrix (case a) or B a sparse star matrix (case b).
     * Each kernel can be sparse or dense, the kernels are ordered as follows:
     *    0:  \f$ M^{-1} ( K^\xi )^T \f$
     *    1:  \f$ M^{-1} ( K^\eta )^T \f$
     *    2:  \f$ M^{-1} ( K^\zeta )^T f$
     *    3:  \f$ A^* \vee B^* \vee C^* \f
     *
     * Remark: The stiffness matrix \f$ M \f$ is diagonal.
     *         The ordering of the global matrix kernels is identical to the ordering of the addresses.
     * 
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * In the case of the time integrator this means prefetching the transposed stiffness matrices (multiplied by the inverse mass matrices)
     * or star matrices of the upcoming operations in the recursive computation of the time derivatives.
     * The last operation,
     * \f[
     *   \left( \frac{\partial^{j_{\max-1}}}{\partial t^{j_{\max-1}} Q_k \right) C_k^*
     * \f],
     * already prefetches the stiffness, unknowns and time integrated unknowns matrix of the upcoming time integration of the next element.
     *
     * @param i_A left/transposed stiffness matrix (case a) or unknowns matrix (case b).
     * @param i_B right/unknowns matrix (case a) or star matrix (case b).
     * @param io_C result matrix.
     * @param i_nonZeroBlockSize size of the non-zero blocks (see function documentation of computeBoundaryIntegal)
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[4])( double *i_A,         double *i_B,         double *io_C,       int i_nonZeroBlockSize,
                                double *i_APrefetch, double *i_BPrefetch, double *i_CPrefetch                         );

    /**
     * Addresses of the transposed stiffness matrices (multiplied by the inverse mass matrix).
     *    0:  \f$ M^{-1} ( K^\xi )^T \f$
     *    1:  \f$ M^{-1} ( K^\eta )^T \f$
     *    2:  \f$ M^{-1} ( K^\zeta )^T f$
     *
     *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks.
     **/ 
    double** m_stiffnessMatrixPointers;

    /**
     * Sets the function pointer for the matrix kernel of the specified matrix, which can be handled sparse or dense.
     *
     * @param i_id id of the matrix.
     * @param i_sparse true if the matrix is handled sparse, false if dense.
     **/
    void setUpMatrixKernel( unsigned int i_id,
                            bool i_sparse );

  public:
    /**
     * Constructor, which initializes the time integrator according to the matrix setup in the given XML-file.
     *
     * @param i_matrixReader XML matrix reader.
     * @param i_memoryManager memory manager of SeisSol.
     **/
    TimeIntegrator( const seissol::XmlParser                   &i_matrixReader,
                    const seissol::initializers::MemoryManager &i_memoryManager );

    /**
     * Computes the time derivatives.
     *   Part of the local time stepping scheme.
     *   Use the function computeTimeIntegral to compute integrals from derivatives.
     *
     * TODO: Switch to compressed storage scheme.
     *
     * @param i_unknowns unknowns of the current time step \f$ t^\text{cell} \f$ for which the time derivatives \f$ \frac{\partial^j}{\partial t^j} \f$ will be computed.
     * @param i_stiffnessMatrices TODO: see description above.
     * @param i_aStar sparse star matrix \f$ A^*_k \f$
     * @param i_bStar sparse star matrix \f$ B^*_k \f$
     * @param i_cStar sparse star matrix \f$ C^*_k \f$
     * @param o_timeDerivatives time derivatives used in the time integration.
     **/
    void computeTimeDerivatives( const double i_unknowns[NUMBEROFUNKNOWNS],
                                       double *i_stiffnessMatrices[3],
                                       double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                       double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                       double i_cStar[STARMATRIX_NUMBEROFNONZEROS], 
                                       double o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] );

    /*
     * DEPRECATED: Fall back code, which uses internal flux matrices.
     */    
    void computeTimeDerivatives( const double i_unknowns[NUMBEROFUNKNOWNS],
                                       double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                       double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                       double i_cStar[STARMATRIX_NUMBEROFNONZEROS], 
                                       double o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] ) {
      computeTimeDerivatives( i_unknowns,
                              m_stiffnessMatrixPointers,
                              i_aStar,
                              i_bStar,
                              i_cStar, 
                              o_timeDerivatives );
     }

    /**
     * Evaluates the taylor series expansion at the given evaluation point.
     *
     * @param i_timeDerivatives time derivatives.
     * @param i_expansionPoint expansion point of the taylor series (point in time when the time derivatives were computed from).
     * @param i_evaluationPoint point in time when the taylor series is evaluated.
     * @param o_unknowns unknowns at the time of the evaluation point.
     **/
    void computeTimeEvaluation( const double  i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                const double &i_expansionPoint,
                                const double &i_evaluationPoint,
                                      double  o_unknowns[NUMBEROFUNKNOWNS] );

    /**
     * Computes the time integrated unknowns from previously computed time derivatives.
     *   Part of the local time stepping scheme.
     *
     * TODO: Switch to compressed storage scheme.
     *
     * @param i_timeDerivatives time derivatives.
     * @param i_deltaTLower length \f$ \Delta t^\text{lo} \f$ of the lower time integration interval with respect to the current time \f$ t^\text{cell} \f$ of the cell.
     * @param i_deltaTUpper length \f$ \Delta t^\text{up} \f$ of the upper time integration interval with respect to the current time \f$ t^\text{cell} \f$ of the cell.
     * @param o_timeIntegratedUnknowns time integrated unknowns over the interval: \f$ [ t^\text{cell} + \Delta t^\text{lo}, t^\text{cell} + \Delta t^\text{up} ] \f$ 
     */
    void computeTimeIntegral( const double  i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                              const double &i_deltaTLower,
                              const double &i_deltaTUpper,
                                    double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] );

    /**
     * Computes the time integral.
     *   Part of the global time stepping scheme.
     *
     * @param i_unknowns unknowns of the current time step \f$ t^n \f$, which will be integrated in time to \f$ t^{n+1} = t^n + \Delta t \f$
     * @param i_stiffnessMatrices TODO: see description above.
     * @param i_aStar sparse star matrix \f$ A^*_k \f$
     * @param i_bStar sparse star matrix \f$ B^*_k \f$
     * @param i_cStar sparse star matrix \f$ C^*_k \f$
     * @param i_deltaT length \f$ \Delta t \f$ of the time integration interval \f$ [t^n, t^{n+1}] \f$
     * @param o_timeIntegratedUnknowns time integrated unknowns:
     *        \f[
     *          I(t^{n}, t^{n+1}, Q_k^n) = \sum_{j=0}^{\mathcal{O}-1} \frac{(t^{n+1}-t^n)^{j+1}}{(j+1)!} \frac{\partial^{j}}{\partial t^j} Q_k(t^n)
     *        \f]
     *        , where the time derivatives \f$ \frac{\partial^{j}}{\partial t^j} Q_k(t^n) \f$ are computed by the following recursive scheme:
     *        \f[
     *          \frac{\partial^{j+1}}{\partial t^{j+1}} Q_k =
     *          - M^{-1} \left( (K^\xi)^T   \left( \frac{\partial^{j}}{\partial t^j} Q_k \right) A_k^* +
     *                          (K^\eta)^T  \left( \frac{\partial^{j}}{\partial t^j} Q_k \right) B_k^* +
     *                          (K^\zeta)^T \left( \frac{\partial^{j}}{\partial t^j} Q_k \right) C_k^* \right)
     *        \f]
     **/
    void computeTimeIntegral( const double  i_unknowns[NUMBEROFUNKNOWNS],
                                    double *i_stiffnessMatrices[3],
                                    double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                    double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                    double  i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                              const double &i_deltaT,
                                    double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] );

    /**
     * DEPRECATED: Fall back code, which uses internal flux matrices.
     **/
    void computeTimeIntegral( const double  i_unknowns[NUMBEROFUNKNOWNS],
                                    double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                    double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                    double  i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                              const double &i_deltaT,
                                    double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
      computeTimeIntegral( i_unknowns,
                           m_stiffnessMatrixPointers,
                           i_aStar,
                           i_bStar,
                           i_cStar,
                           i_deltaT,
                           o_timeIntegratedUnknowns );
    }

};

#endif

