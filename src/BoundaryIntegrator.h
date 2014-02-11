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
 * Copyright (c) 2013-2014
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
 * Boundary integration in SeisSol.
 **/

#ifndef BOUNDARYINTEGRATOR_H_
#define BOUNDARYINTEGRATOR_H_

#ifndef NDEBUG
#warning compiling boundary integrator with assertions
#endif

#include <Monitoring/FlopCounter.hpp>

#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <Initializer/XmlParser.hpp>
#include <Initializer/MemoryManager.h>
#include <generated_code/matrix_kernels/flux_matrices_3d.hpp_include>
#include <generated_code/matrix_kernels/dense_matrices.hpp_include>
#include <Initializer/preProcessorMacros.fpp>

#include <utils/logger.h>

namespace seissol {
  namespace kernels {
    class BoundaryIntegrator;
  }
}

/**
 * Boundary integrator, which computes the fluxes over the four tetrahedral faces.
 **/
class seissol::kernels::BoundaryIntegrator {
  // explicit private for unit tests
  private:
    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global flux matrix (case a) or B the flux solver (case b, 52).
     *   Each kernel can be sparse or dense.
     *   The kernels are ordered element local contribution \f$ F^{-, i} \f$ in front, which is 
     *   followed by the flux matrices for the neighbor element contribution \f$ F^{+, i, j, h} \f$
     *    0:  \f$ F^{-, 1} \f$
     *    1 : \f$ F^{-, 2} \f$
     *    2:  \f$ F^{-, 3} \f$
     *    3 : \f$ F^{-, 4} \f$
     *    4:  \f$ F^+{+, 1, 1, 1} \f$
     *    5:  \f$ F^+{+, 1, 1, 2} \f$
     *    6:  \f$ F^+{+, 1, 1, 3} \f$
     *    7:  \f$ F^+{+, 1, 2, 1} \f$
     *    8:  \f$ F^+{+, 1, 1, 2} \f$
     *    9:  \f$ F^+{+, 1, 1, 3} \f$
     *    [..]
     *    51: \f$ F^+{+, 4, 4, 3} \f$
     *    52: \f$ N_{k,i} A_k^+ N_{k,i}^{-1}\f$ or \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
     *
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * In the case of the boundary integrator this means prefetching flux matrices, time integrated unknowns
     * or flux solvers of the upcoming operations for the first 7 summands.
     * The last operation,
     * \f[
     *     I(t^{n}, t^{n+1}, Q_{k(4)}^n) N_{k,4} A_{k(4)}^- N_{k,4}^{-1}
     * \f],
     * already prefetches the stiffness matrix, unknowns and time integrated unknowns of the upcoming next elements volume integration.
     * 
     *
     * @param i_A left/flux matrix (case a) or unknowns matrix (case b).
     * @param i_B right/unknowns matrix (case a) or flux solver (case b).
     * @param i_C result matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[53])( double *i_A,         double *i_B,         double *io_C,
                                 double *i_APrefetch, double *i_BPrefetch, double *i_CPrefetch );
    
    /**
     * Addresses of the global flux matrices (divided by the inverse diagonal mass matrix).
     *   The order of the pointers (not of the chunks in memory!) is identical to the order of the
     * kernels.
     **/ 
    double** m_fluxMatrixPointers;

    /**
     * Sets the function pointer for the matrix kernel of the specified matrix id, which can be handled sparse or dense.
     *
     * @param i_id id of the matrix.
     * @param i_sparse true if the matrix is handled sparse, false if dense.
     **/
    void setUpMatrixKernel( unsigned int i_id,
                            bool i_sparse );

  public:
    /**
     * Constructor, which initializes the solver according to the matrix setup in the given XML-file.
     *
     * @param i_matrixReader XML matrix reader.
     * @param i_memoryManager memory manager of SeisSol.
     **/
    BoundaryIntegrator( const seissol::XmlParser                   &i_matrixReader,
                        const seissol::initializers::MemoryManager &i_memoryManager );

    /**
     * Computes the boundary integral of a single element.
     *
     * @param i_timeIntegratedUnknownsElement time integrated unknowns of the current element \f$k\f$ and next element \f$k+1\f$ over the inteval \f$ [t^n, t^{n+1}] \f$:
     *        \f$ I(t^n, t^{n+1}, Q_k^n) \f$ and \f$ I(t^n, t^{n+1}, Q_{k+1}^n) \f$.
     * @param i_timeIntegratedUnknownsNeighbors time integrated unknowns of the neighbors \f$k(i) \f$ over the inteval \f$ [t^n, t^{n+1}] \f$
     *        \f$ I(t^n, t^{n+1}, Q_k(i)^n) \f$.
     * @param i_boundaryConditions boundary conditions of the face. Absorbing boundaries (5) change the execution of the kernel.
     * @param i_neighboringIndices indices \f$j(i) \in \{0,1,2,3\}\f$ and \f$h(i) \in \{0,1,2\}\f$,
     *        which depend on the face combinations of the current elements faces \f$ i \in
     *        \{0,1,2,3\}\f$ and the neighboring faces \f$ j(i) \f$ and vertex orientation of the
     *        element and neighboring face \f$ h(i) \f$.
     *        A two dimensional array is expected, where the first index i_neighboringIndices[i][]
     *        selects the face of the current element and the second index i_neighboringIndices[][*]
     *        gives \f$ j(i)\f$: i_neighboringIndices[][0] or \f$ h(i) \f$:  i_neighboringIndices[][1].
     * @param i_fluxMatrices TODO: see above
     * @param i_nApNm1 matrices \f$N_{k,i} A_k^+ N_{k,i}^{-1}\f$,
     *        where \f$N_{k,i}^{-1}\f$ and \f$N_{k,i}\f$ account for the forth and back transformation
     *        to reference space relative to the outerpointing normal of face \f$i\f$ in element
     *        \f$k\f$. \f$A_k^+\f$ is the flux contribution of element \f$k\f$, thus assembled using
     *        the positive eigenvalues of element \f$k\f$.
     * @param i_nAmNm1 matrices \f$N_{k,i} A_{k(i)}^- N_{k,i}^{-1}\f$,
     *        where \f$N_{k,i}^{-1}\f$ and \f$N_{k,i}\f$ account for the forth and back
     *        transformation to reference space relative to the outerpointing normal of face $i$ in
     *        element \f$k\f$. \f$A_{k(i)}^-\f$ is the flux contribuation of the \f$i\f$-th face
     *        neighbor of element \f$k\f$, thus assembled using negative eigenvalues of element \f$k(i)\f$.
     * @param io_unknowns unknowns at the previous time step \f$ t^n \f$ (eventually already updated with
     *        the volume contribution), which will be updated by the boundary integral:
     *        \f[
     *          \frac{|S_k|}{|j_k|} M^{-1}
     *          \bigg(
     *            \sum_{i=1}^4
     *              F^{-,i}
     *              I(t^{n}, t^{n+1}, Q_{k}^n)
     *              N_{k,i} A_k^+ N_{k,i}^{-1}
     *            \\
     *            &\hspace{35pt}+\sum_{i=1}^4
     *              F^{+,i,j,h}
     *              I(t^{n}, t^{n+1}, Q_{k(i)}^n)
     *              N_{k,i} A_{k(i)}^- N_{k,i}^{-1}
     *            \bigg)
     *        \f]
     **/
    void computeBoundaryIntegral(       double  i_timeIntegratedUnknownsElement[2][NUMBEROFUNKNOWNS],
                                        double *i_timeIntegratedUnknownsNeighbors[4],
                                  const int     i_boundaryConditions[4],
                                  const int     i_neighboringIndices[4][2],
                                        double *i_fluxMatrices[52],
                                        double  i_nApNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                        double  i_nAmNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                        double  io_unknowns[NUMBEROFUNKNOWNS] );

    /**
     * DEPRECATED: Fall back code, which uses internal flux matrices
     **/
    void computeBoundaryIntegral(       double  i_timeIntegratedUnknownsElement[2][NUMBEROFUNKNOWNS],
                                        double *i_timeIntegratedUnknownsNeighbors[4],
                                  const int     i_boundaryConditions[4],
                                  const int     i_neighboringIndices[4][2],
                                        double  i_nApNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                        double  i_nAmNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                        double  io_unknowns[NUMBEROFUNKNOWNS] ) {
      computeBoundaryIntegral( i_timeIntegratedUnknownsElement,
                               i_timeIntegratedUnknownsNeighbors,
                               i_boundaryConditions,
                               i_neighboringIndices,
                               m_fluxMatrixPointers,
                               i_nApNm1,
                               i_nAmNm1,
                               io_unknowns );
    };
    
};

#endif
