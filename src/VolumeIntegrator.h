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

#ifndef VOLUMEINTEGRATOR_H_
#define VOLUMEINTEGRATOR_H_

#include <Initializer/MemoryManager.h>
#include <Initializer/XmlParser.hpp>
#include <Initializer/preProcessorMacros.fpp>

namespace seissol {
  namespace kernels {
    class VolumeIntegrator;
  }
}

/**
 * Volume integrator, which computes the volume integral for an element.
 **/
class seissol::kernels::VolumeIntegrator {
  // explicit private for unit tests
  private:
    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global stiffness matrix (case a) or B the sparse star matrix (case b).
     * Each kernel can be sparse or dense, the kernels are ordered as follows:
     *    0:  \f$ K^\xi \f$
     *    1:  \f$ K^\eta \f$
     *    2:  \f$ K^\zeta f$
     *    3:  \f$ A^* \vee B^* \vee C^* \f
     *
     * Remark: The ordering of the global matrix kernels is identical to the ordering of the addresses.
     *
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * In the case of the volume integrator this means prefetching the stiffness matrices (multiplied by the inverse mass matrix)
     * or star matrices of the upcoming matrix operation.
     * The last multiplication,
     * \f[
     *   I(t^{n}, t^{n+1}, Q_{k}^n) C^*_k
     * \f],
     * already prefetches the unknowns, time integrated unknowns and first flux matrix of this elements upcoming boundary integration.
     * 
     * @param i_A left/stiffness matrix (case a) or unknowns matrix (case b).
     * @param i_B right/unknowns matrix (case a) or star matrix (case b).
     * @param io_C result matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[4])( double *i_A,         double *i_B,         double *io_C,
                                double *i_APrefetch, double *i_BPrefetch, double *i_CPrefetch );

    /**
     * Addresses of the stiffness matrices (multiplied by the inverse mass matrix).
     *    0:  \f$ M^{-1} K^\xi \f$
     *    1 : \f$ M^{-1} K^\eta \f$
     *    2:  \f$ M^{-1} K^\zeta f$
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
     * Constructor, which initializes the volume integrator according to the matrix setup in the given XML-file.
     *
     * @param i_matrixReader XML matrix reader.
     * @param i_memoryManager memory manager of SeisSol.
     **/
    VolumeIntegrator( const seissol::XmlParser                   &i_matrixReader,
                      const seissol::initializers::MemoryManager &i_memoryManager );

    /**
     * Computes the volume integral.
     *
     * @param i_timeIntegratedUnknowns time integrated unknowns of the element \f$ k \f$ over the interval \f$ [t^n, t^{n+1}] \f$:
     *        \f$ I(t^n, t^{n+1}, Q_k^n) \f$.
     * @param i_stiffnessMatrices TODO: see above
     * @param i_aStar sparse star matrix \f$ A^*_k \f$
     * @param i_bStar sparse star matrix \f$ B^*_k \f$
     * @param i_cStar sparse star matrix \f$ C^*_k \f$
     * @param io_unknowns unknowns of the previous time step \f$ t^n \f$, which will be updated by the volume integral:
     *        \f[
     *            M^{-1} K^\xi   I(t^{n}, t^{n+1}, Q_{k}^n) A^*_k
     *          + M^{-1} K^\eta  I(t^{n}, t^{n+1}, Q_{k}^n) B^*_k
     *          + M^{-1} K^\zeta I(t^{n}, t^{n+1}, Q_{k}^n) C^*_k
     *        \f]
     *
     **/
    void computeVolumeIntegral( double  i_timeIntegratedUnknowns[NUMBEROFUNKNOWNS],
                                double *i_stiffnessMatrices[3],
                                double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                double  i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                double  io_unknowns[NUMBEROFUNKNOWNS] );

    /**
     * DEPRECATED: Fall back code, which uses internal stiffness matrices.
     **/
    void computeVolumeIntegral( double i_timeIntegratedUnknowns[NUMBEROFUNKNOWNS],
                                double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                double i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                double io_unknowns[NUMBEROFUNKNOWNS] ) {
      computeVolumeIntegral( i_timeIntegratedUnknowns,
                             m_stiffnessMatrixPointers,
                             i_aStar,
                             i_bStar,
                             i_cStar,
                             io_unknowns );
    }

};

#endif

