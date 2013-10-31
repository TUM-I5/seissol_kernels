/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2012
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
 * Benchmarking routine for the Cauchy Kovalewski time integration.
 **/
 
#include "seissol_src/Initializer/preProcessorMacros.fpp"
#include "BenchmarkData.hpp"

// preprocessor helping function to convert integer constants to strings
#define CONCAT_HELPER(a,b,c,d) a ## b ## c ## d
#define CONCAT(a,b,c,d) CONCAT_HELPER(a,b,c,d)

// create corresponding C-function name (-> MatrixGen)
#define CKTIMEINTFUNCTIONNAME CONCAT(cauchyKovalewskiTimeIntegration_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS)

// number of iterations over the grid
#define NUMBEROFBENCHMARKITERATIONS 100

// generated Cauchy Kovalewski routine
extern "C" void CKTIMEINTFUNCTIONNAME( const double i_unknowns[NUMBEROFUNKNOWNS],
                                             double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                             double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                             double i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                             double i_kXiDivMT[KXI_NUMBEROFNONZEROS],
                                             double i_kEtaDivMT[KETA_NUMBEROFNONZEROS],
                                             double i_kZetaDivMT[KZETA_NUMBEROFNONZEROS],
                                       const double &i_deltaT,
                                             double o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] );

/**
 * Computes the Cauchy Kovalewski time integration for the first i_problemSize number of elements
 * in the pseude-discretization (see BenchmarkData.hpp for details regarding the discrization).
 *
 * The computation is repeated NUMBEROFBENCHMARKITERATIONS times.
 *   Remark: In contrast to computations in SeisSol no work is done in between -
 *           After the last element at position [i_problemSize-1] is updated
 *           the first element at position [0] follows directly. This increases the cache reuse.
 *
 * @param i_problemSize number of elements for which the Cauchy Kovalewski time integration is computed.
 * @param io_benchmarkData benchmark data used: $Q_k$, $A^*_k$, $B^*_k$, $C^*_k$,
 *                                                     $K_\xi^T * M^{-1}$, $K_\eta^T * M^{-1}$, $K_\zeta^T * M^{-1}$,
 *                                                     $\Delta t$,
 *                                              $\int_{t^n}^{t^{n+1}} Q_k \d t$
 */
inline void cauchyKovalewskiTimeIntegrationKernel( const int i_problemSize,
                                                         BenchmarkData& io_BenchmarkData ) {
  // repeat the time integration NUMBEROFBENCHMARKITERATIONS times.
  for( int l_iteration = 0; l_iteration < NUMBEROFBENCHMARKITERATIONS; l_iteration++) {
    // iterate of the first i_problem size number of elemetn in the pseudo-discrization.
    for( int l_elementIndex = 0; l_elementIndex < i_problemSize; l_elementIndex++ ) {
      // call generated CCKkernel.
      CKTIMEINTFUNCTIONNAME( io_BenchmarkData.discontinousGalerkinElements[l_elementIndex].unknowns,
                             io_BenchmarkData.discontinousGalerkinElements[l_elementIndex].aStarFlat,
                             io_BenchmarkData.discontinousGalerkinElements[l_elementIndex].bStarFlat,
                             io_BenchmarkData.discontinousGalerkinElements[l_elementIndex].cStarFlat,
                             io_BenchmarkData.kXiDivMTFlat,
                             io_BenchmarkData.kEtaDivMTFlat,
                             io_BenchmarkData.kZetaDivMTFlat,
                             io_BenchmarkData.deltaT,
                             io_BenchmarkData.discontinousGalerkinElements[l_elementIndex].timeIntegratedUnknowns );
    }
  }
}
