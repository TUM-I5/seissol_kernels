/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Test suite, which tests the generated dense matrix kernels.
 **/
#include <cstdlib>
#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"
#include <Initializer/MemoryAllocator.h>
#include <matrix_kernels/dense_matrices.hpp_include>
#include "DenseMatrix.hpp"

namespace unit_test {
  class DenseMatrixKernelTestSuite;
}

/**
 * Unit tests for the dense generated matrix kernels.
 *
 *   Sketch (C = A.B)
 *                    B
 *            ****************
 *            *     k x n    *
 *            ****************
 *      A             C
 *   *******  ****************
 *   *     *  *              *
 *   *  m  *  *              *
 *   *  x  *  *     m x n    *
 *   *  k  *  *              *
 *   *     *  *              *
 *   *******  ****************
 **/
class unit_test::DenseMatrixKernelTestSuite: public CxxTest::TestSuite {
  //private:
    DenseMatrix m_denseMatrix;

    /**
     * Executes C += A.B with the generated kernel.
     *
     * @param i_m #(dense rows) of the left matrix
     * @param i_n #(dense columns) of the left matrix
     * @param i_k #(dense columns) of the right matrix
     * @o_a pointer to matrix A.
     * @o_b pointer to matrix B. 
     * @o_c pointer to matrix C.
     **/
    void executeGeneratedMultiplication( int     i_m, int     i_n, int     i_k,
                                         double *i_a, double *i_b, double *io_c ) {
      #include <unit_tests/dense_matrix_kernels.hpp_include>
    }

    /**
     * Tests generated multiplication (C += A.B) against default execution.
     * @param i_m #(dense rows) of the left matrix
     * @param i_n #(dense columns) of the left matrix
     * @param i_k #(dense columns) of the right matrix
     * @param i_add true: C += A.B, false: C = A.B
     **/
    void testGeneratedMultiplication( int i_m, int i_n, int i_k, bool i_add = true ){
      // matrices for operation C = A.B 
      double *l_a, *l_b, *l_c1, *l_c2;

      m_denseMatrix.allocateMemoryAndSetRandomValues( i_m, i_n, i_k, &l_a, &l_b, &l_c1, &l_c2 );

      // do the generated multiplcation
      executeGeneratedMultiplication( i_m, i_n, i_k, l_a, l_b, l_c1 );

      // do the default multiplication
      m_denseMatrix.executeStandardMultiplication( i_m, i_n, i_k, l_a, l_b, l_c2, i_add );

      // check result
      m_denseMatrix.checkResult( i_m * i_n, l_c1, l_c2 );
    }


  public:
    void testDenseFluxMatrixKernels() {
      // generate random seed
      srand(time(NULL));

      // parameters in the matrix operation
      int l_m, l_n, l_k;      

      // iterate over orders
      for( int l_order = 2; l_order <= s_maximumOrder; l_order++) {
        // compute #(basis functions)
        int l_numberOfBasisFunctions = l_order * (l_order+1) * (l_order+2) / 6;

        // setup matrices for dense flux multiplication
        // TODO: add attenuation here
        l_m = l_numberOfBasisFunctions;
        l_n = 9;
        l_k = l_numberOfBasisFunctions;
        
        // test dense flux kernel
        testGeneratedMultiplication(l_m, l_n, l_k, false);
      }
    }

    void testDenseStarMatrixKernels() {
      // generate random seed
      srand(time(NULL));

      // parameters in the matrix operation
      int l_m, l_n, l_k;      

      // iterate over orders
      for( int l_order = 2; l_order <= s_maximumOrder; l_order++) {
        // compute #(basis functions)
        int l_numberOfBasisFunctions = l_order * (l_order+1) * (l_order+2) / 6;

        // setup matrices for dense star multiplication
        // TODO: add attenuation here
        l_m = l_numberOfBasisFunctions;
        l_n = 9;
        l_k = 9;
        
        // test dense star kernel
        testGeneratedMultiplication(l_m, l_n, l_k);
      }
    }

};
