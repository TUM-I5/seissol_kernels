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
 * Test suite, which tests the generated sparse matrix kernels, which operatare on subblocks of the matrix.
 */
#include <cstdlib>
#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"
#include <Initializer/MemoryAllocator.h>
#include "DenseMatrix.hpp"
#include "generated_code/matrix_kernels/stiffness_matrices_3d.hpp_include"
#include "generated_code/matrix_kernels/star_matrices_3d.hpp_include"

namespace unit_test {
  class BlockedMatrixKernelsTestSuite;
}

class unit_test::BlockedMatrixKernelsTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! aligned memory allocation
    seissol::MemoryAllocator m_memoryAllocator;

    //! path to the matrice directory
    std::string m_matricesDirectory;

    //! dense matrix functionality
    DenseMatrix m_denseMatrix;

    /**
     * Computes the number of basis functions from a given degree of the polynomial basis.
     * @param i_basisDegree degree of the polynomial basis.
     * @return number of basis functions.
     **/
    int computeNumberOfBasisFunctions( int i_basisDegree ) {
      int l_numberOfBasisFunctions  = i_basisDegree;
      l_numberOfBasisFunctions  = (l_numberOfBasisFunctions + 1) * (l_numberOfBasisFunctions + 2) * (l_numberOfBasisFunctions + 3);
      l_numberOfBasisFunctions /= 6;

      return l_numberOfBasisFunctions;
    }

    /**
     * Tests the sparse matrices appearing in the time integration for the given degree.
     * @i_basisDegree polynomial degree of the basis.
     **/
    void testTimeIntegrationMatrices( int i_basisDegree ) {
      // set number of variables, TODO: attenuation
      int l_numberOfVariables = 9;

      // compute the number of basis functions
      int l_numberOfBasisFunctions  = computeNumberOfBasisFunctions( i_basisDegree );

      // setup path to xml file
      std::string l_matricesPath = m_matricesDirectory + "matrices_" + std::to_string(l_numberOfBasisFunctions) + ".xml";

      // setup the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );

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

      // read the stiffness matrices
      l_matrixReader.readGlobalMatrices( "stiffness",
                                         l_matrixIds,  l_matrixNames,
                                         l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                         l_matrixRows, l_matrixColumns, l_matrixValues );

      // itearte over all stiffness matrices
      for( int l_matrixIndex = 0; l_matrixIndex < l_matrixNames.size(); l_matrixIndex++ ) {
        // iterate over all orders in the taylor series expansion
        for( int l_degree = i_basisDegree; l_degree > 0; l_degree-- ) {
          // get the non-zero block size of this order
          int l_nonZeroBlockSize = computeNumberOfBasisFunctions( l_degree );

          // get the current matrix name
          std::string l_matrixName = l_matrixNames[l_matrixIndex];

          // continue with transposed stiffness matrices only, which are used in the ADER time integration.
          if( l_matrixName != "kXiT"  &&
              l_matrixName != "kEtaT" &&
              l_matrixName != "kZetaT"   ) {
            continue;
          }

          // matrices for operation C = A.B 
          double *l_a = NULL, *l_b = NULL, *l_c1 = NULL, *l_c2 = NULL;

          // allocate memory and set to random values
          // For simplicity we are reusing dense memory l_a for the sparse matrix
          m_denseMatrix.allocateMemoryAndSetRandomValues(  l_numberOfBasisFunctions,
                                                           l_numberOfVariables,
                                                           l_numberOfBasisFunctions,
                                                           &l_a, &l_b, &l_c1, &l_c2 );

          // copy sparse matrix values to aligned memory
          for( int l_matrixValueIndex = 0; l_matrixValueIndex < l_matrixValues[l_matrixIndex].size(); l_matrixValueIndex++ ) {
            l_a[l_matrixValueIndex] = l_matrixValues[l_matrixIndex][l_matrixValueIndex];
          }

          // do the blocked sparse matrix multiplication
#include "generated_code/unit_tests/time_sparse_matrix_kernels.hpp_include"

          // reset matrix to zero and set dense matrix values
          std::fill( l_a, l_a+l_numberOfBasisFunctions*l_numberOfBasisFunctions, 0 );
          for( int l_matrixValueIndex = 0; l_matrixValueIndex < l_matrixValues[l_matrixIndex].size(); l_matrixValueIndex++ ) {
            // compute memory index, remark: row- and column-counting in the XML-files starts at 1
            unsigned int l_nonZeroIndex = (l_matrixRows[l_matrixIndex][l_matrixValueIndex]   - 1) +
                                          (l_matrixColumns[l_matrixIndex][l_matrixValueIndex] -1) * l_matrixNumberOfRows[l_matrixIndex];
            l_a[l_nonZeroIndex] = l_matrixValues[l_matrixIndex][l_matrixValueIndex];
          }

          // set zero block
          m_denseMatrix.setZeroBlock( l_numberOfBasisFunctions,
                                      l_numberOfVariables,
                                      l_nonZeroBlockSize,
                                      l_b );

          // do the dense multiplication
          m_denseMatrix.executeStandardMultiplication( l_numberOfBasisFunctions,
                                                       l_numberOfVariables,
                                                       l_numberOfBasisFunctions,
                                                       l_a, l_b, l_c2 );
        
          // check the result
          m_denseMatrix.checkResult( l_numberOfBasisFunctions*l_numberOfVariables,
                                     l_c1, l_c2 );
        }
      }
      // TODO: implement
#warning Time integration unit test for star matrices missing.

    }

  public:
    void tearUp() {
      // get matrices directory
      m_matricesDirectory = m_configuration.getMatricesDirectory();
    }

    void tearDown() {
      // free allocated memory
      m_memoryAllocator.freeMemory();
    }

    void testSparseDenseKernels() {
      // test matrices of all degrees
      for( int l_basisDegree = 1; l_basisDegree < 6; l_basisDegree++ ) {
        testTimeIntegrationMatrices( l_basisDegree );
      }
    }
};
