/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Alexander Heinecke (heinecke AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Alexander_Heinecke,_M.Sc.,_M.Sc._with_honors)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
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
 * C++ representation of a sparse matrix.
 **/

#include <string>
#include <cstdlib>
#include <cxxtest/TestSuite.h>

#include "../../SeisSolGen/ReaderCSC.hpp"
#include "generated_code/matrix_kernels/star_matrices_3d.hpp_include"
#include "generated_code/matrix_kernels/stiffness_matrices_3d.hpp_include"
#include "generated_code/matrix_kernels/flux_matrices_3d.hpp_include"

namespace unit_test {
  class SparseMatrix;
}

class unit_test::SparseMatrix {
  public:
    enum MatrixType{
      StarMatrix,
      StiffnessMatrix,
      FluxMatrix,
      TransposedStiffnessMatrix
    };

    //! path to matrix market file, which contains the definition of the matrix.
    std::string pathToMatrixMarketFile;

    //! type of the sparse matrix
    MatrixType matrixType;

    //! array containing row ids of the non-zero elements.
    int* rowIndices;
    //! array containing column ids of the non-zero elements.
    int* columnIndices;
    //! array containing the non-zero elements.
    double* flatColumnMajor;

    //! number of rows.
    int numberOfRows;
    //! number of columns.
    int numberOfColumns;
    //! number of non-zero elements.
    int numberOfNonZeros;

    //! array containing the dense matrix (including zeros) in a colum-wise storage
    double* denseColumnMajor;

    //! numerical tolerance under which a multiplication failed
    double errorTolerance;

    /**
     * Deep copy operator (memory allocation + copy)
     * @param i_sparseMatrix sparse matrix, which will be copied
     */
    void deepCopy( const SparseMatrix &i_sparseMatrix ) {
      // copy variables
      pathToMatrixMarketFile = i_sparseMatrix.pathToMatrixMarketFile;
      matrixType = i_sparseMatrix.matrixType;
      errorTolerance = i_sparseMatrix.errorTolerance;

      // copy integers
      numberOfRows     = i_sparseMatrix.numberOfRows;
      numberOfColumns  = i_sparseMatrix.numberOfColumns;
      numberOfNonZeros = i_sparseMatrix.numberOfNonZeros;

      // allocate memory
      rowIndices              = (int*)    _mm_malloc( numberOfNonZeros*sizeof(int),     64 ); //TODO: use standard malloc
      columnIndices           = (int*)    _mm_malloc( numberOfNonZeros*sizeof(int),     64 ); //TODO: use standard malloc
      flatColumnMajor         = (double*) _mm_malloc( numberOfNonZeros*sizeof(double), 64 ); //TODO: use standard malloc
      denseColumnMajor        = (double*) malloc( numberOfRows*numberOfColumns*sizeof(double) );

      // copy arrays
      std::copy( i_sparseMatrix.rowIndices,      i_sparseMatrix.rowIndices       + numberOfNonZeros, rowIndices                     );
      std::copy( i_sparseMatrix.columnIndices,   i_sparseMatrix.columnIndices    + numberOfNonZeros, columnIndices                  );
      std::copy( i_sparseMatrix.flatColumnMajor, i_sparseMatrix.flatColumnMajor  + numberOfNonZeros, flatColumnMajor                );
      std::copy( i_sparseMatrix.denseColumnMajor,i_sparseMatrix.denseColumnMajor + (numberOfRows*numberOfColumns), denseColumnMajor );
    }

    /**
     * Copies the values of one array to another.
     *
     * @param i_array input array
     * @param numberOfElements number of elements for each array.
     * @param o_array output array
     */
    void copyValues( const double* i_array, int numberOfElements, double* o_array ) {
      for( int l_i = 0; l_i < numberOfElements; l_i++)
        o_array[l_i] = i_array[l_i];
    }

    /**
     * Sets random values to the given array.
     *
     * @param i_numberOfElements number of elements of the array.
     * @param o_array the array.
     */
    static void setRandomValues( int i_numberOfElements, double* o_array ) {
      // set random values
      for(int l_i = 0; l_i < i_numberOfElements; l_i++) {
        o_array[l_i] = ((double)rand()/(double)RAND_MAX)*10.0;
      }
    }

    /**
     * Sets random values to the sparse matrix entries.
     */
    void setRandomValues() {
      //set random values to the defined matrix
      setRandomValues( numberOfNonZeros, flatColumnMajor );
      //copy them to dense representation
      init_colmajor_from_csc( numberOfColumns, denseColumnMajor, numberOfRows,
                              flatColumnMajor, columnIndices, rowIndices);
    }

    /**
     * Creates a dense test matrix and fills it with random values
     *
     * @param i_numberOfRows number of rows.
     * @param i_numberOfColumns number of columns
     */
    void createDenseTestMatrix( int i_numberOfRows, int i_numberOfColumns, double*& o_testMatrix ) {
      // allocate memory
      o_testMatrix = (double*) malloc( i_numberOfRows*i_numberOfColumns*sizeof(double) );

      setRandomValues( (i_numberOfRows*i_numberOfColumns), o_testMatrix );
    }

    //TODO: redundant code -> BenchmarkCSC.cpp
    void init_colmajor_from_csc(int N, double* A, int LDA, double* values, int* colidx, int* rowidx)
    {
        for (int l = 0; l < (LDA*N); l++)
            A[l] = 0.0;

        for (int t = 0; t < N; t++)
        {
            int lcl_colElems = colidx[t+1] - colidx[t];
            for (int z = 0; z < lcl_colElems; z++)
            {
                A[(t*LDA)+rowidx[colidx[t] + z]] = values[colidx[t] + z];
            }
        }
    }

    //TODO: redundant code -> BenchmarkCSC.cpp
    void mul_colmajor(int M, int N, int K, double* A, int LDA, double* B, int LDB, double* C, int LDC)
    {
        for(int n = 0; n < N; n++)
        {
            for(int k = 0; k < K; k++)
            {
                #pragma simd
                for(int m = 0; m < M; m++)
                {
                    C[(LDC*n)+m] += A[(LDA*k)+m] * B[(LDB*n)+k];
                }
            }
        }
    }

    //TODO: redundant code -> Benchmark.hpp
    static void print_colmajor(double* ptr, size_t nRows, size_t nCols)
    {
        std::cout << std::endl;
        for (size_t l = 0; l < nRows; l++)
        {
            for (size_t k = 0; k < nCols; k++)
            {
                std::cout << ptr[(k*nRows)+l] << " ";
            }
           std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    //TODO: redundant code -> Benchmark.hpp
    static double compare_colmajor(double* ptrA, double* ptrB, size_t nElems)
    {
       double error = 0.0;
     for (size_t i = 0; i < nElems; i++)
     {
         if (fabs(ptrA[i]-ptrB[i]) > error)
           error = fabs(ptrA[i]-ptrB[i]) / fmin(ptrA[i], ptrB[i]);
     }
     return error;
    }

  public:
    // constructor
    SparseMatrix( const std::string &i_pathToMatrixMarketFile,
                  MatrixType i_matrixType,
                  double i_errorTolerance = 0.000001 ) {
      // save variables
      pathToMatrixMarketFile = i_pathToMatrixMarketFile;
      matrixType = i_matrixType;
      errorTolerance = i_errorTolerance;

      // Construct a CSC reader
      seissolgen::ReaderCSC l_readerCSC;
      l_readerCSC.parse_file( pathToMatrixMarketFile,
                              rowIndices,
                              columnIndices,
                              flatColumnMajor,
                              numberOfRows,
                              numberOfColumns,
                              numberOfNonZeros );

      // allocate memory for respective dense matrix
      denseColumnMajor = (double*) malloc( numberOfRows*numberOfColumns*sizeof(double) );

      // initialize the dense matrix
      init_colmajor_from_csc( numberOfColumns, denseColumnMajor, numberOfRows,
                              flatColumnMajor, columnIndices, rowIndices);
    }

    // Copy construtor
    SparseMatrix( const SparseMatrix &i_sparseMatrix ) {
      deepCopy( i_sparseMatrix );
    }

    // Copy assignment operator
    SparseMatrix& operator= (const SparseMatrix &i_sparseMatrix) {
      if (this != &i_sparseMatrix)
        deepCopy( i_sparseMatrix );

      // return current object
      return *this;
    }

    // destructor
    ~SparseMatrix() {
      // free allocated memory
      free(denseColumnMajor);
      //TODO: use standard free
      _mm_free(rowIndices); _mm_free(columnIndices); _mm_free(flatColumnMajor);
    }

    /**
     * Test multiplication of generated Matrix kernels for star matrices.
     *
     */
    bool testGeneratedMultiplication() {
      //! dense matrix containing values which are used during the computation (dense*sparse/sparse*dense).
      double* l_testMatrix;

      //! dense matrix contraing the result of the dense verification multiplication.
      double* l_denseResultMatrix;

      //! dense matrix containing the results of the sparse generated kernels.
      double* l_generatedResultMatrix;

      if( matrixType == StarMatrix) {
        //set random values to the star matrix
        setRandomValues( numberOfNonZeros, flatColumnMajor );
 
        //copy them to dense representation
        init_colmajor_from_csc( numberOfColumns, denseColumnMajor, numberOfRows,
                                flatColumnMajor, columnIndices, rowIndices);

        // iterate over different number of basis functions
        for( int l_numberOfBasisFunctions = 1; l_numberOfBasisFunctions <= 120; l_numberOfBasisFunctions++ ) {
          // allocate the test matrix and fill it with random values
          createDenseTestMatrix( l_numberOfBasisFunctions, numberOfRows, l_testMatrix );

          // alocate the result matrices and fill them with random values
          createDenseTestMatrix( l_numberOfBasisFunctions, numberOfRows, l_denseResultMatrix );
          createDenseTestMatrix( l_numberOfBasisFunctions, numberOfRows, l_generatedResultMatrix );

          // set the same values in both arrays
          copyValues(l_denseResultMatrix, l_numberOfBasisFunctions * numberOfRows, l_generatedResultMatrix);

          // check supported numbers of basis functions only
          if( l_numberOfBasisFunctions == 4 || l_numberOfBasisFunctions == 10 || l_numberOfBasisFunctions == 20) {
            // do the dense multiplication (dgemm syntax)
            mul_colmajor( l_numberOfBasisFunctions, numberOfColumns, numberOfRows,
                          l_testMatrix, l_numberOfBasisFunctions,
                          denseColumnMajor, numberOfRows,
                          l_denseResultMatrix, l_numberOfBasisFunctions );

            // do the generated multiplication
            if( l_numberOfBasisFunctions == 4 )
              generatedMatrixMultiplication_volumeStarMatrix_3D_9_4(l_testMatrix, flatColumnMajor, l_generatedResultMatrix);
            else if ( l_numberOfBasisFunctions == 10 )
              generatedMatrixMultiplication_volumeStarMatrix_3D_9_10(l_testMatrix, flatColumnMajor, l_generatedResultMatrix);
            else if ( l_numberOfBasisFunctions == 20 )
              generatedMatrixMultiplication_volumeStarMatrix_3D_9_20(l_testMatrix, flatColumnMajor, l_generatedResultMatrix);
            else
              TS_ASSERT(false);

            // compare the results
            double l_error = compare_colmajor( l_denseResultMatrix, l_generatedResultMatrix, l_numberOfBasisFunctions*numberOfRows );

            TS_ASSERT( std::abs(l_error) < errorTolerance );
          }
          // free temporary memory
          free(l_testMatrix); free(l_denseResultMatrix); free(l_generatedResultMatrix);
        }
      }
      else if( matrixType == StiffnessMatrix || matrixType == FluxMatrix ) {
        // TODO: attenuation
        for( int l_numberOfVariables = 9; l_numberOfVariables <= 9; l_numberOfVariables++ ) {
          // allocate the test matrix and fill it with random values
          createDenseTestMatrix( numberOfRows, l_numberOfVariables, l_testMatrix );

          // alocate the result matrices and fill them with random values
          createDenseTestMatrix( numberOfRows, l_numberOfVariables, l_denseResultMatrix );
          createDenseTestMatrix( numberOfRows, l_numberOfVariables, l_generatedResultMatrix );

          // set the same values in both arrays
          copyValues(l_denseResultMatrix, numberOfRows * l_numberOfVariables, l_generatedResultMatrix);

          // do the dense multiplication (dgemm syntax)
          mul_colmajor( numberOfRows, l_numberOfVariables, numberOfColumns,
                        denseColumnMajor, numberOfRows,
                        l_testMatrix, numberOfRows,
                        l_denseResultMatrix, numberOfRows );

          // do the generated multiplication (depending on the number of basis functions and type of stiffness matrix)
          #include "generated_code/unit_tests/volume_boundary_sparse_matrix_kernels.hpp_include"

          // compare the results
          double l_error = compare_colmajor( l_denseResultMatrix, l_generatedResultMatrix, numberOfRows*l_numberOfVariables );

          TS_ASSERT( std::abs(l_error) < errorTolerance );

          // free temporary memory
          free(l_testMatrix); free(l_denseResultMatrix); free(l_generatedResultMatrix);
        }
      }

        return true;
    }
};
