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
 * Simple dense matrix functionality for unit tests.
 **/

#ifndef DENSEMATRIX_HPP_
#define DENSEMATRIX_HPP_

#include <typedefs.hpp>
#include <Initializer/XmlParser.hpp>
#include <cxxtest/TestSuite.h>

const int s_maximumOrder = 6;
#ifdef DOUBLE_PRECISION
#if CONVERGENCE_ORDER < 6
const long double s_zeroTolerance = 10e-11;
#else
const long double s_zeroTolerance = 10e-09;
#endif
#endif

#ifdef SINGLE_PRECISION
#if CONVERGENCE_ORDER < 6
const long double s_zeroTolerance = 10e-05;
#else
const long double s_zeroTolerance = 10e-03;
#endif
#endif

namespace unit_test {
  class DenseMatrix;
}

/**
 * Simple dense matrix functionalities: Read XML, set-random values, do matrix-multiplications..
 **/
class unit_test::DenseMatrix {
  //private:
    //! vectors, which hold information about our matrices
    std::vector< unsigned int > m_matrixIds;
    std::vector< std::string  > m_matrixNames;
    std::vector< bool         > m_matrixSparsities;
    std::vector< unsigned int > m_matrixNumberOfRows;
    std::vector< unsigned int > m_matrixNumberOfColumns;

    std::vector< std::vector<unsigned int> > m_matrixRows;
    std::vector< std::vector<unsigned int> > m_matrixColumns;
    std::vector< std::vector<real>       > m_matrixValues;

  public:
    /**
     * Fills array with random values.
     *
     * @param i_length length of the array.
     * @param o_array array to fill.
     **/
    void setRandomValues( int i_length, real *o_array ) {
      for( int l_i = 0; l_i < i_length; l_i++) {
        o_array[l_i] = ((real)rand()/(real)RAND_MAX)*0.1;
      }
    }

    /**
     * Sets the all rows, which are behind the given number of non zero rows to zero.
     * @param i_numberOfRows #rows of the matrix.
     * @param i_numberOfColumns #columns of the matrix.
     * @param i_numberOfNonZeroRows #(top rows), which are set not zero.
     * @param io_matrix column-major matrix, which is modified
     **/
    void setZeroBlock( int     i_numberOfRows,
                       int     i_numberOfColumns,
                       int     i_numberOfNonZeroRows,
                       real    *io_matrix ) {
      for( int l_column = 0; l_column < i_numberOfColumns; l_column++ ) {
        for( int l_row = 0; l_row < i_numberOfRows; l_row++ ) {
          // check if this a zero row
          if( l_row >= i_numberOfNonZeroRows ) {
            // compute matrix index
            int l_matrixIndex = l_column * i_numberOfRows + l_row;
            // set entry to zero
            io_matrix[l_matrixIndex] = 0;
          }
        }
      }
    }
#if 0
    /**
     * Allocate memory and set random values for dense matrix multiplication setup.
     *   Remark: C1 and C2 are set to the same random values 
     *
     * @param i_m #(dense rows) of the left matrix
     * @param i_n #(dense columns) of the left matrix
     * @param i_k #(dense columns) of the right matrix
     * @o_a pointer for matrix A.
     * @o_b pointer for matrix B. 
     * @o_c1 pointer for matrix C1.
     * @o_c2 pointer for matrix C2.
     **/
    void allocateMemoryAndSetRandomValues( int     i_m,  int     i_n,  int    i_k,
                                           double **o_a, double **o_b, double **o_c1, double **o_c2 ) {
      // generate chunk sizes
      std::vector<unsigned long long> l_chunkSizes;
      l_chunkSizes.push_back(i_m * i_k); // A
      l_chunkSizes.push_back(i_k * i_n); // B
      l_chunkSizes.push_back(i_m * i_n); // C1
      l_chunkSizes.push_back(i_m * i_n); // C2

      // double precision data and alignment
      size_t l_elementSize = sizeof(double);
      size_t l_alignment = 64;
      
      // return valus of the memory allocator
      void* l_startOfAllocatedMemory = NULL; // where does the allocated memory start
      void** l_alignedMemoryPointers = NULL;  // addresses of aligned memory
    
      // allocate memory
      m_memoryAllocator.allocateMemory( l_chunkSizes, l_elementSize,  l_alignment,
                                       &l_startOfAllocatedMemory, &l_alignedMemoryPointers );

      // set return pointers
      *o_a =  (double*) l_alignedMemoryPointers[0];
      *o_b =  (double*) l_alignedMemoryPointers[1];
      *o_c1 = (double*) l_alignedMemoryPointers[2];
      *o_c2 = (double*) l_alignedMemoryPointers[3];

      // set random values
      setRandomValues( i_m * i_k, *o_a  );
      setRandomValues( i_k * i_n, *o_b  );
      setRandomValues( i_m * i_n, *o_c1 );

      // copy values of C1 to C2
      for( int l_i = 0; l_i < i_m * i_n; l_i++) {
        (*o_c2)[l_i] = (*o_c1)[l_i];
      }
    }
#endif

    /**
     * Copies the non-zeros of the dense matrix to the sparse matrix.
     *   Remark: No additional checking is done. It's assumed that the
     *           sparse pointer is large enough to hold all non-zeros.
     *
     * @param i_numberOfRows #rows of the matrix.
     * @param i_numberOfColumns #columns of the matrix.
     * @param i_denseMatrix dense matrix representation (column-major).
     * @param o_sparseMatrix pointer to nnz of the sparse matrix (column major).
     **/
    void copyDenseToSparse(       int     i_numberOfRows,
                                  int     i_numberOfColumns,
                            const real   *i_denseMatrix,
                                  real   *o_sparseMatrix ) {
     // #(dense elements)
     int l_numberOfDenseElements = i_numberOfRows * i_numberOfColumns;

     // id of the sparse element
     int l_sparseElement = 0;

     // iterate over dense matrix elements
     for( int l_denseElement = 0; l_denseElement < l_numberOfDenseElements; l_denseElement++ ) {
       if( std::abs( i_denseMatrix[ l_denseElement ] ) > 1e-15 ) {
         o_sparseMatrix[ l_sparseElement ] = i_denseMatrix[ l_denseElement ];
         
         // next sparse element
         l_sparseElement++;
       }
     }
    }

    /*
     * Copies a submatrix of A (sizes of B) to B.
     * If B doesn't fit in A zeros are set.
     *
     * @param i_A values of matrix A.
     * @param i_aNumberOfRows number of rows of A.
     * @param i_aNumberOfColumns number of columns of A.
     * @param i_aLeadingDimension leading dimension of A.
     * @param o_B values of matrix B, which will be set.
     * @param i_bNumberOfRows number of rows of matrix B.
     * @param i_bNumberOfColumns number of columns of matrix B.
     * @param i_bLeadingDimension leading dimension of B.
     */
    void copySubMatrix( const real* i_A,
                        const unsigned int i_aNumberOfRows,
                        const unsigned int i_aNumberOfColumns,
                        const unsigned int i_aLeadingDimension,
                              real* o_B,
                        const unsigned int i_bNumberOfRows,
                        const unsigned int i_bNumberOfColumns,
                        const unsigned int i_bLeadingDimension ) {
      // set matrix B to zero
      for( unsigned int l_index = 0; l_index < i_bLeadingDimension*i_bNumberOfColumns; l_index++ ) {
        o_B[l_index] = (real) 0;
      }

      // copy the entries
      for( unsigned int l_column = 0; l_column < std::min( i_aNumberOfColumns, i_bNumberOfColumns ); l_column++ ) {
        for( unsigned int l_row = 0; l_row < std::min( i_aNumberOfRows, i_bNumberOfRows ); l_row++ ) {
          unsigned int l_aIndex = l_column * i_aLeadingDimension + l_row;
          unsigned int l_bIndex = l_column * i_bLeadingDimension + l_row;

          o_B[l_bIndex] = i_A[l_aIndex];
        }
      }
     
    }

    /*
     * Copies a submatrix of A (sizes of B) to B.
     * If B doesn't fit in A zeros are set.
     *
     * @param i_A values of matrix A.
     * @param i_aNumberOfRows number of rows of A.
     * @param i_aNumberOfColumns number of columns of A.
     * @param o_B values of matrix B, which will be set.
     * @param i_bNumberOfRows number of rows of matrix B.
     * @param i_bNumberOfColumns number of columns of matrix B.
     */
    void copySubMatrix( const real* i_A,
                        const unsigned int i_aNumberOfRows,
                        const unsigned int i_aNumberOfColumns,
                              real* o_B,
                        const unsigned int i_bNumberOfRows,
                        const unsigned int i_bNumberOfColumns ) {
      copySubMatrix( i_A,
                     i_aNumberOfRows,
                     i_aNumberOfColumns,
                     i_aNumberOfRows,
                     o_B,
                     i_bNumberOfRows,
                     i_bNumberOfColumns,
                     i_bNumberOfRows
                   );
    }

    /**
     * Executes C += A.B or C = A.B with a simple for-loop.
     *
     * @param i_m #(dense rows) of the left matrix
     * @param i_n #(dense columns) of the right matrix
     * @param i_k #(dense columns) of the left matrix
     * @o_a pointer to matrix A.
     * @o_b pointer to matrix B. 
     * @o_c pointer to matrix C.
     * @i_add true: C += A.B, false C = A.B
     **/
    void executeStandardMultiplication(       int     i_m,       int     i_n, int     i_k,
                                        const real   *i_a, const real   *i_b, real   *io_c,
                                              bool    i_add = true ) {
      // set result matrix to zero first of required
      if( i_add == false ) {
        for( int l_positionC = 0; l_positionC < i_m * i_n; l_positionC++ ) {
          io_c[l_positionC] = 0;
        }
      }

      for( int l_m = 0; l_m < i_m; l_m++ ) {
        for( int l_n = 0; l_n < i_n; l_n++ ) {
          for( int l_k = 0; l_k < i_k; l_k++ ) {
            // index calculation
            int l_positionA = l_k * i_m + l_m;
            int l_positionB = l_n * i_k + l_k;
            int l_positionC = l_n * i_m + l_m;

            // do the multiplication
            io_c[l_positionC] += i_a[l_positionA] * i_b[l_positionB];
          }
        }
      }
    }

    /**
     * Prints a dense matrix.
     *
     * @param i_numberOfRows number of rows.
     * @param i_numberOfColumnes number of columns.
     * @param i_matrix pointer to matrix elements.
     **/
    static void printDenseMatrix( int i_numberOfRows, int i_numberOfColumns, real* i_matrix ) {
      for( int l_row = 0; l_row < i_numberOfRows; l_row++ ) {
        for( int l_column = 0; l_column < i_numberOfColumns; l_column++ ) {
          std::cout << i_matrix[ l_column * i_numberOfRows + l_row ] << " ";
        }
        std::cout << std::endl;
      } 
    }

    /**
     * Checks if each entry in the first array differs less than the zero tolerance from the ones in the second array.
     *
     * @param i_length length of both arrays.
     * @param i_array1 first array.
     * @param i_array2 second array.
     **/
    void checkResult( int i_length, const real* i_array1, const real* i_array2 ) {
      // check every value individually
      for( int l_i = 0; l_i < i_length; l_i++) {
        if( std::abs(i_array2[l_i]) > 10e-5 ) { 
          TS_ASSERT_DELTA( (long double) ( (long double) i_array1[l_i] / (long double) i_array2[l_i]), 1.0, s_zeroTolerance);
        }
        else {
          TS_ASSERT_DELTA( i_array1[l_i], i_array2[l_i], s_zeroTolerance );
        }
      }
    }

    /*
     * Check that the submatrix of A with size B matches B.
     *
     * @param i_A values of matrix A.
     * @param i_aNumberOfRows number of rows of A.
     * @param i_aNumberOfColumns number of columns of A.
     * @param o_B values of matrix B, which will be set.
     * @param i_bNumberOfRows number of rows of matrix B.
     * @param i_bNumberOfColumns number of columns of matrix B
     */
    void checkSubMatrix( const real* i_A,
                         const unsigned int i_aNumberOfRows,
                         const unsigned int i_aNumberOfColumns,
                         const real* i_B,
                         const unsigned int i_bNumberOfRows,
                         const unsigned int i_bNumberOfColumns ) {
      // assert that B is a submatrix of A
      TS_ASSERT_LESS_THAN_EQUALS( i_bNumberOfRows,    i_aNumberOfRows    );
      TS_ASSERT_LESS_THAN_EQUALS( i_bNumberOfColumns, i_aNumberOfColumns );

      for( unsigned int l_column = 0; l_column < i_bNumberOfColumns; l_column++ ) {
        unsigned int l_aStart = l_column*i_aNumberOfRows;
        unsigned int l_bStart = l_column*i_bNumberOfRows;

        // check result column wise
        checkResult( i_bNumberOfRows, &i_A[l_aStart], &i_B[l_bStart] );
      }
    }

    /**
     * Reads the matrices from the XML-file and stores the internally.
     *
     * @param i_pathToXmlFile location of the XML-file, which contains information about the matrics.
     **/
    void readMatrices( std::string  i_pathToXmlFile ) {
      //!xml-parser
      seissol::XmlParser l_matrixReader( i_pathToXmlFile );

      // read global matrices
      l_matrixReader.readGlobalMatrices( "flux",
                                         m_matrixIds, m_matrixNames,
                                         m_matrixNumberOfRows, m_matrixNumberOfColumns, m_matrixSparsities,
                                         m_matrixRows, m_matrixColumns, m_matrixValues );     

      // read flux solver
      l_matrixReader.readLocalMatrixStructure( "fluxSolver",
                                               m_matrixIds,
                                               m_matrixNumberOfRows, m_matrixNumberOfColumns, m_matrixSparsities,
                                               m_matrixRows, m_matrixColumns );
      // synchronize size of values and name (which are not set by the local matrix reader)
      m_matrixValues.push_back( std::vector< real >());
      m_matrixNames.push_back( std::string() );

      // read stiffness matrices
      l_matrixReader.readGlobalMatrices( "stiffness",
                                         m_matrixIds, m_matrixNames,
                                         m_matrixNumberOfRows, m_matrixNumberOfColumns, m_matrixSparsities,
                                         m_matrixRows, m_matrixColumns, m_matrixValues ); 

      // read star matrix
      l_matrixReader.readLocalMatrixStructure( "starMatrix",
                                               m_matrixIds,
                                               m_matrixNumberOfRows, m_matrixNumberOfColumns, m_matrixSparsities,
                                               m_matrixRows, m_matrixColumns );
    }

    /**
     * Initializes a matrix with the specified id.
     *   Global matrices are set to the corresponding values,
     *   local matrices initialized with random values.
     *
     * @param i_id id of the matrix {0, .., 51} -> flux matrix, {52} -> jacobian (flux solver), {53, 54, 55} -> stiffness matrices, {56, 57, 58} -> transposed stiffness matrices, {59} -> jacobian (star matrix)
     * @param i_numberOfRows number of rows.
     * @param i_numberOfColumns number of columns.
     * @param o_matrix matrix, which is initialized.
     **/
    void initializeMatrix( unsigned int i_id,
                           unsigned int i_numberOfRows,
                           unsigned int i_numberOfColumns,
                           real*        o_matrix ) {
      // assert the matrices have been read
      assert( m_matrixIds.size() > 0 );

      // get the position of our matrix
      unsigned int l_position;
      for( l_position = 0; l_position < m_matrixIds.size(); l_position++ ) {
        if( i_id == m_matrixIds[l_position] ) {
          break;
        }
        else {
          // matrix with the given id couldn't be found
          TS_ASSERT_DIFFERS( l_position, m_matrixIds.size() -1 );
        }
      }

      // flux solver -> random values
      if( i_id == 52 ) {
        // assert correct matrix id
        TS_ASSERT_EQUALS( m_matrixIds[l_position], 52 );
        setRandomValues( i_numberOfRows * i_numberOfColumns,
                         o_matrix );
      } 
      // flux matrix, stiffness matrix or star matrix
      else {
        std::fill( o_matrix, o_matrix + i_numberOfRows * i_numberOfColumns, 0 );
        for( unsigned int l_element = 0; l_element < m_matrixRows[l_position].size(); l_element++ ) {
          // assert that we can save the non-zero element
          TS_ASSERT_LESS_THAN_EQUALS( m_matrixRows[l_position][l_element],    i_numberOfRows);
          TS_ASSERT_LESS_THAN_EQUALS( m_matrixColumns[l_position][l_element], i_numberOfColumns);

          unsigned int l_denseIndex = (m_matrixColumns[l_position][l_element]-1) * i_numberOfRows + (m_matrixRows[l_position][l_element]-1);

          // assert this a valid index
          TS_ASSERT_LESS_THAN( l_denseIndex, i_numberOfRows*i_numberOfColumns );

          if( i_id == 59 ) { // star matrix
            o_matrix[ l_denseIndex ] = ((real)rand()/(real)RAND_MAX)*10.0;
          }
          else if( i_id >= 56 && i_id <= 58 ) { // add minus-sign for time kernel
            o_matrix[ l_denseIndex ] = -m_matrixValues[l_position][l_element];
          }
          else{ // flux or stiffness matrix
            o_matrix[ l_denseIndex ] =  m_matrixValues[l_position][l_element];
          }
        }
      }
    }
};

#endif
