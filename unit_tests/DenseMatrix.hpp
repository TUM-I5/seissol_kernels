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
 * Simple dense matrix functionality for unit tests.
 **/

#ifndef DENSEMATRIX_HPP_
#define DENSEMATRIX_HPP_

#include <Initializer/MemoryAllocator.h>
#include <Initializer/XmlParser.hpp>

const int s_maximumOrder = 6;
// TODO: Add a zero tolerance, which matches the size of the matrix: It doesn't make sense to compare values > 10,000.0 with a zero tolerance of 10e-15..
const double s_zeroTolerance = 10e-10;

namespace unit_test {
  class DenseMatrix;
}

/**
 * Simple dense matrix functionalities: Read XML, set-random values, do matrix-multiplications..
 **/
class unit_test::DenseMatrix {
  //private:
    //! aligned memory allocation
    seissol::MemoryAllocator m_memoryAllocator;

    //! vectors, which hold information about our matrices
    std::vector< unsigned int > m_matrixIds;
    std::vector< std::string  > m_matrixNames;
    std::vector< bool         > m_matrixSparsities;
    std::vector< unsigned int > m_matrixNumberOfRows;
    std::vector< unsigned int > m_matrixNumberOfColumns;

    std::vector< std::vector<unsigned int> > m_matrixRows;
    std::vector< std::vector<unsigned int> > m_matrixColumns;
    std::vector< std::vector<double>       > m_matrixValues;

  public:
    /**
     * Fills array with random values.
     *
     * @param i_length length of the array.
     * @param o_array array to fill.
     **/
    void setRandomValues( int i_length, double *o_array ) {
      for( int l_i = 0; l_i < i_length; l_i++) {
        o_array[l_i] = ((double)rand()/(double)RAND_MAX)*10.0;
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
                       double *io_matrix ) {
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
                            const double *i_denseMatrix,
                                  double *o_sparseMatrix ) {
     // #(dense elements)
     int l_numberOfDenseElements = i_numberOfRows * i_numberOfColumns;

     // id of the sparse element
     int l_sparseElement = 0;

     // iterate over dense matrix elements
     for( int l_denseElement = 0; l_denseElement < l_numberOfDenseElements; l_denseElement++ ) {
       if( i_denseMatrix[ l_denseElement ] > 1e-15 ) {
         o_sparseMatrix[ l_sparseElement ] = i_denseMatrix[ l_denseElement ];
         
         // next sparse element
         l_sparseElement++;
       }
     }
    }

    /**
     * Executes C += A.B with a simple for-loop.
     *
     * @param i_m #(dense rows) of the left matrix
     * @param i_n #(dense columns) of the right matrix
     * @param i_k #(dense columns) of the left matrix
     * @o_a pointer to matrix A.
     * @o_b pointer to matrix B. 
     * @o_c pointer to matrix C.
     **/
    void executeStandardMultiplication(       int     i_m,       int     i_n, int     i_k,
                                        const double *i_a, const double *i_b, double *io_c ) {
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
    void printDenseMatrix( int i_numberOfRows, int i_numberOfColumns, double* i_matrix ) {
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
    void checkResult( int i_length, double* i_array1, double* i_array2 ) {
      // check every value individually
      for( int l_i = 0; l_i < i_length; l_i++) {
        TS_ASSERT_DELTA(i_array1[l_i], i_array2[l_i], s_zeroTolerance);
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
      m_matrixValues.push_back( std::vector< double >());
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
                           double*      o_matrix ) {
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

      // assert we expect the correct size
      TS_ASSERT_EQUALS( m_matrixNumberOfRows[l_position],    i_numberOfRows );
      TS_ASSERT_EQUALS( m_matrixNumberOfColumns[l_position], i_numberOfColumns );

      // flux solver -> random values
      if( i_id == 52 ) {
        // assert correct matrix id
        TS_ASSERT_EQUALS( m_matrixIds[l_position], 52 );
        setRandomValues( i_numberOfRows * i_numberOfColumns,
                        o_matrix );
      } 
      // flux matrix, stiffness matrix or star matrix
      else {
        std::fill( o_matrix, o_matrix + (m_matrixNumberOfRows[l_position] * m_matrixNumberOfColumns[l_position]), 0 );
        for( unsigned int l_element = 0; l_element < m_matrixRows[l_position].size(); l_element++ ) {
          unsigned int l_denseIndex = (m_matrixColumns[l_position][l_element]-1) * m_matrixNumberOfRows[l_position] + (m_matrixRows[l_position][l_element]-1);
          if( i_id == 59 ) { // star matrix
            o_matrix[ l_denseIndex ] = ((double)rand()/(double)RAND_MAX)*10.0;
          }
          else{ // flux or stiffness matrix
            o_matrix[ l_denseIndex ] = m_matrixValues[l_position][l_element];
          }
        }
      }
    }
};

#endif
