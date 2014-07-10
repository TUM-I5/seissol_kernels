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
 * Test suite, which tests the aligned memory allocation.
 **/
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <cxxtest/TestSuite.h>
#include <Initializer/MemoryAllocator.h>

namespace unit_test {
  class MemoryAllocatorTestSuite;
}

class unit_test::MemoryAllocatorTestSuite: public CxxTest::TestSuite {
  //private:

    //! how often should each test be repeated?
    int m_numberOfRepeats = 100;

    //! how many structs should be allocated in the corresonding test?
    int m_numberOfStructs = 97;

    //! byte aligmnments to test again
    size_t m_alignments[9] = {4, 8, 16, 32, 64, 128, 256, 512, 1024}; 

    /**
     * Prints a C++ vector
     *
     * @param i_vector std::vector.
     **/
    void printVector(const std::vector<unsigned long long> &i_vector) {
      std::copy(i_vector.begin(), i_vector.end(), std::ostream_iterator<unsigned long long>(std::cout, " "));
      std::cout << std::endl;
    }

    /**
     * Generates a random distribution of chunks.
     *
     * @retrun random #chunks of random size.
     **/
    std::vector<unsigned long long> generateRandomChunks() {
      std::vector<unsigned long long> l_chunks;
    
      // generate maximum 50 chunks
      int l_numberOfChunks = rand() % 50 + 1;

      // generate random chunks ( 84 * 9 == 756)
      for( int l_chunk = 0; l_chunk < l_numberOfChunks; l_chunk++ ) {
        l_chunks.push_back( rand() % 755 + 1 );
      }

      return l_chunks;
    }

    void testSingleChunk( int i_numberOfStructs ) {
      // generate memory allocator class
      seissol::MemoryAllocator l_memoryAllocator;

      // double precision data
      size_t l_unitSize = sizeof(double);
      
      // return valus of the memory allocator
      void* l_startOfAllocatedMemory = NULL; // where does the allocated memory start
      void** l_alignedMemoryPointers = NULL;  // addresses of aligned memory

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_numberOfRepeats; l_repeat++) {
        // set up memory layout
        std::vector<unsigned long long> l_chunkSizes;
        l_chunkSizes.push_back(rand() % 755 + 1); // single chunk

        // iterate of alignments
        for( int l_i = 0; l_i < 9; l_i++) {
          // allocate memory
          l_memoryAllocator.allocateMemory( l_chunkSizes, l_unitSize,  m_alignments[l_i],
                                            &l_startOfAllocatedMemory, &l_alignedMemoryPointers,
                                            i_numberOfStructs );

          // assert the memory has been allocated
          TS_ASSERT( l_startOfAllocatedMemory != NULL );
          TS_ASSERT( l_alignedMemoryPointers  != NULL );

          // assert our pointers are aligned
          for( int l_struct = 0; l_struct < i_numberOfStructs; l_struct++ ) {
            TS_ASSERT( ((unsigned long long) l_alignedMemoryPointers[l_struct]) % m_alignments[l_i] == 0 );
          }
        }
      }
 
     // free memory
     l_memoryAllocator.freeMemory();
    }

    void testMultipleChunks( int i_numberOfStructs ) {
      // generate memory allocator class
      seissol::MemoryAllocator l_memoryAllocator;

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_numberOfRepeats; l_repeat++) {
        // set up memory layout
        std::vector<unsigned long long> l_chunkSizes = generateRandomChunks();

        //printVector(l_chunkSizes);
 
        // iterate of alignments
        for( int l_i = 0; l_i < 9; l_i++) { 
          size_t l_elementSize = sizeof(double); // double precision data

          // return valus of the memory allocator
          void* l_startOfAllocatedMemory = NULL; // where does the allocated memory start
          void** l_alignedMemoryPointers = NULL;  // addresses of aligned memory

          // allocate memory
          l_memoryAllocator.allocateMemory( l_chunkSizes, l_elementSize,  m_alignments[l_i],
                                           &l_startOfAllocatedMemory, &l_alignedMemoryPointers,
                                            i_numberOfStructs );

          // assert the memory has been allocated
          TS_ASSERT( l_startOfAllocatedMemory != NULL );
          TS_ASSERT( l_alignedMemoryPointers  != NULL );

          // iterate over structs
          for( int l_struct = 0; l_struct < i_numberOfStructs; l_struct++ ) {    
            for( unsigned long long l_chunk = 0; l_chunk < l_chunkSizes.size(); l_chunk++ ) {
              // assert our pointers are aligned
              TS_ASSERT( ((unsigned long long) l_alignedMemoryPointers[l_struct * l_chunkSizes.size() + l_chunk]) % m_alignments[l_i] == 0 );
            }
    
            for( unsigned long long l_chunk = 1; l_chunk < l_chunkSizes.size(); l_chunk++ ) {
              // assert the array can hold the chunks
              unsigned long long l_gapSize = ( (unsigned long long) l_alignedMemoryPointers[l_struct * l_chunkSizes.size() + l_chunk])
                                          -( (unsigned long long) l_alignedMemoryPointers[l_struct * l_chunkSizes.size() + l_chunk-1]);
              TS_ASSERT( l_gapSize / (double) l_chunkSizes[l_chunk-1] >= l_elementSize );

              // assert each chunk acrross the structs is of equal size
              if( l_struct > 0 ) {
                unsigned long long l_gapSizeFirstStruct = ( (unsigned long long) l_alignedMemoryPointers[(l_struct-1) * l_chunkSizes.size() + l_chunk])
                                                         -( (unsigned long long) l_alignedMemoryPointers[(l_struct-1) * l_chunkSizes.size() + l_chunk-1]);

                unsigned long long l_gapSizeSecondStruct = ( (unsigned long long) l_alignedMemoryPointers[(l_struct)  * l_chunkSizes.size() + l_chunk])
                                                          -( (unsigned long long) l_alignedMemoryPointers[(l_struct)  * l_chunkSizes.size() + l_chunk-1]);

                TS_ASSERT( l_gapSizeFirstStruct == l_gapSizeSecondStruct );
              }
            }
          }


        }
      }

      l_memoryAllocator.freeMemory();
    }

  public:
    /**
     * Set up unit tests.
     */
    void setUp() {
      // set seed for random number generation
      srand (time(NULL));
    }
    
    /**
     * Test a single struct and a single chunk.
     */
    void testSingleStructSingleChunk() {
      testSingleChunk( 1 );
    }

    /**
     * Test a single struct and a single chunk.
     */
    void testMultipleStructsSingleChunk() {
      testSingleChunk( m_numberOfStructs );
    }

    /**
     * Test single struct and multiple chunks.
     **/
    void testSingleStructMultipleChunks() {
      testMultipleChunks( 1 );
    }

    /**
     * Test multiple structs and multiple chunks.
     **/
    void testMultipleStructMultipleChunks() {
       testMultipleChunks( m_numberOfStructs );
    }
};
