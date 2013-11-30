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
 * Aligned memory allocation.
 **/
#include "MemoryAllocator.h"

seissol::MemoryAllocator::MemoryAllocator() {
}

void seissol::MemoryAllocator::printMemoryAlignment( std::vector< std::vector<unsigned long long> > i_memoryAlignment ) {
  logDebug() << "printing memory alignment per struct";
  for( unsigned long long l_i = 0; l_i < i_memoryAlignment.size(); l_i++ ) {
    logDebug() << i_memoryAlignment[l_i][0] << ", " << i_memoryAlignment[l_i][1];
  }
}

void seissol::MemoryAllocator::allocateMemory( std::vector<unsigned long long> &i_chunkSizes,
                                               size_t                           i_elementSize,
                                               size_t                           i_alignment,
                                               void                           **o_startOfAllocatedMemory,
                                               void                          ***o_alignedMemoryPointers,
                                               unsigned long long               i_numberOfStructs ) {
  // at least one struct should be allocated
  assert( i_numberOfStructs != 0 );

  //! total size of the individual chunks in case of unaligned memory allocations
  unsigned long long l_unalignedMemorySizePerStruct = i_chunkSizes[0];

  //! memory layout of the chunks 
  std::vector< std::vector<unsigned long long> > l_memoryLayout( i_chunkSizes.size(), std::vector<unsigned long long>(2));

  // first entry is aligned
  l_memoryLayout[0][0] = 0;
  l_memoryLayout[0][1] = i_chunkSizes[0] * i_elementSize;

  // compute memory alignment
  for( unsigned long long l_i = 1; l_i < i_chunkSizes.size(); l_i++ ) {
    // add the current chunk
    l_unalignedMemorySizePerStruct += i_chunkSizes[l_i];

    /* 
     * align start of this chunk with the end of previous chunk as reference.
     *
     * <pre>
     *
     *                           last aligned         first possible             
     *                           address in the       start of this                   aligned start
     *                           previos chunk        chunk                           of this chunk
     *                                |                  |                                |
     *                               \_/                \_/                              \_/
     *
     * [..., * , * , * , * , ..., * , * ,  ... , * , * , - , - , - , ..., - , - , - , - , *, *, *, *, *, ...]
     *       _                                       _  |______________________________|
     *      / \                                     / \                 |
     *       |                                       |             unused memory            
     *    aligned start                           end of
     *    of previous                             previous
     *    chunk                                   chunk
     *
     * </pre>
     */
    //! how far does the previous chunk reach into the aligned memory
    unsigned int l_memoryOvershoot = l_memoryLayout[l_i-1][1]%i_alignment;
    
    // get an aligned starting position
    unsigned long long l_start = l_memoryLayout[l_i-1][1] + (i_alignment - l_memoryOvershoot)%i_alignment;

    // compute length of this memory sections, first possible (relative) address of a new memory chunk
    unsigned long long l_end   = l_start + i_chunkSizes[l_i] * i_elementSize;

    // assert the size of this chunk is covered
    assert( (l_end - l_start) == (i_chunkSizes[l_i] * i_elementSize) );

    // assert we don't have a gap, which is larger than the actual alignment size
    assert( l_start - l_memoryLayout[l_i-1][1] < i_alignment);

    // set start and end
    l_memoryLayout[l_i][0] = l_start;
    l_memoryLayout[l_i][1] = l_end;
  }

  if( i_numberOfStructs > 1 ) {
#ifndef NDEBUG
    printMemoryAlignment( l_memoryLayout );
#endif

   /**
    * Increase size of the last chunk to ensure alignment of the following struct.
    *
    * Example for a struct containing chunks of five and two byte size: 
    *
    * <pre>
    *
    *                   first possible         first possible                                                               increase size of last chunk
    *                   non-aligned start      non-aligned start                                                            to align the start of the next struct
    *                   of the second chunk    of the next struct                                                                  / \
    *                           |                   |                                                                             |   |
    *                          \_/                 \_/                                                                           \_/ \_/
    *                  
    * [..., * , * , * , * , * , - , - , - , * , * , - , - , - , - , - , ...]  --->  [..., * , * , * , * , * , - , - , - , * , * , * , * , - , - , - , ...]
    *       _                               _               _                                                                             _
    *      / \                             / \             / \                                                                           / \
    *       |                               |               |                                                                             |
    *  aligned start                   aligned start     aligned start                                                                 aligned start
    *  of the struct                   of the second     of the next                                                                   of the next
    *  and first chunk                 chunk             struct                                                                        struct
    *
    *      |______________________________________|                                     |_______________________________________________|____________...
    *             single struct allocation                                                     first struct (increased in size)           aligned second struct
    *
    * </pre>
    *
    */
    
    logDebug() << "increasing size of last chunk by " << l_memoryLayout.back()[1]%i_alignment << " bytes.";

    //! how far does the last chunk reach into the aligned memory?
    unsigned int l_memoryOvershoot = l_memoryLayout.back()[1]%i_alignment;

    // increase the size of the last chunk
    l_memoryLayout.back()[1]       += (i_alignment - l_memoryOvershoot)%i_alignment;
  }

#ifndef NDEBUG
  printMemoryAlignment( l_memoryLayout );
#endif

  //! total memory sizes per struct
  l_unalignedMemorySizePerStruct *= i_elementSize;
  unsigned long long l_memorySizePerStruct = l_memoryLayout.back()[1];

  // assert we are below 3 megabyte per struct, (increase this value, if more memory is necessary)
  // maximum are the flux matrices: \f$ 84^2 \cdot 52 \cdot 8 = 2935296 byte \f$ for 84 basis functions or 7th order.
  assert( l_memorySizePerStruct < 3e6 );

  // assert we have at leat the same amount of memory than the unaligned case
  assert( l_unalignedMemorySizePerStruct <= l_memorySizePerStruct );
  
  // assert all gaps are covered (upper bound)
#ifndef NDEBUG
  if( i_numberOfStructs == 1 )
    // we have to fill up a gap for each chunk excluding the last one
    assert( l_memorySizePerStruct - l_unalignedMemorySizePerStruct <= (i_chunkSizes.size()-1) * (i_alignment - 1) );
  else
    // we have to fill up a gap for each chunk
    assert( l_memorySizePerStruct - l_unalignedMemorySizePerStruct <= i_chunkSizes.size() * (i_alignment - 1) );
#endif

  //! total memory size for all structs
  unsigned long long l_unalignedMemorySize = l_unalignedMemorySizePerStruct * i_numberOfStructs;
  unsigned long long l_memorySize          = l_memorySizePerStruct          * i_numberOfStructs;

  // allocate memory
  allocateMemory( l_memorySize, i_alignment, o_startOfAllocatedMemory );
  char* l_startOfAlignedMemory = (char*) *o_startOfAllocatedMemory;

  // store address locally
  m_dataMemoryAddresses.push_back( *o_startOfAllocatedMemory );

  // char pointer to the start of memory (char for pointer arithmetics)
  assert( sizeof(char) == 1 ); 

  // generate pointers to beginnings of individual chunks
  *o_alignedMemoryPointers = (void**) calloc( i_chunkSizes.size()*i_numberOfStructs, sizeof(void*) );

  // store address locally
  m_pointerMemoryAddresses.push_back(*o_alignedMemoryPointers);

  // iterate over the number of structs
  for( int l_struct = 0; l_struct < i_numberOfStructs; l_struct++ ) {
    // set a pointer to the start of each chunk
    for( int l_chunk = 0; l_chunk < l_memoryLayout.size(); l_chunk++ ) {
      // start at the address returned by the malloc
      (*o_alignedMemoryPointers)[l_struct * l_memoryLayout.size() + l_chunk]  = (void*) l_startOfAlignedMemory;
      
      // jump over the previous structs
      (*o_alignedMemoryPointers)[l_struct * l_memoryLayout.size() + l_chunk] += l_struct * l_memorySizePerStruct;

      // jump over the previous chunks
      (*o_alignedMemoryPointers)[l_struct * l_memoryLayout.size() + l_chunk] += l_memoryLayout[l_chunk][0];
    }
  }

  // print some information
  logDebug() << "total size of allocated memory: " << l_memorySize << " bytes." << std::endl
             << "this is a overhead (compared to unaligned allocations) of: ~"
             << l_memorySize - l_unalignedMemorySize
             << " bytes (" << (l_memorySize / (double) l_unalignedMemorySize - 1.0)*100 << "%)."
             << std::endl;
}

void seissol::MemoryAllocator::freeMemory() {
  assert( m_dataMemoryAddresses.size() == m_pointerMemoryAddresses.size());

  for( unsigned long long l_i = 0; l_i < m_dataMemoryAddresses.size(); l_i++ ) {
#ifdef __INTEL_COMPILER
    _mm_free( m_dataMemoryAddresses[l_i] );
#else
  //TODO: Add support for more compilers
#error "Error: Only the Intel compiler is supported.
#endif

    free( m_pointerMemoryAddresses[l_i] );
  }

  // reset memory vectors
  m_dataMemoryAddresses.clear();
  m_pointerMemoryAddresses.clear();
}
