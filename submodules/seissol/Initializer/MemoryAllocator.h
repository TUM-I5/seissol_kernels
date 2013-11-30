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

#ifndef MEMORYALLOCATOR_H_
#define MEMORYALLOCATOR_H_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <utils/logger.h>

namespace seissol {
  class MemoryAllocator;
}

/**
 * Allocates aligned memory for dynamic chunk sizes.
 **/
class seissol::MemoryAllocator {
  //private:
    //! holds all memory addresses, which point to data arrays and have been returned by mallocs calling functions of the memory allocator.
    std::vector< void* > m_dataMemoryAddresses;

    //! holds all memory addresses to arrays representing the chunk structure.
    std::vector< void** > m_pointerMemoryAddresses;

    /**
     * Prints the memory alignment of in terms of relative start and ends in bytes.
     *
     * @param i_memoryAlignment memory alignment.
     **/
    void printMemoryAlignment( std::vector< std::vector<unsigned long long> > i_memoryAlignment );

  public:
    MemoryAllocator();

    /**
     * Allocates a single chunk of memory with the given size and alignment.
     * @param i_size size of the chunk in byte.
     * @param i_alignment alignment of the memory chunk in byte.
     * @param o_pointerToMemory pointer, which points to the aligned memory of the given size.
     **/
    inline void allocateMemory( size_t i_size,
                                size_t i_alignment,
                                void** io_pointerToMemory ) {
      // do the malloc
#ifdef __INTEL_COMPILER
     *io_pointerToMemory = _mm_malloc( i_size, i_alignment );
#else
      //TODO: Add support for more compilers
#error "Error: Only the Intel compiler is supported.
#endif
    };

    /**
     * Allocates aligned memory with a given memory layout.
     *
     * For example two chunks of size 5 and 10 of double precision with an 64 byte alignment, would give an
     * aligned pointer to the start of each chunk, with an 24 byte memory gap between the two chunks because
     * of the alignment. In total 5 * 8 + 24 (gap) + 10 * 8 = 144 bytes are occupied by this setting.
     *
     * @param i_chunkSizes individual chunk sizes (entries per chunk) in memory, which are allocated successive and aligned individually.
     * @param i_elementSize size of each entry.
     * @param i_alignment alignment of the memory chunks.
     * @param o_startOfAllocatedMemory pointer to the starting address as stated by malloc(). This needs to be passed in free()
     * @param o_alignedMemoryPointers array of pointers, which point to the starting addresses of the individual chunks. If i_numberOfStructs > 1
     *                                the chunk at position [l_struct * i_chunkSizes.size() + 0] is the first chunk of the struct with id l_struct. 
     * @param i_numberOfStructs number of structures, which are allocated. Each structure has the layout defined in the chunk sizes.
     */
    void allocateMemory( std::vector<unsigned long long> &i_chunkSizes,
                         size_t                           i_elementSize,
                         size_t                           i_alignment,
                         void**                           o_startOfAllocatedMemory,
                         void***                          o_alignedMemoryPointers,
                         unsigned long long               i_numberOfStructs = 1 );

    /**
     * Frees all memory, which was allocated by functions of the MemoryAllocator.
     **/
    void freeMemory();
};

#endif
