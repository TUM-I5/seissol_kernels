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
 * Memory management of SeisSol.
 **/
#include "MemoryManager.h"

seissol::initializers::MemoryManager::MemoryManager( const seissol::XmlParser &i_matrixReader ) {
  // allocate memory for the pointers to the individual matrices
  m_fluxMatrixPointers      = new double*[53];
  m_stiffnessMatrixPointers = new double*[7];

  // initialize global matrices
  initializeGlobalMatrices( i_matrixReader );
}

seissol::initializers::MemoryManager::~MemoryManager() {
  // free members
  delete[] m_fluxMatrixPointers;
  delete[] m_stiffnessMatrixPointers;

  // free memory of the memory allocate
  m_memoryAllocator.freeMemory();
}

void seissol::initializers::MemoryManager::initializeGlobalMatrix(          bool                       i_sparse,
                                                                   unsigned int                        i_numberOfRows,
                                                                   unsigned int                        i_numberOfColumns,
                                                                      const std::vector<unsigned int> &i_rows,
                                                                      const std::vector<unsigned int> &i_columns,
                                                                      const std::vector<double>       &i_values,
                                                                            double*                    o_matrix ) {
  // assert matching dimensions
  assert( i_rows.size()    == i_columns.size() );
  assert( i_columns.size() == i_values.size()  );

  // easy case: just write the values one after another
  if( i_sparse ) {
    for( unsigned int l_entry = 0; l_entry < i_values.size(); l_entry++) {
      o_matrix[l_entry] = i_values[l_entry];
    }
  }
  // dense matrix: set everything to zero and set only nonzeros
  else {
    // set everything to zero
    std::fill( o_matrix, o_matrix+i_numberOfRows*i_numberOfColumns, 0 );

    // iterate over nonzeros
    for( unsigned int l_entry = 0; l_entry < i_values.size(); l_entry++) {
      // index calculation (counting in XML starts at 1)
      unsigned int l_row    = i_rows[l_entry]    - 1;
      unsigned int l_column = i_columns[l_entry] - 1;

      // jump over columns
      unsigned int l_position = l_column * i_numberOfRows;
      // jump over rows
      l_position += l_row; 

      // set the nonzero
      o_matrix[l_position] = i_values[l_entry];
    } 
  }
}

void seissol::initializers::MemoryManager::initializeGlobalMatrices( const seissol::XmlParser &i_matrixReader ) {
  /*
   * read the global matrices
   */
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

  // read the flux matrices
  i_matrixReader.readGlobalMatrices( "flux",
                                     l_matrixIds,  l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );
  // assert we have all flux matrices
  assert( l_matrixIds.size() == 52 );

  // read the stiffness matrices
  i_matrixReader.readGlobalMatrices( "stiffness",
                                     l_matrixIds,  l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );

  // assert we have all stiffness matrices
  assert( l_matrixIds.size() == 58 );

  /*
   * set up memory chunks
   */
  std::vector<unsigned long long> l_chunkSizes;

  for( int l_matrix = 0; l_matrix < l_matrixIds.size(); l_matrix++) {
    // assert our matrix matches compiled number of basis functions
#ifndef NDEBUG
    assert( l_matrixNumberOfRows[l_matrix]    == NUMBEROFBASISFUNCTIONS );
    assert( l_matrixNumberOfColumns[l_matrix] == NUMBEROFBASISFUNCTIONS );
#endif

    if( l_matrixSparsities[l_matrix] ) {
      logDebug() << "matrix " << l_matrixNames[l_matrix] << " is assumed to be sparse.";

      // add a chunk, large enough for nonzeros
      l_chunkSizes.push_back( l_matrixValues[l_matrix].size() );
    }
    else {
      logDebug() << "matrix " << l_matrixNames[l_matrix] << " is assumed to be dense.";

      // add a chunk, large enough for the complete dense matrix
      l_chunkSizes.push_back( l_matrixNumberOfRows[l_matrix] * l_matrixNumberOfColumns[l_matrix] );
    }
  }

  //! pointer to allocated memory
  void *l_startOfAllocatedMemory = NULL;
  void **l_alignedMemoryPointers = NULL;

  // allocate 64 byte aligned memory
  m_memoryAllocator.allocateMemory( l_chunkSizes, sizeof(double),  64,
                                   &l_startOfAllocatedMemory,
                                   &l_alignedMemoryPointers );
  /* 
   * Initialize flux matrices.
   *
   * 1. reorder the pointers according to the ids
   *      assume matrices A, B, C with ids 0, 1, 2 are stored in the order B, A, C
   *      inside the XML-file, than the memory is allocated in the order of the
   *      XML-file:
   * <pre>
   *             B                     A                  C
   *      [*************, ------, **********, -----, ***********]
   *                        _                   _
   *                       /|\                 /|\
   *                        |                   |
   *                    padding for     ________|
   *                   aligned memory
   *                      access
   *
   *
   *      We store memory addresses of the chunks in the order of the Ids, thus access
   *      to the the correct flux matrix is a simple index calculation.
   *
   *             B                     A                  C
   *      [*************, ------, **********, -----, ***********]
   *       _                      _                  _
   *      /|\                    /|\                /|\
   *       |      matrix          |                  |
   *       |      pointers:       |                  |
   *       |         0------------/                  |
   *       \---------1                               |
   *                 2-------------------------------/
   *
   * </pre>
   *
   * 2. Initialize the matrices with proper values.
   *      For dense matrix kernels zeros are added.
   */
  for( int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++) {
    // reorder this matrix
    unsigned int l_matrixId = l_matrixIds[l_fluxMatrix];

    // store address of this flux matrices
    m_fluxMatrixPointers[l_matrixId] = (double*)l_alignedMemoryPointers[l_fluxMatrix];

    // initialize the stiffness matrix
    initializeGlobalMatrix( l_matrixSparsities[l_fluxMatrix],
                            l_matrixNumberOfRows[l_fluxMatrix],
                            l_matrixNumberOfColumns[l_fluxMatrix],
                            l_matrixRows[l_fluxMatrix],
                            l_matrixColumns[l_fluxMatrix],
                            l_matrixValues[l_fluxMatrix],
                            m_fluxMatrixPointers[l_matrixId] );
  }

  /*
   * Initialize stiffness matrices.
   */
  for( int l_stiffnessMatrix = 0; l_stiffnessMatrix < 6; l_stiffnessMatrix++ ) {
    // jump over the 52 flux matrices
    int l_globalMatrix = l_stiffnessMatrix + 52;

    // store address of this stiffness matrix
    m_stiffnessMatrixPointers[l_stiffnessMatrix] = (double*)l_alignedMemoryPointers[l_globalMatrix];

    // initialize the stiffness matrix
    initializeGlobalMatrix( l_matrixSparsities[l_globalMatrix],
                            l_matrixNumberOfRows[l_globalMatrix],
                            l_matrixNumberOfColumns[l_globalMatrix],
                            l_matrixRows[l_globalMatrix],
                            l_matrixColumns[l_globalMatrix],
                            l_matrixValues[l_globalMatrix],
                            m_stiffnessMatrixPointers[l_stiffnessMatrix] );
  }

  // set overlapping pointers for the prefetches
  m_fluxMatrixPointers[52] = m_stiffnessMatrixPointers[0];
  m_stiffnessMatrixPointers[6] = m_fluxMatrixPointers[0];
}

double** seissol::initializers::MemoryManager::getFluxMatrixPointers() const {
  return m_fluxMatrixPointers;
}

double** seissol::initializers::MemoryManager::getStiffnessMatrixPointers() const {
  return m_stiffnessMatrixPointers;
}
