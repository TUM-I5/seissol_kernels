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
 * Boundary integration in SeisSol.
 **/
#include "BoundaryIntegrator.h"

void seissol::BoundaryIntegrator::setUpMatrixKernel( unsigned int i_id,
                                                     bool i_sparse ) {
  // assert we are not out of bounds
  assert( i_id < 53 );

  // include the generated initialization code
#include <generated_code/initialization/boundary_matrix_kernels.hpp_include>

}

seissol::BoundaryIntegrator::BoundaryIntegrator( const seissol::XmlParser &i_matrixReader,
                                                 const seissol::initializers::MemoryManager &i_memoryManager ) {
  /*
   * Initialize flux matrix pointers.
   */
  m_fluxMatrixPointers = i_memoryManager.getFluxMatrixPointers();

  /*
   * Initialize function pointers to the matrix kernels
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

  // assert we have all matrices
  assert( l_matrixIds.size() == 52 );

  for( int l_i = 0; l_i < 52; l_i++) {
    unsigned int l_matrixId = l_matrixIds[l_i];

    setUpMatrixKernel( l_matrixId,
                       l_matrixSparsities[l_i] );
  }

  // read the jacobi matrix
  i_matrixReader.readLocalMatrixStructure( "fluxSolver",
                                           l_matrixIds,
                                           l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                           l_matrixRows, l_matrixColumns );

  // assert we found the jacobi matrix
  assert( l_matrixIds.size() == 53 );
  assert( l_matrixIds.back() == 52 );

  // set up this kernels
  setUpMatrixKernel( l_matrixIds.back(),
                     l_matrixSparsities[l_matrixIds.back()] );
}

void seissol::BoundaryIntegrator::computeBoundaryIntegral(                double i_timeIntegratedUnknownsElement[2][NUMBEROFUNKNOWNS],
                                                                          double i_timeIntegratedUnknownsNeighbors[4][NUMBEROFUNKNOWNS],
                                                           const unsigned int    i_boundaryConditions[4],
                                                           const unsigned int    i_neighboringIndices[4][2],
                                                                          double i_nApNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                                                          double i_nAmNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                                                          double io_unknowns[NUMBEROFUNKNOWNS] ){
  /*
   * Assert alignments, which are assumed in the matrix kernels.
   */
#if NUMBEROFBASISFUNCTIONS == 1
  // 8 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   %  8 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) %  8 == 0 );
  assert( ((uintptr_t)io_unknowns)                       %  8 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 4
  // 32 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   % 32 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) % 32 == 0 );
  assert( ((uintptr_t)io_unknowns)                       % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 10
  // 16 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   % 16 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) % 16 == 0 );
  assert( ((uintptr_t)io_unknowns)                       % 16 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 20
  // 32 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   % 32 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) % 32 == 0 );
  assert( ((uintptr_t)io_unknowns)                       % 32 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 35
  // 8 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   %  8 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) %  8 == 0 );
  assert( ((uintptr_t)io_unknowns)                       %  8 == 0 );
#elif NUMBEROFBASISFUNCTIONS == 56
  // 64 byte alignment
  assert( ((uintptr_t)i_timeIntegratedUnknownsElement)   % 64 == 0 );
  assert( ((uintptr_t)i_timeIntegratedUnknownsNeighbors) % 64 == 0 );
  assert( ((uintptr_t)io_unknowns)                       % 64 == 0 );
#else
#error Preprocessor flag NUMBEROFBASISFUNCTIONS is not in {1, 4, 10, 20, 35, 56}.
#endif

  /*
   * Computation
   */
  //! temporary product (we have to multiply a matrix from the left and the right, 64 byte aligned)
  double l_temporaryProduct[ NUMBEROFUNKNOWNS ] __attribute__((aligned(64)));

  // iterate over local faces \f$i\f$ and compute the elements local contribution and the neighboring elements contribution
  for( unsigned int l_localFace = 0; l_localFace < 4; l_localFace++) {
    /*
     * Compute the neighboring elements flux matrix id.
     */
    // id of the flux matrix (0-3: element local, 4-51: neighboring element)
    unsigned int l_id;

    if( i_boundaryConditions[l_localFace] != 1 ) {
      // derive memory and kernel index
      // TODO: Precompute this to get rid of the if's in the boundary cases.
      l_id = 4                                        // jump over flux matrices \f$ F^{-, i} \f$
            + l_localFace*12                          // jump over index \f$i\f$
            + i_neighboringIndices[l_localFace][0]*3  // jump over index \f$j\f$
            + i_neighboringIndices[l_localFace][1];   // jump over index \f$h\f$

      // assert we have a neighboring index in the case of non-absorbing boundary conditions.
      assert( l_id >= 4 || i_boundaryConditions[l_localFace] == 5 );
    }
    else { // fall back to local matrices in case of free surface boundary conditions
      l_id = l_localFace;
    }

    /*
     * Element local contribution.
     */
    // set temporary product to zero
    memset(l_temporaryProduct, 0, NUMBEROFUNKNOWNS*sizeof(*l_temporaryProduct));

    // compute element local contribution
    m_matrixKernels[l_localFace]( m_fluxMatrixPointers[l_localFace],  i_timeIntegratedUnknownsElement[0],             l_temporaryProduct,
                                  l_temporaryProduct,                 i_nApNm1[l_localFace],                          io_unknowns         ); // prefetches
    
    m_matrixKernels[52](          l_temporaryProduct,                 i_nApNm1[l_localFace],                          io_unknowns,
                                  m_fluxMatrixPointers[l_id],         i_timeIntegratedUnknownsNeighbors[l_localFace], l_temporaryProduct  ); // prefetches

    /*
     * Neighboring elements contribution.
     */
    // no neighboring element contribution in the case of absorbing boundary conditions
    if( i_boundaryConditions[l_localFace] != 5 ) {
      // set temporary product to zero
      memset(l_temporaryProduct, 0, NUMBEROFUNKNOWNS*sizeof(*l_temporaryProduct));    

      // assert we have a valid index.
      assert( l_id < 52 );

      // compute neighboring elements contribution
      m_matrixKernels[l_id](      m_fluxMatrixPointers[l_id],         i_timeIntegratedUnknownsNeighbors[l_localFace], l_temporaryProduct,
                                  l_temporaryProduct,                 i_nAmNm1[l_localFace],                          io_unknowns         ); // prefetches

      /*
       * Compute the id of the global matrix used in the matrix kernel following the next operator.
       * If the matrix operation is part of the boundary intgrator, that's just the element local flux matrix of next face.
       * Otherwise - l_localFace/3 == 1 - the first stiffness matrix of the volume integrator \f$ M^{-1} K^\xi \f$ for the next element \f$k+1\f$ is next.
       * \f$ M^{-1} K^\xi \f$ can be acessed as the 53rd flux matrix (l_upcomingId == 52). 
       */
      int l_upcomingId = l_localFace + 1 + ( l_localFace / 3 ) * 48;

      /*
       * Prefetch either the next flux matrix and the time integrated unknowns of this element or if this is the last operation - local face \f$i\f$ is four:
       * Prefetch the the first stiffness matrix of the volume integrator \f$ M^{-1} K^\xi \f$ and the time integrated unknowns of the next element.
       */
      m_matrixKernels[52](        l_temporaryProduct,                 i_nAmNm1[l_localFace],                          io_unknowns,
                                  m_fluxMatrixPointers[l_upcomingId], i_timeIntegratedUnknownsElement[l_localFace/3], l_temporaryProduct  ); // prefetches
    }
  }
}
