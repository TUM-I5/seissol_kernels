/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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
 * Simple boundary integrator.
 **/

#ifndef SIMPLEBOUNDARYINTEGRATOR_HPP_
#define SIMPLEBOUNDARYINTEGRATOR_HPP_

#include <Initializer/typedefs.hpp>

#include "DenseMatrix.hpp"

namespace unit_test {
  class SimpleBoundaryIntegrator;
};

/**
 * Simple implementation of the boundary integration in SeisSol.
 **/
class unit_test::SimpleBoundaryIntegrator {
  //private:
    //! flux matrices (elements contribution): \f$ F^{-, i} \f$
    real m_fluxMatricesNeg[4][NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_BASIS_FUNCTIONS];

    //! flux matrices (neighboring elements contribution): \f$ F^+{+, i, j, h} \f$
    real m_fluxMatricesPos[48][NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_BASIS_FUNCTIONS];

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;
  public:
    /**
     * Initialize the matrices.
     *
     * @param i_matricesPath path to the matrix xml file.
     **/
    void initialize( std::string &i_matricesPath ) {
      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( i_matricesPath );

      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        m_denseMatrix.initializeMatrix( l_face,
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        m_fluxMatricesNeg[l_face] );
      } 
     
      for( unsigned int l_id = 0; l_id < 48; l_id++) {
        m_denseMatrix.initializeMatrix( 4+l_id,
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        m_fluxMatricesPos[l_id] );
      }
    }
    
    /**
     * Computes a simple boundary integration.
     *
     * @param i_timeIntegratedUnknownsElement time integrated unknowns of the element.
     * @param i_timeIntegratedUnknownsNeighbors time integrated unknwons of the neighboring elements.
     * @param i_faceTypes types of the faces.
     * @param i_neighboringIndices oriantation of the element in relation to the neighboring element (ref. element).
     * @param i_nApNm1DivM flux solvers for elements contribution.
     * @param i_nAmNm1DivM flux solvers for the neighboring elements contributions.
     * @param io_unknowns degrees of freedom of the element, which are updated.
     **/
    void computeBoundaryIntegration(       real          i_timeIntegratedUnknownsElement[NUMBER_OF_DOFS],
                                           real          i_timeIntegratedUnknownsNeighbors[4][NUMBER_OF_DOFS],
                                     const enum faceType i_faceTypes[4],
                                     const int           i_neighboringIndices[4][2],
                                           real          i_nApNm1DivM[4][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                                           real          i_nAmNm1DivM[4][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                                           real          io_unknowns[NUMBER_OF_DOFS] ) {
      //! temporary product (we have to multiply a matrix from the left and the right)
      real l_temporaryProduct[ NUMBER_OF_DOFS ];


      // compute the elements contribution
      for( unsigned int l_localFace = 0; l_localFace < 4; l_localFace++) {
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBER_OF_DOFS, 0 );

        m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
                                                     m_fluxMatricesNeg[l_localFace], i_timeIntegratedUnknownsElement, l_temporaryProduct );


        m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
                                                     l_temporaryProduct, i_nApNm1DivM[l_localFace], io_unknowns );
      }

      // compute the neighboring elements contribution
      for( unsigned int l_localFace = 0; l_localFace < 4; l_localFace++) {
        // no contribution of the neighboring element in case of absorbing boundary conditions
        if( i_faceTypes[l_localFace] != outflow ) {
          // set temporary product to zero
          std::fill( l_temporaryProduct, l_temporaryProduct+NUMBER_OF_DOFS, 0 );

          // compute matrix index
          unsigned l_fluxIndex = l_localFace * 12
                               + i_neighboringIndices[l_localFace][0] * 3
                               + i_neighboringIndices[l_localFace][1];

          TS_ASSERT( l_fluxIndex < 48 );

          m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
                                                       m_fluxMatricesPos[l_fluxIndex], i_timeIntegratedUnknownsNeighbors[l_localFace], l_temporaryProduct );

          m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
                                                       l_temporaryProduct, i_nAmNm1DivM[l_localFace], io_unknowns );
        }
      }
    }
};
#endif
