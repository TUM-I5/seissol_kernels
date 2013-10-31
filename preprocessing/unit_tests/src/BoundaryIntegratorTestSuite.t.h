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
 * Test suite, which tests the flux matrices.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"
#include "Initializer/preProcessorMacros.fpp"

#define private public
#include "../seissol_src/Solver/kernels/BoundaryIntegrator.h"
#undef private
#include "DenseMatrix.hpp"

namespace unit_test {
  class BoundaryIntegratorTestSuite;
}

class unit_test::BoundaryIntegratorTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! path to the source directory of the unit tests
    std::string m_unitTestsSrcDirectory;

    //! path to the matrice directory
    std::string m_matricesDirectory;

    //! #repeats for random unit tests
    int l_numberOfRepeats;

    //! dense matrix functionality
    DenseMatrix m_denseMatrix;

    // simple boundary integration to test again
    void simpleBoundaryIntegration(                double i_timeIntegratedUnknownsElement[NUMBEROFUNKNOWNS],
                                                   double i_timeIntegratedUnknownsNeighbors[4][NUMBEROFUNKNOWNS],
                                    const unsigned int    i_boundaryConditons[4],
                                    const unsigned int    i_neighboringIndices[4][2],
                                                   double i_nApNm1DivM[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                                   double i_nAmNm1DivM[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                                   double i_fluxMatricesNeg[4][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS],
                                                   double i_fluxMatricesPos[48][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS],
                                                   double io_unknowns[NUMBEROFUNKNOWNS] ) {
      //! temporary product (we have to multiply a matrix from the left and the right)
      double l_temporaryProduct[ NUMBEROFUNKNOWNS ];


      // compute the elements contribution
      for( unsigned int l_localFace = 0; l_localFace < 4; l_localFace++) {
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0 );

        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                     i_fluxMatricesNeg[l_localFace], i_timeIntegratedUnknownsElement, l_temporaryProduct );


        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                     l_temporaryProduct, i_nApNm1DivM[l_localFace], io_unknowns );
      }

      // compute the neighboring elements contribution
      for( unsigned int l_localFace = 0; l_localFace < 4; l_localFace++) {
        // no contribution of the neighboring element in case of absorbing boundary conditions
        if( i_boundaryConditons[l_localFace] != 5 ) { 
          // set temporary product to zero
          std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0 );

          // compute matrix index
          unsigned l_fluxIndex = l_localFace * 12
                               + i_neighboringIndices[l_localFace][0] * 3
                               + i_neighboringIndices[l_localFace][1];

          TS_ASSERT( l_fluxIndex < 48 );

          m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                       i_fluxMatricesPos[l_fluxIndex], i_timeIntegratedUnknownsNeighbors[l_localFace], l_temporaryProduct );

          m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                       l_temporaryProduct, i_nAmNm1DivM[l_localFace], io_unknowns );
        }
      }

    }

  public:
    void setUp() {
      // get number of reapeats
      l_numberOfRepeats = m_configuration.getNumberOfRepeats();

      // get unit tests src directory
      m_unitTestsSrcDirectory = m_configuration.getUnitTestsSrcDirectory();

      m_matricesDirectory = m_configuration.getMatricesDirectory();

      // generate random seed
      srand(time(NULL));
    }

    /**
     * Tests the complete kernel, which does the boundary integration.
     **/
    void testBoundaryIntegrator() {
      //! path to matrix set-up
      std::string l_matricesPath = m_matricesDirectory + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      //! matrix of unknowns
      double l_unknownsUnitTests[NUMBEROFUNKNOWNS] = {0};
      double l_unknownsBoundaryIntegrator[NUMBEROFUNKNOWNS]         __attribute__((aligned(64)));

      //! matrices of time integrated unknowns
      double l_timeIntegratedUnknownsElement[NUMBEROFUNKNOWNS]      __attribute__((aligned(64)));
      double l_timeIntegratedUnknownsNeighbors[4][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

      //! flux solvers matrices (positive eigenvalues): \f$ N_{k,i} A_k^+ N_{k,i}^{-1} \f$
      double l_fluxSolversPos[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! flux solvers matrices (negative eigenvalues): \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      double l_fluxSolversNeg[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! flux matrices (elements contribution): \f$ F^{-, i} \f$
      double l_fluxMatricesNeg[4][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS];

      //! flux matrices (neighboring elements contribution): \f$ F^+{+, i, j, h} \f$
      double l_fluxMatricesPos[48][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS];

      //! boundary conditions
      unsigned int l_boundaryConditions[4] = {0, 0, 0, 0};

      //! neighboring indices
      unsigned int l_neighboringIndices[4][2];

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );

      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );

      // construct a boundary integrator
      seissol::BoundaryIntegrator l_boundaryIntegrator( l_matrixReader, l_memoryManager );

      // repeat the test
      for( int l_repeat = 0; l_repeat < l_numberOfRepeats; l_repeat++) {
        //
        // initialize the matrices
        //
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_unknownsUnitTests );

        // synchronize unknowns matrices 
        for( unsigned int l_element = 0; l_element < NUMBEROFUNKNOWNS; l_element++) {
          l_unknownsBoundaryIntegrator[l_element] = l_unknownsUnitTests[l_element];
        }

        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_timeIntegratedUnknownsElement );


        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                         l_timeIntegratedUnknownsNeighbors[l_face] );

          m_denseMatrix.initializeMatrix( 52,
                                          NUMBEROFVARIABLES,
                                          NUMBEROFVARIABLES,
                                          l_fluxSolversPos[l_face] );
          m_denseMatrix.initializeMatrix( 52,
                                          NUMBEROFVARIABLES,
                                          NUMBEROFVARIABLES,
                                          l_fluxSolversNeg[l_face] );

          m_denseMatrix.initializeMatrix( l_face,
                                          NUMBEROFBASISFUNCTIONS,
                                          NUMBEROFBASISFUNCTIONS,
                                          l_fluxMatricesNeg[l_face] );


          l_neighboringIndices[l_face][0] = rand() % 4;
          l_neighboringIndices[l_face][1] = rand() % 3;
        }

        for( unsigned int l_id = 0; l_id < 48; l_id++) {
          m_denseMatrix.initializeMatrix( 4+l_id,
                                          NUMBEROFBASISFUNCTIONS,
                                          NUMBEROFBASISFUNCTIONS,
                                          l_fluxMatricesPos[l_id] );
        }

        // do a simple boundary intagration
        simpleBoundaryIntegration( l_timeIntegratedUnknownsElement,
                                   l_timeIntegratedUnknownsNeighbors,
                                   l_boundaryConditions,
                                   l_neighboringIndices,
                                   l_fluxSolversPos,
                                   l_fluxSolversNeg,
                                   l_fluxMatricesNeg,
                                   l_fluxMatricesPos,
                                   l_unknownsUnitTests );

        // do boundary integration with generated kernels
        l_boundaryIntegrator.computeBoundaryIntegral( l_timeIntegratedUnknownsElement,
                                                      l_timeIntegratedUnknownsNeighbors,
                                                      l_boundaryConditions,
                                                      l_neighboringIndices,
                                                      l_fluxSolversPos,
                                                      l_fluxSolversNeg,
                                                      l_unknownsBoundaryIntegrator );

        // check the results
        m_denseMatrix.checkResult( NUMBEROFUNKNOWNS,l_unknownsUnitTests, l_unknownsBoundaryIntegrator );

        // test absorbing boundary conditions
        for( int l_face = 0; l_face < 4; l_face++ ) {
          // use absorbing boundary conditions at this face
          l_boundaryConditions[l_face] = 5;
          
          // do a simple boundary intagration
          simpleBoundaryIntegration( l_timeIntegratedUnknownsElement,
                                     l_timeIntegratedUnknownsNeighbors,
                                     l_boundaryConditions,
                                     l_neighboringIndices,
                                     l_fluxSolversPos,
                                     l_fluxSolversNeg,
                                     l_fluxMatricesNeg,
                                     l_fluxMatricesPos,
                                     l_unknownsUnitTests );

          // do boundary integration with generated kernels
          l_boundaryIntegrator.computeBoundaryIntegral( l_timeIntegratedUnknownsElement,
                                                        l_timeIntegratedUnknownsNeighbors,
                                                        l_boundaryConditions,
                                                        l_neighboringIndices,
                                                        l_fluxSolversPos,
                                                        l_fluxSolversNeg,
                                                        l_unknownsBoundaryIntegrator );

          // check the results
          m_denseMatrix.checkResult( NUMBEROFUNKNOWNS,l_unknownsUnitTests, l_unknownsBoundaryIntegrator );
        }
      }
    }
};
