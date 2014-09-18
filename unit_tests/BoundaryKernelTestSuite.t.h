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
 * Test suite, which tests the boundary kernel.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"
#include "Initializer/preProcessorMacros.fpp"

#include "../src/Boundary.h"
#include "DenseMatrix.hpp"
#include "SimpleBoundaryIntegrator.hpp"

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
    void testBoundaryKernel() {
      //! path to matrix set-up
      std::string l_matricesPath = m_matricesDirectory + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      //! matrix of degrees of freedom
      double l_degreesOfFreedom[  NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));
      double l_degreesOfFreedomUT[NUMBEROFBASISFUNCTIONS           *NUMBER_OF_QUANTITIES] = {0};

      //! matrices of time integrated unknowns
      double l_timeIntegrated[            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));
      double l_timeIntegratedNeighbors[4][NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));

      double l_timeIntegratedUT[            NUMBEROFBASISFUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));
      double l_timeIntegratedNeighborsUT[4][NUMBEROFBASISFUNCTIONS*NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));

      //! flux solvers matrices (positive eigenvalues): \f$ N_{k,i} A_k^+ N_{k,i}^{-1} \f$
      double l_fluxSolversPos[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! flux solvers matrices (negative eigenvalues): \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      double l_fluxSolversNeg[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! boundary conditions
      int l_boundaryConditions[4] = {0, 0, 0, 0};

      //! neighboring indices
      int l_neighboringIndices[4][2];

      //! flux matrices
      real *l_fluxMatrices[52];
      real l_fluxMatricesData[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*52] __attribute__((aligned(ALIGNMENT)));
      for( unsigned int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++ ) {
        l_fluxMatrices[l_fluxMatrix] = &l_fluxMatricesData[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*l_fluxMatrix];
      }

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // initialize the flux matrices
      // intitialize stiffness matrices
      for( unsigned int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++ ) {
        m_denseMatrix.initializeMatrix( l_fluxMatrix,
                                        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                        l_fluxMatrices[l_fluxMatrix] );
      }

      // construct boundary integrators
      unit_test::SimpleBoundaryIntegrator l_simpleBoundaryIntegrator;
      l_simpleBoundaryIntegrator.initialize( l_matricesPath );

      seissol::kernels::Boundary l_boundaryKernel;

      // repeat the test
      for( int l_repeat = 0; l_repeat < l_numberOfRepeats; l_repeat++) {

        // initialize DOFs
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_degreesOfFreedomUT );

        // synchronize degrees of freedom
        m_denseMatrix.copySubMatrix( l_degreesOfFreedomUT,
                                     NUMBEROFBASISFUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBEROFBASISFUNCTIONS,
                                     l_degreesOfFreedom,
                                     NUMBEROFBASISFUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );

        // intialize time integrated DOFs and flux solvers
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_timeIntegratedUT );


        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                         l_timeIntegratedNeighborsUT[l_face] );

          m_denseMatrix.initializeMatrix( 52,
                                          NUMBEROFVARIABLES,
                                          NUMBEROFVARIABLES,
                                          l_fluxSolversPos[l_face] );
          m_denseMatrix.initializeMatrix( 52,
                                          NUMBEROFVARIABLES,
                                          NUMBEROFVARIABLES,
                                          l_fluxSolversNeg[l_face] );

          // initialize tet orientations
          l_neighboringIndices[l_face][0] = rand() % 4;
          l_neighboringIndices[l_face][1] = rand() % 3;
        }

        // synchronize time integrated DOFs
        m_denseMatrix.copySubMatrix( l_timeIntegratedUT,
                                     NUMBEROFBASISFUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBEROFBASISFUNCTIONS,
                                     l_timeIntegrated,
                                     NUMBEROFBASISFUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          m_denseMatrix.copySubMatrix( l_timeIntegratedNeighborsUT[l_face],
                                       NUMBEROFBASISFUNCTIONS,
                                       NUMBER_OF_QUANTITIES,
                                       NUMBEROFBASISFUNCTIONS,
                                       l_timeIntegratedNeighbors[l_face],
                                       NUMBEROFBASISFUNCTIONS,
                                       NUMBER_OF_QUANTITIES,
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
        }

        // do a simple boundary intagration
        l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUT,
                                                               l_timeIntegratedNeighborsUT,
                                                               l_boundaryConditions,
                                                               l_neighboringIndices,
                                                               l_fluxSolversPos,
                                                               l_fluxSolversNeg,
                                                               l_degreesOfFreedomUT );

        // do boundary integration with generated kernels
        l_boundaryKernel.computeLocalIntegral(     l_boundaryConditions,
                                                   l_fluxMatrices,
                                                   l_timeIntegrated,
                                                   l_fluxSolversPos,
                                                   l_degreesOfFreedom );

        l_boundaryKernel.computeNeighborsIntegral( l_boundaryConditions,
                                                   l_neighboringIndices,
                                                   l_fluxMatrices,
                                                   l_timeIntegratedNeighbors,
                                                   l_fluxSolversNeg,
                                                   l_degreesOfFreedom );

        // check the results
        m_denseMatrix.checkSubMatrix( l_degreesOfFreedom,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_degreesOfFreedomUT,
                                      NUMBEROFBASISFUNCTIONS,
                                      NUMBER_OF_QUANTITIES );

        // test absorbing boundary conditions
        for( int l_face = 0; l_face < 4; l_face++ ) {
          // use absorbing boundary conditions at this face
          l_boundaryConditions[l_face] = 5;
          
          // do a simple boundary intagration
          l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUT,
                                                                 l_timeIntegratedNeighborsUT,
                                                                 l_boundaryConditions,
                                                                 l_neighboringIndices,
                                                                 l_fluxSolversPos,
                                                                 l_fluxSolversNeg,
                                                                 l_degreesOfFreedomUT );

          // do boundary integration with optimized matrix kernels
          l_boundaryKernel.computeLocalIntegral(     l_boundaryConditions,
                                                     l_fluxMatrices,
                                                     l_timeIntegrated,
                                                     l_fluxSolversPos,
                                                     l_degreesOfFreedom );

          l_boundaryKernel.computeNeighborsIntegral( l_boundaryConditions,
                                                     l_neighboringIndices,
                                                     l_fluxMatrices,
                                                     l_timeIntegratedNeighbors,
                                                     l_fluxSolversNeg,
                                                     l_degreesOfFreedom );

          // check the results
          m_denseMatrix.checkSubMatrix( l_degreesOfFreedom,
                                        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                        NUMBER_OF_QUANTITIES,
                                        l_degreesOfFreedomUT,
                                        NUMBEROFBASISFUNCTIONS,
                                        NUMBER_OF_QUANTITIES );
        }

      }

    }
};
