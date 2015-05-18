/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
#include <string>
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

    //! #repeats for random unit tests
    int m_numberOfRepeats;

    //! dense matrix functionality
    DenseMatrix m_denseMatrix;

  public:
    void setUp() {
      // get number of reapeats
      m_numberOfRepeats = m_configuration.getNumberOfRepeats();

      // get unit tests src directory
      m_unitTestsSrcDirectory = m_configuration.getUnitTestsSrcDirectory();

      // generate random seed
      srand(time(NULL));
    }

    /**
     * Tests the complete kernel, which does the boundary integration.
     **/
    void testBoundaryKernel() {
      //! path to matrix set-up
      std::string l_matricesPath = m_configuration.getMatricesFile();

      //! matrix of degrees of freedom
      real l_degreesOfFreedom[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
      real l_degreesOfFreedomUT[NUMBER_OF_DOFS] = {0};

      //! matrices of time integrated unknowns
      real l_timeIntegrated[            NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
      real l_timeIntegratedNeighbors[4][NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
      real* l_timeIntegratedNeighborsPT[4];
      l_timeIntegratedNeighborsPT[0] = l_timeIntegratedNeighbors[0]; l_timeIntegratedNeighborsPT[1] = l_timeIntegratedNeighbors[1];
      l_timeIntegratedNeighborsPT[2] = l_timeIntegratedNeighbors[2]; l_timeIntegratedNeighborsPT[3] = l_timeIntegratedNeighbors[3];

      real l_timeIntegratedUT[            NUMBER_OF_DOFS] __attribute__((aligned(ALIGNMENT)));
      real l_timeIntegratedNeighborsUT[4][NUMBER_OF_DOFS] __attribute__((aligned(ALIGNMENT)));

      //! flux solvers matrices (positive eigenvalues): \f$ N_{k,i} A_k^+ N_{k,i}^{-1} \f$
      real l_fluxSolversPos[4][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];

      //! flux solvers matrices (negative eigenvalues): \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      real l_fluxSolversNeg[4][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];

      //! face types
      enum faceType l_faceTypes[4] = {regular, regular, regular, regular};

      //! neighboring indices
      int l_neighboringIndices[4][2];

      //! flux matrices
      real* l_fluxMatrices[52*2];
      real* l_fluxMatricesData = (real*) _mm_malloc( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*52*2*sizeof(real), ALIGNMENT );

      for( unsigned int l_fluxMatrix = 0; l_fluxMatrix < 52*2; l_fluxMatrix++ ) {
        l_fluxMatrices[l_fluxMatrix] = l_fluxMatricesData + (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*l_fluxMatrix);
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

        if( m_configuration.getSparseSwitch( l_fluxMatrix ) == true ) {
          m_denseMatrix.copyDenseToSparse( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                           NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                           l_fluxMatrices[l_fluxMatrix],
                                           l_fluxMatrices[l_fluxMatrix+52] );
        }
        else {
          m_denseMatrix.copySubMatrix( l_fluxMatrices[l_fluxMatrix],
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                       l_fluxMatrices[l_fluxMatrix+52],
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
        }

      }

      // construct boundary integrators
      unit_test::SimpleBoundaryIntegrator l_simpleBoundaryIntegrator;
      l_simpleBoundaryIntegrator.initialize( l_matricesPath );

      seissol::kernels::Boundary l_boundaryKernel;

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_numberOfRepeats; l_repeat++) {
        // initialize DOFs
        m_denseMatrix.setRandomValues( NUMBER_OF_DOFS,
                                       l_degreesOfFreedomUT );

        // synchronize degrees of freedom
        m_denseMatrix.copySubMatrix( l_degreesOfFreedomUT,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     l_degreesOfFreedom,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );

        // intialize time integrated DOFs and flux solvers
        m_denseMatrix.setRandomValues( NUMBER_OF_DOFS,
                                       l_timeIntegratedUT );


        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          m_denseMatrix.setRandomValues( NUMBER_OF_DOFS,
                                         l_timeIntegratedNeighborsUT[l_face] );

          m_denseMatrix.initializeMatrix( 52,
                                          NUMBER_OF_QUANTITIES,
                                          NUMBER_OF_QUANTITIES,
                                          l_fluxSolversPos[l_face] );
          m_denseMatrix.initializeMatrix( 52,
                                          NUMBER_OF_QUANTITIES,
                                          NUMBER_OF_QUANTITIES,
                                          l_fluxSolversNeg[l_face] );

          // initialize tet orientations
          l_neighboringIndices[l_face][0] = rand() % 4;
          l_neighboringIndices[l_face][1] = rand() % 3;
        }

        // synchronize time integrated DOFs
        m_denseMatrix.copySubMatrix( l_timeIntegratedUT,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     l_timeIntegrated,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          m_denseMatrix.copySubMatrix( l_timeIntegratedNeighborsUT[l_face],
                                       NUMBER_OF_BASIS_FUNCTIONS,
                                       NUMBER_OF_QUANTITIES,
                                       NUMBER_OF_BASIS_FUNCTIONS,
                                       l_timeIntegratedNeighbors[l_face],
                                       NUMBER_OF_BASIS_FUNCTIONS,
                                       NUMBER_OF_QUANTITIES,
                                       NUMBER_OF_ALIGNED_BASIS_FUNCTIONS );
        }

        // do a simple boundary intagration
        l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUT,
                                                               l_timeIntegratedNeighborsUT,
                                                               l_faceTypes,
                                                               l_neighboringIndices,
                                                               l_fluxSolversPos,
                                                               l_fluxSolversNeg,
                                                               l_degreesOfFreedomUT );

        // do boundary integration with generated kernels
        l_boundaryKernel.computeLocalIntegral(     l_faceTypes,
                                                   l_fluxMatrices+52,
                                                   l_timeIntegrated,
                                                   l_fluxSolversPos,
                                                   l_degreesOfFreedom );

        l_boundaryKernel.computeNeighborsIntegral( l_faceTypes,
                                                   l_neighboringIndices,
                                                   l_fluxMatrices+52,
                                                   l_timeIntegratedNeighborsPT,
                                                   l_fluxSolversNeg,
                                                   l_degreesOfFreedom );

        // check the results
        m_denseMatrix.checkSubMatrix( l_degreesOfFreedom,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_degreesOfFreedomUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );

        // test absorbing boundary conditions
        for( int l_face = 0; l_face < 4; l_face++ ) {
          // use absorbing boundary conditions at this face
          l_faceTypes[l_face] = outflow;
          
          // do a simple boundary intagration
          l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUT,
                                                                 l_timeIntegratedNeighborsUT,
                                                                 l_faceTypes,
                                                                 l_neighboringIndices,
                                                                 l_fluxSolversPos,
                                                                 l_fluxSolversNeg,
                                                                 l_degreesOfFreedomUT );

          // do boundary integration with optimized matrix kernels
          l_boundaryKernel.computeLocalIntegral(     l_faceTypes,
                                                     l_fluxMatrices+52,
                                                     l_timeIntegrated,
                                                     l_fluxSolversPos,
                                                     l_degreesOfFreedom );

          l_boundaryKernel.computeNeighborsIntegral( l_faceTypes,
                                                     l_neighboringIndices,
                                                     l_fluxMatrices+52,
                                                     l_timeIntegratedNeighborsPT,
                                                     l_fluxSolversNeg,
                                                     l_degreesOfFreedom );

          // check the results
          m_denseMatrix.checkSubMatrix( l_degreesOfFreedom,
                                        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                        NUMBER_OF_QUANTITIES,
                                        l_degreesOfFreedomUT,
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        NUMBER_OF_QUANTITIES );
        }
      }

    // free memory
    _mm_free( l_fluxMatricesData );
    }

    /**
     * Tests the local flop-counters.
     **/
    void testLocalFlops() {
      unsigned int   l_nonZeroFlops  = 0;
      unsigned int   l_hardwareFlops = 0;
      enum faceType  l_faceTypes[4];

      // volume kernel
      seissol::kernels::Boundary l_boundaryKernel;

      /*
       * sparse <= dense
       */
      // get the flops for a single call
      l_faceTypes[0] = l_faceTypes[1] = l_faceTypes[2] = l_faceTypes[3] = regular;
      l_boundaryKernel.flopsLocalIntegral( l_faceTypes, l_nonZeroFlops, l_hardwareFlops );

      // assert we are at least doing the sparse flops in hardware
      TS_ASSERT_LESS_THAN_EQUALS( l_nonZeroFlops, l_hardwareFlops );

      /*
       * hardware <= dense
       */
      unsigned int l_denseFlops =  (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES * 2) * 4;

      // flux solvers
      l_denseFlops += ( NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * 2 ) * 4;
      TS_ASSERT_LESS_THAN_EQUALS( l_hardwareFlops, l_denseFlops );

      /*
       * dr == 0
       */
      // get the flops for an DR-only setting
      l_faceTypes[0] = l_faceTypes[1] = l_faceTypes[2] = l_faceTypes[3] = dynamicRupture;
      l_boundaryKernel.flopsLocalIntegral( l_faceTypes, l_nonZeroFlops, l_hardwareFlops );

      // assert we don't count any flops as they are done externally
      TS_ASSERT_EQUALS( 0, l_nonZeroFlops );
      TS_ASSERT_EQUALS( 0, l_hardwareFlops );
    }

    /**
     * Tests the neighboring flop-counters.
     **/
    void testNeighboringFlops() {
      unsigned int   l_nonZeroFlops  = 0;
      unsigned int   l_hardwareFlops = 0;
      enum faceType  l_faceTypes[4];
      int   l_neighboringIndices[4][2];

      // volume kernel
      seissol::kernels::Boundary l_boundaryKernel;

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_numberOfRepeats; l_repeat++) {
        // set up neighboring indices
        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          l_neighboringIndices[l_face][0] = rand()%4;
          l_neighboringIndices[l_face][1] = rand()%3;
        }

        /*
         * sparse <= dense
         */
        // get the flops for a single call
        l_faceTypes[0] = l_faceTypes[1] = l_faceTypes[2] = l_faceTypes[3] = regular;

        l_boundaryKernel.flopsNeighborsIntegral( l_faceTypes, l_neighboringIndices, l_nonZeroFlops, l_hardwareFlops );

        // assert we are at least doing the sparse flops in hardware
        TS_ASSERT_LESS_THAN_EQUALS( l_nonZeroFlops, l_hardwareFlops );

        /*
         * hardware <= dense
         */
        unsigned int l_denseFlops =  (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES * 2) * 4;

        // flux solvers
        l_denseFlops += ( NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * 2 ) * 4;
        TS_ASSERT_LESS_THAN_EQUALS( l_hardwareFlops, l_denseFlops );

        /*
         * outflow == 0
         */
        // get the flops for an outflow-only setting
        l_faceTypes[0] = l_faceTypes[1] = l_faceTypes[2] = l_faceTypes[3] = outflow;
        l_boundaryKernel.flopsNeighborsIntegral( l_faceTypes, l_neighboringIndices, l_nonZeroFlops, l_hardwareFlops );

        // assert we don't count any flops
        TS_ASSERT_EQUALS( 0, l_nonZeroFlops );
        TS_ASSERT_EQUALS( 0, l_hardwareFlops );
      }
    }
};
