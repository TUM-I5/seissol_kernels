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
 * Test suite, which tests the time kernel.
 **/

#include <cassert>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"

#include "../src/common.hpp"
#define private public
#include "../src/Time.h"
#undef private
#include "DenseMatrix.hpp"
#include "SimpleTimeIntegrator.hpp"

namespace unit_test {
  class TimeKernelTestSuite;
}

class unit_test::TimeKernelTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! Simple time integration functionality
    unit_test::SimpleTimeIntegrator m_simpleTimeIntegrator;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

    /**
     * Allocated aligned memory for the given data structures.
     *
     * @param o_stiffnessMatrices stiffness matrices.
     * @param o_degreesOfFreedom degrees of freedom.
     * @param o_starMatrices star matrices.
     * @param o_timeDerivatives time derivatives in compressed format.
     * @param o_timeIntegrated time integrated DOFs.
     * @param o_timeExtrapolated time extrapolated DOFs.
     **/
    void allocateMemory( real**& o_stiffnessMatrices,
                         real*&  o_degreesOfFreedom,
                         real*&  o_timeDerivatives,
                         real*&  o_timeIntegrated,
                         real*&  o_timeExtrapolated) {

      unsigned int l_alignedNumberOfBasisFunctions = seissol::kernels::getNumberOfAlignedBasisFunctions();
      unsigned int l_alignedTimeDerivativesSize    = seissol::kernels::Time::getAlignedTimeDerivativesSize();

      // allocate memory for the stiffness matrices
      unsigned l_sizesStiff[2];
      l_sizesStiff[0] = seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 );
      l_sizesStiff[1] = seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER );


      // allocate memory
      o_stiffnessMatrices    = (real**) malloc( 9*sizeof(real*) );
      for( unsigned int l_c = 0; l_c < 6; l_c++ ) {
        o_stiffnessMatrices[l_c] = (real*) _mm_malloc( l_sizesStiff[0]*l_sizesStiff[1]*sizeof(real),                  ALIGNMENT );
      }

      o_degreesOfFreedom     = (real*) _mm_malloc( l_alignedNumberOfBasisFunctions*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );

      o_timeDerivatives      = (real*) _mm_malloc( l_alignedTimeDerivativesSize*NUMBER_OF_QUANTITIES*sizeof(real),    ALIGNMENT );

      o_timeIntegrated       = (real*) _mm_malloc( l_alignedNumberOfBasisFunctions*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );
      o_timeExtrapolated     = (real*) _mm_malloc( l_alignedNumberOfBasisFunctions*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );
    }

    /**
     * Frees memory of the given data structures.
     *
     * @param o_stiffnessMatrices stiffness matrices.
     * @param o_degreesOfFreedom degrees of freedom.
     * @param o_starMatrices star matrices.
     * @param o_timeDerivatives time derivatives in compressed format.
     * @param o_timeIntegrated time integrated DOFs.
     * @param o_timeExtrapolated time extrapolated DOFs.
     **/
    void freeMemory( real**& o_stiffnessMatrices,
                     real*&  o_degreesOfFreedom,
                     real*&  o_timeDerivatives,
                     real*&  o_timeIntegrated,
                     real*&  o_timeExtrapolated ) {
      _mm_free( o_degreesOfFreedom  );

      for( unsigned int l_c = 0; l_c < 6; l_c++ ) {
        _mm_free( o_stiffnessMatrices[l_c] );
      }
      free( o_stiffnessMatrices );

      _mm_free( o_timeDerivatives      );
      _mm_free( o_timeIntegrated       );
      _mm_free( o_timeExtrapolated     );
    }

  public:
    /**
     * Sets up the time kernel.
     **/
    void setUp() {
      // generate random seed
      srand(time(NULL));
    }

    /**
     * Test the time kernel.
     **/
    void testTimeKernel() {
      // allocate memory
      real** l_stiffnessMatrices, *l_degreesOfFreedom, *l_timeDerivatives, *l_timeIntegrated, *l_timeExtrapolated;

      real l_starMatricesDense[3][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
      real l_starMatrices[3][STAR_NNZ];

#if STAR_NNZ == NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES
      real* l_starPointer = l_starMatricesDense[0];
#else
      real* l_starPointer = l_starMatrices[0];
#endif

      allocateMemory( l_stiffnessMatrices,
                      l_degreesOfFreedom,
                      l_timeDerivatives,
                      l_timeIntegrated,
                      l_timeExtrapolated );

      real l_degreesOfFreedomUT[NUMBER_OF_DOFS];
      real l_timeIntegratedUT[  NUMBER_OF_DOFS];
      real l_timeDerivativesUT[ CONVERGENCE_ORDER][NUMBER_OF_DOFS];
      real l_timeExtrapolatedUT[NUMBER_OF_DOFS];

      // location of the matrices
      std::string l_matricesPath = m_configuration.getMatricesFile();

      // read matrices
      m_denseMatrix.readMatrices( l_matricesPath );

      // intitialize stiffness matrices
      for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
        m_denseMatrix.initializeMatrix( 56+l_c,
                                        seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 ),
                                        seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER   ),
                                        l_stiffnessMatrices[l_c] );

        // dense matrix
        if( m_configuration.getSparseSwitch( 56+l_c ) == false ) {
          l_stiffnessMatrices[l_c+6] = l_stiffnessMatrices[l_c];
        }
        // sparse matrix
        else {
          m_denseMatrix.copyDenseToSparse( seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 ),
                                           seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER   ),
                                           l_stiffnessMatrices[l_c],
                                           l_stiffnessMatrices[l_c+3] );

          l_stiffnessMatrices[l_c+6] = l_stiffnessMatrices[l_c+3];
        }
      }

      // initialize a simple time intgrator
      m_simpleTimeIntegrator.initialize( l_matricesPath );

      // create a new time kernel
      seissol::kernels::Time l_timeKernel;

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        /*
         * Set random values for time derivatives
         */
        m_denseMatrix.setRandomValues( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                       l_degreesOfFreedom );

        m_denseMatrix.setRandomValues( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                       l_timeIntegrated );

        // copy values to UT-datastructures
        m_denseMatrix.copySubMatrix( l_degreesOfFreedom,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     l_degreesOfFreedomUT,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES );

        m_denseMatrix.copySubMatrix( l_timeIntegrated,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     l_timeIntegratedUT,
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES );

        for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
          // initialize star matrices
          m_denseMatrix.initializeMatrix( 59,
                                          NUMBER_OF_QUANTITIES,
                                          NUMBER_OF_QUANTITIES,
                                          l_starMatricesDense[l_c] );

          // init sparse if requested
          if( m_configuration.getSparseSwitch(59) == true ) {
            m_denseMatrix.copyDenseToSparse( NUMBER_OF_QUANTITIES,
                                             NUMBER_OF_QUANTITIES,
                                             l_starMatricesDense[l_c],
                                             l_starMatrices[l_c] );
          }
        }

        /*
         * Test hybrid ADER computation: time derivatives and time integrated DOFs in one step.
         */
        // time step width of this cell
        real l_cellDeltaT[2];
        l_cellDeltaT[0] = 0; m_denseMatrix.setRandomValues( 1, &l_cellDeltaT[1] );

        // reference implementation
        m_simpleTimeIntegrator.computeTimeDerivation( l_degreesOfFreedomUT,
                                                      l_starMatricesDense[0],
                                                      l_starMatricesDense[1],
                                                      l_starMatricesDense[2],
                                                      l_timeDerivativesUT );

        m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUT,
                                                       l_cellDeltaT,
                                                       l_timeIntegratedUT );

        // optimized implementation
        l_timeKernel.computeAder(                      l_cellDeltaT[1],
                                                       l_stiffnessMatrices+6,
                                                       l_degreesOfFreedom,
                                  (real(*) [STAR_NNZ]) l_starPointer,
                                                       l_timeIntegrated,
                                                       l_timeDerivatives );

        // check the results
        m_denseMatrix.checkResult( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                   l_degreesOfFreedom,
                                   l_timeDerivatives );

        m_denseMatrix.checkResult( NUMBER_OF_DOFS,
                                   l_degreesOfFreedomUT,
                                   l_timeDerivativesUT[0] );

        for( int l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
          unsigned int l_offset = l_timeKernel.m_derivativesOffsets[l_derivative];
#ifdef SINGLE_PRECISION
          if( !(CONVERGENCE_ORDER == 2 && ALIGNMENT >= 32) && !(CONVERGENCE_ORDER == 3 && ALIGNMENT >= 64) ) {
#endif
#ifdef DOUBLE_PRECISION
          if( !(CONVERGENCE_ORDER == 2 && ALIGNMENT >= 64) ) {
#endif
            m_denseMatrix.checkSubMatrix(  l_timeDerivativesUT[l_derivative],
                                           NUMBER_OF_BASIS_FUNCTIONS,
                                           NUMBER_OF_QUANTITIES,
                                          &l_timeDerivatives[ l_offset ],
                                           l_timeKernel.m_numberOfAlignedBasisFunctions[l_derivative],
                                           NUMBER_OF_QUANTITIES );
          }
          else {
            // 64 byte aligned data extendes the number of basis functions (4) of degree 2
            m_denseMatrix.checkSubMatrix( &l_timeDerivatives[ l_offset ],
                                           l_timeKernel.m_numberOfAlignedBasisFunctions[l_derivative],
                                           NUMBER_OF_QUANTITIES,
                                           l_timeDerivativesUT[l_derivative],
                                           NUMBER_OF_BASIS_FUNCTIONS,
                                           NUMBER_OF_QUANTITIES );
          }
        }

        m_denseMatrix.checkSubMatrix( l_timeIntegrated,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_timeIntegratedUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );


        /*
         * Test time integration from derivatives
         */
        real l_time[3]; real l_deltaT[2];
        m_denseMatrix.setRandomValues( 3, l_time );
        l_time[1] += l_time[0]; l_time[2] += l_time[1];
        l_deltaT[0] = l_time[1] - l_time[0]; l_deltaT[1] = l_time[2] - l_time[0];

        l_timeKernel.computeIntegral( l_time[0],
                                      l_time[1],
                                      l_time[2],
                                      l_timeDerivatives,
                                      l_timeIntegrated );

        m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUT,
                                                       l_deltaT,
                                                       l_timeIntegratedUT );

        // check the results
        m_denseMatrix.checkSubMatrix( l_timeIntegrated,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_timeIntegratedUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );

        /*
         * Test time extrapolation from derivatives
         */
        l_timeKernel.computeExtrapolation( l_time[0],
                                           l_time[1],
                                           l_timeDerivatives,
                                           l_timeExtrapolated );

        m_simpleTimeIntegrator.computeExtrapolation( l_timeDerivativesUT,
                                                     l_deltaT[0],
                                                     l_timeExtrapolatedUT );

        // check the results
        m_denseMatrix.checkSubMatrix( l_timeExtrapolated,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_timeExtrapolatedUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );

        /*
         * Test pure time integration.
         */
        l_cellDeltaT[0] = 0; m_denseMatrix.setRandomValues( 1, &l_cellDeltaT[1] );

        // reference implementation
        m_simpleTimeIntegrator.computeTimeDerivation( l_degreesOfFreedomUT,
                                                      l_starMatricesDense[0],
                                                      l_starMatricesDense[1],
                                                      l_starMatricesDense[2],
                                                      l_timeDerivativesUT );

        m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUT,
                                                       l_cellDeltaT,
                                                       l_timeIntegratedUT );

        // optimized implementation
        l_timeKernel.computeAder(                      l_cellDeltaT[1],
                                                       l_stiffnessMatrices+6,
                                                       l_degreesOfFreedom,
                                  (real(*) [STAR_NNZ]) l_starPointer,
                                                       l_timeIntegrated );

        // check the results
        m_denseMatrix.checkSubMatrix( l_timeIntegrated,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_timeIntegratedUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );
      }

      // free memory
      freeMemory( l_stiffnessMatrices,
                  l_degreesOfFreedom,
                  l_timeDerivatives,
                  l_timeIntegrated,
                  l_timeExtrapolated );
    }

    /*
     * Tests integration of neighboring cells (buffer copy or time integration).
     */
    void testNeighborsIntegration() {
      // time step width of the cell
      real l_timeStepWidth;
      m_denseMatrix.setRandomValues( 1, &l_timeStepWidth );

      // location of the matrices
      std::string l_matricesPath = m_configuration.getMatricesFile();

      // read matrices
      m_denseMatrix.readMatrices( l_matricesPath );

      // initialize a simple time intgrator
      m_simpleTimeIntegrator.initialize( l_matricesPath );

      // create a new time kernel
      seissol::kernels::Time l_timeKernel;

      //! face types
      enum faceType l_faceTypes[4] = {regular, regular, regular, regular};

      // pointers to time buffers/derivatives
      real *l_timeDofs[4];

      // final time integrated values of neighboring cells
      real *l_timeIntegrated[4];
      real  l_timeIntegratedUT[ 4 * NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES ];

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        /*
         * Setup
         */
        // generate a random LTS setup
        char l_ltsSetup = (rand() % 16);

        // iterate over neighbors and allocate memory for either an time integrated buffer or time derivatives
        for( unsigned int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          if( ( l_ltsSetup >> l_neighbor )%2 == 0 ) {
            // allocate memory for a time integrated buffer
            l_timeDofs[l_neighbor] = (real*) _mm_malloc( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );

            // fill with random values
            m_denseMatrix.setRandomValues( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                           l_timeDofs[l_neighbor] );
          }
          else {
            // allocate memory for time derivatives
            l_timeDofs[l_neighbor] = (real*) _mm_malloc( seissol::kernels::Time::getAlignedTimeDerivativesSize()*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );

            // fill with random values
            m_denseMatrix.setRandomValues( seissol::kernels::Time::getAlignedTimeDerivativesSize()*NUMBER_OF_QUANTITIES,
                                           l_timeDofs[l_neighbor] );

            // set entries with zero-padding to zero
            unsigned int l_firstEntry = 0;

            for( unsigned int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
              m_denseMatrix.setZeroBlock( seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order ),
                                          NUMBER_OF_QUANTITIES,
                                          seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER-l_order ),
                                         &l_timeDofs[l_neighbor][l_firstEntry] );

              l_firstEntry += seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order ) * NUMBER_OF_QUANTITIES;
            }
          }
        }

        // current timings of all cells
        double l_currentTime[5];
        m_denseMatrix.setRandomValues( 5, l_currentTime );

        for( int l_neighbor = 1; l_neighbor < 5; l_neighbor++ ) {
          l_currentTime[0] = std::max( l_currentTime[0], l_currentTime[l_neighbor] );
        }

        /*
         * time integration
         */
        // reference implementation
        for( int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          if( ( l_ltsSetup >> l_neighbor )%2 == 0 ) {
            m_denseMatrix.copySubMatrix( l_timeDofs[l_neighbor],
                                         NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                         NUMBER_OF_QUANTITIES,
                                        &l_timeIntegratedUT[l_neighbor*NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES],
                                         NUMBER_OF_BASIS_FUNCTIONS,
                                         NUMBER_OF_QUANTITIES );
          }
          else {
            // copy time derivatives to non-compressed storage scheme
            real l_timeDerivativesUT[CONVERGENCE_ORDER][NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES];

            seissol::kernels::convertAlignedCompressedTimeDerivatives( l_timeDofs[l_neighbor],
                                                                       l_timeDerivativesUT );

            // do the time integration
            real l_deltaT[2];
            l_deltaT[0] = l_currentTime[0]                   - l_currentTime[l_neighbor+1];
            l_deltaT[1] = l_currentTime[0] + l_timeStepWidth - l_currentTime[l_neighbor+1];

            m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUT,
                                                           l_deltaT,
                                                          &l_timeIntegratedUT[l_neighbor*NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES] );
          }
        }

        // optimized implementation
        real l_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
        l_timeKernel.computeIntegrals( l_ltsSetup,
                                       l_faceTypes,
                                       l_currentTime,
                                       l_timeStepWidth,
                                       l_timeDofs,
                                       l_integrationBuffer,
                                       l_timeIntegrated );

        /*
         * check the results
         */
        for( int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          m_denseMatrix.checkSubMatrix( l_timeIntegrated[l_neighbor],
                                        NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                        NUMBER_OF_QUANTITIES,
                                        &l_timeIntegratedUT[l_neighbor*NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES],
                                        NUMBER_OF_BASIS_FUNCTIONS,
                                        NUMBER_OF_QUANTITIES );
        }

        /*
         * Free memory
         */
        for( int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          _mm_free( l_timeDofs[l_neighbor] );
        }
      }
    }

    void testAderFlops() {
      unsigned int l_nonZeroFlops  = 0;
      unsigned int l_hardwareFlops = 0;

      // volume kernel
      seissol::kernels::Time l_timeKernel;

      // get the flops for a single call
      l_timeKernel.flopsAder( l_nonZeroFlops, l_hardwareFlops );

      // assert we are at least doing the sparse flops in hardware
      TS_ASSERT_LESS_THAN_EQUALS( l_nonZeroFlops, l_hardwareFlops );


      // derive dense flops
      unsigned int l_denseFlops = NUMBER_OF_ALIGNED_DOFS;

      for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
        // stiffness matrices
        l_denseFlops += ( seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER - l_derivative     ) *
                          seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER - l_derivative + 1 ) *
                          NUMBER_OF_QUANTITIES * 2
                        ) * 3;
        // star matrices
        l_denseFlops += ( seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER - l_derivative     ) *
                          STAR_NNZ * 2
                        ) * 3;

        // update of the time integrated dofs
        l_denseFlops += seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER - l_derivative ) * NUMBER_OF_QUANTITIES * 2;
      }

      // assert we are not doing more FLOPS in hardware than for dense-only stifness execution
      TS_ASSERT_LESS_THAN_EQUALS( l_hardwareFlops, l_denseFlops );

    }

    /*
     * Tests the derivation of the LTS setup.
     */
    void testLtsSetup() {
      unsigned int  l_localClusterId;
      unsigned int  l_neighboringClusterIds[4];
      unsigned int  l_faceNeighborIds[4];
      enum faceType l_faceTypes[4];
      unsigned short l_ltsSetup;

      /*
       * Global time stepping
       * | - - 0 1 | 1 1 1 1 | 0 0 0 0 |
       */
      l_localClusterId = 4;
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        l_neighboringClusterIds[l_face] = l_localClusterId;
        l_faceNeighborIds[l_face]       = rand()%10000;
        l_faceTypes[l_face]             = regular;
      }
      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 496 );

      /*
       * Global time stepping with dynamic rupture.
       * | - 0 1 1 | 1 1 1 1 | 1 0 0 0 |
       */
      l_faceTypes[3] = dynamicRupture;
      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 1016 );

      /*
       * Local time stepping.
       * | - 1 1 1 | 1 0 1 0 | 0 0 0 1 |
       */
      l_faceTypes[3] = periodic;
      l_neighboringClusterIds[0] = 5;
      l_neighboringClusterIds[2] = 2;
      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 1953 );

      /*
       * Local time stepping with dynamic rupture and boundary conditions.
       * | - 1 1 1 | 1 0 0 0 | 1 0 1 0 |
       */
      l_faceTypes[0]           = outflow;
      l_faceTypes[3]           = dynamicRupture;
      l_neighboringClusterIds[1] = 10;
      l_neighboringClusterIds[2] = 0;
      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 1930 );

      /*
       * Local time stepping with "free surface on buffers".
       * | - 0 1 1 | 0 1 1 1 | 0 0 0 0 |
       */
      l_faceTypes[0] = l_faceTypes[2] = l_faceTypes[3] = regular;
      l_faceTypes[1] = freeSurface;
      l_neighboringClusterIds[0] = l_neighboringClusterIds[1] = l_neighboringClusterIds[2] =  4;
      l_neighboringClusterIds[3] = 3;

      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 880 );

      /*
       * Local time stepping with "free surface on derivatives".
       * | - 0 1 0 | 0 0 1 0 | 0 0 1 0 |
       */

      l_neighboringClusterIds[0] = l_neighboringClusterIds[2] =  3;
      l_ltsSetup = seissol::kernels::Time::getLtsSetup( l_localClusterId, l_neighboringClusterIds, l_faceTypes, l_faceNeighborIds );
      TS_ASSERT_EQUALS( l_ltsSetup, 546 );
    }

    /*
     * Test LTS normalization.
     */
    void testLtsNormalization() {
      unsigned short int l_localSetup;
      unsigned short int l_neighboringSetups[4];

      /*
       * GTS: No changes are required
       */
      l_localSetup = l_neighboringSetups[0] = l_neighboringSetups[1] =
                     l_neighboringSetups[2] = l_neighboringSetups[3] = 240;

      seissol::kernels::Time::normalizeLtsSetup( l_neighboringSetups, l_localSetup );
      TS_ASSERT_EQUALS( l_localSetup, 240 );

      /*
       * LTS with uncritical time relations in the neighbors: no changes are required.
       *
       * l: | - - 0 1 | 1 1 1 1 | 0 0 0 0 | 
       *
       * 1: | - - 1 1 | 0 1 1 0 | 1 0 0 1 |
       * 2: | - - 1 1 | 1 0 0 0 | 0 1 1 1 |
       */
      l_neighboringSetups[1] = 873;
      l_neighboringSetups[2] = 903;

      seissol::kernels::Time::normalizeLtsSetup( l_neighboringSetups, l_localSetup );
      TS_ASSERT_EQUALS( l_localSetup, 240 );

      /*
       * LTS with critical time relation in the neighbor.
       *
       * l: | - - | 0 1 | 1 1 1 1 | 0 0 0 0 |
       *
       * 1: | - 1 | 1 1 | 1 1 1 0 | 0 0 0 0 |
       */
      l_localSetup = l_neighboringSetups[0] =
                     l_neighboringSetups[2] = l_neighboringSetups[3] = 240;

      l_neighboringSetups[1] = 2016;

      seissol::kernels::Time::normalizeLtsSetup( l_neighboringSetups, l_localSetup );
      TS_ASSERT_EQUALS( l_localSetup, 242 );
    }
};
