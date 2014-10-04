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
      l_sizesStiff[1] = l_alignedNumberOfBasisFunctions;

      // allocate memory
      o_stiffnessMatrices    = (real**)    malloc( 3*sizeof(real*) );
      o_stiffnessMatrices[0] = (real*) _mm_malloc( l_sizesStiff[0]*l_sizesStiff[1]*sizeof(real),                      ALIGNMENT );
      o_stiffnessMatrices[1] = (real*) _mm_malloc( l_sizesStiff[0]*l_sizesStiff[1]*sizeof(real),                      ALIGNMENT );
      o_stiffnessMatrices[2] = (real*) _mm_malloc( l_sizesStiff[0]*l_sizesStiff[1]*sizeof(real),                      ALIGNMENT );

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

      _mm_free( o_stiffnessMatrices[0] );
      _mm_free( o_stiffnessMatrices[1] );
      _mm_free( o_stiffnessMatrices[2] );
          free( o_stiffnessMatrices    );

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

      real l_starMatrices[3][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];

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
      // TODO: for now it's all handled dense
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBER_OF_BASIS_FUNCTIONS) + ".xml";

      // read matrices
      m_denseMatrix.readMatrices( l_matricesPath );

      // intitialize stiffness matrices
      for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
        m_denseMatrix.initializeMatrix( 56+l_c,
                                        seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 ),
                                        seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER   ),
                                        l_stiffnessMatrices[l_c] );
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
                                          l_starMatrices[l_c] );
        }

        /*
         * Test time derivatives
         */
        // reference implementation
        m_simpleTimeIntegrator.computeTimeDerivation( l_degreesOfFreedomUT,
                                                      l_starMatrices[0],
                                                      l_starMatrices[1],
                                                      l_starMatrices[2],
                                                      l_timeDerivativesUT );

        // optimized implementation
        l_timeKernel.computeDerivatives( l_stiffnessMatrices,
                                         l_degreesOfFreedom,
                                         l_starMatrices,
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
          m_denseMatrix.checkSubMatrix(  l_timeDerivativesUT[l_derivative],
                                         NUMBER_OF_BASIS_FUNCTIONS,
                                         NUMBER_OF_QUANTITIES,
                                        &l_timeDerivatives[ l_offset ],
                                         l_timeKernel.m_numberOfAlignedBasisFunctions[l_derivative],
                                         NUMBER_OF_QUANTITIES );
        }

        /*
         * Test time integration
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
         * Test time extrapolation
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
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBER_OF_BASIS_FUNCTIONS) + ".xml";

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
      real l_timeIntegrated[4][NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES] __attribute__((aligned(ALIGNMENT)));
      real l_timeIntegratedUT[ 4 * NUMBER_OF_BASIS_FUNCTIONS * NUMBER_OF_QUANTITIES ];

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        /*
         * Setup
         */
        // generate a random LTS setup
        unsigned int l_ltsSetup = ((rand() % 16)<<8);

        // iterate over neighbors and allocate memory for either an time integrated buffer or time derivatives
        for( unsigned int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          if( ( l_ltsSetup >> (8 + l_neighbor) )%2 == 0 ) {
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
        real l_currentTime[5];
        m_denseMatrix.setRandomValues( 5, l_currentTime );

        for( int l_neighbor = 1; l_neighbor < 5; l_neighbor++ ) {
          l_currentTime[0] = std::max( l_currentTime[0], l_currentTime[l_neighbor] );
        }

        /*
         * time integration
         */
        // reference implementation
        for( int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
          if( ( l_ltsSetup >> (8 + l_neighbor) )%2 == 0 ) {
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
        l_timeKernel.computeIntegrals( l_ltsSetup,
                                       l_faceTypes,
                                       l_currentTime,
                                       l_timeStepWidth,
                                       l_timeDofs,
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
};
