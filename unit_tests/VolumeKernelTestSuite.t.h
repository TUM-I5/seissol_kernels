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
 * Test suite, which tests the volume kernel.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"

#include <Initializer/MemoryManager.h>
#include "../src/Time.h"
#include "../src/Volume.h"
#include "DenseMatrix.hpp"
#include "SimpleVolumeIntegrator.hpp"

namespace unit_test {
  class VolumeKernelTestSuite;
}

class unit_test::VolumeKernelTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

    /**
     * Allocated aligned memory for the given data structures.
     *
     * @param o_stiffnessMatrices stiffness matrices.
     * @param o_timeIntegrated time integrated DOFs.
     * @param o_degreesOfFreedom degrees of freedom.
     **/
    void allocateMemory( real**& o_stiffnessMatrices,
                         real*&  o_timeIntegrated,
                         real*&  o_degreesOfFreedom ) {

      unsigned int l_alignedNumberOfBasisFunctions = seissol::kernels::getNumberOfAlignedBasisFunctions();

      // allocate memory for the stiffness matrices
      unsigned l_sizesStiff[2];
      l_sizesStiff[0] = l_alignedNumberOfBasisFunctions;
      l_sizesStiff[1] = seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER-1 );

      // allocate memory
      o_stiffnessMatrices    = (real**) malloc( 9*sizeof(real*) );
      for( unsigned int l_c = 0; l_c < 6; l_c++ ) {
        o_stiffnessMatrices[l_c] = (real*) _mm_malloc( l_sizesStiff[0]*l_sizesStiff[1]*sizeof(real),                  ALIGNMENT );
      }

      o_degreesOfFreedom     = (real*) _mm_malloc( l_alignedNumberOfBasisFunctions*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );

      o_timeIntegrated       = (real*) _mm_malloc( l_alignedNumberOfBasisFunctions*NUMBER_OF_QUANTITIES*sizeof(real), ALIGNMENT );
    }

    /**
     * Frees memory of the given data structures.
     *
     * @param o_stiffnessMatrices stiffness matrices.
     * @param o_timeIntegrated time integrated DOFs.
     * @param o_degreesOfFreedom degrees of freedom.
     **/
    void freeMemory( real**& o_stiffnessMatrices,
                     real*&  o_timeIntegrated,
                     real*&  o_degreesOfFreedom ) {

      for( unsigned int l_c = 0; l_c < 6; l_c++ ) {
        _mm_free( o_stiffnessMatrices[l_c] );
      }
      free( o_stiffnessMatrices );

      _mm_free( o_timeIntegrated );

      _mm_free( o_degreesOfFreedom  );
    }

  public:
    void setUp() {
      // generate random seed
      srand(time(NULL));
    }

    /**
     * Test the volume kernel.
     **/
    void testVolumeKernel() {
      // allocate memory
      real **l_stiffnessMatrices, *l_timeIntegrated, *l_degreesOfFreedom;

      real l_starMatricesDense[3][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
      real l_starMatrices[3][STAR_NNZ];

#if STAR_NNZ == NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES
      real* l_starPointer = l_starMatricesDense[0];
#else
      real* l_starPointer = l_starMatrices[0];
#endif

      allocateMemory( l_stiffnessMatrices,
                      l_timeIntegrated,
                      l_degreesOfFreedom );

      real l_timeIntegratedUT[  NUMBER_OF_DOFS];
      real l_degreesOfFreedomUT[NUMBER_OF_DOFS];

      // setup matrix path
      std::string l_matricesPath = m_configuration.getMatricesFile();

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // set up volume integrators
      seissol::kernels::Volume l_volumeKernel;
      unit_test::SimpleVolumeIntegrator l_simpleVolumeIntegrator;
      l_simpleVolumeIntegrator.initialize( l_matricesPath );

      // intitialize stiffness matrices
      for( unsigned int l_c = 0; l_c < 3; l_c++ ) {
        // init stiffness matrix
        m_denseMatrix.initializeMatrix( 53+l_c,
                                        seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER   ),
                                        seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER-1 ),
                                        l_stiffnessMatrices[l_c] );

        if( m_configuration.getSparseSwitch( 53+l_c ) == false ) {
          l_stiffnessMatrices[l_c+6] = l_stiffnessMatrices[l_c];
        }
        // sparse matrix
        else {
          m_denseMatrix.copyDenseToSparse( seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER   ),
                                           seissol::kernels::getNumberOfBasisFunctions(        CONVERGENCE_ORDER-1 ),
                                           l_stiffnessMatrices[l_c],
                                           l_stiffnessMatrices[l_c+3] );

          l_stiffnessMatrices[l_c+6] = l_stiffnessMatrices[l_c+3];
        }
      }

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        // set random values to DOFs and time integrated DOFs
        m_denseMatrix.setRandomValues( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                       l_degreesOfFreedom );

        m_denseMatrix.setRandomValues( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES,
                                       l_timeIntegrated );

        // set padded values to zero
        m_denseMatrix.setZeroBlock( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                    NUMBER_OF_QUANTITIES,
                                    NUMBER_OF_BASIS_FUNCTIONS,
                                    l_degreesOfFreedom );

        m_denseMatrix.setZeroBlock( NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                    NUMBER_OF_QUANTITIES,
                                    NUMBER_OF_BASIS_FUNCTIONS,
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
         * Test volume integration
         */
        // standard multiplication
        l_simpleVolumeIntegrator.computeVolumeIntegration( l_timeIntegratedUT,
                                                           l_starMatricesDense,
                                                           l_degreesOfFreedomUT );
        // genrated multiplication
        l_volumeKernel.computeIntegral(                      l_stiffnessMatrices+6,
                                                             l_timeIntegrated,
                                        (real(*) [STAR_NNZ]) l_starPointer,
                                                             l_degreesOfFreedom );

        // check the results
        m_denseMatrix.checkSubMatrix( l_degreesOfFreedom,
                                      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES,
                                      l_degreesOfFreedomUT,
                                      NUMBER_OF_BASIS_FUNCTIONS,
                                      NUMBER_OF_QUANTITIES );
      }

      freeMemory( l_stiffnessMatrices,
                  l_timeIntegrated,
                  l_degreesOfFreedom );
    }

    /**
     * Tests the flop-counters.
     **/
    void testIntegralFlops() {
      unsigned int l_nonZeroFlops  = 0;
      unsigned int l_hardwareFlops = 0;

      // volume kernel
      seissol::kernels::Volume l_volumeKernel;

      // get the flops for a single call
      l_volumeKernel.flopsIntegral( l_nonZeroFlops, l_hardwareFlops );

      // assert we are at least doing the sparse flops in hardware
      TS_ASSERT_LESS_THAN_EQUALS( l_nonZeroFlops, l_hardwareFlops );


      // stiffness multiplications
      unsigned int l_denseFlops =  (  NUMBER_OF_ALIGNED_BASIS_FUNCTIONS
                                    * seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER - 1 )
                                    * NUMBER_OF_QUANTITIES * 2) * 3;

      // star matrices
      l_denseFlops += ( STAR_NNZ * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * 2 ) * 3;

      // assert we are not doing more FLOPS in hardware than for dense-only stifness execution
      TS_ASSERT_LESS_THAN_EQUALS( l_hardwareFlops, l_denseFlops );
    }
};
