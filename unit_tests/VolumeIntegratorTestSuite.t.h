/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Test suite, which tests the volume integrator.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"

#include <Initializer/MemoryManager.h>
#include "../src/VolumeIntegrator.h"
#include "DenseMatrix.hpp"
#include "SimpleVolumeIntegrator.hpp"

namespace unit_test {
  class VolumeIntegratorTestSuite;
}

class unit_test::VolumeIntegratorTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

  public:
    void setUp() {
      // generate random seed
      srand(time(NULL));
    }

    /**
     * Tests the complete kernel, which does the volume integration.
     **/
    void testVolumeIntegrator() {
      //! setup matrix path
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );
      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );
      // set up volume integrators
      seissol::kernels::VolumeIntegrator l_volumeIntegrator( l_matrixReader, l_memoryManager );
      unit_test::SimpleVolumeIntegrator l_simpleVolumeIntegrator;
      l_simpleVolumeIntegrator.initialize( l_matricesPath );

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // unknowns: \f$ Q_k \f$ for the volume integration in the unit tests and the kernel
      //            updated by the volume integrated unknowns:
      //            \f$ Q_k += K^\xi I(t^{n}, t^{n+1}, Q_{k}^n) A^* + K^\eta I(t^{n}, t^{n+1}, Q_{k}^n) B^* + K^\zeta I(t^{n}, t^{n+1}, Q_{k}^n) C^* \f$
      double l_unknownsUnitTests[NUMBEROFUNKNOWNS];
      double l_unknownsVolumeIntegrator[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

      // time integrated unknowns I(t^{n}, t^{n+1}, Q_{k}^n)
      double l_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

      // star matrices: \f$ A^*, B^*, C^* \f$
      double l_starMatricesDense[3][NUMBEROFVARIABLES*NUMBEROFVARIABLES];
      double l_starMatricesSparse[3][STARMATRIX_NUMBEROFNONZEROS] __attribute__((aligned(64)));;

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        /*
         * Initialize the matrices.
         */
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_timeIntegratedUnknowns );
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_unknownsUnitTests );

        // synchronize unknowns matrices 
        for( unsigned int l_element = 0; l_element < NUMBEROFUNKNOWNS; l_element++ ) {
          l_unknownsVolumeIntegrator[l_element]  = l_unknownsUnitTests[l_element];
        }

        // iterate over the three coordinates \f$ \xi, \eta, \zeta \f$
        for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) {
          // initialize dense star matrices
          m_denseMatrix.initializeMatrix( 59,
                                          NUMBEROFVARIABLES,
                                          NUMBEROFVARIABLES,
                                          l_starMatricesDense[l_coordinate] );

          // copy non-zeros to sparse star matrix
          m_denseMatrix.copyDenseToSparse( NUMBEROFVARIABLES,
                                           NUMBEROFVARIABLES,
                                           l_starMatricesDense[l_coordinate],
                                           l_starMatricesSparse[l_coordinate] );
        }

        /*
         * Compute the volume integrations
         */
        // standard multiplication
        l_simpleVolumeIntegrator.computeVolumeIntegration( l_timeIntegratedUnknowns,
                                                           l_starMatricesDense,
                                                           l_unknownsUnitTests );

        // genrated multiplication
        l_volumeIntegrator.computeVolumeIntegral( l_timeIntegratedUnknowns,
                                                  l_starMatricesSparse[0],
                                                  l_starMatricesSparse[1],
                                                  l_starMatricesSparse[2],
                                                  l_unknownsVolumeIntegrator );

        /*
         * check the results
         */
        m_denseMatrix.checkResult( NUMBEROFUNKNOWNS,l_unknownsUnitTests, l_unknownsVolumeIntegrator );
      }
    }
};
