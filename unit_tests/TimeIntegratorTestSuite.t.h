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
 * Test suite, which tests the time integrator.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"

#include <Initializer/MemoryManager.h>
#include "../src/TimeIntegrator.h"
#include "DenseMatrix.hpp"
#include "SimpleTimeIntegrator.hpp"

namespace unit_test {
  class TimeIntegratorTestSuite;
}

class unit_test::TimeIntegratorTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! Simple time integration functionality
    unit_test::SimpleTimeIntegrator m_simpleTimeIntegrator;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

  public:
    void setUp() {
      // generate random seed
      srand(time(NULL));
    }

    /**
     * Tests the complete kernel, which does the time integration.
     * Both version are tested: Global (GTS) & local time stepping (LTS).
     **/
    void testTimeIntegrator() {
      //! setup matrix path
#ifdef __MIC__
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".mic.xml";
#else
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";
#endif

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );
      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );
      // set up the time integrator
      seissol::kernels::TimeIntegrator l_timeIntegrator( l_matrixReader, l_memoryManager );

      // stiffness matrices
      double* l_stiffnessMatrices[3];
      for( int l_matrix = 0; l_matrix < 3; l_matrix++ ) {
        l_stiffnessMatrices[l_matrix] = (l_memoryManager.getStiffnessMatrixPointers())[l_matrix+3]; 
      }

      // unknowns of the cell
      double l_unknowns[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));
 
      // time derivatives of the unknowns
      double l_timeDerivativesUnitTests[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));
      double l_timeDerivativesTimeIntegrator[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

      // time integrated unknowns
      double l_timeIntegratedUnknownsUnitTests[NUMBEROFUNKNOWNS] __attribute((aligned(64)));
      double l_timeIntegratedUnknownsTimeIntegrator[NUMBEROFUNKNOWNS] __attribute((aligned(64)));

      // star matrices: \f$ A^*, B^*, C^* \f$
      double l_starMatricesDense[3][NUMBEROFVARIABLES*NUMBEROFVARIABLES];
      double l_starMatricesSparse[3][STARMATRIX_NUMBEROFNONZEROS] __attribute__((aligned(64)));;

      // ranges for the time integration
      // GTS:      0      -- l_deltaT[1]
      // LTS: l_deltaT[0] -- l_deltaT[1]
      double l_deltaT[2];

      // initialize stiffness matrices
      m_simpleTimeIntegrator.initialize( l_matricesPath );

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // repeat the test
      for( int l_repeat = 0; l_repeat < m_configuration.getNumberOfRepeats(); l_repeat++) {
        /*
         * Initialize the matrices.
         */
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_unknowns );

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

        // set random time step widths
        m_denseMatrix.setRandomValues( 2,
                                       l_deltaT );
        l_deltaT[0] = std::abs( l_deltaT[0] );
        l_deltaT[1] = std::abs( l_deltaT[1] )+ 0.5;

        /*
         * LTS: Compute the time derivatives.
         */
        // standard multiplication
        m_simpleTimeIntegrator.computeTimeDerivation( l_unknowns,
                                                      l_starMatricesDense[0],
                                                      l_starMatricesDense[1],
                                                      l_starMatricesDense[2],
                                                      l_timeDerivativesUnitTests );
        // genrated multiplication
        l_timeIntegrator.computeTimeDerivatives( l_unknowns,
                                                 l_stiffnessMatrices,
                                                 l_starMatricesSparse[0],
                                                 l_starMatricesSparse[1],
                                                 l_starMatricesSparse[2],
                                                 l_timeDerivativesTimeIntegrator );

        /*
         * LTS: Compute the time integrals.
         */
        // standard evaluation
        m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUnitTests,
                                                       l_deltaT,
                                                       l_timeIntegratedUnknownsUnitTests );

        // "generated" evaluation
        l_timeIntegrator.computeTimeIntegral( l_timeDerivativesTimeIntegrator,
                                              l_deltaT[0],
                                              l_deltaT[1],
                                              l_timeIntegratedUnknownsTimeIntegrator );

        /*
         * LTS: Check the results.
         */
        // time derivatives
        for( int l_order = 0; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) {
          m_denseMatrix.checkResult( NUMBEROFUNKNOWNS,
                                     l_timeDerivativesUnitTests[l_order],
                                     l_timeDerivativesTimeIntegrator[l_order] );
        }

        // time integrals
        m_denseMatrix.checkResult( NUMBEROFUNKNOWNS,
                                   l_timeIntegratedUnknownsUnitTests,
                                   l_timeIntegratedUnknownsTimeIntegrator );

        /*
         * GTS: Compute time integrals.
         */
        // standard evaluation (uses time derivatives of the LTS functionality)
        l_deltaT[0] = 0;
        m_simpleTimeIntegrator.computeTimeIntegration( l_timeDerivativesUnitTests,
                                                       l_deltaT,
                                                       l_timeIntegratedUnknownsUnitTests );

        // generated multiplication
        l_timeIntegrator.computeTimeIntegral( l_unknowns,
                                              l_stiffnessMatrices,
                                              l_starMatricesSparse[0],
                                              l_starMatricesSparse[1],
                                              l_starMatricesSparse[2],
                                              l_deltaT[1],
                                              l_timeIntegratedUnknownsTimeIntegrator );

        /*
         * GTS: check the results.
         */
        m_denseMatrix.checkResult( NUMBEROFUNKNOWNS, l_timeIntegratedUnknownsUnitTests, l_timeIntegratedUnknownsTimeIntegrator );
      }
    }
};
