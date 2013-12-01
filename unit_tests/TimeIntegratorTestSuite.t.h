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
 * Test suite, which tests the time integrator.
 **/

#include <iostream>
#include <cxxtest/TestSuite.h>
#include "configuration.hpp"

#include <Initializer/MemoryManager.h>
#include "../src/TimeIntegrator.h"
#include "DenseMatrix.hpp"

namespace unit_test {
  class TimeIntegratorTestSuite;
}

class unit_test::TimeIntegratorTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

    /**
     * Simple time derivation to test again.
     *
     * @param i_unknowns unknowns for which the time derivatives are computed.
     * @param i_transposedStiffnessMatrices transposed stiffness matrices (multiplied by the inverse mass matrix):  \f$ M^{-1} (K^\xi)^T, M^{-1} (K^\eta)^T, M^{-1} (K^\zeta)^T \f$
     * @param i_starMatrices  star matrices:
     *        0: \f$ A^*_k \f$
     *        1: \f$ B^*_k \f$
     *        2: \f$ C^*_k \f$
     * @param o_timeDerivatives time derivates.
     **/
    void simpleTimeDerivation( const double i_unknowns[NUMBEROFUNKNOWNS],
                               const double i_transposedStiffnessMatrices[3][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS],
                               const double i_starMatrices[3][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                     double o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] ) {
      // temporary matrix for two-step multiplication
      double l_temporaryProduct[NUMBEROFUNKNOWNS];
      
      // copy unknowns to zeroth derivative
      for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
        o_timeDerivatives[0][l_entry] = i_unknowns[l_entry];
      }

      // iterate over order in time
      for( int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) {
        // set time derivatives of this order to zero
        std::fill( o_timeDerivatives[l_order], o_timeDerivatives[l_order]+NUMBEROFUNKNOWNS, 0);

        // iterate over the three reference coordinates: \f$ \xi, \eta, \zeta \f$
        for( int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) {
          // set temporary product to zero
          std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0);

          // stiffness matrix multiplication
          m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                       i_transposedStiffnessMatrices[ l_coordinate ], o_timeDerivatives[l_order-1], l_temporaryProduct ); 

          // star matrix multiplcation
          m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                       l_temporaryProduct, i_starMatrices[ l_coordinate ], o_timeDerivatives[ l_order ] );
        }
      }
    }

    /**
     * Simple time integration to test again.
     *
     * @param i_timeDerivatives time derivatives of the unknowns in the cell.
     * @param i_deltaT integration boundaries \f$ \Delta t^\text{lo} \f$ and \f$ \Delta t^\text{up} \f$ relative to the time level \f$ \t^\text{cell} \f$ of the cell; integration is performned over the interval \f$ [ t^\text{cell} + \Delta t^\text{lo}, t^\text{cell} + \Delta t^\text{up} ] \f$
     * @param o_timeIntegratedUnknowns unknowns integrated over the specified interval.
     **/
    void simpleTimeIntegration( const double i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                const double i_deltaT[2],
                                      double o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
      // assert integration forward in time
      TS_ASSERT( i_deltaT[1] > i_deltaT[0] );

      // initialization of Taylor series factors
      double l_taylorSeriesFactors[2];
      l_taylorSeriesFactors[0] = i_deltaT[0];
      l_taylorSeriesFactors[1] = i_deltaT[1];
      double l_taylorSeriesFactor = l_taylorSeriesFactors[1] - l_taylorSeriesFactors[0];

      // update time integrated unknowns with zeroth derivatives
      for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
        o_timeIntegratedUnknowns[l_entry] = l_taylorSeriesFactor * i_timeDerivatives[0][l_entry];
      }

      // iterate over order in time
      for( int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) {
        // compute factor of the taylor series
        l_taylorSeriesFactors[0] = -l_taylorSeriesFactors[0] * i_deltaT[0] / double(l_order+1);
        l_taylorSeriesFactors[1] = -l_taylorSeriesFactors[1] * i_deltaT[1] / double(l_order+1);
        l_taylorSeriesFactor = l_taylorSeriesFactors[1] - l_taylorSeriesFactors[0];
       
        // update time integrated unknowns
        for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
          o_timeIntegratedUnknowns[l_entry] += l_taylorSeriesFactor * i_timeDerivatives[l_order][l_entry];
        } 
      }
    }

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
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );
      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );
      // set up the time integrator
      seissol::kernels::TimeIntegrator l_timeIntegrator( l_matrixReader, l_memoryManager );

      // unknowns of the cell
      double l_unknowns[NUMBEROFUNKNOWNS] __attribute__((aligned(64)));
 
      // time derivatives of the unknowns
      double l_timeDerivativesUnitTests[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));
      double l_timeDerivativesTimeIntegrator[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));

      // time integrated unknowns
      double l_timeIntegratedUnknownsUnitTests[NUMBEROFUNKNOWNS] __attribute((aligned(64)));
      double l_timeIntegratedUnknownsTimeIntegrator[NUMBEROFUNKNOWNS] __attribute((aligned(64)));

      // stiffness matrices (multiplied by the inverse of the mass matrix): \f$ M^{-1} K^\xi, M^{-1} K^\eta, M^{-1} K^\zeta \f$
      double l_stiffnessMatricesDense[3][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS];

      // star matrices: \f$ A^*, B^*, C^* \f$
      double l_starMatricesDense[3][NUMBEROFVARIABLES*NUMBEROFVARIABLES];
      double l_starMatricesSparse[3][STARMATRIX_NUMBEROFNONZEROS] __attribute__((aligned(64)));;

      // ranges for the time integration
      // GTS:      0      -- l_deltaT[1]
      // LTS: l_deltaT[0] -- l_deltaT[1]
      double l_deltaT[2];

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
          // initialize unit test transposed stiffness matrix (multiplied by the inverse of the mass matrix)
          m_denseMatrix.initializeMatrix( 56+l_coordinate,
                                          NUMBEROFBASISFUNCTIONS,
                                          NUMBEROFBASISFUNCTIONS,
                                          l_stiffnessMatricesDense[l_coordinate] );

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
        l_deltaT[1] = std::abs( l_deltaT[1] )+ 10;

        /*
         * LTS: Compute the time derivatives.
         */
        // standard multiplication
        simpleTimeDerivation( l_unknowns,
                              l_stiffnessMatricesDense,
                              l_starMatricesDense,
                              l_timeDerivativesUnitTests );
        // genrated multiplication
        l_timeIntegrator.computeTimeDerivatives( l_unknowns,
                                                 l_starMatricesSparse[0],
                                                 l_starMatricesSparse[1],
                                                 l_starMatricesSparse[2],
                                                 l_timeDerivativesTimeIntegrator );

        /*
         * LTS: Compute the time integrals.
         */
        // standard evaluation
        simpleTimeIntegration( l_timeDerivativesUnitTests,
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
        simpleTimeIntegration( l_timeDerivativesUnitTests,
                               l_deltaT,
                               l_timeIntegratedUnknownsUnitTests );

        // generated multiplication
        l_timeIntegrator.computeTimeIntegral( l_unknowns,
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
