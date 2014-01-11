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
#include "../src/BoundaryIntegrator.h"
#undef private
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
    void testBoundaryIntegrator() {
      //! path to matrix set-up
      std::string l_matricesPath = m_matricesDirectory + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      //! matrix of unknowns
      double l_unknownsUnitTests[NUMBEROFUNKNOWNS] = {0};
      double l_unknownsBoundaryIntegrator[NUMBEROFUNKNOWNS]         __attribute__((aligned(64)));

      //! matrices of time integrated unknowns
      double l_timeIntegratedUnknownsElement[2][NUMBEROFUNKNOWNS]   __attribute__((aligned(64)));
      double l_timeIntegratedUnknownsNeighbors[4][NUMBEROFUNKNOWNS] __attribute__((aligned(64)));
      double* l_timeIntegratedUnknownsNeighborPointers[4];
      for( int l_face = 0; l_face < 4; l_face++ ) {
        l_timeIntegratedUnknownsNeighborPointers[l_face] = &l_timeIntegratedUnknownsNeighbors[l_face][0]; 
      }

      //! flux solvers matrices (positive eigenvalues): \f$ N_{k,i} A_k^+ N_{k,i}^{-1} \f$
      double l_fluxSolversPos[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! flux solvers matrices (negative eigenvalues): \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
      double l_fluxSolversNeg[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES];

      //! boundary conditions
      int l_boundaryConditions[4] = {0, 0, 0, 0};

      //! neighboring indices
      int l_neighboringIndices[4][2];

      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );

      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );

      // construct boundary integrators
      unit_test::SimpleBoundaryIntegrator l_simpleBoundaryIntegrator;
      l_simpleBoundaryIntegrator.initialize( l_matricesPath );
      seissol::kernels::BoundaryIntegrator l_boundaryIntegrator( l_matrixReader, l_memoryManager );

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
                                       l_timeIntegratedUnknownsElement[0] );
        m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS,
                                       l_timeIntegratedUnknownsElement[1] );


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

          l_neighboringIndices[l_face][0] = rand() % 4;
          l_neighboringIndices[l_face][1] = rand() % 3;
        }

        // do a simple boundary intagration
        l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUnknownsElement[0],
                                                               l_timeIntegratedUnknownsNeighbors,
                                                               l_boundaryConditions,
                                                               l_neighboringIndices,
                                                               l_fluxSolversPos,
                                                               l_fluxSolversNeg,
                                                               l_unknownsUnitTests );

        // do boundary integration with generated kernels
        l_boundaryIntegrator.computeBoundaryIntegral( l_timeIntegratedUnknownsElement,
                                                      l_timeIntegratedUnknownsNeighborPointers,
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
          l_simpleBoundaryIntegrator.computeBoundaryIntegration( l_timeIntegratedUnknownsElement[0],
                                                                 l_timeIntegratedUnknownsNeighbors,
                                                                 l_boundaryConditions,
                                                                 l_neighboringIndices,
                                                                 l_fluxSolversPos,
                                                                 l_fluxSolversNeg,
                                                                 l_unknownsUnitTests );

          // do boundary integration with generated kernels
          l_boundaryIntegrator.computeBoundaryIntegral( l_timeIntegratedUnknownsElement,
                                                        l_timeIntegratedUnknownsNeighborPointers,
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
