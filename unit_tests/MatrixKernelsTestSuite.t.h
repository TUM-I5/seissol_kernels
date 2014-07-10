/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2012, SeisSol Group
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
 * Test suite, which tests generated, plain matrix kernels.
 **/
#include <vector>
#include <string>
#include <sstream>

#include <cxxtest/TestSuite.h>
#include <utils/logger.h>

#include <seissol/Initializer/MemoryAllocator.h>
#include <seissol/Initializer/MemoryManager.h>
#include "../src/TimeIntegrator.h"

#include "configuration.hpp"
#include "SparseMatrix.hpp"
#include <seissol/Initializer/preProcessorMacros.fpp>

// extern declaration of SeisSol's "old" CauchKovalewski3D routine.
extern "C" void c_bind_cauchy_kovalewski_3d( const double* i_unknowns,
                                               const int    &i_numberOfBasisFunctions,
                                               const int    &i_numberOfQuantities,
                                               const int    &i_orderOfTaylorSeriesExpansion,
                                               const double* i_aStarDivM,
                                               const double* i_bStarDivM,
                                               const double* i_cStarDivM,
                                               const double* i_kXi,
                                               const double* i_kEta,
                                               const double* i_kZeta,
                                               const double    &i_deltaT,
                                                     double* o_timeIntegratedUnknowns );

namespace unit_test {
  class MatrixTestSuite;
}

class unit_test::MatrixTestSuite: public CxxTest::TestSuite {
  //private:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    //! How often should each unit test with random numbers repeated?
    int m_numberOfRepeats;

    //! vector of strings containing star matrices
    std::vector<SparseMatrix> m_starMatrices;

    //! vector of strings containing stiffness matrices
    std::vector<SparseMatrix> m_stiffnessMatrices;

    //! vector of strings containing flux matrices
    std::vector<SparseMatrix> m_fluxMatrices;

    //! location of the sparse matrix files
    std::string m_matricesDirectory;

    //! minimum order to test
    unsigned int m_minimumOrder;

    //! maximum order to test
    unsigned int m_maximumOrder;

    /**
     * Parse stiffness matrices (up to the given order) and star matrices.
     *
     * @param i_maximumOrder maximum order of the stiffness matrices
     */
    void parseMatrices( int i_maximumOrder ) {
      // compute maximum degree
      int l_maximumDegree = i_maximumOrder - 1;

      // currently there's only one star matrix (pure elastics)
      //   Remark: The Star matrix is parsed three time to simulate A*, B* and C*
      m_starMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/starMatrix_3D_csc.mtx", SparseMatrix::StarMatrix) );
      m_starMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/starMatrix_3D_csc.mtx", SparseMatrix::StarMatrix) );
      m_starMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/starMatrix_3D_csc.mtx", SparseMatrix::StarMatrix) );

      // add stiffness matrices for different number of basis functions basis functions
      for(int l_degree = m_minimumOrder-1; l_degree <= l_maximumDegree; l_degree++) {
        // convert integer to string
        std::stringstream l_stringStream;
        l_stringStream << l_degree;
        std::string l_degreeToString = l_stringStream.str();

        // add stiffness matrices for 3 basis functions
        m_stiffnessMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/kXiDivMT_3D_"  +l_degreeToString+"_csc.mtx", SparseMatrix::TransposedStiffnessMatrix) );
        m_stiffnessMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/kEtaDivMT_3D_" +l_degreeToString+"_csc.mtx", SparseMatrix::TransposedStiffnessMatrix) );
        m_stiffnessMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+"/kZetaDivMT_3D_"+l_degreeToString+"_csc.mtx", SparseMatrix::TransposedStiffnessMatrix) );

        // add flux matrices

        // file name of flux matrices
        std::string l_fileName;

        for( int l_faceLocal = 1; l_faceLocal <= 4; l_faceLocal++) {
          // clear stringstream
          l_stringStream.str("");
          l_stringStream.clear();
              
          // create function name string
          l_stringStream << "fM" << l_faceLocal << "DivM_3D_" +l_degreeToString+"_csc.mtx";
          l_fileName = l_stringStream.str();

          // add flux matrix for element local contribution
          m_fluxMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+'/'+l_fileName, SparseMatrix::FluxMatrix) );

          for( int l_neighboringFace = 1; l_neighboringFace <= 4; l_neighboringFace++) {
            for( int l_vertexOrientation = 1; l_vertexOrientation <= 3; l_vertexOrientation++ ) {
              // clear stringstream
              l_stringStream.str("");
              l_stringStream.clear();
              
              // create function name string
              l_stringStream << "fP" << l_faceLocal << l_neighboringFace << l_vertexOrientation << "DivM_3D_" << l_degreeToString + "_csc.mtx";
              std::string l_fileName = l_stringStream.str();

              // add flux matrix for neighboring element contribution
              m_fluxMatrices.push_back( unit_test::SparseMatrix(m_matricesDirectory+'/'+l_fileName, SparseMatrix::FluxMatrix) );
            }
          }
        }

      }
    }

    // reset matrices
    void deleteMatrices() {
      m_starMatrices.clear();
      m_stiffnessMatrices.clear();
      m_fluxMatrices.clear();
    }

  public:
    /**
     * Setup local configuration before every unit test.
     */
    void setUp() {
      // get number of reapeats.
      m_numberOfRepeats = m_configuration.getNumberOfRepeats();

      // get matrices directory
      m_matricesDirectory = m_configuration.getMatricesDirectory();

      // get min. and max. order
      m_minimumOrder = m_configuration.getMinimumOrder();
      m_maximumOrder = m_configuration.getMaximumOrder();

      // parse matrices
      parseMatrices( m_maximumOrder );
    }

    /**
     * Delete matrices after every unit test.
     */
    void tearDown() {
      deleteMatrices();
    }

    /**
     * Unit test for generated matrix kernels containing the star matrix multiplications (dense*sparse).
     */
    void testGeneratedStarMultiplication() {
// TODO: Add call to logger.
#if 0
      std::cout << std::endl << "*** testing star matrix multiplications ***" << std::endl;
#endif
      for( int l_i = 0; l_i < m_starMatrices.size(); l_i++ ) {
// TODO: Add call to logger.
#if 0
        std::cout << "testing: " << (m_starMatrices[l_i]).pathToMatrixMarketFile << std::endl;
#endif

        // repeat the test
        for( int l_numberOfTest = 0; l_numberOfTest < m_numberOfRepeats; l_numberOfTest++ ) {
          TS_ASSERT(m_starMatrices[l_i].testGeneratedMultiplication());
        }
      }
    }

    /**
     * Unit test for matrix kernels containing the stiffness matrix multiplications (sparse*dense).
     */
    void testGeneratedStiffnessMultiplication() {
// TODO: Add call to logger.
#if 0
      std::cout << std::endl << "*** testing stiffness matrix multiplications ***" << std::endl;
#endif
      for( int l_i = 0; l_i < m_stiffnessMatrices.size(); l_i++ ) {
// TODO: Add call to logger.
#if 0
        std::cout << "testing: " << m_stiffnessMatrices[l_i].pathToMatrixMarketFile << std::endl;
#endif

        // repeat the test
        for( int l_numberOfTests = 0; l_numberOfTests < m_numberOfRepeats; l_numberOfTests++ ) {
          TS_ASSERT(m_stiffnessMatrices[l_i].testGeneratedMultiplication());
        }
      }
    }

    /**
     * Unit test for matrix kernels containing the flux matrix multiplications (sparse*dense).
     */
    void testGeneratedFluxMultiplication() {
// TODO: Add call to logger.
#if 0
      std::cout << std::endl << "*** testing flux matrix multiplications ***" << std::endl;
#endif

      for( int l_i = 0; l_i < m_fluxMatrices.size(); l_i++ ) {
// TODO: Add call to logger.
#if 0
        std::cout << "testing: " << m_fluxMatrices[l_i].pathToMatrixMarketFile << std::endl;
#endif

        // repeat the test
        for( int l_numberOfTests = 0; l_numberOfTests < m_numberOfRepeats; l_numberOfTests++ ) {
          TS_ASSERT(m_fluxMatrices[l_i].testGeneratedMultiplication());
        }
      }
    }

    /**
     * Unit test for the Fortran Cauchy Kovalweski kernel, which calls generated C matrix multiplications.
     */
    void testCauchyKovalewskiKernel() {
      int l_numberOfBasisFunctions       = NUMBEROFBASISFUNCTIONS;
      int l_numberOfVariables            = NUMBEROFVARIABLES;
      int l_orderOfTaylorSeriesExpansion = ORDEROFTAYLORSERIESEXPANSION+1;

      //! setup matrix path
      std::string l_matricesPath = m_configuration.getMatricesDirectory() + "/matrices_" + std::to_string(NUMBEROFBASISFUNCTIONS) + ".xml";

      // set up the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );
      // create a new memory manager
      seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );
      // set up the volume integrator
      seissol::kernels::TimeIntegrator l_timeIntegrator( l_matrixReader, l_memoryManager );      

// TODO: Add call to logger.
#if 0
      std::cout
        << std::endl
        << "*** testing Cauchy Kovalewski kernel with (#basis functions, #variables, order of Taylor series expansion): ("
        << l_numberOfBasisFunctions << ", " << l_numberOfVariables << ", " << l_orderOfTaylorSeriesExpansion << ")" << std::endl;
#endif

      int l_stiffnessPosition;
      if (NUMBEROFBASISFUNCTIONS == 4) {
        l_stiffnessPosition = 2 - m_minimumOrder;
      }
      else if (NUMBEROFBASISFUNCTIONS == 10) {
        l_stiffnessPosition = 5 - m_minimumOrder;
      }
      else if (NUMBEROFBASISFUNCTIONS == 20) {
        l_stiffnessPosition = 8 - m_minimumOrder;
      }
      else if (NUMBEROFBASISFUNCTIONS == 35) {
        l_stiffnessPosition = 11 - m_minimumOrder;
      }
      else if (NUMBEROFBASISFUNCTIONS == 56) {
        l_stiffnessPosition = 14 - m_minimumOrder;
      }
      else if (NUMBEROFBASISFUNCTIONS == 84) {
        l_stiffnessPosition = 17 - m_minimumOrder;
      }
      else {
        TS_FAIL("number of basis functions not in {4, 10, 20, 35, 56, 84}.");
        return;
      }

      for( int l_numberOfTests = 0; l_numberOfTests < m_numberOfRepeats; l_numberOfTests++ ) {
        // temporary pointers
        double *l_unknowns;
        double *l_aStarDivMDense, *l_bStarDivMDense, *l_cStarDivMDense;
        double *l_aStarDivMFlat, *l_bStarDivMFlat, *l_cStarDivMFlat;
        double *l_kXiDense, *l_kEtaDense, *l_kZetaDense;
        double *l_kXiFlat, *l_kEtaFlat, *l_kZetaFlat;
        double l_deltaT;
        double *l_timeIntegratedUnknownsCK3D, *l_timeIntegratedUnknownsGen;

        // allocate memory for unknowns and time integrated unknowns (64 byte aligned)
        seissol::MemoryAllocator l_memoryAllocator;
        l_memoryAllocator.allocateMemory( l_numberOfBasisFunctions*l_numberOfVariables*sizeof(double),
                                          64,
                                          (void**) &l_unknowns );
        l_memoryAllocator.allocateMemory( l_numberOfBasisFunctions*l_numberOfVariables*sizeof(double),
                                          64,
                                          (void**) &l_timeIntegratedUnknownsCK3D );
        l_memoryAllocator.allocateMemory( l_numberOfBasisFunctions*l_numberOfVariables*sizeof(double),
                                          64,
                                          (void**) &l_timeIntegratedUnknownsGen );

        // set random values to star matrices, unknowns and time integrated unknowns
        m_starMatrices[0].setRandomValues();
        m_starMatrices[1].setRandomValues();
        m_starMatrices[2].setRandomValues();
        SparseMatrix::setRandomValues( l_numberOfBasisFunctions*l_numberOfVariables, l_unknowns                   );
        SparseMatrix::setRandomValues( l_numberOfBasisFunctions*l_numberOfVariables, l_timeIntegratedUnknownsCK3D );
        SparseMatrix::setRandomValues( l_numberOfBasisFunctions*l_numberOfVariables, l_timeIntegratedUnknownsGen  );


        // assign dense random star and stiffness matrices
        l_aStarDivMDense = m_starMatrices[0].denseColumnMajor;
        l_bStarDivMDense = m_starMatrices[1].denseColumnMajor;
        l_cStarDivMDense = m_starMatrices[2].denseColumnMajor;

        l_kXiDense =   m_stiffnessMatrices[l_stiffnessPosition  ].denseColumnMajor;
        l_kEtaDense =  m_stiffnessMatrices[l_stiffnessPosition+1].denseColumnMajor;
        l_kZetaDense = m_stiffnessMatrices[l_stiffnessPosition+2].denseColumnMajor;

        // set time width constant
        l_deltaT = 0.05;

        // call SeisSol's CK kernel
        c_bind_cauchy_kovalewski_3d( l_unknowns,
                                     l_numberOfBasisFunctions,
                                     l_numberOfVariables,
                                     l_orderOfTaylorSeriesExpansion,
                                     l_aStarDivMDense, l_bStarDivMDense, l_cStarDivMDense,
                                     l_kXiDense, l_kEtaDense, l_kZetaDense,
                                     l_deltaT,
                                     l_timeIntegratedUnknownsCK3D );

        // assign flat column majo matrices
        l_aStarDivMFlat = m_starMatrices[0].flatColumnMajor;
        l_bStarDivMFlat = m_starMatrices[1].flatColumnMajor;
        l_cStarDivMFlat = m_starMatrices[2].flatColumnMajor;

        l_kXiFlat =   m_stiffnessMatrices[l_stiffnessPosition  ].flatColumnMajor;
        l_kEtaFlat =  m_stiffnessMatrices[l_stiffnessPosition+1].flatColumnMajor;
        l_kZetaFlat = m_stiffnessMatrices[l_stiffnessPosition+2].flatColumnMajor;

        // call generated CCKkernel
        l_timeIntegrator.computeTimeIntegral( l_unknowns,
                                              l_aStarDivMFlat,
                                              l_bStarDivMFlat,
                                              l_cStarDivMFlat,
                                              l_deltaT,
                                              l_timeIntegratedUnknownsGen );

        // compute absolute error in infinity norm
        double l_error =  SparseMatrix::compare_colmajor(l_timeIntegratedUnknownsCK3D, l_timeIntegratedUnknownsGen , l_numberOfBasisFunctions * l_numberOfVariables );

        // abort if the error is too large
        if( l_error > 10e-12 ) {
          std::cout << "relative infinity norm error: " << l_error << std::endl;
          std::cout << "-- \"old\" routine --" << std::endl;
          SparseMatrix::print_colmajor(l_timeIntegratedUnknownsCK3D, l_numberOfBasisFunctions, l_numberOfVariables );
          std::cout << "-- generated routine --" << std::endl;
          SparseMatrix::print_colmajor(l_timeIntegratedUnknownsGen, l_numberOfBasisFunctions, l_numberOfVariables );
          TS_FAIL( "relative infinity norm error to large. check relative error! (see result matrices above)" );
        }
      }
    }
};
