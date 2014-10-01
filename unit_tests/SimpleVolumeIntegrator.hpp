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
 * Simple volume integrator.
 **/

#ifndef SIMPLEVOLUMEINTEGRATOR_HPP_
#define SIMPLEVOLUMEINTEGRATOR_HPP_

#include <Initializer/typedefs.hpp>
#include <seissol/Initializer/preProcessorMacros.fpp>

#include "DenseMatrix.hpp"

namespace unit_test {
  class SimpleVolumeIntegrator;
};

/**
 * Simple implementation of the volume integration in SeisSol.
 **/
class unit_test::SimpleVolumeIntegrator {
  //private:
    //! stiffness matrices (multiplied by the inverse of the mass matrix): \f$ M^{-1} K^\xi, M^{-1} K^\eta, M^{-1} K^\zeta \f$
    double m_stiffnessMatrices[3][NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_BASIS_FUNCTIONS];

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;
  public:
    /**
     * Initializes the volume integrator.
     *
     * @param i_matricesPath path to the matrix xml file.
     **/
    void initialize( std::string &i_matricesPath ) {
      // read in the matrices for our unit tests
      m_denseMatrix.readMatrices( i_matricesPath );

        // iterate over the three coordinates \f$ \xi, \eta, \zeta \f$
        for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) {
          // inistialize unit test stiffness matrix (multiplied by the inverse of the mass matrix)
          m_denseMatrix.initializeMatrix( 53+l_coordinate,
                                          NUMBER_OF_BASIS_FUNCTIONS,
                                          NUMBER_OF_BASIS_FUNCTIONS,
                                          m_stiffnessMatrices[l_coordinate] );
        }
    }
 
 
    /**
     * Simple volume integration to test again.
     *
     * @param i_timeIntegratedUnknowns time integrated unknowns of the element \f$ k \f$ over the interval \f$ [t^n, t^{n+1}] \f$:
     *        \f$ I(t^n, t^{n+1}, Q_k^n) \f$.
     * @param i_starMatrices star matrices:
     *        0: \f$ A^*_k \f$
     *        1: \f$ B^*_k \f$
     *        2: \f$ C^*_k \f$
     * @param io_unknowns unknowns of the previous time step \f$ t^n \f$, which will be updated by the volume integral:
     *        \f[
     *            M^{-1} K^\xi   I(t^{n}, t^{n+1}, Q_{k}^n) A^*_k
     *          + M^{-1} K^\eta  I(t^{n}, t^{n+1}, Q_{k}^n) B^*_k
     *          + M^{-1} K^\zeta I(t^{n}, t^{n+1}, Q_{k}^n) C^*_k
     *        \f]
     **/
    void computeVolumeIntegration( const double i_timeIntegratedUnknowns[NUMBER_OF_DOFS],
                                   double const i_starMatrices[3][NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                                   double       io_unknowns[NUMBER_OF_DOFS] ) {
      // temporary matrix for two-step multiplications
      double l_temporaryProduct[NUMBER_OF_DOFS];

      // iterate over the three reference coordinate: \f$ \xi, \eta, \zeta \f$
      for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) {
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBER_OF_DOFS, 0 );

        // stiffness matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
                                                     m_stiffnessMatrices[ l_coordinate ], i_timeIntegratedUnknowns, l_temporaryProduct );

        // star matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
                                                     l_temporaryProduct, i_starMatrices[ l_coordinate ], io_unknowns );
      }
    }
};
#endif
