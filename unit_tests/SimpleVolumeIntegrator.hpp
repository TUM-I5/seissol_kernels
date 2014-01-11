/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer) *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2013-2014
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
 * Simple volume integrator.
 **/

#ifndef SIMPLEVOLUMEINTEGRATOR_HPP_
#define SIMPLEVOLUMEINTEGRATOR_HPP_

#include <Initializer/typedefs.hpp>

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
    double m_stiffnessMatrices[3][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS];

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
                                          NUMBEROFBASISFUNCTIONS,
                                          NUMBEROFBASISFUNCTIONS,
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
    void computeVolumeIntegration( const double i_timeIntegratedUnknowns[NUMBEROFUNKNOWNS],
                                   const double i_starMatrices[3][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                         double io_unknowns[NUMBEROFUNKNOWNS] ) {
      // temporary matrix for two-step multiplications
      double l_temporaryProduct[NUMBEROFUNKNOWNS];

      // iterate over the three reference coordinate: \f$ \xi, \eta, \zeta \f$
      for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) {
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0 );

        // stiffness matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                     m_stiffnessMatrices[ l_coordinate ], i_timeIntegratedUnknowns, l_temporaryProduct );

        // star matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                     l_temporaryProduct, i_starMatrices[ l_coordinate ], io_unknowns );
      }
    }
};
#endif
