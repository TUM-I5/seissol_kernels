/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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
 * Simple time integrator.
 **/

#ifndef SIMPLETIMEINTEGRATOR_HPP_
#define SIMPLETIMEINTEGRATOR_HPP_

#include <Initializer/typedefs.hpp>
#include <seissol/Initializer/preProcessorMacros.fpp>

#include "DenseMatrix.hpp"

namespace unit_test {
  class SimpleTimeIntegrator;
};

/**
 * Simple implementation of the ADER-time integration in SeisSol.
 **/
class unit_test::SimpleTimeIntegrator {
  //private:
    //! transposed stiffness matrices (multiplied by the inverse of the mass matrix): \f$ M^{-1} K^\xi, M^{-1} K^\eta, M^{-1} K^\zeta \f$
    real m_transposedStiffnessMatrices[3][NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS];

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;
  public:
    /**
     * Initializes the time integrator with the given file.
     *
     * @param i_matricesPath path to matrices XML-file.
     **/
    void initialize( std::string &i_matricesPath ) {
      m_denseMatrix.readMatrices( i_matricesPath );

      // iterate over the three coordinates \f$ \xi, \eta, \zeta \f$
      for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) { 
        // initialize the transposed stiffness matrix (multiplied by the inverse of the mass matrix)
        m_denseMatrix.initializeMatrix( 56+l_coordinate,
                                        NUMBEROFBASISFUNCTIONS,
                                        NUMBEROFBASISFUNCTIONS,
                                        m_transposedStiffnessMatrices[l_coordinate] );
      }
    }

    /** 
     * Computes the time derivatives of the given unknowns.
     *
     * @param i_unknowns unknowns for which the time derivatives are computed.
     * @param i_aStar star matrix \f$ A^*_k \f$
     * @param i_bStar star matrix \f$ B^*_k \f$
     * @param i_cStar star matrix \f$ C^*_k \f$
     * @param o_timeDerivatives time derivates.
     **/
    void computeTimeDerivation( const real i_unknowns[NUMBEROFUNKNOWNS],
                                const real i_aStar[NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                const real i_bStar[NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                const real i_cStar[NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                      real o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] ) {
      // temporary matrix for two-step multiplication
      real l_temporaryProduct[NUMBEROFUNKNOWNS];
     
      // copy unknowns to zeroth derivative
      for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) { 
        o_timeDerivatives[0][l_entry] = i_unknowns[l_entry];
      }   

      // iterate over order in time
      for( int l_order = 1; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) { 
        // set time derivatives of this order to zero
        std::fill( o_timeDerivatives[l_order], o_timeDerivatives[l_order+1], 0); 

        /*
         * reference coordinate: \f$ \xi \f$
         */
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0); 

        // stiffness matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                     m_transposedStiffnessMatrices[0], o_timeDerivatives[l_order-1], l_temporaryProduct );

        // star matrix multiplcation
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                     l_temporaryProduct, i_aStar, o_timeDerivatives[l_order] );   
                              
        /*
         * reference coordinate: \f$ \eta \f$
         */
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0); 
        
        // stiffness matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                     m_transposedStiffnessMatrices[1], o_timeDerivatives[l_order-1], l_temporaryProduct ); 

        // star matrix multiplcation
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                     l_temporaryProduct, i_bStar, o_timeDerivatives[l_order] );                                      
        /*
         * reference coordinate: \f$ \zeta \f$
         */
        // set temporary product to zero
        std::fill( l_temporaryProduct, l_temporaryProduct+NUMBEROFUNKNOWNS, 0); 
        
        // stiffness matrix multiplication
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFBASISFUNCTIONS,
                                                     m_transposedStiffnessMatrices[2], o_timeDerivatives[l_order-1], l_temporaryProduct ); 

        // star matrix multiplcation
        m_denseMatrix.executeStandardMultiplication( NUMBEROFBASISFUNCTIONS, NUMBEROFVARIABLES, NUMBEROFVARIABLES,
                                                     l_temporaryProduct, i_cStar, o_timeDerivatives[l_order] );
      } 
    }

    /**
     * Computes the time integrals of the given unknowns.
     *
     * @param i_timeDerivatives time derivatives of the unknowns in the cell.
     * @param i_deltaT integration boundaries \f$ \Delta t^\text{lo} \f$ and \f$ \Delta t^\text{up} \f$ relative to the time level \f$ \t^\text{cell} \f$ of the cell; integration is performned over the interval \f$ [ t^\text{cell} + \Delta t^\text{lo}, t^\text{cell} + \Delta t^\text{up} ] \f$
     * @param o_timeIntegratedUnknowns unknowns integrated over the specified interval.
     **/
    void computeTimeIntegration( const real i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                 const real i_deltaT[2],
                                       real o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
      // assert integration forward in time
      TS_ASSERT( i_deltaT[1] > i_deltaT[0] );

      // initialization of Taylor series scalars
      real l_firstTerm  = (real)  1;
      real l_secondTerm = (real)  1;
      real l_factorial  = (real) -1;
      real l_scalar;

      // reset time integrated deegrees of freeomd with zeroth derivatives
      for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
        o_timeIntegratedUnknowns[l_entry] = 0;
      }

      // iterate over order in time
      for( int l_order = 0; l_order < ORDEROFTAYLORSERIESEXPANSION; l_order++ ) {
        // compute factor of the taylor series
        l_firstTerm  *=  i_deltaT[1];
        l_secondTerm *=  i_deltaT[0];
        l_factorial  *= -(real)(l_order+1);

        l_scalar  = l_firstTerm - l_secondTerm;
        l_scalar /= l_factorial;
       
        // update time integrated unknowns
        for( int l_entry = 0; l_entry < NUMBEROFUNKNOWNS; l_entry++ ) {
          o_timeIntegratedUnknowns[l_entry] += l_scalar * i_timeDerivatives[l_order][l_entry];
        } 
      }
    }

   /**
    * Compute the extrapolation of the DOFs based on the Taylor series expansion. Extrapolation is computed at the expansion point plus a \f$ \Delta_t \f$.
    *
    * @param i_timeDerivatives time derivatives of the Taylor series expansion.
    * @param i_deltaT \f$ \Delta_t \f$ from the expansion point to the extrapolation point.
    * @param o_timeExtrapolated time extrapolated degrees of freedom.
    **/
   void computeExtrapolation( const real i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                                    real i_deltaT,
                                    real o_timeExtrapolated[NUMBEROFUNKNOWNS] ) {
      // copy 0th derivative to extrapolation
      for( unsigned int l_dof = 0; l_dof < NUMBEROFUNKNOWNS; l_dof++ ) {
        o_timeExtrapolated[l_dof] = i_timeDerivatives[0][l_dof];
      }

      real l_scalar = (real) 1;

      // iterate of all derivatives
      for( unsigned int l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
        l_scalar *= -i_deltaT;
        l_scalar /= (real) l_derivative;

        // add the contribution of this derivative
        for( unsigned int l_dof =0; l_dof < NUMBEROFUNKNOWNS; l_dof++ ) {
          o_timeExtrapolated[l_dof] += l_scalar * i_timeDerivatives[l_derivative][l_dof];
        }
      }
   }

};
#endif
