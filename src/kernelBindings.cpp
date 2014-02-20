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
 * Binds the boundary and volume integrator to C for usage in Fortran.
 **/
#ifdef __INTEL_OFFLOAD
#ifdef __MIC__
#define DIRTY_EXCLUDE_ON_MIC
#endif
#endif

#ifndef DIRTY_EXCLUDE_ON_MIC

#include "BoundaryIntegrator.h"
#include "VolumeIntegrator.h"
#include "TimeIntegrator.h"

// setup path of the XML-file, which contains information about the matrices.
// TODO: Add this to input-parameters of SeisSol.
#ifdef __INTEL_OFFLOAD
#define MATRIXXMLFILE "matrices_" STR(NUMBEROFBASISFUNCTIONS) ".offload_mic.xml"
#else
#ifdef __MIC__
#define MATRIXXMLFILE "matrices_" STR(NUMBEROFBASISFUNCTIONS) ".mic.xml"
#else
#define MATRIXXMLFILE "matrices_" STR(NUMBEROFBASISFUNCTIONS) ".xml"
#endif
#endif

// C-function name of the boundary integration.
#define BOUNDARYINTFUNCTIONNAME CONCAT_4( boundaryIntegration_,         NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS )

// C-function name of the volume integration
#define VOLUMEINTFUNCTIONNAME CONCAT_4( volumeIntegration_,             NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS )

// C-function name of the time integration
#define TIMEINTFUNCTIONNAME CONCAT_4( cauchyKovalewskiTimeIntegration_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS )

// C-funcion name of the time derivative computation
#define TIMEDERFUNCTIONNAME CONCAT_4( computeTimeDerivatives_,          NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS )

// C-function name of the time integration, which uses previously computed time derivatives and allows for LTS.
#define LTSTIMEINTFUNCTIONNAME CONCAT_4( computeTimeIntegrationLTS_, NUMBEROFVARIABLES, _, NUMBEROFBASISFUNCTIONS )

// set up the xml-parser
seissol::XmlParser l_matrixReader( MATRIXXMLFILE );

// create a new memory manager
seissol::initializers::MemoryManager l_memoryManager( l_matrixReader );

// set up the boundary integrator
seissol::kernels::BoundaryIntegrator l_boundaryIntegrator( l_matrixReader, l_memoryManager );

// set up the volume integrator
seissol::kernels::VolumeIntegrator   l_volumeIntegrator(   l_matrixReader, l_memoryManager );

// set up the time integrator
seissol::kernels::TimeIntegrator     l_timeIntegrator(     l_matrixReader, l_memoryManager );

// prevent name mangling
extern "C" {
  /**
   * Simple forward of the boundary integration. Details can be found in the BoundaryIntegrator-class.
   **/
  void BOUNDARYINTFUNCTIONNAME(       double  i_timeIntegratedUnknownsElement[2][NUMBEROFUNKNOWNS],
                                      double *i_timeIntegratedUnknownsNeighbors[4],
                                const int     i_boundaryConditions[4],
                                const int     i_neighboringIndices[4][2],
                                      double  i_nApNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                      double  i_nAmNm1[4][NUMBEROFVARIABLES*NUMBEROFVARIABLES],
                                      double  io_unknowns[NUMBEROFUNKNOWNS] ) {
    l_boundaryIntegrator.computeBoundaryIntegral( i_timeIntegratedUnknownsElement,
                                                  i_timeIntegratedUnknownsNeighbors,
                                                  i_boundaryConditions,
                                                  i_neighboringIndices,
                                                  i_nApNm1,
                                                  i_nAmNm1,
                                                  io_unknowns );
  }

  /**
   * Simple forward of the volume integration. Details can be found in the VolumeIntegrator-class.
   **/
  void VOLUMEINTFUNCTIONNAME( double i_timeIntegratedUnknowns[NUMBEROFUNKNOWNS],
                              double i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                              double i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                              double i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                              double io_unknowns[NUMBEROFUNKNOWNS] ) {
    l_volumeIntegrator.computeVolumeIntegral( i_timeIntegratedUnknowns,
                                              i_aStar,
                                              i_bStar,
                                              i_cStar,
                                              io_unknowns );
  }

  /**
   * Simple forward of the time integration. Details can be found in the TimeIntegrator-class.
   **/
  void TIMEINTFUNCTIONNAME( const double  i_unknowns[NUMBEROFUNKNOWNS],
                                  double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                  double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                  double  i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                            const double &i_deltaT,
                                  double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
    l_timeIntegrator.computeTimeIntegral( i_unknowns,
                                          i_aStar,
                                          i_bStar,
                                          i_cStar,
                                          i_deltaT,
                                          o_timeIntegratedUnknowns );
  }

  /**
   * Simple forward of the time derivative computation. Details can be found in the TimeIntegrator-class.
   **/
  void TIMEDERFUNCTIONNAME( const double  i_unknowns[NUMBEROFUNKNOWNS],
                                  double  i_aStar[STARMATRIX_NUMBEROFNONZEROS],
                                  double  i_bStar[STARMATRIX_NUMBEROFNONZEROS],
                                  double  i_cStar[STARMATRIX_NUMBEROFNONZEROS],
                                  double  o_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS] ) {
    l_timeIntegrator.computeTimeDerivatives( i_unknowns,
                                             i_aStar,
                                             i_bStar,
                                             i_cStar,
                                             o_timeDerivatives );
  }

  /**
   * Simple forward of the time integration, which uses prevoiusly computed time derivatives.
   * Details can be found in the TimeIntegrator-class.
   **/
  void LTSTIMEINTFUNCTIONNAME( const double  i_timeDerivatives[ORDEROFTAYLORSERIESEXPANSION][NUMBEROFUNKNOWNS],
                               const double &i_deltaTLower,
                               const double &i_deltaTUpper,
                                     double  o_timeIntegratedUnknowns[NUMBEROFUNKNOWNS] ) {
    l_timeIntegrator.computeTimeIntegral( i_timeDerivatives,
                                          i_deltaTLower,
                                          i_deltaTUpper,
                                          o_timeIntegratedUnknowns );
  }
}
#endif
