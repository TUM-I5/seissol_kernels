/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2012
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
 * Interface between the Cauchy Kovalewski / ADER kernel and benchit.
 **/
 
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "interface.h"

#include <sstream>

#include "cauchyKovalewskiTimeIntegrationKernel.hpp"
#include "seissol_src/Monitoring/SeisSolStopwatch.hpp"


//! used benchmark data
BenchmarkData l_benchmarkData;

//! SeisSol stopwatch
SeisSolStopwatch l_stopWatch;

/* Reads the environment variables used by this kernel. */
/**
 * Read used environmnet variables
 *
 * @param o_benchmarkData benchmark data, where the variables are stored to.
 */
void evaluate_environment(BenchmarkData &o_benchmarkData) {
  int errors = 0;
  char * p = 0;
  
  // read minimum number of elements in the pseudo-discretization.
  p = bi_getenv( "BENCHIT_KERNEL_PROBLEMSIZE_MIN", 0 );
  if ( p == NULL )
    errors++;
  else
    o_benchmarkData.minimumNumberOfElements = atoi( p );
  
  // read maximum number of elements in the pseudo-discretization.
  p = bi_getenv( "BENCHIT_KERNEL_PROBLEMSIZE_MAX", 0 );
  if ( p == NULL )
    errors++;
  else
    o_benchmarkData.maximumNumberOfElements = atoi( p );
  
  // read increment between different mesh sizes.
  p = bi_getenv( "BENCHIT_KERNEL_PROBLEMSIZE_INCREMENT", 0 );
  if ( p == NULL )
    errors++;
  else o_benchmarkData.incrementPerStep = atoi( p );
 
  // aboirt in the case of errors.
  if ( errors > 0 ) {
    std::cerr << "There's at least one environment variable not set!" << std::endl;
    exit( 1 );
  }
  
  // compute number of steps
  o_benchmarkData.numberOfSteps = (int) ( o_benchmarkData.maximumNumberOfElements -  o_benchmarkData.minimumNumberOfElements + 1 ) / o_benchmarkData.incrementPerStep;
  if (( o_benchmarkData.maximumNumberOfElements - o_benchmarkData.minimumNumberOfElements + 1 ) % o_benchmarkData.incrementPerStep != 0)
    o_benchmarkData.incrementPerStep++;
}

/** 
 * The implementation of the bi_getinfo from the BenchIT interface.
 *  Here the infostruct is filled with informations about the
 *  kernel.
 *  @param infostruct pointer to a structure filled with zero's
 */
void bi_getinfo( bi_info * pinfo ) {
  // set up kernel description
  std::stringstream l_kernelDescription;
  l_kernelDescription << "Cauchy Kovalewski / ADER kernel with " << NUMBEROFBASISFUNCTIONS
                      << "basis functions and " << NUMBEROFVARIABLES << " variables";

  (void) memset ( pinfo, 0, sizeof( bi_info ) );

  // get environment variables for the kernel
  evaluate_environment(l_benchmarkData);
  
  // set information
  pinfo->codesequence = bi_strdup( "start kernel; do nothing; " );
  pinfo->kerneldescription = bi_strdup( l_kernelDescription.str().c_str() );
  pinfo->xaxistext = bi_strdup( "number of elements" );
  pinfo->maxproblemsize = l_benchmarkData.numberOfSteps;
  pinfo->num_processes = 1;
  pinfo->num_threads_per_process = 0;
  pinfo->kernel_execs_mpi1 = 0;
  pinfo->kernel_execs_mpi2 = 0;
  pinfo->kernel_execs_pvm = 0;
  pinfo->kernel_execs_omp = 0;
  pinfo->kernel_execs_pthreads = 0;
  
#ifdef __INTEL_COMPILER
  // add cycles functions if the Intel compiler is used.
  pinfo->numfunctions = 2;
#else
  // pure time information in case of all other compilers
  pinfo->numfunctions = 1;
#endif

  // allocate memory for y-axis texts and properties.
  pinfo->yaxistexts = (char**) malloc( pinfo->numfunctions * sizeof( char* ) );
  if ( pinfo->yaxistexts == NULL ) {
    std::cerr << "Allocation of yaxistexts failed." << std::endl << std::flush;
    exit( 127 );
  }
  
  pinfo->outlier_direction_upwards = (int*) malloc( pinfo->numfunctions * sizeof( int ) );
  if ( pinfo->outlier_direction_upwards == NULL ) {
    std::cerr << "Allocation of outlier direction failed." << std::endl << std::flush;
    exit( 127 );
  }
  pinfo->legendtexts = (char**) malloc( pinfo->numfunctions * sizeof( char* ) );
  if ( pinfo->legendtexts == NULL ) {
    std::cerr << "Allocation of legendtexts failed." << std::endl << std::flush;
    exit( 127 );
  }
  pinfo->base_yaxis = (double*) malloc( pinfo->numfunctions * sizeof( double ) );
  if ( pinfo->base_yaxis == NULL ) {
    std::cerr << "Allocation of base yaxis failed." << std::endl << std::flush;
    exit( 127 );
  }
  
  // set up y axis texts and properties.
  pinfo->yaxistexts[0] = bi_strdup( "time in s" );
  pinfo->outlier_direction_upwards[0] = 1;
  pinfo->base_yaxis[0] = 0; // linear scale
  pinfo->legendtexts[0] = bi_strdup( "time in s" );

  // set up second y-axis (CPU-cycles) for the Intel compiler.
#ifdef __INTEL_COMPILER  
  pinfo->yaxistexts[1] = bi_strdup( "cpu cycles" );
  pinfo->outlier_direction_upwards[1] = 1;
  pinfo->base_yaxis[1] = 0; // linear scale
  pinfo->legendtexts[1] = bi_strdup( "cpu cycles" );
#endif

}



/**
 * Implementation of the bi_init of the BenchIT interface.
 */
void* bi_init( int ) {
  // allocate memory for the given elements
  l_benchmarkData.initializeDiscontinousGalerkinElements( l_benchmarkData.maximumNumberOfElements ); 
  
  // print information about the benchmark
  std::cout << std::endl
            << "ADER kernel:" << std::endl
            << "  maximum(#elements)  : " <<  l_benchmarkData.maximumNumberOfElements << std::endl
            << "  minimum(#elements)  : " <<  l_benchmarkData.minimumNumberOfElements << std::endl
            << "  increcement per step: " <<  l_benchmarkData.incrementPerStep        << std::endl
            << "  #steps              : " <<  l_benchmarkData.numberOfSteps           << std::endl;
  
  return (void*) NULL;
}



/** The central function within each kernel. This function
 *  is called for each measurment step seperately.
 *  @param  o_results pointer to a field of doubles, the size of the field depends on the number
 *                    of functions, there are #functions+1 doubles.
 *  @return 0 if the measurment was sucessfull, something else in the case of an error.
 */
int bi_entry( void*, int i_problemSize, double* o_results ) {

#ifdef __INTEL_COMPILER
  //! CPU cycles necessary for kernel execution
  size_t l_cycles;
#endif

  //! execution time of the kernel
  double l_executionTime;

  // calculate real problemsize
  int l_problemSize = l_benchmarkData.minimumNumberOfElements + ( (i_problemSize - 1)  * l_benchmarkData.incrementPerStep );

  /* check wether the pointer to store the results in is valid or not */
  if ( o_results == NULL ) return 1;
  
  // start stop watch
  l_stopWatch.start();
  
#ifdef __INTEL_COMPILER
  // Save current cycle counter
  l_cycles = __rdtsc();
#endif
  
  // execute kernel
  cauchyKovalewskiTimeIntegrationKernel( l_problemSize, l_benchmarkData );
  
  
#ifdef __INTEL_COMPILER
    // set clock cycles
    l_cycles = __rdtsc() - l_cycles;
#endif
  
  // get exectution time
  l_executionTime = l_stopWatch.stop();

  o_results[0] = (double)l_problemSize;
  //o_results[1] = dtime;
  o_results[1] = l_executionTime;
  
#ifdef __INTEL_COMPILER
  o_results[2] = l_cycles;
#endif

  return 0;
}

/** Clean up the memory.
 */
void bi_cleanup( void* ) {}

