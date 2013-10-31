/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Heinecke (heinecke AT in.tum.de)
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
 * Standalone-performance benchmark for the Cauchy Kovalewski / ADER kernel
 **/
 
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sstream>

#include "cauchyKovalewskiTimeIntegrationKernel.hpp"
#include "seissol_src/Monitoring/SeisSolStopwatch.hpp"


//! used benchmark data
BenchmarkData l_benchmarkData;
//! SeisSol stopwatch
SeisSolStopwatch l_stopWatch;

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
  		std::cout << "wrong usage! Use this benchmark as followed:" << std::endl << "  ./driver.exe #elems_start #elems_stop inc" << std::endl;
  		return -1; 
  }

  l_benchmarkData.minimumNumberOfElements = atoi( argv[1] );
  l_benchmarkData.maximumNumberOfElements = atoi( argv[2] );
  l_benchmarkData.incrementPerStep = atoi( argv[3] );
 
  // compute number of steps
  l_benchmarkData.numberOfSteps = (int) ( l_benchmarkData.maximumNumberOfElements - l_benchmarkData.minimumNumberOfElements ) / l_benchmarkData.incrementPerStep;
    
  // allocate memory for the given elements
  l_benchmarkData.initializeDiscontinousGalerkinElements( l_benchmarkData.maximumNumberOfElements ); 
  
#if 0
  // print information about the benchmark
  std::cout << std::endl
            << "ADER kernel:" << std::endl
            << "  maximum(#elements)  : " <<  l_benchmarkData.maximumNumberOfElements << std::endl
            << "  minimum(#elements)  : " <<  l_benchmarkData.minimumNumberOfElements << std::endl
            << "  increcement per step: " <<  l_benchmarkData.incrementPerStep        << std::endl
            << "  #steps              : " <<  l_benchmarkData.numberOfSteps           << std::endl << std::endl;
#endif
            
#ifdef __INTEL_COMPILER
    std::cout << "problemsize,executiontime,cycles" << std::endl; 
#else
    std::cout << "problemsize,executiontime" << std::endl; 
#endif
            
  for (int i_problemSize = 0; i_problemSize <= l_benchmarkData.numberOfSteps; i_problemSize++)
  {
#ifdef __INTEL_COMPILER
    //! CPU cycles necessary for kernel execution
    size_t l_cycles;
#endif

    //! execution time of the kernel
    double l_executionTime;

    // calculate real problemsize
    int l_problemSize = l_benchmarkData.minimumNumberOfElements + ( (i_problemSize)  * l_benchmarkData.incrementPerStep );
  
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

#ifdef __INTEL_COMPILER
    std::cout << (double)l_problemSize << "," << l_executionTime << "," << l_cycles << std::endl;
#else
    std::cout << (double)l_problemSize << "," << l_executionTime << std::endl; 
#endif
  }

  return 0;
}
