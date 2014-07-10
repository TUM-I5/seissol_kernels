/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Counts the floating point operations in SeisSol.
 **/
#ifndef NDEBUG

#include "FlopCounter.hpp"

  // Define the FLOP counter.
  unsigned long long num_flops = 0;

  unsigned long long l_numberOfTimeIntegrationFlops = 0;
  unsigned long long l_numberOfVolumeIntegrationFlops = 0;
  unsigned long long l_numberOfBoundaryIntegrationFlops = 0;

  // prevent name mangling
  extern "C" {
    /**
     * Resets the FLOP counter to zero.
     */
    void resetFlops() {
      num_flops = 0;
    }

    /**
     * Stores the boundary FLOPs and resets the counter to zero.
     */
    void addTimeFlops() {
      l_numberOfTimeIntegrationFlops += num_flops;
      resetFlops();
    }

    /**
     * Stores the volume FLOPs and resets the counter to zero.
     */
    void addVolumeFlops() {
      l_numberOfVolumeIntegrationFlops += num_flops;
      resetFlops();
    }

    /**
     * Stores the boundary FLOPs and resets the counter to zero.
     */
    void addBoundaryFlops() {
      l_numberOfBoundaryIntegrationFlops += num_flops;
      resetFlops();
    }
    
    /**
     * Prints the measured FLOPS.
     */
    void printFlops() {
      logInfo() << "FLOPS - "
                << "time: "    << l_numberOfTimeIntegrationFlops   << ", "
                << "vol: "     << l_numberOfVolumeIntegrationFlops << ", "
                << "bnd: "     << l_numberOfBoundaryIntegrationFlops;
    }
  }
#endif
