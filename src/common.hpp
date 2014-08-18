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
 * Common kernel-level functions
 **/

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include "typedefs.hpp"

namespace seissol {
  namespace kernels {
    /**
     * Gets the number of basis functions for the given convergence order.
     *
     * @param i_convergenceOrder convergence order.
     * @return number of basis funcitons.
     **/
    static unsigned int getNumberOfBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER ) {
      return i_convergenceOrder*(i_convergenceOrder+1)*(i_convergenceOrder+2)/6;
    }

    /**
     * Get the # of basis functions aligned to the given boundaries.
     *
     * @param i_convergenceOrder convergence order.
     * @param i_alignment alignment in bits.
     * @return aligned number of basis functions.
     **/
    static unsigned int getNumberOfAlignedBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER,
                                                          unsigned int i_alignment        = ALIGNMENT ) {
      unsigned int l_nonZeroBits = getNumberOfBasisFunctions( i_convergenceOrder) * sizeof(real);

      // compute corresponding aligned # of basis functions
      unsigned int l_alignedBits = l_nonZeroBits + (i_alignment-(l_nonZeroBits % i_alignment))%i_alignment;
      assert( l_alignedBits % sizeof(real) == 0 );

      return l_alignedBits / sizeof( real );
    }
  }
}

#endif
