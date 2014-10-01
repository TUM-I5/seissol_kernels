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
 * Typedefs, structs and macros for the implementation.
 **/

// define order of taylor series relative to the number of basis functions
#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_BASIS_FUNCTIONS 4

#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_BASIS_FUNCTIONS 10

#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_BASIS_FUNCTIONS 20

#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_BASIS_FUNCTIONS 35

#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_BASIS_FUNCTIONS 56

#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_BASIS_FUNCTIONS 84

#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_BASIS_FUNCTIONS 120

#else
#error Preprocessor flag CONVERGENCE_ORDER is not in {2, 3, 4, 5, 6, 7, 8}.
#endif

// aligned number of basis functions
#if ALIGNMENT == 32

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 4
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 12
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 20
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 36
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 84
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#endif

#elif ALIGNMENT == 64

#if CONVERGENCE_ORDER == 2
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 8
#elif CONVERGENCE_ORDER == 3
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 16
#elif CONVERGENCE_ORDER == 4
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 24
#elif CONVERGENCE_ORDER == 5
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 40
#elif CONVERGENCE_ORDER == 6
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 56
#elif CONVERGENCE_ORDER == 7
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 88
#elif CONVERGENCE_ORDER == 8
#define NUMBER_OF_ALIGNED_BASIS_FUNCTIONS 120
#endif

#else

#error ALIGNMENT is not in {32, 64}.

#endif

// number of entries in dense stiffness and flux matrices.
#define NUMBER_OF_DOFS         (NUMBER_OF_BASIS_FUNCTIONS        *NUMBER_OF_QUANTITIES)
#define NUMBER_OF_ALIGNED_DOFS (NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES)

// use double precision for floating point numbers
typedef double real;
