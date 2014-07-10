#if 0
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2012-2013, SeisSol Group
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
 * Defines preprocessor variables and macros used in SeisSol.
 */
#endif

#if 0
!
! Fortran Logging
!
#endif

#ifndef __STDC__

#ifndef LOGLEVEL0
#define LOGLEVEL0 LOGLEVEL
#endif

#define FORTRAN_STDERR 0
#define FORTRAN_STDOUT 6

#define FORTRAN_LINE_SIZE 1500

#define logDebug(format)   if (LOGLEVEL .ge. 3) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Debug   |'; if (LOGLEVEL .ge. 3) write(FORTRAN_STDERR, format)
#define logInfo(format)    if (LOGLEVEL .ge. 2) write(FORTRAN_STDOUT, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Info    |'; if (LOGLEVEL .ge. 2) write(FORTRAN_STDOUT, format)
#define logWarning(format) if (LOGLEVEL .ge. 1) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Warning |'; if (LOGLEVEL .ge. 1) write(FORTRAN_STDERR, format)
#define logError(format)   if (LOGLEVEL .ge. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', myrank, ' | Error   |'; if (LOGLEVEL .ge. 0) write(FORTRAN_STDERR, format)

#define logDebug0(format)   if (LOGLEVEL0 .ge. 3 .and. myrank .eq. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Debug   |'; if (LOGLEVEL0 .ge. 3 .and. myrank .eq. 0) write(FORTRAN_STDERR, format)
#define logInfo0(format)    if (LOGLEVEL0 .ge. 2 .and. myrank .eq. 0) write(FORTRAN_STDOUT, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Info    |'; if (LOGLEVEL0 .ge. 2 .and. myrank .eq. 0) write(FORTRAN_STDOUT, format)
#define logWarning0(format) if (LOGLEVEL0 .ge. 1 .and. myrank .eq. 0) write(FORTRAN_STDERR, '(A, I8, A)', advance='no') 'Rank: ', 0, ' | Warning |'; if (LOGLEVEL0 .ge. 1 .and. myrank .eq. 0) write(FORTRAN_STDERR, format)

#endif


#if 0
! Manual instrumentation for Scalasca with epik.
! Override function calls if not compiled with EPIK.
#endif

#if defined(EPIK) && !defined(__STDC__)
#include "epik_user.inc"
#else
#define EPIK_FUNC_REG(str)
#define EPIK_FUNC_START()
#define EPIK_FUNC_END()
#define EPIK_USER_REG(id,str)
#define EPIK_USER_START(id)
#define EPIK_USER_END(id)
#endif


#if 0
!
! Generated Kernels
!
#endif

#ifdef GENERATEDKERNELS

#if 0
! preprocessor concatenation
#endif
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define CONCAT_HELPER_4(a,b,c,d) a ## b ## c ## d
#define CONCAT_4(a,b,c,d) CONCAT_HELPER_4(a,b,c,d)
#define CONCAT_HELPER_6(a,b,c,d,e,f) a ## b ## c ## d ## e ## f
#define CONCAT_6(a,b,c,d,e,f) CONCAT_HELPER_6(a,b,c,d, e, f)

#ifndef __STDC__
#include <seissol/generated_code/pre_processor/matrixSizes3D.fpp>
#else
#include <seissol/generated_code/pre_processor/matrixSizes3D.h>
#endif

#if 0
! zero tolerance constant
#endif
#ifndef ZEROTOLERANCE
#define ZEROTOLERANCE 1e-10
#endif

#if 0
! check if the necessary precompiler macros are defined
#endif
#ifndef NUMBEROFVARIABLES
#error Preprocessor flag NUMBEROFVARIABLES not set.
#endif

#ifndef NUMBEROFBASISFUNCTIONS
#error  Preprocessor flag NUMBEROFBASISFUNCTIONS not set.
#endif


#if 0
! check valid preprocessor flags.
! * The number of variables depends on the degree of attenuation.
!   Whereas 9 variables correspond to the elastic wave equations.
! * The number of basis functions depends on the order of the discontinuous Galerkin method.
!   Therfore valid values are given by the formula (O)*(O+1)*(O+2)/6, where O is the order of method.
#endif
#if NUMBEROFVARIABLES != 9
#error Preprocessor flag NUMBEROFVARIABLES is not equal 9 (elastic wave equations).
#endif

#if 0
! define flat array sizes
#endif
#define NUMBEROFUNKNOWNS (NUMBEROFBASISFUNCTIONS*NUMBEROFVARIABLES)

#if 0
// define order of taylor series relative to the number of basis functions
#endif
#if NUMBEROFBASISFUNCTIONS == 1
#define ORDEROFTAYLORSERIESEXPANSION 1

#elif NUMBEROFBASISFUNCTIONS == 4
#define ORDEROFTAYLORSERIESEXPANSION 2

#elif NUMBEROFBASISFUNCTIONS == 10
#define ORDEROFTAYLORSERIESEXPANSION 3

#elif NUMBEROFBASISFUNCTIONS == 20
#define ORDEROFTAYLORSERIESEXPANSION 4

#elif NUMBEROFBASISFUNCTIONS == 35
#define ORDEROFTAYLORSERIESEXPANSION 5

#elif NUMBEROFBASISFUNCTIONS == 56
#define ORDEROFTAYLORSERIESEXPANSION 6

#elif NUMBEROFBASISFUNCTIONS == 84
#define ORDEROFTAYLORSERIESEXPANSION 7

#else
#error Preprocessor flag NUMBEROFBASISFUNCTIONS is not in {1, 4, 10, 20, 35, 56, 84}.
#endif

#if 0
! number of entries in dense stiffness and flux matrices.
#endif
#define NUMBEROFBASISFUNCTIONSSQUARED (NUMBEROFBASISFUNCTIONS*NUMBEROFBASISFUNCTIONS)

#if 0
! number of entries in dense matrices, which compute the numerical flux.
#endif
#define NUMBEROFVARIABLESSQUARED (NUMBEROFVARIABLES*NUMBEROFVARIABLES)


#if 0
! fortran specific variables
#endif
#ifndef __STDC__
! position of the stderr stream
#ifndef STDERR
#define STDERR 0
#endif

#if (NUMBEROFBASISFUNCTIONS !=  1) .and. (NUMBEROFBASISFUNCTIONS !=   4)  .and.\
    (NUMBEROFBASISFUNCTIONS != 10) .and. (NUMBEROFBASISFUNCTIONS !=  20)  .and.\
    (NUMBEROFBASISFUNCTIONS != 35) .and. (NUMBEROFBASISFUNCTIONS !=  56)  .and.\
    (NUMBEROFBASISFUNCTIONS != 84) .and. (NUMBEROFBASISFUNCTIONS != 120)
#error Preprocessor flag NUMBEROFBASISFUNCTIONS is not in {1, 4, 10, 20, 35, 56, 84, 120}.
#endif

#endif

#endif
