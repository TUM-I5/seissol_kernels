! @file
! This file is part of SeisSol.
! 
! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
! @author Alexander Heinecke (heinecke AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Alexander_Heinecke,_M.Sc.,_M.Sc._with_honors)
! 
! @date 2012-11-29 18:54:13.117382
! 
! @section LICENSE
! This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
! 
! According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting the this software software.
! 
! The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
! 
! Copyright (c) 2012
! Technische Universitaet Muenchen
! Department of Informatics
! Chair of Scientific Computing
! http://www5.in.tum.de/
! 
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 
! Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
! Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
!
! @section DESCRIPTION
! 
! Defines matrix sizes in the preprocessing phase.
!   Remark: This file was generated. It should not be modified by hand.
! 
#if NUMBEROFBASISFUNCTIONS == 3
#define FP11_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP11_NUMBEROFNONZEROS 5
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP11_NUMBEROFNONZEROS 14
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP11_NUMBEROFNONZEROS 30
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP11_NUMBEROFNONZEROS 55
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP12_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP12_NUMBEROFNONZEROS 8
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP12_NUMBEROFNONZEROS 31
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP12_NUMBEROFNONZEROS 85
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP12_NUMBEROFNONZEROS 184
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP13_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP13_NUMBEROFNONZEROS 8
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP13_NUMBEROFNONZEROS 31
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP13_NUMBEROFNONZEROS 85
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP13_NUMBEROFNONZEROS 184
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP21_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP21_NUMBEROFNONZEROS 8
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP21_NUMBEROFNONZEROS 31
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP21_NUMBEROFNONZEROS 85
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP21_NUMBEROFNONZEROS 184
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP22_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP22_NUMBEROFNONZEROS 9
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP22_NUMBEROFNONZEROS 36
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP22_NUMBEROFNONZEROS 100
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP22_NUMBEROFNONZEROS 223
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP23_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP23_NUMBEROFNONZEROS 7
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP23_NUMBEROFNONZEROS 26
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP23_NUMBEROFNONZEROS 70
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP23_NUMBEROFNONZEROS 155
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP31_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP31_NUMBEROFNONZEROS 8
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP31_NUMBEROFNONZEROS 31
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP31_NUMBEROFNONZEROS 85
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP31_NUMBEROFNONZEROS 184
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP32_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP32_NUMBEROFNONZEROS 7
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP32_NUMBEROFNONZEROS 26
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP32_NUMBEROFNONZEROS 70
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP32_NUMBEROFNONZEROS 155
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define FP33_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define FP33_NUMBEROFNONZEROS 9
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define FP33_NUMBEROFNONZEROS 36
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define FP33_NUMBEROFNONZEROS 100
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define FP33_NUMBEROFNONZEROS 223
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define KETAT_NUMBEROFNONZEROS 0
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define KETAT_NUMBEROFNONZEROS 2
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define KETAT_NUMBEROFNONZEROS 10
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define KETAT_NUMBEROFNONZEROS 30
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define KETAT_NUMBEROFNONZEROS 70
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define KETA_NUMBEROFNONZEROS 0
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define KETA_NUMBEROFNONZEROS 2
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define KETA_NUMBEROFNONZEROS 10
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define KETA_NUMBEROFNONZEROS 30
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define KETA_NUMBEROFNONZEROS 70
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define KXIT_NUMBEROFNONZEROS 0
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define KXIT_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define KXIT_NUMBEROFNONZEROS 4
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define KXIT_NUMBEROFNONZEROS 13
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define KXIT_NUMBEROFNONZEROS 30
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define KXI_NUMBEROFNONZEROS 0
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define KXI_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define KXI_NUMBEROFNONZEROS 4
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define KXI_NUMBEROFNONZEROS 13
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define KXI_NUMBEROFNONZEROS 30
#endif 

#if NUMBEROFBASISFUNCTIONS == 3
#define M_NUMBEROFNONZEROS 1
#endif 

#if NUMBEROFBASISFUNCTIONS == 6
#define M_NUMBEROFNONZEROS 3
#endif 

#if NUMBEROFBASISFUNCTIONS == 10
#define M_NUMBEROFNONZEROS 6
#endif 

#if NUMBEROFBASISFUNCTIONS == 15
#define M_NUMBEROFNONZEROS 10
#endif 

#if NUMBEROFBASISFUNCTIONS == 21
#define M_NUMBEROFNONZEROS 15
#endif 

