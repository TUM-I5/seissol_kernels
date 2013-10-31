! @file
! This file is part of SeisSol.
!
! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!
! @section LICENSE
! This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
!
! According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
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
! Binds SeisSol's Cauchy Kovalewski procedures to C.

module c_bind_cauchy_kovalewski
  implicit none

  contains

  !
  ! Computes the Cauchy Kovalewski time integration.
  !   This is done by calling the SeisSol's "old" CauchyKovalewski3D routine.
  !
  ! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
  !
  subroutine c_bind_cauchy_kovalewski_3d( i_unknowns,                                                                     &
                                          i_numberOfBasisFunctions, i_numberOfVariables, i_orderOfTaylorSeriesExpansion, &
                                          i_aStarDivM, i_bStarDivM,  i_cStarDivM,                                         &
                                          i_kXi, i_kEta, i_kZeta,                                                         &
                                          i_deltaT,                                                                       &
                                          o_timeIntegratedUnknowns ) bind(c, name='c_bind_cauchy_kovalewski_3d')
    ! SeisSol's "old" CauchyKovalewski3D kernel
    use CauchyKovalewski_mod

    ! Sparse type
    use typesDef

    ! Sparse operators
    use common_operators_mod

    implicit none
    !------------------------------------------------------------------------------------------------------------------------------!
    ! function parameters
    real*8, intent(in), dimension(i_numberOfBasisFunctions, i_numberOfVariables)      :: i_unknowns                                  !< matrix $Q_k(t^n)$, which contains the unkowns at time $t^n$ for a particular element
    integer, intent(in)                                                               :: i_numberOfBasisFunctions, &                 !< number of basis functions (implied by the order of the used method)
                                                                                         i_numberOfVariables, &                      !< number of quantities (depends on the equations, which are solved, i.e. 3D elastics -> 9)
                                                                                         i_orderOfTaylorSeriesExpansion              !< order of the taylor series expansion which implies the number of time derivatives replaced by space derivaties:
                                                                                                                                     ! $\sum_{j=0}^J \frac{(t^{n+1}-t^n)^{j+1}}{(j+1)!} \frac{\partial^{j}}{\partial t^j} Q_k(t^n) =: I(t^{n+1} - t^n)Q_k(t^n)$

    real*8, intent(in), dimension(i_numberOfVariables,     i_numberOfVariables)       :: i_aStarDivM, i_bStarDivM, i_cStarDivM       !< jacobians $A^*_k$, $B^*_k$ and $C*_k$ (already divided by the mass matrix $M$)
    real*8, intent(in), dimension(i_numberOfBasisFunctions, i_numberOfBasisFunctions) :: i_kXi, i_kEta, i_kZeta                      ! stiffness matrices $K^\xi$, $K^eta$ and $K^zeta$
    real*8, intent(in)                                                                :: i_deltaT                                    !< width of time interval used in the integration: $t^{n+1} = t^n + \Delta t$

    real*8, intent(out), dimension(i_numberOfBasisFunctions, i_numberOfVariables)     :: o_timeIntegratedUnknowns                    !< time integrated unknowns: $I(t^{n+1}-t^n}Q_k(t^n) = \int_{t^n}^{t^n+1} Q_k(t) \mathrm{d}t$

    ! local variables
    integer                                                                           :: l_order                                     !< loop variables
    real*8, dimension(i_numberOfBasisFunctions,i_numberOfBasisFunctions,1)            :: l_temporaryStiffnessTensor                  !< temporary tensor for the stiffness matrices (expected by initializing routines)
    real*8, dimension(i_numberOfVariables,i_numberOfVariables,1)                      :: l_temporaryStarTensor                       !< temporary tensor for the star matrices (expected by initializing routines)
    real*8                                                                            :: l_taylorSeriesFactor                        !< factor before every term of the sum in the taylor series expansion: $\frac{(t^{n+1} - t^n)^{j+1}}{(j+1)!}$
    !
    !
    type(tSparseTensor3)                                                             :: l_aStarDivMSp, l_bStarDivMSp, l_cStarDivMSp; ! sparse representations of the jacobians
    type(tSparseTensor3b)                                                            :: l_kXiSp, l_kEtaSp, l_kZetaSp;                ! sparse representations of the jacobians

    type(tSparseTensor3)                                                             :: l_reactionMatrixSp                           ! matrix for attenuation (required by the routine even for elastics, but not used!)
    type(tSparseTensor3b)                                                            :: l_mklmSp                                     ! matrix for attenuation (required by the routine even for elastics, but not used!)
    !
    !------------------------------------------------------------------------------------------------------------------------------!
#if 0
    write(*,*) 'Welcome to Fortran: c_bind_cauchy_kovalewski_3d called'
    write(*,*) '#basisFunctions, #variables: ', i_numberOfBasisFunctions, i_numberOfVariables;
    write(*,*) 'delta T:', i_deltaT
    write(*,*) 'printing unknowns'
    call printDenseMatrix( i_unknowns );
    write(*,*) 'printing time int. unknowns';
    call printDenseMatrix( o_timeIntegratedUnknowns );
    write(*,*) 'printing K^\xi'
    call printDenseMatrix( i_kXi );
    write(*,*) 'printing K^\eta';
    call printDenseMatrix( i_kEta );
    write(*,*) 'printing K^\Zeta';
    call printDenseMatrix( i_kZeta );
    write(*,*) 'printing A_k^*'
    call printDenseMatrix( i_aStarDivM );
    write(*,*) 'printing B_k^*'
    call printDenseMatrix( i_bStarDivM );
    write(*,*) 'printing C_k^*'
    call printDenseMatrix( i_cStarDivM );
#endif
    ! copy stiffness matrices into the temporary tensor and initialize corr. sparse tensors
    l_temporaryStiffnessTensor(:,:,1) = i_kXi;
    call IniSparseTensor3b(l_kXiSp, l_temporaryStiffnessTensor, i_numberOfBasisFunctions, i_numberOfBasisFunctions, 1)
    l_temporaryStiffnessTensor(:,:,1) = i_kEta;
    call IniSparseTensor3b(l_kEtaSp, l_temporaryStiffnessTensor, i_numberOfBasisFunctions, i_numberOfBasisFunctions, 1)
    l_temporaryStiffnessTensor(:,:,1) = i_kZeta;
    call IniSparseTensor3b(l_kZetaSp, l_temporaryStiffnessTensor, i_numberOfBasisFunctions, i_numberOfBasisFunctions, 1)

    ! copy star matrices into the temporary tensor and initialize corr. sparse tensors
    l_temporaryStarTensor(:,:,1) = Transpose(i_aStarDivM);
    call IniSparseTensor3(l_aStarDivMSp, l_temporaryStarTensor, i_numberOfVariables, i_numberOfVariables, 1)
    l_temporaryStarTensor(:,:,1) = Transpose(i_bStarDivM);
    call IniSparseTensor3(l_bStarDivMSp, l_temporaryStarTensor, i_numberOfVariables, i_numberOfVariables, 1)
    l_temporaryStarTensor(:,:,1) = Transpose(i_cStarDivM);
    call IniSparseTensor3(l_cStarDivMSp, l_temporaryStarTensor, i_numberOfVariables, i_numberOfVariables, 1)

    ! initialize required reaction matrices (not used)
    l_temporaryStarTensor = 0;
    call IniSparseTensor3(l_reactionMatrixSp, l_temporaryStarTensor, i_numberOfVariables, i_numberOfVariables, 1)
    call IniSparseTensor3b(l_mklmSp,           l_temporaryStarTensor, i_numberOfVariables, i_numberOfVariables, 1)

    ! call SeisSol's "old" CauchyKovalewski3D kernel
    call cauchyKovalewski3D( TimeIntDof = o_timeIntegratedUnknowns, Dof = i_unknowns, Dt = i_deltaT, &
                             A_Sp = l_aStarDivMSp, B_Sp = l_bStarDivMSp, C_Sp = l_cStarDivMSp, E_Sp = l_reactionMatrixSp, &
                             ADGxi_Sp = l_kXiSp, ADGeta_Sp = l_kEtaSp, ADGzeta_Sp = l_kZetaSp, Mklm_Sp = l_mklmSp, &
                             ReactionTerm = 0, LocDegFr = i_numberOfBasisFunctions, LocDegFrMat = 1, &
                             LocPoly = i_orderOfTaylorSeriesExpansion, nVar = i_numberOfVariables )
    
    ! free memory
    call closeSparseTensor3b( l_kXiSp   )
    call closeSparseTensor3b( l_kEtaSp  )
    call closeSparseTensor3b( l_kZetaSp )
    call closeSparseTensor3( l_aStarDivMSp )
    call closeSparseTensor3( l_bStarDivMSp )
    call closeSparseTensor3( l_cStarDivMSp )
    
    call closeSparseTensor3( l_reactionMatrixSp )
    call closeSparseTensor3b( l_mklmSp )   

#if 0
    write(*,*) 'printing time int. unknowns';
    call printDenseMatrix( o_timeIntegratedUnknowns );
    write(*,*) 'c_bind_cauchy_kovalewski_3d finished.'
#endif
  end subroutine c_bind_cauchy_kovalewski_3d

  !
  ! Prints a given dense matrix.
  !
  ! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
  !
  subroutine printDenseMatrix( i_denseMatrix )
    ! variable declaration
    real*8, dimension(:,:) :: i_denseMatrix;
    integer              :: l_numberOfRows, l_numberOfColumns;

    ! loop variables
    integer :: l_i,l_j;

    l_numberOfRows    = ubound(i_denseMatrix, 1);
    l_numberOfColumns = ubound(i_denseMatrix,2)

#ifdef PRINTMAPLE
    write(*,'(A)',advance='no') 'Matrix ([';
    do l_i = 1,l_numberOfRows
        write(*,'(A)',advance='no') '['
        do l_j = 1,l_numberOfColumns
            write(*,'(1e20.7)', advance='no') i_denseMatrix(l_i,l_j);
            if(l_j .ne. l_numberOfColumns) write(*,'(A)', advance='no') ','
        end do
        write(*,'(A)',advance='no') '],'
    end do
    write(*,'(A)',advance='no') '])';
    write(*,*);
#endif
#ifdef PRINTPATTERN
    do l_i = 1,l_numberOfRows
        do l_j = 1,l_numberOfColumns
            if( abs(l_denseMatrix(l_i,l_j)) < 1e-10 ) then
                write(*,'(A)', advance='no') '-';
            else
                write(*,'(A)', advance='no') '*';
            endif
      end do
        write(*,*);
    end do
#else
    do l_i = 1,l_numberOfRows
      do l_j = 1,l_numberOfColumns
        write(*,'(1e20.7)', advance='no') i_denseMatrix(l_i,l_j);
      end do
        write(*,*);
    end do
#endif
end subroutine printDenseMatrix
end module
