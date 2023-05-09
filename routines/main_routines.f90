
subroutine LU_solver_inplace(A, X, B)
    
    implicit none
    integer :: n, m, i
    real*8, intent(inout) :: A(:,:), B(:,:)
    real*8, intent(out) :: X(:,:)
    real*8 :: Y(size(B,1),size(B,2))
 
    n = size(A, 1)
    m = size(B, 2)
 
    call LU_decomposition_inplace(A)
    call forward_substitution_inplace(A, B, Y)
    call back_substitution(A, Y, X)

    ! write(*,*)'-----------'
    ! do i = 1,n
    !    write(*,*) X(i,:)
    ! enddo
    ! write(*,*)'-----------'

 end subroutine LU_solver_inplace
 
 subroutine LU_solver(A, X, B)
    
    implicit none
    integer :: n, m, i
    real*8, intent(in) :: A(:,:), B(:,:)
    real*8, intent(out) :: X(:,:)
    real*8 :: L(size(A,1),size(A,2)), U(size(A,1),size(A,2)), Y(size(B,1),size(B,2))
 
    n = size(A, 1)
    m = size(B, 2)

    call LU_decomposition(A, L, U)
    call forward_substitution(L, B, Y)
    call back_substitution(U, Y, X)

 end subroutine LU_solver
 

 
subroutine Cholesky_solver(A, X, B)
    
    implicit none
    integer :: n, m, i
    real*8, intent(inout) :: A(:,:)
    real*8, intent(in) :: B(:,:)
    real*8, intent(out) :: X(:,:)
    real*8 :: Y(size(B,1),size(B,2))
 
    n = size(A, 1)
    m = size(B, 2)
 
    call cholenskyDecomp(A)
    call forward_substitution(A, B, Y)
    call back_substitution(A, Y, X)
 
 end subroutine Cholesky_solver
 

! subroutine Jacobi_solver(A, X, B)
!    implicit none 
!    integer :: n, i, j
!    real*8, intent(in)
!    ! É preciso verificar se a matriz A é adequada para o método de Jacobi
!    ! Para cada elemento da diagonal principal, a soma dos outros elementos em linha e em coluna são menores que o elemento da diagonal principal.
!    do i = 1, n
!    enddo
! end subroutine Jacobi_solver


subroutine jacobi_method(A, X, B )
   implicit none
   
   integer :: i,j, iter_max, iter, n,cA, p, flagDiagonalDominant
   real*8 :: tol,sumTmp,sumTmp2
   real*8,dimension(:,:), intent(in):: A, B
   real*8,dimension(:,:), intent(out) :: X
   real*8, dimension(size(X,1),size(X,2)) ::  Xold

   ! Valores iniciais
   n = size(A,1)
   cA = size(A,2)
   p = size(X,2)
   X = 1.0 ! Vetor solução inicial.
   iter_max = 1000
   tol = 1.0e-5
   flagDiagonalDominant = 1
   write(*,*)'--------------'
   ! Verifica se a matriz A é diagonalmente dominante e emite um warning caso não seja.
   do i = 1, n
       if (abs(A(i,i)) <= sum(abs(A(i,:))) - abs(A(i,i))) then
         flagDiagonalDominant = 0
       endif
   end do
   if ( flagDiagonalDominant .eq. 0) then
      write(*,*) "WARNING: Matrix isn't diagonally dominant, convergency not guaranteed. "
   endif
   
   ! Método iterativo de Jacobi
   do iter = 1, iter_max     
       ! Copia X para Xold 
       Xold = X(:,:)
       ! Atualiza o valor de cada elemento de X
       do i = 1, n
            ! Computar somatório fora da diagonal
            sumTmp = 0.0
            do j= 1,cA
               if (i .ne. j)then
                  sumTmp = sumTmp + A(i,j) * Xold(j,1)
               endif
            enddo
            ! Verificar se A é singular.
            if(A(i,i) <= 1.0e-15) then
               write(*,*) 'WARNING: Singular Matrix'
               stop
             endif

            ! Recalcular X(i,1)
           X(i,1) = (B(i,1) - sumTmp) / A(i,i) !  
       end do
      ! if(iter .eq. 80) then
         
      !    stop
      ! endif
       ! Verifica a convergência
      sumTmp = 0
      sumTmp2 = 0
      do i = 1,n
         sumTmp  = sumTmp + (X(i,1) - Xold(i,1))**2
         sumTmp2 = sumTmp2 + X(i,1)**2
      enddo
      sumTmp = sqrt(sumTmp) ; sumTmp2 = sqrt(sumTmp2)
      sumTmp = sumTmp / sumTmp2
       if (sumTmp <= tol) then
           write(*,*) 'Convergência alcançada na iteração', iter
           exit
       endif
   end do
   
   ! ! Imprime a solução
   ! write(*,*) 'Solução:'
   ! do i = 1, n
   !     write(*,*) X(i,1)
   ! end do
   
end subroutine jacobi_method


subroutine Gauss_Seidel_method(A, X, B )
   implicit none
   
   integer :: i,j, iter_max, iter, n,cA, p, flagDiagonalDominant
   real*8 :: tol,sumTmp,sumTmp2
   real*8,dimension(:,:), intent(in):: A, B
   real*8,dimension(:,:), intent(out) :: X
   real*8, dimension(size(X,1),size(X,2)) ::  Xold

   ! Valores iniciais
   n = size(A,1)
   cA = size(A,2)
   p = size(X,2)
   X = 1.0 ! Vetor solução inicial.
   iter_max = 1000
   tol = 1.0e-5
   flagDiagonalDominant = 1
   write(*,*)'--------------'
   ! Verifica se a matriz A é diagonalmente dominante e emite um warning caso não seja.
   do i = 1, n
       if (abs(A(i,i)) <= sum(abs(A(i,:))) - abs(A(i,i))) then
         flagDiagonalDominant = 0
       endif
   end do
   if ( flagDiagonalDominant .eq. 0) then
      write(*,*) "WARNING: Matrix isn't diagonally dominant, convergency not guaranteed. "
   endif
   
   ! Método iterativo de Jacobi
   do iter = 1, iter_max     
       ! Copia X para Xold 
       Xold = X(:,:)
       ! Atualiza o valor de cada elemento de X
       do i = 1, n
            ! Computar somatório fora da diagonal
            sumTmp = 0.0
            do j= 1,cA
               ! Para o método de Gauss - Seidel, se j < i então uso valor já atualizado para X(i), senão utilizo valor da última iteração.
               if (j .lt. i)then
                  sumTmp = sumTmp + A(i,j) * X(j,1)
               else if(j .gt. i) then 
                  sumTmp = sumTmp + A(i,j) * Xold(j,1)
               endif
            enddo
            ! Verificar se A é singular.
            if(A(i,i) <= 1.0e-15) then
               write(*,*) 'WARNING: Singular Matrix'
               stop
             endif

            ! Recalcular X(i,1)
           X(i,1) = (B(i,1) - sumTmp) / A(i,i) !  
       end do
      ! if(iter .eq. 80) then
         
      !    stop
      ! endif
       ! Verifica a convergência
      sumTmp = 0
      sumTmp2 = 0
      do i = 1,n
         sumTmp  = sumTmp + (X(i,1) - Xold(i,1))**2
         sumTmp2 = sumTmp2 + X(i,1)**2
      enddo
      sumTmp = sqrt(sumTmp) ; sumTmp2 = sqrt(sumTmp2)
      sumTmp = sumTmp / sumTmp2
       if (sumTmp <= tol) then
           write(*,*) 'Convergência alcançada na iteração', iter
           exit
       endif
   end do
   
   ! ! Imprime a solução
   ! write(*,*) 'Solução:'
   ! do i = 1, n
   !     write(*,*) X(i,1)
   ! end do
   
end subroutine Gauss_Seidel_method