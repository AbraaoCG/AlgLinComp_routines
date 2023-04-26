
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
    do i = 1,n
        write(*,*) U(i,:)
    enddo
    call back_substitution(U, Y, X)
 
 end subroutine LU_solver
 