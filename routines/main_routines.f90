
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

    write(*,*)'-----------'
    do i = 1,n
       write(*,*) X(i,:)
    enddo
    write(*,*)'-----------'

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

    write(*,*)'-----------'
    do i = 1,n
       write(*,*) X(i,:)
    enddo
    write(*,*)'-----------'
    
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
    ! do i = 1,n
    !     write(*,*) A(i,:)
    ! enddo
    call back_substitution(A, Y, X)
 
 end subroutine Cholesky_solver
 
