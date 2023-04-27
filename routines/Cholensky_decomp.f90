subroutine cholenskyDecomp(A)
    implicit none
    integer :: i, j, k, n, info
    real*8 :: sum
    real*8, intent(inout) :: A(:,:)
  
    n = size(A, 1)
    info = 0
  
    do j = 1, n
      sum = 0.0
      do k = 1, j - 1
        sum = sum + A(j,k)**2
      end do
  
      A(j,j) = sqrt(A(j,j) - sum)
  
      if (abs(A(j,j)) < 1e-12) then
        info = -1
        stop 'WARNING1: A não é uma matriz positiva definida'
      endif
  
      do i = j+1, n
        sum = 0.0
        do k = 1, j - 1
          sum = sum + A(i,k)*A(j,k)
        end do
  
        A(i,j) = (A(i,j) - sum) / A(j,j)
      end do
    end do
  
    do i = 1, n
      do j = i+1, n
        A(i,j) = 0.0
      end do
    end do
  
    if (info == 0) then
      write(*,*) 'Fatorização de Cholesky bem sucedida.'
    endif
  end subroutine cholenskyDecomp