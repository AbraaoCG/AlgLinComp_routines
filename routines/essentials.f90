subroutine matrix_multiplication(A, B, C)
   implicit none
   real*8, dimension(:,:), intent(in) :: A, B
   real*8, dimension(:,:), intent(out) :: C
   integer :: m, n, p, i, j, k

   m = size(A, 1)
   n = size(A, 2)
   p = size(B, 2)

   if (size(B, 1) /= n) then
      print *, "Matrizes não compatíveis para multiplicação"
      stop
   end if

   do i = 1, m
      do j = 1, p
         C(i,j) = 0.0
         do k = 1, n
            C(i,j) = C(i,j) + A(i,k) * B(k,j)
         end do
      end do
   end do

end subroutine matrix_multiplication



subroutine transpose(X, Xt)
   real*8, dimension(:,:) :: X
   real*8, dimension(:,:) :: Xt
   integer :: i, j
   integer :: m, n
   m = size(X, 1)
   n = size(X, 2)
   do i = 1, n
      do j = 1, m
         Xt(i, j) = X(j, i)
      end do
   end do
end subroutine transpose
