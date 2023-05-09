
subroutine LU_decomposition(A, L, U)
   implicit none
   integer :: i, j, k, n
   real*8 :: s
   real*8, intent(in) :: A(:,:)
   real*8, intent(out) :: L(:,:), U(:,:)

   n = size(A, 1)

   L = 0.0
   U = 0.0

   do j = 1, n
      L(j,j) = 1.0
      do i = 1, j
         s = 0
         do k = 1, i - 1
            s = s+ U(k,j) * L(i,k)
         enddo
         U(i,j) = A(i,j) - s
      end do

      do i = j, n
         s = 0
         do  k = 1, j - 1
            s = s+ U(k,j) * L(i,k)
         enddo
         if (U(j,j) /= 0.0d0) then! Verifica Divisão por 0 ( matriz singular )
            L(i,j) = (A(i,j) - s) / U(j,j)
         endif
      end do
   end do

end subroutine LU_decomposition



subroutine LU_decomposition_inplace(A)
   implicit none
   integer :: i, j, k, n, info
   real*8 :: s
   real*8, intent(inout) :: A(:,:)

   n = size(A, 1)

   do j = 1, n ! Iteração em linhas

      do i = 1, j ! Iteração em colunas até a coluna j. (Atualizar Matrix U dentro de A)
         s = 0
         do k = 1, i - 1
            s = s + A(i,k) * A(k,j)
         enddo
         A(i,j) = A(i,j) - s
      end do

      !Se elemento da diagonal principal for 0, então a matriz é singular.
      if (A(j,j) == 0.0d0) then 
         info = j
         return
      endif
      
      ! Iteração em colunas de j até n. (Atualizar Matrix L dentro de A)
      do i = j + 1, n
         s = 0
         do k = 1, j - 1
            s = s + A(i,k) * A(k,j)
         enddo
         A(i,j) = (A(i,j) - s) / A(j,j)
      end do
   end do

   info = 0
end subroutine LU_decomposition_inplace

