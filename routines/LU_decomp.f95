
subroutine back_substitution(U, Y, X)
   implicit none
   real*8, dimension(:,:), intent(in) :: U, Y
   real*8, dimension(:,:), intent(out) :: X
   integer :: m, n, yl, yc, i, j, k
   real*8 :: s

   m = size(U, 1) ! Dimensão da matriz U
   n = size(U, 2) ! Número de colunas de Y
   yl = size(Y, 1)
   yc = size(Y, 2)

   ! Inicializa matriz solução X com zeros
   X = 0.0
   do i = 1, yl
      do j = 1, yc
         X(i,j) = 0.0
      end do
   end do

   do k = 1, yc ! Loop pelas colunas de Y
      do i = m, 1, -1 ! Loop pelas linhas de U e X (de trás para frente)
         s = 0.0
         do j = i+1, m
            s = s + U(i,j) * X(j,k)
         end do
         X(i,k) = (Y(i,k) - s) / U(i,i)
      end do
   end do

end subroutine back_substitution

subroutine forward_substitution(L, B, Y)
   implicit none
   real*8 :: s
   real*8, intent(in) :: L(:,:), B(:,:)
   real*8, intent(out) :: Y(:,:)
   integer :: n, m, i, j, k

   n = size(L,1)  ! Dimensão da matriz L
   m = size(B,2)  ! Número de colunas de B

   ! Inicializa matriz solução Y com zeros
   Y = 0.0d0
   do j = 1, m  ! Loop pelas colunas de B
      do i = 1, n  ! Loop pelas linhas de L e Y
         ! Somatorio dos termos com coeficiente já calculados ( à esquerda do pivo atual i)
         s = 0.0d0
         do k = 1, i-1
            write(*,*) k 
            s = s + L(i,k) * Y(k,j)
         end do
         if (L(i,i) /= 0.0d0) then  ! Verifica Divisão por 0 ( matriz singular )
            Y(i,j) = (B(i,j) - s) / L(i,i)  ! Calculo do coeficiente na linha e coluna atuais
         else
            stop "Matriz singular!"
         end if
      end do
   end do

end subroutine forward_substitution


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
         L(i,j) = (A(i,j) - s) / U(j,j)
      end do
   end do

end subroutine LU_decomposition