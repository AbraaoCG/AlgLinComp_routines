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
         if (U(i,i) /= 0.0d0) then
            X(i,k) = (Y(i,k) - s) / U(i,i)
         else
            stop "WARNING0 :Matriz singular!"
         endif
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
            s = s + L(i,k) * Y(k,j)
         end do
         if (L(i,i) /= 0.0d0) then  ! Verifica Divisão por 0 ( matriz singular )
            Y(i,j) = (B(i,j) - s) / L(i,i)  ! Calculo do coeficiente na linha e coluna atuais
         else
            stop "WARNING0 :Matriz singular!"
         end if
      end do
   end do

end subroutine forward_substitution



subroutine forward_substitution_inplace(L, B, Y)
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
            s = s + L(i,k) * Y(k,j)
         end do
         if (L(i,i) /= 0.0d0) then  ! Verifica Divisão por 0 ( matriz singular )
            !write(*,*)L(i,i)
            Y(i,j) = (B(i,j) - s) / 1 ! Calculo do coeficiente na linha e coluna atuais. Nessa rotina específica para uso quando é realizada uma decomposição LU inplace,
         else                         ! a divisão por L(i,i) não pode ser feita diretamente, pois em A(i,i) é armazenado os dados de U(i,i), enquanto L(i,i) é sempre 1.
            stop "WARNING0 :Matriz singular!"
         end if
      end do
   end do

end subroutine forward_substitution_inplace