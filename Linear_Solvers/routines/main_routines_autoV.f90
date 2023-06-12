subroutine power_method(A, eigenvalue, eigenvector)
  implicit none
  integer :: i, n, iter_max, iter_count
  real*8 :: tol, lambda_old, lambda_new, norm_factor

  real*8, intent(in):: A(:,:)
  real*8, intent(out) :: eigenvalue
  real*8, intent(out) :: eigenvector(:)

  n = size(A, 1)
  iter_max = 10000 ! Número máximo de iterações
  tol = 1.0E-5 ! Tolerância para convergência

  eigenvector = 1.0 / sqrt(real(n))! Inicialização do autovetor

  iter_count = 0
  do while (iter_count < iter_max)
    lambda_old = eigenvalue
    eigenvector = matmul(A, eigenvector) ! Multiplica a matriz pelo autovetor
    eigenvalue = maxval(abs(eigenvector)) ! Fator de normalização
    if ( eigenvalue .lt. 1.0e-12) then
      write(*,*) "WARNING: Singular Matrix, can't find eigenvalue! "
      exit
    endif
    eigenvector = eigenvector / eigenvalue ! Normalização do autovetor

    if (abs(eigenvalue - lambda_old) < tol) exit ! Critério de convergência
    iter_count = iter_count + 1
  end do
  print*,'Convergência alcançada na iteração', iter_count

end subroutine power_method


subroutine jacobi_method(A, eigenvalues, eigenvectors)
  implicit none
  integer :: i, j, p, q, n, iter_max, iter_count, maxOffDiag_i, maxOffDiag_j
  real*8,intent(inout) :: A(:,:)
  real*8, intent(out) :: eigenvalues(:), eigenvectors(:,:)
  real*8 :: tol, c, s, t, sum_offdiag, max_offdiag, theta

  n = size(A, 1)
  iter_max = 10000 ! Nú?mero máximo de iterações
  tol = 1.0E-8 ! Tolerância para convergência
  eigenvectors = 0.0
  do i = 1, n
    eigenvectors(i,i) = 1.0
  end do

  iter_count = 0
  do while (iter_count < iter_max) ! limitando numero maximo de iteracoes ( condicao de parada 1 )
    ! Primeiro, deseja-se achar o maior elemento da matriz A e seu indice.
    max_offdiag = 0.0
    do i = 1, n-1
      do j = i+1, n
        if (abs(A(i,j)) > max_offdiag) then
          max_offdiag = abs(A(i,j))
          p = i
          q = j
        end if
      end do
    end do
   
   ! Se o maior elemento da matrix for menor que a tolerancia, considera-se A diagonalizada.
    if (max_offdiag < tol) exit ! Crité?rio de convergê?ncia

   ! Metodo para diagonalizar A.

    ! c é o cosseno de 'fi' e s é o seno de 'fi'. fi é arcotg(2.0 * A(p,q) / (A(p,p) - A(q,q))
    if (abs(A(p,p) - A(q,q)) < tol) then
      c = sqrt(2.0) / 2.0
      s = sqrt(2.0) / 2.0
    else
      theta = 0.5 * atan(2.0 * A(p,q) / (A(p,p) - A(q,q)))
      c = cos(theta)
      s = sin(theta)
    end if

    ! Sendo P a matriz que incializa como identidade e tem: 
    !  P(p,p) = cos(theta) ; P(p,q) = -sen(theta) ; P(q,q) = cos(theta) ; P(q,p) = sen(theta)
    
    ! Equivalente a multiplicar P transposta por A = A'
    do i = 1, n
      t = A(i,p)
      A(i,p) = c * t + s * A(i,q)
      A(i,q) = -s * t + c * A(i,q)
    end do

    ! Equivalente a multiplicar A' por P = A'' = P^t . A . P
    do j = 1, n
      t = A(p,j)
      A(p,j) = c * t + s * A(q,j)
      A(q,j) = -s * t + c * A(q,j)
    end do

    ! Equivalente a multiplicar X por P
    do i = 1, n
      t = eigenvectors(i,p)
      eigenvectors(i,p) = c * t + s * eigenvectors(i,q)
      eigenvectors(i,q) = -s * t + c * eigenvectors(i,q)
    end do

    iter_count = iter_count + 1
  end do
  print*,'Convergência alcançada na iteração', iter_count

  ! Obtenho autovalores da diagonal de A, que agora foi diagonalizada.
  do i = 1, n
   eigenvalues(i) = A(i,i)
  enddo

  

end subroutine
