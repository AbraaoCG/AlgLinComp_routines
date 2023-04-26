module mod1

    CONTAINS
       include 'routines/essentials.f90'
       include 'routines/LU_decomp.f90'
       include 'routines/main_routines.f90'
    
    end module mod1
    
program LinAlgebra_func
   use mod1
   implicit none
   integer ::n, i
   real*8, allocatable :: A(:,:), Im(:,:), B(:,:)

   n = 3 ! dimensao da matriz
   allocate(A(n,n), Im(n,n), B(n,n))
   ! Criando a matriz A manualmente
   A = reshape([2.0, 1.0, 1.0, &
      4.0, 3.0, 3.0, &
      8.0, 7.0, 9.0], [n,n])

   ! Criando a matriz identidade
   Im = 0.0
   do i = 1, n
      Im(i,i) = 1.0
   end do

   B = 0.0
   call LU_solver(A,B,Im)
   Im = -1.d0
   Call matrix_multiplication(A,B,Im)
   write(*,*)'-----------'
   do i = 1,n
      write(*,*) Im(i,:)
   enddo
end program LinAlgebra_func


