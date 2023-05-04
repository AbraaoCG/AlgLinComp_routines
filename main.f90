module mod1

    CONTAINS
       include 'routines/essentials.f90'
       include 'routines/LU_decomp.f90'
       include 'routines/Cholensky_decomp.f90'
       include 'routines/main_routines.f90'
       
    end module mod1
    
program AXB_Solver
   use mod1
   implicit none
   integer ::n,nb, i, ICOD, io, cb
   character*50 :: filename
   real*8, allocatable :: auxVector(:)
   real*8, allocatable :: A(:,:), X(:,:), B(:,:), Y(:,:), A2(:,:)

   ! Definir número de colunas em cada Vetor B.
   cb = 1
   ICOD = 1

   ! write(*,*) 'Insira nome do arquivo contendo a matrix A: '
   ! read(*,'(A)') filename

   
   filename = 'matrizes/Matriz_A.dat'

   ! Descobrindo número de linhas de A (n).
   n = 0
   OPEN (11,status = 'old', file = filename)
   DO
   READ(11,*,iostat=io)
   IF (io/=0) EXIT
   n = n + 1
   END DO
   rewind(11)

   ! Alocando espaço para A, B e solução X.
  allocate(A(n,n),A2(n,n))
  
  ! Lendo matriz A
  do i = 1,n
   read(11,*)A(i,:)
   A2(i,:) = A(i,:)
  enddo
  CLOSE (11)
  
  write(*,'(A,/,/,A,/,A,/,A,/,A)') 'Insira ICOD:','ICOD = 1 : Solver com decomposição LU ',&
   'ICOD = 2 : Solver com decomposição de Cholesky',&
   'ICOD = 3 : Solver com método iterativo Jacobi',&
   'ICOD = 4 : Solver com método iterativo Gauss-Seidel'
   read(*,*) ICOD

   if(ICOD .eq. 1) then
      call LU_decomposition_inplace(A)
   else if(ICOD .eq. 2) then
      call cholenskyDecomp(A)
   else if(ICOD .eq. 3) then
      print*,'3'
   else if(ICOD .eq. 4) then
      print*,'4'
   endif


   
   
   do while (filename .ne. '-1')
 
      ! write(*,*) 'Insira nome do arquivo contendo B: (Caso deseja fechar o programa, digite -1)'
      ! read(*,'(A)') filename

      filename = 'matrizes/Vetor_B_01.dat'

      if (filename .eq. '-1') then 
         stop
      endif
      
      ! Descobrindo número de linhas de B (nb).
      nb = 0
      OPEN (11, status = 'old', file = filename)

      DO
         READ(11,*,iostat=io)
         IF (io/=0) GO TO 90
         nb = nb + 1
         END DO
      enddo
      90 rewind(11)

      if(nb .eq. n) then
         allocate(B(nb,cb), Y(n,cb), X(n,cb))

         do i = 1,nb
            read(11,*) B(i,:)
            !write(*,*) A(i,:)
         enddo
         
         if ((ICOD .eq. 1 ) ) then
            call forward_substitution_inplace(A, B, Y) ! L . Y = B ;  
            call back_substitution(A, Y, X) ! U . X = Y  
         
         else if ((ICOD .eq. 2 ) ) then
            call forward_substitution(A, B, Y) ! L . Y = B ;  
            call back_substitution(A, Y, X) ! U . X = Y  

         ! else if ((ICOD .eq. 3 ) ) then
         !    call forward_substitution(A, B, Y) ! L . Y = B ;  
         !    call back_substitution(A, Y, X) ! U . X = Y  

         ! else if ((ICOD .eq. 4 ) ) then
         !    call forward_substitution(A, B, Y) ! L . Y = B ;  
         !    call back_substitution(A, Y, X) ! U . X = Y  
         endif
         
         
         
         ! ! Verificar corretude.
         ! do i = 1,n 
         !    write(*,*) B(i,:)
         ! enddo
         ! write(*,*)' ------------------- '
         
         ! CALL matrix_multiplication(A2,X,B)

         ! do i = 1,n 
         !    write(*,*) B(i,:)
         ! enddo


      else
         write(*,*) 'Warning 2 : A matrix and B matrix must have same number of lines.'
      endif

   


   ! ---------------------------- Código de Teste de funcionamento dos algorítimos.

   ! Criando a matriz A manualmente
   ! A = reshape([2.0, 1.0, 1.0, &
   !    4.0, 3.0, 3.0, &
   !    8.0, 7.0, 9.0], [n,n])

   ! A = reshape([4.0, 1.0, 1.0, &             
   ! 1.0, 5.0, 2.0, &            
   !  1.0, 2.0, 6.0], [3,3])


   ! ! Criando a matriz identidade
   ! Im = 0.0
   ! do i = 1, n
   !    Im(i,i) = 1.0
   ! end do

   ! B = 0.0
   ! !call LU_solver(A,B,Im)
   ! call LU_solver(A,B,Im)

   ! Im = -1d0
   ! Call matrix_multiplication(A,B,Im)

   ! write(*,*)'-----------'
   ! do i = 1,n
   !    write(*,*) B(i,:)
   ! enddo
   ! write(*,*)'-----------'

   ! write(*,*)'-----------'
   ! do i = 1,n
   !    write(*,*) Im(i,:)
   ! enddo
   ! write(*,*)'-----------'


    
   ! call Cholesky_solver(A,B,Im)

   ! Im  = -1d0
   ! write(*,*)'-----------'
   ! do i = 1,n
   !    write(*,*) B(i,:)
   ! enddo
   ! write(*,*)'-----------'
   
end program AXB_Solver


