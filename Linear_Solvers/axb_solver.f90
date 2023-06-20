module mod1

    CONTAINS
       include 'routines/essentials.f90'
       include 'routines/LU_decomp.f90'
       include 'routines/Cholesky_decomp.f90'
       include 'routines/main_routines_axb.f90'
       
end module mod1
    
program AXB_Solver
   use mod1
   implicit none
   integer ::n, nb, i, ICOD, io, cb
   character*50 :: filename
   real*8, allocatable :: auxVector(:)
   real*8, allocatable :: A(:,:), X(:,:), B(:,:), Y(:,:)

   ! Definir número de colunas em cada Vetor B.
   cb = 3

   ! Ler caminho do arquivo que contém a matrix A
   write(*,*) 'Insira o caminho ate o arquivo contendo a matrix A: '
   read(*,'(A)') filename

   ! Descobrir número de linhas / colunas de A. ( Supõe-se que seja quadrada).
   n = 0
   OPEN (11,status = 'old', file = filename)
   DO
   READ(11,*,iostat=io)
   IF (io/=0) EXIT
   n = n + 1
   END DO
   rewind(11)

   ! Alocando espaço para A, B e solução X.
  allocate(A(n,n))
  
  ! Lendo matriz A
  do i = 1,n
   read(11,*)A(i,:)
   write(*,*) A(i,:)
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
   endif
   
   do while (filename .ne. '-1')
      write(*,*) 'Insira o caminho ate o arquivo contendo B: (Caso deseja fechar o programa, digite -1)'
      read(*,'(A)') filename
      if (filename .eq. '-1') then 
         stop
      endif
      ! Descobrindo número de linhas de B (nb).
      OPEN (12, status = 'old', file = filename)
      nb = 0
      DO
         READ(12,*,iostat=io)
         IF (io/=0) GO TO 90
         nb = nb + 1
      enddo
      90 rewind(12) 
      if(nb .eq. n) then ! Verifico se número de linhas de B é equivalente ao de A.
         allocate(B(nb,cb), Y(n,cb), X(n,cb)) ! Aloco espaço para matrizes B, Y( intermediária ) e X ( Solução )
         do i = 1,nb
            read(12,*) B(i,:) ! Leitura de B
         enddo
         ! Dependendo do método escolhido, realiza substituição para frente e para trás ou realiza método iterativo.
         if ((ICOD .eq. 1 ) ) then
            call forward_substitution_inplace(A, B, Y) ! L . Y = B ;  
            call back_substitution(A, Y, X) ! U . X = Y  

         else if ((ICOD .eq. 2 ) ) then
            call forward_substitution(A, B, Y) ! L . Y = B ;  
            call back_substitution(A, Y, X) ! U . X = Y  

         else if ((ICOD .eq. 3 ) ) then
            call jacobi_method(A, X, B)

         else if ((ICOD .eq. 4 ) ) then
            call Gauss_Seidel_method(A, X, B)
         endif
         ! Imprimir X
         write(*,*) 'Matrix solucao X: '
         write(*,*) n, cb
         do i = 1,n
            write(*,*) X(i,:)
         enddo

         deallocate(B,Y,X)
      endif
   enddo
      
end program AXB_Solver
