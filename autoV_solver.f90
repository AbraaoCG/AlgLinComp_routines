module mod1

    CONTAINS
       include 'routines/essentials.f90'
       include 'routines/main_routines_autoV.f90'
       
end module mod1
    
program eigen_calculator
   use mod1
   implicit none
   integer ::n, i, ICOD, io
   character*50 :: filename
   real*8, allocatable :: A(:,:), X(:,:), E(:) ! X : matriz com autovetores ; E : vetor com autovalores

   !write(*,*) 'Insira o caminho até o arquivo contendo a matrix A: '
   !read(*,'(A)') filename

   
    !filename = 'matrizes/Matriz_A.dat'

   filename = 'matrizes/exA_3.dat'

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
  allocate(A(n,n))
  
  ! Lendo matriz A
  do i = 1,n
   read(11,*)A(i,:)
  enddo
  CLOSE (11)
  
  write(*,'(A,/,A,/,A)') 'Insira ICOD:','ICOD = 1 : Metodo da potencia ',&
   'ICOD = 2 : Metodo de Jacobi'
   read(*,*) ICOD

   if(ICOD .eq. 1) then
      ! Alocar espaco para 1 autovetor e 1 autovalor nas matrizes criadas X e E.
      allocate(X(n,1) , E(1))
      ! Chamar subrotina com metodo da potencia
      call power_method(A,E(1),X(:,1))
      write(*,'(/,A, f8.4)') 'Autovalor encontrado pelo metodo da potencia: ', E(1)
      write(*,*) 'Autovetor encontrado pelo metodo da potencia: ', X(:,1)

      ! call LU_decomposition_inplace(A)
   else if(ICOD .eq. 2) then
      ! call cholenskyDecomp(A)
      ! Alocar espaco para autovetores e autovalores esperados do metodo de jacobi.
      allocate(X(n,n) , E(n))
      ! Chamar subrotina com metodo de Jacobi.
      call jacobi_method(A,E,X)
      write(*,'(A,/)') 'Segue autovalores de A e, em seguida, matriz com autovetores associados em ordem.'
      write(*,*) E(:)
      write(*,'(/)')
      do i = 1,n
        write(*,*) X(i,:) 
      enddo
   endif


   

   end program