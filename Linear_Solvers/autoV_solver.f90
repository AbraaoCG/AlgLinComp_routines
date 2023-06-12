module mod1

    CONTAINS
       include 'routines/main_routines_autoV.f90'
       
end module mod1
    
program eigen_calculator
   use mod1
   implicit none
   integer ::n, i,j, ICOD,detCOD,io
   character*50 :: filename
   real*8 :: detA
   real*8, allocatable :: A(:,:), X(:,:), E(:)! X : matriz com autovetores ; E : vetor com autovalores

   write(*,*) 'Insira o caminho até o arquivo contendo a matrix A: '
   read(*,'(A)') filename
   ! Descobrindo número de linhas de A (n).
   n = 0
   OPEN (11,status = 'old', file = filename)
   DO
   READ(11,*,iostat=io)
   IF (io/=0) EXIT
   n = n + 1
   END DO
   rewind(11)
   ! Alocando espaço para A
  allocate(A(n,n))
  ! Lendo matriz A
  do i = 1,n
   read(11,*)A(i,:)
  enddo
  CLOSE (11)
  
  write(*,'(A,/,A,/,A)') 'Insira ICOD:','ICOD = 1 : Metodo da potencia ',&
   'ICOD = 2 : Metodo de Jacobi'
   read(*,*) ICOD

   write(*,'(A,/,A)') 'Deseja calculo de determinante?', 'Responda com 0 ou 1. Sim = 1 ; Não = 0 '
   read(*,*) detCOD

   if(ICOD .eq. 1) then
      ! Alocar espaco para 1 autovetor e 1 autovalor nas matrizes criadas X e E.
      allocate(X(n,1) , E(1))
      ! Chamar subrotina com metodo da potencia
      call power_method(A,E(1),X(:,1))
      write(*,'(/,A, f8.4)') 'Autovalor encontrado pelo metodo da potencia: ', E(1)
      write(*,*) 'Autovetor encontrado pelo metodo da potencia: ','[' ,X(:,1), ']'

   else if(ICOD .eq. 2) then
      ! Alocar espaco para autovetores e autovalores esperados do metodo de jacobi.
      allocate(X(n,n) , E(n))
      ! Chamar subrotina com metodo de Jacobi.
      call jacobi_method(A,E,X)
      write(*,'(A,/)') 'Segue autovalores de A e, em seguida, matriz com autovetores associados em ordem.'
      write(*,*),'[ ',E(:) , '  ]'
      write(*,'(/)')
      do i = 1,n
         do j = 1,n
        write(*,'(f8.4,A)', advance = 'no') X(i,j) ,'     '
         enddo
         write(*,'(/)')
      enddo

 
   endif
   if(detCOD .eq. 1) then
      if(ICOD .eq. 1) then
         write(*,*) "WARNING: The power method doesn't allow determinant calculation! Try Jacobi Method instead"
      else if(ICOD .eq. 2) then
         detA = 1
         do i =1,n
            detA = detA * E(i)
         enddo
         write(*,*) 'O determinante de A é ', detA
      endif
   endif
   end program