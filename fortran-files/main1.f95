Program main
use utilidades
implicit none
!-----------------------------------------------------------------
!Inicialización de variables
!variables para la resolución del sistema de ecuaciones:
integer, parameter :: m=2 
integer :: i, n, nmax=100 
real*8 :: x(m), a(m,m), y(m), errf, errx, eps=1d-3
!------------------------------------------------------------------
!Resolución del sistema, promedio y desviacion clases A y D:
!asigna valores iniciales aproximados gráficamente
 x(1)= 10d0 
 x(2)= 10d0 
 n=0 
 !calcula valor de la función para aproximacion inicial 
 call f(x,y,m) 
 !comienza el ciclo de iteraciones 
 do n=1, nmax 
 ! calcula el jacobiano 
  call jacob_f(x,a,m) 
 ! calcula descomposicion LU del jacobiano 
  call factorizar (a,m) 
 ! resuelve el sistema df*dx = -f 
  call sustituir (a,y,m) 
 ! calcula la nueva aproximacion 
  do i=1, m 
    x(i) = x(i) + y(i) 
  end do 
 ! calcula el error relativo de la corrección 
  errx = norma_2(y,m)/(max(eps,norma_2(x,m))) 
 ! calcula la funcion para la nueva aproximacion 
  call f(x,y,m) 
 ! calcula el error de f 
  errf = norma_2(y,m) 
 ! escribe resultados 
  write (*,50) n, (x(i), i=1,m),errx,errf 
  if (max(errx,errf) .lt. eps)   exit 
 end do 
 if (n .gt. nmax) write(*,'("Se ha excedido el numero maximo de iteraciones")') 
50 format (i2,4(1x,f18.12)) 
 !------------------------------------------------------------------------------
 
 






end program
