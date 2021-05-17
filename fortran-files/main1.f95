! A fortran95 program for G95
! By G4

program main
    use utilidades
  implicit none

!-----------------------------------------------------------------
!-----------------------------------------------------------------

!Inicialización de variables:

!variables para la resolución del sistema de ecuaciones:
integer, parameter :: m=2
integer :: i, n, nmax=100
real*8 :: x(m), b(m,m), y(m), errf, errx, eps=1d-3

!Variables para el cálculo EMSR-b:
 real(8) :: A(5,6), v(5), d(320,6), ca(5)

!------------------------------------------------------------------
!------------------------------------------------------------------

!Resolución del sistema, promedio y desviacion clases A y D:

!asigna valores iniciales aproximados gráficamente:
 x(1)= 10d0
 x(2)= 10d0
 n=0

 !calcula valor de la función para aproximacion inicial:
 call f(x,y,m)

 !comienza el ciclo de iteraciones:
 Print*, "Resolucion del sistema no lineal:"
 Print*, "-------------------------------------------------------------------------------"
 Print*,"Iteracion:   Demanda:    Desviacion:        Error demanda:      Error desviacion:  "
 do n=1, nmax

 ! calcula el jacobiano:
  call jacob_f(x,b,m)

 ! calcula descomposicion LU del jacobiano:
  call factorizar (b,m)

 ! resuelve el sistema df*dx = -f
  call sustituir (b,y,m)

 ! calcula la nueva aproximacion:
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
50 format (i2,4(1x,f18.5))

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 !Cálculo de los EMSRs mediante el metodo EMSR-B:

    print*, "-------------------------------------------------------------------------------"
    print*,""
    print *, "Se procede a usar el algoritmo EMSR - b:"
    print*,""

  !Se introducen los valores en la matriz A (media, desviación, etc):
  A(1,:) = (/180d0, x(1), x(2), 0d0, 0d0, 0d0/)
  A(2,:) = (/130d0, 87d0, 8d0, 0d0, 0d0, 0d0/)
  A(3,:) = (/100d0, 89d0, 9d0, 0d0, 0d0, 0d0/)
  A(4,:) = (/80d0, x(1), x(2), 0d0, 0d0, 0d0/)
  A(5,:) = (/40d0, 60d0, 9d0, 0d0, 0d0, 0d0/)

  !Se calculan y muestran las variables conjuntas necesarias para el EMSR-b:
  call valores_emsr(A, 5)
  print *, "Matriz de calculo EMSR-b (A-E): "
  Print*, ""
  Print*, "  Tarifa:  Demanda:  Desviación:  Tmp:  Demanda(tot.):  Desviacion(conj.):"
  call show_array(A, 5, 6)

  !Se calculan y muestran las protecciones mediante newton de cada clase:
  call proteger(A,5,v)
  print *, "Niveles de proteccion (A-E):"
  call show_array(v,5, 1)

  !Se calcula el número de vagones que interesa utilizar, según la función delta:
  print *, "Numero de vagones: ", delta(A, 5, v, dble(1.5), 500, 80, 320)
  print*,""

  !Se calcula y muestra la matriz de EMSRs y tras esta la recomendación inicial de asientos a comercializar:
  call EMSRsProbabilidad(A,5,320,d)
  Print*,"Matriz de EMSRs: "
  print*,""
  print*," Asientos:   Clase A:   Clase B:   Clase C:   Clase D:   Clase E:"
  call show_array(d,320,6)
  call clasificar(d,ca,320,5)
  print*,""
  print*,"Plazas por clase (A-E):"
  call show_array(ca,5,1)

!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------

end program
