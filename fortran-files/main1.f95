program main
    use utilidades
  implicit none

  !-----------------------------------------------------------------
!Inicialización de variables
!variables para la resolución del sistema de ecuaciones:
integer, parameter :: m=2
integer :: i, n, nmax=100
real*8 :: x(m), b(m,m), y(m), errf, errx, eps=1d-3
!Variables para el cálculo EMSR-b:
 real(8) :: A(5,6), v(5)
 !Variables de cosas para Fernando
 d(320,6), ca(5)
!------------------------------------------------------------------
!Resolución del sistema, promedio y desviacion clases A y D:
!asigna valores iniciales aproximados gráficamente
 x(1)= 10d0
 x(2)= 10d0
 n=0
 !calcula valor de la función para aproximacion inicial
 call f(x,y,m)
 !comienza el ciclo de iteraciones
 Print*, "Resolucion del sistema no lineal:"
 Print*, "-------------------------------------------------------------------------------"
 Print*,"Iteracion:   Demanda:    Desviacion:        Error demanda:      Error desviacion:  "
 do n=1, nmax
 ! calcula el jacobiano
  call jacob_f(x,b,m)
 ! calcula descomposicion LU del jacobiano
  call factorizar (b,m)
 ! resuelve el sistema df*dx = -f
  call sustituir (b,y,m)
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
50 format (i2,4(1x,f18.5))
 !------------------------------------------------------------------------------
 !Cálculo de los EMSRs mediante el metodo EMSR-B:
    Print*, "-------------------------------------------------------------------------------"
  print *, "Se procede a usar el algoritmo EMSR - b:"
    Print*,""
  A(1,:) = (/180d0, x(1), x(2), 0d0, 0d0, 0d0/)
  A(2,:) = (/130d0, 87d0, 8d0, 0d0, 0d0, 0d0/)
  A(3,:) = (/100d0, 89d0, 9d0, 0d0, 0d0, 0d0/)
  A(4,:) = (/80d0, x(1), x(2), 0d0, 0d0, 0d0/)
  A(5,:) = (/40d0, 60d0, 9d0, 0d0, 0d0, 0d0/)

  call valores_emsr(A, 5)

  call proteger(A,5,v)

  print *, "MATRIZ VALORES EMSR: "
  Print*, ""
  Print*, "  Tarifa:  Demanda:  Desviación:  Tmp:  Demanda(tot.):  Desviacion(conj.):"

  call show_array(A, 5, 6)

  print *, "NIVELES DE PROTECCION"

  call show_array(v,5, 1)

  print *, "Vagones: ", delta(A, 5, v, dble(1.5), 500, 80, 320)
  print*,""

  print*,"Plazas por clase:"
  call asientos(v,5,a)
  !Plazas en orden descendente de clase:
  call show_array(v,5,1)

!Cosas para Fernando:
!---------------------------------------------------------------
call EMSRsProbabilidad(A,5,320,d)
call show_array(d,320,6)
call clasificar(d,ca,320,5)
print*,""
print*,"Plazas por clase:"
call show_array(ca,5,1)
end program
