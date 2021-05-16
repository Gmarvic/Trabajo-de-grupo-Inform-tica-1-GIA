module utilidades
    contains
!-------------------------------------------------------------------------------
!Subrutinas para el cálculo del sistema no lineal:

!Norma euclídea de un vector de dimensión n:
  function norma_2(v,n)
         implicit none
         integer :: n, k
         real*8 :: v(n), norma_2, s
         s = 0d0
         do k = 1, n
          s = s + abs(v(k))**2
         end do
           norma_2 = sqrt(s)
  end function

!Subrutina de cálculo de valores de las funciones del sistema:
subroutine f(x,y,m)
       implicit none
       integer,intent(in) :: m
       real*8, intent(inout) :: x(m), y(m)
       !Funciones del sistema:
        y(1) = -(exp(x(1)/7)+2-x(2))
        y(2) = -(0.1*x(1)**3-20*x(1)+20-x(2))
end subroutine

!Subrutina con la Jacobiana calculada mediante diferencias finitas:
subroutine jacob_f(x,a,m)
       implicit none
          integer,intent(in) :: m
          real*8,intent(inout) :: x(m), a(m,m)
          real*8 :: h=1d-6,z(m),y1(m),y2(m)
          ! valor de la función en la iteracción dada
          call f(x,y1,m)
          ! calculo parciales respecto a x1
           z(1)=x(1)+h
           z(2)=x(2)
           call f(z,y2,m)
           a(1,1) = (-y2(1)+y1(1))/h
           a(2,1) = (-y2(2)+y1(2))/h
          ! calculo parciales respecto a x2
           z(1)=x(1)
           z(2)=x(2)+h
           call f(z,y2,m)
           a(1,2) = (-y2(1)+y1(1))/h
           a(2,2) = (-y2(2)+y1(2))/h
       end subroutine

!subrutinas de LU:

!factorización en lower y upper:
subroutine factorizar (a,n)
       implicit none
       ! variables y argumentos
          integer :: k, i, h, j, n
          real*8 :: a(n,n), s
          !procedemos con la factorizacion
            do k = 1, n
               do i = k, n
                 s = 0
                 do h =1, k-1
                    s = s + a(i,h) * a(h,k)
                 end do
               a(i,k) = a(i,k) - s
               end do
               do j = k+1, n
                 s = 0
                 do h = 1, k-1
                   s = s + a(k,h) * a(h,j)
                 end do
                 a(k,j) = (1. / a(k,k)) * (a(k,j) - s)
               end do
            end do
       end subroutine

!resolución del sistema:
subroutine sustituir (a,b,n)
            implicit none
            ! variables
            integer :: n, k, h
            real*8 :: a(n,n), b(n), s
            ! resolución:
              do k = 1, n
                 s = 0
                 do h = 1, k-1
                 s = s + a(k,h) * b(h)
                 end do
                 b(k) = (b(k) - s) / a(k,k)
              end do
              do k = n, 1, -1
                 s= 0
                 do h = k+1, n
                    s = s + a(k,h) * b(h)
                 end do
              b(k) = b(k) - s
              end do
      end subroutine

!-------------------------------------------------------------------------------
!Subrutinas de cálculo de los EMSRS y el número de vagones:

    ! Muestra matrices de reales (m,n) por pantalla
  subroutine show_array(A, m, n)
    implicit none
    integer, intent(in) :: m, n
    real(8), intent(in) :: A(m,n)
    CHARACTER(LEN=80) :: String, fstring, temp
    integer :: i

    String = "(A3,"//trim(str(n))//"F10.5, A3)"
    fstring = "(A3,A"//trim(str(n*10))//", A3)"

    temp = ""
    do i = 1, n
      temp = temp//" "
    end do

    write(*, fstring) "_", temp,"_ "
    do i = 1, m-1
      write(*,String) "| ", A(i,:), "|"
    end do
    write(*,String) "|_", A(m,:), " _|"
    print *,
  end subroutine

  ! Convierte un entero en una cadena de characteres
  character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str


  ! funcion distribución gauss, x, media, desviación, f(x)
  function gauss(x, mu, sigma) result(f)
    implicit none
    real(8), intent(in) :: x, mu, sigma
    real(8) :: f

    f = exp(-0.5 * ((x-mu)/sigma)**2)/(sigma*sqrt(2*acos(-1.0)))
  end function


  ! calcula la función de probabilidad acumulada, según t, para una distribución de gauss
  function probabilidad(t, mu, sigma) result(P)
    implicit none
    real(8), intent(in) :: t, mu, sigma
    real(8) :: P, s

    ! beware NaN
    call integralCDF(t, mu, sigma, 1000, s)
    P = 0.5 + s
    if (isnan(P)) P = dble(1)
  end function


  ! Según el método del trapecio se integra la probabilidad compuesta en función a mu y sigma
  ! Función de probabilidad acumulada, CDF
  ! límite de integración t, distribución(mu, sigma), número de intervalos n, solución s
  subroutine integralCDF(t,mu,sigma,n,s)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, mu, sigma
    real(8), intent(inout) :: s

    real(8) :: k, x, h
    integer :: i

    h = (t-mu)/dble(n)

    k = 0d0
    do i = 1, n-1
      x = mu + h*dble(i)
      k = k + gauss(x, mu, sigma)
    end do
    s = 0.5*h*(gauss(t, mu, sigma) + gauss(mu, mu, sigma) + 2.0*k)
  end subroutine



  ! Esta subrutina calcula la demanda agregada, tarifa media ponderada y la desviación conjunta, según el algoritmo EMSR-b, desde y hacia la tabla A(filas, 6), ordenada descendentemente en función de las tarifas de las clases.
  ! La tabla A() se usará para calcular los niveles de protección para cada clase.

  subroutine valores_emsr(A, filas)
    implicit none
    integer, intent(in) :: filas
    real(8), intent(inout) :: A(filas, 6)

    integer :: i, k
    real(8) :: s

    do i = 1, filas

      !Cálculo de la demanda agregada:
      s = 0.d0
      do k = 1, i
        s = s + A(k, 2)
      end do

      A(i,5) = s

      !Cálculo de la tarifa media ponderada
      s = 0.d0
      do k = 1, i
        s = s + A(k, 1)*A(k, 2)

      end do

      A(i,4) = s/A(i,5)

      !Cálculo de la desviación conjunta:
      s = 0.d0
      do k = 1, i
        s = s + A(k, 3)*A(k, 3)
      end do

      A(i, 6) = sqrt(s)
    end do
  end subroutine


  ! aplica la regla de Littlewood; compara cada clase al conjunto restante según el algoritmo EMSRb;
  ! recibe como parámetros la tabla de los valores emsr calculados en la subrutina valores_emsr, el número de filas (clases) que contiene la tabla, y el vector solución v que contiene los niveles de protección para cada clase comparado con el resto.
  ! véase:
  ! https://youtu.be/mZY4CU05PLw
  subroutine proteger (A, filas, v)
    implicit none
    integer, intent(in) :: filas
    real(8), intent(in) :: A(filas, 6)
    real(8), intent(inout) :: v(filas)

    integer :: i, j
    real(8) :: s, x, tol, dx

    v = 0.d0
    ! se define la tolerancia
    tol = dble(0.0001)
    do i = 1, filas-1

      ! solución inicial centrada en mu siempre garantiza convergencia
      ! sólo puede existir una solución
      x = A(i,5)

      ! resolución iterativa según la linealización del complementario cdf, quasi método de Newton

      do j = 1, 50
        ! sólo puede tolerar valores entre (0,1). Si no se cumple esta condición y/o no se encuentra solución después de 50 iteraciones, devuelve por defecto 0.d0
        s = A(i+1, 1)/A(i, 4) - 1 + probabilidad(x, A(i,5), A(i,6))

        if (abs(s) < tol) then
          v(i) = x
          exit
        end if

        dx = -gauss(x, A(i,5), A(i,6))
        x = s/dx + x
      end do
    end do
  end subroutine


  ! función delta, toma los datos de las tarifas de la tabla A, los niveles de protección, v, las tasas, el número de plazas por vagón, número total de asientos, n, resultado d, delta(theta)

  ! la función espera las tarifas de cara al público y no los ingresos (i.e. espera las tarifas, sin los costes añadidos)

  function delta(A, filas, v, tasapax, tasavag, nplazas, n) result(nvagones)
    implicit none

    integer, intent(in) :: filas, n, nplazas, tasavag
    real(8), intent(inout) :: v(filas)
    real(8), intent(in) :: tasapax, A(filas, 6)

    real(8) :: d(n, 2)


    integer :: i, k, clase, nvagones
    real(8) :: s, maxvalue

    maxvalue = 0.d0

    do i = 1, filas
      if (v(i) == 0.d0) v(i) = n
    end do

    s = 0.d0
    clase = 1

    do i = 1, n

      !
      if (i > v(clase)) clase = clase + 1
      if (ceiling(i/dble(nplazas)) - ceiling((i-1)/dble(nplazas)) == 1) s = s - tasavag

      s = s + (A(clase, 1) - tasapax)*(1 - probabilidad(dble(i), A(clase,5), A(clase,3)))

      d(i,1) = i
      d(i, 2) = s

      if (d(i, 2) > maxvalue) then
        k = i
        maxvalue = d(i, 2)
      end if
    end do

    nvagones = ceiling(k/dble(nplazas))

  end function



  ! devuelve un valor normal, según el promedio, la desviación y el valor aleatorio lineal r, el cuál es devuelto como correspondiente

  subroutine randNormal(mu, sigma, r)
    implicit none
    real(8), intent(in) :: mu, sigma
    real(8), intent(inout) :: r

    integer :: j
    real(8) :: s, x, tol, dx

    ! se define la tolerancia
    tol = dble(0.000001)
    x = mu
    do j = 1, 50
      s = r - 1 + probabilidad(x, mu, sigma)

      if (abs(s) < tol) then
        if (x < 0.d0) then
          r = 0.d0
        else
          r = x
        end if
        exit
      end if

      dx = -gauss(x, mu, sigma)
      x = s/dx + x
    end do

  end subroutine

  !Calcula en número de plazas por clase, en base a las protecciones y al número de vagones.

  subroutine asientos(v,n,a)
    implicit none

    integer,intent(in)::n
    real*8,intent(inout)::v(n)
    real*8,intent(in)::a(n,6)

    integer::i=0,m, filas
    real*8::s=0

    filas=80*delta(A, 5, v, dble(1.5), 500, 80, 320)

    do i=1,n
        !Se redondea la protección al entero más alto
        v(i)=int(v(i))+1
        s=s+v(i)

        if (s>=filas)then
        v(i)=dble(filas)-s+v(i)
        do m=i+1,n
            v(m)=0d0
        end do
        exit
        end if

    end do

  end subroutine

!subrutinas de cosas para Fernando (solución probabilística):
!--------------------------------------------------------------------
subroutine EMSRsProbabilidad(A, filas, n, d)
    implicit none

    integer, intent(in) :: filas, n
    real(8), intent(inout) :: d(n, 6)
    real(8), intent(in) ::  A(filas, 6)

    integer :: i,clase
    real(8) :: s

    s = 0.d0
  do clase=1, 5
    do i = 1, n
      s = (A(clase, 4))*(1 - probabilidad(dble(i), A(clase,5), A(clase,6)))

      d(i,1) = i
      d(i, clase+1) = s
    end do
  end do

  end subroutine

  subroutine write_array(A, m, n)
    implicit none
    integer, intent(in) :: m, n
    real(8), intent(in) :: A(m,n)
    CHARACTER(LEN=80) :: String, fstring, temp
    integer :: i

    String = "(A3,"//trim(str(n))//"F10.5, A3)"
    fstring = "(A3,A"//trim(str(n*10))//", A3)"

    temp = ""
    do i = 1, n
      temp = temp//" "
    end do

    write(10, fstring) "_", temp,"_ "
    do i = 1, m-1
      write(10,String) "| ", A(i,:), "|"
    end do
    write(10,String) "|_", A(m,:), " _|"
    print *,
  end subroutine

  !Clasifica los EMSRs para calcular probabilísticamente las protecciones
subroutine clasificar(b,ca, filas, n)
    implicit none
    integer,intent(in)::filas, n
    real*8,intent(in)::b(filas,6)
    real*8,intent(inout):: ca(n)
    integer*8:: clase, i
    real*8::s
    i=0
    clase=1
    s=0

    do

    if (clase-1>0) then
        s=s+ca(clase-1)
        if (s>=filas)then
        ca(clase-1)=dble(filas)-s+ca(clase-1)
        exit
        end if
    end if

    if (s>=filas)then
        ca(clase-1)=dble(filas)-s
        exit
    end if

    clase=clase+1

    do
        i=i+1
      if (b(i,clase)<b(i,clase+1))then
        ca(clase-1)=i-1
        exit
      end if
    end Do

    end Do

  end subroutine

end module
