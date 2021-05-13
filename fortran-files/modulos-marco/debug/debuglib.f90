module lib
  implicit none

contains

  ! Muestra matrices de reales (m,n) por pantalla
  subroutine show_array(A, m, n)
    implicit none
    integer, intent(in) :: m, n
    real(8), intent(in) :: A(m,n)
    CHARACTER(LEN=80) :: String, fstring, temp
    integer :: i

    String = "(A3,"//trim(str(n))//"F7.2, A3)"
    fstring = "(A3,A"//trim(str(n*7))//", A3)"

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


  ! Según el método de (simpson)[DEPRECATED] trapecio se integra la probabilidad compuesta en función a mu y sigma
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



  ! Esta subrutina calcula los valores emsr relevantes, según el algoritmo EMSR-b, desde y hacia la tabla A(filas, 6), ordenada descendentemente en función de las tarifas de las clases.
  ! La tabla A() se usará para calcular los niveles de protección para cada clase.

  subroutine valores_emsr(A, filas)
    implicit none
    integer, intent(in) :: filas
    real(8), intent(inout) :: A(filas, 6)

    integer :: i, k
    real(8) :: s

    do i = 1, filas

      s = 0.d0
      do k = 1, i
        s = s + A(k, 2)
      end do

      A(i,5) = s

      s = 0.d0
      do k = 1, i
        ! TODO ingreso en vez de tarifa
        s = s + A(k, 1)*A(k, 2)

      end do

      A(i,4) = s/A(i,5)

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


  ! la subrutina montecarlo llenará n trenes según Wik17, de las clases bajas a las altas aplicando los niveles de protección, v, con los valores de la tabla;
  ! cada fila de q es un tren, con cada columna representando cada producto, y la cantidad siendo el número de pasajeros. dinero es un array de los ingresos para cada tren. nvagones son el número de vagones con los que el tren está configurado. tasas, y número de productos, clases.

  ! toma valores aleatorios lineales de seed.dat. Son pseudo-aleatorios. Se pueden generar también con otros métodos; xn+1 = l*xn*(1-xn)

  subroutine montecarlo(tabla, v, q, dinero, n, nvagones, tasavag, tasapax, clases)
    implicit none

    integer, intent(in) :: n, nvagones, tasavag, clases
    real(8), intent(in) :: tasapax, tabla(clases, 6), v(clases)
    real(8), intent(inout) :: q(n,clases), dinero(n)

    integer :: i, j, k, libres, nplazas
    real(8) :: t

    nplazas = 80

    open (2, file = 'seed.dat', status = 'old')
    do i = 1, n
      do j = 1, clases
        read(2,*) q(i,j)
      end do
    end do
    close(2)

    do i = 1, n
      do j = 1, clases
        call randNormal(tabla(j, 2), tabla(j, 3), q(i,j))
        if (q(i,j) < 0) q(i,j) = 0.d0
      end do
    end do

    do i = 1, n
      dinero(i) = 0.d0
      libres = nvagones*nplazas

      do j = clases, 1, -1
        t = 0.d0
        do k = j+1, clases
          t = t + q(i, k)
        end do

        if (j == 1) then
          libres = nvagones*nplazas - t
        else
          libres = nvagones*nplazas - v(j-1) - t
        end if

        if (q(i, j) > libres) then
          q(i, j) = libres
          if (libres < 0) q(i, j) = 0.d0
        end if
        dinero(i) = dinero(i) + q(i, j)*(tabla(j, 1))
      end do
    end do

  end subroutine




end module
