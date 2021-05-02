module lib
  implicit none

contains

  ! Muestra matrices de reales (m,n) por pantalla
  subroutine show_array(A, m, n)
    implicit none
    integer, intent(in) :: m, n
    real(8), intent(in) :: A(m,n)

    integer :: i

    do i = 1, m
      !print *, A(i,:)
    end do

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

    call integralCDF(t, mu, sigma, 1000, s)
    P = 0.5 + s
  end function


  ! Según el método de simpson se integra la probabilidad compuesta en función a mu y sigma
  ! Función de probabilidad acumulada, CDF
  ! límite de integración t, distribución(mu, sigma), número de intervalos n, solución s
  subroutine integralCDF(t,mu,sigma,n,s)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, mu, sigma
    real(8), intent(inout) :: s

    real(8) :: s1, s2, h, x
    integer :: i

    h = (t-mu)/dble(n)

    s1 = 0d0
    do i = 1, n-1,2
      x = mu + h*dble(i)
      s1 = s1 + gauss(x, mu, sigma)
    end do

    s2 = 0d0
    do i = 2, n-2, 2
      x = mu + h*dble(i)
      s2 = s2 + gauss(x, mu, sigma)
    end do

    s = (h/dble(3))*(gauss(mu, mu, sigma) + gauss(t, mu, sigma) + dble(4)*s1 + dble(2)*s2)
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
  !TODO !print aid
  subroutine proteger (A, filas, v)
    implicit none
    integer, intent(in) :: filas
    real(8), intent(in) :: A(filas, 6)
    real(8), intent(inout) :: v(filas)

    integer :: i, j
    real(8) :: s, x, tol, dx

    v = 0.d0
    ! se define la tolerancia
    tol = dble(0.0000000001)
    do i = 1, filas-1

      ! solución inicial centrada en mu siempre garantiza convergencia
      x = A(i,5)

      ! resolución iterativa según la linealización del complementario cdf, quasi método de Newton

      do j = 1, 50
        s = A(i+1, 1)/A(i, 4) - 1 + probabilidad(x, A(i,5), A(i,6))

        if (abs(s) < tol) then
          v(i) = x
          !print *, i, v(i), j

          !print *, A(i+1, 1)/A(i, 4)
          !print *, 1 - probabilidad(v(i), A(i,5), A(i,6))
          exit
        end if

        dx = -gauss(x, A(i,5), A(i,6))
        x = s/dx + x
      end do
    end do


    ! Ghost class TODO
    x = A(filas, 5)
    do j = 1, 100
      ! límite superior de probabilidad de ocupar theta asientos (0.1)
      s = 0.1 - 1 + probabilidad(x, A(filas,5), A(filas,6))

      if (abs(s) < tol) then
        v(filas) = x
        !print *, "GHOST CLASS"
        !print *, filas, v(filas), j

        !print *, 0.1
        !print *, 1 - probabilidad(v(filas), A(filas,5), A(filas,6))
        exit
      end if

      dx = -gauss(x, A(filas,5), A(filas,6))
      x = s/dx + x
    end do
  end subroutine


  subroutine vagones(A, filas, v, array, n, prob)
    implicit none
    integer, intent(in) :: n, filas
    ! TODO inout
    real(8), intent(inout) :: A(filas, 6), v(filas)
    real(8), intent(inout) :: array(n, 2), prob(n,2)

    integer :: nplazas
    integer :: tasavagon
    real(8) :: tasapax

    integer :: i, clase
    real(8) :: s, temp


    ! TODO
    v(filas) = dble(321)


    nplazas = 80
    tasavagon = 0
    tasapax = 0

    clase = 1

    s = 0.d0

    do i = 1, n

      if (i > v(clase)) clase = clase + 1

      if (ceiling(i/dble(nplazas)) - ceiling((i-1)/dble(nplazas)) == 1) then
        s = s - tasavagon
      end if

      temp = A(clase,1) - tasapax
      !print*, temp
      temp = temp*prob(i, 2)
      !print*, temp
      s = s + temp !A(clase, 1) - tasapax

      array(i,1) = i
      array(i, 2) = s

      !print *, i, s, clase

    end do



  end subroutine






end module
