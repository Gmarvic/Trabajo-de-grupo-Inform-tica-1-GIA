Module utilidades

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


  function trapecio(a,b,n,ro,mu)

    ! calcula la integral de la probabilidad

    ! a y b, extremos; n, número de subintervalos; ro, desviación; mu, demanda.

    implicit none

    integer, intent(in) :: n

    real*8, intent(in) :: a, b, ro, mu

    real*8 :: trapecio

    real*8 :: s, h, x

    integer :: i

    h=(b-a)/dble(n)

    ! obtenemos aproximación de la integral

    s=0d0

    do i=1, n-1

      x=a+h*dble(i)

      s=s+f(x,ro,mu)

    end do

    trapecio=h*(0.5*f(a,ro,mu)+0.5*f(b,ro,mu)+s)

  end function

  function f(x,ro,mu)

    ! función para el cálculo de probabilidades

    implicit none

    real*8, intent (in) :: x, ro, mu

    real*8 :: f

    real*8, parameter:: pi=acos(-1d0)

    f= (1d0/(ro*sqrt(2d0*pi)))*exp((-(x-mu)**2d0)/(2d0*(ro**2d0)))

  end function

  subroutine escribir_matriz(a,n,m)

    !  a, matriz; n, filas; m, columnas

    implicit none

    integer :: i, j, n, m

    real*8 :: a(n,m)

    do i = 1 , n

      write (*, 100) (a(i, j), j=1,m)

      100 format (20 (f12.7,1x))

    end do

  end subroutine
  end module
