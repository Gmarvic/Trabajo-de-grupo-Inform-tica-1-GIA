Module utilidades 

contains 

 

subroutine emsr() 

¡Cálculo de 350*5 Emsr,s  

end Subroutine 

 

subroutine clasificar() 

end subroutine 

 

function trapecio(a,b,n,ro,mu) 

¡calcula la integral de la probabilidad 

  	¡a y b, extremos; n, número de subintervalos; ro, desviación; mu, demanda. 

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

¡función para el cálculo de probabilidades 

implicit none 

real*8, intent (in) :: x, ro, mu 

real*8 :: f 

real*8, parameter:: pi=acos(-1d0) 

f= (1d0/(ro*sqrt(2d0*pi)))*exp((-(x-mu)**2d0)/(2d0*(ro**2d0))) 

end function 

 

function norma(vector,n) 

¡norma euclídea de un vector de dimensión n 

integer, intent(in) :: n 

real*8, intent (in) :: vector(n) 

real*8 :: norma 

integer :: i 

norma=0d0 

do i=1, n 

norma=norma + vector(i)**2 

end do 

norma=sqrt(norma) 

end function 

 

subroutine escribir_matriz(a,n,m) 

¡ a, matriz; n, filas; m, columnas  

implicit none 

integer :: i, j, n, m 

real*8 :: a(n,m) 

do i = 1 , n 

write (*, 100) (a(i, j), j=1,m) 

100 format (20 (f12.7,1x)) 

end do 

end subroutine 

 

¡modulos de LU: 

 

subroutine factorizar (a,n) 

¡a, matriz; n, dimensión 

¡almacena en a la solución 

implicit none 

integer, intent(in) :: n 

real*8, intent(inout) :: a(n,n) 

integer :: k, i, h, j 

real*8 :: s 

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

a(k,j) = (1 / a(k,k)) * (a(k,j) - s) 

end do 

end do 

end subroutine 

 

subroutine sustituir (a,b,n) 

implicit none 

integer, intent(in) :: n 

real*8, intent(inout) :: a(n,n), b(n) 

integer :: k, h 

real*8 :: s 

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

 

End module 
