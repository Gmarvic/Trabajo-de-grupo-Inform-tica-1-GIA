program main
  use lib
  implicit none

  integer :: n, i, j, clase
  real(8) :: t, mu, sigma, h, arrayA(320, 2), arrayB(320, 2), A(5,6), v(5), B(4,6), tarifas(6,2)

  n = 320
  mu = 264
  sigma = 20

  h = 1.5*3*2*sigma/dble(n)



  ! do i = 1, n
  !   t = i
  !
  !
  !
  !   arrayA(i+1, 1) = t
  !   arrayA(i+1, 2) = 1 - probabilidad(t, mu, sigma)
  !
  !   !print *, i, arrayA(i+1, 1), arrayA(i+1, 2)
  !
  ! end do

  ! datos de la tabla de valores EMSR b

  A(1,:) = (/180, 14, 9, 0, 0, 0/)
  A(2,:) = (/130, 87, 8, 0, 0, 0/)
  A(3,:) = (/100, 89, 9, 0, 0, 0/)
  A(4,:) = (/80, 14, 9, 0, 0, 0/)
  A(5,:) = (/40, 60, 9, 0, 0, 0/)


  ! prueba de congruencia

  ! referencia:
  ! https://youtu.be/mZY4CU05PLw

  B(1,:) = (/500, 16, 5, 0, 0, 0/)
  B(2,:) = (/420, 44, 15, 0, 0, 0/)
  B(3,:) = (/290, 35, 11, 0, 0, 0/)
  B(4,:) = (/125, 0, 0, 0, 0, 0/)



  call valores_emsr(A, 5)

  call proteger(A,5,v)


  !print *, "MATRIZ VALORES EMSR: "

  call show_array(A, 5, 6)

  !print *, "NIVELES DE PROTECCION"

  call show_array(v,5, 1)


  v(5) = dble(321)
  clase = 1

  do i = 1, n
    if (i > v(clase)) clase = clase + 1
    t = i

    mu = A(clase, 5)
    sigma = A(clase, 3)

    arrayA(i, 1) = t
    arrayA(i, 2) = 1 - probabilidad(t, mu, sigma)

    !print *, i, arrayA(i+1, 1), arrayA(i+1, 2)

  end do



  ! arrayA(1,2) = 0.92

  call vagones(A, 5, v, arrayB, 320, arrayA)

  !
  !
  !
  ! do i = 1, 5
  !   tarifas(i,1) = A(i,1)
  !   tarifas(i,2) = A(i,4)
  ! end do
  !
  ! tarifas(6,1) = 9.990
  !
  !
  !
  ! do j = 1, 5
  !
  ! sigma = A(j, 6)
  ! mu = A(j, 5)
  !
  ! h = 3*2*sigma/dble(n)
  !
  ! do i = 0, n
  !   t = -3*sigma + i*h + mu
  !
  !
  !
  !   arrayA(i+1, 1) = t
  !   arrayA(i+1, 2) = -tarifas(j+1,1)/tarifas(j,2) + 1 - probabilidad(t, mu, sigma) !gauss(mu,sigma,t)!probabilidad(t, mu, sigma)
  !
  !   ! !print *, i, arrayA(i+1, 1), arrayA(i+1, 2)
  !
  ! end do
  !
  ! open(j, file = "data"//trim(str(j+10))//".dat", status = 'new')
  ! do i=1,n+1
  !    write(j,*) arrayA(i,1), arrayA(i,2)
  ! end do
  !
  ! close(j)
  !
  !
  !
  ! end do



  ! !print *, probabilidad(dble(0), dble(1), dble(2))

  ! call random_number(h)
  ! !print *, h
  !
  ! call random_number(h)
  ! !print *, h

  ! output data into a file

  ! do i = 1, n
  !   arrayB(i, 2) = arrayB(i,2)**(1/arrayA(i,2))!arrayB(i,2)*(arrayA(i,2)**(0.2))
  ! end do






   open(1, file = 'delta4.dat', status = 'new')
   do i=1,320
      write(1,*) arrayB(i,1), arrayB(i,2)
   end do

   close(1)


end program
