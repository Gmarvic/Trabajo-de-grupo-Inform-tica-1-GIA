program main
  use lib
  implicit none

  real(8) :: r, a(100), b(100), c(100), d(100), e(100), pasajerostotal(100), t, dinero(100)
  integer :: i, n
  n = 100


  r = 0.3

  call randNormal(dble(60), dble(9), r)


  ! opening the file for reading
  open (2, file = 'seedA.dat', status = 'old')
  do i = 1,n
    read(2,*) a(i)
  end do
  close(2)

  open (2, file = 'seedB.dat', status = 'old')
  do i = 1,n
    read(2,*) b(i)
  end do
  close(2)

  open (2, file = 'seedC.dat', status = 'old')
  do i = 1,n
    read(2,*) c(i)
  end do
  close(2)
  open (2, file = 'seedD.dat', status = 'old')
  do i = 1,n
    read(2,*) d(i)
  end do
  close(2)
  open (2, file = 'seedE.dat', status = 'old')
  do i = 1,n
    read(2,*) e(i)
  end do
  close(2)


  do i = 1, n
    call randNormal(dble(14), dble(9), a(i))
    call randNormal(dble(87), dble(8), b(i))
    call randNormal(dble(89), dble(9), c(i))
    call randNormal(dble(14), dble(9), d(i))
    call randNormal(dble(60), dble(9), e(i))

    pasajerostotal(i) = a(i)+b(i)+c(i)+d(i)+e(i)
    ! print *, "Pasajeros en total: ", pasajerostotal(i)

    dinero(i) = a(i)*180 + b(i)*130 + c(i)*100 + d(i)*80 + e(i)*40
    dinero(i) = dinero(i) - 1.5*pasajerostotal(i) - 500*ceiling(pasajerostotal(i)/dble(80))
    print *, dinero(i), pasajerostotal(i)
  end do

  t = 0.d0
  do i = 1, n
    t = t + pasajerostotal(i)
  end do

  print *, "PASAJEROS DE MEDIA: ", t/dble(n)


  t = 0.d0
  do i = 1, n
    t = t + dinero(i)
  end do

  print *, "DINERO DE MEDIA: ", t/dble(n)


  do i = 1, n
    dinero(i) = 0.d0
    if (e(i) < 4*80 - 211) then
      dinero(i) = e(i)*(40 - 1.5)
    else
      e(i) = (4*80 - 211)
      dinero(i) = e(i)*(40 - 1.5)
    end if

    if (d(i) + e(i) < 4*80 - 183) then
      dinero(i) = dinero(i) + d(i)*(80 - 1.5)
    else
      d(i) = (4*80 - e(i) - 183)
      dinero(i) = dinero(i) + d(i)*(80 - 1.5)
    end if

    if (c(i) + d(i) + e(i) < 4*80 - 94) then
      dinero(i) = dinero(i) + c(i)*(100 - 1.5)
    else
      c(i) = (4*80 - e(i) - d(i) - 94)
      dinero(i) = dinero(i) + c(i)*(100 - 1.5)
    end if

    if (b(i) + c(i) + d(i) + e(i) < 4*80 - 9) then
      dinero(i) = dinero(i) + b(i)*(130 - 1.5)
    else
      b(i) = (4*80 - e(i) - d(i) - c(i) - 9)
      dinero(i) = dinero(i) + b(i)*(130 - 1.5)
    end if


    if (a(i) + b(i) + c(i) + d(i) + e(i) < 4*80) then
      dinero(i) = dinero(i) + a(i)*(180 - 1.5)
    else
      a(i) = (4*80 - e(i) - d(i) - c(i) - b(i))
      dinero(i) = dinero(i) + a(i)*(180 - 1.5)
    end if

    dinero(i) = dinero(i) - 500*4

    pasajerostotal(i) = a(i) + b(i) + c(i) + d(i) + e(i)

    print *, a(i), b(i), c(i), d(i), e(i)
    print *, dinero(i), pasajerostotal(i)

  end do

  t = 0.d0
  do i = 1, n
    t = t + pasajerostotal(i)
  end do

  print *, "PASAJEROS DE MEDIA: ", t/dble(n)


  t = 0.d0
  do i = 1, n
    t = t + dinero(i)
  end do
  print *, "DINERO DE MEDIA: ", t/dble(n)


  t = 0.d0
  do i = 1, n
    t = t + a(i)
  end do
  print *, "A: ", t/dble(n)
  t = 0.d0
  do i = 1, n
    t = t + b(i)
  end do
  print *, "B: ", t/dble(n)
  t = 0.d0
  do i = 1, n
    t = t + c(i)
  end do
  print *, "C: ", t/dble(n)
  t = 0.d0
  do i = 1, n
    t = t + d(i)
  end do
  print *, "D: ", t/dble(n)
  t = 0.d0
  do i = 1, n
    t = t + e(i)
  end do
  print *, "E: ", t/dble(n)


end program
