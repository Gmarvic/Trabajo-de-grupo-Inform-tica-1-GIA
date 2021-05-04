program main
  use lib
  implicit none

  real(8) :: tabla(5,6), q(100, 5), pasajerostotal(100), t, dinero(100), tasapax, v(5)
  integer :: i, n, nvagones, tasavag, nplazas, j, clases, libres, k
  n = 100
  nplazas = 80
  clases = 5

  print *, "nvagones?"
  read *, nvagones

  print *, "tasa pasajero?"
  read *, tasapax

  print *, "tasa vagon?"
  read *, tasavag

  t = tasavag/dble(nplazas) - tasapax

  tabla(1,:) = (/dble(180) - t, dble(14), dble(9), dble(0), dble(0), dble(0)/)
  tabla(2,:) = (/dble(130) - t, dble(87), dble(8), dble(0), dble(0), dble(0)/)
  tabla(3,:) = (/dble(100) - t, dble(89), dble(9), dble(0), dble(0), dble(0)/)
  tabla(4,:) = (/dble(80) - t, dble(14), dble(9), dble(0), dble(0), dble(0)/)
  tabla(5,:) = (/dble(40) - t, dble(60), dble(9), dble(0), dble(0), dble(0)/)



  call valores_emsr(tabla, 5)
  call proteger(tabla,5,v)
  print *, "MATRIZ VALORES EMSR: "

  call show_array(tabla, 5, 6)

  print *, "NIVELES DE PROTECCION"

  do i = 1, 5
    if (v(i) == 0.d0) v(i) = tabla(i, 5)+tabla(i,6)
  end do
  call show_array(v,5, 1)






  ! opening the file for reading
  open (2, file = 'seed.dat', status = 'old')
  do i = 1,n
    do j = 1, clases
      read(2,*) q(i,j)
    end do
  end do
  close(2)


  do i = 1, n
    do j = 1, clases
      call randNormal(tabla(j, 2), tabla(j, 3), q(i,j))
    end do
  end do

  do i = 1, n
    t = 0.d0
    do j = 1, clases
      t = t + q(i, j)
    end do
    pasajerostotal(i) =  t
    ! print *, "Pasajeros en total: ", pasajerostotal(i)

    t = 0.d0
    do j = 1, clases
      t = t + q(i, j)*(tabla(j, 1))
    end do
    dinero(i) = t
    ! print *, dinero(i), pasajerostotal(i)
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
    libres = nvagones*nplazas
    do j = clases-1, 1, -1
      t = 0.d0
      do k = j+2, clases
        t = t + q(i, k)
      end do

      libres = nvagones*nplazas - v(j) - t

      if (q(i, j+1) > libres) then
        q(i, j+1) = libres
        if (libres < 0) q(i, j+1) = 0.d0
      end if

      dinero(i) = dinero(i) + q(i, j+1)*(tabla(j+1, 1))
    end do

    t = 0.d0
    do k = 2, clases
      t = t + q(i, k)
    end do

    libres = nvagones*nplazas - t

    if (q(i, 1) > libres) then
      q(i, 1) = libres
      if (libres < 0) q(i, 1) = 0.d0
    end if
    dinero(i) = dinero(i) + q(i, 1)*tabla(1, 1)

    t = 0.d0

    do j = 1, clases
      t = t + q(i, j)
    end do
    pasajerostotal(i) = t


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

  do i = 1, clases
    t = 0.d0
    do j = 1, n
      t = t + q(j, i)
    end do
    print *, i, ":", t/dble(n)
  end do




end program
