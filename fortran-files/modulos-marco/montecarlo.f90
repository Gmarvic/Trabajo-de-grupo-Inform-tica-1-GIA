program main
  use lib
  implicit none

  real(8) :: tabla(5,6), q(200, 5), pasajerostotal(200), t, dinero(200), tasapax, v(5), pasajerosmedia, factor, factorant
  integer :: i, n, nvagones, tasavag, nplazas, j, clases, libres, k, tarifas(5), l
  n = 200
  nplazas = 80
  clases = 5

  print *, "nvagones?"
  read *, nvagones

  print *, "tasa pasajero?"
  read *, tasapax

  print *, "tasa vagon?"
  read *, tasavag

  factor = 1.0
  do l = 1, 10
  t = tasavag/dble(nplazas)/factor + tasapax
  ! t = 0.d0

  tabla(1,:) = (/dble(180) - t, dble(14), dble(9), dble(0), dble(0), dble(0)/)
  tabla(2,:) = (/dble(130) - t, dble(87), dble(8), dble(0), dble(0), dble(0)/)
  tabla(3,:) = (/dble(100) - t, dble(89), dble(9), dble(0), dble(0), dble(0)/)
  tabla(4,:) = (/dble(80) - t, dble(14), dble(9), dble(0), dble(0), dble(0)/)
  tabla(5,:) = (/dble(40) - t, dble(60), dble(9), dble(0), dble(0), dble(0)/)

  tarifas(:) = (/180, 130, 100, 80, 40/)




  call valores_emsr(tabla, 5)
  call proteger(tabla,5,v)
  print *, "MATRIZ VALORES EMSR: "

  call show_array(tabla, 5, 6)

  print *, "NIVELES DE PROTECCION"

  do i = 1, 5
    if (v(i) == 0.d0) v(i) = tabla(i, 5)+tabla(i,6)

    ! v(i) = tabla(i, 5)
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



  pasajerosmedia = t/dble(n)
  print *, "PASAJEROS DE MEDIA: ", pasajerosmedia


  t = 0.d0
  do i = 1, n
    t = t + dinero(i)
  end do
  print *, "INGRESOS MEDIOS", t/dble(n)

  factorant = factor
  factor = 1/(dble(nplazas*nvagones)/pasajerosmedia)
  print *, "factor de correct: ", factor

  do i = 1, clases
    t = 0.d0
    do j = 1, n
      t = t + q(j, i)
    end do
    print *, i, ":", t/dble(n)
  end do

  print *, "--------------------", l, "   --------------------"
  print *, "#vagones : ", nvagones, " tasapax : ", tasapax, " tasa vagones : ", tasavag
  print *, "--------------------------------------------------------"

  if (abs(factor-factorant) < 0.0001) exit
  end do

end program
