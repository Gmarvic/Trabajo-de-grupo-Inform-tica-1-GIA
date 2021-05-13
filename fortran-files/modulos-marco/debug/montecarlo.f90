program main
  use lib
  implicit none

  real(8) :: tabla(5,6), q(200, 5), pasajerostotal(200), t, dinero(200), tasapax, v(5), pasajerosmedia
  real(8) :: factor, factorant
  real(8) ::  arrayA(320, 2), w, arrayB(320, 2)
  integer :: i, n, nvagones, tasavag, nplazas, j, clases, libres, k, l, index, limit

  character(len=1) :: genfiles

  n = 200
  nplazas = 80
  clases = 5



  print *, "nvagones?"
  read *, limit



  print *, "tasa pasajero?"
  read *, tasapax

  print *, "tasa vagon?"
  read *, tasavag


  print *, "generar archivos? [y/n]"
  read *, genfiles



  do index = 1, limit

    nvagones = index
    factor = 1.0
    do l = 1, 10
      t = tasavag/dble(nplazas)/factor + tasapax
      ! t = 0.d0

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

        ! v(i) = tabla(i, 5)
      end do

      call show_array(v,5, 1)

      call montecarlo(tabla, v, q, dinero, n, nvagones, tasavag, tasapax, clases)


      do i = 1, n
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

    arrayB = delta(tabla, 5, v, tasapax, tasavag, nplazas, nplazas*nvagones)



    ! DELTA
    ! l = 1
    ! v(5) = 321
    !
    ! do i = 1, 320
    !   if (i > v(l)) l = l + 1
    !   t = i
    !
    !   arrayA(i, 1) = t
    !   arrayA(i, 2) = (1 - probabilidad(t, tabla(l,5), tabla(l,6)) + 1 - probabilidad(t,tabla(5,5), tabla(5, 6)))/2
    !   arrayA(i, 2) = 1 - probabilidad(t, tabla(l,5), tabla(l,6))
    !
    !   !print *, i, arrayA(i+1, 1), arrayA(i+1, 2)
    !
    ! end do
    !
    ! call vagones(tabla, 5, v, arrayB, 320, arrayA, tasavag)
    !
    !
    ! l = 1
    ! open(10, file = 'acc2.dat', status = 'new')
    ! do i = 1, 320
    !   write(10,*) i, arrayA(i, 2)
    !   l = l + 1
    ! end do
    ! close(1)





    if (genfiles == 'y') then
      ! monte carlo fill
      open(nvagones+2, file = "s"//trim(str(nvagones))//".dat")

      do l = 1 + (index - 1)*50, (index - 1)*50 + 50
        pasajerostotal = 0.d0
        arrayA = 0.d0
        t = 1.0
        w = 0.d0
        do j = clases, 1, -1
          if (q(l, j) > 1.0) then
            do i = int(t), int(q(l, j) + t)
              arrayA(i,1) = i
              !incluir vagones implicitamente o a 80
              arrayA(i, 2) = w + tabla(j, 1)
              w = arrayA(i, 2)
              ! print *, w
              pasajerostotal(l) = i
            end do
            t = t + int(q(l,j))
          end if
        end do


        ! do i=1, int(pasajerostotal(l))
        !   if (arrayA(i,1) > 0.d0) write(nvagones+2,*) arrayA(i,1), arrayA(i,2)
        ! end do
        write(nvagones+2,*) pasajerostotal(l), arrayA(int(pasajerostotal(l)),2)
        write(nvagones+2,*)
      end do
      close(nvagones+2)
    end if

  end do

  if (genfiles == 'y') then
    open(50, file = "deltastar7.dat")

    do i = 1, nplazas*nvagones
      write(50,*) i, arrayB(i, 2)
    end do

    close(50)
  end if

end program
