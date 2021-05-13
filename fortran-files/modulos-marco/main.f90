program main
  use lib
  implicit none

  real(8) :: A(5,6), v(5)


  A(1,:) = (/180, 14, 9, 0, 0, 0/)
  A(2,:) = (/130, 87, 8, 0, 0, 0/)
  A(3,:) = (/100, 89, 9, 0, 0, 0/)
  A(4,:) = (/80, 14, 9, 0, 0, 0/)
  A(5,:) = (/40, 60, 9, 0, 0, 0/)




  call valores_emsr(A, 5)

  call proteger(A,5,v)


  print *, "MATRIZ VALORES EMSR: "

  call show_array(A, 5, 6)

  print *, "NIVELES DE PROTECCION"

  call show_array(v,5, 1)



  print *, "Vagones: ", delta(A, 5, v, dble(1.5), 500, 80, 320)


end program
