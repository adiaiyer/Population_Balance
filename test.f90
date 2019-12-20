program test
implicit none

integer, parameter :: npcon=20
Real(kind=8),DIMENSION(npcon) :: Q_src,settling_vel
Real(kind=8) :: temp_qsrc
integer :: ip


  open(21,file="vel_sintef.dat",ACTION= "read")
  open(20,file="n_init_for_LES_flux.txt",ACTION= "read")
  do ip=1,npcon
    read(20,'(E12.5)') temp_qsrc
    read(21,*)settling_vel(ip)
    Q_src(ip) =  (temp_qsrc)
  enddo
 close(21) 
  close(20)
  print * , Q_src

  print *, settling_vel
end program test


