subroutine energy (kea)
use types,only:rprec
use param
use sim_param,only:u,v,w,ke
implicit none

logical, parameter :: DEBUG = .false.

integer::jx,jy,jz
integer :: jz_min

real(kind=rprec)::kea,denom,temp_w
$if ($MPI)
  real (rprec) :: ke_global
$endif

  jz_min = 1

kea=0._rprec
denom=nx*ny*(nz-jz_min)
do jz=jz_min,nz-1
do jy=1,ny
do jx=1,nx
   temp_w=.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
                                !--MPI: assumes w(jz=0) is in sync here
   ke(jx,jy,jz)=.5_rprec*(u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
   kea=kea+ke(jx,jy,jz)/denom

end do
end do
end do

$if ($MPI)

  call mpi_reduce (kea, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then  !--note its rank here, not coord
    kea = ke_global/nproc
    write (13, *) jt_total, kea
  end if

$else

  write (13, *) jt_total, kea

$endif

end subroutine energy
