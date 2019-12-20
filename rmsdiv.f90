subroutine rmsdiv(rms)
! actually, this is NOT the rms divergence of velocity, its like an
! l_1 norm or something.
use types,only:rprec
use param
use sim_param, only : du=>dudx, dv=>dvdy, dw=>dwdz
use debug_mod
implicit none
integer::jx,jy,jz
integer :: jz_max
real(kind=rprec)::rms
real(rprec),dimension(nx,ny,nz-1) :: rms_test=0._rprec
logical, parameter :: DEBUG = .false.
character (64) :: fname

$if ($MPI)
  real (rprec) :: rms_global
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  jz_max = nz-1
else
  !jz_max = nz
  jz_max = nz-1
end if

rms=0._rprec
do jz=1,jz_max
do jy=1,ny
do jx=1,nx
   rms=rms+abs(du(jx,jy,jz)+dv(jx,jy,jz)+dw(jx,jy,jz))
     rms_test(jx,jy,jz)= abs(du(jx,jy,jz)+dv(jx,jy,jz)+dw(jx,jy,jz)) 
end do
end do
end do
rms=rms/(nx*ny*(jz_max))
!if(coord.eq.41) then !write(*,*)"max_rms",maxval(rms_test)
!do jz=1,nz-1
!  write(*,*) jz,maxval(rms_test(xps,:,jz)),maxloc(rms_test(xps,:,jz)) , coord
!enddo
!endif

$if ($MPI)
  call mpi_reduce (rms, rms_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then
    rms = rms_global/nproc
    !write (*, *) 'rms_global = ', rms_global/nproc
  end if
  !if (rank == 0) rms = rms_global/nproc  !--its rank here, not coord
$endif

if (DEBUG) call DEBUG_write (du(1:nx, 1:ny, 1:nz) + dv(1:nx, 1:ny, 1:nz) +  &
                             dw(1:nx, 1:ny, 1:nz), 'rmsdiv')

end subroutine rmsdiv
