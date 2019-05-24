!--provides divt, jz=1:nz-1
subroutine divstress_uv (divt, tx, ty, tz)
use types,only:rprec
use param,only:ld,ny,nz, BOGUS, VERBOSE, coord !DY coord added by Di Yang for debugging
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(out)::divt
real (rprec), dimension (ld, ny, $lbz:nz), intent (in) :: tx, ty, tz
! sc: we should be able to save some memory here!
! do this by writing a divergence subroutine--then do not store derivs 
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dtxdx,dtydy, dtzdz

logical, parameter :: DEBUG = .true.

!DY Added by Di Yang for debugging
integer (rprec)::jx,jy,jz
!DY end here

if (VERBOSE) write (*, *) 'started divstress_uv'
 
!if (DEBUG) then
!  if (abs (dtxdx(ld, 1, 1) - 0._rprec) > 0.0001_rprec) then
!    write (*, *) 'dtxdx (ld, 1, 1) is not 0'
!    stop
!  end if
!end if

!DY Added by Di Yang for debugging
!!$if(coord.eq.0) then
!!$   do jz=1,1
!!$      do jy=1,ny
!!$         do jx=1,ld
!!$            write(100,*) "In divstress_uv: i,j,k=",jx,jy,jz
!!$            write(100,*) "In divstress_uv: tx,ty,tz=",tx(jx,jy,jz),ty(jx,jy,jz),tz(jx,jy,jz)
!!$         enddo
!!$      enddo
!!$   enddo
!!$endif
!DY end here

! compute stress gradients      
!--MPI: tx 1:nz-1 => dtxdx 1:nz-1
call ddx(dtxdx, tx)  !--really should replace with ddxy (save an fft)
!$if ($MPI)
!  dtdx(:, :, 0) = BOGUS
!$endif
!dtxdx(:, :, nz) = BOGUS

!DY Added by Di Yang for debugging
!!$if(coord.eq.0) then
!!$   do jz=1,1
!!$      do jy=1,ny
!!$         do jx=1,ld
!!$            write(101,*) "In divstress_uv: i,j,k=",jx,jy,jz
!!$            write(101,*) "In divstress_uv: tx,ty,tz=",tx(jx,jy,jz),ty(jx,jy,jz),tz(jx,jy,jz)
!!$            write(101,*) "In divsress_uv: dtxdx,dtydy=",dtxdx(jx,jy,jz),dtydy(jx,jy,jz)
!!$         enddo
!!$      enddo
!!$   enddo
!!$endif
!DY end here

!--MPI: ty 1:nz-1 => dtdy 1:nz-1
call ddy(dtydy, ty)
!$if ($MPI)
!  dtdy(:, :, 0) = BOGUS
!$endif
!dtydy(:, :, nz) = BOGUS

!DY Added by Di Yang for debugging
!!$if(coord.eq.0) then
!!$   do jz=1,1
!!$      do jy=1,ny
!!$         do jx=1,ld
!!$            write(102,*) "In divstress_uv: i,j,k=",jx,jy,jz
!!$            write(102,*) "In divstress_uv: tx,ty,tz=",tx(jx,jy,jz),ty(jx,jy,jz),tz(jx,jy,jz)
!!$            write(102,*) "In divsress_uv: dtxdx,dtydy=",dtxdx(jx,jy,jz),dtydy(jx,jy,jz)
!!$         enddo
!!$      enddo
!!$   enddo
!!$endif
!DY end here

!--MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(dtzdz, tz)
!$if ($MPI)
!  dtzdz(:, :, 0) = BOGUS
!$endif
!if (USE_MPI .and. coord < nproc-1) dtzdz(:, :, nz) = BOGUS

!--MPI following comment only true at bottom process
! the following gives bad results...but it seems like i the
! issue should be taken care of somewhere
! need to correct wall level, since tz will be on uv-node there
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

!DY Added by Di Yang for debugging
!!$if(coord.eq.0) then
!!$   do jz=1,1
!!$      do jy=1,ny
!!$         do jx=1,ld
!!$            write(103,*) "In divstress_uv: i,j,k=",jx,jy,jz
!!$            write(103,*) "In divstress_uv: tx,ty,tz=",tx(jx,jy,jz),ty(jx,jy,jz),tz(jx,jy,jz)
!!$            write(103,*) "In divsress_uv: dtxdx,dtydy,dtzdz=",dtxdx(jx,jy,jz),dtydy(jx,jy,jz),dtzdz(jx,jy,jz)
!!$         enddo
!!$      enddo
!!$   enddo
!!$endif
!DY end here

!--only 1:nz-1 are valid
divt(:, :, 1:nz-1) = dtxdx(:, :, 1:nz-1) + dtydy(:, :, 1:nz-1) +  &
                     dtzdz(:, :, 1:nz-1)

!--Set ld-1, ld to 0 (or could do BOGUS)
divt(ld-1:ld, :, 1:nz-1) = 0._rprec

$if ($MPI)
  divt(:, :, 0) = BOGUS
$endif
divt(:, :, nz) = BOGUS

if (VERBOSE) write (*, *) 'finished divstress_uv'

end subroutine divstress_uv
