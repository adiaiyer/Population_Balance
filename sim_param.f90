module sim_param
use types,only:rprec
use param,only:ld,ny,nz,lh, path,npcon
implicit none

save
public

!character(*),parameter::path='./'  !--moved into param

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !  some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real (rprec), dimension (ld, ny, $lbz:nz) :: u, v, w
real (rprec), dimension (ld, ny, $lbz:nz) :: ke, dragx,dragy,dragz
real (rprec), dimension (ld, ny, $lbz:nz) :: ke_temp
real (rprec), dimension (ld, ny, $lbz:nz) :: dudx=0._rprec, dudy=0._rprec, dudz=0._rprec,    &
                                             dvdx=0._rprec, dvdy=0._rprec, dvdz=0._rprec,    &
                                             dwdx=0._rprec, dwdy=0._rprec, dwdz=0_rprec,    &
                                             RHSx, RHSy, RHSz,    &
                                             RHSx_f, RHSy_f, RHSz_f
real (rprec), dimension (ld, ny, $lbz:nz) :: dudt,dvdt,dwdt

real (rprec), dimension (ld, ny, nz) :: dpdx=0._rprec,  &
                                        dpdy=0._rprec,  &
                                        dpdz=0._rprec
real (rprec), dimension (ld, ny, $lbz:nz) :: dkedx, dkedy, dkedz
real (rprec), dimension (ld, ny, $lbz:nz) :: txx, txy, tyy
real (rprec), dimension (ld, ny, $lbz:nz) :: txz, tyz, tzz

real(kind=rprec),dimension(ld,ny,0:nz)::p

! Added for scalars
!real(kind=rprec),dimension(ld,ny,nz)::theta,q
real(kind=rprec),dimension(ld,ny,$lbz:nz)::theta,q

real (rprec), dimension (ld, ny, $lbz:nz) :: divtx, divty, divtz


! Added for pollen - pollen concentration in grains/m3
! Chamecki - 08/01/2006
real(kind=rprec),dimension(ld,ny,$lbz:nz,npcon)::pcon
real(rprec),dimension(ld,ny,nz) :: fxa=0._rprec, &
                                   fya=0._rprec, &
                                   fza=0._rprec
end module sim_param
