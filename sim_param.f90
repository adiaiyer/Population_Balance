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

logical :: sim_param_initialized = .false.

real (rprec), dimension (:, :, :),allocatable :: u, v, w
real (rprec), dimension (:, :, :),allocatable :: ke, ke_temp, dragx, dragy, dragz
real (rprec), dimension (:, :, :),allocatable :: dudx, dudy, dudz,   &
                                                 dvdx, dvdy, dvdz,   &
                                                 dwdx, dwdy, dwdz,    &
                                                 RHSx, RHSy, RHSz,    &
                                                 RHSx_f, RHSy_f, RHSz_f
real (rprec), dimension (:, :, :),allocatable :: dudt,dvdt,dwdt

real (rprec), dimension (:, :, :),allocatable :: dpdx, dpdy, dpdz
real (rprec), dimension (:, :, :),allocatable :: dkedx, dkedy, dkedz
real (rprec), dimension (:, :, :),allocatable :: txx, txy, tyy
real (rprec), dimension (:, :, :),allocatable :: txz, tyz, tzz
real (rprec), dimension (:, :, :),allocatable :: divtx, divty, divtz
real(kind=rprec),dimension(:,:,:),allocatable :: p

! Added for scalars
real(kind=rprec),dimension(:, :, :),allocatable::theta,q



! Added for pollen - pollen concentration in grains/m3 or # droplets/m^3
! Chamecki - 08/01/2006 , Aiyer - 06/05/2020
real(kind=rprec),dimension(:, :, :,:),allocatable::pcon

contains 
        
subroutine sim_param_init ()

       implicit none

allocate ( u(ld, ny, $lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, $lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, $lbz:nz) ); w = 0.0_rprec
allocate ( dudx(ld, ny, $lbz:nz) ); dudx = 0.0_rprec
allocate ( dudy(ld, ny, $lbz:nz) ); dudy = 0.0_rprec
allocate ( dudz(ld, ny, $lbz:nz) ); dudz = 0.0_rprec
allocate ( dvdx(ld, ny, $lbz:nz) ); dvdx = 0.0_rprec
allocate ( dvdy(ld, ny, $lbz:nz) ); dvdy = 0.0_rprec
allocate ( dvdz(ld, ny, $lbz:nz) ); dvdz = 0.0_rprec
allocate ( dwdx(ld, ny, $lbz:nz) ); dwdx = 0.0_rprec
allocate ( dwdy(ld, ny, $lbz:nz) ); dwdy = 0.0_rprec
allocate ( dwdz(ld, ny, $lbz:nz) ); dwdz = 0.0_rprec
allocate ( RHSx(ld, ny, $lbz:nz) ); RHSx = 0.0_rprec
allocate ( RHSy(ld, ny, $lbz:nz) ); RHSy = 0.0_rprec
allocate ( RHSz(ld, ny, $lbz:nz) ); RHSz = 0.0_rprec
allocate ( RHSx_f(ld, ny, $lbz:nz) ); RHSx_f = 0.0_rprec
allocate ( RHSy_f(ld, ny, $lbz:nz) ); RHSy_f = 0.0_rprec
allocate ( RHSz_f(ld, ny, $lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, $lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, $lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, $lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, $lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, $lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, $lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, $lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, $lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, $lbz:nz) ); divtz = 0.0_rprec

!! Added for scalars  - AA 06/2020
allocate ( theta(ld, ny, $lbz:nz) ); theta = 0.0_rprec
allocate ( q(ld, ny, $lbz:nz) ); q = 0.0_rprec

allocate ( pcon(ld, ny, $lbz:nz, npcon) ); pcon = 0.0_rprec
sim_param_initialized = .true.

end subroutine sim_param_init


end module sim_param
