!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module stat_defs
!*******************************************************************************
use types, only : rprec
use param, only : nx, ny, nz, lh

save
public


type plane_t
    integer :: istart
    real(rprec) :: ldiff
end type plane_t

type zplane_t
    integer :: istart, coord
    real(rprec) :: ldiff
end type zplane_t

type spectra_t
    real(rprec), dimension(:), allocatable :: power
    integer :: istart, coord
    real(rprec) :: ldiff
end type spectra_t

real(rprec) :: spectra_total_time


!  Sums performed over time
type tavg_t
    real(rprec) :: u, v, w, u_w, v_w, w_uv, p
    real(rprec) :: txx, tyy, tzz, txy, txz, tyz
    real(rprec) :: u2, v2, w2, uv, uw, vw, fx, fy, fz
     real(rprec) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
end type tavg_t

$if ($PPCON)
type tavg_pcon_t
        real(rprec) :: pcon
        real(rprec) :: pcon2
end type tavg_pcon_t

type rs_pcon_t
        real(rprec) :: pconp2
end type rs_pcon_t

$endif

type rs_t
    real(rprec) :: up2, vp2, wp2, upvp, upwp, vpwp
end type rs_t

$if ($PPSGS) 
type tavg_sgs_t
    real(rprec) :: cs_opt2, Nu_t
end type tavg_sgs_t
$endif

$if ($PPBUDGET) 
type tavg_budget_t
    real(rprec) :: p
    real(rprec) :: uu, vv, ww, uv, uw, vw
    real(rprec) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(rprec) :: dpdx, dpdy, dpdz
    real(rprec) :: ududx, ududy, ududz, udvdx, udvdy, udvdz, udwdx, udwdy, udwdz
    real(rprec) :: vdudx, vdudy, vdudz, vdvdx, vdvdy, vdvdz, vdwdx, vdwdy, vdwdz
    real(rprec) :: wdudx, wdudy, wdudz, wdvdx, wdvdy, wdvdz, wdwdx, wdwdy, wdwdz
    real(rprec) :: uududx, uvdudy, uwdudz, uudvdx, uvdvdy, uwdvdz, uudwdx, uvdwdy, uwdwdz
    real(rprec) :: vududx, vvdudy, vwdudz, vudvdx, vvdvdy, vwdvdz, vudwdx, vvdwdy, vwdwdz
    real(rprec) :: wududx, wvdudy, wwdudz, wudvdx, wvdvdy, wwdvdz, wudwdx, wvdwdy, wwdwdz
    real(rprec) :: uxux, uyuy, uzuz, vxvx, vyvy, vzvz, wxwx, wywy, wzwz
    real(rprec) :: uxvx, uyvy, uzvz, uxwx, uywy, uzwz, vxwx, vywy, vzwz
    ! real(rprec) :: uyvx, uzwx, vzwy
    real(rprec) :: udpdx, udpdy, udpdz, vdpdx, vdpdy, vdpdz, wdpdx, wdpdy, wdpdz
    real(rprec) :: pdudx, pdudy, pdudz, pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz
    real(rprec) :: lapu, lapv, lapw
    real(rprec) :: ulapu, ulapv, ulapw, vlapu, vlapv, vlapw, wlapu, wlapv, wlapw
end type tavg_budget_t

type budget_t
    real(rprec) :: advxx, advyy, advzz, advxy, advxz, advyz
    real(rprec) :: tflucxx, tflucyy, tfluczz, tflucxy, tflucxz, tflucyz
    real(rprec) :: tpresxx, tpresyy, tpreszz, tpresxy, tpresxz, tpresyz
    real(rprec) :: pstrainxx, pstrainyy, pstrainzz, pstrainxy, pstrainxz, pstrainyz
    real(rprec) :: tviscxx, tviscyy, tvisczz, tviscxy, tviscxz, tviscyz
    real(rprec) :: prodxx, prodyy, prodzz, prodxy, prodxz, prodyz
    real(rprec) :: pdissxx, pdissyy, pdisszz, pdissxy, pdissxz, pdissyz
end type budget_t
$endif


! Create types for outputting data (instantaneous or averaged)
type(plane_t), allocatable, dimension(:) :: xplane, yplane
type(zplane_t), allocatable, dimension(:) :: zplane

type(tavg_t), allocatable, dimension(:,:,:) :: tavg
type(tavg_t), allocatable, dimension(:) :: tavg_zplane

$if ($PPCON)
type(tavg_pcon_t), allocatable, dimension (:,:,:,:) :: tavg_pcon
type(rs_pcon_t), allocatable, dimension (:,:,:,:) :: rs_pcon
$endif

type(rs_t), allocatable, dimension(:,:,:) :: rs
type(rs_t), allocatable, dimension(:) :: rs_zplane, cnpy_zplane

$if ($PPSGS) 
type(tavg_sgs_t), allocatable, dimension(:,:,:) :: tavg_sgs
$endif

$if ($PPBUDGET) 
type(tavg_budget_t), allocatable, dimension(:,:,:) :: tavg_budget
type(budget_t), allocatable, dimension(:,:,:) :: budget
$endif


real(rprec) :: tavg_total_time


! Time between calls of tavg_compute, built by summing dt
real(rprec) :: tavg_dt
! Switch for determining if time averaging has been initialized
logical :: tavg_initialized = .false.

contains

!///////////////////////////////////////////////////////////////////////////////
!/// Spectral RS operators
!///////////////////////////////////////////////////////////////////////////////

!*******************************************************************************
function rs_compute( a , lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_t), dimension(:,:,lbz2:), intent(in) :: a
type(rs_t), allocatable, dimension(:,:,:) :: c

integer :: ubx, uby, ubz

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)

allocate(c(ubx,uby,lbz2:ubz))

c % up2 = a % u2 - a % u * a % u
c % vp2 = a % v2 - a % v * a % v
c % wp2 = a % w2 - a % w * a % w
c % upvp = a % uv - a % u * a % v
!! using u_w and v_w below instead of u and v ensures that the Reynolds
!! stresses are on the same grid as the squared velocities (i.e., w-grid)
c % upwp = a % uw - a % u_w * a % w   !!pj
c % vpwp = a % vw - a % v_w * a % w   !!pj

end function rs_compute

$if ($PPCON)
function rs_compute_pcon( a , lbz2) result(c)
!*******************************************************************************
implicit none
integer, intent(in) :: lbz2
type(tavg_pcon_t), dimension(:,:,lbz2:,:), intent(in) :: a
type(rs_pcon_t), allocatable, dimension(:,:,:,:) :: c

integer :: ubx, uby, ubz,ubpcon

ubx=ubound(a,1)
uby=ubound(a,2)
ubz=ubound(a,3)
ubpcon = ubound(a,4)
allocate(c(ubx,uby,lbz2:ubz,ubpcon))

c % pconp2 = a % pcon2 - a%pcon * a%pcon

end function rs_compute_pcon

$endif

!#ifdef PPBUDGET
!*******************************************************************************
!function budget_compute( a, b, lbz2) result(c)
!!*******************************************************************************
!use param, only: nu_molec
!implicit none
!integer, intent(in) :: lbz2
!type(tavg_budget_t), dimension(:,:,lbz2:), intent(in) :: a
!type(tavg_t), dimension(:,:,lbz2:), intent(in) :: b
!type(budget_t), allocatable, dimension(:,:,:) :: c
!
!integer :: ubx, uby, ubz, i, j, k
!real(rprec) :: Rxx, Ryy, Rzz, Rxy, Rxz, Ryz
!real(rprec) :: Cududx,Cududy, Cududz, Cudvdx, Cudvdy, Cudvdz, Cudwdx, Cudwdy, Cudwdz
!real(rprec) :: Cvdudx,Cvdudy,Cvdudz, Cvdvdx, Cvdvdy, Cvdvdz, Cvdwdx, Cvdwdy, Cvdwdz 
!real(rprec) :: Cwdudx,Cwdudy,Cwdudz, Cwdvdx, Cwdvdy, Cwdvdz, Cwdwdx, Cwdwdy, Cwdwdz 
!real(rprec) :: Cuududx,Cuvdudy,Cuwdudz,Cuudvdx,Cuvdvdy,Cuwdvdz
!real(rprec) :: Cuudwdx, Cuvdwdy, Cuwdwdz
!real(rprec) :: Cvududx, Cvvdudy, Cvwdudz, Cvudvdx, Cvvdvdy, Cvwdvdz
!real(rprec) :: Cvudwdx, Cvvdwdy, Cvwdwdz
!real(rprec) :: Cwududx, Cwvdudy, Cwwdudz, Cwudvdx, Cwvdvdy, Cwwdvdz
!real(rprec) :: Cwudwdx, Cwvdwdy, Cwwdwdz
!real(rprec) :: Cudpdx,Cudpdy, Cudpdz, Cvdpdx, Cvdpdy, Cvdpdz, Cwdpdx, Cwdpdy, Cwdpdz
!real(rprec) :: Cpdudx,Cpdudy, Cpdudz, Cpdvdx, Cpdvdy, Cpdvdz, Cpdwdx, Cpdwdy, Cpdwdz
!real(rprec) :: Cuxux, Cuyuy, Cuzuz, Cvxvx, Cvyvy, Cvzvz, Cwxwx, Cwywy, Cwzwz
!real(rprec) :: Cuxvx, Cuyvy, Cuzvz, Cuxwx, Cuywy, Cuzwz, Cvxwx, Cvywy, Cvzwz
!real(rprec) :: Culapu,Culapv, Culapw, Cvlapu, Cvlapv, Cvlapw, Cwlapu, Cwlapv, Cwlapw
!
!ubx=ubound(a,1)
!uby=ubound(a,2)
!ubz=ubound(a,3)
!
!allocate(c(ubx,uby,lbz2:ubz))
!
!! ----------------------------- Reynolds Stress Budget -----------------------------
!
!do i = 1, ubx
!do j = 1, uby
!do k = 1, ubz
!
!! Compute intermediate terms used in the budget
!! Velocity-Velocity correlation, Reynolds stress
!! Similar to rs_compute, except all are on w-grid
!Rxx = a(i,j,k) % uu - b(i,j,k) % u_w * b(i,j,k) % u_w
!Ryy = a(i,j,k) % vv - b(i,j,k) % v_w * b(i,j,k) % v_w
!Rzz = a(i,j,k) % ww - b(i,j,k) % w   * b(i,j,k) % w
!Rxy = a(i,j,k) % uv - b(i,j,k) % u_w * b(i,j,k) % v_w
!Rxz = a(i,j,k) % uw - b(i,j,k) % u_w * b(i,j,k) % w
!Ryz = a(i,j,k) % vw - b(i,j,k) % v_w * b(i,j,k) % w
!
!! Velocity-Velocity Gradient correlation, ui*dujdxk
!Cududx = a(i,j,k)%ududx - b(i,j,k)%u_w * a(i,j,k)%dudx
!Cududy = a(i,j,k)%ududy - b(i,j,k)%u_w * a(i,j,k)%dudy
!Cududz = a(i,j,k)%ududz - b(i,j,k)%u_w * a(i,j,k)%dudz
!Cudvdx = a(i,j,k)%udvdx - b(i,j,k)%u_w * a(i,j,k)%dvdx
!Cudvdy = a(i,j,k)%udvdy - b(i,j,k)%u_w * a(i,j,k)%dvdy
!Cudvdz = a(i,j,k)%udvdz - b(i,j,k)%u_w * a(i,j,k)%dvdz
!Cudwdx = a(i,j,k)%udwdx - b(i,j,k)%u_w * a(i,j,k)%dwdx
!Cudwdy = a(i,j,k)%udwdy - b(i,j,k)%u_w * a(i,j,k)%dwdy
!Cudwdz = a(i,j,k)%udwdz - b(i,j,k)%u_w * a(i,j,k)%dwdz
!
!Cvdudx = a(i,j,k)%vdudx - b(i,j,k)%v_w * a(i,j,k)%dudx
!Cvdudy = a(i,j,k)%vdudy - b(i,j,k)%v_w * a(i,j,k)%dudy
!Cvdudz = a(i,j,k)%vdudz - b(i,j,k)%v_w * a(i,j,k)%dudz
!Cvdvdx = a(i,j,k)%vdvdx - b(i,j,k)%v_w * a(i,j,k)%dvdx
!Cvdvdy = a(i,j,k)%vdvdy - b(i,j,k)%v_w * a(i,j,k)%dvdy
!Cvdvdz = a(i,j,k)%vdvdz - b(i,j,k)%v_w * a(i,j,k)%dvdz
!Cvdwdx = a(i,j,k)%vdwdx - b(i,j,k)%v_w * a(i,j,k)%dwdx
!Cvdwdy = a(i,j,k)%vdwdy - b(i,j,k)%v_w * a(i,j,k)%dwdy
!Cvdwdz = a(i,j,k)%vdwdz - b(i,j,k)%v_w * a(i,j,k)%dwdz
!
!Cwdudx = a(i,j,k)%wdudx - b(i,j,k)%w * a(i,j,k)%dudx
!Cwdudy = a(i,j,k)%wdudy - b(i,j,k)%w * a(i,j,k)%dudy
!Cwdudz = a(i,j,k)%wdudz - b(i,j,k)%w * a(i,j,k)%dudz
!Cwdvdx = a(i,j,k)%wdvdx - b(i,j,k)%w * a(i,j,k)%dvdx
!Cwdvdy = a(i,j,k)%wdvdy - b(i,j,k)%w * a(i,j,k)%dvdy
!Cwdvdz = a(i,j,k)%wdvdz - b(i,j,k)%w * a(i,j,k)%dvdz
!Cwdwdx = a(i,j,k)%wdwdx - b(i,j,k)%w * a(i,j,k)%dwdx
!Cwdwdy = a(i,j,k)%wdwdy - b(i,j,k)%w * a(i,j,k)%dwdy
!Cwdwdz = a(i,j,k)%wdwdz - b(i,j,k)%w * a(i,j,k)%dwdz
!
!! vel-vel-velGrad triple correlation, ui*uk*dujdxk
!Cuududx = a(i,j,k)%uududx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dudx -        &
!    Rxx*a(i,j,k)%dudx - Cududx*b(i,j,k)%u_w - Cududx*b(i,j,k)%u_w
!Cuvdudy = a(i,j,k)%uvdudy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dudy -        &
!    Rxy*a(i,j,k)%dudy - Cududy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%u_w
!Cuwdudz = a(i,j,k)%uwdudz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dudz   -        &
!    Rxz*a(i,j,k)%dudz - Cududz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%u_w
!Cuudvdx = a(i,j,k)%uudvdx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dvdx -        &
!    Rxx*a(i,j,k)%dvdx - Cudvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%u_w
!Cuvdvdy = a(i,j,k)%uvdvdy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dvdy -        &
!    Rxy*a(i,j,k)%dvdy - Cudvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%u_w
!Cuwdvdz = a(i,j,k)%uwdvdz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dvdz   -        &
!    Rxz*a(i,j,k)%dvdz - Cudvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%u_w
!Cuudwdx = a(i,j,k)%uudwdx - b(i,j,k)%u_w*b(i,j,k)%u_w*a(i,j,k)%dwdx -        &
!    Rxx*a(i,j,k)%dwdx - Cudwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%u_w
!Cuvdwdy = a(i,j,k)%uvdwdy - b(i,j,k)%u_w*b(i,j,k)%v_w*a(i,j,k)%dwdy -        &
!    Rxy*a(i,j,k)%dwdy - Cudwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%u_w
!Cuwdwdz = a(i,j,k)%uwdwdz - b(i,j,k)%u_w*b(i,j,k)%w*a(i,j,k)%dwdz   -        &
!    Rxz*a(i,j,k)%dwdz - Cudwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%u_w
!
!Cvududx = a(i,j,k)%vududx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dudx -        &
!    Rxy*a(i,j,k)%dudx - Cvdudx*b(i,j,k)%u_w - Cududx*b(i,j,k)%v_w
!Cvvdudy = a(i,j,k)%vvdudy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dudy -        &
!    Ryy*a(i,j,k)%dudy - Cvdudy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%v_w
!Cvwdudz = a(i,j,k)%vwdudz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dudz   -        &
!    Ryz*a(i,j,k)%dudz - Cvdudz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%v_w
!Cvudvdx = a(i,j,k)%vudvdx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dvdx -        &
!    Rxy*a(i,j,k)%dvdx - Cvdvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%v_w
!Cvvdvdy = a(i,j,k)%vvdvdy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dvdy -        &
!    Ryy*a(i,j,k)%dvdy - Cvdvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%v_w
!Cvwdvdz = a(i,j,k)%vwdvdz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dvdz   -        &
!    Ryz*a(i,j,k)%dvdz - Cvdvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%v_w
!Cvudwdx = a(i,j,k)%vudwdx - b(i,j,k)%v_w*b(i,j,k)%u_w*a(i,j,k)%dwdx -        &
!    Rxy*a(i,j,k)%dwdx - Cvdwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%v_w
!Cvvdwdy = a(i,j,k)%vvdwdy - b(i,j,k)%v_w*b(i,j,k)%v_w*a(i,j,k)%dwdy -        &
!    Ryy*a(i,j,k)%dwdy - Cvdwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%v_w
!Cvwdwdz = a(i,j,k)%vwdwdz - b(i,j,k)%v_w*b(i,j,k)%w*a(i,j,k)%dwdz   -        &
!    Ryz*a(i,j,k)%dwdz - Cvdwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%v_w
!
!Cwududx = a(i,j,k)%wududx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dudx -          &
!    Rxz*a(i,j,k)%dudx - Cwdudx*b(i,j,k)%u_w - Cududx*b(i,j,k)%w
!Cwvdudy = a(i,j,k)%wvdudy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dudy -          &
!    Ryz*a(i,j,k)%dudy - Cwdudy*b(i,j,k)%v_w - Cvdudy*b(i,j,k)%w
!Cwwdudz = a(i,j,k)%wwdudz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dudz   -          &
!    Rzz*a(i,j,k)%dudz - Cwdudz*b(i,j,k)%w   - Cwdudz*b(i,j,k)%w
!Cwudvdx = a(i,j,k)%wudvdx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dvdx -          &
!    Rxz*a(i,j,k)%dvdx - Cwdvdx*b(i,j,k)%u_w - Cudvdx*b(i,j,k)%w
!Cwvdvdy = a(i,j,k)%wvdvdy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dvdy -          &
!    Ryz*a(i,j,k)%dvdy - Cwdvdy*b(i,j,k)%v_w - Cvdvdy*b(i,j,k)%w
!Cwwdvdz = a(i,j,k)%wwdvdz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dvdz   -          &
!    Rzz*a(i,j,k)%dvdz - Cwdvdz*b(i,j,k)%w   - Cwdvdz*b(i,j,k)%w
!Cwudwdx = a(i,j,k)%wudwdx - b(i,j,k)%w*b(i,j,k)%u_w*a(i,j,k)%dwdx -          &
!    Rxz*a(i,j,k)%dwdx - Cwdwdx*b(i,j,k)%u_w - Cudwdx*b(i,j,k)%w
!Cwvdwdy = a(i,j,k)%wvdwdy - b(i,j,k)%w*b(i,j,k)%v_w*a(i,j,k)%dwdy -          &
!    Ryz*a(i,j,k)%dwdy - Cwdwdy*b(i,j,k)%v_w - Cvdwdy*b(i,j,k)%w
!Cwwdwdz = a(i,j,k)%wwdwdz - b(i,j,k)%w*b(i,j,k)%w*a(i,j,k)%dwdz   -          &
!    Rzz*a(i,j,k)%dwdz - Cwdwdz*b(i,j,k)%w   - Cwdwdz*b(i,j,k)%w
!
!! vel-presGrad correlation, ui*dpdxj
!Cudpdx = a(i,j,k)%udpdx - b(i,j,k)%u_w * a(i,j,k)%dpdx
!Cudpdy = a(i,j,k)%udpdy - b(i,j,k)%u_w * a(i,j,k)%dpdy
!Cudpdz = a(i,j,k)%udpdz - b(i,j,k)%u_w * a(i,j,k)%dpdz
!Cvdpdx = a(i,j,k)%vdpdx - b(i,j,k)%v_w * a(i,j,k)%dpdx
!Cvdpdy = a(i,j,k)%vdpdy - b(i,j,k)%v_w * a(i,j,k)%dpdy
!Cvdpdz = a(i,j,k)%vdpdz - b(i,j,k)%v_w * a(i,j,k)%dpdz
!Cwdpdx = a(i,j,k)%wdpdx - b(i,j,k)%w * a(i,j,k)%dpdx
!Cwdpdy = a(i,j,k)%wdpdy - b(i,j,k)%w * a(i,j,k)%dpdy
!Cwdpdz = a(i,j,k)%wdpdz - b(i,j,k)%w * a(i,j,k)%dpdz
!
!! pres-velGrad correlation, p*duidxj
!Cpdudx = a(i,j,k)%pdudx - a(i,j,k)%p * a(i,j,k)%dudx
!Cpdudy = a(i,j,k)%pdudy - a(i,j,k)%p * a(i,j,k)%dudy
!Cpdudz = a(i,j,k)%pdudz - a(i,j,k)%p * a(i,j,k)%dudz
!Cpdvdx = a(i,j,k)%pdvdx - a(i,j,k)%p * a(i,j,k)%dvdx
!Cpdvdy = a(i,j,k)%pdvdy - a(i,j,k)%p * a(i,j,k)%dvdy
!Cpdvdz = a(i,j,k)%pdvdz - a(i,j,k)%p * a(i,j,k)%dvdz
!Cpdwdx = a(i,j,k)%pdwdx - a(i,j,k)%p * a(i,j,k)%dwdx
!Cpdwdy = a(i,j,k)%pdwdy - a(i,j,k)%p * a(i,j,k)%dwdy
!Cpdwdz = a(i,j,k)%pdwdz - a(i,j,k)%p * a(i,j,k)%dwdz
!
!! velGrad-velGrad correlation, duidxk*dujdxk, i=j
!Cuxux = a(i,j,k)%uxux - (a(i,j,k)%dudx**2)
!Cuyuy = a(i,j,k)%uyuy - (a(i,j,k)%dudy**2)
!Cuzuz = a(i,j,k)%uzuz - (a(i,j,k)%dudz**2)
!Cvxvx = a(i,j,k)%vxvx - (a(i,j,k)%dvdx**2)
!Cvyvy = a(i,j,k)%vyvy - (a(i,j,k)%dvdy**2)
!Cvzvz = a(i,j,k)%vzvz - (a(i,j,k)%dvdz**2)
!Cwxwx = a(i,j,k)%wxwx - (a(i,j,k)%dwdx**2)
!Cwywy = a(i,j,k)%wywy - (a(i,j,k)%dwdy**2)
!Cwzwz = a(i,j,k)%wzwz - (a(i,j,k)%dwdz**2)
!
!! velGrad-velGrad correlation, duidxk*dujdxk, i/=j
!Cuxvx = a(i,j,k)%uxvx - a(i,j,k)%dudx * a(i,j,k)%dvdx
!Cuyvy = a(i,j,k)%uyvy - a(i,j,k)%dudy * a(i,j,k)%dvdy
!Cuzvz = a(i,j,k)%uzvz - a(i,j,k)%dudz * a(i,j,k)%dvdz
!Cuxwx = a(i,j,k)%uxwx - a(i,j,k)%dudx * a(i,j,k)%dwdx
!Cuywy = a(i,j,k)%uywy - a(i,j,k)%dudy * a(i,j,k)%dwdy
!Cuzwz = a(i,j,k)%uzwz - a(i,j,k)%dudz * a(i,j,k)%dwdz
!Cvxwx = a(i,j,k)%vxwx - a(i,j,k)%dvdx * a(i,j,k)%dwdx
!Cvywy = a(i,j,k)%vywy - a(i,j,k)%dvdy * a(i,j,k)%dwdy
!Cvzwz = a(i,j,k)%vzwz - a(i,j,k)%dvdz * a(i,j,k)%dwdz
!
!! vel-Laplacian correlation, ui*lap(uj)
!Culapu = a(i,j,k)%ulapu - b(i,j,k)%u_w * a(i,j,k)%lapu
!Culapv = a(i,j,k)%ulapv - b(i,j,k)%u_w * a(i,j,k)%lapv
!Culapw = a(i,j,k)%ulapw - b(i,j,k)%u_w * a(i,j,k)%lapw
!Cvlapu = a(i,j,k)%vlapu - b(i,j,k)%v_w * a(i,j,k)%lapu
!Cvlapv = a(i,j,k)%vlapv - b(i,j,k)%v_w * a(i,j,k)%lapv
!Cvlapw = a(i,j,k)%vlapw - b(i,j,k)%v_w * a(i,j,k)%lapw
!Cwlapu = a(i,j,k)%wlapu - b(i,j,k)%w * a(i,j,k)%lapu
!Cwlapv = a(i,j,k)%wlapv - b(i,j,k)%w * a(i,j,k)%lapv
!Cwlapw = a(i,j,k)%wlapw - b(i,j,k)%w * a(i,j,k)%lapw
!
!! Now compute terms to be outputted
!! Advection, Uk*dRijdxk
!c(i,j,k) % advxx = 2.0_rprec*( b(i,j,k)%u_w*Cududx + b(i,j,k)%v_w*Cududy + b(i,j,k)%w*Cududz)
!c(i,j,k) % advyy = 2.0_rprec*( b(i,j,k)%u_w*Cvdvdx + b(i,j,k)%v_w*Cvdvdy + b(i,j,k)%w*Cvdvdz)
!c(i,j,k) % advzz = 2.0_rprec*( b(i,j,k)%u_w*Cwdwdx + b(i,j,k)%v_w*Cwdwdy + b(i,j,k)%w*Cwdwdz)
!c(i,j,k) % advxy = (b(i,j,k)%u_w*Cudvdx + b(i,j,k)%v_w*Cudvdy + b(i,j,k)%w*Cudvdz)            &
!    + (b(i,j,k)%u_w*Cvdudx + b(i,j,k)%v_w*Cvdudy + b(i,j,k)%w*Cvdudz)
!c(i,j,k) % advxz = (b(i,j,k)%u_w*Cudwdx + b(i,j,k)%v_w*Cudwdy + b(i,j,k)%w*Cudwdz)            &
!    + (b(i,j,k)%u_w*Cwdudx + b(i,j,k)%v_w*Cwdudy + b(i,j,k)%w*Cwdudz)
!c(i,j,k) % advyz = (b(i,j,k)%u_w*Cvdwdx + b(i,j,k)%v_w*Cvdwdy + b(i,j,k)%w*Cvdwdz)            &
!    + (b(i,j,k)%u_w*Cwdvdx + b(i,j,k)%v_w*Cwdvdy + b(i,j,k)%w*Cwdvdz)
!
!! Transport by fluctuations, d(ui*uj*uk)dxk = ui*uk*dujdxk + uj*uk*duidxk
!c(i,j,k) % tflucxx = 2.0_rprec*(Cuududx + Cuvdudy + Cuwdudz)
!c(i,j,k) % tflucyy = 2.0_rprec*(Cvudvdx + Cvvdvdy + Cvwdvdz)
!c(i,j,k) % tfluczz = 2.0_rprec*(Cwudwdx + Cwvdwdy + Cwwdwdz)
!c(i,j,k) % tflucxy = (Cuudvdx + Cuvdvdy + Cuwdvdz) + (Cvududx + Cvvdudy + Cvwdudz)
!c(i,j,k) % tflucxz = (Cuudwdx + Cuvdwdy + Cuwdwdz) + (Cwududx + Cwvdudy + Cwwdudz)
!c(i,j,k) % tflucyz = (Cvudwdx + Cvvdwdy + Cvwdwdz) + (Cwudvdx + Cwvdvdy + Cwwdvdz)
!
!! Transport by pressure-velocity, d(p*(uj*Iik+ui*Ijk))dxk = d(p*uj)dxi + d(p*ui)dxj
!c(i,j,k) % tpresxx = 2.0_rprec*(Cpdudx + Cudpdx)
!c(i,j,k) % tpresyy = 2.0_rprec*(Cpdvdy + Cvdpdy)
!c(i,j,k) % tpreszz = 2.0_rprec*(Cpdwdz + Cwdpdz)
!c(i,j,k) % tpresxy = (Cpdvdx + Cvdpdx) + (Cpdudy + Cudpdy)
!c(i,j,k) % tpresxz = (Cpdwdx + Cwdpdx) + (Cpdudz + Cudpdz)
!c(i,j,k) % tpresyz = (Cpdwdy + Cwdpdy) + (Cpdvdz + Cvdpdz)
!
!! Pressure-Strain, 2*p*sij = p*(dujdxi + duidxj)
!c(i,j,k) % pstrainxx = 2.0_rprec*Cpdudx
!c(i,j,k) % pstrainyy = 2.0_rprec*Cpdvdy
!c(i,j,k) % pstrainzz = 2.0_rprec*Cpdwdz
!c(i,j,k) % pstrainxy = Cpdudy + Cpdvdx
!c(i,j,k) % pstrainxz = Cpdudz + Cpdwdx
!c(i,j,k) % pstrainyz = Cpdvdz + Cpdwdy
!
!! Production Rate, -(Rik*dujdxk + Rjk*duidxk)
!c(i,j,k) % prodxx = -2.0_rprec*(Rxx*a(i,j,k)%dudx + Rxy*a(i,j,k)%dudy + Rxz*a(i,j,k)%dudz)
!c(i,j,k) % prodyy = -2.0_rprec*(Rxy*a(i,j,k)%dvdx + Ryy*a(i,j,k)%dvdy + Ryz*a(i,j,k)%dvdz)
!c(i,j,k) % prodzz = -2.0_rprec*(Rxz*a(i,j,k)%dwdx + Ryz*a(i,j,k)%dwdy + Rzz*a(i,j,k)%dwdz)
!c(i,j,k) % prodxy = -(Rxx*a(i,j,k)%dvdx + Rxy*a(i,j,k)%dvdy + Rxz*a(i,j,k)%dvdz +                 &
!    Rxy*a(i,j,k)%dudx + Ryy*a(i,j,k)%dudy + Ryz*a(i,j,k)%dudz)
!c(i,j,k) % prodxz = -(Rxx*a(i,j,k)%dwdx + Rxy*a(i,j,k)%dwdy + Rxz*a(i,j,k)%dwdz +                 &
!    Rxz*a(i,j,k)%dudx + Ryz*a(i,j,k)%dudy + Rzz*a(i,j,k)%dudz)
!c(i,j,k) % prodyz = -(Rxy*a(i,j,k)%dwdx + Ryy*a(i,j,k)%dwdy + Ryz*a(i,j,k)%dwdz +                 &
!    Rxz*a(i,j,k)%dvdx + Ryz*a(i,j,k)%dvdy + Rzz*a(i,j,k)%dvdz)
!
!! Pseudo-dissipation, 2*nu*duidxk*dujdxk
!c(i,j,k) % pdissxx = 2.0_rprec*nu_molec*(Cuxux + Cuyuy + Cuzuz)
!c(i,j,k) % pdissyy = 2.0_rprec*nu_molec*(Cvxvx + Cvyvy + Cvzvz)
!c(i,j,k) % pdisszz = 2.0_rprec*nu_molec*(Cwxwx + Cwywy + Cwzwz)
!c(i,j,k) % pdissxy = 2.0_rprec*nu_molec*(Cuxvx + Cuyvy + Cuzvz)
!c(i,j,k) % pdissxz = 2.0_rprec*nu_molec*(Cuxwx + Cuywy + Cuzwz)
!c(i,j,k) % pdissyz = 2.0_rprec*nu_molec*(Cvxwx + Cvywy + Cvzwz)
!
!! Transport by viscous diffusion
!! Cuilapuj already multiplied by viscosity because divtxj was used
!! Change sign of Cuilapuj because of how divtxj is set up
!c(i,j,k) % tviscxx = c(i,j,k) % pdissxx - 2.0_rprec*Culapu
!c(i,j,k) % tviscyy = c(i,j,k) % pdissyy - 2.0_rprec*Cvlapv
!c(i,j,k) % tvisczz = c(i,j,k) % pdisszz - 2.0_rprec*Cwlapw
!c(i,j,k) % tviscxy = c(i,j,k) % pdissxy - (Culapv + Cvlapu)
!c(i,j,k) % tviscxz = c(i,j,k) % pdissxz - (Culapw + Cwlapu)
!c(i,j,k) % tviscyz = c(i,j,k) % pdissyz - (Cvlapw + Cwlapv)
!
!! ----------------------- Mean kinetic energy balance --------------------------
!
!end do 
!end do
!end do
!
!! ----------------------------- Prepare Output ---------------------------------
!! Move all terms to the RHS
!c % advxx = - c % advxx
!c % advyy = - c % advyy
!c % advzz = - c % advzz
!c % advxy = - c % advxy
!c % advxz = - c % advxz
!c % advyz = - c % advyz
!
!c % tflucxx = - c % tflucxx
!c % tflucyy = - c % tflucyy
!c % tfluczz = - c % tfluczz
!c % tflucxy = - c % tflucxy
!c % tflucxz = - c % tflucxz
!c % tflucyz = - c % tflucyz
!
!c % tpresxx = - c % tpresxx
!c % tpresyy = - c % tpresyy
!c % tpreszz = - c % tpreszz
!c % tpresxy = - c % tpresxy
!c % tpresxz = - c % tpresxz
!c % tpresyz = - c % tpresyz
!
!c % pdissxx = - c % pdissxx
!c % pdissyy = - c % pdissyy
!c % pdisszz = - c % pdisszz
!c % pdissxy = - c % pdissxy
!c % pdissxz = - c % pdissxz
!c % pdissyz = - c % pdissyz
!
!end function budget_compute
!#endif


end module stat_defs
