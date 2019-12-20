! c = - (u X vort)
!...Compute rotational convective term in physical space  (X-mom eq.)
!...For staggered grid
!...u2, du2, and du4 are on W nodes, rest on UVP nodes
!...(Fixed k=Nz plane on 1/24)
!...Corrected near wall on 4/14/96 {w(DZ/2)=0.5*w(DZ)}
! sc: better to have 0.25*w(dz), since really parabolic.
!--provides cx, cy, cz at 1:nz-1
subroutine convec (cx,cy,cz)
use types,only:rprec
use param
use sim_param, only : u1=>u, u2=>v, u3=>w, du1d2=>dudy, du1d3=>dudz,   &
                      du2d1=>dvdx, du2d3=>dvdz, du3d1=>dwdx, du3d2=>dwdy
use fft
use debug_mod
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real (rprec), dimension (ld, ny, $lbz:nz), intent (out) :: cx, cy, cz

logical, parameter :: DEBUG = .false.

integer::jz
integer :: jz_min, jz_max

!--save forces heap storage
real(kind=rprec), save, dimension(ld_big,ny2,nz)::cc_big
!real(kind=rprec),dimension(ld_big,ny2,nz)::cc_big
!--save forces heap storage
real (rprec), save, dimension (ld_big, ny2, $lbz:nz) :: u1_big, u2_big, u3_big
!real (rprec), dimension (ld_big, ny2, $lbz:nz) :: u1_big, u2_big, u3_big
!--MPI: only u1_big(0:nz-1), u2_big(0:nz-1), u3_big(1:nz) are used
!--save forces heap storage 
real (rprec), save, dimension (ld_big, ny2, nz) :: vort1_big, vort2_big,  &
                                                   vort3_big
!real (rprec), dimension (ld_big, ny2, nz) :: vort1_big, vort2_big, vort3_big
!--MPI: only vort1_big(1:nz), vort2_big(1:nz), vort3_big(1:nz-1) are used
real(kind=rprec)::ignore_me,const,factor

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
real(kind=rprec),dimension($lbz:nz) :: ust
real(rprec) :: z
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (VERBOSE) write (*, *) 'started convec'

!...Recall dudz, and dvdz are on UVP node for k=1 only
!...So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! sc: it would be nice to NOT have to loop through slices here
! Loop through horizontal slices

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
!DY Start here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(OCEAN_FLAG .and. STOKES_FLAG) then
   do jz=$lbz,nz
      $if ($MPI)
      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
      $else
      z=(real(jz)-0.5_rprec)*dz
      $endif
      ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z)
   enddo
endif
!+++++++++++++++
!DY End here
!+++++++++++++++

const=1._rprec/(nx*ny)
!$omp parallel do default(shared) private(jz)		
do jz = $lbz, nz
  !--MPI: u1_big, u2_big needed at jz = 0, u3_big not needed though
  !--MPI: could get u{1,2}_big
! use cx,cy,cz for temp storage here! 
   cx(:,:,jz)=const*u1(:,:,jz)
   cy(:,:,jz)=const*u2(:,:,jz)
   cz(:,:,jz)=const*u3(:,:,jz)
! do forward fft on normal-size arrays
   call rfftwnd_f77_one_real_to_complex(forw,cx(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,cy(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,cz(:,:,jz),ignore_me)
! zero pad: padd takes care of the oddballs
   call padd(u1_big(:,:,jz),cx(:,:,jz))
   call padd(u2_big(:,:,jz),cy(:,:,jz))
   call padd(u3_big(:,:,jz),cz(:,:,jz))
! Back to physical space
! the normalization should be ok...
   call rfftwnd_f77_one_complex_to_real(back_big,u1_big(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back_big,u2_big(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back_big,u3_big(:,:,jz),ignore_me)
end do

do jz = 1, nz
! now do the same, but with the vorticity!
!--if du1d3, du2d3 are on u-nodes for jz=1, then we need a special
!  definition of the vorticity in that case which also interpolates
!  du3d1, du3d2 to the u-node at jz=1
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
        (jz == 1) ) then
     !--du3d2(jz=1) should be 0, so we could use this
     cx(:, :, 1) = const * ( 0.5_rprec * (du3d2(:, :, 1) +     &
                                           du3d2(:, :, 2)) -   &
                              du2d3(:, :, 1) )
     !--du3d1(jz=1) should be 0, so we could use this
     cy(:, :, 1) = const * ( du1d3(:, :, 1) -                &
                              0.5_rprec * (du3d1(:, :, 1) +  &
                                           du3d1(:, :, 2)) )
   else
     cx(:,:,jz)=const*(du3d2(:,:,jz)-du2d3(:,:,jz))
     cy(:,:,jz)=const*(du1d3(:,:,jz)-du3d1(:,:,jz))
   end if
   cz(:,:,jz)=const*(du2d1(:,:,jz)-du1d2(:,:,jz))

   call rfftwnd_f77_one_real_to_complex(forw,cx(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,cy(:,:,jz),ignore_me)
   call rfftwnd_f77_one_real_to_complex(forw,cz(:,:,jz),ignore_me)
   call padd(vort1_big(:,:,jz),cx(:,:,jz))
   call padd(vort2_big(:,:,jz),cy(:,:,jz))
   call padd(vort3_big(:,:,jz),cz(:,:,jz))

! Back to physical space
! the normalization should be ok...
   call rfftwnd_f77_one_complex_to_real(back_big,vort1_big(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back_big,vort2_big(:,:,jz),ignore_me)
   call rfftwnd_f77_one_complex_to_real(back_big,vort3_big(:,:,jz),ignore_me)
end do
!$omp end parallel do

! CX
! redefinition of const
const=1._rprec/(nx2*ny2)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! the cc's contain the normalization factor for the upcoming fft's
  cc_big(:,:,1)=const*(u2_big(:,:,1)*(-vort3_big(:,:,1))&
       +0.5_rprec*u3_big(:,:,2)*(vort2_big(:,:,2)))
  !--vort2(jz=1) is located on uvp-node        ^  try with 1 (experimental)
  !--the 0.5 * u3(:,:,2) is the interpolation of u3 to the first uvp node
  !  above the wall (could arguably be 0.25 * u3(:,:,2))

  jz_min = 2
else
  jz_min = 1
end if

!$omp parallel do default(shared) private(jz)
do jz=jz_min,nz-1
   cc_big(:,:,jz)=const*(u2_big(:,:,jz)*(-vort3_big(:,:,jz))&
        +0.5_rprec*(u3_big(:,:,jz+1)*(vort2_big(:,:,jz+1))&
        +u3_big(:,:,jz)*(vort2_big(:,:,jz))))
end do
!$omp end parallel do

! Loop through horizontal slices
!$omp parallel do default(shared) private(jz)	
do jz=1,nz-1
   call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),ignore_me)
! un-zero pad
! note: cc_big is going into cx!!!!
   call unpadd(cx(:,:,jz),cc_big(:,:,jz))
! Back to physical space
   call rfftwnd_f77_one_complex_to_real(back,cx(:,:,jz),ignore_me)
end do
!$omp end parallel do

! CY
! const should be 1./(nx2*ny2) here

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! the cc's contain the normalization factor for the upcoming fft's
!DY  cc_big(:,:,1)=const*(u1_big(:,:,1)*(vort3_big(:,:,1))&
!DY       +0.5_rprec*u3_big(:,:,2)*(-vort1_big(:,:,2)))
  if(OCEAN_FLAG .and. STOKES_FLAG) then
    cc_big(:,:,1)=const*((u1_big(:,:,1)+ust(1))*(vort3_big(:,:,1))&
         +0.5_rprec*u3_big(:,:,2)*(-vort1_big(:,:,2)))
  else
    cc_big(:,:,1)=const*(u1_big(:,:,1)*(vort3_big(:,:,1))&
         +0.5_rprec*u3_big(:,:,2)*(-vort1_big(:,:,2)))
  endif
  !--vort1(jz=1) is uvp-node                    ^ try with 1 (experimental)
  !--the 0.5 * u3(:,:,2) is the interpolation of u3 to the first uvp node
  !  above the wall (could arguably be 0.25 * u3(:,:,2))

  jz_min = 2
else
  jz_min = 1
end if
!++++++++++++++++
!DY End here
!++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!$omp parallel do default(shared) private(jz)
do jz = jz_min, nz - 1
!DY   cc_big(:,:,jz)=const*(u1_big(:,:,jz)*(vort3_big(:,:,jz))&
!DY        +0.5_rprec*(u3_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
!DY        +u3_big(:,:,jz)*(-vort1_big(:,:,jz))))
  if(OCEAN_FLAG .and. STOKES_FLAG) then
    cc_big(:,:,jz)=const*((u1_big(:,:,jz)+ust(jz))*(vort3_big(:,:,jz))&
        +0.5_rprec*(u3_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +u3_big(:,:,jz)*(-vort1_big(:,:,jz))))
  else
    cc_big(:,:,jz)=const*(u1_big(:,:,jz)*(vort3_big(:,:,jz))&
        +0.5_rprec*(u3_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +u3_big(:,:,jz)*(-vort1_big(:,:,jz))))
  endif
end do
!$omp end parallel do
!++++++++++++++++
!DY End here
!++++++++++++++++

!$omp parallel do default(shared) private(jz)		
do jz=1,nz-1
   call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),ignore_me)
 ! un-zero pad
! note: cc_big is going into cy!!!!
   call unpadd(cy(:,:,jz),cc_big(:,:,jz))

! Back to physical space
   call rfftwnd_f77_one_complex_to_real(back,cy(:,:,jz),ignore_me)
end do
!$omp end parallel do

! CZ

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! There is no convective acceleration of w at wall or at top.
  cc_big(:,:,1)=0._rprec

  jz_min = 2
else
  jz_min = 1
end if

!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!
!  jz_max = nz - 1
!else
!  jz_max = nz
!end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang for ocean LES with Stokes drift
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!$omp parallel do default(shared) private(jz)
do jz=jz_min, nz - 1
!DY   cc_big(:,:,jz)=const*0.5_rprec*(&
!DY        (u1_big(:,:,jz)+u1_big(:,:,jz-1))*(-vort2_big(:,:,jz))&
!DY        +(u2_big(:,:,jz)+u2_big(:,:,jz-1))*(vort1_big(:,:,jz))&
!DY         )
  if(OCEAN_FLAG .and. STOKES_FLAG) then
    cc_big(:,:,jz)=const*0.5_rprec*(&
        (u1_big(:,:,jz)+u1_big(:,:,jz-1)+ust(jz)+ust(jz-1))*(-vort2_big(:,:,jz))&
        +(u2_big(:,:,jz)+u2_big(:,:,jz-1))*(vort1_big(:,:,jz))&
         )
  else
    cc_big(:,:,jz)=const*0.5_rprec*(&
        (u1_big(:,:,jz)+u1_big(:,:,jz-1))*(-vort2_big(:,:,jz))&
        +(u2_big(:,:,jz)+u2_big(:,:,jz-1))*(vort1_big(:,:,jz))&
         )
  endif
end do
!$omp end parallel do
!++++++++++++++++
!DY End here
!++++++++++++++++

! Loop through horizontal slices
!$omp parallel do default(shared) private(jz)		
do jz=1,nz - 1
   call rfftwnd_f77_one_real_to_complex(forw_big,cc_big(:,:,jz),ignore_me)

! un-zero pad
! note: cc_big is going into cz!!!!
   call unpadd(cz(:,:,jz),cc_big(:,:,jz))

! Back to physical space
   call rfftwnd_f77_one_complex_to_real(back,cz(:,:,jz),ignore_me)
end do
!$omp end parallel do

$if ($MPI)
  cx(:, :, 0) = BOGUS
  cy(:, :, 0) = BOGUS
  cz(: ,:, 0) = BOGUS
$endif

!--top level is not valid
cx(:, :, nz) = BOGUS
cy(:, :, nz) = BOGUS
cz(:, :, nz) = BOGUS

if (DEBUG) call DEBUG_write (cx(:, :, 1:nz), 'convec.z.cx')
if (DEBUG) call DEBUG_write (cy(:, :, 1:nz), 'convec.z.cy')
if (DEBUG) call DEBUG_write (cz(:, :, 1:nz), 'convec.z.cz')

if (VERBOSE) write (*, *) 'finished convec'

end subroutine convec
