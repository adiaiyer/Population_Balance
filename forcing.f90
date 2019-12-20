! modified this so that is just calculated the force--it does not do the
! time advancement
subroutine forcing ()
!subroutine forcing(jt)
use types,only:rprec
use param
use sim_param
use immersedbc
$if ($TREES)
  use trees
$endif
$if ($LVLSET)
  use level_set
$endif
implicit none

!integer,intent(in)::jt
integer::px,py,lx,ly,lz
integer :: jx,jy,jz,i,isize,jsize

real (rprec) :: Rx, Ry, Rz 

! start calculation of body forces (fx,fy,fz)
! 'force' is the mean pressure gradient
if (use_bldg) then
   do i=1,n_bldg
     px=bldg_pts(1,i)
     py=bldg_pts(2,i)
     lx=bldg_pts(3,i)
     ly=bldg_pts(4,i)
     lz=bldg_pts(5,i)
     do jz=1,lz
     do jy=py,py+ly
     do jx=px,px+lx

       ! forces after pressure update
       Rx = -tadv1*dpdx(jx,jy,jz)
       Ry = -tadv1*dpdy(jx,jy,jz)
       Rz = -tadv1*dpdz(jx,jy,jz)

       fx(jx,jy,jz) = ((u_des(jx,jy,jz)-u(jx,jy,jz))/dt - Rx)
       fy(jx,jy,jz) = ((v_des(jx,jy,jz)-v(jx,jy,jz))/dt - Ry)
       fz(jx,jy,jz) = ((w_des(jx,jy,jz)-w(jx,jy,jz))/dt - Rz)

     end do
     end do
     end do
   end do
   ! end calculation of forces
endif

$if ($TREES)
  call add_tree_force ()
$endif

$if ($LVLSET)
  call add_level_set_force ()
$endif

end subroutine forcing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_cond ()
use types, only : rprec
use param, only : face_avg, jx_s, nx, ny, nz, pi, read_inflow_file,  &
                  passive_scalar, x_relax
use sim_param, only : u, v,w, theta
use io, only : inflow_read
implicit none

integer :: jx, jy, jz
integer :: i, isize, jsize
integer :: jx_0, jx_s1

real (rprec) :: factor

!---------------------------------------------------------------------

!--read from file
if (read_inflow_file) then  !--read vel inflow @ jx = jx_s from file
  call inflow_read ()  !--this sets u, v, w at (jx_s,:,:)
else
  u(jx_s, :, :) = face_avg
  v(jx_s, :, :) = 0._rprec
  w(jx_s, :, :) = 0._rprec
end if

!--just an experiment
jx_s1 = modulo (jx_s + 1 - 1, nx) + 1
u(jx_s1, :, :) = u(jx_s, :, :)
v(jx_s1, :, :) = v(jx_s, :, :)
w(jx_s1, :, :) = w(jx_s, :, :)

jx_0 = jx_s - x_relax  !--velocity at this point is not forced here 
  
do i = 1, x_relax-1

  jx = jx_s - i

  factor = 0.5_rprec * (1._rprec - cos (pi * real (i, rprec) / x_relax))
  !factor = real (i, rprec) / x_relax

  u(jx, 1:ny, 1:nz) = u(jx_s, 1:ny, 1:nz) +                              &
                      factor * (u(jx_0, 1:ny, 1:nz) - u(jx_s, 1:ny, 1:nz))
  v(jx, 1:ny, 1:nz) = v(jx_s, 1:ny, 1:nz) +                              &
                      factor * (v(jx_0, 1:ny, 1:nz) - v(jx_s, 1:ny, 1:nz))
  w(jx, 1:ny, 1:nz) = w(jx_s, 1:ny, 1:nz) +                              &
                      factor * (w(jx_0, 1:ny, 1:nz) - w(jx_s, 1:ny, 1:nz))

  if (passive_scalar) then
     theta(jx, 1:ny, 1:nz) = 0._rprec +                                  &
                             factor * (theta(jx_0, 1:ny, 1:nz) - 0._rprec)
  end if

end do

end subroutine inflow_cond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--provides u, v, w at 1:nz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine project ()
use param
use sim_param
use immersedbc
use scalars_module, only:inflow_pollen
implicit none

logical, parameter :: DEBUG = .false.

integer :: jx, jy, jz
integer :: jz_min, jz_max

real (rprec) :: RHS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for oil plume with momentum flux
real (rprec) :: w0_avg
integer :: il,ir,jl,jr
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!---------------------------------------------------------------------

do jz = 1, nz - 1
  do jy = 1, ny
    do jx = 1, nx
      RHS = -tadv1 * dpdx(jx, jy, jz)
      u(jx, jy, jz) = (u(jx, jy, jz) + dt * (RHS + fx(jx, jy, jz)))
      RHS = -tadv1 * dpdy(jx, jy, jz)
      v(jx, jy, jz) = (v(jx, jy, jz) + dt * (RHS + fy(jx, jy, jz)))

      !if (DEBUG) then
      !  if ( isnan (u(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in u at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !  if ( isnan (v(jx, jy, jz)) ) then
      !    write (*, *) $str($context_doc)
      !    write (*, *) 'nan in v at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

    end do
  end do
end do

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  jz_min = 2
else
  jz_min = 1
end if

do jz = jz_min, nz - 1
  do jy = 1, ny
    do jx = 1, nx
      RHS = -tadv1 * dpdz(jx, jy, jz)
      w(jx, jy, jz) = (w(jx, jy, jz) + dt * (RHS + fz(jx, jy, jz)))

      !if (DEBUG) then
      !  if ( isnan (w(jx, jy, jz)) ) then
      !    write (*, *) 'nan in w at (jx, jy, jz) = ', jx, jy, jz
      !    stop
      !  end if
      !end if

    end do
  end do
end do

if(inflow)call inflow_cond

! Zero inflow condition for pollen
! Chamecki 08/21/2006
if(inflow_PCon) call inflow_pollen()

!--left this stuff last, so BCs are still enforced, no matter what
!  inflow_cond does
$if ($MPI)
  !--send velocity info down & recv velocity info from above
  call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 1,  &
                     u(1, 1, nz), ld*ny, MPI_RPREC, up, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 2,  &
                     v(1, 1, nz), ld*ny, MPI_RPREC, up, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 3,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, 3,   &
                     comm, status, ierr)                     
$endif

!--enfore bc at top
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then

  if (force_top_bot .and. inflow) then
    u(:, :, nz) = face_avg
    v(:, :, nz) = 0._rprec
  else if (use_mean_p_force) then
    ! no-stress top
    u(:,:,nz)=u(:,:,nz-1)
    ! no-stress top
    v(:,:,nz)=v(:,:,nz-1)
  end if

  w(:, :, nz)=0._rprec

end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Changed by Di Yang for no-slip seabed for ocean simulation
if(OCEAN_FLAG) then
   ! Use bottom bc for ocean simulation
   ! For ocean simulations, the domain is upside down.
   ! jz=nz_tot-1 is the first grid above the bottom boundary in the physical space.
   if(U_PLUME_FLAG) then
      IF (pointsource .AND. (jt_total>=ini_src .AND. jt_total<end_src)) THEN
         il=xps-(int(b0_plume/dx)+5)
         ir=xps+(int(b0_plume/dx)+5)
         jl=yps-(int(b0_plume/dy)+5)
         jr=yps+(int(b0_plume/dy)+5)
         if(il.lt.1.or.ir.gt.nx.or.jl.lt.1.or.jr.gt.ny) then
            print*, "Bottom oil plume reaches the boundary!"
            print*, "(il,ir,jl,jr)=",il,ir,jl,jr
            stop
         endif
         if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
            do jy=jl,jr
               do jx=il,ir
                  w(jx,jy,nz)=-ubm*exp(-(((jx-xps)*dx)**2+((jy-yps)*dy)**2)/(2._rprec*b0_plume**2)) &
                       *(1._rprec-exp(-100._rprec*dt*(jt_total-ini_src)))
               enddo
            enddo
            w0_avg=sum(w(1:nx,1:ny,nz))/float(nx*ny)
            print*, "w_plume_center=",w(xps,yps,nz)
            print*, "w0_avg=",w0_avg
            do jy=1,ny
               do jx=1,nx
                  w(jx,jy,nz)=w(jx,jy,nz)-w0_avg
               enddo
            enddo
         endif
      ENDIF
   endif
   if(SEABED_FLAG) then
      if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
         u(:,:,nz)=0._rprec
         v(:,:,nz)=0._rprec
         w(:, :, nz)=0._rprec
         print*, "In forcing: no-slip seabed!"
      endif
   endif
endif
!DY End here
!++++++++++++++++++++++++++++++++++++++++++

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! just a test
  if (lbc_mom == 'stress free') then
    if (force_top_bot) then
      u(:, :, 1) = face_avg
      v(:, :, 1) = 0._rprec
    else
      u(:, :, 1) = u(:, :, 2)
      v(:, :, 1) = v(:, :, 2)
    end if
  end if

  w(:, :, 1)=0._rprec

end if

end subroutine project
