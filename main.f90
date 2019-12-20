program main

use types,only:rprec
use param
use sim_param
use io,only:openfiles,output_loop,output_final,inflow_write,avg_stats
!use output_slice, only: output_slice_loop
use fft
use immersedbc
use test_filtermodule
use topbc,only:setsponge,sponge
use bottombc,only:num_patch,avgpatch,ustar_avg
use scalars_module,only:beta_scal,obukhov,theta_all_in_one,RHS_T,RHS_PCon,RHS_Tf,pollen_all_in_one2,scale_us,deposition,Real_dep,P_surf_flux,P_surf_flux_dep,beta_pcon
use scalars_module2,only:patch_or_remote,timestep_conditions
$if ($LVLSET)
  use level_set
$endif
$if ($TREES)
  use trees
$endif
!--just for debugging
use debug_mod

implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

integer,parameter::wbase=1  !--controls the frequency of screen diagnostics
integer :: jx,jy,jz
integer :: ip
integer :: jt_diurnal

logical, parameter :: DEBUG = .false.,DEBUG2 = .false.

character (64) :: fnamev

real(kind=rprec) rmsdivvel,kea,kestor,testu   
real(kind=rprec),dimension(nz)::u_ndim
real(kind=rprec),dimension(ld,ny,nz)::S_hat,Pr_0
real(kind=rprec)::const,tt
real (rprec) :: force,z,a_leaf
real (rprec) :: ug_time_factor,ug_period1,ug_period2,ug_period3,ug_period4
real(kind=rprec),dimension(2)::timestep_vars

$if ($MPI)
  integer :: np, coords(1)
$endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
real(kind=rprec),dimension($lbz:nz) :: ust
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for oil plume with momentum flux
real (rprec) :: w0_avg
integer :: il,ir,jl,jr
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!--Initialize computational environment (MPI or not)
$if ($MPI)

  !--check for consistent preprocessor & param.f90 definitions of 
  !  MPI and $MPI
  if (.not. USE_MPI) then
    write (*, *) 'inconsistent use of USE_MPI and $MPI'
    stop
  end if

  call mpi_init (ierr)
  call mpi_comm_size (MPI_COMM_WORLD, np, ierr)
  call mpi_comm_rank (MPI_COMM_WORLD, global_rank, ierr)

  !--check if run-time number of processes agrees with nproc parameter
  if (np /= nproc) then
    write (*, *) 'runtime number of procs = ', np,  &
                 ' not equal to nproc = ', nproc
    stop
!+++++ Added by DY
!+++++ To be identical to LES_JHU
!  else
!     nproc = np
!+++++ Change by DY end here
  end if

  !--set up a 1d cartesian topology
!+++++ Changed by DY
!+++++ To be identical to LES_JHU
  call mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
                        .true., comm, ierr)
!  call mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
!                        .false., comm, ierr)
!+++++ Change by DY end here
  !--slight problem here for ghost layers:
  !  u-node info needs to be shifted up to proc w/ rank "up",
  !  w-node info needs to be shifted down to proc w/ rank "down"
  call mpi_cart_shift (comm, 0, 1, down, up, ierr)
  call mpi_comm_rank (comm, rank, ierr)
  call mpi_cart_coords (comm, rank, 1, coords, ierr)
  coord = coords(1)  !--use coord (NOT rank) to determine global position

  !--rank->coord and coord->rank conversions
!+++++ Changed by DY
!+++++ To be identical to LES_JHU
!!$  do ip = 0, np-1
!!$    call mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
!!$    call mpi_cart_coords (comm, ip, 1, coord_of_rank(ip), ierr)
!!$  end do
  allocate( rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1) )
  do ip = 0, nproc-1
     call mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
     call mpi_cart_coords (comm, ip, 1, coords, ierr)
     coord_of_rank(ip) = coords(1)
  end do
!+++++ Change by DY end here

  write (*, *) 'Hello! from process with coord = ', coord

  !--set the MPI_RPREC variable
  if (rprec == kind (1.e0)) then
    MPI_RPREC = MPI_REAL
    MPI_CPREC = MPI_COMPLEX
  else if (rprec == kind (1.d0)) then
    MPI_RPREC = MPI_DOUBLE_PRECISION
    MPI_CPREC = MPI_DOUBLE_COMPLEX
  else
    write (*, *) 'error defining MPI_RPREC/MPI_CPREC'
    stop
  end if

$else

  if (nproc /= 1) then
    write (*, *) 'nproc /=1 for non-MPI run is an error'
    stop
  end if
  if (USE_MPI) then
    write (*, *) 'inconsistent use of USE_MPI and $MPI'
    stop
  end if

$endif

tt=0

!--Initialize surface boundary conditions
  call patch_or_remote ()

!--I.C.s
  call initial()

! formulate the fft plans--may want to use FFTW_USE_WISDOM
! initialize the kx,ky arrays
  call init_fft()

! Open output files      
  call openfiles()
  print *,'Starting from time = ',jt_total

  !DY Added by Di Yang for releasing dye at later stage
  if(PCon_FLAG.and.DYE_RELEASE_DELAY) then
     if(jt_total.le.ini_dye) then
        PCon(:,:,:,ip_dye) = 0._rprec
        RHS_PCon(:,:,:,ip_dye) = 0._rprec
        print *, 'Ereasing existing dye concentration'
        print *, 'jt_total,ini_dye=',jt_total,ini_dye
     end if
   end if

!--initialize test filter
!--this is used for lower BC, even if no dynamic model
!TSC standard dynamic or Lagrangian
  call test_filter_init (2._rprec * filter_size, G_test)

if (model == 3 .or. model == 5) then  !--scale dependent dynamic
  call test_filter_init (4._rprec * filter_size, G_test_test)
end if

if (ubc == 1) then
    call setsponge()
!    print *,'sponge value calculated for damping layer'
else
    sponge=0._rprec
end if

!--Print computational environment on screen
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  print *, 'Number of timesteps', nsteps
  print *, 'dt and Cs = ', dt, cs
  print *, 'Nx, Ny, Nz = ', nx, ny, nz
  if (USE_MPI) print *, 'Number of processes = ', nproc
  print *, 'Number of patches = ', num_patch
  !print *, 'sampling stats every ', c_count, ' timesteps'
  !print *, 'writing stats every ', p_count, ' timesteps'
  if (molec) print*, 'molecular viscosity (dimensional) ', nu_molec
end if

!--MPI: u,v,w should be set at jz = 0:nz before getting here, except
!  bottom process which is BOGUS (starts at 1)
! time Loop
do jt=1,nsteps   

  if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num

  ! advance total time
  tt=tt+dt

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for oil plume with non-zero momentum flux
!DY Start here
  if(OCEAN_FLAG.and.U_PLUME_FLAG) then
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
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! save previous time's right-hand-sides for Adams-Bashforth Integration
  !  (In subroutine "STEP" use first order time advancement on first time 
  !  step).
     RHSx_f = RHSx
     RHSy_f = RHSy
     RHSz_f = RHSz

  ! Call obukhov to calculate the MO functions !!
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) call obukhov()

  !--no MPI here yet
  if (use_bldg) then
    call building_interp (u, v, w, .04_rprec, 3)
    call building_interp (dudx, dvdx, dwdx, .04_rprec, 3)
    call building_interp (dudy, dvdy, dwdy, .04_rprec, 3)
  end if

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,B4_filt:u,v,w',jz,u(4,4,jz),v(4,4,jz),w(4,4,jz)
        end do
     end if
  endif
  ! kill oddballs and calculate horizontal derivatives
  !--MPI: all these arrays are 0:nz
  !--MPI: on exit, u,dudx,dudy,v,dvdx,dvdy,w,dwdx,dwdy valid on jz=0:nz,
  !  except on bottom process (0 level set to BOGUS, starts at 1)
  call filt_da (u, dudx, dudy)
  call filt_da (v, dvdx, dvdy)
  call filt_da (w, dwdx, dwdy)

  ! finite differences
  !--MPI: on exit of ddz_uv, have dudz, dvdz at 1:nz, except
  !  bottom process has 2:nz
  call ddz_uv(dudz,u)
  call ddz_uv(dvdz,v)
  !--MPI: on exit of ddz_w, have dwdz at 0:nz-1, except top process
  !  has 0:nz, and bottom process has 1:nz-1
  call ddz_w(dwdz,w)

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_filt:u,dudx,dudy,dudz',jz,u(6,4,jz),dudx(6,4,jz),dudy(6,4,jz),&
                dudz(6,4,jz),coord
        end do
     endif
  end if

  !TS calculate wall stress and calculate derivatives at wall
  if (dns_bc) then
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress_dns ()
     end if
  else
  !TS "impose" wall stress and calculate derivatives at wall
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress ()  !--provides txz, tyz, dudz, dvdz at jz=1
                           !--MPI: bottom process only
     end if
     if(use_bldg) call walldudx_building
  end if

  ! compute turbulent viscosity (const.)
  if (dns_bc .and. molec) then
    call dns_stress(txx,txy,txz,tyy,tyz,tzz)
  else
    !--MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz
    call sgs_stag()
    $if ($TREES)
      call add_tree_BC ()
    $endif
  end if
  if(use_bldg)then
     call wallstress_building(txy,txz,tyz)
     call building_mask(u,v,w)
  endif

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_sgs:txx,tyy,txy,txz,tyz,coord',jz,txx(6,4,jz),tyy(6,4,jz),txy(6,4,jz),&
                txz(6,4,jz),tyz(6,4,jz),coord
        end do
     endif
  end if

  !xx----VK -ADDED FOR SCALARS !! --xxxxxx
  !TS ADD passive_scalar

  if(S_FLAG.and.(jt.GE.SCAL_INIT))  then
    call theta_all_in_one
  else
  beta_scal=0._rprec
  end if
  !xx----VK -ADDED FOR SCALARS !! --xxxxxx

!++++++++++++++++++++
!DY For debugging
!  goto 666
!DY End here
!++++++++++++++++++++

  ! Time advance for pollen
  ! Chamecki - 08/01/2006

  beta_pcon=0._rprec
!DY Initialize beta_pcon to be zero
!DY IF active_pcon=.true. beta_pcon will be calculated in pollen_all_in_one or pollen_all_in_one2

  IF (PCon_FLAG .AND. (jt.GE.PCon_init)) THEN
    
    ! Pseudo-spectral
!!$    IF (PCon_scheme==1) CALL pollen_all_in_one
     IF( PCon_scheme==1) then
        print*, "Pseudo-spectral method for particle, stop!"
        stop
     ENDIF
    
    ! QUICK and SMART
    IF (PCon_scheme==2 .OR. PCon_scheme==3) CALL pollen_all_in_one2

  END IF

!++++++++++++++++++++
!DY For debugging
!666 continue
!DY End here
!++++++++++++++++++++

  $if ($MPI)
     !--exchange ghost-node information for tij
     !--send stuff up to ghost nodes
     !--move this into sgs_stag?
     call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 6,   &
                        tzz(:, :, 0), ld*ny, MPI_RPREC, down, 6,  &
                        comm, status, ierr)
  $endif

  ! compute divergence of SGS shear stresses     
  ! note: the divt's and the diagonal elements of t are equivalenced!
  !--actually, they are not equivalenced in this version

  !--provides divtz 1:nz-1
  call divstress_uv(divtx, txx, txy, txz)
  call divstress_uv(divty, txy, tyy, tyz)
  !--provides divtz 1:nz-1, except 1:nz at top process
  call divstress_w(divtz, txz, tyz, tzz)

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_divstress:tzz,divtx,divty,divtz,coord',jz,tzz(6,4,jz),divtx(6,4,jz),divty(6,4,jz),&
                divtz(6,4,jz),coord
        end do
     endif
  end if

  if (VERBOSE) write (*, *) 'main about to call convec'

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,B4_convec:RHSx,RHSy,RHSz,coord',jz,RHSx(4,4,jz),RHSy(4,4,jz),&
                RHSz(4,4,jz),coord
        end do
     end if
  endif
  !--provides RHS{x,y,z} 1:nz-1
  call convec(RHSx,RHSy,RHSz)

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_convec:RHSx,RHSy,RHSz,coord',jz,RHSx(4,4,jz),RHSy(4,4,jz),&
                RHSz(4,4,jz),coord
        end do
     endif
  end if

  if (use_bldg) call building_mask (u, v, w)

  ! Compute preliminary RHS matrices for pressure calculation
  RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

  dudt(:,:,1:nz-1)=-divtx(:,:,1:nz-1)
  dvdt(:,:,1:nz-1)=-divty(:,:,1:nz-1)
  dwdt(:,:,1:nz-1)=-divtz(:,:,1:nz-1)

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+divT:RHSx,RHSy,RHSz,coord',jz,RHSx(4,4,jz),RHSy(4,4,jz),&
                RHSz(4,4,jz),coord
        end do
     endif
 end if

  if (S_FLAG .and. (.not.passive_scalar)) then
    !--add buoyancy term...only valid for theta
    RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) + beta_scal(:, :, 1:nz-1)
  end if

!+++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY Add buoyancy term due to droplets
  if(active_pcon) then
     RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) + beta_pcon(:, :, 1:nz-1)
  endif
!DY End here
!+++++++++++++++++++++++++++++++++++++++++

  if ((jan_diurnal_run) .and. (ug_diurnal_variation)) then !perform inc. and dec. of geos. wind

       if ((jt_total+1)*dt*z_i/u_star .le. 1000._rprec) then
            jt_diurnal=jt_total+1
       else
            jt_diurnal=jt_total+1-nint(1000._rprec/(dt*z_i/u_star))
       end if
       ! create four ug parameters: ug_period1=ug;ug_period2,ug_period3,ug_period4
       if (nz_tot .ge. 192) then
          ug_period1=ug*u_star;ug_period2=6._rprec
          ug_period3=4._rprec;ug_period4=ug*u_star
       else
          ug_period1=ug*u_star;ug_period2=6._rprec
          ug_period3=6._rprec;ug_period4=ug*u_star
       end if 

       if ((jt_diurnal*dt*z_i/u_star/3600._rprec).GT.1.66_rprec .AND. &
          (jt_diurnal*dt*z_i/u_star/3600._rprec).LE.3.33_rprec) then
           ug_time_factor=(ug_period1+(ug_period2-ug_period1)/(3.33_rprec-1.66_rprec)&
                         *(jt_diurnal*dt*z_i/u_star/3600._rprec - 1.66_rprec))/u_star/ug
       elseif ((jt_diurnal*dt*z_i/u_star/3600._rprec).GT.3.33_rprec .AND. &
              (jt_diurnal*dt*z_i/u_star/3600._rprec).LE.15._rprec) then
           ug_time_factor=(ug_period2+(ug_period3-ug_period2)/(15._rprec-3.33_rprec)&
                         *(jt_diurnal*dt*z_i/u_star/3600._rprec - 3.33_rprec))/u_star/ug
       elseif ((jt_diurnal*dt*z_i/u_star/3600._rprec).GT.15._rprec .AND. &
              (jt_diurnal*dt*z_i/u_star/3600._rprec).LT.17._rprec) then
           ug_time_factor=(ug_period3+(ug_period4-ug_period3)/(17._rprec-15._rprec)&
                         *(jt_diurnal*dt*z_i/u_star/3600._rprec - 15._rprec))/u_star/ug
       else
           ug_time_factor=ug_period4/u_star/ug
       end if
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       if (jt .eq. 1) open(60,file=path//'output/ug_time.out',position="append") 
       write (60,'(i6.6,1x,e12.6)') jt_total+1,ug_time_factor*ug*u_star
       print *,'jt_diurnal,ug_time_factor*ug',jt_diurnal,ug_time_factor*ug*u_star
     end if
  else
       ug_time_factor=1._rprec
       if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
         if (jt .eq. nsteps) print *,'Used const Ug = ',jt,ug_time_factor*ug*u_star
       end if
  end if !end if clause for jan_diurnal_run

  if (coriolis_forcing) then
    ! This is to put in the coriolis forcing using coriol,ug and vg as
    ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
    ! Note that ug and vg are non-dimensional (using u_star in param.f90)
    RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) +                 &
                         coriol * v(:, :, 1:nz-1) - ug_time_factor*coriol * vg
    RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) -                 &
                         coriol * u(:, :, 1:nz-1) + ug_time_factor*coriol * ug

    dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+coriol*v(:,:,1:nz-1)-coriol*vg
    dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)-coriol*u(:,:,1:nz-1)+coriol*ug
  end if

!++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation
  if (OCEAN_FLAG) then
     ! Enforce coriolis forcing for ocean simulation without geostrophic wind vector
     if(STOKES_FLAG) then
        do jz=$lbz,nz
           $if ($MPI)
           z = (coord*(nz-1) + jz - 0.5_rprec) * dz
           $else
           z=(real(jz)-0.5_rprec)*dz
           $endif
           ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z)
        enddo
        do jz=1,nz-1
           do jy=1,ny
              do jx=1,nx
                 RHSx(jx,jy,jz) = RHSx(jx,jy,jz) +coriol * v(jx,jy,jz)
                 RHSy(jx,jy,jz) = RHSy(jx,jy,jz) -coriol * (u(jx,jy,jz)+ust(jz))
                 dudt(jx,jy,jz)=dudt(jx,jy,jz)+coriol*v(jx,jy,jz)
                 dvdt(jx,jy,jz)=dvdt(jx,jy,jz)-coriol*(u(jx,jy,jz)+ust(jz))
              enddo
           enddo
        enddo
     else
        RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) +coriol * v(:, :, 1:nz-1)
        RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) -coriol * u(:, :, 1:nz-1)
        dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+coriol*v(:,:,1:nz-1)
        dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)-coriol*u(:,:,1:nz-1)
     endif
  end if
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++

  !XXXXXX%%%%%  Add damping terms to the momentum RHS %%%XXXXXXXXXXXX
  if (ubc==1 .and. damping_method==2) then !add damping terms to the momentum RHS
      do jz=1,nz-1 
      RHSx(1:nx,1:ny,jz)=RHSx(1:nx,1:ny,jz)-0.5*(sponge(jz)+sponge(jz+1))*&
                   (u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/(nx*ny))
      RHSy(1:nx,1:ny,jz)=RHSy(1:nx,1:ny,jz)-0.5*(sponge(jz)+sponge(jz+1))*&
                   (v(1:nx,1:ny,jz)-sum(v(1:nx,1:ny,jz))/(nx*ny))
      RHSz(1:nx,1:ny,jz)=RHSz(1:nx,1:ny,jz)-0.5*sponge(jz)*&
                   (w(1:nx,1:ny,jz)-sum(w(1:nx,1:ny,jz))/(nx*ny))
      end do
  elseif (ubc==1 .and. damping_method==1) then
      do jz=1,nz-1
      RHSz(1:nx,1:ny,jz)=RHSz(1:nx,1:ny,jz)-sponge(jz)*&
                   w(1:nx,1:ny,jz)
      end do
  end if 
  !XXXXXX%%%%%  Sponge/dampling block ends %%%%%%%%%%%%%%%XXXXXXXXXXXX

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_coriol:RHSx,RHSy,RHSz,coord',jz,RHSx(4,4,jz),RHSy(4,4,jz),&
                RHSz(4,4,jz),coord
        end do
     endif
  end if
   !--calculate u^(*) (intermediate vel field)
   !  at this stage, p, dpdx_i are from previous time step
   !  (assumes old dpdx has NOT been added to RHSx_f, etc)
   !  we add force (mean press forcing) here so that u^(*) is as close
   !  to the final velocity as possible

  if (use_mean_p_force) then
     force = mean_p_force
  else
     force = 0._rprec
  end if

  if ((jt == 1) .and. (.not. initu)) then
    ! if initu, then this is read from the initialization file
    ! else for the first step put RHS_f=RHS
    !--i.e. at first step, take an Euler step
    RHSx_f=RHSx
    RHSy_f=RHSy
    RHSz_f=RHSz
  end if

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+divT:RHSx,RHSx_f,RHSy,RHSy_f,RHSz,RHSz_f,coord',jz,RHSx(4,4,jz),&
                RHSx_f(4,4,jz),RHSy(4,4,jz),RHSy_f(4,4,jz),RHSz(4,4,jz),RHSz_f(4,4,jz),coord
        end do
     endif
  end if

  ! Added for the Ragweed field
  IF (rag06 .AND. PCon_FLAG) force=force*(scale_us**2)

  !--only 1:nz-1 are valid
  IF (use_res_canopy) THEN

    dragx=0._rprec
    dragy=0._rprec
    dragz=0._rprec

    IF (use_field_scale) THEN

      do jz=1,nz-1

        $if ($MPI)
          z = (coord*(nz-1) + jz - 0.5_rprec) * dz
        $else
          z=(jz-.5_rprec)*dz
        $endif

        if (z.le.h_leaf1) then
          a_leaf = a_leaf1
        else if (z.le.h_leaf2) then
          a_leaf = a_leaf2
        else if (z.le.h_leaf3) then
          a_leaf = a_leaf3
        else if (z.le.h_leaf4) then
          a_leaf = a_leaf4
        else if (z.le.h_leaf5) then
          a_leaf = a_leaf5
        else if (z.le.h_leaf6) then
          a_leaf = a_leaf6
        else if (z.le.h_leaf7) then
          a_leaf = a_leaf7
        else if (z.le.h_canopy) then
          a_leaf = a_leaf7*(h_canopy-z)/(h_canopy-h_leaf7)
        end if
        a_leaf = a_leaf*LAIp/LAIw

        if (z.le.h_canopy) then
          do jy=1,ny
          do jx=1,nx
            dragx(jx,jy,jz)=-Cd*a_leaf*u(jx,jy,jz)* &
                            (2._rprec*ke(jx,jy,jz))**0.5_rprec
            dragy(jx,jy,jz)=-Cd*a_leaf*v(jx,jy,jz)* &
                            (2._rprec*ke(jx,jy,jz))**0.5_rprec
            dragz(jx,jy,jz)=-Cd*a_leaf*w(jx,jy,jz)* &
                            (2._rprec*ke(jx,jy,jz))**0.5_rprec
          end do
          end do
        end if

      end do

    END IF

  IF (use_force_angle) THEN
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                           &
                      dt * ( tadv1 * RHSx(:, :, 1:nz-1) +         &
                             tadv2 * RHSx_f(:, :, 1:nz-1) +       &
                             dragx(:,:,1:nz-1) + force*COS(force_angle) )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                             tadv2 * RHSy_f(:, :, 1:nz-1) + &
                             dragy(:,:,1:nz-1) + force*SIN(force_angle) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                             tadv2 * RHSz_f(:, :, 1:nz-1) + &
                             dragz(:,:,1:nz-1) )

    dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+dragx(:,:,1:nz-1)+mean_p_force*COS(force_angle)
    dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)+dragy(:,:,1:nz-1)+mean_p_force*SIN(force_angle)
    dwdt(:,:,1:nz-1)=dwdt(:,:,1:nz-1)+dragz(:,:,1:nz-1)
  ELSE
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                           &
                      dt * ( tadv1 * RHSx(:, :, 1:nz-1) +         &
                             tadv2 * RHSx_f(:, :, 1:nz-1) +       &
                             dragx(:,:,1:nz-1) + force )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                             tadv2 * RHSy_f(:, :, 1:nz-1) + &
                             dragy(:,:,1:nz-1) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                             tadv2 * RHSz_f(:, :, 1:nz-1) + &
                             dragz(:,:,1:nz-1) )

    dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+dragx(:,:,1:nz-1)+mean_p_force
    dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)+dragy(:,:,1:nz-1)
    dwdt(:,:,1:nz-1)=dwdt(:,:,1:nz-1)+dragz(:,:,1:nz-1)
  END IF

  ELSE

  IF (use_force_angle) THEN
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                           &
                      dt * ( tadv1 * RHSx(:, :, 1:nz-1) +         &
                             tadv2 * RHSx_f(:, :, 1:nz-1) + force*COS(force_angle) )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                             tadv2 * RHSy_f(:, :, 1:nz-1) + force*SIN(force_angle) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                             tadv2 * RHSz_f(:, :, 1:nz-1) )

    dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+mean_p_force*COS(force_angle)
    dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)+mean_p_force*SIN(force_angle)
  ELSE
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) +                           &
                      dt * ( tadv1 * RHSx(:, :, 1:nz-1) +         &
                             tadv2 * RHSx_f(:, :, 1:nz-1) + force )
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSy(:, :, 1:nz-1) +  &
                             tadv2 * RHSy_f(:, :, 1:nz-1) )
    w(:, :, 1:nz-1) = w(:, :, 1:nz-1) +                    &
                      dt * ( tadv1 * RHSz(:, :, 1:nz-1) +  &
                             tadv2 * RHSz_f(:, :, 1:nz-1) )

    dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+mean_p_force
  END IF

  END IF

  $if ($MPI)
    !--after this point, u,v,w at jz = 0 are not useful, until updated
    u(:, :, 0) = BOGUS
    v(:, :, 0) = BOGUS
    w(:, :, 0) = BOGUS
  $endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for oil plume with non-zero momentum flux
!DY Start here
!!$    if(OCEAN_FLAG.and.U_PLUME_FLAG) then
!!$       il=xps-(int(b0_plume/dx)+1)
!!$       ir=xps+(int(b0_plume/dx)+1)
!!$       jl=yps-(int(b0_plume/dy)+1)
!!$       jr=yps+(int(b0_plume/dy)+1)
!!$       if(il.lt.1.or.ir.gt.nx.or.jl.lt.1.or.jr.gt.ny) then
!!$          print*, "Bottom oil plume reaches the boundary!"
!!$          print*, "(il,ir,jl,jr)=",il,ir,jl,jr
!!$          stop
!!$       endif
!!$       if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$          do jy=jl,jr
!!$             do jx=il,ir
!!$                w(jx,jy,nz)=-ubm*exp(-(((jx-xps)*dx)**2+((jy-yps)*dy)**2)/b0_plume**2)
!!$             enddo
!!$          enddo
!!$          w0_avg=sum(w(1:nx,1:ny,nz))/float(nx*ny)
!!$          do jy=1,ny
!!$             do jx=1,nx
!!$                w(jx,jy,nz)=w(jx,jy,nz)-w0_avg
!!$             enddo
!!$          enddo
!!$       endif
!!$    endif
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+RHS:u,v,w,coord',jz,u(4,4,jz),&
                v(4,4,jz),w(4,4,jz),coord
        end do
     endif
  end if

  if (DEBUG) write (*, *) 'main: before press_stag'

  !--solve Poisson equation for pressure
  !--do we ever need p itself, or only its gradient? -> probably
  !  do not need to store p
  !call press_stag (p, dpdx, dpdy)
  !--provides p, dpdx, dpdy at 0:nz-1
  call press_stag_array (p, dpdx, dpdy)

  if (DEBUG) write (*, *) 'main: after press_stag'

  !--calculate dpdz here
  !--careful, p is not dimensioned the same as the others
  dpdz(1:nx, 1:ny, 1:nz-1) = (p(1:nx, 1:ny, 1:nz-1) -   &
                              p(1:nx, 1:ny, 0:nz-2)) / dz
  dpdz(:, :, nz) = BOGUS

  !--if really wanted to, could avoid storing pressure gradients
  !  just add them directly to RHS in press_stag
  RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) - dpdx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) - dpdy(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) - dpdz(:, :, 1:nz-1)

  dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)-dpdx(:,:,1:nz-1)
  dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)-dpdy(:,:,1:nz-1)
  dwdt(:,:,1:nz-1)=dwdt(:,:,1:nz-1)-dpdz(:,:,1:nz-1)

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+press_stag:p,dpdx,dpdy,dpdz,RHSx,RHSy,RHSz,coord',jz,p(4,4,jz),&
                dpdx(4,4,jz),dpdy(4,4,jz),dpdz(4,4,jz),RHSx(4,4,jz),RHSy(4,4,jz),RHSz(4,4,jz),coord
        end do
     endif
  end if

  call forcing ()
  
  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+forcing:u,v,w,coord',jz,u(4,4,jz),&
                v(4,4,jz),w(4,4,jz),coord
        end do
     endif
  end if

  !   provides u, v, w at 1:nz 
  call project ()
  
  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_+project:u,v,w,coord',jz,u(4,4,jz),&
                v(4,4,jz),w(4,4,jz),coord
        end do
     endif
  end if

  $if ($MPI)
    !--exchange ghost-node information
    !--send stuff up to ghost layer (0) (for z-derivs)
    !--nz should already be in sync with 1 level: done in project()
    call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)

    call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)

    call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       w(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)
  $endif

  if (DEBUG2) then
     if(coord.eq.0) then
        do jz=0,nz
           print *,'jz,AFTER_sync_endJT:u,v,w,coord',jz,u(4,4,jz),&
                v(4,4,jz),w(4,4,jz),coord
        end do
     endif
  end if
  !--MPI: at this point, have u, v, w at 0:nz

  call energy(kea)

  $if ($MPI)
    if (coord == 0) then
      ke(:,:,0)=0._rprec
    endif
    call mpi_sendrecv (ke(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       ke(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)
    call mpi_sendrecv (ke(1, 1, 1), ld*ny, MPI_RPREC, down, 2,  &
                       ke(1, 1, nz), ld*ny, MPI_RPREC, up, 2,   &
                       comm, status, ierr)
  $endif

  ke_temp=ke

  call filt_da(ke_temp,dkedx,dkedy)
  call ddz_uv(dkedz,ke)

  dudt(:,:,1:nz-1)=dudt(:,:,1:nz-1)+dkedx(:,:,1:nz-1)
  dvdt(:,:,1:nz-1)=dvdt(:,:,1:nz-1)+dkedy(:,:,1:nz-1)
  dwdt(:,:,1:nz-1)=dwdt(:,:,1:nz-1)+dkedz(:,:,1:nz-1)

  $if ($MPI)
    call mpi_sendrecv (dudt(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       dudt(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)

    call mpi_sendrecv (dvdt(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       dvdt(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)

    call mpi_sendrecv (dwdt(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       dwdt(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)
  $endif

  call avg_stats ()  !--only does something once every n_avg_stats steps

  if (modulo (jt, wbase) == 0) then

    call rmsdiv (rmsdivvel)
    call timestep_conditions(timestep_vars(1),timestep_vars(2))
    ! timestep_vars(1) is CFL and timestep_vars(2) is viscous stability limit

    call output_loop ()
    !  call output_slice_loop()

    if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then

      if ((S_FLAG) .or. (coriolis_forcing)) then
        write (6, 7778) wt_s, S_FLAG, patch_flag, remote_flag, &
                        coriolis_forcing, ug*u_star,lbc
        write (6, 7779) jt_total, dt, jt_total*(dt*z_i/u_star),rmsdivvel, kea, &
                      timestep_vars(1),timestep_vars(2)
      else
        write (6, 7777) jt_total, dt, rmsdivvel, kea, timestep_vars(1),timestep_vars(2)

      end if
    end if

  end if
   
7777 format ('jt,dt,divu,ke:',1x,i6.6,3(1x,e9.4))
7778 format ('wt_s(K-m/s),Scalars,patch_flag,remote_flag,&
             &coriolis,Ug(m/s),lbc:',(f7.3,1x,L2,1x,i2,1x,i2,1x,L2,1x,f7.3,1x,i2))
7779 format ('jt,dt,time(s),divu,ke,CFL,visc_stab:',1x,i6.6,6(1x,e9.4))


  if (write_inflow_file) call inflow_write ()
                              !--for creating inflow_BC file

  if (DEBUG) write (*, *) $str($context_doc), ' reached line ', $line_num
  
  
!  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  ! Added this for smaller output for video
!   IF (output_video) THEN
!    IF ((jt_total>=video_start) .AND. (jt_total<video_end) .AND. modulo(jt_total,video_freq)==0) THEN
!      WRITE (fnamev, '(A,I6.6,A)') path // 'output/vel_short', jt_total, '.out'
!      OPEN(123,file=fnamev,form='UNFORMATTED')
!      WRITE(123)REAL(u(:, :, 1:nz),4),REAL(v(:, :, 1:nz),4),REAL(w(:, :, 1:nz),4), &
!                REAL(PCon(:,:,1:nz),4),REAL(deposition(:,:),4),REAL(Real_dep(:,:),4),&
!                REAL(ustar_avg(:,:),4)

       ! Only 2D outputs for now
!       WRITE(123)REAL(u(:,:,1),4),REAL(v(:,:,1),4),REAL(w(:,:,1),4),REAL(PCon(:,:,1),4), &
!                 REAL(u(:,ny/2,1:nz),4),REAL(v(:,ny/2,1:nz),4),REAL(w(:,ny/2,1:nz),4),REAL(PCon(:,ny/2,1:nz),4), &
!                 REAL(deposition(:,:),4),REAL(Real_dep(:,:),4),REAL(ustar_avg(:,:),4), &
!                 REAL(P_surf_flux(:,:),4),REAL(P_surf_flux_dep(:,:),4)
!       WRITE(123)REAL(u(:,:,1:nz),4),REAL(v(:,:,1:nz),4),REAL(w(:,:,1:nz),4),REAL(PCon(:,:,1:nz),4), &
!                 REAL(deposition(:,:),4),REAL(Real_dep(:,:),4),REAL(ustar_avg(:,:),4), &
!                 REAL(P_surf_flux(:,:),4),REAL(P_surf_flux_dep(:,:),4)
!     CLOSE(123)
!    END IF
!   END IF
!  endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for oil plume with non-zero momentum flux
!DY Start here
!!$    if(OCEAN_FLAG.and.U_PLUME_FLAG) then
!!$       il=xps-(int(b0_plume/dx)+1)
!!$       ir=xps+(int(b0_plume/dx)+1)
!!$       jl=yps-(int(b0_plume/dy)+1)
!!$       jr=yps+(int(b0_plume/dy)+1)
!!$       if(il.lt.1.or.ir.gt.nx.or.jl.lt.1.or.jr.gt.ny) then
!!$          print*, "Bottom oil plume reaches the boundary!"
!!$          print*, "(il,ir,jl,jr)=",il,ir,jl,jr
!!$          stop
!!$       endif
!!$       if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$          do jy=jl,jr
!!$             do jx=il,ir
!!$                w(jx,jy,nz)=-ubm*exp(-(((jx-xps)*dx)**2+((jy-yps)*dy)**2)/b0_plume**2)
!!$                print*, jx,jy,w(jx,jy,nz)
!!$             enddo
!!$          enddo
!!$          w0_avg=sum(w(1:nx,1:ny,nz))/float(nx*ny)
!!$          print*, "w0_avg=",w0_avg
!!$          do jy=1,ny
!!$             do jx=1,nx
!!$                w(jx,jy,nz)=w(jx,jy,nz)-w0_avg
!!$             enddo
!!$          enddo
!!$          open(2000000+jt_total)
!!$          write(2000000+jt_total,*) 'variables=x,y,w'
!!$          write(2000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' f=point'
!!$          do jy=1,ny
!!$             do jx=1,nx
!!$                write(2000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,w(jx,jy,nz)
!!$             enddo
!!$          enddo
!!$505       format(4e12.4)
!!$          close(2000000+jt_total)
!!$       endif
!!$    endif
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end do  !--end time loop

if (output) call output_final (jt)

$if ($MPI)
  call mpi_finalize (ierr)
$endif

end program main
