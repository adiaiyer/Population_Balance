module scalars_module2
use types,only:rprec
use param !, jt_global => jt  !--rename to avoid name clashes
                            !--could also modify all routines to access jt
                            !  from param module, not argument list
use bottombc !makes obukhov functions available
use sim_param,only:u,v,w,theta,q,path,PCon
use sgsmodule,only: Nu_t
use scalars_module,only: L,wstar,dTdz,dqdz,sgs_t3,sgs_q3 , &
                         DPCondz,sgs_PCon3,res_PCon3,Kc_t,Cs2Sc, &
                         sgs_PCon1
implicit none

!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
real(kind=rprec):: T_s_min, T_s_max,z_o_min,z_o_max
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
! Part II of the scalar files .. also look at scalars_module.f90
!contains subroutines:
! ic_scal --- initialize velocity fields and scalars
! patch_or_remote - initialize surface boundary conditions
! scalar_in - read surface data from an external data file (remote-sensed)
! append_zeros_string - append zeros in front of a string - called by scalar_in (NOT USED !!)
! scalar_slice - outputs time-averaged x-z slices of scalar variables
! controlled by c_count & p_count from param; Uses file unit numbers from (36-47)
! obukhov_slice - outputs the obukhov variables (phi,psi, L,wstar);
! called from scalar_slice and toggled by parameter obukhov_output (specified in scalars_module.f)
! DIURNAL_FORCING - sets the surface temperature equal to the mean temperature from a data file via a 
! diurnal forcing over time changing the surface temperature at each time step
! 
! Authored by Vijayant Kumar
! Last modified - June 16, 2005
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

contains

!!$!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!$subroutine ic_scal()
!!$!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!!$!c...Log profile that is modified to flatten at z=z_i
!!$!c .. Modified by Vijayant to put scalars back
!!$! Last modified April 14, 2004
!!$
!!$!use types,only:rprec
!!$!use sim_param,only:u,v,w,theta,q
!!$!use bottombc
!!$
!!$implicit none
!!$real(kind=rprec),dimension(nz)::ubar,vbar,wbar
!!$real(kind=rprec)::rms, noise, arg, arg2,theta_mean,ustar_ratio
!!$real(kind=rprec)::z,w_star,T_star,q_star,ran3
!!$real(kind=rprec)::z_turb_limit,perturb_height_factor,z_inv
!!$integer::jx,jy,jz,seed,jz_abs
!!$!+++++++++++++++++++++++++++++++++++++
!!$!DY Added by Di Yang for ocean les
!!$real(kind=rprec):: dep_ek
!!$complex(kind=rprec)::gamma,Ubarcomp
!!$!DY End here
!!$!+++++++++++++++++++++++++++++++++++++
!!$
!!$$if ($MPI)
!!$  integer :: sendcounts(nproc)
!!$  integer :: recvcounts(nproc)
!!$  integer :: displs(nproc)
!!$$endif
!!$
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot
!!$
!!$!      if (patch_flag .eq. 1) theta_mean=theta_s1
!!$!      if (patch_flag .eq. 1) theta_mean=T_init
!!$
!!$if (lbc .eq. 0) then
!!$!           theta_mean=T_s_min-2._rprec !T_s_min is dimensional while T_s is non-dimensional
!!$   theta_mean=T_s_min-0.001*T_s_min !T_s_min is dimensional while T_s is non-dimensional
!!$   print *,'theta_mean = ',theta_mean
!!$else
!!$   theta_mean=T_init
!!$   print *,'theta_mean = ',theta_mean
!!$end if
!!$
!!$if (wt_s .lt. 0._rprec) then
!!$   perturb_height_factor=0.10_rprec
!!$   z_inv=0.10_rprec*z_i
!!$else
!!$   perturb_height_factor=0.3_rprec
!!$   z_inv=0.57_rprec*z_i
!!$!++++++++++++++++++++++++++++++++++++++++++++++++
!!$!DY Added by Di Yang for ocean simulation
!!$!DY For z_i=300 m, z_inv=100 m.
!!$   if(OCEAN_FLAG) then
!!$      perturb_height_factor=0.16666667_rprec
!!$      z_inv=0.166666667_rprec*z_i
!!$   endif
!!$!DY End here
!!$!++++++++++++++++++++++++++++++++++++++++++++++++
!!$end if
!!$z_turb_limit=perturb_height_factor*z_i
!!$      
!!$if (wt_s .eq. 0.0_rprec) then
!!$! Compute the values of w_star etc. using the default value of
!!$! wt_s = 0.06
!!$   w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
!!$! w_star is of O(1) with z_i=500 and wt_s=0.06
!!$   T_star=0.06_rprec/w_star
!!$   q_star=T_star
!!$else
!!$   w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
!!$   T_star=wt_s/w_star
!!$   q_star=T_star
!!$end if
!!$!      print *,'w_star,T_star,q_star,seed',w_star,T_star,q_star,seed
!!$        
!!$!      zo_avg=sum(zo)/float(nx*ny)
!!$
!!$$if ($MPI)
!!$print *,'Modified Log Profile for IC for coord = ',coord
!!$$else
!!$print *,'Modified Log Profile for IC'
!!$$endif
!!$do jz=1,nz
!!$   $if ($MPI)
!!$   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$   $else
!!$   z=(real(jz)-0.5_rprec)*dz
!!$   $endif
!!$!        z=(real(jz)-0.5_rprec)*dz*z_i
!!$!c IC in equilibrium with rough surface (rough dominates in effective zo)
!!$!        arg2=z/(sum(zo)/float(nx*ny))
!!$   ustar_ratio = ug*vonk/log(z_inv/zo1)
!!$   arg2=z/(zo1/z_i)
!!$   arg=ustar_ratio*(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
!!$   if (coriolis_forcing) then
!!$      ubar(jz)=ug
!!$      vbar(jz)=vg
!!$      wbar(jz)=0._rprec
!!$! Note that ug and vg have already been non-dimensionalized in param.f90
!!$!        ubar(jz)=arg/30._rprec
!!$   else
!!$      ubar(jz)=arg
!!$      vbar(jz)=0._rprec
!!$      wbar(jz)=0._rprec
!!$   end if
!!$!C sc: I changed around some parenthesis here
!!$!        if (z.gt.(0.6_rprec*z_i)) then
!!$!        if (z.gt.(z_turb_limit)) then
!!$!        print *, 'Statement executed for the scalars'
!!$!        ubar(jz)=ubar(jz-1)
!!$!        end if
!!$   if (z.gt.z_inv/z_i) then
!!$      ubar(jz)=ug
!!$   endif
!!$!++++++++++++++++++++++++++++++++++++++++++++++
!!$!DY Added by Di Yang for ocean simulation
!!$   if(OCEAN_FLAG .and. (.not. STOKES_FLAG)) then
!!$      dep_ek = sqrt(2._rprec*nu_ek/(coriol*u_star/z_i))/z_i 
!!$      ! coriol was normalied in param.f90, nu_ek is dimensional.
!!$!           ubar(jz)=-(1._rprec/vonk)*log(arg2)+(1._rprec/vonk)*log((nz_tot-0.5_rprec)*dz/(zo1/z_i))
!!$!           call Ekman_layer()
!!$      ubar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * cos(-z/dep_ek-pi/4._rprec))
!!$      vbar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * sin(-z/dep_ek-pi/4._rprec))
!!$!      write(100,*) z,ubar(jz),vbar(jz)
!!$!      print*, "nu_ek,coriol,dep_ek=",nu_ek,coriol,dep_ek
!!$   else if(OCEAN_FLAG .and. STOKES_FLAG) then
!!$      gamma = cmplx(0._rprec,1._rprec)*coriol*U_stokes/(4._rprec*wavenm_w**2*(nu_ek/(u_star*z_i)) &
!!$           -cmplx(0._rprec,1._rprec)*coriol)
!!$      Ubarcomp = cmplx(1._rprec,-1._rprec)/sqrt(2._rprec*coriol*(nu_ek/(u_star*z_i))) &
!!$           *(1._rprec-2._rprec*wavenm_w*(nu_ek/(u_star*z_i))*gamma)*exp(cmplx(1._rprec,1._rprec) &
!!$           /sqrt(2._rprec)*sqrt(coriol/(nu_ek/(u_star*z_i)))*(-z))+gamma*exp(2._rprec*wavenm_w*(-z))
!!$      ubar(jz) = real(Ubarcomp)
!!$      vbar(jz) = aimag(Ubarcomp)
!!$   endif
!!$!DY End here
!!$!++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$end do
!!$
!!$!      do jz=1,nz
!!$!       print *,'k, ubar:',jz,ubar(jz)
!!$!      end do
!!$
!!$!+++++++++++++++++++++++++++++++++++++++++
!!$!DY Changed by Di Yang for ocean LES
!!$rms = 3._rprec
!!$!rms = 0.00002_rprec
!!$!DY End here
!!$!+++++++++++++++++++++++++++++++++++++++++
!!$do jz=1,nz
!!$   $if ($MPI)
!!$   jz_abs = coord * (nz-1) + jz
!!$   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
!!$   $else
!!$   jz_abs = jz
!!$   z = (jz-.5_rprec) * dz * z_i
!!$   $endif
!!$   seed = -80 - jz_abs  !--trying to make consistent init for MPI
!!$   do jy=1,ny
!!$      do jx=1,nx
!!$!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!!$!c...Taking std dev of vel as 1 at all heights
!!$
!!$!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!!$!c u should also put L_z=z_i i.e. the inversion layer height should
!!$!c be equal to the height of the domain and in that case the second
!!$!c part of the subsequent if loop will never execute. This is
!!$!c ensured by putting an OR statement in the if clause, which makes 
!!$!c sure that only the first part of if block is executed and not the
!!$!c block after else
!!$
!!$!            z=(real(jz)-0.5_rprec)*dz*z_i
!!$         if (z .LE. z_turb_limit) then
!!$!       if (z .LE. z_i) then
!!$            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$            u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
!!$            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$            v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
!!$            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$            w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
!!$            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$            theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
!!$            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$            q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
!!$         else
!!$            u(jx,jy,jz)=ubar(jz)
!!$            v(jx,jy,jz)=vbar(jz)
!!$            w(jx,jy,jz)=wbar(jz)
!!$!         if (wt_s .ge. 0._rprec) then
!!$            if ((z .gt. z_turb_limit) .and. (z .le. z_inv)) then
!!$               theta(jx,jy,jz)=theta_mean/T_scale
!!$            else
!!$               theta(jx,jy,jz)=(theta_mean+(z-z_inv)*inv_strength)/T_scale
!!$!                  theta(jx,jy,jz)=(theta_mean+(z-z_turb_limit)*inv_strength)/T_scale
!!$            end if
!!$!         else
!!$!         theta(jx,jy,jz)=(theta_mean+(z-z_turb_limit)*inv_strength)/T_scale
!!$!         end if 
!!$            q(jx,jy,jz)=q_mix
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         u(jx,jy,jz)=noise*w_star/u_star*0.01+ubar(jz)
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         v(jx,jy,jz)=noise*w_star/u_star*0.01+vbar(jz)
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         w(jx,jy,jz)=noise*w_star/u_star*0.01
!!$!         theta(jx,jy,jz)=(theta_mean+(z-z_i)*inv_strength)/T_scale
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         q(jx,jy,jz)=q_mix+10._rprec*noise*(1-z/z_i)*q_star
!!$         end if
!!$      end do
!!$   end do
!!$end do
!!$
!!$  !...BC for W
!!$if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!$   w(1:nx, 1:ny, 1) = 0._rprec
!!$end if
!!$if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$   w(1:nx, 1:ny, nz) = 0._rprec
!!$endif
!!$
!!$  !...BC for U, V & T
!!$if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$   u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
!!$   v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
!!$!    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+inv_strength/T_scale*z_i*dz
!!$   theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
!!$end if
!!$
!!$$if ($MPI)
!!$sendcounts = size (u(:,:,1:nz))
!!$recvcounts = size (u(:,:,1:nz-1))
!!$displs = coord_of_rank * recvcounts
!!$call mpi_gatherv (u(1,1,1), sendcounts, MPI_RPREC,&
!!$     u_tot(1,1,1), sendcounts, displs,       &
!!$     MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$call mpi_gatherv (v(1,1,1), sendcounts, MPI_RPREC,&
!!$     v_tot(1,1,1), sendcounts, displs,       &
!!$     MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$call mpi_gatherv (w(1,1,1), sendcounts, MPI_RPREC,&
!!$     w_tot(1,1,1), sendcounts, displs,       &
!!$     MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
!!$     theta_tot(1,1,1), sendcounts, displs,       &
!!$     MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$$endif
!!$
!!$!VK Display the mean vertical profiles of the initialized variables on the
!!$!screen
!!$if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!$   open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
!!$
!!$   if (.not. USE_MPI) then
!!$      do jz=1,nz
!!$         z = (jz - 0.5_rprec) * dz
!!$         write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$              float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$              (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!!$         write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$              float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$              (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!!$      end do
!!$   else
!!$      do jz=1,nz_tot
!!$         z = (jz - 0.5_rprec) * dz
!!$         write(6,7781) jz,z,(sum(u_tot(:,:,jz))/float(nx*ny))*u_star,(sum(v_tot(:,:,jz))/&
!!$              float(nx*ny))*u_star,(sum(w_tot(:,:,jz))/float(nx*ny))*u_star,&
!!$              (sum(theta_tot(:,:,jz))/float(nx*ny))*T_scale
!!$         write(44,7781) jz,z,(sum(u_tot(:,:,jz))/float(nx*ny))*u_star,(sum(v_tot(:,:,jz))/&
!!$              float(nx*ny))*u_star,(sum(w_tot(:,:,jz))/float(nx*ny))*u_star,&
!!$              (sum(theta_tot(:,:,jz))/float(nx*ny))*T_scale
!!$      end do
!!$   endif
!!$   
!!$   close(44)
!!$7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
!!$endif
!!$
!!$end subroutine ic_scal




!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal_v2()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

!use types,only:rprec
!use sim_param,only:u,v,w,theta,q
!use bottombc

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar,thetabar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean,ustar_ratio
real(kind=rprec)::z,w_star,T_star,q_star,ran3
real(kind=rprec)::z_turb_limit,perturb_height_factor,z_inv
integer::jx,jy,jz,seed,jz_abs
!+++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean les
real(kind=rprec):: dep_ek
complex(kind=rprec)::gamma,Ubarcomp
!DY End here
!+++++++++++++++++++++++++++++++++++++

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif

!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
!!$real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot

real (rprec), dimension(1:nz_tot) :: ubar_tot
real (rprec), dimension(1:nz_tot) :: vbar_tot
real (rprec), dimension(1:nz_tot) :: wbar_tot
real (rprec), dimension(1:nz_tot) :: thetabar_tot

!      if (patch_flag .eq. 1) theta_mean=theta_s1
!      if (patch_flag .eq. 1) theta_mean=T_init

if (lbc .eq. 0) then
!           theta_mean=T_s_min-2._rprec !T_s_min is dimensional while T_s is non-dimensional
   theta_mean=T_s_min-0.001*T_s_min !T_s_min is dimensional while T_s is non-dimensional
   print *,'theta_mean = ',theta_mean
else
   theta_mean=T_init
   print *,'theta_mean = ',theta_mean
end if

if (wt_s .lt. 0._rprec) then
   perturb_height_factor=0.10_rprec
   z_inv=0.10_rprec*z_i
else
   perturb_height_factor=0.3_rprec
   z_inv=0.57_rprec*z_i
!++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation
!DY For z_i=300 m, z_inv=100 m.
   if(OCEAN_FLAG) then
      perturb_height_factor=0.1_rprec
      z_inv=0.1_rprec*z_i
   endif
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++
end if
z_turb_limit=perturb_height_factor*z_i
      
if (wt_s .eq. 0.0_rprec) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
   w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
   T_star=0.06_rprec/w_star
   q_star=T_star
else
   w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
   T_star=wt_s/w_star
   q_star=T_star
end if
!      print *,'w_star,T_star,q_star,seed',w_star,T_star,q_star,seed
        
!      zo_avg=sum(zo)/float(nx*ny)

$if ($MPI)
print *,'Modified Log Profile for IC for coord = ',coord
$else
print *,'Modified Log Profile for IC'
$endif
do jz=1,nz
   $if ($MPI)
   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
   $else
   z=(real(jz)-0.5_rprec)*dz
   $endif
!        z=(real(jz)-0.5_rprec)*dz*z_i
!c IC in equilibrium with rough surface (rough dominates in effective zo)
!        arg2=z/(sum(zo)/float(nx*ny))
   ustar_ratio = ug*vonk/log(z_inv/zo1)
   arg2=z/(zo1/z_i)
   arg=ustar_ratio*(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
   if (coriolis_forcing) then
      ubar(jz)=ug
      vbar(jz)=vg
      wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
!        ubar(jz)=arg/30._rprec
   else
      ubar(jz)=arg
      vbar(jz)=0._rprec
      wbar(jz)=0._rprec
   end if
!C sc: I changed around some parenthesis here
!        if (z.gt.(0.6_rprec*z_i)) then
!        if (z.gt.(z_turb_limit)) then
!        print *, 'Statement executed for the scalars'
!        ubar(jz)=ubar(jz-1)
!        end if
   if (z.gt.z_inv/z_i) then
      ubar(jz)=ug
   endif
!++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation
   if(OCEAN_FLAG .and. (.not. STOKES_FLAG)) then
!!$      if(CROSSFLOW_FLAG.and.coriol<1.e-10_rprec) then
!!$         ubar(jz) = ucross
!!$         vbar(jz) = 0._rprec
!!$      endif
      dep_ek = sqrt(2._rprec*nu_ek/(coriol*u_star/z_i))/z_i
      ! coriol was normalied in param.f90, nu_ek is dimensional.
!           ubar(jz)=-(1._rprec/vonk)*log(arg2)+(1._rprec/vonk)*log((nz_tot-0.5_rprec)*dz/(zo1/z_i))
!           call Ekman_layer()
      ubar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * cos(-z/dep_ek-pi/4._rprec))
      vbar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * sin(-z/dep_ek-pi/4._rprec))
!      write(100,*) z,ubar(jz),vbar(jz)
!      print*, "nu_ek,coriol,dep_ek=",nu_ek,coriol,dep_ek
   else if(OCEAN_FLAG .and. STOKES_FLAG) then
!!$      if(coriol<1.e-10_rprec) then
!!$         ubar(jz) = 0._rprec
!!$         vbar(jz) = 0._rprec
!!$      endif
      gamma = cmplx(0._rprec,1._rprec)*coriol*U_stokes/(4._rprec*wavenm_w**2*(nu_ek/(u_star*z_i)) &
           -cmplx(0._rprec,1._rprec)*coriol)
      Ubarcomp = cmplx(1._rprec,-1._rprec)/sqrt(2._rprec*coriol*(nu_ek/(u_star*z_i))) &
           *(1._rprec-2._rprec*wavenm_w*(nu_ek/(u_star*z_i))*gamma)*exp(cmplx(1._rprec,1._rprec) &
           /sqrt(2._rprec)*sqrt(coriol/(nu_ek/(u_star*z_i)))*(-z))+gamma*exp(2._rprec*wavenm_w*(-z))
      ubar(jz) = real(Ubarcomp)
      vbar(jz) = aimag(Ubarcomp)
   endif
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++

end do

!      do jz=1,nz
!       print *,'k, ubar:',jz,ubar(jz)
!      end do

!+++++++++++++++++++++++++++++++++++++++++
!DY Changed by Di Yang for ocean LES
rms = 0.1_rprec
!rms = 0.00002_rprec
!DY End here
!+++++++++++++++++++++++++++++++++++++++++
do jz=1,nz
   $if ($MPI)
   jz_abs = coord * (nz-1) + jz
   z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
   $else
   jz_abs = jz
   z = (jz-.5_rprec) * dz * z_i
   $endif
   seed = -80 - jz_abs  !--trying to make consistent init for MPI
   do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

!            z=(real(jz)-0.5_rprec)*dz*z_i
         if (z .LE. z_turb_limit) then
!       if (z .LE. z_i) then
            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
            u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
            v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
            w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
            noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
            theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
            q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
         else
            u(jx,jy,jz)=ubar(jz)
            v(jx,jy,jz)=vbar(jz)
            w(jx,jy,jz)=wbar(jz)
!         if (wt_s .ge. 0._rprec) then
            if ((z .gt. z_turb_limit) .and. (z .le. z_inv)) then
               !DY Changed by Di Yang for Deepwater Horizon
               theta(jx,jy,jz)=theta_mean/T_scale
!               theta(jx,jy,jz)=(theta_mean-z_i**2._rprec*inv_strength)/T_scale
            else
               theta(jx,jy,jz)=(theta_mean+(z-z_inv)*inv_strength)/T_scale
               !DY for quadratic stratification
              ! theta(jx,jy,jz)=(theta_mean-(z_i-(z-z_inv))**2._rprec*inv_strength)/T_scale
!                  theta(jx,jy,jz)=(theta_mean+(z-z_turb_limit)*inv_strength)/T_scale
            end if
!         else
!         theta(jx,jy,jz)=(theta_mean+(z-z_turb_limit)*inv_strength)/T_scale
!         end if 
            q(jx,jy,jz)=q_mix
!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!         u(jx,jy,jz)=noise*w_star/u_star*0.01+ubar(jz)
!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!         v(jx,jy,jz)=noise*w_star/u_star*0.01+vbar(jz)
!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!         w(jx,jy,jz)=noise*w_star/u_star*0.01
!         theta(jx,jy,jz)=(theta_mean+(z-z_i)*inv_strength)/T_scale
!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!         q(jx,jy,jz)=q_mix+10._rprec*noise*(1-z/z_i)*q_star
         end if
      end do
   end do
end do

  !...BC for W
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   w(1:nx, 1:ny, 1) = 0._rprec
end if
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   w(1:nx, 1:ny, nz) = 0._rprec
endif

  !...BC for U, V & T
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
   v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
!    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+inv_strength/T_scale*z_i*dz
   theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
end if

$if ($MPI)
sendcounts = size (ubar(1:nz))
recvcounts = size (ubar(1:nz-1))
displs = coord_of_rank * recvcounts

do jz=1,nz
   ubar(jz) = sum(u(:,:,jz))/float(nx*ny)
   vbar(jz) = sum(v(:,:,jz))/float(nx*ny)
   wbar(jz) = sum(w(:,:,jz))/float(nx*ny)
   thetabar(jz) = sum(theta(:,:,jz))/float(nx*ny)
enddo

call mpi_gatherv (ubar(1), sendcounts, MPI_RPREC,&
     ubar_tot(1), sendcounts, displs,       &
     MPI_RPREC, rank_of_coord(0), comm, ierr)
call mpi_gatherv (vbar(1), sendcounts, MPI_RPREC,&
     vbar_tot(1), sendcounts, displs,       &
     MPI_RPREC, rank_of_coord(0), comm, ierr)
call mpi_gatherv (wbar(1), sendcounts, MPI_RPREC,&
     wbar_tot(1), sendcounts, displs,       &
     MPI_RPREC, rank_of_coord(0), comm, ierr)
call mpi_gatherv (thetabar(1), sendcounts, MPI_RPREC,&
     thetabar_tot(1), sendcounts, displs,       &
     MPI_RPREC, rank_of_coord(0), comm, ierr)
$endif

!VK Display the mean vertical profiles of the initialized variables on the
!screen
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")

   if (.not. USE_MPI) then
      do jz=1,nz
         z = (jz - 0.5_rprec) * dz
         write(6,7781) jz,z,ubar(jz)*u_star,vbar(jz)*u_star,wbar(jz)*u_star,&
              thetabar(jz)*T_scale
         write(44,7781) jz,z,ubar(jz)*u_star,vbar(jz)*u_star,wbar(jz)*u_star,&
              thetabar(jz)*T_scale
      end do
   else
      do jz=1,nz_tot
         z = (jz - 0.5_rprec) * dz
         write(6,7781) jz,z,ubar_tot(jz)*u_star,vbar_tot(jz)*u_star,wbar_tot(jz)*u_star,&
              thetabar_tot(jz)*T_scale
         write(44,7781) jz,z,ubar_tot(jz)*u_star,vbar_tot(jz)*u_star,wbar_tot(jz)*u_star,&
              thetabar_tot(jz)*T_scale
      end do
   endif
   
   close(44)
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
endif

end subroutine ic_scal_v2





!!$!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$!DY Star here
!!$subroutine Ekman_layer()
!!$!DY Added by Di Yang
!!$!DY Theoretical solution for laminar Ekman layer in the upper ocean
!!$
!!$use param
!!$
!!$implicit none
!!$integer :: jz
!!$real(kind=rprec):: dep_ek,z
!!$real(kind=rprec),dimension(nz)::ubar,vbar,wbar
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$!Calculate Ekman layer depth
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$
!!$dep_ek = sqrt(2._rprec*nu_ek/coriol)/z_i
!!$
!!$!~end here
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$!Calculate Ekman layer velocity
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$
!!$do jz=1,nz
!!$   $if ($MPI)
!!$   z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$   $else
!!$   z = (real(jz)-0.5_rprec)*dz
!!$   $endif
!!$   ubar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * cos(-z/dep_ek-pi/4._rprec))
!!$   vbar(jz) = sqrt(2._rprec)/coriol/dep_ek * exp(-z/dep_ek) * (1._rprec * sin(-z/dep_ek-pi/4._rprec))
!!$!The stress has been normalized by ustar, and thus equals to 1.
!!$enddo
!!$
!!$!~end here
!!$
!!$end subroutine Ekman_layer





!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine patch_or_remote()
!use param
!use bottombc
!use scalars_module
implicit none
real(kind=rprec):: T_s_init_mean,z_o_init_mean
integer :: i,jy
!z_os already defined in scalars_module
!zo,T_s and q_s defined in bottombc
! April 20, 2004 - so far contains both patch_or_remote
! and scalar_in subroutine

if (patch_flag .eq. 1) then
print *, 'Assigning values to patches'
!call patches(zo,T_s,q_s2d,patch,patchnum)
 call patches()

!z_os already defined in scalars_module
!zo,T_s and q_s defined in bottombc
z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)


! Commented by Ying Pan 03/20/2012
! zo_PCon set in bottomc.f90 for all cases
!  ! Added for pollen
!  ! Chamecki - 08/11/2006
!  IF (rag06 .or. fieldsize .or. fieldsize_stripe) THEN
!    ! for this case it has already been set in bottombc.f90
!  ELSE
!    zo_PCon(:,:)=z_ref/z_i
!  END IF
!
!  PRINT*,'zo_PCon assigned'


! Calculate minimum and maximum values of T_s and zo for use in
! ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if
!c sets temperature field and roughness field at the bottom , x-y plane
!c Added by Vijayant
else if (remote_flag .eq. 1) then
print *, 'Assigning remote-sensed values to the surface'
   call scalar_in() ! Updated T_s and zo loaded from bottombc
   
! Calculate minimum and maximum values of T_s and zo for use in
! ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if

   T_s_init_mean=sum(T_s)/float(nx*ny)
!   z_o_init_mean=sum(zo)/float(nx*ny)
! perform logarithmic average for zo
   z_o_init_mean=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

if (remote_to_patch_flag) then
   call remote_to_patch(T_s,1)
   call remote_to_patch(zo,2)
end if

if (remote_homog_flag .eq. 1) then
print *,'Homogeinizing the remote b.c.s'
   T_s=T_s_init_mean
   zo=z_o_init_mean
end if


!   print *,zo(8:12,8:12),sum(zo)/float(nx*ny)

! Non-dimensionalize T
   T_s=(T_s)/T_scale
   zo=zo/z_i
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!open(unit=56,file=path//'interp_data/remote_temp.txt',status="unknown")
!open(unit=57,file=path//'interp_data/remote_z_m.txt',status="unknown")
!      do jz=1,ny
!          read(56,5168) (T_s(i,jz),i=1,nx)
!          read(57,5168) (zo(i,jz),i=1,nx)
!      enddo
!close(56);close(57)

!print *, T_s(5,5),zo(5,5)
!   T_s_init_mean=sum(T_s)/float(nx*ny)
!   z_o_init_mean=sum(zo)/float(nx*ny)

!if (remote_homog_flag .eq. 1) then
!print *,'Homogeinizing the remote b.c.s'
!   T_s=T_s_init_mean
!   zo=z_o_init_mean
!end if
!print *, T_s(5,5),zo(5,5)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! z_os is the scalar roughness length. Divide momentum
! roughness length by 10 to get scalar roughness length
! data (look in the scalar_in routine above)
z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)


      if (GABLS_test) then
        zo(:,:)=0.1_rprec/z_i
        z_os(:,:)=zo(:,:)
        print *,'setting zo and zo_s for GABLS case = 0.1 m'
      end if
end if

! Write surface bc input to a file
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
open(unit=57,file=path//'output/surface_bc_input.txt',status="unknown",position="append")
do jy=1,ny
write(57,5168) (T_s(i,jy),i=1,nx)
end do
do jy=1,ny
write(57,5168) (zo(i,jy),i=1,nx)
end do
close(57)
end if

5168     format(1400(E14.5))

end subroutine patch_or_remote


!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
subroutine scalar_in()
!subroutine scalar_in(T_s,zo)
!subroutine scalar_in(T_s,zo,suffix2,suffix3)
!c This reads in the scalar input from a interpolated scalar
!c file and assigns it to the x-y plane at z=0 i.e. at the ground
!c This for the moment is used to read the momentum roughness, z0
!c and temperature from the USDA remote-sensed data set. The data 
!c can be interpolated using bilinear/cubic/spline interpolation (MATLAB)
!c for the grid size(nx,ny)
!c Authored by Vijayant Kumar
!c Last modified on April 11th, 2004

!use param
use bottombc,only:T_s,zo !Load the variables from bottombc and update in here
implicit none
integer:: ii,jj
character(len=6):: suffix2,suffix3
      
       write(suffix2,'(i6.6)') nx ! create a string from nx   
!       call append_zeros_string(suffix2,nx) ! Append leading zeros to string
       
      if (coarse_grain_flag) then
      write(suffix3,'(i6.6)') stencil_pts    
!      call append_zeros_string(suffix3,stencil_pts) 

 
open(unit=77,file='../interp_data/coarse_grained/interp_temp_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
open(unit=78,file='../interp_data/coarse_grained/interp_z_m_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
       
     
print *,'interp_temp_cg_'//suffix2(4:6)//'pts_'&
//suffix3(4:6)//'.out loaded from scalar_in.f'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      T_s=T_s+273.15_rprec ! Convert T to Kelvin
      close(77)
      close(78)

       else

!open(unit=77,file='./interp_data/interp_temp_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
!open(unit=78,file='./interp_data/interp_z_m_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')

open(unit=77,file=path//'interp_data/interp_temp_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
open(unit=78,file=path//'interp_data/interp_z_m_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
print *,path//'interp_data/interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'
       

      do jj=1,ny
!          do ii=1,nx
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      print *,'Mean T_s (K) = ',sum(T_s)/float(nx*ny)
      print *,'Mean zo (m) = ',sum(zo)/float(nx*ny)
      close(77)
      close(78)
          T_s=T_s+273.15_rprec ! Convert T to Kelvin
      end if
5169     format(1400(E16.10))
end subroutine scalar_in

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!subroutine diurnal_forcing
!!use scalars_module
!implicit none
!integer:: ind
!real(kind=rprec),dimension(no_days*86400,1):: T_s_vector
!
!if ((patch_flag .EQ. 1) .and. (lbc .eq. 0)) then
!     open(unit=56,file='./interp_data/surf_temp_UTAH_2days.txt')
!     do ind=1,nsteps-SCAL_init+1
!        read(56,*)  T_s_vector(ind,1)
!     end do
!     close(56)
!     T_s=T_s_vector(jt-SCAL_init+1,1)/T_scale
!end if
!
!end subroutine diurnal_forcing

   


!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxx-appends character zeros to a character string--XXXX
!!! NOT NEEDED ANYMORE - May 5, 2004

!!subroutine append_zeros_string(string_in,number_in)
!!integer :: number_in
!!character(len=6) string_in

!!if (number_in<100) then
!!    if (number_in<10) then
!!        string_in(4:6)='00'//string_in(6:6)
!!    else
!!        string_in(4:6)='0'//string_in(5:6)
!!    end if
!!end if

!!end subroutine append_zeros_string

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx-------scalar output subroutine-----XXXXXXXXXXXXXXXXXXXXX

subroutine scalar_slice()
!subroutine scalar_slice(w,t,q,sgs_t3,sgs_q3,dTdz,beta_scal,Nu_t,jt)
!c This is exactly the same like the subroutine avgslice with the
!c only difference being that it averages the scalar variables
!c to find the y-averaged instantaneous x-z slices of variables
!c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
!c It also outputs the average covariance between wt and wq
!use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3 
!use output_slice,only: collocate_MPI_averages
implicit none
!$if ($MPI)
!  integer :: recvcounts(nproc)
!  integer :: displs(nproc)
!$endif
integer:: i, j, k
!!integer, parameter :: nz_tot = (nz - 1) * nproc + 1
real(kind=rprec),dimension(nx,nz-1),save:: atheta,t2,q2,asgs_t3,awt
real(kind=rprec),dimension(nx,nz-1),save:: adTdz,anu_t,t3,var_t
!real(kind=rprec),dimension(nx,nz-1),save:: aq,adqdz,awq,q2,asgs_q3
real(kind=rprec):: ttheta1,tt2,tsgst,twt,tdTdz,arg1,arg2,fr
real(kind=rprec):: tnu_t,tt3
!real(kind=rprec):: tq1,tq2,tsgsq,twq,tdqdz
real(kind=rprec),dimension(:,:),allocatable:: avg_scalar_out

fr=(1._rprec/float(p_count))*float(c_count)

!if (jt .EQ. SCAL_init) then
if (jt .EQ. c_count) then
atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec
anu_t=0._rprec;t3=0._rprec;var_t=0._rprec
!aq=0._rprec;q2=0._rprec;asgs_q3=0._rprec;awq=0._rprec;adqdz=0._rprec
end if

do k=1,nz-1
do i=1,nx
ttheta1=0._rprec;tt2=0._rprec;tsgst=0._rprec;twt=0._rprec;tdTdz=0._rprec
tnu_t=0._rprec;tt3=0._rprec
!tq1=0._rprec;tq2=0._rprec;tsgsq=0._rprec;twq=0._rprec;tdqdz=0._rprec

do j=1,ny  
    ttheta1=ttheta1+theta(i,j,k)
!    tt2=tt2+theta(i,j,k)*theta(i,j,k)
    tsgst=tsgst+sgs_t3(i,j,k)
    tdTdz=tdTdz+dTdz(i,j,k)
    tnu_t=tnu_t+Nu_t(i,j,k)
!    tt3=tt3+theta(i,j,k)*theta(i,j,k)*theta(i,j,k)
     if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec
!         arg2=0._rprec
      else  
         arg1=(theta(i,j,k)+theta(i,j,k-1))/2.
!         arg2=(q(i,j,k)+q(i,j,k-1))/2.
      end if

    twt=twt+w(i,j,k)*arg1

!    tq1=tq1+q(i,j,k)
!    tq2=tq2+q(i,j,k)*q(i,j,k)
!    tsgsq=tsgsq+sgs_q3(i,j,k)
!    tdqdz=tdqdz+dqdz(i,j,k)
!    twq=twq+w(i,j,k)*arg2
end do

var_t(i,k)=var_t(i,k)+fr*sum((theta(1:nx,1:ny,k)-sum(theta(1:nx,1:ny,k))/(nx*ny))**2)/(nx*ny)

atheta(i,k)=atheta(i,k)+(fr)*ttheta1/ny
!t2(i,k)=t2(i,k)+(fr)*tt2/ny
asgs_t3(i,k)=asgs_t3(i,k)+(fr)*tsgst/ny
awt(i,k)=awt(i,k)+(fr)*twt/ny
adTdz(i,k)=adTdz(i,k)+(fr)*tdTdz/ny
anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t/ny
!t3(i,k)=t3(i,k)+(fr)*tt3/ny
! Humidity variables
!aq(i,k)=aq(i,k)+(fr)*tq1/ny
!q2(i,k)=q2(i,k)+(fr)*tq2/ny
!asgs_q3(i,k)=asgs_q3(i,k)+(fr)*tsgsq/ny
!awq(i,k)=awq(i,k)+(fr)*twq/ny
!adqdz(i,k)=adqdz(i,k)+(fr)*tdqdz/ny
! Humidity variables
end do
end do
      
if (mod(jt,p_count)==0) then
  allocate(avg_scalar_out(1:nx,1:nz_tot-1));
  call collocate_MPI_averages_N(atheta,avg_scalar_out,35,'theta')
  call collocate_MPI_averages_N(t2,avg_scalar_out,36,'t2')
  call collocate_MPI_averages_N(asgs_t3,avg_scalar_out,37,'sgs_t3')
  call collocate_MPI_averages_N(awt,avg_scalar_out,38,'wt')
  call collocate_MPI_averages_N(adTdz,avg_scalar_out,39,'dTdz')
!  call collocate_MPI_averages_N(aq,avg_scalar_out,40,'q')
!  call collocate_MPI_averages_N(q2,avg_scalar_out,41,'q2')
!  call collocate_MPI_averages_N(adqdz,avg_scalar_out,42,'dqdz')
!  call collocate_MPI_averages_N(asgs_q3,avg_scalar_out,43,'sgs_q3')
!  call collocate_MPI_averages_N(awq,avg_scalar_out,44,'wq')
  call collocate_MPI_averages_N(anu_t,avg_scalar_out,45,'Nu_t')
  call collocate_MPI_averages_N(t3,avg_scalar_out,46,'t3')
  call collocate_MPI_averages_N(var_t,avg_scalar_out,47,'var_t');deallocate(avg_scalar_out)

atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec;
anu_t=0._rprec;t3=0._rprec;
!aq=0._rprec;q2=0._rprec;asgs_q3=0._rprec;awq=0._rprec;adqdz=0._rprec
var_t=0._rprec
end if

 5168     format(1400(E14.5))
!if (obukhov_output==1) then
!    call obukhov_slice !output obukhov variables ...
!end if

end subroutine scalar_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
!integer, parameter:: lbound_loc=average_dim_num*(nx-nz+1)+nz-1
!integer, parameter:: ubound_loc=average_dim_num*(nz-2)+1
!integer, parameter:: lbound_global=average_dim_num*(nx-nz_tot+1)+nz_tot-1
!integer, parameter:: ubound_global=average_dim_num*(nz_tot-2)+1

character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!!real(kind=rprec),dimension(nx,nz-1)::avg_var_proc
!!real(kind=rprec),dimension(nx,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  

  open(file_ind,file=trim(local_filename),status="unknown",position="append")

      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
           do ind1=1,nx
             if (ABS(avg_var_tot_domain(ind1,ind2)) .LT. TINYS) then
               avg_var_tot_domain(ind1,ind2) = 0._rprec
             endif
           end do
         end do
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         do ind1=1,nz_tot-1
           if (ABS(avg_var_tot_domain(ind1,1)) .LT. TINYS) then
             avg_var_tot_domain(ind1,1) = 0._rprec
           endif
         end do
         write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

 5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
!integer, parameter:: lbound_loc=average_dim_num*(nx-nz+1)+nz-1
!integer, parameter:: ubound_loc=average_dim_num*(nz-2)+1
!integer, parameter:: lbound_global=average_dim_num*(nx-nz_tot+1)+nz_tot-1
!integer, parameter:: ubound_global=average_dim_num*(nz_tot-2)+1

real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!!real(kind=rprec),dimension(nx,nz-1)::avg_var_proc
!!real(kind=rprec),dimension(nx,nz_tot-1)::avg_var_tot_domain
  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

!  open(file_ind,file=path//'output/aver_u.out',status="unknown",position="append")
  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

 5168     format(1400(E14.5))
end subroutine collocate_MPI_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!subroutine obukhov_VIJ_slice()
!use scalars_module,only: L,wstar,sgs_t3,jt
!!use scalars_module,only: L,wstar,ustar_avg,sgs_t3
!use sim_param,only: path
!implicit none
!integer:: i, k
!!integer, parameter:: output_dim_flag=1-average_dim_num/2
!integer, parameter:: dim2_OBU=(ny-1)*(1-average_dim_num/2)+1
!!! The size of the output variable is setup such that when average_dim_num(param.f90)
!!! The size of the output variable is setup such that when average_dim_num(param.f90)
!!! is set to 2(output z,1 arrays), the variables being here originally of size nx,ny
!!! are output as nx,1 whereas if average_dim_num is 1 (output x,z arrays), the variables
!!! here are output as (nx,ny) with the averaging just in time
!real(kind=rprec),dimension(nx,dim2_OBU),save:: aphi_m,apsi_m,aphi_h,&
!                                               apsi_h,aL,awstar,austar,awT_sfc
!real(kind=rprec):: fr
!
!
!fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
!!c ------------------------VK----------------------------------
!if (jt .LE. SCAL_init) then
!aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
!aL=0._rprec;awstar=0._rprec;austar=0._rprec;awT_sfc=0._rprec
!end if
!!c ------------------------VK----------------------------------
!
!if (average_dim_num .eq. 2) then 
!    aphi_m(1:nx,1)=aphi_m(1:nx,1)+(fr)*sum(phi_m,dim=2)
!    apsi_m(1:nx,1)=apsi_m(1:nx,1)+(fr)*sum(psi_m,dim=2)
!    aphi_h(1:nx,1)=aphi_h(1:nx,1)+(fr)*sum(phi_h,dim=2)
!    apsi_h(1:nx,1)=apsi_h(1:nx,1)+(fr)*sum(psi_h,dim=2)
!    aL(1:nx,1)=aL(1:nx,1)+(fr)*sum(L,dim=2)
!    awstar(1:nx,1)=awstar(1:nx,1)+(fr)*sum(wstar,dim=2)
!    austar(1:nx,1)=austar(1:nx,1)+(fr)*sum(ustar_avg,dim=2)
!    awT_sfc(1:nx,1)=awT_sfc(1:nx,1)+(fr)*sum(sgs_t3(1:nx,1:ny,1),dim=2)
!else if (average_dim_num .eq. 1) then
!    aphi_m(1:nx,:)=aphi_m(1:nx,:)+(fr)*phi_m(1:nx,:)
!    apsi_m(1:nx,:)=apsi_m(1:nx,1:ny)+(fr)*psi_m
!    aphi_h(1:nx,1:ny)=aphi_h(1:nx,1:ny)+(fr)*phi_h
!    apsi_h(1:nx,1:ny)=apsi_h(1:nx,1:ny)+(fr)*psi_h
!    aL(1:nx,1:ny)=aL(1:nx,1:ny)+(fr)*L
!    awstar(1:nx,1:ny)=awstar(1:nx,1:ny)+(fr)*wstar
!    austar(1:nx,1:ny)=austar(1:nx,1:ny)+(fr)*ustar_avg
!    awT_sfc(1:nx,1:ny)=awT_sfc(1:nx,1:ny)+(fr)*sgs_t3(1:nx,1:ny,1)
!end if
!if (mod(jt,p_count)==0) then
!!cprint *, 'outputting averaged data to files'
!open(48,file=path//'output/aver_phi_m.out',&
!status="unknown",position="append")
!open(49,file=path//'output/aver_psi_m.out',&
!status="unknown",position="append")
!open(50,file=path//'output/aver_phi_h.out',&
!status="unknown",position="append")
!open(52,file=path//'output/aver_psi_h.out',&
!status="unknown",position="append")
!open(53,file=path//'output/aver_L.out',&
!status="unknown",position="append")
!open(54,file=path//'output/aver_wstar.out',&
!status="unknown",position="append")
!open(55,file=path//'output/aver_ustar.out',&
!status="unknown",position="append")
!open(56,file=path//'output/wT_sfc_slice.out',&
!status="unknown",position="append")
!
!if (average_dim_num .eq. 2) then 
!        write(48,5168) jt*dt,(aphi_m(i,1),i=1,nx)
!        write(49,5168) jt*dt,(apsi_m(i,1),i=1,nx)
!        write(50,5168) jt*dt,(aphi_h(i,1),i=1,nx)
!        write(52,5168) jt*dt,(apsi_h(i,1),i=1,nx)
!        write(53,5168) jt*dt,(aL(i,1),i=1,nx)
!        write(54,5168) jt*dt,(awstar(i,1),i=1,nx)
!        write(55,5168) jt*dt,(austar(i,1),i=1,nx)
!        write(56,5168) jt*dt,(awT_sfc(i,1),i=1,nx)
!elseif (average_dim_num .eq. 1) then
!    do k=1,ny
!        write(48,5168) jt*dt,(aphi_m(i,k),i=1,nx)
!        write(49,5168) jt*dt,(apsi_m(i,k),i=1,nx)
!        write(50,5168) jt*dt,(aphi_h(i,k),i=1,nx)
!        write(52,5168) jt*dt,(apsi_h(i,k),i=1,nx)
!        write(53,5168) jt*dt,(aL(i,k),i=1,nx)
!        write(54,5168) jt*dt,(awstar(i,k),i=1,nx)
!        write(55,5168) jt*dt,(austar(i,k),i=1,nx)
!        write(56,5168) jt*dt,(awT_sfc(i,k),i=1,nx)
!    end do
!end if
!aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
!aL=0._rprec;awstar=0._rprec;austar=0._rprec;awT_sfc=0._rprec
!end if

!close(48);close(49);close(50);close(51);close(52);close(53)
!close(54);close(55);close(56)
! 5168     format(1400(E14.5))

!end subroutine obukhov_VIJ_slice

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!subroutine obukhov_slice
!use scalars_module,only: L,wstar,sgs_t3,jt
!!use scalars_module,only: L,wstar,ustar_avg,sgs_t3
!use sim_param,only: path
!implicit none
!integer:: i, k
!!real(kind=rprec),dimension(nx,nz):: phi_m,psi_m,phi_h,psi_h
!!real(kind=rprec),dimension(nx,nz):: L_fin,wstar_fin
!real(kind=rprec),dimension(nx,nz),save:: aphi_m,apsi_m,aphi_h,apsi_h,aL,awstar,austar,awT_sfc
!real(kind=rprec):: fr
!
!
!fr=(1._rprec/float(p_count))*float(c_count)
!!c ------------------------VK----------------------------------
!if (jt .LE. SCAL_init) then
!aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
!aL=0._rprec;awstar=0._rprec;austar=0._rprec;awT_sfc=0._rprec
!end if
!!c ------------------------VK----------------------------------
!
!do k=1,ny
!do i=1,nx
!aphi_m(i,k)=aphi_m(i,k)+(fr)*phi_m(i,k)
!apsi_m(i,k)=apsi_m(i,k)+(fr)*psi_m(i,k)
!aphi_h(i,k)=aphi_h(i,k)+(fr)*phi_h(i,k)
!apsi_h(i,k)=apsi_h(i,k)+(fr)*psi_h(i,k)
!aL(i,k)=aL(i,k)+(fr)*L(i,k)
!awstar(i,k)=awstar(i,k)+(fr)*wstar(i,k)
!austar(i,k)=austar(i,k)+(fr)*ustar_avg(i,k)
!awT_sfc(i,k)=awT_sfc(i,k)+(fr)*sgs_t3(i,k,1)
!end do
!end do
!
!if (mod(jt,p_count)==0) then
!!cprint *, 'outputting averaged data to files'
!open(48,file=path//'output/aver_phi_m.out',&
!status="unknown",position="append")
!open(49,file=path//'output/aver_psi_m.out',&
!status="unknown",position="append")
!open(50,file=path//'output/aver_phi_h.out',&
!status="unknown",position="append")
!open(52,file=path//'output/aver_psi_h.out',&
!status="unknown",position="append")
!open(53,file=path//'output/aver_L.out',&
!status="unknown",position="append")
!open(54,file=path//'output/aver_wstar.out',&
!status="unknown",position="append")
!open(55,file=path//'output/aver_ustar.out',&
!status="unknown",position="append")
!open(56,file=path//'output/wT_sfc_slice.out',&
!status="unknown",position="append")
!
!do k=1,ny
!write(48,5168) jt*dt,(aphi_m(i,k),i=1,nx)
!write(49,5168) jt*dt,(apsi_m(i,k),i=1,nx)
!write(50,5168) jt*dt,(aphi_h(i,k),i=1,nx)
!write(52,5168) jt*dt,(apsi_h(i,k),i=1,nx)
!write(53,5168) jt*dt,(aL(i,k),i=1,nx)
!write(54,5168) jt*dt,(awstar(i,k),i=1,nx)
!write(55,5168) jt*dt,(austar(i,k),i=1,nx)
!write(56,5168) jt*dt,(awT_sfc(i,k),i=1,nx)
!end do
!
!aphi_m=0._rprec;apsi_m=0._rprec;aphi_h=0._rprec;apsi_h=0._rprec
!aL=0._rprec;awstar=0._rprec;austar=0._rprec;awT_sfc=0._rprec
!end if
!
!close(48);close(49);close(50);close(51);close(52);close(53)
!close(54);close(55);close(56)
! 5168     format(1400(E14.5))
!
!end subroutine obukhov_slice
!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine timestep_conditions(CFL,visc_stab)
! This subroutine computes CFl and viscous stability and is called every wbase timesetps
! from main.f90
implicit none
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
  real(kind=rprec),dimension(4,nproc)::max_var_tot_domain
$endif

real(kind=rprec) :: delta, u_res_max, nu_max
real(kind=rprec),dimension(1,4) :: max_vels
real(kind=rprec),intent(out) :: CFL, visc_stab
 
$if ($MPI)
  recvcounts = size (max_vels)
  displs = coord_of_rank * recvcounts 
  max_vels(1,1)=maxval(u(1:nx,1:ny,1:nz-1))
  max_vels(1,2)=maxval(v(1:nx,1:ny,1:nz-1))
  max_vels(1,3)=maxval(w(1:nx,1:ny,1:nz-1))
  max_vels(1,4)=maxval(Nu_t(1:nx,1:ny,1:nz-1))
!  print *,'coord = ',rank,'max_u = ',max_vels(1,1),'max_v = ',max_vels (1,2),'mxx_w = ',max_vels(1,3),'max_Nu_t = ',max_vels(1,4)
  call mpi_gatherv (max_vels(1,1), size(max_vels), MPI_RPREC,                &
                    max_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
!  if (rank == 0) then
!    print *,'After mpi_gatherv: ','max_u = ',maxval(max_var_tot_domain(1,:)),'max_v = ',maxval(max_var_tot_domain(2,:)),'max_w = ',maxval(max_var_tot_domain(3,:)),'max_Nu_t = ',maxval(max_var_tot_domain(4,:))
!  endif
  u_res_max=sqrt(maxval(max_var_tot_domain(1,:))**2+maxval(max_var_tot_domain(2,:))**2+&
  maxval(max_var_tot_domain(3,:))**2)
  nu_max=maxval(max_var_tot_domain(4,:))
$else
  u_res_max = sqrt(maxval(u(1:nx,1:ny,1:nz-1)**2+v(1:nx,1:ny,1:nz-1)**2+&
  w(1:nx,1:ny,1:nz-1)**2)) ! don't bother with interp here
  nu_max=maxval(Nu_t(1:nx,1:ny,1:nz-1))
$endif  

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    delta = min(dx, dy, dz)
    CFL = u_res_max*dt/delta
    visc_stab = dt*nu_max/delta**2
end if

end subroutine timestep_conditions

!!$!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!$!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!$subroutine ic_scal_GABLS()
!!$!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!!$!c...Log profile that is modified to flatten at z=z_i
!!$!c .. Modified by Vijayant to put scalars back
!!$! Last modified April 14, 2004
!!$
!!$!use types,only:rprec
!!$!use sim_param,only:u,v,w,theta,q
!!$!use bottombc
!!$
!!$implicit none
!!$real(kind=rprec),dimension(nz)::ubar,vbar,wbar
!!$real(kind=rprec)::rms, noise, arg, arg2,theta_mean
!!$real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
!!$integer::jx,jy,jz,seed,jz_abs
!!$
!!$!      if (patch_flag .eq. 1) theta_mean=theta_s1
!!$!      if (patch_flag .eq. 1) theta_mean=T_init
!!$       theta_mean=T_init
!!$       perturb_height_factor=50._rprec/z_i
!!$      
!!$      if (wt_s .eq. 0.0_rprec) then
!!$! Compute the values of w_star etc. using the default value of
!!$! wt_s = 0.06
!!$      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
!!$! w_star is of O(1) with z_i=500 and wt_s=0.06
!!$      T_star=0.06_rprec/w_star
!!$      q_star=T_star
!!$      else
!!$      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
!!$      T_star=wt_s/w_star
!!$      q_star=T_star
!!$      end if
!!$!      print *,'w_star,T_star,q_star,seed',w_star,T_star,q_star,seed
!!$        
!!$!      zo_avg=sum(zo)/float(nx*ny)
!!$
!!$         $if ($MPI)
!!$            print *,'Modified Log Profile for IC for coord = ',coord
!!$         $else
!!$            print *,'Modified Log Profile for IC'
!!$         $endif
!!$       do jz=1,nz
!!$         $if ($MPI)
!!$            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$         $else
!!$            z=(real(jz)-0.5_rprec)*dz
!!$         $endif
!!$!        z=(real(jz)-0.5_rprec)*dz*z_i
!!$!c IC in equilibrium with rough surface (rough dominates in effective zo)
!!$        arg2=z/(sum(zo)/float(nx*ny))
!!$        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
!!$        if (coriolis_forcing) then
!!$        ubar(jz)=ug
!!$        vbar(jz)=vg
!!$        wbar(jz)=0._rprec
!!$! Note that ug and vg have already been non-dimensionalized in param.f90
!!$!        ubar(jz)=arg/30._rprec
!!$        else
!!$        ubar(jz)=arg
!!$        vbar(jz)=0._rprec
!!$        wbar(jz)=0._rprec
!!$        end if
!!$!C sc: I changed around some parenthesis here
!!$        if (z.gt.(perturb_height_factor*z_i)) then
!!$!        print *, 'Statement executed for the scalars'
!!$        ubar(jz)=ubar(jz-1)
!!$        end if
!!$       end do
!!$
!!$
!!$!      do jz=1,nz
!!$!       print *,'k, ubar:',jz,ubar(jz)
!!$!      end do
!!$
!!$  rms = 3._rprec
!!$  do jz=1,nz
!!$  $if ($MPI)
!!$    jz_abs = coord * (nz-1) + jz
!!$    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
!!$  $else
!!$    jz_abs = jz
!!$    z = (jz-.5_rprec) * dz * z_i
!!$  $endif
!!$  seed = -80 - jz_abs  !--trying to make consistent init for MPI
!!$    do jy=1,ny
!!$      do jx=1,nx
!!$!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!!$!c...Taking std dev of vel as 1 at all heights
!!$
!!$!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!!$!c u should also put L_z=z_i i.e. the inversion layer height should
!!$!c be equal to the height of the domain and in that case the second
!!$!c part of the subsequent if loop will never execute. This is
!!$!c ensured by putting an OR statement in the if clause, which makes 
!!$!c sure that only the first part of if block is executed and not the
!!$!c block after else
!!$
!!$!            z=(real(jz)-0.5_rprec)*dz*z_i
!!$       if (z.le.perturb_height_factor*z_i) then
!!$!         u(jx,jy,jz)=ubar(jz)
!!$!         v(jx,jy,jz)=vbar(jz)
!!$!         w(jx,jy,jz)=wbar(jz)
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         if (jz .eq. 1) ubar(1)=0._rprec
!!$         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz) !noise
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
!!$         q(jx,jy,jz)=q_mix
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
!!$       else
!!$         u(jx,jy,jz)=ubar(jz)
!!$         v(jx,jy,jz)=vbar(jz)
!!$         w(jx,jy,jz)=wbar(jz)
!!$         theta(jx,jy,jz)=(theta_mean+(z-perturb_height_factor*z_i)*inv_strength)/T_scale
!!$         q(jx,jy,jz)=q_mix
!!$        end if
!!$       end do
!!$     end do
!!$  end do
!!$
!!$
!!$  !...BC for W
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!$    w(1:nx, 1:ny, 1) = 0._rprec
!!$  end if
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$    w(1:nx, 1:ny, nz) = 0._rprec
!!$  endif
!!$
!!$  !...BC for U, V & T
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
!!$    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
!!$!    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+inv_strength/T_scale*z_i*dz
!!$    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
!!$  end if
!!$
!!$!VK Display the mean vertical profiles of the initialized variables on the
!!$!screen
!!$open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
!!$do jz=1,nz
!!$     $if ($MPI)
!!$       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$     $else
!!$       z = (jz - 0.5_rprec) * dz
!!$     $endif
!!$     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!!$     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!!$end do
!!$close(44)
!!$7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
!!$
!!$
!!$end subroutine ic_scal_GABLS

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine budget_TKE_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To compute the budget of TKE, we need the following terms:
! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
! calculated/outputted in avgslice.f90 and scalar_slice.f90
! So, the rest 4 terms will be computed/outputted here
! Similarly for temperature variance budget, we need
! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from 
! scalar_slice.f90.. so we just need to compute term 2 here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use sim_param,only:path,u,v,w,theta,p,txz,tyz
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:dissip
implicit none
integer::i,j,k
real(kind=rprec),dimension(nx,nz-1),save::awu2,awv2,awT2,auv,awp
real(kind=rprec),dimension(nx,nz-1),save::autau13,avtau23,adissip
real(kind=rprec)::twu2,twv2,twt2,tuv,twp,arg1,arg2,arg3,arg4,arg5,arg6,fr
real(kind=rprec)::arg7
real(kind=rprec)::tutau13,tvtau23,tdissip
real(kind=rprec),dimension(:,:),allocatable::avg_budget_out
real(kind=rprec),dimension(nz-1)::ubar_profile,vbar_profile,Tbar_profile

fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)

do k=1,nz-1
ubar_profile(k)=sum(u(1:nx,1:ny,k))/(nx*ny)
vbar_profile(k)=sum(v(1:nx,1:ny,k))/(nx*ny)
Tbar_profile(k)=sum(theta(1:nx,1:ny,k))/(nx*ny)
end do

do k=1,Nz-1
do i=1,Nx
   twu2=0._rprec;twv2=0._rprec;twT2=0._rprec;tuv=0._rprec;twp=0._rprec;
   tutau13=0._rprec;tvtau23=0._rprec;tdissip=0._rprec

   do j=1,Ny
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg5=0._rprec;arg6=0._rprec
      else  
!         arg1=(u(i,j,k)*u(i,j,k)+u(i,j,k-1)*u(i,j,k-1))/2.
!         arg2=(v(i,j,k)*v(i,j,k)+v(i,j,k-1)*v(i,j,k-1))/2.
!         arg3=(theta(i,j,k)*theta(i,j,k)+theta(i,j,k-1)*theta(i,j,k-1))/2.
         arg4=(p(i,j,k)+p(i,j,k-1))/2.
         arg5=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
         arg6=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
         arg7=(theta(i,j,k)-Tbar_profile(k)+theta(i,j,k-1)-Tbar_profile(k-1))/2.
      end if
      twu2=twu2+w(i,j,k)*arg5*arg5 !directly computes <w(u')^2>
      twv2=twv2+w(i,j,k)*arg6*arg6 !directly computes <w(v')^2>
! Also note that since <w> = 0, there is no need to calculate it as
! <w^3> = <w'^3> and we are outputting <w^3> in avgslice
      twT2=twT2+w(i,j,k)*arg7
      twp=twp+w(i,j,k)*arg4
! <u'v'> is not as simple as <u'w'> since <u> .ne. 0 whereas <w>=0
! therefore, we directly calculate <u'v'> here
      tuv=tuv+(u(i,j,k)-ubar_profile(k))*(v(i,j,k)-vbar_profile(k))
      tutau13=tutau13+arg5*txz(i,j,k) !computes SGS transport of TKE i.e. <u'\tau_{13}>
      tvtau23=tvtau23+arg6*tyz(i,j,k) !computes SGS transport of TKE i.e. <v'\tau_{13}>
      tdissip=tdissip+dissip(i,j,k) ! outputs dissip calculated in sgs stag..
   end do
   awu2(i,k)=awu2(i,k)+(fr)*twu2/Ny
   awv2(i,k)=awv2(i,k)+(fr)*twv2/Ny
   awT2(i,k)=awT2(i,k)+(fr)*twT2/Ny
   awp(i,k)=awp(i,k)+(fr)*twp/Ny
   auv(i,k)=auv(i,k)+(fr)*tuv/Ny
   autau13(i,k)=autau13(i,k)+(fr)*tutau13/Ny
   avtau23(i,k)=avtau23(i,k)+(fr)*tvtau23/Ny
   adissip(i,k)=adissip(i,k)+(fr)*tdissip/Ny
end do
end do

if (mod(jt,p_count)==0) then
  allocate(avg_budget_out(1:nx,1:(nz_tot-1)));
  call collocate_MPI_averages_N(awu2,avg_budget_out,61,'awu2')
  call collocate_MPI_averages_N(awv2,avg_budget_out,62,'awv2')
  call collocate_MPI_averages_N(awT2,avg_budget_out,63,'awT2')
  call collocate_MPI_averages_N(awp,avg_budget_out,64,'awp')
  call collocate_MPI_averages_N(auv,avg_budget_out,65,'auv')
  call collocate_MPI_averages_N(autau13,avg_budget_out,66,'autau13')
  call collocate_MPI_averages_N(avtau23,avg_budget_out,67,'avtau23')
  call collocate_MPI_averages_N(adissip,avg_budget_out,68,'adissip');deallocate(avg_budget_out)
!VK Zero out the outputted averages !!
  awu2=0._rprec;awv2=0._rprec;awT2=0._rprec;awp=0._rprec;auv=0._rprec
  autau13=0._rprec;avtau23=0._rprec;adissip=0._rprec
end if

end subroutine budget_TKE_scalar

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine budget_TKE_scalar
!! To compute the budget of TKE, we need the following terms:
!! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
!! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
!! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
!! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
!! calculated/outputted in avgslice.f90 and scalar_slice.f90
!! So, the rest 4 terms will be computed/outputted here
!! Similarly for temperature variance budget, we need
!! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from 
!! scalar_slice.f90.. so we just need to compute term 2 here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use sim_param,only:path,u,v,w,theta,p,txz,tyz
!use param,only:dz,p_count,c_count,jt
!use sgsmodule,only:dissip
!implicit none
!integer::i,j,k
!real(kind=rprec),dimension(nx,nz-1),save::awu2,awv2,awT2,auv,awp
!real(kind=rprec),dimension(nx,nz-1),save::autau13,avtau23,adissip
!real(kind=rprec)::twu2,twv2,twt2,tuv,twp,arg1,arg2,arg3,arg4,arg5,arg6,fr
!real(kind=rprec)::tutau13,tvtau23,tdissip
!real(kind=rprec),dimension(:,:),allocatable::avg_out
!real(kind=rprec),dimension(nz-1)::ubar_profile,vbar_profile
!
!fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
!
!do k=1,nz-1
!ubar_profile(k)=sum(u(1:nx,1:ny,k))/(nx*ny)
!vbar_profile(k)=sum(v(1:nx,1:ny,k))/(nx*ny)
!end do
!
!do k=1,Nz-1
!do i=1,Nx
!   twu2=0._rprec;twv2=0._rprec;twT2=0._rprec;tuv=0._rprec;twp=0._rprec;
!   tutau13=0._rprec;tvtau23=0._rprec;tdissip=0._rprec
!
!   do j=1,Ny
!      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
!         arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg5=0._rprec;arg6=0._rprec
!      else  
!         arg1=(u(i,j,k)*u(i,j,k)+u(i,j,k-1)*u(i,j,k-1))/2.
!         arg2=(v(i,j,k)*v(i,j,k)+v(i,j,k-1)*v(i,j,k-1))/2.
!         arg3=(theta(i,j,k)*theta(i,j,k)+theta(i,j,k-1)*theta(i,j,k-1))/2.
!         arg4=(p(i,j,k)+p(i,j,k-1))/2.
!         arg5=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
!         arg6=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
!      end if
!      twu2=twu2+w(i,j,k)*arg1
!      twv2=twv2+w(i,j,k)*arg2
!      twT2=twT2+w(i,j,k)*arg3
!      twp=twp+w(i,j,k)*arg4
!      tuv=tuv+u(i,j,k)*v(i,j,k)
!      tutau13=tutau13+arg5*txz(i,j,k)
!      tvtau23=tvtau23+arg6*tyz(i,j,k)
!      tdissip=tdissip+dissip(i,j,k)
!   end do
!   awu2(i,k)=awu2(i,k)+(fr)*twu2/Ny
!   awv2(i,k)=awv2(i,k)+(fr)*twv2/Ny
!   awT2(i,k)=awT2(i,k)+(fr)*twT2/Ny
!   awp(i,k)=awp(i,k)+(fr)*twp/Ny
!   auv(i,k)=auv(i,k)+(fr)*tuv/Ny
!   autau13(i,k)=autau13(i,k)+(fr)*tutau13/Ny
!   avtau23(i,k)=avtau23(i,k)+(fr)*tvtau23/Ny
!   adissip(i,k)=adissip(i,k)+(fr)*tdissip/Ny
!end do
!end do
!
!if (mod(jt,p_count)==0) then
!  allocate(avg_out(1:nx,1:(nz_tot-1)));
!  call collocate_MPI_averages_N(awu2,avg_out,44,'awu2')
!  call collocate_MPI_averages_N(awv2,avg_out,44,'awv2')
!  call collocate_MPI_averages_N(awT2,avg_out,44,'awT2')
!  call collocate_MPI_averages_N(awp,avg_out,44,'awp')
!  call collocate_MPI_averages_N(auv,avg_out,44,'auv')
!  call collocate_MPI_averages_N(autau13,avg_out,44,'autau13')
!  call collocate_MPI_averages_N(avtau23,avg_out,44,'avtau23')
!  call collocate_MPI_averages_N(adissip,avg_out,44,'adissip');deallocate(avg_out)
!!VK Zero out the outputted averages !!
!  awu2=0._rprec;awv2=0._rprec;awT2=0._rprec;awp=0._rprec;auv=0._rprec
!  autau13=0._rprec;avtau23=0._rprec;adissip=0._rprec
!end if
!
!end subroutine budget_TKE_scalar
!
!!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!$subroutine ic_scal_GABLS_diurnal()
!!$!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!!$!c...Log profile that is modified to flatten at z=z_i
!!$!c .. Modified by Vijayant to put scalars back
!!$! Last modified April 14, 2004
!!$
!!$!use types,only:rprec
!!$!use sim_param,only:u,v,w,theta,q
!!$!use bottombc
!!$
!!$implicit none
!!$real(kind=rprec),dimension(nz)::ubar,vbar,wbar,T_bar,tke_bar,tke_sum
!!$real(kind=rprec)::rms, noise, arg, arg2,theta_mean
!!$real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
!!$real(kind=rprec),dimension(:,:,:),allocatable::wtemp
!!$integer::jx,jy,jz,seed,jz_abs,ii
!!$character (64) :: fname, temp
!!$
!!$tke_bar=0._rprec !initialize mean TKE profile as 0
!!$tke_sum=0._rprec !initialize mean TKE profile as 0
!!$
!!$!      if (patch_flag .eq. 1) theta_mean=theta_s1
!!$!      if (patch_flag .eq. 1) theta_mean=T_init
!!$       theta_mean=T_init
!!$       perturb_height_factor=800._rprec/z_i
!!$      
!!$      if (wt_s .eq. 0.0_rprec) then
!!$      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
!!$! w_star is of O(1) with z_i=500 and wt_s=0.06
!!$      T_star=0.06_rprec/w_star
!!$      q_star=T_star
!!$      else
!!$      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
!!$      T_star=wt_s/w_star
!!$      q_star=T_star
!!$      end if
!!$
!!$         $if ($MPI)
!!$            print *,'Modified Log Profile for IC for coord = ',coord
!!$         $else
!!$            print *,'Modified Log Profile for IC'
!!$         $endif
!!$       do jz=1,nz
!!$         $if ($MPI)
!!$            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$         $else
!!$            z=(real(jz)-0.5_rprec)*dz
!!$         $endif
!!$!        z=(real(jz)-0.5_rprec)*dz*z_i
!!$!c IC in equilibrium with rough surface (rough dominates in effective zo)
!!$        arg2=z/(sum(zo)/float(nx*ny))
!!$        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
!!$        if (coriolis_forcing) then
!!$        ubar(jz)=ug
!!$        vbar(jz)=vg
!!$        wbar(jz)=0._rprec
!!$! Note that ug and vg have already been non-dimensionalized in param.f90
!!$!        ubar(jz)=arg/30._rprec
!!$        else
!!$        ubar(jz)=arg
!!$        vbar(jz)=0._rprec
!!$        wbar(jz)=0._rprec
!!$        end if
!!$!C sc: I changed around some parenthesis here
!!$        if (z.gt.(perturb_height_factor*z_i)) then
!!$!        print *, 'Statement executed for the scalars'
!!$        ubar(jz)=ubar(jz-1)
!!$        end if
!!$       end do
!!$       
!!$!       open(unit=44,file=path//'Tbar_diurnal_GABLS.dat',status='unknown')
!!$!       read(44,5169) (T_bar(ii),ii=1,nz)
!!$5169     format(1400(E16.10))
!!$!       print *,'T_bar',T_bar
!!$
!!$!      do jz=1,nz
!!$!       print *,'k, ubar:',jz,ubar(jz)
!!$!      end do
!!$
!!$!  rms = 3._rprec
!!$  do jz=1,nz
!!$  $if ($MPI)
!!$    jz_abs = coord * (nz-1) + jz
!!$    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
!!$  $else
!!$    jz_abs = jz
!!$    z = (jz-.5_rprec) * dz * z_i
!!$  $endif
!!$  seed = -80 - jz_abs  !--trying to make consistent init for MPI
!!$  if (z .lt. 800._rprec) tke_bar(jz)=0.5_rprec*(1-z/800._rprec)
!!$
!!$  rms=(2._rprec*tke_bar(jz)/3._rprec)**0.5_rprec !2/3(1/2(u^2))=1/3(u^2) for each (u,v,w)
!!$    do jy=1,ny
!!$      do jx=1,nx
!!$!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!!$!c...Taking std dev of vel as 1 at all heights
!!$
!!$!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!!$!c u should also put L_z=z_i i.e. the inversion layer height should
!!$!c be equal to the height of the domain and in that case the second
!!$!c part of the subsequent if loop will never execute. This is
!!$!c ensured by putting an OR statement in the if clause, which makes 
!!$!c sure that only the first part of if block is executed and not the
!!$!c block after else
!!$
!!$!            z=(real(jz)-0.5_rprec)*dz*z_i
!!$       if (z.le.perturb_height_factor*z_i) then
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz) !noise
!!$         u(jx,jy,jz)=noise/u_star+ubar(jz) !noise
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$         v(jx,jy,jz)=noise/u_star+vbar(jz) !noise
!!$!         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
!!$         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
!!$         w(jx,jy,jz)=noise/u_star+wbar(jz) !noise
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
!!$         
!!$         call T_pos_gabls(theta_mean,z)
!!$!         theta_mean=T_bar(jz)
!!$         theta(jx,jy,jz)=theta_mean/T_scale
!!$         q(jx,jy,jz)=q_mix
!!$!         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
!!$!         q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
!!$       else
!!$         u(jx,jy,jz)=ubar(jz)
!!$         v(jx,jy,jz)=vbar(jz)
!!$         w(jx,jy,jz)=wbar(jz)
!!$         call T_pos_gabls(theta_mean,z)
!!$!         theta_mean=T_bar(jz)
!!$         theta(jx,jy,jz)=theta_mean/T_scale
!!$!         theta(jx,jy,jz)=(theta_mean+(z-perturb_height_factor*z_i)*inv_strength)/T_scale
!!$         q(jx,jy,jz)=q_mix
!!$        end if
!!$       end do
!!$     end do
!!$  end do
!!$
!!$
!!$  !...BC for W
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!!$    w(1:nx, 1:ny, 1) = 0._rprec
!!$  end if
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$    w(1:nx, 1:ny, nz) = 0._rprec
!!$  endif
!!$
!!$  !...BC for U, V & T
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
!!$    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
!!$!    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+inv_strength/T_scale*z_i*dz
!!$    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
!!$  end if
!!$
!!$!calculate initial TKE profile
!!$allocate(wtemp(ld,ny,nz)) !wtemp allocated
!!$wtemp=0._rprec
!!$wtemp(:,:,1:nz-1)=0.5_rprec*(w(:,:,1:nz-1)+w(:,:,2:nz));
!!$  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!!$   wtemp(:,:,nz)=w(:,:,nz);
!!$  end if
!!$do jz=1,nz
!!$   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
!!$   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((v(1:nx,1:ny,jz)-sum(v(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
!!$!   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((w(1:nx,1:ny,jz)-sum(w(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
!!$   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((wtemp(1:nx,1:ny,jz)-sum(wtemp(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
!!$end do
!!$deallocate(wtemp)
!!$
!!$!VK Display the mean vertical profiles of the initialized variables on the
!!$!screen
!!$    write(fname,'(A,i6.6,A)')path//'output/init_profiles.dat'
!!$    $if ($MPI)
!!$      write (temp, '(".c",i0)') coord
!!$      fname = trim (fname) // temp
!!$    $endif
!!$open(unit=44,file=fname,status="unknown",position="append")
!!$do jz=1,nz
!!$     $if ($MPI)
!!$       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!!$     $else
!!$       z = (jz - 0.5_rprec) * dz
!!$     $endif
!!$     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,tke_sum(jz)*u_star*u_star
!!$     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!!$     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!!$     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,real(tke_sum(jz))*u_star*u_star
!!$end do
!!$close(44)
!!$
!!$    write(fname,'(A,i6.6,A)')path//'output/init_profiles_3d.bin'
!!$    $if ($MPI)
!!$      write (temp, '(".c",i0)') coord
!!$      fname = trim (fname) // temp
!!$    $endif
!!$open(unit=44,file=fname,form="unformatted")
!!$write(44) real(u),real(v),real(w),real(theta);close(44)
!!$
!!$7781 format('jz, z, ubar, vbar, wbar,T_bar,TKE:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F16.10))
!!$!print *,'tke_bar',tke_bar
!!$
!!$end subroutine ic_scal_GABLS_diurnal
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine T_pos_gabls(T1,z_loc)
implicit none
real(kind=rprec):: T1,z_loc
if (z_loc<=200._rprec) then
    T1=288._rprec-z_loc*(288._rprec-286._rprec)/200._rprec;
else if (z_loc>200._rprec .AND. z_loc<=850._rprec) then
    T1=286._rprec
else if (z_loc>850._rprec .AND. z_loc<=1000._rprec) then
    T1=286._rprec+(z_loc-850._rprec)*(292._rprec-286._rprec)/(1000._rprec-850._rprec);
else if (z_loc>1000._rprec .AND. z_loc<=2000._rprec) then
    T1=292._rprec+(z_loc-1000._rprec)*(300._rprec-292._rprec)/(2000._rprec-1000._rprec);
else if (z_loc>2000._rprec .AND. z_loc<=3500._rprec) then
    T1=300._rprec+(z_loc-2000._rprec)*(310._rprec-300._rprec)/(3500._rprec-2000._rprec);
else if (z_loc>3500._rprec .AND. z_loc<=4000._rprec) then
    T1=310._rprec+(z_loc-3500._rprec)*(312._rprec-310._rprec)/(4000._rprec-3500._rprec);
else
    print *,'z not contained in the if block range !!!'
    stop
end if
end subroutine T_pos_gabls
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX






subroutine pollen_slice()
!c This is exactly the same like the subroutine avgslice with the
!c only difference being that it averages the scalar variables
!c to find the y-averaged instantaneous x-z slices of variables
!c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
!c It also outputs the average covariance between wt and wq
!use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3 
!use output_slice,only: collocate_MPI_averages

  implicit none

  integer:: i, j, k

  real(kind=rprec),dimension(nx,nz-1),save:: atheta,t2,q2,asgs_t3,awt
  real(kind=rprec),dimension(nx,nz-1),save:: adTdz,anu_t,t3,var_t,asc_t,acs2sc_t
  real(kind=rprec),dimension(nx,nz-1),save:: asgs_PCon1,auPCon

  real(kind=rprec):: ttheta1,tt2,tsgst,twt,tdTdz,arg1,arg2,fr
  real(kind=rprec):: tnu_t,tt3,tsc_t,tcs2sc_t
  real(kind=rprec):: tsgs_PCon1,tuPCon

  real(kind=rprec),dimension(:,:),allocatable:: avg_scalar_out

  
  fr=(1._rprec/float(p_count))*float(c_count)

  if (jt .EQ. c_count) then
    atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec
    anu_t=0._rprec;t3=0._rprec;var_t=0._rprec;asc_t=0._rprec;acs2sc_t=0._rprec
    asgs_PCon1=0._rprec;auPCon=0._rprec
  end if

  do k=1,nz-1
    do i=1,nx
      ttheta1=0._rprec;tt2=0._rprec;tsgst=0._rprec;twt=0._rprec;tdTdz=0._rprec
      tnu_t=0._rprec;tt3=0._rprec;tsc_t=0._rprec;tcs2sc_t=0._rprec
      tsgs_PCon1=0._rprec;tuPCon=0._rprec

      do j=1,ny  
        ttheta1=ttheta1+PCon(i,j,k,npcon)
        tsgst=tsgst+sgs_PCon3(i,j,k)
        tdTdz=tdTdz+dPCondz(i,j,k)
        ! tnu_t=tnu_t+Nu_t(i,j,k)
        tsc_t=tsc_t+Kc_t(i,j,k)
        tcs2sc_t=tcs2sc_t+Cs2Sc(i,j,k)
	tt2=tt2+PCon(i,j,k,npcon)*PCon(i,j,k,npcon)
        tsgs_PCon1=tsgs_PCon1+sgs_PCon1(i,j,k)
        tuPCon=tuPCon+u(i,j,k)*PCon(i,j,k,npcon)

        IF (PCon_scheme==1) THEN
	
  	  if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
            arg1=0._rprec
          else  
            arg1=(PCon(i,j,k,npcon)+PCon(i,j,k-1,npcon))/2.
          end if 
          twt=twt+w(i,j,k)*arg1
	
	ELSE
	  twt=twt+res_PCon3(i,j,k)
	END IF

      end do


      var_t(i,k)=var_t(i,k)+fr*sum((PCon(1:nx,1:ny,k,npcon)-sum(PCon(1:nx,1:ny,k,npcon))/(nx*ny))**2)/(nx*ny)

      atheta(i,k)=atheta(i,k)+(fr)*ttheta1/ny
      asgs_t3(i,k)=asgs_t3(i,k)+(fr)*tsgst/ny
      awt(i,k)=awt(i,k)+(fr)*twt/ny
      adTdz(i,k)=adTdz(i,k)+(fr)*tdTdz/ny
      ! anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t/ny
      asc_t(i,k)=asc_t(i,k)+(fr)*tsc_t/ny
      acs2sc_t(i,k)=acs2sc_t(i,k)+(fr)*tcs2sc_t/ny
      t2(i,k)=t2(i,k)+(fr)*tt2/ny
      asgs_PCon1(i,k)=asgs_PCon1(i,k)+(fr)*tsgs_PCon1/ny
      auPCon(i,k)=auPCon(i,k)+(fr)*tuPCon/ny

    end do
  end do
      
  if (mod(jt,p_count)==0) then
    allocate(avg_scalar_out(1:nx,1:nz_tot-1));
    call collocate_MPI_averages_N(atheta,avg_scalar_out,35,'PCon')
    call collocate_MPI_averages_N(t2,avg_scalar_out,36,'PCon2')
    call collocate_MPI_averages_N(asgs_t3,avg_scalar_out,37,'sgs_PCon3')
    call collocate_MPI_averages_N(awt,avg_scalar_out,38,'wPCon')
    call collocate_MPI_averages_N(adTdz,avg_scalar_out,39,'dPCondz')
    call collocate_MPI_averages_N(asgs_PCon1,avg_scalar_out,40,'sgs_PCon1')
    call collocate_MPI_averages_N(auPCon,avg_scalar_out,41,'uPCon')
    call collocate_MPI_averages_N(t3,avg_scalar_out,46,'PCon3')
    call collocate_MPI_averages_N(var_t,avg_scalar_out,47,'var_PCon');
    call collocate_MPI_averages_N(asc_t,avg_scalar_out,48,'Kc_sgs');
    call collocate_MPI_averages_N(acs2sc_t,avg_scalar_out,49,'Cs2Sc_sgs');
    deallocate(avg_scalar_out)

   atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec;
   anu_t=0._rprec;t3=0._rprec;
   var_t=0._rprec;asc_t=0._rprec;acs2sc_t=0._rprec;
   asgs_PCon1=0._rprec;auPCon=0._rprec;
  end if

5168     format(1400(E14.5))

end subroutine pollen_slice


end module scalars_module2
