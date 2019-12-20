subroutine initial()
use types,only:rprec
use param
use sim_param,only:path,u,v,w,RHSx,RHSy,RHSz,theta,q,PCon,ke
!+++++++++++++++++++++++++++
!DY Changed by Di Yang
!DY use sgsmodule , only : Cs_opt2, Cs_opt2_avg, F_LM, F_MM, F_QN, F_NN 
use sgsmodule , only : Cs_opt2, Cs_opt2_avg, F_LM, F_MM, F_QN, F_NN, F_KX,F_XX,F_KX2,F_XX2
!DY End here
!+++++++++++++++++++++++++++
!use bottombc,only:psi_m ! added by VK
use scalars_module,only:RHS_T,sgs_t3,psi_m,RHS_PCon,deposition,Real_dep,Kc_t ! added by VK
use scalars_module2,only:ic_scal,ic_scal_gabls,ic_scal_GABLS_diurnal ! added by VK
! VK -label 122 assigned to vel_sc.out for reading input files in case of
! scalars
!!!!XXXXXXXXXX--------Added by Vijayant----XXXXXXX!!!!!
use immersedbc,only:fx,fy,fz,u_des,v_des,w_des,n_bldg,bldg_pts
use io,only:mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

!real(kind=rprec),dimension(ld,ny,nz)::crap
logical, parameter :: use_add_random = .false.
character (64) :: fname, temp
real(kind=rprec),dimension(ld,ny,nz)::crap
real(kind=rprec)::z
real(kind=rprec)::temp_w

integer::i,jz,jx,jy, ipcon

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif

real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSx_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSy_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSz_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot
real (rprec), dimension(ld, ny, 1:nz_tot, npcon) :: PCon_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: Cs_opt2_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: RHS_T_tot
real(kind=rprec),dimension(ld,ny,nz_tot,npcon) :: RHS_PCon_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: Kc_t_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_KX_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_XX_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_KX2_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_XX2_tot

!real (rprec), dimension(ld, ny, 1:nz_tot, 3) :: tmp1,tmp2

Cs_opt2_avg=0._rprec
fx=0._rprec;fy=0._rprec;fz=0._rprec
u_des=0._rprec;v_des=0._rprec;w_des=0._rprec
mean_u=0._rprec;mean_u2=0._rprec;mean_v=0._rprec;mean_v2=0._rprec
mean_w=0._rprec;mean_w2=0._rprec

!VK Modified so that when running scalars, the file is named vel_sc.out
!VK else called vel.out

if (S_FLAG) then
    fname = path // 'vel_sc.out'
else
    fname = path // 'vel.out'
end if

!$if ($MPI)
!  write (temp, '(".c",i0)') coord
!  fname = trim (fname) // temp
!$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  open(11,file=fname,form='unformatted')
endif

if(initu)then
  if(initsc) then
  
    print *,'Reading initial velocity and temperature from file'

    if((.not. USE_MPI)) then

    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),theta(:,:,1:nz), &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz),RHSz(:, :, 1:nz),         &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m
      case (2:3)
        read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),    &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),         &
                 RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
                 
      case (4)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2 
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),  &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),       &
                    RHS_T(:,:,1:nz),sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM,  &
                    crap, crap
        end if
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),  &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),       &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM, &
                    F_QN, F_NN
        end if
      case default
        write (*, *) 'initial: invalid model number'
    end select

    else if ((USE_MPI .and. coord == 0)) then

!++++++++++++++++++++++++++++++++++++++++++++++
!DY Changed by Di Yang
!DY Because this is not consistent with the output format in subroutine checkpoint_final (in io.f90).
!DY This causes the simulation crash during restart.
!DY The new format is copied directly from subroutine checkpoint_final.
!DY The current version is incomplete, mainly because the "if...else" structure and logical conditions are inconsistent
!DY with those in subroutine checkpoint_final.  Use the current version for testing only.

!DY    select case (model)
!DY      case (1)
!DY        read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot), &
!DY                  RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot),RHSz_tot(:, :, 1:nz_tot),         &
!DY                  RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m
!DY      case (2:3)
!DY        read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot),    &
!DY                 RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),         &
!DY                 RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m, Cs_opt2_tot

!DY      case (4)
!DY        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
!DY          read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot), &
!DY                    RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),      &
!DY                    RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m, Cs_opt2_tot
!DY        else
!DY          read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot),  &
!DY                    RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),       &
!DY                    RHS_T_tot(:,:,1:nz_tot),sgs_t3(:,:,1), psi_m, Cs_opt2_tot, F_LM_tot, F_MM_tot
!DY        end if
!DY      case (5)
!DY        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
!DY          read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot), &
!DY                    RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),      &
!DY                    RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m, Cs_opt2_tot
!DY        else
!DY          read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),theta_tot(:,:,1:nz_tot),  &
!DY                    RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),       &
!DY                    RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m, Cs_opt2_tot, F_LM_tot, F_MM_tot, &
!DY                    F_QN_tot, F_NN_tot
!DY        end if
!DY      case default
!DY        write (*, *) 'initial: invalid model number'
!DY    end select
       if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
          read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot),   &
               PCon_tot(:,:,1:nz_tot,1:npcon), RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot), &
               RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m,RHS_PCon_tot(:,:,1:nz_tot,1:npcon), &
               Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot, &
               deposition(:,:),Real_dep,Kc_t_tot,F_KX_tot,F_XX_tot,F_KX2_tot,F_XX2_tot
!!$          read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot),   &
!!$               tmp1(:,:,1:nz_tot,1:3), RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot), &
!!$               RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m,tmp2(:,:,1:nz_tot,1:3), &
!!$               Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot, &
!!$               deposition(:,:),Real_dep,Kc_t_tot,F_KX_tot,F_XX_tot,F_KX2_tot,F_XX2_tot
!          do ipcon=1,npcon
!             do jz=1,nz_tot
!                do jy=1,ny
!                   do jx=1,nx
!                      PCon_tot(jx,jy,jz,ipcon)=tmp1(jx,jy,jz,ipcon)
!                      RHS_PCon_tot(jx,jy,jz,ipcon)=tmp2(jx,jy,jz,ipcon)
!                   enddo
!                enddo
!             enddo
!          enddo
       elseif (S_FLAG) then !WITH SCALARS
          read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot),   &
               RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),          &
               RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m,Cs_opt2_tot, F_LM_tot, F_MM_tot, &
               F_QN_tot, F_NN_tot
       elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
          read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), PCon_tot(:,:,1:nz_tot,1:npcon),   &
               RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),          &
               RHS_PCon_tot(:,:,1:nz_tot,1:npcon), Cs_opt2_tot, F_LM_tot, F_MM_tot,F_QN_tot, F_NN_tot,            &
               deposition(:,:),Real_dep,Kc_t_tot,F_KX_tot,F_XX_tot,F_KX2_tot,F_XX2_tot
       else ! No SCALARS
          read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot),           &
               RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),  &
               Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
       end if
!DY Changes by Di Yang end here
!+++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for testing
!  write(100,*) 'variables=x,y,z,PCon'
!       open(200)
!       write(200,*) 'zone t="restart" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
!       do jz=1,nz_tot
!          do jy=1,ny
!             do jx=1,nx
!                write(200,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot,PCon_tot(jx,jy,jz), &
!                     u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz),theta_tot(jx,jy,jz)
!             enddo
!          enddo
!       enddo
!505    format(8e12.4)
!       close(200)
!DY End here
!++++++++++++++++++++++++++++++++++++++++++

    endif    

!    END IF


    
!TS INITIALIZE THE ZERO CONCENTRATION FIELD IF jt_total=0
    if(passive_scalar)then
       open(1,file=path//'run')
       read(1,*)i
       close(1)
       if(i.eq.0)then
          theta=0._rprec;RHS_T=0._rprec
       endif
    endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then    
open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
    if (.not. USE_MPI) then
    do jz=1,nz
            z=(real(jz)-0.5_rprec)*dz
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
    enddo
    else
    do jz=1,nz_tot
            z=(real(jz)-0.5_rprec)*dz
     write(6,7781) jz,z,(sum(u_tot(:,:,jz))/float(nx*ny))*u_star,(sum(v_tot(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w_tot(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta_tot(:,:,jz))/float(nx*ny))*T_scale
     write(44,7781) jz,z,(sum(u_tot(:,:,jz))/float(nx*ny))*u_star,(sum(v_tot(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w_tot(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta_tot(:,:,jz))/float(nx*ny))*T_scale
    end do
    endif
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
close(44)
endif

  else

    print *,'Reading initial velocity field from file'

    if((.not. USE_MPI)) then

    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),       &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz)
      case (2:4)
        read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2
        else
          read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2, F_LM, F_MM, F_QN, F_NN
        end if
      case default
        write (*, *) 'initial: invalid model number'
    end select

    else if ((USE_MPI .and. coord == 0)) then

    select case (model)
      case (1)
        read (11) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot),       &
                  RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot)
      case (2:4)
        read(11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),             &
                 RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),  &
                 Cs_opt2_tot
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),             &
                    RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),  &
                    Cs_opt2_tot
        else
          read(11) u_tot(:, :, 1:nz_tot),v_tot(:, :, 1:nz_tot),w_tot(:, :, 1:nz_tot),             &
                   RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),  &
                   Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
        end if
      case default
        write (*, *) 'initial: invalid model number'
    end select

    endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord ==0)) then
    if (.not. USE_MPI) then
    do jz=1,nz
      write(6,7780) jz,sum(u(:,:,jz))/float(nx*ny),sum(v(:,:,jz))/&
           float(nx*ny),sum(w(:,:,jz))/float(nx*ny)
    end do
    endif
    if (USE_MPI .and. coord == 0) then
    do jz=1,nz_tot
      write(6,7780) jz,sum(u_tot(:,:,jz))/float(nx*ny),sum(v_tot(:,:,jz))/&
           float(nx*ny),sum(w_tot(:,:,jz))/float(nx*ny)
    end do
    endif
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))
endif

  end if
$if ($MPI)
  sendcounts = size (u(:,:,1:nz-1))
  recvcounts = size (u(:,:,1:nz))
  displs = coord_of_rank * sendcounts
  call mpi_scatterv (u_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    u(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (v_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    v(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (w_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    w(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (RHSx_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    RHSx(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (RHSy_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    RHSy(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (RHSz_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    RHSz(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (Cs_opt2_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    Cs_opt2(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_LM_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_LM(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_MM_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_MM(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_QN_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_QN(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_NN_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_NN(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  if (S_FLAG) then
  call mpi_scatterv (theta_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    theta(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (RHS_T_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    RHS_T(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif

!+++++++++++++++++++++++++++++++
!DY Added by Di Yang
  if (PCon_FLAG) then
     do ipcon=1,npcon
        call mpi_scatterv (PCon_tot(1,1,1,ipcon), recvcounts, displs,MPI_RPREC,&
             PCon(1,1,1,ipcon), recvcounts,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
        call mpi_scatterv (RHS_PCon_tot(1,1,1,ipcon), recvcounts, displs,MPI_RPREC,&
             RHS_PCon(1,1,1,ipcon), recvcounts,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  call mpi_scatterv (Kc_t_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    Kc_t(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  sendcounts = size (F_KX(:,:,1:nz-1))
  recvcounts = size (F_KX(:,:,1:nz))
  call mpi_scatterv (F_KX_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_KX(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_XX_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_XX(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_KX2_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_KX2(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_scatterv (F_XX2_tot(1,1,1), recvcounts, displs,MPI_RPREC,&
                    F_XX2(1,1,1), recvcounts,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  
  endif
!DY End here
!++++++++++++++++++++++++++++++++

$endif
!    do jz=1,nz
!      write(6,7790) jz,sum(u(:,:,jz))/float(nx*ny),sum(v(:,:,jz))/&
!           float(nx*ny),sum(w(:,:,jz))/float(nx*ny)
!    end do
!7790 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))

else
  if (dns_bc) then
     print*, 'Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    print*, 'Creating initial fields'
    if (S_FLAG) then
       if (GABLS_test) then
         print *, 'Creating initial velocity & scalar fields for the GABLS test case'
         call ic_scal_GABLS()
       else if (GABLS_diurnal_test) then
         print *, 'Creating initial velocity & scalar fields for the GABLS diurnal case'
         call ic_scal_GABLS_diurnal()
       else
         print *, 'Creating initial velocity & scalar fields'
         call ic_scal()
       end if
    else
       print*, 'Creating initial velocity field'
       call ic()
    end if
  end if
end if

! bldg stuff
if (use_bldg) then
   open(1,file=path//'bldg.dat')
   read(1,*) n_bldg
   allocate(bldg_pts(5,n_bldg))
   do i=1,n_bldg
      read(1,*) bldg_pts(1:5,i)
      if(bldg_pts(5,i).ge.nz)then
         write(*,*)"lz should be less than nz"
         stop
      end if
   end do
   close(1)
end if

$if ($MPI)
  !--synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz 
  call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                     u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                     v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                     w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                     comm, status, ierr)
  call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
                     u(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 5,  &
                     v(1, 1, nz), ld*ny, MPI_RPREC, up, 5,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 6,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, 6,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, 7,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, 7,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, 8,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, 8,   &
                     comm, status, ierr)
!++++++++++++++++++++++++
!DY Added by Di Yang
  if (PCon_FLAG) then
     do ipcon=1,npcon
        call mpi_sendrecv (PCon(1, 1, nz-1, ipcon), ld*ny, MPI_RPREC, up, 1,  &
             PCon(1, 1, 0, ipcon), ld*ny, MPI_RPREC, down, 1,   &
             comm, status, ierr)
     enddo
  endif
!DY End here
!+++++++++++++++++++++++++++
$endif

if (USE_MPI .and. coord == 0) then
  !--set 0-level velocities to BOGUS
  u(:, :, $lbz) = BOGUS
  v(:, :, $lbz) = BOGUS
  w(:, :, $lbz) = BOGUS
  theta(:, :, $lbz) = BOGUS
end if

! kinetic energy
do jz=1,nz-1
do jy=1,ny
do jx=1,nx

   temp_w=.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
   ke(jx,jy,jz)=.5_rprec*(u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)

end do
end do
end do

! Always initialize pollen concentration with zeros
! Chamecki - 08/01/2006
IF (PCon_FLAG .AND. .NOT.(initPCon)) THEN
  PRINT*,'Pollen profile initialized with zeros'
  PCon(:,:,:,:)=0._rprec
  if (USE_MPI .and. coord == 0) then
   PCon(:,:,$lbz,:) = BOGUS
  end if
END IF

!DY+++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY Force negative particle concentration to be zero
!DY Start here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

!if(PCon_FLAG) then
!   do ipcon=1,npcon
!      do jz=$lbz,nz
!         do jy=1,ny
!            do jx=1,ld
!               if(PCon(jx,jy,jz,ipcon).lt.0.) PCon(jx,jy,jz,ipcon)=0.
!            enddo
!         enddo
!      enddo
!   enddo
!endif

!DY+++++++++++++
!DY End here
!DY+++++++++++++


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_random ()
implicit none

real (rprec), parameter :: rms = 0.2_rprec

integer :: i, j, k
integer :: seed

real (rprec) :: noise
real (rprec) :: ran3

!---------------------------------------------------------------------

seed = -80

do k = 1, nz
  do j = 1, ny
    do i = 1, nx
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      u(i, j, k) = u(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      v(i, j, k) = v(i, j, k) + noise
      noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
      w(i, j, k) = w(i, j, k) + noise
    end do
  end do
end do

end subroutine add_random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine initial
