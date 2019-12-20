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
use scalars_module2,only:ic_scal_v2!,ic_scal_gabls,ic_scal_GABLS_diurnal ! added by VK
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

Cs_opt2_avg=0._rprec
fx=0._rprec;fy=0._rprec;fz=0._rprec
u_des=0._rprec;v_des=0._rprec;w_des=0._rprec
mean_u=0._rprec;mean_u2=0._rprec;mean_v=0._rprec;mean_v2=0._rprec
mean_w=0._rprec;mean_w2=0._rprec

if(initu)then
   
   if(S_FLAG.and.PCon_FLAG) then
      write(fname,'(A,I4.4,A)') path//'restart/vel_sc_pcon_',2000+coord,'.out'
   else if(S_FLAG) then
      write(fname,'(A,I4.4,A)') path//'restart/vel_sc_',2000+coord,'.out'
   else if(PCon_FLAG) then
      write(fname,'(A,I4.4,A)') path//'restart/vel_pcon_',2000+coord,'.out'
   else
      write(fname,'(A,I4.4,A)') path//'restart/vel_',2000+coord,'.out'
   endif

   open(unit=2000+coord,file=fname,form="unformatted")
   
!   if(initsc) then
      
   if (S_FLAG .AND. PCon_FLAG) then !WITH SCALARS AND POLLENS
      read (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
           PCon(:,:,1:nz,1:npcon), RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz), &
           RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m,RHS_PCon(:,:,1:nz,1:npcon), &
           Cs_opt2, F_LM, F_MM, F_QN, F_NN, deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
      print *,'Reading initial velocity, temperature, and particle concentration from file'
   elseif (S_FLAG) then !WITH SCALARS
      read (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
           RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
           RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m,Cs_opt2, F_LM, F_MM, F_QN, F_NN
      print *,'Reading initial velocity and temperature from file'
   elseif (PCon_FLAG) then ! Pollen
      read (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), PCon(:,:,1:nz,1:npcon),   &
           RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
           RHS_PCon(:,:,1:nz,1:npcon), Cs_opt2, F_LM, F_MM,F_QN, F_NN,            &
           deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
      print *,'Reading initial velocity and particle concentration from file'
   else ! No SCALARS
      read (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
           RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
           Cs_opt2, F_LM, F_MM, F_QN, F_NN
      print *,'Reading initial velocity from file'
   end if

   close(2000+coord)
   
!   end if

else
   if (dns_bc) then
      print*, 'Creating initial velocity field with DNS BCs'
      call ic_dns()
   else
      print*, 'Creating initial fields'
      if (S_FLAG) then
         if (GABLS_test) then
            print *, 'Creating initial velocity & scalar fields for the GABLS test case'
!!$            call ic_scal_GABLS()
            print*, "Subroutine ic_scal_GABLS is commented out to save memory"
            stop
         else if (GABLS_diurnal_test) then
            print *, 'Creating initial velocity & scalar fields for the GABLS diurnal case'
!!$            call ic_scal_GABLS_diurnal()
            print*, "Subroutine ic_scal_GABLS_diurnal is commented out to save memory"
         else
            print *, 'Creating initial velocity & scalar fields'
            call ic_scal_v2()
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
if (PCon_FLAG) then
   do ipcon=1,npcon
      call mpi_sendrecv (PCon(1, 1, nz-1, ipcon), ld*ny, MPI_RPREC, up, 1,  &
           PCon(1, 1, 0, ipcon), ld*ny, MPI_RPREC, down, 1,   &
           comm, status, ierr)
   enddo
endif
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
