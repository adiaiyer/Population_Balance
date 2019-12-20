! put everything onto w-nodes, follow original version
!--provides txx, txy, tyy, tzz for jz=1:nz-1; txz, tyz for 1:nz
subroutine sgs_stag ()
use types,only:rprec
use param
use sim_param,only: u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,  &
                    txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:u_lag,v_lag,w_lag,Cs_opt2,Nu_t,magS
!------------------------- Vij Comment begins---------------
! 04/14/2004 - Added Nu_t to the list of variables from sgsmodule
! 05/19/2004 - Replaced all references to visc by Nu_t; deleted local 
!              declaration of visc
! 05/19/2004 - Replace all .5 by 0.5
!-------------------------Vij Comment ends ------------------
use immersedbc,only:building_mask,building_interp
use test_filtermodule,only:filter_size
!+++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for no-slip seabed
use bottombc,only:zo,psi_m,phi_m,ustar_avg,d0
use test_filtermodule
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++
$if ($TREES)
  use trees
$endif
implicit none
real(kind=rprec),dimension(nz)::l,ziko,zz
real (rprec), dimension (ld, ny, nz) :: S11, S12, S22, S33, S13, S23
!real(kind=rprec),dimension(ld,ny,nz):: dissip
real(kind=rprec)::ux, uy, uz, vx, vy, vz, wx, wy, wz
real(kind=rprec),dimension(ld,ny) :: txzp, tyzp,S
real(kind=rprec) :: delta, nu, const
integer::jx,jy,jz
integer :: jz_min, jz_max

!+++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for no-slip seabed
real(kind=rprec),dimension(nx,ny)::u_avg,u_aver,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec),parameter::zo_seabed=0.0001_rprec/z_i
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++

if (VERBOSE) write (*, *) 'started sgs_stag'

delta=filter_size*(dx*dy*dz)**(1._rprec/3._rprec) ! nondimensional
! Cs is Smagorinsky's constant. l is a filter size (non-dim.)  

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! save the wall level
  txzp(:,:)=txz(:,:,1)
  tyzp(:,:)=tyz(:,:,1)

  ! calculate Sij on w-nodes
  ! calculate |S| on w-nodes

  ! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
  do jy=1,ny
  do jx=1,nx              
     ux=dudx(jx,jy,1)  ! uvp-node
     uy=dudy(jx,jy,1)  ! uvp-node
     uz=dudz(jx,jy,1)  ! uvp-node
     vx=dvdx(jx,jy,1)  ! uvp-node
     vy=dvdy(jx,jy,1)  ! uvp-node
     vz=dvdz(jx,jy,1)  ! uvp-node 
     ! special case
     wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
     wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
     wz=dwdz(jx,jy,1)  ! uvp-node
     S11(jx,jy,1)=ux          ! uvp-node
     S12(jx,jy,1)=0.5_rprec*(uy+vx) ! uvp-node
     ! taken care of with wall stress routine
     S13(jx,jy,1)=0.5_rprec*(uz+wx) ! uvp
     S22(jx,jy,1)=vy          ! uvp-node
     ! taken care of with wall stress routine 
     S23(jx,jy,1)=0.5_rprec*(vz+wy) ! uvp
     S33(jx,jy,1)=wz          ! uvp-node
  end do
  end do

  jz_min = 2

else

  jz_min = 1

end if

$if ($MPI)
  !--this is only required b/c of the unnatural placement of all strains
  !  onto w-nodes, be careful not to overwrite nz on top process with garbage
  !--dwdz(jz=0) is already known, except at bottom process (OK)
  !call mpi_sendrecv (dwdz(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
  !                   dwdz(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
  !                   comm, status, ierr)
  call mpi_sendrecv (dwdz(1, 1, 1), ld*ny, MPI_RPREC, down, 2,  &
                     dwdz(1, 1, nz), ld*ny, MPI_RPREC, up, 2,   &
                     comm, status, ierr)
$endif

! calculate derivatives/strain on w-nodes
!$vvohmygod parallel do default(shared) private(jx,jy,jz)	
!--in MPI version, calculating up to nz saves some interprocess exchange
!  later but note that dwdz is not provided w/o some communication
!  (unless its the top process) 
do jz=jz_min, nz
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11(jx,jy,jz)=ux          ! w-node
   S12(jx,jy,jz)=0.5_rprec*(uy+vx) ! w-node
   S13(jx,jy,jz)=0.5_rprec*(uz+wx) ! w-node
   S22(jx,jy,jz)=vy          ! w-node
   S23(jx,jy,jz)=0.5_rprec*(vz+wy) ! w-node
   S33(jx,jy,jz)=wz          ! w-node
end do
end do
end do
!$ffohmygod end parallel do

! This part computes the average velocity during cs_count times steps
! This is used with the lagrangian model only
if (model == 4 .OR. model==5) then
   u_lag = u_lag+u
   v_lag = v_lag+v
   w_lag = w_lag+w
end if

if (sgs) then!ref01
   if((model == 1))then  !For traditional Smagorinsky   ref02


     $if ($TREES)

       l = delta
       call trees_Cs ()

     $else

       ! Define parameters (Co and nn) for wallfunction
       Cs_opt2 = Co**2  ! constant coefficient

       if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
         do jz = 1, nz

           if (jz == 1) then
              zz(jz) = (jz - 0.5_rprec) * dz
           else  ! w-nodes
              zz(jz) = (jz - 1) * dz
           end if

           ! z's nondimensional, l here is on w-nodes, except at bottom
           l(jz) = ( Co**(nnn)*(vonk*zz(jz))**(-nnn) +  &
                     (delta)**(-nnn) )**(-1._rprec/nnn)
         end do
       else
         do jz = 1, nz
           zz(jz) = ((jz - 1) + coord * (nz - 1)) * dz  !--w-nodes

           ! z's nondimensional, l here is on w-nodes, except at bottom
           l(jz) = ( Co**(nnn)*(vonk*zz(jz))**(-nnn) +  &
                     (delta)**(-nnn) )**(-1._rprec/nnn)

         end do 
       end if

     $endif

   else ! for dynamic procedures: initialization  !cont ref02
      l = delta  ! constant equal to delta
      if ((jt == 1) .and. (inilag)) then !ref05
            print *,'CS_opt2 initialiazed'
            Cs_opt2 = 0.03_rprec

! make sure that the "node conventions" in these dynamic models
! match with those done here
      elseif(((jt.GE.DYN_init).OR.(initu)).AND.(mod(jt,cs_count)==0)) then!cont ref05
         if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
             if (jt ==DYN_init .OR. jt .eq. cs_count) print *,'running dynamic model = ',model
         end if
         if (model == 2) then  ! Standard dynamic model !ref06
           !ziko = Cs_opt2(5,5,:)  !--what does this do?
           !--provides ziko 1:nz
           call std_dynamic(ziko,S11,S12,S13,S22,S23,S33)
           forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)
         else if (model==3) then !cont ref06
           ! Plane average dynamic  continue ref 05
           !ziko = Cs_opt2(5,5,:)  !--what does this do?
           call scaledep_dynamic(ziko,u,v,w,S11,S12,S13,S22,S23,S33)
           forall (jz = 1:nz) Cs_opt2(:, :, jz) = ziko(jz)
         else if (model==4) then ! Lagrangian SS !cont ref06
           call lagrange_Ssim(S11,S12,S13,S22,S23,S33)
         elseif (model==5) then ! LagSD !cont ref06
           call lagrange_Sdep(S11,S12,S13,S22,S23,S33)
         end if !end ref06
      end if !end ref05
      ! for if jt == 1 and 2 else ifs end ref 5
   end if !(for model =1 or else) end ref02   
end if ! for sgs-if end ref01

!           do jz = 1, nz
!           print *,'jt,Cs_opt2',jt,sum(Cs_opt2(1:nx,1:ny,jz))/float(nx*ny)
!           end do

! define |S| and viscosity on w-nodes (on uvp node for jz=1)
!$comp parallel do default(shared) private(jx,jy,jz)
!--MPI: going up to nz saves communication
do jz = 1, nz
do jy=1,ny
do jx=1,nx
   S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +&
        S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +&
        S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
   Nu_t(jx,jy,jz)=S(jx,jy)*Cs_opt2(jx,jy,jz)*(l(jz)**2)
   !DY********************************************************
   !DY Commented out by Di Yang for higher resolution simulations
   !DY dissip(jx,jy,jz)=-(S(jx,jy)**3)*Cs_opt2(jx,jy,jz)*(l(jz)**2)
   !DY********************************************************
 
   ! Store |S|
   magS(jx,jy,jz)=S(jx,jy)
   
end do
end do

!if (mod(jt,100) == 0) then
!      if(jz==1.or.jz==nz/4.or.jz==nz/2) then
!         write(90)real(jt*dt),real(jz),real(Cs_opt2(1:NX,1:NY,jz))
!         write(93)real(jt*dt),real(jz),real(Nu_t(1:NX,1:NY,jz))
!         write(94)real(jt*dt),real(jz),real(dissip(1:NX,1:NY,jz))
!      end if
!end if
end do

!	if (mod(jt,50)==0) pause
!	pause
!$comp end parallel do

! evaluate nu_t= c_s^2 l^2 |S|
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! bottom
  !$comp parallel do default(shared) private(jx,jy,const)
  do jy=1,ny
  do jx=1,nx
     nu = 0._rprec
     const = 0._rprec
     if (sgs) then
        const = Nu_t(jx,jy,1)
     end if
     if (molec) then
        nu = (nu_molec/(u_star*z_i)) ! dim/less 
     end if
  ! everything on right node here
     txx(jx,jy,1) = -2._rprec*(const+nu)*S11(jx,jy,1)
     txy(jx,jy,1) = -2._rprec*(const+nu)*S12(jx,jy,1)
     tyy(jx,jy,1) = -2._rprec*(const+nu)*S22(jx,jy,1)
     tzz(jx,jy,1) = -2._rprec*(const+nu)*S33(jx,jy,1)
  end do
  end do
  !$comp end parallel do

  jz_min = 2  !--wall level set by wallstress routine

else

  jz_min = 1

end if

! middle
!$comp parallel do default(shared) private(jx,jy,jz,const)	
!--note that nz level is available for the interpolations at (jz+1) here
!--txx, tyy, tzz, txy not needed at nz
do jz=jz_min, nz - 1
do jy=1,ny
do jx=1,nx
   nu=0._rprec
   const=0._rprec
   if (sgs) then
      const=0.5_rprec*(Nu_t(jx,jy,jz) + Nu_t(jx,jy,jz+1))
   end if
   if (molec) then
      nu=(nu_molec/(u_star*z_i)) ! dim/less
   end if
   txx(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S11(jx,jy,jz) + S11(jx,jy,jz+1))
   txy(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S12(jx,jy,jz) + S12(jx,jy,jz+1))
   tyy(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S22(jx,jy,jz) + S22(jx,jy,jz+1))
   tzz(jx,jy,jz)=-2._rprec*(const+nu)*&
        0.5_rprec*(S33(jx,jy,jz) + S33(jx,jy,jz+1))
end do
end do
end do
!$comp end parallel do

!$comp parallel do default(shared) private(jx,jy,jz,const)	
do jz = jz_min, nz - 1
do jy=1,ny
do jx=1,nx
   nu = 0._rprec
   const = 0._rprec
   if (sgs) const=Nu_t(jx,jy,jz)
   if (molec) nu=(nu_molec/(u_star*z_i)) ! dim/less
   txz(jx,jy,jz)=-2._rprec*(const + nu) * S13(jx,jy,jz)
   tyz(jx,jy,jz)=-2._rprec*(const + nu) * S23(jx,jy,jz)
end do
end do
end do
!$comp end parallel do

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  ! don't need to do this, unless we accidentally mess with the wall level
  txz(:,:,1)=txzp(:,:)
  tyz(:,:,1)=tyzp(:,:)
end if

$if ($MPI)
  !--recv information for top nodes: txy, txz only
  !--other components not needed, since no equation is solved there
  !  (values are just copied)
  call mpi_sendrecv (txz(1, 1, 1), ld*ny, MPI_RPREC, down, 3,  &
                     txz(1, 1, nz), ld*ny, MPI_RPREC, up, 3,   &
                     comm, status, ierr)
  call mpi_sendrecv (tyz(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
                     tyz(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
                     comm, status, ierr)

  !--also set 0-layer to bogus values
  txx(:, :, 0) = BOGUS
  txy(:, :, 0) = BOGUS
  txz(:, :, 0) = BOGUS
  tyy(:, :, 0) = BOGUS
  tyz(:, :, 0) = BOGUS
  tzz(:, :, 0) = BOGUS  !--tzz is updated in main (move to here)

  txx(:, :, nz) = BOGUS
  txy(:, :, nz) = BOGUS
  tyy(:, :, nz) = BOGUS
  tzz(:, :, nz) = BOGUS
  
$endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Changed by Di Yang for no-slip seabed for ocean simulation
!DY In ocean simulations, the domain is actually upside down.
!DY jz=nz_tot-1 is the first grid above the bottom boundary in the physical space.
if(OCEAN_FLAG.and.SEABED_FLAG) then
   if(coord.eq.nproc-1) then
      u1=u(:,:,nz-1)
      v1=v(:,:,nz-1)
      call test_filter(u1,G_test)
      call test_filter(v1,G_test)
      denom=log(0.5_rprec*dz/zo_seabed)
      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
      ustar_avg=u_avg*vonk/denom
      do jy=1,ny
         do jx=1,nx
            const=(ustar_avg(jx,jy)**2) /u_avg(jx,jy)
            ! const has oppisite sign as the bc at jz=1
            txz(jx,jy,nz)=const *u1(jx,jy)
            tyz(jx,jy,nz)=const *v1(jx,jy)
            !this is as in Moeng 84
            dudz(jx,jy,nz)=-ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,nz-1)/u_avg(jx,jy)
            dvdz(jx,jy,nz)=-ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,nz-1)/u_avg(jx,jy)
            dudz(jx,jy,nz)=merge(0._rprec,dudz(jx,jy,nz),u(jx,jy,nz-1).eq.0._rprec)
            dvdz(jx,jy,nz)=merge(0._rprec,dvdz(jx,jy,nz),v(jx,jy,nz-1).eq.0._rprec)
         end do
      end do
!!$      print*, "In sgs_stag: no-slip seabed!"
!!$      print*, "txz(xps,yps,nz)=",txz(xps,yps,nz)
   endif
   goto 999 ! use wall function in subroutine wallstress for txz and tyz
   ! For ocean simulations, the domain is upside down.
   ! jz=nz_tot-1 is the first grid above the bottom boundary in the physical space.
elseif ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  txz(:,:,nz)=0._rprec
  tyz(:,:,nz)=0._rprec
end if
999 continue
!DY End here
!*****************************************************************

if (VERBOSE) write (*, *) 'finished sgs_stag'

end subroutine sgs_stag
