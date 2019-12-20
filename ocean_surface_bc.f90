module bottombc
use types,only:rprec
use param,only:nx,ny,ld,z_ref,rag06,fieldsize,fieldsize_stripe,field_xs,field_xe,field_ys,field_ye
implicit none
! num_patch= numbr of patches, zo?=surface roughness for the patches types
! ptypes=number of patches types to be used, usually we use 2
integer,parameter::num_patch=1,ptypes=2
!real(kind=rprec),parameter::zo1=.1_rprec,zo2=.46_rprec
real(kind=rprec),parameter::zo1=0.15_rprec,zo2=0.15_rprec
real(kind=rprec),parameter::do1=0.00_rprec,do2=0.75_rprec
integer,parameter::square=0
real(kind=rprec),parameter::zo_out=.15_rprec, zo_square=.15_rprec
integer,parameter::Nx_sq=5,Ny_sq=3
!real(kind=rprec),parameter::t_scale=300._rprec
!real(kind=rprec),parameter::theta_s1=300.15_rprec,theta_s2=310.15_rprec
real(kind=rprec),parameter::theta_s1=288._rprec,theta_s2=320.15_rprec
!real(kind=rprec),parameter::theta_s1=310.15_rprec,theta_s2=300.15_rprec
!real(kind=rprec),parameter::theta_s1=305.15_rprec,theta_s2=300.15_rprec

real(kind=rprec),dimension(nx,ny)::ustar_avg
! This is the common local ustar computed in obukhov and used in wallstress.f90
real(kind=rprec),dimension(nx,ny)::zo,T_s,q_s
real(kind=rprec),dimension(nx,ny)::z_os
!TS add for non-neutral case
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m,phi_h,psi_h
!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90) 
!VK for non-neutral scenario
integer,dimension(nx,ny)::patch
integer,dimension(num_patch)::patchnum


! Variables added for pollen
! Chamecki - 08/04/2006
real(kind=rprec),dimension(nx,ny)        :: zo_PCon            ! Roughness for pollen concentration
real(kind=rprec),dimension(nx,ny)        :: PCon_s             ! Surface concentration of pollen

real(kind=rprec),dimension(nx,ny)        :: d0                 ! Displacement height


contains
subroutine patches()
!VK This assigns momentum roughness, temperature and wetness 
!VK for the different patches
!VK and fills the lookup tables patch and patchnum
!VK This is called from routines patch_or_remote.f90 depending
!VK whether to use remotely-sensed data or patches
!use types,only:rprec
use param
implicit none

integer::i, j, patchtype, begini, endi, type
!integer,dimension(ptypes),intent(inout)::patchnum
!integer,dimension(nx,ny),intent(inout)::patch

! sc: 
! without this, some compilers may give total junk
! this was cause the u_avg to be WAY off, and so txz was WAY off, etc.
patchnum=0
!
if (square==1) then!for square patch
   zo(:,:)=zo_out/z_i
   patch(:,:)= 1
   zo(Nx/2-(Nx_sq-1)/2:Nx/2+(Nx_sq-1)/2,&
        Ny/2-(Ny_sq-1)/2:Ny/2+(Ny_sq-1)/2) = zo_square/z_i
   patch(Nx/2-(Nx_sq-1)/2:Nx/2+(Nx_sq-1)/2,&
        Ny/2-(Ny_sq-1)/2:Ny/2+(Ny_sq-1)/2) = 2
!TS   patchnum(2)=Nx_sq*Ny_sq
!TS   patchnum(1)=Nx*Ny-patchnum(2)
   patchnum=Nx*Ny
else!for regular patches
   do j=1,ny
      endi = 0
      type=1
      do patchtype=1,num_patch
         begini=endi+1
         endi=(nx*patchtype)/num_patch
         do i=begini,endi
            if (type.eq.1) then
               zo(i,j)=zo1/z_i
               T_s(i,j)=theta_s1/t_scale
               q_s(i,j)=q_s1
	       

	       ! Added for pollen
	       PCon_s(i,j)=PCon_sfc

	       
               patch(i,j)=type
            end if
            if(type.eq.2) then
               zo(i,j)=zo2/z_i
               T_s(i,j)=theta_s2/t_scale
               q_s(i,j)=q_s2
               patch(i,j)=type
            end if
            patchnum(type)=patchnum(type)+1
          end do
          if (type.eq.1) then
             type=2
          elseif (type.eq.2) then
             type=1
          end if
       end do       
      end do
end if

  

  ! This is to setup the zo for the simulation of the ragweed field experiment
  IF (rag06) THEN
  
    ! Upwind field
    zo(1:6,:)=zo1/z_i
    d0(1:6,:)=do1/z_i
    PCon_s(1:6,:)=0._rprec
    
    ! Ragweed stripe
    zo(7:10,:)=zo2/z_i
    d0(7:10,:)=do2/z_i
    PCon_s(7:10,:)=PCon_sfc
    
    ! Downwind field
    zo(11:100,:)=zo1/z_i
    d0(11:100,:)=do1/z_i    
    PCon_s(11:100,:)=0._rprec
    
    ! Correct the sides of the stripe
    zo(:,1:25)=zo1/z_i
    d0(:,1:25)=do1/z_i
    PCon_s(:,1:25)=0._rprec
    zo(:,46:50)=zo1/z_i
    d0(:,46:50)=do1/z_i
    PCon_s(:,46:50)=0._rprec
    

    ! Set reference height for pollen (usually set in scalarr_module2.f90)
    zo_PCon(:,:)=zo(:,:)+d0(:,:)
  
  ELSE IF (fieldsize) THEN
  
    ! Homogeneous
    zo(:,:)=zo1/z_i
    d0(:,:)=do1/z_i
    PCon_s(:,:)=0._rprec


    ! Ragweed field 4
    zo(6:37,35:66)=zo2/z_i
    d0(6:37,35:66)=do2/z_i
    PCon_s(6:37,35:66)=PCon_sfc
    

    ! Set reference height for pollen (usually set in scalarr_module2.f90)
    zo_PCon(:,:)=zo(:,:)+d0(:,:)
  
  
  ELSE IF (fieldsize_stripe) THEN
  
    ! Homogeneous
    zo(:,:)=zo1/z_i
    d0(:,:)=do1/z_i
    PCon_s(:,:)=0._rprec


    ! Ragweed field 4
    zo(field_xs:field_xe,:)=zo2/z_i
    d0(field_xs:field_xe,:)=do2/z_i
    PCon_s(field_xs:field_xe,:)=PCon_sfc
    

    ! Set reference height for pollen (usually set in scalarr_module2.f90)
    zo_PCon(:,:)=zo(:,:)+d0(:,:)
  
  ELSE
  
    ! Homogeneous
    zo(:,:)=zo1/z_i
    d0(:,:)=do1/z_i
    PCon_s(:,:)=PCon_sfc
  
    ! Set reference height for pollen (usually set in scalarr_module2.f90)
    zo_PCon(:,:)=zo(:,:)+d0(:,:)

  END IF

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  OPEN(1,FILE=path//'output/check_patches.out',STATUS="unknown",POSITION="append")
  DO i=1,nx
    WRITE(1,'(201F18.8)')i*dx,(zo(i,j),j=1,ny)
  END DO
  CLOSE(1)
  end if



end subroutine patches
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
subroutine remote_to_patch(T_s_in,T_or_z)
implicit none
integer:: ii,jj,kk,T_or_z
real(kind=rprec),dimension(nx,ny):: T_s_in,dummy
real(kind=rprec):: sigma_multiplier,crap1,crap2,crap3
real(kind=rprec):: patch_cold, patch_hot

sigma_multiplier=1._rprec
dummy=0._rprec
if (T_or_z == 1) then
    crap1=sum(T_s_in)/float(nx*ny)
elseif (T_or_z == 2) then
    crap1=exp(sum(dlog(T_s_in))/float(nx*ny))
else
    print *,'Wrong choice of T_or_z in remote_to_patch().. STOPPING'
    stop
end if

crap2=sqrt(sum((T_s_in-crap1)**2)/float(nx*ny))
!crap2=sqrt(sum((T_s_in-crap1)**2._rprec)/(real(nx)*real(ny)))

patch_hot=crap1+sigma_multiplier*crap2;
patch_cold=crap1-sigma_multiplier*crap2;

print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
!First do the patch business for temperature
if ((patch_hot .lt. 0._rprec) .OR. (patch_cold .lt. 0._rprec)) then
print *,'Hot & Cold patch calculation yields negative T_s'
!print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
print *,'Trying sigma_multiplier = 0.75'
     if (patch_cold < 0._rprec) then
        sigma_multiplier=0.75_rprec
        patch_cold=crap1-sigma_multiplier*crap2;
        patch_hot=crap1+sigma_multiplier*crap2;
        crap3=0.5_rprec*(patch_cold+patch_hot)
        print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
        if (patch_cold < 0._rprec) then
           print *,'sigma = 0.75 FAILED, Trying sigma_= 0.5'
           sigma_multiplier=0.5_rprec
           patch_cold=crap1-sigma_multiplier*crap2;
           patch_hot=crap1+sigma_multiplier*crap2;
           crap3=0.5_rprec*(patch_cold+patch_hot)
           print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
           if (patch_cold < 0._rprec) then
              print *,'sigma = 0.5 FAILED, STOPPING NOW...'
              print *,'This message is from the subroutine patch_to_remote in scalars_module2.f90.'
           end if
        end if
     end if
end if
! Assign to patches
if (T_or_z .eq. 1) then
  dummy(1:nx/2,:)=patch_hot
  dummy(nx/2+1:nx,:)=patch_cold
   print *,'2 patch T_s field generated: patch_hot,patch_cold = ',patch_hot,patch_cold
elseif (T_or_z .eq. 2) then
! Assign cold patch to the first patch and hot patch to the second one
! Note that this is exactly opposite to what we do for temperature as
! in REALITY, hot surfaces have low roughness and vice-versa
  dummy(1:nx/2,:)=patch_cold
  dummy(nx/2+1:nx,:)=patch_hot
   print *,'2 patch roughnesss field generated: patch_smooth,patch_rough = ',patch_cold,patch_hot
end if

T_s_in=dummy

end subroutine remote_to_patch

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX

subroutine avgpatch(u_avg)
! computes the averaged value of a variable (at the wall) over a patch
! and assigns it to an nx X ny array

! sc: 
! note: this is inefficient: should calculate directly to u_avg
! --no maybe its good
!use types,only:rprec
use sim_param,only:u
!use bottombc,only:ptypes
implicit none
real(kind=rprec),dimension(nx,ny),intent(inout)::u_avg
!integer,dimension(ptypes),intent(in)::patchnum
!integer,dimension(nx,ny),intent(in)::patch
integer::i,j,k
real(kind=rprec),dimension(ptypes)::temp
!real(kind=rprec),dimension(ld,ny,1),intent(in)::u

temp=0._rprec
do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      temp(patch(i,j))=temp(patch(i,j))+u(i,j,1)/real(patchnum(patch(i,j)))
   end if
end do
end do
end do

do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      u_avg(i,j)=real(temp(k))
   end if
end do
end do
end do
end subroutine avgpatch

end module bottombc
