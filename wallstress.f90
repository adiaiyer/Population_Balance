! For use with staggered grid LES
! JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
subroutine wallstress ()
use types,only:rprec
!++++++++++++++++++++++++++++++++++++++++++++++
!DY Changed by Di Yang for ocean simulation
use param,only:jt,nsteps,dz,ld,lh,nx,ny,nz,vonk,lbc_mom,WALL_HOMOG_FLAG,coord,OCEAN_FLAG,nu_molec,nu_sgs,u_star,z_i,SEABED_FLAG,nproc
!DY coord added by Di Yang for debugging
!DY OCEAN_FLAG added by Di Yang for ocean simulation
use sgsmodule,only:Nu_t !DY Added by Di Yang for calculation of dudz(:,:,1)
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use bottombc,only:zo,psi_m,phi_m,ustar_avg,d0
use test_filtermodule
implicit none
integer::jx,jy
!real(kind=rprec),dimension(nx,ny)::ustar,u_avg,u_aver,denom
real(kind=rprec),dimension(nx,ny)::u_avg,u_aver,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec)::const
!+++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation
!real(kind=rprec)::
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation

if(OCEAN_FLAG) then

   if(coord.eq.0) then
      do jy=1,ny
         do jx=1,nx
            txz(jx,jy,1)=1. !DY Everything has been normalized by u_star
            tyz(jx,jy,1)=0.
!!            if(abs(Nu_t(jx,jy,1)).lt.1.e-10) then
!!!               print*, "Zero value of Nu_t for BC: in subroutine wallstress!"
!!!               stop
!!               dudz(jx,jy,1)=-txz(jx,jy,1)/(Nu_t(jx,jy,1)+nu_molec)
!!            else
!!               dudz(jx,jy,1)=-txz(jx,jy,1)/Nu_t(jx,jy,1)
!!            endif
!!$            if(Nu_t(jx,jy,1)*u_star*z_i.gt.1.e-3) then
!!$               dudz(jx,jy,1)=-txz(jx,jy,1)/(Nu_t(jx,jy,1)+nu_molec)
!!$            else if(dudz(jx,jy,2).lt.0.) then
!!$               dudz(jx,jy,1)=dudz(jx,jy,2)
!!$            else
!!$               dudz(jx,jy,1)=-txz(jx,jy,1)/nu_sgs
!!$            endif
            dudz(jx,jy,1)=dudz(jx,jy,2)
            dvdz(jx,jy,1)=0.
         enddo
      enddo
   endif

!!$!DY Added for no-slip seabed condition
!!$!DY In ocean simulations, the domain is actually upside down.
!!$!DY jz=nz_tot-1 is the first grid above the bottom boundary in the physical space.
!!$   print*, "checkpoint1: wallstress",seabed_flag,nproc
!!$   if(SEABED_FLAG) then
!!$      print*, "checkpoint2: wallstress"
!!$      if(coord.eq.nproc-1) then
!!$         print*, "checkpoint3: wallstress"
!!$         u1=u(:,:,nz-1)
!!$         v1=v(:,:,nz-1)
!!$         call test_filter(u1,G_test)
!!$         call test_filter(v1,G_test)
!!$         denom=log((0.5_rprec*dz-d0)/zo)
!!$         u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
!!$         ustar_avg=u_avg*vonk/denom
!!$         do jy=1,ny
!!$            do jx=1,nx
!!$               const=-(ustar_avg(jx,jy)**2) /u_avg(jx,jy)
!!$               txz(jx,jy,nz)=const *u1(jx,jy)
!!$               tyz(jx,jy,nz)=const *v1(jx,jy)
!!$               !this is as in Moeng 84
!!$               dudz(jx,jy,nz)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,nz-1)/u_avg(jx,jy)
!!$               dvdz(jx,jy,nz)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,nz-1)/u_avg(jx,jy)
!!$               dudz(jx,jy,nz)=merge(0._rprec,dudz(jx,jy,nz),u(jx,jy,nz-1).eq.0._rprec)
!!$               dvdz(jx,jy,nz)=merge(0._rprec,dvdz(jx,jy,nz),v(jx,jy,nz-1).eq.0._rprec)
!!$            end do
!!$         end do
!!$         print*, "In wallstress: no-slip seabed!"
!!$         print*, "txz(10,10,nz)=",txz(10,10,nz)
!!$      endif
!!$   endif
         
   goto 999 !DY Skip all the calculation for ABL

endif

!DY End here
!+++++++++++++++++++++++++++++++++++++++++++

select case (lbc_mom)

  case ('wall')

    u1=u(:,:,1)
    v1=v(:,:,1)

!DY Added by Di Yang for debugging
!!$if(coord.eq.0) then
!!$   do jy=1,ny
!!$      do jx=1,ld
!!$         write(201,*) "In divstress_uv: i,j=",jx,jy
!!$         write(201,*) "In divstress_uv: u1,v1=",u1(jx,jy),v1(jx,jy)
!!$      enddo
!!$   enddo
!!$endif
!DY end here

!    if ((patch_flag .eq. 1) .AND. (num_patch .eq. 1)) then
    if (WALL_HOMOG_FLAG) then
!!  calculate u_star in the average sense !!

!      denom=log(0.5_rprec*dz/zo)-sum(psi_m)/float(nx*ny)
      denom=log((0.5_rprec*dz-d0)/zo)-sum(psi_m)/float(nx*ny)

      u_avg=sum(sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2))/float(nx*ny)

!DY Added by Di Yang for debugging
!!$      if(coord.eq.0) then
!!$         do jy=1,ny
!!$            do jx=1,ld
!!$               write(202,*) "In divstress_uv: i,j=",jx,jy
!!$               write(202,*) "In divstress_uv: u1,v1,u_avg,psi_m=",u1(jx,jy),v1(jx,jy),u_avg(jx,jy),psi_m(jx,jy)
!!$            enddo
!!$         enddo
!!$      endif
!DY end here

      if (jt .eq. nsteps) print *,'USED WALL HOMOG conds in wallstress'
    else
      call test_filter(u1,G_test)
      call test_filter(v1,G_test)

!      denom=log(0.5_rprec*dz/zo)-psi_m
      denom=log((0.5_rprec*dz-d0)/zo)-psi_m

      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)

!DY Added by Di Yang for debugging
!!$      if(coord.eq.0) then
!!$         do jy=1,ny
!!$            do jx=1,ld
!!$               write(203,*) "In wallstress: i,j=",jx,jy
!!$               write(203,*) "In wallstress: d0,zo=",d0(jx,jy),zo(jx,jy)
!!$               write(203,*) "In wallstress: denom,u_avg,psi_m=",denom(jx,jy),u_avg(jx,jy),psi_m(jx,jy)
!!$               write(203,*) "In wallstress: u1,v1=",u1(jx,jy),v1(jx,jy)
!!$            enddo
!!$         enddo
!!$      endif
!DY end here

    end if
!    ustar=u_avg*vonk/denom
    ustar_avg=u_avg*vonk/denom

    do jy=1,ny
    do jx=1,nx
!       const=-(ustar(jx,jy)**2) /u_avg(jx,jy)
       const=-(ustar_avg(jx,jy)**2) /u_avg(jx,jy)
       txz(jx,jy,1)=const *u1(jx,jy)
       tyz(jx,jy,1)=const *v1(jx,jy)
    !this is as in Moeng 84

!       dudz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
!           *phi_m(jx,jy)
!       dvdz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
!           *phi_m(jx,jy)
       dudz(jx,jy,1)=ustar_avg(jx,jy)/((0.5_rprec*dz-d0(jx,jy))*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)
       dvdz(jx,jy,1)=ustar_avg(jx,jy)/((0.5_rprec*dz-d0(jx,jy))*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)


       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
    end do
    end do

!DY Added by Di Yang for debugging
!!$    if(coord.eq.0) then
!!$       do jy=1,ny
!!$          do jx=1,ld
!!$             write(204,*) "In divstress_uv: i,j=",jx,jy
!!$             write(204,*) "In wallstress: ustar_avg,u_avg,denom=",ustar_avg(jx,jy),u_avg(jx,jy),denom(jx,jy)
!!$             write(204,*) "In divstress_uv: u1,v1,u_avg=",u1(jx,jy),v1(jx,jy),u_avg(jx,jy)
!!$             write(204,*) "In divstress_uv: txz,tyz,dudz,dvdz=",txz(jx,jy,1),tyz(jx,jy,1),dudz(jx,jy,1),dvdz(jx,jy,1)
!!$          enddo
!!$       enddo
!!$    endif
!DY end here

  case ('stress free')

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case default

    write (*, *) 'invalid lbc_mom'

end select

999 continue

end subroutine wallstress
