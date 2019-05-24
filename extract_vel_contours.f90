module types
implicit none

public
integer, parameter :: rprec = kind (1.d0)
end module types

!      #######################################################################
!      #######################################################################
!      ##                         postprocessing file for velocity          ##   
!      ##                         and concentration data written by         ##
!      ##                         lesgo. The file is written as binary      ##
!      ##                         data for each processor at every          ##
!      ##                         "base" timesteps.                         ##
!      ##                                                                   ##
!      #######################################################################

PROGRAM POSTPROCESSING
use types,only:rprec
  implicit none

!Parameters
  INTEGER,PARAMETER                           :: nproc=96
  REAL(rprec),PARAMETER                       :: TINYS=1e-24
  INTEGER,PARAMETER                           :: ny=150
  INTEGER,PARAMETER                           :: nz_tot=385
  INTEGER,PARAMETER                           :: nz = (385-1)/nproc +1
  INTEGER,PARAMETER                           :: nx=384
  INTEGER,PARAMETER                           :: nsteps=200
  INTEGER,PARAMETER                           :: base=200
  INTEGER,PARAMETER                           :: npcon=15
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =1
  INTEGER,PARAMETER                          :: jt_end = 30
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2
  REAL(rprec),PARAMETER                      :: L_y = 0.78125
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0001
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,dmax
  REAL(rprec),dimension(:,:,:),allocatable   ::u,v,w,avg_u,avg_v,avg_w,dissip,lpcon
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon0,breakage_freq0,Re0
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,avg_Pcon!,lpcon
  CHARACTER(len=64)                           :: file,path,plane
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,plane_cut,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329

  dir ='output_data/'
    print *,"hello"
   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),dissip0(ld,ny,1:nz),pcon0(nx,ny,1:nz,npcon))
   allocate(u(nx,ny,1:nz_tot),v(nx,ny,1:nz_tot),w(nx,ny,1:nz_tot),dissip(nx,ny,1:nz_tot),&
          avg_w(ld,ny,1:nz_tot),pcon(nx,ny,1:nz_tot,npcon),lpcon(nx,ny,1:nz_tot))
allocate(dmax(nx,ny,nz_tot))  
avg_w=0._rprec  
!  do jt=jt_start,jt_end

 path = 'post_proc_output/av_375_650/'
   write(file,'(A,A)')TRIM(PATH), "av_dissip_375_650.out"
!   write(file,'(A)') "Pcon_30000.out"
   open(20,FILE=file,form='unformatted')
!   read(20) PCon(1:nx,1:ny,1:nz_tot,:)
   read(20) dissip(:,:,1:nz_tot)
  close(20)
   write(file,'(A,A)')TRIM(PATH), "av_pcon_375_650.out"
!   write(file,'(A)') "Pcon_30000.out"
   open(20,FILE=file,form='unformatted')
!   read(20) PCon(1:nx,1:ny,1:nz_tot,:)
   read(20) Pcon(:,:,1:nz_tot,:)
  close(20)
!print * , minval(dissip)

dissip = -dissip*u_star**3

dmax = dissip**(-0.4_rprec)*(1.9e-2_rprec/1018._rprec)**(0.6_rprec)
!    counter = jt-jt_start+1
!    jt=30
! !   counter=1
!    jt_total=jt*base
!    do ip=1,nproc
!       write(file,'(A,A,I6.6,A,I2.2,A)')TRIM(dir),'vel_pcon_0',jt_total,'_20',ip-1,'.out'
!       open(unit=2000+ip,FILE=file,form='unformatted')
!       read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,:),dissip0(:,:,1:nz) 
!       nzs=(ip-1)*(nz-1)
!       u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
!       v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
!       w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
!       PCon(:,:,nzs+1:nzs+nz,:)=Pcon0(:,:,1:nz,:)
!       dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)      
!       close(2000+ip)
!    enddo
!    avg_w = avg_w +w/nt
!    write(*,*) maxloc(-dissip(:,:,:)*u_star**3),maxval(-dissip*u_star**3)
!    
!!    write(*,*) "at zps"
!!    do jx=xps-1,xps+1
!!    do jy=yps-1,yps+1
!!       write(*,*) -dissip(jx,jy,zps)*u_star**3,jx,jy
!!    enddo
!!    enddo
!!       do jx=xps-1,xps+1
!!    do jy=yps-1,yps+1
!!       write(*,*) -dissip(jx,jy,zps+1)*u_star**3,jx,jy
!!    enddo
!!    enddo
!!
!!       write(*,*) "at zps+1"
!!    do jx=xps-1,xps+1
!!    do jy=yps-1,yps+1
!!       write(*,*) -dissip(jx,jy,zps+2)*u_star**3,jx,jy
!!    enddo
!!    enddo
!
    lpcon = 1e-20_rprec
    where(dissip .ge. 1e-20_rprec) lpcon = log10(dissip)

    where(pcon .le. 1e-2)
      pcon=1e-10_rprec
    elsewhere
      pcon=1._rprec
    endwhere
    
    print *, maxval(lpcon)
 !  write(*,*) maxval(pcon(:,:,:,1)),maxloc(pcon(:,:,:,1))
!!    print *,"Opened files,timestep=",jt_total
!
!!   enddo
!  write(*,*) maxval(-avg_w*u_star),maxloc(-avg_w*u_star)
505 format(18e15.3)
!   path='post_proc_output/output_plane/'

   plane = 'yps'
write(file,'(A,A,A,A)')TRIM(path),'av_dissip_masked_dmax',TRIM(plane),'.dat'
open(20,File=file,status='unknown')
write(20,*) 'variables=x,z,log(eps),dmax'
write(20,*)'zone t="',1,'"i=',nx,'k=',nz_tot,'f=point'
do jz=1,nz_tot
  do jx=1,nx
    write(20,505)(jx-1)*L_x/nx,-(jz-1)*z_i/nz_tot,lpcon(jx,yps,jz)*pcon(jx,yps,jz,5),dmax(jx,yps,jz)
  enddo
enddo

!!Select plane
!   plane_cut=yps
!   write(file,'(A,A,A,A,I6.6,A)')TRIM(path),'vel_noocean',TRIM(plane),'_',jt_total,'.dat'
!   open(unit=2000+counter,FILE=file,status='unknown')
!   write(2000+counter,*) 'variables=x,z,avg_u,avg_v,avg_w,dissip'
!   write(2000+counter,*) 'zone t="',581+counter,'"i=',nx,'k=',nz_tot,'f=point'
!   do jz=1,nz_tot
!     do jx=1,nx
!       write(2000+counter,505)(jx-1)*L_x/nx,-(jz-1)*z_i/nz_tot,u(jx,plane_cut,jz),&
!                 v(jx,plane_cut,jz)*u_star,-w(jx,plane_cut,jz)*u_star,-dissip(jx,plane_cut,jz)*u_star**3
!     enddo
!   enddo
!   close(2000+counter)
!
!
!   write(file,'(A,A,A,A,I6.6,A)')TRIM(path),'pcon_',TRIM(plane),'_',jt_total,'.dat'
!   open(unit=2000+counter,FILE=file,status='unknown')
!   write(2000+counter,*) 'variables=x,z,p1,p7,p12,p15'
!   write(2000+counter,*) 'zone t="',581+counter,'"i=',nx,'k=',nz_tot,'f=point'
!   do jz=1,nz_tot
!     do jx=1,nx
!       write(2000+counter,505)(jx-1)*L_x/nx,-(jz-1)*z_i/nz_tot,lpcon(jx,plane_cut,jz,1),lpcon(jx,plane_cut,jz,7),&
!                              lpcon(jx,plane_cut,jz,12),lpcon(jx,plane_cut,jz,15)
!     enddo
!   enddo
!   close(2000+counter)
!  enddo
deallocate(u0,v0,w0,u,v,w,pcon0,pcon,dissip)
END PROGRAM POSTPROCESSING










