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
  INTEGER,PARAMETER                           :: nproc=48
  REAL(rprec),PARAMETER                       :: TINYS=1e-24
  INTEGER,PARAMETER                           :: ny=288
  INTEGER,PARAMETER                           :: nz_tot=385
  INTEGER,PARAMETER                           :: nz = (nz_tot-1)/nproc +1
  INTEGER,PARAMETER                           :: nx=288
  INTEGER,PARAMETER                           :: nsteps=200
  INTEGER,PARAMETER                           :: base=5
  INTEGER,PARAMETER                           :: npcon=1
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =1
  INTEGER,PARAMETER                          :: jt_end = 4
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=1.
  REAL(rprec),PARAMETER                      :: L_y = 1.
  REAL(rprec),PARAMETER                      :: u_star = 1.
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  real(rprec)                                :: junk
  REAL(rprec),dimension(:,:,:),allocatable   ::u0,v0,w0,dissip0,dissip,avg_dissip
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,avg_u,avg_v,avg_w,u2,v2,w2
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon0,breakage_freq0,Re0,Pcon20
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,avg_Pcon,pcon2
  CHARACTER(len=64)                           :: file
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2
  integer(kind=8)  :: reclen
  dir ='output/'
!    print *,"hello"
   allocate(u0(nx,ny,1:nz),v0(nx,ny,1:nz),w0(nx,ny,1:nz))
   allocate(u(nx,ny,1:nz_tot),v(nx,ny,1:nz_tot),w(nx,ny,1:nz_tot))
   allocate(Pcon0(nx,ny,1:nz,npcon),Pcon(nx,ny,1:nz_tot,npcon))
   allocate(Pcon20(nx,ny,1:nz,npcon),Pcon2(nx,ny,1:nz_tot,npcon))
   allocate(u2(nx,ny,1:nz_tot),v2(nx,ny,1:nz_tot),w2(nx,ny,1:nz_tot))
  inquire(iolength=reclen) u0
    do ip=1,nproc
!ip =42
       write(file,'(A,A,A,I3.3,A)')TRIM(dir),'veluv_avg','.c',ip-1,'.bin'
       open(unit=2000+ip,FILE=file,form='unformatted',access='direct',action='read',recl=reclen)
       read(2000+ip,rec=1) u0(:,:,1:nz)
       read(2000+ip,rec=2) v0(:,:,1:nz)
       read(2000+ip,rec=3) w0(:,:,1:nz)
       nzs=(ip-1)*(nz-1)
       u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
       v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
       w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
      close(2000+ip)
!      print *,ip, minval(w(xps,yps,:)), maxval(w(xps,yps,:))
    enddo

    do ip=1,nproc
!ip =42
       write(file,'(A,A,A,I3.3,A)')TRIM(dir),'vel2_avg','.c',ip-1,'.bin'
       open(unit=2000+ip,FILE=file,form='unformatted',access='direct',action='read',recl=reclen)
       read(2000+ip,rec=1) u0(:,:,1:nz)
       read(2000+ip,rec=2) v0(:,:,1:nz)
       read(2000+ip,rec=3) w0(:,:,1:nz)
       nzs=(ip-1)*(nz-1)
       u2(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
       v2(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
       w2(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
      close(2000+ip)
!      print *,ip, minval(w(xps,yps,:)), maxval(w(xps,yps,:))
    enddo
    do ip=1,nproc
!ip =42
       write(file,'(A,A,A,I3.3,A)')TRIM(dir),'pcon_avg','.c',ip-1,'.bin'
       open(unit=2000+ip,FILE=file,form='unformatted',access='direct',action='read',recl=reclen)
       do ipcon = 1,npcon
       read(2000+ip,rec=ipcon) pcon0(:,:,1:nz,ipcon)
       enddo
       nzs=(ip-1)*(nz-1)
      Pcon(:,:,nzs+1:nzs+nz,:) = Pcon0(:,:,1:nz,:)
       close(2000+ip)
!      print *,ip, minval(w(xps,yps,:)), maxval(w(xps,yps,:))
    enddo
 
    do ip=1,nproc
!ip =42
       write(file,'(A,A,A,I3.3,A)')TRIM(dir),'pcon2_avg','.c',ip-1,'.bin'
       open(unit=2000+ip,FILE=file,form='unformatted',access='direct',action='read',recl=reclen)
       do ipcon = 1,npcon
       read(2000+ip,rec=ipcon) pcon20(:,:,1:nz,ipcon)
       enddo
       nzs=(ip-1)*(nz-1)
      Pcon2(:,:,nzs+1:nzs+nz,:) = Pcon20(:,:,1:nz,:)
       close(2000+ip)
!      print *,ip, minval(w(xps,yps,:)), maxval(w(xps,yps,:))
    enddo
   print * , minval(u(xps,yps,:)),maxval(u(xps,yps,:))
   print * , minval(v(xps,yps,:)),maxval(v(xps,yps,:))
   print * , minval(w(xps,yps,:)),maxval(w(xps,yps,:))
   print * , minval(pcon(xps,yps,:,1)),maxval(pcon(xps,yps,:,1))

   print * , minval(u2(xps,yps,:)),maxval(u2(xps,yps,:))
   print * , minval(v2(xps,yps,:)),maxval(v2(xps,yps,:))
   print * , minval(w2(xps,yps,:)),maxval(w2(xps,yps,:))
   print * , minval(pcon2(xps,yps,:,1)),maxval(pcon2(xps,yps,:,1))
   print * , (u(xps,yps,100)),(u(xps,yps,300))
   print * , (v(xps,yps,100)),(v(xps,yps,300))
   print * , (w(xps,yps,100)),(w(xps,yps,300))
   505 format(18e15.7)

!inquire(iolength=reclen) u(1:nx,1:ny,1:nz_tot)
!!print *, reclen
!open(unit=20,file='tavg_uvw.out',access="direct",form='unformatted',recl=reclen)
!write(20,rec=1) u(1:nx,1:ny,1:nz_tot)
!write(20,rec=2) v(1:nx,1:ny,1:nz_tot)
!write(20,rec=3) w(1:nx,1:ny,1:nz_tot)
!close(20)

deallocate(u,v,w,u0,v0,w0)

END PROGRAM POSTPROCESSING










