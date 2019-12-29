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
  INTEGER,PARAMETER                          :: jt_start =2
  INTEGER,PARAMETER                          :: jt_end = 12
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
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,avg_u,avg_v,avg_w
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon0,breakage_freq0,Re0
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,avg_Pcon
  CHARACTER(len=64)                           :: file
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2
  integer(kind=8)  :: reclen
  dir ='output_data/'
  allocate(Pcon0(ld,ny,1:nz,npcon))
!    print *,"hello"
   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz))
   allocate(u(ld,ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot),avg_u(ld,ny,1:nz_tot),avg_v(ld,ny,1:nz_tot))
   allocate(avg_w(ld,ny,1:nz_tot))
   allocate(dissip0(ld,ny,1:nz),dissip(ld,ny,1:nz_tot),avg_dissip(ld,ny,1:nz_tot))
  avg_w=0._rprec
  avg_u=0._rprec
  avg_v=0._rprec
   avg_dissip=0._rprec
  avg_w=0._rprec
  do jt=jt_start,jt_end
    counter = jt-jt_start+1
 !   counter=1
    jt_total=jt*base
    do ip=1,nproc
       write(file,'(A,A,I6.6,A,I2.2,A)')TRIM(dir),'vel_pcon_0',jt_total,'_20',ip-1,'.out'
       open(unit=2000+ip,FILE=file,form='unformatted')
       read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,:),dissip0(:,:,1:nz)
       nzs=(ip-1)*(nz-1)
       u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
       v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
       w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
       dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
       close(2000+ip)
    enddo
  print *,"Opened files,timestep=",jt_total
 ! print *,minloc(w(xps,yps,:)), minval(w(xps,yps,:))!,maxval(w(xps,yps,:))
  avg_u = avg_u +u/nt    
    avg_v = avg_v +v/nt
     avg_w = avg_w +w/nt
    avg_dissip = avg_dissip + dissip/nt
  enddo

  print *,minval(avg_u(xps,yps,:)),maxval(avg_u(xps,yps,:))
  print *,minval(avg_v(xps,yps,:)),maxval(avg_v(xps,yps,:))
  print *,minval(avg_w(xps,yps,:)),maxval(avg_w(xps,yps,:))
505 format(18e15.7)
!inquire(iolength=reclen) avg_w(1:nx,1:ny,1:nz_tot)
!open(unit=20,file='av_uvweps.out',access='direct',recl=reclen,form='unformatted',action='read')
!read(20,rec=3) avg_w
!read(20,rec=4)avg_dissip
!close(20)
!
deallocate(u,v,w,u0,v0,w0)
deallocate(Pcon0)
deallocate(dissip0,dissip)

!inquire(iolength=reclen) avg_dissip(1:nx,1:ny,1:nz_tot)
!
!open(unit=20,file='av_uvweps.out',access="direct",form='unformatted',recl=reclen)
!write(20,rec=1) avg_u(1:nx,1:ny,1:nz_tot)
!write(20,rec=2) avg_v(1:nx,1:ny,1:nz_tot)
!write(20,rec=3) avg_w(1:nx,1:ny,1:nz_tot)
!!write(20,rec=4) avg_dissip(1:nx,1:ny,1:nz_tot)
!close(20)

!open(unit=20,file='center_vel.dat',form='formatted',status='unknown')
!write(20,*) 'z','U_c'
!do jz=1,nz_tot
! write(20,505) (jz-1)*z_i/nz_tot,avg_w(xps,yps,jz)
!enddo
!close(20)
!
!open(unit=20,file='center_dissip.dat',form='formatted',status='unknown')
!write(20,*) 'z','dissip'
!do jz=1,nz_tot
! write(20,505) (jz-1)*z_i/nz_tot,avg_dissip(xps,yps,jz)
!enddo
!close(20)
!deallocate(avg_u,avg_v,avg_w,avg_dissip)

END PROGRAM POSTPROCESSING










