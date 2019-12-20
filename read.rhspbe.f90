
module types
implicit none

public
integer, parameter :: rprec = kind (1.d0)
end module types
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
  INTEGER,PARAMETER                          :: jt_start = 150
  INTEGER,PARAMETER                          :: jt_end = 200
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2_rprec
  REAL(rprec),PARAMETER                      :: L_y = 0.78125_rprec
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,avg_dissip,tot
  REAL(rprec),dimension(:,:,:),allocatable::u,v,w,avg_u,avg_v,avg_w,dissip,rhs_v
  REAL(rprec),dimension(:,:,:,:),allocatable::Pcon0,br0,Re0,rhs_pbe0,rhs
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,avg_Pcon,Rhs_pbe,tot_vol
  CHARACTER(len=64)                           :: file
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=330
  real,dimension(npcon)                       ::dsd, diameter,vol
  real(rprec) :: q_src,dx=l_x/nx,dy=L_y/ny,dz=z_i/nz_tot 

real(rprec),parameter::pi=3.1415926535897932384626433_rprec
  dir ='output_data/'
   allocate(Pcon0(ld,ny,1:nz,npcon),Pcon(ld,ny,1:nz_tot,npcon),Rhs_pbe0(ld,ny,1:nz,npcon),Rhs_pbe(ld,ny,1:nz_tot,npcon))
   allocate(u0(ld,ny,1:nz),tot_vol(ld,ny,nz_tot,npcon),rhs(ld,ny,nz_tot,npcon),rhs_v(ld,ny,nz_tot))
   allocate(dissip0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),br0(ld,ny,1:nz,npcon),Re0(ld,ny,1:nz,npcon))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)
 write(*,*) diameter
open(10,file="rhs_pbe_q.dat",position="append")
q_src =(6.047e4_rprec/(dx*z_i*dy*z_i*dz*z_i))/(1_rprec/(z_i/u_star))!*dx*dy*dz*pi/6_rprec*diameter(npcon)**3
write(*,*) q_src
jt = 300
!!  do  jt = 50,100
   jt_total = jt*base
        do ip=1,nproc
      write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(dir),'vel_pcon_00',jt_total,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip) u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz), &
                  br0(:,:,1:nz,1:npcon),Re0(:,:,1:nz,:),Rhs_pbe0(:,:,1:nz,1:npcon)
      nzs=(ip-1)*(nz-1)
      Rhs_pbe(:,:,nzs+1:nzs+nz,1:npcon) = Rhs_pbe0(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
  print *,"Opened files,timestep=",jt_total
 
  do ip=1,npcon
     rhs(:,:,:,ip) = pi/6_rprec*Rhs_pbe(:,:,:,ip)*diameter(ip)**3_rprec
   enddo
  write(*,*) maxval(rhs),minval(rhs)
!write(*,*) (rhs(xps,yps,zps,ip),ip=1,npcon)
!write(*,*) sum(rhs(xps,yps,zps,:))
rhs_v = sum(rhs(:,:,:,1:npcon),4)
write(*,*) size(rhs_v)
  write(*,*) rhs_v(xps,yps,zps),minloc(rhs_v),minval(rhs_v),sum(rhs_v(2:320,2:149,20:383)),maxval(rhs_v)
do ip=1,npcon
dsd(ip) = sum(rhs(20:80,60:80,100:340,ip))*dx*dy*dz
enddo
!do jz=1,nz_tot
!do jy=1,ny
!do jx=1,nx
!
!dsd = dsd + rhs_v(jx,jy,jz)
!enddo
!enddo
!enddo
!write(*,*) dsd
!vol = pi/6_rprec*dsd*diameter**3
!write(*,*) sum(vol)

write(10,*) (dsd(ip),ip=1,npcon),sum(dsd),dsd(npcon)+q_src
!enddo
close(10)
505 format(18e15.7)
!!open(unit=22,FILE='rhsv_40000.dat')
!!write(22,*) 'variables=x,y,z,rhs_v'
!!write(22,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!!      do jz=1,nz_tot
!!          do jy=1,ny
!!            do jx=1,nx
!!               write(22,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,rhs_v(jx,jy,jz)
!!                
!!            enddo
!!          enddo
!!       enddo
!!close(22)
!


END PROGRAM POSTPROCESSING
