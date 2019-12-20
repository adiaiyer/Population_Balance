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
  INTEGER,PARAMETER                           :: base=100
  INTEGER,PARAMETER                           :: npcon=15
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =300
  INTEGER,PARAMETER                          :: jt_end = 400
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2._rprec
  REAL(rprec),PARAMETER                      :: L_y = 0.78125_rprec
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,dissip
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,breakage_freq,Re,tot_vol,rhs1,rhs2
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(10,npcon)               :: dsd,rhs_source,rhs_source2
  CHARACTER(len=64)                           :: file,path,path1
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,i,jx,jy,jz,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  real                                        :: vol 
real(rprec)                       ::               dx=L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,d50
  real(rprec),dimension(:,:),allocatable            ::M_tot,M_y,M_z,y_cm,z_cm,tot
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
!  dir ='output_data/'
allocate(Pcon(nx,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon))
allocate(rhs1(nx,ny,1:nz_tot,npcon),rhs2(nx,ny,1:nz_tot,npcon))
allocate(dissip(nx,ny,1:nz_tot))

open(21,File='diameter.dat')
do ip=1,npcon
  read(21,*) diameter(ip)
enddo
close(21)
path1 = 'post_proc_output/' 


write(file,'(A,A)')TRIM(PATH1), "av_rhs_pcon_300_400.out"
!write(file,'(A)') "Pcon_30000.out"
open(20,FILE=file,form='unformatted')

read(20) PCon(:,:,1:nz_tot,:),rhs1(:,:,1:nz_tot,:),rhs2(:,:,1:nz_tot,:),dissip(:,:,1:nz_tot)
!read(20) u(:,:,1:nz_tot),v(:,:,1:nz_tot),w(:,:,1:nz_tot),PCon(:,:,1:nz_tot,:),&
!          dissip(:,:,1:nz_tot)
close(20)
write(*,*) minval(rhs1(:,:,:,15)), minval(rhs1(:,:,:,13))

do ip=1,npcon
  tot_vol(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!  write(*,*) maxval(tot_vol(:,:,:,ip))
enddo
   

path ='post_proc_output/'

505 format(18e15.7)
d_half = (diameter(2:npcon)+ diameter(1:npcon-1))/2_rprec
d_d(1) = (d_half(1)-diameter(1))*2
d_d(npcon) = diameter(npcon)-d_half(npcon-1)
d_d(2:npcon-1) = d_half(2:npcon-1) - d_half(1:npcon-2)
d_d = d_d*1.e6_rprec
jy_start = yps-3
jy_end   = yps+3
jz_start = 168
jz_end   = 178
counter=1

do i=1,5
 jx = 96 + 48*(i-1)
do ip=1,npcon
 dsd(i,ip) = sum(pcon(jx,40:110,:,ip))/size(pcon(jx,40:110,:,ip))
enddo
enddo



write(file,'(A,A)')TRIM(path), "dsd_y2z.dat"
open(21,FILE=file,status='unknown')
   do i=1,5
  write(21,505) (40-1)*l_y/ny,(110-1)*l_y/ny,(96+(i-1)*48)*L_x/nx ,(dsd(i,ip)/d_d(ip),ip=1,npcon)   
   enddo
close(20)

deallocate(Pcon,tot_vol)
deallocate(rhs1,rhs2)
deallocate(dissip)
!deallocate(u,v,w,dissip,Re,tot_vol,tot,d50,breakage_freq)
END PROGRAM POSTPROCESSING



