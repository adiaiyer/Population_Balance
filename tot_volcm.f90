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
 REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,dissip
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,breakage_freq,pcon0,tot_vol
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(10,npcon)               :: dsd
  CHARACTER(len=64)                           :: file,path,path1
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  real                                        :: vol 
real(rprec)                       ::               dx=L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,d50
  real(rprec),dimension(:),allocatable        :: M_tot,M_y,M_z,y_cm,z_cm
  real(rprec),dimension(:,:),allocatable      :: tot
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
  dir ='output_data/'
!   allocate(Pcon(nx,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon))
!   allocate(u(nx,1:ny,1:nz_tot),v(nx,ny,1:nz_tot),w(nx,ny,1:nz_tot))
   allocate(Pcon(ld,ny,1:nz_tot,npcon),tot_vol(ld,ny,1:nz_tot,npcon))
   allocate(u(ld,1:ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot))


   allocate(Pcon0(ld,ny,1:nz,npcon))
   allocate(dissip(ld,ny,1:nz_tot))

allocate(d50(310-xps+1),tot(310-jx+1,npcon-1))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)

allocate(weight_x(nx),weight_y(ny),weight_z(nz_tot))



 jt=200
      jt_total=jt*base
    do ip=1,nproc
      write(file,'(A,A,I6.6,A,I3.3,A)')TRIM(dir),'vel_pcon_0',jt_total,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon)
      nzs=(ip-1)*(nz-1)
!      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
!      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
!      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
!      dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
      Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
!      breakage_freq(:,:,nzs+1:nzs+nz,1:npcon) = breakage_freq0(:,:,1:nz,1:npcon)
!!      Re(:,:,nzs+1:nzs+nz,1:npcon) = Re0(:,:,1:nz,1:npcon)
!      Rhs1(:,:,nzs+1:nzs+nz,1:npcon )= Rhs10(:,:,1:nz,1:npcon)
!      Rhs2(:,:,nzs+1:nzs+nz,1:npcon) = Rhs20(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
  print *,"Opened files,timestep=",jt_total




  weight_y = dy/2
  weight_z = dz/2
  weight_y(2:ny-1) = dy
  weight_z(2:nz_tot-1) = dz
!  weight_x = dx/2
!  weight_x(2:ax-1) = dx

! path1 = 'post_proc_output/' 
!   write(file,'(A,A)')TRIM(PATH1), "av_vel_pcon_300_452.out"
!!!   write(file,'(A)') "Pcon_30000.out"
!   open(20,FILE=file,form='unformatted')
!!   read(20) PCon(:,:,1:nz_tot,:)
!   read(20) u(:,:,1:nz_tot),v(:,:,1:nz_tot),w(:,:,1:nz_tot),PCon(:,:,1:nz_tot,:)!,&
!! !           dissip(:,:,1:nz_tot)
!  close(20)
!!   write(*,*) maxval(pcon(:,:,:,1)), maxloc(pcon(:,:,:,1))
!  
     do ip=1,npcon
     tot_vol(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!     write(*,*) maxval(tot_vol(:,:,:,ip))
     enddo
   
   allocate(M_tot(310-xps+1),M_y(310-xps+1),M_z(310-xps+1),y_cm(310-xps+1),z_cm(310-xps+1))
   M_tot=0._rprec
   M_y = 0._rprec 
   M_z = 0._rprec
   y_cm= 0._rprec 
   z_cm= 0._rprec 
    do jx = xps,310      
    do jy= 1,ny
    do jz = 1,nz_tot
        M_tot(jx-xps+1) = M_tot(jx-xps+1) + weight_y(jy)*weight_z(jz)*sum(tot_vol(jx,jy,jz,:))
        M_y(jx-xps+1) =  M_y(jx-xps+1) + (jy-1)*L_y/ny*weight_y(jy)*weight_z(jz)*sum(tot_vol(jx,jy,jz,:))
        M_z(jx-xps+1) =  M_z(jx-xps+1) + (jz-1)*z_i/nz_tot*weight_y(jy)*weight_z(jz)*sum(tot_vol(jx,jy,jz,:))
 
    enddo
    enddo
 
    enddo
    y_cm = M_y/M_tot - (yps-1)*L_y/ny
    z_cm = M_z/M_tot    

path ='post_proc_output/pcon40k/'
!
505 format(18e15.7)
write(file,'(A,A)') TRIM(PATH) , 'y_cm_tot.dat'
open(unit=22,FILE=file,status='unknown')

write(file,'(A,A)') TRIM(PATH) , 'z_cm_tot.dat'
open(unit=23,FILE=file,status='unknown')
write(23,*)'variables=x,z_cm'
write(23,*) 'zone t="',1,'"i=',310-xps+1,'f=point'

write(22,*)'variables=x,y_cm'
write(22,*) 'zone t="',1,'"i=',310-xps+1,'f=point'
   
do jx = xps,310
  write(22,505) (jx-1)*l_x/nx,y_cm(jx-xps+1)
 
  write(23,505) (jx-1)*l_x/nx,(nz_tot-1._rprec)/nz_tot-z_cm(jx-xps+1)
enddo

close(22)


close(23)

!   jx_start = 150
!   jx_end   = 240
!   jy_start = yps-3
!   jy_end   = yps+3
!   jz_start = 168
!   jz_end   = 178
!
!   
!do jx =xps,310
!tot(jx-xps+1,1)= sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1))&
!                  /sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!    do ip=2,npcon-1
!      tot(jx-xps+1,ip) = tot(jx-xps+1,ip-1) +sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,ip))&
!                                       /sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!        if(tot(jx-xps+1,ip) .ge. 0.5_rprec)  then
!          d50(jx-xps+1) = diameter(ip-1) + (diameter(ip)-diameter(ip-1))*(0.5-tot(jx-xps+1,ip-1))&
!                                             /(tot(jx-xps+1,ip)-tot(jx-xps+1,ip-1))
!          exit
!        endif 
!    enddo
!enddo
!
!write(file,'(A,A,I3.3,A)') TRIM(PATH) , "d50_", int((jz_end-5)*100/384),".dat"
!open(unit=24,FILE=file,status='unknown')
!write(24,*)'variables=x,d50'
!write(24,*) 'zone t="',1,'"i=',310-xps+1,'f=point'
!   
!do jx = xps,310
!  write(24,505) (jx-1)*l_x/nx,d50(jx-xps+1)
!enddo
!
!close(24)
!
!
deallocate(u,v,w,pcon,tot_vol,tot,d50)!,breakage_freq)
END PROGRAM POSTPROCESSING
