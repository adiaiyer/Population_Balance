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
  INTEGER,PARAMETER                           :: ny=288
  INTEGER,PARAMETER                           :: nz_tot=385
  INTEGER,PARAMETER                           :: nz = (385-1)/nproc +1
  INTEGER,PARAMETER                           :: nx=288
  INTEGER,PARAMETER                           :: nsteps=200
  INTEGER,PARAMETER                           :: base=350
  INTEGER,PARAMETER                           :: npcon=20
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =215
  INTEGER,PARAMETER                          :: jt_end = 460
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = 2.5/nproc
  REAL(rprec),PARAMETER                      :: L_x=1._rprec
  REAL(rprec),PARAMETER                      :: L_y = 1_rprec
  REAL(rprec),PARAMETER                      :: u_star = 1
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0001
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,dissip,u0,v0,w0,dissip0
  REAL(rprec),dimension(:,:,:,:),allocatable ::  Pcon,breakage_freq,Re,tot_vol,tot_area
  Real(rprec),dimension(:,:,:,:),allocatable :: Pcon0,avg_pcon
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(npcon)               :: dsd
  CHARACTER(len=64)                           ::  file,path,path1i,format_spec,fname
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,reclen,jt,jx,jy,jz,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  real                                        :: vol,half_sum
real(rprec)                       ::               d32,dx=L_x/nx,dy=L_y/ny,dz=L_z/z_i/(nz-1)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,d50
  real(rprec),dimension(:,:),allocatable            ::M_tot,M_y,M_z,y_cm,z_cm,tot
  logical :: file_exists
  real(rprec),parameter::pi=3.1415926535897932384626433_rprec
  dir ='output_data/'
!   allocate(Pcon(nx,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon))
!   allocate(u(1:nx,1:ny,1:nz_tot),v(nx,ny,1:nz_tot),w(nx,ny,1:nz_tot))
!   allocate(breakage_freq(nx,ny,1:nz_tot,npcon),Re(nx,ny,1:nz_tot,npcon))
!   allocate(dissip(nx,ny,1:nz_tot))
format_spec = '(I6,2F14.7)'
   allocate(Pcon0(ld,ny,1:nz,npcon),Pcon(ld,ny,1:nz_tot,npcon),tot_area(nx,ny,1:nz_tot,npcon),tot_vol(ld,ny,1:nz_tot,npcon))
    print *,"hello"
   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),dissip0(ld,ny,1:nz),&
  u(ld,ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot),dissip(ld,ny,1:nz_tot),avg_pcon(nx,ny,1:nz_tot,npcon))
allocate(tot(310-jx+1,npcon-1))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)

allocate(weight_x(nx),weight_y(ny),weight_z(nz_tot))

  weight_y = dy/2
  weight_z = dz/2
  weight_y(2:ny-1) = dy
  weight_z(2:nz_tot-1) = dz
!  weight_x = dx/2
!  weight_x(2:ax-1) = dx

 path = 'post_proc_output/' 
!   write(file,'(A)') "pcon_050000.out"
open(10,file='tot_vol.dat')
!jt= 460
!    do jt=jt_start,jt_end
      jt_total=88000
    do ip=1,nproc
      write(file,'(A,A,I6.6,A,I3.3,A)')TRIM(dir),'vel_pcon_0',jt_total,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip) u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz)!, &
!                  breakage_freq0(:,:,1:nz,1:npcon),Re0(:,:,1:nz,1:npcon),rhs10(:,:,1:nz,1:npcon),rhs20(:,:,1:nz,1:npcon)
      nzs=(ip-1)*(nz-1)
      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)      
      dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
      Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
!  print *,"Opened files,timestep=",jt_total


!inquire(iolength=reclen)u(1:nx,1:ny,1:nz_tot)
!write(fname,'(A,A)')TRIM(PATH), "avg_pcon_342_469.out"
!inquire(FILE = fname,EXIST=file_exists)
!
!if(file_exists) then
!        open(10,file=fname,access='direct',RECL=reclen,form='unformatted')
!        do ip = 1,npcon
!                read(10,rec=ip) pcon(1:nx,1:ny,1:nz_tot,ip)
!        enddo
!        close(10)
!endif

!!   write(file,'(A)') "Pcon_30000.out"
!   open(20,FILE=file,form='unformatted')
!   read(20) PCon(1:nx,1:ny,1:nz_tot,:)
! !  read(20) u(:,:,1:nz_tot),v(:,:,1:nz_tot),w(:,:,1:nz_tot),PCon(:,:,1:nz_tot,:),&
! !           dissip(:,:,1:nz_tot)
!  close(20)

!  print *, maxval(-dissip(55,:,:))*u_star**3,maxloc(-dissip(55,:,:)),&
!            maxval(-dissip(:,:,290:310))*u_star**3,maxloc(-dissip(:,:,290:310)),&
!            maxval(-dissip(145,:,168))*u_star**3,maxloc(-dissip(145,:,168))
!!   write(*,*) maxval(pcon(:,:,:,1)), maxloc(pcon(:,:,:,1))
!  
     do ip=1,npcon
     tot_vol(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!     write(*,*) maxval(tot_vol(:,:,:,ip))
     enddo
     half_sum = sum(tot_vol(1:nx,:,150:zps+2,:))
 !   write(10,format_spec) jt_total, sum(tot_vol),half_sum
!    enddo    
!   close(10)
!   allocate(M_tot(310-xps+1,npcon),M_y(310-xps+1,npcon),M_z(310-xps+1,npcon),y_cm(310-xps+1,npcon),z_cm(310-xps+1,npcon))
!   M_tot=0._rprec
!   M_y = 0._rprec 
!   M_z = 0._rprec
!   y_cm= 0._rprec 
!   z_cm= 0._rprec 
!    do jx = xps,310      
!    do jy= 1,ny
!    do jz = 1,nz_tot
!        M_tot(jx-xps+1,:) = M_tot(jx-xps+1,:) + weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
!        M_y(jx-xps+1,:) =  M_y(jx-xps+1,:) + (jy-1)*L_y/ny*weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
!        M_z(jx-xps+1,:) =  M_z(jx-xps+1,:) + (jz-1)*z_i/nz_tot*weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
! 
!    enddo
!    enddo
! 
!    enddo
!    y_cm = M_y/M_tot - (yps-1)*L_y/ny
!    z_cm = M_z/M_tot    
!
!path ='post_proc_output/'
!!
505 format(23e15.7)
!write(file,'(A,A)') TRIM(PATH) , 'y_cm.dat'
!open(unit=22,FILE=file,status='unknown')
!
!write(file,'(A,A)') TRIM(PATH) , 'z_cm.dat'
!open(unit=23,FILE=file,status='unknown')
!write(23,*)'variables=x,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15'
!write(23,*) 'zone t="',1,'"i=',310-xps+1,'f=point'
!
!write(22,*)'variables=x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15'
!write(22,*) 'zone t="',1,'"i=',310-xps+1,'f=point'
!   
!do jx = xps,310
!  write(22,505) (jx-1)*l_x/nx,(y_cm(jx-xps+1,ip),ip=1,npcon)
! 
!  write(23,505) (jx-1)*l_x/nx,(-z_cm(jx-xps+1,ip),ip=1,npcon)
!enddo
!
!close(22)
!
!
!close(23)
!
!  ! write(*,*)  sum(tot_vol(1:250,:,zps-4,:))*maxval(-w(:,:,zps-4)*u_star), sum(tot_vol(1:250,:,zps-8,:))&
!!*maxval(-w(:,:,zps-8)*u_star) 
   d_half = (diameter(2:npcon)+ diameter(1:npcon-1))/2_rprec
   d_d(1) = (d_half(1)-diameter(1))*2
   d_d(npcon) = (diameter(npcon)-d_half(npcon-1))*2
   d_d(2:npcon-1) = d_half(2:npcon-1) - d_half(1:npcon-2)
   d_d = d_d*1.e6_rprec
   jx_start = xps-2
   jx_end   = xps+2
   jy_start = yps-2
   jy_end   = yps+2
!   jz_start = 168
!   jz_end   = 178

  do ip=1,npcon
!    tot_vol(:,:,:,ip) = pi/6_rprec*pcon(1:nx,:,:,ip)*diameter(ip)**3_rprec
    tot_area(:,:,:,ip) = pi*pcon(1:nx,:,:,ip)*diameter(ip)**2_rprec
  enddo
do  jz=30,300,5 

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

!   jx_start = 166
!   jx_end   = 167
!   jy_start = yps
!   jy_end   = yps
!   jz_start = 172
!   jz_end   = 172
!   counter=1
!
!do jz = 150,300,20
!     d32 =  6 * sum(tot_vol(xps-2:xps+2,yps-2:yps+2,jz,1:npcon)) / &
!                           sum(tot_area(xps-2:xps+2,yps-2:yps+2,jz,1:npcon))
     do ip=1,npcon
      
       dsd(ip) =sum(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))&
                      /size(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))
     enddo
   write(file,'(A,A,A)')TRIM(path), "dsd_z_centerline_2",".dat"
   
!   open(unit=21,FILE='d32_centerline_2.dat',status='unknown',position='append',action='write')
!   write(21,505) (zps-jz)*2.5/384,d32

 !   close(21)
   open(20,FILE=file,status='unknown',position = 'append',action = 'write')
     write(20,505)  (zps-jz)*2.5*100/385,(dsd(ip)/d_d(ip),ip=1,npcon)
! 
   enddo
   close(20)
!do jx =xps,310
!tot(jx-xps+1,1)=sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!    do ip=2,npcon-1
!      tot(jx-xps+1,ip) = tot(jx-xps+1,ip-1)+&
!               sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,ip))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!        if (tot(jx-xps+1,ip) .ge. 0.5_rprec)  then
!          d50(jx-xps+1) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(jx-xps+1,ip-1))/&
!                                                    (tot(jx-xps+1,ip)-tot(jx-xps+1,ip-1))
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
!     do ip=1,npcon
!
!       dsd(ip) = sum(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))&
!                         /size(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))
!     enddo
!
!
!   write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(path), "dsd_z_",jt_total,'_', int((jz)*100/384),".dat"
!
!     open(20,FILE=file)
!     write(20,505) (dsd(ip)/d_d(ip),ip=1,npcon)
!close(20)
!
!jz =270
!
!     do ip=1,npcon
!
!       dsd(ip) = sum(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))&
!                         /size(pcon(jx_start:jx_end,jy_start:jy_end,jz,ip))
!     enddo
!
!
!
!   write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(path), "dsd_z_",jt_total,'_', int((jz)*100/384),".dat"
!     open(20,FILE=file)
!     write(20,505) (dsd(ip)/d_d(ip),ip=1,npcon)
!close(20)
!!do jx =xps,310
!!tot(jx-xps+1,1)=sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!!    do ip=2,npcon-1
!!      tot(jx-xps+1,ip) = tot(jx-xps+1,ip-1) &
!!             +sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,ip))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!!        if (tot(jx-xps+1,ip) .ge. 0.5_rprec)  then
!!          d50(jx-xps+1) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(jx-xps+1,ip-1))&
!!                                               /(tot(jx-xps+1,ip)-tot(jx-xps+1,ip-1))
!!          exit
!!        endif
!!    enddo
!!enddo
!!
!!write(file,'(A,A,I3.3,A)') TRIM(PATH) , "d50_", int((jz_end-5)*100/384),".dat"
!!open(unit=24,FILE=file,status='unknown')
!!write(24,*)'variables=x,d50'
!!write(24,*) 'zone t="',1,'"i=',310-xps+1,'f=point'
!!
!!do jx = xps,310
!!  write(24,505) (jx-1)*l_x/nx,d50(jx-xps+1)
!!enddo
!!
!!close(24)
!
!
!
!
!
!!write(file,'(A,A)')TRIM(PATH),"av_vel.dat"
!!open(21,FILE=file)
!!write(21,*) 'variables=x,y,z,u,v,w'
!!write(21,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!!do jz = 1,nz_tot
!!do jy=1,ny
!!do jx=1,nx
!!
!!  write(21,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,u(jx,jy,jz)*u_star,v(jx,jy,jz)*u_star,&
!!       -w(jx,jy,jz)*u_star
!!enddo
!!enddo
!!enddo
!
!close(21)
!write(file,'(A,A)')TRIM(PATH),"pcon_90400.dat"
!open(21,FILE=file)
!write(21,*) 'variables=x,y,z,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15'
!write(21,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!! jy=yps
!      do jz=1,nz_tot
!          do jy=1,ny
!            do jx=1,nx
!             write(21,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,(pcon(jx,jy,jz,ip),ip=1,npcon)
!           enddo
!         enddo
!        enddo
!
!enddo
!inquire(iolength=reclen)u(1:nx,1:ny,1:nz_tot)
!open(30,file='pcon_161000.out',form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(30,rec=ip) Pcon(1:nx,1:ny,1:nz_tot,ip)
!enddo
!
!close(30)
!deallocate(u,v,w,dissip)
END PROGRAM POSTPROCESSING



!REAL FUNCTION vol_int(intg,ix_l,ix_r,iy_l,iy_r,iz_l,iz_r)
!
!use types, only : rprec
!use param, only : dx,dy,dz
!implicit none
!
!real(rprec), dimesnion(:,:,:,:),allocatable:: intg
!integer ::           ix_l,ix_r,iy_l,iy_r,iz_l,iz_r
!!INEGRATE along X
!integer,parameter:: ax = ix_r-ix_l+1,&
!                    ay = iy_r-iy_l+1,&
!                    az = iz_r-iz_l+1
!real(rprec),dimension(ax,ay) :: weight = dx/2
!real(rprec),dimension(az) :: weight_z = dz/2
!
!weight_xy(2:ix_r-1,2:iy_r-1) = dx
!weight_z(2:iz_r-1)  = dz
!
!do jz=1,az
!do jy=1,ay
!do jx=1,ax
!
! vol_int = vol_int+ weight_xy(jx,jy)*weight(jz)*intg(jx,jy,jz)
!
!enddo
!enddo
!enddo
!
!RETURN
!
!END
!
!
!
!REAL FUNCTION surface_int(intg,il_1,ir_1,il_2,ir_2,plane)
!
!use types, only : rprec
!use param, only : dx,dy,dz
!
!integer :: plane,il_1,ir_1,il_2,ir_2
!integer,parameter:: n1 = ir_1-il_1+1,&
!                    n2 = ir_2-il_2+1
!real(rprec), dimension(n1,n2) :: intg
!!plane =1 x
!!plane =2 y
!!plane =3 z
!
!SELECT CASE (plane)
!  CASE(1) 
!    weight_1 = dy/2
!    weight_2 = dz/2
!    weight_1(2:n1-1) = dy
!    weight_2(2:n1-1) = dz
!  CASE(2)
!    weight_1 = dx/2
!    weight_2 = dz/2
!    weight_1(2:n1-1) = dx
!    weight_2(2:n2-1) = dz
!
!  CASE(3)
!    weight_1 = dx/2
!    weight_2 = dy/2
!    weight_1(2:n1-1) = dx
!    weight_2(2:n2-1) = dy
!  
!  END SELECT
!
!weight = matmul(weight_1,weight_2)
!
!surface_int = sum(weight*intg)
!
!
!END

