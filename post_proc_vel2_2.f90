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
  INTEGER,PARAMETER                           :: npcon=20
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =300
  INTEGER,PARAMETER                          :: jt_end = 400
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2._rprec
  REAL(rprec),PARAMETER                      :: L_y = 0.78125_rprec
  REAL(rprec),PARAMETER                      :: u_star = 1.
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,dissip
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,breakage_freq,Re,tot_vol,pcon0
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(:,:),allocatable               :: dsd
  CHARACTER(len=64)                           :: file,path,path1
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  integer ::                       count_n=3
  real                                        :: vol,rat 
real(rprec)                       ::               dx=L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,d50
  real(rprec),dimension(:,:),allocatable            ::M_tot,M_y,M_z,y_cm,z_cm,tot
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
  dir ='output_data/'
   allocate(Pcon(ld,ny,1:nz_tot,npcon),tot_vol(ld,ny,1:nz_tot,npcon))
   allocate(u(ld,1:ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot))
   allocate(breakage_freq(ld,ny,1:nz_tot,npcon),Re(ld,ny,1:nz_tot,npcon))
   allocate(dissip(ld,ny,1:nz_tot))
allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),pcon0(ld,ny,1:nz,npcon),dissip0(ld,ny,1:nz))
allocate(d50(310-xps+1),tot(310-jx+1,npcon-1))

diameter(1) = 20e-6_rprec
diameter(npcon) = 1200e-6_rprec
rat = (60._rprec)**(1._rprec/19._rprec)
!print *, d_d
do ip=2,npcon-1
   diameter(ip) = diameter(1)*rat**(ip-1)
enddo
!   open(21,File='diameter.dat')
!   do ip=1,npcon
!     read(21,*) diameter(ip)
!   enddo
!   close(21)

allocate(weight_x(nx),weight_y(ny),weight_z(nz_tot))
allocate(dsd(count_n,npcon))
  weight_y = dy/2
  weight_z = dz/2
  weight_y(2:ny-1) = dy
  weight_z(2:nz_tot-1) = dz
!  weight_x = dx/2
!  weight_x(2:ax-1) = dx

! path = 'post_proc_output/av_375_650/' 
!   write(file,'(A,A)')TRIM(PATH), "av_rhs_pcon_0375_0650.out"
!   write(file,'(A)') "Pcon_30000.out"
!   open(20,FILE=file,form='unformatted')
!   read(20) PCon(1:nx,1:ny,1:nz_tot,:)
 !  read(20) u(:,:,1:nz_tot),v(:,:,1:nz_tot),w(:,:,1:nz_tot),PCon(:,:,1:nz_tot,:),&
 !           dissip(:,:,1:nz_tot)
 ! close(20)


    jt=40
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
       PCon(:,:,nzs+1:nzs+nz,:)=Pcon0(:,:,1:nz,:)
       dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
       close(2000+ip)
    enddo



!  print *, maxval(-dissip(55,:,:))*u_star**3,maxloc(-dissip(55,:,:)),&
!            maxval(-dissip(:,:,290:310))*u_star**3,maxloc(-dissip(:,:,290:310)),&
!            maxval(-dissip(145,:,168))*u_star**3,maxloc(-dissip(145,:,168))
!!   write(*,*) maxval(pcon(:,:,:,1)), maxloc(pcon(:,:,:,1))
!  
     do ip=1,npcon
     tot_vol(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!!     write(*,*) maxval(tot_vol(:,:,:,ip))
     enddo

505 format(23e15.7e2)

!open(unit=22,file='pcon_init.dat',status='unknown')
!write(22,505) ((tot_vol(59,yps,271,ip)*(1.9e-3_rprec/60_rprec)/sum(tot_vol(59,yps,271,:)))*&
!                 (6._rprec/pi/diameter(ip)**3_rprec),ip=1,npcon)
!close(22)
   
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

path ='post_proc_output/'
!
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
!  write(23,505) (jx-1)*l_x/nx,((nz_tot-1._rprec)/nz_tot-z_cm(jx-xps+1,ip),ip=1,npcon)
!enddo
!
!close(22)
!
!
!close(23)

  ! write(*,*)  sum(tot_vol(1:250,:,zps-4,:))*maxval(-w(:,:,zps-4)*u_star), sum(tot_vol(1:250,:,zps-8,:))&
!*maxval(-w(:,:,zps-8)*u_star) 
   d_half = (diameter(2:npcon)+ diameter(1:npcon-1))/2_rprec
   d_d(1) = (d_half(1)-diameter(1))*2
   d_d(npcon) = diameter(npcon)-diameter(npcon-1)
   d_d(2:npcon-1) = d_half(2:npcon-1) - d_half(1:npcon-2)
   d_d = d_d*1.e6_rprec
!   jx_start = 150
!   jx_end   = 170
!   jy_start = yps-3
!   jy_end   = yps+3
!   jz_start = 168
!   jz_end   = 178

  
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


!   jx_start = 150
   jx_start = 59
!   jx_end   = 170
   jy_start = yps-3
   jy_end   = yps+3
!   jz_start = 168
!   jz_end   = 178
   jz = 271
   counter=1
   print *, pcon(jx_start,yps,jz,1)
!   print *, d_d
!   do jx=jx_start,jx_end,10
     do ip=1,npcon
      
       dsd(counter,ip) =(pcon(jx_start,yps,jz,ip))!&
!                      /size(pcon(jx,yps,jz_start:jz,ip))
   enddo
     
!       write(*,*) dsd(counter,1), dsd(counter,4) , dsd(counter,8) ,pcon(jx_start,yps,jz_end-4,1)
 !  counter=counter+1
 !  enddo
!path = 'post_proc_output/av_186_301/'
  ! write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz_end-5)*100/384),".dat"
   write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz-1)*100_rprec/nz_tot),".dat"
   
   open(20,FILE=file)
   do counter=1,count_n
     write(20,505) (jx_start+(counter-1)*10)*L_x/nx ,sum(pi/6._rprec*dsd(counter,:)*diameter**3),&
                   (dsd(counter,ip)/d_d(ip),ip=1,npcon)
   enddo
 
   close(20)
!   counter=1
!!   jz_start = 188
! !  jz_end   = 198
!    jz = 193
!!do jx =xps,310
!!tot(jx-xps+1,1)=sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!!    do ip=2,npcon-1
!!      tot(jx-xps+1,ip) = tot(jx-xps+1,ip-1)+&
!!               sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,ip))/sum(tot_vol(jx,jy_start:jy_end,jz_start:jz_end,1:npcon-1))
!!        if (tot(jx-xps+1,ip) .ge. 0.5_rprec)  then
!!          d50(jx-xps+1) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(jx-xps+1,ip-1))/&
!!                                                    (tot(jx-xps+1,ip)-tot(jx-xps+1,ip-1))
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
!!
!
!
!   do jx=jx_start,jx_end,10
!     do ip=1,npcon
!
!       dsd(counter,ip) = (pcon(jx,yps,jz,ip))!&
!!                         /size(pcon(jx,yps,jz,ip))
!     enddo
!
!   counter=counter+1
!   enddo
!
!!   write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz_end-5)*100/384),".dat"
!   write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz-1)*100._rprec/nz_tot),".dat"
!
!     open(20,FILE=file)
!   do counter=1,count_n
!     write(20,505) (jx_start+(counter-1)*10)*L_x/nx,sum(pi/6._rprec*dsd(counter,:)*diameter**3),&
!                  (dsd(counter,ip)/d_d(ip),ip=1,npcon)
!   enddo
!close(20)
!
!   counter=1
!!   jz_start = 158
!!   jz_end   = 168
!    jz = 163
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
!
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
!   do jx=jx_start,jx_end,10
!     do ip=1,npcon
!
!       dsd(counter,ip) = (pcon(jx,yps,jz,ip))!&
!                 !        /size(pcon(jx,yps,jz,ip))
!     enddo
!
!   counter=counter+1
!   enddo
!
!!   write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz_end-5)*100/384),".dat"
!   write(file,'(A,A,I3.3,A)')TRIM(path), "dsd_z_", int((jz-1)*100_rprec/nz_tot),".dat"
!
!     open(20,FILE=file)
!   do counter=1,count_n
!     write(20,505) (jx_start+(counter-1)*10)*L_x/nx,sum(pi/6._rprec*dsd(counter,:)*diameter**3),&
!                  (dsd(counter,ip)/d_d(ip),ip=1,npcon)
!   enddo
!close(20)

!write(file,'(A,A)')TRIM(PATH),"av_vel.dat"
!open(21,FILE=file)
!write(21,*) 'variables=x,y,z,u,v,w'
!write(21,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!do jz = 1,nz_tot
!do jy=1,ny
!do jx=1,nx
!
!  write(21,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,u(jx,jy,jz)*u_star,v(jx,jy,jz)*u_star,&
!       -w(jx,jy,jz)*u_star
!enddo
!enddo
!enddo
!
!close(21)
deallocate(u,v,w,dissip,Re,tot_vol,tot,breakage_freq)
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
!close(21)
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

