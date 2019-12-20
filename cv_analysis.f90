module types
implicit none

public
integer, parameter :: rprec = kind (1.d0)
end module types








PROGRAM CV_ANALYSIS

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
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125_rprec 
  REAL(rprec),PARAMETER                      :: ucross =0.15_rprec/u_star
  REAL(rprec),PARAMETER                      :: dt = 0.0002*u_star/z_i
!  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,dissip
 REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,tot_vol
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,breakage_freq,Pcon0,tot_vol1
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(3,npcon)               :: dsd
  CHARACTER(len=64)                           :: file,path
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        ix_l,ix_r,iy_l,iy_r,iz_l,iz_r,ax,ay,az
  real(rprec)                        ::vol_int=0._rprec,vol_int2=0._rprec,u_av,srf_int_t=0_rprec,&
                                       srf_int_l=0._rprec,srf_int_r=0._rprec,srf_int_b=0._rprec
  real(rprec)                       ::               dx =L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z   
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
real(rprec) ::                Q_src,residual



  dir ='output_data/'
   allocate(Pcon(ld,ny,nz_tot,npcon),tot_vol1(ld,ny,nz_tot,npcon))
   allocate(u(ld,ny,nz_tot),v(ld,ny,nz_tot),w(ld,ny,nz_tot))
   allocate(breakage_freq(nx,ny,1:nz_tot,npcon),Pcon0(ld,ny,nz,npcon))
   allocate(dissip(nx,ny,1:nz_tot),u0(ld,ny,nz),v0(ld,ny,nz),w0(ld,ny,nz))
   allocate(tot_vol(nx,ny,nz_tot))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)

q_src =  (6.047e4_rprec/(dx*z_i*dy*z_i*dz*z_i))/(1_rprec/(z_i/u_star))  
ix_l = 160
ix_r = 210
iy_l = 25
iy_r = 125
iz_l = 100
iz_r = 280


  ax = ix_r-ix_l+1
  ay = iy_r-iy_l+1
  az = iz_r-iz_l+1

allocate(weight_x(ax),weight_y(ay),weight_z(az))

  weight_y = dy/2
  weight_z = dz/2
  weight_y(2:ay-1) = dy
  weight_z(2:az-1) = dz
  weight_x = dx/2
  weight_x(2:ax-1) = dx


jt = 202
    do ip=1,nproc
      write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(dir),'vel_pcon_00',jt*base,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip) u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon)!,dissip0(:,:,1:nz), &
      nzs=(ip-1)*(nz-1)
      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)      
      Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
  print *,"Opened files,timestep=",jt*base
!write(file,'(A)') "av_vel_pcon_186_301.out"
!   open(20,FILE=file,form='unformatted')
!   read(20)u(1:nx,1:ny,1:nz_tot),v(1:nx,1:ny,1:nz_tot),w(1:nx,1:ny,1:nz_tot),PCon(1:nx,1:ny,1:nz_tot,1:npcon),&
!            dissip(1:nx,1:ny,1:nz_tot)
!   close(20)


     do ip=1,npcon
     tot_vol1(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!     write(*,*) maxval(tot_vol1(:,:,:,ip))
     enddo
!     write(*,*) maxval(tot_vol1),maxval(Pcon(:,:,:,10))
   tot_vol = sum(tot_vol1(1:nx,1:ny,1:nz_tot,:),4)
   write(*,*)"size = ", maxval(tot_vol)
write(*,*) "min u", minval(ucross+u)
  u_av = ucross + sum(u(ix_r,iy_l:iy_r,iz_l:iz_r))/size(u(ix_r,iy_l:iy_r,iz_l:iz_r))
  do jz=iz_l,iz_r
  do jy = iy_l,iy_r
  
    srf_int_r = srf_int_r +   weight_y(jy-iy_l+1)*weight_z(jz-iz_l+1)*tot_vol(ix_r,jy,jz)*(ucross+u(ix_r,jy,jz))
!     write(*,*) tot_vol(ix_r,jy,jz)*(ucross+u(ix_r,jy,jz)*u_star)
    srf_int_l = srf_int_l +   weight_y(jy-iy_l+1)*weight_z(jz-iz_l+1)*tot_vol(ix_r,jy,jz)*(-(ucross+u(ix_l,jy,jz)))

  enddo
  enddo

write(*,*)"weight" , maxval(weight_y),maxval(weight_z)
  write(*,*) "surf_right =" , srf_int_r
  do jy =iy_l,iy_r
  do jx = ix_l,ix_r
  
    srf_int_t = srf_int_t +weight_y(jy-iy_l+1)*weight_x(jx-ix_l+1)*tot_vol(jx,jy,iz_l)*(-w(jx,jy,iz_l))
   
    srf_int_b = srf_int_b +weight_y(jy-iy_l+1)*weight_x(jx-ix_l+1)*tot_vol(jx,jy,iz_r)*(w(jx,jy,iz_r))
enddo
enddo
  write(*,*) "surface_top=" ,srf_int_t
  write(*,*) "maxw" , maxval(-w(:,:,iz_r)),minval(-w(:,:,iz_r))
do jz=iz_l,iz_r
do jy=iy_l,iy_r
do jx=ix_l,ix_r

 vol_int = vol_int+ weight_x(jx-ix_l+1)*weight_y(jy-iy_l+1)*weight_z(jz-iz_l+1)*tot_vol(jx,jy,jz)
enddo
enddo
enddo

write(*,*) "vol_int",vol_int
    do ip=1,nproc
      write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(dir),'vel_pcon_00',(jt+1)*base,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon)!,dissip0(:,:,1:nz),&
!                  breakage_freq0(:,:,1:nz,1:npcon), Re0(:,:,1:nz,:)
      nzs=(ip-1)*(nz-1)
      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
      Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
   write(*,*)"Opened file", (jt+1)*base
     do ip=1,npcon
     tot_vol1(:,:,:,ip) = pi/6_rprec*pcon(:,:,:,ip)*diameter(ip)**3_rprec
!     write(*,*) maxval(tot_vol(:,:,:,ip))
     enddo
 tot_vol = sum(tot_vol1(1:nx,1:ny,1:nz_tot,:),4)
     do jz=iz_l,iz_r
do jy=iy_l,iy_r
do jx=ix_l,ix_r

 vol_int2 = vol_int2+weight_x(jx-ix_l+1)*weight_y(jy-iy_l+1)*weight_z(jz-iz_l+1)*tot_vol(jx,jy,jz)
enddo
enddo
enddo
write (*,*) "vol_2=" ,vol_int2
residual = (vol_int2-vol_int) + dt*base*(srf_int_t + srf_int_r +srf_int_b+srf_int_l) !- pi/6._rprec*diameter(npcon)**3*Q_src*dx*dy*dz)

write(*,*) "surface",  (srf_int_t + srf_int_r +srf_int_b+srf_int_l) 

open(10,file="residual.dat")
write(10,*)ix_l,ix_r,iy_l,iy_r,iz_l,iz_r,vol_int2,vol_int,srf_int_t,srf_int_r,srf_int_l,residual
close(10)

deallocate(u,v,w,pcon,dissip)

END PROGRAM CV_ANALYSIS



!REAL FUNCTION vol_int(intg,ix_l,ix_r,iy_l,iy_r,iz_l,iz_r)
!
!use types, only : rprec
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

