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
  INTEGER,PARAMETER                          :: jt_start =500
  INTEGER,PARAMETER                          :: jt_end = 650
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2._rprec
  REAL(rprec),PARAMETER                      :: L_y = 0.78125_rprec
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,dissip
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,breakage_freq,pcon0,tot_vol,tot_area
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(10,npcon)               :: dsd
  CHARACTER(len=64)                           :: file,path,path1
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,i,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  real                                        :: vol,t_vol
real(rprec)                       ::               dx=L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,tot
  real(rprec),dimension(:,:),allocatable::  M_tot,M_y,M_z,y_cm,z_cm,d50,d32
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
 integer,parameter :: n=6
  dir ='output_data/'

  allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),pcon0(ld,ny,1:nz,npcon),dissip0(ld,ny,1:nz)) 
   allocate(Pcon(ld,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon),tot_area(nx,ny,1:nz_tot,npcon))
   allocate(dissip(ld,ny,1:nz_tot))

   allocate(d50(n,jt_end-jt_start+1),tot(npcon-1))
   allocate(d32(n,jt_end-jt_start+1))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)
path ='post_proc_output/'
!
505 format(18e15.7)

!  open(unit=1001,FILE="d50_140_295_z231.dat",status='unknown',position="append")
!  write(1001,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
!  open(unit=1002,FILE="d50_140_295_z188.dat",status='unknown',position="append")
!  write(1002,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
!  open(unit=1003,FILE="d50_140_295_z168.dat",status='unknown',position="append")
!  write(1003,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
!  open(unit=1004,FILE="d50_140_295_z208.dat",status='unknown',position="append")
!  write(1004,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 

  open(unit=1001,FILE="d32_140_295_z231.dat",status='unknown',position="append")
  write(1001,*)'variables=t,d32_1,d32_2,d32_3,d32_4,d32_5,d32_6'
  open(unit=1002,FILE="d32_140_295_z188.dat",status='unknown',position="append")
  write(1002,*)'variables=t,d32_1,d32_2,d32_3,d32_4,d32_5,d32_6'
  open(unit=1003,FILE="d32_140_295_z168.dat",status='unknown',position="append")
  write(1003,*)'variables=t,d32_1,d32_2,d32_3,d32_4,d32_5,d32_6'
  open(unit=1004,FILE="d32_140_295_z208.dat",status='unknown',position="append")
  write(1004,*)'variables=t,d32_1,d32_2,d32_3,d32_4,d32_5,d32_6'


do jt=jt_start,jt_end
!  jt = 134
  jt_total=jt*base
  do ip=1,nproc
    write(file,'(A,A,I6.6,A,I3.3,A)')TRIM(dir),'vel_pcon_0',jt_total,'_2',ip-1,'.out'
    open(unit=2000+ip,FILE=file,form='unformatted')
    read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz)
    nzs=(ip-1)*(nz-1)
    dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
    Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
    close(2000+ip)
  enddo
  print *,"Opened files,timestep=",jt_total

  do ip=1,npcon
    tot_vol(:,:,:,ip) = pi/6_rprec*pcon(1:nx,:,:,ip)*diameter(ip)**3_rprec
    tot_area(:,:,:,ip) = pi*pcon(1:nx,:,:,ip)*diameter(ip)**2_rprec
  enddo

    jy=yps  
    jz =231
  do i=1,n

    jx=145 + 30*(i-1)
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        print *,ip, tot(ip)
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(i,jt-jt_start+1) = diameter(ip-1) + (diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1)) 
!!          write(*,*)i,ip
!          exit
!        endif 
         d32(i,jt-jt_start+1) =6*sum(tot_vol(jx,jy,jz,1:npcon))/&
                 sum(tot_area(jx,jy,jz,1:npcon))
!      enddo
  
  enddo
!  write(file,'(A,A,I3.3,A)') TRIM(PATH) ,"d50_145_220_z_",int((jz-1)*100/384),".dat"
!  open(unit=1000+n,FILE=file,status='unknown',position="append")
!  write(1000+n,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
  write(1001,505) (jt)*base*dt_dim,(d32(i,jt-jt_start+1),i=1,n)
!  close(1000+n)

  jz=188
  do i =1,n

    jx=145 + 25*(i-1)
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(i,jt-jt_start+1) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!          exit
!        endif
      d32(i,jt-jt_start+1) =6*sum(tot_vol(jx,jy,jz,1:npcon))/&
                 sum(tot_area(jx,jy,jz,1:npcon))
 !     enddo

  enddo
  write(1002,505) (jt)*base*dt_dim,(d32(i,jt-jt_start+1),i=1,n)
  jz=168
  do i =1,n

    jx=145 + 25*(i-1)
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(i,jt-jt_start+1) = diameter(ip-1)+(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!          exit
!        endif
       d32(i,jt-jt_start+1) =6*sum(tot_vol(jx,jy,jz,1:npcon))/&
                 sum(tot_area(jx,jy,jz,1:npcon))
  !    enddo

  enddo
  write(1003,505) (jt)*base*dt_dim,(d32(i,jt-jt_start+1),i=1,n)


jz=208
  do i =1,n

    jx=145 + 25*(i-1)
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(i,jt-jt_start+1) = diameter(ip-1)+(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!          exit
!        endif
      d32(i,jt-jt_start+1) =6*sum(tot_vol(jx,jy,jz,1:npcon))/&
                 sum(tot_area(jx,jy,jz,1:npcon))
   !   enddo

  enddo
  write(1004,505) (jt)*base*dt_dim,(d32(i,jt-jt_start+1),i=1,n)


enddo
  close(1001)
  close(1002)
  close(1003)
  close(1004)



deallocate(u0,v0,w0,dissip0,pcon0,tot_vol,tot,d50,pcon,dissip)
END PROGRAM POSTPROCESSING




