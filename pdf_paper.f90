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
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,dissip,tot_area
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,breakage_freq0,Re0,pcon0,tot_vol,&
                                                rhs10,rhs1
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
  real(rprec),dimension(:,:),allocatable::M_tot,M_y,M_z,y_cm,z_cm,d50,rhs_source
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
 integer,parameter :: n=6
  dir ='output_data/'

  allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),pcon0(ld,ny,1:nz,npcon),dissip0(ld,ny,1:nz),& 
            breakage_freq0(ld,ny,1:nz,npcon),Re0(ld,ny,1:nz,npcon),Rhs10(ld,ny,1:nz,npcon)) 

   allocate(Pcon(ld,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon),rhs1(ld,ny,1:nz_tot,npcon))
   allocate(dissip(ld,ny,1:nz_tot),tot_area(nx,ny,1:nz_tot))

   allocate(d50(n,jt_end-jt_start+1),tot(npcon-1),rhs_source(n,npcon))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)
path ='post_proc_output/'
!
505 format(18e15.7)

  open(unit=1001,FILE="d50_140_295_z231.dat",status='unknown',position="append")
  write(1001,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
  open(unit=1002,FILE="d50_140_295_z188.dat",status='unknown',position="append")
  write(1002,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
  open(unit=1003,FILE="d50_140_295_z168.dat",status='unknown',position="append")
  write(1003,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 
  open(unit=1004,FILE="d50_140_295_z208.dat",status='unknown',position="append")
  write(1004,*)'variables=t,d50_1,d50_2,d50_3,d50_4,d50_5,d50_6' 

  open(unit=1005,FILE="ar_eps_140_295_z231.dat",status='unknown',position="append")
  write(1005,*)'variables=t,ar_1,ar_2,ar_3,ar_4,ar_5,ar_6,e_1,e_2,e_3,e_4,e_5,e_6'
  open(unit=1006,FILE="ar_eps_140_295_z188.dat",status='unknown',position="append")
  write(1006,*)'variables=t,ar_1,ar_2,ar_3,ar_4,ar_5,ar_6,e_1,e_2,e_3,e_4,e_5,e_6'
  open(unit=1007,FILE="ar_eps_140_295_z168.dat",status='unknown',position="append")
  write(1007,*)'variables=t,ar_1,ar_2,ar_3,ar_4,ar_5,ar_6,e_1,e_2,e_3,e_4,e_5,e_6'
  open(unit=1008,FILE="ar_eps_140_295_z208.dat",status='unknown',position="append")
  write(1008,*)'variables=t,ar_1,ar_2,ar_3,ar_4,ar_5,ar_6,e_1,e_2,e_3,e_4,e_5,e_6'


  open(unit=1009,FILE="rhs_140_295_z231.dat",status='unknown',position="append")
  write(1009,*)'variables=t,x,rhs_1,rhs_2,rhs_3,rhs_4,rhs_5,rhs_6...15' 
  open(unit=1010,FILE="rhs_140_295_z188.dat",status='unknown',position="append")
  write(1010,*)'variables=t,x,rhs_1,rhs_2,rhs_3,rhs_4,rhs_5,rhs_6...15' 
  open(unit=1011,FILE="rhs_140_295_z168.dat",status='unknown',position="append")
  write(1011,*)'variables=t,x,rhs_1,rhs_2,rhs_3,rhs_4,rhs_5,rhs_6...15' 
  open(unit=1012,FILE="rhs_140_295_z208.dat",status='unknown',position="append")
  write(1012,*)'variables=t,x,rhs_1,rhs_2,rhs_3,rhs_4,rhs_5,rhs_6...15' 


do jt=jt_start,jt_end

  jt_total=jt*base
  do ip=1,nproc
    write(file,'(A,A,I6.6,A,I3.3,A)')TRIM(dir),'vel_pcon_0',jt_total,'_2',ip-1,'.out'
    open(unit=2000+ip,FILE=file,form='unformatted')
    read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz),&
                 breakage_freq0(:,:,1:nz,1:npcon),Re0(:,:,1:nz,1:npcon),rhs10(:,:,1:nz,1:npcon)
    nzs=(ip-1)*(nz-1)
    rhs1(:,:,nzs+1:nzs+nz,1:npcon) = rhs10(:,:,1:nz,1:npcon)
    dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
    Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
    close(2000+ip)
  enddo
  print *,"Opened files,timestep=",jt_total

  do ip=1,npcon
    tot_vol(:,:,:,ip) = pi/6_rprec*pcon(1:nx,:,:,ip)*diameter(ip)**3_rprec
    tot_area(:,:,:) = pi*pcon(1:nx,:,:,ip)*diameter(ip)**2_rprec
  enddo

  jy=yps  
  jz =231
  do i=1,n
    jx=145 + 30*(i-1)
!!!!!!!!!!!!!!! RHS PDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    rhs_source = 0._rprec
    rhs_source(i,1:npcon) =rhs1(jx,jy,jz,1:npcon)/pcon(jx,jy,jz,1:npcon)
    write(1009,505) jt*base*dt_dim,jx*l_x/nx,(rhs_source(i,ip),ip=1,npcon)
    print *, "write success"
!!!!!!!!!!!!!!! D50 PDF  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
    tot = 0._rprec
    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
    do ip=2,npcon-1
      tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
      print *,ip, tot(ip)
      if(tot(ip) .ge. 0.5_rprec)  then
        d50(i,jt-jt_start+1) = diameter(ip-1) + (diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
                                             /(tot(ip)-tot(ip-1)) 
        exit
      endif 
    enddo
  enddo

  write(1001,505) (jt)*base*dt_dim,(d50(i,jt-jt_start+1),i=1,n)
  write(1005,505) jt*base*dt_dim,(tot_area(jx,jy,jz),jx=140,295,30),(-dissip(jx,jy,jz),jx=140,295,30)
  jz=188
  do i =1,n

    jx=145 + 25*(i-1)
!!!!!!!!!!!!!!! RHS PDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    rhs_source = 0._rprec
    rhs_source(i,1:npcon) =rhs1(jx,jy,jz,1:npcon)/pcon(jx,jy,jz,1:npcon)
    write(1010,505) jt*base*dt_dim,jx*l_x/nx,(rhs_source(i,ip),ip=1,npcon)

!!!!!!!!!!!!!!! D50 PDF  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
    tot = 0._rprec
    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
      do ip=2,npcon-1
        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
        if(tot(ip) .ge. 0.5_rprec)  then
          d50(i,jt-jt_start+1) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
                                               /(tot(ip)-tot(ip-1))
          exit
        endif
      enddo

  enddo 

  write(1006,505) jt*base*dt_dim,(tot_area(jx,jy,jz),jx=140,295,30),(-dissip(jx,jy,jz),jx=140,295,30)
  write(1002,505) (jt)*base*dt_dim,(d50(i,jt-jt_start+1),i=1,n)

  jz=168
  do i =1,n

    jx=145 + 25*(i-1)
!!!!!!!!!!!!!!! RHS PDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    rhs_source = 0._rprec
    rhs_source(i,1:npcon) =rhs1(jx,jy,jz,1:npcon)/pcon(jx,jy,jz,1:npcon)
    write(1011,505) jt*base*dt_dim,jx*l_x/nx,(rhs_source(i,ip),ip=1,npcon)

!!!!!!!!!!!!!!! D50 PDF  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
    tot = 0._rprec
    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
      do ip=2,npcon-1
        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
        if(tot(ip) .ge. 0.5_rprec)  then
          d50(i,jt-jt_start+1) = diameter(ip-1)+(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
                                               /(tot(ip)-tot(ip-1))
          exit
        endif
      enddo

  enddo
  write(1003,505) (jt)*base*dt_dim,(d50(i,jt-jt_start+1),i=1,n)

  write(1007,505) jt*base*dt_dim,(tot_area(jx,jy,jz),jx=140,295,30),(-dissip(jx,jy,jz),jx=140,295,30)

  jz=208
  do i =1,n

    jx=145 + 25*(i-1)
!!!!!!!!!!!!!!! RHS PDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    rhs_source = 0._rprec
    rhs_source(i,1:npcon) =rhs1(jx,jy,jz,1:npcon)/pcon(jx,jy,jz,1:npcon)
    write(1012,505) jt*base*dt_dim,jx*l_x/nx,(rhs_source(i,ip),ip=1,npcon)

!!!!!!!!!!!!!!! D50 PDF  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
    tot = 0._rprec
    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
      do ip=2,npcon-1
        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
        if(tot(ip) .ge. 0.5_rprec)  then
          d50(i,jt-jt_start+1) = diameter(ip-1)+(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
                                               /(tot(ip)-tot(ip-1))
          exit
        endif
      enddo

  enddo
  write(1004,505) (jt)*base*dt_dim,(d50(i,jt-jt_start+1),i=1,n)
  write(1008,505) jt*base*dt_dim,(tot_area(jx,jy,jz),jx=140,295,30),(-dissip(jx,jy,jz),jx=140,295,30)

enddo
  close(1001)
  close(1002)
  close(1003)
  close(1004)
  close(1005)
  close(1006)
  close(1007)
  close(1008)
  close(1009)
  close(1010)
  close(1011)
  close(1012)



deallocate(u0,v0,w0,dissip0,pcon0,tot_vol,tot,d50,pcon,dissip,rhs1,tot_Area,Re0,breakage_freq0)
END PROGRAM POSTPROCESSING




