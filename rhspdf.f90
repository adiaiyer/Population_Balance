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
  INTEGER,PARAMETER                          :: jt_start =100
  INTEGER,PARAMETER                          :: jt_end = 460
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
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,breakage_freq0,pcon0,tot_vol,rhs1,rhs10,re0
  Real(rprec),dimension(npcon)               :: diameter,d_d
  Real(rprec),dimension(npcon-1)               :: d_half
  Real(rprec),dimension(10,npcon)               :: dsd,rhs_source
  CHARACTER(len=64)                           :: file,path,path1
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,i,jt_total,nzs,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2,zps=329
  integer::                        jy_start,jy_end,jz_start,jz_end,jx_start,jx_end
  real                                        :: vol,t_vol
real(rprec)                       ::               dx=L_x/nx,dy=L_y/ny,dz=z_i/(nz_tot)
  real(rprec),dimension(:),allocatable        :: weight_x,weight_y,weight_z,tot
  real(rprec),dimension(:,:),allocatable::  M_tot,M_y,M_z,y_cm,z_cm,d50
real(rprec),parameter::pi=3.1415926535897932384626433_rprec
 integer,parameter :: n=4
  dir ='output_data/'

  allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),pcon0(ld,ny,1:nz,npcon),dissip0(ld,ny,1:nz),&
            breakage_freq0(ld,ny,1:nz,npcon),Re0(ld,ny,1:nz,npcon),Rhs10(ld,ny,1:nz,npcon)) 
   allocate(Pcon(ld,ny,1:nz_tot,npcon),tot_vol(nx,ny,1:nz_tot,npcon),Rhs1(ld,ny,1:nz_tot,npcon))
   allocate(dissip(ld,ny,1:nz_tot))

   allocate(d50(n,jt_end-jt_start+1),tot(npcon-1))
   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
   close(21)
path ='post_proc_output/rhs_pdf/'
!
505 format(18e15.7)

do jt=jt_start,jt_end
  jt_total=jt*base
 do ip=1,nproc
      write(file,'(A,A,I5.5,A,I3.3,A)')TRIM(dir),'vel_pcon_00',jt_total,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip)u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz),&
                  breakage_freq0(:,:,1:nz,1:npcon),Re0(:,:,1:nz,1:npcon),rhs10(:,:,1:nz,1:npcon)
      nzs=(ip-1)*(nz-1)
      pcon(:,:,nzs+1:nzs+nz,1:npcon )= pcon0(:,:,1:nz,1:npcon)
      Rhs1(:,:,nzs+1:nzs+nz,1:npcon )= Rhs10(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo


  print *,"Opened files,timestep=",jt_total

  do ip=1,npcon
    tot_vol(:,:,:,ip) = pi/6_rprec*pcon(1:nx,:,:,ip)*diameter(ip)**3_rprec
  enddo

jx_start = 150
jx_end   = 240
jy_start = yps-3
jy_end   = yps+3
jz_start = 168
jz_end   = 178
counter=1

  do jx=jx_start,jx_end,10
    do ip=1,npcon

      rhs_source(counter,ip) =sum(rhs1(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip)/&
                                   pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))&
                               /size(pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))

    enddo
    write(file,'(A,A,I3.3,A)') TRIM(PATH) , "rhs",int((jz_end-5)*100/384),".dat"
    open(unit=1000+counter,FILE=file,status='unknown',access='append')
    write(1000+counter,505) (jt-1)*base*dt_dim,(jx-1)*l_x/nx,(rhs_source(counter,ip),ip=1,npcon)
    close(1000+counter)
    counter=counter+1
 
  enddo
  jz_start = 188
  jz_end   = 198
  counter=1

  do jx=jx_start,jx_end,10
    do ip=1,npcon

      rhs_source(counter,ip)=sum(rhs1(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip)/&
                                   pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))&
                               /size(pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))

    enddo
    write(file,'(A,A,I3.3,A)') TRIM(PATH) , "rhs",int((jz_end-5)*100/384),".dat"
    open(unit=1000+counter,FILE=file,status='unknown',access='append')
    write(1000+counter,505)(jt-1)*base*dt_dim,(jx-1)*l_x/nx,(rhs_source(counter,ip),ip=1,npcon)
    close(1000+counter)
    counter=counter+1

  enddo
  counter=1
  jz_start = 158
  jz_end   = 168
  do jx=jx_start,jx_end,10
    do ip=1,npcon

      rhs_source(counter,ip)=sum(rhs1(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip)/&
                                   pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))&
                               /size(pcon(jx:jx+10,jy_start:jy_end,jz_start:jz_end,ip))

    enddo
    write(file,'(A,A,I3.3,A)') TRIM(PATH) , "rhs",int((jz_end-5)*100/384),".dat"
    open(unit=1000+counter,FILE=file,status='unknown',access='append')
    write(1000+counter,505)(jt-1)*base*dt_dim,(jx-1)*l_x/nx,(rhs_source(counter,ip),ip=1,npcon)
    close(1000+counter)
    counter=counter+1

  enddo
enddo



deallocate(u0,v0,w0,dissip0,pcon0,tot_vol,tot,rhs1,pcon,dissip)
END PROGRAM POSTPROCESSING




