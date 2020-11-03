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
  INTEGER,PARAMETER                           :: nz = (nz_tot-1)/nproc +1
  INTEGER,PARAMETER                           :: nx=288
  INTEGER,PARAMETER                           :: nsteps=200
  INTEGER,PARAMETER                           :: base=5
  INTEGER,PARAMETER                           :: npcon=1
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =1
  INTEGER,PARAMETER                          :: jt_end = 4
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=2.5_rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=1.
  REAL(rprec),PARAMETER                      :: L_y = 1.
  REAL(rprec),PARAMETER                      :: u_star = 1.
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  real(rprec)                                :: junk
  REAL(rprec),dimension(:,:,:),allocatable   :: txx0,tyy0,tzz0,txz0,tyz0,txy0
  REAL(rprec),dimension(:,:,:),allocatable   :: txx,tyy,tzz,txz,tyz,txy
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon0,breakage_freq0,Re0,Pcon20
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,avg_Pcon,pcon2
  CHARACTER(len=64)                           :: filename
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2
  integer(kind=8)  :: reclen
  dir ='output/'
!    print *,"hello"
   allocate(txx0(nx,ny,1:nz),tyy0(nx,ny,1:nz),tzz0(nx,ny,1:nz))
   allocate(txz0(nx,ny,1:nz),txy0(nx,ny,1:nz),tyz0(nx,ny,1:nz))
   allocate(txx(nx,ny,1:nz_tot),tyy(nx,ny,1:nz_tot),tzz(nx,ny,1:nz_tot))
   allocate(txz(nx,ny,1:nz_tot),txy(nx,ny,1:nz_tot),tyz(nx,ny,1:nz_tot))

  inquire(iolength=reclen) txx0
    do ip=1,nproc
!ip =42
       write(filename,'(A,A,A,I3.3,A)')TRIM(dir),'tau_avg','.c',ip-1,'.bin'
       open(unit=2000+ip,FILE=filename,form='unformatted',access='direct',action='read',recl=reclen)
       read(2000+ip,rec=1) txx0(:,:,1:nz)
       read(2000+ip,rec=2) txy0(:,:,1:nz)
       read(2000+ip,rec=3) tyy0(:,:,1:nz)
       read(2000+ip,rec=4) txz0(:,:,1:nz)
       read(2000+ip,rec=5) tyz0(:,:,1:nz)
       read(2000+ip,rec=6) tzz0(:,:,1:nz)
       nzs=(ip-1)*(nz-1)
       txx(:,:,nzs+1:nzs+nz) = txx0(:,:,1:nz)
       tyy(:,:,nzs+1:nzs+nz) = tyy0(:,:,1:nz)
       tzz(:,:,nzs+1:nzs+nz) = tzz0(:,:,1:nz)
       txy(:,:,nzs+1:nzs+nz) = txy0(:,:,1:nz)
       txz(:,:,nzs+1:nzs+nz) = txz0(:,:,1:nz)
       tyz(:,:,nzs+1:nzs+nz) = tyz0(:,:,1:nz)
      close(2000+ip)


    enddo
  dir = 'post_proc_output/'

  write(filename,'(A,A)') TRIM(dir), 'subgrid_stress.out'
  inquire(iolength=reclen) txx
    open(unit =20, file = filename, form = 'unformatted', &
            access = 'direct', recl = reclen)

    write(20,rec=1) txx(1:nx,:,:)
    write(20,rec=2) tyy(1:nx,:,:)
    write(20,rec=3) tzz(1:nx,:,:)
    write(20,rec=4) txy(1:nx,:,:)
    write(20,rec=5) txz(1:nx,:,:)
    write(20,rec=6) tyz(1:nx,:,:)

    close(20)

 !   505 format(18e15.7)


END PROGRAM POSTPROCESSING










