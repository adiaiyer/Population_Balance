PROGRAM READ_DATA
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
  INTEGER,PARAMETER                          :: jt_start =1
  INTEGER,PARAMETER                          :: jt_end = 30
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2
  REAL(rprec),PARAMETER                      :: L_y = 0.78125
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0001
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0
  REAL(rprec),dimension(:,:,:),allocatable   :: u,v,w,avg_u,avg_v,avg_w
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon0,breakage_freq0,Re0
  REAL(rprec),dimension(:,:,:,:),allocatable :: Pcon,avg_Pcon
  CHARACTER(len=64)                           :: file,path,plane
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,plane_cut,counter
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2

  dir ='output_data/'
    print *,"hello"
   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),pcon0(ld,ny,1:nz,npcon))
   allocate(u(ld,ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot),avg_w(ld,ny,1:nz_tot),pcon(ld,ny,1:nz_tot,npcon))
