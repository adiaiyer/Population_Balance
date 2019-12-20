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
  INTEGER,PARAMETER                           :: npcon=16
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start =400
  INTEGER,PARAMETER                          :: jt_end = 600
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2
  REAL(rprec),PARAMETER                      :: L_y = 0.78125
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  INTEGER::ip,jt,jx,jy,jz,jt_total,nzs,ipcon
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2


 REAL(rprec),dimension(:,:,:),allocatable   ::u,v,w,dissip

 REAL(rprec),dimension(:,:,:,:),allocatable   ::breakage_freq,Re,pcon
  CHARACTER(len=64)                           :: file
  CHARACTER(len=64)                           :: dir
505 format(18e15.7)
allocate(pcon(ld,ny,1:nz_tot,npcon))!breakage_freq(ld,ny,1:nz_tot,npcon),Re(ld,ny,1:nz_tot,npcon))
open(unit=20,file='Pcon_10000.out',form='unformatted')
read(20)pcon(:,:,1:nz_tot,1:npcon-1) 
close(20)
open(20,file='Pcon_10000.dat')
write(20,*)'variables=x,y,z,br1,br2,re3,re4,re5,re6,re7,re8,re9,re10,re11,re12,re13,re14,re15'
write(20,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
      do jz=1,nz_tot
          do jy=1,ny
            do jx=1,nx
               write(20,505)(jx-1)*L_x/nx,-(jz-1)*z_i/nz_tot,(pcon(jx,yps,jz,ipcon),ipcon=1,npcon-1)
            enddo
          enddo
       enddo
close(20)
deallocate(pcon)

END PROGRAM POSTPROCESSING
