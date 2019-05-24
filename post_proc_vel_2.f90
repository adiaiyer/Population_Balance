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
  real(rprec)                                :: junk
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0
  REAL(rprec),dimension(:,:,:),allocatable   ::u,v,w,avg_u,avg_v,avg_w,dissip,avg_dissip
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon0,breakage_freq0,Re0
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,avg_Pcon,breakage_freq,avg_breakage_freq
  CHARACTER(len=64)                           :: file
  CHARACTER(len=64)                           :: dir,path
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon
  logical                                     :: exist
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2

  dir ='output_data/'
!   allocate(Pcon0(ld,ny,1:nz,npcon),Pcon(ld,ny,1:nz_tot,npcon),avg_pcon(ld,ny,1:nz_tot,npcon))
   allocate(dissip0(ld,ny,1:nz),dissip(ld,ny,1:nz_tot),avg_dissip(ld,ny,1:nz_tot))
   allocate(breakage_freq0(ld,ny,1:nz,npcon),breakage_freq(ld,ny,1:nz_tot,npcon),&
           avg_breakage_freq(ld,ny,1:nz_tot,npcon))
   print *,"hello"
!   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz))
 !  allocate(u(ld,ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot),avg_u(ld,ny,1:nz_tot),avg_v(ld,ny,1:nz_tot))
!   allocate(avg_w(ld,ny,1:nz_tot))
!  avg_w=0._rprec
!  avg_u=0._rprec
!  avg_v=0._rprec
!  avg_pcon=0._rprec
  avg_dissip=0._rprec
  avg_breakage_freq =0._rprec
!  do jt=jt_start,jt_end
   jt=200
      jt_total=jt*base
    do ip=1,nproc
      write(file,'(A,A,I5.5,A,I2.2,A)')TRIM(dir),'vel_pcon_00',jt_total,'_20',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
!      read(2000+ip) u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon-1)
       read(2000+ip) junk,junk,junk,junk,breakage_freq0(:,:,1:nz,1:npcon-1),junk,dissip0(:,:,1:nz)
      nzs=(ip-1)*(nz-1)
!      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
!      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
!      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)
       dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)      
!      Pcon(:,:,nzs+1:nzs+nz,1:npcon-1) = Pcon0(:,:,1:nz,1:npcon-1)
       breakage_freq(:,:,nzs+1:nzs+nz,1:npcon-1) =breakage_freq0(:,:,1:nz,1:npcon-1)
      close(2000+ip)
      write(*,*) maxval(abs(dissip)),ip,nzs+nz
   enddo
  print *,"Opened files,timestep=",jt_total
!   avg_u = avg_u +u/nt    
!    avg_v = avg_v +v/nt
!     avg_w = avg_w +w/nt
!    avg_pcon = avg_pcon +pcon/nt
 avg_dissip=avg_dissip+(dissip**2)/nt
  avg_breakage_freq = avg_breakage_freq + breakage_freq/nt    


!  write(file,'(A)') 'max_vel_proc.txt'
!  open(unit = 30, FILE=file,status="unknown",position="append",action="write")
!  write(30,*)jt_total,maxval(-w*u_star)
 ! close(30)
!  enddo
!  print *, "w= ", maxval(-w*u_star)
avg_dissip = sqrt(avg_dissip)
505 format(18e15.7)
!deallocate(Pcon0,Pcon)
deallocate(dissip0,breakage_freq0,dissip,breakage_freq)
!open(unit=20,FILE='av_vel.out',form='unformatted')
!open(unit=21,FILE='P1_15_10000.dat')
!open(unit=25,FILE='rise_vel_source.dat')
!write(25,*) 'variables=x,z,w'
!write(25,*) 'zone t="',1,'"i=',nx,'k=',nz_tot,'f=point'
!  do jz=1,nz_tot
!    do jx=1,nx
!    write(25,505)(jx-1)*L_x/nx, -(jz-1)*z_i/nz_tot,-avg_w(jx,yps,jz)*u_star
!    enddo
!  enddo
!close(25)


!write(20) avg_u(:,:,1:nz_tot),avg_v(:,:,1:nz_tot),avg_w(:,:,1:nz_tot)
!close(20)
!open(unit=20,FILE='vel_10.dat')
!!open(unit=21,FILE='Pcon_10.dat')
!deallocate(u,v,w,u0,v0,w0)
!write(20,*) 'variables=x,y,z,avg_u,avg_v,avg_w'
!write(20,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!      do jz=1,nz_tot
!          do jy=1,ny
!            do jx=1,nx
!               write(20,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,avg_u(jx,jy,jz),&
!                 avg_v(jx,jy,jz), avg_w(jx,jy,jz)
!            enddo
!          enddo
!       enddo
!close(20)
!deallocate(avg_u,avg_v,avg_w)
!write(21,*) 'variables=x,y,z,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15'
!write(21,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!      do jz=1,nz_tot
!          do jy=1,ny
!            do jx=1,nx
!             write(21,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,(avg_pcon(jx,jy,jz,ipcon)*z_i**3,ipcon=1,npcon-1)
!           enddo
!          enddo
!        enddo
!deallocate(avg_pcon)
!!
!close(21)
!path='output_av/'
!write(file,'(A)')'av_dissiproot_br_freq.out'
!open(unit=20,form='unformatted')
!write(20)avg_dissip(:,:,1:nz_tot)*u_star**3/z_i,avg_breakage_freq(:,:,1:nz_tot,:)*u_star/z_i
!close(20)

deallocate(avg_dissip,avg_breakage_freq)

END PROGRAM POSTPROCESSING










