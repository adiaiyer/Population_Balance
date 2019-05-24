module types
implicit none

public
integer, parameter :: rprec = kind (1.d0)
end module types

subroutine calculate_integral(Pcon,i_start,i_end,j_start,j_end,k_start,k_end,dsd_pcon)

integer, intent(in) :: j_start,j_end,k_start,k_end,i_start,i_end
integer             :: ipcon,npcon=15
real,dimension(:,:,:,:),allocatable :: Pcon,Pcon_vol
real,dimension(:), allocatable,intent(out) :: dsd_pcon

allocate(PCon(ld,ny,1:nz_tot,npcon),Pcon_vol(i_end-i_start+1,j_end-j_start+1,k_end-k_start+1,npcon))

Pcon_vol = Pcon(i_start:i_end,j_start:j_end,k_start:k_end,:)

allocate(dsd_pcon(npcon))
do ipcon=1,npcon
  dsd_pcon(ipcon) = sum(Pcon_vol(:,:,:,ipcon))
enddo

end subroutine calculate_integral


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
  INTEGER,PARAMETER                          :: jt_start = 375
  INTEGER,PARAMETER                          :: jt_end = 650
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = z_i/nproc
  REAL(rprec),PARAMETER                      :: L_x=2
  REAL(rprec),PARAMETER                      :: L_y = 0.78125
  REAL(rprec),PARAMETER                      :: u_star = 0.00006125
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0002
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,avg_dissip
  REAL(rprec),dimension(:,:,:),allocatable   ::u,v,w,avg_u,avg_v,avg_w,dissip
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon0,breakage_freq0,Re0,avg_rhs1,avg_rhs2
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,avg_Pcon,breakage_freq,tot_vol,rhs1,rhs2,&
                                               rhs10,rhs20,Re
  CHARACTER(len=64)                           :: file,path,filename
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon
  logical                                     :: exist
  real(rprec)                                        :: t_vol,dy=l_y/ny,dz=z_i/nz_tot
  INTEGER,PARAMETER                           :: xps=nx/8,yps=ny/2
  real(rprec),dimension(:),allocatable::diameter,tot,d50,weight_x,weight_y,weight_z
  real(rprec),dimension(:,:),allocatable               ::M_tot,M_y,M_z,z_cm,y_cm
  
 real(rprec),parameter::pi=3.1415926535897932384626433_rprec


  dir ='output_data/'



   allocate(Pcon0(ld,ny,1:nz,npcon),Pcon(ld,ny,1:nz_tot,npcon),avg_pcon(ld,ny,1:nz_tot,npcon),breakage_freq(ld,ny,1:nz_tot,npcon))
    print *,"hello"
   allocate(u0(ld,ny,1:nz),v0(ld,ny,1:nz),w0(ld,ny,1:nz),breakage_freq0(ld,ny,1:nz,npcon),dissip0(ld,ny,1:nz),&
           Re0(ld,ny,1:nz,npcon),Rhs10(ld,ny,1:nz,npcon),Rhs20(ld,ny,1:nz,npcon),&
  u(ld,ny,1:nz_tot),v(ld,ny,1:nz_tot),w(ld,ny,1:nz_tot),avg_u(ld,ny,1:nz_tot),avg_v(ld,ny,1:nz_tot))
    allocate(dissip(ld,ny,1:nz_tot),Rhs1(ld,ny,1:nz_tot,npcon),Rhs2(ld,ny,1:nz_tot,npcon))
   allocate(avg_w(ld,ny,1:nz_tot),avg_dissip(ld,ny,1:nz_tot),avg_rhs1(ld,ny,1:nz_tot,npcon),avg_rhs2(ld,ny,1:nz_tot,npcon),&
            tot_vol(ld,ny,nz_tot,npcon))
   allocate(tot(npcon-1),d50(nx),diameter(npcon))


allocate(weight_x(nx),weight_y(ny),weight_z(nz_tot))

  weight_y = dy/2
  weight_z = dz/2
  weight_y(2:ny-1) = dy
  weight_z(2:nz_tot-1) = dz

  avg_w=0._rprec
  avg_u=0._rprec
  avg_v=0._rprec
  avg_pcon=0._rprec
  avg_dissip = 0._rprec
  avg_rhs1 = 0._rprec
  avg_rhs2 = 0._rprec
  do jt=jt_start,jt_end
!   jt=399
      jt_total=jt*base
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
!      breakage_freq(:,:,nzs+1:nzs+nz,1:npcon) = breakage_freq0(:,:,1:nz,1:npcon)
!!      Re(:,:,nzs+1:nzs+nz,1:npcon) = Re0(:,:,1:nz,1:npcon)
!      Rhs1(:,:,nzs+1:nzs+nz,1:npcon )= Rhs10(:,:,1:nz,1:npcon)
!      Rhs2(:,:,nzs+1:nzs+nz,1:npcon) = Rhs20(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
  print *,"Opened files,timestep=",jt_total
!  print *,pcon(295,yps,168,1),rhs1(295,yps,168,1)


!  write(*,*) maxval(Pcon(:,:,:,15))
   avg_u = avg_u +u/nt    
   avg_v = avg_v +v/nt
   avg_w = avg_w +w/nt
   avg_pcon = avg_pcon +pcon/nt
   avg_dissip = avg_dissip + dissip/nt
!   avg_rhs1 = avg_rhs1+ rhs1/nt
!   avg_rhs2 = avg_rhs2 + rhs2/nt
enddo

path = './post_proc_output/av_375_650/'
write(filename,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH),'av_dissip_',jt_start, '_',jt_end,'.out'
open(unit=12,Form='unformatted',File = filename)
write(12) avg_dissip(1:nx,1:ny,1:nz_tot)
close(12)

write(filename,'(A,A,I3.3,A,I3.3,A)')TRIM(PATH), 'av_u_',jt_start, '_',jt_end,'.out'
open(unit=12,Form='unformatted',File = filename)
write(12) avg_u(1:nx,1:ny,1:nz_tot)
close(12)

write(filename,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH),'av_v_',jt_start, '_',jt_end,'.out'
open(unit=12,Form='unformatted',File = filename)
write(12) avg_v(1:nx,1:ny,1:nz_tot)
close(12)

write(filename,'(A,A,I3.3,A,I3.3,A)')TRIM(PATH), 'av_w_',jt_start, '_',jt_end,'.out'
open(unit=12,Form='unformatted',File = filename)
write(12) avg_w(1:nx,1:ny,1:nz_tot)
close(12)

write(filename,'(A,A,I3.3,A,I3.3,A)')TRIM(PATH), 'av_pcon_',jt_start, '_',jt_end,'.out'
open(unit=12,Form='unformatted',File = filename)
write(12) avg_pcon(1:nx,1:ny,1:nz_tot,:)
close(12)





!!! Now that you calculated avg pcon, calculate d50 diameters right here.
!   open(21,File='diameter.dat')
!   do ip=1,npcon
!     read(21,*) diameter(ip)
!   enddo
!  close(21)
!
!  do ip=1,npcon
!    tot_vol(:,:,:,ip) = pi/6_rprec*avg_pcon(:,:,:,ip)*diameter(ip)**3_rprec
!  enddo
!allocate(M_tot(310-xps+1,npcon),M_y(310-xps+1,npcon),M_z(310-xps+1,npcon),y_cm(310-xps+1,npcon),z_cm(310-xps+1,npcon))
!   M_tot=0._rprec
!   M_y = 0._rprec
!   M_z = 0._rprec
!   y_cm= 0._rprec
!   z_cm= 0._rprec
!    do jx = xps,310
!    do jy= 1,ny
!    do jz = 1,nz_tot
!        M_tot(jx-xps+1,:) = M_tot(jx-xps+1,:) +weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
!        M_y(jx-xps+1,:) =  M_y(jx-xps+1,:) +(jy-1)*L_y/ny*weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
!        M_z(jx-xps+1,:) =  M_z(jx-xps+1,:) +(jz-1)*z_i/nz_tot*weight_y(jy)*weight_z(jz)*tot_vol(jx,jy,jz,:)
!
!    enddo
!    enddo
!
!    enddo
!    y_cm = M_y/M_tot - (yps-1)*L_y/ny
!    z_cm = M_z/M_tot
!
!path ='./'
!!
!505 format(18e15.7)
!write(file,'(A)')  'y_cm.dat'
!open(unit=22,FILE=file,status='unknown')
!
!write(file,'(A)')  'z_cm.dat'
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
!  write(23,505) (jx-1)*l_x/nx,(384._rprec/nz_tot-z_cm(jx-xps+1,ip),ip=1,npcon)
!enddo
!
!close(22)
!
!
!close(23)
!
!    jy=yps
!    jz =168
!
!
!  write(file,'(A,I2.2,A)')'d50_av_z_', int((384-jz-1)*100/385),'.dat'
!  open(34,FILE=file,status='unknown',position='append')
!  write(34,*) 'x,d50'
!  do jx = 48,5*nx/6
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        print *,ip, tot(ip)
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(jx) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!!          write(*,*)i,ip
!          exit
!        endif
!      enddo
!    write(34,'(E15.7E2,E15.7E2)') (jx-1)*l_x/nx,d50(jx)
!  enddo
!
!  close(34)
!  jz =188
!  write(file,'(A,I2.2,A)')'d50_av_z_', int((384-jz-1)*100/385),'.dat'
!  open(34,FILE=file,status='unknown',position='append')
!  write(34,*) 'x,d50'
!
!  do jx = 48,5*nx/6
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        print *,ip, tot(ip)
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(jx) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!!          write(*,*)i,ip
!          exit
!        endif
!      enddo
!    write(34,'(E15.7E2,E15.7E2)') (jx-1)*l_x/nx,d50(jx)
!  enddo
!  close(34)
!
!  jz =208
!  write(file,'(A,I2.2,A)')'d50_av_z_', int((384-jz-1)*100/385),'.dat'
!  open(34,FILE=file,status='unknown',position='append')
!  write(34,*) 'x,d50'
!  do jx = 48,5*nx/6
!    t_vol = sum(tot_vol(jx,jy,jz,1:npcon-1))
!    tot = 0._rprec
!    tot(1)= tot_vol(jx,yps,jz,1)/t_vol
!      do ip=2,npcon-1
!        tot(ip) = tot(ip-1) +tot_vol(jx,jy,jz,ip)/t_vol
!        print *,ip, tot(ip)
!        if(tot(ip) .ge. 0.5_rprec)  then
!          d50(jx) = diameter(ip-1) +(diameter(ip)-diameter(ip-1))*(0.5-tot(ip-1))&
!                                               /(tot(ip)-tot(ip-1))
!!          write(*,*)i,ip
!          exit
!        endif
!      enddo
!    write(34,'(E15.7E2,E15.7E2)') (jx-1)*l_x/nx,d50(jx)
!
!  enddo
!close(34)
!where(avg_pcon .LE. 1e-24_rprec)
!   avg_pcon = 1e-24_rprec
!endwhere
!
!
! avg_pcon = log10(avg_pcon)
!
!
!jx =200
!
!  write(file,'(A,A)')'av_conc','.dat'
!  open(35,FILE=file,status='unknown',position='append')
!  write(35,*) 'z,p1,p7,p12,p15'
!do jz=1,nz_tot
! write(35,505) (385._rprec-jz)/nz_tot,avg_pcon(jx,yps,jz,1),avg_pcon(jx,yps,jz,7),avg_pcon(jx,yps,jz,12),avg_pcon(jx,yps,jz,15)
!enddo
!close(35)
!open(unit=20,FILE="Pcon_30000.out",form="unformatted")
!write(20) PCon(1:nx,1:ny,1:nz_tot,1:npcon)
!close(20)
!
deallocate(Pcon0)
deallocate(breakage_freq0)!breakage_freq)
deallocate(u,v,w,u0,v0,w0,Re0,dissip0)!dissip)
!deallocate(rhs10,rhs20,rhs1,rhs2)
write(file,'(A,I4.4,A,I4.4,A)')'av_rhs_pcon_',jt_start,'_',jt_end,'.out'
open(unit=20,FILE=file,form='unformatted')
!write(20) avg_u(1:nx,1:ny,1:nz_tot),avg_v(1:nx,1:ny,1:nz_tot),avg_w(1:nx,1:ny,1:nz_tot)
write(20)avg_PCon(1:nx,1:ny,1:nz_tot,1:npcon)!,avg_rhs1(1:nx,1:ny,1:nz_tot,1:npcon),avg_rhs2(1:nx,1:ny,1:nz_tot,1:npcon),&
!          avg_dissip(1:nx,1:ny,1:nz_tot)
close(20)
!
!
!!open(unit=21,FILE='Pcon_8000.dat')
!!write(22,*) 'variables=x,y,z,avg_u,avg_v,avg_w,dissip'
!!write(22,*) 'zone t="',1,'"i=',nx,'j=',ny,'k=',nz_tot,'f=point'
!!      do jz=1,nz_tot
!!          do jy=1,ny
!!            do jx=1,nx
!!               write(22,505)(jx-1)*L_x/nx,(jy-1)*L_y/ny,-(jz-1)*z_i/nz_tot,avg_u(jx,jy,jz)*u_star,&
!!                 avg_v(jx,jy,jz)*u_star,-avg_w(jx,jy,jz)*u_star,avg_dissip(jx,jy,jz)*(u_star**3)/z_i
!!            enddo
!!          enddo
!!       enddo
!!close(22)
deallocate(avg_w,avg_dissip)
!
!write(21,*) 'variables=x,z,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15'
!write(21,*) 'zone t="',1,'"i=',nx,'k=',nz_tot,'f=point'
! jy=yps
!      do jz=1,nz_tot
!    !      do jy=1,ny
!            do jx=1,nx
!             write(21,505)(jx-1)*L_x/nx,-(jz-1)*z_i/nz_tot,(pcon(jx,jy,jz,ipcon)*z_i**3,ipcon=1,npcon)
!           enddo
!   !      enddo
!        enddo
deallocate(avg_pcon,pcon,avg_rhs1,avg_rhs2)

END PROGRAM POSTPROCESSING


