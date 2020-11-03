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
  INTEGER,PARAMETER                           :: ny=288
  INTEGER,PARAMETER                           :: nz_tot=385
  INTEGER,PARAMETER                           :: nz = (385-1)/nproc +1
  INTEGER,PARAMETER                           :: nx=288
  INTEGER,PARAMETER                           :: nsteps=200
  INTEGER,PARAMETER                           :: base=350
  INTEGER,PARAMETER                           :: npcon=20
  INTEGER,PARAMETER                          :: lh=nx/2+1,ld=2*lh
  INTEGER,PARAMETER                          :: jt_start = 330
  INTEGER,PARAMETER                          :: jt_end = 485
  INTEGER,PARAMETER                           :: nt=jt_end-jt_start+1
  REAL(rprec),PARAMETER                      :: z_i=1._rprec
  REAL(rprec),PARAMETER                      :: L_z = 2.5_rprec/nproc
  REAL(rprec),PARAMETER                      :: L_x=1
  REAL(rprec),PARAMETER                      :: L_y = 1
  REAL(rprec),PARAMETER                      :: u_star = 1
  REAL(rprec),PARAMETER                      :: dt_dim = 0.0001
  REAL(rprec),PARAMETER                      :: dt = dt_dim*u_star/z_i
  REAL(rprec),PARAMETER                      :: T_scale = 300.
  REAL(rprec),dimension(:,:,:),allocatable   :: u0,v0,w0,dissip0,avg_dissip
  REAL(rprec),dimension(:,:,:),allocatable   ::u,v,w,avg_u,avg_v,avg_w,dissip
  REAL (rprec), dimension(:,:,:),allocatable :: u2,v2,w2,avg_u2,avg_v2,avg_w2
  REAL (rprec), dimension(:,:,:),allocatable :: uw2,vw2,uv,avg_uw2,avg_vw2,avg_uv
  REAL (rprec), dimension(:,:,:),allocatable :: uw,vw,avg_uw,avg_vw,C,w_uv,avg_w_uv
  REAL (rprec), dimension(:,:,:),allocatable :: d32,tot_ar,tot_vol,tvol2,avg_tvol2
  REAL (rprec), dimension(:,:,:),allocatable ::mu_d32,mu_tot_ar,sigma_d32,sigma_tot_ar
  REAL (rprec), dimension(:,:,:),allocatable ::mu_d32old,mu_tot_arold
  REAL (rprec), dimension(:,:,:),allocatable ::mu_epsold,sigma_eps
  REAL(rprec),dimension(:,:,:,:),allocatable ::ts,mu_ts,mu_tsold,sigma_ts
  REAL(rprec),dimension(:,:,:,:),allocatable ::uc,vc,avg_uc,avg_vc,wc,avg_wc
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon0,breakage_freq0,Re0,avg_rhs1,avg_rhs2
  REAL(rprec),dimension(:,:,:,:),allocatable ::Pcon,avg_Pcon,breakage_freq,rhs1,rhs2,&
                                               rhs10,rhs20,Re,pcon2,avg_pcon2
  CHARACTER(len=64)                           :: file,path
  CHARACTER(len=64)                           :: dir
  INTEGER                                     ::ip,jt,jx,jy,jz,jt_total,nzs,ipcon,reclen,counter
  logical                                     :: file_exists,initialized=.false.
  real(rprec)                                        ::  t_vol,dy=l_y/ny,dz=L_z/z_i/(nz-1)
  INTEGER,PARAMETER                           :: xps=nx/2,yps=ny/2
  real(rprec),dimension(:),allocatable::diameter,tot,d50,weight_x,weight_y,weight_z
  real(rprec),dimension(:,:),allocatable               ::M_tot,M_y,M_z,z_cm,y_cm
 integer,dimension(:,:,:),allocatable :: pdf_counter
integer,dimension(:,:,:,:),allocatable :: ts_counter

 real(rprec),parameter::pi=3.1415926535897932384626433_rprec
character (128) :: fname_vel,fname_vel2,fname_pcon,fname_pcon2,fname

  dir ='output_data/'



!allocate(pdf_counter(nx,ny,1:nz_tot))
!allocate(ts_counter(nx,ny,1:nz_tot,npcon))


!allocate(Re0(ld,ny,1:nz,npcon))
!allocate(breakage_freq0(ld,ny,1:nz,npcon))

!!! Pcon allocate
allocate(Pcon0(ld,ny,1:nz,npcon))
allocate(Pcon(ld,ny,1:nz_tot,npcon))
allocate(avg_pcon(nx,ny,1:nz_tot,npcon))
!allocate(Pcon2(nx,ny,1:nz_tot,npcon))
!allocate(avg_pcon2(nx,ny,1:nz_tot,npcon))   

!allocate(rhs10(ld,ny,1:nz,npcon))
!allocate(rhs1(ld,ny,1:nz_tot,npcon))
!allocate(avg_rhs1(nx,ny,1:nz_tot,npcon))
!! Pcon fluxes
!allocate(avg_uc(nx,ny,1:nz_tot,npcon))
!allocate(avg_vc(nx,ny,1:nz_tot,npcon))
!allocate(avg_wc(nx,ny,1:nz_tot,npcon))
!allocate(uc(nx,ny,1:nz_tot,npcon))
!allocate(vc(nx,ny,1:nz_tot,npcon))
!allocate(wc(nx,ny,1:nz_tot,npcon))

!! Velocity field

!file read
allocate(u0(ld,ny,1:nz))
allocate(v0(ld,ny,1:nz))
allocate(w0(ld,ny,1:nz))
allocate(dissip0(ld,ny,1:nz))

! total vel
allocate(u(ld,ny,1:nz_tot))
allocate(v(ld,ny,1:nz_tot))
allocate(w(ld,ny,1:nz_tot))
!allocate(uw(ld,ny,1:nz_tot))
!allocate(vw(ld,ny,1:nz_tot))
!allocate(w_uv(ld,ny,1:nz_tot))
allocate(dissip(ld,ny,1:nz_tot))
!!avg_vel
allocate(avg_u(nx,ny,1:nz_tot))
allocate(avg_v(nx,ny,1:nz_tot))
allocate(avg_w(nx,ny,1:nz_tot))
!allocate(avg_uw(nx,ny,1:nz_tot))
!allocate(avg_vw(nx,ny,1:nz_tot))
!allocate(avg_w_uv(nx,ny,1:nz_tot))
allocate(avg_dissip(nx,ny,1:nz_tot))

!!vel 2

!allocate(u2(nx,ny,1:nz_tot))
!allocate(v2(nx,ny,1:nz_tot))
!allocate(w2(nx,ny,1:nz_tot))
!allocate(uw2(nx,ny,1:nz_tot))
!allocate(vw2(nx,ny,1:nz_tot))
!allocate(uv(nx,ny,1:nz_tot))

!! avg vel 2
!allocate(avg_u2(nx,ny,1:nz_tot))
!allocate(avg_v2(nx,ny,1:nz_tot))
!allocate(avg_w2(nx,ny,1:nz_tot))
!allocate(avg_uw2(nx,ny,1:nz_tot))
!allocate(avg_vw2(nx,ny,1:nz_tot))
!allocate(avg_uv(nx,ny,1:nz_tot))


allocate(diameter(npcon))

!allocate(d32(nx,ny,1:nz_tot))
!allocate(mu_d32(nx,ny,1:nz_tot))
!allocate(mu_d32old(nx,ny,1:nz_tot))
!allocate(sigma_d32(nx,ny,1:nz_tot))
!
!
!allocate(mu_epsold(nx,ny,1:nz_tot))
!allocate(sigma_eps(nx,ny,1:nz_tot))
!
!allocate(ts(nx,ny,1:nz_tot,npcon))
!allocate(mu_ts(nx,ny,1:nz_tot,npcon))
!allocate(mu_tsold(nx,ny,1:nz_tot,npcon))
!allocate(sigma_ts(nx,ny,1:nz_tot,npcon))
!
!
!allocate(tot_ar(nx,ny,1:nz_tot))
!allocate(mu_tot_ar(nx,ny,1:nz_tot))
!allocate(mu_tot_arold(nx,ny,1:nz_tot))
!allocate(sigma_tot_ar(nx,ny,1:nz_tot))

!allocate(tot_vol(nx,ny,1:nz_tot))
!allocate(tvol2(nx,ny,1:nz_tot))
!allocate(avg_tvol2(nx,ny,1:nz_tot))

!allocate(C(nx,ny,1:nz_tot))



path = './post_proc_output/'

inquire(iolength=reclen) u(1:nx,1:ny,1:nz_tot)
write(fname,'(A,A)')TRIM(PATH), "avg_vel_470_500.out"
inquire(FILE = fname,EXIST=file_exists)

avg_pcon= 0._rprec
!avg_pcon2 = 0._rprec
!avg_rhs1 = 0._rprec
!tvol2 = 0._rprec
!avg_tvol2 = 0._rprec
!avg_uc= 0._rprec
!avg_vc= 0._rprec
!avg_wc= 0._rprec



avg_u= 0._rprec
avg_v= 0._rprec
avg_w= 0._rprec
!avg_uw= 0._rprec
!avg_vw= 0._rprec
!avg_w_uv = 0._rprec


!avg_u2= 0._rprec
!avg_v2= 0._rprec
!avg_w2= 0._rprec
!avg_uw2= 0._rprec
!avg_vw2= 0._rprec
!avg_uv= 0._rprec


!mu_d32 = 0._rprec
!mu_tot_ar = 0._rprec
!sigma_d32 = 0._rprec
!sigma_tot_ar = 0._rprec
!
!sigma_eps = 0._rprec
!
!mu_ts = 0._rprec
!sigma_ts = 0._rprec
!C =0._rprec
!tot_ar = 0._rprec
!tot_vol = 0._rprec
!ts = 0._rprec
!ts_counter = 1
!pdf_counter = 1
If (initialized) then

        if(file_exists) then
                open(10,file=fname,access='direct',RECL=reclen,form='unformatted')
                read(10,rec=1)avg_u(1:nx,1:ny,1:nz_tot)
                read(10,rec=2)avg_v(1:nx,1:ny,1:nz_tot)
                read(10,rec=3)avg_w(1:nx,1:ny,1:nz_tot)
                close(10)
        else
                avg_w=0._rprec
                avg_u=0._rprec
                avg_v=0._rprec
        endif
        
        
        write(fname,'(A,A)')TRIM(PATH), "avg_velrms2_470_500.out"
        inquire(FILE = fname,EXIST=file_exists)
        
        if(file_exists) then
                open(10,file=fname,access='direct',RECL=reclen,form='unformatted')
                read(10,rec=1)avg_u2(1:nx,1:ny,1:nz_tot)
                read(10,rec=2)avg_v2(1:nx,1:ny,1:nz_tot)
                read(10,rec=3)avg_w2(1:nx,1:ny,1:nz_tot)
                close(10)
        else
                avg_w2=0._rprec
                avg_u2=0._rprec
                avg_v2=0._rprec
        endif
        
        
        write(fname,'(A,A)')TRIM(PATH), "avg_pcon_470_500.out"
        inquire(FILE = fname,EXIST=file_exists)
        
        if(file_exists) then
                open(10,file=fname,access='direct',RECL=reclen,form='unformatted')
                do ip = 1,npcon
                        read(10,rec=ip) avg_pcon(1:nx,1:ny,1:nz_tot,ip)
                enddo
                close(10)
        else
                avg_pcon=0._rprec
        endif
        

        write(fname,'(A,A)')TRIM(PATH), "avg_pcon_2470_500.out"
        inquire(FILE = fname,EXIST=file_exists)
        
        if(file_exists) then
                open(10,file=fname,access='direct',RECL=reclen,form='unformatted')
                do ip = 1,npcon
                        read(10,rec=ip) avg_pcon2(1:nx,1:ny,1:nz_tot,ip)
                enddo
                close(10)
        else
                avg_pcon2 = 0._rprec
        endif
endif


   open(21,File='diameter.dat')
   do ip=1,npcon
     read(21,*) diameter(ip)
   enddo
  close(21)
  avg_dissip = 0._rprec
!  avg_rhs1 = 0._rprec
!  avg_rhs2 = 0._rprec
counter = 0
  do jt=jt_start,jt_end
!  jt=375
      jt_total=jt*350
    do ip=1,nproc
      write(file,'(A,A,I6.6,A,I3.3,A)')TRIM(dir),'vel_pcon_0',jt_total,'_2',ip-1,'.out'
      open(unit=2000+ip,FILE=file,form='unformatted')
      read(2000+ip) u0(:,:,1:nz),v0(:,:,1:nz),w0(:,:,1:nz),Pcon0(:,:,1:nz,1:npcon),dissip0(:,:,1:nz), &
                  breakage_freq0(:,:,1:nz,1:npcon),Re0(:,:,1:nz,1:npcon),rhs10(:,:,1:nz,1:npcon)
      nzs=(ip-1)*(nz-1)
      u(:,:,nzs+1:nzs+nz) = u0(:,:,1:nz)
      v(:,:,nzs+1:nzs+nz) = v0(:,:,1:nz)
      w(:,:,nzs+1:nzs+nz) = w0(:,:,1:nz)      
      dissip(:,:,nzs+1:nzs+nz) = dissip0(:,:,1:nz)
      Pcon(:,:,nzs+1:nzs+nz,1:npcon) = Pcon0(:,:,1:nz,1:npcon)
!      breakage_freq(:,:,nzs+1:nzs+nz,1:npcon) = breakage_freq0(:,:,1:nz,1:npcon)
!      Re(:,:,nzs+1:nzs+nz,1:npcon) = Re0(:,:,1:nz,1:npcon)
!      Rhs1(:,:,nzs+1:nzs+nz,1:npcon )= Rhs10(:,:,1:nz,1:npcon)
!      Rhs2(:,:,nzs+1:nzs+nz,1:npcon) = Rhs20(:,:,1:nz,1:npcon)
      close(2000+ip)
   enddo
  print *,"Opened files,timestep=",jt_total
!  print *,pcon(295,yps,168,1),rhs1(295,yps,168,1)


   counter = counter +1
!   uw(:,:,2:nz_tot) = (u(:,:,2:nz_tot) + u(:,:,1:nz_tot-1))/2 
!   vw(:,:,2:nz_tot) = (v(:,:,2:nz_tot) + v(:,:,1:nz_tot-1))/2
!   w_uv(:,:,1:nz_tot-1) = (w(:,:,1:nz_tot-1) + w(:,:,2:nz_tot))*0.5_rprec
!   u2 = u(1:nx,:,:)**2 ! on uv-grid
!   v2 = v(1:nx,:,:)**2 ! on uv-grid
!   w2 = w(1:nx,:,:)**2 ! on w-grid
!
!   uw2 = uw(1:nx,:,:)*w(1:nx,:,:)! on w-grid
!   vw2 = vw(1:nx,:,:)*w(1:nx,:,:)! on w-grid
!   
!   uv = u(1:nx,:,:)*v(1:nx,:,:) ! on uv-grid

!tot_ar = 0._rprec
!tot_vol = 0._rprec
!   do ip =1,npcon
!           uc(:,:,:,ip) = u(1:nx,:,:)*pcon(1:nx,:,:,ip) * pi/6 * diameter(ip)**3
!           vc(:,:,:,ip) = v(1:nx,:,:)*pcon(1:nx,:,:,ip) * pi/6 * diameter(ip)**3
!           wc(:,:,:,ip) = w_uv(1:nx,:,:)*pcon(1:nx,:,:,ip)  * pi/6 * diameter(ip)**3
!           pcon2(:,:,:,ip) = (pcon(1:nx,:,:,ip)  * pi/6 * diameter(ip)**3)**2
!           tot_ar = tot_ar + pcon(1:nx,:,:,ip) * pi * diameter(ip)**2
!           tot_vol = tot_vol + pcon(1:nx,:,:,ip)  * pi/6 * diameter(ip)**3
!   enddo
!  where (tot_vol .ge. 1.413e-8_rprec ) 
!        d32 = 6 * tot_vol/tot_ar
!        mu_d32old = mu_d32
!        mu_d32  = mu_d32  + (d32 - mu_d32)/pdf_counter
!        sigma_d32 = sigma_d32 + (d32 - mu_d32old)*(d32-mu_d32)
!        pdf_counter = pdf_counter + 1 
!endwhere
!         do ip = 1,npcon
!         where(pcon(1:nx,:,:,ip) .ge. maxval(pcon(1:nx,:,:,ip))*1.e-4_rprec)
!                ts(:,:,:,ip) = rhs1(1:nx,:,:,ip)/pcon(1:nx,:,:,ip)
!                mu_tsold(:,:,:,ip)= mu_ts(:,:,:,ip)
!                mu_ts(:,:,:,ip)  = mu_ts(:,:,:,ip)  + (ts(:,:,:,ip)- mu_ts(:,:,:,ip))/ts_counter(:,:,:,ip)
!                sigma_ts(:,:,:,ip) = sigma_ts(:,:,:,ip) + (ts(:,:,:,ip) - mu_tsold(:,:,:,ip))&
!                *(ts(:,:,:,ip)-mu_ts(:,:,:,ip))
!                ts_counter(:,:,:,ip) = ts_counter(:,:,:,ip)+1
!        endwhere
!        enddo
!  mu_epsold = avg_dissip
!  avg_dissip  = avg_dissip  + (dissip(1:nx,:,:) - avg_dissip)/counter
!  sigma_eps = sigma_eps + (dissip(1:nx,:,:) - mu_epsold)*(dissip(1:nx,:,:)-avg_dissip)
!   tvol2 = tot_vol**2
!   mu_tot_arold = mu_tot_ar
!   mu_tot_ar = mu_tot_ar + (tot_ar - mu_tot_ar)/counter
!   sigma_tot_ar = sigma_tot_ar + (tot_ar - mu_tot_arold)*(tot_ar-mu_tot_ar)


!  write(*,*) maxval(Pcon(:,:,:,15))

   avg_u = avg_u +u(1:nx,:,:)/nt    
   avg_v = avg_v +v(1:nx,:,:)/nt
   avg_w = avg_w +w(1:nx,:,:)/nt
!
!   avg_uw = avg_uw + uw(1:nx,:,:)/nt  ! w grid
!   avg_vw = avg_vw + vw(1:nx,:,:)/nt  ! w grid
!   avg_w_uv = avg_w_uv + w_uv(1:nx,:,:)/nt ! uv grid
   avg_pcon = avg_pcon +pcon(1:nx,:,:,:)/nt

!   avg_uw2 = avg_uw2 + uw2/nt ! w grid
!   avg_vw2 = avg_vw2 + vw2/nt ! w grid
!   avg_uv  = avg_uv  +  uv/nt
!
!   avg_u2 = avg_u2 +u2/nt    
!   avg_v2 = avg_v2 +v2/nt
!   avg_w2 = avg_w2 +w2/nt
!   avg_uc = avg_uc + uc/nt ! uv grid
!   avg_vc = avg_vc + vc/nt ! uv grid
!   avg_wc = avg_wc + wc/nt ! uv grid
!   avg_pcon2 = avg_pcon2 + pcon2/nt
!   avg_tvol2 = avg_tvol2 + tvol2/nt
!   avg_rhs1 = avg_rhs1+ rhs1/nt
!   avg_rhs2 = avg_rhs2 + rhs2/nt

print * ,'yolo'
enddo

!sigma_d32 = sigma_d32/(pdf_counter)
!sigma_tot_ar = sigma_tot_ar/(counter-1)
!sigma_eps = sigma_eps/(counter-1)
!sigma_ts = sigma_ts/(ts_counter)

!do ip = 1,npcon
!
!           C = C + avg_pcon(1:nx,:,:,ip)  * pi/6 * diameter(ip)**3
!enddo


!!open(unit=10,FILE='restart_av.out',form='unformatted')
open(unit=10,FILE='restart_vel_pcon.out',form='unformatted')
write(10)avg_u,avg_v,avg_w,avg_pcon
close(10)

inquire(iolength=reclen) u(1:nx,1:ny,1:nz_tot)
write(fname_vel,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_vel_',jt_start,'_',jt_end,'.out'
open(unit=10,FILE=fname_vel,form='unformatted',access='direct',recl=reclen)
write(10,rec=1) avg_u(1:nx,1:ny,1:nz_tot)
write(10,rec=2) avg_v(1:nx,1:ny,1:nz_tot)
write(10,rec=3) avg_w(1:nx,1:ny,1:nz_tot)
close(10)
!
!
!write(fname_vel2,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_velrms_',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_vel2,form='unformatted',access='direct',recl=reclen)
!write(10,rec=1) avg_u2(1:nx,1:ny,1:nz_tot) - (avg_u(1:nx,1:ny,1:nz_tot))**2
!write(10,rec=2) avg_v2(1:nx,1:ny,1:nz_tot) - (avg_v(1:nx,1:ny,1:nz_tot))**2
!write(10,rec=3) avg_w2(1:nx,1:ny,1:nz_tot) - (avg_w(1:nx,1:ny,1:nz_tot))**2
!close(10)
!
!write(fname_vel2,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_rs_cross',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_vel2,form='unformatted',access='direct',recl=reclen)
!write(10,rec=1) avg_uw2(1:nx,1:ny,1:nz_tot) - avg_uw * avg_w
!write(10,rec=2) avg_vw2(1:nx,1:ny,1:nz_tot) - avg_vw * avg_w
!write(10,rec=3) avg_uv(1:nx,1:ny,1:nz_tot)  - avg_u  * avg_v
!close(10)

write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_pcon_',jt_start,'_',jt_end,'.out'
open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
do ip = 1,npcon
        write(10,rec=ip) avg_pcon(1:nx,1:ny,1:nz_tot,ip)
enddo
close(10)

!write(fname_pcon2,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_pcon_rms',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon2,form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(10,rec=ip) avg_pcon2(1:nx,1:ny,1:nz_tot,ip)- (avg_pcon(1:nx,1:ny,1:nz_tot,ip) * pi/6 * diameter(ip)**3)**2
!enddo
!close(10)

!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_pconx_flux_',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(10,rec=ip) avg_uc(1:nx,1:ny,1:nz_tot,ip) - avg_u * avg_pcon(:,:,:,ip) * pi/6 * diameter(ip) **3
!enddo
!close(10)
!
!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_pcony_flux_',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(10,rec=ip) avg_vc(1:nx,1:ny,1:nz_tot,ip) -avg_v * avg_pcon(:,:,:,ip) * pi/6 * diameter(ip) **3  
!enddo
!close(10)
!
!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_C_rms',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!        write(10,rec=1) avg_tvol2 - C**2
!close(10)


!write(fname_vel,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'pdf_stats_2',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_vel,form='unformatted',access='direct',recl=reclen)
!write(10,rec=1) mu_d32
!write(10,rec=2) sigma_d32
!write(10,rec=3) mu_tot_ar
!write(10,rec=4) sigma_tot_ar
!close(10)
!
!write(fname_vel,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'pdf_eps_',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_vel,form='unformatted',access='direct',recl=reclen)
!write(10,rec=1) avg_dissip
!write(10,rec=2) sigma_eps
!close(10)
!
!
!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'pdf_muts_2',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(10,rec=ip) mu_ts(1:nx,1:ny,1:nz_tot,ip)
!enddo
!close(10)
!
!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'pdf_sigmats_2',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!do ip = 1,npcon
!        write(10,rec=ip) sigma_ts(1:nx,1:ny,1:nz_tot,ip)
!enddo
!close(10)
!write(fname_pcon,'(A,A,I3.3,A,I3.3,A)') TRIM(PATH), 'avg_C_flux_tot_',jt_start,'_',jt_end,'.out'
!open(unit=10,FILE=fname_pcon,form='unformatted',access='direct',recl=reclen)
!write(10,rec=1) sum(avg_uc(1:nx,1:ny,1:nz_tot,:),4) - avg_u * C
!write(10,rec=2) sum(avg_vc(1:nx,1:ny,1:nz_tot,:),4) - avg_v * C
!write(10,rec=3) sum(avg_wc(1:nx,1:ny,1:nz_tot,:),4) - avg_w_uv * C
!close(10)


deallocate(Pcon0)
deallocate(Pcon)
deallocate(avg_pcon)
!deallocate(Pcon2)
!deallocate(avg_pcon2)   
!deallocate(C)
!deallocate(avg_uc)
!deallocate(avg_vc)
!deallocate(avg_wc)
!deallocate(uc)
!deallocate(vc)
!deallocate(wc)
deallocate(u0)
deallocate(v0)
deallocate(w0)
deallocate(dissip0)
deallocate(u)
deallocate(v)
deallocate(w)
!deallocate(uw)
!deallocate(vw)
!deallocate(avg_u)
!deallocate(avg_v)
!deallocate(avg_w)
!deallocate(avg_uw)
!deallocate(avg_vw)
!deallocate(u2)
!deallocate(v2)
!deallocate(w2)
!deallocate(uw2)
!deallocate(vw2)
!deallocate(uv)
!deallocate(avg_u2)
!deallocate(avg_v2)
!deallocate(avg_w2)
!deallocate(avg_uw2)
!deallocate(avg_vw2)
!deallocate(avg_uv)
!deallocate(diameter)
!deallocate(d32)
!deallocate(mu_d32)
!deallocate(mu_d32old)
!deallocate(sigma_d32)
!deallocate(tot_ar)
!deallocate(mu_tot_ar)
!deallocate(mu_tot_arold)
!deallocate(sigma_tot_ar)
!deallocate(tot_vol)

END PROGRAM POSTPROCESSING


