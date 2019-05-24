subroutine readfiles()

use param
use types,only:rprec
use sim_param,only : u,v,w,theta,PCon
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX,F_XX,F_KX2,F_XX2
use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,psi_m,deposition,Real_dep,Kc_t

implicit none

INTEGER   :: itrun,ntrun=1,ip
character (64) :: fname
INTEGER , PARAMETER :: nproc=48

!!!!!!! DATA VARIABLES!!!!!!!!
real (rprec), dimension (ld, ny, $lbz:nz) :: u0, v0, w0,theta0
real (rprec), dimesnion (ld,ny,$lbz:nz,npcon) :: Pcon0
real (rprec), dimension(ld, ny, 1:nz_tot, npcon) :: PCon_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot



   do itrun=1,ntrun
 
      !-----------------
      ! read in data
      !-----------------
 
      idfile = ntstart-ntrun*ntd + (irun-1)*ntrun*ntd  + ntd*itrun
      write(cfile,'(a,i7.7,a)') 'vel_sc',idfile
 
      do ip=1,nproc
         if(ntstart+ntd*ntrun*(irun-1).lt.1e5) then
            write(fname,'(a,i7.7,a,i3.3,a)')
cpath(irun)//'vel_sc',idfile,'_',ip-1,'.out'
         else if(ntstart+ntd*ntrun*(irun-1).lt.1e6) then
            write(fname,'(a,i7.7,a,i3.3,a)')
cpath2(irun)//'vel_sc',idfile,'_',ip-1,'.out'
         endif
         open(unit=2000+ip,file=fname,form='unformatted')
 
         print*, 'Open file successfully!',ip-1
 
         read (2000+ip) u0(:,:,1:nz0), v0(:,:,1:nz0), w0(:,:,1:nz0), &
              theta0(:,:,1:nz0), PCon0(:,:,1:nz0,1:npcon), &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
 
         nzs=(ip-1)*(nz0-1)
         u(:,:,nzs+1:nzs+nz0) = u0(:,:,1:nz0)
         v(:,:,nzs+1:nzs+nz0) = v0(:,:,1:nz0)
         w(:,:,nzs+1:nzs+nz0) = w0(:,:,1:nz0)
         theta(:,:,nzs+1:nzs+nz0) = theta0(:,:,1:nz0)
         PCon(:,:,nzs+1:nzs+nz0,1:npcon) =
max(PCon0(:,:,1:nz0,1:npcon),1.e-24_rprec)
 
         close(2000+ip)
 
      enddo
 
      print*, 'File '//cfile//' readin completely!'
 
      !-------------
      ! end here
      !-------------
 
   enddo
