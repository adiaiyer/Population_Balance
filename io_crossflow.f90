module io
use types,only:rprec
use param, only : ld, nx, ny, nz, jx_s, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc,      &
                  average_dim_num, nz_tot,jt_total,p_count,dt,z_i,u_star
implicit none
private
public :: openfiles,output_loop,output_final,                   &
          mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,         &
          inflow_read, inflow_write, avg_stats,calculate_avgs,   &
          create_folder_for_output
public :: average_dim_select_flag, dim1_size, dim2_size,        &
          dim1_global, dim2_global, collocate_MPI_averages,     &
          compute_avg_var,post_spec,log_file_create

!integer, parameter :: nz_tot = (nz - 1) * nproc + 1

integer,parameter::num_hour_out=1
!integer,parameter::base=nint(num_hour_out*3600/(dt*z_i/u_star)),nwrite=base
integer,parameter::base=100,nwrite=base
!integer,parameter::base=100,nwrite=base
integer,parameter:: avg_base=1
!!!!  io_spec=.true. output plan-averaged spectrum
logical,parameter:: io_spec=.true.,output_fields_3d_flag=.false.
!integer,parameter::spec_write_freqz=p_count/3, fields_3d_write_freqz=p_count*6
!integer,parameter::spec_write_freqz=base/num_hour_out, fields_3d_write_freqz=p_count*6
integer,parameter::spec_write_freqz=3000, fields_3d_write_freqz=p_count*6
!integer,parameter::spec_write_start=1,spec_write_end=nint(86400._rprec/(dt*z_i/u_star))
!integer,parameter::spec_write_start=base*6,spec_write_end=base*12
integer,parameter::spec_write_start=1,spec_write_end=1000000
!integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! slice_inst sets the value of the y node for the chosen x-z inst slice
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimensio:n of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process 
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical,parameter:: inst_slice_flag=.false.
logical,parameter:: suppress_out =.true. !AA condition to suppres output
integer,parameter:: num_vars=4 ! works currently only for u,v,w,T due to the size difference in Cs
integer,parameter:: slice_inst=(nz_tot-1)/2, inst_slice_freqz=5, inst_array_lim=200
!real(kind=rprec),dimension(ld,nz+1,num_vars*inst_array_lim)::inst_slice_array
!integer:: inst_slice_counter


logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

integer, parameter :: n_avg_stats = 100
                      !--interval for updates in avg_stats
character (*), parameter :: end_hdr_avg = '# end header'


!! --------------------------------------------------------------------
!! The following block defines parameters for use in avgslice and scalar_slice
!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

!!!!  time_spec>0 output time series spectrum (need additional calcu.)
integer,parameter::time_spec=0
integer::n_obs, jt_total_init
integer,allocatable::obs_pt(:,:)

!!!!  io_mean=.true. output small domain time-averaged velocity
logical,parameter::io_mean=.false.
integer,parameter::jx_pls=1,jx_ple=nx,width=ny/2-1
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width+1
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz)::&
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

!!!!  io_lambda2
logical,parameter::io_lambda2=.false.
real(kind=rprec),dimension(nx,ny,nz)::lam2
CHARACTER(len=255)::cwd,makedirectory,folder,conc,vel,Re_data,Rhs_data,output1,&
                      brk_freq,all_files,all_files2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine openfiles()
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

integer::i

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    read (1, *) jt_total
    jt_total_init=jt_total 
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0
    jt_total_init=jt_total 
  end if

end if

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
  open(13,file=path//'output/check_ke.out',status="unknown",position="append")
end if

if(time_spec.gt.0)then
  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
  if(jt_total.eq.0)rewind(15)
endif

if(io_mean)then
  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
  if(jt_total.eq.0)then
    rewind(51)
    write(51)jx_pls,jx_ple,jy_pls,jy_ple
  endif
endif

fCS1plan = path // 'output/CS1plan.out'
fCS2plan = path // 'output/CS2plan.out'
fCS4plan = path // 'output/CS4plan.out'
fVISCplan = path // 'output/VISCplan.out'
fDISSplan = path // 'output/DISSplan.out'
fCS1Vplan = path // 'output/CS1Vplan.out'
fCS2Vplan = path // 'output/CS2Vplan.out'
fCS4Vplan = path // 'output/CS4Vplan.out'

$if ($MPI)
  !--append coordinate identifiers
  write (temp, '(".c",i0)') coord
  fCS1plan = trim (fCS1plan) // temp
  fCS2plan = trim (fCS2plan) // temp
  fCS4plan = trim (fCS4plan) // temp
  fVISCplan = trim (fVISCplan) // temp
  fDISSplan = trim (fDISSplan) // temp
  fCS1Vplan = trim (fCS1Vplan) // temp
  fCS2Vplan = trim (fCS2Vplan) // temp
  fCS4Vplan = trim (fCS4Vplan) // temp
$endif

!open (90, file=fCS1plan, form='unformatted')
!open (91, file=fCS2plan, form='unformatted')
!open (92, file=fCS4plan, form='unformatted')
!open (93, file=fVISCplan, form='unformatted')
!open (94, file=fDISSplan, form='unformatted')
!open (95, file=fCS1Vplan, form='unformatted')
!open (96, file=fCS2Vplan, form='unformatted')
!open (97, file=fCS4Vplan, form='unformatted')

!open(98,file="avg_Pconjet.dat",status='unknown',position='append')

if(time_spec.gt.0)then
open(1,file=path//'obs.pt')
read(1,*)n_obs
allocate(obs_pt(1:2,n_obs))
do i=1,n_obs
read(1,*)obs_pt(1:2,i)
enddo
close(1)
endif

end subroutine openfiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop
use param,only:output,dt,c_count,S_FLAG,SCAL_init,jt,jan_diurnal_run,PCon_FLAG,PCon_init,Pcon_count
use sim_param,only:path,u,v,w,dudz,dudx,p,&
     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:Cs_opt2
use scalars_module,only:sgs_t3
use scalars_module2,only:scalar_slice,budget_TKE_scalar,pollen_slice
implicit none
!integer,intent(in)::jt
!real(kind=rprec),dimension(ld,ny,nz)::usp,vsp
real(kind=rprec),dimension(nz)::u_ndim
character(len=20)::req

character (64) :: fname, temp

integer::jx,jy,jz
integer:: fields_3d_write_freqz_temp

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension(ld,$lbz:nz,num_vars*inst_array_lim),save::inst_slice_array
integer,save:: inst_slice_counter

jt_total=jt_total+1
call calculate_mean

!print *,'Entering output loop'
!!$!+++++++++++++++++++++++++++++++++++++
!!$!DY Added by Di Yang
!!$!DY Output oil inflow condition for ENDLESS
!!$call output_plane_oil_flux
!!$!+++++++++++++++++++++++++++++++++++++

if ((use_avgslice) .and. (mod(jt,c_count)==0)) then
         call avgslice()
       if ((S_FLAG) .and. (jt.GE.SCAL_init)) then
         call scalar_slice() ! Uses file unit numbers (36-47)
       endif
       if ((PCon_FLAG) .and. (jt.GE.PCon_init)) then
         call pollen_slice()
       endif
end if

!AA Added By Aditya Aiyer
!Call subroutine to compute average over z
!if ((mod(jt_total,1)==0)) then
!  print *, 'Entering average'
!  call calculate_avgs
!endif
!AA End here


if (output) then
  if (mod(jt_total,base)==0) then
!    print *,'writing to output files'
    if (S_FLAG .OR. PCon_FLAG) then
!*****************************
!DY Modified by Di Yang
!DY Allow larger jt_total
!DY Start here
!       print*, "io.f90 checkpoint1"
!DY    write (fname, '(a,i6.6,a)') path // 'output/vel_sc', jt_total, '.out'
       write (fname, '(a,i7.7,a,i3.3,a)') path // 'output/vel_sc', jt_total, '_', coord, '.out'
    else
!DY    write (fname, '(a,i6.6,a)') path // 'output/vel', jt_total, '.out'
       write (fname, '(a,i7.7,a,i3.3,a)') path // 'output/vel', jt_total, '_', coord, '.out'
    end if

!AA Changed checkpoint_v2 to checkpoint
!    call checkpoint_v2 (1000+coord)
        !call checkpoint (1000+coord)
!AA End of Change

!    print*, "io.f90 checkpoint4"
 !   close(1000+coord)

!DY End here
!****************************
 call output_vel_conc_binary
    !+++++++++++++++++++++++++++++++++++++
    !DY Added by Di Yang
    !DY Output oil inflow condition for ENDLESS
!    call output_plane_oil_flux
    !+++++++++++++++++++++++++++++++++++++
  end if
end if


  if ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then 
!  if ((io_spec) .and. mod(jt_total,spec_write_freqz)==0) then 
    if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
  end if

  if (time_spec.gt.0) call timeseries_spec

end subroutine output_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_avg_var(avg_var,raw_var,output_dim)
use param,only:nx,ny,nz,c_count,p_count
integer :: output_dim
!real(kind=rprec),dimension(:,:)::avg_var
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var
!real(kind=rprec),dimension(output_dim*(nx-1)+1,nz-1)::avg_var
real(kind=rprec),dimension(1:nx,1:ny,1:nz-1):: raw_var
real(kind=rprec):: fr
     
fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
select case (output_dim)
  case(1) !average over y and t
      avg_var=avg_var+fr*sum(raw_var(1:nx,1:ny,1:nz-1),dim=2)/ny
  case(2) ! average over x,y and t
       PRINT*,'fix io.f90 line 385'
!      avg_var(:,1)=avg_var(:,1)+fr*sum(sum(raw_var(1:nx,1:ny,1:nz-1),dim=1),dim=2)/(nx*ny)
  
end select
end subroutine compute_avg_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
!integer, parameter:: lbound_loc=average_dim_num*(nx-nz+1)+nz-1
!integer, parameter:: ubound_loc=average_dim_num*(nz-2)+1
!integer, parameter:: lbound_global=average_dim_num*(nx-nz_tot+1)+nz_tot-1
!integer, parameter:: ubound_global=average_dim_num*(nz_tot-2)+1

real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!!real(kind=rprec),dimension(nx,nz-1)::avg_var_proc
!!real(kind=rprec),dimension(nx,nz_tot-1)::avg_var_tot_domain

!local_filename=path//'output/aver_'//filename_str//'.out'

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  

      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

 5168     format(1400(E14.5))
end subroutine collocate_MPI_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
!integer, parameter:: lbound_loc=average_dim_num*(nx-nz+1)+nz-1
!integer, parameter:: ubound_loc=average_dim_num*(nz-2)+1
!integer, parameter:: lbound_global=average_dim_num*(nx-nz_tot+1)+nz_tot-1
!integer, parameter:: ubound_global=average_dim_num*(nz_tot-2)+1

character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!!real(kind=rprec),dimension(nx,nz-1)::avg_var_proc
!!real(kind=rprec),dimension(nx,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  

  open(file_ind,file=trim(local_filename),status="unknown",position="append")

      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
           do ind1=1,nx
             if (ABS(avg_var_tot_domain(ind1,ind2)) .LT. TINYS) then
               avg_var_tot_domain(ind1,ind2) = 0._rprec
             endif
           end do
         end do
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         do ind1=1,nz_tot-1
           if (ABS(avg_var_tot_domain(ind1,1)) .LT. TINYS) then
             avg_var_tot_domain(ind1,1) = 0._rprec
           endif
         end do
         write(file_ind,5168) jt_total*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

 5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine avgslice
use sim_param,only:path,u,v,w,dudz,dvdz, txx, txz, tyy, tyz, tzz, p,dudt,dvdt,dwdt
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
!use scalars_module2,only:collocate_MPI_averages
implicit none
!!integer, parameter :: nz_tot = (nz - 1) * nproc + 1
integer::i,j,k
real(kind=rprec),dimension(nx,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs
real(kind=rprec),dimension(nx,nz-1),save::adudz,advdz,aCs_Ssim,abeta_sgs,abetaclip_sgs
real(kind=rprec),dimension(nx,nz-1),save::atxx,atxz,atyy,atyz,atzz
real(kind=rprec),dimension(nx,nz-1),save::u3,v3,w3
real(kind=rprec),dimension(nx,nz-1),save::adudt,advdt,adwdt
real(kind=rprec)::tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,&
     tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr, arg1, arg2,tu3,tv3,tw3,   &
     tdudt,tdvdt,tdwdt     
real(kind=rprec)::tCs_Ssim
real(kind=rprec),dimension(:,:),allocatable::avg_out


!print *,'u from io',u(1:2,1:2,1:2)
fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
do k=1,Nz-1
do i=1,Nx
   tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
   ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
   ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
   tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
   tCs=0._rprec;tCs_Ssim=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec
   tdudt=0._rprec;tdvdt=0._rprec;tdwdt=0._rprec

   do j=1,Ny
      tu1=tu1+u(i,j,k)
      tv1=tv1+v(i,j,k)
      tw1=tw1+w(i,j,k)
      tp1=tp1+p(i,j,k)
      ttxx=ttxx+txx(i,j,k)
      ttxz=ttxz+txz(i,j,k)
      ttyy=ttyy+tyy(i,j,k)
      ttyz=ttyz+tyz(i,j,k)
      ttzz=ttzz+tzz(i,j,k)
      tdudz=tdudz+dudz(i,j,k)
      tdvdz=tdvdz+dvdz(i,j,k)
      tu2=tu2+u(i,j,k)*u(i,j,k)
      tv2=tv2+v(i,j,k)*v(i,j,k)
      tw2=tw2+w(i,j,k)*w(i,j,k)
      tp2=tp2+p(i,j,k)*p(i,j,k)
      tCs=tCs+sqrt(Cs_opt2(i,j,k))
      tCs_Ssim=tCs_Ssim+sqrt(Cs_Ssim(i,j,k))
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec
         arg2=0._rprec
      else  
         arg1=(u(i,j,k)+u(i,j,k-1))/2.
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
      end if
      tuw=tuw+w(i,j,k)*arg1
      tvw=tvw+w(i,j,k)*arg2
      tu3=tu3+u(i,j,k)*u(i,j,k)*u(i,j,k)
      tv3=tv3+v(i,j,k)*v(i,j,k)*v(i,j,k)
      tw3=tw3+w(i,j,k)*w(i,j,k)*w(i,j,k)
      tdudt=tdudt+dudt(i,j,k)
      tdvdt=tdvdt+dvdt(i,j,k)
      tdwdt=tdwdt+dwdt(i,j,k)
   end do
   au(i,k)=au(i,k)+(fr)*tu1/Ny
   av(i,k)=av(i,k)+(fr)*tv1/Ny
   aw(i,k)=aw(i,k)+(fr)*tw1/Ny
   ap(i,k)=ap(i,k)+(fr)*tp1/Ny
   adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
   advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
   u2(i,k)=u2(i,k)+(fr)*tu2/Ny
   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
   p2(i,k)=p2(i,k)+fr*tp2/Ny
   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
   aCs_Ssim(i,k)=aCs_Ssim(i,k)+fr*tCs_Ssim/Ny
   abeta_sgs(i,k)=abeta_sgs(i,k)+fr*Beta_avg(k)
   abetaclip_sgs(i,k)=abetaclip_sgs(i,k)+fr*Betaclip_avg(k)
   u3(i,k)=u3(i,k)+(fr)*tu3/Ny
   v3(i,k)=v3(i,k)+(fr)*tv3/Ny
   w3(i,k)=w3(i,k)+(fr)*tw3/Ny
   adudt(i,k)=adudt(i,k)+(fr)*tdudt/Ny
   advdt(i,k)=advdt(i,k)+(fr)*tdvdt/Ny
   adwdt(i,k)=adwdt(i,k)+(fr)*tdwdt/Ny
end do
end do

!if(mod(jt,p_count)==0.and.jt.gt.2000) then
if (mod(jt,p_count)==0) then

  allocate(avg_out(1:nx,1:(nz_tot-1)));
call collocate_MPI_averages_N(au,avg_out,20,'u')
call collocate_MPI_averages_N(av,avg_out,21,'v')
call collocate_MPI_averages_N(aw,avg_out,22,'w')
call collocate_MPI_averages_N(ap,avg_out,23,'p')
call collocate_MPI_averages_N(u2,avg_out,24,'u2')
call collocate_MPI_averages_N(v2,avg_out,25,'v2')
call collocate_MPI_averages_N(w2,avg_out,26,'w2')
call collocate_MPI_averages_N(p2,avg_out,32,'p2')
call collocate_MPI_averages_N(atxx,avg_out,27,'txx')
call collocate_MPI_averages_N(atxz,avg_out,28,'txz')
call collocate_MPI_averages_N(atyy,avg_out,29,'tyy')
call collocate_MPI_averages_N(atyz,avg_out,30,'tyz')
call collocate_MPI_averages_N(atzz,avg_out,31,'tzz')
call collocate_MPI_averages_N(auw,avg_out,33,'uw')
call collocate_MPI_averages_N(avw,avg_out,34,'vw')
call collocate_MPI_averages_N(aCs,avg_out,35,'Cs')
call collocate_MPI_averages_N(adudz,avg_out,36,'dudz')
call collocate_MPI_averages_N(advdz,avg_out,37,'dvdz')
call collocate_MPI_averages_N(aCs_Ssim,avg_out,38,'Cs_Ssim')
call collocate_MPI_averages_N(abeta_sgs,avg_out,39,'beta_sgs')
call collocate_MPI_averages_N(abetaclip_sgs,avg_out,40,'betaclip_sgs');
call collocate_MPI_averages_N(u3,avg_out,41,'u3')
call collocate_MPI_averages_N(v3,avg_out,42,'v3')
call collocate_MPI_averages_N(w3,avg_out,43,'w3')
call collocate_MPI_averages_N(adudt,avg_out,44,'dudt')
call collocate_MPI_averages_N(advdt,avg_out,45,'dvdt')
call collocate_MPI_averages_N(adwdt,avg_out,46,'dwdt');deallocate(avg_out)

!VK Zero out the outputted averages !!
   au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
   w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
   atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
   adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
   abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
   adudt=0._rprec;advdt=0._rprec;adwdt=0._rprec;
end if
end subroutine avgslice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint_final (lun)
!use param, only : ld,ny,nz,nz_tot,S_FLAG,PCon_FLAG
use param
use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta, PCon
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX,F_XX,F_KX2,F_XX2
use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,psi_m,deposition,Real_dep,Kc_t
implicit none

integer, intent (in) :: lun

integer:: ipcon

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif

real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSx_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSy_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: RHSz_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot
real (rprec), dimension(ld, ny, 1:nz_tot, npcon) :: PCon_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: Cs_opt2_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: RHS_T_tot
real(kind=rprec),dimension(ld,ny,nz_tot,npcon) :: RHS_PCon_tot
real(kind=rprec),dimension(ld,ny,nz_tot) :: Kc_t_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_KX_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_XX_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_KX2_tot
real(kind=rprec),dimension(nx,ny,nz_tot) :: F_XX2_tot

!---------------------------------------------------------------------
$if ($MPI)
  sendcounts = size (u(:,:,1:nz))
  recvcounts = size (u(:,:,1:nz-1))
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (u(1,1,1), sendcounts, MPI_RPREC,&
                    u_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (v(1,1,1), sendcounts, MPI_RPREC,&
                    v_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (w(1,1,1), sendcounts, MPI_RPREC,&
                    w_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (RHSx(1,1,1), sendcounts, MPI_RPREC,&
                    RHSx_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (RHSy(1,1,1), sendcounts, MPI_RPREC,&
                    RHSy_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (RHSz(1,1,1), sendcounts, MPI_RPREC,&
                    RHSz_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (Cs_opt2(1,1,1), sendcounts, MPI_RPREC,&
                    Cs_opt2_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_LM(1,1,1), sendcounts, MPI_RPREC,&
                    F_LM_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_MM(1,1,1), sendcounts, MPI_RPREC,&
                    F_MM_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_QN(1,1,1), sendcounts, MPI_RPREC,&
                    F_QN_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_NN(1,1,1), sendcounts, MPI_RPREC,&
                    F_NN_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  if (S_FLAG) then
  call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
                    theta_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (RHS_T(1,1,1), sendcounts, MPI_RPREC,&
                    RHS_T_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif
  if (PCon_FLAG) then
     do ipcon=1,npcon
        call mpi_gatherv (PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
             PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
        call mpi_gatherv (RHS_PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
             RHS_PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  call mpi_gatherv (Kc_t(1,1,1), sendcounts, MPI_RPREC,&
                    Kc_t_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  sendcounts = size (F_KX(:,:,1:nz))
  recvcounts = size (F_KX(:,:,1:nz-1))
  call mpi_gatherv (F_KX(1,1,1), sendcounts, MPI_RPREC,&
                    F_KX_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_XX(1,1,1), sendcounts, MPI_RPREC,&
                    F_XX_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_KX2(1,1,1), sendcounts, MPI_RPREC,&
                    F_KX2_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (F_XX2(1,1,1), sendcounts, MPI_RPREC,&
                    F_XX2_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif
$endif

if((.not. USE_MPI)) then 

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
              PCon(:,:,1:nz,1:npcon), RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz), &
              RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, RHS_PCon(:,:,1:nz,1:npcon), &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN, &
              deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
elseif (S_FLAG) then !WITH SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
              RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m,Cs_opt2, F_LM, F_MM, &
              F_QN, F_NN
elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), PCon(:,:,1:nz,1:npcon),   &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
              RHS_PCon(:,:,1:nz,1:npcon), Cs_opt2, F_LM, F_MM,F_QN, F_NN,            &
              deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
else ! No SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
              RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
              Cs_opt2, F_LM, F_MM, F_QN, F_NN
end if

else if ((USE_MPI .and. coord == 0)) then

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot),   &
              PCon_tot(:,:,1:nz_tot,1:npcon), RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot), &
              RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m,RHS_PCon_tot(:,:,1:nz_tot,1:npcon), &
              Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot, &
              deposition(:,:),Real_dep,Kc_t_tot,F_KX_tot,F_XX_tot,F_KX2_tot,F_XX2_tot
elseif (S_FLAG) then !WITH SCALARS
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot),   &
              RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),          &
              RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:,1), psi_m,Cs_opt2_tot, F_LM_tot, F_MM_tot, &
              F_QN_tot, F_NN_tot
elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), PCon_tot(:,:,1:nz_tot,1:npcon),   &
              RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),          &
              RHS_PCon_tot(:,:,1:nz_tot,1:npcon), Cs_opt2_tot, F_LM_tot, F_MM_tot,F_QN_tot, F_NN_tot,            &
	      deposition(:,:),Real_dep,Kc_t_tot,F_KX_tot,F_XX_tot,F_KX2_tot,F_XX2_tot
else ! No SCALARS
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot),           &
              RHSx_tot(:, :, 1:nz_tot), RHSy_tot(:, :, 1:nz_tot), RHSz_tot(:, :, 1:nz_tot),  &
              Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
end if

endif

end subroutine checkpoint_final

!!AA Code output data in binary

subroutine output_vel_conc_binary

use param
use sim_param, only : u,v,w,PCon
use sgsmodule,only : Cs_opt2, F_LM,F_MM,F_QN,F_NN,F_KX,F_XX,F_KX2,F_XX2
use scalars_module, only : deposition,RHS_pbe,breakage_freq,fitl,Re
use bottombc, only : u_star_avg
use sgsmodule
implicit none

character (64) :: fname

if (S_FLAG .and. Pcon_FLAG) then
  write(fname,'(A,I7.7,A,I4.4,A)')path//'restart/vel_pcon_',jt_total,'_',2000+coord,'.out'
endif

open(unit=2000+coord,file=fname,form='unformatted')

if(S_FLAG .and. PCON_FLAG) then
  write (2000+coord)u(:,:,:1:nz),v(:,:,:,1:nz),w(:,:,1:nz),Pcon(:,:,1:nz,1:npcon-1),&
        breakage_freq(:,:,1:nz,1:npcon-1),Re(:,:,1:nz,1:npcon-1),dissipation(:,:,1:nz)
endif
close(2000+coord)

end subroutine output_vel_conc_binary

!!AA END here


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY The original version writes data into a single file, which requires large memory for high resolution
!DY This version writes data on each cpu into separate files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine checkpoint_final_v2
!use param, only : ld,ny,nz,nz_tot,S_FLAG,PCon_FLAG
use param
use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta, PCon
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX,F_XX,F_KX2,F_XX2
use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,psi_m,deposition,Real_dep,Kc_t
implicit none

integer:: ipcon

character (64) :: fname

if(S_FLAG.and.PCon_FLAG) then
   write(fname,'(A,I4.4,A)') path//'restart/vel_sc_pcon_',2000+coord,'.out'
else if(S_FLAG) then
   write(fname,'(A,I4.4,A)') path//'restart/vel_sc_',2000+coord,'.out'
else if(PCon_FLAG) then
   write(fname,'(A,I4.4,A)') path//'restart/vel_pcon_',2000+coord,'.out'
else
   write(fname,'(A,I4.4,A)') path//'restart/vel_',2000+coord,'.out'
endif

open(unit=2000+coord,file=fname,form='unformatted')

if (S_FLAG .AND. PCon_FLAG) then !WITH SCALARS AND POLLENS
   write (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
        PCon(:,:,1:nz,1:npcon), RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz), &
        RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m,RHS_PCon(:,:,1:nz,1:npcon), &
        Cs_opt2, F_LM, F_MM, F_QN, F_NN, deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
elseif (S_FLAG) then !WITH SCALARS
   write (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
        RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
        RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m,Cs_opt2, F_LM, F_MM, F_QN, F_NN
elseif (PCon_FLAG) then ! Pollen
   write (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), PCon(:,:,1:nz,1:npcon),   &
        RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),          &
        RHS_PCon(:,:,1:nz,1:npcon), Cs_opt2, F_LM, F_MM,F_QN, F_NN,            &
        deposition(:,:),Real_dep,Kc_t,F_KX,F_XX,F_KX2,F_XX2
else ! No SCALARS
   write (2000+coord) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),           &
        RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
        Cs_opt2, F_LM, F_MM, F_QN, F_NN
end if

close(2000+coord)

end subroutine checkpoint_final_v2
!+++++++++++++++++++++++++++++++++++++++++++++++
!DY Subroutine checkpoint_final_v2 end here
!+++++++++++++++++++++++++++++++++++++++++++++++





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)
use param
use sim_param, only : u, v, w, theta, PCon
use scalars_module, only : deposition,Real_dep,P_surf_flux,P_surf_flux_dep
use bottombc, only: ustar_avg
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean LES with Stokes drift
use sgsmodule
!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

integer, intent (in) :: lun
!DY Added by Di Yang for testing
integer jx,jy,jz,ipcon
!DY End here
!AA DEfine Str

character (LEN=50) :: str

!AA End Here

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif

real (rprec), dimension(ld, ny, 1:nz_tot) :: u_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: v_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: w_tot
real (rprec), dimension(ld, ny, 1:nz_tot) :: theta_tot
real (rprec), dimension(ld, ny, 1:nz_tot, npcon) :: PCon_tot
!real (rprec), dimension(ld, ny, 1:nz_tot) :: Cs_opt2_tot
!real (rprec), dimension(ld, ny, 1:nz_tot) :: Nu_t_tot
!logical,parameter :: suppress_out
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot

!---------------------------------------------------------------------
$if ($MPI) !LV1
  sendcounts = size (u(:,:,1:nz))
  recvcounts = size (u(:,:,1:nz-1))
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (u(1,1,1), sendcounts, MPI_RPREC,&
                    u_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (v(1,1,1), sendcounts, MPI_RPREC,&
                    v_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (w(1,1,1), sendcounts, MPI_RPREC,&
                    w_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$  call mpi_gatherv (Cs_opt2(1,1,1), sendcounts, MPI_RPREC,&
!!$                    Cs_opt2_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$  call mpi_gatherv (Nu_t(1,1,1), sendcounts, MPI_RPREC,&
!!$                    Nu_t_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$
!!$  call mpi_gatherv (F_LM(1,1,1), sendcounts, MPI_RPREC,&
!!$                    F_LM_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$  call mpi_gatherv (F_MM(1,1,1), sendcounts, MPI_RPREC,&
!!$                    F_MM_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$  call mpi_gatherv (F_QN(1,1,1), sendcounts, MPI_RPREC,&
!!$                    F_QN_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!$  call mpi_gatherv (F_NN(1,1,1), sendcounts, MPI_RPREC,&
!!$                    F_NN_tot(1,1,1), sendcounts, displs,       &
!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)

  if (S_FLAG) then  !LV2
  call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
                    theta_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif  !LV2
  if (PCon_FLAG) then !LV2
     do ipcon=1,npcon
        call mpi_gatherv (PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
             PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  endif !LV2
$endif  !LV1

!open(5,file = "hello.dat")
!write(5,*) "Hey i am in the checkpoint"
!close(5)

!AA Changed .not.
if((.not. USE_MPI)) then  !LV1
!if((USE_MPI)) then
!AA Change Ends 
if (S_FLAG .AND. PCon_FLAG) then !LV2 ! With scalars and pollen - Ying 11/01/2010
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
              PCon(:,:,1:nz,1:npcon), deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
elseif (S_FLAG) then !WITH SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz)
elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), PCon(:,:,1:nz,1:npcon),   &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
else ! No SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz)
end if  !LV2

else if ((USE_MPI .and. coord == 0)) then

if (S_FLAG .AND. PCon_FLAG) then !LV2 ! With scalars and pollen - Ying 11/01/2010
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
              theta_tot(:,:,1:nz_tot), PCon_tot(:,:,1:nz_tot,1:npcon), &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!DY Added by Di Yang for testing
!  write(100,*) 'variables=x,y,z,PCon'
!!$  open(1000000+jt_total)
!!$!  write(1000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,PCon3,,u,v,w,theta,Cs_opt2,Nu_t'
!!$  write(1000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
!!$  write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
!!$  do jz=1,nz_tot
!!$     do jy=1,ny
!!$        do jx=1,nx
!!$           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!!$!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
!!$                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), &
!!$                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
!!$                theta_tot(jx,jy,jz),Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  enddo
505 format(12e15.7)
!!$  close(1000000+jt_total)
!!$  open(3000000+jt_total)
!!$!  write(1000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,PCon3,,u,v,w,theta,Cs_opt2,Nu_t'
!!$  write(3000000+jt_total,*) 'variables=x,y,z,LM,MM,QN,NN'
!!$  write(3000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
!!$  do jz=1,nz_tot
!!$     do jy=1,ny
!!$        do jx=1,nx
!!$           write(3000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!!$                F_LM_tot(jx,jy,jz),F_MM_tot(jx,jy,jz),F_QN_tot(jx,jy,jz),F_NN_tot(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  enddo
!!$  close(3000000+jt_total)

!AA Commenting rest of the file opens
  
!if(.not. OCEAN_FLAG) !Ocean flag is true so this won't open
  open(11000000+jt_total)
!!$  write(11000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(11000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta'
  write(11000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',1,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=ny/2,ny/2
        do jx=1,nx
           write(11000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2),&
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(11000000+jt_total)
!endif



!AA Change begins here
     write(str,'(A,I6.6,A)') 'Oilslice.',jt_total,'.out'
     open(3,FILE=TRIM(str))
     write(3,*) 'variables =x,y,z,PCon1,PCon2,PCon3,w'
     write(3,*) 'zone t="',jt_total,'"i=',nx,' j=',ny,' k=',nz_tot,'f =point'
     if(OCEAN_FLAG) then  !LV3
        do jz=1,nz_tot
           do jy=1,ny
              do jx=1,nx
                 write(3,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,&
                   !                PCon_tot(jx,jz,1),PCon_tot(jx,jz,2),PCon_tot(jx,jz,3),
                   !                &
                   PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), PCon_tot(jx,jy,jz,3),&
                   -w_tot(jx,jy,jz)
              enddo
           enddo
        enddo
     else
        do jz=1,nz_tot
           do jx=1,nx
              write(3,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot,&
                   !                PCon_tot(jx,jz,1),PCon_tot(jx,jz,2),PCon_tot(jx,jz,3),
                   !                &
                   PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                   w_tot(jx,jy,jz) 
           enddo
        enddo
     endif  !LV3
     close(3)
!AA Ends Here


 ! if ((.not. suppress_out)) !LVa
!AA Make sure this file does'nt open
!if(.not. OCEAN_FLAG)
  open(12000000+jt_total)
!!$  write(12000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(12000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,PCon3,u,v,w,theta'
  write(12000000+jt_total,*) 'zone t="',jt_total,'" i=',1,' j=',ny,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=1,ny
        do jx=nx/2,nx/2
           write(12000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), PCon_tot(jx,jy,jz,3),&
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(12000000+jt_total)
  open(13000000+jt_total)
!!$  write(13000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(13000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,PCon3,u,v,w,theta'
  write(13000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',1,' f=point'
  do jz=1,1
     do jy=1,ny
        do jx=1,nx
           write(13000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), PCon_tot(jx,jy,jz,3),&
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(13000000+jt_total)
!endif
!!$  open(31000000+jt_total)
!!$  write(31000000+jt_total,*) 'variables=x,y,z,LM,MM,QN,NN'
!!$  write(31000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',1,' k=',nz_tot,' f=point'
!!$  do jz=1,nz_tot
!!$     do jy=ny/2,ny/2
!!$        do jx=1,nx
!!$           write(31000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!!$                F_LM_tot(jx,jy,jz),F_MM_tot(jx,jy,jz),F_QN_tot(jx,jy,jz),F_NN_tot(jx,jy,jz)
!!$        enddo
!!$     enddo
!!$  enddo
!!$  close(31000000+jt_total)
!DY End here
elseif (S_FLAG) then !WITH SCALARS
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
              theta_tot(:,:,1:nz_tot)
!DY Added by Di Yang for testing
!  write(100,*) 'variables=x,y,z,PCon'
  open(1000000+jt_total)
!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta,Cs_opt2,Nu_t'
  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta'
  write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=1,ny
        do jx=1,nx
           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot,&
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz),theta_tot(jx,jy,jz)!, &
!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(1000000+jt_total)
!DY End here
elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
              PCon_tot(:,:,1:nz_tot,1:npcon),   &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!DY Added by Di Yang for testing
!  write(100,*) 'variables=x,y,z,PCon'
  open(1000000+jt_total)
!  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,Pcon2,Pcon3,Cs_opt2,Nu_t'
!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,Cs_opt2,Nu_t'
  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,PCon3'
  write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=1,ny
        do jx=1,nx
           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), PCon_tot(jx,jy,jz,3)!, &
!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(1000000+jt_total)
!DY End here
else ! No SCALARS
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot)

!endif !LVa

endif  !LV2

endif   !LV1

end subroutine checkpoint



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY The original version writes data into a single file, which requires large memory for high resolution
!DY This version writes data on each cpu into separate files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine checkpoint_v2 (lun)
!use param
!use sim_param, only : u, v, w, theta, PCon
!use scalars_module, only :deposition,Real_dep,P_surf_flux,P_surf_flux_dep,RHS_pbe,Re,breakage_freq,Pcon_neg
!use bottombc, only: ustar_avg
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!DY Added by Di Yang for ocean LES with Stokes drift
!use sgsmodule
!!DY End here
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!implicit none
!
!integer, intent (in) :: lun
!!DY Added by Di Yang for testing
!integer jx,jy,jz,ipcon,time_tot
!!DY End here
!
!$if ($MPI)
!  integer :: sendcounts(nproc),sendcounts_yz(nproc),sendcounts_neg(nproc)
!  integer :: recvcounts(nproc),recvcounts_yz(nproc),recvcounts_neg(nproc)
!  integer :: displs(nproc),displs_yz(nproc),displs_neg(nproc)
!$endif
!LOGICAL,PARAMETER :: YZ_PLANE=.true.
!character (LEN=255) :: str,str1,str2,str3,str4,str5,conc_xz,conc_yz,vel_xz,vel_yz
!real (rprec), dimension(ld,ny,1:nz_tot) :: Pcon_neg_tot
!real (rprec), dimension(ld, 1:nz) :: u1,v1,w1,theta1,dissip1
!real (rprec), dimension(ld, 1:nz, npcon) :: PCon1,RHS_pbe1,Re_1,breakage_freq1
!real (rprec), dimension(ld, 1:nz_tot) :: u_tot
!real (rprec), dimension(ld, 1:nz_tot) :: v_tot
!real (rprec), dimension(ld, 1:nz_tot) :: w_tot
!real (rprec), dimension(ld, 1:nz_tot) :: theta_tot
!real (rprec), dimension(ld, 1:nz_tot) :: dissip_tot
!real (rprec), dimension(ld, 1:nz_tot, npcon) :: PCon_tot,RHS_pbe_tot,Re_tot,breakage_freq_tot
!!real (rprec), dimension(ld, ny, 1:nz_tot) :: Cs_opt2_tot
!!real (rprec), dimension(ld, ny, 1:nz_tot) :: Nu_t_tot
!!AA
!real (rprec), dimension(ny,1:nz) :: u1_yz,v1_yz,w1_yz,theta1_yz,dissip1_yz
!real (rprec), dimension(ny,1:nz, npcon) :: PCon1_yz,RHS_pbe1_yz,Re_1_yz,breakage_freq1_yz 
!real (rprec), dimension(ny,1:nz_tot) :: u_tot_yz,v_tot_yz,w_tot_yz,theta_tot_yz,dissip_tot_yz
!real (rprec), dimension(ny,1:nz_tot,npcon)::Pcon_tot_yz,RHS_pbe_tot_yz,Re_tot_yz,breakage_freq_tot_yz
!logical, parameter :: check_neg=.false.
!!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
!!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
!!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
!!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot
!
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!DY Added by Di Yang for ocean LES with Stokes drift
!real(kind=rprec),dimension($lbz:nz) :: ust
!real (rprec) :: z
!!DY End here
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!DY Added by Di Yang for ocean LES with Stokes drift
!!DY Start here
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!if(OCEAN_FLAG .and. STOKES_FLAG) then !LV1
!   do jz=$lbz,nz
!      $if ($MPI)
!      z = (coord*(nz-1) + jz - 0.5_rprec) * dz
!      $else
!      z=(real(jz)-0.5_rprec)*dz
!      $endif
!      ust(jz)=U_stokes*exp(-2._rprec*wavenm_w*z) + ucross
!   enddo
!else
!   ust = 0._rprec
!endif !LV1
!!+++++++++++++++ 
!!DY End here
!!+++++++++++++++
!
!!---------------------------------------------------------------------
!$if ($MPI) !LV1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!Changed by Di Yang to include Stokes drift and uniform crossflow 
!!DY  u1 = u(:,yps,1:nz)
!do jz=1,nz
!   do jx=1,nx
!      u1(jx,jz) = u(jx,yps,jz) + ust(jz)
!   enddo
!enddo
!
!!AA add yz plane stuff
!if(yz_plane) then
!  do jz=1,nz
!    do jy=1,nx
!      u1_yz(jy,jz) = u(xps,jy,jz) +ust(jz)
!    enddo
!  enddo
!endif
!!AA end
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  v1 = v(:,yps,1:nz)
!  w1 = w(:,yps,1:nz)
!  theta1 = theta(:,yps,1:nz)
!  PCon1 = PCon(:,yps,1:nz,:)
!!AA Added RHS_PBE and Dissipation
!  RHS_pbe1 = RHS_pbe(:,yps,1:nz,:)
!  dissip1 = dissip(:,yps,1:nz)
!  Re_1 = Re(:,yps,1:nz,:)
!  breakage_freq1 = breakage_freq(:,yps,1:nz,:)
!!AA End here
!
!!AA Add variables for yz plane
!if(yz_plane) then
!  v1_yz = v(xps,:,1:nz)
!  w1_yz = w(xps,:,1:nz)
!  theta1_yz = theta(xps,:,1:nz)
!  PCon1_yz = PCon(xps,:,1:nz,:)
!!AA Added RHS_PBE and Dissipation
!  RHS_pbe1_yz = RHS_pbe(xps,:,1:nz,:)
!  dissip1_yz = dissip(xps,:,1:nz)
!  Re_1_yz = Re(xps,:,1:nz,:)
!  breakage_freq1_yz = breakage_freq(xps,:,1:nz,:)
!endif
!!AA End Here
!!write(*,*) "reached here"
!  sendcounts = size (u1(:,1:nz))
!  recvcounts = size (u1(:,1:nz-1))
!  displs = coord_of_rank * recvcounts
!!  write(*,*) sendcounts
!!  write(*,*) displs
!  call mpi_gatherv (u1(1,1), sendcounts, MPI_RPREC,&
!                    u_tot(1,1), sendcounts, displs, &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  call mpi_gatherv (v1(1,1), sendcounts, MPI_RPREC,&
!                    v_tot(1,1), sendcounts, displs,&
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  call mpi_gatherv (w1(1,1), sendcounts, MPI_RPREC,&
!                    w_tot(1,1), sendcounts, displs,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$  call mpi_gatherv (Cs_opt2(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    Cs_opt2_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$  call mpi_gatherv (Nu_t(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    Nu_t_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$
!!!$  call mpi_gatherv (F_LM(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    F_LM_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$  call mpi_gatherv (F_MM(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    F_MM_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$  call mpi_gatherv (F_QN(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    F_QN_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!!!$  call mpi_gatherv (F_NN(1,1,1), sendcounts, MPI_RPREC,&
!!!$                    F_NN_tot(1,1,1), sendcounts, displs,       &
!!!$                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!  if (S_FLAG) then  !LV2
!  call mpi_gatherv (theta1(1,1), sendcounts, MPI_RPREC,&
!                    theta_tot(1,1), sendcounts, displs,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  endif !LV2
!  if (PCon_FLAG) then !LV2
!     do ipcon=1,npcon
!        call mpi_gatherv (PCon1(1,1,ipcon), sendcounts, MPI_RPREC,&
!             PCon_tot(1,1,ipcon), sendcounts, displs,       &
!             MPI_RPREC, rank_of_coord(0), comm, ierr)
!     enddo
!  endif !LV2
!  if (PBE_FLAG) then
!
!    call mpi_gatherv (Pcon_neg(1,1,1), sendcounts_neg, MPI_RPREC,&
!                    Pcon_neg_tot(1,1,1), sendcounts_neg, displs_neg,&
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!    do ipcon=1,npcon
!      call mpi_gatherv (RHS_pbe1(1,1,ipcon), sendcounts, MPI_RPREC,&
!           RHS_pbe_tot(1,1,ipcon), sendcounts, displs,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!      
!       call mpi_gatherv (Re_1(1,1,ipcon), sendcounts, MPI_RPREC,&
!           Re_tot(1,1,ipcon), sendcounts, displs,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!      call mpi_gatherv (breakage_freq1(1,1,ipcon), sendcounts, MPI_RPREC,&
!           breakage_freq_tot(1,1,ipcon), sendcounts, displs,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!     enddo
!     call mpi_gatherv (dissip1(1,1), sendcounts, MPI_RPREC,&
!                    dissip_tot(1,1), sendcounts, displs,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!     
!  endif
!!write(*,*) "finished mpi for xz plane"
!!AA MPI stuff for the yz plane data
!if(yz_plane) then
!
!!  write(*,*) "inside yz plane loop"
!  sendcounts_yz = size (u1_yz(:,1:nz))
!  recvcounts_yz = size (u1_yz(:,1:nz-1))
!  displs_yz = coord_of_rank * recvcounts_yz
!  call mpi_gatherv (u1_yz(1,1), sendcounts_yz, MPI_RPREC,&
!                    u_tot_yz(1,1), sendcounts_yz, displs_yz,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  call mpi_gatherv (v1_yz(1,1), sendcounts_yz, MPI_RPREC,&
!                    v_tot_yz(1,1), sendcounts_yz, displs_yz,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  call mpi_gatherv (w1_yz(1,1), sendcounts_yz, MPI_RPREC,&
!                    w_tot_yz(1,1), sendcounts_yz, displs_yz,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!  
!
!if (S_FLAG) then  !LV2
!  call mpi_gatherv (theta1_yz(1,1), sendcounts_yz, MPI_RPREC,&
!                    theta_tot_yz(1,1), sendcounts_yz, displs_yz,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!  endif !LV2
!  if (PCon_FLAG) then !LV2
!     do ipcon=1,npcon
!        call mpi_gatherv (PCon1_yz(1,1,ipcon), sendcounts_yz, MPI_RPREC,&
!             PCon_tot_yz(1,1,ipcon), sendcounts_yz, displs_yz,       &
!             MPI_RPREC, rank_of_coord(0), comm, ierr)
!     enddo
!  endif !LV2
!  if (PBE_FLAG) then
!    do ipcon=1,npcon
!      call mpi_gatherv (RHS_pbe1_yz(1,1,ipcon), sendcounts_yz, MPI_RPREC,&
!           RHS_pbe_tot_yz(1,1,ipcon), sendcounts_yz, displs_yz,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!       call mpi_gatherv (Re_1_yz(1,1,ipcon), sendcounts_yz, MPI_RPREC,&
!           Re_tot_yz(1,1,ipcon), sendcounts_yz, displs_yz,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!      call mpi_gatherv (breakage_freq1_yz(1,1,ipcon), sendcounts_yz, MPI_RPREC,&
!           breakage_freq_tot_yz(1,1,ipcon), sendcounts_yz, displs_yz,             &
!           MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!     enddo
!     call mpi_gatherv (dissip1_yz(1,1), sendcounts_yz, MPI_RPREC,&
!                    dissip_tot_yz(1,1), sendcounts_yz, displs_yz,       &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!  endif
!endif
!!AA End here
!$endif !LV1
!
!if((.not. USE_MPI)) then !LV1
!
!if (S_FLAG .AND. PCon_FLAG) then ! LV2 With scalars and pollen - Ying 11/01/2010
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz),   &
!              PCon(:,:,1:nz,1:npcon), deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!elseif (S_FLAG) then !WITH SCALARS
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), theta(:,:,1:nz)
!elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), PCon(:,:,1:nz,1:npcon),   &
!              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!else ! No SCALARS
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz)
!end if !LV2
!
!else if (USE_MPI) then
!
!if (S_FLAG .AND. PCon_FLAG) then !LV2 ! With scalars and pollen - Ying 11/01/2010
!!!$  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
!!!$              theta(:,:,1:nz), PCon(:,:,1:nz,1:npcon), &
!!!$              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!!!$              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!time_tot = nsteps*dt_dim
!
!505 format(12e15.7)
!call getcwd(cwd)
!
!!write(*,*) TRIM(cwd)
!
!!write(folder,'(A,A,A,I2.2,A,I2.2,A)')TRim(cwd),'/','output_',time_tot,'_',npcon-1,'drops'
!!write(folder,'(A,A,A,I2.2,A,I4.4,A)')TRim(cwd),'/','output_3dsmooth',time_tot,'_',int(C1),'forcing'
!!write(vel,'(A,A,A,I2.2,A)')TRIM(folder),'/','vel_',time_tot,'/'
!!write(conc,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','conc_',time_tot,'_',npcon-1,'drops/'
!!write(Re_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Re_data_',time_tot,'_',npcon-1,'drops/'
!!write(brk_freq,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','brk_freq_',time_tot,'_',npcon-1,'drops/'
!!write(Rhs_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Rhs_',time_tot,'_',npcon-1,'drops/'
!
!
!write(folder,'(A,I3.3,A,I2.2,A)')'output_',time_tot,'_',npcon-1,'_c600'
!write(vel,'(A,A,A,I3.3,A)')TRIM(folder),'/','vel_',time_tot,'/'
!write(conc,'(A,A,A,I3.3,A)')TRIM(folder),'/','conc_',time_tot,'/'
!write(Re_data,'(A,A,A,I3.3,A)') TRIM(folder),'/','Re_data_',time_tot,'/'
!write(brk_freq,'(A,A,A,I3.3,A)')TRIM(folder),'/','brk_freq_',time_tot,'/'
!write(Rhs_data,'(A,A,A,I3.3,A)')TRIM(folder),'/','Rhs_',time_tot,'/'
!
!
!!write(*,*) TRIM(vel)
!!write(*,*) Trim(conc)
!
!!write(*,*) LEN(conc)
!
!write(vel_xz,'(A,A,I6.6,A)') TRIM(vel),'vel_',jt_total,'.dat'
!write(conc_xz,'(A,A,I6.6,A)') TRIM(conc),'conc_',jt_total,'.dat'
!!write(*,*) vel_xz
!!write(*,*) TRIM(conc_xz)
!  if(coord.eq.0) then !LV3
!     open(51,file=TRIM(conc_xz))
!     open(52,file=TRIM(vel_xz))
!!AA Added Pcon4 and 5
!     write(51,*) 'variables=x,z,PCon1,PCon2,PCon3,PCon4,PCon5,PCon6,PCon7,PCon8,PCon9,PCon10'
!     write(51,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!     write(52,*) 'variables=x,z,u,v,w,theta,dissipation'
!     write(52,*) 'zone t="',jt_total,'" i=',nx,'k=',nz_tot,'f=point'
!     if(yz_plane) then
!       write(conc_yz,'(A,A,I6.6,A)') TRIM(conc),'concyz_',jt_total,'.dat'
!       write(vel_yz,'(A,A,I6.6,A)') TRIM(vel),'vel_yz_',jt_total,'.dat'
!       open(53,file=TRIM(conc_yz))
!       open(54,file=TRIM(vel_yz))
!       write(53,*) 'variables=y,z,PCon1,PCon2,PCon3,PCon4,PCon5,PCon6,PCon7,PCon8,PCon9,PCon10'
!       write(53,*) 'zone t="',jt_total,'" i=',ny,' k=',nz_tot,'f=point'
!       write(54,*) 'variables=y,z,u,v,w,theta,dissipation'
!       write(54,*) 'zone t="',jt_total,'"i=',ny,'k=',nz_tot,'f=point'
!     endif
!    if(PBE_FLAG) then  !LV4
!     
!        do jz=1,nz_tot
!           do jx=1,nx
!              write(51,505)(jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot,(Pcon_tot(jx,jz,ipcon),ipcon=1,npcon-1)
!           enddo
!        enddo
!        if(yz_plane) then
!        do jz=1,nz_tot
!           do jy=1,ny
!              write(53,505)(jy-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,(Pcon_tot_yz(jy,jz,ipcon),ipcon=1,npcon-1)
!           enddo
!        enddo
!        endif
!     endif
!     if(OCEAN_FLAG) then
!        do jz=1,nz_tot
!          do jx=1,nx
!              write(52,505)(jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot,&
!                            u_tot(jx,jz),v_tot(jx,jz),-w_tot(jx,jz), &
!                   2._rprec-theta_tot(jx,jz),dissip_tot(jx,jz)!,Cs_opt2_tot(jx,jz),Nu_t_tot(jx,jz)
!           enddo
!        enddo
!        if(yz_plane) then
!        do jz=1,nz_tot
!          do jy=1,ny
!              write(54,505)(jy-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,&
!                            u_tot_yz(jy,jz),v_tot_yz(jy,jz),-w_tot_yz(jy,jz), &
!                   2._rprec-theta_tot_yz(jy,jz),dissip_tot_yz(jy,jz)!,Cs_opt2_tot(jx,jz),Nu_t_tot(jx,jz)
!           enddo
!        enddo
!        endif
!     endif
!!     else
! !       do jz=1,nz_tot
!  !         do jx=1,nx
!   !           write(11000000+jt_total,505)(jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot, &
!      !             !                PCon_tot(jx,jz,1),PCon_tot(jx,jz,2),PCon_tot(jx,jz,3), &
!     !              PCon_tot(jx,jz,1), PCon_tot(jx,jz,2), &
!    !               u_tot(jx,jz),v_tot(jx,jz),w_tot(jx,jz), &
!   !                theta_tot(jx,jz)!,Cs_opt2_tot(jx,jz),Nu_t_tot(jx,jz)
!  !         enddo
! !       enddo
! !    endif !LV4
!!AA Chang ends here    
!    close(51)
!    close(52)
!    close(53)
!    close(54)
!
!  endif !LV3
!
!!DY End here
!elseif (S_FLAG) then !WITH SCALARS
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
!              theta(:,:,1:nz)
!!DY Added by Di Yang for testing
!  if(coord.eq.0) then !LV3
!!  write(100,*) 'variables=x,y,z,PCon'
!     open(1000000+jt_total)
!!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta,Cs_opt2,Nu_t'
!     write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta'
!     write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!     do jz=1,nz_tot
!        do jx=1,nx
!           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot,&
!                u_tot(jx,jz),v_tot(jx,jz),w_tot(jx,jz),theta_tot(jx,jz)!, &
!!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
!        enddo
!     enddo
!     close(1000000+jt_total)
!  endif !LV3
!!DY End here
!elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
!              PCon(:,:,1:nz,1:npcon),   &
!              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!!DY Added by Di Yang for testing
!  if(coord.eq.0) then   !LV3
!!  write(100,*) 'variables=x,y,z,PCon'
!     open(1000000+jt_total)
!!  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,Pcon2,Pcon3,Cs_opt2,Nu_t'
!!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,Cs_opt2,Nu_t'
!     write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,PCon3'
!     write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!     do jz=1,nz_tot
!        do jx=1,nx
!           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot, &
!!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
!                PCon_tot(jx,jz,1), PCon_tot(jx,jz,2), PCon_tot(jx,jz,3)!, &
!!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
!        enddo
!     enddo
!     close(1000000+jt_total)
!  endif !LV3
!!DY End here
!else ! No SCALARS
!  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz)
!end if !LV2
!
!
!if (PBE_FLAG) then !LV2
!
!
!  if (coord.eq.0) then !LV3
!    write(str,'(A,A,I6.6,A)') TRIM(Rhs_data),'RHS_pbe.',jt_total,'.dat'
!    open(15,FILE=TRIM(str))
!    write(15,*)'variables=x,z,Rhs1,Rhs2,Rhs3,Rhs4,Rhs5,Rhs6,Rhs7,Rhs8,Rhs9,Rhs10'
!    write(15,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jx=1,nx
!        write(15,505)(jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot,(RHS_PBE_TOT(JX,JZ,ipcon),IPCON=1,NPCON-1)  
!      enddo
!    enddo
!    write(str1,'(A,A,I6.6,A)') TRIM(Re_data),'Parameters_',jt_total,'.dat'
!    open(104,FILE=TRIM(str1))
!    write(104,*)'variables=x,z,Re1,Re2,Re3,Re4,Re5,Re6,Re7,Re8,Re9,Re10'
!    write(104,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jx=1,nx
!        write(104,505)(jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot,(RE_TOT(JX,JZ,ipcon),IPCON=1,NPCON-1)
!         !   Re_tot(jx,jz,1),Re_tot(jx,jz,2),Re_tot(jx,jz,3),Re_tot(jx,jz,4),Re_tot(jx,jz,5),&
!         !   Re_tot(jx,jz,6),Re_tot(jx,jz,7),Re_tot(jx,jz,8),Re_tot(jx,jz,9),Re_tot(jx,jz,10)  
!      enddo
!    enddo
!    write(str2,'(A,A,I6.6,A)') TRIM(BRK_FREQ),'brk_freq_',jt_total,'.dat'
!    open(105,FILE=TRIM(str2))
!    write(105,*)'variables=x,z,br1,br2,br3,br4,br5,br6,br7,br8,br9,br10'
!    write(105,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jx=1,nx
!        write(105,505)(jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot,(BREAKAGE_FREQ_TOT(JX,JZ,ipcon),IPCON=1,NPCON-1)
!!            breakage_freq_tot(jx,jz,1),breakage_freq_tot(jx,jz,2),breakage_freq_tot(jx,jz,3),breakage_freq_tot(jx,jz,4),breakage_freq_tot(jx,jz,5),&
!!            breakage_freq_tot(jx,jz,6),breakage_freq_tot(jx,jz,7),breakage_freq_tot(jx,jz,8),breakage_freq_tot(jx,jz,9),breakage_freq_tot(jx,jz,10)
!      enddo
!    enddo
!
!   
!  endif  !LV3
!  close(104)
!  close(105)
!  close(15)
!
!
!
!
! if (coord.eq.0) then !LV3
!!    write(str3,'(A,I6.6,A)') 'RHS_pbe_yz.',jt_total,'.dat'
!!    open(106,FILE=TRIM(str3))
!!    write(106,*)'variables=y,z,Rhs1,Rhs2,Rhs3,Rhs4,Rhs5,Rhs6,Rhs7,Rhs8,Rhs9,Rhs10'
!!    write(106,*) 'zone t="',jt_total,'" i=',ny,' k=',nz_tot,' f=point'
!!    do jz=1,nz_tot
!!      do jx=1,ny
!!        write(106,505) (jx-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot, &
!!        Rhs_pbe_tot_yz(jx,jz,1),RHS_pbe_tot_yz(jx,jz,2),RHS_pbe_tot_yz(jx,jz,3),RHS_pbe_tot_yz(jx,jz,4),RHS_pbe_tot_yz(jx,jz,5),&
!!        RHS_pbe_tot_yz(jx,jz,6),RHS_pbe_tot_yz(jx,jz,7),RHS_pbe_tot_yz(jx,jz,8),RHS_pbe_tot_yz(jx,jz,9),RHS_pbe_tot_yz(jx,jz,10)
!!      enddo
!!    enddo
!!AA
!if(yz_plane) then
!    write(str4,'(A,A,I6.6,A)') TRIM(Re_data), 'Parameters_yz_',jt_total,'.dat'
!    open(107,FILE=TRIM(str4))
!    write(107,*)'variables=y,z,Re1,Re2,Re3,Re4,Re5,Re6,Re7,Re8,Re9,Re10'
!    write(107,*) 'zone t="',jt_total,'" i=',ny,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jx=1,ny
!        write(107,505) (jx-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,&
!            Re_tot_yz(jx,jz,1),Re_tot_yz(jx,jz,2),Re_tot_yz(jx,jz,3),Re_tot_yz(jx,jz,4),Re_tot_yz(jx,jz,5),&
!            Re_tot_yz(jx,jz,6),Re_tot_yz(jx,jz,7),Re_tot_yz(jx,jz,8),Re_tot_yz(jx,jz,9),Re_tot_yz(jx,jz,10)
!      enddo
!    enddo
!    write(str5,'(A,A,I6.6,A)') Trim(brk_freq),'brk_freq_yz_',jt_total,'.dat'
!    open(108,FILE=TRIM(str5))
!    write(108,*)'variables=y,z,br1,br2,br3,br4,br5,br6,br7,br8,br9,br10'
!    write(108,*) 'zone t="',jt_total,'" i=',ny,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jx=1,ny
!        write(108,505) (jx-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,&
!            breakage_freq_tot_yz(jx,jz,1),breakage_freq_tot_yz(jx,jz,2),breakage_freq_tot_yz(jx,jz,3),breakage_freq_tot_yz(jx,jz,4),breakage_freq_tot_yz(jx,jz,5),&
!            breakage_freq_tot_yz(jx,jz,6),breakage_freq_tot_yz(jx,jz,7),breakage_freq_tot_yz(jx,jz,8),breakage_freq_tot_yz(jx,jz,9),breakage_freq_tot_yz(jx,jz,10)
!      enddo
!    enddo
!endif
!!AA
!  endif  !LV3
!  if (check_neg) then !LV4
!    write(str3,'(A,I6.6,A)') 'pcon_neg_',jt_total,'.dat'
!    open(106,file=TRIM(str3))
!    write(106,*) 'variables = x,y,z,Pcon_neg'
!    write(106,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
!    do jz=1,nz_tot
!      do jy=1,ny
!        do jx=1,nx
!          write(106,*),(jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,Pcon_neg_tot(jx,jy,jz)
!        enddo
!      enddo
!    enddo
! ! close(105)
!    close(106)
!  endif !LV4
!  close(107)
!  close(108)
!
!endif !Lv2
!
!
!
!endif !LV1
!
!end subroutine checkpoint_v2
!
!AA added subroutine for averages



subroutine calculate_avgs()
use param
use sim_param, only :u, PCon
use scalars_module, only:RHS_pbe,Re,breakage_freq
use bottombc, only: ustar_avg
use sgsmodule
integer jx,jy,jz,ipcon

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
real (rprec), dimension(ld,1:nz, npcon) :: PCon1,PCon2,PCon3,breakage_freq1,Re_1,tau_c1,RHS_pbe1
real (rprec), dimension(ld,1:nz_tot, npcon) ::PCon_tot1,PCon_tot2,PCon_tot3,breakage_freq_tot,Re_tot,tau_c_tot,&
RHS_pbe_tot
real (rprec), dimension(ld,ny,1:nz_tot) :: Pcon_neg_tot
!real (rprec), dimension(1:nz) :: dissip1,dissip2,dissip3
!real (rprec), dimension(1:nz) :: dissip_tot1,dissip_tot2,dissip_tot3
real (rprec), dimension(ld,1:nz) :: u1,dissip1
real (rprec), dimension(ld,1:nz_tot) :: dissip_tot
real (rprec), dimension(npcon) :: sum1=0._rprec,sum2=0._rprec,sum3=0._rprec
character (LEN=50) :: str
  u1 = u(:,yps,1:nz)
  PCon1 = PCon(:,yps,1:nz,:)
  breakage_freq1 = breakage_freq(:,yps,1:nz,:)
  dissip1 = dissip(:,yps,1:nz)
  Re_1 = Re(:,yps,1:nz,:)
  RHS_pbe1 = RHS_pbe(:,yps,1:nz,:)
!  tau_c1 = tau_c(:,yps,1:nz,:)
!  PCon2 = PCon(,yps,1:nz,:)
!  PCon3 = PCon(ld,yps,1:nz,:)
!  dissip1 = dissip(xps,yps,1:nz)
!  dissip2 = dissip(xps+4,yps,1:nz)
!  dissip3 = dissip(xps-4,yps,1:nz)
  sendcounts = size (u1(:,1:nz))
  recvcounts = size (u1(:,1:nz-1))
  displs = coord_of_rank * recvcounts

!  sendcounts_neg = size(u(:,:,1:nz))
!  recvcounts_neg = size(u(:,:,1:nz-1))
!  displs = coord_of_rank * recvcounts
!AA Gather negative conc array

!  call mpi_gatherv (Pcon_neg(1,1,1), sendcounts_neg, MPI_RPREC,&  
!                    Pcon_neg_tot(1,1,1), sendcounts_neg, displs_neg,&
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)


  call mpi_gatherv (dissip1(1,1), sendcounts, MPI_RPREC,&
                    dissip_tot(1,1), sendcounts, displs, &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!  call mpi_gatherv (dissip2(1), sendcounts, MPI_RPREC,&
!                    dissip_tot2(1), sendcounts, displs, &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!  call mpi_gatherv (dissip3(1), sendcounts, MPI_RPREC,&
!                    dissip_tot3(1), sendcounts, displs, &
!                    MPI_RPREC, rank_of_coord(0), comm, ierr)
!


  if (PCon_FLAG) then !LV2
     do ipcon=1,npcon
        call mpi_gatherv (PCon1(1,1,ipcon), sendcounts, MPI_RPREC,&
             PCon_tot1(1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)

        call mpi_gatherv (breakage_freq1(1,1,ipcon), sendcounts, MPI_RPREC,&
           breakage_freq_tot(1,1,ipcon), sendcounts, displs,             &
           MPI_RPREC, rank_of_coord(0), comm, ierr)
        call mpi_gatherv (Re_1(1,1,ipcon), sendcounts, MPI_RPREC,&
             Re_tot(1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
        call mpi_gatherv (RHS_pbe1(1,1,ipcon), sendcounts, MPI_RPREC,&
             RHS_pbe_tot(1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
!        call mpi_gatherv (tau_c1(1,1,ipcon), sendcounts, MPI_RPREC,&
!             tau_c_tot(1,1,ipcon), sendcounts, displs,       &
!             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  endif !LV2
606 format(12e15.7)
do ipcon=1,npcon
  do jz=1,nz_tot
    sum1(ipcon) = sum1(ipcon) + PCon_tot1(xps,jz,ipcon)
    sum2(ipcon) = sum2(ipcon) + PCon_tot1(xps+4,jz,ipcon)
    sum3(ipcon) = sum3(ipcon) + Pcon_tot1(xps-4,jz,ipcon)
  enddo
enddo
sum1 = sum1/nz_tot
sum2 = sum2/nz_tot
sum3 = sum3/nz_tot
!  sum2 = sum(PCon_tot1(xps+4,:,:),Dim=2)/nz_tot
!  sum3 = sum(PCon_tot1(xps-4,:,:),Dim=2)/nz_tot
if (coord .EQ. 0) then
  open(98,file="Pconjetz4.dat",status='unknown',position='append')
  write(98,606) (87-1)*z_i/nz_tot , (PCon_tot1(xps,87,ipcon),ipcon=1,npcon-1)

  open(97,file="Pconjetz7.dat",status='unknown',position='append')
  write(97,606) (158-1)*z_i/nz_tot , (PCon_tot1(xps,158,ipcon),ipcon=1,npcon-1)

  open(96,file="Pconjetz5.dat",status='unknown',position='append')
  write(96,606) (108-1)*z_i/nz_tot , (PCon_tot1(xps,108,ipcon),ipcon=1,npcon-1)

  open(60,file = "Breakjetz4.dat",status='unknown',position='append')
  write(60,606) (87-1)*z_i/nz_tot , (breakage_freq_tot(xps,87,ipcon),ipcon=1,npcon-1)

  open(61,file = "Breakjetz7.dat",status='unknown',position='append')
  write(61,606) (158-1)*z_i/nz_tot , (breakage_freq_tot(xps,158,ipcon),ipcon=1,npcon-1)

  open(62,file = "Breakjetz5.dat",status='unknown',position='append')
  write(62,606) (108-1)*z_i/nz_tot ,(breakage_freq_tot(xps,108,ipcon),ipcon=1,npcon-1)


   open(93,file="dissipz4.dat",status='unknown',position='append')
  write(93,606) (87-1)*z_i/nz_tot , dissip_tot(xps,87)

  open(94,file="dissipz7.dat",status='unknown',position='append')
  write(94,606) (158-1)*z_i/nz_tot , dissip_tot(xps,158)

  open(95,file="dissipz5.dat",status='unknown',position='append')
  write(95,606) (108-1)*z_i/nz_tot , dissip_tot(xps,108)

  open(90,file = "Rez4.dat",status='unknown',position='append')
  write(90,606) (87-1)*z_i/nz_tot ,(Re_tot(xps,87,ipcon),ipcon=1,npcon-1)

  open(91,file = "Rez7.dat",status='unknown',position='append')
  write(91,606) (158-1)*z_i/nz_tot ,(Re_tot(xps,158,ipcon),ipcon=1,npcon-1)

  open(92,file = "Rez5.dat",status='unknown',position='append')
  write(92,606) (108-1)*z_i/nz_tot,(Re_tot(xps,108,ipcon),ipcon=1,npcon-1)

 ! write(str,'(A,I6.6,A)') 'pcon_neg_',jt_total,'.dat'
 ! open(105,file=TRIM(str))
 ! write(105,*) 'variables = x,y,z,Pcon_neg'
 ! write(105,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
 ! do jz=1,nz_tot
!    do jy=1,ny
!      do jx=1,nx
!        write(105,*),(jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,-(jz-1)*z_i/nz_tot,Pcon_neg_tot(jx,jy,jz)
!      enddo
!    enddo
!  enddo
!  close(105)
  close(90)
   close(91)
 close(92)
 close(93)
 close(94)
 close(95)

  close(60)
  close(61)
  close(62)
  close(96)
  close(97)
  close(98)
endif
end subroutine calculate_avgs



!+++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY Output the oil flux on a plane near the bottom of OML
!DY For oil plume inflow condition to ENDLESS
!+++++++++++++++++++++++++++++++++++++++++++++++
subroutine output_plane_oil_flux
use param
use sim_param, only : w, PCon

implicit none

integer jx,jy,jz,ipcon

integer, parameter :: jz_p1 = 9, jz_p2 = 18
integer, parameter :: jx_min = 65, jx_max = 172
integer, parameter :: jy_min = 1, jy_max = ny
integer :: jz_p1_local, jz_p2_local, icpu_p1, icpu_p2

icpu_p1 = int(jz_p1/(nz-1))
icpu_p2 = int(jz_p2/(nz-1))
jz_p1_local = jz_p1 - icpu_p1*(nz-1)+1
jz_p2_local = jz_p2 - icpu_p2*(nz-1)+1

if(coord.eq.0) then
   print*, icpu_p1,icpu_p2,jz_p1_local,jz_p2_local
endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == icpu_p1)) then
   open(31000000+jt_total)
   write(31000000+jt_total,*) 'variables=x,y,z,PCon2,w'
   write(31000000+jt_total,*) 'zone t="',jt_total,'" i=',jx_max-jx_min+1,' j=',jy_max-jy_min+1,' k=',1,' f=point'
   do jy = jy_min, jy_max
      do jx = jx_min, jx_max
         write(31000000+jt_total,*) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz_p1-1)*z_i/nz_tot, &
              PCon(jx,jy,jz_p1_local,2), (-w(jx,jy,jz_p1_local)-w(jx,jy,jz_p1_local+1))/2._rprec
      enddo
   enddo
   close(31000000+jt_total)
endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == icpu_p2)) then
   open(32000000+jt_total)
   write(32000000+jt_total,*) 'variables=x,y,z,PCon2,w'
   write(32000000+jt_total,*) 'zone t="',jt_total,'" i=',jx_max-jx_min+1,' j=',jy_max-jy_min+1,' k=',1,' f=point'
   do jy = jy_min, jy_max
      do jx = jx_min, jx_max
         write(32000000+jt_total,*) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz_p2-1)*z_i/nz_tot, &
              PCon(jx,jy,jz_p2_local,2), (-w(jx,jy,jz_p2_local)-w(jx,jy,jz_p2_local+1))/2._rprec
      enddo
   enddo
   close(32000000+jt_total)
endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   open(41000000+jt_total)
   write(41000000+jt_total,*) 'variables=x,y,z,PCon2,w'
   write(41000000+jt_total,*) 'zone t="',jt_total,'" i=',jx_max-jx_min+1,' j=',jy_max-jy_min+1,' k=',1,' f=point'
   do jy = jy_min, jy_max
      do jx = jx_min, jx_max
         write(41000000+jt_total,*) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(1-1)*z_i/nz_tot, &
              PCon(jx,jy,1,2), w(jx,jy,1)
      enddo
   enddo
   close(41000000+jt_total)
   open(42000000+jt_total)
   write(42000000+jt_total,*) 'variables=x,y,z,PCon2,w'
   write(42000000+jt_total,*) 'zone t="',jt_total,'" i=',jx_max-jx_min+1,' j=',jy_max-jy_min+1,' k=',1,' f=point'
   do jy = jy_min, jy_max
      do jx = jx_min, jx_max
         write(42000000+jt_total,*) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(2-1)*z_i/nz_tot, &
              PCon(jx,jy,2,2), w(jx,jy,2)
      enddo
   enddo
   close(42000000+jt_total)
   open(43000000+jt_total)
   write(43000000+jt_total,*) 'variables=x,y,z,PCon2,w'
   write(43000000+jt_total,*) 'zone t="',jt_total,'" i=',jx_max-jx_min+1,' j=',jy_max-jy_min+1,' k=',1,' f=point'
   do jy = jy_min, jy_max
      do jx = jx_min, jx_max
         write(43000000+jt_total,*) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(3-1)*z_i/nz_tot, &
              PCon(jx,jy,3,2), w(jx,jy,3)
      enddo
   enddo
   close(43000000+jt_total)
endif

end subroutine output_plane_oil_flux



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun

integer, parameter :: lun_default = 11

integer::jx,jy,jz
integer :: lun

logical :: opn

!---------------------------------------------------------------------

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Modified by Di Yang
!DY Original version requires large memory for high-resolution case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!DY if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

!DY inquire (unit=lun, opened=opn)

!DY if (.not. opn) then
!DY  write (*, *) 'output_final: lun=', lun, ' is not open'
!DY  stop
!DY end if

!DY rewind (lun)

!DY endif

!DY call checkpoint_final (lun)
call checkpoint_final_v2
!+++++++++++++++
!DY End here
!+++++++++++++++

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
close (lun)
endif

if ((cumulative_time) .and. (lun == lun_default)) then
  !--only do this for true final output, not intermediate recording
  open (1, file=fcumulative_time)
!  if ((DYN_init .ne. 1) .AND. (GABLS_diurnal_test)) jt_total=jt_total_init
  write (1, *) jt_total
  close (1)
end if

end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_lambda2_out
use sim_param,only:path
implicit none
character(len=24)::fname
call lambda2()
write(fname,'(A13,i6.6,A4)')path//'output/lam-',jt_total,'.out'
open(1,file=fname,form='unformatted')
write(1)nx,ny,nz
write(1)real(lam2)
close(1)
end subroutine io_lambda2_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lambda2()
use types,only:rprec
use sim_param,only:u,v,w,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use param,only:dx,dy,dz
implicit none
!TSreal(kind=rprec),dimension(nx)::dpdx,dpdy,dpdz,aeu,awu,anu,asu,afu,abu
!TSreal(kind=rprec)::fne,fnn,fnf,S11,S22,S33,S12,O12,S13,O13,S23,O23
real(kind=rprec)::S11,S22,S33,S12,O12,S13,O13,S23,O23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz
integer::jx,jy,jz
! following used for eispack call...
integer::neis,nmeis,matzeis,ierreis,iv1eis(3)
double precision,dimension(3,3)::aeis,zeis
double precision,dimension(3)::wreis,wieis,fv1eis
double precision::ave

! assignments for eispack call
neis=3
nmeis=3
matzeis=0
ierreis=0
!TSfne=1._rprec/dx
!TSfnn=1._rprec/dy
!TSfnf=1._rprec/dz
!TSprint*,fne,fnn,fnf
lam2=0._rprec

! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
jz=1
do jy=1,ny
do jx=1,nx              
   ux=dudx(jx,jy,1)  ! uvp-node
   uy=dudy(jx,jy,1)  ! uvp-node
   uz=dudz(jx,jy,1)  ! uvp-node
   vx=dvdx(jx,jy,1)  ! uvp-node
   vy=dvdy(jx,jy,1)  ! uvp-node
   vz=dvdz(jx,jy,1)  ! uvp-node 
! special case
   wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
   wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
   wz=dwdz(jx,jy,1)  ! uvp-node
   S11=ux          ! uvp-node
   S12=0.5_rprec*(uy+vx) ! uvp-node
! taken care of with wall stress routine
   S13=0.5_rprec*(uz+wx) ! uvp
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! uvp-node
! taken care of with wall stress routine 
   S23=0.5_rprec*(vz+wy) ! uvp
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! uvp-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
  write (*, *) 'rg temporarily removed, sorry'; stop
  !call rg(nmeis,neis,aeis,wreis,wieis,matzeis,zeis,iv1eis,fv1eis,ierreis)
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
! calculate derivatives/strain on w-nodes
do jz=2,nz-1  
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11=ux          ! w-node
   S12=0.5_rprec*(uy+vx) ! w-node
   S13=0.5_rprec*(uz+wx) ! w-node
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! w-node
   S23=0.5_rprec*(vz+wy) ! w-node
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! w-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
   !call rg(nmeis,neis,aeis,wreis,wieis,matzeis,zeis,iv1eis,fv1eis,ierreis)
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
end do

print*,'minmax',minval(lam2),maxval(lam2)
end subroutine lambda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_mean_out
implicit none
write(51)real(mean_u),real(mean_u2),real(mean_v),real(mean_v2),&
     real(mean_w),real(mean_w2)
!--the nz/4*3 stuff has to go
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
end subroutine io_mean_out

subroutine calculate_mean
use sim_param,only:u,v,w
use sgsmodule,only:Cs_opt2,Cs_opt2_avg
implicit none
Cs_opt2_avg(:,:,:)=Cs_opt2_avg(:,:,:)+Cs_opt2(:,:,:)/nwrite
!TS
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
end subroutine calculate_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timeseries_spec
use sim_param,only:u,v,w,theta
implicit none
integer::jx,jy,jz,i
if(mod(jt_total,time_spec)==0.and.jt_total.gt.2000)then
jx=NX/8
jy=NY/2+1
jz=NZ/2
!TSwrite(15)real(u(jx+NX/24*2,jy,jz:jz+3)),real(u(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(u(jx+NX/24*6,jy,jz:jz+3)),&
!TS     real(v(jx+NX/24*2,jy,jz:jz+3)),real(v(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(v(jx+NX/24*6,jy,jz:jz+3)),&
!TS     real(w(jx+NX/24*2,jy,jz:jz+3)),real(w(jx+NX/24*4,jy,jz:jz+3)),&
!TS     real(w(jx+NX/24*6,jy,jz:jz+3))
endif
end subroutine timeseries_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine post_spec(jt_local)
use sim_param,only:path,u,v,w,theta,PCon
use param
use fft
implicit none
real(kind=rprec),dimension(nx/2,nz)::spectra_u,spectra_v,spectra_w,&
     spectra_theta
real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
integer,intent(in)::jt_local
integer::k,jz!,z

real(kind=rprec),dimension(nz-1)::z   ! Ying changed z to array (09/30/2010)
                                      ! Chamecki changed z from integer to real (08/04/2006)
real(kind=rprec),dimension(nz_tot-1)::z_tot

!character(len=25)::fname1,fname2,fname3,fname4
character(len=64)::fname1,fname2,fname3,fname4
$if ($MPI)
  $define $lbz 0
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$else
  $define $lbz 1
$endif

$if ($MPI)
  recvcounts = size (z)
  displs = coord_of_rank * recvcounts
  do jz=1,nz-1
    z(jz) = (coord*(nz-1)+jz-0.5_rprec)*dz*z_i
  enddo
  call mpi_gatherv (z(1), size (z), MPI_RPREC,&
                    z_tot(1), recvcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
$else
  do jz=1,nz-1
    z(jz)=(jz-0.5_rprec)*dz*z_i
    z_tot(jz)=z(jz)
  enddo
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  write(fname1,'(a,a)') path//'output/spec_x','.dat'
  open(82,file=fname1,form='formatted')
  do jz = 1,nz_tot-1
    write(82,*) (real(kx(k,1)/z_i*z_tot(jz)),k=1,nx/2)  
  enddo
  close(82)
endif

!write(fname1,'(a,a)') path//'output/spec_x','.dat'
!open(82,file=fname1,form='formatted')
do jz=1,nz-1
!do jz=$lbz,nz-1
!   z=(jz-0.5_rprec)*dz*z_i
!   write(82,*) (real(2*pi/L_x*k/z_i*z),k=1,nx/2-1)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! kx=(2pi/L_x)*(0:nx/2) => L_x is already non-dim by z_i
!! => k_x=[(2pi/L_x)*(0:nx/2)/z_i] and 
!! kz is given by = k_x*z => kz=[(2pi/L_x)*(0:nx/2)*([1:nz]-0.5)*dz]
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   write(82,*) (real(kx(k,1)/z_i*z),k=1,nx/2)
! Calculate spectrum now !! for each slice in MPI case !   
!   call spectrum(u(:, :, jz), spectra_u(:,jz))
!   call spectrum(v(:, :, jz), spectra_v(:,jz))
!   call spectrum(w(:, :, jz), spectra_w(:,jz))
!   call spectrum(theta(:,:,jz),spectra_theta(:,jz))
   call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
   call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
   call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
   call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))
   
   ! Replace temperature spectra by Pollen concentration
   ! Chamecki - 08/10/2006
   IF (PCon_FLAG) call spectrum(PCon(:, :, jz,npcon), spectra_uvwT(4,:,jz))
   
enddo
!   close(82)
!   print *,'spectra_sample_U',spectra_u(:,4)
!   print *,'spectra_sample_U2',spectra_uvwT(1,:,4)
!   print *,'spectra_sample_V',spectra_v(:,4)
!   print *,'spectra_sample_V2',spectra_uvwT(2,:,4)
!   print *,'spectra_sample_W',spectra_w(:,4)
!   print *,'spectra_sample_W2',spectra_uvwT(3,:,4)
!   print *,'spectra_sample_T',spectra_theta(:,4)
!   print *,'spectra_sample_T2',spectra_uvwT(4,:,4)
$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
                    spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
$else
  spectra_uvwT_tot=spectra_uvwT
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
   write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
   open(83,file=fname1,form='unformatted')

!write(fname1,'(a,i6.6,a)')path//'output/spec_u',jt_local,'.dat'
!open(1,file=fname1,form='formatted')
!write(fname2,'(A,i6.6,A)')path//'output/spec_v',jt_local,'.dat'
!open(2,file=fname2,form='formatted')
!write(fname3,'(A,i6.6,A)')path//'output/spec_w',jt_local,'.dat'
!open(3,file=fname3,form='formatted')
!write(fname4,'(A,i6.6,A)')path//'output/spec_t',jt_local,'.dat'
!open(4,file=fname4,form='formatted')

   write(83) real(spectra_uvwT_tot(:,1:nx/2,:))
   close(83)
!do jz=1,nz
!   write(1,*)real(spectra_u(2:nx/2,jz))
!   write(2,*)real(spectra_v(2:nx/2,jz))
!   write(3,*)real(spectra_w(2:nx/2,jz))
!   write(4,*)real(spectra_theta(2:nx/2,jz))
!enddo
!close(1);close(2);close(3);close(4)
end if

end subroutine post_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum(u, spec)
use fft
implicit none      
real(kind=rprec),dimension(ld,ny),intent(in)::u
real(kind=rprec),dimension(nx/2),intent(out)::spec  !--assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8, save :: plan
logical, save :: init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

! initialize
spec(:)=0._rprec
do jy=1,ny
   vel_r(:)= u(1:nx,jy)/real(nx,kind=rprec)
! check this normaliztion-part of forward
! call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)
! compute magnitudes
! the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
   spec(1)=spec(1)+0.5*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
!        print *,'k,vel,spec',k,vel_c(k),spec(k)
   end do

   !--assume Nyquist is 0
   !spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
end do
spec(:)=spec(:)/real(Ny,kind=rprec) ! for average over Ny
end subroutine spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_stats ()
use param
use sim_param, only : u, v, w, txz, theta
use sgsmodule
use fft, only : kx
implicit none

!--choose naming convention that does not conflict with qpost
character (*), parameter :: fubar_avg = 'output/ubar-avg_stats.dat'
character (*), parameter :: fvbar_avg = 'output/vbar-avg_stats.dat'
character (*), parameter :: fupr2bar_avg = 'output/upr2bar-avg_stats.dat'
character (*), parameter :: fstressbar_avg = 'output/stressbar-avg_stats.dat'
character (*), parameter :: fEozbar_avg = 'output/Eozbar-avg_stats.dat'
character (*), parameter :: fthetabar_avg = 'output/thetabar-avg_stats.dat'
character (*), parameter :: fCs2bar_avg = 'output/Cs2bar-avg_stats.dat'
character (*), parameter :: fNutbar_avg = 'output/Nutbar-avg_stats.dat'

!integer, parameter :: nz_tot = (nz - 1) * nproc + 1
integer, parameter :: hdr_len = 256

logical, parameter :: DEBUG = .false.

character (hdr_len) :: Eozbar_hdr

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer, save :: n_ubar_avg
integer, save :: n_vbar_avg
integer, save :: n_upr2bar_avg
integer, save :: n_stressbar_avg
integer, save :: n_Eozbar_avg
integer, save :: n_thetabar_avg
integer, save :: n_Cs2bar_avg
integer, save :: n_Nutbar_avg
integer :: jz

logical, save :: init = .false.

real (rprec) :: z
real (rprec) :: zu(1, nz_tot-1)
real (rprec) :: kz_z(2, nx/2)
real (rprec), save :: ubar_avg(1, nz_tot-1)  !--1 is <u>
real (rprec), save :: vbar_avg(1, nz_tot-1)  !--1 is <v>
real (rprec), save :: upr2bar_avg(3, nz_tot-1)  !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_avg(3, nz_tot-1)
                      !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_avg(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec), save :: thetabar_avg(1, nz_tot-1)  !--1 is <theta>
real (rprec), save :: Cs2bar_avg(1, nz_tot-1)
real (rprec), save :: Nutbar_avg(1, nz_tot-1)
!--tot is a temp for current stats at nz_tot size
real (rprec), save :: ubar_tot(1, nz_tot-1)  !--1 is <u>
real (rprec), save :: vbar_tot(1, nz_tot-1)  !--1 is <v>
real (rprec), save :: upr2bar_tot(3, nz_tot-1)  !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_tot(3, nz_tot-1)
                      !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_tot(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec), save :: thetabar_tot(1, nz_tot-1)  !--1 is <theta>
real (rprec), save :: Cs2bar_tot(1, nz_tot-1)
real (rprec), save :: Nutbar_tot(1, nz_tot-1)
real (rprec) :: upr(nx, ny), vpr(nx, ny), wpr(nx, ny)
real (rprec) :: ubar(nz-1), vbar(nz-1), wbar(nz-1), thetabar(nz-1)
real (rprec) :: Cs2bar(nz-1), Nutbar(nz-1)
real (rprec) :: upr2bar(3, nz-1)
real (rprec) :: stressbar(3, nz-1)
real (rprec) :: Eozbar(nx/2, nz-1)
real (rprec) :: ubar0,vbar0

!---------------------------------------------------------------------

!--check whether or not to actually do anything
!--motivation for doing this way is that it cleans up interface in main
if (modulo (jt, n_avg_stats) /= 0) goto 001  !--do nothing, exit cleanly

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  if (.not. init) then  !--initialization

    call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
    call init_avg (fvbar_avg, 1, vbar_avg, n_vbar_avg)
    call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
    call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg) 
    do jz = 1, nz-2
      call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                     leaveopn='yes')
    end do
    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)
    call init_avg (fthetabar_avg, 1, thetabar_avg, n_thetabar_avg)
    call init_avg (fCs2bar_avg, 1, Cs2bar_avg, n_Cs2bar_avg)
    call init_avg (fNutbar_avg, 1, Nutbar_avg, n_Nutbar_avg)

    init = .true.

  end if
end if

do jz = 1, nz-1

  ubar(jz) = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
  vbar(jz) = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
  wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)
  thetabar(jz) = sum (theta(1:nx, 1:ny, jz)) / (nx * ny)
  Cs2bar(jz) = sum (Cs_opt2(1:nx, 1:ny, jz)) / (nx * ny)
  Nutbar(jz) = sum (Nu_t(1:nx, 1:ny, jz)) / (nx * ny)

!+++++++++++++++++++++++++++++++++++
!DY Added by Di Yang, start here
enddo

ubar0 = sum (u(1:nx, 1:ny, 0)) / (nx * ny)
vbar0 = sum (v(1:nx, 1:ny, 0)) / (nx * ny)

do jz = 1, nz-1
!DY End here
!+++++++++++++++++++++++++++++++++++

   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
        (jz == 1) ) then
      upr = 0._rprec
      vpr = 0._rprec
      wpr = 0._rprec
   else
    !--see qpost for u/w-node interpolation
    !--convention will be to put up, vp, wp on w-nodes
!+++++++++++++++++++++++++++
!DY Modified by Di Yang
!DY In the old version, there is no physical value for ubar(0) for when using MPI
!DY      upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
!DY           u(1:nx, 1:ny, jz-1) - ubar(jz-1))
!DY      vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
!DY           v(1:nx, 1:ny, jz-1) - vbar(jz-1))
      if(jz==1) then
         upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
              u(1:nx, 1:ny, jz-1) - ubar0)
         vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
              v(1:nx, 1:ny, jz-1) - vbar0)
      else
         upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
              u(1:nx, 1:ny, jz-1) - ubar(jz-1))
         vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
              v(1:nx, 1:ny, jz-1) - vbar(jz-1))
      endif
      wpr = w(1:nx, 1:ny, jz) - wbar(jz)
   end if
 
  upr2bar(1, jz) = sum (upr**2) / (nx * ny)
  upr2bar(2, jz) = sum (vpr**2) / (nx * ny)
  upr2bar(3, jz) = sum (wpr**2) / (nx * ny)

  stressbar(1, jz) = sum (upr * wpr) / (nx * ny) 
  stressbar(2, jz) = sum (txz(1:nx, 1:ny, jz)) / (nx * ny)
  stressbar(3, jz) = sum (stressbar(1:2, jz))

  !--energy spectra
  call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
  z = (jz - 0.5_rprec) * dz
  Eozbar(:, jz) = Eozbar(:, jz) / z

end do

!--collect current stats into nz_tot sized arrays
$if ($MPI)

  if (DEBUG) then
    write (*, *) coord, ': ubar(1) = ', ubar(1)
  end if

  recvcounts = size (ubar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (ubar(1), size (ubar), MPI_RPREC,                &
                    ubar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (vbar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (vbar(1), size (vbar), MPI_RPREC,                &
                    vbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (upr2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (upr2bar(1, 1), size (upr2bar), MPI_RPREC,          &
                    upr2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)  

  recvcounts = size (stressbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (stressbar(1, 1), size (stressbar), MPI_RPREC,        &
                    stressbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Eozbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Eozbar(1, 1), size (Eozbar), MPI_RPREC,              &
                    Eozbar_tot(1, 1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (thetabar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (thetabar(1), size (thetabar), MPI_RPREC,                &
                    thetabar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Cs2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Cs2bar(1), size (Cs2bar), MPI_RPREC,                &
                    Cs2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Nutbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Nutbar(1), size (Nutbar), MPI_RPREC,                &
                    Nutbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

$else

  ubar_tot(1, :) = ubar

  vbar_tot(1, :) = vbar

  upr2bar_tot = upr2bar

  stressbar_tot = stressbar

  Eozbar_tot(1, :, :) = Eozbar

  thetabar_tot(1, :) = thetabar

  Cs2bar_tot(1, :) = Cs2bar

  Nutbar_tot(1, :) = Nutbar

$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !--calculation of cumulative average stats
  ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
  n_ubar_avg = n_ubar_avg + 1

  vbar_avg = (n_vbar_avg * vbar_avg + vbar_tot) / (n_vbar_avg + 1)
  n_vbar_avg = n_vbar_avg + 1

  upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
                (n_upr2bar_avg + 1)
  n_upr2bar_avg = n_upr2bar_avg + 1

  stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
                  (n_stressbar_avg + 1)
  n_stressbar_avg = n_stressbar_avg + 1

  Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
  n_Eozbar_avg = n_Eozbar_avg + 1

  thetabar_avg = (n_thetabar_avg * thetabar_avg + thetabar_tot) / (n_thetabar_avg + 1)
  n_thetabar_avg = n_thetabar_avg + 1

  Cs2bar_avg = (n_Cs2bar_avg * Cs2bar_avg + Cs2bar_tot) / (n_Cs2bar_avg + 1)
  n_Cs2bar_avg = n_Cs2bar_avg + 1

  Nutbar_avg = (n_Nutbar_avg * Nutbar_avg + Nutbar_tot) / (n_Nutbar_avg + 1)
  n_Nutbar_avg = n_Nutbar_avg + 1

  !--prepare list of z-coordinates
  forall (jz=1:nz_tot-1) zu(1, jz) = (jz - 0.5_rprec) * dz

  !--prepare  header, optional

  !--write out to file
  call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
  call write_avg (fvbar_avg, n_vbar_avg, zu, vbar_avg)
  call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
  call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)
  call write_avg (fthetabar_avg, n_thetabar_avg, zu, thetabar_avg)
  call write_avg (fCs2bar_avg, n_Cs2bar_avg, zu, Cs2bar_avg)
  call write_avg (fNutbar_avg, n_Nutbar_avg, zu, Nutbar_avg)

  !--this is a bit awkward: maybe just have another routine to do it right
  Eozbar_hdr = 'zone' !--this is for tecplot... 
  kz_z(1, :) = kx(1:nx/2, 1) * zu(1, 1)
  kz_z(2, :) = zu(1, 1)
  call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, 1),  &
                  hdr=Eozbar_hdr) 

  do jz = 2, nz_tot - 1

    kz_z(1, :) = kx(1:nx/2, 1) * zu(1, jz)
    kz_z(2, :) = zu(1, jz)
  
    call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, jz),  &
                    hdr=Eozbar_hdr, position='append') 

  end do

end if

001 continue  !--exit cleanly

end subroutine avg_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_ccol  !--num. coord columns: x, y, etc.

real (rprec), intent (out) :: a_avg(:, :)
integer, intent (out) :: n_avg

character (*), optional, intent (in) :: leaveopn

character (128) :: buff

logical :: exst, opn

integer :: j

real (rprec) :: z(n_ccol)

!---------------------------------------------------------------------

inquire (file=file_avg, exist=exst, opened=opn)

if (exst) then

  if (.not. opn) then
    open (1, file=file_avg)

    read (1, '(a)') buff

    if (buff(1:1) == '#') then
      read (buff(2:), *) n_avg
    else
      write (*, *) 'avg_stats: error'
      write (*, *) trim (file_avg), ' does not have expected format on line 1'
      stop  !--need to replace all stops with nice mpi exits
    end if
  end if

  !--skip data header lines here
  do
    read (1, '(a)') buff
    if (trim (buff) == trim (end_hdr_avg)) exit
  end do

  do j = 1, size (a_avg, 2)
    read (1, *) z, a_avg(:, j)  !--z is just placeholder here
  end do

  if (present (leaveopn)) then
    if (leaveopn /= 'yes') close (1)  !--case sensitive here
  else
    close (1)
  end if

else

  n_avg = 0
  a_avg = 0._rprec

end if

end subroutine init_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_avg
real (rprec), intent (in) :: x(:, :)  !--coord columns for x, y, etc
real (rprec), intent (in) :: a_avg(:, :)

character (*), optional, intent (in) :: hdr
character (*), optional, intent (in) :: position

character (64) :: r_fmt, fmt
character (32) :: posn

integer :: j

!---------------------------------------------------------------------

!--check sizes compatible
if (size (x, 2) /= size (a_avg, 2)) then
  write (*, *) 'write_avg: error with sizes of x, a_avg'
  stop
end if

if (present (position)) then
  posn = position
else
  posn = 'rewind'
end if

open (1, file=file_avg, position=posn)

if (trim (posn) /= 'append') then  !--case sensitive
  write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
end if

if (present (hdr)) then
  !--write data header, if present
  write (1, '(a)') trim (hdr)
end if

!--write something to indicate end of header, always do this
write (1, '(a)') end_hdr_avg

!--define output format
write (r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                           '.', precision (1._rprec)
write (fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                         '(1x,', trim (r_fmt), '))'

!--write to file
do j = 1, size (a_avg, 2)
  write (1, fmt) x(:, j), a_avg(:, j)
end do

close (1)

end subroutine write_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt, jx_s, jt_start_write
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

integer, save :: rec = 0
integer :: iolen

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

if (.not. initialized) then

  inquire (unit=lun, opened=opn)
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  inquire (iolength=iolen) u(jx_s, :, :), v(jx_s, :, :), w(jx_s, :, :)
  open (unit=lun, file=inflow_file, access='direct', recl=iolen)

  !--rec should be 0 at this point
  initialized = .true.

end if

if (jt == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt >= jt_start_write) then
  rec = rec + 1
  write (unit=lun, rec=rec) u(jx_s, :, :), v(jx_s, :, :), w(jx_s, :, :)
end if

end subroutine inflow_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_read ()
use param, only : ny, nz, pi, nsteps, jt
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .true.

logical, parameter :: DEBUG = .false.
logical, save :: init_DEBUG = .false.
integer :: jy, jz

integer :: i
integer :: iolen
integer, save :: rec = 0
integer, save :: nrec

logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt
$if ($MPI)
  real (rprec) :: u_tmp(ny, 0:nz), v_tmp(ny, 0:nz), w_tmp(ny, 1:nz+1)
$else
  real (rprec) :: u_tmp(ny, nz), v_tmp(ny, nz), w_tmp(ny, nz)
$endif

!---------------------------------------------------------------------

if (.not. initialized) then

  inquire (file=inflow_file, exist=exst)
  if (.not. exst) then
    write (*, *) sub // ': ', inflow_file, ' does not exist'
    stop
  end if

  inquire (unit=lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    !--figure out the record length
    inquire (iolength=iolen) u(jx_s, :, :), v(jx_s, :, :), w(jx_s, :, :)

    !--figure out the number of records
    call len_da_file (inflow_file, iolen, nrec)

    write (*, *) sub // ': number of records = ', nrec

    !--check minimum length requirement
    if (2 * (l_blend + 2) >= nrec) then
      write (*, *) sub // ': ', inflow_file, 'is too short to recycle'
      stop
    end if

    open (unit=lun, file=inflow_file, access='direct', recl=iolen)
  end if

  initialized = .true.

end if

rec = rec + 1

read (unit=lun, rec=rec) u(jx_s, :, :), v(jx_s, :, :), w(jx_s, :, :)

if ((recycle) .and. (nrec - rec <= l_blend)) then

  i = nrec - rec
  wgt = 0.5_rprec * ( 1._rprec + cos (pi * real (i, rprec) / l_blend) )

  read (unit=lun, rec=l_blend-i+1) u_tmp, v_tmp, w_tmp

  u(jx_s, :, :) = wgt * u(jx_s, :, :) + (1._rprec - wgt) * u_tmp
  v(jx_s, :, :) = wgt * v(jx_s, :, :) + (1._rprec - wgt) * v_tmp
  w(jx_s, :, :) = wgt * w(jx_s, :, :) + (1._rprec - wgt) * w_tmp

  if (rec == nrec) rec = l_blend + 1  !--reset the record counter
    !--the first read above just read the last record;
    !  the second read above just read record l_blend+1
    !  so, on the next pass, the read below should read record l_blend + 2

end if

if (DEBUG) then  !--write out slices as an ascii time series

  if (.not. init_DEBUG) then

    inquire (unit=88, exist=exst, opened=opn)

    if (exst .and. .not. opn) then
      open (unit=88, file='inflow_read_debug.dat')
      write (88, '(a)') 'variables = "y" "z" "t" "u" "v" "w"'
      write (88, '(3(a,i0))') 'zone, f=point, i= ', ny, ', j= ', nz,  &
                              ', k= ', nsteps
    else
      write (*, *) sub // ': problem opening debug file'
      stop
    end if

    init_DEBUG = .true.

  end if

  do jz = 1, nz
    do jy = 1, ny
      write (88, '(3(1x,i0),3(1x,es12.5))') jy, jz, jt, u(jx_s,jy,jz),  &
                                            v(jx_s,jy,jz), w(jx_s,jy,jz)
    end do
  end do
end if

end subroutine inflow_read


subroutine create_folder_for_output

use param
implicit none
  
!  CHARACTER(len=255)::cwd,makedirectory,folder,conc,vel,Re_data,Rhs_data,output1,&
!                      brk_freq,all_files,all_files2
  integer :: time_tot
  time_tot = nsteps*dt_dim
  CALL getcwd(cwd)
!  write(folder,'(A,A,A,I2.2,A,I2,A)')TRim(cwd),'/','output_',time_tot,'_',npcon-1,'drops'
!  write(*,*) Trim(folder)
!  write(makedirectory,'(A,A,A,A)') 'mkdir'
!  write(folder,'(A,A,A,I2.2,A,I2.2,A)')TRim(cwd),'/','output_',time_tot,'_',npcon-1,'drops' 
!write(folder,'(A,A,A,I2.2,A,I4.4,A)')TRim(cwd),'/','output_3dsmooth',time_tot,'_',int(C1),'forcing'
!write(vel,'(A,A,A,I2.2,A)')TRIM(folder),'/','vel_',time_tot,'n/'
!write(conc,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','conc_',time_tot,'_',npcon-1,'drops/'
!write(Re_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Re_data_',time_tot,'_',npcon-1,'drops/'
!write(brk_freq,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','brk_freq_',time_tot,'_',npcon-1,'drops/'
!write(Rhs_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Rhs_',time_tot,'_',npcon-1,'drops/' 

!  write(makedirectory,'(A,A,A,A)') 'mkdir'
 ! write(folder,'(A,A,A,I2.2,A,I4.4,A)')TRim(cwd),'/','output_3dsmooth',time_tot,'_',int(C1),'forcing'
 ! write(output1,'(A,X,A)') Trim(makedirectory),Trim(folder)
!  write(vel,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','vel_',time_tot,'_',npcon-1,'drops/'
!  write(conc,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','conc_',time_tot,'_',npcon-1,'drops/'
!  write(Re_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Re_data_',time_tot,'_',npcon-1,'drops/'
!  write(brk_freq,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','brk_freq_',time_tot,'_',npcon-1,'drops/'
!  write(Rhs_data,'(A,A,A,I2.2,A,I2.2,A)')TRIM(folder),'/','Rhs_',time_tot,'_',npcon-1,'drops/'

 write(makedirectory,'(A,A,A,A)') 'mkdir'
!  write(folder,'(A,A,A,I2.2,A,I2.2,A)')TRim(cwd),'/','output_',time_tot,'_',npcon-1,'drops' 
  write(folder,'(A,I3.3,A,I2.2,A)')'output_',time_tot,'_',npcon-1,'_c600'
  write(output1,'(A,X,A)') Trim(makedirectory),Trim(folder)
  write(vel,'(A,A,A,I3.3,A)')TRIM(folder),'/','vel_',time_tot,'/'
  write(conc,'(A,A,A,I3.3,A)')TRIM(folder),'/','conc_',time_tot,'/'
  write(Re_data,'(A,A,A,I3.3,A)')TRIM(folder),'/','Re_data_',time_tot,'/'
  write(brk_freq,'(A,A,A,I3.3,A)')TRIM(folder),'/','brk_freq_',time_tot,'/'
  write(Rhs_data,'(A,A,A,I3.3,A)')TRIM(folder),'/','Rhs_',time_tot,'/'


!  write(vel,'(A,A,A,I2.2,A,I2,A)') TRIM(folder),'/','vel_',time_tot,'_',npcon-1,'drops'
!  write(conc,'(A,A,A,I2.2,A,I2,A)') TRIM(folder),'/','conc_',time_tot,'_',npcon-1,'drops'
!  write(Re_data,'(A,A,A,I1,A,I2,A)') TRIM(folder),'/', 'Re_data_',time_tot,'_',npcon-1,'drops'
!!  write(brk_freq,'(A,A,A,I1,A,I2,A)') TRIM(folder),'/','brk_freq_',time_tot,'_',npcon-1,'drops'
!  write(Rhs_data,'(A,A,A,I1,A,I2,A)') TRIM(folder),'/','Rhs_',time_tot,'_',npcon-1,'drops'
  write(all_files,'(A,X,A,X,A,X,A,X,A)') Trim(makedirectory),Trim(vel),Trim(conc),trim(brk_freq)
  write(all_files2,'(A,X,A,X,A,X,A,X,A)') Trim(makedirectory),Trim(Re_data),Trim(Rhs_data)
  call system(output1)
  call system(all_files)
  call system(all_files2)

end subroutine create_folder_for_output



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine len_da_file(fname, lenrec, length)
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
  return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file

subroutine log_file_create

use param
implicit none

open(120,file="log.dat")
write(120,*) "System Parameters used in this simulation"
write(120,'(A,2X,I6.6,2X,3I3.3)') "number of timesteps,nx,ny,nz = ",nsteps,nx,ny,nz
write(120,'(A,4F6.3)') "dt,dx,dy,dz", dt_dim,dx,dy,dz
write(120,'(A,4E15.7)') "Surface tension , Oh_ratio, dmax , kb", sigma,mu_ratio,dmax,k_b
write(120,*) "Random noise in both x and y direction"
close(120)

end subroutine log_file_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module io
