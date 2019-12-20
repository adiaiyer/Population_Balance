module io
use types,only:rprec
use param, only : ld, nx, ny, nz, jx_s, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc,      &
                  average_dim_num, nz_tot,jt_total,p_count,dt,z_i,u_star
implicit none
private
public :: openfiles,output_loop,output_final,                   &
          mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,         &
          inflow_read, inflow_write, avg_stats
public :: average_dim_select_flag, dim1_size, dim2_size,        &
          dim1_global, dim2_global, collocate_MPI_averages,     &
          compute_avg_var,post_spec

!integer, parameter :: nz_tot = (nz - 1) * nproc + 1

integer,parameter::num_hour_out=1
!integer,parameter::base=nint(num_hour_out*3600/(dt*z_i/u_star)),nwrite=base
!integer,parameter::base=10,nwrite=base
integer,parameter::base=250,nwrite=base

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
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process 
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical,parameter:: inst_slice_flag=.false.
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
use param,only:output,dt,c_count,S_FLAG,SCAL_init,jt,jan_diurnal_run,PCon_FLAG,PCon_init
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


if ((use_avgslice) .and. (mod(jt,c_count)==0)) then
         call avgslice()
       if ((S_FLAG) .and. (jt.GE.SCAL_init)) then
         call scalar_slice() ! Uses file unit numbers (36-47)
       endif
       if ((PCon_FLAG) .and. (jt.GE.PCon_init)) then
         call pollen_slice()
       endif
end if

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

!DY if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!DY     open(1,file=fname,form='unformatted')
!DY end if

!DY     call checkpoint (1)

!DY if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!DY     close(1)
!DY end if

!    print*, "io.f90 checkpoint2"
    open(1000+coord,file=fname,form='unformatted')
!    print*, "io.f90 checkpoint3"
    call checkpoint_v2 (1000+coord)
!    print*, "io.f90 checkpoint4"
    close(1000+coord)

!DY End here
!****************************

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

!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot

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

  if (S_FLAG) then 
  call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
                    theta_tot(1,1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif
  if (PCon_FLAG) then
     do ipcon=1,npcon
        call mpi_gatherv (PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
             PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  endif
$endif

if((.not. USE_MPI)) then

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
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
end if

else if ((USE_MPI .and. coord == 0)) then

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
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
505 format(12e12.4)
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
  open(11000000+jt_total)
!!$  write(11000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(11000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta'
  write(11000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',1,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=ny/2,ny/2
        do jx=1,nx
           write(11000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), &
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(11000000+jt_total)
  open(12000000+jt_total)
!!$  write(12000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(12000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta'
  write(12000000+jt_total,*) 'zone t="',jt_total,'" i=',1,' j=',ny,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=1,ny
        do jx=nx/2,nx/2
           write(12000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), &
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(12000000+jt_total)
  open(13000000+jt_total)
!!$  write(13000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
  write(13000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta'
  write(13000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',1,' f=point'
  do jz=1,1
     do jy=1,ny
        do jx=1,nx
           write(13000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), &
                u_tot(jx,jy,jz),v_tot(jx,jy,jz),w_tot(jx,jy,jz), &
                theta_tot(jx,jy,jz)!,Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(13000000+jt_total)
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
  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2'
  write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' j=',ny,' k=',nz_tot,' f=point'
  do jz=1,nz_tot
     do jy=1,ny
        do jx=1,nx
           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jy-1)*L_y*z_i/ny,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2)!, &
!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
  enddo
  close(1000000+jt_total)
!DY End here
else ! No SCALARS
  write (lun) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot)
end if

endif

end subroutine checkpoint



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang
!DY The original version writes data into a single file, which requires large memory for high resolution
!DY This version writes data on each cpu into separate files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint_v2 (lun)
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

$if ($MPI)
  integer :: sendcounts(nproc)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif

real (rprec), dimension(ld, 1:nz) :: u1,v1,w1,theta1
real (rprec), dimension(ld, 1:nz, npcon) :: PCon1
real (rprec), dimension(ld, 1:nz_tot) :: u_tot
real (rprec), dimension(ld, 1:nz_tot) :: v_tot
real (rprec), dimension(ld, 1:nz_tot) :: w_tot
real (rprec), dimension(ld, 1:nz_tot) :: theta_tot
real (rprec), dimension(ld, 1:nz_tot, npcon) :: PCon_tot
!real (rprec), dimension(ld, ny, 1:nz_tot) :: Cs_opt2_tot
!real (rprec), dimension(ld, ny, 1:nz_tot) :: Nu_t_tot

!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_LM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_MM_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_QN_tot
!!$real(kind=rprec),dimension(ld,ny,nz_tot) :: F_NN_tot

!---------------------------------------------------------------------
$if ($MPI)
  u1 = u(:,yps,1:nz)
  v1 = v(:,yps,1:nz)
  w1 = w(:,yps,1:nz)
  theta1 = theta(:,yps,1:nz)
  PCon1 = PCon(:,yps,1:nz,:)
  sendcounts = size (u1(:,1:nz))
  recvcounts = size (u1(:,1:nz-1))
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (u1(1,1), sendcounts, MPI_RPREC,&
                    u_tot(1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (v1(1,1), sendcounts, MPI_RPREC,&
                    v_tot(1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  call mpi_gatherv (w1(1,1), sendcounts, MPI_RPREC,&
                    w_tot(1,1), sendcounts, displs,       &
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

  if (S_FLAG) then 
  call mpi_gatherv (theta1(1,1), sendcounts, MPI_RPREC,&
                    theta_tot(1,1), sendcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
  endif
  if (PCon_FLAG) then
     do ipcon=1,npcon
        call mpi_gatherv (PCon1(1,1,ipcon), sendcounts, MPI_RPREC,&
             PCon_tot(1,1,ipcon), sendcounts, displs,       &
             MPI_RPREC, rank_of_coord(0), comm, ierr)
     enddo
  endif
$endif

if((.not. USE_MPI)) then

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
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
end if

else if (USE_MPI) then

if (S_FLAG .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
              theta(:,:,1:nz), PCon(:,:,1:nz,1:npcon), &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)

505 format(12e12.4)

  if(coord.eq.0) then
     open(11000000+jt_total)
!!$  write(11000000+jt_total,*) 'variables=x,y,z,PCon1,PCon2,u,v,w,theta,Cs_opt2,Nu_t'
     write(11000000+jt_total,*) 'variables=x,z,PCon1,PCon2,u,v,w,theta'
     write(11000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
     if(OCEAN_FLAG) then
        do jz=1,nz_tot
           do jx=1,nx
              write(11000000+jt_total,505) (jx-1)*L_x*z_i/nx,-(jz-1)*z_i/nz_tot, &
                   !                PCon_tot(jx,jz,1),PCon_tot(jx,jz,2),PCon_tot(jx,jz,3), &
                   PCon_tot(jx,jz,1), PCon_tot(jx,jz,2), &
                   u_tot(jx,jz),v_tot(jx,jz),-w_tot(jx,jz), &
                   2._rprec-theta_tot(jx,jz)!,Cs_opt2_tot(jx,jz),Nu_t_tot(jx,jz)
           enddo
        enddo
     else
        do jz=1,nz_tot
           do jx=1,nx
              write(11000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot, &
                   !                PCon_tot(jx,jz,1),PCon_tot(jx,jz,2),PCon_tot(jx,jz,3), &
                   PCon_tot(jx,jz,1), PCon_tot(jx,jz,2), &
                   u_tot(jx,jz),v_tot(jx,jz),w_tot(jx,jz), &
                   theta_tot(jx,jz)!,Cs_opt2_tot(jx,jz),Nu_t_tot(jx,jz)
           enddo
        enddo
     endif
     close(11000000+jt_total)
  endif

!DY End here
elseif (S_FLAG) then !WITH SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
              theta(:,:,1:nz)
!DY Added by Di Yang for testing
  if(coord.eq.0) then
!  write(100,*) 'variables=x,y,z,PCon'
     open(1000000+jt_total)
!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta,Cs_opt2,Nu_t'
     write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,theta'
     write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
     do jz=1,nz_tot
        do jx=1,nx
           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot,&
                u_tot(jx,jz),v_tot(jx,jz),w_tot(jx,jz),theta_tot(jx,jz)!, &
!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
     close(1000000+jt_total)
  endif
!DY End here
elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz), &
              PCon(:,:,1:nz,1:npcon),   &
              deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
              P_surf_flux(:,:), P_surf_flux_dep(:,:)
!DY Added by Di Yang for testing
  if(coord.eq.0) then
!  write(100,*) 'variables=x,y,z,PCon'
     open(1000000+jt_total)
!  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,Pcon2,Pcon3,Cs_opt2,Nu_t'
!!$  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,Cs_opt2,Nu_t'
     write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2'
     write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nx,' k=',nz_tot,' f=point'
     do jz=1,nz_tot
        do jx=1,nx
           write(1000000+jt_total,505) (jx-1)*L_x*z_i/nx,(jz-1)*z_i/nz_tot, &
!                PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
                PCon_tot(jx,jz,1), PCon_tot(jx,jz,2)!, &
!                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
        enddo
     enddo
     close(1000000+jt_total)
  endif
!DY End here
else ! No SCALARS
  write (lun) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz)
end if

endif

end subroutine checkpoint_v2



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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module io
