module param
use types,only:rprec
$if ($MPI)
  use mpi
$endif
implicit none

save

private rprec  !--this is dumb.
public

!--mpi stuff
$if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 32
$else
  $define $MPI_LOGICAL .false.
  $define $NPROC 1
$endif

logical, parameter :: USE_MPI = $MPI_LOGICAL

$undefine $MPI_LOGICAL

$if ($MPI)
  integer :: status(MPI_STATUS_SIZE)
$endif

character (*), parameter :: path = './'
!character (*), parameter :: path = './test_remote/homog/'
!character (*), parameter :: path = './MPI_neutral_test/np_2/'

!--this stuff must defined, even if not using MPI
integer, parameter :: nproc = $NPROC  !--this must be 1 if no MPI
integer :: ierr
integer :: comm
integer :: up, down
integer :: global_rank
integer :: MPI_RPREC, MPI_CPREC
integer :: rank = -1   !--init to bogus (so its defined, even if no MPI)
integer :: coord = -1  !--same here
!DY Start 1
!DY integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
integer, allocatable, dimension(:) :: rank_of_coord, coord_of_rank
!DY End 1
!--end mpi stuff

logical, parameter :: VERBOSE = .false.  !--prints small stuff to screen
                      !--use DEBUG to write lots of data/files

integer,parameter:: nx=150,ny=150,nz=(280-1)/nproc + 1
!integer,parameter:: nx=32,ny=32,nz=(32-1)/nproc + 1
integer,parameter:: nz_tot=(nz - 1) * nproc + 1 
integer,parameter:: nx2=3*nx/2,ny2=3*ny/2
integer,parameter:: lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

real (rprec), parameter :: BOGUS = -1234567890._rprec
real (rprec), parameter :: TINYS = 1.e-30_rprec

real(rprec),parameter::pi=3.1415926535897932384626433_rprec
real(kind=rprec),parameter::L_x=2.*0.4222222_rprec,L_y=2.*0.4222222_rprec
!real(kind=rprec),parameter::L_x=10._rprec,L_y=5._rprec
real(rprec),parameter::z_i=0.9_rprec, L_z=z_i/nproc
                            !--L_z is not nondimensionalized by z_i yet
! set the aspect ratio of the box, already nondimensional
real(rprec),parameter::dz=L_z/z_i/(nz-1)
real(rprec),parameter::dx=L_x/nx,dy=L_y/ny

! dt_dim=0.1  => 1 hour  = 36000 steps
!                3 hours = 108000 steps
integer,parameter::nsteps=5000 !simulation steps
!integer,parameter::nsteps=10 !simulation steps

!--Coriolis stuff; coriol=non-dim coriolis parameter,
! ug=horiz geostrophic vel, vg=transverse geostrophic vel
logical,parameter::coriolis_forcing=.false.
!real(rprec),parameter::ug_dim=3._rprec,vg_dim=-9._rprec 
real(rprec),parameter::ug_dim=8._rprec,vg_dim=0._rprec
!real(rprec),parameter::coriol=1.39E-04*z_i/u_star,      &
!                       ug=u_star/u_star,vg=0._rprec/u_star

! u_star=0.45 if coriolis_forcing=.FALSE. and =ug if coriolis_forcing=.TRUE.
real(rprec),parameter::u_star=0.00006125_rprec,Pr=0.4_rprec
!real(rprec),parameter::u_star=sqrt(ug_dim**2+vg_dim**2),Pr=.4_rprec

real(rprec),parameter::dt_dim=0.001_rprec !dimensional time step in seconds
real(rprec),parameter::dt=dt_dim*u_star/z_i


!real(rprec),parameter::coriol=9.125E-05*z_i/u_star! for latitude = 38.74N
!real(rprec),parameter::coriol=1.3942923E-04*z_i/u_star !for latitude = 37.6N
real(rprec),parameter::coriol=7.0E-05*z_i/u_star !for latitude = 28.73N
real(rprec),parameter::ug=ug_dim/u_star,vg=vg_dim/u_star

real(rprec),parameter::vonk=.4_rprec

integer,parameter::c_count=5,p_count=1000 !p_count=6000 => 10 min.
integer, parameter :: cs_count = 1  !--tsteps between dynamic Cs updates
logical,parameter::output=.true.
logical, parameter :: use_avgslice = .true.
! The parameter average_dim_num controls the number of dimensions averaged
! while outputting files from avgslice .. 
! Possible values of average_dim_num: 2 => average over x,y,t ; 1 => average over y,t
integer, parameter :: average_dim_num = 1

! Output short binary for postprocess and video (includes velocity field, pollen concentration and deposition)
logical, parameter :: output_video=.TRUE.
integer, parameter :: video_start=6000
integer, parameter :: video_end=600000
integer, parameter :: video_freq=1000

! nu_molec is dimensional m^2/s
real(rprec),parameter::nu_molec=1.14e-5_rprec

logical,parameter::use_bldg=.false.
logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.

!Model type: 1->Smagorinsky; 2->Dynamic; 3->Scale dependent
!            4->Lagrangian scale-sim   5-> Lagragian scale-dep
!Models type: 1->static prandtl, 2->Dynamic
!Cs is the Smagorinsky Constant
!Co and nnn are used in the mason model for smagorisky coeff
integer,parameter::model=5,models=1,nnn=2,BETA_FLAG=1

! calc_lag_Ssim enables calculation of a Cs based on lag Ssim
! using BETA = 1 while using lag sdep for the real Cs. The Cs calc using
! calc_lag_ssim is only used for plotting purposes
! CAUTION : use only with model=5 and only when required (e.g. diurnal runs)
!logical,parameter::calc_lag_Ssim=.true.
logical,parameter::calc_lag_Ssim=.false.

real(kind=rprec),parameter::Co=0.16_rprec
real(kind=rprec),parameter::cs=0.2_rprec

!Test filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::ifilter=1

! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
! damping method = 1: use Nieuwstadt's method, = 2: use Raleigh damping method
!integer,parameter::ubc=1,damping_method=1
integer,parameter::ubc=0,damping_method=1

character (*), parameter :: lbc_mom = 'wall'
                            !--'wall', 'stress free'
                            
! jx_s:bufferzone starting point,x_relax:buffer zone length
logical,parameter::inflow=.false.
integer,parameter::jx_s=nx,x_relax=8
real(kind=rprec),parameter::face_avg=1._rprec
logical, parameter :: read_inflow_file = .false.

logical, parameter :: write_inflow_file = .false.
                      !--records at position jx_s
integer, parameter :: jt_start_write = 170000

! forcing along top and bottom bdrys
! if inflow is true and force_top_bot is true, then the top & bottom
! velocities are forced to the inflow velocity
logical, parameter :: force_top_bot = .false.

logical, parameter :: use_mean_p_force = .false.
real (rprec), parameter :: mean_p_force = 1._rprec * z_i/L_z/nproc
                                          !--usually just z_i/L_z

! Constant angle for the pressure forcing (no angle if use_force_angle=.FALSE. )
logical, parameter :: use_force_angle = .FALSE. 
real (rprec), parameter :: force_angle = -pi/6._rprec

logical, parameter :: use_trees = .false.

! Resolved canopy
logical, parameter :: use_res_canopy = .FALSE.
logical, parameter :: use_field_scale = .FALSE.
logical, parameter :: use_plant_scale = .FALSE.

! canopy height
real (rprec), parameter :: h_canopy = 2.67/z_i

! drag coefficient
real (rprec), parameter :: Cd = 0.17

! leaf area index
real (rprec), parameter :: LAIp = 6.0
real (rprec), parameter :: LAIw = 2.9

! leaf area density
real (rprec), parameter :: a_leaf1 = 0.99
real (rprec), parameter :: a_leaf2 = 2.17
real (rprec), parameter :: a_leaf3 = 3.93
real (rprec), parameter :: a_leaf4 = 4.32
real (rprec), parameter :: a_leaf5 = 4.68
real (rprec), parameter :: a_leaf6 = 4.35
real (rprec), parameter :: a_leaf7 = 0.79

! height corresponding to leaf area density
real (rprec), parameter :: h_leaf1 = 0.33/z_i
real (rprec), parameter :: h_leaf2 = 0.66/z_i
real (rprec), parameter :: h_leaf3 = 1.0/z_i
real (rprec), parameter :: h_leaf4 = 1.33/z_i
real (rprec), parameter :: h_leaf5 = 1.66/z_i
real (rprec), parameter :: h_leaf6 = 2.0/z_i
real (rprec), parameter :: h_leaf7 = 2.33/z_i

integer :: jt,jt_total  ! global time-step counter

! time advance parameters (AB2)
real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1

!------xxxxxxxxx--SCALARS_PARAMETERS--xxxxxxxxx---------------
! S_FLAG=1 for Theta and q, =0 for no scalars
!logical,parameter::S_FLAG=.TRUE.,coupling_flag=.FALSE.,mo_flag=.TRUE.
logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.FALSE.,WALL_HOMOG_FLAG=.FALSE.
!logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.TRUE.,WALL_HOMOG_FLAG=.TRUE.
!integer,parameter::DYN_init=1800, SCAL_init=4001
!integer,parameter::DYN_init=1500, SCAL_init=3000
!integer,parameter::DYN_init=3000, SCAL_init=3000
integer,parameter::DYN_init=1000, SCAL_init=1
!integer,parameter::DYN_init=1, SCAL_init=1

! lbc=0: prescribed surface temperature, lbc=1 prescribed surface flux 
integer,parameter :: lbc=1, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0
logical,parameter :: remote_to_patch_flag=.FALSE.! create a 2 patch T_s field using the remote-sensed data
! The corresponding remote_to_patch subroutine is in bottombc.f90
integer,parameter :: diurnal_forcing_flag=0, no_days=1
logical,parameter :: jan_diurnal_run=.false.,ug_diurnal_variation=.false.
logical,parameter :: GABLS_diurnal_test=.FALSE.
!integer,parameter::patch_flag=1, remote_flag=0, time_start=0
! initu=.TRUE. & initsc=.FALSE read velocity fields from a binary file
! initu=.TRUE. & initsc=.TRUE. read velocity & scalar fields from a binary file
! initu=.FALSE. & S_FLAG=.TRUE. initialize velocity & scalar fields using ran
! initu=.FALSE. & S_FLAG=.FALSE. initialize velocity fields using ran
!--initu = true to read from a file; false to create with random noise
!logical, parameter :: initu = .false.
!--initlag = true to initialize cs, FLM & FMM; false to read from vel.out
!logical, parameter :: inilag = .false.
logical,parameter :: initu=.false.,initsc=.false.,inilag=.true.,interp=.TRUE.
!logical,parameter :: initu=.FALSE.,initsc=.FALSE.,inilag=.TRUE.,interp=.TRUE.
!logical,parameter :: initu=.TRUE.,initsc=.FALSE.,inilag=.TRUE.,interp=.FALSE.
! Added a new parameter - passive_scalar for passive scalars with bldngs
logical,parameter :: passive_scalar=.false.,GABLS_test=.false.
logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
logical,parameter :: check_dt=.TRUE.
integer,parameter :: stencil_pts=4
logical,parameter :: coarse_grain_flag=.FALSE.
!inversion strength (K/m)
!real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0.010_rprec,dTdz_top=0._rprec
real(kind=rprec),parameter::g=9.81_rprec, inv_strength=255._rprec,dTdz_top=inv_strength
real(kind=rprec),parameter :: q_s1=11._rprec,q_s2=13._rprec,q_mix=12._rprec,q_top=12._rprec,wq_s=1.77e-10_rprec
!real(kind=rprec),parameter :: theta_s1=310.15_rprec,theta_s2=300.15_rprec
real(kind=rprec),parameter :: theta_top=300._rprec,T_scale=300._rprec,wt_s=0.0_rprec,T_init=283.15_rprec
!real(kind=rprec),parameter :: theta_top=300._rprec,T_scale=300._rprec,wt_s=0.2_rprec,T_init=300._rprec
!real(kind=rprec),parameter :: theta_top=300._rprec,T_scale=300._rprec,wt_s=0.2_rprec,T_init=300._rprec

real(kind=rprec),parameter::cap_thick=80._rprec, z_decay=1._rprec




!------xxxxxxxxx--POLLEN_PARAMETERS--xxxxxxxxx---------------

!  note always set 
!  lbc=1

! Flag to include pollen concentration in the simulation
! PCon_FLAG=.TRUE.   - include pollen
! PCon_FLAG=.FALSE.  - do not include pollen
logical,parameter :: PCon_FLAG=.TRUE.

! Flag to initialize pollen concentration
! initPCon=.TRUE.  - read pollen concentration from vel_sc.out
! initPCon=.FALSE. - initialize pollen concentration with zero
logical,parameter :: initPCon=.TRUE.

! Scale for pollen concentration (in grains/m3 or g/m3 or kg/m3)
real(kind=rprec),parameter :: PCon_scale=1._rprec

! Numerical scheme for pollen concentration
! nxc,nyc,nzc are the number of points to be allocated in the pollen arrays
! 1 - pseudo-spectral (same as velocity field)
! 2 - QUICK finite volumes
! 3 - SMART finite volumes
integer,parameter :: PCon_scheme=3

! Flag for periodic or inflow/outflow BC (not valid for PCon_scheme=1)
! periodicbcx=.TRUE.   - periodic BC in x
! periodicbcx=.FALSE.  - inflow BC at x=0 and outflow at x=Lx
logical,parameter :: periodicbcx=.FALSE.
! periodicbcy=.TRUE.   - periodic BC in y
! periodicbcy=.FALSE.  - inflow BC at y=Ly and outflow at y=0
logical,parameter :: periodicbcy=.FALSE.


! Time step to start evolving pollen concentration field
integer,parameter :: PCon_init=1

! Flag to specify bottom bounday condition
! lbcp=0  - prescribed surface pollen concentration
! lbcp=1  - prescribed pollen surface flux
integer,parameter :: lbcp=1
! Prescribed surface pollen concentration
real(kind=rprec),parameter :: PCon_sfc=500._rprec/PCon_scale
! Prescribed pollen surface flux
real(kind=rprec),parameter :: sfc_flux=0._rprec/(PCon_scale*u_star)


! Schmidt number for pollen concentration
! SGS Schmidt number
real(kind=rprec),parameter :: Sc=0.4_rprec
! Turbulent Schimdt number for boundary condition
real(kind=rprec),parameter :: Sc_boun=0.95_rprec

! Scalar Eddy diffusivity model:
! 1-> Constant Schmidt number; 
! 2-> Plane Averaged Scale Invariant (PASI);
! 3-> Plane Averaged Scale Dependent (PASD);
! 4-> Lagrangian Averaged Scale Invariant (LASI);
! 5-> Lagragian Averaged Scale Dependent (LASD);
integer,parameter :: model_sc=1

! Lagrangian initialization
! FALSE - Do NOT read initialization for the Lagrangian variables
! TRUE  - Read initialization for the Lagrangian variables
! Note: If inilag=.TRUE. or initsc=.FALSE. then inilag_sc is set to false by default
logical,parameter :: inilag_sc=.FALSE.


! Reference height were concentrations are specified (similar to z_o for momentum)
real(kind=rprec),parameter :: z_ref=100._rprec

! Flag to include gravitational settling
! settling=.TRUE.   - include settling
! settling=.FALSE.  - do not include settling
logical,parameter :: settling=.TRUE.

! Number of droplet/bubble sizes
integer,parameter :: npcon=2

! Settling velocity to be used
! Here is the base value for a set of particle sizes
! settling_vel=settling_vel0*2._rprec**(ipcon-1)
real(kind=rprec),parameter :: settling_vel0=0.06_rprec/u_star
!real(kind=rprec),parameter :: settling_vel0=0.0027_rprec/u_star
real(kind=rprec),parameter,dimension(npcon) :: settling_vel=(/0.06_rprec/u_star,0.0_rprec/u_star/)
!real(kind=rprec),parameter :: settling_vel(1)=0.1722_rprec/u_star
!real(kind=rprec),parameter :: settling_vel(2)=0.0112_rprec/u_star


! Flag to include particle inertia (uses the Eulerian equilibrium of Shotorban and Balachandar (2007))
! PCon_acc=.TRUE.  - include inertial
! PCon_acc=.FALSE. - do not include inertial
logical,parameter :: PCon_acc=.TRUE.

! Flag to include buffer zone for zero concentration inflow
! inflow_PCon=.TRUE.   - zero concentration inflow (do not modify velocity field)
! inflow_PCon=.FALSE.  - periodic concentration inflow
logical,parameter::inflow_PCon=.FALSE.
! Point where zero concentration ir reached and size of buffer
integer,parameter::jx_p=nx,x_relaxp=8



! For point source within the domain
logical,parameter :: pointsource=.TRUE.
! For single point source
logical,parameter :: single_point=.TRUE.
! x, y and z location of the source (node #)
integer,parameter :: xps=nx/2
integer,parameter :: yps=ny/2
!integer,parameter :: zps=70
integer,parameter :: zps=nz_tot-26
! Source strength in grains/m3 per second
!real(kind=rprec),parameter :: Q_src=(0.11875_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
!real(kind=rprec),parameter :: Q_src=(1.0_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
!real(kind=rprec),parameter:: Q_src=(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
!real(kind=rprec),parameter :: Q_src=(0.05_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
!     (2.28e9_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))/)
real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(2.1e-6_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
     (6.45e-9_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))/)
! timesteps to begin and end point source
integer,parameter :: ini_src=1
integer,parameter :: end_src=1200000


! For the field experiment rag06
logical,parameter :: rag06=.FALSE.
logical,parameter :: fieldsize=.FALSE.
logical,parameter :: fieldsize_stripe=.TRUE.

! Starting and ending grid points of the field
! for fieldsize_stripe=.TRUE., only field_xs and field_xe are used
integer,parameter :: field_xs=9  ! Starting x grid point
integer,parameter :: field_xe=11 ! Ending x grid point
integer,parameter :: field_ys=49  ! Starting y grid point
integer,parameter :: field_ye=51 ! Ending y grid point

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for buoyancy force due to oil droplet concentration
logical,parameter :: active_pcon=.true.
!DY active_pcon=.true.: add droplet buoyancy term to vertical momentum equation
!DY active_pcon=.false.: passive droplet, no feedback to flow

!real (rprec), parameter :: densratio_pcon=859.8708_rprec/1031._rprec
!real (rprec), parameter :: densratio_pcon=99.4598_rprec/1030._rprec
real (kind=rprec), parameter,dimension(npcon) :: densratio_pcon=(/1.4_rprec/1000._rprec,1000._rprec/1000._rprec/)
!DY densratio is the density ratio between droplet and continuous phase fluid
!DY Here, the value is for oil/water or gas/water

real (rprec), parameter :: V_pcon0=1._rprec/1000._rprec
!real (rprec), parameter :: V_pcon0=4.189e-9_rprec/z_i**3
!real (kind=rprec), parameter,dimension(npcon) :: V_pcon=(/4.189e-9_rprec,6.545e-11_rprec/)
real (kind=rprec), parameter,dimension(npcon) :: V_pcon=(/1._rprec/1.4_rprec,1._rprec/1000._rprec/)

!DY V_pcon is the volume of particles per PCon unit
!DY If PCon is for gram/m^3, V_pcon is 1/rho_p
!DY If PCon is for droplets/m^3, V_pcon is the volume per droplet

!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for sgs droplet acceleration
logical,parameter :: pconsgs_acc_flag=.false.
!DY pconsgs_acc_flag=.true.: add droplet acceleration to PCon sgs flux
!DY pconsgs_acc_flag=.false: no sgs acceleration in PCon sgs flux

integer,parameter:: model_psgsacc=1
!DY model_psgsacc=1: based on acceleration
!DY model_psgsacc=2: based on velocity gradient tensor
!DY model_psgsacc=1: based on pressure Hessian

real,parameter:: Cs_pa=1.0e4_rprec
!DY Cs_pa: prescribed model coefficient for droplet sgs flux due to acceleration

!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for ocean simulation
logical,parameter :: OCEAN_FLAG=.true.
!DY For ocean simulation:
!DY S_FLAG=.true.
!DY coriolis_forcing=.false. This turns off ug,vg. But coriolis is included in main.f90
!DY use_mean_p_force = .false.  Surface stress will be balanced by coriolis force due to Ekman Spiral.

real (rprec), parameter :: nu_ek=1.16e-2_rprec
!DY nu_ek is the eddy viscosity for the Ekman layer initial condition

real (rprec), parameter :: nu_sgs=1.0e-3_rprec
!DY nu_sgs is the subgrid-scale eddy viscosity, prescribed for SGS model test

!DY For sea water: density=1030 kg/m^3, viscosity=1.08e-3 Pa*s
!DY For crude oil at Gulf of Mexico: density=950 kg/m^3, rising velocity=(1030-950)*9.8*d^2/(18*1.08e-3)
!DY If d=5e-4 m, rising velocity=0.01 m/s
!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for wave-induced Stokes drift effect
logical,parameter :: STOKES_FLAG=.true.
!DY Only for ocean simulation: i.e. OCEAN_FLAG=.true.

real (rprec), parameter :: amp_w = 0.004_rprec/z_i
real (rprec), parameter :: lambda_w = 0.3_rprec/z_i
real (rprec), parameter :: wavenm_w = 2._rprec * pi / lambda_w
real (rprec), parameter :: omega_w = sqrt(g*z_i/u_star**2*wavenm_w)
real (rprec), parameter :: U_stokes = omega_w * wavenm_w * amp_w**2
!DY amp_w is the wave amplitude
!DY lambda_w is the wavelength
!DY wavenm_w is the wavenumber
!DY omega_w is the angular frequency
!DY U_stokes is the magnitude of the wave-induced Stokes drift

!DY End here
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!DY Added by Di Yang for deep horizon oil plume with momentum flux
!DY Only for ocean simulation: i.e. OCEAN_FLAG=.true.
!DY Consider Gaussian profile here

logical,parameter :: U_PLUME_FLAG=.false.

real (rprec), parameter :: ubm=0.1186_rprec/u_star
!real (rprec), parameter :: ubm=11.86_rprec/u_star
!DY ubm is the centerline velocity for the inner plume at the bottom boundary

real (rprec), parameter :: b0_plume=0.007_rprec/z_i
!DY b0_plume is the initial width of the inner plume at the bottom boundary

logical,parameter :: GAUSSIAN_SOURCE_2d_FLAG=.true.
!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y)
integer,parameter :: N_GAU = 10._rprec
!DY power coefficient for Gaussian distribution
!DY For regular Gaussian: N_GAU = 2
!DY For super Gaussian: N_GAU has large velue, e.g. N_GAU = 10


logical,parameter :: GAUSSIAN_SOURCE_3d_FLAG=.false.
!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y) and h0_plume in (z)
real (rprec), parameter :: h0_plume=10._rprec/z_i
!DY ho_plume: standard diviation of vertical Gaussian

logical,parameter :: SEABED_FLAG=.true.
!DY if SEABED_FLAG=.true.: no-slip bc for seabed
!DY if SEABED_FLAG=.false.: free-slip bc for seabed

integer, parameter :: iboussinesq = 3
!DY equation of state: rho = rho_0 * ( 1 - alpha_t * theta ) 
!DY iboussinesq = 1: beta = - g * rho / rho_0
!DY iboussinesq = 2: beta = g * ( rho_bar - rho ) / rho_bar
!DY iboussinesq = 3: beta = g * ( rho_0 - rho ) / rho_0 
!DY These options for Boussinesq approximation only work with OCEAN_FLAG

logical,parameter :: PCON_INLET_X=.true., PCON_INLET_Y=.true.

!DY BUBBLE_BURST_FLAG  - out flux surface condition for bubble burst at sea surface
logical,parameter :: BUBBLE_BURST_FLAG = .true.
integer,parameter :: ip_bubble_min=1, ip_bubble_max=npcon
!DY PCon from ip_bubble_min to ip_bubble_max are for bubbles with bursting surface condition
!DY PCon out of this range still use

!DY For delayed dye release
logical,parameter :: DYE_RELEASE_DELAY = .false.
integer,parameter :: ip_dye = 2 ! PCon with ipcon=ip_dye represents the dye
integer,parameter :: ini_dye = 100000 
! When restartting the run with jt_total<=ini_dye, erase existing dye concentration

!DY End here
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module param
