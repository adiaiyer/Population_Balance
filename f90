Binary file lesgo-mpi matches
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(0.11875_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(1.0_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter:: Q_src=(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(0.05_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
param_crossflow_Coriol3D.f90:real(kind=rprec),parameter,dimension(npcon) :: Q_src=(/(0._rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)),&
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(2.1e-6_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(6.448e-2_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
param_crossflow_Coriol3D.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y)
param_crossflow_Coriol3D.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y) and h0_plume in (z)
param.f90:!real(kind=rprec),parameter :: Q_src=(0.11875_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param.f90:!real(kind=rprec),parameter :: Q_src=(1.0_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param.f90:!real(kind=rprec),parameter:: Q_src=(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param.f90:!real(kind=rprec),parameter :: Q_src=(0.05_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
param.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
param.f90:real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(6.448_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
param.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y)
param.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y) and h0_plume in (z)
scalars_module_backup.f90:      released=((jt_total-ini_src+1)*dt*Q_src*dx*dy*dz)
scalars_module_backup.f90:      released=((end_src-ini_src)*dt*Q_src*dx*dy*dz)
scalars_module_backup.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module_backup.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src/4
scalars_module_backup.f90:!    RHS(xps-1,yps,zps)=RHS(xps-1,yps,zps)+Q_src/8
scalars_module_backup.f90:!    RHS(xps+1,yps,zps)=RHS(xps+1,yps,zps)+Q_src/8
scalars_module_backup.f90:!    RHS(xps,yps-1,zps)=RHS(xps,yps-1,zps)+Q_src/8
scalars_module_backup.f90:!    RHS(xps,yps+1,zps)=RHS(xps,yps+1,zps)+Q_src/8
scalars_module_backup.f90:!    RHS(xps-1,yps-1,zps)=RHS(xps-1,yps-1,zps)+Q_src/16
scalars_module_backup.f90:!    RHS(xps-1,yps+1,zps)=RHS(xps-1,yps+1,zps)+Q_src/16
scalars_module_backup.f90:!    RHS(xps+1,yps-1,zps)=RHS(xps+1,yps-1,zps)+Q_src/16
scalars_module_backup.f90:!    RHS(xps+1,yps+1,zps)=RHS(xps+1,yps+1,zps)+Q_src/16
scalars_module_backup.f90:    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src/4/2
scalars_module_backup.f90:    RHS(xps-1,yps,zps)=RHS(xps-1,yps,zps)+Q_src/8/2
scalars_module_backup.f90:    RHS(xps+1,yps,zps)=RHS(xps+1,yps,zps)+Q_src/8/2
scalars_module_backup.f90:    RHS(xps,yps-1,zps)=RHS(xps,yps-1,zps)+Q_src/8/2
scalars_module_backup.f90:    RHS(xps,yps+1,zps)=RHS(xps,yps+1,zps)+Q_src/8/2
scalars_module_backup.f90:    RHS(xps-1,yps-1,zps)=RHS(xps-1,yps-1,zps)+Q_src/16/2
scalars_module_backup.f90:    RHS(xps-1,yps+1,zps)=RHS(xps-1,yps+1,zps)+Q_src/16/2
scalars_module_backup.f90:    RHS(xps+1,yps-1,zps)=RHS(xps+1,yps-1,zps)+Q_src/16/2
scalars_module_backup.f90:    RHS(xps+1,yps+1,zps)=RHS(xps+1,yps+1,zps)+Q_src/16/2
scalars_module_backup.f90:    RHS(xps,yps,zps-1)=RHS(xps,yps,zps-1)+Q_src/4/4
scalars_module_backup.f90:    RHS(xps-1,yps,zps-1)=RHS(xps-1,yps,zps-1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps+1,yps,zps-1)=RHS(xps+1,yps,zps-1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps,yps-1,zps-1)=RHS(xps,yps-1,zps-1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps,yps+1,zps-1)=RHS(xps,yps+1,zps-1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps-1,yps-1,zps-1)=RHS(xps-1,yps-1,zps-1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps-1,yps+1,zps-1)=RHS(xps-1,yps+1,zps-1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps+1,yps-1,zps-1)=RHS(xps+1,yps-1,zps-1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps+1,yps+1,zps-1)=RHS(xps+1,yps+1,zps-1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+Q_src/4/4
scalars_module_backup.f90:    RHS(xps-1,yps,zps+1)=RHS(xps-1,yps,zps+1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps+1,yps,zps+1)=RHS(xps+1,yps,zps+1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps,yps-1,zps+1)=RHS(xps,yps-1,zps+1)+Q_src/8/4
scalars_module_backup.f90:    RHS(xps,yps+1,zps+1)=RHS(xps,yps+1,zps+1)+Q_src/8/4 
scalars_module_backup.f90:    RHS(xps-1,yps-1,zps+1)=RHS(xps-1,yps-1,zps+1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps-1,yps+1,zps+1)=RHS(xps-1,yps+1,zps+1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps+1,yps-1,zps+1)=RHS(xps+1,yps-1,zps+1)+Q_src/16/4
scalars_module_backup.f90:    RHS(xps+1,yps+1,zps+1)=RHS(xps+1,yps+1,zps+1)+Q_src/16/4
scalars_module_backup.f90:        RHS(i,j,zps)=RHS(i,j,zps)+Q_src
scalars_module_backup.f90:      released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src*dx*dy*dz) &
scalars_module_backup.f90:              -0.5_rprec*((jt_total-ini_src)*dt*Q_src*dx*dy*dz)
scalars_module_backup.f90:      released=((end_src-ini_src)*dt*Q_src*dx*dy*dz)
scalars_module_backup.f90:           RHS(xps,yps,zps_local)=RHS(xps,yps,zps_local)+(0.75_rprec)*Q_src
scalars_module_backup.f90:           RHS(xps,yps,zps_local+1)=RHS(xps,yps,zps_local+1)+(0.25_rprec)*Q_src
scalars_module_backup.f90:                 RHS(i,j,zps_local)=RHS(i,j,zps_local)+Q_src
scalars_module_backup.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module_backup.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
scalars_module_backup.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
scalars_module_backup.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
scalars_module_BC.f90:          released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src(ipcon)*dx*dy*dz) &
scalars_module_BC.f90:               -0.5_rprec*((jt_total-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_BC.f90:          released=((end_src-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_BC.f90:            +(0.75_rprec)*Q_src(ipcon)
scalars_module_BC.f90:            +(0.25_rprec)*Q_src(ipcon)
scalars_module_BC.f90:                 RHS(i,j,zps_local,ipcon)=RHS(i,j,zps_local,ipcon)+Q_src(ipcon)
scalars_module_BC.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module_BC.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
scalars_module_BC.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
scalars_module_BC.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
scalars_module_crossflow.f90:          released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src(ipcon)*dx*dy*dz) &
scalars_module_crossflow.f90:               -0.5_rprec*((jt_total-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_crossflow.f90:          released=((end_src-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_crossflow.f90:                    RHS(jx,jy,zps_local,ipcon)=RHS(jx,jy,zps_local,ipcon)+source_fluc*(0.75_rprec)*Q_src(ipcon) & 
scalars_module_crossflow.f90:                    RHS(jx,jy,zps_local+1,ipcon)=RHS(jx,jy,zps_local+1,ipcon)+source_fluc*(0.25_rprec)*Q_src(ipcon) & 
scalars_module_crossflow.f90:              RHS(xps,yps,zps_local,ipcon)=RHS(xps,yps,zps_local,ipcon)+(0.75_rprec)*Q_src(ipcon)*source_fluc
scalars_module_crossflow.f90:              RHS(xps,yps,zps_local+1,ipcon)=RHS(xps,yps,zps_local+1,ipcon)+(0.25_rprec)*Q_src(ipcon)*source_fluc
scalars_module_crossflow.f90:                 RHS(i,j,zps_local,ipcon)=RHS(i,j,zps_local,ipcon)+Q_src(ipcon)
scalars_module_crossflow.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module_crossflow.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
scalars_module_crossflow.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
scalars_module_crossflow.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
scalars_module_crossflow_mod.f90:          released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src(ipcon)*dx*dy*dz) &
scalars_module_crossflow_mod.f90:               -0.5_rprec*((jt_total-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_crossflow_mod.f90:          released=((end_src-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module_crossflow_mod.f90:                    RHS(jx,jy,zps_local,ipcon)=RHS(jx,jy,zps_local,ipcon)+source_fluc*(0.75_rprec)*Q_src(ipcon) & 
scalars_module_crossflow_mod.f90:                    RHS(jx,jy,zps_local+1,ipcon)=RHS(jx,jy,zps_local+1,ipcon)+source_fluc*(0.25_rprec)*Q_src(ipcon) & 
scalars_module_crossflow_mod.f90:              RHS(xps,yps,zps_local,ipcon)=RHS(xps,yps,zps_local,ipcon)+(0.75_rprec)*Q_src(ipcon)*source_fluc
scalars_module_crossflow_mod.f90:              RHS(xps,yps,zps_local+1,ipcon)=RHS(xps,yps,zps_local+1,ipcon)+(0.25_rprec)*Q_src(ipcon)*source_fluc
scalars_module_crossflow_mod.f90:                 RHS(i,j,zps_local,ipcon)=RHS(i,j,zps_local,ipcon)+Q_src(ipcon)
scalars_module_crossflow_mod.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module_crossflow_mod.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
scalars_module_crossflow_mod.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
scalars_module_crossflow_mod.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
scalars_module.f90:          released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src(ipcon)*dx*dy*dz) &
scalars_module.f90:               -0.5_rprec*((jt_total-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module.f90:          released=((end_src-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
scalars_module.f90:                    RHS(jx,jy,zps_local,ipcon)=RHS(jx,jy,zps_local,ipcon)+source_fluc*(0.75_rprec)*Q_src(ipcon) & 
scalars_module.f90:                    RHS(jx,jy,zps_local+1,ipcon)=RHS(jx,jy,zps_local+1,ipcon)+source_fluc*(0.25_rprec)*Q_src(ipcon) & 
scalars_module.f90:              RHS(xps,yps,zps_local,ipcon)=RHS(xps,yps,zps_local,ipcon)+(0.75_rprec)*Q_src(ipcon)*source_fluc
scalars_module.f90:              RHS(xps,yps,zps_local+1,ipcon)=RHS(xps,yps,zps_local+1,ipcon)+(0.25_rprec)*Q_src(ipcon)*source_fluc
scalars_module.f90:                 RHS(i,j,zps_local,ipcon)=RHS(i,j,zps_local,ipcon)+Q_src(ipcon)
scalars_module.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
scalars_module.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
scalars_module.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
scalars_module.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(0.11875_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(1.0_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter:: Q_src=(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter :: Q_src=(0.05_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star))
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(3.74e7_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
t.param_crossflow_Coriol3D.f90:real(kind=rprec),parameter,dimension(npcon) :: Q_src=(/(2.1e-6_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)),&
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(2.1e-6_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
t.param_crossflow_Coriol3D.f90:!real(kind=rprec),parameter,dimension(npcon):: Q_src=(/(6.448e-2_rprec/(dx*z_i*dy*z_i*dz*z_i))/(PCon_scale/(z_i/u_star)), &
t.param_crossflow_Coriol3D.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y)
t.param_crossflow_Coriol3D.f90:!DY if GAUSSIAN_SOURCE_FLAG=.true.: smoothly distribute Q_src by Gaussian filter with b0_plume in (x,y) and h0_plume in (z)
t.scalars_module_crossflow_mod.f90:          released=1.5_rprec*((jt_total-ini_src+1)*dt*Q_src(ipcon)*dx*dy*dz) &
t.scalars_module_crossflow_mod.f90:               -0.5_rprec*((jt_total-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
t.scalars_module_crossflow_mod.f90:          released=((end_src-ini_src)*dt*Q_src(ipcon)*dx*dy*dz)
t.scalars_module_crossflow_mod.f90:                    RHS(jx,jy,zps_local,ipcon)=RHS(jx,jy,zps_local,ipcon)+source_fluc*(0.75_rprec)*Q_src(ipcon) & 
t.scalars_module_crossflow_mod.f90:                    RHS(jx,jy,zps_local+1,ipcon)=RHS(jx,jy,zps_local+1,ipcon)+source_fluc*(0.25_rprec)*Q_src(ipcon) & 
t.scalars_module_crossflow_mod.f90:              RHS(xps,yps,zps_local,ipcon)=RHS(xps,yps,zps_local,ipcon)+(0.75_rprec)*Q_src(ipcon)*source_fluc
t.scalars_module_crossflow_mod.f90:              RHS(xps,yps,zps_local+1,ipcon)=RHS(xps,yps,zps_local+1,ipcon)+(0.25_rprec)*Q_src(ipcon)*source_fluc
t.scalars_module_crossflow_mod.f90:                 RHS(i,j,zps_local,ipcon)=RHS(i,j,zps_local,ipcon)+Q_src(ipcon)
t.scalars_module_crossflow_mod.f90:!    RHS(xps,yps,zps)=RHS(xps,yps,zps)+Q_src
t.scalars_module_crossflow_mod.f90:!DY           RHS(xps,yps,zps)=RHS(xps,yps,zps)+(0.75_rprec)*Q_src
t.scalars_module_crossflow_mod.f90:!DY           RHS(xps,yps,zps+1)=RHS(xps,yps,zps+1)+(0.25_rprec)*Q_src
t.scalars_module_crossflow_mod.f90:!DY                 RHS(i,j,zps)=RHS(i,j,zps)+Q_src
