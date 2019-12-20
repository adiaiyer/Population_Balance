! this is the w-node version
!--provides Cs_opt2 1:nz
!--MPI: requires u,v on 0:nz, except bottom node 1:nz
subroutine lagrange_Ssim(S11,S12,S13,S22,S23,S33)
! standard dynamic model to calculate the Smagorinsky coefficient
! this is done layer-by-layer to save memory
! everything is done to be on uv-nodes
! -note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary. mem requirement ~ n^2, not n^3
use types,only:rprec
use param
use sim_param,only:u,v,w
use sgsmodule,only:F_LM,F_MM,Beta,Cs_opt2,opftime,Beta_avg,Betaclip_avg
use test_filtermodule
use immersedbc,only:n_bldg,bldg_pts,building_interp
$if ($TREES)
  use trees, only : trees_lag_dyn, trees_Cs_lag_dyn
$endif
implicit none

real (rprec), dimension(ld,ny,nz) :: S11,S12,S13,S22,S23,S33

logical, parameter :: DEBUG = .false.

integer :: jx,jy,jz
integer :: i, px, py, lx, ly, lz
real (rprec) :: w_nod(ld, ny)
!real(kind=rprec), dimension(ld,ny,nz+1) :: w_nod
real(kind=rprec) :: u_rms,v_rms,uw_rms,w_rms,u_av,v_av,w_av
real(kind=rprec), dimension(nx,ny) :: u_pr,v_pr,w_pr
real(kind=rprec) :: fractus

real(kind=rprec), dimension(ld,ny) :: L11,L12,L13,L22,L23,L33
real(kind=rprec), dimension(ld,ny) :: M11,M12,M13,M22,M23,M33
real(kind=rprec), dimension(ld,ny) :: fourbeta

real(kind=rprec), dimension(nz) :: LMvert,MMvert

real(kind=rprec), dimension(ld,ny) :: LM,MM,Tn,epsi,dumfac

real(kind=rprec), dimension(ld,ny) :: S_bar,S11_bar,S12_bar,&
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar, S_S12_bar,&
     S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar

real(kind=rprec), dimension(ld,ny) :: u_bar,v_bar,w_bar
!     hat is tilda, bar, hat, three filterings
real(kind=rprec), dimension(ld,ny) :: S
real(kind=rprec) :: delta,const
real(kind=rprec) :: lagran_dt,opftdelta,powcoeff

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real (rprec), dimension(ld, ny, $lbz:nz) :: u_temp, v_temp

logical, parameter :: write_output = .false.

logical, save :: F_LM_MM_init = .false.

!---------------------------------------------------------------------

if (DEBUG) write (*, *) 'started lagrange_Ssim'
!     do jz = 1, nz
!       print *,'jt,Cs_opt2,F_LM,F_MM',jt,sum(Cs_opt2(1:nx,1:ny,jz))/float(nx*ny)&
!      ,sum(F_LM(1:nx,1:ny,jz))/float(nx*ny),sum(F_MM(1:nx,1:ny,jz))/float(nx*ny)
!     end do


delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)
!TS Add the opftdelta
!TS FOR BUILDING (opftime=1*1.5)
opftdelta = opftime*delta
powcoeff = -1._rprec/8._rprec

u_temp=u
v_temp=v

if(use_bldg)then
  !--these calls smooth u_temp, v_temp, w, Sij
  call building_interp(u_temp,v_temp,w,.25_rprec,1)
  call building_interp(S11,S12,S13,.25_rprec,1)
  call building_interp(S22,S23,S33,.25_rprec,1)

   !--zero F_LM in & near bldgs
   do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   do jz=1,lz
   do jy=py-1,py+ly+1
   do jx=px-1,px+lx+1
      F_LM(jx,jy,jz)=0._rprec
   enddo
   enddo
   !--Neumann condition on F_MM around buildings
   F_MM(px,py+1:py+ly-1,jz)=F_MM(px-1,py+1:py+ly-1,jz)
   F_MM(px+lx,py+1:py+ly-1,jz)=F_MM(px+ly+1,py+1:py+ly-1,jz)
   F_MM(px:px+lx,py,jz)=F_MM(px:px+lx,py-1,jz)
   F_MM(px:px+lx,py+ly,jz)=F_MM(px:px+lx,py+ly+1,jz)
   enddo
   enddo
endif

$if ($TREES)
  call trees_lag_dyn (S11, S12, S13, S22, S23, S33)  !--Beta treatment in here
$else
  !      Beta(:,:,1) = (1.d0-0.65d0*Exp(-0.7d0*0.5d0*dz/delta))
  !	Beta(ld,:,1) = 0.0
  !	Beta(ld-1,:,1) = 0.0
  !	do jz= 2,nz
  !	Beta(:,:,jz) = (1.d0-0.65d0*Exp(-0.7d0*(jz-1)*dz/delta))
  !	Beta(ld,:,jz) = 0.0
  !	Beta(ld-1,:,jz) = 0.0			
  !  	end do
  Beta = 1._rprec
$endif

call interpolag_Ssim()

!! (this is the time step used in the lagrangian computations)
lagran_dt=dt*real(cs_count,kind=rprec)
fractus= 1._rprec/real(ny*nx,kind=rprec)
const = 2._rprec*delta**2

!--need to get rid of this-it is waste of mem/cpu
!do jz =1,nz-1
!! w_nod is the value of w interpolated on uvp nodes
!   w_nod(:,:,jz) = (w(:,:,jz)+w(:,:,jz+1))*.5_rprec
!end do
!w_nod(:,:,nz) = w(:,:,nz-1)  !--this makes NO sense

do jz = 1,nz
! using L_ij as temp storage here
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
        (jz == 1) ) then
     !!! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
     ! put on uvp-nodes
     u_bar(:,:) = u_temp(:,:,1)! first uv-node
     v_bar(:,:) = v_temp(:,:,1)! first uv-node
     w_bar(:,:) = .25_rprec*w(:,:,2)! first uv-node
   else  ! w-nodes
     u_bar(:,:) = .5_rprec*(u_temp(:,:,jz) + u_temp(:,:,jz-1))
     v_bar(:,:) = .5_rprec*(v_temp(:,:,jz) + v_temp(:,:,jz-1))
     w_bar(:,:) = w(:,:,jz)
   end if
   L11=u_bar*u_bar
   L12=u_bar*v_bar
   L13=u_bar*w_bar
   L23=v_bar*w_bar
   L22=v_bar*v_bar
   L33=w_bar*w_bar

   call test_filter(u_bar,G_test)
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)
   call test_filter(L11,G_test)  ! in-place filtering
   L11 = L11 - u_bar*u_bar
   call test_filter(L12,G_test)
   L12 = L12 - u_bar*v_bar
   call test_filter(L13,G_test)
   L13 = L13 - u_bar*w_bar
   call test_filter(L22,G_test)
   L22 = L22 - v_bar*v_bar
   call test_filter(L23,G_test)
   L23 = L23 - v_bar*w_bar
   call test_filter(L33,G_test)
   L33 = L33 - w_bar*w_bar
! calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2+S22(:,:,jz)**2+S33(:,:,jz)**2+&
        2._rprec*(S12(:,:,jz)**2+S13(:,:,jz)**2+S23(:,:,jz)**2)))
! S_ij already on w-nodes
   S11_bar(:,:) = S11(:,:,jz)  
   S12_bar(:,:) = S12(:,:,jz)  
   S13_bar(:,:) = S13(:,:,jz)  
   S22_bar(:,:) = S22(:,:,jz)  
   S23_bar(:,:) = S23(:,:,jz)  
   S33_bar(:,:) = S33(:,:,jz)

   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   S_S11_bar(:,:) = S(:,:)*S11(:,:,jz)
   S_S12_bar(:,:) = S(:,:)*S12(:,:,jz)
   S_S13_bar(:,:) = S(:,:)*S13(:,:,jz)
   S_S22_bar(:,:) = S(:,:)*S22(:,:,jz)
   S_S23_bar(:,:) = S(:,:)*S23(:,:,jz)
   S_S33_bar(:,:) = S(:,:)*S33(:,:,jz)

   call test_filter(S_S11_bar,G_test)
   call test_filter(S_S12_bar,G_test)
   call test_filter(S_S13_bar,G_test)
   call test_filter(S_S22_bar,G_test)
   call test_filter(S_S23_bar,G_test)
   call test_filter(S_S33_bar,G_test)     

! now put beta back into M_ij
   fourbeta=4._rprec*Beta(:,:,jz)

   M11 = const*(S_S11_bar - fourbeta*S_bar*S11_bar)
   M12 = const*(S_S12_bar - fourbeta*S_bar*S12_bar)
   M13 = const*(S_S13_bar - fourbeta*S_bar*S13_bar)
   M22 = const*(S_S22_bar - fourbeta*S_bar*S22_bar)
   M23 = const*(S_S23_bar - fourbeta*S_bar*S23_bar)
   M33 = const*(S_S33_bar - fourbeta*S_bar*S33_bar)

   LM=L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
   MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)

! not an efficient initialization
   if (inilag) then
     if ((.not. F_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
       print *,'F_MM and F_LM initialized' 
       F_MM (:,:,jz) = MM
       F_LM (:,:,jz) = 0.025_rprec*MM
       F_MM(ld-1:ld,:,jz)=1._rprec
       F_LM(ld-1:ld,:,jz)=1._rprec

       if (jz == nz) F_LM_MM_init = .true.
     end if
   end if


    if(inflow)then  !--may need to change this
       Tn=merge(.1_rprec*const*S**2,MM,MM.le..1_rprec*const*S**2)
       MM=Tn
       LM(jx_s-x_relax+1:jx_s,1:ny)=0._rprec
       F_LM(jx_s-x_relax+1:jx_s,1:ny,jz)=0._rprec
    endif

   Tn =max(real(F_LM(:,:,jz)*F_MM(:,:,jz)),real(1E-32))
   Tn = opftdelta*(Tn**powcoeff)
   dumfac = lagran_dt/Tn
   epsi = dumfac / (1.0_rprec+dumfac)

   F_LM(:,:,jz)=(epsi*LM + (1.0_rprec-epsi)*F_LM(:,:,jz))
   F_MM(:,:,jz)=(epsi*MM + (1.0_rprec-epsi)*F_MM(:,:,jz))
   F_LM(:,:,jz)= max(real(1E-32),real(F_LM(:,:,jz)))

   Cs_opt2(:,:,jz) = F_LM(:,:,jz)/F_MM(:,:,jz)
   Cs_opt2(ld,:,jz) = 0._rprec
   Cs_opt2(ld-1,:,jz) = 0._rprec
   Cs_opt2(:,:,jz)= max(real(1E-32),real(Cs_opt2(:,:,jz)))
   Beta_avg(jz)=sum(Beta(1:nx,1:ny,jz))
   Betaclip_avg(jz)=sum(Beta(1:nx,1:ny,jz))

   if (mod(jt,200) == 0) then
      LMvert(jz) = sum(sum(LM,DIM=1),DIM=1)/ny/nx
      MMvert(jz) = sum(sum(MM,DIM=1),DIM=1)/ny/nx
   end if

   ! Note: these are not the rms but the standard dev
   ! Square these values to get the rms
   !--you mean take the square-root?
   if ((write_output) .and. (mod(jt,10)==0 .AND. jt>1)) then
!	all the folowing is on UVP nodes
      u_av  = sum(sum(u(:,:,jz),DIM=1),DIM=1)*fractus
      v_av  = sum(sum(v(:,:,jz),DIM=1),DIM=1)*fractus
      if (jz == nz) then
        w_nod(:, :) = w(:, :, nz)
      else
        w_nod(:, :) = 0.5_rprec * (w(:, :, jz) + w(:, :, jz+1))
      end if
      w_av  = sum(sum(w_nod(:,:),DIM=1),DIM=1)*fractus
      do jx=1,nx
         u_pr (jx,:) = u(jx,:,jz)-u_av
         v_pr (jx,:) = v(jx,:,jz)-v_av
         w_pr (jx,:) = w_nod(jx,:)-w_av
      end do

      u_rms = sqrt(sum(sum((u_pr*u_pr),DIM=1),DIM=1)*fractus)
      v_rms = sqrt(sum(sum((v_pr*v_pr),DIM=1),DIM=1)*fractus)
      w_rms = sqrt(sum(sum((w_pr*w_pr),DIM=1),DIM=1)*fractus)
      uw_rms=- sum(sum((u_pr*w_pr),DIM=1),DIM=1)*fractus

      !--format needs redo here
      write(96,*) jt*dz,u_av,u_rms,v_rms,w_rms,uw_rms

   end if
! this ends the main jz=1-nz loop          
end do

!TS
if(use_bldg)then
!if(.false.)then
   do i=1,n_bldg
   px=bldg_pts(1,i)
   py=bldg_pts(2,i)
   lx=bldg_pts(3,i)
   ly=bldg_pts(4,i)
   lz=bldg_pts(5,i)
   Cs_opt2(px:px+lx,py:py+ly,1:lz)=1.E-24
   enddo
!   do jx=px+lx-2,px+lx+4
!   print*,jx,real(Cs_opt2(jx,py-2:py+ly+2,8))
!   enddo
!   do jx=px+lx-2,px+lx+4
!   print*,jx,real(Cs_opt2(jx,py-2:py+ly+2,14))
!   enddo
!   print*,'maxval',maxval(Cs_opt2(px-5:px+lx+10,py-10:py+ly+10,4:20))
endif

$if ($TREES)
  !--set Cs to epsilon inside trees
  call trees_Cs_lag_dyn ()
$endif

if (DEBUG) write (*, *) 'finished lagrange_Ssim'

end subroutine lagrange_Ssim
