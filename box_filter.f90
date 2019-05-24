!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                             box_filter_3D.f90                              ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/25/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is a 3D box filter to be used with the finite volume
!              discretization of the scalar field. The filter is performed at scale
!              2*delta
!
!    ################################################################################
!
!    NOTE:  The vector size (ld:ny+1,nz) was chosen to accommodate all the vectors
!           that will be filtered. The parameters fx, fy and fz indicate the real 
!           size of the specific vector being filtered.
!		 
!    ################################################################################
!

SUBROUTINE box_filter_3D(f_c,f_3Dh,fx,fy,fz)

! Modules
  USE types,only:rprec
  USE param,only:ld,nx,ny,nz

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),INTENT(in)   :: f_c(ld,ny+1,nz)              ! Original field
  REAL(kind=rprec),INTENT(inout):: f_3Dh(ld,ny+1,nz)            ! 3D filtered field
  INTEGER                       :: fx,fy,fz                     ! Valid length of original field


!
! BEGINNING CODE   
!


  ! This is for the boundary nodes (for now, they are not filtered at all)
  f_3Dh=f_c
  f_3Dh=0._rprec

  ! Treat all internal nodes
  f_3Dh(2:fx-1,2:fy-1,2:fz-1)=((1._rprec/8._rprec)*f_c(2:fx-1,2:fy-1,2:fz-1)+ &
			       (1._rprec/16._rprec)*(f_c(1:fx-2,2:fy-1,2:fz-1)+f_c(3:fx,2:fy-1,2:fz-1)+ &
						     f_c(2:fx-1,1:fy-2,2:fz-1)+f_c(2:fx-1,3:fy,2:fz-1)+ &
						     f_c(2:fx-1,2:fy-1,1:fz-2)+f_c(2:fx-1,2:fy-1,3:fz))+ &
			       (1._rprec/32._rprec)*(f_c(1:fx-2,1:fy-2,2:fz-1)+f_c(1:fx-2,3:fy,2:fz-1)+ &
						     f_c(3:fx,1:fy-2,2:fz-1)+f_c(3:fx,3:fy,2:fz-1)+ &
						     f_c(1:fx-2,2:fy-1,1:fz-2)+f_c(3:fx,2:fy-1,1:fz-2)+ &
						     f_c(2:fx-1,1:fy-2,1:fz-2)+f_c(2:fx-1,3:fy,1:fz-2)+ &
						     f_c(1:fx-2,2:fy-1,3:fz)+f_c(3:fx,2:fy-1,3:fz)+ &
						     f_c(2:fx-1,1:fy-2,3:fz)+f_c(2:fx-1,3:fy,3:fz))+ &
			       (1._rprec/64._rprec)*(f_c(1:fx-2,1:fy-2,1:fz-2)+f_c(1:fx-2,3:fy,1:fz-2)+ &
						     f_c(3:fx,1:fy-2,1:fz-2)+f_c(3:fx,3:fy,1:fz-2)+ &
						     f_c(1:fx-2,1:fy-2,3:fz)+f_c(1:fx-2,3:fy,3:fz)+ &
						     f_c(3:fx,1:fy-2,3:fz)+f_c(3:fx,3:fy,3:fz)))

  ! Treat boundary nodes
  ! For now the boundary nodes are left unfiltered  
  
END SUBROUTINE box_filter_3D





!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                          box_filter_3D_4delta.f90                          ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                06/06/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is a 3D box filter to be used with the finite volume
!              discretization of the scalar field. The filter is performed at scale
!              4*delta
!
!    ################################################################################
!
!    NOTE:  The vector size (ld:ny+1,nz) was chosen to accommodate all the vectors
!           that will be filtered. The parameters fx, fy and fz indicate the real 
!           size of the specific vector being filtered.
!		 
!    ################################################################################
!

SUBROUTINE box_filter_3D_4delta(f_c,f_3Dh,fx,fy,fz)

! Modules
  USE types,only:rprec
  USE param,only:ld,nx,ny,nz

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),INTENT(in)   :: f_c(ld,ny+1,nz)              ! Original field
  REAL(kind=rprec),INTENT(inout):: f_3Dh(ld,ny+1,nz)            ! 3D filtered field
  INTEGER                       :: fx,fy,fz                     ! Valid length of original field


!
! BEGINNING CODE   
!


  ! This is for the boundary nodes (for now, they are not filtered at all)
  f_3Dh=f_c
  f_3Dh=0._rprec

  ! Treat all internal nodes
  f_3Dh(3:fx-2,3:fy-2,3:fz-2)=((1._rprec/64._rprec)*(f_c(3:fx-2,3:fy-2,3:fz-2)+ &
                                                     f_c(2:fx-3,3:fy-2,3:fz-2)+f_c(4:fx-1,3:fy-2,3:fz-2)+ &
                                                     f_c(3:fx-2,2:fy-3,3:fz-2)+f_c(3:fx-2,4:fy-1,3:fz-2)+ &
                                                     f_c(2:fx-3,2:fy-3,3:fz-2)+f_c(2:fx-3,4:fy-1,3:fz-2)+ &
                                                     f_c(4:fx-1,2:fy-3,3:fz-2)+f_c(4:fx-1,4:fy-1,3:fz-2)+ &
                                                     f_c(3:fx-2,3:fy-2,2:fz-3)+f_c(3:fx-2,3:fy-2,4:fz-1)+ &
                                                     f_c(2:fx-3,3:fy-2,2:fz-3)+f_c(4:fx-1,3:fy-2,2:fz-3)+ &
                                                     f_c(3:fx-2,2:fy-3,2:fz-3)+f_c(3:fx-2,4:fy-1,2:fz-3)+ &
                                                     f_c(2:fx-3,3:fy-2,4:fz-1)+f_c(4:fx-1,3:fy-2,4:fz-1)+ &
                                                     f_c(3:fx-2,2:fy-3,4:fz-1)+f_c(3:fx-2,4:fy-1,4:fz-1)+ &
                                                     f_c(2:fx-3,2:fy-3,2:fz-3)+f_c(2:fx-3,4:fy-1,2:fz-3)+ &
                                                     f_c(4:fx-1,2:fy-3,2:fz-3)+f_c(4:fx-1,4:fy-1,2:fz-3)+ &
                                                     f_c(2:fx-3,2:fy-3,4:fz-1)+f_c(2:fx-3,4:fy-1,4:fz-1)+ &
                                                     f_c(4:fx-1,2:fy-3,4:fz-1)+f_c(4:fx-1,4:fy-1,4:fz-1))+ &
                              (1._rprec/128._rprec)*(f_c(3:fx-2,3:fy-2,5:fz)+ &
                                                     f_c(2:fx-3,3:fy-2,5:fz)+f_c(4:fx-1,3:fy-2,5:fz)+ &
                                                     f_c(3:fx-2,2:fy-3,5:fz)+f_c(3:fx-2,4:fy-1,5:fz)+ &
                                                     f_c(2:fx-3,2:fy-3,5:fz)+f_c(2:fx-3,4:fy-1,5:fz)+ &
                                                     f_c(4:fx-1,2:fy-3,5:fz)+f_c(4:fx-1,4:fy-1,5:fz)+ &
						     f_c(3:fx-2,3:fy-2,5:fz)+ &
                                                     f_c(2:fx-3,3:fy-2,1:fz-4)+f_c(4:fx-1,3:fy-2,1:fz-4)+ &
                                                     f_c(3:fx-2,2:fy-3,1:fz-4)+f_c(3:fx-2,4:fy-1,1:fz-4)+ &
                                                     f_c(2:fx-3,2:fy-3,1:fz-4)+f_c(2:fx-3,4:fy-1,1:fz-4)+ &
                                                     f_c(4:fx-1,2:fy-3,1:fz-4)+f_c(4:fx-1,4:fy-1,1:fz-4)+ &						     
						     f_c(1:fx-4,3:fy-2,3:fz-2)+f_c(5:fx,3:fy-2,3:fz-2)+ &
						     f_c(1:fx-4,2:fy-3,3:fz-2)+f_c(5:fx,2:fy-3,3:fz-2)+ &
						     f_c(1:fx-4,4:fy-1,3:fz-2)+f_c(5:fx,4:fy-1,3:fz-2)+ &
						     f_c(3:fx-2,1:fy-4,3:fz-2)+f_c(3:fx-2,5:fy,3:fz-2)+ &
						     f_c(2:fx-3,1:fy-4,3:fz-2)+f_c(2:fx-3,5:fy,3:fz-2)+ &
						     f_c(4:fx-1,1:fy-4,3:fz-2)+f_c(4:fx-1,5:fy,3:fz-2)+ &
						     f_c(1:fx-4,3:fy-2,2:fz-3)+f_c(5:fx,3:fy-2,2:fz-3)+ &
						     f_c(1:fx-4,2:fy-3,2:fz-3)+f_c(5:fx,2:fy-3,2:fz-3)+ &
						     f_c(1:fx-4,4:fy-1,2:fz-3)+f_c(5:fx,4:fy-1,2:fz-3)+ &
						     f_c(3:fx-2,1:fy-4,2:fz-3)+f_c(3:fx-2,5:fy,2:fz-3)+ &
						     f_c(2:fx-3,1:fy-4,2:fz-3)+f_c(2:fx-3,5:fy,2:fz-3)+ &
						     f_c(4:fx-1,1:fy-4,2:fz-3)+f_c(4:fx-1,5:fy,2:fz-3)+ &
						     f_c(1:fx-4,3:fy-2,4:fz-1)+f_c(5:fx,3:fy-2,4:fz-1)+ &
						     f_c(1:fx-4,2:fy-3,4:fz-1)+f_c(5:fx,2:fy-3,4:fz-1)+ &
						     f_c(1:fx-4,4:fy-1,4:fz-1)+f_c(5:fx,4:fy-1,4:fz-1)+ &
						     f_c(3:fx-2,1:fy-4,4:fz-1)+f_c(3:fx-2,5:fy,4:fz-1)+ &
						     f_c(2:fx-3,1:fy-4,4:fz-1)+f_c(2:fx-3,5:fy,4:fz-1)+ &
						     f_c(4:fx-1,1:fy-4,4:fz-1)+f_c(4:fx-1,5:fy,4:fz-1))+ &
                              (1._rprec/256._rprec)*(f_c(1:fx-4,1:fy-4,3:fz-2)+f_c(5:fx,1:fy-4,3:fz-2)+ &
			                             f_c(1:fx-4,5:fy,3:fz-2)+f_c(5:fx,5:fy,3:fz-2)+ &
						     f_c(1:fx-4,1:fy-4,2:fz-3)+f_c(5:fx,1:fy-4,2:fz-3)+ &
			                             f_c(1:fx-4,5:fy,2:fz-3)+f_c(5:fx,5:fy,2:fz-3)+ &
						     f_c(1:fx-4,1:fy-4,4:fz-1)+f_c(5:fx,1:fy-4,4:fz-1)+ &
			                             f_c(1:fx-4,5:fy,4:fz-1)+f_c(5:fx,5:fy,4:fz-1)+ &
						     f_c(1:fx-4,3:fy-2,1:fz-4)+f_c(5:fx,3:fy-2,1:fz-4)+ &
						     f_c(1:fx-4,2:fy-3,1:fz-4)+f_c(5:fx,2:fy-3,1:fz-4)+ &
						     f_c(1:fx-4,4:fy-1,1:fz-4)+f_c(5:fx,4:fy-1,1:fz-4)+ &
						     f_c(3:fx-2,1:fy-4,1:fz-4)+f_c(3:fx-2,5:fy,1:fz-4)+ &
						     f_c(2:fx-3,1:fy-4,1:fz-4)+f_c(2:fx-3,5:fy,1:fz-4)+ &
						     f_c(4:fx-1,1:fy-4,1:fz-4)+f_c(4:fx-1,5:fy,1:fz-4)+ &
						     f_c(1:fx-4,3:fy-2,5:fz)+f_c(5:fx,3:fy-2,5:fz)+ &
						     f_c(1:fx-4,2:fy-3,5:fz)+f_c(5:fx,2:fy-3,5:fz)+ &
						     f_c(1:fx-4,4:fy-1,5:fz)+f_c(5:fx,4:fy-1,5:fz)+ &
						     f_c(3:fx-2,1:fy-4,5:fz)+f_c(3:fx-2,5:fy,5:fz)+ &
						     f_c(2:fx-3,1:fy-4,5:fz)+f_c(2:fx-3,5:fy,5:fz)+ &
						     f_c(4:fx-1,1:fy-4,5:fz)+f_c(4:fx-1,5:fy,5:fz))+ &
                              (1._rprec/512._rprec)*(f_c(1:fx-4,1:fy-4,1:fz-4)+f_c(5:fx,1:fy-4,1:fz-4)+ &
			                             f_c(1:fx-4,5:fy,1:fz-4)+f_c(5:fx,5:fy,1:fz-4)+ &						     
						     f_c(1:fx-4,1:fy-4,5:fz)+f_c(5:fx,1:fy-4,5:fz)+ &
			                             f_c(1:fx-4,5:fy,5:fz)+f_c(5:fx,5:fy,5:fz)))

  ! Treat boundary nodes
  ! For now the boundary nodes are left unfiltered  
  
END SUBROUTINE box_filter_3D_4delta





!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                           box_filter_2D_xy.f90                             ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                05/25/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is a horizontal 2D box filter to be used with the finite 
!              volume discretization of the scalar field.
!
!    ################################################################################
!
!    NOTE:  The vector size (ld:ny+1,nz) was chosen to accommodate all the vectors
!           that will be filtered. The parameters fx, fy and fz indicate the real 
!           size of the specific vector being filtered.
!		 
!    ################################################################################
!

SUBROUTINE box_filter_2D_xy(f_c,f_2Dxy,fx,fy,fz)

! Modules
  USE types,only:rprec
  USE param,only:ld,nx,ny,nz

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),INTENT(in)   :: f_c(ld,ny+1,nz)              ! Original field
  REAL(kind=rprec),INTENT(inout):: f_2Dxy(ld,ny+1,nz)           ! 2D filtered field
  INTEGER                       :: fx,fy,fz                     ! Valid length of original field


!
! BEGINNING CODE   
!


  ! This is for the boundary nodes (for now, they are not filtered at all)
  f_2Dxy=f_c

  ! Treat all internal nodes
  f_2Dxy(2:fx-1,2:fy-1,1:fz)=((1._rprec/4._rprec)*f_c(2:fx-1,2:fy-1,1:fz)+ &
                              (1._rprec/8._rprec)*(f_c(1:fx-2,2:fy-1,1:fz)+f_c(3:fx,2:fy-1,1:fz)+ &
                                                   f_c(2:fx-1,1:fy-2,1:fz)+f_c(2:fx-1,3:fy,1:fz))+ &
                              (1._rprec/16._rprec)*(f_c(1:fx-2,1:fy-2,1:fz)+f_c(1:fx-2,3:fy,1:fz)+ &
                                                    f_c(3:fx,1:fy-2,1:fz)+f_c(3:fx,3:fy,1:fz)))

  ! Treat boundary nodes
  ! For now the boundary nodes are left unfiltered  
  
END SUBROUTINE box_filter_2D_xy





!
!    ################################################################################
!    ################################################################################
!    ##                                                                            ## 
!    ##                        box_filter_2D_xy_4delta.f90                         ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                                                            ##
!    ##                        The Johns Hopkins University                        ##
!    ##                                                                            ##
!    ##                                06/07/2007                                  ##
!    ##                                                                            ##
!    ################################################################################
!    ################################################################################
!
!    ################################################################################
!
!    PURPOSE:  This routine is a horizontal 2D box filter to be used with the finite 
!              volume discretization of the scalar field. The filter is performed at
!              scale 4*delta
!
!    ################################################################################
!
!    NOTE:  The vector size (ld:ny+1,nz) was chosen to accommodate all the vectors
!           that will be filtered. The parameters fx, fy and fz indicate the real 
!           size of the specific vector being filtered.
!		 
!    ################################################################################
!

SUBROUTINE box_filter_2D_xy_4delta(f_c,f_2Dxy,fx,fy,fz)

! Modules
  USE types,only:rprec
  USE param,only:ld,nx,ny,nz

! No implicit variables
  IMPLICIT NONE

! Main variables
  REAL(kind=rprec),INTENT(in)   :: f_c(ld,ny+1,nz)              ! Original field
  REAL(kind=rprec),INTENT(inout):: f_2Dxy(ld,ny+1,nz)           ! 2D filtered field
  INTEGER                       :: fx,fy,fz                     ! Valid length of original field


!
! BEGINNING CODE   
!


  ! This is for the boundary nodes (for now, they are not filtered at all)
  f_2Dxy=f_c

  ! Treat all internal nodes
  f_2Dxy(3:fx-2,3:fy-2,1:fz)=((1._rprec/16._rprec)*(f_c(3:fx-2,3:fy-2,1:fz)+ &
                                                    f_c(2:fx-3,3:fy-2,1:fz)+f_c(4:fx-1,3:fy-2,1:fz)+ &
                                                    f_c(3:fx-2,2:fy-3,1:fz)+f_c(3:fx-2,4:fy-1,1:fz)+ &
                                                    f_c(2:fx-3,2:fy-3,1:fz)+f_c(2:fx-3,4:fy-1,1:fz)+ &
                                                    f_c(4:fx-1,2:fy-3,1:fz)+f_c(4:fx-1,4:fy-1,1:fz))+ &
			      (1._rprec/32._rprec)*(f_c(1:fx-4,3:fy-2,1:fz)+f_c(5:fx,3:fy-2,1:fz)+ &
			                            f_c(1:fx-4,2:fy-3,1:fz)+f_c(5:fx,2:fy-3,1:fz)+ &
			                            f_c(1:fx-4,4:fy-1,1:fz)+f_c(5:fx,4:fy-1,1:fz)+ &					    
						    f_c(3:fx-2,1:fy-4,1:fz)+f_c(3:fx-2,5:fy,1:fz)+ &
			                            f_c(2:fx-3,1:fy-4,1:fz)+f_c(2:fx-3,5:fy,1:fz)+ &
			                            f_c(4:fx-1,1:fy-4,1:fz)+f_c(4:fx-1,5:fy,1:fz))+ &
			      (1._rprec/64._rprec)*(f_c(1:fx-4,1:fy-4,1:fz)+f_c(5:fx,1:fy-4,1:fz)+ &
						    f_c(1:fx-4,5:fy,1:fz)+f_c(5:fx,5:fy,1:fz)))

  ! Treat boundary nodes
  ! For now the boundary nodes are left unfiltered  
  
END SUBROUTINE box_filter_2D_xy_4delta
