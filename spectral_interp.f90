
! 
!   These routines perform spectral interpolations in x, y or both.
!   They are meant to be used to obtain velocity components at element faces for the
! finite-volume solution of scalar dispersion.
!


SUBROUTINE spect_interp_x(f_c, f_int_x)

  USE types,only:rprec
  USE param,only:lh,nx,ny,nz,dx
  USE fft
  IMPLICIT NONE
  INTEGER::jz

  ! Complex variable for fftw
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(in)   ::f_c
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(inout)::f_int_x
  REAL(kind=rprec)::const, ignore_me

  ! Normalization constant
  const=1._rprec/(nx*ny)
  
  ! Loop over vertical levels
  DO jz=1,nz

    ! Normalization (temporay storage in f_int_x)
    f_int_x(:,:,jz)=const*f_c(:,:,jz)

    ! FFT
    CALL rfftwnd_f77_one_real_to_complex(forw,f_int_x(:,:,jz),ignore_me)

    f_int_x(lh,:,jz)=0._rprec
    f_int_x(:,ny/2+1,jz)=0._rprec
    
    ! Shift in the x-direction
    f_int_x(:,:,jz)=EXP((dx/2._rprec)*eye*kx(:,:))*f_int_x(:,:,jz)

    ! Inverse FFT
    CALL rfftwnd_f77_one_complex_to_real(back,f_int_x(:,:,jz),ignore_me)
    
  END DO
  
END SUBROUTINE spect_interp_x



SUBROUTINE spect_interp_y(f_c, f_int_y)

  USE types,only:rprec
  USE param,only:lh,nx,ny,nz,dy
  USE fft
  IMPLICIT NONE
  INTEGER::jz

  ! Complex variable for fftw
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(in)   ::f_c
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(inout)::f_int_y
  REAL(kind=rprec)::const, ignore_me

  ! Normalization constant
  const=1._rprec/(nx*ny)
  
  ! Loop over vertical levels
  DO jz=1,nz

    ! Normalization (temporay storage in f_int_y)
    f_int_y(:,:,jz)=const*f_c(:,:,jz)

    ! FFT
    CALL rfftwnd_f77_one_real_to_complex(forw,f_int_y(:,:,jz),ignore_me)

    f_int_y(lh,:,jz)=0._rprec
    f_int_y(:,ny/2+1,jz)=0._rprec
    
    ! Shift in the y-direction
    f_int_y(:,:,jz)=EXP((dy/2._rprec)*eye*ky(:,:))*f_int_y(:,:,jz)

    ! Inverse FFT
    CALL rfftwnd_f77_one_complex_to_real(back,f_int_y(:,:,jz),ignore_me)
    
  END DO
  
END SUBROUTINE spect_interp_y



SUBROUTINE spect_interp_xy(f_c, f_int_xy)

  USE types,only:rprec
  USE param,only:lh,nx,ny,nz,dx,dy
  USE fft
  IMPLICIT NONE
  INTEGER::jz

  ! Complex variable for fftw
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(in)   ::f_c
  COMPLEX(kind=rprec),DIMENSION(lh,ny,nz),INTENT(inout)::f_int_xy
  REAL(kind=rprec)::const, ignore_me

  ! Normalization constant
  const=1._rprec/(nx*ny)
  
  ! Loop over vertical levels
  DO jz=1,nz

    ! Normalization (temporay storage in f_int_xy)
    f_int_xy(:,:,jz)=const*f_c(:,:,jz)

    ! FFT
    CALL rfftwnd_f77_one_real_to_complex(forw,f_int_xy(:,:,jz),ignore_me)

    f_int_xy(lh,:,jz)=0._rprec
    f_int_xy(:,ny/2+1,jz)=0._rprec
    
    ! Shift in the x-direction
    f_int_xy(:,:,jz)=EXP((dx/2._rprec)*eye*kx(:,:))*f_int_xy(:,:,jz)
    
    ! Shift in the y-direction
    f_int_xy(:,:,jz)=EXP((dy/2._rprec)*eye*ky(:,:))*f_int_xy(:,:,jz)

    ! Inverse FFT
    CALL rfftwnd_f77_one_complex_to_real(back,f_int_xy(:,:,jz),ignore_me)
    
  END DO
  
END SUBROUTINE spect_interp_xy
