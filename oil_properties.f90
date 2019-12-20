module oil_prop

  use types,only:rprec
  use param,only:rho_c,rho_d,mu_c,g,npcon
  $IF ($MPI)
    use mpi
  $endif
  
  implicit none


  
  real(rprec),dimension(npcon)  :: diameter,A,w_0
  integer :: k
  A = 0.15_rprec*(rho_c*diameter/nu_c)**(0.687_rprec)
  w_0 = (rho_c-rho_d)/(18._rprec*mu_c)*g*diameter**2
  
  diameter (1)=dmin
  diameter (npcon) = dmax 
  do ipcon=2,npcon-1
    diameter(ipcon) = dmin*d_diameter**(ipcon-1)
  enddo   
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rise_vel()  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
  
    real(rprec) :: x,x0,fx
    integer :: iters
  
    debug=.true.
  
    x0 = 0.2_rprec
    do k=1,npcon
    call solve(func,func_p,x_0,x,iters,debug,k)
  
    print 11, x, iters
11  format('solver returns x = ', e22.15, ' after', i3, ' iterations')
    fx = func(x)
    print 12, fx
12  format('the value of f(x) is ', e22.15)  
    end subroutine rise_vel
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine solve(f,fp,x0,x,iters,debug,ipcon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    use types, only:rprec
    !! Estimate the zero of f(x) using Newton's method
    !  Input :
    !!   f: The function to find root
    !!  fp: function returning derivative of f
    !!  x0 : initial guess
    !!  debug: logical, prints iterations if debug=.true. 
    !!  return:
    !!  x : the estimate of the zero if converges
    !!  iters : no. of iterations
    
    implicit none
    
    real(rprec), intent(in) :: x0
    real(rprec), external :: f,fp
    logical, intent(in) :: debug
    real(rprec), intent(out) :: x
    integer, intent(in) :: ipcon
    integer, intent(out) :: iters
    
    ! Declaration of variables
    
    real(rprec) :: deltax,fx,fxprime
    integer :: i
    
    ! initial guess
    
    x=x0
    
    if(debug) then
      print 14, x
    
    14 format('Initial guess: x=', e22.15)
    endif
    
    do i=1,maxiter
    
      fx = f(x,ipcon)
      fxprime = fp(x,ipcon)
    
    
      if(abs(fx) < tol) then
        exit
      endif
    
      deltax = fx/fxprime
    
      x = x-deltax
    
      if(debug) then
        print 12,k,x  
    12  format('After', i3, ' iterations, x = ', e22.15)
      endif
    
    enddo
    
    
    if(k>maxiter) then
    
      fx = f(x,ipcon)
      if(abs(fx) >tol) then
       print *, '** did not converege'
       endif
      endif
    
     iters = k-1
    
    end subroutine solve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(rprec) FUNCTION func(w,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Define function to find root
 
    implicit none
    real(rprec), intent(in) :: w
    real(rprec), intent(out) :: func
    integer, intent(in) :: k
  
    func = w - w_0(k)*(1._rprec+A(k)*w**(0.687_rprec))**(-1._rprec)
    
    end 
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(rprec) FUNCTION func_p(w,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate derivative of function
  
    implicit none
    real(rprec), intent(in) ::w
    real(rprec), intent(out) :: func_p
    integer, intent(in) :: k
    func_p = 1._rprec + 0.687_rprec*w_0(k)*A(k)*(1+A(k)*w**(0.687_rprec))*w**(-0.313_rprec)
  
    end 
