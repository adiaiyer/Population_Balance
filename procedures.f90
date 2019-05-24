module procedures
use types,only:rprec
implicit none

contains


 function integrate(x, y, a, b) result(r)
    !! This function constructs a piecewise cubic Hermitian
    !interpolation of an array y(x) based on
    !! discrete numerical data, and subsequently evaluates the integral
    !of the interpolation in the
    !! range (a,b). Note that the mesh spacing of x does not necessarily
    !have to be uniform.
    real(rprec), intent(in)  :: x(:)        !! Variable x
    real(rprec), intent(in)  :: y(size(x))  !! Function y(x)
    real(rprec), intent(in)  :: a           !! Left endpoint
    real(rprec), intent(in)  :: b           !! Right endpoint
    real(rprec)              :: r !! Integral ∫y(x)·dx

    real(rprec), external    :: dpchqa
    real(rprec)              :: d(size(x))
    integer               :: err

    ! Create a PCHIP interpolation of the input data
!    call dpchez(size(x), x, y, d, .false., 0, 0, err)

    ! Integrate the interpolation in the provided range
 !   r = dpchqa(size(x), x, y, d, a, b, err)
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate

  end function
end module

