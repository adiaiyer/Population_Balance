module outflow_fringe
use types,only:rprec
use param,only:nx
implicit none
real(kind=rprec),dimension(nx)::fringe

contains
subroutine setfringe()
use param
implicit none
real(kind=rprec):: xl_fringe, a, b, x

integer :: i

xl_fringe = (ix_fringe-1)*dx
a = L_x - xl_fringe
b = xl_fringe / ( L_x - xl_fringe )

fringe = 0._rprec

do i=ix_fringe+1,nx
   x = (i-1)*dx
   fringe(i) = 0.5_rprec * ( 1._rprec - cos(pi*(x/a-b)))
enddo

open(unit=100)
do i=1,nx
   write(100,*) i,fringe(i),1._rprec-fringe(i)
enddo
close(100)

end subroutine setfringe
end module outflow_fringe
