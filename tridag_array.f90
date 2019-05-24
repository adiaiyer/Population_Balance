subroutine tridag_array(a, b, c, r, u)
! sc: modified for complex!
use types,only:rprec
use param
implicit none

real(kind=rprec),dimension(lh, ny, nz+1),intent(in):: a, b, c
complex(kind=rprec),dimension(lh, ny, nz+1),intent(in) :: r
complex(kind=rprec),dimension(lh, ny, nz+1),intent(out):: u

integer, parameter :: n = nz+1

logical, parameter :: DEBUG = .false.

character (64) :: fmt
integer::jx, jy, j, j_min, j_max
real(kind=rprec)::bet(lh, ny)
real(kind=rprec),dimension(lh, ny, nz+1)::gam

!--want to skip ny/2+1 and 1, 1

$if ($MPI)
  !--wait for c(:,:,1), bet(:,:), u(:,:,1) from "down"
  !--may want irecv here with a wait at the end
  call mpi_recv (c(1, 1, 1), lh*ny, MPI_RPREC, down, 1, comm, status, ierr)
  call mpi_recv (bet(1, 1), lh*ny, MPI_RPREC, down, 2, comm, status, ierr)
  call mpi_recv (u(1, 1, 1), lh*ny, MPI_CPREC, down, 3, comm, status, ierr)
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  do jy = 1, ny
    do jx = 1, lh-1

      if (b(jx, jy, 1) == 0._rprec) then
        write (*, *) 'tridag_array: rewrite eqs, jx, jy= ', jx, jy
        stop
      end if

    end do
  end do

  bet = b(:, :, 1)
  u(:, :, 1) = r(:, :, 1) / bet

  j_min = 1  !--this is only for backward pass
else
  j_min = 2  !--this is only for backward pass
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  j_max = n
else
  j_max = n-1
end if

do j = 2, j_max

  do jy = 1, ny

    if (jy == ny/2+1) cycle

    do jx = 1, lh-1

      if (jx*jy == 1) cycle

      gam(jx, jy, j) = c(jx, jy, j-1) / bet(jx, jy)
      bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

      if (bet(jx, jy) == 0._rprec) then
        write (*, *) 'tridag_array failed at jx,jy,j=', jx, jy, j
        write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j),  &
                     c(jx, jy, j), gam(jx, jy, j), bet(jx, jy)
        stop
      end if

      u(jx, jy, j) = (r(jx, jy, j) - a(jx, jy, j) * u(jx, jy, j-1)) /  &
                     bet(jx, jy)

    end do

  end do

  if (DEBUG) then
     fmt = '(i0,a,i0,1x,"(",es12.5,", ",es12.5,")")'
     write (*, fmt) coord, ': P1: j, u(2,2,j) = ', j, u(2, 2, j)
     fmt = '(i0,a,i0,1x,es12.5)'
     write (*, fmt) coord, ': P1: j, gam(2,2,j) = ', j, gam(2, 2, j)
     write (*, fmt) coord, ': P1: j, bet(2,2) = ', j, bet(2, 2)
     write (*, fmt) coord, ': P1: j, a(2,2,j) = ', j, a(2, 2, j)
     write (*, fmt) coord, ': P1: j, b(2,2,j) = ', j, b(2, 2, j)
     write (*, fmt) coord, ': P1: j, c(2,2,j) = ', j, c(2, 2, j)
     fmt = '(i0,a,i0,1x,"(",es12.5,", ",es12.5,")")'
     write (*, fmt) coord, ': P1: j, r(2,2,j) = ', j, r(2, 2, j)
  end if

end do

$if ($MPI)
  !--send c(n-1), bet, u(n-1) to "up"
  !--may not want blocking sends here
  call mpi_send (c(1, 1, n-1), lh*ny, MPI_RPREC, up, 1, comm, ierr)
  call mpi_send (bet(1, 1), lh*ny, MPI_RPREC, up, 2, comm, ierr)
  call mpi_send (u(1, 1, n-1), lh*ny, MPI_CPREC, up, 3, comm, ierr)

  !--wait for u(n), gam(n) from "up"
  call mpi_recv (u(1, 1, n), lh*ny, MPI_CPREC, up, 4, comm, status, ierr)
  call mpi_recv (gam(1, 1, n), lh*ny, MPI_RPREC, up, 5, comm, status, ierr)
$endif

if (DEBUG) then
  fmt = '(i0,a,1x,"(",es12.5,", ",es12.5,")")'
  write (*, fmt) coord, ': P2: u(2,2,n) = ', u(2, 2, n)
  write (*, fmt) coord, ': P2: u(2,2,n-1) = ', u(2, 2, n-1)
  fmt = '(i0,a,1x,es12.5)'
  write (*, fmt) coord, ': P2: gam(2,2,n) = ', gam(2, 2, n)
end if

!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!  j_max = n-1
!else
!  j_max = n
!end if

do j = n-1, j_min, -1

  if (DEBUG) then
    fmt = '(i0,a,i0,1x,"(",es12.5,", ",es12.5,")")'
    write (*, fmt) coord, ': P2: j, u_i(2,2,j) = ', j, u(2, 2, j)
    write (*, fmt) coord, ': P2: j, u(2,2,j+1) = ', j, u(2, 2, j+1)
    fmt = '(i0,a,i0,1x,es12.5)'
    write (*, fmt) coord, ': P2: j, gam(2,2,j+1) = ', j, gam(2, 2, j+1)
  end if

  !--intend on removing cycle statements/repl with something faster
  do jy = 1, ny

    if (jy == ny/2+1) cycle

    do jx = 1, lh-1

      if (jx*jy == 1) cycle

      u(jx, jy, j) = u(jx, jy, j) - gam(jx, jy, j+1) * u(jx, jy, j+1)

    end do

  end do

  if (DEBUG) then
    fmt = '(i0,a,i0,"(",es12.5,", ",es12.5,")")'
    write (*, fmt) coord, ': P2: j, u_f(2,2,j) = ', j, u(2, 2, j)
  end if

end do

$if ($MPI)
  !--send u(2), gam(2) to "down"
  call mpi_send (u(1, 1, 2), lh*ny, MPI_CPREC, down, 4, comm, ierr)
  call mpi_send (gam(1, 1, 2), lh*ny, MPI_RPREC, down, 5, comm, ierr)
$endif

end subroutine tridag_array
