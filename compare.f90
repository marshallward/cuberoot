use iso_fortran_env, only : real64, real128
use cubes
implicit none

integer :: k, n
real(kind=real128) :: q
real, allocatable :: x(:)
integer :: i

real, allocatable :: nbits(:)

integer :: io_unit

!print *, "Resolution (as k in 2**-k)?"
!read (*,*) k
k = 17
n = 1 + 7 * 2**(k-3)
allocate(x(n))
x = [(i*2.**-k, i=2**(k-3), 2**k)]

!n = 1+10**7
!allocate(x(n))
!x = [(i*1e-3, i=0,n)]

!print '(*(a26))', "x", "Quad Newton", "x**(1./3.)", "Double Newton", "Newton No-div"
allocate(nbits(8))
nbits(:) = 0

! Absolute error with respect to a *rounded* quadratic precision number.
open(file='err_quad.txt', newunit=io_unit)
do i = 1,n
  q = cuberoot_newton_quad(real(x(i), kind=real128))
  write(io_unit, '(*(es24.16,2x))') &
    x(i), &
    real(q, kind=real64), &

    ! Absolute error of x**(1/3) relative to quad
    real(abs(x(i)**(1./3.) - q), kind=real64), &
    real(abs(cuberoot_newton(x(i)) - q), kind=real64), &
    real(abs(cuberoot_newton_nodiv(x(i)) - q), kind=real64), &
    real(abs(cuberoot_halley(x(i)) - q), kind=real64), &
    real(abs(cuberoot_halley_nodiv(x(i)) - q), kind=real64), &
    real(abs(cuberoot_lagny(x(i)) - q), kind=real64), &
    real(abs(cbrt_ac(x(i)) - q), kind=real64), &
    real(abs(cuberoot_final(x(i)) - q), kind=real64)

  ! Count the total error (to be converted into bits)
  nbits(1) = nbits(1) + abs(x(i)**(1./3.) - real(q))
  nbits(2) = nbits(2) + abs(cuberoot_newton(x(i)) - real(q))
  nbits(3) = nbits(3) + abs(cuberoot_newton_nodiv(x(i)) - real(q))
  nbits(4) = nbits(4) + abs(cuberoot_halley(x(i)) - real(q))
  nbits(5) = nbits(5) + abs(cuberoot_halley_nodiv(x(i)) - real(q))
  nbits(6) = nbits(6) + abs(cuberoot_lagny(x(i)) - real(q))
  nbits(7) = nbits(7) + abs(cbrt_ac(x(i)) - real(q))
  nbits(8) = nbits(8) + abs(cuberoot_final(x(i)) - real(q))
enddo

! Divide by the ULP error (more or less)
print *, "Total ULP errors (...?)"
print '(*(a14))', "x**(1./3.)", "Newton", "Nwtn-Nodiv", "Halley", &
    "Hly-Nodiv", "Lagny/Leroy", "AC", "Final"
print '(*(i14))', int(nbits(:) / 1.1102230246251565E-16)

close(io_unit)

! Absolute error of x^3 - a = 0
open(file='err_root_abs.txt', newunit=io_unit)
do i = 1,n
  q = cuberoot_newton_quad(real(x(i), kind=real128))
  !print '(*(es24.16,2x))', &
  write(io_unit, '(*(es24.16,2x))') &
    x(i), &
    real(q, kind=real64), &

    abs((x(i)**(1./3.))**3 - x(i)), &
    abs(cuberoot_newton(x(i))**3 - x(i)), &
    abs(cuberoot_newton_nodiv(x(i))**3 - x(i)), &
    abs(cuberoot_halley(x(i))**3 - x(i)), &
    abs(cuberoot_halley_nodiv(x(i))**3 - x(i)), &
    abs(cuberoot_lagny(x(i))**3 - x(i)), &
    abs(cuberoot_lagny(x(i))**3 - x(i))
enddo
close(io_unit)

! Relative error of x^3 - a = 0
open(file='err_root_rel.txt', newunit=io_unit)
do i = 1,n
  q = cuberoot_newton_quad(real(x(i), kind=real128))
  !print '(*(es24.16,2x))', &
  write(io_unit, '(*(es24.16,2x))') &
    x(i), &
    real(q, kind=real64), &

    ! Relative error of x^3 - a = 0
    abs((x(i)**(1./3.))**3 - x(i))/x(i), &
    abs(cuberoot_newton(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_newton_nodiv(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_halley(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_halley_nodiv(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_lagny(x(i))**3 - x(i))/x(i), &
    abs(cbrt_ac(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_final(x(i))**3 - x(i))/x(i)
enddo
close(io_unit)

end
