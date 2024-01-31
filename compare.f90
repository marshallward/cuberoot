use iso_fortran_env, only : real64, real128
use cubes
implicit none

integer :: k, n
real(kind=real128) :: q
real, allocatable :: x(:)
integer :: i

integer :: io_unit

!print *, "Resolution (as k in 2**-k)?"
!read (*,*) k
k = 17
n = 1 + 7 * 2**(k-3)
allocate(x(n))

x = [(i*2.**-k, i=2**(k-3), 2**k)]

!print '(*(a26))', "x", "Quad Newton", "x**(1./3.)", "Double Newton", "Newton No-div"

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
    real(abs(cuberoot_lagny(x(i)) - q), kind=real64)
enddo
close(io_unit)

! Absolute error of x^3 - a = 0
open(file='err_root_abs.txt', newunit=io_unit)
do i = 1,n
  q = cuberoot_newton_quad(real(x(i), kind=real128))
  !print '(*(es24.16,2x))', &
  write(io_unit, '(*(es24.16,2x))') &
    x(i), &
    real(q, kind=real64), &

    abs((x(i)**(1./3.))**3 - x(i))/x(i), &
    abs(cuberoot_newton(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_newton_nodiv(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_halley(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_halley_nodiv(x(i))**3 - x(i))/x(i), &
    abs(cuberoot_lagny(x(i))**3 - x(i))/x(i)
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
    abs(cuberoot_lagny(x(i))**3 - x(i))/x(i)
enddo
close(io_unit)

end
