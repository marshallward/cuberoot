use iso_fortran_env, only : real64, real128
use cubes
implicit none

integer :: k, n
real(kind=real128) :: q
real, allocatable :: x(:)
integer :: i

!print *, "Resolution (as k in 2**-k)?"
!read (*,*) k
k = 17
n = 1 + 7 * 2**(k-3)
allocate(x(n))

x = [(i*2.**-k, i=2**(k-3), 2**k)]

!print '(*(a26))', "x", "Quad Newton", "x**(1./3.)", "Double Newton", "Newton No-div"
do i = 1,n
  q = cuberoot_newton_quad(real(x(i), kind=real128))
  print '(*(es24.16,2x))', &
    x(i), &
    real(q, kind=real64), &

    ! Absolute error of x**(1/3) relative to quad
    real(abs(x(i)**(1./3.) - q), kind=real64), &
    real(abs(cuberoot_newton(x(i)) - q), kind=real64), &
    real(abs(cuberoot_newton_nodiv(x(i)) - q), kind=real64), &
    real(abs(cuberoot_halley(x(i)) - q), kind=real64), &
    real(abs(cuberoot_halley_nodiv(x(i)) - q), kind=real64), &
    real(abs(cuberoot_final(x(i)) - q), kind=real64)

    ! Absolute error of x^3 - a = 0
    !abs((x(i)**(1./3.))**3 - x(i)), &
    !abs(cuberoot_newton(x(i))**3 - x(i)), &
    !abs(cuberoot_newton_nodiv(x(i))**3 - x(i)), &
    !abs(cuberoot_halley(x(i))**3 - x(i)), &
    !abs(cuberoot_halley_nodiv(x(i))**3 - x(i)), &
    !abs(cuberoot_final(x(i))**3 - x(i))

    ! Relative error of x^3 - a = 0
    !abs((x(i)**(1./3.))**3 - x(i))/x(i), &
    !abs(cuberoot_newton(x(i))**3 - x(i))/x(i), &
    !abs(cuberoot_newton_nodiv(x(i))**3 - x(i))/x(i), &
    !abs(cuberoot_halley(x(i))**3 - x(i))/x(i), &
    !abs(cuberoot_halley_nodiv(x(i))**3 - x(i))/x(i), &
    !abs(cuberoot_final(x(i))**3 - x(i))/x(i)
enddo

end
