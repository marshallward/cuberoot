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
    x(i)**(1./3.), &
    !cuberoot_newton(x(i)), &
    !abs(x(i)**(1./3.) - cuberoot_newton(x(i)))
    !cuberoot_newton_nodiv(x(i)), &
    cuberoot_halley(x(i)), &
    abs(x(i)**(1./3.) - cuberoot_halley(x(i)))
    !cuberoot_halley_nodiv(x(i)), &
    !cuberoot_lagny(x(i))
enddo

end
