use iso_fortran_env, only : real64, real128
use cubes
implicit none

!integer, parameter :: k = 22
!integer, parameter :: n = 1 + 7 * 2**(k-3)
integer :: k, n
real(kind=real128) :: q
!real :: x(n), r(n)
real, allocatable :: x(:), r(:)
integer :: i, j, m

real :: r1, r2

integer, parameter :: niter = 20
integer :: count_rate, count_max, c1, c2
real :: clock_rate

!read(*,*) k
k = 24
n = 1 + 7 * 2**(k-3)

allocate(x(n), r(n))

x = [(i*2.**-k, i=2**(k-3), 2**k)]

! Set up clock
call system_clock(count_rate=count_rate, count_max=count_max)
clock_rate = real(count_rate)

r = x**(1./3.)
call system_clock(count=c1)
do i = 1, niter
  r = x**(1./3.)
end do
call system_clock(count=c2)
print *, "x**(1./3.):", (c2 - c1) / clock_rate / niter

r = cuberoot_newton(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_newton(x)
end do
call system_clock(count=c2)
print *, "Newton:", (c2 - c1) / clock_rate / niter

r = cuberoot_newton_nodiv(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_newton_nodiv(x)
end do
call system_clock(count=c2)
print *, "Newton No_div:", (c2 - c1) / clock_rate / niter

r = cuberoot_halley(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_halley(x)
end do
call system_clock(count=c2)
print *, "Halley:", (c2 - c1) / clock_rate / niter

r = cuberoot_halley_nodiv(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_halley_nodiv(x)
end do
call system_clock(count=c2)
print *, "Halley no_div:", (c2 - c1) / clock_rate / niter

r = cuberoot_lagny(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_lagny(x)
end do
call system_clock(count=c2)
print *, "Lagny:", (c2 - c1) / clock_rate / niter

r = cbrt_ac(x)
call system_clock(count=c1)
do i = 1, niter
  r = cbrt_ac(x)
end do
call system_clock(count=c2)
print *, "AC:", (c2 - c1) / clock_rate / niter

r = cuberoot_final(x)
call system_clock(count=c1)
do i = 1, niter
  r = cuberoot_final(x)
end do
call system_clock(count=c2)
print *, "Final:", (c2 - c1) / clock_rate / niter

end
