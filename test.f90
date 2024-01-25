use cubes
implicit none
real(8) :: y_e(2), y(2), x
integer :: i

y = [0.125, 1.]
y_e = [0.5, 1.]

do i = 1,2
  x = y(i)
  print *, x, &
      x**(1._8/3._8) - y_e(i), &
      cuberoot_newton(x) - y_e(i), &
      cuberoot_newton_nodiv(x) - y_e(i), &
      cuberoot_halley(x) - y_e(i), &
      cuberoot_halley_nodiv(x) - y_e(i), &
      cuberoot_final(x) - y_e(i)
end do

end
