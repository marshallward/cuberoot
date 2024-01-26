use cubes
implicit none
real(8) :: y_e(13), y(13), x
integer :: i

y = [(2.**(3*i), i=-6,6)]
y_e = [(2.**i, i=-6,6)]

do i = 1,size(y)
  x = y(i)
  print *, x, &
      !x**(1._8/3._8) - y_e(i), &
      !cuberoot_newton(x) - y_e(i), &
      !cuberoot_newton_nodiv(x) - y_e(i), &
      !cuberoot_halley(x) - y_e(i), &
      !cuberoot_halley_nodiv(x) - y_e(i), &
      !cuberoot_final(x) - y_e(i)
      x**(1._8/3._8), &
      cuberoot_newton(x), &
      cuberoot_newton_nodiv(x), &
      cuberoot_halley(x), &
      cuberoot_halley_nodiv(x)
      !cuberoot_final(x)
end do

end
