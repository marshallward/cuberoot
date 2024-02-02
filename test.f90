use cubes
implicit none
real(8) :: y_e(13), y(13), x
integer :: i

y = [(2.**(3*i), i=-6,5), 1.0 - 0.5*epsilon(0.5)]
y_e = [(2.**i, i=-6,5), 1.0]

print '(*(a24))', "x", "Quad Newton", "x**(1./3.)", "Double Newton", "Newton No-div"
do i = 1,size(y)
  x = y(i)
  if (i == size(y)) print *, "y = 1 - 0.5 ULP"
  print '(*(es24.16))', &
      x, &
      !x**(1._8/3._8) - y_e(i), &
      !cuberoot_newton(x) - y_e(i), &
      !cuberoot_newton_nodiv(x) - y_e(i), &
      !cuberoot_halley(x) - y_e(i), &
      !cuberoot_halley_nodiv(x) - y_e(i), &
      !cuberoot_final(x) - y_e(i)
      real(cuberoot_newton_quad(real(x,16)),8), &
      x**(1._8/3._8), &
      cuberoot_newton(x), &
      cuberoot_newton_nodiv(x), &
      cuberoot_halley(x), &
      cuberoot_halley_nodiv(x), &
      cuberoot_lagny(x)
end do

end
