module cubes

use iso_fortran_env, only : real64, real128
implicit none

! Hallberg proposes starting with (3/8)^(1/3) for Newton iteration.
! This equates the relative error of the first estimate x1 for a=0.125 and 1.
! That is, (x1 - 1/2) / x1 = (x1' - 1) / x1'.
! This is less relevant to Halley solvers, but 0.7 seems a good choice overall.
real, parameter :: s = (3./8.)**(1./3.)

! This guess represents an upper bound, albeit biased towards higher values.
!real, parameter :: s = 1.

contains

! Compute the cube root in quadrature precision
! ... but only to double precision (as print* shows, at least).
elemental function cuberoot_newton_quad(x) result(root)
  real(kind=real128), intent(in) :: x
  real(kind=real128) :: root

  real(kind=real128) :: r
  integer :: i

  if (x == 0.) then
    root = 0.
  else
    root = 1._real128
    do i = 1, 6
      r = root
      root = ((2._real128)*r**3 + x) / ((3._real128)*r**2)
    enddo
  endif
end function cuberoot_newton_quad


! Six iteration Newton iteration solver for r**3 - x = 0
elemental function cuberoot_newton(x) result(r)
  real, intent(in) :: x
  real :: r
  !real, parameter :: s = (3./8.)**(1./3.)

  integer :: i

  if (x == 0.) then
    r = x
  else
    ! Implicitly initialize with r = s, followed one iteration.
    r = (2.*(s**3) + x) / (3.*(s**2))

    ! Do not simplify!  The form r = r - f/f' minimizes noise around the ULP.
    do i = 1, 5
      r = r - (r**3 - x) / (3.*(r**2))
    enddo
  endif
end function cuberoot_newton


! Three Halley + one Newton iterative solver for r**3 - x = 0
elemental function cuberoot_halley(x) result(r)
  real, intent(in) :: x
  real :: r

  integer :: i

  if (x == 0.) then
    r = 0.
  else
    ! Implicit initialization of r = s followed by one Halley iteration
    r = s * (s**3 + 2. * x) / (2.*(s**3) + x)

    ! This simplified form is faster than r = r - 2 f f' / (2f'*f' - f f'')
    do i = 1, 2
      r = r * (r**3 + 2.*x) / (2.*(r**3) + x)
    enddo
  endif

  ! Finalize with Newton iteration to minimize ULP noise.
  r = r - (r**3 - x) / (3.*(r**2))
end function cuberoot_halley


!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot_newton_nodiv(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  integer :: ex_3 ! One third of the exponent part of x, used to rescale x to get a.
  integer :: itt

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    ex_3 = ceiling(exponent(x) / 3.)
    ! Here asx is in the range of 0.125 <= asx < 1.0
    asx = scale(abs(x), -3*ex_3)
    !asx = x

    ! This first estimate is one iteration of Newton's method with a starting guess of s.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near s.
    num = 2. * s**3 + asx
    den = 3. * s**2

    ! Iteratively determine Root = asx**1/3 using Newton's method, noting that in this case Newton's
    ! method converges monotonically from above and needs no bounding.  For the range of asx from
    ! 0.125 to 1.0 with the first guess used above, 6 iterations suffice to converge to roundoff.

    do itt=1,4
      ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
      ! equivalently as Root = (2.0*Root**2 + asx) / (3.0 * Root**2).
      ! Keeping the estimates in a fractional form Root = num / den allows this calculation with
      ! fewer (or no) real divisions during the iterations before doing a single real division
      ! at the end, and it is therefore more computationally efficient.

      num_prev = num ; den_prev = den
      num = 2.0 * (num_prev**3) + asx * (den_prev**3)
      den = 3.0 * (den_prev * (num_prev**2))
    enddo

    root = num / den

    ! Finalize with Newton
    root = root - (root**3 - asx) / (3. * (root**2))

    root = sign(scale(root, ex_3), x)
  endif
end function cuberoot_newton_nodiv


!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot_halley_nodiv(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  integer :: ex_3 ! One third of the exponent part of x, used to rescale x to get a.
  integer :: itt

  real :: r

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    !!ex_3 = ceiling(exponent(x) / 3.)
    !!! Here asx is in the range of 0.125 <= asx < 1.0
    !!asx = scale(abs(x), -3*ex_3)
    asx = x

    ! Iteratively determine Root = asx**1/3 using Halley's method and then Newton's method, noting
    ! that in this case Newton's method and Halley's menthod both converge monotonically from above
    ! and need no bounding.

    !   This first estimate is one iteration of Halley's method with a starting guess of 1.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near 1.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.
    num = s * (s**3 + 2.*x)
    den = 2. * (s**3) + x

    do itt = 1, 2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num
      den_prev = den

      num = num_prev * (num_prev**3 + 2. * asx * (den_prev**3))
      den = den_prev * (2. * (num_prev**3) + asx * (den_prev**3))
    enddo
    root = num / den

    ! Finalize with complete form
    root = root - (root**3 - asx) / (3. * (root**2))

    !root = sign(scale(root_asx, ex_3), x)
  endif

end function cuberoot_halley_nodiv


!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot_final(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  integer :: ex_3 ! One third of the exponent part of x, used to rescale x to get a.
  integer :: itt

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    !ex_3 = ceiling(exponent(x) / 3.)
    !! Here asx is in the range of 0.125 <= asx < 1.0
    !asx = scale(abs(x), -3*ex_3)
    asx = x

    ! Iteratively determine Root = asx**1/3 using Halley's method and then Newton's method, noting
    ! that in this case Newton's method and Halley's menthod both converge monotonically from above
    ! and need no bounding.  Halley's method is slightly mroe complicated that Newton's method, but
    ! converges in a third fewer iterations.  The combiation used here saves 2 iterations over just
    ! using Newton's method.

    !   This first estimate is one iteration of Halley's method with a starting guess of 1.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near 1.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.
    !num = 0.5 + asx
    !den = 1. + 0.5 * asx

    ! Initialize with 0.7071...
    !num = 0.707106
    !den = 1.0

    ! Explicitly apply the first step
    num = 0.707106 * (0.707106**3 + 2. * asx)
    den = 2. * (0.707106**3) + asx

    ! Equivalent to:  root_asx = (1.0 + 2.0*asx) / (2.0 + asx)

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den
      num = num_prev * (num_prev**3 + 2. * asx * (den_prev**3))
      den = den_prev * (2. * (num_prev**3) + asx * (den_prev**3))
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo

    ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
    ! equivalently as Root = (2.0*Root**3 + asx) / (3.0 * Root**2).

    ! For asx in the range of 0.125 to 1., this iteration of Newton's method gives answers that
    ! are close to being within roundoff of the true solution.
    !root_asx = (2. * (num**3) + asx * (den**3)) / (3. * (den * (num**2)))
    ! Equivalent to:  root_asx = (2.0*root_asx**3 + asx) / (3.0*root_asx**2)

    ! Or simply just divide as in cuberoot_halley_nodiv
    root_asx = num / den

    ! One final iteration of Newton's method with the tradional correction form polishes
    ! up the root and gives a solution that is within the last bit of the true solution.
    root_asx = root_asx - (root_asx**3 - asx) / (3. * (root_asx**2))

    !root = sign(scale(root_asx, ex_3), x)
    root = root_asx
  endif
end function cuberoot_final

end
