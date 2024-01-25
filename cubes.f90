module cubes

use iso_fortran_env, only : real64, real128
implicit none

contains

elemental function cuberoot_newton_quad(x) result(root)
  real(kind=real128), intent(in) :: x
  real(kind=real128) :: root

  real(kind=real128) :: r
  integer :: i

  if (x == 0.) then
    root = 0.
  else
    root = 1._real128
    do i=1,6
      r = root
      root = ((2._real128)*r**3 + x) / ((3._real128)*r**2)
    enddo
  endif
end function cuberoot_newton_quad


! Six iteration Newton iteration solver for r**3 - x = 0
elemental function cuberoot_newton(x) result(r)
  real, intent(in) :: x
  real :: r

  integer :: i

  if (x == 0.) then
    r = x
  else
    ! Implicitly initialize with r = 1, followed one iteration.
    ! This appears to be faster than explicit r = 1 with six iterations.
    r = (2. + x) / 3.

    ! Do not simplify!  The form r = r - f/f' minimizes noise around the ULP.
    do i = 1, 5
      r = r - (r**3 - x) / (3.*r**2)
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
    ! Implicit initialization of r = 1 followed by one Halley iteration
    r = (1. + 2.*x) / (2. + x)

    ! The simplified form is faster than r = r - f f'' / (2f'*f' - f f'')
    do i = 1, 2
      r = r * (r**3 + 2.*x) / (2.*r**3 + x)
    enddo
  endif

  ! Finalize with Newton iteration to minimize ULP noise.
  r = r - (r**3 - x) / (3.*r**2)
end function cuberoot_halley


!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot_nodiv(x) result(root)
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
  real, parameter :: den_min = 2.**(minexponent(1.) / 4 + 4)  ! A value of den that triggers rescaling [C]
  real, parameter :: den_max = 2.**(maxexponent(1.) / 4 - 2)  ! A value of den that triggers rescaling [C]
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

    ! This first estimate is one iteration of Newton's method with a starting guess of 1.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near 1.
    num = 2.0 + asx
    den = 3.0
    ! Iteratively determine Root = asx**1/3 using Newton's method, noting that in this case Newton's
    ! method converges monotonically from above and needs no bounding.  For the range of asx from
    ! 0.125 to 1.0 with the first guess used above, 6 iterations suffice to converge to roundoff.

    !do itt=1,9
    do itt=1,6
      ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
      ! equivalently as Root = (2.0*Root**2 + asx) / (3.0 * Root**2).
      ! Keeping the estimates in a fractional form Root = num / den allows this calculation with
      ! fewer (or no) real divisions during the iterations before doing a single real division
      ! at the end, and it is therefore more computationally efficient.

      num_prev = num ; den_prev = den
      num = 2.0 * num_prev**3 + asx * den_prev**3
      den = 3.0 * (den_prev * num_prev**2)

      !if ((num * den_prev == num_prev * den) .or. (itt == 9)) then
      !  !   If successive estimates of root are identical, this is a converged solution.
      !  root_asx = num / den
      !  exit
      !elseif (num * den_prev > num_prev * den) then
      !  !   If the estimates are increasing, this also indicates convergence, but for a more subtle
      !  ! reason.  Because Newton's method converges monotonically from above (at least for infinite
      !  ! precision math), the only reason why this estimate could increase is if the iterations
      !  ! have converged to a roundoff-level limit cycle around an irrational or otherwise
      !  ! unrepresentable solution, with values only changing in the last bit or two.  If so, we
      !  ! should stop iterating and accept the one of the current or previous solutions, both of
      !  ! which will be within numerical roundoff of the true solution.
      !  root_asx = num / den
      !  ! Pick the more accurate of the last two iterations.
      !  ! Given that both of the two previous iterations are within roundoff of the true
      !  ! solution, this next step might be overkill.
      !  if ( abs(den_prev**3*root_asx**3 - den_prev**3*asx) > abs(num_prev**3 - den_prev**3*asx) ) then
      !    ! The previous iteration was slightly more accurate, so use that for root_asx.
      !    root_asx = num_prev / den_prev
      !  endif
      !  exit
      !endif

      ! Because successive estimates of the numerator and denominator tend to be the cube of their
      ! predecessors, the numerator and denominator need to be rescaled by division when they get
      ! too large or small to avoid overflow or underflow in the convergence test below.
      if ((den > den_max) .or. (den < den_min)) then
        num = scale(num, -exponent(den))
        den = scale(den, -exponent(den))
      endif

    enddo

    !root = sign(scale(root_asx, ex_3), x)
    !root = root_asx

    root = num / den
  endif

end function cuberoot_nodiv


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
  real, parameter :: den_min = 2.**(minexponent(1.) / 4 + 4)  ! A value of den that triggers rescaling [C]
  real, parameter :: den_max = 2.**(maxexponent(1.) / 4 - 2)  ! A value of den that triggers rescaling [C]
  logical :: converged
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
    num = 1.0 + 2.0*asx
    den = 2.0 + asx
    !converged = .false.

    do itt=1,3
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den
      num = num_prev * (num_prev**3 + 2.0 * asx * den_prev**3)
      den = den_prev * (2.0 * num_prev**3 + asx * den_prev**3)

      !if (num * den_prev == num_prev * den) then
      !  converged = .true.
      !  root_asx = num / den
      !  exit
      !endif
    enddo

    !! For the range of asx from 0.125 to 1.0 with the first guess of 1.0 and 3 iterations with
    !! Halley's method, 2 more iterations with Newton's method suffice to converge within roundoff.
    !if (.not.converged) then ; do itt=1,4

    !  ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
    !  ! equivalently as Root = (2.0*Root**3 + asx) / (3.0 * Root**2).
    !  num_prev = num ; den_prev = den
    !  num = 2.0 * num_prev**3 + asx * den_prev**3
    !  den = 3.0 * (den_prev * num_prev**2)

    !  ! Because successive estimates of the numerator and denominator tend to be the cube of their
    !  ! predecessors, the numerator and denominator need to be rescaled when they get too large or
    !  ! small to avoid overflow or underflow in the convergence test below.
    !  if ((den > den_max) .or. (den < den_min)) then
    !    num = scale(num, -exponent(den))
    !    den = scale(den, -exponent(den))
    !  endif

    !  if ((num * den_prev == num_prev * den) .or. (num**3 == asx * den**3) .or. (itt == 4)) then
    !    !   If successive estimates of root are identical, this is a converged solution.
    !    root_asx = num / den
    !    exit
    !  elseif (num * den_prev > num_prev * den) then
    !    !   If the estimates are increasing, this also indicates convergence, but for a more subtle
    !    ! reason.  Because Newton's method converges monotonically from above (at least for infinite
    !    ! precision math), the only reason why this estimate could increase is if the iterations
    !    ! have converged to a roundoff-level limit cycle around an irrational or otherwise
    !    ! unrepresentable solution, with values only changing in the last bit or two.  If so, we
    !    ! should stop iterating and accept the one of the current or previous solutions, both of
    !    ! which will be within numerical roundoff of the true solution.
    !    root_asx = num / den
    !    ! Pick the more accurate of the last two iterations.
    !    ! Given that both of the two previous iterations are within roundoff of the true
    !    ! solution, this next step might be overkill.
    !    if ( abs(den_prev**3*root_asx**3 - den_prev**3*asx) > abs(num_prev**3 - den_prev**3*asx) ) then
    !      ! The previous iteration was slightly more accurate, so use that for root_asx.
    !      root_asx = num_prev / den_prev
    !    endif
    !    exit
    !  endif

    !enddo ; endif

    !root = sign(scale(root_asx, ex_3), x)
    root = num / den


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
    num = 0.5 + asx
    den = 1.0 + 0.5*asx
    ! Equivalent to:  root_asx = (1.0 + 2.0*asx) / (2.0 + asx)

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den
      num = num_prev * (num_prev**3 + 2.0 * asx * den_prev**3)
      den = den_prev * (2.0 * num_prev**3 + asx * den_prev**3)
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo
    ! At this point the error in the root is better than 1 part in 2e7.

    ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
    ! equivalently as Root = (2.0*Root**3 + asx) / (3.0 * Root**2).

    ! For asx in the range of 0.125 to 1., this iteration of Newton's method gives answers that
    ! are close to being within roundoff of the true solution.
    root_asx = (2.0 * num**3 + asx * den**3) / ( 3.0 * (den * num**2) )
    ! Equivalent to:  root_asx = (2.0*root_asx**3 + asx) / (3.0*root_asx**2)

    ! One final iteration of Newton's method with the tradional correction form polishes
    ! up the root and gives a solution that is within the last bit of the true solution.
    root_asx = root_asx - (root_asx**3 - asx) / (3.0 * root_asx**2)

    !root = sign(scale(root_asx, ex_3), x)
    root = root_asx
  endif
end function cuberoot_final

end
