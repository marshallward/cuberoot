module cubes

use iso_fortran_env, only : real64, real128
use iso_fortran_env, only : int64
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


!> Rescale `a` to the range [0.125, 1) while preserving its fractional term.
pure subroutine rescale(a, x, e_a, s_a)
  real(kind=real64), intent(in) :: a
    !< The real parameter to be rescaled.
  real(kind=real64), intent(out) :: x
    !< The rescaled value of `a`
  integer(kind=int64), intent(out) :: e_a
    !< The biased exponent of `a`
  integer(kind=int64), intent(out) :: s_a
    !< The sign bit of `a`

  integer, parameter :: bias = 1023
    !< The double precision exponent offset

  integer(kind=int64) :: xb
    !< A floating point number, bit-packed as an integer
  integer(kind=int64) :: e_scaled
    !< The new rescaled exponent of `a` (i.e. the exponent of `x`)

  ! Pack bits of `a` into `xb` and extract its exponent and sign
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, 63, 1)
  e_a = ibits(xb, 52, 11)

  ! Decompose the exponent as `e = modulo(e,3) + 3*(e/3)` and extract the
  ! rescaled exponent, shifted to to select the triple (-3,-2,-1).
  e_scaled = modulo(e_a, 3) - 3 + bias

  ! Insert the new 11-bit exponent into `xb`, while also setting the sign bit
  ! to zero, ensuring that `xb` is always positive.
  call mvbits(e_scaled, 0, 12, xb, 52)

  ! Transfer the final modified value to `x`
  x = transfer(xb, 1._real64)
end subroutine rescale


!> Descale a real number to its original base, and apply the cube root to the
!! remaining exponent.
pure function descale_cbrt(x, e_a, s_a) result(r)
  real(kind=real64), intent(in) :: x
    !< Cube root of the rescaled value, which was rescaled to [0.125, 1.0)
  integer(kind=real64), intent(in) :: e_a
    !< Exponent of the original value to be cube rooted
  integer(kind=real64), intent(in) :: s_a
    !< Sign bit of the original value to be cube rooted
  real(kind=real64) :: r
    !< Restored vale with the cube root applied to its exponent

  integer, parameter :: bias = 1023
    !< The double precision exponent offset

  integer(kind=real64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=real64) :: e_r
    ! Exponent of the descaled value

  ! Extract the exponent of the rescaled value, in {-3, -2, -1}
  xb = transfer(x, 1_8)
  e_r = ibits(xb, 52, 11)

  ! Apply the cube root to the old exponent (after removing its bias) and add
  ! to the rescaled exponent.  Correct the previous -3 with a +1.
  e_r = e_r + (e_a/3 - bias/3 + 1)

  ! Apply the corrected exponent and sign and convert back to double precision.
  call mvbits(e_r, 0, 11, xb, 52)
  call mvbits(s_a, 0, 1, xb, 63)
  r = transfer(xb, 1._8)
end function descale_cbrt


! Six iteration Newton iteration solver for r**3 - x = 0
elemental function cuberoot_newton(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a
  integer :: i

  if (a == 0.) then
    r = a
  else
    ! Rescale to 0.125 < x < 1
    call rescale(a, x, e_a, s_a)

    ! Implicitly initialize with r = s, followed one iteration.
    r = (2.*(s**3) + x) / (3.*(s**2))

    ! Do not simplify!  The form r = r - f/f' minimizes noise around the ULP.
    do i = 1, 5
      r = r - (r**3 - x) / (3.*(r**2))
    enddo

    ! Unscale and apply cuberoot to exponent
    r = descale_cbrt(r, e_a, s_a)
  endif
end function cuberoot_newton


! Three Halley + one Newton iterative solver for r**3 - x = 0
elemental function cuberoot_halley(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a
  integer :: i

  if (a == 0.) then
    r = 0.
  else
    ! Rescale to 0.125 < x < 1
    call rescale(a, x, e_a, s_a)

    ! Implicit initialization of r = s followed by one Halley iteration
    r = s * (s**3 + 2. * x) / (2.*(s**3) + x)

    ! This simplified form is faster than r = r - 2 f f' / (2f'*f' - f f'')
    do i = 1, 2
      r = r * (r**3 + 2.*x) / (2.*(r**3) + x)
    enddo

    ! Finalize with Newton iteration to minimize ULP noise.
    r = r - (r**3 - x) / (3.*(r**2))

    ! Unscale and apply cuberoot to exponent
    r = descale_cbrt(r, e_a, s_a)
  endif
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
  integer(8) :: xb, e, e_s, e_r, e_new
  integer(8) :: a_s

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale(x, asx, e_s, a_s)

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

    ! Descale exponent and take its cube root
    root = descale_cbrt(root, e_s, a_s)
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
  integer(8) :: xb, e, e_s, e_r, e_new
  integer(8) :: a_s

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale(x, asx, e_s, a_s)

    ! Iteratively determine Root = asx**1/3 using Halley's method and then Newton's method, noting
    ! that in this case Newton's method and Halley's menthod both converge monotonically from above
    ! and need no bounding.

    !   This first estimate is one iteration of Halley's method with a starting guess of 1.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near 1.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.
    num = s * (s**3 + 2.*asx)
    den = 2. * (s**3) + asx

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

    ! Descale exponent and take its cube root
    root = descale_cbrt(root, e_s, a_s)
  endif
end function cuberoot_halley_nodiv


! Used in the Principia project, but so far seems no better than the others.
! I may be missing something.
elemental function cuberoot_lagny(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a
  integer :: i

  if (a == 0.) then
    r = a
  else
    call rescale(a, x, e_a, s_a)

    ! Implicitly initialize with r = s, followed one iteration.
    r = (2.*(s**3) + x) / (3.*(s**2))

    ! NOTE: They promote Kahan initialization, which I have not seriously
    ! explored yet.  (I tried but have yet to see much benefit.)

    do i = 1, 3
      r = r + r*(x - r**3)/(2*r**3 + x)

      ! Other apparently more accurate implementations, but they're rather slow
      !r = (sqrt(3.) * r**2 + sqrt(4.*x*r - r**4))*sqrt(1./12.)/r
      !r = 0.5*r + sqrt(0.25*r**2 + (x-r**3)/(3.*r))
      !r = k*r + sqrt(l*r**2 + (x - r**3)/(m*r))
    enddo

    r = descale_cbrt(r, e_a, s_a)
  endif
end function cuberoot_lagny


end
