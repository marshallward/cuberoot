module cubes

use, intrinsic :: iso_fortran_env, only : int64
use, intrinsic :: iso_fortran_env, only : real64, real128

use, intrinsic :: iso_c_binding, only : c_double

use, intrinsic :: ieee_arithmetic, only : ieee_fma

implicit none

! Hallberg proposes starting with (3/8)^(1/3) for Newton iteration.
! This equates the relative error of the first estimate x1 for a=0.125 and 1.
! That is, (x1 - 1/2) / x1 = (x1' - 1) / x1'.
! This is less relevant to Halley solvers, but 0.7 seems a good choice overall.
real, parameter :: s = (3./8.)**(1./3.)

! This guess represents an upper bound, albeit biased towards higher values.
!real, parameter :: s = 1.

! Floating point model, if bit layout from high to low is (sign, exp, frac)

integer, parameter :: bias = maxexponent(1.) - 1
  !< The double precision exponent offset
integer, parameter :: signbit = storage_size(1.) - 1
  !< Position of sign bit
integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
  !< Bit size of exponent
integer, parameter :: expbit = signbit - explen
  !< Position of lowest exponent bit
integer, parameter :: fraclen = expbit
  !< Length of fractional part

interface
  pure function cbrt_ac_c(a) result(r) bind(c, name="cbrt_ac")
    import :: c_double

    real(kind=c_double), value, intent(in) :: a
    real(kind=c_double) :: r
  end function cbrt_ac_c
end interface

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
    do i = 1, 20
      r = root
      root = ((2._real128)*r**3 + x) / ((3._real128)*r**2)
    enddo
  endif
end function cuberoot_newton_quad


!> Rescale `a` to the range [0.125, 1) and compute its cube-root exponent.
pure subroutine rescale_cbrt(a, x, e_r, s_a)
  real, intent(in) :: a
    !< The real parameter to be rescaled for cube root
  real, intent(out) :: x
    !< The rescaled value of a
  integer(kind=int64), intent(out) :: e_r
    !< Cube root of the exponent of the rescaling of `a`
  integer(kind=int64), intent(out) :: s_a
    !< The sign bit of a

  integer(kind=int64) :: xb
    ! Floating point value of a, bit-packed as an integer
  integer(kind=int64) :: e_a
    ! Unscaled exponent of a
  integer(kind=int64) :: e_x
    ! Exponent of x
  integer(kind=int64) :: e_div, e_mod
    ! Quotient and remainder of e in e = 3*(e/3) + modulo(e,3).

  ! Pack bits of a into xb and extract its exponent and sign.
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen) - bias
  ! Compute terms of exponent decomposition e = 3*(e/3) + modulo(e,3).
  ! (Fortran division is round-to-zero, so we must emulate floor division.)
  e_mod = modulo(e_a, 3_int64)
  e_div = (e_a - e_mod)/3

  ! Our scaling decomposes e_a into e = {3*(e/3) + 3} + {modulo(e,3) - 3}.

  ! The first term is a perfect cube, whose cube root is computed below.
  e_r = e_div + 1

  ! The second term ensures that x is shifted to [0.125, 1).
  e_x = e_mod - 3

  ! Insert the new 11-bit exponent into xb and write to x and extend the
  ! bitcount to 12, so that the sign bit is zero and x is always positive.
  call mvbits(e_x + bias, 0, explen + 1, xb, fraclen)
  x = transfer(xb, 1.)
end subroutine rescale_cbrt


!> Undo the rescaling of a real number back to its original base.
pure function descale(x, e_a, s_a) result(a)
  real, intent(in) :: x
    !< The rescaled value which is to be restored.
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the unscaled value
  real :: a
    !< Restored value with the corrected exponent and sign

  integer(kind=int64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=int64) :: e_x
    ! Biased exponent of x

  ! Apply the corrected exponent and sign to x.
  xb = transfer(x, 1_8)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1.)
end function descale


! Six iteration Newton iteration solver for r**3 - x = 0
elemental function cuberoot_newton(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a, xb
  integer :: i

  if (a == 0.) then
    r = a
  else
    ! Rescale to 0.125 < x < 1
    call rescale_cbrt(a, x, e_a, s_a)

    ! Implicitly initialize with r = s, followed one iteration.
    r = (2.*(s**3) + x) / (3.*(s**2))

    ! Do not simplify!  The form r = r - f/f' minimizes noise around the ULP.
    do i = 1, 5
      r = r - (r**3 - x) / (3.*(r**2))
    enddo

    ! Unscale and apply cuberoot to exponent
    r = descale(r, e_a, s_a)
  endif
end function cuberoot_newton


! Three Halley + one Newton iterative solver for r**3 - x = 0
elemental function cuberoot_halley(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a, xb
  integer :: i

  if (a == 0.) then
    r = 0.
  else
    ! Rescale to 0.125 < x < 1
    call rescale_cbrt(a, x, e_a, s_a)

    ! Implicit initialization of r = s followed by one Halley iteration
    r = s * (s**3 + 2. * x) / (2.*(s**3) + x)

    ! Round the result to 17-bits (Principia says 17, but is it not 16 bits?)
    ! (This is really important btw!!  Not sure why yet...)
    xb = transfer(r, 1_int64)
    xb = iand(xb, z'fffffff000000000')
    r = transfer(xb, 1._8)

    ! This simplified form is faster than r = r - 2 f f' / (2f'*f' - f f'')
    do i = 1, 2
      r = r * (r**3 + 2.*x) / (2.*(r**3) + x)
    enddo

    ! Finalize with Newton iteration to minimize ULP noise.
    r = r - (r**3 - x) / (3.*(r**2))

    ! Unscale and apply cuberoot to exponent
    r = descale(r, e_a, s_a)
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
    call rescale_cbrt(x, asx, e_s, a_s)

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
    root = descale(root, e_s, a_s)
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
    call rescale_cbrt(x, asx, e_s, a_s)

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
    root = descale(root, e_s, a_s)
  endif
end function cuberoot_halley_nodiv


! Used in the Principia project, but so far seems no better than the others.
! I may be missing something.
elemental function cuberoot_lagny_old(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a
  integer :: i

  if (a == 0.) then
    r = a
  else
    call rescale_cbrt(a, x, e_a, s_a)

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

    r = descale(r, e_a, s_a)
  endif
end function cuberoot_lagny_old


elemental function cuberoot_lagny(a) result(r)
  real, intent(in) :: a
  real :: r

  real, parameter :: G = 0.10096781215580288

  ! The actual parameters in the formula; not used in Principia
  real, parameter :: la = sqrt(3.)
  real, parameter :: lb = 4.
  real, parameter :: lc = 1./sqrt(12.)

  !! Steve Canon's optimized corrections?
  !! These do not seem to even converge... am I using them incorrectly?
  !real, parameter :: la = 1.7329127602084620
  !real, parameter :: lb = 4.0029873779316976
  !real, parameter :: lc = 0.23020325578115797

  real :: x
  integer(int64) :: e_a, s_a
  real :: r2, r3, r4, x2, x3
  real :: num, den, dr
  integer :: i

  integer(int64) :: xb

  ! Kahan initialization
  integer(int64) :: Yb, Cb, Qb

  ! Bit correction
  real :: r0, r1, rt
  logical :: may_correct, do_correct
  real :: cb_a, cb_b
  ! From 0x1.7C73DBBD9FA60p-66
  real, parameter :: tau = 2.0140991445038357E-20

  if (a == 0.) then
    r = a
  else
    call rescale_cbrt(a, x, e_a, s_a)

    x2 = x*x
    x3 = x2*x

    ! Use s=0.7 Halley estimate for now
    r = (2.*(s**3) + x) / (3.*(s**2))

    !!! TODO: Use Kahan estimate?  Is this OK?

    ! So AFAIK, the Kahan estimate is about the [eta,8*eta) range, and the
    ! optimal range is something like 0.1009...  But we do [0.125,1) so
    ! I am probably overthinking the whole thing.

    ! The follow turn out to be *TERRIBLE*.  Would like to see how Kahan or
    ! Canon stating values fare...
    !r = (2. + x) / 3.
    !r = 1

    ! 2. Lagny's irrational method
    r2 = r * r
    r4 = r2 * r2
    r = (la * r2 + sqrt(lb * x * r - r4)) * (lc / r)

    !! This can replace the Lagny iteration.  It is more accurate, but a bit
    !! slower.  The same effect can be obtained by 2 Lagny irrational
    !! iterations.  Mostly this fixed up problem near x=1.
    !do i=1,2
    !r = r * (r**3 + 2.*x) / (2.*(r**3) + x)
    !enddo

    ! Point being: There may be some opportunity to speedup the first stage.

    ! 3. Round the result to "17 bits" (Principia says 17, but looks like 16?)
    ! (This is really important btw!!  Not sure why yet...)
    xb = transfer(r, 1_int64)
    xb = iand(xb, z'fffffff000000000')
    r = transfer(xb, 1._real64)

    ! Fifth order Lagny-Schröder
    r2 = r*r
    r3 = r2*r
    num = (r3 - x)*((((10.*r3) + 16.*x)*r3) + x2)
    den = r2*((((15.*r3) + (51.*x))*r3) + (15*(x2)))

    dr = num/den

    r = r - dr

    !! This is supposed to fix the final bit, but currently does not work.

    !! Correct the bit?
    !r0 = r - dr
    !r1 = r - r0 - dr
    !rt = r0 + 2*r1

    !! (First check... dont really get that yet)
    !may_correct = (abs(0.5*(rt - r0) - r1) <= tau * r0) .and. (rt /= r0)

    !! (Then fix the bit if needed)
    !if (may_correct) then
    !  cb_a = min(r0, rt)
    !  cb_b = 0.5 * abs(r0 - rt)
    !  ! Replace with CbrtOneBit...
    !  do_correct = cbrt_one_bit(x, cb_a, cb_b)
    !  if (do_correct) then
    !    r = max(r0, rt)
    !  else
    !    r = min(r0, rt)
    !  endif
    !else
    !  r = r0
    !endif

    r = descale(r, e_a, s_a)
  endif
end function cuberoot_lagny

pure function cbrt_one_bit(x, a, b) result(do_correct)
  real, intent(in) :: x
  real, intent(in) :: a
  real, intent(in) :: b
  logical :: do_correct

  ! This is supposed to be a rather lengthy function which sweeps through and
  ! does bit-for-bit calculation up to the 54th bit.  But I am nowhere near
  ! ready to try and implement something like this.

  do_correct = .false.
end function cbrt_one_bit

!----

elemental function cbrt_ac(a) result(r)
  real, intent(in) :: a
  real :: r

  real(kind=c_double) :: a_c
  real(kind=c_double) :: r_c

  a_c = real(a, kind=c_double)
  r_c = cbrt_ac_c(a_c)
  r = real(r_c, kind=real64)
end function cbrt_ac


! Carbon copy of cbrt_ac (I hope...)
!
! This is literally just a Fortran implementation of this amazing solver which
! was posted on Stack Overflow.  (Or at least that's the plan.)
!
! Someday I hope this semi-anonymous person receives the proper accolades, but
! for now please see the citations below.
!
! URL:  https://stackoverflow.com/a/73354137/317172
! User: https://stackoverflow.com/users/2439725/wim
!
elemental function cuberoot_final(a) result(r)
  real, intent(in) :: a
  real :: r

  real :: x
  integer(int64) :: e_a, s_a, xb
  integer :: i

  real :: q
  real :: xqq, q2, r2, r2_h, r2_l
  real :: num, den

  if (a == 0.) then
    r = a
  else
    ! cbrt_ac does not do a full rescale, just one 0.125 rescale.  How does it
    ! get away with this??

    ! Rescale to 0.125 < x < 1
    call rescale_cbrt(a, x, e_a, s_a)
    !x = a

    ! Solve for a**(-2./3.), then transform to a**(1./3.)

    ! Initialize with fancy integer tricks.
    ! TODO: Needs explanation
    xb = transfer(x, 1_int64)
    xb = int(z'6A8EB53800000000', int64) - 2 * (xb / 3)
    q = transfer(xb, 1._real64)

    ! Newton iteration of f(q) = q**2 - x**(-3)
    xqq = (x * q) * q
    q = q + (1./3.) * (xqq * (-xqq) + q)

    ! Repeat, with an altered coefficient (?)
    xqq = (x * q) * q
    q = q + 0.33523333333 * (xqq * (-xqq) + q)

    ! Repeat
    xqq = (x * q) * q
    q = q + (1./3.) * (xqq * (-xqq) + q)

    ! Now we use the q = x**(-2/3) estimate to solve for x**(1/3).
    ! The product r is our estimate of x**(1/3).
    r = q * x

    ! Solve for f(r) = r - x**3, again using Newton iteration, but approximate
    ! f(r)/f'(r) = (x - r**3)/(3*r2) as q*(1/3)*(x - r**3).
    ! (How?  r ~= x**(1/3) and q ~= x**(-2/3), so 1/r**2 ~= q.)
    ! (Why?  I think it avoids a division...?)

    ! We want to maximize the residual in x - r**3, so split r**2 = r2_l + r2_h
    ! to preserve the FMA residual if available.
    r2_h = r * r

    ! r2_l = (r * r) - r2_h
    ! Every compiler wants to substitute r*r with r2_h so r2_l needs an
    ! explicit FMA operation.
    r2_l = ieee_fma(r, r, -r2_h)

    ! The rest of the FMAs complete without any trouble.
    r = r + q * (1./3.) * (-r * r2_l + (-r * r2_h + x))

    ! Author says r approximates x**(1/3) within 0.50002 ULP.
    ! Time to finish the job!

    !! Apply one final Halley iteration to reduce it below 0.5 ULP.
    !r2_h = r*r
    !r2_l = r*r - r2_h

    !num = -r * r2_l + (-r * r2_h + a)
    !den = 3. * a - 2.*num
    !r = r + r * num / den

    ! Unscale and apply cuberoot to exponent
    r = descale(r, e_a, s_a)
  endif
end function cuberoot_final

end
