Cube root solvers
=================

This repo contains several implementations of cube root solvers in Fortran.

Requirements
------------

This first requirement is vital, and the primary motivation for this
investigation.  The others are desirable, but can be sacrified if necessary.

1. Solutions must be power-of-two dimensionally invariant.

   That is, ``2^(-N) * cuberoot(2**(3*N) x) = x`` without bit roundoff.

   This permits us to continue our dimensional testing in expressions with cube
   root (and is the primary motivation of this work.)

   This is addressed by scaling all solutions to a value between 1/8 and 1,
   then unscaling the result.

2. Certain values should be mathematically exact.

   We would *really* like ``cuberoot(0.125) = 0.5`` and ``cuberoot(1.0) =
   1.0``.  Probably others as well, but these are rather important, partly in
   support of dimensional scaling.

3. Solutions should be accurate within one ULP.

   We do not expect mathematical equivalence (for now), but do expect accuracy
   comparable to ``x**(1./3.)``, and ideally within bit roundoff of
   quadratic precision ("real(16)").

4. Solutions must be independent of external dependencies, such as math
   libraries.

   ``x**(1./3.)`` is not an option because the ``**`` operator is ambiguous and
   may depend on an external math library, such as the C standard library's
   ``libm``.  Our solution must be implemented in primitive mathematical
   operations.

5. Solutions must be bit-reproducible across environments, compilers, and
   platforms.

   This is achieved by enforcing an order of operations, and not relying on
   other transcendentals, such as ``exp`` and ``log``.  IEEE-754 requires
   division to be bit-reproducible, but vendors may not necessarily comply.

6. Solutions should be as fast or faster than the intrinsic ``**3`` operator.

   Our solutions are faster than GNU ``libm``, but slower than Intel's
   ``cbrt``.  We may be able to improve this in the future.


Methods
=======

All methods here currently use some iterative solver of ``x**3 - a = 0``.

(TODO: Proper explanations; for now, look at the code.)

* Intrinsic: ``x**(1./3.)``

* Newton: ``r = r - f /f'`` (x6)

* Halley ``r = r - 2 f f' / (f'**2 - f f'')`` (x3 + 1 Newton)

* "No-division" Newton: Numerator and denominator are computed separately,
  with a final division at the end.  (4x no-div, 1 regular Newton)

* "No-division Halley: Same, with Halley's method (3x Halley, 1x std Newton)

There is also a "Final Halley" test which is a placeholder for testing.


Results
=======

Errors
------

TODO


Timings
-------

TODO


Comments
========

The overall structure is the same for all of these solutions.

1. Select an appropriate initial value.  ``r = 0.7`` appears to be a strong
   choice for all solvers.  There is some mathematical motivation for
   refinement, with varying impact.

2. Pre-compute the first iteration.  For example, don't do ``root = 1.``.
   Instead, do this for Newton iteration::

      root = (2. + x) / 3.

   This can be about 20% faster than explicitly computing the first iteration.

3. Regardless of method, finalize with an unsimplified Newton iteration::

      root = root - (root**3 - x) / (3. * (root**2))

   Something about the "root = root + correct" form cleans up the final few
   bits.  (Needs a proper mathematical explanation... TODO?)

4. Take care with exponents in order of operations.  For truly baffling
   reasons, Intel Fortran does not regard the exponent as an independent
   operator with highest priority in expressions like this::

      a * b**3

   and will happily compute it as::

      (a * b) * b**2

   if it sees some trivial optimization (even if it is no better than ``a *
   (b * b**2)``.

   For bit equivalence across compilers, wrap all exponents in parentheses.::

      a * (b**3)

* Loop iterations were selected to produce double precision results.  These
  would need to be tuned for single or quadruple precision.

* The "no-division" methods should be used with some care.  While the ratio
  will converge, there is nothing constraining the magnitudes of these values,
  and they may grow beyond the limits of double precision.  This is not a
  problem within three or four iterations, but often needs to "renormalize"
  after about six iterations.  For example, one can do something like::

      num = num / den ; den = 1

* Earlier versions of these functions included various convergence tests in
  order to avoid redundant iterations.  But it was found that the tests
  themselves exceeded the cost of checking for convergence, and it was faster
  to simply run a fixed number of iterations.

* Starting from ``root = 1`` is viable, but it can produce inflated errors near
  0.125.  It also causes ``cuberoot(0.125)`` to have an error in 1 ULP.
  Starting near 0.7 seems to fix this issue and significally reduce most errors
  around 0.125.


Summary
=======

All methods seem capable of achieving the required goals.  Every method is
accurate and competitively fast.  There is no "wrong" choice.

The fastest method was the "no-division" Halley method with a final Newton
iteration.

None of the methods were able to
exactly produce results from quadratic precision, but all were equivalent
within one ULP.  (Note that ``x**(1./3.)`` was also only exact within one ULP).
