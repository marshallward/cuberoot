/* Accurate double precision cube root function cbrt_ac() and surrounding test code */

/* Solver found at Stack Overflow
 * https://stackoverflow.com/a/73354137/317172
 * User: https://stackoverflow.com/users/2439725/wim
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <quadmath.h>

/* gcc -O3 -Wall -m64 -std=c99 -march=skylake  cbrt_tst.c -lm -lquadmath */
/* Assumptions:
                - rounding mode: round nearest 
                - safe math: with gcc don't use -ffast-math, with icc enable -fp-model=precise
                - ieee-754 floating point math, higher intermediate floating point precision should not be used in the calculations, i.e. FLT_EVAL_METHOD == 0, which is default with the gcc x86-64 compiler.
*/    

uint64_t d2i (double x);
double i2d (uint64_t i);

double cbrt_ac(double z){                            // Accurate cube root (double)
    double a, y, r, r2_h, r2_l, y_a2y4, ayy, diff, diff3, denom;
    uint64_t ai, ai23, aim23;
    int issmall;
    
    a = fabs(z);
    issmall = (a <  0.015625);                       // Scale large, small and/or subnormal numbers to avoid underflow, overflow or subnormal numbers
    a = issmall ? (a * 0x1.0p+210) : (a * 0.125);
    ai = d2i(a);
    if ((ai >= 0x7FF0000000000000ull) || (z == 0.0)){    // Inf, 0.0 and NaN
        r = z + z;
    }
    else
    { 
        ai23 = 2 * (ai/3);                           // Integer division. The compiler, with suitable optimization level, should generate a much more efficient multiplication by 0xAAAAAAAAAAAAAAAB
        aim23 = 0x6A8EB53800000000ull - ai23;        // This uses a similar idea as the "fast inverse square root" approximation, see https://en.wikipedia.org/wiki/Fast_inverse_square_root
        y = i2d(aim23);                              // y is an approximation of a^(-2/3)

        ayy = (a * y) * y;                           // First Newton iteration for f(y)=a^2-y^-3 to calculate a better approximation y=a^(-2/3)
        y_a2y4 = fma(ayy, -ayy, y);
        y = fma(y_a2y4, 0.33333333333333333333, y);

        ayy = (a * y) * y;                           // Second Newton iteration
        y_a2y4 = fma(ayy, -ayy, y);
        y = fma(y_a2y4, 0.33523333333, y);           // This is a small modification to the exact Newton parameter 1/3 which gives slightly better results

        ayy = (a * y) * y;                           // Third Newton iteration
        y_a2y4 = fma(ayy, -ayy, y);
        y = fma(y_a2y4, 0.33333333333333333333, y);

        r = y * a;                                   // Now r = y * a is an approximation of a^(1/3), because y approximates a^(-2/3).
        r2_h = r * r;                                // Compute one pseudo Newton step with g(r)=a-r^3, but instead of dividing by f'(r)=3r^2 we multiply with 
                                                     // the approximation 0.3333...*y (division is usually a relatively expensive operation)
        r2_l = fma(r, r, -r2_h);                     // For better accuracy we split r*r=r^2 as r^2=r2_h+r2_l exactly.
        diff = fma(r2_h, -r, a);                     // Compute diff=a-r^3 accurately: diff=(a-r*r2_h)-r*r2_l with two fma instructions
        diff = fma(r2_l, -r, diff);               
        diff3 = diff * 0.33333333333333333333; 
        r = fma(diff3, y, r);                        // Now r approximates a^(1/3) within about 0.50002 ulp
                                                     
        r2_h = r * r;                                // One final Halley iteration
        r2_l = fma(r, r, -r2_h);
        diff = fma(r2_h, -r, a);
        diff = fma(r2_l, -r, diff);
        denom = fma(a, 3.0, -2.0 * diff);     
        r = fma(diff/denom, r, r);

        r = issmall ? (r * 0x1.0p-70) : (r * 2.0);   // Undo scaling
        r = copysign(r, z);
    }
    return r;                                     
}

union dbl_uint64{
    double d;
    uint64_t i;
};

uint64_t d2i (double x){                 // Bitwise copy (type-punning) from double to uint64_t (no conversion)
    union dbl_uint64 tmp ;               // With C++ use memcopy instead of this "union trick"   
    tmp.d = x;
    return tmp.i;
}

double i2d (uint64_t i){                 // Bitwise copy (type-punning) from uint64_t to double (no conversion)
    union dbl_uint64 tmp ;               // With C++ use memcopy instead of this "union trick" 
    tmp.i = i;
    return tmp.d;
}



