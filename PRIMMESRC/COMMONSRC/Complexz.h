/* Slight modification of the SuperLU header file
 * to avoid conflicts with f2c and g2c libraries 
 * A. Stathopoulos, Sept 15, 2006                     */
/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#ifndef COMPLEX_H
#define COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { 
        double r, i; 
} Complex_Z;

/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add_primme(c, a, b) { (c).r = (a).r + (b).r; \
                         (c).i = (a).i + (b).i; }

/* Complex Subtraction c = a - b */
#define z_sub_primme(c, a, b) { (c).r = (a).r - (b).r; \
                         (c).i = (a).i - (b).i; }

/* Complex-Double Multiplication */
#define zd_mult_primme(c, a, b) { (c).r = (a).r * (b); \
                           (c).i = (a).i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult_primme(c, a, b) { \
        double cr, ci; \
        cr = (a).r * (b).r - (a).i * (b).i; \
        ci = (a).i * (b).r + (a).r * (b).i; \
        (c).r = cr; \
        (c).i = ci; \
    }

/* Complex conjugate-Complex Multiplication */
#define zconjz_mult_primme(c, a, b) { \
        double cr, ci; \
        cr = (a).r * (b).r + (a).i * (b).i; \
        ci = (a).r * (b).i - (a).i * (b).r; \
        (c).r = cr; \
        (c).i = ci; \
    }

/* Complex equality testing */
#define z_eq_primme(a, b)  ( (a).r == (b).r && (a).i == (b).i )


/* Prototypes for functions in dcomplex.c */
void   z_div_primme(Complex_Z *, Complex_Z *, Complex_Z *);
double z_abs_primme(Complex_Z);     /* exact sqrt(r^2+i^2) */
double z_abs1_primme(Complex_Z);    /* approximate  |r|+|i| */
void   z_exp_primme(Complex_Z *, Complex_Z *);
void   d_cnjg_primme(Complex_Z *r, Complex_Z *z);
double d_imag_primme(Complex_Z *);

#ifdef __cplusplus
}
#endif

#endif  /* COMPLEX_H */
