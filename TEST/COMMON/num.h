
#ifndef NUM_H
#define NUM_H

#ifdef USE_DOUBLECOMPLEX
#  include "../../PRIMMESRC/ZSRC/numerical_z.h"
#  include <complex.h>
#  undef I
#  define IMAGINARY _Complex_I
#  define PRIMME_NUM complex double
#  define REAL_PART(x) (creal(x))
#  define REAL_PARTZ(x) ((x).r)
#  define CONJ(x) (conj(x))
#  define SUF(NAME) NAME ## _zprimme
#  define SUFX(NAME) Num_z ## NAME ## _zprimme
#  define PREFIX(NAME) z ## NAME
#  define COMPLEXZ(X) ((Complex_Z*)(X))
   static inline Complex_Z COMPLEXZV(PRIMME_NUM x) {Complex_Z z={creal(x), cimag(x)}; return z;}
#else
#  include "../../PRIMMESRC/DSRC/numerical_d.h"
#  define IMAGINARY 0.0
#  define PRIMME_NUM double
#  define REAL_PART(x) (x)
#  define REAL_PARTZ(x) (x)
#  define CONJ(x) (x)
#  define SUF(NAME) NAME ## _dprimme
#  define SUFX(NAME) Num_d ## NAME ## _dprimme
#  define PREFIX(NAME) d ## NAME
#  define COMPLEXZ(X) (X)
#  define COMPLEXZV(X) (X)
#endif
#define MACHINE_EPSILON 1.11e-16

#endif
