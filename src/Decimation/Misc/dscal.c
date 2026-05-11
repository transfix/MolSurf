/* dscal.f -- translated by f2c (version of 3 February 1990  3:36:42).
   You must link the resulting object file with the libraries:
        -lF77 -lI77 -lm -lc   (in that order)
*/

#include <Decimation/Misc/f2c.h>

/* ----------------------------------------------------------------------- */
/*  IMSL Name:  DSCAL (Double precision version) */

/*  Computer:   VAXUNX/DOUBLE */

/*  Revised:    August 9, 1986 */

/*  Purpose:    Multiply a vector by a scalar, y = ay, both double */
/*              precision. */

/*  Usage:      CALL DSCAL (N, DA, DX, INCX) */

/*  Arguments: */
/*     N      - Length of vector X.  (Input) */
/*     DA     - Double precision scalar.  (Input) */
/*     DX     - Double precision vector of length N*INCX.  (Input/Output) */
/*              DSCAL replaces X(I) with DA*X(I) for I=1,...,N. X(I) */
/*              refers to a specific element of DX. See INCX argument */
/*              description. */
/*     INCX   - Displacement between elements of DX.  (Input) */
/*              X(I) is defined to be DX(1+(I-1)*INCX). INCX must be */
/*              greater than zero. */

/*  GAMS:       D1a6 */

/*  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations */
/*              STAT/LIBRARY Mathematical Support */

/*  Copyright:  1986 by IMSL, Inc.  All Rights Reserved. */

/*  Warranty:   IMSL warrants only that IMSL testing has been applied */
/*              to this code.  No other warranty, expressed or implied, */
/*              is applicable. */

/* ----------------------------------------------------------------------- */

int dscal_(integer *n, doublereal *da, double *dx, integer *incx) {
  /* System generated locals */
  integer i_1, i_2;

  /* Local variables */
  static integer i, m, ns, mp1;

  /* Parameter adjustments */
  --dx;

  /* Function Body */
  /*                                  SPECIFICATIONS FOR ARGUMENTS */
  /*                                  SPECIFICATIONS FOR LOCAL VARIABLES */

  if (*n > 0) {
    if (*incx != 1) {
      /*                                  CODE FOR INCREMENTS NOT EQUAL
      TO 1. */
      ns = *n * *incx;
      i_1 = ns;
      i_2 = *incx;
      for (i = 1; i_2 < 0 ? i >= i_1 : i <= i_1; i += i_2) {
        dx[i] = *da * dx[i];
        /* L10: */
      }
    } else {
      /*                                  CODE FOR INCREMENTS EQUAL TO
      1. */
      /*                                  CLEAN-UP LOOP SO REMAINING
      VECTOR */
      /*                                  LENGTH IS A MULTIPLE OF 5. */
      m = *n - *n / 5 * 5;
      i_2 = m;
      for (i = 1; i <= i_2; ++i) {
        dx[i] = *da * dx[i];
        /* L30: */
      }
      mp1 = m + 1;
      i_2 = *n;
      for (i = mp1; i <= i_2; i += 5) {
        dx[i] = *da * dx[i];
        dx[i + 1] = *da * dx[i + 1];
        dx[i + 2] = *da * dx[i + 2];
        dx[i + 3] = *da * dx[i + 3];
        dx[i + 4] = *da * dx[i + 4];
        /* L40: */
      }
    }
  }
  return 0;
} /* dscal_ */
