
/* dsvdc.f -- translated by f2c (version of 3 February 1990  3:36:42).
   You must link the resulting object file with the libraries:
        -lF77 -lI77 -lm -lc   (in that order)
*/

#include <math.h>
#include <Decimation/Misc/f2c.h>

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b44 = -1.;
static doublereal c_b141 = 1.;

int dsvdc_(doublereal *x, integer *ldx, integer *n, integer *p, doublereal *s,
           doublereal *e, doublereal *u, integer *ldu, doublereal *v,
           integer *ldv, doublereal *work, integer *job, integer *info)
/*
doublereal *x;
integer *ldx, *n, *p;
doublereal *s, *e, *u;
integer *ldu;
doublereal *v;
integer *ldv;
doublereal *work;
integer *job, *info;
*/
{
  /* System generated locals */
  integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i_1, i_2, i_3;
  doublereal d_1, d_2, d_3, d_4, d_5, d_6, d_7;

  /* Builtin functions */
  extern double d_sign(doublereal *, doublereal *);
  extern double sqrt(double);

  /* Local variables */
  static integer kase;
  extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
                          integer *);
  static integer jobu, iter;
  extern drot_(integer *, doublereal *, integer *, doublereal *, integer *,
               doublereal *, doublereal *);
  static doublereal test;
  extern doublereal dnrm2_(integer *, doublereal *, integer *);
  static integer nctp1;
  static doublereal b, c;
  static integer nrtp1;
  static doublereal f, g;
  static integer i, j, k, l, m;
  static doublereal t, scale;
  extern int dscal_(integer *, doublereal *, doublereal *, integer *);
  static doublereal shift;
  extern int dswap_(integer *, doublereal *, integer *, doublereal *,
                    integer *);
  extern int drotg_(doublereal *, doublereal *, doublereal *, doublereal *);
  static integer maxit;
  extern int daxpy_(integer *, doublereal *, doublereal *, integer *,
                    doublereal *, integer *);
  static logical wantu, wantv;
  static doublereal t1, ztest, el;
  static integer kk;
  static doublereal cs;
  static integer ll, mm, ls;
  static doublereal sl;
  static integer lu;
  static doublereal sm, sn;
  static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
  static doublereal emm1, smm1;

  /* Parameter adjustments */
  x_dim1 = *ldx;
  x_offset = x_dim1 + 1;
  x -= x_offset;
  --s;
  --e;
  u_dim1 = *ldu;
  u_offset = u_dim1 + 1;
  u -= u_offset;
  v_dim1 = *ldv;
  v_offset = v_dim1 + 1;
  v -= v_offset;
  --work;

  /* Function Body */

  /*     dsvdc is a subroutine to reduce a double precision nxp matrix x */
  /*     by orthogonal transformations u and v to diagonal form.  the */
  /*     diagonal elements s(i) are the singular values of x.  the */
  /*     columns of u are the corresponding left singular vectors, */
  /*     and the columns of v the right singular vectors. */

  /*     on entry */

  /*         x         double precision(ldx,p), where ldx.ge.n. */
  /*                   x contains the matrix whose singular value */
  /*                   decomposition is to be computed.  x is */
  /*                   destroyed by dsvdc. */

  /*         ldx       integer. */
  /*                   ldx is the leading dimension of the array x. */

  /*         n         integer. */
  /*                   n is the number of rows of the matrix x. */

  /*         p         integer. */
  /*                   p is the number of columns of the matrix x. */

  /*         ldu       integer. */
  /*                   ldu is the leading dimension of the array u. */
  /*                   (see below). */

  /*         ldv       integer. */
  /*                   ldv is the leading dimension of the array v. */
  /*                   (see below). */

  /*         work      double precision(n). */
  /*                   work is a scratch array. */

  /*         job       integer. */
  /*                   job controls the computation of the singular */
  /*                   vectors.  it has the decimal expansion ab */
  /*                   with the following meaning */

  /*                        a.eq.0    do not compute the left singular */
  /*                                  vectors. */
  /*                        a.eq.1    return the n left singular vectors */
  /*                                  in u. */
  /*                        a.ge.2    return the first min(n,p) singular */
  /*                                  vectors in u. */
  /*                        b.eq.0    do not compute the right singular */
  /*                                  vectors. */
  /*                        b.eq.1    return the right singular vectors */
  /*                                  in v. */

  /*     on return */

  /*         s         double precision(mm), where mm=min(n+1,p). */
  /*                   the first min(n,p) entries of s contain the */
  /*                   singular values of x arranged in descending */
  /*                   order of magnitude. */

  /*         e         double precision(p), */
  /*                   e ordinarily contains zeros.  however see the */
  /*                   discussion of info for exceptions. */

  /*         u         double precision(ldu,k), where ldu.ge.n.  if */
  /*                                   joba.eq.1 then k.eq.n, if joba.ge.2
   */
  /*                                   then k.eq.min(n,p). */
  /*                   u contains the matrix of left singular vectors. */
  /*                   u is not referenced if joba.eq.0.  if n.le.p */
  /*                   or if joba.eq.2, then u may be identified with x */
  /*                   in the subroutine call. */

  /*         v         double precision(ldv,p), where ldv.ge.p. */
  /*                   v contains the matrix of right singular vectors. */
  /*                   v is not referenced if job.eq.0.  if p.le.n, */
  /*                   then v may be identified with x in the */
  /*                   subroutine call. */

  /*         info      integer. */
  /*                   the singular values (and their corresponding */
  /*                   singular vectors) s(info+1),s(info+2),...,s(m) */
  /*                   are correct (here m=min(n,p)).  thus if */
  /*                   info.eq.0, all the singular values and their */
  /*                   vectors are correct.  in any event, the matrix */
  /*                   b = trans(u)*x*v is the bidiagonal matrix */
  /*                   with the elements of s on its diagonal and the */
  /*                   elements of e on its super-diagonal (trans(u) */
  /*                   is the transpose of u).  thus the singular */
  /*                   values of x and b are the same. */

  /*     linpack. this version dated 08/14/78 . */
  /*              correction made to shift 2/84. */
  /*     g.w. stewart, university of maryland, argonne national lab. */

  /*     dsvdc uses the following functions and subprograms. */

  /*     external drot */
  /*     blas daxpy,ddot,dscal,dswap,dnrm2,drotg */
  /*     fortran dabs,dmax1,max0,min0,mod,dsqrt */

  /*     internal variables */

  /*     set the maximum number of iterations. */

  maxit = 150;

  /*     determine what is to be computed. */

  wantu = FALSE_;
  wantv = FALSE_;
  jobu = *job % 100 / 10;
  ncu = *n;
  if (jobu > 1) {
    ncu = min(*n, *p);
  }
  if (jobu != 0) {
    wantu = TRUE_;
  }
  if (*job % 10 != 0) {
    wantv = TRUE_;
  }

  /*     reduce x to bidiagonal form, storing the diagonal elements */
  /*     in s and the super-diagonal elements in e. */

  *info = 0;
  /* Computing MAX */
  i_1 = *n - 1;
  nct = min(*p, i_1);
  /* Computing MAX */
  /* Computing MAX */
  i_3 = *p - 2;
  i_1 = 0, i_2 = min(*n, i_3);
  nrt = max(i_2, i_1);
  lu = max(nct, nrt);
  if (lu < 1) {
    goto L170;
  }
  i_1 = lu;
  for (l = 1; l <= i_1; ++l) {
    lp1 = l + 1;
    if (l > nct) {
      goto L20;
    }

    /*           compute the transformation for the l-th column and */
    /*           place the l-th diagonal in s(l). */

    i_2 = *n - l + 1;
    s[l] = dnrm2_(&i_2, &x[l + l * x_dim1], &c__1);
    if (s[l] == 0.) {
      goto L10;
    }
    if (x[l + l * x_dim1] != 0.) {
      s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
    }
    i_2 = *n - l + 1;
    d_1 = 1. / s[l];
    dscal_(&i_2, &d_1, &x[l + l * x_dim1], &c__1);
    x[l + l * x_dim1] += 1.;
  L10:
    s[l] = -s[l];
  L20:
    if (*p < lp1) {
      goto L50;
    }
    i_2 = *p;
    for (j = lp1; j <= i_2; ++j) {
      if (l > nct) {
        goto L30;
      }
      if (s[l] == 0.) {
        goto L30;
      }

      /*              apply the transformation. */

      i_3 = *n - l + 1;
      t = -ddot_(&i_3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1) /
          x[l + l * x_dim1];
      i_3 = *n - l + 1;
      daxpy_(&i_3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &c__1);
    L30:

      /*           place the l-th row of x into  e for the */
      /*           subsequent calculation of the row transformation. */

      e[j] = x[l + j * x_dim1];
      /* L40: */
    }
  L50:
    if (!wantu || l > nct) {
      goto L70;
    }

    /*           place the transformation in u for subsequent back */
    /*           multiplication. */

    i_2 = *n;
    for (i = l; i <= i_2; ++i) {
      u[i + l * u_dim1] = x[i + l * x_dim1];
      /* L60: */
    }
  L70:
    if (l > nrt) {
      goto L150;
    }

    /*           compute the l-th row transformation and place the */
    /*           l-th super-diagonal in e(l). */

    i_2 = *p - l;
    e[l] = dnrm2_(&i_2, &e[lp1], &c__1);
    if (e[l] == 0.) {
      goto L80;
    }
    if (e[lp1] != 0.) {
      e[l] = d_sign(&e[l], &e[lp1]);
    }
    i_2 = *p - l;
    d_1 = 1. / e[l];
    dscal_(&i_2, &d_1, &e[lp1], &c__1);
    e[lp1] += 1.;
  L80:
    e[l] = -e[l];
    if (lp1 > *n || e[l] == 0.) {
      goto L120;
    }

    /*              apply the transformation. */

    i_2 = *n;
    for (i = lp1; i <= i_2; ++i) {
      work[i] = 0.;
      /* L90: */
    }
    i_2 = *p;
    for (j = lp1; j <= i_2; ++j) {
      i_3 = *n - l;
      daxpy_(&i_3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &c__1);
      /* L100: */
    }
    i_2 = *p;
    for (j = lp1; j <= i_2; ++j) {
      i_3 = *n - l;
      d_1 = -e[j] / e[lp1];
      daxpy_(&i_3, &d_1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &c__1);
      /* L110: */
    }
  L120:
    if (!wantv) {
      goto L140;
    }

    /*              place the transformation in v for subsequent */
    /*              back multiplication. */

    i_2 = *p;
    for (i = lp1; i <= i_2; ++i) {
      v[i + l * v_dim1] = e[i];
      /* L130: */
    }
  L140:
  L150:
      /* L160: */
      ;
  }
L170:

  /*     set up the final bidiagonal matrix or order m. */

  /* Computing MAX */
  i_1 = *p, i_2 = *n + 1;
  m = min(i_2, i_1);
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;
  if (nct < *p) {
    s[nctp1] = x[nctp1 + nctp1 * x_dim1];
  }
  if (*n < m) {
    s[m] = 0.;
  }
  if (nrtp1 < m) {
    e[nrtp1] = x[nrtp1 + m * x_dim1];
  }
  e[m] = 0.;

  /*     if required, generate u. */

  if (!wantu) {
    goto L300;
  }
  if (ncu < nctp1) {
    goto L200;
  }
  i_1 = ncu;
  for (j = nctp1; j <= i_1; ++j) {
    i_2 = *n;
    for (i = 1; i <= i_2; ++i) {
      u[i + j * u_dim1] = 0.;
      /* L180: */
    }
    u[j + j * u_dim1] = 1.;
    /* L190: */
  }
L200:
  if (nct < 1) {
    goto L290;
  }
  i_1 = nct;
  for (ll = 1; ll <= i_1; ++ll) {
    l = nct - ll + 1;
    if (s[l] == 0.) {
      goto L250;
    }
    lp1 = l + 1;
    if (ncu < lp1) {
      goto L220;
    }
    i_2 = ncu;
    for (j = lp1; j <= i_2; ++j) {
      i_3 = *n - l + 1;
      t = -ddot_(&i_3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &c__1) /
          u[l + l * u_dim1];
      i_3 = *n - l + 1;
      daxpy_(&i_3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &c__1);
      /* L210: */
    }
  L220:
    i_2 = *n - l + 1;
    dscal_(&i_2, &c_b44, &u[l + l * u_dim1], &c__1);
    u[l + l * u_dim1] += 1.;
    lm1 = l - 1;
    if (lm1 < 1) {
      goto L240;
    }
    i_2 = lm1;
    for (i = 1; i <= i_2; ++i) {
      u[i + l * u_dim1] = 0.;
      /* L230: */
    }
  L240:
    goto L270;
  L250:
    i_2 = *n;
    for (i = 1; i <= i_2; ++i) {
      u[i + l * u_dim1] = 0.;
      /* L260: */
    }
    u[l + l * u_dim1] = 1.;
  L270:
      /* L280: */
      ;
  }
L290:
L300:

  /*     if it is required, generate v. */

  if (!wantv) {
    goto L350;
  }
  i_1 = *p;
  for (ll = 1; ll <= i_1; ++ll) {
    l = *p - ll + 1;
    lp1 = l + 1;
    if (l > nrt) {
      goto L320;
    }
    if (e[l] == 0.) {
      goto L320;
    }
    i_2 = *p;
    for (j = lp1; j <= i_2; ++j) {
      i_3 = *p - l;
      t = -ddot_(&i_3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * v_dim1],
                 &c__1) /
          v[lp1 + l * v_dim1];
      i_3 = *p - l;
      daxpy_(&i_3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * v_dim1],
             &c__1);
      /* L310: */
    }
  L320:
    i_2 = *p;
    for (i = 1; i <= i_2; ++i) {
      v[i + l * v_dim1] = 0.;
      /* L330: */
    }
    v[l + l * v_dim1] = 1.;
    /* L340: */
  }
L350:

  /*     main iteration loop for the singular values. */

  mm = m;
  iter = 0;
L360:

  /*        quit if all the singular values have been found. */

  /*     ...exit */
  if (m == 0) {
    goto L620;
  }

  /*        if too many iterations have been performed, set */
  /*        flag and return. */

  if (iter < maxit) {
    goto L370;
  }
  *info = m;
  /*     ......exit */
  goto L620;
L370:

  /*        this section of the program inspects for */
  /*        negligible elements in the s and e arrays.  on */
  /*        completion the variables kase and l are set as follows. */

  /*           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m */
  /*           kase = 2     if s(l) is negligible and l.lt.m */
  /*           kase = 3     if e(l-1) is negligible, l.lt.m, and */
  /*                        s(l), ..., s(m) are not negligible (qr step). */

  /*           kase = 4     if e(m-1) is negligible (convergence). */

  i_1 = m;
  for (ll = 1; ll <= i_1; ++ll) {
    l = m - ll;
    /*        ...exit */
    if (l == 0) {
      goto L400;
    }
    test = (d_1 = s[l], abs(d_1)) + (d_2 = s[l + 1], abs(d_2));
    ztest = test + (d_1 = e[l], abs(d_1));
    if (ztest != test) {
      goto L380;
    }
    e[l] = 0.;
    /*        ......exit */
    goto L400;
  L380:
      /* L390: */
      ;
  }
L400:
  if (l != m - 1) {
    goto L410;
  }
  kase = 4;
  goto L480;
L410:
  lp1 = l + 1;
  mp1 = m + 1;
  i_1 = mp1;
  for (lls = lp1; lls <= i_1; ++lls) {
    ls = m - lls + lp1;
    /*           ...exit */
    if (ls == l) {
      goto L440;
    }
    test = 0.;
    if (ls != m) {
      test += (d_1 = e[ls], abs(d_1));
    }
    if (ls != l + 1) {
      test += (d_1 = e[ls - 1], abs(d_1));
    }
    ztest = test + (d_1 = s[ls], abs(d_1));
    if (ztest != test) {
      goto L420;
    }
    s[ls] = 0.;
    /*           ......exit */
    goto L440;
  L420:
      /* L430: */
      ;
  }
L440:
  if (ls != l) {
    goto L450;
  }
  kase = 3;
  goto L470;
L450:
  if (ls != m) {
    goto L460;
  }
  kase = 1;
  goto L470;
L460:
  kase = 2;
  l = ls;
L470:
L480:
  ++l;

  /*        perform the task indicated by kase. */

  switch (kase) {
  case 1:
    goto L490;
  case 2:
    goto L520;
  case 3:
    goto L540;
  case 4:
    goto L570;
  }

  /*        deflate negligible s(m). */

L490:
  mm1 = m - 1;
  f = e[m - 1];
  e[m - 1] = 0.;
  i_1 = mm1;
  for (kk = l; kk <= i_1; ++kk) {
    k = mm1 - kk + l;
    t1 = s[k];
    drotg_(&t1, &f, &cs, &sn);
    s[k] = t1;
    if (k == l) {
      goto L500;
    }
    f = -sn * e[k - 1];
    e[k - 1] = cs * e[k - 1];
  L500:
    if (wantv) {
      drot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &cs, &sn);
    }
    /* L510: */
  }
  goto L610;

  /*        split at negligible s(l). */

L520:
  f = e[l - 1];
  e[l - 1] = 0.;
  i_1 = m;
  for (k = l; k <= i_1; ++k) {
    t1 = s[k];
    drotg_(&t1, &f, &cs, &sn);
    s[k] = t1;
    f = -sn * e[k];
    e[k] = cs * e[k];
    if (wantu) {
      drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &c__1, &cs,
            &sn);
    }
    /* L530: */
  }
  goto L610;

  /*        perform one qr step. */

L540:

  /*           calculate the shift. */

  /* Computing MAX */
  d_6 = (d_1 = s[m], abs(d_1)), d_7 = (d_2 = s[m - 1], abs(d_2)),
  d_6 = max(d_7, d_6), d_7 = (d_3 = e[m - 1], abs(d_3)), d_6 = max(d_7, d_6),
  d_7 = (d_4 = s[l], abs(d_4)), d_6 = max(d_7, d_6),
  d_7 = (d_5 = e[l], abs(d_5));
  scale = max(d_7, d_6);
  sm = s[m] / scale;
  smm1 = s[m - 1] / scale;
  emm1 = e[m - 1] / scale;
  sl = s[l] / scale;
  el = e[l] / scale;
  /* Computing 2nd power */
  d_1 = emm1;
  b = ((smm1 + sm) * (smm1 - sm) + d_1 * d_1) / 2.;
  /* Computing 2nd power */
  d_1 = sm * emm1;
  c = d_1 * d_1;
  shift = 0.;
  if (b == 0. && c == 0.) {
    goto L550;
  }
  /* Computing 2nd power */
  d_1 = b;
  shift = sqrt(d_1 * d_1 + c);
  if (b < 0.) {
    shift = -shift;
  }
  shift = c / (b + shift);
L550:
  f = (sl + sm) * (sl - sm) + shift;
  g = sl * el;

  /*           chase zeros. */

  mm1 = m - 1;
  i_1 = mm1;
  for (k = l; k <= i_1; ++k) {
    drotg_(&f, &g, &cs, &sn);
    if (k != l) {
      e[k - 1] = f;
    }
    f = cs * s[k] + sn * e[k];
    e[k] = cs * e[k] - sn * s[k];
    g = sn * s[k + 1];
    s[k + 1] = cs * s[k + 1];
    if (wantv) {
      drot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &c__1, &cs,
            &sn);
    }
    drotg_(&f, &g, &cs, &sn);
    s[k] = f;
    f = cs * e[k] + sn * s[k + 1];
    s[k + 1] = -sn * e[k] + cs * s[k + 1];
    g = sn * e[k + 1];
    e[k + 1] = cs * e[k + 1];
    if (wantu && k < *n) {
      drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &c__1, &cs,
            &sn);
    }
    /* L560: */
  }
  e[m - 1] = f;
  ++iter;
  goto L610;

  /*        convergence. */

L570:

  /*           make the singular value  positive. */

  if (s[l] >= 0.) {
    goto L580;
  }
  s[l] = -s[l];
  if (wantv) {
    dscal_(p, &c_b44, &v[l * v_dim1 + 1], &c__1);
  }
L580:

  /*           order the singular value. */

L590:
  if (l == mm) {
    goto L600;
  }
  /*           ...exit */
  if (s[l] >= s[l + 1]) {
    goto L600;
  }
  t = s[l];
  s[l] = s[l + 1];
  s[l + 1] = t;
  if (wantv && l < *p) {
    dswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
  }
  if (wantu && l < *n) {
    dswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
  }
  ++l;
  goto L590;
L600:
  iter = 0;
  --m;
L610:
  goto L360;
L620:
  return 0;
} /* dsvdc_ */

int dswap_(integer *n, doublereal *dx, integer *incx, doublereal *dy,
           integer *incy) {
  /* System generated locals */
  integer i_1;

  /* Local variables */
  static integer i, m;
  static doublereal dtemp;
  static integer ix, iy, mp1;

  /* Parameter adjustments */
  --dx;
  --dy;

  /* Function Body */

  /*     interchanges two vectors. */
  /*     uses unrolled loops for increments equal one. */
  /*     jack dongarra, linpack, 3/11/78. */

  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  /*       code for unequal increments or equal increments not equal */
  /*         to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) {
    dtemp = dx[ix];
    dx[ix] = dy[iy];
    dy[iy] = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*       code for both increments equal to 1 */

  /*       clean-up loop */

L20:
  m = *n % 3;
  if (m == 0) {
    goto L40;
  }
  i_1 = m;
  for (i = 1; i <= i_1; ++i) {
    dtemp = dx[i];
    dx[i] = dy[i];
    dy[i] = dtemp;
    /* L30: */
  }
  if (*n < 3) {
    return 0;
  }
L40:
  mp1 = m + 1;
  i_1 = *n;
  for (i = mp1; i <= i_1; i += 3) {
    dtemp = dx[i];
    dx[i] = dy[i];
    dy[i] = dtemp;
    dtemp = dx[i + 1];
    dx[i + 1] = dy[i + 1];
    dy[i + 1] = dtemp;
    dtemp = dx[i + 2];
    dx[i + 2] = dy[i + 2];
    dy[i + 2] = dtemp;
    /* L50: */
  }
  return 0;
} /* dswap_ */

int daxpy_(integer *n, doublereal *da, doublereal *dx, integer *incx,
           doublereal *dy, integer *incy) {
  /* System generated locals */
  integer i_1;

  /* Local variables */
  static integer i, m, ix, iy, mp1;

  /* Parameter adjustments */
  --dx;
  --dy;

  /* Function Body */

  /*     constant times a vector plus a vector. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */

  if (*n <= 0) {
    return 0;
  }
  if (*da == 0.) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  /*        code for unequal increments or equal increments */
  /*          not equal to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) {
    dy[iy] += *da * dx[ix];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*        code for both increments equal to 1 */

  /*        clean-up loop */

L20:
  m = *n % 4;
  if (m == 0) {
    goto L40;
  }
  i_1 = m;
  for (i = 1; i <= i_1; ++i) {
    dy[i] += *da * dx[i];
    /* L30: */
  }
  if (*n < 4) {
    return 0;
  }
L40:
  mp1 = m + 1;
  i_1 = *n;
  for (i = mp1; i <= i_1; i += 4) {
    dy[i] += *da * dx[i];
    dy[i + 1] += *da * dx[i + 1];
    dy[i + 2] += *da * dx[i + 2];
    dy[i + 3] += *da * dx[i + 3];
    /* L50: */
  }
  return 0;
} /* daxpy_ */

doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy,
                 integer *incy) {
  /* System generated locals */
  integer i_1;
  doublereal ret_val;

  /* Local variables */
  static integer i, m;
  static doublereal dtemp;
  static integer ix, iy, mp1;

  /* Parameter adjustments */
  --dx;
  --dy;

  /* Function Body */

  /*     forms the dot product of two vectors. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */

  ret_val = 0.;
  dtemp = 0.;
  if (*n <= 0) {
    return ret_val;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  /*        code for unequal increments or equal increments */
  /*          not equal to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) {
    dtemp += dx[ix] * dy[iy];
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  ret_val = dtemp;
  return ret_val;

  /*        code for both increments equal to 1 */

  /*        clean-up loop */

L20:
  m = *n % 5;
  if (m == 0) {
    goto L40;
  }
  i_1 = m;
  for (i = 1; i <= i_1; ++i) {
    dtemp += dx[i] * dy[i];
    /* L30: */
  }
  if (*n < 5) {
    goto L60;
  }
L40:
  mp1 = m + 1;
  i_1 = *n;
  for (i = mp1; i <= i_1; i += 5) {
    dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] +
            dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] +
            dx[i + 4] * dy[i + 4];
    /* L50: */
  }
L60:
  ret_val = dtemp;
  return ret_val;
} /* ddot_ */

doublereal dnrm2_(integer *n, doublereal *dx, integer *incx) {
  /* Initialized data */
  static doublereal zero = 0.;
  static doublereal one = 1.;
  static doublereal cutlo = 8.232e-11;
  static doublereal cuthi = 1.304e19;

  /* System generated locals */
  integer i_1, i_2;
  doublereal ret_val, d_1;

  /* Builtin functions */
  double sqrt(double);

  /* Local variables */
  static doublereal xmax;
  static integer next, i, j, nn;
  static doublereal hitest, sum;

  /* Parameter adjustments */
  --dx;

  /* Function Body */

  /*     euclidean norm of the n-vector stored in dx() with storage */
  /*     increment incx . */
  /*     if    n .le. 0 return with result = 0. */
  /*     if n .ge. 1 then incx must be .ge. 1 */

  /*           c.l.lawson, 1978 jan 08 */

  /*     four phase method     using two built-in constants that are */
  /*     hopefully applicable to all machines. */
  /*         cutlo = maximum of  dsqrt(u/eps)  over all known machines. */
  /*         cuthi = minimum of  dsqrt(v)      over all known machines. */
  /*     where */
  /*         eps = smallest no. such that eps + 1. .gt. 1. */
  /*         u   = smallest positive no.   (underflow limit) */
  /*         v   = largest  no.            (overflow  limit) */

  /*     brief outline of algorithm.. */

  /*     phase 1    scans zero components. */
  /*     move to phase 2 when a component is nonzero and .le. cutlo */
  /*     move to phase 3 when a component is .gt. cutlo */
  /*     move to phase 4 when a component is .ge. cuthi/m */
  /*     where m = n for x() real and m = 2*n for complex. */

  /*     values for cutlo and cuthi.. */
  /*     from the environmental parameters listed in the imsl converter */
  /*     document the limiting values are as follows.. */
  /*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
   */
  /*                   univac and dec at 2**(-103) */
  /*                   thus cutlo = 2**(-51) = 4.44089e-16 */
  /*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
  /*                   thus cuthi = 2**(63.5) = 1.30438e19 */
  /*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
  /*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
  /*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
  /*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
  /*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

  if (*n > 0) {
    goto L10;
  }
  ret_val = zero;
  goto L300;

L10:
  next = 0;
  sum = zero;
  nn = *n * *incx;
  /*                                                 begin main loop */
  i = 1;
L20:
  switch (next) {
  case 0:
    goto L30;
  case 1:
    goto L50;
  case 2:
    goto L70;
  case 3:
    goto L110;
  }
L30:
  if ((d_1 = dx[i], abs(d_1)) > cutlo) {
    goto L85;
  }
  next = 1;
  xmax = zero;

  /*                        phase 1.  sum is zero */

L50:
  if (dx[i] == zero) {
    goto L200;
  }
  if ((d_1 = dx[i], abs(d_1)) > cutlo) {
    goto L85;
  }

  /*                                prepare for phase 2. */
  next = 2;
  goto L105;

  /*                                prepare for phase 4. */

L100:
  i = j;
  next = 3;
  sum = sum / dx[i] / dx[i];
L105:
  xmax = (d_1 = dx[i], abs(d_1));
  goto L115;

  /*                   phase 2.  sum is small. */
  /*                             scale to avoid destructive underflow. */

L70:
  if ((d_1 = dx[i], abs(d_1)) > cutlo) {
    goto L75;
  }

  /*                     common code for phases 2 and 4. */
  /*                     in phase 4 sum is large.  scale to avoid overflow.
   */

L110:
  if ((d_1 = dx[i], abs(d_1)) <= xmax) {
    goto L115;
  }
  /* Computing 2nd power */
  d_1 = xmax / dx[i];
  sum = one + sum * (d_1 * d_1);
  xmax = (d_1 = dx[i], abs(d_1));
  goto L200;

L115:
  /* Computing 2nd power */
  d_1 = dx[i] / xmax;
  sum += d_1 * d_1;
  goto L200;

  /*                  prepare for phase 3. */

L75:
  sum = sum * xmax * xmax;

  /*     for real or d.p. set hitest = cuthi/n */
  /*     for complex      set hitest = cuthi/(2*n) */

L85:
  hitest = cuthi / (real)(*n);

  /*                   phase 3.  sum is mid-range.  no scaling. */

  i_1 = nn;
  i_2 = *incx;
  for (j = i; i_2 < 0 ? j >= i_1 : j <= i_1; j += i_2) {
    if ((d_1 = dx[j], abs(d_1)) >= hitest) {
      goto L100;
    }
    /* L95: */
    /* Computing 2nd power */
    d_1 = dx[j];
    sum += d_1 * d_1;
  }
  ret_val = sqrt(sum);
  goto L300;

L200:
  i += *incx;
  if (i <= nn) {
    goto L20;
  }

  /*              end of main loop. */

  /*              compute square root and adjust for scaling. */

  ret_val = xmax * sqrt(sum);
L300:
  return ret_val;
} /* dnrm2_ */

int drot_(integer *n, doublereal *dx, integer *incx, doublereal *dy,
          integer *incy, doublereal *c, doublereal *s) {
  /* System generated locals */
  integer i_1;

  /* Local variables */
  static integer i;
  static doublereal dtemp;
  static integer ix, iy;

  /* Parameter adjustments */
  --dx;
  --dy;

  /* Function Body */

  /*     applies a plane rotation. */
  /*     jack dongarra, linpack, 3/11/78. */

  if (*n <= 0) {
    return 0;
  }
  if (*incx == 1 && *incy == 1) {
    goto L20;
  }

  /*       code for unequal increments or equal increments not equal */
  /*         to 1 */

  ix = 1;
  iy = 1;
  if (*incx < 0) {
    ix = (-(*n) + 1) * *incx + 1;
  }
  if (*incy < 0) {
    iy = (-(*n) + 1) * *incy + 1;
  }
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) {
    dtemp = *c * dx[ix] + *s * dy[iy];
    dy[iy] = *c * dy[iy] - *s * dx[ix];
    dx[ix] = dtemp;
    ix += *incx;
    iy += *incy;
    /* L10: */
  }
  return 0;

  /*       code for both increments equal to 1 */

L20:
  i_1 = *n;
  for (i = 1; i <= i_1; ++i) {
    dtemp = *c * dx[i] + *s * dy[i];
    dy[i] = *c * dy[i] - *s * dx[i];
    dx[i] = dtemp;
    /* L30: */
  }
  return 0;
} /* drot_ */

int drotg_(doublereal *da, doublereal *db, doublereal *c, doublereal *s) {
  /* System generated locals */
  doublereal d_1, d_2;

  /* Builtin functions */
  double sqrt(double), d_sign(doublereal *, doublereal *);

  /* Local variables */
  static doublereal r, scale, z, roe;

  /*     construct givens plane rotation. */
  /*     jack dongarra, linpack, 3/11/78. */

  roe = *db;
  if (abs(*da) > abs(*db)) {
    roe = *da;
  }
  scale = abs(*da) + abs(*db);
  if (scale != 0.) {
    goto L10;
  }
  *c = 1.;
  *s = 0.;
  r = 0.;
  goto L20;
L10:
  /* Computing 2nd power */
  d_1 = *da / scale;
  /* Computing 2nd power */
  d_2 = *db / scale;
  r = scale * sqrt(d_1 * d_1 + d_2 * d_2);
  r = d_sign(&c_b141, &roe) * r;
  *c = *da / r;
  *s = *db / r;
L20:
  z = 1.;
  if (abs(*da) > abs(*db)) {
    z = *s;
  }
  if (abs(*db) >= abs(*da) && *c != 0.) {
    z = 1. / *c;
  }
  *da = r;
  *db = z;
  return 0;
} /* drotg_ */

double d_sign(doublereal *a, doublereal *b) {
  double x;
  x = (*a >= 0 ? *a : -*a);
  return (*b >= 0 ? x : -x);
}

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
