#if 0
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Box.h>
#include <X11/Shell.h>
#include <X11/Xaw/Dialog.h>
#endif

#include <stdio.h>
#include <math.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <assert.h>
#include <Decimation/Misc/misc.h>

void Exchangerowcolumn();

/*************************************************************************/
/* GaussInverse                                                          */
/* - Inverse the matrix A by Gauss elimination methods                   */
/*************************************************************************/
void gaussinverse(a, n, eps, message) int n, *message;
double *a, eps;
/*   n        -integer, the dimension of the two dimensional matrix a[n][n].
     a        -pointer, point to a matrix to be inversed and also the
                        computing result
     eps      -controllor, that control the computation by testing the
                        pivot element
     message, -pointer, point to an integer, if pivot element less than eps
                        message = 0,
                        otherwise message = 1.
*/
{
  double max;
  int k, ik, jk, i, j, *z;

  z = (int *)malloc(2 * n * sizeof(int));
  *message = 1;
  for (k = 0; k < n; k++) {
    max = 0.0;
    for (i = k; i < n; i++)
      for (j = k; j < n; j++)
        if (fabs(*(a + i * n + j)) > max) {
          ik = i;
          jk = j;
          max = fabs(*(a + i * n + j));
        };
    if (max < eps || max == 0.0) {
      *message = 0;
      printf("The matrix in gaussinverse is singular, %f\n", max);
      return;
    };
    max = 1.0 / *(a + ik * n + jk);
    *(a + ik * n + jk) = 1.0;
    z[2 * k] = ik;
    z[2 * k + 1] = jk;
    Exchangerowcolumn(a, n, k, ik, jk);
    for (j = 0; j < n; j++)
      *(a + k * n + j) = *(a + k * n + j) * max;
    for (i = 0; i < n; i++)
      if (i != k) {
        max = *(a + i * n + k);
        *(a + i * n + k) = 0.0;
        for (j = 0; j < n; j++)
          *(a + i * n + j) = *(a + i * n + j) - max * *(a + k * n + j);
      }
  };
  for (k = n - 2; k > -1; k--) {
    ik = z[2 * k + 1];
    jk = z[2 * k];
    Exchangerowcolumn(a, n, k, ik, jk);
  }
  free(z);
}

/*************************************************************************/
/* gaussInverse3                                                          */
/* - Inverse 3 * 3  matrix A by Gauss elimination methods                   */
/*************************************************************************/
void gaussinverse3(a, eps, message) int *message;
double *a, eps;
{
  double max;
  int n, k, ik, jk, i, j, z[6];

  n = 3;
  *message = 1;
  for (k = 0; k < n; k++) {
    max = 0.0;
    for (i = k; i < n; i++)
      for (j = k; j < n; j++)
        if (fabs(*(a + i * n + j)) > max) {
          ik = i;
          jk = j;
          max = fabs(*(a + i * n + j));
        };
    if (max < eps || max == 0.0) {
      *message = 0;
      printf("The matrix in gaussinverse3 is singular, %lf\n", max);
      return;
    };
    max = 1.0 / *(a + ik * n + jk);
    *(a + ik * n + jk) = 1.0;
    z[2 * k] = ik;
    z[2 * k + 1] = jk;
    Exchangerowcolumn(a, n, k, ik, jk);
    for (j = 0; j < n; j++)
      *(a + k * n + j) = *(a + k * n + j) * max;
    for (i = 0; i < n; i++)
      if (i != k) {
        max = *(a + i * n + k);
        *(a + i * n + k) = 0.0;
        for (j = 0; j < n; j++)
          *(a + i * n + j) = *(a + i * n + j) - max * *(a + k * n + j);
      }
  };
  for (k = n - 2; k > -1; k--) {
    ik = z[2 * k + 1];
    jk = z[2 * k];
    Exchangerowcolumn(a, n, k, ik, jk);
  }
}

/*************************************************************************/
/* Exchangerowcolumn                                                     */
/* - Exchange two rows and two columns of the matrix A                   */
/*************************************************************************/
void Exchangerowcolumn(a, n, k, ik, jk) int n, k, ik, jk;
double *a;
/*  exchange the (k,ik) rows
                 (k,ik) columns.
*/
{
  double b;
  int j;
  /*  exchange the (k,ik) rows */
  if (ik != k) {
    for (j = 0; j < n; j++) {
      b = *(a + ik * n + j);
      *(a + ik * n + j) = *(a + k * n + j);
      *(a + k * n + j) = b;
    }
  }
  /*  exchange the (k,ik) column  */
  if (jk != k) {
    for (j = 0; j < n; j++) {
      b = *(a + j * n + jk);
      *(a + j * n + jk) = *(a + j * n + k);
      *(a + j * n + k) = b;
    }
  }
}
/*************************************************************************/
/* MatrixTimesMatrix                                                     */
/* times a matrix a with matrix v, the result is x. i.e., x = a*v        */
/*************************************************************************/
void MatrixTimesMatrix(a, v, x, m, n, l) int m, n, l;
double *a; /* pointe to a (m,n) matrix               */
double *v; /* pointe to a (n,l) matrix               */
double *x; /* pointe to a (m,l) matrix               */
{
  int i, j, k;
  double sum;
  for (k = 0; k < l; k++)
    for (i = 0; i < m; i++) {
      sum = 0.0;
      for (j = 0; j < n; j++)
        sum = sum + *(a + i * n + j) * *(v + j * l + k);
      *(x + i * l + k) = sum;
    }
}
/*************************************************************************/
/* LeastSquareSolver                                                      */
/* Solve the equation ax = b in the least square sense                   */
/*************************************************************************/
int LeastSquareSolver(a, m, n, b, eps, solution, error)
double a[];        /* the given matrix whose element is ordered in  */
                   /*   a[i*n +j].                                  */
int m, n;          /* the number of row and column of the matrix a  */
double b[];        /* the right-hand side                           */
double eps;        /* the singularity controller, the function      */
                   /*  return 0 if fails, otherwise return 1        */
double solution[]; /* the solution of the equation                  */
double *error;     /* the least square error                        */
{
  int i, j, k, message;
  double *ata, *atb, sum;

  if ((ata = (double *)malloc(sizeof(double) * n * n)) == NULL) {
    printf("LeastSquareSolver:  malloc() fails - ata\n");
    return (0);
  }
  if ((atb = (double *)malloc(sizeof(double) * n)) == NULL) {
    printf("LeastSquareSolver:  malloc() fails - atb\n");
    return (0);
  }

  *error = 10000.0;
  for (i = 0; i < n; i++)
    solution[i] = 0.0;
  for (k = 0; k < n; k++)
    for (i = 0; i < n; i++) {
      sum = 0.0;
      for (j = 0; j < m; j++)
        sum = sum + *(a + j * n + i) * *(a + j * n + k);
      *(ata + i * n + k) = sum;
    }
  for (i = 0; i < n; i++) {
    sum = 0.0;
    for (j = 0; j < m; j++)
      sum = sum + *(a + j * n + i) * *(b + j);
    *(atb + i) = sum;
  }

  gaussinverse(ata, n, eps, &message);
  if (message == 0) {
    return 0;
  }

  MatrixTimesMatrix(ata, atb, solution, n, n, 1);

  *error = 0.0;
  for (i = 0; i < m; i++) {
    sum = *(b + i);
    for (j = 0; j < n; j++)
      sum = sum - *(a + i * n + j) * solution[j];
    *error = *error + sum * sum;
  }

  free(ata);
  free(atb);
  return 1;
}

/*************************************************************************/
/* QRdecomposition- QR decomposition of a matrix a = QR                  */
/*        Input:    a is m*n matrix stored in row-wise                   */
/*        output:   a is m*n matrix stored in row-wise, a = R            */
/*                  q is m*m matrix stored in row-wise                   */
/*        return:   == 0, fail;   == 1 Ok.                               */
/*************************************************************************/
int QRdecomposition(a, m, n, q)
int m, n;
double a[], q[];
{
  int i, j, k, l, nn, p, jj;
  double u, alpha, w, t;
  if (m < n) {
    printf("fail:row<column!\n");
    return (0);
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++) {
      l = i * m + j;
      q[l] = 0.0;
      if (i == j)
        q[l] = 1.0;
    }
  }
  nn = n;
  if (m == n)
    nn = m - 1;
  for (k = 0; k < nn; k++) {
    u = 0.0;
    l = k * n + k;
    for (i = k; i < m; i++) {
      w = fabs(a[i * n + k]);
      if (w > u)
        u = w;
    }
    alpha = 0.0;
    for (i = k; i < m; i++) {
      t = a[i * n + k] / u;
      alpha += t * t;
    }
    if (a[l] > 0.0)
      u = -u;
    alpha = u * sqrt(alpha);
    if (fabs(alpha) + 1.0 == 1.0) {
      printf("fail:column-linear-dependent!\n");
      return (0);
    }
    u = sqrt(2.0 * alpha * (alpha - a[l]));
    if ((u + 1.0) != 1.0) {
      a[l] = (a[l] - alpha) / u;
      for (i = k + 1; i < m; i++) {
        p = i * n + k;
        a[p] /= u;
      }
      for (j = 0; j < m; j++) {
        t = 0.0;
        for (jj = k; jj < m; jj++)
          t += a[jj * n + k] * q[jj * m + j];
        for (i = k; i < m; i++) {
          p = i * m + j;
          q[p] -= 2.0 * t * a[i * n + k];
        }
      }
      for (j = k + 1; j < n; j++) {
        t = 0.0;
        for (jj = k; jj < m; jj++)
          t += a[jj * n + k] * a[jj * n + j];
        for (i = k; i < m; i++) {
          p = i * n + j;
          a[p] -= 2.0 * t * a[i * n + k];
        }
      }
      a[l] = alpha;
      for (i = k + 1; i < m; i++)
        a[i * n + k] = 0.0;
    }
  }
  for (i = 0; i < m - 1; i++) {
    for (j = i + 1; j < m; j++) {
      p = i * m + j;
      l = j * m + i;
      t = q[p];
      q[p] = q[l];
      q[l] = t;
    }
  }
  return (1);
}

/*************************************************************************/
/* MinimalLeastSquare -- sloving Ax = b with ||x|| = min                 */
/*                       by SVD decompositiom, the matrix is stored in   */
/*                       column wise, where A is a n*p matrix            */
/*   A = U D V^T.  The out put is U,V. The elements is arranged in column*/
/*************************************************************************/
void MinimalLeastSquare(A, ldx, n, p, x, b) double A[], x[], b[];
int ldx;  /* ldx >= n -- the leading dimension of A     */
int n, p; /* (n,p)-- the size of the matrix A           */
{
  double *u, *v, *s, *e, *work;
  int ldu, ldv, job, info, mm, m;
  int i, j, k;

  ldu = n; /* ldu >= n             */
  ldv = p; /* ldv >= p             */
  job = 11;

  u = (double *)malloc(n * n * sizeof(double));
  v = (double *)malloc(p * p * sizeof(double));
  mm = n + 1;
  if (mm < p)
    mm = p;
  s = (double *)malloc(mm * sizeof(double));
  e = (double *)malloc(p * sizeof(double));
  work = (double *)malloc(n * sizeof(double));

  dsvdc_(A, &ldx, &n, &p, s, e, u, &ldu, v, &ldv, work, &job, &info);

  m = n;
  if (n > p)
    m = p;

  /* determine nonzero singular values of x */
  k = 0;
  for (i = 0; i < m; i++) {
    if (s[i] / s[0] < 0.00001)
      break;
    k = k + 1;
  }
  /* compute (u^T b)_k and Y_k              */
  for (i = 0; i < k; i++) {
    work[i] = 0.0;
    for (j = 0; j < n; j++) {
      work[i] = work[i] + u[i * n + j] * b[j];
    }
    work[i] = work[i] / s[i]; /* work is Y_k   */
  }

  /* compute X = V_k V_k                    */
  for (i = 0; i < p; i++) {
    x[i] = 0.0;
    for (j = 0; j < k; j++) {
      x[i] = x[i] + v[i + j * p] * work[j];
    }
  }

  free(u);
  free(v);
  free(s);
  free(e);
  free(work);
}

/*----------------------------------------------------------------------------
 * BGaussInverse -- Inverse the matrix A by Gauss elimination methods
 * Return:    0 - if pivot element less than eps
 *            1 - otherwise
 *            The matrix is first normalized by divide the maximal element.
 *---------------------------------------------------------------------------*/
int Normalized_GaussInverse(a, n, eps)
int n;     /* dimension of the two dimensional matrix a[n][n]      */
double *a, /* matrix to be inversed and also the computing result  */
    eps;   /* control the computation by testing the pivot element */
{
  double max, maxele;
  int k, ik, jk, i, j, *z;

  maxele = 0.0;
  for (k = 0; k < n * n; k++) {
    if (maxele < fabs(a[k]))
      maxele = fabs(a[k]);
  }
  maxele = 1.0 / maxele;
  if (maxele > 1.0) {
    for (k = 0; k < n * n; k++) {
      a[k] = a[k] * maxele;
    }
  }

  z = (int *)malloc(2 * n * sizeof(int));
  for (k = 0; k < n; k++) {
    max = 0.0;
    for (i = k; i < n; i++)
      for (j = k; j < n; j++)
        if (fabs(*(a + i * n + j)) > max) {
          ik = i;
          jk = j;
          max = fabs(*(a + i * n + j));
        }

    if (max < eps)
      return 0;

    max = 1.0 / *(a + ik * n + jk);
    *(a + ik * n + jk) = 1.0;
    z[2 * k] = ik;
    z[2 * k + 1] = jk;
    Exchangerowcolumn(a, n, k, ik, jk);

    for (j = 0; j < n; j++)
      *(a + k * n + j) = *(a + k * n + j) * max;

    for (i = 0; i < n; i++)
      if (i != k) {
        max = *(a + i * n + k);
        *(a + i * n + k) = 0.0;
        for (j = 0; j < n; j++)
          *(a + i * n + j) = *(a + i * n + j) - max * *(a + k * n + j);
      }
  }

  for (k = n - 2; k > -1; k--) {
    ik = z[2 * k + 1];
    jk = z[2 * k];
    Exchangerowcolumn(a, n, k, ik, jk);
  }

  if (maxele > 1.0) {
    for (k = 0; k < n * n; k++) {
      a[k] = a[k] * maxele;
    }
  }

  free(z);
  return 1;
}
