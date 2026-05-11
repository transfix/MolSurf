/*****************************************************************************/
/*****************************************************************************/
/**                                                                         **/
/** This SHASTRA software is not in the Public Domain. It is distributed on **/
/** a person to person basis, solely for educational use and permission is  **/
/** NOT granted for its transfer to anyone or for its use in any commercial **/
/** product.  There is NO warranty on the available software and neither    **/
/** Purdue University nor the Applied Algebra and Geometry group directed   **/
/** by C.  Bajaj accept responsibility for the consequences of its use.     **/
/**                                                                         **/
/*****************************************************************************/
/*****************************************************************************/

/*
 * BezierMath.c  --  The library  of the Bezier related functions.
 *      created by:     Xu Guoliang, Lin Youming
 *      from:           Sep., 29, 1993
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <Decimation/Bezier/BezierMath.h>

static char errmsg[80]; /* error message */

/* prototypes for internal functions */
static float TensorCastel_k(/* float *bern, float *t, int *n, int k */);
static float TensorCastel_1(/* float *bern, float *t, int *n */);
static float TensorCastel_2(/* float *bern, float *t, int *n */);
static float TensorCastel_3(/* float *bern, float *t, int *n */);
static float TensorCastel_4(/* float *bern, float *t, int *n */);
static float BaryCastel_k(/* float *bern, float *t, int *n, int k */);
static float BaryCastel_1(/* float *bern, float *t, int *n */);
static float BaryCastel_2(/* float *bern, float *t, int *n */);
static float BaryCastel_3(/* float *bern, float *t, int *n */);
static float BaryCastel_4(/* float *bern, float *t, int *n */);
static void ExchangeRowColumn(/* double *a,int n,int k,int ik,int jk */);

/*---------------------------------------------------------------------------
 * BernPolyCreate ---
 *     create a new Bernstein polynomial. The degrees and
 *     dimension is create. The space for coefficients and
 *     domain are also create. Specification of domain and
 *     coefficients are option. They can be set by other
 *     routines. A zero poly will be created if coefficients
 *     are not specified.
 *---------------------------------------------------------------------------*/
BernPoly BernPolyCreate(barydim, tensordim, barydegree, tensordegree, coeff,
                        simplex, hypercube)
int barydim, tensordim; /* (barydim+tensordim) is the real dim        */
int barydegree;         /* degree of the barycentric part             */
int tensordegree[];     /* degree of the tensor part that has
                           length tensordim                           */
float coeff[];          /* coefficients          */
float *simplex;         /* [barydim+1][barydim], vertices of simplex */
float hypercube[][2];   /* [tensordim][2], knotes of tensor part     */
{
  BernPoly newpoly;
  int i, max;

  errmsg[0] = '\0';

  /* check dimension */
  if ((barydim < 0) || (tensordim < 0) || (barydim + tensordim) <= 0) {
    strcpy(errmsg, "Dimension error in creating poly");
    return (NULL);
  }

  /* check degree of barycentric part  */
  if ((barydim > 0) && (barydegree < 0)) {
    strcpy(errmsg, "Barycentric degree error in creating poly");
    return (NULL);
  }

  /* check degree of tensor part  */
  if ((tensordim > 0) && tensordegree != NULL)
    for (i = 0; i < tensordim; i++)
      if (tensordegree[i] < 0) {
        strcpy(errmsg, "Tensor degree error in creating poly");
        return (NULL);
      }

  newpoly = (BernPoly)malloc(sizeof(struct BernPoly));

  newpoly->barydim = barydim;
  newpoly->tensordim = tensordim;
  newpoly->barydegree = barydegree;

  /* allocate memory for simplex, and put number in */
  if (barydim == 0) /* no barycentric part */
    newpoly->simplex = NULL;
  else {
    newpoly->simplex = (OneVertex *)malloc((barydim + 1) * sizeof(OneVertex));
    for (i = 0; i <= barydim; i++)
      (newpoly->simplex[i]).vertex = (float *)malloc(barydim * sizeof(float));

    if (simplex != NULL)
      if (!BernPolySetSimplex(newpoly, simplex)) {
        return (NULL);
      }
  }

  /* allocate memory for tensordegree and hypercube, and put number in */
  if (tensordim == 0) { /* no tensor part */
    newpoly->tensordegree = NULL;
    newpoly->hypercube = NULL;
  } else {
    newpoly->tensordegree = (int *)malloc(tensordim * sizeof(int));
    if (tensordegree == NULL) {
      for (i = 0; i < tensordim; i++)
        newpoly->tensordegree[i] = 0;
    } else {
      for (i = 0; i < tensordim; i++)
        newpoly->tensordegree[i] = tensordegree[i];
    }

    newpoly->hypercube = (Knote *)malloc(tensordim * sizeof(Knote));
    for (i = 0; i < tensordim; i++) /* for Bern poly only */
      newpoly->hypercube[i].order = 2;

    for (i = 0; i < tensordim; i++)
      (newpoly->hypercube[i]).knote =
          (float *)malloc((newpoly->hypercube[i].order) * sizeof(float));

    if (hypercube != NULL)
      if (!BernPolySetHypercube(newpoly, hypercube)) {
        BernPolyFree(newpoly);
        return (NULL);
      }
  }

  /* allocate memory for coefficients, and put number in */
  max = BernPolyNumCoeff(newpoly);
  newpoly->coeff = (float *)malloc(max * sizeof(float));
  BernPolySetCoeff(newpoly, coeff);

  return (newpoly);
}

/*---------------------------------------------------------------------------
 * BernPolyNumCoeff --
 *     return the number of coefficients of a given poly
 *---------------------------------------------------------------------------*/
int BernPolyNumCoeff(poly)
BernPoly poly;
{
  int max, i;

  errmsg[0] = '\0';
  max = BCombination(BernPolyBaryDegree(poly), BernPolyBaryDim(poly));

  if (BernPolyTensorDim(poly) > 0)
    for (i = 0; i < BernPolyTensorDim(poly); i++)
      max = max * (poly->tensordegree[i] + 1);

  return (max);
}

/*---------------------------------------------------------------------------
 * BernPolySetCoeff --
 *     set coefficient for polynomial
 *---------------------------------------------------------------------------*/
void BernPolySetCoeff(poly, coeff) BernPoly poly;
float *coeff; /* must be in order */
{
  int i, num_coeff;

  errmsg[0] = '\0';
  num_coeff = BernPolyNumCoeff(poly);
  if (coeff == NULL) {
    for (i = 0; i < num_coeff; i++)
      BernPolyCoeff(poly, i) = 0.0;
  } else {
    for (i = 0; i < num_coeff; i++)
      BernPolyCoeff(poly, i) = coeff[i];
  }
}

/*---------------------------------------------------------------------------
 * BernPolySetSimplex --
 *     set the domain for barycentric part
 *--------------------------------------------------------------------------*/
int BernPolySetSimplex(poly, simplex)
BernPoly poly;
float *simplex; /* (barydim+1)x(barydim), vertices of simplex */
{
  int i, j, n;
  double *matrix;

  errmsg[0] = '\0';
  if (BernPolyBaryDim(poly) == 0)
    return 1;

  if (simplex == NULL) {
    sprintf(errmsg, "No simplex given");
    return 0;
  }

  /* check if the simplex is affine independent */
  n = BernPolyBaryDim(poly) + 1;
  if ((matrix = (double *)malloc(n * n * sizeof(double))) == NULL) {
    strcpy(errmsg, "Can not allocate memory");
    return 0;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < (n - 1); i++)
      matrix[i * n + j] = simplex[j * BernPolyBaryDim(poly) + i];

    matrix[(n - 1) * n + j] = 1.0;
  }

  if (!BGaussInverse(matrix, n, 0.0000000001)) {
    /*	if ( !BGaussInverse(matrix, n, 0.000) ) {*/ /* 1998/9/12 changed to
                                                     allow degenerate simplex */
    strcpy(errmsg, "The simplex is not affine independent");
    return 0;
  }

  /* set value of simplex */
  for (i = 0; i <= BernPolyBaryDim(poly); i++)
    for (j = 0; j < BernPolyBaryDim(poly); j++)
      BernPolySimplex(poly, i, j) = simplex[i * BernPolyBaryDim(poly) + j];

  return 1;
}

/*---------------------------------------------------------------------------
 * BernPolySetHypercube --
 *     set the domain for tensor part
 *--------------------------------------------------------------------------*/
int BernPolySetHypercube(poly, hypercube)
BernPoly poly;
float hypercube[][2];
{
  register int i;

  errmsg[0] = '\0';
  if (BernPolyTensorDim(poly) == 0)
    return 1;

  if (hypercube == NULL) {
    sprintf(errmsg, "No hypercube given");
    return 0;
  }

  for (i = 0; i < BernPolyTensorDim(poly); i++)
    if ((hypercube[i][1] - hypercube[i][0]) < 0.0001) {
      strcpy(errmsg, "Invalid hypercube");
      return 0;
    }

  for (i = 0; i < BernPolyTensorDim(poly); i++) {
    BernPolyHypercube(poly, i, 0) = hypercube[i][0];
    BernPolyHypercube(poly, i, 1) = hypercube[i][1];
  }

  return 1;
}

/*---------------------------------------------------------------------------
 * BernZeroPolyCreate --
 *     create a zero poly which has same dimension and domain
 *     as the poly given, and has zero degree
 *--------------------------------------------------------------------------*/
BernPoly BernZeroPolyCreate(poly)
BernPoly poly;
{
  BernPoly new;
  register int i, j;

  new = BernPolyCreate(BernPolyBaryDim(poly), BernPolyTensorDim(poly), 0, NULL,
                       NULL, NULL, NULL);
  if (new == NULL)
    return NULL;

  for (i = 0; i <= BernPolyBaryDim(poly); i++)
    for (j = 0; j < BernPolyBaryDim(poly); j++)
      BernPolySimplex(new, i, j) = BernPolySimplex(poly, i, j);

  for (i = 0; i < BernPolyTensorDim(poly); i++) {
    BernPolyHypercube(new, i, 0) = BernPolyHypercube(poly, i, 0);
    BernPolyHypercube(new, i, 1) = BernPolyHypercube(poly, i, 1);
  }

  return new;
}

/*---------------------------------------------------------------------------
 * BernPolyCopy -- Copy a poly that is the same as given poly
 *---------------------------------------------------------------------------*/
BernPoly BernPolyCopy(poly)
BernPoly poly;
{
  BernPoly p;
  int i, j, max;

  errmsg[0] = '\0';
  p = BernPolyCreate(BernPolyBaryDim(poly), BernPolyTensorDim(poly),
                     BernPolyBaryDegree(poly), poly->tensordegree, NULL, NULL,
                     NULL);
  if (p == NULL)
    return (NULL);

  /* copy coefficients     */
  max = BernPolyNumCoeff(poly);
  for (i = 0; i < max; i++)
    BernPolyCoeff(p, i) = BernPolyCoeff(poly, i);

  /* copy the domain  simplex*/
  for (i = 0; i <= BernPolyBaryDim(poly); i++)
    for (j = 0; j < BernPolyBaryDim(poly); j++)
      BernPolySimplex(p, i, j) = BernPolySimplex(poly, i, j);

  /* copy the domain hypercub*/
  for (i = 0; i < BernPolyTensorDim(poly); i++) {
    BernPolyHypercube(p, i, 0) = BernPolyHypercube(poly, i, 0);
    BernPolyHypercube(p, i, 1) = BernPolyHypercube(poly, i, 1);
  }

  return (p);
}

/*---------------------------------------------------------------------------
 * BernPolyCheckDef --
 *     check if two polynomials have same dimension and same domain
 *--------------------------------------------------------------------------*/
int BernPolyCheckDef(poly1, poly2)
BernPoly poly1, poly2;
{
  register int i, j;
  double d, err;

  errmsg[0] = '\0';

  /* check dimension */
  if ((BernPolyBaryDim(poly1) != BernPolyBaryDim(poly2)) ||
      (BernPolyTensorDim(poly1) != BernPolyTensorDim(poly2))) {
    strcpy(errmsg, "The dimensions of two polys are not match!");
    return 0;
  }

  /* check domain simplex */
  if (BernPolyBaryDim(poly1) > 0) {
    err = 0.0;
    for (i = 0; i <= BernPolyBaryDim(poly1); i++)
      for (j = 0; j < BernPolyBaryDim(poly1); j++) {
        d = BernPolySimplex(poly1, i, j) - BernPolySimplex(poly2, i, j);
        err = err + d * d;
      }

    if (sqrt(err) > 0.00001) {
      sprintf(errmsg, "The simplex of two polys are not match!");
      return 0;
    }
  }

  /* check domain hypercube */
  if (BernPolyTensorDim(poly1) > 0) {
    err = 0.0;
    for (i = 0; i < BernPolyTensorDim(poly1); i++) {
      d = BernPolyHypercube(poly1, i, 0) - BernPolyHypercube(poly2, i, 0);
      err = err + d * d;
      d = BernPolyHypercube(poly1, i, 1) - BernPolyHypercube(poly2, i, 1);
      err = err + d * d;
    }

    if (sqrt(err) > 0.00001) {
      sprintf(errmsg, "The hypercube of two polys are not match!");
      return 0;
    }
  }

  return 1;
}

/*---------------------------------------------------------------------------
 * IsBernZeroPoly -- check if given poly is a zero poly
 *---------------------------------------------------------------------------*/
int IsBernZeroPoly(poly)
BernPoly poly;
{
  int i, max;

  errmsg[0] = '\0';
  max = BernPolyNumCoeff(poly);
  for (i = 0; i < max; i++)
    if (BernPolyCoeff(poly, i) != 0.0)
      return 0;

  /* for zero poly, change degree to zero */
  BernPolyBaryDegree(poly) = 0;
  for (i = 0; i < BernPolyTensorDim(poly); i++)
    BernPolyTensorDegree(poly, i) = 0;

  return 1;
}

/*---------------------------------------------------------------------------
 * BernPolyAdd --
 *     sum two polynomial and the result is a new created polynomial
 *---------------------------------------------------------------------------*/
BernPoly BernPolyAdd(poly1, poly2)
BernPoly poly1, poly2; /* the given polynomials                     */
{
  BernPoly p1, p2, temp;
  int max, i, k;

  errmsg[0] = '\0';

  if (!BernPolyCheckDef(poly1, poly2))
    return NULL;

  /* poly1 is a zero poly */
  if (IsBernZeroPoly(poly1))
    return BernPolyCopy(poly2);

  /* poly2 is a zero poly */
  if (IsBernZeroPoly(poly2))
    return BernPolyCopy(poly1);

  p1 = BernPolyCopy(poly1);
  p2 = BernPolyCopy(poly2);
  if ((p1 == NULL) || (p2 == NULL))
    return NULL;

  /* elevation of barycentric part */
  if (BernPolyBaryDegree(p1) > BernPolyBaryDegree(p2))
    for (i = BernPolyBaryDegree(p2); i < BernPolyBaryDegree(p1); i++) {
      temp = BernPolyDegElevate(p2, 1);
      BernPolyFree(p2);
      p2 = temp;
      if (p2 == NULL)
        return NULL;
    }

  if (BernPolyBaryDegree(p1) < BernPolyBaryDegree(p2))
    for (i = BernPolyBaryDegree(p1); i < BernPolyBaryDegree(p2); i++) {
      temp = BernPolyDegElevate(p1, 1);
      BernPolyFree(p1);
      p1 = temp;
      if (p1 == NULL)
        return NULL;
    }

  /* elevation of tensor part */
  for (k = 0; k < BernPolyTensorDim(poly1); k++) {
    if (BernPolyTensorDegree(p1, k) > BernPolyTensorDegree(p2, k))
      for (i = BernPolyTensorDegree(p2, k); i < BernPolyTensorDegree(p1, k);
           i++) {
        temp = BernPolyDegElevate(p2, BernPolyBaryDim(p2) + k + 1);
        BernPolyFree(p2);
        p2 = temp;
        if (p2 == NULL)
          return NULL;
      }

    if (BernPolyTensorDegree(p1, k) < BernPolyTensorDegree(p2, k))
      for (i = BernPolyTensorDegree(p1, k); i < BernPolyTensorDegree(p2, k);
           i++) {
        temp = BernPolyDegElevate(p1, BernPolyBaryDim(p1) + k + 1);
        BernPolyFree(p1);
        p1 = temp;
        if (p1 == NULL)
          return NULL;
      }
  }

  /* add the coefficients */
  max = BernPolyNumCoeff(p1);
  for (i = 0; i < max; i++)
    BernPolyCoeff(p1, i) = BernPolyCoeff(p1, i) + BernPolyCoeff(p2, i);

  BernPolyFree(p2);
  return (p1);
}

/*---------------------------------------------------------------------------
 * BernPolySub -- the difference of  two polynomials and the result is
 *                a new created polynomial
 *---------------------------------------------------------------------------*/
BernPoly BernPolySub(poly1, poly2)
BernPoly poly1, poly2; /* the given polynomials    */
{
  BernPoly p, sum;
  int max, i;

  errmsg[0] = '\0';

  p = BernPolyCopy(poly2);
  if (p == NULL)
    return NULL;

  max = BernPolyNumCoeff(p);
  for (i = 0; i < max; i++) {
    BernPolyCoeff(p, i) = -BernPolyCoeff(p, i);
  }

  sum = BernPolyAdd(poly1, p);
  BernPolyFree(p);
  return (sum);
}

/*---------------------------------------------------------------------------
 * BernPolyAddConst -- a polynomial add a constant the result is a new
 *                     created polynomial that has the same form as the
 *                     first poly
 *---------------------------------------------------------------------------*/
BernPoly BernPolyAddConst(poly, constant)
BernPoly poly;  /* the given polynomial                     */
float constant; /* the given constant                       */
{
  BernPoly p;
  int max, i;

  errmsg[0] = '\0';
  p = BernPolyCopy(poly);
  if (p == NULL)
    return (NULL);

  max = BernPolyNumCoeff(p);
  for (i = 0; i < max; i++)
    BernPolyCoeff(p, i) = constant + BernPolyCoeff(p, i);

  return (p);
}

/*---------------------------------------------------------------------------
 * BernPolyMultConst -- a polynomial times a constant the result is a new
 *                      created polynomial that has the same form as the
 *                      first poly
 *---------------------------------------------------------------------------*/
BernPoly BernPolyMultConst(poly, constant)
BernPoly poly;  /* the given polynomial                 */
float constant; /* the given constant                   */
{
  BernPoly p;
  int max, i;

  errmsg[0] = '\0';
  p = BernPolyCopy(poly);
  if (p == NULL)
    return (NULL);

  max = BernPolyNumCoeff(p);
  for (i = 0; i < max; i++)
    BernPolyCoeff(p, i) = constant * BernPolyCoeff(p, i);

  return (p);
}

/*---------------------------------------------------------------------------
 * BernPolyDegElevate -- Elevate a polynomial in a given direction. A new poly
 *                       is created that  has higher degree by 1 in the given
 *                       direction
 *---------------------------------------------------------------------------*/
BernPoly BernPolyDegElevate(poly, dir)
BernPoly poly; /* the given polynomial                            */
int dir;       /* the dirction of the elevation that >= 1 and
                  <= dimension of the poly. If the dir is not in
                                              this range, a NULL poly is returned.
                                          If dir <= d1, the dimension of barycentric
                  part,       then the total degree is incresed by 1          */
{
  int d1, d2, int0;
  int i, j, k, ii, i1, i2, max1, max2, max3, max;
  double sum;

  errmsg[0] = '\0';
  d1 = BernPolyBaryDim(poly);
  d2 = BernPolyTensorDim(poly);

  if (dir < 1 || dir > d1 + d2) {
    strcpy(errmsg, "Invalid dirction in degree elevation!");
    return NULL;
  }

  max1 = BCombination(BernPolyBaryDegree(poly), d1);

  max2 = 1;
  for (i = 0; i < d2; i++)
    max2 = max2 * (BernPolyTensorDegree(poly, i) + 1);

  /* elevate barycentric part*/
  if (dir <= d1) {
    BernPoly p;
    int index[10], index1[10];

    p = BernPolyCreate(d1, d2, BernPolyBaryDegree(poly) + 1, poly->tensordegree,
                       NULL, NULL, NULL);
    if (p == NULL)
      return NULL;

    /* form the new coefficients */
    max3 = BCombination(poly->barydegree + 1, d1);
    max = max2 * max3; /* max is the total number of coeffs    */

    for (i = 0; i < max; i++) {
      i1 = i / max2;      /* index in barycentric part of new poly */
      i2 = i - i1 * max2; /* index in tensor part of new poly      */
      sum = 0.0;
      BaryInvIndexMap(i1, index, poly->barydegree + 1, d1);
      ii = poly->barydegree + 1;
      for (k = 0; k < d1; k++)
        ii = ii - index[k];
      index[d1] = ii;
      for (j = 0; j <= d1; j++) {
        for (k = 0; k <= d1; k++)
          index1[k] = index[k];

        if (index1[j] > 0) {
          index1[j] = index1[j] - 1;
          int0 = BaryIndexMap(index1, poly->barydegree, d1);
          sum = sum + index[j] * poly->coeff[int0 * max2 + i2];
        }
      }
      p->coeff[i] = sum / (poly->barydegree + 1);
    }

    /* copy the domain  simplex*/
    for (i = 0; i <= d1; i++)
      for (j = 0; j < d1; j++)
        (p->simplex[i]).vertex[j] = (poly->simplex[i]).vertex[j];

    /* copy the domain hypercub*/
    for (i = 0; i < d2; i++) {
      (p->hypercube[i]).knote[0] = (poly->hypercube[i]).knote[0];
      (p->hypercube[i]).knote[1] = (poly->hypercube[i]).knote[1];
    }

    return (p); /* the barycentric part is finished      */
  }

  /* elevate tensor part*/
  if (dir > d1) {
    BernPoly p;
    int degree2[10], i, j, dir1;
    int intsecond, intfirst, index2[10];

    for (i = 0; i < d2; i++)
      degree2[i] = poly->tensordegree[i];

    dir1 = dir - d1;
    degree2[dir1 - 1] =
        degree2[dir1 - 1] + 1; /* degree2 is the degree of
                                                                          tensor
                                  part of the new poly */

    max3 = 1;
    for (i = 0; i < d2; i++)
      max3 = max3 * (degree2[i] + 1);

    max = max1 * max3; /* max is the total number of coeffs */

    /* create a poly           */
    p = BernPolyCreate(d1, d2, poly->barydegree, degree2, NULL, NULL, NULL);
    if (p == NULL)
      return NULL;

    /* produce the coefficients */
    for (i = 0; i < max; i++) {
      p->coeff[i] = 0.0;
      i1 = i / max3;      /* index in barycentric part of new poly */
      i2 = i - i1 * max3; /* index in tensor part of new poly      */

      /* 1. produce index needed */
      TensorInvIndexMap(i2, index2, degree2, d2);
      intsecond = TensorIndexMap(index2, poly->tensordegree, d2);
      if (index2[dir1 - 1] < degree2[dir1 - 1]) {
        intsecond = TensorIndexMap(index2, poly->tensordegree, d2);
        p->coeff[i] = (1.0 - (float)index2[dir1 - 1] / degree2[dir1 - 1]) *
                      poly->coeff[i1 * max2 + intsecond];
      }

      if (index2[dir1 - 1] > 0) {
        index2[dir1 - 1] = index2[dir1 - 1] - 1;
        intfirst = TensorIndexMap(index2, poly->tensordegree, d2);
        p->coeff[i] = p->coeff[i] + (1.0 + index2[dir1 - 1]) /
                                        degree2[dir1 - 1] *
                                        poly->coeff[i1 * max2 + intfirst];
      }
    }

    /* copy the domain  simplex*/
    for (i = 0; i <= d1; i++) {
      for (j = 0; j < d1; j++)
        (p->simplex[i]).vertex[j] = (poly->simplex[i]).vertex[j];
    }

    /* copy the domain hypercub*/
    for (i = 0; i < d2; i++) {
      (p->hypercube[i]).knote[0] = (poly->hypercube[i]).knote[0];
      (p->hypercube[i]).knote[1] = (poly->hypercube[i]).knote[1];
    }

    return (p); /* tensor part is finished  */
  }
  return (NULL);
}

/*---------------------------------------------------------------------------
 * BernPolyDerivative-- create a new poly that is the partial derivative at
 *                      the given direction for a given poly
 *---------------------------------------------------------------------------*/
BernPoly BernPolyDerivative(poly, dir)
BernPoly poly; /* the given polynomial                      */
int dir;       /* the dirction of the partial that >= 1 and
                  <= dimension of the poly.                 */
{
  int d1, d2;
  int i, i1, i2, max1, max2, max3, max;

  errmsg[0] = '\0';
  d1 = BernPolyBaryDim(poly);
  d2 = BernPolyTensorDim(poly);

  if (dir < 1 || dir > d1 + d2) {
    strcpy(errmsg, "Invalid dirction for partial derivation!");
    return NULL;
  }

  max1 = BCombination(BernPolyBaryDegree(poly), d1);

  max2 = 1;
  for (i = 0; i < d2; i++)
    max2 = max2 * (BernPolyTensorDegree(poly, i) + 1);

  /* partial for tensor part   */
  if (dir > d1) {
    int i, j, dir1, intsecond, intfirst;
    int degree2[10], index2[10];
    BernPoly p;

    for (i = 0; i < d2; i++)
      degree2[i] = poly->tensordegree[i];

    dir1 = dir - d1;

    /* if degree is zero, return zero poly*/
    if (degree2[dir1 - 1] == 0)
      return BernZeroPolyCreate(poly);

    degree2[dir1 - 1] = degree2[dir1 - 1] - 1;
    /* degree2 is partial's degree of tensor
           part of the new poly                  */
    max3 = 1;
    for (i = 0; i < d2; i++)
      max3 = max3 * (degree2[i] + 1);

    max = max1 * max3; /* max is the total number of coeffs */

    /* create a poly           */
    p = BernPolyCreate(d1, d2, poly->barydegree, degree2, NULL, NULL, NULL);
    if (p == NULL)
      return NULL;

    /* produce the coefficients*/
    for (i = 0; i < max; i++) {
      i1 = i / max3;      /* index in barycentric part of new poly */
      i2 = i - i1 * max3; /* index in tensor part of new poly      */

      /* 1. produce index needed */
      TensorInvIndexMap(i2, index2, degree2, d2);
      intsecond = TensorIndexMap(index2, poly->tensordegree, d2);

      index2[dir1 - 1] = index2[dir1 - 1] + 1;
      intfirst = TensorIndexMap(index2, poly->tensordegree, d2);

      /* 2. divided difference   */
      p->coeff[i] = poly->tensordegree[dir1 - 1] *
                    (poly->coeff[i1 * max2 + intfirst] -
                     poly->coeff[i1 * max2 + intsecond]) /
                    ((poly->hypercube[dir1 - 1]).knote[1] -
                     (poly->hypercube[dir1 - 1]).knote[0]);
    }

    /* copy the domain  simplex*/
    for (i = 0; i <= d1; i++)
      for (j = 0; j < d1; j++)
        (p->simplex[i]).vertex[j] = (poly->simplex[i]).vertex[j];

    /* copy the domain hypercub*/
    for (i = 0; i < d2; i++) {
      (p->hypercube[i]).knote[0] = (poly->hypercube[i]).knote[0];
      (p->hypercube[i]).knote[1] = (poly->hypercube[i]).knote[1];
    }

    return (p); /* tensor part is finished       */
  }

  /* partial for  barycentric form*/
  if (dir <= d1) {
    int j, *index, int0, int1;
    float sum;
    double *matrix, dq;
    BernPoly p;

    /* return zero if degree is zero*/
    if (poly->barydegree == 0)
      return BernZeroPolyCreate(poly);

    matrix = (double *)malloc(d1 * d1 * sizeof(double));
    index = (int *)malloc(d1 * d1 * sizeof(int));

    /* form the matrix              */
    for (i = 0; i < d1; i++)
      for (j = 0; j < d1; j++)
        matrix[i * d1 + j] =
            ((poly->simplex[j + 1]).vertex[i] - (poly->simplex[0]).vertex[i]);

    /* inverse the matrix           */
    if (!BGaussInverse(matrix, d1, 0.0000000001)) {
      strcpy(errmsg, "The simplex domain is degenerted!");
      return NULL;
    }

    max3 = BCombination(poly->barydegree - 1, d1);
    max = max2 * max3; /* max is the total number of coeffs  */

    /* create a poly             */
    p = BernPolyCreate(d1, d2, poly->barydegree - 1, poly->tensordegree, NULL,
                       NULL, NULL);
    if (p == NULL)
      return NULL;

    /* produce the coefficients  */
    for (i = 0; i < max; i++) {
      i1 = i / max2;      /* index in barycentric part of new poly */
      i2 = i - i1 * max2; /* index in tensor part of new poly      */
      sum = 0.0;
      for (j = 0; j < d1; j++) {
        BaryInvIndexMap(i1, index + j * d1, poly->barydegree - 1, d1);
        int0 = BaryIndexMap(index + j * d1, poly->barydegree, d1);
        /*
                                index[j*d1+dir-1] = index[j*d1+dir-1] + 1;
        */
        index[j * d1 + j] = index[j * d1 + j] + 1; /* 1995,6,6 changed by xu*/
        int1 = BaryIndexMap(index + j * d1, poly->barydegree, d1);
        dq = poly->barydegree *
             (poly->coeff[int1 * max2 + i2] - poly->coeff[int0 * max2 + i2]);
        /*
                                sum = sum + dq*matrix[(dir-1)*d1 + j];
        */
        sum = sum + dq * matrix[(dir - 1) + j * d1]; /* 1995,6,6 changed by xu*/
      }
      p->coeff[i] = sum;
    }

    /* copy the domain  simplex*/
    for (i = 0; i <= d1; i++)
      for (j = 0; j < d1; j++)
        (p->simplex[i]).vertex[j] = (poly->simplex[i]).vertex[j];

    /* copy the domain hypercub*/
    for (i = 0; i < d2; i++) {
      (p->hypercube[i]).knote[0] = (poly->hypercube[i]).knote[0];
      (p->hypercube[i]).knote[1] = (poly->hypercube[i]).knote[1];
    }

    free(matrix);
    free(index);
    return (p); /* the barycentric part is finished      */
  }
  return (NULL);
}

/*---------------------------------------------------------------------------
 * BernPolyEvaluate --- Evaluate a Bernstein polynimail at a given point.
 *---------------------------------------------------------------------------*/
float BernPolyEvaluate(poly, point)
BernPoly poly; /* the given polynomial                           */
float *point;  /* the evaluating point. The point is in the form
                           (t_1,...,t_d1,s_1,...,s_d2)
                  where t_i is bary part and assume
                 t_0 = 1-t_1-...-t_d1 and s_i is the tensor part */
{
  float *b, result;
  int max, max1, index[4], i, j, k, l, in;

  errmsg[0] = '\0';
  max = BCombination(poly->barydegree, poly->barydim);

  b = (float *)malloc(max * sizeof(float));

  /* Barycentric dimension  0  */
  if (poly->barydim == 0)
    result =
        TensorCastel_k(poly->coeff, point, poly->tensordegree, poly->tensordim);

  /* Barycentric dimension  1  */
  if (poly->barydim == 1) {
    max1 = 1;
    for (i = 0; i < poly->tensordim; i++)
      max1 = max1 * (poly->tensordegree[i] + 1);

    for (i = 0; i <= poly->barydegree; i++)
      b[i] = TensorCastel_k(poly->coeff + i * max1, point + 1,
                            poly->tensordegree, poly->tensordim);

    result = BaryCastel_1(b, point, poly->barydegree);
  }

  /* Barycentric dimension  2  */
  if (poly->barydim == 2) {
    max1 = 1;
    for (i = 0; i < poly->tensordim; i++)
      max1 = max1 * (poly->tensordegree[i] + 1);

    for (i = 0; i <= poly->barydegree; i++) {
      index[0] = i;
      for (j = 0; j <= poly->barydegree - i; j++) {
        index[1] = j;
        in = BaryIndexMap(index, poly->barydegree, poly->barydim);
        b[in] = TensorCastel_k(poly->coeff + in * max1, point + 2,
                               poly->tensordegree, poly->tensordim);
      }
    }

    result = BaryCastel_2(b, point, poly->barydegree);
  }

  /* Barycentric dimension  3  */
  if (poly->barydim == 3) {
    max1 = 1;
    for (i = 0; i < poly->tensordim; i++)
      max1 = max1 * (poly->tensordegree[i] + 1);

    for (i = 0; i <= poly->barydegree; i++) {
      index[0] = i;
      for (j = 0; j <= poly->barydegree - i; j++) {
        index[1] = j;
        for (k = 0; k <= poly->barydegree - i - j; k++) {
          index[2] = k;
          in = BaryIndexMap(index, poly->barydegree, poly->barydim);
          b[in] = TensorCastel_k(poly->coeff + in * max1, point + 3,
                                 poly->tensordegree, poly->tensordim);
        }
      }
    }

    result = BaryCastel_3(b, point, poly->barydegree);
  }

  /* Barycentric dimension  4  */
  if (poly->barydim == 4) {
    max1 = 1;
    for (i = 0; i < poly->tensordim; i++)
      max1 = max1 * (poly->tensordegree[i] + 1);

    for (i = 0; i <= poly->barydegree; i++) {
      index[0] = i;
      for (j = 0; j <= poly->barydegree - i; j++) {
        index[1] = j;
        for (k = 0; k <= poly->barydegree - i - j; k++) {
          index[2] = k;
          for (l = 0; l <= poly->barydegree - i - j - k; l++) {
            index[3] = l;
            in = BaryIndexMap(index, poly->barydegree, poly->barydim);
            b[in] = TensorCastel_k(poly->coeff + in * max1, point + 4,
                                   poly->tensordegree, poly->tensordim);
          }
        }
      }
    }

    result = BaryCastel_4(b, point, poly->barydegree);
  }

  if (poly->barydim > 4) {
    sprintf(errmsg,
            "No functions for evaluating Bernstein polynomial with dim = %d",
            k);

    result = 0.0;
  }

  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * BernPolyFree -- release everything
 *---------------------------------------------------------------------------*/
void BernPolyFree(poly) BernPoly poly;
{
  register int i;

  if (poly == NULL)
    return;

  if (poly->coeff != NULL)
    free(poly->coeff);

  if ((BernPolyBaryDim(poly) > 0) && (poly->simplex != NULL)) {
    for (i = 0; i <= BernPolyBaryDim(poly); i++)
      free((poly->simplex[i]).vertex);

    free(poly->simplex);
  }

  if ((BernPolyTensorDim(poly) > 0) && (poly->hypercube != NULL)) {
    free(poly->tensordegree);
    for (i = 0; i < BernPolyTensorDim(poly); i++)
      free((poly->hypercube[i]).knote);

    free(poly->hypercube);
  }

  free(poly);
}

/*---------------------------------------------------------------------------
 * BCombination -- compute combination number n + k choose k
 *---------------------------------------------------------------------------*/
int BCombination(n, k)
int n, k;
{
  int result, d, i;

  result = 1;
  d = 1;
  for (i = 1; i <= k; i++) {
    result = result * (n + i);
    d = d * i;
  }

  return (result / d);
}

/*---------------------------------------------------------------------------
 * TensorCastel_k -- compute k dimensional Bernstein polynomial in tensor form
 *                   by Casteljau algorithm.
 *---------------------------------------------------------------------------*/
static float TensorCastel_k(bern, t, n, k)
float *bern; /* the coefficients of the given Bernsein        */
float *t;    /* the evaluation point (t_1,t_2,...,t_k)        */
int *n;      /* the degree of the Bernsein (n_1,n_2,..., n_k) */
int k;       /* dimension k >= 0                              */
{
  float result;

  switch (k) {
  case 0:
    result = bern[0];
    break;

  case 1:
    result = TensorCastel_1(bern, t, n);
    break;

  case 2:
    result = TensorCastel_2(bern, t, n);
    break;

  case 3:
    result = TensorCastel_3(bern, t, n);
    break;

  case 4:
    result = TensorCastel_4(bern, t, n);
    break;

  default:
    sprintf(errmsg,
            "No functions for evaluating dimension %d Bernstein polynomial", k);
    result = 0.0;
    break;
  }

  return (result);
}

/*---------------------------------------------------------------------------
 * TensorCastel_1 -- compute one dimensional Bernstein polynomial by Casteljau
 *                   algorithm.
 *---------------------------------------------------------------------------*/
static float TensorCastel_1(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point                      */
int *n;      /* the degree of the Bernsein                */
{
  float *b, result;
  int r, i;

  b = (float *)malloc((n[0] + 1) * sizeof(float));
  for (i = 0; i <= n[0]; i++)
    b[i] = bern[i];

  for (r = 1; r <= n[0]; r++)
    for (i = 0; i <= n[0] - r; i++)
      b[i] = (1.0 - t[0]) * b[i] + t[0] * b[i + 1];

  result = b[0];
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * TensorCastel_2 -- compute 2 dimensional Bernstein polynomial by Casteljau
 *                   algorithm.
 *---------------------------------------------------------------------------*/
static float TensorCastel_2(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2)            */
int *n;      /* the degree of the Bernsein (n_1,n_2)      */
{
  float *b, result;
  int i, max;

  max = 1;
  for (i = 0; i < 2; i++)
    max = max * (n[i] + 1);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  for (i = 0; i <= n[0]; i++)
    b[i] = TensorCastel_1(b + i * (n[1] + 1), t + 1, n + 1);

  result = TensorCastel_1(b, t, n);
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * TensorCastel_3 -- compute 3 dimensional Bernstein polynomial by Casteljau
 *                   algorithm.
 *---------------------------------------------------------------------------*/
static float TensorCastel_3(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2,t_3)        */
int *n;      /* the degree of the Bernsein (n_1,n_2,n_3)  */
{
  float *b, result;
  int i, j, n1, n2, max;

  max = 1;
  for (i = 0; i < 3; i++)
    max = max * (n[i] + 1);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  n1 = n[1] + 1;
  n2 = n[2] + 1;
  for (i = 0; i <= n[0]; i++)
    for (j = 0; j <= n[1]; j++)
      b[i * n1 + j] = TensorCastel_1(b + i * n1 * n2 + j * n2, t + 2, n + 2);

  result = TensorCastel_2(b, t, n);
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * TensorCastel_4 -- compute 4 dimensional Bernstein polynomial by Casteljau
 *                   algorithm.
 *---------------------------------------------------------------------------*/
static float TensorCastel_4(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2,...,t_4)    */
int *n;      /* the degree of the Bernsein (n_1,n_2,...,
                 n_4)                                     */
{
  float *b, result;
  int i, j, k, n1, n2, n3, max;

  max = 1;
  for (i = 0; i < 4; i++)
    max = max * (n[i] + 1);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  n1 = n[1] + 1;
  n2 = n[2] + 1;
  n3 = n[3] + 1;
  for (i = 0; i <= n[0]; i++)
    for (j = 0; j <= n[1]; j++)
      for (k = 0; k <= n[2]; k++)
        b[i * n1 * n2 + j * n2 + k] = TensorCastel_1(
            b + i * n3 * n2 * n1 + j * n2 * n3 + k * n3, t + 3, n + 3);

  result = TensorCastel_3(b, t, n);
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * BaryCastel_k -- compute k dimensional Bernstein polynomial by Casteljau
 *                 algorithm.
 *---------------------------------------------------------------------------*/
static float BaryCastel_k(bern, t, n, k)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2,...,t_k)    */
             /* assume t_0 = 1 - t_1 -...-t_k             */
int n;       /* the degree of the Bernsein                */
int k;       /* dimension k > 0                           */
{
  float result;

  switch (k) {
  case 0:
    result = bern[0];
    break;

  case 1:
    result = BaryCastel_1(bern, t, n);
    break;

  case 2:
    result = BaryCastel_2(bern, t, n);
    break;

  case 3:
    result = BaryCastel_3(bern, t, n);
    break;

  case 4:
    result = BaryCastel_4(bern, t, n);
    break;

  default:
    sprintf(errmsg,
            "No functions for evaluating dimension %d Bernstein polynomial", k);
    result = 0.0;
    break;
  }

  return (result);
}

/*---------------------------------------------------------------------------
 * BaryCastel_4 -- compute 4 dimensional Bernstein polynomial in bartcentric
 *                 form by Casteljau algorithm.
 *---------------------------------------------------------------------------*/
static float BaryCastel_4(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2,...,t_4)    */
int n;       /* the degree of the Bernsein  n             */
{
  float *b, result;
  int i, j, k, l, max, r, s, index[4], index1[4], index2[4], index3[4],
      index4[4];

  max = BCombination(n, 4);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  for (r = 1; r <= n; r++)
    for (i = 0; i <= n - r; i++) {
      index[0] = i;
      for (j = 0; j <= n - r - i; j++) {
        index[1] = j;
        for (k = 0; k <= n - r - i - j; k++) {
          index[2] = k;
          for (l = 0; l <= n - r - i - j - k; l++) {
            index[3] = l;
            for (s = 0; s < 4; s++) {
              index1[s] = index[s];
              index2[s] = index[s];
              index3[s] = index[s];
              index4[s] = index[s];
            }

            index1[0] = index1[0] + 1;
            index2[1] = index2[1] + 1;
            index3[2] = index3[2] + 1;
            index4[3] = index4[3] + 1;

            b[BaryIndexMap(index, n - r, 4)] =
                t[0] * b[BaryIndexMap(index1, n - r + 1, 4)] +
                t[1] * b[BaryIndexMap(index2, n - r + 1, 4)] +
                t[2] * b[BaryIndexMap(index3, n - r + 1, 4)] +
                t[3] * b[BaryIndexMap(index4, n - r + 1, 4)] +
                (1.0 - t[0] - t[1] - t[2] - t[3]) *
                    b[BaryIndexMap(index, n - r + 1, 4)];
          }
        }
      }
    }

  result = b[0];
  free(b);
  return (result);
}
/*---------------------------------------------------------------------------
 * BaryCastel_3 -- compute 3 dimensional Bernstein polynomial in bartcentric
 *                 form by Casteljau algorithm.
 *---------------------------------------------------------------------------*/
static float BaryCastel_3(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2,t_3)        */
int n;       /* the degree of the Bernsein  n             */
{
  float *b, result;
  int i, j, k, max, r, s, index[4], index1[4], index2[4], index3[4];

  max = BCombination(n, 3);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  for (r = 1; r <= n; r++)
    for (i = 0; i <= n - r; i++) {
      index[0] = i;
      for (j = 0; j <= n - r - i; j++) {
        index[1] = j;
        for (k = 0; k <= n - r - i - j; k++) {
          index[2] = k;
          for (s = 0; s < 3; s++) {
            index1[s] = index[s];
            index2[s] = index[s];
            index3[s] = index[s];
          }
          index1[0] = index1[0] + 1;
          index2[1] = index2[1] + 1;
          index3[2] = index3[2] + 1;
          b[BaryIndexMap(index, n - r, 3)] =
              t[0] * b[BaryIndexMap(index1, n - r + 1, 3)] +
              t[1] * b[BaryIndexMap(index2, n - r + 1, 3)] +
              t[2] * b[BaryIndexMap(index3, n - r + 1, 3)] +
              (1.0 - t[0] - t[1] - t[2]) * b[BaryIndexMap(index, n - r + 1, 3)];
        }
      }
    }

  result = b[0];
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * BaryCastel_2 -- compute 2 dimensional Bernstein polynomial in bartcentric
 *                 form by Casteljau algorithm.
 *---------------------------------------------------------------------------*/
static float BaryCastel_2(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1,t_2)            */
int n;       /* the degree of the Bernsein  n             */
{
  float *b, result;
  int i, j, max, r, s, index[4], index1[4], index2[4];

  max = BCombination(n, 2);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  for (r = 1; r <= n; r++)
    for (i = 0; i <= n - r; i++) {
      index[0] = i;
      for (j = 0; j <= n - r - i; j++) {
        index[1] = j;
        for (s = 0; s < 2; s++) {
          index1[s] = index[s];
          index2[s] = index[s];
        }
        index1[0] = index1[0] + 1;
        index2[1] = index2[1] + 1;
        b[BaryIndexMap(index, n - r, 2)] =
            t[0] * b[BaryIndexMap(index1, n - r + 1, 2)] +
            t[1] * b[BaryIndexMap(index2, n - r + 1, 2)] +
            (1.0 - t[0] - t[1]) * b[BaryIndexMap(index, n - r + 1, 2)];
      }
    }

  result = b[0];
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * BaryCastel_1 -- compute 1 dimensional Bernstein polynomial in bartcentric
 *                 form by Casteljau algorithm.
 *---------------------------------------------------------------------------*/
static float BaryCastel_1(bern, t, n)
float *bern; /* the coefficients of the given Bernsein    */
float *t;    /* the evaluation point (t_1)                */
int n;       /* the degree of the Bernsein  n             */
{
  float *b, result;
  int i, max, r, s, index[4], index1[4];

  max = BCombination(n, 1);

  b = (float *)malloc(max * sizeof(float));
  for (i = 0; i < max; i++)
    b[i] = bern[i];

  for (r = 1; r <= n; r++)
    for (i = 0; i <= n - r; i++) {
      index[0] = i;
      for (s = 0; s < 1; s++)
        index1[s] = index[s];
      index1[0] = index1[0] + 1;
      b[BaryIndexMap(index, n - r, 1)] =
          t[0] * b[BaryIndexMap(index1, n - r + 1, 1)] +
          (1.0 - t[0]) * b[BaryIndexMap(index, n - r + 1, 1)];
    }

  result = b[0];
  free(b);
  return (result);
}

/*---------------------------------------------------------------------------
 * BaryIndexMap -- map an index (i1,i2,...,id) to an integer, where the index
 *                 is in barycentric form
 *---------------------------------------------------------------------------*/
int BaryIndexMap(index, m, d)
int index[]; /* the given index                  */
int m;       /* the degree                       */
int d;       /* the dimension                    */
{
  int i, j, k, sum;

  /*	printf("BaryIndexMap:  index = ");
          for (i=0; i < d; i++)  printf("%d  ", index[i]);
          printf("\n\tdegree = %d\n\tdimension = %d\n", m, d);
  */

  k = 0;
  sum = 0;
  for (i = 1; i <= d; i++) {
    for (j = 1; j <= index[i - 1]; j++) {
      k = k + BCombination(m - sum - j + 1, d - i);
    }

    sum = sum + index[i - 1];
  }

  return (k);
}

/*---------------------------------------------------------------------------
 * TensorIndexMap -- map an index (i1,i2,...,id) to an integer
 *---------------------------------------------------------------------------*/
int TensorIndexMap(index, m, d)
int index[]; /* the given index                     */
int m[];     /* the degree that has length d        */
int d;       /* the dimension                       */
{
  int i, j, pro, sum;

  /*	printf("TensorIndexMap:  index = ");
          for (i=0; i < d; i++)  printf("%d  ", index[i]);
          printf("\n\tdegree = ");
          for (i=0; i < d; i++)  printf("%d  ", index[i]);
          printf("\n\tdimension = %d\n", d);
  */
  sum = 0;
  for (i = 1; i <= d; i++) {
    pro = 1;
    for (j = i + 1; j <= d; j++)
      pro = pro * (m[j - 1] + 1);

    sum = sum + index[i - 1] * pro;
  }

  return (sum);
}

/*---------------------------------------------------------------------------
 * TensorInvIndexMap -- map an integer to a tensor form  index (i1,i2,...,id)
 *---------------------------------------------------------------------------*/
void TensorInvIndexMap(ninteger, index, m,
                       d) int ninteger; /* the given integer */
int index[]; /* to be determined index                  */
int m[];     /* degree that has length d                */
int d;       /* dimension                               */
{
  int i, j, sum1, sum;

  sum = ninteger;
  for (i = 1; i <= d; i++) {
    j = d + 1 - i - 1;
    sum1 = sum / (m[j] + 1);
    index[j] = sum - sum1 * (m[j] + 1);
    sum = sum1;
  }
}

/*---------------------------------------------------------------------------
 * BaryInvIndexMap -- map an integer to an barycentric form  index
 *                    (i1,i2,...,id) by using a stupid way--searching
 *---------------------------------------------------------------------------*/
void BaryInvIndexMap(ninteger, index, m,
                     d) int ninteger; /* the given integer */
int index[]; /* to be determined index                  */
int m;       /* degree                                  */
int d;       /* dimension                               */
{
  int i, j, k, l;

  switch (d) {
  case 0:
    break;

  case 1:
    index[0] = ninteger;
    break;

  case 2:
    for (i = 0; i <= m; i++) {
      index[0] = i;
      for (j = 0; j <= m - i; j++) {
        index[1] = j;
        if (BaryIndexMap(index, m, d) == ninteger)
          return;
      }
    }
    break;

  case 3:
    for (i = 0; i <= m; i++) {
      index[0] = i;
      for (j = 0; j <= m - i; j++) {
        index[1] = j;
        for (k = 0; k <= m - i - j; k++) {
          index[2] = k;
          if (BaryIndexMap(index, m, d) == ninteger)
            return;
        }
      }
    }
    break;

  case 4:
    for (i = 0; i <= m; i++) {
      index[0] = i;
      for (j = 0; j <= m - i; j++) {
        index[1] = j;
        for (k = 0; k <= m - i - j; k++) {
          index[2] = k;
          for (l = 0; l <= m - i - j - k; l++) {
            index[3] = l;
            if (BaryIndexMap(index, m, d) == ninteger)
              return;
          }
        }
      }
    }
    break;

  default:
    index[0] = -1;
    sprintf(errmsg, "Invalid dimension %d", d);
    break;
  }
}

/*---------------------------------------------------------------------------
 * BernPolyIndexMap -- map an index:
 *                        (i0, i1, i2, ..., id1, j1, ..., jd2)
 *                     to an integer
 *---------------------------------------------------------------------------*/
int BernPolyIndexMap(poly, index)
BernPoly poly;
int index[]; /* (i0, i1, i2, ..., id1, j1, ..., jd2) */
{
  int i, barymap, tensormap, max, *indexBary, *indexTensor;

  if (BernPolyBaryDim(poly) == 0) {
    indexBary = index;
  } else {
    indexBary = index + 1;
  }
  barymap =
      BaryIndexMap(indexBary, BernPolyBaryDegree(poly), BernPolyBaryDim(poly));

  indexTensor = indexBary + BernPolyBaryDim(poly);
  tensormap = TensorIndexMap(indexTensor, &BernPolyTensorDegree(poly, 0),
                             BernPolyTensorDim(poly));

  max = 1;
  if (barymap > 0) {
    for (i = 0; i < BernPolyTensorDim(poly); i++)
      max = max * (BernPolyTensorDegree(poly, i) + 1);
  }

  return (barymap * max + tensormap);
}

/*---------------------------------------------------------------------------
 * BernPolyInvIndexMap -- map an integer to an index:
 *                        (i0, i1, i2, ..., id1, j1, ..., jd2)
 *---------------------------------------------------------------------------*/
void BernPolyInvIndexMap(poly, m, index) BernPoly poly;
int m, index[];
{
  int i, max, mBary, mTensor, *indexTensor;

  max = 1;
  for (i = 0; i < BernPolyTensorDim(poly); i++) {
    max = max * (BernPolyTensorDegree(poly, i) + 1);
  }

  mBary = m / max;
  mTensor = m - mBary * max;

  if (BernPolyBaryDim(poly) == 0) {
    indexTensor = index;
  } else {
    indexTensor = index + BernPolyBaryDim(poly) + 1;
    BaryInvIndexMap(mBary, index + 1, BernPolyBaryDegree(poly),
                    BernPolyBaryDim(poly));
    index[0] = BernPolyBaryDegree(poly);
    for (i = 1; i <= BernPolyBaryDim(poly); i++)
      index[0] -= index[i];
  }

  TensorInvIndexMap(mTensor, indexTensor, &BernPolyTensorDegree(poly, 0),
                    BernPolyTensorDim(poly));
}

/*----------------------------------------------------------------------------
 * BGaussInverse -- Inverse the matrix A by Gauss elimination methods
 * Return:    0 - if pivot element less than eps
 *            1 - otherwise
 *---------------------------------------------------------------------------*/
int BGaussInverse(a, n, eps)
int n;     /* dimension of the two dimensional matrix a[n][n]      */
double *a, /* matrix to be inversed and also the computing result  */
    eps;   /* control the computation by testing the pivot element */
{
  double max;
  int k, ik, jk, i, j, *z;

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
    ExchangeRowColumn(a, n, k, ik, jk);

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
    ExchangeRowColumn(a, n, k, ik, jk);
  }

  free(z);
  return 1;
}

/*----------------------------------------------------------------------------
 * ExchangeRowColumn - Exchange two rows and two columns of the matrix A
 *---------------------------------------------------------------------------*/
static void ExchangeRowColumn(a, n, k, ik, jk) double *a;
int n, k, ik, jk; /*  exchange the (k,ik) rows, (k,jk) columns  */
{
  double b;
  int j;

  /*  exchange the (k,ik) rows */
  if (ik != k)
    for (j = 0; j < n; j++) {
      b = *(a + ik * n + j);
      *(a + ik * n + j) = *(a + k * n + j);
      *(a + k * n + j) = b;
    }

  /*  exchange the (k,jk) column  */
  if (jk != k)
    for (j = 0; j < n; j++) {
      b = *(a + j * n + jk);
      *(a + j * n + jk) = *(a + j * n + k);
      *(a + j * n + k) = b;
    }
}

/*----------------------------------------------------------------------------
 * ErrorMessage -- return error message
 *---------------------------------------------------------------------------*/
char *ErrorMessage() { return errmsg; }

/*----------------------------------------------------------------------------
 * PrintErrorMessage -- produce error message on the standard error
 *---------------------------------------------------------------------------*/
void PrintErrorMessage() { fprintf(stderr, "%s\n", errmsg); }

/*----------------------------------------------------------------------------
 * BernPolySetCoeffByIndex -- Set a coefficient to a value by given index:
 *                            (i0, i1, i2, ..., id1, j1, ..., jd2)
 *---------------------------------------------------------------------------*/
void BernPolySetCoeffByIndex(poly, index, c) BernPoly poly; /* poly */
int *index; /* index of barycentric part and tensor part */
float c;    /* value of coefficient */
{
  int i;

  i = BernPolyIndexMap(poly, index);
  BernPolyCoeff(poly, i) = c;
}

/*----------------------------------------------------------------------------
 * BernPolyGetCoeffByIndex -- Get the value of a coefficient by given index:
 *                            (i0, i1, i2, ..., id1, j1, ..., jd2)
 *---------------------------------------------------------------------------*/
float BernPolyGetCoeffByIndex(poly, index)
BernPoly poly; /* poly */
int *index;    /* index of barycentric part and tensor part */
{
  int i;

  i = BernPolyIndexMap(poly, index);
  return BernPolyCoeff(poly, i);
}

/*----------------------------------------------------------------------------
 * BernPolyBaryBaseEvaluate -- Evaluate the base function of barycentric
 *                             BB form
 *---------------------------------------------------------------------------*/
float BernPolyBaryBaseEvaluate(var, index, degree, dim)
float var[]; /* variable                                      */
int index[]; /* the power on the var                          */
int degree;  /* degree of the base                            */
int dim;     /* dimension of the base                         */
{
  float result, coeff;
  int i, j;

  result = 1.0;
  for (i = 1; i <= degree; i++) {
    result = result * i;
  }

  coeff = 1.0;
  for (i = 0; i <= dim; i++) {
    for (j = 1; j <= index[i]; j++) {
      result = result * var[i];
      coeff = coeff * j;
    }
  }
  result = result / coeff;
  return (result);
}
