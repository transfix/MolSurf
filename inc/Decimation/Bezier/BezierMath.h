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
 * BezierMath.h  --  The .h  file of the Bezier related functions. All these
 *                   functions deal with mixed form Bernstein base polynomial.
 *                   See the Programmer Manual for the details of the mixed
 *                   form Bernstein base polynomial and related algorithms.
 *                       The dimension could be any in theory, but here we
 *                   assume  it is less than 5. However, it is not difficicult
 *                   to increase the dimension if necessary.
 *                       The C code for these functions are in the file
 *                   BezerLib.c
 *      created by:     Xu Guoliang, Lin Youming
 *      from:           Sep., 29, 1993
 */

#ifndef BEZIER_H
#define BEZIER_H

/* one vertex, for barycentric part */
typedef struct {
  float *vertex;
} OneVertex;

/* for tensor part */
typedef struct {
  int order;    /* (order+1) is the number of base of B-spline,
             if order = tensordegree, the spline
                                     degenerate to Bezier polynomial           */
  float *knote; /* the knotes that has length order          */
} Knote;

typedef struct BernPoly {
  int barydim;        /* dimension for barycentric part       */
  int tensordim;      /* dimension for tensor part            */
  int barydegree;     /* degree of barycentric part           */
  int *tensordegree;  /* a list of degree for tensor part     */
  float *coeff;       /* coefficients                         */
  OneVertex *simplex; /* vertices of the barycentric part that
                 has length (barydim+1)               */
  Knote *hypercube;   /* knotes of tensor part that has
             length tensordim                     */
} *BernPoly;

/* main functions     */
BernPoly BernZeroPolyCreate(/* BernPoly */);
BernPoly BernPolyCreate(/* int, int, int *, float *, float *, float * */);
int BernPolyNumCoeff(/* BernPoly */);
void BernPolySetCoeff(/* BernPoly, float [] */);
int BernPolySetSimplex(/* BernPoly, float [] */);
int BernPolySetHypercube(/* BernPoly, float [] */);
BernPoly BernPolyCopy(/* BernPoly */);
int BernPolyCheckDef(/* BernPoly, BernPoly */);
int IsBernZeroPoly(/* BernPoly */);
BernPoly BernPolyAdd(/* BernPoly, BernPoly */);
BernPoly BernPolySub(/* BernPoly, BernPoly */);
BernPoly BernPolyAddConst(/* BernPoly, float */);
BernPoly BernPolyMultConst(/* BernPoly, float */);
BernPoly BernPolyDegElevate(/* BernPoly, int */);
BernPoly BernPolyDerivative(/* BernPoly, int */);
float BernPolyEvaluate(/* BernPoly, float [] */);
float BernPolyBaryBaseEvaluate(/* float [], int [], int, int */);
void BernPolyFree(/* BernPoly */);
void BernPolySetCoeffByIndex(/* BernPoly, int [], float */);
float BernPolyGetCoeffByIndex(/* BernPoly, int [] */);
int BernPolyIndexMap(/* BernPoly, int [] */);
void BernPolyInvIndexMap(/* BernPoly, int, int [] */);

/* index map function, not recommended */
int BaryIndexMap(/* int [], int, int */);
int TensorIndexMap(/* int [], int, int */);
void TensorInvIndexMap(/* int, int [], int, int */);
void BaryInvIndexMap(/* int, int [], int, int */);

/* error function */
char *ErrorMessage();
void PrintErrorMessage();

/* some utility functions may be used */
int BCombination(/* int, int */);
int BGaussInverse(/* double *, int, double */);

/* macros */
#define BernPolyBaryDim(poly) ((poly)->barydim)
#define BernPolyTensorDim(poly) ((poly)->tensordim)
#define BernPolyBaryDegree(poly) ((poly)->barydegree)
#define BernPolyTensorDegree(poly, i) ((poly)->tensordegree[i])
#define BernPolyCoeff(poly, i) ((poly)->coeff[i])
#define BernPolySimplex(poly, i, j) (((poly)->simplex[i]).vertex[j])
#define BernPolyHypercube(poly, i, j) (((poly)->hypercube[i]).knote[j])

#endif
