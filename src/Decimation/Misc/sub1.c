#include <stdio.h>
#include <math.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <Decimation/Misc/misc.h>
#include <Decimation/Bezier/BezierMath.h>

static int IsophoteIntersectPoint();
extern int ISOPHOTE;
extern float ISOPHOTE_LENGTH;

void NormalizationOf3Vector(n) float n[];
{
  float result;
  result = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (result == 0.0) {
    printf("the length of the vector is zero\n");
    return;
  }
  n[0] = n[0] / result;
  n[1] = n[1] / result;
  n[2] = n[2] / result;
}

void DNormalizationOf3Vector(vector) double vector[3];
{
  double length;

  length = sqrt(vector[0] * vector[0] + vector[1] * vector[1] +
                vector[2] * vector[2]);
  if (length == 0.0) {
    printf("the length of the vector is zero\n");
    return;
  }
  vector[0] = vector[0] / length;
  vector[1] = vector[1] / length;
  vector[2] = vector[2] / length;
}

/*------------------------------------------------------------------------
 * DotProduct -- compute dot product
 *-----------------------------------------------------------------------*/
float DotProduct(p1, p2)
float p1[], p2[];
{
  float result;

  result = p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];

  /* normalize                     */
  /*
  result = result/sqrt((p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]));
  result = result/sqrt((p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]));
  */
  return (result);
}

/*------------------------------------------------------------------------
 * LongDotProduct -- compute dot product
 *-----------------------------------------------------------------------*/
float LongDotProduct(p1, p2, n)
float p1[], p2[];
int n;
{
  int i;
  float result;

  result = 0.0;

  for (i = 0; i < n; i++) {
    result = result + p1[i] * p2[i];
  }

  return (result);
}

/*************************************************************************/
/* InnerProduct                                                          */
/*************************************************************************/
double InnerProduct(u, v)
double u[], v[];
{
  double w;
  w = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  return (w);
}

/*************************************************************************/
/* CrossProduct                                                          */
/* For the given two vectors U and V, compute their unit normal  cross   */
/* product vector  W with  the given direction tan                       */
/*************************************************************************/
void CrossProduct(U, V, tan, W) double U[], V[]; /* the given two vectors */
double tan[]; /* the given direction            */
double W[];   /* the cross-product              */
{
  double Normal, InnerPro;
  /* cross product               */
  W[0] = U[1] * V[2] - U[2] * V[1];
  W[1] = U[2] * V[0] - U[0] * V[2];
  W[2] = U[0] * V[1] - U[1] * V[0];
  /* normal                      */
  Normal = W[0] * W[0] + W[1] * W[1] + W[2] * W[2];
  Normal = sqrt(Normal);

  /* inner product of W and tan  */
  InnerPro = W[0] * tan[0] + W[1] * tan[1] + W[2] * tan[2];
  if (InnerPro < 0.0)
    Normal = -Normal;

  W[0] = W[0] / Normal;
  W[1] = W[1] / Normal;
  W[2] = W[2] / Normal;
}

/*************************************************************************/
/* CrossProductF                                                         */
/* For the given two vectors U and V, compute their unit normal  cross   */
/* product vector  W with  the given direction tan                       */
/*************************************************************************/
void CrossProductF(U, V, tan, W) float U[], V[]; /* the given two vectors */
float tan[]; /* the given direction            */
float W[];   /* the cross-product              */
{
  float Normal, InnerPro;
  /* cross product               */
  W[0] = U[1] * V[2] - U[2] * V[1];
  W[1] = U[2] * V[0] - U[0] * V[2];
  W[2] = U[0] * V[1] - U[1] * V[0];
  /* normal                      */
  Normal = W[0] * W[0] + W[1] * W[1] + W[2] * W[2];
  Normal = sqrt(Normal);

  /* inner product of W and tan  */
  InnerPro = W[0] * tan[0] + W[1] * tan[1] + W[2] * tan[2];
  if (InnerPro < 0.0)
    Normal = -Normal;

  W[0] = W[0] / Normal;
  W[1] = W[1] / Normal;
  W[2] = W[2] / Normal;
}

/*************************************************************************/
/* CrossProduct1                                                         */
/* For the given two vectors U and V, compute their unit normal  cross   */
/* product vector  W                                                     */
/*************************************************************************/
void CrossProduct1(U, V, W) double U[], V[]; /* the given two vectors */
double W[]; /* the cross-product              */
{
  double Normal;
  /* cross product               */
  W[0] = U[1] * V[2] - U[2] * V[1];
  W[1] = U[2] * V[0] - U[0] * V[2];
  W[2] = U[0] * V[1] - U[1] * V[0];
  /* normal                      */
  Normal = W[0] * W[0] + W[1] * W[1] + W[2] * W[2];
  Normal = sqrt(Normal);

  W[0] = W[0] / Normal;
  W[1] = W[1] / Normal;
  W[2] = W[2] / Normal;
}

/* ---------------------------------------------------------------------------*/
void CrossProduct2(u, v, w) float u[], v[], w[];
{
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - v[2] * u[0];
  w[2] = u[0] * v[1] - v[0] * u[1];
}

/* ---------------------------------------------------------------------------*/
void CrossProduct2D(u, v, w) double u[], v[], w[];
{
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - v[2] * u[0];
  w[2] = u[0] * v[1] - v[0] * u[1];
}

/* ---------------------------------------------------------------------------*/
/* ProjectVectorToPlane -- Project a vector to palne defined by               */
/*                (p-point)^T normal = 0                                      */
/* where vector is the given point and normal is the given normal             */
/* provector is the result of the projection                                  */
/* ---------------------------------------------------------------------------*/
void ProjectVectorToPlane(vector, normal, point, provector) float vector[],
    normal[], point[], provector[];
{
  double matrix[3][3], nx, ny, nz;
  float b, c;
  int i;

  nx = normal[0];
  ny = normal[1];
  nz = normal[2];
  RotateMatrix_z(nx, ny, nz, matrix);
  b = 0.0;
  c = 0.0;
  for (i = 0; i < 3; i++) {
    b = b + vector[i] * matrix[i][0];
    c = c + vector[i] * matrix[i][1];
  }

  for (i = 0; i < 3; i++) {
    provector[i] = b * matrix[i][0] + c * matrix[i][1];
  }
}

/* ---------------------------------------------------------------------------*/
/* MaxOfOneVarFunc -- Compute maximum value of one variable function by       */
/*                    search  and refinement method                           */
/* ---------------------------------------------------------------------------*/
void MaxOfOneVarFunc(float (*F)(float), float a, float b, int n, float epsilon,
                     float *xmax, float *Fmax)
/*float (*F)();      */ /* F is the one variable function               */
/*float a,b;         */ /* [a, b] is interval on where the max is taken */
/*int   n;           */ /* The searched number of intervals             */
/*float epsilon;     */ /* error control parameter                      */
/*float *xmax,*Fmax; */ /* The maximal point and maxmal value           */
{
  int i;
  float xl, xm, xr, Fl, Fm, Fr, x, Fx, dx;

  if (a > b) {
    printf("Invalid interval [a,b] in MaxOfOneVarFunc: a, b = %f,%f\n", a, b);
    return;
  }
  dx = b - a;
  if (dx < epsilon) {
    *xmax = 0.5 * (a + b);
    *Fmax = (*F)(*xmax);
    return;
  }
  dx = dx / n;
  xl = a;
  xm = a + dx;
  Fl = (*F)(a);
  *Fmax = Fl;
  *xmax = a;
  Fm = (*F)(a + dx);

  for (i = 2; i <= n; i++) {
    xr = a + i * dx;
    Fr = (*F)(xr);
    if (Fm >= Fl && Fm >= Fr) {

      RefineConvex(F, xl, xm, xr, Fl, Fm, Fr, epsilon, &x, &Fx);
      /*
             MaxOfOneVarFunc(F, xl, xr, n, epsilon, &x,&Fx);
      */
      if (*Fmax < Fx) {
        *Fmax = Fx;
        *xmax = x;
      }
    }
    xl = xm;
    Fl = Fm;
    xm = xr;
    Fm = Fr;
  }
  if (*Fmax < Fr) {
    *Fmax = Fr;
    *xmax = b;
  }
}

/* ---------------------------------------------------------------------------*/
/* RefineConvex -- Compute maximum value of one variable convex function by   */
/*                 three point interpolation  refinement method               */
/* ---------------------------------------------------------------------------*/
void RefineConvex(float (*F)(float), float xl, float xm, float xr, float Fl,
                  float Fm, float Fr, float epsilon, float *xmax, float *Fmax)
/*float (*F)();      */ /* F is the one variable function               */
/*float xl,xm,xr;    */ /* Three points xl < xm < xr                    */
/*float Fl,Fm,Fr;    */ /* The function value on the given three points */
/*float epsilon;     */ /* error control parameter                      */
/*float *xmax,*Fmax; */ /* The maximal point and maxmal value           */
{
  float x, Fx;
  x = ThreePointInterpMax(xl, xm, xr, Fl, Fm, Fr);
  if (x < xl || x > xr) {
    x = xl + 0.618 * (xr - xl);
  }
  Fx = (*F)(x);
  if (x > xm) {
    if (fabs(x - xm) < epsilon || fabs(xr - x) < epsilon) {
      *xmax = x;
      *Fmax = Fx;
      return;
    }
    if (Fx < Fm) {
      RefineConvex(F, xl, xm, x, Fl, Fm, Fx, epsilon, xmax, Fmax);
      return;
    } else {
      RefineConvex(F, xm, x, xr, Fm, Fx, Fr, epsilon, xmax, Fmax);
      return;
    }
  } else {
    if (fabs(x - xl) < epsilon || fabs(xm - x) < epsilon) {
      *xmax = x;
      *Fmax = Fx;
      return;
    }
    if (Fx < Fm) {
      RefineConvex(F, x, xm, xr, Fx, Fm, Fr, epsilon, xmax, Fmax);
      return;
    } else {
      RefineConvex(F, xl, x, xm, Fl, Fx, Fm, epsilon, xmax, Fmax);
      return;
    }
  }
}
/* ---------------------------------------------------------------------------*/
/* ThreePointInterpMax -- Compute maximum value of a poly of degree 2 defined */
/*                        interpolating three points                          */
/* ---------------------------------------------------------------------------*/
float ThreePointInterpMax(float xl, float xm, float xr, float Fl, float Fm,
                          float Fr)
/*float xl,xm,xr;*/ /* Three points xl < xm < xr                    */
/*float Fl,Fm,Fr;*/ /* The function value on the given three points */
{
  float x, af1, af2, bt;

  af1 = (Fm - Fl) / (xm - xl);
  af2 = (Fr - Fl) / (xr - xl);
  bt = (af2 - af1) / (xr - xm);
  if (fabs(bt) < 0.00001) { /* linear function                             */
    if (Fl >= Fr) {
      return (xl);
    } else {
      return (xr);
    }
  }
  x = 0.5 * (xm + xl - af1 / bt);
  /*
  printf("xl,xm,xr,x = %f,%f,%f,%f,%f\n",xl,xm,xr,x,bt);
  */
  return (x);
}

/* ---------------------------------------------------------------------------*/
/* Determinant -- compute the determinant of three vectors [v1,v2,v3] in R^3  */
/* ---------------------------------------------------------------------------*/
float Determinant(v1, v2, v3)
float v1[], v2[], v3[];
{
  float result;

  result = v1[0] * v2[1] * v3[2] + v1[1] * v2[2] * v3[0] +
           v1[2] * v2[0] * v3[1] - v1[2] * v2[1] * v3[0] -
           v1[1] * v2[0] * v3[2] - v1[0] * v2[2] * v3[1];
  return (result);
}

/* ---------------------------------------------------------------------------*/
/* DDeterminant-- compute the determinant of three vectors [v1,v2,v3] in R^3  */
/* ---------------------------------------------------------------------------*/
double DDeterminant(v1, v2, v3)
double v1[], v2[], v3[];
{
  double result;

  result = v1[0] * v2[1] * v3[2] + v1[1] * v2[2] * v3[0] +
           v1[2] * v2[0] * v3[1] - v1[2] * v2[1] * v3[0] -
           v1[1] * v2[0] * v3[2] - v1[0] * v2[2] * v3[1];
  return (result);
}

/* ---------------------------------------------------------------------------*/
/* DDeterminant_Matrix -- compute the determinant of  a matrix M              */
/* ---------------------------------------------------------------------------*/
double DDeterminant_Matrix(M)
double M[];
{
  double result;

  result = M[0] * M[4] * M[8] + M[3] * M[7] * M[2] + M[6] * M[1] * M[5] -
           M[6] * M[4] * M[2] - M[3] * M[1] * M[8] - M[0] * M[7] * M[5];
  return (result);
}

/* ---------------------------------------------------------------------------*/
/* BernPolyFormalBartial2D -- compute the partial derivative of Barycentric   */
/*                            form  Bernstein poly in the formal form         */
/* ---------------------------------------------------------------------------*/
void BernPolyFormalBartial2D(poly, dpoly, degree, dir) BernPoly poly, dpoly;
int degree;
int dir; /* dir = 0,1,2                                  */
{
  int index[3], i, j, k;
  float coeff;

  if (degree == 0) {
    return;
  }

  for (i = 0; i <= degree - 1; i++) {
    for (j = 0; j <= degree - 1 - i; j++) {
      k = degree - 1 - i - j;
      index[0] = i;
      index[1] = j;
      index[2] = k;
      index[dir] = index[dir] + 1;
      coeff = BernPolyGetCoeffByIndex(poly, index);
      coeff = degree * coeff;
      index[dir] = index[dir] - 1;
      BernPolySetCoeffByIndex(dpoly, index, coeff);
    }
  }
}

/* ---------------------------------------------------------------------------*/
/* IsophoteIntersectPoint -- Compute a isophote point on a edge with normal   */
/* ---------------------------------------------------------------------------*/
int IsophoteIntersectPoint(p1, p2, n1, n2, L, c, p, n)
float p1[], p2[], n1[], n2[], L[], c, p[], n[];
{
  int i, j, number;
  float n21[3], ln1, ln2, c2;
  double A, B, C, root[2];

  c2 = c * c;
  for (i = 0; i < 3; i++) {
    n21[i] = n2[i] - n1[i];
  }
  ln2 = DotProduct(L, n2);
  ln1 = DotProduct(L, n1);

  A = (ln2 - ln1) * (ln2 - ln1) - DotProduct(n21, n21) * c2;
  B = ln1 * (ln2 - ln1) - c2 * DotProduct(n1, n21);
  B = B + B;
  C = ln1 * ln1 - c2 * DotProduct(n1, n1);

  FindRoots2(A, B, C, root, &number);
  if (number > 0) {
    for (i = 0; i < number; i++) {
      if (root[i] >= -0.000001 && root[i] <= 1.000001) {
        ln1 = root[i];
        for (j = 0; j < 3; j++) {
          p[j] = (1.0 - ln1) * p1[j] + ln1 * p2[j];
          n[j] = (1.0 - ln1) * n1[j] + ln1 * n2[j];
        }
        return (1);
      }
    }
  }

  return (0);
}

void DrawPolygon_Line(k, vertex, normal) int k;
Vertex *vertex;
Normal *normal;
{
  int i, l;
  float p1[3], p2[3];

  // DrawPolygon(k, vertex, normal);

  for (i = 0; i < k; i++) {
    l = (i + 1) % k;

    p1[0] = vertex[i][0] + 0.005 * normal[i][0];
    p1[1] = vertex[i][1] + 0.005 * normal[i][1];
    p1[2] = vertex[i][2] + 0.005 * normal[i][2];

    p2[0] = vertex[l][0] + 0.005 * normal[l][0];
    p2[1] = vertex[l][1] + 0.005 * normal[l][1];
    p2[2] = vertex[l][2] + 0.005 * normal[l][2];
    // LineAbs(p1[0],p1[1],p1[2], p2[0],p2[1],p2[2]);
  }
  for (i = 0; i < k; i++) {
    l = (i + 1) % k;
    p1[0] = vertex[i][0] - 0.005 * normal[i][0];
    p1[1] = vertex[i][1] - 0.005 * normal[i][1];
    p1[2] = vertex[i][2] - 0.005 * normal[i][2];

    p2[0] = vertex[l][0] - 0.005 * normal[l][0];
    p2[1] = vertex[l][1] - 0.005 * normal[l][1];
    p2[2] = vertex[l][2] - 0.005 * normal[l][2];
    // LineAbs(p1[0],p1[1],p1[2], p2[0],p2[1],p2[2]);
  }
}

/* ---------------------------------------------------------------------------*/
/* XuDrawPolygon -- Compute a isophote point on a edge with normal            */
/* ---------------------------------------------------------------------------*/
void XuDrawPolygon(k, vertex, normal) Vertex *vertex;
Normal *normal;
int k;
{
  float p12[3], p13[3], p23[3], n12[3], n13[3], n23[3], L[3], c, p1[3], p2[3],
      p3[3], n1[3], n2[3], n3[3], d;
  int i12, i13, i23, i;

  if (ISOPHOTE == 0) {
    // DrawPolygon(k, vertex, normal);
    return;
  }
  if (k != 3) {
    return;
  }

  if (ISOPHOTE == 1) {
    for (i = 0; i < 3; i++) {
      p1[i] = vertex[0][i];
      p2[i] = vertex[1][i];
      p3[i] = vertex[2][i];
      n1[i] = normal[0][i];
      n2[i] = normal[1][i];
      n3[i] = normal[2][i];
    }

    L[0] = 0.57735027;
    L[1] = 0.57735027;
    L[2] = 0.57735027;

    for (i = 1; i < 10; i++) {
      c = i / 9.0;

      i12 = IsophoteIntersectPoint(p1, p2, n1, n2, L, c, p12, n12);
      i13 = IsophoteIntersectPoint(p1, p3, n1, n3, L, c, p13, n13);
      i23 = IsophoteIntersectPoint(p2, p3, n2, n3, L, c, p23, n23);

      /*
      printf("i12, i13, i23 = %d, %d, %d\n", i12, i13, i23);
  */

      d = 0.01;
      if (i12 == 1 && i13 == 1) {
        // LineAbs(p12[0]+d*n12[0],p12[1]+d*n12[1],p12[2]+d*n12[2],
        //         p13[0]+d*n13[0],p13[1]+d*n13[1],p13[2]+d*n13[2]);

        // LineAbs(p12[0]-d*n12[0],p12[1]-d*n12[1],p12[2]-d*n12[2],
        //         p13[0]-d*n13[0],p13[1]-d*n13[1],p13[2]-d*n13[2]);

        ISOPHOTE_LENGTH =
            ISOPHOTE_LENGTH + sqrt((p12[0] - p13[0]) * (p12[0] - p13[0]) +
                                   (p12[1] - p13[1]) * (p12[1] - p13[1]) +
                                   (p12[2] - p13[2]) * (p12[2] - p13[2]));
      }
      if (i12 == 1 && i23 == 1) {
        // LineAbs(p12[0]+d*n12[0],p12[1]+d*n12[1],p12[2]+d*n12[2],
        //         p23[0]+d*n23[0],p23[1]+d*n23[1],p23[2]+d*n23[2]);

        // LineAbs(p12[0]-d*n12[0],p12[1]-d*n12[1],p12[2]-d*n12[2],
        //         p23[0]-d*n23[0],p23[1]-d*n23[1],p23[2]-d*n23[2]);

        ISOPHOTE_LENGTH =
            ISOPHOTE_LENGTH + sqrt((p12[0] - p23[0]) * (p12[0] - p23[0]) +
                                   (p12[1] - p23[1]) * (p12[1] - p23[1]) +
                                   (p12[2] - p23[2]) * (p12[2] - p23[2]));
      }
      if (i13 == 1 && i23 == 1) {
        // LineAbs(p13[0]+d*n13[0],p13[1]+d*n13[1],p13[2]+d*n13[2],
        //         p23[0]+d*n23[0],p23[1]+d*n23[1],p23[2]+d*n23[2]);

        // LineAbs(p13[0]-d*n13[0],p13[1]-d*n13[1],p13[2]-d*n13[2],
        //         p23[0]-d*n23[0],p23[1]-d*n23[1],p23[2]-d*n23[2]);

        ISOPHOTE_LENGTH =
            ISOPHOTE_LENGTH + sqrt((p13[0] - p23[0]) * (p13[0] - p23[0]) +
                                   (p13[1] - p23[1]) * (p13[1] - p23[1]) +
                                   (p13[2] - p23[2]) * (p13[2] - p23[2]));
      }
    } /*end i loop            */
  } /*end if ISOPHOTE == 1  */
}

/*************************************************************************/
/* RotateMatrix -- Compute a rotate matrix M that rotate  vector n1 to n2*/
/*                 along the vector n3, that is, M * n1 = n2             */
/*                 where n1, n2, n3  are normalized vectors and n3 is    */
/*                 perpendicular to n1 and n2.                           */
/*************************************************************************/
void RotateMatrix(n1, n2, n3, matrix) double n1[], n2[],
    n3[];            /* the  given vectors                 */
double matrix[3][3]; /* the rotating matrix                */
{
  int i;
  double matrix1[3][3], matrix2[3][3], n13[3], n23[3];

  CrossProduct1(n1, n3, n13);
  CrossProduct1(n2, n3, n23);

  matrix1[0][0] = n2[0];
  matrix1[1][0] = n2[1];
  matrix1[2][0] = n2[2];

  matrix1[0][1] = n23[0];
  matrix1[1][1] = n23[1];
  matrix1[2][1] = n23[2];

  matrix1[0][2] = n3[0];
  matrix1[1][2] = n3[1];
  matrix1[2][2] = n3[2];

  matrix2[0][0] = n1[0];
  matrix2[0][1] = n1[1];
  matrix2[0][2] = n1[2];

  matrix2[1][0] = n13[0];
  matrix2[1][1] = n13[1];
  matrix2[1][2] = n13[2];

  matrix2[2][0] = n3[0];
  matrix2[2][1] = n3[1];
  matrix2[2][2] = n3[2];

  for (i = 0; i < 3; i++) {
    matrix[i][0] = matrix1[i][0] * matrix2[0][0] +
                   matrix1[i][1] * matrix2[1][0] +
                   matrix1[i][2] * matrix2[2][0];

    matrix[i][1] = matrix1[i][0] * matrix2[0][1] +
                   matrix1[i][1] * matrix2[1][1] +
                   matrix1[i][2] * matrix2[2][1];

    matrix[i][2] = matrix1[i][0] * matrix2[0][2] +
                   matrix1[i][1] * matrix2[1][2] +
                   matrix1[i][2] * matrix2[2][2];
  }
}

/*************************************************************************/
/* RotateMatrix_z -- Rotate the vector (nx,ny,nz)^T to be z-axis         */
/*                 That is, matrix^T*(nx,ny,nz)^T = [0,0, |n|]           */
/*************************************************************************/
void RotateMatrix_z(nx, ny, nz, matrix) double nx, ny,
    nz;              /* the  normal  at the given point    */
double matrix[3][3]; /* the rotating matrix                */
{
  double c1, c2, s1, s2, normal, normalz;

  normal = sqrt(nx * nx + ny * ny + nz * nz);
  normalz = sqrt(nx * nx + ny * ny);
  c1 = nz / normal;
  c2 = -ny / normalz;
  s1 = normalz / normal;
  s2 = nx / normalz;

  if (normalz < 0.001) {
    matrix[0][0] = 1.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = 0.0;
    matrix[1][1] = 1.0;
    matrix[1][2] = 0.0;
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 1.0;
  }
  if (normalz >= 0.001) {
    matrix[0][0] = c2;
    matrix[0][1] = -c1 * s2;
    matrix[0][2] = s1 * s2;
    matrix[1][0] = s2;
    matrix[1][1] = c1 * c2;
    matrix[1][2] = -s1 * c2;
    matrix[2][0] = 0.0;
    matrix[2][1] = s1;
    matrix[2][2] = c1;
  }
}

/*-----------------------------------------------------------------------------
 * TriangleNormal--compute the normal of a triangle, by right-handed rule
 *---------------------------------------------------------------------------*/
void TriangleNormal(p1, p2, p3, normal) float p1[], p2[],
    p3[];       /* vertices of the triangle                */
float normal[]; /* the computed normal                     */
{
  float x1, y1, z1, x2, y2, z2;

  x1 = p1[0] - p3[0];
  y1 = p1[1] - p3[1];
  z1 = p1[2] - p3[2];

  x2 = p2[0] - p3[0];
  y2 = p2[1] - p3[1];
  z2 = p2[2] - p3[2];

  normal[0] = y1 * z2 - y2 * z1;
  normal[1] = x2 * z1 - x1 * z2;
  normal[2] = x1 * y2 - x2 * y1;
  x1 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
  x1 = sqrt(x1);
  if (x1 != 0.0) {
    normal[0] = normal[0] / x1;
    normal[1] = normal[1] / x1;
    normal[2] = normal[2] / x1;
  } else {
    printf("normal = %f , %f, %f\n", normal[0], normal[1], normal[2]);
  }
}

/*---------------------------------------------------------------------------
SortingInteger -- Sorting an integer  array from small to big
----------------------------------------------------------------------------*/
void SortingInteger(point, index, number) int point[];
int index[];
int number;
{
  int n, temp, gap, i, j;

  n = number;
  for (i = 0; i < n; i++) {
    index[i] = i;
  }

  for (gap = n / 2; gap > 0; gap /= 2)
    for (i = gap; i < n; i++) {
      for (j = i - gap; j >= 0 && point[index[j]] > point[index[j + gap]];
           j -= gap) {
        temp = index[j];
        index[j] = index[j + gap];
        index[j + gap] = temp;
      }
    }
}

/*---------------------------------------------------------------------------
SortingDouble -- Sorting a double array from small to big
----------------------------------------------------------------------------*/
void SortingDouble(point, index, number) double point[];
int index[];
int number;
{
  int n, temp, gap, i, j;
  n = number;
  /*
  if (n > 256) n = 256;
  */
  for (i = 0; i < n; i++) {
    index[i] = i;
  }
  for (gap = n / 2; gap > 0; gap /= 2)
    for (i = gap; i < n; i++) {
      for (j = i - gap; j >= 0 && point[index[j]] > point[index[j + gap]];
           j -= gap) {
        temp = index[j];
        index[j] = index[j + gap];
        index[j + gap] = temp;
      }
    }
}

/*---------------------------------------------------------------------------
ReadHowmanyComponents -- Determine how mant component has of a file
----------------------------------------------------------------------------*/
int ReadHowmanyComponents(filename)
char *filename;
{
  char str[256];
  FILE *fp;
  float f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;
  int nvert, nface, component;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return (0);
  }

  fgets(str, 256, fp);
  sscanf(str, "%d %d", &nvert, &nface);

  fgets(str, 256, fp);
  component = sscanf(str, "%f %f %f %f %f %f %f %f %f %f %f %f\n", &f0, &f1,
                     &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11);

  fclose(fp);
  return (component);
}

/*---------------------------------------------------------------------------
DistanceOfTwoLines -- compute the distance of two lines in the space
                      defined by p1(1-s) + p2 s , q1(1-t) + q2 t
----------------------------------------------------------------------------*/
float DistanceOfTwoLines(p1, p2, q1, q2, s, t)
float p1[], p2[], q1[], q2[], *s, *t;
{
  int i;
  float result, pp[3], qq[3], pq[3], a, b, c, b1, b2;

  for (i = 0; i < 3; i++) {
    pp[i] = p2[i] - p1[i];
    qq[i] = q1[i] - q2[i];
    pq[i] = p1[i] - q1[i];
  }
  a = DotProduct(pp, pp);
  b = DotProduct(pp, qq);
  c = DotProduct(qq, qq);

  b1 = -DotProduct(pq, pp);
  b2 = -DotProduct(pq, qq);

  *t = a * c - b * b;

  if (fabs(*t) < 0.000001) { /* two line parallel */
    *t = 0.0;
    *s = b1 / a;
    result = DotProduct(pq, pq) + *s * *s * a - 2 * *s * b1;
  } else {

    *s = (b1 * c - b2 * b) / (*t);
    *t = (b2 * a - b1 * b) / (*t);

    result = DotProduct(pq, pq) + *s * *s * a + *t * *t * c - 2 * *s * b1 -
             2 * *t * b2 + 2 * *s * *t * b;
  }
  return (sqrt(result));
}

/*---------------------------------------------------------------------------
SizeOfTriangle  -- the size of a triangle, defined by the average of the
                   lengths of its three edges
----------------------------------------------------------------------------*/
float SizeOfTriangle(v1, v2, v3)
float v1[], v2[], v3[];
{
  float result, v12[3], v13[3], v23[3];

  v12[0] = v1[0] - v2[0];
  v12[1] = v1[1] - v2[1];
  v12[2] = v1[2] - v2[2];

  v13[0] = v1[0] - v3[0];
  v13[1] = v1[1] - v3[1];
  v13[2] = v1[2] - v3[2];

  v23[0] = v2[0] - v3[0];
  v23[1] = v2[1] - v3[1];
  v23[2] = v2[2] - v3[2];

  result = sqrt(DotProduct(v12, v12)) + sqrt(DotProduct(v13, v13)) +
           sqrt(DotProduct(v23, v23));

  result = result / 3.0;
  return (result);
}

/*---------------------------------------------------------------------------
AreaOfTriangle  -- the area of a triangle
----------------------------------------------------------------------------*/
float AreaOfTriangle(v1, v2, v3)
float v1[], v2[], v3[];
{
  float a, b, c, p, result, v12[3], v13[3], v23[3];

  v12[0] = v1[0] - v2[0];
  v12[1] = v1[1] - v2[1];
  v12[2] = v1[2] - v2[2];

  v13[0] = v1[0] - v3[0];
  v13[1] = v1[1] - v3[1];
  v13[2] = v1[2] - v3[2];

  v23[0] = v2[0] - v3[0];
  v23[1] = v2[1] - v3[1];
  v23[2] = v2[2] - v3[2];

  a = sqrt(DotProduct(v12, v12));
  b = sqrt(DotProduct(v13, v13));
  c = sqrt(DotProduct(v23, v23));

  p = 0.5 * (a + b + c);
  p = p * (p - a) * (p - b) * (p - c);
  if (p < 0.0)
    p = 0.0;

  result = sqrt(p);
  return (result);
}

/*---------------------------------------------------------------------------
LinearSystemOrder2  -- Solving a linear system of order 2
                       [ c1  c2 ]                 [y1]
                       [ c3  c4 ]  is the matrix  [y2] - is the right-handed
----------------------------------------------------------------------------*/
int LinearSystemOrder2(float c1, float c2, float c3, float c4, float y1,
                       float y2, float *s, float *t) {
  float det;

  det = c1 * c4 - c2 * c3;

  if (fabs(det) < 0.000001) {
    return (0);
  }

  *s = (y1 * c4 - y2 * c2) / det;
  *t = (c1 * y2 - c3 * y1) / det;
  return (1);
}

/*---------------------------------------------------------------------------
DLinearSystemOrder2  -- Solving a linear system of order 2
                       [ c1  c2 ]                 [y1]
                       [ c3  c4 ]  is the matrix  [y2] - is the right-handed
----------------------------------------------------------------------------*/
int DLinearSystemOrder2(c1, c2, c3, c4, y1, y2, s, t)
double c1, c2, c3, c4, y1, y2, *s, *t;
{
  double det;

  det = c1 * c4 - c2 * c3;

  if (fabs(det) < 0.000000000001) {
    return (0);
  }

  *s = (y1 * c4 - y2 * c2) / det;
  *t = (c1 * y2 - c3 * y1) / det;
  return (1);
}

/*---------------------------------------------------------------------------
Angle_Of_Three_Points -- compute the angle of three points p1,p2,p3. That
                         is the angle between p1-p2 and p3-p2
----------------------------------------------------------------------------*/
float Angle_Of_Three_Points(p1, p2, p3)
float p1[], p2[], p3[];
{
  float p12[3], p32[3], result, length1, length3;

  p12[0] = p1[0] - p2[0];
  p12[1] = p1[1] - p2[1];
  p12[2] = p1[2] - p2[2];

  p32[0] = p3[0] - p2[0];
  p32[1] = p3[1] - p2[1];
  p32[2] = p3[2] - p2[2];

  length1 = DotProduct(p12, p12);
  if (length1 < 0.0000001)
    return (0.0);
  length3 = DotProduct(p32, p32);
  if (length3 < 0.0000001)
    return (0.0);

  length1 = sqrt(length1);
  length3 = sqrt(length3);

  result = DotProduct(p12, p32) / (length1 * length3);

  if (result > 1.0)
    result = 1.0;
  if (result < -1.0)
    result = -1.0;

  result = acos(result);

  return (result);
}

/*---------------------------------------------------------------------------
Edge_Normal_Of_Triangle -- compute the normal of triangle [p1p2p3] on the
                           edge [p1p2], points to p3
----------------------------------------------------------------------------*/
void Edge_Normal_Of_Triangle(p1, p2, p3, normal) float p1[], p2[], p3[],
    normal[];
{
  float p12[3], p32[3], af, length1;

  p12[0] = p1[0] - p2[0];
  p12[1] = p1[1] - p2[1];
  p12[2] = p1[2] - p2[2];

  p32[0] = p3[0] - p2[0];
  p32[1] = p3[1] - p2[1];
  p32[2] = p3[2] - p2[2];

  length1 = DotProduct(p12, p12);
  if (length1 < 0.0000001) {
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    return;
  }

  af = DotProduct(p12, p32) / length1;
  normal[0] = p32[0] - af * p12[0];
  normal[1] = p32[1] - af * p12[1];
  normal[2] = p32[2] - af * p12[2];

  length1 = sqrt(DotProduct(normal, normal));

  if (length1 < 0.0000001) {
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    return;
  }
  normal[0] = normal[0] / length1;
  normal[1] = normal[1] / length1;
  normal[2] = normal[2] / length1;
}

/*---------------------------------------------------------------------------
Barycentric_Coord_3D_Triangle -- compute the barycentric coordinate of a
                                 point in term a 3D triangle in space.
----------------------------------------------------------------------------*/
void Barycentric_Coord_3D_Triangle(p, p1, p2, p3, b123) float p[], p1[], p2[],
    p3[], b123[];
{
  int i, info;
  float pp3[3], pp13[3], pp23[3];
  float c1, c2, c3, c4, y1, y2;

  for (i = 0; i < 3; i++) {
    pp3[i] = p[i] - p3[i];
    pp13[i] = p1[i] - p3[i];
    pp23[i] = p2[i] - p3[i];
  }
  c1 = DotProduct(pp13, pp13);
  c2 = DotProduct(pp23, pp13);
  c3 = c2;
  c4 = DotProduct(pp23, pp23);
  y1 = DotProduct(pp3, pp13);
  y2 = DotProduct(pp3, pp23);
  if (c1 > c4) {
    c2 = c2 / c1;
    c3 = c3 / c1;
    c4 = c4 / c1;
    y1 = y1 / c1;
    y2 = y2 / c1;
    c1 = 1.0;
  } else {
    c1 = c1 / c4;
    c2 = c2 / c4;
    c3 = c3 / c4;
    y1 = y1 / c4;
    y2 = y2 / c4;
    c4 = 1.0;
  }

  info = LinearSystemOrder2(c1, c2, c3, c4, y1, y2, b123, b123 + 1);
  if (info == 0) {
    printf("The triangle is degenerate, in Barycentric_Coord_3D_Triangle\n");
    b123[0] = 0.3333333;
    b123[1] = 0.3333333;
    b123[2] = 0.3333333;
    printf("c1,c2,c3,c4,y1,y2 = %f,%f, %f, %f,%f, %f\n", c1, c2, c3, c4, y1,
           y2);
    printf("b123 = %f,%f, %f\n", b123[0], b123[1], b123[2]);
    printf("p = %f,%f,%f\n", p[0], p[1], p[2]);
    return;
  }
  b123[2] = 1.0 - b123[0] - b123[1];
}

/*---------------------------------------------------------------------------
  OneD_degree_elevate1 -- Degree elevation once for one dimension poly
----------------------------------------------------------------------------*/
void OneD_degree_elevate1(degree, incoeff, outcoeff) int degree;
float incoeff[], outcoeff[];
{
  int i, degree1;

  degree1 = degree + 1;

  outcoeff[0] = incoeff[0];
  outcoeff[degree1] = incoeff[degree];

  for (i = 1; i <= degree; i++) {
    outcoeff[i] = i * incoeff[i - 1] + (degree1 - i) * incoeff[i];
    outcoeff[i] = outcoeff[i] / degree1;
  }
}

/*---------------------------------------------------------------------------
  OneD_degree_elevate -- Degree elevation  for one dimension poly
----------------------------------------------------------------------------*/
void OneD_degree_elevate(indegree, incoeff, outdegree, outcoeff) int indegree,
    outdegree;
float incoeff[], outcoeff[];
{
  int i, j;
  float *tem;

  if (outdegree <= indegree) {
    for (i = 0; i <= indegree; i++) {
      outcoeff[i] = incoeff[i];
    }
    return;
  }

  tem = (float *)malloc((outdegree + 1) * sizeof(float));

  for (j = 0; j <= indegree; j++) {
    tem[j] = incoeff[j];
  }

  for (i = indegree; i < outdegree; i++) {
    OneD_degree_elevate1(i, tem, outcoeff);
    if (i < outdegree - 1) {
      for (j = 0; j <= i + 1; j++) {
        tem[j] = outcoeff[j];
      }
    }
  }

  free(tem);
}

/*---------------------------------------------------------------------------
  Intersection_Two_Lines  -- Compute the intersection point of two lines
                             (1-s)*p1 + s*p2,  (1-t)*q1 + t*q2
----------------------------------------------------------------------------*/
int Intersection_Two_Lines(p1, p2, q1, q2, s, t, point)
float p1[], p2[], q1[], q2[], *s, *t, point[];
{
  float c1, c2, c3, c4, y1, y2;
  int info;

  c1 = p2[0] - p1[0];
  c2 = q1[0] - q2[0];
  c3 = p2[1] - p1[1];
  c4 = q1[1] - q2[1];
  y1 = q1[0] - p1[0];
  y2 = q1[1] - p1[1];

  info = LinearSystemOrder2(c1, c2, c3, c4, y1, y2, s, t);
  if (info == 1) {
    point[0] = (1.0 - *s) * p1[0] + *s * p2[0];
    point[1] = (1.0 - *s) * p1[1] + *s * p2[1];
  }

  return (info);
}

/*---------------------------------------------------------------------------
Distance_Two_points -- Distance of two 2D points
----------------------------------------------------------------------------*/
float Distance_Two_points(p1, p2)
float p1[], p2[];
{
  float result;

  result = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                (p1[1] - p2[1]) * (p1[1] - p2[1]));

  return (result);
}

/*---------------------------------------------------------------------------
Distance_Two_points -- Distance of two 3D points
----------------------------------------------------------------------------*/
float Distance_Two_points3D(p1, p2)
float p1[], p2[];
{
  float result;

  result = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                (p1[1] - p2[1]) * (p1[1] - p2[1]) +
                (p1[2] - p2[2]) * (p1[2] - p2[2]));

  return (result);
}

/* -------------------------------------------------------------------------
Search in a sequence, assume the the sequence is in the increasing order
by bisection approach for n > 4.  For n <= 4, using linear search. The
return is the index of the find object in the sequence.
----------------------------------------------------------------------------*/
int SearchingInSequence(int sequence[], int n, int object) {
  int i, left, right, middle;

  if (n <= 4) {
    /* linear search   */
    for (i = 0; i < n; i++) {
      if (sequence[i] == object)
        return (i);
    }
    return (-1);
  }

  left = 0;
  right = n - 1;

  /* bisection search  */
  while (right >= left) {
    middle = (left + right) >> 1;
    if (sequence[middle] == object)
      return (middle);
    if (sequence[middle] < object)
      left = middle + 1;
    else
      right = middle - 1;
  }
  return (-1);
}

/* ---------------------------------------------------------------------------
 */
/* OutCenter -- Compute the out center of a triangle */
/* ---------------------------------------------------------------------------
 */
int CircumCenter(p1, p2, p3, center, barycord)
float p1[], p2[], p3[], center[], barycord[];
{
  double matrix[9], b[3], solution[3];
  float p13[3], p12[3], ap13[3], ap12[3];
  int i, info;

  for (i = 0; i < 3; i++) {
    p13[i] = p1[i] - p3[i];
    p12[i] = p1[i] - p2[i];
    ap13[i] = 0.5 * (p1[i] + p3[i]);
    ap12[i] = 0.5 * (p1[i] + p2[i]);
    matrix[6 + i] = 1.0;
  }
  matrix[0] = DotProduct(p13, p1);
  matrix[1] = DotProduct(p13, p2);
  matrix[2] = DotProduct(p13, p3);
  matrix[3] = DotProduct(p12, p1);
  matrix[4] = DotProduct(p12, p2);
  matrix[5] = DotProduct(p12, p3);
  b[0] = DotProduct(p13, ap13);
  b[1] = DotProduct(p12, ap12);
  ;
  b[2] = 1.0;

  gaussinverse(matrix, 3, 0.000000000001, &info);

  if (info == 0) {
    printf("The three point is colinear and the matrix  is singular!\n");
    return (0);
  }

  MatrixTimesMatrix(matrix, b, solution, 3, 3, 1);
  barycord[0] = solution[0];
  barycord[1] = solution[1];
  barycord[2] = solution[2];

  for (i = 0; i < 3; i++) {
    center[i] = p1[i] * barycord[0] + p2[i] * barycord[1] + p3[i] * barycord[2];
  }

  return (1);
}
