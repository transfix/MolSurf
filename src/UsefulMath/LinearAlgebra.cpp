/*
  Copyright 2011 The University of Texas at Austin

        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MolSurf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <UsefulMath/LinearAlgebra.h>

using CCVOpenGLMath::LinearAlgebra;
using CCVOpenGLMath::Vector;

LinearAlgebra::LinearAlgebra() {}

LinearAlgebra::~LinearAlgebra() {}

bool LinearAlgebra::getCylinderFit(int n, double *x, double *y, double *z,
                                   Vector *p1, Vector *p2, double *radius) {
  // 1: Get the best fit Lxy on xy projection of points
  // 2: Get the best fit Lxz on xz projection of points
  // 3: Combine Lxy, Lxz to get the axis of cylinder
  // 4. Define the centroid to be a point on the axis
  // 5. Let radius of cyinder be average of radii of above two fits
  // 6. Some how get the extent of the axis
  double m1, m2, c1, c2, radius1, radius2;
  // step 1:
  if (!leastSquares(n, x, y, &m1, &c1, &radius1)) {
    return false;
  }
  // step 2:
  if (!leastSquares(n, x, z, &m2, &c2, &radius2)) {
    return false;
  }
  // step 3:
  // sin t = sqrt( m^2 / ( 1 + m^2 ) )
  // cos t = sqrt( 1 / ( 1 + m^2 ) )
  // The dcs of line 1 are norm( cos t1, sin t1, 0,              0 )
  // The dcs of line 2 are norm( cos t2, 0,              sin t2, 0 )
  // The DCs are norm( cost1+cos t2, sint1, sin t2, 0 )
  double sin_t1 = sqrt(m1 * m1 / (1.0 + m1 * m1));
  double cos_t1 = sqrt(1 / (1.0 + m1 * m1));
  if (m1 < 0) {
    sin_t1 = -1 * sin_t1;
  }
  double sin_t2 = sqrt(m2 * m2 / (1.0 + m2 * m2));
  double cos_t2 = sqrt(1 / (1.0 + m2 * m2));
  if (m2 < 0) {
    sin_t2 = -1 * sin_t2;
  }
  Vector DCs =
      Vector((float)(cos_t1 + cos_t2), (float)sin_t1, (float)sin_t2, 0);
  DCs.normalize();
  // step 4:
  double x0, y0, z0; // point on axis ( say centroid )
  if (!mean(x, n, &x0)) {
    return false;
  }
  if (!mean(y, n, &y0)) {
    return false;
  }
  if (!mean(z, n, &z0)) {
    return false;
  }
  // step 5:
  *radius = (radius1 + radius2) * 0.5;
  // step 6:
  // minLen = maxLen = 0
  // for each point P
  //             get norm ( p-p0 )
  //             compute cos t = dot prod ( axis dcs with prev result )
  //             len = len( ( p-p0 ) * cos t )
  //             update minLen, maxLen with len and sign of cos t
  //     end
  // Endpoints are center + max Len, center - minLen
  double minLen = 0, maxLen = 0;
  int i;
  for (i = 0; i < n; i++) {
    Vector diff((float)(x[i] - x0), (float)(y[i] - y0), (float)(z[i] - z0), 0);
    Vector normdiff = diff;
    normdiff.normalize();
    double cos_t = DCs.dot(normdiff);
    Vector proj(diff * (float)cos_t);
    double len = proj.norm();
    if (cos_t < 0) {
      len *= -1;
    }
    if (len < minLen) {
      minLen = len;
    }
    if (len > maxLen) {
      maxLen = len;
    }
  }
  // Endpoints are center + max Len, center - minLen
  p1->set((float)(x0 + minLen * DCs[0]), (float)(y0 + minLen * DCs[1]),
          (float)(z0 + minLen * DCs[2]), 1);
  p2->set((float)(x0 + maxLen * DCs[0]), (float)(y0 + maxLen * DCs[1]),
          (float)(z0 + maxLen * DCs[2]), 1);
  return true;
}

bool LinearAlgebra::leastSquares(int n, double *x, double *y, double *m,
                                 double *c, double *radius) {
  if (!x) {
    return false;
  }
  if (!y) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  double xMean = 0, yMean = 0, xSumSquare = 0, ySumSquare = 0, xDotY = 0;
  if (!mean(x, n, &xMean)) {
    return false;
  }
  if (!mean(y, n, &yMean)) {
    return false;
  }
  if (!sumOfSquares(x, n, &xSumSquare)) {
    return false;
  }
  if (!sumOfSquares(y, n, &ySumSquare)) {
    return false;
  }
  if (!dotProduct(x, y, n, &xDotY)) {
    return false;
  }
  double denomM = n * xMean * yMean - xDotY;
  // hope it wont happen ! if( (denomM < 0.000000001) && (denomM > -0.000000001
  // ) ) return false; // slope is 90 degrees, not really a error !
  double M = 0.5 *
             (ySumSquare - n * yMean * yMean - xSumSquare + n * xMean * xMean) /
             denomM;
  // quadratic has two solutions, get the least residual as the answer
  double m1 = -M + sqrt(M * M + 1);
  double m2 = -M - sqrt(M * M + 1);
  double c1 = yMean - m1 * xMean;
  double c2 = yMean - m2 * xMean;
  double residual1 = 0;
  double residual2 = 0;
  int i;
  double dist1 = 0, dist2 = 0;
  for (i = 0; i < n; i++) {
    dist1 += fabs(y[i] - c1 - m1 * x[i]) / sqrt(1 + m1 * m1);
    dist2 += fabs(y[i] - c2 - m2 * x[i]) / sqrt(1 + m2 * m2);
    residual1 +=
        (y[i] - c1 - m1 * x[i]) * (y[i] - c1 - m1 * x[i]) / (1 + m1 * m1);
    residual2 +=
        (y[i] - c2 - m2 * x[i]) * (y[i] - c2 - m2 * x[i]) / (1 + m2 * m2);
  }
  dist1 = dist1 / (double)n;
  dist2 = dist2 / (double)n;
  if (residual1 < residual2) {
    *m = m1;
    *c = c1;
    *radius = dist1;
  } else {
    *m = m2;
    *c = c2;
    *radius = dist2;
  }
  return true;
}

bool LinearAlgebra::mean(double *x, int n, double *mean) {
  if (!x) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  double sum = 0;
  if (!summation(x, n, &sum)) {
    return false;
  }
  *mean = sum / (double)n;
  return true;
}

bool LinearAlgebra::summation(double *x, int n, double *sum) {
  if (!x) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  int i;
  *sum = 0;
  for (i = 0; i < n; i++) {
    *sum += x[i];
  }
  return true;
}

bool LinearAlgebra::sumOfSquares(double *x, int n, double *sumSquare) {
  if (!x) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  int i;
  *sumSquare = 0;
  for (i = 0; i < n; i++) {
    *sumSquare += x[i] * x[i];
  }
  return true;
}

bool LinearAlgebra::dotProduct(double *x, double *y, int n, double *dotProd) {
  if (!x) {
    return false;
  }
  if (!y) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  if (!dotProd) {
    return false;
  }
  int i;
  *dotProd = 0;
  for (i = 0; i < n; i++) {
    *dotProd += x[i] * y[i];
  }
  return true;
}

bool LinearAlgebra::discretize(double *x, int n, double posVal, double negVal) {
  if (!x) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  int i;
  for (i = 0; i < n; i++) {
    if (x[i] < 0) {
      x[i] = negVal;
    }
    if (x[i] > 0) {
      x[i] = posVal;
    }
  }
  return true;
}

bool LinearAlgebra::correlate(double *x, double *y, int n,
                              double *correlationCoefficient) {
  if (!x) {
    return false;
  }
  if (!y) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  if (!correlationCoefficient) {
    return false;
  }
  double dotProd = 0;
  double norm1 = 0;
  double norm2 = 0;
  if (!dotProduct(x, y, n, &dotProd)) {
    return false;
  }
  if (!dotProduct(x, x, n, &norm1)) {
    return false;
  }
  if (!dotProduct(y, y, n, &norm2)) {
    return false;
  }
  norm1 = sqrt(norm1);
  norm2 = sqrt(norm2);
  double denom = norm1 * norm2;
  if (denom < 0.0000000000000001) {
    denom = 0.0000000000000001;
  }
  *correlationCoefficient = dotProd / denom;
  return true;
}

//  if usePositiveValues is true, take all positive x values and correlate with
//  corresponding y values.
bool LinearAlgebra::selectivelyCorrelate(double *x, int rangeToCorrelate,
                                         double *y, int n,
                                         double *correlationCoefficient) {
  if (!x) {
    return false;
  }
  if (!y) {
    return false;
  }
  if (n < 1) {
    return false;
  }
  if (!correlationCoefficient) {
    return false;
  }
  // find how many significant values are there in x.
  int nUseful = 0;
  if (rangeToCorrelate == 0) {
    nUseful = n;
  } else {
    int i;
    for (i = 0; i < n; i++) {
      if ((rangeToCorrelate > 0) && (x[i] > 0)) {
        nUseful++;
      } else if ((rangeToCorrelate < 0) && (x[i] < 0)) {
        nUseful++;
      }
    }
  }
  if (nUseful <= 0) {
    return false;
  }
  double *xUseful = new double[nUseful];
  double *yUseful = new double[nUseful];
  // fill up above arrays
  {
    int c = 0;
    int i;
    for (i = 0; i < n; i++) {
      if (rangeToCorrelate == 0) {
        xUseful[c] = x[i];
        yUseful[c] = y[i];
        c++;
      }
      if ((rangeToCorrelate > 0) && (x[i] > 0)) {
        xUseful[c] = x[i];
        yUseful[c] = y[i];
        c++;
      } else if ((rangeToCorrelate < 0) && (x[i] < 0)) {
        xUseful[c] = x[i];
        yUseful[c] = y[i];
        c++;
      }
    }
  }
  // correllate
  *correlationCoefficient = 0;
  if (!correlate(xUseful, yUseful, nUseful, correlationCoefficient)) {
    return false;
  }
  return true;
}

bool LinearAlgebra::solveSystem(double a11, double a12, double a13, double a21,
                                double a22, double a23, double a31, double a32,
                                double a33, double b1, double b2, double b3,
                                double *x, double *y, double *z) {
  double determinant = a11 * (a22 * a33 - a32 * a23) -
                       a12 * (a21 * a33 - a31 * a23) +
                       a13 * (a21 * a32 - a31 * a22);
  double xdeterminant = b1 * (a22 * a33 - a32 * a23) -
                        a12 * (b2 * a33 - b3 * a23) +
                        a13 * (b2 * a32 - b3 * a22);
  double ydeterminant = a11 * (b2 * a33 - b3 * a23) -
                        b1 * (a21 * a33 - a31 * a23) +
                        a13 * (a21 * b3 - a31 * b2);
  double zdeterminant = a11 * (a22 * b3 - a32 * b2) -
                        a12 * (a21 * b3 - a31 * b2) +
                        b1 * (a21 * a32 - a31 * a22);
  if (fabs(determinant) < 1e-10) {
    return false;
  }
  if (fabs(xdeterminant) < 1e-10) {
    return false;
  }
  if (fabs(ydeterminant) < 1e-10) {
    return false;
  }
  if (fabs(zdeterminant) < 1e-10) {
    return false;
  }
  *x = xdeterminant / determinant;
  *y = ydeterminant / determinant;
  *z = zdeterminant / determinant;
  return true;
}

bool LinearAlgebra::solve2x2System(double a11, double a12, double b1,
                                   double a21, double a22, double b2, double *x,
                                   double *y) {
  double d = a11 * a22 - a12 * a21;
  if (fabs(d) < 1e-10) {
    return false;
  }
  *x = (b1 * a22 - b2 * a12) / d;
  *y = (b2 * a11 - b1 * a21) / d;
  return true;
}

bool LinearAlgebra::solveDependentEquations(double c11, double c12, double c13,
                                            double c21, double c22, double c23,
                                            double c31, double c32, double c33,
                                            double *vec) {
  // try making x as 1, find y and z.
  {
    vec[0] = 1;
    if (solve2x2System(c12, c13, -c11, c22, c23, -c21, &(vec[1]), &(vec[2]))) {
      return true;
    }
    if (solve2x2System(c12, c13, -c11, c32, c33, -c31, &(vec[1]), &(vec[2]))) {
      return true;
    }
    if (solve2x2System(c22, c23, -c21, c32, c33, -c31, &(vec[1]), &(vec[2]))) {
      return true;
    }
  }
  // try making y as 1, find x and z.
  {
    vec[1] = 1;
    if (solve2x2System(c11, c13, -c12, c21, c23, -c22, &(vec[0]), &(vec[2]))) {
      return true;
    }
    if (solve2x2System(c11, c13, -c12, c31, c33, -c32, &(vec[0]), &(vec[2]))) {
      return true;
    }
    if (solve2x2System(c21, c23, -c22, c31, c33, -c32, &(vec[0]), &(vec[2]))) {
      return true;
    }
  }
  // try making z as 1, find x and y.
  {
    vec[2] = 1;
    if (solve2x2System(c11, c12, -c13, c21, c22, -c23, &(vec[0]), &(vec[1]))) {
      return true;
    }
    if (solve2x2System(c11, c12, -c13, c31, c32, -c33, &(vec[0]), &(vec[1]))) {
      return true;
    }
    if (solve2x2System(c21, c22, -c23, c31, c32, -c33, &(vec[0]), &(vec[1]))) {
      return true;
    }
  }
  return false;
}

bool LinearAlgebra::solveEigenSystem(double c11, double c12, double c13,
                                     double c21, double c22, double c23,
                                     double c31, double c32, double c33,
                                     double *k1Vec, double *k2Vec, double g1,
                                     double g2) {
  if (!solveDependentEquations(c11 - g1, c12, c13, c21, c22 - g1, c23, c31, c32,
                               c33 - g1, k1Vec)) {
    return false;
  }
  if (!solveDependentEquations(c11 - g2, c12, c13, c21, c22 - g2, c23, c31, c32,
                               c33 - g2, k2Vec)) {
    return false;
  }
  return true;
}

// Symmetric Householder reduction to tridiagonal form.
void LinearAlgebra::tridiagonalHouseholder(double mat[3][3], double val[3],
                                           double val2[3]) {
  //  This is derived from the Algol procedures tred2 by
  //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.

  for (int i = 0; i < 3; i++)
    val[i] = mat[2][i];

  // Householder reduction to tridiagonal form.
  for (int i = 2; i > 0; i--) {
    // Scale to avoid under/overflow.
    double scale = 0, h = 0;

    for (int j = 0; j < i; j++)
      scale += fabs(val[j]);

    if (scale == 0) {
      val2[i] = val[i - 1];

      for (int j = 0; j < i; j++) {
        val[j] = mat[i - 1][j];
        mat[i][j] = mat[j][i] = 0;
      }
    } else {
      // Generate Householder vector.
      for (int j = 0; j < i; j++) {
        val[j] = val[j] / scale;
        h += val[j] * val[j];
      }

      double f = val[i - 1], g = sqrt(h);

      if (f > 0)
        g = -g;

      val2[i] = scale * g;
      h -= f * g;
      val[i - 1] = f - g;

      for (int j = 0; j < i; j++)
        val2[j] = 0;

      // Apply similarity transformation to remaining columns.
      for (int j = 0; j < i; j++) {
        f = val[j];
        mat[j][i] = f;
        g = val2[j] + mat[j][j] * f;
        for (int k = j + 1; k <= i - 1; k++) {
          g += mat[k][j] * val[k];
          val2[k] += mat[k][j] * f;
        }
        val2[j] = g;
      }

      f = 0;
      for (int j = 0; j < i; j++) {
        val2[j] = val2[j] / h;
        f += val[j] * val2[j];
      }

      for (int j = 0; j < i; j++)
        val2[j] = val2[j] - 0.5 * f * val[j] / h;

      for (int j = 0; j < i; j++) {
        f = val[j];
        g = val2[j];

        for (int k = j; k <= i - 1; k++)
          mat[k][j] -= (f * val2[k] + g * val[k]);

        val[j] = mat[i - 1][j];
        mat[i][j] = 0;
      }
    }

    val[i] = h;
  }

  // Accumulate transformations.
  for (int i = 0; i < 2; i++) {
    mat[2][i] = mat[i][i];
    mat[i][i] = 1;

    double h = val[i + 1];

    if (h != 0) {
      for (int j = 0; j <= i; j++)
        val[j] = mat[j][i + 1] / h;

      for (int j = 0; j <= i; j++) {
        float g = 0;

        for (int k = 0; k <= i; k++)
          g += mat[k][i + 1] * mat[k][j];

        for (int k = 0; k <= i; k++)
          mat[k][j] -= g * val[k];
      }
    }

    for (int j = 0; j <= i; j++)
      mat[j][i + 1] = 0;
  }

  for (int j = 0; j < 3; j++) {
    val[j] = mat[2][j];
    mat[2][j] = 0;
  }

  mat[2][2] = 1;
  val2[0] = 0;
}

// Symmetric tridiagonal QL algorithm.
void LinearAlgebra::tridiagonalQL(double mat[3][3], double val[3],
                                  double val2[3]) {
  //  This is derived from the Algol procedures tql2, by
  //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.

  for (int i = 1; i < 3; i++)
    val2[i - 1] = val2[i];

  val2[2] = 0;

  double f = 0, tst1 = 0, eps = pow(2.0, -52.0);

  for (int l = 0; l < 3; l++) {
    // Find small subdiagonal element
    if (tst1 < fabs(val[l]) + fabs(val2[l]))
      tst1 = fabs(val[l]) + fabs(val2[l]);

    int m = l;

    while (m < 3) {
      if (fabs(val2[m]) <= eps * tst1)
        break;
      m++;
    }

    // If m == l, val[ l ] is an eigenvalue,
    // otherwise, iterate.
    if (m > l) {
      int iter = 0;

      do {
        iter++; // (Could check iteration count here.)

        // Compute implicit shift
        double g = val[l], p = (val[l + 1] - g) / (2.0 * val2[l]);
        double r = sqrt(p * p + 1);

        if (p < 0)
          r = -r;

        val[l] = val2[l] / (p + r);
        val[l + 1] = val2[l] * (p + r);

        double dl1 = val[l + 1], h = g - val[l];

        for (int i = l + 2; i < 3; i++)
          val[i] = val[i] - h;

        f = f + h;

        // Implicit QL transformation.
        p = val[m];
        double c = 1, c2 = 1, c3 = 1;
        double el1 = val2[l + 1];
        double s = 0, s2 = 0;

        for (int i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * val2[i];
          h = c * p;
          r = sqrt(p * p + val2[i] * val2[i]);
          val2[i + 1] = s * r;
          s = val2[i] / r;
          c = p / r;
          p = c * val[i] - s * g;
          val[i + 1] = h + s * (c * g + s * val[i]);

          // Accumulate transformation.
          for (int j = 0; j < 3; j++) {
            h = mat[j][i + 1];
            mat[j][i + 1] = s * mat[j][i] + c * h;
            mat[j][i] = c * mat[j][i] - s * h;
          }
        }

        p = -s * s2 * c3 * el1 * val2[l] / dl1;
        val2[l] = s * p;
        val[l] = c * p;

        // Check for convergence.
      } while (fabs(val2[l]) > eps * tst1);
    }

    val[l] = val[l] + f;
    if (val[l] < 0)
      val[l] = 0;
    val2[l] = 0;
  }

  // Sort eigenvalues and corresponding vectors.
  for (int i = 0; i < 2; i++) {
    int k = i;
    double p = val[i];

    for (int j = i + 1; j < 3; j++)
      if (val[j] > p) {
        k = j;
        p = val[j];
      }

    if (k != i) {
      val[k] = val[i];
      val[i] = p;

      for (int j = 0; j < 3; j++) {
        p = mat[j][i];
        mat[j][i] = mat[j][k];
        mat[j][k] = p;
      }
    }
  }
}

void LinearAlgebra::eigenSymmetric3x3(const Matrix &mat,
                                      std::vector<Vector> &vec,
                                      std::vector<float> &val) {
  double mat2[3][3];
  double val2[3], val3[3];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mat2[i][j] = mat.get(j, i);

  /*  for ( int i = 0; i < 3; i++ )
      {
       for ( int j = 0; j < 3; j++ )
         printf( "%9f ", mat2[ i ][ j ] );

       printf( "\n" );
      }  */

  tridiagonalHouseholder(mat2, val2, val3);

  tridiagonalQL(mat2, val2, val3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      vec[i][j] = mat2[j][i];

    vec[i][3] = 0;
    val[i] = val2[i];
  }
}
