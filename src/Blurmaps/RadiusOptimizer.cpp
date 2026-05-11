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
#include <Blurmaps/RadiusOptimizer.h>

using namespace CVCUtility;

RadiusOptimizer::RadiusOptimizer() {
  numb_vert = 0;
  vert = 0;
  m_Initialized = false;
}

RadiusOptimizer::~RadiusOptimizer() { delete[] vert; }

bool RadiusOptimizer::init() {
  FILE *sphereFP = fopen("./Blurmaps/unitSphere.raw", "r");
  if (!sphereFP) {
    m_Initialized = false;
    return false;
  }
  int temp;
  if (fscanf(sphereFP, "%d %d\n", &numb_vert, &temp) != 2)
    error("Error reading from unitSphere.raw!");
  vert = new double[numb_vert * 3];
  // read sphere data
  int i;
  for (i = 0; i < numb_vert; i++) {
    double x, y, z;
    if (fscanf(sphereFP, "%lf %lf %lf\n", &x, &y, &z) != 3) {
      m_Initialized = false;
      fclose(sphereFP);
      return false;
    }
    vert[3 * i + 0] = x;
    vert[3 * i + 1] = y;
    vert[3 * i + 2] = z;
  }
  m_Initialized = true;
  fclose(sphereFP);
  return true;
}

bool RadiusOptimizer::getExtents(double *centers, double *radii,
                                 int numberOfAtoms, double *min, double *max,
                                 double *minRadius, double *maxRadius,
                                 double *avgRadius) {
  for (int i = 0; i < 3; i++) {
    min[i] = centers[3 * 0 + i] - radii[0];
    max[i] = centers[3 * 0 + i] + radii[0];
  }
  // read local atoms
  double radiusSum = 0.0;
  *minRadius = radii[0];
  *maxRadius = radii[0];
  for (int c = 0; c < numberOfAtoms; c++) {
    double r = radii[c];
    radiusSum += r;
    if (*minRadius > r) {
      *minRadius = r;
    }
    if (*maxRadius < r) {
      *maxRadius = r;
    }
    for (int i = 0; i < 3; i++) {
      if (centers[3 * c + i] - r < min[i]) {
        min[i] = centers[3 * c + i] - r;
      }
      if (centers[3 * c + i] + r > max[i]) {
        max[i] = centers[3 * c + i] + r;
      }
    }
  }
  *avgRadius = radiusSum / numberOfAtoms;
  for (int i = 0; i < 3; i++) {
    min[i] -= 0.5;
    max[i] += 0.5;
  }
  return true;
}

bool RadiusOptimizer::classifyVolume(double *centers, double *radii,
                                     int numberOfAtoms, double *minExt,
                                     double *maxExt, double *span,
                                     unsigned char *volume, int size) {
  int l;
  for (l = 0; l < numberOfAtoms; l++) {
    int imin[3];
    int imax[3];
    {
      int i;
      for (i = 0; i < 3; i++) {
        double w = centers[l * 3 + i] - minExt[i];
        int iw = (int)(w / span[i]);
        double r = radii[l];
        int ir = (int)(r / span[i] + 4);
        imin[i] = iw - ir;
        if (imin[i] < 0) {
          imin[i] = 0;
        }
        imax[i] = iw + ir;
        if (imax[i] >= size) {
          imax[i] = size - 1;
        }
      }
    }
    int i, j, k;
    for (i = imin[0]; i <= imax[0]; i++) {
      int isize = i * size;
      double x = minExt[0] + i * span[0] - centers[l * 3 + 0];
      for (j = imin[1]; j <= imax[1]; j++) {
        int jsize = size * (isize + j);
        double y = minExt[1] + j * span[1] - centers[l * 3 + 1];
        for (k = imin[2]; k <= imax[2]; k++) {
          int ksize = jsize + k;
          double z = minExt[2] + k * span[2] - centers[l * 3 + 2];
          double dist = x * x + y * y + z * z;
          if (dist < radii[l]) {
            volume[ksize] = 1;
          }
        }
      }
    }
  } // end l loop
  return true;
}

bool RadiusOptimizer::getOptimizedRadii(double *centers, double *radii,
                                        int numberOfAtoms) {
  if (!m_Initialized || !vert || !centers || !radii || numberOfAtoms < 1) {
    return false;
  }
  double minExt[3];
  double maxExt[3];
  double avgRadius;
  double minRadius;
  double maxRadius;
  if (!getExtents(centers, radii, numberOfAtoms, minExt, maxExt, &minRadius,
                  &maxRadius, &avgRadius)) {
    return false;
  }
  int size = 256;
  unsigned char *volume =
      (unsigned char *)malloc(size * size * size * sizeof(unsigned char));
  memset(volume, '0', size * size * size);
  double span[3];
  {
    int i;
    for (i = 0; i < 3; i++) {
      span[i] = (maxExt[i] - minExt[i]) / (size - 1.0);
    }
  }
  if (!classifyVolume(centers, radii, numberOfAtoms, minExt, maxExt, span,
                      volume, size)) {
    delete[] volume;
    return false;
  }
  double max_mini = 1.0;
  int adj = 0;
  double *adjust = (double *)malloc(numberOfAtoms * sizeof(double));
  {
    int i;
    for (i = 0; i < numberOfAtoms; i++) {
      adjust[i] = radii[i] * radii[i];
    }
  }
  double *adjust_local = (double *)malloc(numberOfAtoms * sizeof(double));
  while (fabs(max_mini) > 0.0001 && adj < 11) {
    printf("\tIteration number %d\n", adj);
    // for (adj = 0; adj < 11; adj++) {
    double maxi_dis = 0.0;
    double mean_dis = 0.0;
    double mean_integral = 0.0;
    double mean_dis_local = 0.0;
    int total = 0;
    double total_area = 0.0;
    int total_local = 0;
    int j = 0;
    double mini_max = 1000.0;
    double mean_max = 0.0;
    double max_mini = 0.0;
    int l;
    for (l = 0; l < numberOfAtoms; l++) {
      double cx = centers[l * 3 + 0];
      double cy = centers[l * 3 + 1];
      double cz = centers[l * 3 + 2];
      double r = radii[l];
      double mini_local = 1000.0;
      double local_max = 0.0;
      adjust_local[l] = adjust[l];
      int i;
      for (i = 0; i < numb_vert; i++) {
        int i3 = i * 3;
        double x = cx + r * vert[i3];
        double y = cy + r * vert[i3 + 1];
        double z = cz + r * vert[i3 + 2];
        // find a box
        int ixa = (int)((x - minExt[0]) / span[0]);
        int ixb = ixa + 1;
        int iya = (int)((y - minExt[1]) / span[1]);
        int iyb = iya + 1;
        int iza = (int)((z - minExt[2]) / span[2]);
        int izb = iza + 1;
        // printf("ixa, iya, iza = %d, %d, %d, %d, %d, %d, %d\n", ixa, iya, iza,
        // ixb, iyb, izb, (ixb*size+iyb)*size + izb);
        double sum = volume[(ixa * size + iya) * size + iza] +
                     volume[(ixa * size + iya) * size + izb] +
                     volume[(ixa * size + iyb) * size + iza] +
                     volume[(ixa * size + iyb) * size + izb] +
                     volume[(ixb * size + iya) * size + iza] +
                     volume[(ixb * size + iya) * size + izb] +
                     volume[(ixb * size + iyb) * size + iza] +
                     volume[(ixb * size + iyb) * size + izb];
        if (sum < 8) {
          // printf("sum = %d\n", sum);
          // solving one dimensional equations
          double error = 10000.0;
          double t = 0.0;
          int k = 0;
          j = j + 1;
          int classify = 1;
          while (error > 0.0001 && fabs(t) < 10.0 && k < 15 && classify == 1 &&
                 t > -r) {
            double p[3];
            p[0] = x + t * vert[i3];
            p[1] = y + t * vert[i3 + 1];
            p[2] = z + t * vert[i3 + 2];
            double C = -0.5;
            double func, grad[3];
            Func_Gradient_Sum_Gauss(centers, numberOfAtoms, p, C, &func, grad,
                                    &classify, l, adjust);
            if (classify == 1) {
              error = func / (vert[i3] * grad[0] + vert[i3 + 1] * grad[1] +
                              vert[i3 + 2] * grad[2]);
              t = t - error;
              error = fabs(error);
              k = k + 1;
            }
          }
          if (error <= 0.0001 && classify == 1 && t > -r) {
            if (mini_local > t) {
              mini_local = t;
            }
            if (maxi_dis < t) {
              maxi_dis = t;
            }
            if (local_max < fabs(t)) {
              local_max = fabs(t);
            }
            mean_dis = mean_dis + t;
            mean_integral = mean_integral + t * r * r;
            total = total + 1;
            total_area = total_area + r * r;
          }
        }
      } // end i loop
      if (fabs(mini_local) < 999.0) {
        // if (maxi_dis <  mini_local) maxi_dis = mini_local;
        if (fabs(max_mini) < fabs(mini_local) && adjust[l] > 0.0) {
          max_mini = fabs(mini_local);
        }
        mean_dis_local = mean_dis_local + fabs(mini_local);
        mean_max = mean_max + local_max;
        if (mini_max > local_max) {
          mini_max = local_max;
        }
        total_local = total_local + 1;
        adjust_local[l] =
            adjust[l] + r * r - (mini_local + r) * (mini_local + r);
      }
    } // end l loop
    {
      int c;
      for (c = 0; c < numberOfAtoms; c++) {
        double r = radii[c];
        adjust[c] = adjust_local[c];
        if (adjust[c] > r * r) {
          adjust[c] = r * r;
        }
        if (adjust[c] < 0.0) {
          adjust[c] = 0.0;
        }
      }
    }
    adj = adj + 1;
  } // end while
  {
    int c;
    for (c = 0; c < numberOfAtoms; c++) {
      radii[c] = sqrt(adjust[c]);
    }
  }
  delete[] adjust;
  delete[] adjust_local;
  delete[] volume;
  return true;
}

void RadiusOptimizer::Func_Gradient_Sum_Gauss(double *centers, int numb_atom,
                                              double *p, double C, double *func,
                                              double *grad, int *classify,
                                              int curAtom, double *radii) {
  int l;
  double w, e, x, y, z, result, this_var, CC, C2;
  result = 0.0;
  grad[0] = 0.0;
  grad[1] = 0.0;
  grad[2] = 0.0;
  *classify = 1;
  x = p[0] - centers[curAtom * 3 + 0];
  y = p[1] - centers[curAtom * 3 + 1];
  z = p[2] - centers[curAtom * 3 + 2];
  CC = C;
  C2 = CC + CC;
  // CC = C/atoms[l5+3];  // original
  w = radii[curAtom];
  w = CC * (x * x + y * y + z * z - w);
  // printf("adjust = %f\n", adjust[mo]);
  this_var = exp(w);
  for (l = 0; l < numb_atom; l++) {
    x = p[0] - centers[curAtom * 3 + 0];
    y = p[1] - centers[curAtom * 3 + 1];
    z = p[2] - centers[curAtom * 3 + 2];
    w = CC * (x * x + y * y + z * z - w);
    if (w > -14.0) {
      e = exp(w);
      if (e > this_var) {
        *classify = 0;
        return;
      }
      result = result + e;
      e = C2 * e;
      grad[0] = grad[0] + e * x;
      grad[1] = grad[1] + e * y;
      grad[2] = grad[2] + e * z;
    }
  }
  *func = result - 1.0;
  //*func = result - 0.5;
}

void RadiusOptimizer::Curvature_Sum_Gauss(double *atoms, int numb_atom,
                                          double *p, double C, int *classify,
                                          int mo, double *maxk) {
  int i, l, l5;
  double w, e, x, y, z, this_var, CC, C2, g, w1, w2;
  double fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz;
  double t1[3], t2[3], t11[3], t12[3], t22[3];
  double phx, phy, phxx, phxy, phyy;
  double g11, g12, g22, g111, g112, g122, g211, g212, g222;
  double first, second, g1122, g1212, K, H3[3], H2, H, k1, k2;
  fx = 0.0;
  fy = 0.0;
  fz = 0.0;
  fxx = 0.0;
  fyy = 0.0;
  fzz = 0.0;
  fxy = 0.0;
  fxz = 0.0;
  fyz = 0.0;
  *classify = 1;
  l5 = mo * 6;
  x = p[0] - atoms[l5];
  y = p[1] - atoms[l5 + 1];
  z = p[2] - atoms[l5 + 2];
  CC = C;
  C2 = CC + CC;
  // CC = C/atoms[l5+3];  // original
  w = CC * (x * x + y * y + z * z - atoms[l5 + 3]);
  this_var = exp(w);
  for (l = 0; l < numb_atom; l++) {
    l5 = l * 6;
    x = p[0] - atoms[l5];
    y = p[1] - atoms[l5 + 1];
    z = p[2] - atoms[l5 + 2];
    w = CC * (x * x + y * y + z * z - atoms[l5 + 3]);
    if (w > -14.0) {
      e = exp(w);
      if (e > this_var) {
        *classify = 0;
        return;
      }
      g = C2 * e;
      fx = fx + g * x;
      fy = fy + g * y;
      fz = fz + g * z;
      fxx = fxx + g * (x * x * C2 + 1.0);
      fyy = fyy + g * (y * y * C2 + 1.0);
      fzz = fzz + g * (z * z * C2 + 1.0);
      fxy = fxy + g * y * C2 * x;
      fxz = fxz + g * z * C2 * x;
      fyz = fyz + g * z * C2 * y;
    }
  } // end l loop
  // printf("fx, fy, fz, fxx, fxy, fyz = %f, %f, %f, %f, %f, %f,%f,%f,%f\n", fx,
  // fy, fz, fxx, fyy,fzz, fxy, fxz,fyz);
  if (fabs(fz) >= fabs(fx) && fabs(fz) >= fabs(fy)) {
    phx = -fx / fz;
    phy = -fy / fz;
    phxx = -(fxx + 2.0 * fxz * phx + fzz * phx * phx) / fz;
    phyy = -(fyy + 2.0 * fyz * phy + fzz * phy * phy) / fz;
    phxy = -(fxy + fxz * phy + fyz * phx + fzz * phx * phy) / fz;
    t1[0] = 1.0;
    t1[1] = 0.0;
    t1[2] = phx;
    t2[0] = 0.0;
    t2[1] = 1.0;
    t2[2] = phy;
    t11[0] = 0.0;
    t11[1] = 0.0;
    t11[2] = phxx;
    t22[0] = 0.0;
    t22[1] = 0.0;
    t22[2] = phyy;
    t12[0] = 0.0;
    t12[1] = 0.0;
    t12[2] = phxy;
  }
  if (fabs(fy) >= fabs(fx) && fabs(fy) >= fabs(fz)) {
    phx = -fx / fy;
    phy = -fz / fy;
    phxx = -(fxx + 2.0 * fxy * phx + fyy * phx * phx) / fy;
    phyy = -(fzz + 2.0 * fyz * phy + fyy * phy * phy) / fy;
    phxy = -(fxz + fxy * phy + fyz * phx + fyy * phx * phy) / fy;
    t1[0] = 1.0;
    t1[1] = phx;
    t1[2] = 0.0;
    t2[0] = 0.0;
    t2[1] = phy;
    t2[2] = 1.0;
    t11[0] = 0.0;
    t11[1] = phxx;
    t11[2] = 0.0;
    t22[0] = 0.0;
    t22[1] = phyy;
    t22[2] = 0.0;
    t12[0] = 0.0;
    t12[1] = phxy;
    t12[2] = 0.0;
  }
  if (fabs(fx) >= fabs(fy) && fabs(fx) >= fabs(fz)) {
    // printf("xxxxxxxxxx\n");
    phx = -fz / fx;
    phy = -fy / fx;
    phxx = -(fzz + 2.0 * fxz * phx + fxx * phx * phx) / fx;
    phyy = -(fyy + 2.0 * fxy * phy + fxx * phy * phy) / fx;
    phxy = -(fyz + fxz * phy + fxy * phx + fxx * phx * phy) / fx;
    t1[0] = phx;
    t1[1] = 0.0;
    t1[2] = 1.0;
    t2[0] = phy;
    t2[1] = 1.0;
    t2[2] = 0.0;
    t11[0] = phxx;
    t11[1] = 0.0;
    t11[2] = 0.0;
    t22[0] = phyy;
    t22[1] = 0.0;
    t22[2] = 0.0;
    t12[0] = phxy;
    t12[1] = 0.0;
    t12[2] = 0.0;
  }
  g11 = DotProduct(t1, t1);
  g12 = DotProduct(t1, t2);
  g22 = DotProduct(t2, t2);
  g = g11 * g22 - g12 * g12;
  // printf("g11, g12, g22, g = %f, %f, %f, %f\n", g11, g12, g22, g);
  g111 = DotProduct(t1, t11);
  g112 = DotProduct(t1, t12);
  g122 = DotProduct(t1, t22);
  g211 = DotProduct(t2, t11);
  g212 = DotProduct(t2, t12);
  g222 = DotProduct(t2, t22);
  // mean curvature
  first = g22 * g111 + g11 * g122 - 2 * g12 * g112;
  second = g22 * g211 + g11 * g222 - 2 * g12 * g212;
  w1 = (first * g22 - second * g12) / g;
  w2 = (-first * g12 + second * g11) / g;
  for (i = 0; i < 3; i++) {
    H3[i] = (g22 * t11[i] + g11 * t22[i] - 2 * g12 * t12[i] - w1 * t1[i] -
             w2 * t2[i]) /
            (g + g);
  }
  if (H3[0] * fx + H3[1] * fy + H3[2] * fz < 0.0) {
    H2 = DotProduct(H3, H3);
    H = sqrt(H2);
    g1122 = DotProduct(t11, t22);
    g1212 = DotProduct(t12, t12);
    first = g1122 - g1212;
    second =
        (g111 * g22 - g211 * g12) * g122 + (-g111 * g12 + g211 * g11) * g222 -
        (g112 * g22 - g212 * g12) * g112 - (-g112 * g12 + g212 * g11) * g212;
    K = (first - second / g) / g;
    // printf("Gauss = %f, Mean = %f, direction = %f, p = %f, %f, %f\n",K, H,
    // H3[0]*fx+H3[1]*fy+H3[2]*fz, p[0],p[1],p[2]);
    //  compute principal curvature
    w = H2 - K;
    if (w < 0.0) {
      printf("w = %f\n", w);
      w = 0.0;
    }
    w = sqrt(w);
    k1 = fabs(H + w);
    k2 = fabs(H - w);
    // printf("k1 = %f,k2 = %f, %f, %f\n", k1, k2, k1*k2, 0.5*(k1+k2));
    *maxk = k1;
    if (*maxk < k2) {
      *maxk = k2;
    }
    // printf("k1 = %f,k2 = %f, %f, %f, maxk = %f\n", k1, k2, k1*k2,
    // 0.5*(k1+k2), *maxk);
  }
}

void RadiusOptimizer::Spherical_Fit(double *x, double *y, double *z, int npts,
                                    double *radius) {
  int i, j, k, info;
  double r;
  double matrixlarge[650], matrix4[16], leftlarge[162], left4[4], solution[4],
      sum;
  for (i = 0; i < npts; i++) {
    matrixlarge[4 * i] = x[i] + x[i];
    matrixlarge[4 * i + 1] = y[i] + y[i];
    matrixlarge[4 * i + 2] = z[i] + z[i];
    matrixlarge[4 * i + 3] = 1.0;
    leftlarge[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
    // printf("leftlarge[i] = %f\n", leftlarge[i]);
  }
  for (i = 0; i < 4; i++) {
    sum = 0.0;
    for (k = 0; k < npts; k++) {
      sum = sum + matrixlarge[4 * k + i] * leftlarge[k];
    }
    left4[i] = sum;
    for (j = 0; j < 4; j++) {
      sum = 0.0;
      for (k = 0; k < npts; k++) {
        sum = sum + matrixlarge[4 * k + i] * matrixlarge[4 * k + j];
      }
      matrix4[i * 4 + j] = sum;
    }
  }
  gaussinverse(matrix4, 4, 0.000L, &info);
  if (info == 0) {
    printf("The matrix  is singular!\n");
    return;
  }
  MatrixTimesMatrix(matrix4, left4, solution, 4, 4, 1);
  r = solution[0] * solution[0] + solution[1] * solution[1] +
      solution[2] * solution[2] + solution[3];
  r = sqrt(r);
  *radius = r;
}

/*************************************************************************/
/* GaussInverse                                                          */
/* - Inverse the matrix A by Gauss elimination methods                   */
/*   matrix is instored in the order a[0][0], a[0][1], ....              */
/*************************************************************************/
void RadiusOptimizer::gaussinverse(double *a, int n, double eps, int *message)
// int             n, *message;
// double     *a, eps;
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
    max = 0.0L;
    for (i = k; i < n; i++)
      for (j = k; j < n; j++)
        if (fabs(*(a + i * n + j)) > max) {
          ik = i;
          jk = j;
          max = fabs(*(a + i * n + j));
        };
    if (max < eps || max == 0.0L) {
      *message = 0;
      printf("The matrix in gaussinverse is singular, %e\n", max);
      return;
    };
    max = 1.0 / *(a + ik * n + jk);
    *(a + ik * n + jk) = 1.0L;
    z[2 * k] = ik;
    z[2 * k + 1] = jk;
    Exchangerowcolumn(a, n, k, ik, jk);
    for (j = 0; j < n; j++) {
      *(a + k * n + j) = *(a + k * n + j) * max;
    }
    for (i = 0; i < n; i++)
      if (i != k) {
        max = *(a + i * n + k);
        *(a + i * n + k) = 0.0L;
        for (j = 0; j < n; j++) {
          *(a + i * n + j) = *(a + i * n + j) - max * *(a + k * n + j);
        }
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
/* Exchangerowcolumn                                                     */
/* - Exchange two rows and two columns of the matrix A                   */
/*************************************************************************/
void RadiusOptimizer::Exchangerowcolumn(double *a, int n, int k, int ik, int jk)
// int             n,k,ik,jk;
// double          *a;
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

/* MatrixTimesMatrix                                                     */
/* times a matrix a with matrix v, the result is x. i.e., x = a*v        */
void RadiusOptimizer::MatrixTimesMatrix(double *a, double *v, double *x, int m,
                                        int n, int l)
// int         m,n,l;
// double      *a;            /* pointe to a (m,n) matrix               */
// double      *v;            /* pointe to a (n,l) matrix               */
// double      *x;            /* pointe to a (m,l) matrix               */
{
  int i, j, k;
  double sum;
  for (k = 0; k < l; k++)
    for (i = 0; i < m; i++) {
      sum = 0.0L;
      for (j = 0; j < n; j++) {
        sum = sum + *(a + i * n + j) * *(v + j * l + k);
      }
      *(x + i * l + k) = sum;
    }
}

// DotProduct -- compute dot product
double RadiusOptimizer::DotProduct(double *p1, double *p2)
// double   p1[],p2[];
{
  double result;
  result = p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
  /* normalize                     */
  /*
  result = result/sqrt((p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]));
  result = result/sqrt((p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]));
  */
  return (result);
}
