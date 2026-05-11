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
 * root.c  --  Find roots for polynomials
 *      Xu, guoliang     April, 1993, Modified Nov. 1996
 */

#include <stdio.h>
#include <math.h>

/*---------------------------------------------------------------------------
FindRoots1 -- compute degree 1 polynomial roots
----------------------------------------------------------------------------*/
void FindRoots1(a, b, root, number) double a, b, root[];
int *number;
{
  if (fabs(a) < 0.000000001) {
    *number = 0;
    return;
  }
  root[0] = -b / a;
  *number = 1;
}
/*---------------------------------------------------------------------------
FindRoots2 -- compute degree 2 polynomial roots of the equation
              a x^2 + b x + c = 0
----------------------------------------------------------------------------*/
void FindRoots2(a, b, c, root, number) double a, b, c, root[];
int *number;
{
  double det;

  if (fabs(a) < 0.000000001) {
    FindRoots1(b, c, root, number);
    return;
  }

  det = b * b - 4 * a * c;
  if (fabs(det) < 0.0000001)
    det = 0.0;
  if (det < 0.0) {
    *number = 0;
    return;
  }
  *number = 2;
  det = sqrt(det);
  root[0] = (-b + det) / (a + a);
  root[1] = (-b - det) / (a + a);
}

/*---------------------------------------------------------------------------
FindRoots2More -- compute degree 2 polynomial roots of the equation
              a x^2 + b x + c = 0, try to find more roots
----------------------------------------------------------------------------*/
void FindRoots2More(a, b, c, root, number) double a, b, c, root[];
int *number;
{
  double det;

  if (fabs(a) < 0.000000001) {
    FindRoots1(b, c, root, number);
    return;
  }

  det = b * b - 4 * a * c;
  printf("det = %e \n", det);
  /*
  if (fabs(det) < 0.0001)  det = 0.0;
  if (det < 0.0) {
     *number = 0;
     return;
  }
  */

  /* changed for bob  */
  if (det <= -0.0001) {
    *number = 0;
    return;
  }
  if (det < 0.0 && det > -0.0001) {
    det = 0.0;
  }

  *number = 2;
  det = sqrt(det);
  root[0] = (-b + det) / (a + a);
  root[1] = (-b - det) / (a + a);
}

/*---------------------------------------------------------------------------
FindRoots3 -- compute degree 3 polynomial ax^3 + bx^2 + cx + d = 0 roots
----------------------------------------------------------------------------*/
void FindRoots3(a, b, c, d, root, number) double a, b, c, d, root[];
int *number;
{
  double ba3, a2, b2, a227, b3, ac, abc, p, q;
  double det, q2, p3, r, theta, onethird, p33, sdet, result;
  double abs1, abs2, sign1, sign2;

  if (fabs(a) < 0.000000001) {
    FindRoots2(b, c, d, root, number);
    return;
  }
  *number = 0;
  ba3 = b / (3 * a);
  a2 = a * a;
  b2 = b * b;
  a227 = a2 * 27;
  b3 = b2 * b;
  ac = a * c;
  abc = ac * b;
  p = (3 * ac - b2) / (3 * a2);
  q = (b3 + b3 - 9 * abc + d * a227) / (a * a227);

  onethird = 0.333333333333;
  q2 = 0.5 * q;
  p3 = p / 3.0;
  p33 = p3 * p3 * p3;
  det = q2 * q2 + p33;
  if (det >= 0.0) {
    sdet = sqrt(det);
    abs1 = -q2 + sdet;
    abs2 = -q2 - sdet;
    if (abs1 >= 0)
      sign1 = 1.0;
    else
      sign1 = -1.0;
    if (abs2 >= 0)
      sign2 = 1.0;
    else
      sign2 = -1.0;
    root[*number] = sign1 * pow(fabs(abs1), onethird) +
                    sign2 * pow(fabs(abs2), onethird) - ba3;
    ;
    *number = *number + 1;
    if (det == 0.0) {
      if (-q2 >= 0)
        sign1 = 1.0;
      else
        sign1 = -1.0;
      root[*number] = -sign1 * pow(fabs(-q2), onethird) - ba3;
      ;
      root[*number + 1] = root[*number];
      *number = *number + 2;
    }
    return;
  }
  r = sqrt(-p33);
  theta = acos(-q2 / r) * onethird;
  q2 = 2 * sqrt(-p3);
  root[*number] = q2 * cos(theta) - ba3;
  ;
  *number = *number + 1;
  root[*number] = q2 * cos(theta + 2.094395102) - ba3;
  ;
  *number = *number + 1;
  root[*number] = q2 * cos(theta + 4.188790205) - ba3;
  ;
  *number = *number + 1;
}

/*---------------------------------------------------------------------------
FindRoots4 -- compute all the real roots of a degree 4 polynomial equation
              ax^4 + bx^3 + cx^2 + dx + e = 0
----------------------------------------------------------------------------*/
void FindRoots4(a, b, c, d, e, root, number) double a, b, c, d, e, root[];
int *number;
{
  double aa, bb, cc, dd, ee, y, ybc, b4c, byd, ra, roots[4], max;
  int flag;

  max = fabs(a);
  if (fabs(b) > max)
    max = fabs(b);
  if (fabs(c) > max)
    max = fabs(c);
  if (fabs(d) > max)
    max = fabs(d);
  if (fabs(e) > max)
    max = fabs(e);

  aa = a / max;
  bb = b / max;
  cc = c / max;
  dd = d / max;
  ee = e / max;

  if (fabs(aa) < 0.000001) {
    FindRoots3(bb, cc, dd, ee, root, number);
    return;
  }
  bb = bb / aa;
  cc = cc / aa;
  dd = dd / aa;
  ee = ee / aa;

  if (fabs(ee) < 0.000001) {
    FindRoots3(1.0, bb, cc, dd, root, number);
    root[*number] = 0.0;
    *number = *number + 1;
    return;
  }

  b4c = bb * bb - 4 * cc;
  FindRoots3(8.0, -4 * cc, 2 * bb * dd - 8 * ee, -ee * b4c - dd * dd, root,
             number);
  y = root[0];

  ybc = 8 * y + b4c;
  if (ybc < 0.0)
    ybc = -ybc;
  ybc = sqrt(ybc);
  byd = bb * y - dd;

  if (ybc <= 0.0) {
    ra = 0.0;
  } else {
    ra = byd / ybc;
  }

  *number = 0;
  FindRoots2More(1.0, 0.5 * (bb + ybc), y + ra, root, &flag);
  if (flag != 0) {
    *number = 2;
  }
  FindRoots2More(1.0, 0.5 * (bb - ybc), y - ra, roots, &flag);
  if (flag != 0) {
    root[*number] = roots[0];
    root[*number + 1] = roots[1];
    *number = *number + 2;
  }
}

/*---------------------------------------------------------------------------
BernsteinRoot -- find roots of  Bernstein polynimial of degree 1, or 2 or
                 3 or 4  in [0,1].
----------------------------------------------------------------------------*/
void BernsteinRoot(weight, roots, number, degree) double weight[];
double roots[];
int *number, degree;
{
  double a, b, c, d, e, x, r1, r2, root[2];
  int i, k, flag;

  *number = 0;
  if (degree == 1) {
    if (weight[0] == 0.0 && weight[1] == 0.0) {
      *number = 1000;
      return;
    }
    if (weight[0] == weight[1]) {
      return;
    }
    x = weight[0] / (weight[0] - weight[1]);
    if (x >= 0.0 && x <= 1.0) {
      roots[*number] = x;
      *number = *number + 1;
    }
    return;
  }
  if (degree == 2) {
    if (weight[0] == 0.0 && weight[1] == 0.0 && weight[2] == 0.0) {
      *number = 1000;
      return;
    }
    if (weight[0] == 0.0 && weight[2] == 0.0) {
      roots[0] = 0.0;
      roots[1] = 1.0;
      *number = 2;
      return;
    }

    if (weight[0] == 0.0 && weight[1] == 0.0) {
      roots[0] = 0.0;
      *number = 1;
      return;
    }
    if (weight[0] == 0.0) {
      roots[0] = 0.0;
      *number = 1;
      x = 2 * weight[1] / (2 * weight[1] - weight[2]);
      if (x >= 0.0 && x <= 1.0) {
        roots[*number] = x;
        *number = *number + 1;
      }
      return;
    }

    if (weight[2] == 0.0 && weight[1] == 0.0) {
      roots[0] = 1.0;
      *number = 1;
      return;
    }

    if (weight[2] == 0.0) {
      roots[*number] = 1.0;
      *number = *number + 1;
      x = weight[0] / (weight[0] - 2 * weight[1]);
      if (x >= 0.0 && x <= 1.0) {
        roots[*number] = x;
        *number = *number + 1;
      }
      return;
    }
    a = weight[2];
    b = 2 * weight[1];
    c = weight[0];
    FindRoots2(a, b, c, root, &flag);
    r1 = root[0];
    r2 = root[1];
    if (flag == 0)
      return;
    r1 = r1 / (1.0 + r1);
    r2 = r2 / (1.0 + r2);
    if (r1 >= 0.0 && r1 <= 1.0) {
      roots[*number] = r1;
      *number = *number + 1;
    }
    if (r2 >= 0.0 && r2 <= 1.0) {
      roots[*number] = r2;
      *number = *number + 1;
    }
    return;
  }
  if (degree == 3) {
    if (weight[0] == 0.0 && weight[1] == 0.0 && weight[2] == 0.0 &&
        weight[3] == 0.0) {
      *number = 1000;
      return;
    }

    if (weight[0] == 0.0 && weight[3] == 0.0) {
      roots[*number] = 0.0;
      *number = *number + 1;
      roots[*number] = 1.0;
      *number = *number + 1;
      x = weight[1] / (weight[1] - weight[2]);
      if (x >= 0.0 && x <= 1.0) {
        roots[*number] = x;
        *number = *number + 1;
      }
      return;
    }

    if (weight[0] == 0.0 && weight[1] == 0.0 && weight[2] == 0.0) {
      roots[0] = 0.0;
      *number = 1;
      return;
    }

    if (weight[0] == 0.0 && weight[1] == 0.0) {
      roots[0] = 0.0;
      *number = 1;
      x = 3 * weight[2] / (3 * weight[2] - weight[3]);
      if (x >= 0.0 && x <= 1.0) {
        roots[*number] = x;
        *number = *number + 1;
      }
      return;
    }

    if (weight[0] == 0.0) {
      roots[*number] = 0.0;
      *number = *number + 1;
      a = weight[3];
      b = 3 * weight[2];
      c = 3 * weight[1];
      FindRoots2(a, b, c, root, &flag);
      r1 = root[0];
      r2 = root[1];
      if (flag == 0)
        return;
      r1 = r1 / (1.0 + r1);
      r2 = r2 / (1.0 + r2);
      if (r1 >= 0.0 && r1 <= 1.0) {
        roots[*number] = r1;
        *number = *number + 1;
      }
      if (r2 >= 0.0 && r2 <= 1.0) {
        roots[*number] = r2;
        *number = *number + 1;
      }
      return;
    }

    if (weight[3] == 0.0 && weight[1] == 0.0 && weight[2] == 0.0) {
      roots[0] = 1.0;
      *number = 1;
      return;
    }

    if (weight[3] == 0.0 && weight[2] == 0.0) {
      roots[0] = 1.0;
      *number = 1;
      x = weight[0] / (weight[0] - 3 * weight[1]);
      if (x >= 0.0 && x <= 1.0) {
        roots[*number] = x;
        *number = *number + 1;
      }
      return;
    }

    if (weight[3] == 0.0) {
      roots[*number] = 1.0;
      *number = *number + 1;
      a = 3 * weight[2];
      b = 3 * weight[1];
      c = weight[0];
      FindRoots2(a, b, c, root, &flag);
      r1 = root[0];
      r2 = root[1];
      if (flag == 0)
        return;
      r1 = r1 / (1.0 + r1);
      r2 = r2 / (1.0 + r2);
      if (r1 >= 0.0 && r1 <= 1.0) {
        roots[*number] = r1;
        *number = *number + 1;
      }
      if (r2 >= 0.0 && r2 <= 1.0) {
        roots[*number] = r2;
        *number = *number + 1;
      }
      return;
    }

    a = weight[3];
    b = 3 * weight[2];
    c = 3 * weight[1];
    d = weight[0];
    FindRoots3(a, b, c, d, roots, &k);
    for (i = 0; i < k; i++) {
      roots[i] = roots[i] / (1.0 + roots[i]);
      if (roots[i] >= 0.0 && roots[i] <= 1.0) {
        roots[*number] = roots[i];
        *number = *number + 1;
      }
    }
    return;
  }

  if (degree == 4) {
    /* here need to check if the coefficients is zero   */
    /* especially a should not be zero                  */
    a = weight[4];
    b = 4 * weight[3];
    c = 6 * weight[2];
    d = 4 * weight[1];
    e = weight[0];
    FindRoots4(a, b, c, d, e, roots, &k);
    for (i = 0; i < k; i++) {
      roots[i] = roots[i] / (1.0 + roots[i]);
      if (roots[i] >= 0.0 && roots[i] <= 1.0) {
        roots[*number] = roots[i];
        *number = *number + 1;
      }
    }
    return;
  }
}
/*---------------------------------------------------------------------------
BernsteinRootHigh -- find roots of  Bernstein polynimial of any degree
                     in [0,1].
----------------------------------------------------------------------------*/
void BernsteinRootHigh(weight, roots, number, degree) double weight[];
double roots[];
int *number, degree;
{
  int i;
  struct complex {
    double real, imag;
  } complex_root[10];
  double cni, weight1[6], work[24];
  int iwork[6];

  cni = 1.0;
  for (i = 0; i <= degree; i++) {
    weight1[degree - i] = cni * weight[i];
    cni = cni * (degree - i) / (i + 1.0);
  }

  /*
  dzplrc (&degree, weight1, complex_root);

  *number = 0;
  for (i = 0; i < degree; i++ ) {
      if (complex_root[i].real >= 0.0 &&
          fabs(complex_root[i].imag) < 0.0001) {
          roots[*number] = complex_root[i].real;
          *number = *number + 1;
      }
  }
  */
  *number = sturm(weight1, degree, 0.0, 1000.0, 0.0001, work, iwork);
  for (i = 0; i < *number; i++) {
    /*
    roots[i] = roots[i]/(1.0 + roots[i]);
    */
    roots[i] = weight1[i] / (1.0 + weight1[i]);
  }
}
