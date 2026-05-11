#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <X11/StringDefs.h>
#include <Decimation/Control/cmds.h>
#include <Decimation/Bezier/BezierMath.h>
#include <Decimation/Misc/misc.h>

static void Draw1_8Sphere();
/*---------------------------------------------------------------------------
 * DrawLevelShpere -- Draw sphere with 8*4^level faces
----------------------------------------------------------------------------*/
void DrawLevelSphere(float center[], float r, int level)
/*
float   center[];
float   r;
int     level;
*/
{
  float p1[3], p2[3], p3[3], p4[3], p5[3], p6[3];
  int i, k;

  for (i = 0; i < 3; i++) {
    p1[i] = center[i];
    p2[i] = center[i];
    p3[i] = center[i];
    p4[i] = center[i];
    p5[i] = center[i];
    p6[i] = center[i];
  }
  p1[0] = p1[0] + r;
  p2[0] = p2[0] - r;
  p3[1] = p3[1] + r;
  p4[1] = p4[1] - r;
  p5[2] = p5[2] + r;
  p6[2] = p6[2] - r;

  k = level;
  Draw1_8Sphere(p1, p3, p5, center, r, &k);
  k = level;
  Draw1_8Sphere(p1, p3, p6, center, r, &k);
  k = level;
  Draw1_8Sphere(p2, p3, p5, center, r, &k);
  k = level;
  Draw1_8Sphere(p2, p3, p6, center, r, &k);
  k = level;
  Draw1_8Sphere(p1, p4, p5, center, r, &k);
  k = level;
  Draw1_8Sphere(p1, p4, p6, center, r, &k);
  k = level;
  Draw1_8Sphere(p2, p4, p5, center, r, &k);
  k = level;
  Draw1_8Sphere(p2, p4, p6, center, r, &k);
}
/* -------------------------------------------------------------------------*/
void Draw1_8Sphere(p1, p2, p3, center, r, k) float p1[], p2[], p3[], center[],
    r;
int *k;
{
  Vertex *vertices;
  Normal *normals;
  float p4[3], p5[3], p6[3], x, y, z;
  int i, kk;

  *k = *k - 1;
  for (i = 0; i < 3; i++) {
    p4[i] = p2[i] + p3[i] - 2 * center[i];
    p5[i] = p1[i] + p3[i] - 2 * center[i];
    p6[i] = p2[i] + p1[i] - 2 * center[i];
  }
  x = 0.0;
  y = 0.0;
  z = 0.0;
  for (i = 0; i < 3; i++) {
    x = p4[i] * p4[i] + x;
    y = p5[i] * p5[i] + y;
    z = p6[i] * p6[i] + z;
  }
  x = sqrt(x);
  y = sqrt(y);
  z = sqrt(z);
  for (i = 0; i < 3; i++) {
    p4[i] = r * p4[i] / x + center[i];
    p5[i] = r * p5[i] / y + center[i];
    p6[i] = r * p6[i] / z + center[i];
  }

  if (*k >= 0) {
    kk = *k;
    Draw1_8Sphere(p1, p5, p6, center, r, &kk);
    kk = *k;
    Draw1_8Sphere(p4, p5, p6, center, r, &kk);
    kk = *k;
    Draw1_8Sphere(p4, p5, p3, center, r, &kk);
    kk = *k;
    Draw1_8Sphere(p2, p4, p6, center, r, &kk);
  }
  if (*k < 0) {
    vertices = (Vertex *)malloc(3 * sizeof(Vertex));
    normals = (Normal *)malloc(3 * sizeof(Normal));
    for (i = 0; i < 3; i++) {
      vertices[0][i] = p1[i];
      vertices[1][i] = p2[i];
      vertices[2][i] = p3[i];
      normals[0][i] = (p1[i] - center[i]) / r;
      normals[1][i] = (p2[i] - center[i]) / r;
      normals[2][i] = (p3[i] - center[i]) / r;
    }
    DrawPolygon(3, vertices, normals);
  }

  /*
  if (*k <= 0 ) {
  vertices = (Vertex *) malloc (3 * sizeof (Vertex));
  normals = (Normal *) malloc (3 * sizeof (Normal));
  for (i = 0; i < 3; i++ ) {
      vertices [0][i] = p1[i];
      vertices [1][i] = p5[i];
      vertices [2][i] = p6[i];
      normals [0][i] = (p1[i] - center[i])/r;
      normals [1][i] = (p5[i] - center[i])/r;
      normals [2][i] = (p6[i] - center[i])/r;
  }
  DrawPolygon_Line(3, vertices, normals);
  vertices = (Vertex *) malloc (3 * sizeof (Vertex));
  normals = (Normal *) malloc (3 * sizeof (Normal));
  for (i = 0; i < 3; i++ ) {
      vertices [0][i] = p2[i];
      vertices [1][i] = p4[i];
      vertices [2][i] = p6[i];
      normals [0][i] = (p2[i] - center[i])/r;
      normals [1][i] = (p4[i] - center[i])/r;
      normals [2][i] = (p6[i] - center[i])/r;
  }
  DrawPolygon_Line(3, vertices, normals);
  vertices = (Vertex *) malloc (3 * sizeof (Vertex));
  normals = (Normal *) malloc (3 * sizeof (Normal));
  for (i = 0; i < 3; i++ ) {
      vertices [0][i] = p4[i];
      vertices [1][i] = p5[i];
      vertices [2][i] = p6[i];
      normals [0][i] = (p4[i] - center[i])/r;
      normals [1][i] = (p5[i] - center[i])/r;
      normals [2][i] = (p6[i] - center[i])/r;
  }
  DrawPolygon_Line(3, vertices, normals);
  vertices = (Vertex *) malloc (3 * sizeof (Vertex));
  normals = (Normal *) malloc (3 * sizeof (Normal));
  for (i = 0; i < 3; i++ ) {
      vertices [0][i] = p3[i];
      vertices [1][i] = p4[i];
      vertices [2][i] = p5[i];
      normals [0][i] = (p3[i] - center[i])/r;
      normals [1][i] = (p4[i] - center[i])/r;
      normals [2][i] = (p5[i] - center[i])/r;
  }
  DrawPolygon_Line(3, vertices, normals);
  }
  */
}

/*-----------------------------------------------------------------------------
 * DrawTube -- Draw a tube with end point p0 p1 and radius r
 *---------------------------------------------------------------------------*/
void DrawTube(p0, p1, r, k) float p0[], p1[], r;
int k;
{
  float c1, c2, s1, s2, p01[3], length, length1, theta, theta1, theta2;
  float q1[3], q2[3], nq1[3], nq2[3], x1, x2, y1, y2, nx1, nx2, ny1, ny2;
  int i, j;
  Vertex *vertices;
  Normal *normals;

  length = 0.0;
  for (i = 0; i < 3; i++) {
    p01[i] = p0[i] - p1[i];
    length = length + p01[i] * p01[i];
  }
  length = sqrt(length);
  length1 = sqrt(p01[0] * p01[0] + p01[1] * p01[1]);
  c1 = 1.0;
  s1 = 0.0;
  c2 = 1.0;
  s2 = 0.0;
  if (length1 > 0) {
    c1 = p01[2] / length;
    s1 = length1 / length;
    c2 = -p01[1] / length1;
    s2 = p01[0] / length1;
  }
  theta = 6.2831853 / k;
  for (j = 0; j < k; j++) {
    theta1 = j * theta;
    theta2 = theta1 + theta;
    nx1 = cos(theta1);
    ny1 = sin(theta1);
    nx2 = cos(theta2);
    ny2 = sin(theta2);

    x1 = nx1 * r;
    y1 = ny1 * r;
    x2 = nx2 * r;
    y2 = ny2 * r;
    q1[0] = c2 * x1 - c1 * s2 * y1;
    q1[1] = s2 * x1 + c1 * c2 * y1;
    q1[2] = s1 * y1;

    q2[0] = c2 * x2 - c1 * s2 * y2;
    q2[1] = s2 * x2 + c1 * c2 * y2;
    q2[2] = s1 * y2;

    nq1[0] = c2 * nx1 - c1 * s2 * ny1;
    nq1[1] = s2 * nx1 + c1 * c2 * ny1;
    nq1[2] = s1 * ny1;

    nq2[0] = c2 * nx2 - c1 * s2 * ny2;
    nq2[1] = s2 * nx2 + c1 * c2 * ny2;
    nq2[2] = s1 * ny2;

    vertices = (Vertex *)malloc(4 * sizeof(Vertex));
    normals = (Normal *)malloc(4 * sizeof(Normal));
    for (i = 0; i < 3; i++) {
      vertices[0][i] = q1[i] + p0[i];
      vertices[1][i] = q2[i] + p0[i];
      vertices[2][i] = q2[i] + p1[i];
      vertices[3][i] = q1[i] + p1[i];

      normals[0][i] = nq1[i];
      normals[1][i] = nq2[i];
      normals[2][i] = nq2[i];
      normals[3][i] = nq1[i];
    }
    DrawPolygon(4, vertices, normals);
  }
}

/*-----------------------------------------------------------------------------
 * DrawCone -- Draw a cone with end point p1, dirction p0-p1, radius r
 *             and hight h
 *---------------------------------------------------------------------------*/
void DrawCone(p0, p1, r, h, k) float p0[], p1[], r, h;
int k;
{
  float c1, c2, s1, s2, p01[3], length, length1, theta, theta1, theta2;
  float q1[3], q2[3], nq1[3], nq2[3], n[3], x1, x2, y1, y2, nx1, nx2, ny1, ny2,
      nz, nzl;
  int i, j;
  Vertex *vertices;
  Normal *normals;

  length = 0.0;
  for (i = 0; i < 3; i++) {
    p01[i] = p0[i] - p1[i];
    length = length + p01[i] * p01[i];
  }
  length = sqrt(length);
  length1 = sqrt(p01[0] * p01[0] + p01[1] * p01[1]);
  c1 = 1.0;
  s1 = 0.0;
  c2 = 1.0;
  s2 = 0.0;
  if (length1 > 0) {
    c1 = p01[2] / length;
    s1 = length1 / length;
    c2 = -p01[1] / length1;
    s2 = p01[0] / length1;
  }
  n[0] = s1 * s2;
  n[1] = -s1 * c2;
  n[2] = c1;

  nz = r / sqrt(r * r + h * h);
  nzl = sqrt(1.0 + r * r / (r * r + h * h));
  theta = 6.2831853 / k;
  for (j = 0; j < k; j++) {
    theta1 = j * theta;
    theta2 = theta1 + theta;
    nx1 = cos(theta1);
    ny1 = sin(theta1);
    nx2 = cos(theta2);
    ny2 = sin(theta2);

    x1 = nx1 * r;
    y1 = ny1 * r;
    x2 = nx2 * r;
    y2 = ny2 * r;
    q1[0] = c2 * x1 - c1 * s2 * y1;
    q1[1] = s2 * x1 + c1 * c2 * y1;
    q1[2] = s1 * y1;

    q2[0] = c2 * x2 - c1 * s2 * y2;
    q2[1] = s2 * x2 + c1 * c2 * y2;
    q2[2] = s1 * y2;

    nx1 = nx1 / nzl;
    ny1 = ny1 / nzl;
    nx2 = nx2 / nzl;
    ny2 = ny2 / nzl;
    nz = nz / nzl;

    nq1[0] = c2 * nx1 - c1 * s2 * ny1 + s1 * s2 * nz;
    nq1[1] = s2 * nx1 + c1 * c2 * ny1 - s1 * c2 * nz;
    nq1[2] = s1 * ny1 + c1 * nz;

    nq2[0] = c2 * nx2 - c1 * s2 * ny2 + s1 * s2 * nz;
    nq2[1] = s2 * nx2 + c1 * c2 * ny2 - s1 * c2 * nz;
    nq2[2] = s1 * ny2 + c1 * nz;

    vertices = (Vertex *)malloc(3 * sizeof(Vertex));
    normals = (Normal *)malloc(3 * sizeof(Normal));
    for (i = 0; i < 3; i++) {
      vertices[0][i] = q1[i] + p1[i];
      vertices[1][i] = q2[i] + p1[i];
      vertices[2][i] = ((h + length) * p1[i] - h * p0[i]) / length;

      normals[0][i] = nq1[i];
      normals[1][i] = nq2[i];
      normals[2][i] = 0.0;
    }
    DrawPolygon(3, vertices, normals);
    vertices = (Vertex *)malloc(3 * sizeof(Vertex));
    normals = (Normal *)malloc(3 * sizeof(Normal));
    for (i = 0; i < 3; i++) {
      vertices[0][i] = q1[i] + p1[i];
      vertices[1][i] = q2[i] + p1[i];
      vertices[2][i] = p1[i];

      normals[0][i] = n[i];
      normals[1][i] = n[i];
      normals[2][i] = n[i];
    }
    DrawPolygon(3, vertices, normals);
  }
}

/*-----------------------------------------------------------------------------
 * DrawCircle -- Draw a cone with center p, radius r
 *               and resolution n
 *---------------------------------------------------------------------------*/
void DrawCircle(p, r, n) float p[], r;
int n;
{
  int i;
  float x0, x1, y0, y1, ct, cti;

  x0 = p[0] + r;
  y0 = p[1];
  ct = 6.283185307 / n;
  for (i = 1; i < n + 1; i++) {
    cti = ct * i;
    x1 = p[0] + r * cos(cti);
    y1 = p[1] + r * sin(cti);

    LineAbs(x0, y0, 0.0, x1, y1, 0.0);
    x0 = x1;
    y0 = y1;
  }
}
