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
#include <SphericalPatchIntersections/SphericalPatch.h>

SphericalPatch::SphericalPatch() {}

SphericalPatch::~SphericalPatch() {}

double SphericalPatch::euclideanDist3D(double p1[3], double p2[3]) {
  return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
              (p1[1] - p2[1]) * (p1[1] - p2[1]) +
              (p1[2] - p2[2]) * (p1[2] - p2[2]));
}

// We should ensure that center != point before coming in to this function
double SphericalPatch::intersectCircle(double center[2], double radius,
                                       double point[2],
                                       double intersection[2]) {
  double temp0 = center[0] - point[0];
  double temp1 = center[1] - point[1];
  double distToCenter = sqrt((temp0) * (temp0) + (temp1) * (temp1));
  intersection[0] = (-temp0) / distToCenter * radius + center[0];
  intersection[1] = (-temp1) / distToCenter * radius + center[1];
  return fabs(radius - distToCenter);
}

bool SphericalPatch::intersectionInCube(double intersection[3], double bmin[3],
                                        double bmax[3]) {
  if (intersection[0] < bmin[0]) {
    return false;
  }
  if (intersection[0] > bmax[0]) {
    return false;
  }
  if (intersection[1] < bmin[1]) {
    return false;
  }
  if (intersection[1] > bmax[1]) {
    return false;
  }
  if (intersection[2] < bmin[2]) {
    return false;
  }
  if (intersection[2] > bmax[2]) {
    return false;
  }
  return true;
}

bool SphericalPatch::intersectionInCube2D(double intersection[2], double minx,
                                          double maxx, double miny,
                                          double maxy) {
  if (intersection[0] < minx) {
    return false;
  }
  if (intersection[0] > maxx) {
    return false;
  }
  if (intersection[1] < miny) {
    return false;
  }
  if (intersection[1] > maxy) {
    return false;
  }
  return true;
}

double SphericalPatch::pointToCircle(int dim0, int dim1, int dim2,
                                     double sphereToCircle, double r,
                                     double c[3], double point[3],
                                     double bmin[3], double bmax[3],
                                     double circleDim0Coordinate) {
  double intersection[2];
  double circleCenter[2];
  double circleRadius;
  double pointOnPlane[2];
  if (sphereToCircle > r) {
    return -1;
  }
  circleCenter[0] = c[dim1];
  circleCenter[1] = c[dim2];
  circleRadius = sqrt(r * r - sphereToCircle * sphereToCircle);
  pointOnPlane[0] = point[dim1];
  pointOnPlane[1] = point[dim2];
  double cur2dDist =
      intersectCircle(circleCenter, circleRadius, pointOnPlane, intersection);
  if (intersectionInCube2D(intersection, bmin[dim1], bmax[dim1], bmin[dim2],
                           bmax[dim2])) {
    return sqrt((point[dim0] - circleDim0Coordinate) *
                    (point[dim0] - circleDim0Coordinate) +
                cur2dDist * cur2dDist);
  }
  return -1;
}

double SphericalPatch::pointToSphere(double r, double c[3], double point[3],
                                     double bmin[3], double bmax[3]) {
  double d = euclideanDist3D(point, c);
  if (d == 0) {
    return r;
  }
  // if ray sphere intersection is inside cube, return distance from point to
  // sphere.
  double intersection[3];
  intersection[0] = (point[0] - c[0]) * r / d + c[0];
  intersection[1] = (point[1] - c[1]) * r / d + c[1];
  intersection[2] = (point[2] - c[2]) * r / d + c[2];
  if (intersectionInCube(intersection, bmin, bmax)) {
    return fabs(r - d);
  }
  return -1;
}

double SphericalPatch::pointToLineSegment(double r, double c[3],
                                          double point[3], double bmin[3],
                                          double bmax[3], double vertex[3],
                                          int rayD) {
  // will do a ray sphere intersection test
  // B = 2 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc))
  // C = (X0 - Xc)^2 + (Y0 - Yc)^2 + (Z0 - Zc)^2 - Sr^2
  double B = (vertex[rayD] - c[rayD]);
  double C = (vertex[0] - c[0]) * (vertex[0] - c[0]) +
             (vertex[1] - c[1]) * (vertex[1] - c[1]) +
             (vertex[2] - c[2]) * (vertex[2] - c[2]) - r * r;
  double det = 4 * (B * B - C);
  if (det < 0) {
    return -1;
  }
  double sqrtdet2 = sqrt(det) / 2;
  double t1 = (-B - sqrtdet2);
  double t2 = (-B + sqrtdet2);
  ;
  double intersection[3];
  intersection[0] = vertex[0];
  intersection[1] = vertex[1];
  intersection[2] = vertex[2];
  double curDist = -1;
  intersection[rayD] = vertex[rayD] + t1;
  if (intersection[rayD] >= bmin[rayD] && intersection[rayD] <= bmax[rayD]) {
    curDist = euclideanDist3D(point, intersection);
  }
  intersection[rayD] = vertex[rayD] + t2;
  if (intersection[rayD] >= bmin[rayD] && intersection[rayD] <= bmax[rayD]) {
    double d = euclideanDist3D(point, intersection);
    if (curDist > d) {
      curDist = d;
    }
  }
  return curDist;
}

double SphericalPatch::getClosestDist(double r, double c[3], double point[3],
                                      double bmin[3], double bmax[3]) {
  if (c[0] == point[0] && c[0] == point[0] && c[0] == point[0]) {
    return r;
  }
  double dist = -1;
  dist = pointToSphere(r, c, point, bmin, bmax);
  if (dist >= 0) {
    return dist;
  }
  // need to check each face. get circle if any. Find min dist. If any inside,
  // return that. circle on xmin
  double curDist;
  {
    int dim0 = 0, dim1 = 1, dim2 = 2;
    double sphereToCircle = fabs(c[dim0] - bmin[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmin[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // circle on xmax
  {
    int dim0 = 0, dim1 = 1, dim2 = 2;
    double sphereToCircle = fabs(c[dim0] - bmax[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmax[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // circle on ymin
  {
    int dim0 = 1, dim1 = 0, dim2 = 2;
    double sphereToCircle = fabs(c[dim0] - bmin[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmin[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // circle on ymax
  {
    int dim0 = 1, dim1 = 0, dim2 = 2;
    double sphereToCircle = fabs(c[dim0] - bmax[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmax[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // circle on zmin
  {
    int dim0 = 1, dim1 = 2, dim2 = 0;
    double sphereToCircle = fabs(c[dim0] - bmin[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmin[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // circle on zmax
  {
    int dim0 = 1, dim1 = 2, dim2 = 0;
    double sphereToCircle = fabs(c[dim0] - bmax[dim0]);
    curDist = pointToCircle(dim0, dim1, dim2, sphereToCircle, r, c, point, bmin,
                            bmax, bmax[dim0]);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // intersect with all 12 edges
  //
  double vertex[3];
  // xmin, ymin, zmin 1 0 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmin[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 0);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymax, zmin 1 0 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmax[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 0);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymax, zmax 1 0 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmax[1];
    vertex[2] = bmax[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 0);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymax, zmax 1 0 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmin[1];
    vertex[2] = bmax[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 0);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymin, zmin 0 1 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmin[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 1);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmax, ymin, zmin 0 1 0
  {
    vertex[0] = bmax[0];
    vertex[1] = bmin[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 1);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmax, ymin, zmax 0 1 0
  {
    vertex[0] = bmax[0];
    vertex[1] = bmin[1];
    vertex[2] = bmax[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 1);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymin, zmax 0 1 0
  {
    vertex[0] = bmin[0];
    vertex[1] = bmin[1];
    vertex[2] = bmax[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 1);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymin, zmin 0 0 1
  {
    vertex[0] = bmin[0];
    vertex[1] = bmin[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 2);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmax, ymin, zmin 0 0 1
  {
    vertex[0] = bmax[0];
    vertex[1] = bmin[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 2);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmax, ymax, zmin 0 0 1
  {
    vertex[0] = bmax[0];
    vertex[1] = bmax[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 2);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  // xmin, ymax, zmin 0 0 1
  {
    vertex[0] = bmin[0];
    vertex[1] = bmax[1];
    vertex[2] = bmin[2];
    curDist = pointToLineSegment(r, c, point, bmin, bmax, vertex, 2);
    if ((dist < 0) || (dist > curDist && curDist >= 0)) {
      dist = curDist;
    }
  }
  return dist;
}
