/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <DynamicAdaptiveGrid/GridPoint.h>

using namespace DynamicAdaptiveGrid;

int GridPoint::gridPointIDGenerator = 0;

GridPoint::GridPoint() : Point(0, 0, 0) {
  functionValue = 0;
  state = OUT_POINT;
  hasVal = false;
  generateGridPointID();

  paircount = 0;
  triplecount = 0;
  pairprob = 0;
  tripleprob = 0;
  pairtestcount = 0;
  tripletestcount = 0;

  defaultMode = 2;

  markAsNoVal();
}

GridPoint::GridPoint(const GridPoint &gp) : Point(gp.x, gp.y, gp.z) {
  functionValue = gp.functionValue;
  state = gp.state;
  hasVal = gp.hasVal;
  generateGridPointID();

  paircount = 0;
  triplecount = 0;
  pairprob = 0;
  tripleprob = 0;
  pairtestcount = 0;
  tripletestcount = 0;

  defaultMode = 2;

  markAsNoVal();
}

GridPoint::GridPoint(double a, double b, double c) : Point(a, b, c) {
  functionValue = 0;
  state = OUT_POINT;
  hasVal = false;
  generateGridPointID();

  paircount = 0;
  triplecount = 0;
  pairprob = 0;
  tripleprob = 0;
  pairtestcount = 0;
  tripletestcount = 0;

  defaultMode = 2;

  markAsNoVal();
}

GridPoint::GridPoint(double a, double b, double c, double F) : Point(a, b, c) {
  functionValue = F;
  state = OUT_POINT;
  hasVal = false;
  generateGridPointID();

  paircount = 0;
  triplecount = 0;
  pairprob = 0;
  tripleprob = 0;
  pairtestcount = 0;
  tripletestcount = 0;

  defaultMode = 2;

  markAsNoVal();
}

GridPoint::GridPoint(double a, double b, double c, double F, int s)
    : Point(a, b, c) {
  functionValue = F;
  state = s;
  hasVal = false;
  generateGridPointID();

  paircount = 0;
  triplecount = 0;
  pairprob = 0;
  tripleprob = 0;
  pairtestcount = 0;
  tripletestcount = 0;

  defaultMode = 2;

  markAsNoVal();
}

void GridPoint::generateGridPointID() { gridPointID = gridPointIDGenerator++; }

bool GridPoint::hasValue() { return hasVal; }

int GridPoint::getState() { return state; }

double GridPoint::getFunctionValue(int mode) {
  if (mode == 0) {
    if (hasVDW)
      return functionValue;
    else {
      vdwAtPoint();
      return functionValue;
    }
  } else if (mode == 1) {
    if (hasSES)
      return functionValue;
    else {
      sesAtPoint();
      return functionValue;
    }
  } else if (mode == 2) {
    if (hasGauss)
      return functionValue;
    else {
      gaussAtPoint();
      return functionValue;
    }
  }
}

double GridPoint::getFunctionValue() { return getFunctionValue(defaultMode); }

int GridPoint::getGridPointID() { return gridPointID; }

void GridPoint::setState(int s) {
  if (s < UNDEFINED_POINT)
    state = s;
  else
    state = OUT_POINT;
}

void GridPoint::setFunctionValue(double f) { functionValue = f; }

void GridPoint::addToFunctionValue(double f) { functionValue += f; }

void GridPoint::setPosition(double a, double b, double c) {
  x = a;
  y = b;
  z = c;
}

void GridPoint::setGridCells(vector<Point *> gcs) {
  gridCells.clear();

  for (int i = 0; i < gcs.size(); i++)
    gridCells.push_back((GridCell *)gcs[i]);
}

vector<GridCell *> GridPoint::getGridCells() { return gridCells; }

GridCell *GridPoint::getGridCell(int i) {
  if (i < gridCells.size())
    return gridCells[i];
  return NULL;
}

int GridPoint::getNumGridCells() { return gridCells.size(); }

void GridPoint::vdwAtPoint() {
  vector<DPG::Point *> atms;
  exposedAtomPG->range(this, solventRad + 3.0, atms);

  double mydist = -50.0;
  bool inside = false;
  int size = atms.size();
  double distsq;
  double dx, dy, dz, r, dist;

  for (int i = 0; i < size; i++) {
    DPG::Ball *a = (DPG::Ball *)atms[i];
    dx = a->getX() - x;
    dy = a->getY() - y;
    dz = a->getZ() - z;
    r = a->getRadius() + solventRad;

    distsq = dx * dx + dy * dy + dz * dz;

    dist = r - sqrt(distsq);

    if (dist > mydist)
      mydist = dist;
  }

  functionValue = mydist;
  hasVal = true;
}

void GridPoint::sesAtPoint() {
  vector<DPG::Point *> atms;
  exposedAtomPG->range(this, solventRad + 3.0, atms);

  int size = atms.size();

  double sdf = 50.0;
  double dist;

  double singleSDF = 50.0;
  double pairSDF = 50.0;
  double tripleSDF = 50.0;

  bool insideVDW = false;
  bool inside = false;
  bool inCone = false;
  bool inTet = false;

  bool *insideAtom = new bool[size];
  bool **pairIntersect = new bool *[size];
  int i, j;
  double dx, dy, dz, dist2;
  double R1, R2, R3, R;
  DPG::Ball *b = new DPG::Ball();
  b->setRadius(solventRad);

  for (i = 0; i < size; i++) {
    pairIntersect[i] = new bool[size];

    for (j = 0; j < size; j++)
      pairIntersect[i][j] = false;

    pairIntersect[i][i] = true;
  }

  for (i = 0; i < size; i++) {
    DPG::Ball *a = (DPG::Ball *)atms[i];

    dx = a->getX() - x;
    dy = a->getY() - y;
    dz = a->getZ() - z;
    R = a->getRadius() + solventRad;

    dist2 = sqrt(dx * dx + dy * dy + dz * dz);
    dist = R - dist2;

    if (dist > 0) // ignore points outside SAS.
    {
      if (dist < singleSDF) {
        double ratio = R / dist2;

        b->setX(a->getX() + dx * ratio);
        b->setY(a->getY() + dy * ratio);
        b->setZ(a->getZ() + dz * ratio);

        if (!exposedAtomPG->doesOverlap(b, 2.0))
          singleSDF = dist;
      }
      inside = true;
      insideAtom[i] = true;

      if (dist > solventRad) // inside vdw
      {
        insideVDW = true;
      }
    }
  }

  // Checking pairwise intersection cases. Finding closest distance to the SAS.

  if (inside) {
    for (int i = 0; i < size; i++) {
      for (int j = i + 1; j < size; j++) {
        DPG::Ball *a1 = (DPG::Ball *)atms[i];
        DPG::Ball *a2 = (DPG::Ball *)atms[j];
        R1 = a1->getRadius() + solventRad;
        R2 = a2->getRadius() + solventRad;
        R = R1 + R2;

        pairtestcount++;

        if (a1->distsq(a2) < R * R) {
          pairIntersect[i][j] = true;
          pairIntersect[j][i] = true;
          double tempsdf;

          if (getSDFForIntersectingSpheres(a1, a2, x, y, z, solventRad,
                                           tempsdf)) {

            if (tempsdf > R1 && tempsdf > R2) {
              cout << "In cone, but distance is too large " << tempsdf << endl;
              pairprob++;
            }

            else if (tempsdf < pairSDF) {
              pairSDF = tempsdf;
              inCone = true;
            }

            paircount++;
          }

        } // end if double intersect

      } // end for j

    } // end for i

  } // end if inside

  // Checking triple intersection cases. Finding closest distance to the SAS.
  /*
          if(inside)
          {
                  for(int i=0; i<size; i++)
                  {
                          for(int j=i+1; j<size; j++)
                          {
                                  for(int k=j+1;k<size;k++)
                                  {
                                          DPG::Ball* a1 = (DPG::Ball*)atms[i];
                                          DPG::Ball* a2 = (DPG::Ball*)atms[j];
                                          DPG::Ball* a3 = (DPG::Ball*)atms[k];
                                          R1 = a1->getRadius() + solventRad;
                                          R2 = a2->getRadius() + solventRad;
                                          R3 = a3->getRadius() + solventRad;
                                          R = R1 + R2;

                                          double tempsdf;


                                          if(pairIntersect[i][j] &&
     pairIntersect[j][k] && pairIntersect[i][k])	// all pairs intersect.
     so possible triple intersection
                                          {
                                                  tripletestcount++;

                                                  if(getPointToCircleDistance(a1,
     a2, a3->getX(), a3->getY(), a3->getZ(), solventRad) <= R3)
                                                  {
                                                          if(getSDFForTripleIntersect(a1,
     a2, a3, x, y, z, solventRad, tempsdf))
                                                          {
                                                                  if(tempsdf >
     R1  && tempsdf > R2  && tempsdf > R3)
                                                                  {
                                                                          tripleprob++;
                                                                          cout<<"Problem
     in SDF: "<<tempsdf<<endl;
                                                                  }

                                                                  else
     if(tripleSDF > tempsdf)
                                                                  {
                                                                          tripleSDF
     = tempsdf; inTet = true;
                                                                  }

                                                                  triplecount++;
                                                          }
                                                  }
                                          }

                                  }	//end for k

                          }// end for j

                  }// end for i

          }//end if inside
  */
  if (inside && sdf > singleSDF)
    sdf = singleSDF;
  if (inCone && sdf > pairSDF)
    sdf = pairSDF;
  if (inTet && sdf > tripleSDF)
    sdf = tripleSDF;

  if (!inside) {
    sdf = 0.0;
  }

  functionValue = sdf;
  hasSES = true;
  hasVal = true;

  delete[] insideAtom;
  for (int i = 0; i < size; i++)
    delete[] pairIntersect[i];
  delete[] pairIntersect;
}

double GridPoint::getPointToCircleDistance(
    DPG::Ball *a1, DPG::Ball *a2, double x, double y, double z,
    double solventRad) // circle specified by two intersecting spheres
{
  float R1 = a1->getRadius() + solventRad;
  float R2 = a2->getRadius() + solventRad;

  float dx = a1->getX() - a2->getX();
  float dy = a1->getY() - a2->getY();
  float dz = a1->getZ() - a2->getZ();

  float distSq = (dx * dx + dy * dy + dz * dz);
  float dist = sqrt(distSq);

  float R1Sq = R1 * R1;
  float R2Sq = R2 * R2;

  float vecC2C1x = dx / dist; // normalized vector pointing from center of atom
                              // 2 to center of atom 1
  float vecC2C1y = dy / dist;
  float vecC2C1z = dz / dist;

  float A = vecC2C1x; // equation of plane containing the circle
  float B = vecC2C1y;
  float C = vecC2C1z;
  float D =
      ((a2->getX() * a2->getX() - a1->getX() * a1->getX()) +
       (a2->getY() * a2->getY() - a1->getY() * a1->getY()) +
       (a2->getZ() * a2->getZ() - a1->getZ() * a1->getZ()) + (R1Sq - R2Sq)) /
      (2 * dist);

  float pointPlaneDist = A * x + B * y + C * z + D;
  float projectedX = x - A * pointPlaneDist;
  float projectedY = y - B * pointPlaneDist;
  float projectedZ = z - C * pointPlaneDist;

  float distanceFromCenter1ToCircle =
      (R1Sq - R2Sq + distSq) /
      (2 *
       dist); // solution of R1^2 - d1^2 = R2^2 - d2^2 = r^2 and d1+d2 = dist
  float distanceFromCenter2ToCircle = dist - distanceFromCenter1ToCircle;

  float circleCenterX = (a1->getX() * distanceFromCenter2ToCircle +
                         a2->getX() * distanceFromCenter1ToCircle) /
                        dist;
  float circleCenterY = (a1->getY() * distanceFromCenter2ToCircle +
                         a2->getY() * distanceFromCenter1ToCircle) /
                        dist;
  float circleCenterZ = (a1->getZ() * distanceFromCenter2ToCircle +
                         a2->getZ() * distanceFromCenter1ToCircle) /
                        dist;

  float circleRad =
      sqrt(R1Sq - distanceFromCenter1ToCircle * distanceFromCenter1ToCircle);

  float onPlaneDirFromCCX = projectedX - circleCenterX;
  float onPlaneDirFromCCY = projectedY - circleCenterY;
  float onPlaneDirFromCCZ = projectedZ - circleCenterZ;

  float onPlaneDistanceToCircleCenter =
      sqrt(onPlaneDirFromCCX * onPlaneDirFromCCX +
           onPlaneDirFromCCY * onPlaneDirFromCCY +
           onPlaneDirFromCCZ * onPlaneDirFromCCZ);
  float ratio =
      circleRad /
      onPlaneDistanceToCircleCenter; // possible numerical problem when the
                                     // projection is close to the center. have
                                     // to fix this.

  float closestPointOnCircleX = circleCenterX + ratio * onPlaneDirFromCCX;
  float closestPointOnCircleY = circleCenterY + ratio * onPlaneDirFromCCY;
  float closestPointOnCircleZ = circleCenterZ + ratio * onPlaneDirFromCCZ;

  float dPtX = x - closestPointOnCircleX;
  float dPtY = y - closestPointOnCircleY;
  float dPtZ = z - closestPointOnCircleZ;

  return sqrt(dPtX * dPtX + dPtY * dPtY + dPtZ * dPtZ);
}

bool GridPoint::getSDFForIntersectingSpheres(DPG::Ball *a1, DPG::Ball *a2,
                                             double x, double y, double z,
                                             double solventRad, double &sdf) {
  float R1 = a1->getRadius() + solventRad;
  float R2 = a2->getRadius() + solventRad;

  float dx = a1->getX() - a2->getX();
  float dy = a1->getY() - a2->getY();
  float dz = a1->getZ() - a2->getZ();

  float dx1 = a1->getX() - x;
  float dy1 = a1->getY() - y;
  float dz1 = a1->getZ() - z;

  float dx2 = a2->getX() - x;
  float dy2 = a2->getY() - y;
  float dz2 = a2->getZ() - z;

  float distSq = (dx * dx + dy * dy + dz * dz);
  float dist = sqrt(distSq);

  float R1Sq = R1 * R1;
  float R2Sq = R2 * R2;

  float distanceFromCenter1ToCircle =
      (R1Sq - R2Sq + distSq) /
      (2 *
       dist); // solution of R1^2 - d1^2 = R2^2 - d2^2 = r^2 and d1+d2 = dist
  float distanceFromCenter2ToCircle = dist - distanceFromCenter1ToCircle;

  float cosTheta1 =
      distanceFromCenter1ToCircle / R1; // half-angles subtended from the
                                        // centers to the circle of intersection
  float cosTheta2 = distanceFromCenter2ToCircle / R2;

  float vecC2C1x = dx / dist; // normalized vector pointing from center of atom
                              // 2 to center of atom 1
  float vecC2C1y = dy / dist;
  float vecC2C1z = dz / dist;

  float dist1Sq = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  float dist1 = sqrt(dist1Sq);

  float vecPtC1x = dx1 / dist1; // normalized vector pointing from the given
                                // Point to center of atom 1
  float vecPtC1y = dy1 / dist1;
  float vecPtC1z = dz1 / dist1;

  float dist2Sq = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  float dist2 = sqrt(dist2Sq);

  float vecPtC2x = dx2 / dist2; // normalized vector pointing from the given
                                // point to center of atom 2
  float vecPtC2y = dy2 / dist2;
  float vecPtC2z = dz2 / dist2;

  float cosThetaPt1 =
      vecPtC1x * vecC2C1x + vecPtC1y * vecC2C1y +
      vecPtC1z *
          vecC2C1z; // angle between the axis and the line connecting Pt and C1
  float cosThetaPt2 = -(
      vecPtC2x * vecC2C1x + vecPtC2y * vecC2C1y +
      vecPtC2z *
          vecC2C1z); // angle between the axis and the line connecting Pt and C2

  if (cosThetaPt1 > cosTheta1 &&
      cosThetaPt2 >
          cosTheta2) // inside both cones. SDF is distance from Pt to circle
  {

    float A = vecC2C1x; // equation of plane containing the circle
    float B = vecC2C1y;
    float C = vecC2C1z;
    float D =
        ((a2->getX() * a2->getX() - a1->getX() * a1->getX()) +
         (a2->getY() * a2->getY() - a1->getY() * a1->getY()) +
         (a2->getZ() * a2->getZ() - a1->getZ() * a1->getZ()) + (R1Sq - R2Sq)) /
        (2 * dist);

    float pointPlaneDist = A * x + B * y + C * z + D;
    float projectedX = x - A * pointPlaneDist;
    float projectedY = y - B * pointPlaneDist;
    float projectedZ = z - C * pointPlaneDist;

    float circleCenterX = (a1->getX() * distanceFromCenter2ToCircle +
                           a2->getX() * distanceFromCenter1ToCircle) /
                          dist;
    float circleCenterY = (a1->getY() * distanceFromCenter2ToCircle +
                           a2->getY() * distanceFromCenter1ToCircle) /
                          dist;
    float circleCenterZ = (a1->getZ() * distanceFromCenter2ToCircle +
                           a2->getZ() * distanceFromCenter1ToCircle) /
                          dist;

    float circleRad =
        sqrt(R1Sq - distanceFromCenter1ToCircle * distanceFromCenter1ToCircle);

    float onPlaneDirFromCCX = projectedX - circleCenterX;
    float onPlaneDirFromCCY = projectedY - circleCenterY;
    float onPlaneDirFromCCZ = projectedZ - circleCenterZ;

    float onPlaneDistanceToCircleCenter =
        sqrt(onPlaneDirFromCCX * onPlaneDirFromCCX +
             onPlaneDirFromCCY * onPlaneDirFromCCY +
             onPlaneDirFromCCZ * onPlaneDirFromCCZ);
    float ratio =
        circleRad /
        onPlaneDistanceToCircleCenter; // possible numerical problem when the
                                       // projection is close to the center.
                                       // have to fix this.

    float closestPointOnCircleX = circleCenterX + ratio * onPlaneDirFromCCX;
    float closestPointOnCircleY = circleCenterY + ratio * onPlaneDirFromCCY;
    float closestPointOnCircleZ = circleCenterZ + ratio * onPlaneDirFromCCZ;

    DPG::Ball *b = new DPG::Ball(closestPointOnCircleX, closestPointOnCircleY,
                                 closestPointOnCircleZ, solventRad);

    if (!exposedAtomPG->doesOverlap(b, 2.0)) {
      float dPtX = x - closestPointOnCircleX;
      float dPtY = y - closestPointOnCircleY;
      float dPtZ = z - closestPointOnCircleZ;

      sdf = sqrt(dPtX * dPtX + dPtY * dPtY + dPtZ * dPtZ);

      return true; // reentrant case found
    }
  }

  return false;
}

bool GridPoint::getSDFForTripleIntersect(DPG::Ball *a1, DPG::Ball *a2,
                                         DPG::Ball *a3, double x, double y,
                                         double z, double solventRad,
                                         double &sdf) {
  float R1 = a1->getRadius() + solventRad;
  float R2 = a2->getRadius() + solventRad;
  float R3 = a3->getRadius() + solventRad;

  // computing distances and the planes. Two planes are enough
  float dist12, dist32;
  float plane12[4], plane32[4];

  float dx, dy, dz;

  dx = a1->getX() - a2->getX();
  dy = a1->getY() - a2->getY();
  dz = a1->getZ() - a2->getZ();
  dist12 = sqrt(dx * dx + dy * dy + dz * dz);

  plane12[0] = dx / dist12;
  plane12[1] = dy / dist12;
  plane12[2] = dz / dist12;
  plane12[3] =
      (a2->getX() * a2->getX() - a1->getX() * a1->getX() +
       a2->getY() * a2->getY() - a1->getY() * a1->getY() +
       a2->getZ() * a2->getZ() - a1->getZ() * a1->getZ() + R1 * R1 - R2 * R2) /
      (2.0 * dist12);

  dx = a3->getX() - a2->getX();
  dy = a3->getY() - a2->getY();
  dz = a3->getZ() - a2->getZ();
  dist32 = sqrt(dx * dx + dy * dy + dz * dz);

  plane32[0] = dx / dist32;
  plane32[1] = dy / dist32;
  plane32[2] = dz / dist32;
  plane32[3] =
      (a2->getX() * a2->getX() - a3->getX() * a3->getX() +
       a2->getY() * a2->getY() - a3->getY() * a3->getY() +
       a2->getZ() * a2->getZ() - a3->getZ() * a3->getZ() + R3 * R3 - R2 * R2) /
      (2.0 * dist32);

  // intersecting the planes to get a line
  float cross[3];
  cross[0] = plane12[1] * plane32[2] - plane12[2] * plane32[1];
  cross[1] = plane12[2] * plane32[0] - plane12[0] * plane32[2];
  cross[2] = plane12[0] * plane32[1] - plane12[1] * plane32[0];

  float dot[3];
  dot[0] = plane12[0] * plane32[0];
  dot[1] = plane12[1] * plane32[1];
  dot[2] = plane12[2] * plane32[2];

  float pointOnPlane[3];

  if (abs(cross[0]) >= abs(cross[1]) &&
      abs(cross[0]) >= abs(cross[2])) // set x coordinate to 0
  {
    pointOnPlane[0] = 0.0;
    pointOnPlane[1] =
        (plane32[3] * plane12[2] - plane12[3] * plane32[2]) / cross[0];
    pointOnPlane[2] =
        (plane12[3] * plane32[1] - plane32[3] * plane12[1]) / cross[0];
  } else if (abs(cross[1]) >= abs(cross[0]) &&
             abs(cross[1]) >= abs(cross[2])) // set y coordinate to 0
  {
    pointOnPlane[0] =
        (plane12[3] * plane32[2] - plane32[3] * plane12[2]) / cross[1];
    pointOnPlane[1] = 0.0;
    pointOnPlane[2] =
        (plane32[3] * plane12[0] - plane12[3] * plane32[0]) / cross[1];
  } else if (abs(cross[2]) >= abs(cross[1]) &&
             abs(cross[2]) >= abs(cross[0])) // set z coordinate to 0
  {
    pointOnPlane[0] =
        (plane32[3] * plane12[1] - plane12[3] * plane32[1]) / cross[2];
    pointOnPlane[1] =
        (plane12[3] * plane32[0] - plane32[3] * plane12[0]) / cross[2];
    pointOnPlane[2] = 0.0;
  }

  // The parametric equation of the line is pointOnPlane + t*cross, i.e x(t) =
  // pointOnPlane[0] + t*cross[0] ... plug this into the equation of a sphere,
  // and solve for t t^2 (i^2 + j^2 + k^2) + 2t( (x0-x1)i + (y0-y1)j + (z0-z1)k
  // ) + ( (x0-x1)^2 + (y0-y1)^2 + (z0-z1)^2 - r1^2) = 0 where, cross = [i,j,k];
  // pointOnPlane = [x0,y0,z0]; a1.xyz = [x1,y1,z1]; R1 = r1;

  float x01 = pointOnPlane[0] - a1->getX();
  float y01 = pointOnPlane[1] - a1->getY();
  float z01 = pointOnPlane[2] - a1->getZ();

  float a = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
  float b = 2.0 * (x01 * cross[0] + y01 * cross[1] + dz * cross[2]);
  float c = x01 * x01 + y01 * y01 + z01 * z01 - R1 * R1;

  float det = b * b - 4.0 * a * c;
  if (det < 0) {
    return false;
  } else
    det = sqrt(det);

  float t1 = (-b + det) / (2.0 * a);
  float t2 = (-b - det) / (2.0 * a);

  float point1[3] = {pointOnPlane[0] + t1 * cross[0],
                     pointOnPlane[1] + t1 * cross[1],
                     pointOnPlane[2] + t1 * cross[2]};
  float point2[3] = {pointOnPlane[0] + t2 * cross[0],
                     pointOnPlane[1] + t2 * cross[1],
                     pointOnPlane[2] + t2 * cross[2]};

  float givenPoint[3] = {x, y, z};
  float dist1 = 10.0;
  float dist2 = 10.0;
  bool sdfcomputed = false;

  DPG::Ball *b1 = new DPG::Ball(point1[0], point1[1], point1[2], solventRad);
  DPG::Ball *b2 = new DPG::Ball(point2[0], point2[1], point2[2], solventRad);

  if (!exposedAtomPG->doesOverlap(b1, 2.0)) {
    {
      float dx1 = point1[0] - givenPoint[0];
      float dy1 = point1[1] - givenPoint[1];
      float dz1 = point1[2] - givenPoint[2];
      dist1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
      sdfcomputed = true;
    }
  }
  if (!exposedAtomPG->doesOverlap(b2, 2.0)) {
    {
      float dx2 = point2[0] - givenPoint[0];
      float dy2 = point2[1] - givenPoint[1];
      float dz2 = point2[2] - givenPoint[2];
      dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
      sdfcomputed = true;
    }
  }

  sdf = dist1 < dist2 ? dist1 : dist2;

  return sdfcomputed;
}

void GridPoint::gaussAtPoint() {
  vector<DPG::Point *> atms;
  exposedAtomPG->range(this, solventRad + 3.0, atms);

  double val = 0.0;
  double blobby = -2.3;
  int size = atms.size();
  double dx, dy, dz, distsq, r, expval;

  for (int i = 0; i < size; i++) {
    DPG::Ball *a = (DPG::Ball *)atms[i];
    dx = a->getX() - x;
    dy = a->getY() - y;
    dz = a->getZ() - z;

    distsq = dx * dx + dy * dy + dz * dz;
    r = a->getRadius();

    expval = blobby * distsq / (r * r) - blobby;
    val += (exp(expval));
  }

  functionValue = val;
  hasGauss = true;
  hasVal = true;
}

void GridPoint::gaussAtPointSingleAtom(PDBParser::Atom *at, bool sub) {
  double blobby = -2.3;

  double ax = at->getX();
  double ay = at->getY();
  double az = at->getZ();
  double ar = at->getRadius();

  double distsq =
      (x - ax) * (x - ax) + (y - ay) * (y - ay) + (z - az) * (z - az);
  double expval = blobby * distsq / (ar * ar) - blobby;
  double val = (exp(expval));

  if (sub)
    functionValue -= val;
  else
    functionValue += val;
  hasVal = true;
}
