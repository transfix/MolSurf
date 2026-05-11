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

#ifndef _PG_H_
#define _PG_H_

#include <vector>
#include <DPG/1D-IntegerRange.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#define INITHASHSIZE 51

namespace DPG {

class Point {
protected:
  double x;
  double y;
  double z;

public:
  int id;

  Point() {
    x = 0;
    y = 0;
    z = 0;
  }

  Point(const Point &p) {
    x = p.x;
    y = p.y;
    z = p.z;
  }

  Point(double a, double b, double c) {
    x = a;
    y = b;
    z = c;
  }

  ~Point() {}

  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }

  void setX(double a) { x = a; }
  void setY(double a) { y = a; }
  void setZ(double a) { z = a; }

  void translate(double a, double b, double c) {
    x += a;
    y += b;
    z += c;
  }
  void translate(double a) {
    x += a;
    y += a;
    z += a;
  }

  bool operator==(const Point &p) { return (x == p.x && y == p.y && z == p.z); }

  double distsq(Point *b) {
    double res;
    double dx, dy, dz;
    dx = x - b->x;
    dy = y - b->y;
    dz = z - b->z;
    res = dx * dx + dy * dy + dz * dz;
    return res;
  }

  double distance(Point *b) { return sqrt(distsq(b)); }

  string toString() // used to generate key for hashing
  {
    ostringstream out("");

    double xmod, ymod, zmod;
    if (x > 0.0)
      xmod = x * 10.0 + 0.5;
    else
      xmod = x * 10.0 - 0.5;

    if (y > 0.0)
      ymod = y * 10.0 + 0.5;
    else
      ymod = y * 10.0 - 0.5;

    if (z > 0.0)
      zmod = z * 10.0 + 0.5;
    else
      zmod = z * 10.0 - 0.5;

    out << (int)(xmod) << " ";

    out << (int)(ymod) << " ";

    out << (int)(zmod);

    return out.str();
  }
};

class WeightedPoint : public Point {
protected:
  double weight;

public:
  WeightedPoint() : Point() { weight = 0; }

  WeightedPoint(const WeightedPoint &p) {
    x = p.x;
    y = p.y;
    z = p.z;
    weight = p.weight;
  }

  WeightedPoint(double w) : Point() { weight = w; }

  WeightedPoint(double a, double b, double c, double w) : Point(a, b, c) {
    weight = w;
  }

  double getWeight() { return weight; }
  void setWeight(double w) { weight = w; }

  double powerDistance(Point *b) { return sqrt(distsq(b) - weight); }

  bool operator==(const WeightedPoint &p) {
    return (x == p.x && y == p.y && z == p.z && weight == p.weight);
  }
};

class Ball : public WeightedPoint {
protected:
  double radius;

public:
  Ball() : WeightedPoint() { radius = 0; }

  Ball(const Ball &p) {
    x = p.x;
    y = p.y;
    z = p.z;
    weight = p.weight;
    radius = p.radius;
  }

  Ball(double r) : WeightedPoint(r * r) { radius = r; }

  Ball(double a, double b, double c, double r) : WeightedPoint(a, b, c, r * r) {
    radius = r;
  }

  double getRadius() { return radius; }
  void setRadius(double r) { radius = r; }

  double surfaceDistance(Ball *b) // computes distance from surface of one ball
                                  // to the surface of the other
  {
    return distance(b) - radius - b->radius;
  }

  bool doesIntersect(Ball *b) {
    double dist = surfaceDistance(b);
    if (dist > 0)
      return false;
    return true;
  }

  bool operator==(const Ball &p) {
    return (x == p.x && y == p.y && z == p.z && radius == p.radius);
  }
};

class gridcell {
public:
  struct cellID {
    int x;
    int y;
    int z;
  };

  cellID ID;
  vector<Point *> balls;

  gridcell(int a, int b, int c) {
    ID.x = a;
    ID.y = b;
    ID.z = c;
  }

  ~gridcell() {
    // cout<<"Inside balls destructor"<<endl;balls.clear();
    balls.clear();
  }
};

class line {
public:
  struct lineID {
    int y;
    int z;
  };

  lineID ID;
  IntegerRange<gridcell *> *RR;

  line(int b, int c) {
    ID.y = b;
    ID.z = c;
    RR = new IntegerRange<gridcell *>();
  }

  ~line() {
    // cout<<"Inside line destructor"<<endl;
  }
};

class plane {
public:
  struct planeID {
    int z;
  };

  planeID ID;
  IntegerRange<line *> *RR;

  plane(int c) {
    ID.z = c;
    RR = new IntegerRange<line *>();
  }

  ~plane() {
    // cout<<"Inside plane destructor"<<endl;
  }
};

class grid {
public:
  IntegerRange<plane *> *RR;

  grid() { RR = new IntegerRange<plane *>(); }

  ~grid() {
    // cout<<"Inside grid destructor"<<endl;
  }
};

class PG {
private:
  grid g;

  double minx, miny, minz, maxx, maxy, maxz; // bounding box

  double rmax;
  double DIM;       // Dimension of each grid in DPG. DIM should be O(r_{max}).
  double TRANSLATE; // The underlying range reporting only works with positive
                    // integers. So, to ensure all grids and points have
                    // positive coordinates we translate all points by a fixed
                    // margin. It is generally good to be conservative and
                    // translate more than the most negative input coordinate

  void computeTranslate() {
    double min = minx;
    if (miny < min)
      min = miny;
    if (minz < min)
      min = minz;

    double xDim = maxx - minx;
    double yDim = maxy - miny;
    double zDim = maxz - minz;

    double maxDim = xDim;
    if (maxDim < yDim)
      maxDim = yDim;
    if (maxDim < zDim)
      maxDim = zDim;

    TRANSLATE = maxDim;
    if (min < 0)
      TRANSLATE -= min;
  }

public:
  PG(double D, double xlate,
     double r) // use when r_{max} and xlate is already known
  {
    DIM = D;
    rmax = r;
    TRANSLATE = xlate;
  }

  PG(double mnx, double mny, double mnz, double mxx, double mxy, double mxz,
     double divisionsize,
     int inithashsize =
         INITHASHSIZE) // useful when bounding box is known.
                       // but the actual coordinates are not known yet.
  {
    maxx = mxx;
    maxy = mxy;
    maxz = mxz;
    minx = mnx;
    miny = mny;
    minz = mnz;

    computeTranslate();

    if (divisionsize != 0)
      DIM = divisionsize;
    else
      DIM = 3.5;

    rmax = 1.5;
  }

  PG(vector<Point *>
         *alist) // dpg itself will compute the bounding box and fix the
                 // TRANSLATE and r_{max} values. and adds the points
  {
    int numPoints = alist->size();

    Point *p = alist->at(0);

    double x, y, z;
    x = p->getX();
    y = p->getY();
    z = p->getZ();

    minx = x;
    miny = y;
    minz = z;
    maxx = x;
    maxy = y;
    maxz = z;

    for (int i = 0; i < numPoints; i++) {
      Point *p = alist->at(i);
      x = p->getX();
      y = p->getY();
      z = p->getZ();

      if (x < minx)
        minx = x;
      if (y < miny)
        miny = y;
      if (z < minz)
        minz = z;

      if (x > minx)
        maxx = x;
      if (y > maxy)
        maxy = y;
      if (z > maxz)
        maxz = z;
    }

    computeTranslate();

    DIM = 3.5;
    rmax = 1.5;

    addPoints(alist);
  }

  void boundingbox(double &mnx, double &mny, double &mnz, double &mxx,
                   double &mxy, double &mxz) {
    mnx = this->minx;
    mny = this->miny;
    mnz = this->minz;
    mxx = this->maxx;
    mxy = this->maxy;
    mxz = this->maxz;
  }

  /* get methods */

  double getdivsize() { return DIM; }
  double getRmax() { return rmax; }

  /* Supported Queries*/

  bool range(Point *q, double delta,
             vector<Point *>
                 &result); // returns all points within a distance delta from q
  vector<Point *> range(Point *q, double delta);
  bool doesOverlap(Ball *q, double delta);
  bool pointsWithinRange(Point *q,
                         double delta); // returns true if there are points
                                        // within a distance delta from q
  bool
  lightRange(Point *q, double delta,
             vector<Point *> &
                 result); // returns all points within a distance delta+DIM from
                          // q. this avoids pairwise distance computations.
  vector<Point *> lightRange(Point *q, double delta);

  /* Supported updates*/

  void addPoint(Point *a);
  void removePoint(Point *a);

  void addPoints(vector<Point *> *alist) {
    int i, size;
    size = alist->size();

    for (i = 0; i < size; i++)
      addPoint(alist->at(i));
  }

  void removePoints(vector<Point *> *alist) {
    int i, size;
    size = alist->size();

    for (i = 0; i < size; i++)
      removePoint(alist->at(i));
  }

  void move(Point *oldP, Point *newP) {
    removePoint(oldP);
    addPoint(newP);
  }
}; // end definition of DPG

} // namespace DPG

#endif
