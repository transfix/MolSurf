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
#ifndef DISTANCE_TRANSFORM_V2_
#define DISTANCE_TRANSFORM_V2_

#include <Utility/utility.h>
#include <SignDistanceFunction_v2/FaceVertSet3D.h>
#include <SignDistanceFunction_v2/RawivParser.h>
#include <SignDistanceFunction_v2/reg3data.h>
#include <SignDistanceFunction_v2/mtxlib.h>

class DistanceTransformCell {
public:
  DistanceTransformCell() {}

  ~DistanceTransformCell() {}

  bool hasTriangle() const { return triList.length() > 0; }

  dynamic_array<int> triList;
};

class DistanceTransform {
  Reg3Data<float> *p_Data;
  DistanceTransformCell *p_DistanceTransformCells;
  FaceVertSet3D *p_Surf;

public:
  // sx, sy, sz: scale factor of the bounding volume
  DistanceTransform(FaceVertSet3D &fvs, int dim[3], float Distance = 20,
                    float sx = 2, float sy = 2, float sz = 2);

  // define the grid explicitly
  DistanceTransform(FaceVertSet3D &fvs, int dim[3], float bbox[3],
                    float spacing[3]);

  DistanceTransform(FaceVertSet3D &fvs, const Reg3Data<float> &reg3);

  ~DistanceTransform(void);

  void transform();

  void writeRawiv(const char *fname) { RawivParser().write(*p_Data, fname); }

  void transform1D(int len, float f[], float d[], int parent[], float span);
  const static float MAX_FLOAT;
  bool intersectDistanceTransformCell(const Point3f &v0, const Vector3f &norm,
                                      int i, int j, int k);
  bool nearSurface(int i, int j, int k);
  float computeNearDistance(int ix, int iy, int iz, Point3f &pnt);
  double distance2Triangle(const Point3f &pnt, int nt, Point3f &ne);
  double distance2Edge(const Point3f &pnt, const Point3f &v1, const Point3f &v2,
                       Point3f &ne);
  bool pointInTriangle(const Point3f &res, const Point3f &v0, const Point3f &v1,
                       const Point3f &v2, const Vector3f &norm);
  int nearestPlane(const dynamic_array<int> &triList, const Point3f &vert);
  double rayTriangleIntersection(int nt, const Point3f &begin,
                                 const Point3f &end);
  // return -1, if point is inside, return 1 if point is outside
  int inOrOut(const dynamic_array<int> &triList, const Point3f &vert,
              const Point3f &nearPnt);
  // return 1 if they intersect
  bool TriangleCubeIntersection(int nt, const Point3f &lower,
                                const Point3f &upper);
  const Reg3Data<float> &getReg3Data() { return *p_Data; }
  const Reg3Data<float> *getData() { return p_Data; }

private:
  void init();
};

#endif
