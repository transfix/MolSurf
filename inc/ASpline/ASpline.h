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
#ifndef ASPLINE_MAIN_H
#define ASPLINE_MAIN_H

#include <Utility/utility.h>
#include <ASpline/ASplinePrimitives.h>
#include <ASpline/quadrature.h>
#include <Geometry/Geometry.h>
#include <LBIE_lib/Geoframe.h> //Duplicates functionality we need

const double boundingBoxSize = 1000000000;

namespace ASPLINE {
class Vertex;
class Face;
class ASpline {
public:
  ASpline();
  ASpline(FILE *fp);           // From a .raw type file
  ASpline(Geoframe *geoframe); // From a surface represented as a Geoframe
  ASpline(Geometry *goemetry); // From a surface reprseneted as a Geometry
  Geoframe *getGeoframe();     // Coerce back into a Geoframe
  Geometry *getGeometry();     // Coerce back into a Geometry
  // Craig: move these three following functions here from /MolSurf
  void saveNewRawncFile(FILE *fp, int resolution);
  void printRawnWithQuadrature(const char *fname, int resolution,
                               int numQuadPts);
  void printQuadratureResults(const char *fname);
  void addPoint(Vertex *point);
  void addFacet(Face *facet);
  void Triangle_Bezier_Patch(int resolution);
  void Triangle_Bezier_Patch_Nodes(const char *type, int numOfPts);
  void Triangle_Bezier_Patch_With_Quadrature(const char *type, int numOfPts,
                                             int resolution);
  void Triangle_Patch_Volume();
  // private: // These are used already so just leave them public
  std::vector<Vertex *> verticies;
  std::vector<Face *> faces;
  int numbpts, numbtris;
  double min[3], max[3]; // used for the bounding box
  bool normalFlag, colorFlag;
  double volume;

private:
  void Cubic_Bezier_Coeffs(double *p1, double *p2, double *p3, double *n1,
                           double *n2, double *n3, double *b300, double *b210,
                           double *b120, double *b030, double *b021,
                           double *b012, double *b003, double *b102,
                           double *b201);
  void Simple_Bezier_B111(double *p1, double *p2, double *p3, double *n1,
                          double *n2, double *n3, double *b300, double *b210,
                          double *b120, double *b030, double *b201,
                          double *b021, double *b111);
  void Triangle_Patch(int t, double *p1, double *p2, double *p3, double *n1,
                      double *n2, double *n3, double *b300, double *b210,
                      double *b120, double *b030, double *b021, double *b012,
                      double *b003, double *b102, double *b201, double *b111_1,
                      double *b111_2, double *b111_3, int resolution);
  void Triangle_Patch_Nodes(Quadrature *my_quadr, int i, double *p1, double *p2,
                            double *p3, double *n1, double *n2, double *n3,
                            double *b300, double *b210, double *b120,
                            double *b030, double *b021, double *b012,
                            double *b003, double *b102, double *b201,
                            double *b111_1, double *b111_2, double *b111_3);
  void Triangle_Patch_With_Quadrature(Quadrature *my_quadr, int t, double *p1,
                                      double *p2, double *p3, double *n1,
                                      double *n2, double *n3, double *b300,
                                      double *b210, double *b120, double *b030,
                                      double *b021, double *b012, double *b003,
                                      double *b102, double *b201,
                                      double *b111_1, double *b111_2,
                                      double *b111_3, int resolution);
  double Patch_Volume(Quadrature *my_quadr, int i, double *p1, double *p2,
                      double *p3, double *n1, double *n2, double *n3,
                      double *b300, double *b210, double *b120, double *b030,
                      double *b021, double *b012, double *b003, double *b102,
                      double *b201, double *b111_1, double *b111_2,
                      double *b111_3);
  double BB2lmd(double *b300, double *b210, double *b120, double *b030,
                double *b021, double *b012, double *b003, double *b102,
                double *b201, double *b111_1, double *b111_2, double *b111_3,
                double b1, double b2, double b3);
  double Bezier_LambdaPolyZeroC1(double *F1, double *N1, double *N2, double *N3,
                                 double *D1, double *D2, double *D3);
  double Bezier_ThreeLambdaPolyValueC1(double *F1, double *N1, double *N2,
                                       double *N3, double *D1, double *D2,
                                       double *D3, double lmd);
  double Bezier_ThreeLambdaDervPolyValueC1(double *F1, double *N1, double *N2,
                                           double *N3, double *D1, double *D2,
                                           double *D3, double lmd);
  void BB_Gradient(double *p1_lmd, double *p2_lmd, double *p3_lmd, double *n1,
                   double *n2, double *n3, double *b300, double *b210,
                   double *b120, double *b030, double *b021, double *b012,
                   double *b003, double *b102, double *b201, double *b111_1,
                   double *b111_2, double *b111_3, double b1, double b2,
                   double b3, double lmd, double *n);
  double DotProduct(double *v1, double *v2);
  void CrossProduct(double *v1, double *v2, double *v);
  double Determinant(double *v1, double *v2, double *v3);
  double Triangle_area(double *v1, double *v2, double *v3);
  void PolyDerivative(double *p, int n, double *dp);
  void PolyProduct(double *p1, int m, double *p2, int n, double *p3);
  double HornerEvaluate(double *p, int n, double t);
  void BoundingBox();
  int computeTriIndex(int i, int j, int resolution);
};
} // namespace ASPLINE

#endif
