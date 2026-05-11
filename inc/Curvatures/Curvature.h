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
#ifndef CCV_CURVATURE_H
#define CCV_CURVATURE_H

#include <Utility/utility.h>
#include <UsefulMath/LinearAlgebra.h>

class Curvature {
public:
  Curvature(int numberOfPoints, float *points, double *HandK, double *normals,
            double *k1Vector, double *k2Vector);
  virtual ~Curvature();
  bool getCurvatures();
  bool write(const char *filename);
  static bool read(const char *filename, int *numberOfPoints, double **HandK,
                   double **normals, double **k1Vector, double **k2Vector);

protected:
  virtual void evalCurvature(double *phiX, double *phiY, double *phiZ,
                             double *phiXX, double *phiYY, double *phiZZ,
                             double *phiXY, double *phiXZ, double *phiYZ,
                             double x, double y, double z) = 0;
  static bool read2Values(FILE *fp, int numberOfPoints, double **values);
  static bool read3Values(FILE *fp, int numberOfPoints, double **values);
  int m_NumberOfPoints;
  float *m_Points;
  double *m_HandK;
  double *m_Normals;
  double *m_K1Vector;
  double *m_K2Vector;
  bool m_Initialized;
};

#endif
