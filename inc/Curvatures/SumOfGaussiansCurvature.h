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
#ifndef CCV_SUM_OF_GAUSSIANS_CURVATURE_H
#define CCV_SUM_OF_GAUSSIANS_CURVATURE_H

#include <Curvatures/Curvature.h>
#include <Curvatures/CurvaturesGridVoxel.h>

class CurvaturesGridVoxel;

class SumOfGaussiansCurvature : public Curvature {
public:
  SumOfGaussiansCurvature(int numberOfGaussians, double *gaussianCenters,
                          int numberOfGridDivisions, double maxFunctionError,
                          double blobbiness, int numberOfPoints, float *points,
                          double *HandK, double *normals, double *k1Vector,
                          double *k2Vector);
  virtual ~SumOfGaussiansCurvature();
  bool initialize();

protected:
  void getMinMax();
  void getIndices(double x, double y, double z, int *gridx, int *gridy,
                  int *gridz);
  void populateGrid(double overlapExtent);
  void createGrid();
  inline void evalCurvature(double *phiX, double *phiY, double *phiZ,
                            double *phiXX, double *phiYY, double *phiZZ,
                            double *phiXY, double *phiXZ, double *phiYZ,
                            double x, double y, double z);
  inline int getIndex(double x, double y, double z);
  inline void evalCurvatureDueToOneAtom(
      double *phiX, double *phiY, double *phiZ, double *phiXX, double *phiYY,
      double *phiZZ, double *phiXY, double *phiXZ, double *phiYZ, double xc,
      double yc, double zc, double rad, double x, double y, double z);
  double m_Min[3], m_Max[3];
  CurvaturesGridVoxel *m_CurvaturesGridVoxels;
  int m_NumberOfGaussians;
  double *m_GaussianCenters;
  int m_NumberOfGridDivisions;
  double m_MaxFunctionError;
  double m_Blobbiness;
};

#endif
