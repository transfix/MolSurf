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
#ifndef _HLEVELSET_H_
#define _HLEVELSET_H_

#include <boost/tuple/tuple.hpp>
#include <Geometry/Geometry.h>
#include <vector>

class SimpleVolumeData;
namespace PDBParser {
class GroupOfAtoms;
class Atom;
}; // namespace PDBParser

class HLevelSet {
public:
  HLevelSet();
  virtual ~HLevelSet();
  SimpleVolumeData *
  getHigherOrderLevelSetSurface(PDBParser::GroupOfAtoms *molecule,
                                unsigned int *dim, int nIterations,
                                double weight);
  bool computeFunction_ajrkN(float *func_h, float *func_phi, unsigned int *dim,
                             float *minExtent, float *maxExtent);
  bool getIntersectionPoints(const std::vector<float> vertex_Positions,
                             std::vector<float> vertex_PointColors,
                             int *newobject, int numobject, unsigned int *dim);
  bool computeGeometryNormal(Geometry *geo, unsigned int *dim);
  bool computeGeometryError(Geometry *geo, unsigned int *dim);
  void computeGeometryDerivative(Geometry *geo, float *coeff,
                                 unsigned int *dim);
  bool fixOrientation(Geometry *geo);
  float *coeff;
  bool getAtomListAndExtent(PDBParser::GroupOfAtoms *molecule,
                            std::vector<PDBParser::Atom *> &atomList,
                            float *minExtent, float *maxExtent);
  bool computeFunction(const std::vector<float> Vertices,
                       const std::vector<float> Radius, float *func_phi,
                       unsigned int *dim, float *minExtent, float *maxExtent,
                       int nIterations, double weight);
  bool computeFunction_Zhang(std::vector<float> vertex_Positions,
                             float *funcvalue, unsigned int *dim, float *minExt,
                             float *maxExt, float *weightc);
  bool computeFunction_Zhang_sdf(float *funcvalue, float *weightc);
  bool computeFunction_Xu_Li(float *vertJPosition, int size, float *gridvalue,
                             float *funcvalue, unsigned int *dim, float *minExt,
                             float *maxExt, int *newobject, int numobject);
  void getBoundingBox(const std::vector<float> Vertices,
                      const std::vector<float> Radius, float min[3],
                      float max[3], int probeRadius, double blobbiness);
  void EvaluateCubicSplineOrder1PartialsAtAnyGivenPoint(float *, int, int, int,
                                                        int, int, int, float *,
                                                        float *);
  void ComputeTensorForGivenPoint(int, int, int, int, int, int, float *,
                                  float *, float *, float *, float *);
  void getLowerCorner(float *p, int &u, int &v, int &w);

private:
  unsigned int Dim[3];
  int End;
  float *Funcvalue;
  float *Funcvalue_bak, *Dfang;
  float Dxyz[3], MinExt[3], MaxExt[3], Dt;
  float minExtent[3], maxExtent[3];
  float dxyz[3];
  float *boundarysign;
};

#endif
