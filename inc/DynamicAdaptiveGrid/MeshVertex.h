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
#ifndef DAG_MESH_VERTEX_H
#define DAG_MESH_VERTEX_H

#include <DPG/PG.h>

namespace DynamicAdaptiveGrid {
enum {
  RAW_VERTEX,
  RAWN_VERTEX,
  RAWNC_VERTEX,
};

class MeshVertex : public DPG::Point {
private:
  static int meshVertexIDGenerator;
  int meshVertexID;
  void generateMeshVertexID();
  void init();
  void initNorm();
  void initCol();

protected:
  int type;
  double normal[3];
  double color[4];
  int numSimplices; // numSimplices = the number of simplices that share this
                    // vertex
  bool active;      // active = true when numSimplices > 0
  int index;
  // arand- created this stuff but didn't implement
  // I just want things to build.
  int ID;
  bool isAct;

public:
  MeshVertex();
  MeshVertex(double x, double y, double z);
  MeshVertex(const Point &p);
  MeshVertex(const MeshVertex &mv);
  MeshVertex(const Point &p, double *norm);
  MeshVertex(const Point &p, double *norm, double *col);
  ~MeshVertex();

  double *getNormal();
  double *Color();
  int getType();
  int getID();
  int getNumSimplices();
  bool isActive();

  void setNormal(double nx, double ny, double nz);
  void setColor(double r, double g, double b, double alpha);
  void setType(int t);
  void setPosition(double a, double b, double c);
  void incSimplices();
  void decSimplices();
};

}; // namespace DynamicAdaptiveGrid

#endif
