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
#ifndef DAG_MESH_SIMPLEX_H
#define DAG_MESH_SIMPLEX_H

#include <DPG/PG.h>
#include <DynamicAdaptiveGrid/MeshVertex.h>

namespace DynamicAdaptiveGrid {
class MeshSimplex {
private:
  static int meshSimplexIDGenerator;
  int meshSimplexID;
  void generateMeshSimplexID();
  void init();

protected:
  vector<MeshVertex *> vertices; // must be in counterclockwise order. vertices
                                 // cannot be removed.
  int numVertices; // numVertices = the number of vertices of the simplex
  bool degenerate; // degenerate = true when numVertices < 3. Further checks for
                   // colinear point will be added later

public:
  MeshSimplex();
  MeshSimplex(const MeshSimplex &mv);
  MeshSimplex(vector<MeshVertex *> verts);
  ~MeshSimplex();

  int getID();
  int getNumVertices();
  vector<MeshVertex *> getVertices();
  bool isDegenerate();

  void addVertices(vector<MeshVertex *> verts);
  void addVertex(MeshVertex *v);
  void removeVertices();
};
}; // namespace DynamicAdaptiveGrid

#endif
