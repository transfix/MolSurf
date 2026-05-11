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
// Classes Vertex and Face

#ifndef GEOMETRYCN_H
#define GEOMETRYCN_H

#include <Utility/utility.h>

namespace COMPUTENORMALS {
class Vertex {
public:
  Vertex() {}
  double Point[3];
  double Normal[3];
  int NumFaces;
  double Color[3];
  double weight;
  int facevalence;
  int edgevalence;
  std::vector<int> faceneighbor;
  std::vector<int> edgeneighbor;
};

class Face {
public:
  Face() {}
  int Index[3];      /* the indexs of a triangle-increase order*/
  int IndexInRAW[3]; /* the indices of a triangle in the origianl raw file */
  std::vector<Vertex *> Node; /* Gaussian integration nodes             */
  int Orien;                  /* 1--orientaioned, othwise no            */
  int AdjTri[3];              /* neighbor triangles index               */
  double Center[3];           /* the center of the triangle             */
  double Normal[3];           /* face normal                            */
  int whichnml[3];            /* Index[0] normal index                  */
  int Tou;                    /* the begin index of the insided points  */
  int Wei;                    /* the end index+1 of the insided points  */
  int State;                  /* 1--removed; 2-- new, 0--old;           */
  std::vector<Face *> subFacets;
};
} // namespace COMPUTENORMALS

#endif
