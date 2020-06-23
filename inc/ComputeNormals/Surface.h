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
#ifndef SURFACECN_H
#define SURFACECN_H

#include <Utility/utility.h>
#include <ComputeNormals/geometry.h>

using namespace COMPUTENORMALS;

namespace COMPUTENORMALS
{
	class Vertex;
	class Face;
	class Surface
	{
		public:
			Surface();
			void addPoint(Vertex* point);
			void addFacet(Face* facet);
			void Produce_Normal(bool item);
			std::vector<Vertex*> m_Vertices;
			std::vector<Face*> m_Faces;
			int numbpts, numbtris;
			double min[3], max[3];
			bool normalFlag, colorFlag;
		private:
			void TriangleNormal(double* p1, double* p2, double* p3, double* normal);
			void ProduceFaceNormal();
			void Valence_Neighbor();
			void NormalByAveraging();
			void LoopNormal(double* points, int n, double* normal);
			void NormalByLoop();
			double DotProduct(double* v1, double* v2);
			void CrossProduct(double* v1, double* v2, double* v);
			void BoundingBox();
			int computeTriIndex(int i, int j, int resolution);
	};
}

#endif
