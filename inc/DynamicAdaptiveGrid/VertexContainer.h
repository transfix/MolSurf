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
#ifndef DAG_VERTEX_CONTAINER_H
#define DAG_VERTEX_CONTAINER_H

#include <DynamicAdaptiveGrid/MeshVertex.h>
#include <iostream>
#include <string>
#include <map>

using namespace std;

namespace DynamicAdaptiveGrid
{
	class VertexContainer
	{
		protected:
			map<string,MeshVertex*> vertexCont;
			map<string,int> vertMap;
//			VertexMap vertMap;

		public:
			VertexContainer();
			~VertexContainer();			

			void addVertex(string, MeshVertex*);
			void removeVertex(string);
			void updateVertex(string, MeshVertex*);
			MeshVertex* getVertex(string);
			void printVerticesAndUpdateMap(std::ofstream *outputFile);
			int getVertexIndex(string);
			int size();
			void clearMap();
			void clearContainer();
	};
	
};

#endif
