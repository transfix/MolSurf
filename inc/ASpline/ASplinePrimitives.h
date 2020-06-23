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
#ifndef ASPLINE_PRIMITIVES_H
#define ASPLINE_PRIMITIVES_H

#include <Utility/utility.h>

namespace ASPLINE
{
	class Vertex
	{
		public:
			Vertex()
			{
			}
			double   point[3];
			double   normal[3];
			double	 color[3];
			double   weight;
			int	 numFaces;

			// Print contents to a file
			void print(FILE* fp)
			{
			  //fprintf(fp, "%f %f %f\n", point[0], point[1], point[2]);
			  fprintf(fp, "%f %f %f %f %f %f %f\n", point[0], point[1], point[2],
				  normal[0], normal[1], normal[2], weight);
			}

			// Copy constructor
			Vertex(const Vertex& v)
			{
				point[0] = v.point[0];
				point[1] = v.point[1];
				point[2] = v.point[2];
				normal[0] = v.normal[0];
				normal[1] = v.normal[1];
				normal[2] = v.normal[2];
				color[0] = v.color[0];
				color[1] = v.color[1];
				color[2] = v.color[2];
				weight = v.weight;
			}

			//Constructor
			Vertex(double x, double y, double z, double nx, double ny, double nz, double red, double green, double blue, double weight_ = 1)
			{
				point[0] = x;
				point[1] = y;
				point[2] = z;
				normal[0] = nx;
				normal[1] = ny;
				normal[2] = nz;
				color[0] = red;
				color[0] = green;
				color[0] = blue;
				weight = weight_;
			}


	};

	class Face
	{
		public:
			Face()
			{
			}
			int index[3];      		/* the indexs of a triangle-increase order*/
			int indexInRAW[3]; 		/* the indices of a triangle in the origianl raw file */
			std::vector<Vertex*> node;	/* Gaussian integration nodes             */
			int orien;        		/* 1--orientaioned, othwise no            */
			int adjTri[3];    		/* neighbor triangles index               */
			double center[3];		/* the center of the triangle             */
			double normal[3];		/* face normal                            */
			int whichnml[3];		/* index[0] normal index                  */
			int tou;			/* the begin index of the insided points  */
			int wei;			/* the end index+1 of the insided points  */
			int state;			/* 1--removed; 2-- new, 0--old;           */
			std::vector<Face*> subFacets;

			//Print contents to a file
			void print(FILE* fp)
			{
				for(int i=0; i<node.size(); i++)
				{
					node[i]->print(fp);
				}
			}

			//Copy constructor
			Face(const Face& f)
			{
				index[0] = f.index[0];
				index[1] = f.index[1];
				index[2] = f.index[2];
				indexInRAW[0] = f.indexInRAW[0];
				indexInRAW[1] = f.indexInRAW[1];
				indexInRAW[2] = f.indexInRAW[2];
				node = f.node; //rely on stl's copy constructor here
				orien = f.orien;
				adjTri[0] = f.adjTri[0];
				adjTri[1] = f.adjTri[1];
				adjTri[2] = f.adjTri[2];
				whichnml[0] = f.whichnml[0];
				whichnml[1] = f.whichnml[1];
				whichnml[2] = f.whichnml[2];
				center[0] = f.center[0];
				center[1] = f.center[1];
				center[2] = f.center[2];
				tou = f.tou;
				wei = f.wei;
				state = f.state;
				subFacets = f.subFacets; //rely on stl's copy constructor here
			}


			// Assign indexInRAW to index such that index[0] <= index[1] <= index[2]
			// Inline for now but should be moved to Face.h and Face.cpp eventually
			Face(int index0, int index1, int index2)
			{
				indexInRAW[0] = index0;
				indexInRAW[1] = index1;
				indexInRAW[2] = index2;
				adjTri[0] = -1;
				adjTri[1] = -1;
				adjTri[2] = -1;
				whichnml[0] = -1;
				whichnml[1] = -1;
				whichnml[2] = -1;
				if(indexInRAW[0] < indexInRAW[1] && indexInRAW[1] < indexInRAW[2])
				{
					index[0] = indexInRAW[0];
					index[1] = indexInRAW[1];
					index[2] = indexInRAW[2];
				}
				else if(indexInRAW[0] < indexInRAW[2] && indexInRAW[2] < indexInRAW[1])
				{
					index[0] = indexInRAW[0];
					index[1] = indexInRAW[2];
					index[2] = indexInRAW[1];
				}
				else if(indexInRAW[1] < indexInRAW[2] && indexInRAW[2] < indexInRAW[0])
				{
					index[0] = indexInRAW[1];
					index[1] = indexInRAW[2];
					index[2] = indexInRAW[0];
				}
				else if(indexInRAW[1] < indexInRAW[0] && indexInRAW[0] < indexInRAW[2])
				{
					index[0] = indexInRAW[1];
					index[1] = indexInRAW[0];
					index[2] = indexInRAW[2];
				}
				else if(indexInRAW[2] < indexInRAW[1] && indexInRAW[1] < indexInRAW[0])
				{
					index[0] = indexInRAW[2];
					index[1] = indexInRAW[1];
					index[2] = indexInRAW[0];
				}
				else if(indexInRAW[2] < indexInRAW[0] && indexInRAW[0] < indexInRAW[1])
				{
					index[0] = indexInRAW[2];
					index[1] = indexInRAW[0];
					index[2] = indexInRAW[1];
				}
			}
	};
}
#endif
