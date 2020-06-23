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
#ifndef _Little_Geoframe_H__
#define _Little_Geoframe_H__

#include <AdaptiveIsocontouring/comp.h>
#include <Utility/utility.h>

#define GRAD

class LittleGeoframe
{
	public:
		LittleGeoframe();
		~LittleGeoframe();
		void Clear()
		{
			delete [] quads;
			delete [] normals;
			delete [] triangles;
			delete [] verts;
			numtris=0;
			numverts=0;
			numquads=0;
			vsize=10000;
			tsize=10000;
			qsize=10000;
			qsize=10000;
			verts  = (float(*)[3])malloc(sizeof(float[3]) * vsize);
			normals  = (float(*)[3])malloc(sizeof(float[3]) * tsize);
			triangles   = (unsigned int (*)[3])malloc(sizeof(unsigned int[3]) * tsize);
			quads = (unsigned int (*)[4])malloc(sizeof(unsigned int[4]) * qsize);
			quad_render = (bool*) malloc(sizeof(bool) * qsize);
			colors = (float(*)[4])malloc(sizeof(float[4]) * qsize);
		}
		void Clear2()
		{
			numtris = 0;
			numverts = 0;
			numquads=0;
		}
		int addQuad(unsigned int* v, int num);
		int addTri(unsigned int v1, unsigned int v2, unsigned int v3);
		int addVert(float v_pos[3], float norm[3], float c[4]);
		int getNumTris(void)
		{
			return numtris;
		}
		int getNumVerts(void)
		{
			return numverts;
		}
		int getNumQuads(void)
		{
			return numquads;
		}
		void setNumQuads(int i) { numquads = i; }
		void setNumVerts(int i) { numverts = i;}
		void setNumTris(int i) { numtris = i;}
		void setQuadRender(int quad_ind, bool render);
		void calculatenormals();
		void computeNormal(int i0, int i1, int i2, float* out);
		void calculateTriangleNormal(float* norm, unsigned int c);
		void calculateQuadNormals();
		void calculateExtents();
		int getQuad(int i, int j) { return quads[i][j]; }

	private:
		int numverts;
		int numtris;
		int numquads;
		int tsize,vsize,qsize;
		float(*verts)[3];
		float(*normals)[3];
		float(*colors)[4];
		unsigned int (*triangles)[3];
		unsigned int (*quads)[4];
		bool (*quad_render);
		double biggestDim;
		double centerx, centery, centerz;
};

#endif
