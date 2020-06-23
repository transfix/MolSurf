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
// Craig: moved to LittleGeoframe to avoid collision with LBIE's larger Geoframe class

#include <AdaptiveIsocontouring/Geoframe.h>

extern void normalize(float* vect, int size);

LittleGeoframe::LittleGeoframe()
{
	normals = 0;
	numverts = 0;
	numtris = 0;
	numquads=0;
	vsize=10000;
	tsize=10000;
	qsize=10000;
	qsize=10000;
	verts  = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	normals  = (float(*)[3])malloc(sizeof(float[3]) * tsize);
	triangles   = (unsigned int (*)[3])malloc(sizeof(unsigned int[3]) * tsize);
	quads = (unsigned int (*)[4])malloc(sizeof(unsigned int[4]) * qsize);
	quad_render = (bool*)malloc(sizeof(bool) * qsize);
	colors = (float(*)[4])malloc(sizeof(float[4]) * qsize);
}

LittleGeoframe::~LittleGeoframe()
{
	delete [] verts;
	delete [] quads;
	delete [] quad_render;
	delete [] triangles;
	delete [] normals;
}

void cross(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void LittleGeoframe::calculateTriangleNormal(float* norm, unsigned int c)
{
	float v1[3], v2[3];
	int vert;
	vert = triangles[c][0];
	v1[0] = v2[0] = -verts[vert][0];
	v1[1] = v2[1] = -verts[vert][1];
	v1[2] = v2[2] = -verts[vert][2];
	vert = triangles[c][1];
	v1[0] += verts[vert][0];
	v1[1] += verts[vert][1];
	v1[2] += verts[vert][2];
	vert = triangles[c][2];
	v2[0] += verts[vert][0];
	v2[1] += verts[vert][1];
	v2[2] += verts[vert][2];
	cross(norm, v1, v2);
}

void LittleGeoframe::computeNormal(int i0, int i1, int i2, float* norm)
{
	float* v0 = verts[i0];
	float* v1 = verts[i1];
	float* v2 = verts[i2];
	float t0[3];
	float t1[3];
	for(int j=0; j<3; j++)
	{
		t0[j] = verts[i0][j] - verts[i1][j];
		t1[j] = verts[i2][j] - verts[i1][j];
	}
	float ret[3];
	cross(ret, t1, t0);
	copyVect(ret, norm, 3);
}

void LittleGeoframe::calculateQuadNormals()
{
	for(int j=0; j<numquads; j++)
	{
		if(quads[j][0] == quads[j][1] || quads[j][0] == quads[j][2])
		{
			computeNormal(quads[j][1], quads[j][2], quads[j][3], normals[j]);
			cout<<"done 1"<<endl;
		}
		else if(quads[j][1] == quads[j][2])
		{
			computeNormal(quads[j][0], quads[j][2], quads[j][3], normals[j]);
			cout<<"done 2"<<endl;
		}
		else
		{
			computeNormal(quads[j][0], quads[j][1], quads[j][2], normals[j]);
		}
	}
}

int LittleGeoframe::addQuad(unsigned int* v , int num)
{
	assert(num==3 || num==4);
	if(numquads >= qsize)
	{
		qsize<<=1;
		quads = (unsigned int (*)[4])realloc(quads, sizeof(unsigned int[4]) * qsize);
		quad_render = (bool*)realloc(quad_render, sizeof(bool) * qsize);
		normals = (float(*)[3])realloc(normals,sizeof(float[3])*qsize);
	}
	if(num==4)
	{
		quads[numquads][0] = v[0];
		quads[numquads][1] = v[1];
		quads[numquads][2] = v[2];
		quads[numquads][3] = v[3];
		return numquads++;
	}
	else if(num==3)
	{
		triangles[numtris][0] = v[0];
		triangles[numtris][1] = v[1];
		triangles[numtris][2] = v[2];
		return numtris++;
	}
}

int LittleGeoframe::addTri(unsigned int v1,unsigned int v2,unsigned int v3)
{
	if(numtris+1 >= tsize)
	{
		tsize<<=1;
		triangles = (unsigned int (*)[3])realloc(triangles, sizeof(unsigned int[3]) * tsize);
	}
	triangles[numtris][0] = v1;
	triangles[numtris][1] = v2;
	triangles[numtris][2] = v3;
	return numtris++;
}

int LittleGeoframe::addVert(float v_pos[3], float norm[3], float c[4])
{
	int i;
	if(numverts+1 > vsize)
	{
		vsize<<=1;
		verts = (float(*)[3])realloc(verts,sizeof(float[3])*vsize);
		colors = (float(*)[4])realloc(colors,sizeof(float[4])*vsize);
	}
	for(i=0; i<3; i++)
	{
		verts[numverts][i]=v_pos[i];
	}
	for(i=0; i<4; i++)
	{
		colors[numverts][i] = c[i];
	}
	return numverts++;
}


void LittleGeoframe::calculatenormals()
{
	float normal[3];
	float len;
	// for each triangle
	for(int c=0; c<numtris; c++)
	{
		calculateTriangleNormal(normal, c);
		normals[c][0] = normal[0];
		normals[c][1] = normal[1];
		normals[c][2] = normal[2];
	}
	// normalize the vectors
	for(int vert=0; vert<numtris; vert++)
	{
		len = (float) sqrt(normals[vert][0] * normals[vert][0] +
						   normals[vert][1] * normals[vert][1] +
						   normals[vert][2] * normals[vert][2]);
		normals[vert][0]/=len;
		normals[vert][1]/=len;
		normals[vert][2]/=len;
	}
}

void LittleGeoframe::calculateExtents()
{
	float max_x, min_x;
	float max_y, min_y;
	float max_z, min_z;
	float value;
	for(int c=0; c<numverts; c++)
	{
		if(c==0)
		{
			max_x = min_x = verts[c][0];
			max_y = min_y = verts[c][1];
			max_z = min_z = verts[c][2];
		}
		else
		{
			value = verts[c][0];
			max_x = (value>max_x?value:max_x);
			min_x = (value<min_x?value:min_x);
			value = verts[c][1];
			max_y = (value>max_y?value:max_y);
			min_y = (value<min_y?value:min_y);
			value = verts[c][2];
			max_z = (value>max_z?value:max_z);
			min_z = (value<min_z?value:min_z);
		}
	}
	biggestDim = (max_y-min_y>max_x-min_x?max_y-min_y:max_x-min_x);
	biggestDim = (max_z-min_z>biggestDim?max_z-min_z:biggestDim);
	centerx = (max_x+min_x)/2.0;
	centery = (max_y+min_y)/2.0;
	centerz = (max_z+min_z)/2.0;
}


void LittleGeoframe::setQuadRender(int quad_ind, bool render)
{
	quad_render[quad_ind] = render;
}

