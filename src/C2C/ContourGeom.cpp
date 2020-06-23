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
// ContourGeom - Class for a contour surface
// Copyright (c) 1998 Emilio Camahort, Dan Schikore

#include <C2C/ContourGeom.h>

// ContourGeom() - basic constructor
ContourGeom::ContourGeom()
{
	nvert = 0;
	ntri = 0;
	vsize =  500;
	tsize = 1000;
	vert  = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	vnorm = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	vcol  = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	tri   = (u_int(*)[3])malloc(sizeof(u_int[3]) * tsize);
	pvset = new set<VPosition, LtPos>();
	edgeset = new set<EdgeIndex, LtEdge>();
}

// copy constructor: It copies the vert, norm and triangle array
ContourGeom::ContourGeom(const ContourGeom& con3d)
{
	nvert = con3d.nvert;
	ntri = con3d.ntri;
	vsize = con3d.vsize;
	tsize = con3d.tsize;
	vert = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	vnorm = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	vcol = (float(*)[3])malloc(sizeof(float[3]) * vsize);
	tri   = (u_int(*)[3])malloc(sizeof(u_int[3]) * tsize);
	memcpy(vert, con3d.vert, sizeof(float[3])*vsize);
	memcpy(vnorm, con3d.vnorm, sizeof(float[3])*vsize);
	memcpy(vcol, con3d.vcol, sizeof(float[3])*vsize);
	memcpy(tri, con3d.tri, sizeof(u_int[3])*tsize);
	pvset = new set<VPosition, LtPos>();
	edgeset = new set<EdgeIndex, LtEdge>();
}

// ~ContourGeom() - free allocated memory
ContourGeom::~ContourGeom()
{
	free(vert);
	free(vnorm);
	free(tri);
	free(vcol);
	delete pvset;
	delete edgeset;
}

// AddVert() - add a vertex with the given (unit) normal
int
ContourGeom::AddVert(float x,float y,float z, float nx,float ny,float nz, float cr,float cg,float cb)
{
	int n = nvert++;
	if(nvert > vsize)
	{
		vsize<<=1;
		vert  = (float(*)[3])realloc(vert, sizeof(float[3]) * vsize);
		vnorm = (float(*)[3])realloc(vnorm, sizeof(float[3]) * vsize);
		vcol = (float(*)[3])realloc(vcol, sizeof(float[3]) * vsize);
	}
	vert[n][0] = x;
	vert[n][1] = y;
	vert[n][2] = z;
	vnorm[n][0] = nx;
	vnorm[n][1] = ny;
	vnorm[n][2] = nz;
	vcol[n][0] = cr;
	vcol[n][1] = cg;
	vcol[n][2] = cb;
	return(n);
}

// AddVertUnique() - add a vertex with the given (unit) normal
int ContourGeom::AddVertUnique(float x, float y, float z, float nx, float ny, float nz, float cr, float cg, float cb)
{
	VPosition vtx;
	vtx.idx = nvert;
	vtx.cord[0] = x;
	vtx.cord[1] = y;
	vtx.cord[2] = z;
	set<VPosition, LtPos>::iterator it = pvset->find(vtx);
	if(it != pvset->end())
	{
		return it->idx;
	}
	pvset->insert(vtx);
	return(AddVert(x,y,z, nx,ny,nz, cr,cg,cb));
}

int ContourGeom::AddVertUnique(float p[3], float n[3], float c[3], EdgeIndex ei)
{
	ei.idx = nvert;
	set<EdgeIndex, LtEdge>::iterator it = edgeset->find(ei);
	if(it != edgeset->end())
	{
		return it->idx;
	}
	edgeset->insert(ei);
	return (AddVert(p, n, c));
}
// AddTri() - add a triangle indexed by it's 3 vertices
int ContourGeom::AddTri(u_int v1, u_int v2, u_int v3)
{
	int n = ntri++;
	if(ntri > tsize)
	{
		tsize<<=1;
		tri = (u_int(*)[3])realloc(tri, sizeof(u_int[3]) * tsize);
	}
	tri[n][0] = v1;
	tri[n][1] = v2;
	tri[n][2] = v3;
	return(n);
}

void ContourGeom::merge(ContourGeom* con)
{
	int _ntri = ntri;
	int _nvert = nvert;
	ntri += con->ntri;
	nvert += con->nvert;
	if(nvert > vsize)
	{
		vsize = nvert * 2;
		vert = (float(*)[3])realloc(vert, sizeof(float[3]) * vsize);
		vnorm = (float(*)[3])realloc(vnorm, sizeof(float[3]) * vsize);
		vcol = (float(*)[3])realloc(vcol, sizeof(float[3]) * vsize);
	}
	memcpy(vert + _nvert, con->vert, sizeof(float[3])*con->nvert);
	memcpy(vnorm + _nvert, con->vnorm, sizeof(float[3])*con->nvert);
	memcpy(vcol + _nvert, con->vcol, sizeof(float[3])*con->nvert);
	if(ntri > tsize)
	{
		tsize = ntri * 2;
		tri = (unsigned int (*)[3])realloc(tri, sizeof(int[3]) * tsize);
	}
	memcpy(tri + _ntri, con->tri, sizeof(int[3])*con->ntri);
	// adjust indices of added triangles
	for(int i = _ntri; i < ntri; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			tri[i][j] += _nvert;
		}
	}
}
