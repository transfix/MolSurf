/*
  Copyright 2011 The University of Texas at Austin

	Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <DynamicAdaptiveGrid/MeshSimplex.h>

using namespace DynamicAdaptiveGrid;

int MeshSimplex::meshSimplexIDGenerator = 0;

MeshSimplex::MeshSimplex()
{
	init();
	generateMeshSimplexID();
}

MeshSimplex::MeshSimplex(const MeshSimplex& ms)
{
	numVertices = ms.numVertices;
	degenerate = ms.degenerate;	
	vertices.clear();
	addVertices(ms.vertices);
	generateMeshSimplexID();
}

MeshSimplex::MeshSimplex(vector<MeshVertex *> verts)
{
	init();
	addVertices(verts);
	generateMeshSimplexID();
}

MeshSimplex::~MeshSimplex()
{

}

void MeshSimplex::generateMeshSimplexID()
{
	meshSimplexID = meshSimplexIDGenerator++;
}

void MeshSimplex::init()
{
	numVertices = 0;
	degenerate = true;
	vertices.clear();
}

int MeshSimplex::getID()
{
	return meshSimplexID;
}

int MeshSimplex::getNumVertices()
{
	return numVertices;
}

vector<MeshVertex *> MeshSimplex::getVertices()
{
	return vertices;
}

void MeshSimplex::removeVertices()
{
	for(int i=0;i<numVertices;i++)
	{
		vertices[i]->decSimplices();
	}
	vertices.clear();
}

bool MeshSimplex::isDegenerate()
{
	return degenerate;
}

void MeshSimplex::addVertices(vector<MeshVertex *> verts)
{
	numVertices = verts.size();
	for(int i=0;i<numVertices;i++)
		addVertex(verts[i]);
}

void MeshSimplex::addVertex(MeshVertex * v)
{
	vertices.push_back(v);
	numVertices++;
	v->incSimplices();
}

