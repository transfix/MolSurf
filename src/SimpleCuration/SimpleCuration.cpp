/*
  Copyright 2011 The University of Texas at Austin

	Authors: Alex Rand <arand@ices.utexas.edu>
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
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>


#include <GeometryFileTypes/GeometryLoader.h>
#include <Geometry/Geometry.h>

using namespace std;


//geometry->AllocateTris(numverts, triangleIndices.size() / 3);
Geometry * RemoveInteriorPockets(Geometry * inGeom) {

  //cout << "Begin Simple Curation" << endl;

  int nVert = inGeom->m_NumTriVerts;
  int nTri = inGeom->m_NumTris;

  vector<set<int> > nbs;

  nbs.resize(nVert);
  for (int i=0; i<nTri; i++) {
    int a = inGeom->m_Tris[3*i];
    int b = inGeom->m_Tris[3*i+1];
    int c = inGeom->m_Tris[3*i+2];

    nbs[a].insert(b);
    nbs[a].insert(c);
    nbs[b].insert(c);
    nbs[b].insert(a);
    nbs[c].insert(a);
    nbs[c].insert(b);

  }
  

  vector<int> compSizes;

  int * comp = new int[nVert];

  for (int i=0; i<nVert; i++) {
    comp[i] = -1;
  }
  
  int start = 0;

  int iComp=0;
  
  while (start < nVert) {
    if (comp[start] == -1) {
      // do the bfs
      set<int> toDo;
      set<int> found;
      
      toDo.insert(start);
      found.insert(start);
      
      while (toDo.size() > 0) {
	set<int>::iterator it = toDo.begin();
	int next = *it;
	toDo.erase(it);

	comp[next] = iComp;
	
	for (it = nbs[next].begin(); it!= nbs[next].end(); it++) {
	  int another = *it;
	  if (found.find(another) == found.end()) {
	    toDo.insert(another);
	    found.insert(another);
	  }
	}

      }

      //cout << "Component " << iComp << " has size " << found.size() << endl;

      compSizes.push_back(found.size());

      iComp++;
    } 
    start++;
  }
 
  // identify the largest component...
  int maxComp = 0;
  for (int i=0; i<compSizes.size(); i++) {
    if (compSizes[i] > compSizes[maxComp]) {
      maxComp = i;
    }
  }
  
  //cout << "Large component: " << maxComp << endl;

  // create the new geometry
  int mVert = 0;
  for (int i=0; i<nVert; i++) {
    if (comp[i] == maxComp) {
      comp[i] = mVert;
      mVert++;
    } else {
      comp[i] = -1;
    }
  }

  int mTri = 0;
  for (int i=0; i<nTri; i++) {
    int a = inGeom->m_Tris[3*i];
    if (comp[a] >= 0) {
      mTri++;
    }
  }

  Geometry* outGeom = new Geometry;
  outGeom->AllocateTris(mVert, mTri);

  if (inGeom->m_TriVertColorsTransparent) {
    outGeom->AllocateTriVertColors();
  }

  int j=0;
  for (int i=0; i<nVert; i++) {  
    if (comp[i] >= 0) {

      outGeom->m_TriVerts[3*j+0] = inGeom->m_TriVerts[3*i+0];
      outGeom->m_TriVerts[3*j+1] = inGeom->m_TriVerts[3*i+1];
      outGeom->m_TriVerts[3*j+2] = inGeom->m_TriVerts[3*i+2];
      // do the normals too if available...
      if (inGeom->m_TriVertNormals) {
	outGeom->m_TriVertNormals[3*j+0] = inGeom->m_TriVertNormals[3*i+0];
	outGeom->m_TriVertNormals[3*j+1] = inGeom->m_TriVertNormals[3*i+1];
	outGeom->m_TriVertNormals[3*j+2] = inGeom->m_TriVertNormals[3*i+2];
      }
      if (inGeom->m_TriVertColorsTransparent) {
	outGeom->m_TriVertColorsTransparent[4*j+0] = inGeom->m_TriVertColorsTransparent[4*i+0];
	outGeom->m_TriVertColorsTransparent[4*j+1] = inGeom->m_TriVertColorsTransparent[4*i+1];
	outGeom->m_TriVertColorsTransparent[4*j+2] = inGeom->m_TriVertColorsTransparent[4*i+2];
	outGeom->m_TriVertColorsTransparent[4*j+3] = inGeom->m_TriVertColorsTransparent[4*i+3];
      }

      j++;
    }
  }


  j = 0;
  // copy the triangles...
  for (int i=0; i<nTri; i++) {
    int a = inGeom->m_Tris[3*i];
    int b = inGeom->m_Tris[3*i+1];
    int c = inGeom->m_Tris[3*i+2];

    if (comp[a] >= 0) {
      outGeom->m_Tris[3*j] = comp[a];
      outGeom->m_Tris[3*j+1] = comp[b];
      outGeom->m_Tris[3*j+2] = comp[c];
      j++;
    }
  }

  delete [] comp;
  return outGeom;

}
