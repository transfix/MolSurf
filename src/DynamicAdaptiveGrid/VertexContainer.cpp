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
#include <DynamicAdaptiveGrid/VertexContainer.h>
#include <iostream>
#include <fstream>

using namespace DynamicAdaptiveGrid;

VertexContainer::VertexContainer() {}

VertexContainer::~VertexContainer() {
  vertexCont.clear();
  vertMap.clear();
}

void VertexContainer::addVertex(string str, MeshVertex *mv) {
  vertexCont[str] = mv;
}

void VertexContainer::removeVertex(string str) { vertexCont.erase(str); }

void VertexContainer::updateVertex(string str, MeshVertex *mv) {
  vertexCont.erase(str);
  vertexCont[str] = mv;
}

MeshVertex *VertexContainer::getVertex(string str) {
  if (vertexCont.find(str) == vertexCont.end())
    return NULL;
  return vertexCont[str];
}

void VertexContainer::clearMap() { vertMap.clear(); }

void VertexContainer::clearContainer() { vertexCont.clear(); }

int VertexContainer::getVertexIndex(string str) {
  if (vertMap.find(str) == vertMap.end())
    return -1;
  return vertMap[str];
}

int VertexContainer::size() { return vertexCont.size(); }

void VertexContainer::printVerticesAndUpdateMap(std::ofstream *outputFile) {
  vertMap.clear();

  //      cout<<"Found "<<vertexCont.size()<<" vertices"<<endl;

  map<string, MeshVertex *>::iterator vertStart = vertexCont.begin();
  map<string, MeshVertex *>::iterator vertEnd = vertexCont.end();

  map<string, MeshVertex *>::iterator vertIt;
  int index = 0;

  for (vertIt = vertStart; vertIt != vertEnd; vertIt++) {
    MeshVertex *temp = (*vertIt).second;
    (*outputFile) << temp->getX() << " " << temp->getY() << " " << temp->getZ()
                  << endl;
    vertMap[temp->toString()] = index++;
  }
}
