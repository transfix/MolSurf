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
#include <DynamicAdaptiveGrid/MeshVertex.h>

using namespace DynamicAdaptiveGrid;

int MeshVertex::meshVertexIDGenerator = 0;

MeshVertex::MeshVertex() : Point(0, 0, 0) {
  init();
  initNorm();
  initCol();
  generateMeshVertexID();
}

MeshVertex::MeshVertex(double x, double y, double z) : Point(x, y, z) {
  init();
  initNorm();
  initCol();
  generateMeshVertexID();
}

MeshVertex::MeshVertex(const Point &p) : Point(p.getX(), p.getY(), p.getZ()) {
  init();
  initNorm();
  initCol();
  generateMeshVertexID();
}

MeshVertex::MeshVertex(const MeshVertex &mv)
    : Point(mv.getX(), mv.getY(), mv.getZ()) {
  numSimplices = mv.numSimplices;
  active = mv.active;
  type = mv.type;
  setNormal(mv.normal[0], mv.normal[1], mv.normal[2]);
  setColor(mv.color[0], mv.color[1], mv.color[2], mv.color[3]);
  generateMeshVertexID();
}

MeshVertex::MeshVertex(const Point &p, double *norm)
    : Point(p.getX(), p.getY(), p.getZ()) {
  init();
  setNormal(norm[0], norm[1], norm[2]);
  initCol();
  generateMeshVertexID();
}

MeshVertex::MeshVertex(const Point &p, double *norm, double *col)
    : Point(p.getX(), p.getY(), p.getZ()) {
  init();
  setNormal(norm[0], norm[1], norm[2]);
  setColor(col[0], col[1], col[2], col[3]);
  generateMeshVertexID();
}

MeshVertex::~MeshVertex() {}

void MeshVertex::init() {
  numSimplices = 0;
  type = RAW_VERTEX;
  active = false;
}

void MeshVertex::initNorm() {
  normal[0] = 0;
  normal[1] = 0;
  normal[2] = 0;
}

void MeshVertex::initCol() {
  color[0] = 0;
  color[1] = 0;
  color[2] = 0;
  color[3] = 0;
}

void MeshVertex::generateMeshVertexID() {
  meshVertexID = meshVertexIDGenerator++;
}

double *MeshVertex::getNormal() { return normal; }

double *MeshVertex::Color() { return color; }

int MeshVertex::getType() { return type; }

int MeshVertex::getID() { return ID; }

int MeshVertex::getNumSimplices() { return numSimplices; }

bool MeshVertex::isActive() { return isAct; }

void MeshVertex::setNormal(double nx, double ny, double nz) {
  normal[0] = nx;
  normal[1] = ny;
  normal[2] = nz;

  if (type < RAWN_VERTEX)
    type = RAWN_VERTEX;
}

void MeshVertex::setColor(double r, double g, double b, double alpha) {
  color[0] = r;
  color[1] = g;
  color[2] = b;
  color[3] = alpha;

  if (type >= RAWN_VERTEX)
    type = RAWNC_VERTEX;
}

void MeshVertex::setType(int t) { type = t; }

void MeshVertex::setPosition(double a, double b, double c) {
  x = a;
  y = b;
  z = c;
}

void MeshVertex::incSimplices() {
  numSimplices++;
  active = true;
}

void MeshVertex::decSimplices() {
  if (active)
    numSimplices--;
  if (numSimplices < 1)
    active = false;
}
