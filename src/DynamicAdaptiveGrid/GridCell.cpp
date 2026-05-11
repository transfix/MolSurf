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
#include <DynamicAdaptiveGrid/GridCell.h>
#include <DynamicAdaptiveGrid/MarchingCubeTables.h>

using namespace DynamicAdaptiveGrid;

int GridCell::gridCellIDGenerator = 0;

GridCell::GridCell() : Point(0, 0, 0) {
  functionValue = 0;
  state = OUT_CELL;
  cellDimension = 0.7;
  initContour();
  generateGridCellID();
}

GridCell::GridCell(const GridCell &gc) : Point(gc.x, gc.y, gc.z) {
  functionValue = gc.functionValue;
  state = gc.state;
  cellDimension = gc.cellDimension;

  initContour();
  setContour(gc);

  gridPoints.clear();
  for (int i = 0; i < gc.gridPoints.size(); i++)
    gridPoints.push_back(gc.gridPoints[i]);

  generateGridCellID();
}

GridCell::GridCell(double a, double b, double c) : Point(a, b, c) {
  functionValue = 0;
  state = OUT_CELL;
  cellDimension = 0.7;

  initContour();
  generateGridCellID();
}

GridCell::GridCell(double a, double b, double c, double F) : Point(a, b, c) {
  functionValue = F;
  state = OUT_CELL;
  cellDimension = 0.7;

  initContour();
  generateGridCellID();
}

GridCell::GridCell(double a, double b, double c, double F, int s)
    : Point(a, b, c) {
  functionValue = F;
  setState(s);
  cellDimension = 0.7;

  initContour();
  generateGridCellID();
}

GridCell::GridCell(double a, double b, double c, double F, int s, double d)
    : Point(a, b, c) {
  functionValue = F;
  setState(s);
  cellDimension = d;

  initContour();
  generateGridCellID();
}

GridCell::GridCell(double a, double b, double c, double F, int s, double d,
                   vector<GridPoint *> gps) {
  functionValue = F;
  setState(s);
  cellDimension = d;

  initContour();

  setGridPoints(gps);

  generateGridCellID();
}

GridCell::~GridCell() {
  gridPoints.clear();

  if (contourState == STABLE_CONTOUR)
    destroyContour();

  simplices.clear();
}

void GridCell::generateGridCellID() { gridCellID = gridCellIDGenerator++; }

int GridCell::getState() { return state; }

int GridCell::getContourState() { return contourState; }

double GridCell::getFunctionValue() { return functionValue; }

int GridCell::getGridCellID() { return gridCellID; }

double GridCell::getCellDimension() { return cellDimension; }

vector<GridPoint *> GridCell::getGridPoints() { return gridPoints; }

GridPoint *GridCell::getGridPoint(int i) {
  if (i < gridPoints.size())
    return gridPoints[i];
  return NULL;
}

int GridCell::getNumGridPoints() { return gridPoints.size(); }

int GridCell::getNumSimplices() { return simplices.size(); }

void GridCell::setState(int s) {
  if (s < UNDEFINED_CELL)
    state = s;
  else
    state = UNDEFINED_CELL;
}

void GridCell::setFunctionValue(double f) { functionValue = f; }

void GridCell::setCellDimension(double d) { cellDimension = d; }

void GridCell::setGridPoints(vector<GridPoint *> gps) {
  int size = gps.size();

  gridPoints.clear();

  for (int i = 0; i < size; i++) {
    gridPoints.push_back(gps[i]);
  }

  for (int iEdge = 0; iEdge < 12; iEdge++) {
    GridPoint *a = gridPoints[a2iEdgeConnection[iEdge][0]];
    GridPoint *b = gridPoints[a2iEdgeConnection[iEdge][1]];

    edgeNames[iEdge] = a->toString() + b->toString();
  }
}

void GridCell::addGridPoint(
    GridPoint *gp) // adds a gridpoint maintaining the order of vertices used in
                   // the marchingcubes
{
  //      cout<<gp->getX()<<" "<<gp->getY()<<" "<<gp->getZ()<<"   ";

  bool xGreater = gp->getX() > x;
  bool yGreater = gp->getY() > y;
  bool zGreater = gp->getZ() > z;

  int flag = 0;

  if (xGreater)
    flag += 1;
  if (yGreater)
    flag += 2;
  if (zGreater)
    flag += 4;

  int index;

  switch (flag) {
  case 0:
    index = 0;
    break;
  case 1:
    index = 1;
    break;
  case 2:
    index = 3;
    break;
  case 3:
    index = 2;
    break;
  case 4:
    index = 4;
    break;
  case 5:
    index = 5;
    break;
  case 6:
    index = 7;
    break;
  case 7:
    index = 6;
    break;
  default:
    index = 0;
    break;
  }

  //      cout<<index<<endl;

  gridPoints[index] = gp;
}

void GridCell::setPosition(double a, double b, double c) {
  x = a;
  y = b;
  z = c;
}

void GridCell::markForUpdate() { contourState = MARKED_FOR_UPDATE_CONTOUR; }

void GridCell::markForCreate() { contourState = MARKED_FOR_CREATE_CONTOUR; }

void GridCell::markForDelete() { contourState = MARKED_FOR_DELETE_CONTOUR; }

int GridCell::classifyCell() {
  bool vdwFound = false;
  bool sasFound = false;
  bool outFound = false;

  //      cout<<state<<" ";

  for (int i = 0; i < gridPoints.size(); i++) {
    int st = gridPoints[i]->getState();

    if (st == VDW_POINT)
      vdwFound = true;
    else if (st == SAS_POINT)
      sasFound = true;
    else if (st == OUT_POINT)
      outFound = true;

    //              cout<<st<<" ";
  }

  //      cout<<endl;

  if (vdwFound && !sasFound && !outFound)
    state = BURRIED_CELL;

  else if (vdwFound && sasFound && !outFound)
    state = VDW_CELL;

  else if (!vdwFound && sasFound && !outFound)
    state = BAND_CELL;

  else if (!vdwFound && sasFound && outFound)
    state = SAS_CELL;

  else if (!vdwFound && !sasFound && outFound)
    state = OUT_CELL;

  else {
    state = UNDEFINED_CELL;
    //              cout<<vdwFound<<" "<<sasFound<<" "<<outFound<<endl;
    //              cout<<x<< " " << y << " " << z <<endl;
  }

  //      cout<<state<<endl;

  return state;
}

int GridCell::classifyCell(vector<GridPoint *> gps) {
  gridPoints.clear();

  for (int i = 0; i < gps.size(); i++)
    gridPoints.push_back((GridPoint *)gps[i]);

  return classifyCell();
}

bool GridCell::computeContour() { return computeContour(0.0); }

bool GridCell::computeContour(double isoValue) {
  vmarchcube(isoValue, false);
  contourState = STABLE_CONTOUR;
  return true;
}

void GridCell::initContour() {
  numSimplices = 0;
  numVertices = 0;
  contourState = NO_CONTOUR;
}

void GridCell::setContour(
    const GridCell &gc) // we are not going to create new simplices. simply copy
                        // the pointers
{
  if (numSimplices != 0)
    destroyContour();

  for (int i = 0; i < 12; i++)
    vertices[i] = gc.vertices[i];

  numSimplices = gc.simplices.size();

  for (int i = 0; i < numSimplices; i++)
    simplices.push_back(gc.simplices[i]);

  contourState = STABLE_CONTOUR;
}

void GridCell::destroyContour() {
  for (int i = 0; i < simplices.size(); i++) {
    simplices[i]->removeVertices();
  }
  simplices.clear();

  for (int i = 0; i < numVertices;
       i++) // Vertices can be shared with other cells. So though it is safe to
            // clear the vector,
  { // but it is safe to delete vertices only if they are no longer active (not
    // shared with other cells)
    MeshVertex *v = vertices[i];

    if (v) {
      if (!v->isActive()) {
        vertCon->removeVertex(edgeNames[i]);
      }
      delete v;
    }

    vertices[i] = NULL;
  }

  initContour();
}

void GridCell::setVertexContainer(VertexContainer *vc) { vertCon = vc; }

float GridCell::fGetOffset(const float &fValue1, const float &fValue2,
                           const float &fValueDesired) {
  double fDelta = fValue2 - fValue1;
  if (fDelta == 0.0) {
    return 0.5;
  }
  return (fValueDesired - fValue1) / fDelta;
}

void GridCell::vmarchcube(
    float iso, bool updateMode) // if updateMode is true, then existing
                                // meshVertices will be overwritten
{
  int iCorner, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
  float fOffset;

  // Find which vertices are inside of the surface and which are outside

  int numVertices = gridPoints.size() < 8 ? gridPoints.size() : 8;

  int iVertexFlags = 0; // bitwise flags to mark vertices which are inside (1),
                        // or outside (0) the iso-surface

  for (int i = 0; i < numVertices; i++) {
    //              cout<<gridPoints[i]->getFunctionValue()<<" ";
    if (gridPoints[i]->getFunctionValue() <= iso)
      iVertexFlags |= 1 << i;
  }
  //      cout<<endl;

  // Find which edges are intersected by the surface

  iEdgeFlags = aiCubeEdgeFlags[iVertexFlags]; // Another set of bitwise flags to
                                              // mark edges which are
                                              // intersected by the iso-surface

  // If the cube is entirely inside or outside of the surface, then there will
  // be no intersections

  if (iEdgeFlags == 0) {
    //              cout<<"Not intersected"<<endl;
    return;
  }

  // Find the point of intersection of the surface with each edge
  // Then find the normal to the surface at those points

  for (iEdge = 0; iEdge < 12; iEdge++) {
    if (iEdgeFlags & (1 << iEdge)) {
      GridPoint *a = gridPoints[a2iEdgeConnection[iEdge][0]];
      GridPoint *b = gridPoints[a2iEdgeConnection[iEdge][1]];

      fOffset = fGetOffset(a->getFunctionValue(), b->getFunctionValue(), iso);

      string edgeName = edgeNames[iEdge];

      //                      cout<<a2iEdgeConnection[iEdge][0]<<"
      //                      "<<a2iEdgeConnection[iEdge][1]<<endl;
      //                      cout<<edgeName<<endl;

      if (vertCon->getVertex(edgeName) == NULL || updateMode) {
        // add the vert...
        //                              cout<<"Adding new vertex"<<endl;
        double xTemp = (1 - fOffset) * a->getX() + fOffset * b->getX();
        double yTemp = (1 - fOffset) * a->getY() + fOffset * b->getY();
        double zTemp = (1 - fOffset) * a->getZ() + fOffset * b->getZ();

        if (!updateMode) {
          MeshVertex *mv = new MeshVertex(xTemp, yTemp, zTemp);
          vertCon->addVertex(edgeName, mv);
        } else {
          vertCon->getVertex(edgeName)->setPosition(xTemp, yTemp, zTemp);
        }
        numVertices++;
      }

      vertices[iEdge] = vertCon->getVertex(edgeName);
    }
  }
  //      cout<<"Added vertexes"<<endl;

  // Add the triangles that were found.  There can be up to five per cube
  for (iTriangle = 0; iTriangle < 5; iTriangle++) {
    if (a2iTriangleConnectionTable[iVertexFlags][3 * iTriangle] < 0)
      break;

    MeshSimplex *ms = new MeshSimplex();
    for (iCorner = 0; iCorner < 3; iCorner++) {
      int iVertex =
          a2iTriangleConnectionTable[iVertexFlags][3 * iTriangle + iCorner];
      ms->addVertex(vertices[iVertex]);
      numSimplices++;
      //                      cout<<"Adding triangle"<<endl;
    }
    simplices.push_back(ms);
  }
  //      cout<<"Added triangles"<<endl;
}

void GridCell::printMesh(std::ofstream *outputFile) {
  int numSimplices = simplices.size();
  int numverts;

  for (int i = 0; i < numSimplices; i++) {
    vector<MeshVertex *> verts = simplices[i]->getVertices();

    numverts = verts.size();

    if (numverts != 3)
      cout << "Error: bad simplex" << endl;

    for (int j = 0; j < numverts; j++) {
      (*outputFile) << vertCon->getVertexIndex(verts[j]->toString()) << " ";
    }
    (*outputFile) << endl;
  }
}
