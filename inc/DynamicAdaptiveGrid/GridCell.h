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
#ifndef DAG_GRIDCELL_H
#define DAG_GRIDCELL_H

#include <DPG/PG.h>
#include <DynamicAdaptiveGrid/GridPoint.h>
#include <DynamicAdaptiveGrid/MeshSimplex.h>
#include <DynamicAdaptiveGrid/VertexContainer.h>
#include <iostream>
#include <fstream>

namespace DynamicAdaptiveGrid {
enum { BURRIED_CELL, VDW_CELL, BAND_CELL, SAS_CELL, OUT_CELL, UNDEFINED_CELL };

enum {
  STABLE_CONTOUR,
  MARKED_FOR_DELETE_CONTOUR,
  MARKED_FOR_CREATE_CONTOUR,
  MARKED_FOR_UPDATE_CONTOUR,
  NO_CONTOUR
};

class GridPoint;

class GridCell : public DPG::Point {
private:
  static int gridCellIDGenerator;
  int gridCellID;
  void generateGridCellID();

  void initContour();
  void setContour(const GridCell &gc);

  void sortGridPoints();
  float fGetOffset(const float &fValue1, const float &fValue2,
                   const float &fValueDesired);
  void vmarchcube(float iso, bool updateMode);

protected:
  double functionValue;
  int state;
  int contourState;
  double cellDimension;
  vector<GridPoint *>
      gridPoints; // given the cell dimension and position, it is easy to find
                  // the location of the gridpoints. However, storing pointers
                  // to the gridpoint objects allow us to lookup/set the
                  // function value and state of the gridpoints. It is useful
                  // for updating the classification of the cell.

  int numSimplices;
  vector<MeshSimplex *> simplices;

  int numVertices;
  MeshVertex *vertices[12];
  string edgeNames[12];
  VertexContainer *vertCon; // global map of edge->mesh vertices. Used to check
                            // for duplicacy.

public:
  GridCell();
  GridCell(const GridCell &gc);
  GridCell(double a, double b, double c);
  GridCell(double a, double b, double c, double F);
  GridCell(double a, double b, double c, double F, int s);
  GridCell(double a, double b, double c, double F, int s, double d);
  GridCell(double a, double b, double c, double F, int s, double d,
           vector<GridPoint *> gps);

  ~GridCell();

  int getState();
  int getContourState();
  double getFunctionValue();
  int getGridCellID();
  double getCellDimension();
  vector<GridPoint *> getGridPoints();
  GridPoint *getGridPoint(int i);
  int getNumGridPoints();
  int getNumSimplices();

  void setState(int s);
  void markForDelete();
  void markForCreate();
  void markForUpdate();
  void setFunctionValue(double f);
  void setCellDimension(double d);
  void setGridPoints(vector<GridPoint *> gps);
  void addGridPoint(GridPoint *gp);
  void setVertexContainer(VertexContainer *vc);

  void setPosition(double a, double b, double c);

  int classifyCell();
  int classifyCell(vector<GridPoint *> gps);

  void destroyContour();
  bool computeContour();
  bool computeContour(double isoValue);

  void printMesh(std::ofstream *outputFile);
};

}; // namespace DynamicAdaptiveGrid

#endif
