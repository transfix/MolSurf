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
#ifndef DAG_H
#define DAG_H

#include <DynamicAdaptiveGrid/MeshVertex.h>
#include <DynamicAdaptiveGrid/MeshSimplex.h>
#include <DynamicAdaptiveGrid/VertexContainer.h>
#include <DynamicAdaptiveGrid/GridPoint.h>
#include <DynamicAdaptiveGrid/GridCell.h>
#include <PDBParser/Atom.h>
#include <DPG/PG.h>

namespace DynamicAdaptiveGrid {
class DynamicAdaptiveGrids {
private:
  map<string, GridCell *> gridCellMap;
  map<string, GridPoint *> gridPointMap;

  vector<GridCell *> gridCells;
  vector<GridPoint *> gridPoints;

  map<string, GridCell *> newlyUpdated;
  vector<GridCell *> toBeDeleted;
  vector<GridCell *> toBeUpdated;
  vector<GridCell *> toBeCreated;

  DPG::PG *allAtomPG;
  DPG::PG *exposedAtomPG;
  DPG::PG *gridPointPG;
  DPG::PG *exposedGridPointPG;
  DPG::PG *gridCellPG;
  DPG::PG *exposedGridCellPG;

  //                      GridPoint*** gridPointArray;
  //                      GridCell*** gridCellArray;

  int mode;
  double gridSpacing;
  double solventRadius;
  double dpgCellSpacing;
  double dpgxlate;
  double isovalue;
  double minx, miny, minz, maxx, maxy, maxz;
  double xDim, yDim, zDim;
  int xCells, yCells, zCells;
  int numExposedAtoms;
  int numExposedCells;

  // void createGrid();
  // void destroyGrid();
  void init();
  void cleanup();

  // Finds gridpoints inside the bounding sphere. Creates new gridpoints if
  // necessary.
  void findGridPointsInBall(double x, double y, double z, double d,
                            vector<GridPoint *> &gps);

  void updateMinMax(PDBParser::Atom *a);
  void addSingleAtom(PDBParser::Atom *a);
  void removeSingleAtom(PDBParser::Atom *a);

  // recompute the contours inside updated cells
  bool updateContour();

  // marked cells are updated based on new states of its gridpoints.
  bool updateCells();

  // If any gridpoint is not already created, then it is created and marked as
  // outside the SAS. Then updates its own state.
  bool updateCell(GridCell *gc);

  // marks cells around a newly updated gridpoint for update
  bool markNeighborCells(GridPoint *gp);

protected:
  VertexContainer vertCon;

public:
  DynamicAdaptiveGrids(double grid_spacing, double pg_spacing, int mode,
                       double solvent_rad, double isovalue, double minx,
                       double miny, double minz, double maxx, double maxy,
                       double maxz);
  ~DynamicAdaptiveGrids();

  bool addAtom(PDBParser::Atom *a);
  bool removeAtom(PDBParser::Atom *a);
  bool moveAtom(PDBParser::Atom *from, PDBParser::Atom *to);

  bool addAtoms(vector<PDBParser::Atom *> &atms);
  bool removeAtoms(vector<PDBParser::Atom *> &atms);

  bool printSurface(string surfaceFileName);
};
}; // namespace DynamicAdaptiveGrid

#endif
