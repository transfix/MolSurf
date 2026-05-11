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
#include <InterfaceStats/OutlineInterface.h>
#include <algorithm>

OutlineInterface::OutlineInterface(string sf1, string sf2, double iw, bool rm,
                                   bool v) {
  verbose = v;

  receptorSurfaceFileName = sf1;
  ligandSurfaceFileName = sf2;

  interfaceWidth = iw;

  rawcMode = rm;
  filesRead = false;
  initialize();
}

OutlineInterface::~OutlineInterface() { cleanup(); }

void OutlineInterface::cleanup() {
  receptorPoints.clear();
  receptorPointColors.clear();

  ligandPoints.clear();
  ligandPointsXformed.clear();

  receptorTriangles.clear();
  ligandTriangles.clear();

  if (pgL)
    delete pgL;

  if (receptorPointsOnInterface)
    delete receptorPointsOnInterface;
  if (receptorPointsOnInterfaceBoundary)
    delete receptorPointsOnInterfaceBoundary;
}

bool OutlineInterface::initialize() {
  filesRead = false;

  FILE *receptorSurfaceFile = fopen(receptorSurfaceFileName.c_str(), "rt");
  if (!receptorSurfaceFile) {
    cout << "could not open surface file " << receptorSurfaceFileName.c_str()
         << endl;
    return false;
  }

  FILE *ligandSurfaceFile = fopen(ligandSurfaceFileName.c_str(), "rt");
  if (!ligandSurfaceFile) {
    cout << "could not open surface file " << ligandSurfaceFileName.c_str()
         << endl;
    fclose(receptorSurfaceFile);
    return false;
  }

  double x, y, z, nx, ny, nz;
  int index1, index2, index3;
  double dx, dy, dz;
  double a, b, c;
  double s;
  double area;

  if (fscanf(receptorSurfaceFile, "%d %d", &numVertR, &numTriangR) != 2) {
    cout << "Wrong format in receptor surface file\n";
    return false;
  }

  if (verbose)
    cout << "Receptor has: " << numVertR << " vertices and " << numTriangR
         << " triangles\n";

  receptorPointsOnInterfaceBoundary = new bool[numVertR];
  receptorPointsOnInterface = new bool[numVertR];

  for (int i = 0; i < numVertR; i++) {
    receptorPointsOnInterface[i] = false;
    receptorPointsOnInterfaceBoundary[i] = false;

    vertexcolor vc;

    if (rawcMode) {
      if (fscanf(receptorSurfaceFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &z,
                 &(vc.r), &(vc.g), &(vc.b)) != 6) {
        cout << "Wrong format in receptor surface file. At vertex no: " << i
             << "\n";
        cout << "x = " << x << ", y = " << y << ", z = " << z
             << ", r = " << vc.r << ", g = " << vc.g << ", b = " << vc.b
             << "\n";
        return false;
      }
    } else {
      if (fscanf(receptorSurfaceFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x,
                 &y, &z, &nx, &ny, &nz, &(vc.r), &(vc.g), &(vc.b)) != 9) {
        cout << "Wrong format in receptor surface file. At vertex no: " << i
             << "\n";
        cout << "x = " << x << ", y = " << y << ", z = " << z
             << ", r = " << vc.r << ", g = " << vc.g << ", b = " << vc.b
             << "\n";
        return false;
      }
    }

    DPG::WeightedPoint *p = new DPG::WeightedPoint(x, y, z, 0.0);
    receptorPoints.push_back(p);
    receptorPointColors.push_back(vc);
  }

  if (verbose)
    cout << "Points read from receptor surface and DPG created\n";

  for (int i = 0; i < numTriangR; i++) {
    if (fscanf(receptorSurfaceFile, "%d %d %d", &index1, &index2, &index3) !=
        3) {
      cout << "Wrong format in receptor surface file\n";
      return false;
    }

    triangle tr;
    tr.vertex[0] = index1;
    tr.vertex[1] = index2;
    tr.vertex[2] = index3;

    receptorTriangles.push_back(tr);
  }

  fclose(receptorSurfaceFile);

  if (verbose)
    cout << "Triangles read from receptor surface\n";

  if (fscanf(ligandSurfaceFile, "%d %d", &numVertL, &numTriangL) != 2) {
    cout << "Wrong format in ligand surface file\n";
    return false;
  }

  if (verbose)
    cout << "Ligand as: " << numVertL << " vertices and " << numTriangL
         << " triangles\n";

  for (int i = 0; i < numVertL; i++) {
    if (fscanf(ligandSurfaceFile, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx,
               &ny, &nz) != 6) {
      cout << "Wrong format in ligand surface file. At vertex no: " << i
           << "\n";
      cout << "x = " << x << ", y = " << y << ", z = " << z << ", nx = " << nx
           << ", ny = " << ny << ", nz = " << nz << "\n";
      return false;
    }

    DPG::Point *p = new DPG::Point(x, y, z);
    ligandPoints.push_back(p);
  }

  if (verbose)
    cout << "Points read from ligand surface\n";

  fclose(ligandSurfaceFile);

  return true;
}

bool OutlineInterface::detectInterfacePoints(Matrix xform) {
  pgL = new DPG::PG(10.0, 1000.0, 3.0);
  ligandPointsXformed.clear();

  double x, y, z;
  double newx, newy, newz;

  for (int i = 0; i < numVertL; i++) {
    x = ligandPoints[i]->getX();
    y = ligandPoints[i]->getY();
    z = ligandPoints[i]->getZ();

    newx = xform.get(0, 0) * x + xform.get(0, 1) * y + xform.get(0, 2) * z +
           xform.get(0, 3);
    newy = xform.get(1, 0) * x + xform.get(1, 1) * y + xform.get(1, 2) * z +
           xform.get(1, 3);
    newz = xform.get(2, 0) * x + xform.get(2, 1) * y + xform.get(2, 2) * z +
           xform.get(2, 3);

    DPG::Point *p = new DPG::Point(newx, newy, newz);
    pgL->addPoint(p);
    ligandPointsXformed.push_back(p);
  }

  if (verbose)
    cout << "DPG for ligand created\n";

  for (int i = 0; i < numVertR; i++) {
    receptorPointsOnInterfaceBoundary[i] = false;

    if (pgL->pointsWithinRange(receptorPoints[i], interfaceWidth))
      receptorPointsOnInterface[i] = true;
    else
      receptorPointsOnInterface[i] = false;
  }

  if (verbose)
    cout << "Receptor's interface with ligand computed\n";

  int vert1, vert2, vert3;
  for (int i = 0; i < numTriangR; i++) {
    vert1 = receptorTriangles[i].vertex[0];
    vert2 = receptorTriangles[i].vertex[1];
    vert3 = receptorTriangles[i].vertex[2];

    if (vert1 < 0 || vert2 < 0 || vert3 < 0 || vert1 >= numVertR ||
        vert2 >= numVertR || vert3 >= numVertR) {
      cout << "Error: Some vertex indices are out of bounds\n";
    } else {
      if ((receptorPointsOnInterface[vert1] ||
           receptorPointsOnInterface[vert2] ||
           receptorPointsOnInterface[vert3]) &&
          !(receptorPointsOnInterface[vert1] &&
            receptorPointsOnInterface[vert2] &&
            receptorPointsOnInterface[vert3])) {
        receptorPointsOnInterfaceBoundary[vert1] = true;
        receptorPointsOnInterfaceBoundary[vert2] = true;
        receptorPointsOnInterfaceBoundary[vert3] = true;
      }
    }
  }

  if (verbose)
    cout << "Boundary vertices detected\n";

  return true;
}

bool OutlineInterface::computeOutline(Matrix xform, string outputFName,
                                      vertexcolor outlineColor) {
  outlineFileName = outputFName;

  if (detectInterfacePoints(xform))
    return prepareBindingSiteOutlineSurfaceFile(outlineColor);
}

bool OutlineInterface::prepareBindingSiteOutlineSurfaceFile(
    vertexcolor outlineColor) {
  FILE *bindingSitesFile = fopen(outlineFileName.c_str(), "wt");
  if (!bindingSitesFile) {
    cout << "could not create surface file " << outlineFileName.c_str() << endl;
    return false;
  }

  // preparing overlap surface for receptor
  fprintf(bindingSitesFile, "%d %d\n", numVertR, numTriangR);

  for (int i = 0; i < numVertR; i++) {
    fprintf(bindingSitesFile, "%lf %lf %lf", receptorPoints[i]->getX(),
            receptorPoints[i]->getY(), receptorPoints[i]->getZ());

    if (receptorPointsOnInterfaceBoundary[i]) {
      fprintf(bindingSitesFile, " %lf %lf %lf\n", outlineColor.r,
              outlineColor.g, outlineColor.b);
    } else {
      fprintf(bindingSitesFile, " %lf %lf %lf\n", receptorPointColors[i].r,
              receptorPointColors[i].g, receptorPointColors[i].b);
    }
  }

  if (verbose)
    cout << "Points written to file\n";

  for (int i = 0; i < numTriangR; i++) {
    fprintf(bindingSitesFile, "%d %d %d\n", receptorTriangles[i].vertex[0],
            receptorTriangles[i].vertex[1], receptorTriangles[i].vertex[2]);
  }

  if (verbose)
    cout << "Triangles written to file\n";

  fclose(bindingSitesFile);

  return true;
}
