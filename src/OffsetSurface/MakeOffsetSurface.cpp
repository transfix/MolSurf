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
#include <OffsetSurface/MakeOffsetSurface.h>
#include <GeometryFileTypes/GeometryLoader.h>

#include <cstdlib>
#include <cmath>

Geometry *MakeOffsetSurface(Geometry *surf, double dist) {

  cout << surf->m_NumTris << " " << surf->m_NumTriVerts << endl;

  for (int i = 0; i < surf->m_NumTriVerts; i++) {
    // cout << surf->m_TriVerts[3*i] << " "
    //  << surf->m_TriVerts[3*i+1] << " "
    //  << surf->m_TriVerts[3*i+2] << endl;

    double norm = sqrt(
        surf->m_TriVertNormals[3 * i] * surf->m_TriVertNormals[3 * i] +
        surf->m_TriVertNormals[3 * i + 1] * surf->m_TriVertNormals[3 * i + 1] +
        surf->m_TriVertNormals[3 * i + 2] * surf->m_TriVertNormals[3 * i + 2]);

    surf->m_TriVerts[3 * i] += surf->m_TriVertNormals[3 * i] * dist / norm;
    surf->m_TriVerts[3 * i + 1] +=
        surf->m_TriVertNormals[3 * i + 1] * dist / norm;
    surf->m_TriVerts[3 * i + 2] +=
        surf->m_TriVertNormals[3 * i + 2] * dist / norm;
  }

  return surf;
}

void usageOffsetSurface() {
  cout << endl
       << "Usage: MolSurf -offsetSurface <rawn input file> <rawn output file> "
          "<offset distance>"
       << endl
       << "    Create an offset mesh." << endl;
}
bool offsetSurface(int argc, char **argv) {
  if (argc != 5) {
    usageOffsetSurface();
    return false;
  }

  Geometry *geometry = GeometryLoader().loadFile(argv[2]);
  double distance = atof(argv[4]);
  geometry = MakeOffsetSurface(geometry, distance);

  GeometryLoader().saveFile(argv[3], geometry);
  delete geometry;

  return true;
}
