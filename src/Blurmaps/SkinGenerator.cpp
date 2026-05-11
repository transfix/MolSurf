/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit
 * Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is
 * free software; you can redistribute it and/or
  modify it under the terms of
 * the GNU Lesser General Public
  License version 2.1 as published by the Free
 * Software Foundation.


  MolSurf is distributed in the hope that it will be
 * useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of

 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General
 * Public License for more details.

  You should have received a copy of the
 * GNU General Public License
  along with MolSurf; if not, write to the Free
 * Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
*/

#include <Blurmaps/SkinGenerator.h>
#include <PDBParser/FlattenGOA.h>
#include <DPG/PG.h>

SkinGenerator::SkinGenerator(PDBParser::GroupOfAtoms *molecule, string surf,
                             double pr, double fbw, double cd,
                             PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  // BlurMapsDataManager::flattenGOA(molecule, m_AtomList, collectionData, 0, 0,
  // 0, radiusType, PDBParser::ATOM_TYPE, false);
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);

  surfaceFileName = surf;
  probeRadius = pr;
  floatingBandWidth = fbw;
  clashDistance = cd;
}

SkinGenerator::~SkinGenerator() { m_AtomList.clear(); }

bool SkinGenerator::populate(string xyzFileName) {
  FILE *surfaceFile, *xyzFile;

  int numVert, numTriang;

  double x, y, z, nx, ny, nz;

  double dist = probeRadius + floatingBandWidth;

  string xyzrFileName;
  xyzrFileName.assign(xyzFileName, 0, xyzFileName.length() - 3);
  xyzrFileName += "xyzr";
  FILE *xyzrFile;

  vector<DPG::Point *> atomlist;

  surfaceFile = fopen(surfaceFileName.c_str(), "rt");

  if (!surfaceFile) {
    cout << "could not open surface file" << endl;
    return false;
  }

  else {
    xyzFile = fopen(xyzFileName.c_str(), "wt");

    if (!xyzFile) {
      cout << "could not create xyz file" << endl;
      fclose(surfaceFile);
      return false;
    }

    xyzrFile = fopen(xyzrFileName.c_str(), "wt");

    if (!xyzrFile) {
      cout << "could not create xyzr file" << endl;
      fclose(surfaceFile);
      return false;
    }

    // Adding grown skin using the surface

    int retval = 0;
    retval = fscanf(surfaceFile, "%d", &numVert);
    retval = fscanf(surfaceFile, "%d", &numTriang);

    DPG::PG *pg = new DPG::PG(10.0, 1000.0, 3.0);
    vector<DPG::Point *> results;

    int yesCount = 0;
    int noCount = 0;

    for (int i = 0; i < numVert; i++) {
      retval = fscanf(surfaceFile, "%lf", &x);
      retval = fscanf(surfaceFile, "%lf", &y);
      retval = fscanf(surfaceFile, "%lf", &z);
      retval = fscanf(surfaceFile, "%lf", &nx);
      retval = fscanf(surfaceFile, "%lf", &ny);
      retval = fscanf(surfaceFile, "%lf", &nz);

      x = x + nx * dist;
      y = y + ny * dist;
      z = z + nz * dist;

      DPG::Point *p = new DPG::Point(x, y, z);

      results.clear();
      results = pg->range(p, clashDistance);

      if (results.size() == 0) {
        yesCount++;

        pg->addPoint(p);

        atomlist.push_back(p);
      } else
        noCount++;
    }

    cout << "Accepted " << yesCount << " skin molecules" << endl;
    cout << "Rejected " << noCount << " skin molecules" << endl;

    // printing

    int numAtoms = m_AtomList.size();

    fprintf(xyzFile, "%d\n", numAtoms + yesCount);

    for (int i = 0; i < numAtoms; i++) {
      fprintf(xyzFile, "I %lf %lf %lf %lf\n", m_AtomList[i]->m_Position[0],
              m_AtomList[i]->m_Position[1], m_AtomList[i]->m_Position[2],
              m_AtomList[i]->getCharge());
      fprintf(xyzrFile, "%lf %lf %lf %lf\n", m_AtomList[i]->m_Position[0],
              m_AtomList[i]->m_Position[1], m_AtomList[i]->m_Position[2],
              m_AtomList[i]->getRadius());
    }

    for (int i = 0; i < yesCount; i++) {
      fprintf(xyzFile, "E %lf %lf %lf %lf\n", atomlist[i]->getX(),
              atomlist[i]->getY(), atomlist[i]->getZ(), 0.0);
      fprintf(xyzrFile, "%lf %lf %lf %lf\n", atomlist[i]->getX(),
              atomlist[i]->getY(), atomlist[i]->getZ(), 1.1);
    }

    fclose(surfaceFile);
    fclose(xyzFile);
    fclose(xyzrFile);
    return true;
  }
}

bool SkinGenerator::populate(string xyzFileName, string xyzrFileName) {
  populate(xyzFileName);
}

bool SkinGenerator::populate(string xyzFileName, string xyzrFileName,
                             string f2dFileName) {
  populate(xyzFileName);
}
