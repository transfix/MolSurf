/*
  Copyright 2011 The University of Texas at Austin

        Author: Zhang Qin <zqyork@ices.utexas.edu>
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
#include <FastMolSurf/fastMolSurf.h>
#include <PDBParser/FlattenGOA.h>

using namespace CVCUtility;
using namespace PDBParser;

// TODO: add these functions to the FastMoleSurface class
//       move commandline parsing to MolSurf/main
//       add functionality to MolSurfAPI
//       test the code and add tests to the test scripts...

float VolumeExcludedFunction(float dis, float radius, float width) {
  if (width >= radius) {
    width = radius;
  }
  if (dis <= radius - width) {
    return 0;
  } else if (dis >= radius + width) {
    return 1;
  } else if (width > 0) {
    return (0.5 + 0.75 * (dis - radius) / width -
            0.25 * (dis - radius) * (dis - radius) * (dis - radius) /
                (width * width * width));
  } else {
    error("Width is less than 0!");
  }
}

bool getAtomListAndExtent(GroupOfAtoms *molecule,
                          vector<PDBParser::Atom *> &atomList, float *minExt,
                          float *maxExt) {
  if (!molecule) {
    return false;
  }
  GroupOfAtoms::RADIUS_TYPE radiusType = GroupOfAtoms::VDW_RADIUS;
  // vector <PDBParser::Atom*> atomList;
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  FlattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);
  float stapnt[3] = {0.0, 0.0, 0.0};
  float endpnt[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    minExt[i] = 0.0;
    maxExt[i] = 0.0;
  }
  for (int i = 0; i < atomList.size(); i++) {
    float radius = atomList[i]->getRadius();
    for (int j = 0; j < 3; j++) {
      stapnt[j] = (atomList[i]->m_Position[j] - radius) < stapnt[j]
                      ? atomList[i]->m_Position[j] - radius
                      : stapnt[j];
      endpnt[j] = (atomList[i]->m_Position[j] + radius) > endpnt[j]
                      ? atomList[i]->m_Position[j] + radius
                      : endpnt[j];
    }
  }
  float scale = 1.2;
  for (int j = 0; j < 3; j++) {
    minExt[j] =
        (endpnt[j] + stapnt[j]) / 2.0 - (endpnt[j] - stapnt[j]) * scale / 2.0;
    maxExt[j] = minExt[j] + (endpnt[j] - stapnt[j]) * scale;
  }
  return true;
}

void getVolumeData(float *data, unsigned int *dim,
                   vector<PDBParser::Atom *> &atomList, float *minExt, float dx,
                   bool SAS) {
  double probeRadius = 1.4;
  float radius;
  for (int i = 0; i < atomList.size(); i++) {
    radius = atomList[i]->getRadius();
    if (SAS) {
      radius += probeRadius;
    }
    float width = 1.0;
    // int range1 = (int)((radius-width)/dx);
    int range = (int)((radius + width) / dx) + 2;
    int a[3];
    for (int j = 0; j < 3; j++) {
      a[j] = (int)((atomList[i]->m_Position[j] - minExt[j]) / dx);
    }
    int u[3];
    float p[3];
    int index;
    for (u[0] = a[0] - range; u[0] < a[0] + range; u[0]++) {
      if (u[0] < 0) {
        continue;
      }
      for (u[1] = a[1] - range; u[1] < a[1] + range; u[1]++) {
        if (u[1] < 0) {
          continue;
        }
        for (u[2] = a[2] - range; u[2] < a[2] + range; u[2]++) {
          if (u[2] < 0) {
            continue;
          }
          index = u[0] + u[1] * dim[0] + u[2] * dim[0] * dim[1];
          float distance = 0.0;
          for (int j = 0; j < 3; j++) {
            p[j] = minExt[j] + u[j] * dx;
            distance += (p[j] - atomList[i]->m_Position[j]) *
                        (p[j] - atomList[i]->m_Position[j]);
          }
          distance = sqrt(distance);
          // cout<<"function value: " << VolumeExcludedFunction(distance,
          // radius, width)<<endl;
          data[index] *= VolumeExcludedFunction(distance, radius, width);
        }
      }
    }
  }
}

// interface needs to be cleaned up to match everything else...
int constructFastMoleSurface(int argc, char *argv[]) {
  /*
        if(argc<7)
        {
                cout<<"Usage: ./FastMoleSurface <input.pdb> <dim0> <dim1> <dim2>
     <bool SAS> <output.raw>" <<endl; return 1;
        }
  */
  GroupOfAtoms *molecule = new GroupOfAtoms();
  GOALoader *goaload = new GOALoader();
  molecule = goaload->loadFile(argv[2], false);
  vector<PDBParser::Atom *> atomList;
  float minExt[3], maxExt[3];
  unsigned int dim[3];
  dim[0] = atoi(argv[3]);
  dim[1] = atoi(argv[4]);
  dim[2] = atoi(argv[5]);
  if (!getAtomListAndExtent(molecule, atomList, minExt, maxExt)) {
    return 0;
  }
  cout << "Molecule " << argv[2] << " has " << atomList.size() << " atoms."
       << endl;
  float interval[3];
  for (int i = 0; i < 3; i++) {
    interval[i] = maxExt[i] - minExt[i];
  }
  float maxInterval;
  int maxDim;
  maxInterval = Max_Of_Three(interval[0], interval[1], interval[2]);
  maxDim = Max_Of_Three(dim[0], dim[1], dim[2]);
  float dx = maxInterval / (float)maxDim;
  for (int i = 0; i < 3; i++) {
    dim[i] = (int)(interval[i] / dx) + 1;
    maxExt[i] = minExt[i] + (dim[i] - 1) * dx;
  }
  float *data = new float[dim[0] * dim[1] * dim[2]];
  for (int i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    data[i] = 1.0;
  }
  bool SAS = false;
  if (strcmp(argv[6], "true") == 0) {
    SAS = true;
  }
  getVolumeData(data, dim, atomList, minExt, dx, SAS);
  SimpleVolumeData *sdata = new SimpleVolumeData(dim);
  sdata->setDimensions(dim);
  sdata->setNumberOfVariables(1);
  sdata->setData(0, data);
  sdata->setType(0, SimpleVolumeData::FLOAT);
  sdata->setName(0, "Solvent Excluded Function Method");
  sdata->setMinExtent(minExt);
  sdata->setMaxExtent(maxExt);
  Geometry *geome = 0;
  geome = SimpleVolumeDataIsocontourer::getIsocontour(sdata, 0.5);
  GeometryLoader().saveFile(argv[7], geome);
  delete geome;
  delete sdata; // sdata has deleted the data for me. It shouldn't work like
                // this.
  cout << "done!" << endl;
  return 0;
}
