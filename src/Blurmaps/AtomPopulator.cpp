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
#include <Blurmaps/AtomPopulator.h>
#include <PDBParser/FlattenGOA.h>

AtomPopulator::AtomPopulator(PDBParser::GroupOfAtoms *molecule, int dim1,
                             int dim2, int dim3, double probeRadius,
                             PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  m_ProbeRadius = probeRadius;
  m_floatingBandWidth = 0.0;
  m_RadiusType = radiusType;
  m_SasInteriorVertex = 0;

  // "flatten" the GOA
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  // BlurMapsDataManager::flattenGOA(molecule, m_AtomList, collectionData, 0, 0,
  // 0, radiusType, PDBParser::ATOM_TYPE, false);
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);
  // compute the bounding box
  {
    m_Min[0] = m_Min[1] = m_Min[2] = 0.;
    m_Max[0] = m_Max[1] = m_Max[2] = 0.;
    BlurMapsDataManager::getBoundingBox(m_AtomList, m_Min, m_Max, m_RadiusType,
                                        m_ProbeRadius * 2 + 1, 0);
    m_Dim[0] = dim1;
    m_Dim[1] = dim2;
    m_Dim[2] = dim3;
    // origin,
    m_Orig[0] = m_Min[0];
    m_Orig[1] = m_Min[1];
    m_Orig[2] = m_Min[2];
    // and span
    m_Span[0] = (m_Max[0] - m_Min[0]) / (float)(dim1 - 1);
    m_Span[1] = (m_Max[1] - m_Min[1]) / (float)(dim2 - 1);
    m_Span[2] = (m_Max[2] - m_Min[2]) / (float)(dim3 - 1);
  }
}

AtomPopulator::AtomPopulator(PDBParser::GroupOfAtoms *molecule, int dim1,
                             int dim2, int dim3, double probeRadius,
                             double floatingBandWidth,
                             PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  m_ProbeRadius = probeRadius;
  m_floatingBandWidth = floatingBandWidth;
  m_RadiusType = radiusType;
  m_SasInteriorVertex = 0;

  // "flatten" the GOA
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  // BlurMapsDataManager::flattenGOA(molecule, m_AtomList, collectionData, 0, 0,
  // 0, radiusType, PDBParser::ATOM_TYPE, false);
  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);
  // compute the bounding box
  {
    m_Min[0] = m_Min[1] = m_Min[2] = 0.;
    m_Max[0] = m_Max[1] = m_Max[2] = 0.;
    BlurMapsDataManager::getBoundingBox(m_AtomList, m_Min, m_Max, m_RadiusType,
                                        m_ProbeRadius * 2 + 1, 0);
    m_Dim[0] = dim1;
    m_Dim[1] = dim2;
    m_Dim[2] = dim3;
    // origin,
    m_Orig[0] = m_Min[0];
    m_Orig[1] = m_Min[1];
    m_Orig[2] = m_Min[2];
    // and span
    m_Span[0] = (m_Max[0] - m_Min[0]) / (float)(dim1 - 1);
    m_Span[1] = (m_Max[1] - m_Min[1]) / (float)(dim2 - 1);
    m_Span[2] = (m_Max[2] - m_Min[2]) / (float)(dim3 - 1);
  }
}

AtomPopulator::~AtomPopulator() {
  if (m_SasInteriorVertex) {
    delete[] m_SasInteriorVertex;
  }
  m_SasInteriorVertex = 0;
  if (m_SasBoundaryCells) {
    delete[] m_SasBoundaryCells;
  }
  m_SasBoundaryCells = 0;
  {
    int n = m_AtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_AtomList[i];
    }
    m_AtomList.clear();
  }
  {
    int n = m_SkinAtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_SkinAtomList[i];
    }
    m_SkinAtomList.clear();
  }
}

bool AtomPopulator::addAtom(int cellx, int celly, int cellz, bool *atomPresent,
                            int numSpheres) {
  // this is to help in double convolution. Comment out otherwise.
  //     if( numSpheres > 2 ) return true;
  //     return false;
  int iRadiusExtent = m_ProbeRadius / m_Span[0];
  int jRadiusExtent = m_ProbeRadius / m_Span[1];
  int kRadiusExtent = m_ProbeRadius / m_Span[2];
  int i, j, k;
  int mini, minj, mink;
  int maxi, maxj, maxk;
  mini = cellx - iRadiusExtent;
  if (mini < 0) {
    mini = 0;
  }
  minj = celly - jRadiusExtent;
  if (minj < 0) {
    minj = 0;
  }
  mink = cellz - kRadiusExtent;
  if (mink < 0) {
    mink = 0;
  }
  maxi = cellx + iRadiusExtent;
  if (maxi >= m_Dim[0] - 1) {
    maxi = m_Dim[0] - 1 - 1;
  }
  maxj = celly + jRadiusExtent;
  if (maxj >= m_Dim[1] - 1) {
    maxj = m_Dim[1] - 1 - 1;
  }
  maxk = cellz + kRadiusExtent;
  if (maxk >= m_Dim[2] - 1) {
    maxk = m_Dim[2] - 1 - 1;
  }
  bool stericClash = false;
  for (i = mini; i <= maxi; i++) {
    for (j = minj; j <= maxj; j++) {
      for (k = mink; k <= maxk; k++) {
        int c = k * (m_Dim[0] - 1) * (m_Dim[1] - 1) + j * (m_Dim[0] - 1) + i;
        if (atomPresent[c]) {
          stericClash = true;
        }
        {}
      }
    }
  }
  return !stericClash;
}

double lookUpCharge(const PDBParser::Atom &atom) {
  if (atom.m_ElementIndex >= 0 &&
      atom.m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[atom.m_ElementIndex].charge;
  }
  return 0; // unknown atom. Just say neutral and return 0.
}

bool AtomPopulator::printSkin(const char *filename, bool printMiscFiles) {
  if (printMiscFiles) {
    char skinfname[256];
    strcpy(skinfname, filename);
    strcat(skinfname, "_pseudo_skin.pdb");
    FILE *fp = fopen(skinfname, "w");
    if (!fp) {
      return false;
    }
    int i;
    for (i = 0; i < m_SkinAtomList.size(); i++) {
      int atomNum = i % 100000;
      PDBParser::Atom *at = m_SkinAtomList[i];
      fprintf(fp, "ATOM  ");       // "ATOM "
      fprintf(fp, "%5d", atomNum); // serial
      fprintf(fp, " ");            // char 12 is blank
      const char *atom_name = at->getName();
      if (strlen(atom_name) >= 4) {
        fprintf(fp, "%-4.4s", at->getName());
      } else {
        fprintf(fp, " %-3.3s", at->getName());
      }
      fprintf(fp, " ");                         // altLoc
      fprintf(fp, "%3s", at->getResidueName()); // resName
      fprintf(fp, " ");                         // char 21 is blank
      fprintf(fp, "%c", at->getChain());        // chainID
      fprintf(fp, "%4d", at->getResidueNum());  // resSeq
      fprintf(fp, " ");                         // iCode
      fprintf(fp, "   ");                       // chars 28, 29, 30 are blank
      fprintf(fp, "%8.3f", at->m_Position[0]);  // x
      fprintf(fp, "%8.3f", at->m_Position[1]);  // y
      fprintf(fp, "%8.3f", at->m_Position[2]);  // z
      fprintf(fp, "%6.2f", at->getRadius());    // occupancy
      fprintf(fp, "      ");                    // tempFactor
      fprintf(fp, "    ");                      // segID
      fprintf(fp, "      ");                    // chars 67-72 are blank
      fprintf(fp, "XX");
      fprintf(fp, "  "); // charge
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  {
    // write docking input file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      return false;
    }
    fprintf(fp, "%lu\n", m_AtomList.size() + m_SkinAtomList.size());
    int i;
    for (i = 0; i < m_AtomList.size(); i++) {
      PDBParser::Atom *at = m_AtomList[i];
      fprintf(fp, "I     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0],
              at->m_Position[1], at->m_Position[2], at->getCharge());
    }
    for (i = 0; i < m_SkinAtomList.size(); i++) {
      PDBParser::Atom *at = m_SkinAtomList[i];
      fprintf(fp, "E     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0],
              at->m_Position[1], at->m_Position[2], at->getCharge());
    }
    fclose(fp);
  }
  return true;
}

bool AtomPopulator::getSasInterior() {
  if (!m_SasInteriorVertex) {
    return false;
  }
  {
    int c;
    for (c = 0; c < m_Dim[0] * m_Dim[1] * m_Dim[2]; c++) {
      m_SasInteriorVertex[c] = 0;
    }
  }
  unsigned int zoff = 0;
  int i;
  for (i = 0; i < m_AtomList.size(); i++) {
    PDBParser::Atom *at = m_AtomList[i];
    float radius;
    radius = at->getRadius(m_RadiusType);
    radius += m_ProbeRadius;
    double c[3], maxRad;
    const float eps = 1.0e-3f;
    unsigned int amax[3], amin[3], j, k, l;
    // calculate the maximum radius for this atom
    maxRad = radius + 1;
    // compute the dataset coordinates of the atom's center
    for (j = 0; j < 3; j++) {
      c[j] = (at->m_Position[j] - m_Orig[j]) / m_Span[j];
      c[j] = ((c[j] - floor(c[j])) >= 0.5) ? ceil(c[j]) : floor(c[j]);
    }
    // then compute the bounding box of the atom (maxRad^3)
    for (j = 0; j < 3; j++) {
      int tempAMin = (int)(c[j] - (maxRad / m_Span[j]) - 1);
      tempAMin = (tempAMin < 0) ? 0 : tempAMin;
      amin[j] = tempAMin;
      int tempAMax = (int)(c[j] + (maxRad / m_Span[j]) + 1);
      tempAMax = (tempAMax > m_Dim[j]) ? m_Dim[j] : tempAMax;
      amax[j] = tempAMax;
    }
    c[0] = at->m_Position[0]; // m_Orig[0] + c[0]*m_Span[0];
    c[1] = at->m_Position[1]; // m_Orig[1] + c[1]*m_Span[1];
    c[2] = at->m_Position[2]; // m_Orig[2] + c[2]*m_Span[2];
    // begin blurring kernel
    for (l = amin[2]; l < amax[2]; l++) {
      for (k = amin[1]; k < amax[1]; k++) {
        for (j = amin[0]; j < amax[0]; j++) {
          int n = j + k * m_Dim[0] + (l - zoff) * m_Dim[0] * m_Dim[1];
          if (n < 0 || n >= (m_Dim[0] * m_Dim[1] * m_Dim[2])) {
            printf("Out of range index");
          }
          float pnt[3];
          pnt[0] = m_Orig[0] + j * m_Span[0];
          pnt[1] = m_Orig[1] + k * m_Span[1];
          pnt[2] = m_Orig[2] + l * m_Span[2];
          double r = (c[0] - pnt[0]) * (c[0] - pnt[0]) +
                     (c[1] - pnt[1]) * (c[1] - pnt[1]) +
                     (c[2] - pnt[2]) * (c[2] - pnt[2]);
          if (r < radius * radius) {
            m_SasInteriorVertex[n]++;
          }
        }
      }
    }
    // end blurring kernel
  }
  return true;
}

bool AtomPopulator::getSasBoundary(unsigned char *sasBoundary) {
  if (!sasBoundary || !m_SasInteriorVertex) {
    return false;
  }
  int c = 0;
  int i, j, k;
  for (i = 0; i < m_Dim[0]; i++) {
    for (j = 0; j < m_Dim[1]; j++) {
      for (k = 0; k < m_Dim[2]; k++) {
        sasBoundary[c] = 0;
        // if is boundary, make boundary vertex true;
        if (!m_SasInteriorVertex[c]) {
          if (isBoundary(i, j, k, m_SasInteriorVertex)) {
            sasBoundary[c] = m_SasInteriorVertex[c];
          }
        }
        c++;
      }
    }
  }
  return true;
}

bool AtomPopulator::isBoundaryCell(unsigned char *sasCellInterior) {
  bool interiorFlag, exteriorFlag;
  interiorFlag = exteriorFlag = false;
  int i;
  for (i = 0; i < 8; i++) {
    if (sasCellInterior[i]) {
      interiorFlag = true;
    } else {
      exteriorFlag = true;
    }
  }
  return (interiorFlag & exteriorFlag);
}

unsigned char AtomPopulator::getMaxInArray(unsigned char *sasCellInterior,
                                           int n) {
  if (!sasCellInterior || n < 1) {
    return 0;
  }
  unsigned char maxVal = sasCellInterior[0];
  int i;
  for (i = 1; i < n; i++) {
    if (sasCellInterior[i] > maxVal) {
      maxVal = sasCellInterior[i];
    }
  }
  return maxVal;
}

bool AtomPopulator::getSasBoundaryCells() {
  if (!m_SasBoundaryCells || !m_SasInteriorVertex) {
    return false;
  }
  int numberCell, i, j;
  numberCell = (m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1);
  int cellVertexArray[8];
  unsigned char sasCellInterior[8];
  for (i = 0; i < numberCell; i++) {
    m_SasBoundaryCells[i] = 0;
    TrilinearGrid::getCellVertices(i, cellVertexArray, m_Dim);
    for (j = 0; j < 8; j++) {
      sasCellInterior[j] = m_SasInteriorVertex[cellVertexArray[j]];
    }
    if (isBoundaryCell(sasCellInterior)) {
      m_SasBoundaryCells[i] = getMaxInArray(sasCellInterior, 8);
    }
  }
  return true;
}

bool AtomPopulator::getBoundaryCells() {
  // get the SAS boundary cells
  unsigned char *sasBoundary =
      new unsigned char[m_Dim[0] * m_Dim[1] * m_Dim[2]];
  {
    if (!m_SasInteriorVertex) {
      m_SasInteriorVertex = new unsigned char[m_Dim[0] * m_Dim[1] * m_Dim[2]];
    }
    if (!getSasInterior()) {
      delete[] m_SasInteriorVertex;
      m_SasInteriorVertex = 0;
      delete[] sasBoundary;
      sasBoundary = 0;
      return false;
    }
    if (!getSasBoundary(sasBoundary)) {
      delete[] m_SasInteriorVertex;
      m_SasInteriorVertex = 0;
      delete[] sasBoundary;
      sasBoundary = 0;
      return false;
    }
    if (!getSasBoundaryCells()) {
      delete[] m_SasInteriorVertex;
      m_SasInteriorVertex = 0;
      delete[] sasBoundary;
      sasBoundary = 0;
      return false;
    }
  }
  return true;
}

bool AtomPopulator::populateAtomList() {
  m_SasBoundaryCells =
      new unsigned char[(m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1)];
  if (!getBoundaryCells()) {
    delete[] m_SasBoundaryCells;
    m_SasBoundaryCells = 0;
  }
  //////////////////////////////////////////////////////
  //
  //  1. Create boundary cells array
  //  2. Mix them up
  //  3. Create a atom present array. Could reuse boundary cells array if
  //  needed.
  //  4. Choose one cell
  //  5.     Depending on hueristic, add atom or not.
  //
  //////////////////////////////////////////////////////
  // 1. Create boundary cells array
  int numBoundaryCells = 0;
  int *boundaryCellIndexArray = 0;
  {
    int c = 0;
    int numberCell = (m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1);
    for (c = 0; c < numberCell; c++) {
      if (m_SasBoundaryCells[c]) {
        numBoundaryCells++;
      }
    }
  }
  {
    if (numBoundaryCells == 0) {
      return false;
    }
    boundaryCellIndexArray = new int[numBoundaryCells];
  }
  {
    int c = 0;
    int numberCell = (m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1);
    int count = 0;
    for (c = 0; c < numberCell; c++) {
      if (m_SasBoundaryCells[c]) {
        boundaryCellIndexArray[count] = c;
        count++;
      }
    }
  }
  // 2. Mix them up
  {
    int c;
    for (c = 0; c < numBoundaryCells; c++) {
      double r = (double)rand() / ((double)(RAND_MAX));
      int index = r * numBoundaryCells;
      if (index >= numBoundaryCells) {
        index = numBoundaryCells - 1;
      }
      if (index < 0) {
        index = 0;
      }
      // swap current and the random index!
      int temp = boundaryCellIndexArray[c];
      boundaryCellIndexArray[c] = boundaryCellIndexArray[index];
      boundaryCellIndexArray[index] = temp;
    }
  }
  // 3. Create a sphere present array and init it to false
  bool *atomPresent =
      new bool[(m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1)];
  {
    int c = 0;
    int numberCell = (m_Dim[0] - 1) * (m_Dim[1] - 1) * (m_Dim[2] - 1);
    for (c = 0; c < numberCell; c++) {
      atomPresent[c] = false;
    }
  }
  // 3. See if you want to add it depending on a hueristic
  int numadded = 0;
  int numskipped = 0;
  {
    int c;
    for (c = 0; c < numBoundaryCells; c++) {
      int cellIndex = boundaryCellIndexArray[c];
      int numSpheres = (int)(m_SasBoundaryCells[cellIndex]);
      int x, y, z;
      TrilinearGrid::cell2xyz(cellIndex, x, y, z, m_Dim);
      if (addAtom(x, y, z, atomPresent, numSpheres)) {
        PDBParser::Atom *atom = new PDBParser::Atom();
        float pnt[3];
        pnt[0] = m_Orig[0] + (x + 0.5) * m_Span[0];
        pnt[1] = m_Orig[1] + (y + 0.5) * m_Span[1];
        pnt[2] = m_Orig[2] + (z + 0.5) * m_Span[2];
        atom->m_Position[0] = pnt[0];
        atom->m_Position[1] = pnt[1];
        atom->m_Position[2] = pnt[2];
        atom->setRadius(m_ProbeRadius - m_floatingBandWidth);
        m_SkinAtomList.push_back(atom);
        atomPresent[cellIndex] = true;
        numadded++;
      }
      // else
      //    numskipped++;
    }
  }
  return true;
}

bool AtomPopulator::populate(const char *optFileName, bool printMiscFiles) {
  if (!populateAtomList()) {
    int n = m_SkinAtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_SkinAtomList[i];
    }
    m_SkinAtomList.clear();
    return false;
  }
  int n = m_SkinAtomList.size();
  if (!printSkin(optFileName, printMiscFiles)) {
    return false;
  }
  {
    int n = m_SkinAtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_SkinAtomList[i];
    }
    m_SkinAtomList.clear();
  }
  return true;
}

bool AtomPopulator::populate(PDBParser::GroupOfAtoms *molecule) {
  if (!populateAtomList()) {
    int n = m_SkinAtomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete m_SkinAtomList[i];
    }
    m_SkinAtomList.clear();
    return false;
  }
  PDBParser::GroupOfAtoms *firstResidue =
      molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
  int n = m_SkinAtomList.size();
  int i;
  for (i = 0; i < n; i++) {
    PDBParser::Atom *atom = m_SkinAtomList[i];
    if (!atom) {
      continue;
    }
    firstResidue->addAtom(
        atom->getAtomNum(), atom->getElementIndex(), atom->getResidueNum(),
        (char *)atom->getName(), atom->getResidueId(), atom->getChain(),
        atom->m_Position[0], atom->m_Position[1], atom->m_Position[2],
        atom->getRadius(), atom->getCharge(), false);
    delete m_SkinAtomList[i];
  }
  m_SkinAtomList.clear();
  return true;
}
