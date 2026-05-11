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

#include <Blurmaps/F2dGenerator.h>
#include <PDBParser/FlattenGOA.h>
#include <DPG/PG.h>

F2DGenerator::F2DGenerator(PDBParser::GroupOfAtoms *molecule, string xyzf,
                           string f2df, bool r,
                           PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                           double sar) {
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }

  FlattenGOA(molecule, m_AtomList, collectionData, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);

  xyzFileName = xyzf;
  outputFileName = f2df;
  receptor = r;
  skinAtomRad = sar;
}

F2DGenerator::~F2DGenerator() { m_AtomList.clear(); }

double F2DGenerator::diff(double a, double b) {
  double c = 0;

  if (a > b)
    c = a - b;
  else if (b > a)
    c = b - a;

  return c;
}

void F2DGenerator::printAtom(int i, bool internal) {
  int atomNum = i % 100000;
  PDBParser::Atom *at = m_AtomList[i];

  // writing f2d
  fprintf(f2dFile, "ATOM   ");          // "ATOM "
  fprintf(f2dFile, "%5d", atomNum + 1); // serial
  fprintf(f2dFile, " ");                // char 12 is blank
  const char *atom_name = at->getName();
  if (strlen(atom_name) >= 4) {
    fprintf(f2dFile, "%-4.4s", at->getName());
  } else {
    fprintf(f2dFile, " %-3.3s", at->getName());
  }
  fprintf(f2dFile, " ");                         // altLoc
  fprintf(f2dFile, "%3s", at->getResidueName()); // resName
  fprintf(f2dFile, " ");                         // char 21 is blank
  fprintf(f2dFile, " ");                         // chain is blank
  //      fprintf(f2dFile, "%c", at->getChain());         // chainID
  fprintf(f2dFile, "%4d", at->getResidueNum()); // resSeq
  fprintf(f2dFile, " ");                        // iCode
  fprintf(f2dFile, "   ");                      // chars 28, 29, 30 are blank
  fprintf(f2dFile, "%8.3f", at->m_Position[0]); // x
  fprintf(f2dFile, "%8.3f", at->m_Position[1]); // y
  fprintf(f2dFile, "%8.3f", at->m_Position[2]); // z
  fprintf(f2dFile, " ");
  fprintf(f2dFile, "% 6.4f", at->getCharge()); // charge
  fprintf(f2dFile, "%7.4f", at->getRadius());  // occupancy
  fprintf(f2dFile, " ");
  if (internal)
    fprintf(f2dFile, "I");
  else
    fprintf(f2dFile, "E");
  fprintf(f2dFile, "\n");

  // writing pqr
  fprintf(newpqrFile, "ATOM  ");           // "ATOM "
  fprintf(newpqrFile, "%5d", atomNum + 1); // serial
  fprintf(newpqrFile, " ");                // char 12 is blank
  if (strlen(atom_name) >= 4) {
    fprintf(newpqrFile, "%-4.4s", at->getName());
  } else {
    fprintf(newpqrFile, " %-3.3s", at->getName());
  }
  fprintf(newpqrFile, " ");                         // altLoc
  fprintf(newpqrFile, "%3s", at->getResidueName()); // resName
  fprintf(newpqrFile, " ");                         // char 21 is blank
  fprintf(newpqrFile, "%c", at->getChain());        // chainID
  fprintf(newpqrFile, "%4d", at->getResidueNum());  // resSeq
  fprintf(newpqrFile, " ");                         // iCode
  fprintf(newpqrFile, "   ");                      // chars 28, 29, 30 are blank
  fprintf(newpqrFile, "%8.3f", at->m_Position[0]); // x
  fprintf(newpqrFile, "%8.3f", at->m_Position[1]); // y
  fprintf(newpqrFile, "%8.3f", at->m_Position[2]); // z
  fprintf(newpqrFile, " ");
  fprintf(newpqrFile, "% 6.4f", at->getCharge()); // charge
  fprintf(newpqrFile, "%7.4f", at->getRadius());  // occupancy
  fprintf(newpqrFile, "\n");
}

void F2DGenerator::printSkinAtom(double x, double y, double z, int i) {
  int atomNum = i % 100000;
  PDBParser::Atom *at = m_AtomList[i];

  // writing f2d
  fprintf(f2dFile, "HETATM ");          // "ATOM "
  fprintf(f2dFile, "%5d", atomNum + 1); // serial
  fprintf(f2dFile, "  O   HOH  ");
  fprintf(f2dFile, "%4d", numres); // resSeq
  fprintf(f2dFile, " ");           // iCode
  fprintf(f2dFile, "   ");         // chars 28, 29, 30 are blank
  fprintf(f2dFile, "%8.3f", x);    // x
  fprintf(f2dFile, "%8.3f", y);    // y
  fprintf(f2dFile, "%8.3f", z);    // z
  fprintf(f2dFile, "  0.0000 %5.4f E", skinAtomRad);
  fprintf(f2dFile, "\n");

  // writing pqr
  fprintf(newpqrFile, "HETATM ");          // "ATOM "
  fprintf(newpqrFile, "%5d", atomNum + 1); // serial
  fprintf(newpqrFile, "  O   HOH 9");
  fprintf(newpqrFile, "%4d", numres); // resSeq
  fprintf(newpqrFile, " ");           // iCode
  fprintf(newpqrFile, "   ");         // chars 28, 29, 30 are blank
  fprintf(newpqrFile, "%8.3f", x);    // x
  fprintf(newpqrFile, "%8.3f", y);    // y
  fprintf(newpqrFile, "%8.3f", z);    // z
  fprintf(newpqrFile, "  0.0000 %5.4f E", skinAtomRad);
  fprintf(newpqrFile, "\n");
}

bool F2DGenerator::generate() {
  string newpqrFileName;
  newpqrFileName.assign(xyzFileName, 0, xyzFileName.length() - 4);
  newpqrFileName += "_f2d.pqr";

  FILE *xyzFile = fopen(xyzFileName.c_str(), "rt");

  if (!xyzFile) {
    cout << "could not open xyz file" << endl;
    return false;
  }

  else {
    f2dFile = fopen(outputFileName.c_str(), "wt");

    if (!f2dFile) {
      cout << "could not create f2d file" << endl;
      fclose(xyzFile);
      return false;
    }

    newpqrFile = fopen(newpqrFileName.c_str(), "wt");

    if (!newpqrFile) {
      cout << "could not create pqr file" << endl;
      fclose(xyzFile);
      fclose(f2dFile);
      return false;
    }

    // insert all internal atom coords into DPG
    int numXYZ;
    double x, y, z, q;

    DPG::PG *pg = new DPG::PG(10.0, 1000.0, 3.0);
    vector<DPG::Point *> skinAtomList;

    int retval = fscanf(xyzFile, "%d", &numXYZ);

    int internal = 0;
    int external = 0;

    for (int i = 0; i < numXYZ; i++) {
      char c[10];
      retval = fscanf(xyzFile, "%s", c);
      retval = fscanf(xyzFile, "%lf", &x);
      retval = fscanf(xyzFile, "%lf", &y);
      retval = fscanf(xyzFile, "%lf", &z);
      retval = fscanf(xyzFile, "%lf", &q);

      DPG::Point *p = new DPG::Point(x, y, z);

      string intext(c);

      if (intext.compare("I") == 0) {
        pg->addPoint(p);
        internal++;
      } else if (intext.compare("E") == 0 && receptor) {
        skinAtomList.push_back(p);
        external++;
      }
    }

    // printing internal and external atoms from pqr

    int numAtoms = m_AtomList.size();

    vector<DPG::Point *> results;
    vector<int> externalIndices;

    double close = 0.001;

    for (int i = 0; i < numAtoms; i++) {
      DPG::Point *p = new DPG::Point(m_AtomList[i]->m_Position[0],
                                     m_AtomList[i]->m_Position[1],
                                     m_AtomList[i]->m_Position[2]);

      bool isInternal = false;

      results.clear();
      results = pg->range(p, 1.0);

      int size = results.size();
      if (size > 0) {
        for (int j = 0; j < size; j++) {
          double X, Y, Z;
          X = results[j]->getX();
          Y = results[j]->getY();
          Z = results[j]->getZ();

          if (diff(p->getX(), X) < close && diff(p->getY(), Y) < close &&
              diff(p->getZ(), Z) < close) {
            isInternal =
                true; // found corresponding atoms in internal atoms list
            printAtom(i, isInternal);
          }
        }
      }

      if (!isInternal)
        externalIndices.push_back(i);
    }

    // If it is ligand, adding skin atoms
    if (!receptor) {
      int numSkinAtoms = externalIndices.size();
      for (int i = 0; i < numSkinAtoms; i++) {
        printAtom(externalIndices[i], false);
      }
    }

    // If it is receptor, adding grown skin atoms

    if (receptor) {
      numres = m_AtomList[numAtoms - 1]->getResidueNum() + 1;
      int numSkinAtoms = skinAtomList.size();

      for (int i = 0; i < numSkinAtoms; i++) {
        printSkinAtom(skinAtomList[i]->getX(), skinAtomList[i]->getY(),
                      skinAtomList[i]->getZ(), numAtoms + i);
      }
    }

    fclose(f2dFile);
    fclose(newpqrFile);
  }
}
