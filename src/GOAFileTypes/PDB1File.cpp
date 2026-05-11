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
// #include <Blurmaps/BlurMapsDataManager.h>
#include <GOAFileTypes/PDB1File.h>
#include <math.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/GOAFileIO.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/parserPDBtoGOA.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

PDB1File PDB1File::ms_PDB1FileRepresentative;

PDB1File::PDB1File() {}

PDB1File::~PDB1File() {}

PDBParser::GroupOfAtoms *PDB1File::loadFile(const string &fileName,
                                            bool deleteWater) {
  PDBParser::parserPDBtoGOA *p = new PDBParser::parserPDBtoGOA();
  PDBParser::GroupOfAtoms *molecule =
      p->parsePDBModels(fileName.c_str(), deleteWater);
  delete p;
  return molecule;
}

bool PDB1File::checkType(const string &fileName) { return false; }

bool PDB1File::saveFile(PDBParser::GroupOfAtoms *molecule,
                        const string &fileName, unsigned int level,
                        CCVOpenGLMath::Matrix *transformation) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  bool ret = PDBParser::writeGOA2PDB(fp, molecule, level, 0);
  fclose(fp);
  return ret;
}

GOAFileType *PDB1File::getRepresentative() {
  return &ms_PDB1FileRepresentative;
}
