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
#ifndef CCV_PDB1_FILE_H
#define CCV_PDB1_FILE_H

#include <GOAFileTypes/GOAFileType.h>

namespace PDB1Parser {
class GroupOfAtoms;
};

class PDB1File : public GOAFileType {
public:
  virtual ~PDB1File();
  virtual PDBParser::GroupOfAtoms *loadFile(const string &fileName,
                                            bool deleteWater = false);
  virtual bool checkType(const string &fileName);
  virtual bool saveFile(PDBParser::GroupOfAtoms *groupOfAtoms,
                        const string &fileName, unsigned int level = 0,
                        CCVOpenGLMath::Matrix *transformation = 0);
  virtual string extension() { return "pdb1"; };
  virtual string filter() { return "PDB1 files (*.pdb1)"; };
  static PDB1File ms_PDB1FileRepresentative;
  static GOAFileType *getRepresentative();

protected:
  PDB1File();
};

#endif
