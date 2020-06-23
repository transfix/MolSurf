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
#ifndef CCV_GOA_FILE_TYPE_H
#define CCV_GOA_FILE_TYPE_H

#include <stdio.h>
#include <string>
#include <UsefulMath/Matrix.h>

using std::string;

namespace PDBParser
{
	class GroupOfAtoms;
};

class GOAFileType
{
	public:
		GOAFileType();
		virtual ~GOAFileType();
		virtual PDBParser::GroupOfAtoms* loadFile(const string& fileName, bool deleteWater = false) = 0;
		virtual bool checkType(const string& fileName) = 0;
		virtual bool saveFile(PDBParser::GroupOfAtoms* groupOfAtoms, const string& fileName, unsigned int level=0, CCVOpenGLMath::Matrix* transformation = 0) = 0;
		virtual string extension() = 0;
		virtual string filter() = 0;
};

#endif
