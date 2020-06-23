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
#ifndef CCV_PQR_FILE_H
#define CCV_PQR_FILE_H

#include <GOAFileTypes/GOAFileType.h>

namespace PDBParser
{
	class GroupOfAtoms;
};

class PQRFile : public GOAFileType
{
	public:
		virtual ~PQRFile();
		virtual PDBParser::GroupOfAtoms* loadFile(const string& fileName, bool deleteWater = false);
		virtual bool checkType(const string& fileName);
		virtual bool saveFile(PDBParser::GroupOfAtoms* groupOfAtoms, const string& fileName, unsigned int level=0, CCVOpenGLMath::Matrix* transformation = 0);
		virtual string extension()
		{
			return "pqr";
		};
		virtual string filter()
		{
			return "PQR files (*.pqr)";
		};
		static PQRFile ms_PQRFileRepresentative;
		static GOAFileType* getRepresentative();
	protected:
		PQRFile();
		void setChainColors(PDBParser::GroupOfAtoms* chain, int n);
	private:
		void getString(char* line, char* str, int pos, int len);
		void getChar(char* line, char* str, int pos);
		void getInt(char* line, int* i, int pos, int len);
		void getFloat(char* line, float* d, int pos, int len);
};

#endif
