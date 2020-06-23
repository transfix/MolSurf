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
#ifndef CCV_VOLUME_LOADER_H
#define CCV_VOLUME_LOADER_H

#include <Utility/utility.h>
#include <VolumeFileTypes/RawIVFile.h>
#include <VolumeFileTypes/RawVFile.h>
#include <VolumeFileTypes/DXFile.h>
#include <VolumeFileTypes/MRCFile.h>
#include <VolumeFileTypes/VolumeFileType.h>

class VolumeFileType;
class SimpleVolumeData;

class VolumeLoader
{
	public:
		VolumeLoader();
		virtual ~VolumeLoader();
		bool saveFile(const string& fileName, const string& selectedFilter, SimpleVolumeData* simpleVolumeData, unsigned int variable=0);
		SimpleVolumeData* loadFile(const string& fileName);
		bool saveFile(const string& fileName, SimpleVolumeData* simpleVolumeData);
		string getLoadFilterString();
		string getSaveFilterString();
		bool isValidExtension(string extension);
	protected:
		string getAllExtensions();
		bool endsWith(string str, string substr);
		SimpleVolumeData* tryAll(const string& fileName);
		void addVolumeFileType(VolumeFileType* type);
		std::map<string, VolumeFileType*> m_ExtensionMap;
		std::map<string, VolumeFileType*> m_FilterMap;
};

#endif
