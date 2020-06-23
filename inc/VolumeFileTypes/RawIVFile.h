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
#ifndef CCV_RAWIV_FILE_H
#define CCV_RAWIV_FILE_H

#include <Utility/utility.h>
#include <VolumeFileTypes/VolumeFileType.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <ByteOrder/ByteSwapping.h>

#if defined(_LARGEFILE_SOURCE)
# define FOPEN fopen64
#else
# define FOPEN fopen
#endif

class RawIVFile
 : public VolumeFileType
{
	public:
		virtual ~RawIVFile();
		virtual SimpleVolumeData* loadFile(const string& fileName);
		virtual bool checkType(const string& fileName);
		virtual bool saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable=0);
		virtual string extension()
		{
			return "rawiv";
		};
		virtual string filter()
		{
			return "RawIV files (*.rawiv)";
		};
		static RawIVFile ms_RawIVFileRepresentative;
		static VolumeFileType* getRepresentative();
	protected:
		RawIVFile();
};

#endif
