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
#include <VolumeFileTypes/RawIVFile.h>

using namespace CVCUtility;

RawIVFile RawIVFile::ms_RawIVFileRepresentative;

RawIVFile::RawIVFile()
{
}

RawIVFile::~RawIVFile()
{
}

SimpleVolumeData* RawIVFile::loadFile(const string& fileName)
{
	SimpleVolumeData* simpleVolumeData = new SimpleVolumeData(128, 128, 128);
	FILE* fp=FOPEN(fileName.c_str(), "rb");
	unsigned int dims[3];
	unsigned int temp_numverts, temp_numcells; //for large files, these are meaningless, and anyway not used.
	float minExt[3], maxExt[3], orig[3], span[3];
	// check to make sure the file exists
	if(!fp)
	{
		printf("Error: could not open file %s\n", fileName.c_str());
		delete simpleVolumeData;
		simpleVolumeData = 0;
		return 0;
	}
	// read the header
	freadSafely(minExt, 3, sizeof(float), fp);
	freadSafely(maxExt, 3, sizeof(float), fp);
	freadSafely(&temp_numverts, 1, sizeof(unsigned int), fp);
	freadSafely(&temp_numcells, 1, sizeof(unsigned int), fp);
	freadSafely(dims, 3, sizeof(unsigned int), fp);

	freadSafely(orig, 3, sizeof(float), fp);
	freadSafely(span, 3, sizeof(float), fp);

	if(isLittleEndian())
	{
		swapByteOrder(minExt, 3);
		swapByteOrder(maxExt, 3);
		swapByteOrder(&temp_numverts, 1);
		swapByteOrder(&temp_numcells, 1);
		swapByteOrder(dims, 3);
		swapByteOrder(orig, 3);
		swapByteOrder(span, 3);
	}

	// find out how large the data is
	Q_ULLONG dataStart = ftell(fp), dataSize;
	fseek(fp, 0, SEEK_END);
	dataSize = ftell(fp) - dataStart;
	fseek(fp, dataStart, SEEK_SET);
	// a small sanity check to make sure this file is valid
	int testSize = dims[0]*dims[1]*dims[2];
	if(!testSize)
	{
		error("File " + fileName + " has a 0-sized dimension!");
	}
	if(dataSize % testSize)
	{
		error("Error: rawiv file " + fileName  + " header dimensions don't match file size");
	}

	// call some set...() functions
	simpleVolumeData->setNumberOfVariables(1);
	simpleVolumeData->setDimensions(dims);
	simpleVolumeData->setMinExtent(minExt);
	simpleVolumeData->setMaxExtent(maxExt);
	// figure out the data type
	switch(dataSize / ((Q_ULLONG)dims[0]*dims[1]*dims[2]))
	{
		case 1:
			simpleVolumeData->setType(0, SimpleVolumeData::UCHAR);
			break;
		case 2:
			simpleVolumeData->setType(0, SimpleVolumeData::USHORT);
			break;
		case 4:
			simpleVolumeData->setType(0, SimpleVolumeData::FLOAT);
			break;
		default:
			simpleVolumeData->setType(0, SimpleVolumeData::NO_TYPE);
			break;
	}
	// allocate space for the data
	unsigned char* data = (unsigned char*)malloc((Q_ULLONG)sizeof(unsigned char)*dataSize);
	// read the data
	freadSafely(data, dataSize, 1, fp);
	// swap the byte order if needed
	if(isLittleEndian())
	{
		switch(simpleVolumeData->getTypeSize(0))
		{
			case 1:
				break;
			case 2:
				swapByteOrder((unsigned short*)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			case 4:
				swapByteOrder((float*)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			case 8:
				swapByteOrder((double*)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			default:
				break;
		}
	}
	simpleVolumeData->setData(0, data);
	// close the file
	fclose(fp);
	return simpleVolumeData;
}

bool RawIVFile::checkType(const string& fileName)
{
	return false;
}

bool RawIVFile::saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable)
{
	// no data? bail.
	if(!simpleVolumeData->getData(variable))
	{
		return false;
	}
	// unsupported datatype for the requested variable? bail.
	if(simpleVolumeData->getType(variable) == SimpleVolumeData::ULONG || simpleVolumeData->getType(variable) == SimpleVolumeData::DOUBLE)
	{
		return false;
	}
	FILE* fp = FOPEN(fileName.c_str(), "wb");
	char header[68];
	// failed to open the file? bail.
	if(!fp)
	{
		return false;
	}
	// make the data big endian
	simpleVolumeData->makeVariablesBigEndian();
	// create the header
	simpleVolumeData->createRawIVHeader(header);
	// write the header
	fwrite(header, sizeof(header), 1, fp);
	// write the data
	fwrite(simpleVolumeData->getData(variable), simpleVolumeData->getDataSize(variable), 1, fp);
	// make the data native endian (poorly named function)
	simpleVolumeData->makeVariablesBigEndian();
	// clean up
	fclose(fp);
	return true;
}

VolumeFileType* RawIVFile::getRepresentative()
{
	return &ms_RawIVFileRepresentative;
}
