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
#include <VolumeFileTypes/RawVFile.h>

using namespace CVCUtility;

RawVFile RawVFile::ms_RawVFileRepresentative;

RawVFile::RawVFile() {}

RawVFile::~RawVFile() {}

SimpleVolumeData *RawVFile::loadFile(const string &fileName) {
  SimpleVolumeData *simpleVolumeData = new SimpleVolumeData(128, 128, 128);
  FILE *fp = FOPEN(fileName.c_str(), "rb");
  unsigned int i, magic, dims[3], numTimeSteps, numVariables;
  float minExt[4], maxExt[4];
  // check to make sure the file exists
  if (!fp) {
    printf("Error: could not open file %s\n", fileName.c_str());
    delete simpleVolumeData;
    simpleVolumeData = 0;
    return 0;
  }
  // read the header
  // read the magic number
  freadSafely(&magic, 1, sizeof(unsigned int), fp);
  if (isLittleEndian()) {
    swapByteOrder(&magic, 1);
  }
  if (magic != 0xBAADBEEF) {
    // doh! this isn't a RawV file
    printf("Error: file format is not recognized. (expecting RawV)\n");
    fclose(fp);
    return 0;
  }
  // the dimensions, # of timesteps, # of variables, and min/max extents
  freadSafely(dims, 3, sizeof(unsigned int), fp);
  freadSafely(&numTimeSteps, 1, sizeof(unsigned int), fp);
  freadSafely(&numVariables, 1, sizeof(unsigned int), fp);
  freadSafely(minExt, 4, sizeof(float), fp);
  freadSafely(maxExt, 4, sizeof(float), fp);
  if (isLittleEndian()) {
    swapByteOrder(dims, 3);
    swapByteOrder(&numTimeSteps, 1);
    swapByteOrder(&numVariables, 1);
    swapByteOrder(minExt, 4);
    swapByteOrder(maxExt, 4);
  }
  // call a bunch of set...() functions
  // setNumberOfVariables() has the convenient side effect of cleaning up
  // any data that this instance might already be holding. (effectively wiping
  // the slate clean)
  simpleVolumeData->setNumberOfVariables(numVariables);
  // the rest ...
  simpleVolumeData->setDimensions(dims);
  simpleVolumeData->setMinExtent(minExt);
  simpleVolumeData->setMaxExtent(maxExt);
  // finish reading the header (variable names and types)
  for (i = 0; i < numVariables; i++) {
    char name[64];
    unsigned char type;
    // read the type and the name of a variable
    freadSafely(&type, 1, sizeof(unsigned char), fp);
    freadSafely(name, 64, sizeof(char), fp);
    name[63] = '\0'; // just in case
    // call set...() functions
    switch (type) {
    case 1:
      simpleVolumeData->setType(i, SimpleVolumeData::UCHAR);
      break;
    case 2:
      simpleVolumeData->setType(i, SimpleVolumeData::USHORT);
      break;
    case 3:
      simpleVolumeData->setType(i, SimpleVolumeData::ULONG);
      break;
    case 4:
      simpleVolumeData->setType(i, SimpleVolumeData::FLOAT);
      break;
    case 5:
      simpleVolumeData->setType(i, SimpleVolumeData::DOUBLE);
      break;
    default:
      break;
    }
    simpleVolumeData->setName(i, name);
  }
  // finally, read the data
  for (i = 0; i < numVariables; i++) {
    // now that all the meta-data is filled in, allocate space for the variables
    unsigned char *data = new unsigned char[simpleVolumeData->getDataSize(i)];
    // read the data for this variables first (and only?) timestep
    freadSafely(data, simpleVolumeData->getDataSize(i), 1, fp);
    // swap the byte order if needed
    if (isLittleEndian()) {
      switch (simpleVolumeData->getTypeSize(i)) {
      case 1:
        break;
      case 2:
        swapByteOrder((unsigned short *)data, dims[0] * dims[1] * dims[2]);
        break;
      case 4:
        swapByteOrder((float *)data, dims[0] * dims[1] * dims[2]);
        break;
      case 8:
        swapByteOrder((double *)data, dims[0] * dims[1] * dims[2]);
        break;
      default:
        break;
      }
    }
    simpleVolumeData->setData(i, data);
    // seek past any remaining timesteps
    if (numTimeSteps > 1) {
      fseek(fp, simpleVolumeData->getDataSize(i) * (numTimeSteps - 1),
            SEEK_CUR);
    }
  }
  // close the file
  fclose(fp);
  return simpleVolumeData;
}

bool RawVFile::checkType(const string &fileName) { return false; }

bool RawVFile::saveFile(SimpleVolumeData *simpleVolumeData,
                        const string &fileName, unsigned int variable) {
  // no data? bail.
  if (!simpleVolumeData->getData(variable)) {
    return false;
  }
  // unsupported datatype for the requested variable? bail.
  if (simpleVolumeData->getType(variable) == SimpleVolumeData::ULONG ||
      simpleVolumeData->getType(variable) == SimpleVolumeData::DOUBLE) {
    return false;
  }
  FILE *fp = FOPEN(fileName.c_str(), "wb");
  char *header = NULL;
  unsigned int hsize = 0;
  // failed to open the file? bail.
  if (!fp) {
    return false;
  }
  // make the data big endian
  simpleVolumeData->makeVariablesBigEndian();
  // create the header
  simpleVolumeData->createRawVHeader(&header, &hsize);
  // write the header
  fwrite(header, hsize, 1, fp);
  // write the data
  for (unsigned int i = 0; i < simpleVolumeData->getNumberOfVariables(); i++) {
    fwrite(simpleVolumeData->getData(i), simpleVolumeData->getDataSize(i), 1,
           fp);
  }
  // make the data native endian (poorly named function)
  simpleVolumeData->makeVariablesBigEndian();
  // clean up
  delete[] header;
  fclose(fp);
  return true;
}

VolumeFileType *RawVFile::getRepresentative() {
  return &ms_RawVFileRepresentative;
}
