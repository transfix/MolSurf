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
#ifndef CCV_DX_FILE_H
#define CCV_DX_FILE_H

#include <Utility/utility.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <ByteOrder/ByteSwapping.h>
#include <VolumeFileTypes/VolumeFileType.h>

class DXFile : public VolumeFileType {
public:
  virtual ~DXFile();
  virtual SimpleVolumeData *loadFile(const string &fileName);
  virtual bool checkType(const string &fileName);
  virtual bool saveFile(SimpleVolumeData *simpleVolumeData,
                        const string &fileName, unsigned int variable = 0);
  virtual string extension() { return "dx"; };
  virtual string filter() { return "DX files (*.dx)"; };
  static DXFile ms_DXFileRepresentative;
  static VolumeFileType *getRepresentative();

protected:
  DXFile();
  bool isCommentOrEmpty(const char *line);
  void tryToGetDimensions(bool *dimsSet, unsigned int *dims, const char *line);
  void tryToGetOrigin(bool *originSet, float *origin, const char *line);
  void tryToGetSpans(int *deltasFound, float *span, const char *line);
  void tryToGetDataHeader(bool *dataHeaderFound, unsigned int *dims,
                          const char *line);
};

#endif
