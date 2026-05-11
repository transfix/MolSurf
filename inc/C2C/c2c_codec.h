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
#ifndef C2C_CODEC_H
#define C2C_CODEC_H

#include <C2C/ContourGeom.h>
#include <C2C/arithlib/libdjvu++/ByteStream.h>

// decoding
#include <C2C/DiskIO.h>
#include <C2C/FileC2CBuffer.h>
#include <C2C/StreamC2CBuffer.h>
#include <C2C/Decoder.h>
// encoding
#include <C2C/RawSliceFactory.h>
#include <C2C/RawVSliceFactory.h>
#include <C2C/BlockFactory.h>
#include <C2C/CompCon.h>

ContourGeom *decodeC2CFile(const char *fileName, bool &color);
ContourGeom *decodeC2CBuffer(void *data, int size, unsigned char type,
                             bool &color);
void encodeC2CFile(const char *inFile, const char *outFile, unsigned char type,
                   float isoval);
void writeC2CFile(void *data, unsigned char *red, unsigned char *green,
                  unsigned char *blue, unsigned char type, const char *outFile,
                  float isoval, int dim[3], float orig[3], float span[3]);
ByteStream *encodeC2CBuffer(void *data, unsigned char *red,
                            unsigned char *green, unsigned char *blue,
                            unsigned char type, float isoval, int dim[3],
                            float orig[3], float span[3]);

#endif
