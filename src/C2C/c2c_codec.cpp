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
#include <C2C/c2c_codec.h>

ContourGeom* decodeC2CFile(const char* fileName, bool& color)
{
	C2CBuffer<unsigned char> *cbuf;
	C2CBuffer<unsigned short> *sbuf;
	C2CBuffer<float> *fbuf;
	Decoder<float> *fdec;
	Decoder<unsigned short> *sdec;
	Decoder<unsigned char> *cdec;
	DiskIO* input = new BufferedIO(fileName);
	ContourGeom* con = NULL;
	int datatype;
	unsigned char c;
	// read the datatype
	input->open();
	input->get(&c, 1);
	//printf("decoding c2c data...\n");
	datatype = c;
	// tell the caller if there is color or not
	if(datatype > 2)
	{
		color = true;
	}
	else
	{
		color = false;
	}
	// decode the file
	switch(datatype)
	{
		case 0:
		case 3:
			{
				cbuf = new FileC2CBuffer<unsigned char>(input);
				cdec = new Decoder<unsigned char>(cbuf, color);
				con = cdec->constructCon();
				delete cdec;
				delete cbuf;
				break;
			}
		case 1:
		case 4:
			{
				sbuf = new FileC2CBuffer<unsigned short>(input);
				sdec = new Decoder<unsigned short>(sbuf, color);
				con = sdec->constructCon();
				delete sdec;
				delete sbuf;
				break;
			}
		case 2:
		case 5:
			{
				fbuf = new FileC2CBuffer<float>(input);
				fdec = new Decoder<float>(fbuf, color);
				con = fdec->constructCon();
				delete fdec;
				delete fbuf;
				break;
			}
		default:
			break;
	}
	input->close();
	delete input;
	return con;
}

ContourGeom* decodeC2CBuffer(void* data, int size, unsigned char type, bool& color)
{
	MemoryByteStream* mem = new MemoryByteStream(data, size);
	C2CBuffer<unsigned char> *cbuf;
	C2CBuffer<unsigned short> *sbuf;
	C2CBuffer<float> *fbuf;
	Decoder<unsigned char> *cdec;
	Decoder<unsigned short> *sdec;
	Decoder<float> *fdec;
	ContourGeom* con=NULL;
	// tell the caller if there is color or not
	if(type > 2)
	{
		color = true;
	}
	else
	{
		color = false;
	}
	// buffer -> contour
	switch((int)type)
	{
		case 0:
		case 3:
			{
				cbuf = new StreamC2CBuffer<unsigned char>(mem);
				cdec = new Decoder<unsigned char>(cbuf, color);
				con = cdec->constructCon();
				delete cdec;
				delete cbuf;
				break;
			}
		case 1:
		case 4:
			{
				sbuf = new StreamC2CBuffer<unsigned short>(mem);
				sdec = new Decoder<unsigned short>(sbuf, color);
				con = sdec->constructCon();
				delete sdec;
				delete sbuf;
				break;
			}
		case 2:
		case 5:
			{
				fbuf = new StreamC2CBuffer<float>(mem);
				fdec = new Decoder<float>(fbuf, color);
				con = fdec->constructCon();
				delete fdec;
				delete fbuf;
				break;
			}
		default:
			printf("unknown datatype %d\n", type);
			delete mem; // mem is not deleted for us in this case
			break;
	}
	// mem is deleted by the StreamC2CBuffer destructor
	mem = NULL;
	return con;
}

// type
// 0 - unsigned char rawiv
// 1 - unsigned short rawiv
// 2 - float rawiv
// 3 - unsigned char rawv
// 4 - unsigned short rawv
// 5 - float rawv
void encodeC2CFile(const char* inFile, const char* outFile, unsigned char type, float isoval)
{
	SliceFactory<unsigned char> *cfac;
	SliceFactory<unsigned short> *sfac;
	SliceFactory<float> *ffac;
	CompCon<unsigned char> *ccon;
	CompCon<unsigned short> *scon;
	CompCon<float> *fcon;
	switch((int)type)
	{
			// colorless models
		case 0:
			{
				cfac = new RawSliceFactory<unsigned char>(inFile);
				ccon = new CompCon<unsigned char>(cfac);
				ccon->setOutputFile(outFile, 0);
				ccon->marchingCubes(isoval);
				delete ccon;
				delete cfac;
				break;
			}
		case 1:
			{
				sfac = new RawSliceFactory<unsigned short>(inFile);
				scon = new CompCon<unsigned short>(sfac);
				scon->setOutputFile(outFile, 1);
				scon->marchingCubes(isoval);
				delete scon;
				delete sfac;
				break;
			}
		case 2:
			{
				ffac = new RawSliceFactory<float>(inFile);
				fcon = new CompCon<float>(ffac);
				fcon->setOutputFile(outFile, 2);
				fcon->marchingCubes(isoval);
				delete fcon;
				delete ffac;
				break;
			}
			// color models
		case 3:
			{
				cfac = new RawVSliceFactory<unsigned char>(inFile);
				ccon = new CompCon<unsigned char>(cfac);
				ccon->setOutputFile(outFile, 3);
				ccon->marchingCubes(isoval);
				delete ccon;
				delete cfac;
				break;
			}
		case 4:
			{
				sfac = new RawVSliceFactory<unsigned short>(inFile);
				scon = new CompCon<unsigned short>(sfac);
				scon->setOutputFile(outFile, 4);
				scon->marchingCubes(isoval);
				delete scon;
				delete sfac;
				break;
			}
		case 5:
			{
				ffac = new RawVSliceFactory<float>(inFile);
				fcon = new CompCon<float>(ffac);
				fcon->setOutputFile(outFile, 5);
				fcon->marchingCubes(isoval);
				delete fcon;
				delete ffac;
				break;
			}
	}
}

void writeC2CFile(void* data, unsigned char* red, unsigned char* green, unsigned char* blue, unsigned char type, const char* outFile, float isoval, int dim[3], float orig[3], float span[3])
{
	BlockFactory<unsigned char> *cfac;
	BlockFactory<unsigned short> *sfac;
	BlockFactory<float> *ffac;
	CompCon<unsigned char> *ccon;
	CompCon<unsigned short> *scon;
	CompCon<float> *fcon;
	switch((int)type)
	{
		case 0:
		case 3:
			{
				cfac = new BlockFactory<unsigned char>((unsigned char*)data,red,green,blue,dim,orig,span);
				ccon = new CompCon<unsigned char>(cfac);
				ccon->setOutputFile(outFile, type);
				ccon->marchingCubes(isoval);
				delete ccon;
				delete cfac;
				break;
			}
		case 1:
		case 4:
			{
				sfac = new BlockFactory<unsigned short>((unsigned short*)data,red,green,blue, dim,orig,span);
				scon = new CompCon<unsigned short>(sfac);
				scon->setOutputFile(outFile, type);
				scon->marchingCubes(isoval);
				delete scon;
				delete sfac;
				break;
			}
		case 2:
		case 5:
			{
				ffac = new BlockFactory<float>((float*)data,red,green,blue,dim,orig,span);
				fcon = new CompCon<float>(ffac);
				fcon->setOutputFile(outFile, type);
				fcon->marchingCubes(isoval);
				delete fcon;
				delete ffac;
				break;
			}
	}
}

ByteStream* encodeC2CBuffer(void* data, unsigned char* red, unsigned char* green, unsigned char* blue, unsigned char type, float isoval, int dim[3], float orig[3], float span[3])
{
	MemoryByteStream* ret = new MemoryByteStream();
	BlockFactory<unsigned char> *cfac;
	BlockFactory<unsigned short> *sfac;
	BlockFactory<float> *ffac;
	CompCon<unsigned char> *ccon;
	CompCon<unsigned short> *scon;
	CompCon<float> *fcon;
	switch((int)type)
	{
		case 0:
		case 3:
			{
				cfac = new BlockFactory<unsigned char>((unsigned char*)data,red,green,blue, dim,orig,span);
				ccon = new CompCon<unsigned char>(cfac);
				ccon->setOutputStream(ret);
				ccon->marchingCubes(isoval);
				delete ccon;
				delete cfac;
				break;
			}
		case 1:
		case 4:
			{
				sfac = new BlockFactory<unsigned short>((unsigned short*)data,red,green,blue, dim,orig,span);
				scon = new CompCon<unsigned short>(sfac);
				scon->setOutputStream(ret);
				scon->marchingCubes(isoval);
				delete scon;
				delete sfac;
				break;
			}
		case 2:
		case 5:
			{
				ffac = new BlockFactory<float>((float*)data,red,green,blue,dim,orig,span);
				fcon = new CompCon<float>(ffac);
				fcon->setOutputStream(ret);
				fcon->marchingCubes(isoval);
				delete fcon;
				delete ffac;
				break;
			}
		default:
			break;
	}
	return ret;
}
