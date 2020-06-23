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
#ifndef STREAM_C2C_BUFFER_H
#define STREAM_C2C_BUFFER_H

#include <Utility/utility.h>
#include <C2C/arithlib/encode.h>
#include <C2C/arithlib/libdjvu++/ByteStream.h>
#include <C2C/C2CBuffer.h>
#include <C2C/util.h>

template <class T>
class StreamC2CBuffer : public C2CBuffer<T>
{

	public:
		// Construct a buffer from a ByteStream
		StreamC2CBuffer(ByteStream* io);

		~StreamC2CBuffer();

		virtual int getSlice(T* data, unsigned char* colors, unsigned char* bits);

		virtual int getLayer(Cell* &cells);

		virtual void getDimension(int dim[3], float orig[3], float span[3]);

		virtual int currentSlice()
		{
			return nslice;
		}

		virtual int currentLayer()
		{
			return nlayer;
		}

		virtual float getIsovalue()
		{
			return val;
		}
	private:
		inline void readHeader();

		ByteStream* m_io;
		int m_dim[3];
		float m_orig[3], m_span[3];
		int nlayer, nslice;
		float val;
};

template <class T>
StreamC2CBuffer<T>::StreamC2CBuffer(ByteStream* io) : m_io(io)
{
	readHeader();
}

template <class T>
StreamC2CBuffer<T>::~StreamC2CBuffer()
{
	delete m_io;
}

template <class T>
inline void StreamC2CBuffer<T>::readHeader()
{
	m_io->read(m_orig, sizeof(float)*3);
	m_io->read(m_span, sizeof(float)*3);
	m_dim[0] = m_io->read32();
	m_dim[1] = m_io->read32();
	m_dim[2] = m_io->read32();
	m_io->read(&val, sizeof(float));
	nlayer = 0;
	nslice = 0;
}


template <class T>
void  StreamC2CBuffer<T>::getDimension(int dim[3], float orig[3], float span[3])
{
	for(int i = 0; i < 3; i++)
	{
		dim[i] = m_dim[i];
		orig[i] = m_orig[i];
		span[i] = m_span[i];
	}
}

template <class T>
int StreamC2CBuffer<T>::getSlice(T* data, unsigned char* colors,
								 unsigned char* bits)
{
	memset(data, 0, sizeof(T)*m_dim[0]*m_dim[1]);
	memset(bits, 0, m_dim[0]*m_dim[1]);
	if(colors)
	{
		memset(colors, 0, 3*m_dim[0]*m_dim[1]);
	}
	int nv;
	nv = m_io->read32();
#ifdef _DEBUG
	printf("num of verts = %d\n", nv);
#endif
	if(nv == 0)
	{
		return 0;
	}
	// read the bitmap
	BitBuffer* bbuf = readBitBuffer(m_io);
#ifdef ZP_CODEC
	bbuf->zp_decode();
#else
	bbuf->arith_decode();
#endif
	u_char(*vtrs)[2];
	vtrs = (u_char(*)[2])malloc(sizeof(u_char[2])*nv);
	int count = 0;
	for(int j = 0; j < m_dim[1]; j++)
	{
		for(int i = 0; i < m_dim[0]; i++)
		{
			BIT bit = bbuf->get_a_bit();
			if(bit == 1)
			{
				vtrs[count][0] = i;
				vtrs[count][1] = j;
				count++;
			}
		}
	}
	delete bbuf;
#ifdef _DEBUG
	assert(count == nv);
#endif
	// read function values
	BitBuffer* vbuf = readBitBuffer(m_io);
	T* tmp = new T[nv];
	decode_vals(vbuf, tmp, nv, val);
	for(int i = 0; i < nv; i++)
	{
		int n = vtrs[i][0] + vtrs[i][1]*m_dim[0];
		data[n] = tmp[i];
		bits[n] = 1;
	}
	delete[] tmp;
	delete vbuf;
	// read the colors
	if(colors)
	{
		unsigned char* ctmp = new unsigned char [3*nv];
		m_io->read(ctmp, 3*nv);
		for(int i=0; i < nv; i++)
		{
			int n = vtrs[i][0] + vtrs[i][1]*m_dim[0];
			colors[n*3+0] = ctmp[i*3+0];
			colors[n*3+1] = ctmp[i*3+1];
			colors[n*3+2] = ctmp[i*3+2];
		}
		delete [] ctmp;
	}
	free(vtrs);
	return count;
}

template <class T>
int StreamC2CBuffer<T>::getLayer(Cell* &cells)
{
	int nc;
	nc = m_io->read32();
#ifdef _DEBUG
	printf("num of cells = %d\n", nc);
#endif
	cells = new Cell[nc];
	if(nc == 0)
	{
		return nc;
	}
	int count = 0;
	BitBuffer* bbuf = readBitBuffer(m_io);
#ifdef ZP_CODEC
	bbuf->zp_decode();
#else
	bbuf->arith_decode();
#endif
	for(int j = 0; j < m_dim[1]-1; j++)
	{
		for(int i = 0; i < m_dim[0]-1; i++)
		{
			BIT bit = bbuf->get_a_bit();
			if(bit == 1)
			{
				cells[count].iy = j;
				cells[count].ix = i;
				count++;
			}
		}
	}
#ifdef _DEBUG
	assert(count == nc);
#endif
	delete bbuf;
	return nc;
}
#endif
