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
#ifndef FILE_C2C_BUFFER_H
#define FILE_C2C_BUFFER_H

#include <Utility/utility.h>
#include <C2C/BufferedIO.h>
#include <C2C/C2CBuffer.h>
#include <C2C/DiskIO.h>
#include <C2C/util.h>

template <class T> class FileC2CBuffer : public C2CBuffer<T> {

public:
  /// Construct a buffer from file IO
  FileC2CBuffer(DiskIO *io);

  FileC2CBuffer(const char *fname);

  ~FileC2CBuffer();

  virtual int getSlice(T *data, unsigned char *colors, unsigned char *bits);

  virtual int getLayer(Cell *&cells);

  virtual void getDimension(int dim[3], float orig[3], float span[3]);

  virtual int currentSlice() { return nslice; }

  virtual int currentLayer() { return nlayer; }

  virtual float getIsovalue() { return val; }

private:
  inline void readHeader();

  DiskIO *m_io;
  int m_dim[3];
  float m_orig[3], m_span[3];
  int nlayer, nslice;
  float val;
  char *in_file;
};

template <class T> FileC2CBuffer<T>::FileC2CBuffer(DiskIO *io) : m_io(io) {
  in_file = NULL;
  readHeader();
}

template <class T> FileC2CBuffer<T>::FileC2CBuffer(const char *fname) {
#ifdef WIN32
  in_file = _strdup(fname);
#else
  in_file = strdup(fname);
#endif
  m_io = new BufferedIO(fname, DiskIO::READ);
  m_io->open();
  // read the first byte of data type
  // Constructor from DiskIO* doesn't read this byte
  // XXX:prok - This will cause some serious trouble for color c2c files.
  // I don't know of any code that uses this constructor.
  unsigned char t;
  m_io->get(&t, 1);
  readHeader();
}

template <class T> inline void FileC2CBuffer<T>::readHeader() {
  m_io->get(m_orig, 3);
  m_io->get(m_span, 3);
  m_io->get(m_dim, 3);
  m_io->get(&val, 1);
  nlayer = 0;
  nslice = 0;
}

template <class T> FileC2CBuffer<T>::~FileC2CBuffer() {
  if (in_file != NULL) {
    free(in_file);
    m_io->close();
    delete m_io;
  }
}

template <class T>
void FileC2CBuffer<T>::getDimension(int dim[3], float orig[3], float span[3]) {
  for (int i = 0; i < 3; i++) {
    dim[i] = m_dim[i];
    orig[i] = m_orig[i];
    span[i] = m_span[i];
  }
}

template <class T>
int FileC2CBuffer<T>::getSlice(T *data, unsigned char *colors,
                               unsigned char *bits) {
  memset(data, 0, sizeof(T) * m_dim[0] * m_dim[1]);
  memset(bits, 0, m_dim[0] * m_dim[1]);
  if (colors) {
    memset(colors, 0, 3 * m_dim[0] * m_dim[1]);
  }
  int nv;
  m_io->get(&nv, 1);
#ifdef _DEBUG
  printf("num of verts = %d\n", nv);
#endif
  if (nv == 0) {
    return 0;
  }
  // read the bitmap
  BitBuffer *bbuf = readBitBuffer(m_io);
#ifdef ZP_CODEC
  bbuf->zp_decode();
#else
  bbuf->arith_decode();
#endif
  u_char(*vtrs)[2];
  vtrs = (u_char(*)[2])malloc(sizeof(u_char[2]) * nv);
  int count = 0;
  for (int j = 0; j < m_dim[1]; j++) {
    for (int i = 0; i < m_dim[0]; i++) {
      BIT bit = bbuf->get_a_bit();
      if (bit == 1) {
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
  BitBuffer *vbuf = readBitBuffer(m_io);
  T *tmp = new T[nv];
  decode_vals(vbuf, tmp, nv, val);
  for (int i = 0; i < nv; i++) {
    int n = vtrs[i][0] + vtrs[i][1] * m_dim[0];
    data[n] = tmp[i];
    bits[n] = 1;
  }
  delete[] tmp;
  delete vbuf;
  // read color values (if they're there)
  if (colors) {
    unsigned char *ctmp = new unsigned char[3 * nv];
    m_io->get(ctmp, 3 * nv);
    for (int i = 0; i < nv; i++) {
      int n = vtrs[i][0] + vtrs[i][1] * m_dim[0];
      colors[n * 3 + 0] = ctmp[i * 3 + 0];
      colors[n * 3 + 1] = ctmp[i * 3 + 1];
      colors[n * 3 + 2] = ctmp[i * 3 + 2];
    }
    delete[] ctmp;
  }
  free(vtrs);
  return count;
}

template <class T> int FileC2CBuffer<T>::getLayer(Cell *&cells) {
  int nc;
  m_io->get(&nc, 1);
#ifdef _DEBUG
  printf("num of cells = %d\n", nc);
#endif
  cells = new Cell[nc];
  if (nc == 0) {
    return nc;
  }
  int count = 0;
  BitBuffer *bbuf = readBitBuffer(m_io);
#ifdef ZP_CODEC
  bbuf->zp_decode();
#else
  bbuf->arith_decode();
#endif
  for (int j = 0; j < m_dim[1] - 1; j++) {
    for (int i = 0; i < m_dim[0] - 1; i++) {
      BIT bit = bbuf->get_a_bit();
      if (bit == 1) {
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
