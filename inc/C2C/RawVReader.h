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
// rawvreader.h: this class reads slices of a (4 variable RGBA) rawV file
// (based on rawreader.h)

#ifndef C2C_RAWV_READER_H
#define C2C_RAWV_READER_H

#include <Utility/utility.h>
#include <C2C/DiskIO.h>
#include <C2C/SliceCache.h>
#include <C2C/Slice.h>

template <class T> class RawVReader {

public:
  RawVReader(DiskIO *io, int dim[3], unsigned char types[4], int timesteps,
             int datastart);
  ~RawVReader() {}

  // get the nth slice
  // return NULL if n is greater than the last slice number
  Slice<T> *getSlice(int n);

  /// slice size
  int size() const;

private:
  int typesize(unsigned char type) {
    switch (type) {
    case 1:
      return 1; // unsigned char
    case 2:
      return 2; // unsigned short
    case 3:
      return 4; // unsigned int/long
    case 4:
      return 4; // float
    case 5:
      return 8; // double
    }
    return -1; // we really don't want to get here
  }

  DiskIO *m_io;
  int m_dim[3];
  int m_varstarts[4];
  unsigned char m_types[4];
  int next; // current slice
};

template <class T>
RawVReader<T>::RawVReader(DiskIO *io, int dim[3], unsigned char types[4],
                          int timesteps, int datastart)
    : m_io(io) {
  int i;
  for (i = 0; i < 3; i++) {
    m_dim[i] = dim[i];
  }
  for (i = 0; i < 4; i++) {
    m_types[i] = types[i];
  }
  for (i = 0; i < 4; i++) {
    if (i == 0) {
      m_varstarts[0] = datastart;
    } else
      m_varstarts[i] = m_varstarts[i - 1] + dim[0] * dim[1] * dim[2] *
                                                typesize(types[i - 1]) *
                                                timesteps;
  }
  next = 0;
}

template <class T> Slice<T> *RawVReader<T>::getSlice(int n) {
  if (n >= m_dim[2]) {
    return NULL;
  }
  next = n;
  T *data = new T[m_dim[0] * m_dim[1]];
  unsigned char *red = new unsigned char[m_dim[0] * m_dim[1]],
                *green = new unsigned char[m_dim[0] * m_dim[1]],
                *blue = new unsigned char[m_dim[0] * m_dim[1]];
  // seek to the start of the red variable
  m_io->seek(m_varstarts[0] + (n * m_dim[0] * m_dim[1]), DiskIO::FROM_START);
  // read the red data
  m_io->get(red, m_dim[0] * m_dim[1]);
  // seek to the start of the red variable
  m_io->seek(m_varstarts[1] + (n * m_dim[0] * m_dim[1]), DiskIO::FROM_START);
  // read the green data
  m_io->get(green, m_dim[0] * m_dim[1]);
  // seek to the start of the red variable
  m_io->seek(m_varstarts[2] + (n * m_dim[0] * m_dim[1]), DiskIO::FROM_START);
  // read the blue data
  m_io->get(blue, m_dim[0] * m_dim[1]);
  // seek to the start of the data variable
  m_io->seek(m_varstarts[3] + (n * m_dim[0] * m_dim[1] * typesize(m_types[3])),
             DiskIO::FROM_START);
  // read the data
  m_io->get(data, m_dim[0] * m_dim[1]);
  Slice<T> *slc = new Slice<T>(data, red, green, blue, m_dim);
  next++;
  delete[] data;
  delete[] red;
  delete[] green;
  delete[] blue;
  return slc;
}

template <class T> int RawVReader<T>::size() const {
  int ret = 0;
  for (int i = 0; i < 4; i++) {
    switch (m_types[i]) {
    case 1:
      // unsigned char
      // assumed to be 1 byte
      ret += m_dim[0] * m_dim[1] * sizeof(unsigned char);
      break;
    case 2:
      // unsigned short
      // assumed to be 2 bytes
      ret += m_dim[0] * m_dim[1] * sizeof(unsigned short);
      break;
    case 3:
      // unsigned int/long
      // assumed to be 4 bytes
      ret += m_dim[0] * m_dim[1] * sizeof(unsigned int);
      break;
    case 4:
      // float
      // assumed to be 4 bytes
      ret += m_dim[0] * m_dim[1] * sizeof(float);
      break;
    case 5:
      // double
      // assumed to be 8 bytes
      ret += m_dim[0] * m_dim[1] * sizeof(double);
      break;
    default:
      break;
    }
  }
  return ret;
}

#endif
