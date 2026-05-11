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
#ifndef RAW_SLICE_FAC_H
#define RAW_SLICE_FAC_H

#include <Utility/utility.h>
#include <C2C/BufferedIO.h>
#include <C2C/RawReader.h>
#include <C2C/SliceCache.h>
#include <C2C/SliceFactory.h>

// This class reads slices from a rawiv file
template <class T> class RawSliceFactory : virtual public SliceFactory<T> {
public:
  RawSliceFactory(const char *rawfile);

  virtual ~RawSliceFactory();

  virtual Slice<T> *getNextSlice();

  virtual void reset() { next = 0; }

  virtual int getNumSlices() { return m_dim[2]; }

  virtual void getDimension(int dim[3], float orig[3], float span[3]) {
    for (int i = 0; i < 3; i++) {
      dim[i] = m_dim[i];
      orig[i] = m_orig[i];
      span[i] = m_span[i];
    }
  }

  virtual int current() { return next; }

private:
  void readRawHeader();

  void loadInitPages();

  int m_dim[3];
  float m_orig[3], m_span[3];
  int next; // the number of the next available slice
  SliceCache<T> *cache;
  BufferedIO *m_io;
  RawReader<T> *reader;
};

// Alas, template functions must stay in the header.

template <class T> RawSliceFactory<T>::RawSliceFactory(const char *rawfile) {
  m_io = new BufferedIO(rawfile, DiskIO::READ);
  m_io->open();
  readRawHeader();
  reader = new RawReader<T>(m_io, m_dim);
  cache = new SliceCache<T>(m_dim[0] * m_dim[1] * m_dim[2] * sizeof(T));
  next = 0;
  loadInitPages();
}

template <class T> RawSliceFactory<T>::~RawSliceFactory() {
  delete cache;
  delete reader;
  m_io->close();
  delete m_io;
}

template <class T> Slice<T> *RawSliceFactory<T>::getNextSlice() {
  if (cache->contains(next)) {
    Slice<T> *slc = cache->get(next);
    next++;
    return slc;
  }
  return reader->getSlice(next++);
}

template <class T> void RawSliceFactory<T>::readRawHeader() {
  float extent[6]; // seems not useful anywhere
  int nv, nc;      // may not be correct anyway
  // retrieve header info
  m_io->get(extent, 6);
  m_io->get(&nv, 1);
  m_io->get(&nc, 1);
  m_io->get(m_dim, 3);
  m_io->get(m_orig, 3);
  m_io->get(m_span, 3);
}

template <class T> void RawSliceFactory<T>::loadInitPages() {
  int sz = m_dim[0] * m_dim[1] * sizeof(T);
  int np = (cache->capacity() < sz) ? 1 : cache->capacity() / sz;
  for (int i = 0; i < np && i < m_dim[2]; i++) {
    Slice<T> *slc = reader->getSlice(i);
    cache->add(i, slc);
  }
}
#endif
