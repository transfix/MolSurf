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
#ifndef COMP_CON_H
#define COMP_CON_H

#include <Utility/utility.h>
#include <C2C/arithlib/libdjvu++/ByteStream.h>
#include <C2C/BufferedIO.h>
#include <C2C/cubes.h>
#include <C2C/SingleConStep.h>
#include <C2C/SliceFactory.h>
#include <C2C/Slice.h>
#include <C2C/util.h>

class Computable {
public:
  virtual ~Computable() {}

  virtual void marchingCubes(float isoval) = 0;

  virtual void setOutputFile(const char *fname, unsigned char t) = 0;

  virtual void setOutputStream(ByteStream *_stream) = 0;

  virtual void reset() = 0;
};

template <class T> class CompCon : public Computable {
public:
  // Construct a CompCon from a slice factory
  CompCon(SliceFactory<T> *fac);

  virtual ~CompCon();

  void setOutputFile(const char *fname, unsigned char t = 0);

  void setOutputStream(ByteStream *_stream);

  // get one new slice
  bool getNewSlice();

  // march through cells between current two slices
  void MarchingSlice(float isoval);

  // get the vertex values of one cell
  void getCellValues(int i, int j, T data[8]);

  /// mark end vertices for edge e of cell (i, j)
  void markEdge(int i, int j, int e);

  // march through all layers of cells
  void marchingCubes(float isoval);

  // reset the states for next extraction
  void reset();

  // accessors
  Slice<T> *firstSlice() { return first; }

  Slice<T> *secondSlice() { return next; }

  int getNSlices() { return nslices; }

protected:
  void writeHeader(float isoval);

  void encodeSlice(int, u_char *, float); // write out data in encoded format

  void runLength(int, u_char *); // run-length encoding of cell indices

  void en_arithIndex(int, u_char *); // arithmetic encoding cell indices

private:
  int dx, dy, nslices;    // x, y, z dimension of the data
  int count;              // current slice
  Slice<T> *first, *next; // two adjacent slices to compute isosurface
  float orig[3], span[3];
  SliceFactory<T> *m_fac;

  // wflag: false: return a singleconstep, true: write to file or a stream
  bool wflag;
  // out_flag: 0: write to a diskio, 1: write to a stream
  int out_flag;

  SingleConStep<T> *scon;
  ByteStream *stream;
  DiskIO *out;
};

template <class T> CompCon<T>::CompCon(SliceFactory<T> *fac) : m_fac(fac) {
  int dim[3];
  m_fac->getDimension(dim, orig, span);
  dx = dim[0];
  dy = dim[1];
  nslices = dim[2];
  count = 0;
  first = m_fac->getNextSlice();
  count++;
  next = m_fac->getNextSlice();
  count++;
  wflag = false;
  scon = NULL;
  out = NULL;
  stream = NULL;
}

template <class T> CompCon<T>::~CompCon() {
  if (wflag && out != NULL) {
    out->close(false);
    delete out;
  }
  if (first != NULL) {
    delete first;
  }
  if (next != NULL) {
    delete next;
  }
}

template <class T> void CompCon<T>::reset() {
  count = 0;
  if (first != NULL) {
    delete first;
  }
  if (next != NULL) {
    delete next;
  }
  m_fac->reset();
  first = m_fac->getNextSlice();
  count++;
  next = m_fac->getNextSlice();
  count++;
}

template <class T>
void CompCon<T>::setOutputFile(const char *fname, unsigned char t) {
  if (out != NULL) {
    out->close(false);
    delete out;
  }
  out = new BufferedIO(fname, DiskIO::WRITE);
  bool r = out->open();
#ifdef _DEBUG
  if (!r) {
    printf("cannot open %s to write\n", fname);
    exit(1);
  }
#endif
  wflag = true;
  // output the data type as the first byte
  out->put(&t, 1);
  out_flag = 0;
}

template <class T> void CompCon<T>::setOutputStream(ByteStream *_stream) {
  stream = _stream;
  wflag = true;
  out_flag = 1;
}

template <class T> bool CompCon<T>::getNewSlice() {
  // printf("get new slice %d\n", count);
  if (count >= nslices) {
    // fprintf(stderr, "no more slices\n");
    return false;
  }
  if (first != NULL) {
    delete first; // remove the old slice from memory
  }
  first = next; // change the new slice to old
  if (count < nslices) {
    next = m_fac->getNextSlice();
  }
  count++;
  return true;
}

template <class T> void CompCon<T>::getCellValues(int i, int j, T data[8]) {
  int offset = i + j * dx;
  data[0] = first->verts[offset].getValue();
  data[1] = first->verts[offset + 1].getValue();
  data[2] = next->verts[offset + 1].getValue();
  data[3] = next->verts[offset].getValue();
  data[4] = first->verts[offset + dx].getValue();
  data[5] = first->verts[offset + dx + 1].getValue();
  data[6] = next->verts[offset + dx + 1].getValue();
  data[7] = next->verts[offset + dx].getValue();
}

template <class T>
void CompCon<T>::encodeSlice(int count_c, u_char *cells, float isoval) {
#ifdef _DEBUG
  printf("# of intersected cells = %d at layer %d\n", count_c, count - 1);
#endif
  switch (out_flag) {
  case 0:
    first->writeOut(out, isoval);
    out->put(&count_c, 1);
    break;
  case 1:
    first->writeOut(stream, isoval);
    stream->write32(count_c);
    break;
  default:
    return;
  }
  en_arithIndex(count_c, cells);
}

// arithmetic encoding the intersected cells
template <class T> void CompCon<T>::en_arithIndex(int count_c, u_char *cells) {
  if (count_c == 0) {
    return;
  }
  BIT bit = 0;
  BitBuffer *bbuf = new BitBuffer();
  BIT *bits = new BIT[(dy - 1) * (dx - 1)];
#ifdef SQUARE_ZAG
  // arrange bits in continuous order
  for (int j = 0; j < dy - 1; j++) {
    for (int i = 0; i < dx - 1; i++) {
      int n = i + j * (dx - 1);
      if (j % 2 == 0) {
        if (cubeedges[cells[n]][0] > 0) {
          bits[n] = 1;
        } else {
          bits[n] = 0;
        }
      } else {
        int m = (dx - 1) * (j + 1) - i - 1;
        if (cubeedges[cells[m]][0] > 0) {
          bits[n] = 1;
        } else {
          bits[n] = 0;
        }
      }
    }
  }
#elif defined(TRIANGLE_ZAG)
  int i = 0, j = 0, n = 0;
  int dir = 0; // ZAG direction
  while (i < (dx - 1) && j < (dy - 1)) {
    int m = i + j * (dx - 1);
    if (cubeedges[cells[m]][0] > 0) {
      bits[n++] = 1;
    } else {
      bits[n++] = 0;
    }
    if (j == 0 && dir == 0) {
      if (i < dx - 2) {
        i++;
      } else {
        j++;
      }
      dir = (dir + 1) % 2;
      continue;
    }
    if (i == 0 && dir == 1) {
      if (j < dx - 2) {
        j++;
      } else {
        i++;
      }
      dir = (dir + 1) % 2;
      continue;
    }
    if (i == dx - 2 && dir == 0) {
      j++;
      dir = (dir + 1) % 2;
      continue;
    }
    if (j == dx - 2 && dir == 1) {
      i++;
      dir = (dir + 1) % 2;
      continue;
    }
    if (dir == 0) {
      i++;
      j--;
      continue;
    }
    if (dir == 1) {
      i--;
      j++;
      continue;
    }
  }
#else
  for (int j = 0; j < dy - 1; j++) {
    for (int i = 0; i < dx - 1; i++) {
      int n = i + j * (dx - 1);
      bits[n] = (cubeedges[cells[n]][0] > 0) ? 1 : 0;
    }
  }
#endif
  bbuf->put_bits((dy - 1) * (dx - 1), bits);
  delete[] bits;
#ifdef ZP_CODEC
  bbuf->zp_encode();
#else
  bbuf->arith_encode();
#endif
  switch (out_flag) {
  case 0:
    writeBitBuffer(bbuf, out);
    break;
  case 1:
    writeBitBuffer(bbuf, stream);
    break;
  }
#ifdef _DEBUG
  printf("layer bitmap uses %d compressed bytes\n", bbuf->getNumBytes());
#endif
  delete bbuf;
}

// run-length encoding of the intersected cells
template <class T> void CompCon<T>::runLength(int count_c, u_char *cells) {
  if (count_c == 0) {
    return;
  }
  int nstrip = 0;
  u_char len = 0;
  vector<u_char> vec;
  u_char startx, starty; // starting cell of a strip
  bool stripOn = false;
  for (int j = 0; j < dy - 1; j++) {
    for (int i = 0; i < dx - 1; i++) {
      int n = i + j * (dx - 1);
      if (!stripOn && cubeedges[cells[n]][0] > 0) // a strip starts
      {
        startx = i;
        starty = j;
        nstrip++;
        len = 1;
        stripOn = true;
      } else if (stripOn && cubeedges[cells[n]][0] > 0 && len < 255) {
        // a strip continues
        len++;
        // here is an assumption that len <= 256, will fix it later
      } else if (stripOn && (cubeedges[cells[n]][0] == 0 || len == 255)) {
        // a strip ends
        vec.push_back(startx);
        vec.push_back(starty);
        vec.push_back(len);
        stripOn = false;
        len = 0;
      } else {
        // do nothing
      }
    }
  }
  if (stripOn && len > 0) {
    vec.push_back(startx);
    vec.push_back(starty);
    vec.push_back(len);
  }
  printf("num of strips = %d\n", nstrip);
  assert(vec.size() == nstrip * 3);
  for (int i = 0; i < nstrip * 3; i++) {
    switch (out_flag) {
    case 0:
      out->put(&vec[i], 1);
      break;
    case 1:
      stream->write8(vec[i]);
      break;
    }
  }
}

template <class T> void CompCon<T>::MarchingSlice(float isoval) {
  // byte array representing the configuration of cells
  u_char *cells = new u_char[(dx - 1) * (dy - 1)];
  memset(cells, 0, sizeof(u_char) * (dx - 1) * (dy - 1));
  int i, j, count_c = 0;
  T data[8];
  for (j = 0; j < dy - 1; j++) {
    for (i = 0; i < dx - 1; i++) {
      int n = i + j * (dx - 1);
      getCellValues(i, j, data);
      if ((float)data[0] < isoval) {
        cells[n] |= 0x01;
      }
      if ((float)data[1] < isoval) {
        cells[n] |= 0x02;
      }
      if ((float)data[2] < isoval) {
        cells[n] |= 0x04;
      }
      if ((float)data[3] < isoval) {
        cells[n] |= 0x08;
      }
      if ((float)data[4] < isoval) {
        cells[n] |= 0x10;
      }
      if ((float)data[5] < isoval) {
        cells[n] |= 0x20;
      }
      if ((float)data[6] < isoval) {
        cells[n] |= 0x40;
      }
      if ((float)data[7] < isoval) {
        cells[n] |= 0x80;
      }
      int code = cells[n];
      // increment the # of intersected cell
      if (cubeedges[code][0] > 0) {
        count_c++;
      }
      // mark the end vertices of intersected edges
      for (int e = 0; e < cubeedges[code][0]; e++) {
        markEdge(i, j, cubeedges[code][1 + e]);
      }
    }
  }
  if (wflag) {
    encodeSlice(count_c, cells, isoval);
  } else {
    scon->addSlice(*first);
    Layer *lay = new Layer(dx, dy, cells);
    scon->addLayer(*lay);
    delete lay;
  }
  delete[] cells;
}

template <class T> void CompCon<T>::writeHeader(float isoval) {
  switch (out_flag) {
  case 0:
    out->put(orig, 3);
    out->put(span, 3);
    out->put(&dx, 1);
    out->put(&dy, 1);
    out->put(&nslices, 1);
    out->put(&isoval, 1);
    break;
  case 1:
    stream->write(orig, sizeof(float) * 3);
    stream->write(span, sizeof(float) * 3);
    stream->write32(dx);
    stream->write32(dy);
    stream->write32(nslices);
    stream->write(&isoval, sizeof(float));
    break;
  default:
    break;
  };
}

template <class T> void CompCon<T>::marchingCubes(float isoval) {
  if (wflag) {
    // write out the global mesh information
    writeHeader(isoval);
  }
  for (int i = 0; i < nslices - 1; i++) {
    MarchingSlice(isoval);
    getNewSlice();
  }
  if (wflag) {
    switch (out_flag) {
    case 0:
      next->writeOut(out, isoval);
      break;
    case 1:
      next->writeOut(stream, isoval);
      break;
    }
  } else {
    scon->addSlice(*next);
  }
}

template <class T> void CompCon<T>::markEdge(int i, int j, int e) {
  int ix, iy;
  ix = i + edgeinfo[e].di;
  iy = j + edgeinfo[e].dj;
  if (edgeinfo[e].dk == 0) {
    first->setBit(ix, iy, edgedir[edgeinfo[e].dir]);
    switch (edgeinfo[e].dir) {
    case 0:
      first->setBit(ix + 1, iy, NEGX);
      break;
    case 1:
      first->setBit(ix, iy + 1, NEGY);
      break;
    case 2:
      next->setBit(ix, iy, NEGZ);
    }
  } else if (edgeinfo[e].dk == 1) {
    next->setBit(ix, iy, edgedir[edgeinfo[e].dir]);
    switch (edgeinfo[e].dir) {
    case 0:
      next->setBit(ix + 1, iy, NEGX);
      break;
    case 1:
      next->setBit(ix, iy + 1, NEGY);
      break;
    }
  }
}

#endif
