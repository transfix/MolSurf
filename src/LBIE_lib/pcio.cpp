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
#include <LBIE_lib/pcio.h>
#include <assert.h>

size_t getFloat(float *flts, size_t n, FILE *fp) {
  unsigned char *pb = new unsigned char[n * 4];
  unsigned char *pf = (unsigned char *)flts;
  size_t nbytes = fread(pb, 1, 4 * n, fp);
  // swap the byte order
  if (nbytes == n * 4) {
    for (size_t i = 0; i < n; i++) {
      pf[4 * i] = pb[4 * i + 3];
      pf[4 * i + 1] = pb[4 * i + 2];
      pf[4 * i + 2] = pb[4 * i + 1];
      pf[4 * i + 3] = pb[4 * i];
    }
  }
  delete pb;
  return nbytes;
}

size_t putFloat(float *flts, size_t n, FILE *fp) {
  unsigned char *pb = new unsigned char[n * 4];
  unsigned char *pf = (unsigned char *)flts;
  for (size_t i = 0; i < n; i++) {
    pb[4 * i] = pf[4 * i + 3];
    pb[4 * i + 1] = pf[4 * i + 2];
    pb[4 * i + 2] = pf[4 * i + 1];
    pb[4 * i + 3] = pf[4 * i];
  }
  size_t nbytes = fwrite(pb, 1, 4 * n, fp);
  delete pb;
  return nbytes;
}

size_t getInt(int *Ints, size_t n, FILE *fp) {
  unsigned char *pb = new unsigned char[4 * n];
  unsigned char *pf = (unsigned char *)Ints;
  int nbytes = fread(pb, 1, 4 * n, fp);
  for (size_t i = 0; i < n; i++) {
    pf[4 * i] = pb[4 * i + 3];
    pf[4 * i + 1] = pb[4 * i + 2];
    pf[4 * i + 2] = pb[4 * i + 1];
    pf[4 * i + 3] = pb[4 * i];
  }
  delete pb;
  return nbytes;
}

size_t getShort(short *shts, size_t n, FILE *fp) {
  unsigned char *pb = new unsigned char[n * 2];
  unsigned char *ps = (unsigned char *)shts;
  size_t nbytes = fread(pb, sizeof(unsigned char), n * 2, fp);
  // swap the byte order
  if (nbytes == n * 2) {
    for (size_t i = 0; i < n; i++) {
      ps[2 * i] = pb[2 * i + 1];
      ps[2 * i + 1] = pb[2 * i];
    }
  }
  delete pb;
  return nbytes;
}

size_t getUnChar(unsigned char *chars, size_t n, FILE *fp) {
  unsigned char *pb = new unsigned char[n];
  unsigned char *ps = (unsigned char *)chars;
  size_t nbytes = fread(pb, sizeof(unsigned char), n, fp);
  if (nbytes == n) {
    for (size_t i = 0; i < n; i++) {
      ps[i] = pb[i];
    }
  }
  delete pb;
  return nbytes;
}
