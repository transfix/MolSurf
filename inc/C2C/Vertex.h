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
#ifndef VERTEX_H
#define VERTEX_H

#include <Utility/utility.h>
#include <C2C/arithlib/unitypes.h>

const u_char POSX = 0x01;
const u_char POSY = 0x02;
const u_char POSZ = 0x04;
const u_char NEGX = 0x08;
const u_char NEGY = 0x10;
const u_char NEGZ = 0x20;

static u_char edgedir[6] = {POSX, POSY, POSZ, NEGX, NEGY, NEGZ};

template <class T> class Slice;

template <class T> class Vertex {
private:
  T val;
  u_char code;

public:
  Vertex() {
    val = 0;
    code = 0;
  }
  Vertex(const T &x) : val(x) { code = 0; }
  Vertex(const Vertex<T> &vtx) : val(vtx.val) { code = vtx.code; }
  ~Vertex() {}

  Vertex<T> &operator=(const Vertex<T> &vtx);
  void setValue(const T &x, u_char c = 0);
  void setBit(u_char);
  bool isUsed() const;
  const T &getValue() { return val; }

  // reset the code bits
  void reset() { code = 0; }

  friend class Slice<T>;
};

template <class T> inline void Vertex<T>::setBit(u_char direct) {
  code |= direct;
}

template <class T> inline bool Vertex<T>::isUsed() const {
  return (code != 0x0);
}

template <class T> Vertex<T> &Vertex<T>::operator=(const Vertex<T> &vx) {
  if (this != &vx) {
    val = vx.val;
    code = vx.code;
  }
  return *this;
}

template <class T> inline void Vertex<T>::setValue(const T &x, u_char c) {
  val = x;
  code = c;
}

#endif
