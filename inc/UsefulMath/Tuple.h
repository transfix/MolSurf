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
#if !defined(AFX_TUPLE_H__5AD4C604_B71A_4924_941A_15A0955C4E4E__INCLUDED_)
#define AFX_TUPLE_H__5AD4C604_B71A_4924_941A_15A0955C4E4E__INCLUDED_

#include <Utility/utility.h>

namespace CCVOpenGLMath {
class Tuple {
public:
  Tuple(float x, float y, float z, float w);
  Tuple();
  virtual ~Tuple();
  Tuple(const Tuple &copy);
  Tuple &operator=(const Tuple &copy);
  void print() const;
  Tuple &set(float x, float y, float z, float w);
  Tuple &set(float *array);
  Tuple &set(const Tuple &copy);
  float &operator[](unsigned int i);
  const float &operator[](unsigned int i) const;

protected:
  float p[4];
};
}; // namespace CCVOpenGLMath

#endif
