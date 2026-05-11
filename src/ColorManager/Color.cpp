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
#include <ColorManager/Color.h>

using namespace ColorManager;

Color::Color(float r, float g, float b, float a) : Tuple(r, g, b, a) {}

Color::Color() : Tuple() {}

Color::Color(float *array) { set(array); }

Color::~Color() {}

Color::Color(const Color &copy) : Tuple(copy) {}

Color &Color::operator=(const Color &copy) {
  if (this != &copy) {
    set(copy);
  }
  return *this;
}

Color &Color::set(float r, float g, float b, float a) {
  Tuple::set(r, g, b, a);
  return *this;
}

Color &Color::set(float *array) {
  Tuple::set(array);
  return *this;
}

Color &Color::set(const Color &copy) {
  Tuple::set(copy);
  return *this;
}

int Color::getMaxRGBIndex(const Color color) {
  if (color[0] > color[1] && color[0] > color[2]) {
    return 0;
  }
  if (color[1] > color[0] && color[1] > color[2]) {
    return 1;
  }
  return 2;
}
