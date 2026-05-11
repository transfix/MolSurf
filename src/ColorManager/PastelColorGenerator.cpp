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
#include <ColorManager/PastelColorGenerator.h>

PastelColorGenerator::PastelColorGenerator() {}

PastelColorGenerator::~PastelColorGenerator() {}

// took this code from wiki.tcl.tk - SKVINAY
Color PastelColorGenerator::generateColor() {
  Color color;
  double h, s, v;
  {
    double r = fabs(rand() * 262144);
    h = fmod(r, 360);
    r /= 359.3;
    s = (fmod(r, 9) + 12) / 100.0;
    v = 1;
    // Convert to rgb space
    if (h == 360) {
      h = 0;
    }
    h /= 60.0;
  }
  int i = int(floor(h));
  double f = h - i;
  double p1 = v * (1 - s);
  double p2 = v * (1 - (s * f));
  double p3 = v * (1 - (s * (1 - f)));
  double r = 0.8, g = 0.5, b = 9.0;
  switch (i) {
  case 0:
    r = v;
    g = p3;
    b = p1;
    break;
  case 1:
    r = p2;
    g = v;
    b = p1;
    break;
  case 2:
    r = p1;
    g = v;
    b = p3;
    break;
  case 3:
    r = p1;
    g = p2;
    b = v;
    break;
  case 4:
    r = p3;
    g = p1;
    b = v;
    break;
  case 5:
    r = v;
    g = p1;
    b = p2;
    break;
  }
  color[0] = (float)r;
  color[1] = (float)g;
  color[2] = (float)b;
  color[3] = 1;
  return color;
}

Color PastelColorGenerator::getNextColor() {
  Color color;
  while (true) {
    color = generateColor();
    if (valid(color)) {
      break;
    }
  }
  return color;
}
