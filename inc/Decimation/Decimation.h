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
#ifndef DECIMATE_H
#define DECIMATE_H

#include <LBIE_lib/Geoframe.h>

class Decimation {
private:
  double minx, miny, minz, maxx, maxy, maxz;

public:
  Decimation() {}
  void loadGeoframe(Geoframe &geo);
  void saveGeoframe(Geoframe &geo);
  void Decimate(const char *in_file, const char *out_file, int niter,
                float controllor, float variation, float innerangle);
  void DecimatePercentage(const char *in_file, const char *out_file,
                          float percentage, float control, float vari,
                          float inangle);
};

#endif
