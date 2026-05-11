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
#ifndef __LIBLBIE_H__
#define __LIBLBIE_H__

#include <iostream>
#include <math.h>
#include <LBIE_lib/octree.h>
#include <LBIE_lib/Geoframe.h>

#define DEFAULT_ERR 1.2501f
#define DEFAULT_ERR_IN 0.0001f
#define DEFAULT_IVAL -0.5001f
#define DEFAULT_IVAL_IN -9.5001f

#define NO_IMPROVE 0
#define GEO_FLOW 1
#define EDGE_CONTRACT 2
#define JOE_LIU 3
#define MINIMAL_VOL 4
#define OPTIMIZATION 5

namespace LBIE {

class LBIE_Mesher {

public:
  LBIE_Mesher(float in_iso_outer = DEFAULT_IVAL,
              float in_iso_inner = DEFAULT_IVAL_IN,
              float in_error_outer = DEFAULT_ERR,
              float in_error_inner = DEFAULT_ERR_IN,
              MESHTYPE in_meshtype = /*SINGLE*/ TRI,
              int in_improvemethod = NO_IMPROVE); // Meshing constructor

  void run(const char *input, const char *output, float isovalue);
  virtual ~LBIE_Mesher();
  void fileOpen(const char *);
  void setMeshType(MESHTYPE);
  void setMesh();
  void errorChange(float);
  void errorChange_in(float);
  void isovalueChange(float);
  void isovalueChange_in(float);
  MESHTYPE loadMesh(const char *input);
  int saveMesh(const char *output);
  void geoImprove();
  void geoImprove(float delta_t);

protected:
  float l_err;
  float l_err_in;
  float dist_level;
  int fopen_flag;
  char fname_buf[100];
  float isovalue;
  float isovalue_in;
  int dual_flag;
  double biggestDim;
  double centerx, centery, centerz;
  int improve_method;
  MESHTYPE meshtype;
  float iso_inner, iso_outer, inner_err_tol, outer_err_tol;
  Geoframe geoframe;
  Octree *oc;
};
} // namespace LBIE
#endif
