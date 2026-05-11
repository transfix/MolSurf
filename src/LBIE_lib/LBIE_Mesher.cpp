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
#include <Utility/utility.h>
#include <LBIE_lib/LBIE_Mesher.h>
#include <LBIE_lib/Geoframe.h>

using namespace CVCUtility;
using namespace LBIE;

// moved to utility...
/*
bool endsWith(string str, string substr)
{
        if(str.length() < substr.length())
        {
                return false;
        }
        if(str.substr(str.length() - substr.length(), substr.length()) ==
substr)
        {
                return true;
        }
        return false;
}
*/

LBIE_Mesher::LBIE_Mesher(float in_iso_outer, float in_iso_inner,
                         float in_error_outer, float in_error_inner,
                         MESHTYPE in_meshtype, int in_improvemethod) {
  fopen_flag = 0;
  dual_flag = 1;
  // initializing varables.
  iso_outer = in_iso_outer;
  iso_inner = in_iso_inner;
  outer_err_tol = in_error_outer;
  inner_err_tol = in_error_inner;
  improve_method = in_improvemethod;
  setMeshType(in_meshtype);
  oc = new Octree();
}

// Craig: Can we get this to return a Geometry* ?
// Molecule* LBIE_Mesher::run(PDBParser::GroupOfAtoms* molecule, float isovalue)
void LBIE_Mesher::run(const char *input, const char *output, float isovalue) {

  if (endsWith(string(input), "rawiv")) {
    fileOpen(input);
    iso_outer = isovalue; // setting the mesh type

    setMesh();
    // set the appropriate isovalues and error values
    if (iso_outer != DEFAULT_IVAL) {
      isovalueChange(iso_outer);
    }

    if (outer_err_tol != DEFAULT_ERR) {
      errorChange(outer_err_tol);
    }
    if (meshtype == TRI2 || meshtype == TET2) {
      if (iso_inner != DEFAULT_IVAL_IN) {
        isovalueChange_in(iso_inner);
      }
      if (inner_err_tol != DEFAULT_ERR_IN) {
        errorChange_in(inner_err_tol);
      }
    }
  } else if (endsWith(string(input), "raw") ||
             endsWith(string(input), "rawc") ||
             endsWith(string(input), "rawn")) {
    improve_method = GEO_FLOW;
    setMeshType(loadMesh(input));
    geoImprove(isovalue); // this is actually the delta_t parameter
                          //  for quality improvement
  }
  saveMesh(output);
}

LBIE_Mesher::~LBIE_Mesher() {}

void LBIE_Mesher::setMeshType(MESHTYPE mymeshtype) {

  // cout << "SetMeshType " << mymeshtype << endl;

  meshtype = mymeshtype;
  if (mymeshtype == TRI2 || mymeshtype == TET2) {
    dual_flag = 2;
  }
}

void LBIE_Mesher::setMesh() {
  oc->setMeshType(meshtype);
  if (dual_flag == 1) // single isovalue
  {
    oc->collapse();
    oc->compute_qef();
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
  } else if (dual_flag == 2) // interval volume
  {
    oc->collapse_interval();
    oc->compute_qef_interval();
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
  }
}

void LBIE_Mesher::errorChange(float err) {
  if (err <= 0) {
    err = 0.0000001f;
  }
  l_err = err;
  std::cout << "Inside errorChange. outer_err_tol is " << err << "\n";
  geoframe = Geoframe();
  geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
  geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
  if (dual_flag == 1) {
    l_err = err;
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  } else if (dual_flag == 2) {
    l_err = err;
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  }
}

void LBIE_Mesher::errorChange_in(float err) {
  if (err <= 0) {
    err = 0.0000001f;
  }
  l_err = err;
  std::cout << "Inside errorChange_in. inner_err_tol is " << err << "\n";
  geoframe = Geoframe();
  geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
  geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
  if (dual_flag == 1) {
    l_err = err;
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  } else if (dual_flag == 2) {
    l_err_in = err;
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  }
}

void LBIE_Mesher::isovalueChange(float myisovalue) {
  std::cout << "Inside isovalueChange. iso_outer is " << myisovalue << "\n";
  // geoframe = Geoframe(); // arand: this causes crashes in refactored code
  geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
  geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
  oc->set_isovalue(myisovalue);
  if (dual_flag == 1) // single isovalue
  {
    oc->collapse();
    oc->compute_qef();
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  } else if (dual_flag == 2) // interval volume
  {
    oc->collapse_interval();
    oc->compute_qef_interval();
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  }
}

void LBIE_Mesher::isovalueChange_in(float myisovalue_in) {
  std::cout << "Inside isovalueChange_in. iso_inner is " << myisovalue_in
            << "\n";
  // geoframe = Geoframe(); // arand: this causes crashes in refactored code
  geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
  geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
  oc->set_isovalue_in(myisovalue_in);
  dual_flag = 2;
  if (dual_flag == 1) {
    oc->collapse();
    oc->compute_qef();
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  } else if (dual_flag == 2) {
    oc->collapse_interval();
    oc->compute_qef_interval();
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  }
}

void LBIE_Mesher::fileOpen(const char *filename) {
  fopen_flag = 1;
  // geoframe = Geoframe();
  cout << "11 - bound_sign: " << geoframe.bound_sign
       << endl; // allocated somewhere
  cout << "11 - bound_sign[0]: " << geoframe.bound_sign[0]
       << endl; // can't dereference

  strcpy(fname_buf, filename);
  dist_level = 20;
  l_err = DEFAULT_ERR;
  l_err_in = DEFAULT_ERR_IN;
  oc->set_isovalue(DEFAULT_IVAL);
  oc->set_isovalue_in(DEFAULT_IVAL_IN);

  if (dual_flag == 1) // single isovalue
  {

    oc->Octree_loadVolume(filename);
    geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
    geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
    oc->collapse();
    oc->compute_qef();
    oc->traverse_qef(l_err);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  } else if (dual_flag == 2) // interval volume
  {
    oc->Octree_loadVolume(filename);
    geoframe.setSpan(oc->spans[0], oc->spans[1], oc->spans[2]);
    geoframe.setMin(oc->minext[0], oc->minext[1], oc->minext[2]);
    oc->collapse_interval();
    oc->compute_qef_interval();
    oc->traverse_qef_interval(l_err, l_err_in);
    oc->mesh_extract(geoframe, l_err);
    oc->quality_improve(geoframe, improve_method);
  }
  geoframe.calculateExtents();
  biggestDim = geoframe.biggestDim;
  centerx = geoframe.centerx;
  centery = geoframe.centery;
  centerz = geoframe.centerz;
}

int LBIE_Mesher::saveMesh(const char *output) {
  geoframe.write_raw(output, meshtype);
}

MESHTYPE LBIE_Mesher::loadMesh(const char *input) {
  // geoframe = Geoframe(input);
  geoframe.LoadFile(input);
  geoframe.setSpan(1.0, 1.0, 1.0);
  geoframe.setMin(0.0, 0.0, 0.0);
  return (geoframe.getMeshType());
}

void LBIE_Mesher::geoImprove() { geoImprove(0.01f); }

void LBIE_Mesher::geoImprove(float delta_t) {
  oc->setMeshType(meshtype);
  oc->quality_improve(geoframe, improve_method, delta_t);
}
