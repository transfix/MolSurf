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
#ifndef CVC_TEXMOL_POCKET_POCKET_H
#define CVC_TEXMOL_POCKET_POCKET_H

#include <Utility/utility.h>
#include <Blurmaps/MolecularSurface.h>
#include <GeometryFileTypes/GeometryLoader.h>
#include <Geometry/Geometry.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/parserPDBtoGOA.h>
#include <SignDistanceFunction_v2/DistanceTransform.h>
#include <SignDistanceFunction_v2/Geom3DParser.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>
#include <VolumeFileTypes/VolumeLoader.h>

class Geometry;
namespace PDBParser {
class GroupOfAtoms;
}

class Pocket {
public:
  Pocket();
  virtual ~Pocket();
  void constructPockets(PDBParser::GroupOfAtoms *pAtoms, string fileName,
                        float distance, Geometry **pocketGeometry, int size);
  void constructPockets(Geometry *geometry, string fileName, float distance,
                        Geometry **pocketGeometry, int size);
};

#endif
