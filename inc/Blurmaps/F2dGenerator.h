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
#ifndef CCV_F2D_GENERATOR_H
#define CCV_F2D_GENERATOR_H

#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/GOABlur.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <UsefulMath/TrilinearGrid.h>

using CCVOpenGLMath::TrilinearGrid;

namespace PDBParser {
class GroupOfAtoms;
class Atom;
}; // namespace PDBParser

class F2DGenerator {
public:
  F2DGenerator(PDBParser::GroupOfAtoms *molecule, string xyzFileName,
               string outputFileName, bool receptor,
               PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, double sar);

  virtual ~F2DGenerator();
  bool generate();

protected:
  vector<PDBParser::Atom *> m_AtomList;
  bool receptor;
  string xyzFileName;
  string outputFileName;

private:
  double diff(double a, double b);
  void printAtom(int index, bool internal);
  void printSkinAtom(double x, double y, double z, int i);

  double skinAtomRad;
  int numres;
  FILE *f2dFile;
  FILE *newpqrFile;
};

#endif
