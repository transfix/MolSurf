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
#include <Blurmaps/MolecularSurface.h>

MolecularSurface::MolecularSurface() {}

MolecularSurface::~MolecularSurface() {}

Geometry *
MolecularSurface::getMolecularSurface(PDBParser::GroupOfAtoms *molecule,
                                      int size) {
  return getMolecularSurface(molecule, size, 0);
}

Geometry *
MolecularSurface::getMolecularSurface(PDBParser::GroupOfAtoms *molecule,
                                      int size, double enlargedRadius) {
  // cout << "a" << flush;
  PDBParser::Atom *firstAtom;
  PDBParser::GroupOfAtoms *grp =
      molecule->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0]->m_SubGroups[0];
  firstAtom = grp->m_Atoms.at(0);
  bool useDoubleConvolution = false;
  if (useDoubleConvolution) {
    // algorithm
    //
    // Compute a volume with B = -0.5
    // Compute the volume
    // Compute the skin volume with B = -2.3
    // subtract 0.1* skin from volume
    // obtain isocontour with isovalaue = 1;
    double intVolumeBlobbiness = -0.5;
    SimpleVolumeData *interiorSimpleVolumeData = BlurMapsDataManager::getVolume(
        molecule, 0, size, size, size,
        PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, intVolumeBlobbiness,
        PDBParser::GroupOfAtoms::ATOM, 0, 0,
        PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::ATOM_TYPE);
    if (!interiorSimpleVolumeData) {
      return 0;
    }
    double probeRadius = 1.4;
    AtomPopulator *atomPopulator =
        new AtomPopulator(molecule, size, size, size, probeRadius,
                          PDBParser::GroupOfAtoms::VDW_RADIUS);
    PDBParser::GroupOfAtoms *skin =
        PDBParser::GroupOfAtoms::getEmptyHierarchy();
    if (!atomPopulator->populate(skin)) {
      delete skin;
      delete atomPopulator;
      return 0;
    }
    delete atomPopulator;
    double skinVolumeBlobbiness = -5.0;
    SimpleVolumeData *skinSimpleVolumeData = BlurMapsDataManager::getVolume(
        skin, 0, size, size, size, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
        false, skinVolumeBlobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0,
        PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::ATOM_TYPE);
    delete skin;
    if (!skinSimpleVolumeData) {
      delete interiorSimpleVolumeData;
      return 0;
    }
    interiorSimpleVolumeData->addVolume(skinSimpleVolumeData, -0.05, 0.0);
    double isovalue = 3.0;
    return SimpleVolumeDataIsocontourer::getIsocontour(interiorSimpleVolumeData,
                                                       isovalue);
  } else {
    if (!molecule) {
      return NULL;
    }
    SimpleVolumeData *vol = new SimpleVolumeData(size, size, size);
    // unused below? arand commented 5-11-2011
    // SimpleVolumeData* volcolor = new SimpleVolumeData(size, size, size);
    double probeRadius = 1.4;
    SkinRegion2 *skin = new SkinRegion2();
    int depth = 3;
    bool usecolor = false;

    skin->setEnlargedRadius(enlargedRadius);

    if (!skin->getSkinRegion(molecule, size, size, size, probeRadius,
                             PDBParser::GroupOfAtoms::VDW_RADIUS, vol, depth, 0,
                             usecolor, PDBParser::GroupOfAtoms::RESIDUE)) {
      delete vol;
      vol = 0;
      delete skin;
      skin = 0;
      return NULL;
    }
    delete skin;
    Geometry *geometry =
        SimpleVolumeDataIsocontourer::getIsocontour(vol, probeRadius);

    delete vol;
    return geometry;
  }
  // cout << "b" << flush;
}
