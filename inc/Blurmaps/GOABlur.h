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
#ifndef GOA_BLUR_H
#define GOA_BLUR_H

#include <Utility/utility.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/ElementInformation.h>
#include <PDBParser/GOAColor.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/ResidueInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace PDBParser {
class GroupOfAtoms;
class Atom;
}; // namespace PDBParser
class SimpleVolumeData;

/**
 * GOABlur -- A Class for creating blur maps from GroupOfAtoms objects
 */
class GOABlur {
public:
  ///
  GOABlur();
  virtual ~GOABlur();

  /**
   * Set the blobbyness of the blurring kernel
   */
  void setBlobbyness(double b) { mBlobbyness = b; }
  double getBlobbyness() const { return mBlobbyness; }
  /**
   * Compute the ElectronDensity
   * @return in SimpleVolumeData format
   */
  bool generateBlurmap(PDBParser::GroupOfAtoms *molecule,
                       SimpleVolumeData *volumeData,
                       PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                       const char *cmapFile, bool withColor,
                       PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                       PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                       unsigned int level);

  bool generateBlurmap(PDBParser::GroupOfAtoms *molecule,
                       SimpleVolumeData *volumeData, float bbmin[3],
                       float bbmax[3],
                       PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                       const char *cmapFile, bool withColor,
                       PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                       PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                       unsigned int level);

  bool generateBlurmap2(std::vector<PDBParser::Atom *> atomList,
                        SimpleVolumeData *vol,
                        PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                        const char *cmapFile, bool withColor,
                        PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                        PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                        unsigned int level);

protected:
  void blurAtoms(vector<PDBParser::Atom *> &atoms, SimpleVolumeData *vol,
                 PDBParser::GroupOfAtoms::FUNCTIONS densityType, bool withColor,
                 PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);

  double lookUpHydrophobicity(const PDBParser::Atom &atom);
  double lookUpCharge(const PDBParser::Atom &atom);
  float evalDensity(/*const PDBParser::Atom &atom,*/ double ctr[3],
                    float radius, float pnt[3],
                    PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                    double hydrophobicity, double charge);
  double
  getMaxRadiusOfAtomPotential(PDBParser::Atom *atom, double blobbiness,
                              PDBParser::GroupOfAtoms::FUNCTIONS potentialType,
                              double eps);

  std::vector<int>
      m_ColorRanges; // each index is an atom index. If an atom lies
  // between a range, then assign it that ranges color as specified by the next
  // vector
  std::vector<float> m_UserDefinedColors;

private:
  double mBlobbyness;
};

#endif
