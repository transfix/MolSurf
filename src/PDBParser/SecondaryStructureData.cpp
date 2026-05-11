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
#include <PDBParser/SecondaryStructureData.h>

using PDBParser::SecondaryStructureData;

SecondaryStructureData::SecondaryStructureData() {}

SecondaryStructureData::SecondaryStructureData(
    const SecondaryStructureData &secondaryStructureData) {
  m_ResSequenceNumberPairList =
      secondaryStructureData.m_ResSequenceNumberPairList;
  // The geometry is defined as a 7 tuple, containing
  // x1, y1, z1, x2, y2, z2, radius
  m_StrandGeometry = secondaryStructureData.m_StrandGeometry;
  // end for sheets
  // Used for Helices
  m_HelixGeometry = secondaryStructureData.m_HelixGeometry;
}

SecondaryStructureData::~SecondaryStructureData() {}
