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
#include <MultiresolutionSurface/MultiresolutionSurface.h>

MultiresolutionSurface::MultiresolutionSurface(
    PDBParser::GroupOfAtoms *molecule, const Extents *boundary,
    const Extents *subVolume, double innerB, double outerB,
    double innerIsovalue, double outerIsovalue) {
  m_Molecule = molecule;
  m_Boundary = boundary;
  m_SubVolume = subVolume;
  m_InnerB = innerB;
  m_OuterB = outerB;
  m_InnerIsovalue = innerIsovalue;
  m_OuterIsovalue = outerIsovalue;
  updateAtomList();
}

MultiresolutionSurface::~MultiresolutionSurface() {}

bool MultiresolutionSurface::changeInnerB(double innerB) {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_InnerB = innerB;
  printf("New inner b = %lf\n", m_InnerB);
  return true;
}

bool MultiresolutionSurface::changeOuterB(double outerB) {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_OuterB = outerB;
  printf("New outer b = %lf\n", m_OuterB);
  return true;
}

bool MultiresolutionSurface::changeInnerIsovalue(double innerIsovalue) {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_InnerIsovalue = innerIsovalue;
  printf("New innerIsovalue = %lf\n", m_InnerIsovalue);
  return true;
}

bool MultiresolutionSurface::changeOuterIsovalue(double outerIsovalue) {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_OuterIsovalue = outerIsovalue;
  printf("New outerIsovalue = %lf\n", m_OuterIsovalue);
  return true;
}

bool MultiresolutionSurface::updateAtomList() {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_InteriorAtoms.clear(), m_ExteriorAtoms.clear();
  double min[3], max[3];
  min[0] = m_SubVolume->getXMin();
  min[1] = m_SubVolume->getYMin();
  min[2] = m_SubVolume->getZMin();
  max[0] = m_SubVolume->getXMax();
  max[1] = m_SubVolume->getYMax();
  max[2] = m_SubVolume->getZMax();
  m_Molecule->getAtoms(min, max, &m_InteriorAtoms, &m_ExteriorAtoms,
                       PDBParser::GroupOfAtoms::BOTH);
  printf("Number of atoms in =%lu, out = %lu\n", m_InteriorAtoms.size(),
         m_ExteriorAtoms.size() / 4);
  return true;
}

bool MultiresolutionSurface::changeInnerExtent(const Extents *subVolume) {
  if (!m_Molecule || !m_Boundary || !m_SubVolume) {
    return false;
  }
  m_SubVolume = subVolume;
  printf("New extents");
  if (!updateAtomList()) {
    return false;
  }
  return true;
}

Geometry *MultiresolutionSurface::getMultiresSurface() { return m_Geometry; }
