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
#ifndef CCV_MULTI_RESOLUTION_MOLECULAR_SURFACE
#define CCV_MULTI_RESOLUTION_MOLECULAR_SURFACE

#include <Utility/utility.h>
#include <MultiresolutionSurface/Extents.h>
#include <PDBParser/GroupOfAtoms.h>

class Geometry;

class MultiresolutionSurface {
public:
  MultiresolutionSurface(PDBParser::GroupOfAtoms *molecule,
                         const Extents *boundary, const Extents *subVolume,
                         double innerB, double outerB, double innerIsovalue,
                         double outerIsovalue);
  virtual ~MultiresolutionSurface();
  bool changeInnerB(double innerB);
  bool changeOuterB(double outerB);
  bool changeInnerIsovalue(double innerIsovalue);
  bool changeOuterIsovalue(double outerIsovalue);
  bool changeInnerExtent(const Extents *subVolume);
  Geometry *getMultiresSurface();

protected:
  bool updateAtomList();
  PDBParser::GroupOfAtoms *m_Molecule;
  const Extents *m_Boundary;
  const Extents *m_SubVolume;
  double m_InnerB;
  double m_OuterB;
  double m_InnerIsovalue;
  double m_OuterIsovalue;
  Geometry *m_Geometry;
  std::vector<double> m_InteriorAtoms, m_ExteriorAtoms;
};

#endif
