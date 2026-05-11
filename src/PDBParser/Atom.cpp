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
#include <PDBParser/Atom.h>

using PDBParser::Atom;

Atom::Atom() : DPG::Ball(0, 0, 0, PDBParser::elementTable[-1].radius) {
  m_Position[0] = m_Position[1] = m_Position[2] = 0.0;
  m_ResidueNum = 0;
  m_Chain = ' ';
  m_ElementIndex = -1;
  m_IsBoundaryAtom = false;
  m_MutatedAtomInformation = 0;
  m_UserDefinedColor = 0;
  m_interfacepropensity = 0;
  m_donorOrAcceptor = false;
  m_internalHbond = 0;
  // ADD to copy constructer also !
}

Atom::Atom(const Atom &a)
    : DPG::Ball(a.m_Position[0], a.m_Position[1], a.m_Position[2],
                PDBParser::elementTable[a.m_ElementIndex].radius) {
  m_Position[0] = a.m_Position[0];
  m_Position[1] = a.m_Position[1];
  m_Position[2] = a.m_Position[2];
  m_ResidueNum = a.m_ResidueNum;
  m_Chain = a.m_Chain;
  m_ElementIndex = a.m_ElementIndex;
  m_IsBoundaryAtom = a.m_IsBoundaryAtom;
  if (a.m_MutatedAtomInformation) {
    m_MutatedAtomInformation =
        new MutatedAtomInformation(*(a.m_MutatedAtomInformation));
  } else {
    m_MutatedAtomInformation = 0;
  }
  m_UserDefinedColor = a.m_UserDefinedColor;
  m_interfacepropensity = a.m_interfacepropensity;
  m_donorOrAcceptor = a.m_donorOrAcceptor;
  m_internalHbond = a.m_internalHbond;
  m_parent = a.m_parent;
}

Atom &Atom::operator=(const Atom &a) {
  if (this != &a) {
    m_Position[0] = a.m_Position[0];
    m_Position[1] = a.m_Position[1];
    m_Position[2] = a.m_Position[2];
    m_ResidueNum = a.m_ResidueNum;
    m_Chain = a.m_Chain;
    m_ElementIndex = a.m_ElementIndex;
    m_IsBoundaryAtom = a.m_IsBoundaryAtom;
    if (a.m_MutatedAtomInformation) {
      m_MutatedAtomInformation =
          new MutatedAtomInformation(*(a.m_MutatedAtomInformation));
    } else {
      m_MutatedAtomInformation = 0;
    }
    m_UserDefinedColor = a.m_UserDefinedColor;
    m_interfacepropensity = a.m_interfacepropensity;
    m_donorOrAcceptor = a.m_donorOrAcceptor;
    m_internalHbond = a.m_internalHbond;
    m_parent = a.m_parent;
  }
  return *this; // Return ref for multiple assignment
}

Atom::~Atom() {
  delete m_MutatedAtomInformation;
  m_MutatedAtomInformation = 0;
}

float Atom::getRadius() const {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].radius;
  }
  if (m_MutatedAtomInformation) {
    return m_MutatedAtomInformation->m_Radius;
  }

  return PDBParser::defaultInformation.radius;
}

float Atom::getRadius(PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) const {
  if (radiusType == PDBParser::GroupOfAtoms::VDW_RADIUS) {
    return getRadius();
  }
  if (radiusType == PDBParser::GroupOfAtoms::ROUX_RADIUS) {
    if (m_ElementIndex >= 0 &&
        m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
      return PDBParser::elementTable[m_ElementIndex].solventEnlargedRadius;
    }
  }
  return getRadius();
}

float *Atom::getPosition() { return m_Position; }

int Atom::getResidueNum() { return m_ResidueNum; }

int Atom::getAtomNum() { return m_AtomNum; }

const char *Atom::getName() const {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].atomName;
  }
  if (m_MutatedAtomInformation) {
    return m_MutatedAtomInformation->m_Name;
  }
  return PDBParser::defaultInformation.atomName;
}

unsigned char Atom::getResidueId() const {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].residueIndex;
  }
  if (m_MutatedAtomInformation) {
    return m_MutatedAtomInformation->getResidueId();
  }
  return PDBParser::defaultInformation.residueIndex;
}

const char *Atom::getResidueName() const {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].residueName;
  }
  if (m_MutatedAtomInformation) {
    int resId = m_MutatedAtomInformation->getResidueId();
    if (resId >= 0) {
      return residueTable[resId].residueName;
    }
  }
  return PDBParser::defaultInformation.residueName;
}

char Atom::getChain() { return m_Chain; }

int Atom::getElementIndex() { return m_ElementIndex; }

float Atom::getCharge() {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].charge;
  }
  if (m_MutatedAtomInformation) {
    return m_MutatedAtomInformation->m_Charge;
  }
  return PDBParser::defaultInformation.charge;
}

float Atom::getPerAtomHydrophobicity() {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].hydrophobicity;
  }
  return PDBParser::defaultInformation.hydrophobicity;
}

float Atom::getPerResidueHydrophobicity() {
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    return PDBParser::elementTable[m_ElementIndex].perResidueHydrophobicity;
  }
  return PDBParser::defaultInformation.perResidueHydrophobicity;
}

// WARNING: this will change the atom's status to mutated.
void Atom::setRadius(const float r) {
  if (!m_MutatedAtomInformation) {
    m_MutatedAtomInformation = new MutatedAtomInformation(
        getRadius(), getName(), getResidueId(), getCharge());
  }
  m_ElementIndex = -1; // no longer related to table of elements.
  m_MutatedAtomInformation->m_Radius = r;
  radius = r;
}

void Atom::setPosition(const float *position) {
  m_Position[0] = position[0];
  m_Position[1] = position[1];
  m_Position[2] = position[2];
  x = m_Position[0];
  y = m_Position[1];
  z = m_Position[2];
}

void Atom::setPosition(float a, float b, float c) {
  m_Position[0] = a;
  m_Position[1] = b;
  m_Position[2] = c;
  x = m_Position[0];
  y = m_Position[1];
  z = m_Position[2];
}

void Atom::setElementIndex(int elementIndex) {
  m_ElementIndex = elementIndex;
  radius = PDBParser::elementTable[m_ElementIndex].radius;
}

// WARNING: this will change the atom's status to mutated.
void Atom::setName(char *name) {
  if (!m_MutatedAtomInformation) {
    m_MutatedAtomInformation = new MutatedAtomInformation(
        getRadius(), getName(), getResidueId(), getCharge());
  }
  m_ElementIndex = -1; // no longer related to table of elements.
  radius = PDBParser::elementTable[m_ElementIndex].radius;
  strcpy(m_MutatedAtomInformation->m_Name, name);
}

void Atom::setResidueNum(int resNum) { m_ResidueNum = resNum; }

void Atom::setAtomNum(int atNum) { m_AtomNum = atNum; }

// WARNING: this will change the atom's status to mutated.
void Atom::setResidueId(unsigned char residueId) {
  if (!m_MutatedAtomInformation) {
    m_MutatedAtomInformation = new MutatedAtomInformation(
        getRadius(), getName(), getResidueId(), getCharge());
  }
  m_ElementIndex = -1; // no longer related to table of elements.
  radius = PDBParser::elementTable[m_ElementIndex].radius;
  m_MutatedAtomInformation->setResidueId(residueId);
}

void Atom::setChain(char chain) { m_Chain = chain; }

void Atom::setCharge(float charge) {
  if (!m_MutatedAtomInformation) {
    m_MutatedAtomInformation = new MutatedAtomInformation(
        getRadius(), getName(), getResidueId(), getCharge());
  }
  m_ElementIndex = -1; // no longer related to table of elements.
  radius = PDBParser::elementTable[m_ElementIndex].radius;
  m_MutatedAtomInformation->m_Charge = charge;
}

bool Atom::isBoundaryAtom() {
  return m_IsBoundaryAtom;
  const char *n = getName();
  if (strcmp(n, " C  ") == 0) {
    return true;
  }
  if (strcmp(n, " CA ") == 0) {
    return true;
  }
  if (strcmp(n, " N  ") == 0) {
    return true;
  }
  return false;
}

void Atom::setBoundaryFlag(const bool flag) { m_IsBoundaryAtom = flag; }

void Atom::getColor(float *red, float *green, float *blue, float *alpha) {
  if (m_UserDefinedColor) {
    unsigned char aR, aG, aB;
    float aAlpha;
    aAlpha = ((float)(m_UserDefinedColor & 0xff)) / 255.0;
    aR = (unsigned char)((m_UserDefinedColor & 0xff000000) >> 24);
    aG = (unsigned char)((m_UserDefinedColor & 0xff0000) >> 16);
    aB = (unsigned char)((m_UserDefinedColor & 0xff00) >> 8);
    *red = aR / 255.0;
    *green = aG / 255.0;
    *blue = aB / 255.0;
    *alpha = aAlpha;
  }
  // default color: black
  if (m_ElementIndex >= 0 && m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
    *red = PDBParser::elementTable[m_ElementIndex].red;
    *blue = PDBParser::elementTable[m_ElementIndex].blue;
    *green = PDBParser::elementTable[m_ElementIndex].green;
    *alpha = 1.0;
    return;
  }
  if (m_MutatedAtomInformation) {
    *alpha = 1.0;
    if (m_MutatedAtomInformation->getColor(red, green, blue)) {
      return;
    }
  }
  *red = PDBParser::defaultInformation.red;
  *blue = PDBParser::defaultInformation.blue;
  *green = PDBParser::defaultInformation.green;
  *alpha = 1.0;
}

void Atom::getMinMax(double *minX, double *minY, double *minZ, double *maxX,
                     double *maxY, double *maxZ) {
  if (!minX || !minY || !minZ || !maxX || !maxY || !maxZ) {
    return;
  }
  const float r = getRadius();
  float *pos = 0;
  pos = getPosition();
  if (!pos) {
    return;
  }
  *minX = pos[0] - r;
  *minY = pos[1] - r;
  *minZ = pos[2] - r;
  *maxX = pos[0] + r;
  *maxY = pos[1] + r;
  *maxZ = pos[2] + r;
}

void Atom::getMaximumDistanceFromPoint(double xOrigin, double yOrigin,
                                       double zOrigin, double *distance) {
  if (!distance) {
    return;
  }
  float *pos = 0;
  pos = getPosition();
  if (!pos) {
    return;
  }
  double myDist = sqrt((pos[0] - xOrigin) * (pos[0] - xOrigin) +
                       (pos[1] - yOrigin) * (pos[1] - yOrigin) +
                       (pos[2] - zOrigin) * (pos[2] - zOrigin));
  if (*distance < myDist) {
    *distance = myDist;
  }
}

/*  The function returns -1 if there is no intersection  */
/*  in the positive direction. Else, it returns the      */
/*  distance from that point to the ray origin.          */
float Atom::getDistanceFromRay(CCVOpenGLMath::Ray targetVector) {
  // get atom center //
  float *pos = getPosition();
  if (!pos) {
    return -1;
  }
  CCVOpenGLMath::Vector atomCenter(pos[0], pos[1], pos[2], 1);
  // return values from intersection code //
  CCVOpenGLMath::Vector intersectionPoint1;
  CCVOpenGLMath::Vector intersectionPoint2;
  float distance1;
  float distance2;
  // perform ray sphere intersection //
  if (targetVector.intersectSphere(atomCenter, getRadius(), &intersectionPoint1,
                                   &intersectionPoint2, &distance1,
                                   &distance2)) {
    if (distance1 > 0) // should this be >= ? SKVINAY
    {
      return distance1;
    }
    if (distance2 > 0) {
      return distance2;
    }
  }
  // no valid intersections found
  return -1;
}

PDBParser::GroupOfAtoms::TYPE_OF_RESIDUE Atom::getResidueType() const {
  if (getResidueId() == 0) {
    return PDBParser::GroupOfAtoms::UNKNOWN_RESIDUE;
  }
  if (getResidueId() >= PDBParser::RES_TYPE_FIRST_NUCLEOTIDE) {
    return PDBParser::GroupOfAtoms::RNA_RESIDUE;
  }
  return PDBParser::GroupOfAtoms::AMINO_ACID_RESIDUE;
}

double Atom::getFunction(PDBParser::GroupOfAtoms::FUNCTIONS function) {
  if (function >= PDBParser::GroupOfAtoms::NUMBER_OF_FUNCTIONS) {
    return 0;
  }
  switch (function) {
  case PDBParser::GroupOfAtoms::ELECTRON_DENSITY:
    return 1;
  case PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY:
    if (m_ElementIndex >= 0 &&
        m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
      return PDBParser::elementTable[m_ElementIndex].hydrophobicity;
    }
    return 0;
  case PDBParser::GroupOfAtoms::PER_RESIDUE_HYDROPHOBICITY:
    if (m_ElementIndex >= 0 &&
        m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
      int residueIndex = PDBParser::elementTable[m_ElementIndex].residueIndex;
      if (residueIndex < 0 || residueIndex >= TOTAL_RES_TYPE) {
        return 0;
      }
      return PDBParser::residueTable[residueIndex].hydrophobicity;
    }
  default:
    return 0;
  }
}

//    we assume that the transformation does not contain any scale or shear
bool Atom::transform(CCVOpenGLMath::Matrix transformationMatrix) {
  CCVOpenGLMath::Vector pos(m_Position[0], m_Position[1], m_Position[2], 1);
  CCVOpenGLMath::Vector newPos = transformationMatrix * pos;
  if (newPos[3]) {
    x = m_Position[0] = newPos[0] / newPos[3];
    y = m_Position[1] = newPos[1] / newPos[3];
    z = m_Position[2] = newPos[2] / newPos[3];
  } else {
    // what else to do with such points ?
    x = m_Position[0] = newPos[0];
    y = m_Position[1] = newPos[1];
    z = m_Position[2] = newPos[2];
  }
  return true;
}
