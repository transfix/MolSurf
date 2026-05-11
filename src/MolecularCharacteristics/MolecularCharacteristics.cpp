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
#include <MolecularCharacteristics/MolecularCharacteristics.h>

MolecularCharacteristics::MolecularCharacteristics(
    PDBParser::GroupOfAtoms *molecule, SimpleVolumeData *electrostatics,
    int numberOfSpheres, float *activeSite, double *radii,
    char *inputGOAFileName, char *inputElecVolumeFileName) {
  m_Molecule = molecule;
  m_Electrostatics = electrostatics;
  m_NumberOfSpheres = numberOfSpheres;
  m_ActiveSite = activeSite;
  m_Radii = radii;
  strcpy(m_GOAFileName, inputGOAFileName);
  strcpy(m_ElecFileName, inputElecVolumeFileName);
}

MolecularCharacteristics::~MolecularCharacteristics() {}

bool MolecularCharacteristics::getNumberOfAtoms(int *numberOfAtoms) {
  if (!numberOfAtoms) {
    return false;
  }
  return m_Molecule->getNumberOfAtomsRecursive(numberOfAtoms);
}

bool MolecularCharacteristics::getNumberOfChargedAtoms(
    int *numberOfHydrophobicAtoms, int *numberOfHydrophilicAtoms,
    PDBParser::GroupOfAtoms *molecule) {
  if (!numberOfHydrophobicAtoms || !numberOfHydrophilicAtoms || !molecule) {
    return false;
  }
  {
    int i;
    for (i = 0; i < molecule->m_Atoms.size(); i++) {
      PDBParser::Atom *atom = molecule->m_Atoms[i];
      if (atom) {
        double h = 0;
        if (atom->m_ElementIndex >= 0 &&
            atom->m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
          h = PDBParser::elementTable[atom->m_ElementIndex].hydrophobicity;
        }
        if (h > 0) {
          (*numberOfHydrophobicAtoms)++;
        }
        if (h < 0) {
          (*numberOfHydrophilicAtoms)++;
        }
      }
    }
  }
  int i, m;
  m = (int)molecule->m_SubGroups.size();
  for (i = 0; i < m; i++) {
    PDBParser::GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (currSubGroup) {
      if (!getNumberOfChargedAtoms(numberOfHydrophobicAtoms,
                                   numberOfHydrophilicAtoms, currSubGroup)) {
        return false;
      }
    }
  }
  return true;
}

bool MolecularCharacteristics::getNumberOfChargedAtoms(
    int *numberOfHydrophobicAtoms, int *numberOfHydrophilicAtoms, float *point,
    double radius, PDBParser::GroupOfAtoms *molecule) {
  if (!numberOfHydrophobicAtoms || !numberOfHydrophilicAtoms || !molecule) {
    return false;
  }
  {
    int i;
    for (i = 0; i < molecule->m_Atoms.size(); i++) {
      PDBParser::Atom *atom = molecule->m_Atoms[i];
      if (atom) {
        if (dist2(atom->m_Position, point) <= radius) {
          double h = 0;
          if (atom->m_ElementIndex >= 0 &&
              atom->m_ElementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS) {
            h = PDBParser::elementTable[atom->m_ElementIndex].hydrophobicity;
          }
          if (h > 0) {
            (*numberOfHydrophobicAtoms)++;
          }
          if (h < 0) {
            (*numberOfHydrophilicAtoms)++;
          }
        }
      }
    }
  }
  int i, m;
  m = (int)molecule->m_SubGroups.size();
  for (i = 0; i < m; i++) {
    PDBParser::GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (currSubGroup) {
      if (!getNumberOfChargedAtoms(numberOfHydrophobicAtoms,
                                   numberOfHydrophilicAtoms, point, radius,
                                   currSubGroup)) {
        return false;
      }
    }
  }
  return true;
}

double MolecularCharacteristics::dist2(float *p1, float *p2) {
  if (!p1 || !p2) {
    return 0;
  }
  return (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) +
         (p1[2] - p2[2]) * (p1[2] - p2[2]);
}

bool MolecularCharacteristics::getNumberOfAtoms(
    int *numberOfAtoms, float *point, double radius,
    PDBParser::GroupOfAtoms *molecule) {
  if (!numberOfAtoms || !point || (radius <= 0) || !molecule) {
    return false;
  }
  {
    int i;
    for (i = 0; i < molecule->m_Atoms.size(); i++) {
      PDBParser::Atom *atom = molecule->m_Atoms[i];
      if (atom) {
        if (dist2(atom->m_Position, point) <= radius) {
          (*numberOfAtoms)++;
        }
      }
    }
  }
  int i, m;
  m = (int)molecule->m_SubGroups.size();
  for (i = 0; i < m; i++) {
    PDBParser::GroupOfAtoms *currSubGroup = molecule->m_SubGroups.at(i);
    if (currSubGroup) {
      if (!getNumberOfAtoms(numberOfAtoms, point, radius, currSubGroup)) {
        return false;
      }
    }
  }
  return true;
}

double dist(float *p1, float *p2) {
  return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
              (p1[1] - p2[1]) * (p1[1] - p2[1]) +
              (p1[2] - p2[2]) * (p1[2] - p2[2]));
}

// this should be someplace else perhaps ?
bool MolecularCharacteristics::zeroOutVolumeShell(SimpleVolumeData *sData,
                                                  float *center, double radius,
                                                  char *filename) {
  if (!sData || !center || (radius < 0)) {
    return false;
  }
  unsigned int dim[3];
  float orig[3];
  float span[3];
  dim[0] = sData->getWidth();
  dim[1] = sData->getHeight();
  dim[2] = sData->getDepth();
  orig[0] = sData->getMinX();
  orig[1] = sData->getMinY();
  orig[2] = sData->getMinZ();
  span[0] = sData->getSpanX();
  span[1] = sData->getSpanY();
  span[2] = sData->getSpanZ();
  float *data = (float *)sData->getData(0);
  float point[3];
  int i, j, k, c = 0;
  for (i = 0; i < dim[0]; i++) {
    point[0] = orig[0] + i * span[0];
    for (j = 0; j < dim[1]; j++) {
      point[1] = orig[1] + j * span[1];
      for (k = 0; k < dim[2]; k++) {
        point[2] = orig[2] + k * span[2];
        if (dist(point, center) > radius) {
          data[c] = 0;
        }
        c++;
      }
    }
  }
  // create and save the volume
  {
    float max[3];
    max[0] = orig[0] + (dim[0] - 1) * span[0];
    max[1] = orig[1] + (dim[1] - 1) * span[1];
    max[2] = orig[2] + (dim[2] - 1) * span[2];
    float *d = new float[dim[0] * dim[1] * dim[2]];
    {
      int i;
      for (i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
        d[i] = data[i];
      }
    }
    SimpleVolumeData *sData = new SimpleVolumeData(dim);
    sData->setNumberOfVariables(1);
    sData->setData(0, d);
    sData->setType(0, SimpleVolumeData::FLOAT);
    sData->setName(0, "sphereRegion");
    sData->setMinExtent(orig);
    sData->setMaxExtent(max);
    VolumeLoader *vLoader = new VolumeLoader();
    vLoader->saveFile(filename, sData);
    delete vLoader;
    delete sData;
  }
  return true;
}

bool MolecularCharacteristics::getAreaVolume(double *area, double *volume,
                                             double *areaInSphere,
                                             double *volumeInSphere) {
  // blur
  // get area vol
  // in a loop, keep 0'ing out the outer shells,
  //    get area volume
  SimpleVolumeData *sData = 0;
  sData = BlurMapsDataManager::getVolume(
      m_Molecule, 0, 64, 64, 64, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
      false, -2.3, PDBParser::GroupOfAtoms::ATOM, 0, 0,
      PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::GroupOfAtoms::ATOM);
  if (!sData) {
    return false;
  }
  Geometry *geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, 1);
  if (!geometry) {
    delete sData;
    return false;
  }

  int dim[3];
  float orig[3];
  float span[3];
  dim[0] = sData->getWidth();
  dim[1] = sData->getHeight();
  dim[2] = sData->getDepth();
  orig[0] = sData->getMinX();
  orig[1] = sData->getMinY();
  orig[2] = sData->getMinZ();
  span[0] = sData->getSpanX();
  span[1] = sData->getSpanY();
  span[2] = sData->getSpanZ();
  {
    AreaVolume *areaVolumeModule = new AreaVolume();
    double t;
    areaVolumeModule->getVolume((float *)(sData->getData(0)), 1, dim, orig,
                                span, &t, volume);
    delete areaVolumeModule;
  }
  *area = geometry->getArea();
  {
    // get volumes progressively by 0ing out outer shells
    int i;
    for (i = m_NumberOfSpheres - 1; i >= 0; i--) {
      char shellDensFileName[1024];
      sprintf(shellDensFileName, "%s_%d.rawiv", m_GOAFileName, i);
      zeroOutVolumeShell(sData, m_ActiveSite, m_Radii[i], shellDensFileName);
      double t;
      AreaVolume *areaVolumeModule = new AreaVolume();
      areaVolumeModule->getVolume((float *)(sData->getData(0)), 1, dim, orig,
                                  span, &t, &(volumeInSphere[i]));
      delete areaVolumeModule;
    }
  }
  {
    int i;
    for (i = 0; i < m_NumberOfSpheres; i++) {
      areaInSphere[i] = geometry->getArea(m_ActiveSite, m_Radii[i]);
    }
  }
  delete sData;
  return true;
}

bool MolecularCharacteristics::getElecRatio(float *data, int *dim,
                                            double *ratio) {
  double total_neg = 0;
  double total_pos = 0;
  for (int c = 0; c < dim[0] * dim[1] * dim[2]; c++) {
    float v = data[c];
    if (v < 0) {
      total_neg += fabs(v);
    }
    if (v > 0) {
      total_pos += fabs(v);
    }
  }
  if (total_neg != 0) {
    (*ratio) = total_pos / total_neg;
  } else {
    printf("Warning, no neg potential at all\n");
    (*ratio) = 1e10;
  }
  return true;
}

bool MolecularCharacteristics::getElectrostaticCharacteristics(
    double *ratio, double *sphereRatios) {
  if (!m_Electrostatics || !ratio || !sphereRatios) {
    return false;
  }
  int dim[3];
  dim[0] = m_Electrostatics->getWidth();
  dim[1] = m_Electrostatics->getHeight();
  dim[2] = m_Electrostatics->getDepth();
  getElecRatio((float *)(m_Electrostatics->getData(0)), dim, ratio);
  {
    // get volumes progressively by 0ing out outer shells
    int i;
    for (i = m_NumberOfSpheres - 1; i >= 0; i--) {
      char shellElecFileName[1024];
      sprintf(shellElecFileName, "%s_%d.rawiv", m_ElecFileName, i);
      zeroOutVolumeShell(m_Electrostatics, m_ActiveSite, m_Radii[i],
                         shellElecFileName);
      getElecRatio((float *)(m_Electrostatics->getData(0)), dim,
                   &(sphereRatios[i]));
    }
  }
  return true;
}

bool MolecularCharacteristics::getMoleculeCharacteristics(
    int *numberOfAtoms, int *numberOfHydrophobicAtoms,
    int *numberOfHydrophilicAtoms, double *numberOfResidues,
    double *numberOfHydrophobicResidues, double *numberOfHydrophilicResidues,
    double *area, double *volume, int *numberOfAtomsInSphere,
    int *numberOfHydrophobicAtomsInSphere,
    int *numberOfHydrophilicAtomsInSphere, double *numberOfResiduesInSphere,
    double *numberOfHydrophobicResiduesInSphere,
    double *numberOfHydrophilicResiduesInSphere, double *areaInSphere,
    double *volumeInSphere) {
  if (!numberOfAtoms || !numberOfHydrophobicAtoms ||
      !numberOfHydrophilicAtoms || !numberOfResidues ||
      !numberOfHydrophobicResidues || !numberOfHydrophilicResidues || !area ||
      !volume) {
    return false;
  }
  if (m_NumberOfSpheres > 0) {
    if (!m_ActiveSite || !m_Radii || !numberOfAtomsInSphere ||
        !numberOfHydrophobicAtomsInSphere ||
        !numberOfHydrophilicAtomsInSphere || !numberOfResiduesInSphere ||
        !numberOfHydrophobicResiduesInSphere ||
        !numberOfHydrophilicResiduesInSphere || !areaInSphere ||
        !volumeInSphere) {
      return false;
    }
  }
  if (!getNumberOfAtoms(numberOfAtoms)) {
    return false;
  }
  if (!getNumberOfChargedAtoms(numberOfHydrophobicAtoms,
                               numberOfHydrophilicAtoms, m_Molecule)) {
    return false;
  }
  if (m_NumberOfSpheres < 1) {
    return true;
  }
  {
    int i;
    for (i = 0; i < m_NumberOfSpheres; i++) {
      if (!getNumberOfAtoms(&(numberOfAtomsInSphere[i]), m_ActiveSite,
                            m_Radii[i], m_Molecule)) {
        return false;
      }
      if (!getNumberOfChargedAtoms(&(numberOfHydrophobicAtomsInSphere[i]),
                                   &(numberOfHydrophilicAtomsInSphere[i]),
                                   m_ActiveSite, m_Radii[i], m_Molecule)) {
        return false;
      }
    }
  }
  if (!getAreaVolume(area, volume, areaInSphere, volumeInSphere)) {
    return false;
  }
  return true;
}

bool MolecularCharacteristics::appendCharacteristics(const char *filename) {
  if (!filename) {
    return false;
  }
  FILE *fp = fopen(filename, "a");
  if (!fp) {
    return false;
  }
  int numberOfAtoms = 0;
  int numberOfHydrophobicAtoms = 0;
  int numberOfHydrophilicAtoms = 0;
  double numberOfResidues = 0;
  double numberOfHydrophobicResidues = 0;
  double numberOfHydrophilicResidues = 0;
  double area = 0;
  double volume = 0;
  int *numberOfAtomsInSphere = 0;
  int *numberOfHydrophobicAtomsInSphere = 0;
  int *numberOfHydrophilicAtomsInSphere = 0;
  double *numberOfResiduesInSphere = 0;
  double *numberOfHydrophobicResiduesInSphere = 0;
  double *numberOfHydrophilicResiduesInSphere = 0;
  double *areaInSphere = 0;
  double *volumeInSphere = 0;
  if (m_NumberOfSpheres > 0)
    numberOfAtomsInSphere = new int[m_NumberOfSpheres];
  numberOfHydrophobicAtomsInSphere = new int[m_NumberOfSpheres];
  numberOfHydrophilicAtomsInSphere = new int[m_NumberOfSpheres];
  numberOfResiduesInSphere = new double[m_NumberOfSpheres];
  numberOfHydrophobicResiduesInSphere = new double[m_NumberOfSpheres];
  numberOfHydrophilicResiduesInSphere = new double[m_NumberOfSpheres];
  areaInSphere = new double[m_NumberOfSpheres];
  volumeInSphere = new double[m_NumberOfSpheres];
  for (int i = 0; i < m_NumberOfSpheres; i++) {
    numberOfAtomsInSphere[i] = 0;
    numberOfHydrophobicAtomsInSphere[i] = 0;
    numberOfHydrophilicAtomsInSphere[i] = 0;
    numberOfResiduesInSphere[i] = 0;
    numberOfHydrophobicResiduesInSphere[i] = 0;
    numberOfHydrophilicResiduesInSphere[i] = 0;
    areaInSphere[i] = 0;
    volumeInSphere[i] = 0;
  }
  if (!getMoleculeCharacteristics(
          &numberOfAtoms, &numberOfHydrophobicAtoms, &numberOfHydrophilicAtoms,
          &numberOfResidues, &numberOfHydrophobicResidues,
          &numberOfHydrophilicResidues, &area, &volume, numberOfAtomsInSphere,
          numberOfHydrophobicAtomsInSphere, numberOfHydrophilicAtomsInSphere,
          numberOfResiduesInSphere, numberOfHydrophobicResiduesInSphere,
          numberOfHydrophilicResiduesInSphere, areaInSphere, volumeInSphere)) {
    fclose(fp);
    return false;
  }
  double elecRatio;
  double *sphereElecRatios = new double[m_NumberOfSpheres];
  if (m_Electrostatics) {
    if (!getElectrostaticCharacteristics(&elecRatio, sphereElecRatios)) {
      fclose(fp);
      return false;
    }
  }
  fprintf(fp, "%6d %6d %6d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf",
          numberOfAtoms, numberOfHydrophobicAtoms, numberOfHydrophilicAtoms,
          numberOfResidues, numberOfHydrophobicResidues,
          numberOfHydrophilicResidues, area, volume, elecRatio);
  for (int i = 0; i < m_NumberOfSpheres; i++) {
    fprintf(fp, " %6d %6d %6d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf",
            numberOfAtomsInSphere[i], numberOfHydrophobicAtomsInSphere[i],
            numberOfHydrophilicAtomsInSphere[i], numberOfResiduesInSphere[i],
            numberOfHydrophobicResiduesInSphere[i],
            numberOfHydrophilicResiduesInSphere[i], areaInSphere[i],
            volumeInSphere[i], sphereElecRatios[i]);
  }
  fprintf(fp, "\n");
  fclose(fp);
  return true;
}

bool MolecularCharacteristics::appendHeader(const char *filename) {
  if (!filename) {
    return false;
  }
  FILE *fp = fopen(filename, "a");
  if (!fp) {
    return false;
  }
  {
    fprintf(fp, " Atoms  non-p  polar        boo        boo        boo       "
                "Area     Volume pos_by_neg");
    for (int i = 0; i < m_NumberOfSpheres; i++) {
      fprintf(fp, "  Atoms  non-p  polar        boo        boo        boo      "
                  " Area     Volume pos_by_neg");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return true;
}
