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
#include <Blurmaps/BlurMapsDataManager.h>

#include <PDBParser/FlattenGOA.h>

BlurMapsDataManager::BlurMapsDataManager() {}

BlurMapsDataManager::~BlurMapsDataManager() {}
/*
bool BlurMapsDataManager::flattenGOA(PDBParser::GroupOfAtoms* grp,
vector<PDBParser::Atom*> & flatList, PDBParser::CollectionData* collectionData,
        int numberOfTransformations, double* rotations, double* translations,
PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level, bool
colorBySubunits)
{
        if(grp->type == PDBParser::CHAIN_TYPE)
        {
                char cID = (char)grp->resID; // resID is overloaded in GOA to be
chainID for chains if(collectionData)
                {
                        if(collectionData->isTransformationPresent(cID))
                        {
                                int nTrans =
collectionData->getNumberOfTransformations(cID); if(nTrans > 0)
                                {
                                        numberOfTransformations = nTrans;
                                        rotations    = new
double[numberOfTransformations*9]; translations = new
double[numberOfTransformations*3]; if(!collectionData->getRotations(cID,
rotations) || !collectionData->getTranslations(cID, translations))
                                        {
                                                numberOfTransformations = 0;
                                                delete []rotations;
                                                rotations = 0;
                                                delete []translations;
                                                translations = 0;
                                        }
                                }
                        }
                }
        }
        if(grp->type == PDBParser::COLLECTION_TYPE)
        {
                assert(!((grp->m_SubGroups).empty()));
                flattenGOA(grp->m_SubGroups.at(0), flatList, collectionData,
numberOfTransformations, rotations, translations, radiusType, level,
colorBySubunits);
        }
        else
        {
                for(int i=0; i < grp->m_SubGroups.size(); i++)
                {
                        flattenGOA(grp->m_SubGroups.at(i), flatList,
collectionData, numberOfTransformations, rotations, translations, radiusType,
level, colorBySubunits);
                }
        }
        int m = grp->m_Atoms.size();
        for(int i=0; i < m; i++)
        {
                PDBParser::Atom* at = grp->m_Atoms.at(i);
                float x = at->m_Position[0];
                float y = at->m_Position[1];
                float z = at->m_Position[2];
                if(numberOfTransformations > 0 && rotations && translations)
                {
                        for(int j=0; j<numberOfTransformations; j++)
                        {
                                double rot[9];
                                rot[0] = rotations[j*9+0];
                                rot[1] = rotations[j*9+1];
                                rot[2] = rotations[j*9+2];
                                rot[3] = rotations[j*9+3];
                                rot[4] = rotations[j*9+4];
                                rot[5] = rotations[j*9+5];
                                rot[6] = rotations[j*9+6];
                                rot[7] = rotations[j*9+7];
                                rot[8] = rotations[j*9+8];
                                float x1=x*rot[0]+y*rot[1]+z*rot[2];
                                float y1=x*rot[3]+y*rot[4]+z*rot[5];
                                float z1=x*rot[6]+y*rot[7]+z*rot[8];
                                double trans[3];
                                trans[0] = translations[j*3+0];
                                trans[1] = translations[j*3+1];
                                trans[2] = translations[j*3+2];
                                x1 += trans[0];
                                y1 += trans[1];
                                z1 += trans[2];
                                PDBParser::Atom* dupAtom = new
PDBParser::Atom(*at); dupAtom->m_Position[0] = x1; dupAtom->m_Position[1] = y1;
                                dupAtom->m_Position[2] = z1;
                                {
                                        // make only one subunit colored!
                                        if(colorBySubunits)
                                        {
                                                srand(j*7);
                                                double t1 = rand();
                                                double t2 = rand();
                                                unsigned int r,g,b,a;
                                                double red   = rand() /
((double)(RAND_MAX)); double green = rand() / ((double)(RAND_MAX)); double blue
= rand() / ((double)(RAND_MAX)); double alpha = 1.0; r = (unsigned int)(red *
255.0); g = (unsigned int)(green * 255.0); b = (unsigned int)(blue * 255.0); a =
(unsigned int)(alpha * 255.0); dupAtom->m_UserDefinedColor = ((r << 24) &
0xff000000) | ((g << 16) & 0xff0000) | ((b << 8) & 0xff00) | (a & 0xff);
                                        }
                                }
                                flatList.push_back(dupAtom);
                        }
                }
                else
                {
                        PDBParser::Atom* dupAtom = new PDBParser::Atom(*at);
                        dupAtom->m_Position[0] = at->m_Position[0];
                        dupAtom->m_Position[1] = at->m_Position[1];
                        dupAtom->m_Position[2] = at->m_Position[2];
                        flatList.push_back(dupAtom);
                }
        }
        return true;

}
*/

void BlurMapsDataManager::getBoundingBox(
    vector<PDBParser::Atom *> &grp, float min[3], float max[3],
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, int probeRadius,
    double blobbiness) {
  int i, j;
  int na = grp.size();
  float maxRad = 0.0;
  if (na == 0) // no atoms
  {
    for (i = 0; i < 3; i++) {
      min[i] = max[i] = 0.0;
    }
    return;
  }
  for (i = 0; i < 3; i++) {
    PDBParser::Atom *atom = grp[0];
    min[i] = atom->m_Position[i];
    max[i] = atom->m_Position[i];
    float radius = atom->getRadius(radiusType);
    radius += probeRadius;
    float tempRad = radius;
    if (blobbiness > 0.0001) {
      tempRad = radius * sqrt(1.0 + log(1e-2f) / blobbiness);
    }
    if (maxRad < tempRad) {
      maxRad = tempRad;
    }
  }
  for (j = 1; j < na; j++) {
    PDBParser::Atom *atom = grp[j];
    for (i = 0; i < 3; i++) {
      if (atom->m_Position[i] < min[i]) {
        min[i] = atom->m_Position[i];
      }
      if (atom->m_Position[i] > max[i]) {
        max[i] = atom->m_Position[i];
      }
    }
    float radius = atom->getRadius(radiusType);
    radius += probeRadius;
    float tempRad = radius;
    if (blobbiness > 0.0001) {
      tempRad = radius * sqrt(1.0 + log(1e-2f) / blobbiness);
    }
    if (maxRad < tempRad) {
      maxRad = tempRad;
    }
  }
  // adjust the bounding box by 1
  for (i = 0; i < 3; i++) {
    min[i] -= 4 * maxRad;
    max[i] += 4 * maxRad;
  }
}

SimpleVolumeData *BlurMapsDataManager::getMultiLevelBlurredVolume(
    PDBParser::GroupOfAtoms *molecule, const char *volFileName, int dim1,
    int dim2, int dim3, PDBParser::GroupOfAtoms::FUNCTIONS densityType,
    bool withColor, double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
    const char *cmapFile, int gap,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
    CCVOpenGLMath::Matrix *transformation) {
  if (!molecule) {
    return 0;
  }

  SimpleVolumeData *vol = new SimpleVolumeData(dim1, dim2, dim3);
  MultilevelBlur *multiLevelBlur = new MultilevelBlur();
  multiLevelBlur->setBlobbyness(blob);
  if (!multiLevelBlur->generateBlurmap(molecule, vol, densityType, cmapFile,
                                       withColor, colorLevel, radiusType,
                                       level)) {
    return 0;
  }
  if (volFileName && strlen(volFileName)) {
    bool ret = false;
    {
      VolumeLoader *volumeLoader = new VolumeLoader();
      ret = volumeLoader->saveFile(volFileName, vol);
      delete volumeLoader;
    }
  }
  delete multiLevelBlur;
  return vol;
}

SimpleVolumeData *BlurMapsDataManager::getVolume(
    PDBParser::GroupOfAtoms *molecule, const char *volFileName, int dim1,
    int dim2, int dim3, PDBParser::GroupOfAtoms::FUNCTIONS densityType,
    bool withColor, double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
    const char *cmapFile, int gap,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
    CCVOpenGLMath::Matrix *transformation) {
  if (!molecule) {
    return 0;
  }
  if (transformation) {
    molecule->transform(*transformation);
  }
  SimpleVolumeData *vol = new SimpleVolumeData(dim1, dim2, dim3);
  GOABlur *goaBlur = new GOABlur();
  goaBlur->setBlobbyness(blob);
  if (!goaBlur->generateBlurmap(molecule, vol, densityType, cmapFile, withColor,
                                colorLevel, radiusType, level)) {
    return 0;
  }
  if (volFileName && strlen(volFileName)) {
    bool ret = false;
    {
      VolumeLoader *volumeLoader = new VolumeLoader();
      ret = volumeLoader->saveFile(volFileName, vol);
      delete volumeLoader;
    }
  }
  delete goaBlur;
  return vol;
}
SimpleVolumeData *BlurMapsDataManager::getVolume(
    PDBParser::GroupOfAtoms *molecule, const char *volFileName, int dim1,
    int dim2, int dim3, float bbmin[3], float bbmax[3],
    PDBParser::GroupOfAtoms::FUNCTIONS densityType, bool withColor, double blob,
    PDBParser::GroupOfAtoms::GOA_TYPE colorLevel, const char *cmapFile, int gap,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
    CCVOpenGLMath::Matrix *transformation) {
  if (!molecule) {
    return 0;
  }
  if (transformation) {
    molecule->transform(*transformation);
  }
  SimpleVolumeData *vol = new SimpleVolumeData(dim1, dim2, dim3);
  GOABlur *goaBlur = new GOABlur();
  goaBlur->setBlobbyness(blob);
  if (!goaBlur->generateBlurmap(molecule, vol, bbmin, bbmax, densityType,
                                cmapFile, withColor, colorLevel, radiusType,
                                level)) {
    return 0;
  }
  if (volFileName && strlen(volFileName)) {
    bool ret = false;
    {
      VolumeLoader *volumeLoader = new VolumeLoader();
      ret = volumeLoader->saveFile(volFileName, vol);
      delete volumeLoader;
    }
  }
  delete goaBlur;
  return vol;
}

#include <iostream>
using namespace std;

SimpleVolumeData *BlurMapsDataManager::getVolumeFixedResolution(
    PDBParser::GroupOfAtoms *molecule, const char *volFileName,
    double resolution, PDBParser::GroupOfAtoms::FUNCTIONS densityType,
    bool withColor, double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
    const char *cmapFile, int gap,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
    CCVOpenGLMath::Matrix *transformation) {
  if (!molecule) {
    return 0;
  }
  if (transformation) {
    molecule->transform(*transformation);
  }

  int dim1, dim2, dim3;
  float min[3];
  float max[3];

  vector<PDBParser::Atom *> atomList;
  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  bool colorBySubunits = false;
  if (colorLevel == PDBParser::GroupOfAtoms::PROTEIN ||
      colorLevel == PDBParser::GroupOfAtoms::COLLECTION) {
    colorBySubunits = true;
  }
  FlattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, level,
             colorBySubunits);
  {
    // printAtomColors(atomList);
  }
  // compute the bounding box
  min[0] = min[1] = min[2] = 0.;
  max[0] = max[1] = max[2] = 0.;
  getBoundingBox(atomList, min, max, radiusType, 0, blob);

  dim1 = (max[0] - min[0]) / resolution + 3;
  dim2 = (max[1] - min[1]) / resolution + 3;
  dim3 = (max[2] - min[2]) / resolution + 3;

  cout << "arand " << dim1 << " " << dim2 << " " << dim3 << endl;

  SimpleVolumeData *vol = new SimpleVolumeData(dim1, dim2, dim3);
  GOABlur *goaBlur = new GOABlur();
  goaBlur->setBlobbyness(blob);
  if (!goaBlur->generateBlurmap(molecule, vol, densityType, cmapFile, withColor,
                                colorLevel, radiusType, level)) {
    return 0;
  }
  if (volFileName && strlen(volFileName)) {
    bool ret = false;
    {
      VolumeLoader *volumeLoader = new VolumeLoader();
      ret = volumeLoader->saveFile(volFileName, vol);
      delete volumeLoader;
    }
  }
  delete goaBlur;
  return vol;
}

SimpleVolumeData *BlurMapsDataManager::getSkinRegionVolume(
    const char *inputFileName, const char *outputFileName, int dim1, int dim2,
    int dim3, double probeRadius,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, int depth) {
  if (!inputFileName || !outputFileName || (dim1 < 2) || (dim1 < 2) ||
      (dim1 < 2) || (probeRadius <= 0) || (depth < 1)) {
    return 0;
  }
  GOALoader *goaLoader = new GOALoader();
  PDBParser::GroupOfAtoms *molecule = goaLoader->loadFile(inputFileName);
  delete goaLoader;
  if (!molecule) {
    return 0;
  }
  SimpleVolumeData *vol = new SimpleVolumeData(dim1, dim2, dim3);
  // SkinRegion2* skin = new SkinRegion2();
  // if( !skin->getSkinRegion( molecule, dim1, dim2, dim3, probeRadius,
  // radiusType, vol, depth ) )
  SignDistanceMolSurface *molSurf = new SignDistanceMolSurface();

  if (!molSurf->getMolecularSurface(molecule, dim1, dim2, dim3, probeRadius,
                                    radiusType, vol)) {
    delete vol;
    vol = 0;
  }

  delete molSurf;
  delete molecule;
  if (outputFileName && strlen(outputFileName)) {
    bool ret = false;
    {
      VolumeLoader *volumeLoader = new VolumeLoader();
      ret = volumeLoader->saveFile(outputFileName, vol);
      delete volumeLoader;
    }
  }
  return vol;
}

bool BlurMapsDataManager::populateSAS(
    const char *inputFileName, const char *outputFileName, int dim1, int dim2,
    int dim3, double probeRadius, double floatingBandWidth,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  if (!inputFileName || !outputFileName || (dim1 < 2) || (dim1 < 2) ||
      (dim1 < 2) || (probeRadius <= 0)) {
    return false;
  }
  PDBParser::GroupOfAtoms *molecule = 0;
  GOALoader *gLoader = new GOALoader();
  molecule = gLoader->loadFile(inputFileName);
  delete gLoader;
  if (!molecule) {
    return 0;
  }
  AtomPopulator *atomPopulator = new AtomPopulator(
      molecule, dim1, dim2, dim3, probeRadius, floatingBandWidth, radiusType);
  bool ret = atomPopulator->populate(outputFileName);
  delete molecule;
  return ret;
}

bool BlurMapsDataManager::populateSASUsingMesh(
    string pqrFileName, string surfaceFileName, string outputFileName,
    double probeRadius, double floatingBandWidth, double clashDistance,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  PDBParser::GroupOfAtoms *molecule = 0;
  GOALoader *gLoader = new GOALoader();
  molecule = gLoader->loadFile(pqrFileName.c_str());
  delete gLoader;
  if (!molecule) {
    return 0;
  }
  SkinGenerator *skinGenerator =
      new SkinGenerator(molecule, surfaceFileName, probeRadius,
                        floatingBandWidth, clashDistance, radiusType);
  bool ret = skinGenerator->populate(outputFileName);
  delete molecule;
  return ret;
}

bool BlurMapsDataManager::getInterfaceAtoms(
    string pqrFileName1, string pqrFileName2, string outputFileName,
    double interfaceWidth, string atomTypes,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  PDBParser::GroupOfAtoms *molecule1 = 0;
  GOALoader *gLoader1 = new GOALoader();
  molecule1 = gLoader1->loadFile(pqrFileName1.c_str());
  delete gLoader1;
  if (!molecule1) {
    return 0;
  }

  PDBParser::GroupOfAtoms *molecule2 = 0;
  GOALoader *gLoader2 = new GOALoader();
  molecule2 = gLoader2->loadFile(pqrFileName2.c_str());
  delete gLoader2;
  if (!molecule2) {
    return 0;
  }

  InterfaceAtomDetector *interfaceAtomDetector =
      new InterfaceAtomDetector(molecule1, molecule2, outputFileName,
                                interfaceWidth, atomTypes, radiusType);

  bool ret = interfaceAtomDetector->detect();

  delete molecule1;
  delete molecule2;
  return ret;
}

bool BlurMapsDataManager::getInterfaceSurface(string surfaceFileName1,
                                              string surfaceFileName2,
                                              double interfaceWidth,
                                              bool colored, double &area1,
                                              double &area2) {
  InterfaceSurface *interfaceSurface = new InterfaceSurface(
      surfaceFileName1, surfaceFileName2, interfaceWidth, colored);

  double a1 = 0, a2 = 0;

  bool ret = interfaceSurface->detect(a1, a2);

  if (ret) {
    cout << "Area of the interface of the first molecule = " << a1 << " A^2"
         << endl;
    cout << "Area of the interface of the second molecule = " << a2 << " A^2"
         << endl;

    area1 = a1;
    area2 = a2;
  }

  return ret;
}

bool BlurMapsDataManager::generateF2d(
    string pqrFileName, string xyzFileName, string outputFileName,
    bool receptor, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
    double skinAtomRad) {
  PDBParser::GroupOfAtoms *molecule = 0;
  GOALoader *gLoader = new GOALoader();
  molecule = gLoader->loadFile(pqrFileName.c_str());
  delete gLoader;
  if (!molecule) {
    cout << "Molecule not found" << endl;
    return 0;
  }
  F2DGenerator *f2dGenerator = new F2DGenerator(
      molecule, xyzFileName, outputFileName, receptor, radiusType, skinAtomRad);
  bool ret = f2dGenerator->generate();
  delete molecule;
  return ret;
}

double BlurMapsDataManager::getDistanceOfVoxel(int i, int j, int k, int width,
                                               int height, int depth) {
  return (sqrt(((i)-width / 2) * ((i)-width / 2) +
               ((j)-height / 2) * ((j)-height / 2) +
               ((k)-depth / 2) * ((k)-depth / 2)));
}

bool BlurMapsDataManager::outputCurvFiles(
    double *HandK, Geometry *geometry, const char *outputMeanRawSurfaceFileName,
    const char *outputGaussianRawSurfaceFileName,
    const char *curvatureFileName) {
  if (!HandK || !geometry || !outputMeanRawSurfaceFileName ||
      !outputGaussianRawSurfaceFileName || !curvatureFileName) {
    return false;
  }
  if (geometry->m_NumTriVerts < 1) {
    return false;
  }
  double h, k;
  float r, g, b;
  int i;
  delete[] geometry->m_TriVertColorsTransparent;
  geometry->m_TriVertColorsTransparent = new float[geometry->m_NumTriVerts * 4];
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    h = HandK[i * 2 + 0];
    k = HandK[i * 2 + 1];
    r = g = b = 1.0;
    if (h < 0) {
      r = 1.0 + h * 5;
      b = 1.0 + h * 5;
      if (h < -0.2) {
        r = 0.0;
        b = 0.0;
      }
    }
    if (h > 0) {
      g = 1.0 - h * 5;
      b = 1.0 - h * 5;
      if (h > 0.2) {
        g = 0.0;
        b = 0.0;
      }
    }
    geometry->m_TriVertColorsTransparent[4 * i + 0] = r;
    geometry->m_TriVertColorsTransparent[4 * i + 1] = g;
    geometry->m_TriVertColorsTransparent[4 * i + 2] = b;
  }
  GeometryLoader *geometryLoader = new GeometryLoader();
  bool savedMeanCurvatureFile = geometryLoader->saveFile(
      outputMeanRawSurfaceFileName, "Rawnc files (*.rawnc)", geometry);
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    h = HandK[i * 2 + 0];
    k = HandK[i * 2 + 1];
    r = g = b = 1.0;
    if (k < 0) {
      double k2 = -1 * sqrt(-k);
      r = 1.0 + k2 * 5;
      b = 1.0 + k2 * 5;
      if (k2 < -0.2) {
        r = 0.0;
        b = 0.0;
      }
    }
    if (k > 0) {
      double k2 = sqrt(k);
      g = 1.0 - k2 * 5;
      b = 1.0 - k2 * 5;
      if (k2 > 0.2) {
        g = 0.0;
        b = 0.0;
      }
    }
    geometry->m_TriVertColorsTransparent[4 * i + 0] = r;
    geometry->m_TriVertColorsTransparent[4 * i + 1] = g;
    geometry->m_TriVertColorsTransparent[4 * i + 2] = b;
  }
  bool savedGaussianCurvatureFile = geometryLoader->saveFile(
      outputGaussianRawSurfaceFileName, "Rawnc files (*.rawnc)", geometry);
  // print out the either radius > 10A thing
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    h = HandK[i * 2 + 0];
    k = HandK[i * 2 + 1];
    double k1, k2;
    double d = h * h - k;
    if (d < 0) {
      d = 0;
    }
    k1 = h - sqrt(d);
    k2 = h + sqrt(d);
    double r1 = -1.0 / k1;
    double r2 = -1.0 / k2;
    r = g = b = 0.0;
    double minRadius = 8;
    if (r1 < -minRadius || r2 < -minRadius || r1 > minRadius ||
        r2 > minRadius) {
      g = 1.0;
    } else {
      r = 1.0;
    }
    geometry->m_TriVertColorsTransparent[4 * i + 0] = r;
    geometry->m_TriVertColorsTransparent[4 * i + 1] = g;
    geometry->m_TriVertColorsTransparent[4 * i + 2] = b;
  }
  char newname[256];
  strcpy(newname, "curv_either_more_than_10.rawnc");
  bool savedEither10CurvatureFile =
      geometryLoader->saveFile(newname, "Rawnc files (*.rawnc)", geometry);
  // print out the both radius > 10A thing
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    h = HandK[i * 2 + 0];
    k = HandK[i * 2 + 1];
    double k1, k2;
    double d = h * h - k;
    if (d < 0) {
      d = 0;
    }
    k1 = h - sqrt(d);
    k2 = h + sqrt(d);
    double r1 = -1.0 / k1;
    double r2 = -1.0 / k2;
    r = g = b = 0.0;
    if (r1 < -10 && r2 < -10) {
      g = 1.0;
    } else {
      r = 1.0;
    }
    geometry->m_TriVertColorsTransparent[3 * i + 0] = r;
    geometry->m_TriVertColorsTransparent[3 * i + 1] = g;
    geometry->m_TriVertColorsTransparent[3 * i + 2] = b;
  }
  strcpy(newname, "curv_both_more_than_10.rawnc");
  bool savedBoth10CurvatureFile =
      geometryLoader->saveFile(newname, "Rawnc files (*.rawnc)", geometry);
  FILE *fp = fopen("radii.txt", "w");
  // print out colorful thing
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    h = HandK[i * 2 + 0];
    k = HandK[i * 2 + 1];
    double k1, k2;
    double d = h * h - k;
    if (d < 0) {
      d = 0;
    }
    k1 = h - sqrt(d);
    k2 = h + sqrt(d);
    double r1 = -1.0 / k1;
    double r2 = -1.0 / k2;
    fprintf(fp, "%lf %lf\n", r1, r2);
    r = g = b = 0.0;
    if (r1 < -10) {
      r1 = -10;
    }
    if (r1 > 10) {
      r1 = 10;
    }
    r = (10 - r1) / 20;
    g = 1 - r;
    if (r < 0) {
      r = 0;
    }
    if (g < 0) {
      g = 0;
    }
    if (r > 1) {
      r = 1;
    }
    if (g > 1) {
      g = 1;
    }
    geometry->m_TriVertColorsTransparent[3 * i + 0] = r;
    geometry->m_TriVertColorsTransparent[3 * i + 1] = g;
    geometry->m_TriVertColorsTransparent[3 * i + 2] = b;
  }
  fclose(fp);
  strcpy(newname, "curv_colorful.rawnc");
  bool savedcolorfulCurvatureFile =
      geometryLoader->saveFile(newname, "Rawnc files (*.rawnc)", geometry);
  delete geometryLoader;
  // print out the curvatures
  {
    FILE *fpC = 0;
    fpC = fopen(curvatureFileName, "w");
    if (!fpC) {
      return false;
    }
    fprintf(fpC, "%d\n", geometry->m_NumTriVerts);
    for (i = 0; i < geometry->m_NumTriVerts; i++) {
      fprintf(fpC, "%lf %lf\n", HandK[2 * i + 0], HandK[2 * i + 1]);
    }
    fclose(fpC);
  }
  if (savedMeanCurvatureFile && savedGaussianCurvatureFile) {
    return true;
  }
  return false;
}

bool BlurMapsDataManager::getCurvaturesFromIsocontourFile(
    const char *proteinFileName, int dim1, int dim2, int dim3, double blob,
    const char *inputRawSurfaceFileName,
    const char *outputMeanRawSurfaceFileName,
    const char *outputGaussianRawSurfaceFileName, const char *curvatureFileName,
    int numberOfGridDivisions, double maxFunctionError,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level) {
  Geometry *geometry = 0;
  GeometryLoader *geometryLoader = new GeometryLoader();
  geometry = geometryLoader->loadFile(inputRawSurfaceFileName);
  delete geometryLoader;
  if (geometry == 0) {
    return false;
  }
  return getCurvatures(proteinFileName, dim1, dim2, dim3, blob, geometry,
                       outputMeanRawSurfaceFileName,
                       outputGaussianRawSurfaceFileName, curvatureFileName,
                       numberOfGridDivisions, maxFunctionError, radiusType,
                       level);
}

/*************************************************************************/
/*                                                                       */
/*  The type of file being handled here is triangle meshes. The          */
/*  Geometry data structure does support quads, but this function is not */
/*  yet general enough. This should change in the future.                */
/*                                                                       */
/*************************************************************************/
bool BlurMapsDataManager::getCurvatures(
    const char *proteinFileName, int dim1, int dim2, int dim3, double blob,
    Geometry *geometry, const char *outputMeanRawSurfaceFileName,
    const char *outputGaussianRawSurfaceFileName, const char *curvatureFileName,
    int numberOfGridDivisions, double maxFunctionError,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level) {
  if (geometry == 0) {
    return false;
  }
  int numberOfPoints = geometry->m_NumTriVerts;
  float *points = geometry->m_TriVerts;
  if (points == 0) {
    return false;
  }
  if (!proteinFileName) {
    return false;
  }
  GOALoader *gLoader = new GOALoader();
  PDBParser::GroupOfAtoms *molecule = gLoader->loadFile(proteinFileName);
  delete gLoader;
  int numberOfAtoms = 0;
  molecule->getNumberOfAtomsRecursive(&numberOfAtoms);
  if (numberOfAtoms < 1) {
    delete molecule;
    return false;
  }
  double *posRadius = new double[numberOfAtoms * 4];
  int curPosition = 0;
  molecule->getAtomsRecursive(posRadius, numberOfAtoms * 4, &curPosition);
  delete molecule;
  // compute the curvatures
  double *HandK = new double[numberOfPoints * 2];
  double *normals = new double[numberOfPoints * 3];
  double *k1Vec = new double[numberOfPoints * 3];
  double *k2Vec = new double[numberOfPoints * 3];
  SumOfGaussiansCurvature *sumOfGaussiansCurvature =
      new SumOfGaussiansCurvature(
          numberOfAtoms, posRadius, numberOfGridDivisions, maxFunctionError,
          blob, numberOfPoints, points, HandK, normals, k1Vec, k2Vec);
  if (!sumOfGaussiansCurvature->initialize()) {
    delete[] HandK;
    HandK = 0;
    delete[] normals;
    normals = 0;
    delete[] k1Vec;
    k1Vec = 0;
    delete[] k2Vec;
    k2Vec = 0;
    delete[] posRadius;
    posRadius = 0;
    delete sumOfGaussiansCurvature;
    sumOfGaussiansCurvature = 0;
    return false;
  }
  if (!sumOfGaussiansCurvature->getCurvatures()) {
    delete[] HandK;
    HandK = 0;
    delete[] normals;
    normals = 0;
    delete[] k1Vec;
    k1Vec = 0;
    delete[] k2Vec;
    k2Vec = 0;
    delete[] posRadius;
    posRadius = 0;
    delete sumOfGaussiansCurvature;
    sumOfGaussiansCurvature = 0;
    return false;
  }
  char curvfilename[256];
  strcpy(curvfilename, curvatureFileName);
  strcat(curvfilename, ".curv");
  if (!sumOfGaussiansCurvature->write(curvfilename)) {
    delete[] HandK;
    HandK = 0;
    delete[] normals;
    normals = 0;
    delete[] k1Vec;
    k1Vec = 0;
    delete[] k2Vec;
    k2Vec = 0;
    delete[] posRadius;
    posRadius = 0;
    delete sumOfGaussiansCurvature;
    sumOfGaussiansCurvature = 0;
    return false;
  }
  if (!outputCurvFiles(HandK, geometry, outputMeanRawSurfaceFileName,
                       outputGaussianRawSurfaceFileName, curvatureFileName)) {
    delete[] HandK;
    HandK = 0;
    delete[] normals;
    normals = 0;
    delete[] k1Vec;
    k1Vec = 0;
    delete[] k2Vec;
    k2Vec = 0;
    delete[] posRadius;
    posRadius = 0;
    delete sumOfGaussiansCurvature;
    sumOfGaussiansCurvature = 0;
    return false;
  }
  delete[] HandK;
  HandK = 0;
  delete[] normals;
  normals = 0;
  delete[] k1Vec;
  k1Vec = 0;
  delete[] k2Vec;
  k2Vec = 0;
  delete[] posRadius;
  posRadius = 0;
  delete sumOfGaussiansCurvature;
  sumOfGaussiansCurvature = 0;
  return true;
}

Geometry *BlurMapsDataManager::getSurfaceFromPDB(
    PDBParser::GroupOfAtoms *molecule, int dim1, int dim2, int dim3,
    double isovalue, double blobbiness,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, int level) {
  ///////// create the volume ////////////////////
  SimpleVolumeData *sData = 0;
  {
    if (!molecule) {
      return 0;
    }
    sData = BlurMapsDataManager::getVolume(
        molecule, "", dim1, dim2, dim3,
        PDBParser::GroupOfAtoms::ELECTRON_DENSITY, false, blobbiness,
        PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
    if (!sData) {
      return 0;
    }
  }
  ////// create the isosurface /////////////////////
  Geometry *geometry =
      SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
  delete sData;
  return geometry;
}
