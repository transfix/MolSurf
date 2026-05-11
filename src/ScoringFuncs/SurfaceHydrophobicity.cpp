
#include <ScoringFuncs/SurfaceHydrophobicity.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/FlattenGOA.h>
#include <GOAFileTypes/GOALoader.h>

SurfaceHydrophobicity::SurfaceHydrophobicity(string pqrFileName,
                                             string quadFileName) {
  dataReady = false;
  if (!readPQR(pqrFileName))
    printf("could not read pqr\n");
  else if (!readQuad(quadFileName))
    printf("could not read quad\n");
  else {
    computeXlateForPG();
    buildPG();
    dataReady = true;
  }
}

SurfaceHydrophobicity::~SurfaceHydrophobicity() {
  if (dataReady) {
    if (pgPoints)
      delete pgPoints;
    atms.clear();
    quadPoints.clear();
  }
}

double SurfaceHydrophobicity::getSurfaceHydrophobicity() {
  preprocessElementInformationTable(&hydroIndex, &nHydroIndex);

  double totalSurfaceHydrophobicity = 0.0;
  double ip;

  int numAtoms = atms.size();

  vector<DPG::Point *> quadsInRange;

  for (int i = 0; i < atms.size(); i++) {
    ip = getHydrophobicity(atms[i]->getName(), atms[i]->getResidueName(),
                           hydroIndex, nHydroIndex, true, true);

    pgPoints->range(atms[i], atms[i]->getRadius(), quadsInRange);

    int numQuadsInRange = quadsInRange.size();

    for (int j = 0; j < numQuadsInRange; j++) {
      totalSurfaceHydrophobicity +=
          ((DPG::WeightedPoint *)quadsInRange[j])->getWeight() * ip;
    }

    quadsInRange.clear();
  }

  return totalSurfaceHydrophobicity;
}

int SurfaceHydrophobicity::strcmp_nospace(const char *s1, const char *s2) {
  int i = 0, j = 0;

  while (s1[i] && s2[j]) {
    if (!isspace(s1[i]) && !isspace(s2[j])) {
      if (s1[i] == s2[j]) {
        i++;
        j++;
      } else
        return ((int)(s1[i] - s2[j]));
    }

    if (isspace(s1[i]))
      i++;
    if (isspace(s2[j]))
      j++;
  }

  while (isspace(s1[i]))
    i++;
  while (isspace(s2[j]))
    j++;

  return ((int)(s1[i] - s2[j]));
}

float SurfaceHydrophobicity::getHydrophobicity(const char *atomName,
                                               const char *residueName,
                                               int *index, int nIndex,
                                               bool useInterfacePropensity,
                                               bool perResidueHydrophobicity) {
  for (int i = 0; i < nIndex; i++) {
    int j = index[i];

    if (!strcmp_nospace(residueName, PDBParser::elementTable[j].residueName)) {
      int k =
          (i == nIndex - 1) ? PDBParser::MAX_BIOCHEM_ELEMENTS : index[i + 1];

      while (j < k) {
        if (!strcmp_nospace(atomName, PDBParser::elementTable[j].atomName)) {
          if (useInterfacePropensity)
            return ((float)PDBParser::elementTable[j].interfacePropensity);

          if (perResidueHydrophobicity)
            return ((float)PDBParser::elementTable[j].perResidueHydrophobicity);
          else
            return ((float)PDBParser::elementTable[j].hydrophobicity);
        }
        j++;
      }
    }
  }

  return 0;
}

void SurfaceHydrophobicity::preprocessElementInformationTable(int **index,
                                                              int *nIndex) {
  int c = 1;
  int n = sizeof(PDBParser::elementTable) / sizeof(PDBParser::elementTable[0]);

  for (int i = 1; i < n; i++)
    if (strcmp(PDBParser::elementTable[i].residueName,
               PDBParser::elementTable[i - 1].residueName))
      c++;

  (*index) = new int[c];

  *nIndex = c;

  c = 0;
  (*index)[0] = 0;

  for (int i = 1; i < n; i++)
    if (strcmp(PDBParser::elementTable[i].residueName,
               PDBParser::elementTable[i - 1].residueName))
      (*index)[++c] = i;
}

void SurfaceHydrophobicity::computeXlateForPG() {
  double minXYZ;
  DPG::Point *p;
  double x, y, z;

  double minX, minY, minZ;
  double maxX, maxY, maxZ;

  double xDim, yDim, zDim;
  double maxDim;

  int numPoints = quadPoints.size();

  for (int i = 0; i < numPoints; i++) {
    x = quadPoints[i]->getX();
    y = quadPoints[i]->getY();
    z = quadPoints[i]->getZ();

    if (i == 0) {
      minX = maxX = x;
      minY = maxY = y;
      minZ = maxZ = z;
    }

    if (x < minX)
      minX = x;
    else if (x > maxX)
      maxX = x;

    if (y < minY)
      minY = y;
    else if (y > maxY)
      maxY = y;

    if (z < minZ)
      minZ = z;
    else if (z > maxZ)
      maxZ = z;
  }

  minXYZ = minX;
  if (minXYZ > minY)
    minXYZ = minY;
  if (minXYZ > minZ)
    minXYZ = minZ;

  xDim = maxDim = maxX - minX;
  yDim = maxY - minY;
  zDim = maxZ - minZ;

  if (maxDim < yDim)
    maxDim = yDim;
  if (maxDim < zDim)
    maxDim = zDim;

  xlatePG = 0.0;

  if (minXYZ < 0)
    xlatePG = -2.0 * minXYZ;
}

void SurfaceHydrophobicity::buildPG() {
  computeXlateForPG();

  pgPoints = new DPG::PG(3.0, xlatePG, 3.0);

  int numPoints = quadPoints.size();

  for (int i = 0; i < numPoints; i++) {
    pgPoints->addPoint(quadPoints[i]);
  }
}

bool SurfaceHydrophobicity::readPQR(string pqrFileName) {
  PDBParser::GroupOfAtoms *molecule = 0;
  GOALoader *gLoader = new GOALoader();
  molecule = gLoader->loadFile(pqrFileName.c_str());
  delete gLoader;

  if (!molecule) {
    return false;
  }

  PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
  PDBParser::GroupOfAtoms::intToRadiusType(&radiusType, 0);

  PDBParser::CollectionData *collectionData = 0;
  if (molecule->type == PDBParser::COLLECTION_TYPE) {
    collectionData = molecule->m_CollectionData;
  }
  PDBParser::FlattenGOA(molecule, atms, collectionData, 0, 0, 0, radiusType,
                        PDBParser::ATOM_TYPE, false);

  return true;
}

bool SurfaceHydrophobicity::readQuad(string quadFileName) {
  FILE *fp;

  fp = fopen(quadFileName.c_str(), "rt");

  if (fp == NULL) {
    printf("Failed to open quadrature points file (%s)!\n",
           quadFileName.c_str());
    return false;
  }

  double x, y, z, w;
  double nx, ny, nz;

  while (fscanf(fp, (char *)"%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny,
                &nz, &w) == 7) {
    DPG::WeightedPoint *wp = new DPG::WeightedPoint(x, y, z, w);
    quadPoints.push_back(wp);
  }

  fclose(fp);

  return true;
}
