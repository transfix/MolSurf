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
#include <Blurmaps/GOABlur.h>
#include <PDBParser/FlattenGOA.h>

/*  Get the current time in seconds as a double value  */
double getTime() {
#ifdef _WIN32
  time_t ltime;
  time(&ltime);
  return (double)ltime;
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)(t.tv_sec + 1e-6 * t.tv_usec);
#endif
}

GOABlur::GOABlur() {
  mBlobbyness = -2.3442; // this value is from the following paper
  //"Evaluation of Protein Docking Predictions Using Hex 3.1 in CAPRI Rounds 1
  // and 2", D.W. Ritchie (2003) PROTEINS: Struct. Funct. Genet. 52(1), 98-106.
}

GOABlur::~GOABlur() {}

void printAtomColors(vector<PDBParser::Atom *> atomList) {
  FILE *fp = fopen("colors.txt", "w");
  int i;
  for (i = 0; i < atomList.size(); i++) {
    float r, g, b, a;
    atomList[i]->getColor(&r, &g, &b, &a);
    fprintf(fp, "[%c -> %2.1f %2.1f %2.1f %2.1f]\n", atomList[i]->getChain(), r,
            g, b, a);
  }
  fclose(fp);
}

#include <iostream>
using namespace std;

bool GOABlur::generateBlurmap(PDBParser::GroupOfAtoms *molecule,
                              SimpleVolumeData *vol,
                              PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                              const char *cmapFile, bool withColor,
                              PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                              PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                              unsigned int level) {
  float *density = NULL, min[3], max[3];
  unsigned int dim[3];
  vector<PDBParser::Atom *> atomList;
  // find out how big the volume is that we're making
  dim[0] = vol->getWidth();
  dim[1] = vol->getHeight();
#ifdef USE_MPI
  dim[2] = vol->parallelGetDepth();
#else
  dim[2] = vol->getDepth();
#endif
  // do the colormapping
  //
  //    we do one of the following:
  //    1. Use a colormap.
  //    2. Color by ATOM, RESIDUE, SS, CHAIN
  //    3. Color by PROTEIN, COLLECTION. We need to do this when we flatten out
  //    the molecule!
  /*
          {
                  PDBParser::GOAColor* goaColor = new PDBParser::GOAColor();
                  if( cmapFile && strlen(cmapFile))
                  {
                          goaColor->parseColormapFile((char *)cmapFile);
                          goaColor->ApplyColormap(molecule);
                  }
                  else
                          goaColor->ApplyColorsByGOALevel(molecule, colorLevel);
                  delete goaColor;
          }
  */
  // "flatten" the GOA
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
  BlurMapsDataManager::getBoundingBox(atomList, min, max, radiusType, 0,
                                      mBlobbyness);
  //     BlurMapsDataManager::getBoundingBox(atomList, min, max,
  //     radiusType, 1.4, mBlobbyness);  // Wenqi
  vol->setDimensions(dim);
  // allocate memory
  density = new float[dim[0] * dim[1] * dim[2]];
  // you want color with that?
  if (withColor) {
    unsigned char *red = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *green = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *blue = new unsigned char[dim[0] * dim[1] * dim[2]];
    vol->setNumberOfVariables(5);
    vol->setData(0, red);
    vol->setType(0, SimpleVolumeData::UCHAR);
    vol->setName(0, "red");
    vol->setData(1, green);
    vol->setType(1, SimpleVolumeData::UCHAR);
    vol->setName(1, "green");
    vol->setData(2, blue);
    vol->setType(2, SimpleVolumeData::UCHAR);
    vol->setName(2, "blue");
    vol->setData(3, density);
    vol->setType(3, SimpleVolumeData::FLOAT);
    vol->setName(3, "density");
    // a special max density variable for determining colors where atoms
    // overlap
    vol->setData(4, new float[dim[0] * dim[1] * dim[2]]);
    vol->setType(4, SimpleVolumeData::FLOAT);
    vol->setName(4, "max density");
  } else {
    vol->setNumberOfVariables(1);
    vol->setData(0, density);
    vol->setType(0, SimpleVolumeData::FLOAT);
    vol->setName(0, "density");
  }

  // set the min and max
  vol->setMinExtent(min);
  vol->setMaxExtent(max);
  // zero all the data
  vol->setVariablesToZero();
  // blur
  blurAtoms(atomList, vol, densityType, withColor, radiusType);
  {
    int n = atomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete atomList[i];
    }
    atomList.clear();
  }
  // get rid of the maximum density variable
  if (withColor) {
    vol->removeVariable(4);
  }
  return true;
}

bool GOABlur::generateBlurmap(PDBParser::GroupOfAtoms *molecule,
                              SimpleVolumeData *vol, float bbmin[3],
                              float bbmax[3],
                              PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                              const char *cmapFile, bool withColor,
                              PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                              PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                              unsigned int level) {
  float *density = NULL; //, min[3],max[3];
  unsigned int dim[3];
  vector<PDBParser::Atom *> atomList;
  // find out how big the volume is that we're making
  dim[0] = vol->getWidth();
  dim[1] = vol->getHeight();
#ifdef USE_MPI
  dim[2] = vol->parallelGetDepth();
#else
  dim[2] = vol->getDepth();
#endif
  // do the colormapping
  //
  //    we do one of the following:
  //    1. Use a colormap.
  //    2. Color by ATOM, RESIDUE, SS, CHAIN
  //    3. Color by PROTEIN, COLLECTION. We need to do this when we flatten out
  //    the molecule!
  /*
          {
                  PDBParser::GOAColor* goaColor = new PDBParser::GOAColor();
                  if( cmapFile && strlen(cmapFile))
                  {
                          goaColor->parseColormapFile((char *)cmapFile);
                          goaColor->ApplyColormap(molecule);
                  }
                  else
                          goaColor->ApplyColorsByGOALevel(molecule, colorLevel);
                  delete goaColor;
          }
  */
  // "flatten" the GOA
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
  //      min[0] = min[1] = min[2] = 0.;
  //      max[0] = max[1] = max[2] = 0.;
  //      BlurMapsDataManager::getBoundingBox(atomList, min, max, radiusType, 0,
  //      mBlobbyness);
  //     BlurMapsDataManager::getBoundingBox(atomList, min, max,
  //     radiusType, 1.4, mBlobbyness);  // Wenqi
  vol->setDimensions(dim);
  // allocate memory
  density = new float[dim[0] * dim[1] * dim[2]];
  // you want color with that?
  if (withColor) {
    unsigned char *red = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *green = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *blue = new unsigned char[dim[0] * dim[1] * dim[2]];
    vol->setNumberOfVariables(5);
    vol->setData(0, red);
    vol->setType(0, SimpleVolumeData::UCHAR);
    vol->setName(0, "red");
    vol->setData(1, green);
    vol->setType(1, SimpleVolumeData::UCHAR);
    vol->setName(1, "green");
    vol->setData(2, blue);
    vol->setType(2, SimpleVolumeData::UCHAR);
    vol->setName(2, "blue");
    vol->setData(3, density);
    vol->setType(3, SimpleVolumeData::FLOAT);
    vol->setName(3, "density");
    // a special max density variable for determining colors where atoms
    // overlap
    vol->setData(4, new float[dim[0] * dim[1] * dim[2]]);
    vol->setType(4, SimpleVolumeData::FLOAT);
    vol->setName(4, "max density");
  } else {
    vol->setNumberOfVariables(1);
    vol->setData(0, density);
    vol->setType(0, SimpleVolumeData::FLOAT);
    vol->setName(0, "density");
  }

  // set the min and max
  vol->setMinExtent(bbmin);
  vol->setMaxExtent(bbmax);
  // zero all the data
  vol->setVariablesToZero();
  // blur
  blurAtoms(atomList, vol, densityType, withColor, radiusType);
  {
    int n = atomList.size();
    int i;
    for (i = 0; i < n; i++) {
      delete atomList[i];
    }
    atomList.clear();
  }
  // get rid of the maximum density variable
  if (withColor) {
    vol->removeVariable(4);
  }
  return true;
}

bool GOABlur::generateBlurmap2(std::vector<PDBParser::Atom *> atomList,
                               SimpleVolumeData *vol,
                               PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                               const char *cmapFile, bool withColor,
                               PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
                               PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                               unsigned int level) {
  float *density = NULL, min[3], max[3];
  unsigned int dim[3];
  // find out how big the volume is that we're making
  dim[0] = vol->getWidth();
  dim[1] = vol->getHeight();
#ifdef USE_MPI
  dim[2] = vol->parallelGetDepth();
#else
  dim[2] = vol->getDepth();
#endif
  // colormapping should already be done to the atom list
  // compute the bounding box
  min[0] = min[1] = min[2] = 0.;
  max[0] = max[1] = max[2] = 0.;
  BlurMapsDataManager::getBoundingBox(atomList, min, max, radiusType, 0,
                                      mBlobbyness);
  vol->setDimensions(dim);
  // allocate memory
  density = new float[dim[0] * dim[1] * dim[2]];
  // you want color with that?
  if (withColor) {
    unsigned char *red = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *green = new unsigned char[dim[0] * dim[1] * dim[2]],
                  *blue = new unsigned char[dim[0] * dim[1] * dim[2]];
    vol->setNumberOfVariables(5);
    vol->setData(0, red);
    vol->setType(0, SimpleVolumeData::UCHAR);
    vol->setName(0, "red");
    vol->setData(1, green);
    vol->setType(1, SimpleVolumeData::UCHAR);
    vol->setName(1, "green");
    vol->setData(2, blue);
    vol->setType(2, SimpleVolumeData::UCHAR);
    vol->setName(2, "blue");
    vol->setData(3, density);
    vol->setType(3, SimpleVolumeData::FLOAT);
    vol->setName(3, "density");
    // a special max density variable for determining colors where atoms
    // overlap
    vol->setData(4, new float[dim[0] * dim[1] * dim[2]]);
    vol->setType(4, SimpleVolumeData::FLOAT);
    vol->setName(4, "max density");
  } else {
    vol->setNumberOfVariables(1);
    vol->setData(0, density);
    vol->setType(0, SimpleVolumeData::FLOAT);
    vol->setName(0, "density");
  }
  // set the min and max
  vol->setMinExtent(min);
  vol->setMaxExtent(max);
  // zero all the data
  vol->setVariablesToZero();
  // blur
  blurAtoms(atomList, vol, densityType, withColor, radiusType);
  // get rid of the maximum density variable
  if (withColor) {
    vol->removeVariable(4);
  }
  return true;
}

double getElecPotCutoff(double eps) { return 1.0 / (80.0 * eps); }

double GOABlur::getMaxRadiusOfAtomPotential(
    PDBParser::Atom *atom, double blobbiness,
    PDBParser::GroupOfAtoms::FUNCTIONS potentialType, double eps) {
  double maxRad = 0;
  switch (potentialType) {
  case PDBParser::GroupOfAtoms::ELECTRON_DENSITY:
    return atom->getRadius() * sqrt(1.0 + log(eps) / mBlobbyness);
  case PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY:
    return atom->getRadius() * sqrt(1.0 + log(eps) / mBlobbyness);
  case PDBParser::GroupOfAtoms::PER_RESIDUE_HYDROPHOBICITY:
    return atom->getRadius() * sqrt(1.0 + log(eps) / mBlobbyness);
  case PDBParser::GroupOfAtoms::ELECTROSTATIC_POTENTIAL:
    return getElecPotCutoff(eps);
  default:
    return 0;
  }
}

//
////////////////////////////////////////////////////////////////////////////////
// Private Functions
//
//

void GOABlur::blurAtoms(vector<PDBParser::Atom *> &atoms, SimpleVolumeData *vol,
                        PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                        bool withColor,
                        PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  float *dens = 0, *maxDens = 0, orig[3], span[3];
  unsigned char *red = 0, *green = 0, *blue = 0;
  unsigned int dim[3], gdim[3], zoff = 0, i;
#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  // assign data pointers
  if (withColor) {
    dens = (float *)vol->getData(3);
    maxDens = (float *)vol->getData(4);
    red = (unsigned char *)vol->getData(0);
    green = (unsigned char *)vol->getData(1);
    blue = (unsigned char *)vol->getData(2);
  } else {
    dens = (float *)vol->getData(0);
  }
  // get dimension (global and local),
  gdim[0] = dim[0] = vol->getWidth();
  gdim[1] = dim[1] = vol->getHeight();
  gdim[2] = dim[2] = vol->getDepth();
#ifdef USE_MPI
  dim[2] = vol->parallelGetDepth();
  zoff = vol->parallelGetZOffset();
#endif
  // origin,
  orig[0] = vol->getMinX();
  orig[1] = vol->getMinY();
  orig[2] = vol->getMinZ();
  // and span
  span[0] = vol->getSpanX();
  span[1] = vol->getSpanY();
  span[2] = vol->getSpanZ();
  double t1 = getTime();
  // Wenqi
  double rmin = 99999.0;
  for (i = 0; i < atoms.size(); i++) {
    if (rmin > atoms[i]->getRadius()) {
      rmin = atoms[i]->getRadius();
    }
  }
  /* 2010-1-6: arand - removed this blobbyness calculation below
  mBlobbyness = -5.0*log(4)/rmin;
  */
  // printf("size = %lu, min = %f, blobby = %f\n", atoms.size(), rmin,
  // mBlobbyness);
  //  main blurring loop
  for (i = 0; i < atoms.size(); i++) {
    PDBParser::Atom *at = atoms[i];
    float radius = at->getRadius(radiusType);
    //            radius += 1.4;  // Wenqi
    // printf("Radius while blurring: %f\n", radius );
    double hydrophobe = 0;
    double charge = 0;
    double c[3], maxRad;
    const float eps = 1.0e-3f;
    unsigned int amax[3], amin[3], j, k, l;
    // skip this atom if it has an alpha value of 0
    if (densityType == PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY) {
      hydrophobe = at->getPerAtomHydrophobicity();
    }
    if (densityType == PDBParser::GroupOfAtoms::PER_RESIDUE_HYDROPHOBICITY) {
      hydrophobe = at->getPerResidueHydrophobicity();
    }
    charge = at->getCharge();
    maxRad = getMaxRadiusOfAtomPotential(at, mBlobbyness, densityType, eps);
    // compute the dataset coordinates of the atom's center
    for (j = 0; j < 3; j++) {
      c[j] = (at->m_Position[j] - orig[j]) / span[j];
      c[j] = ((c[j] - floor(c[j])) >= 0.5) ? ceil(c[j]) : floor(c[j]);
    }
    // then compute the bounding box of the atom (maxRad^3)
    for (j = 0; j < 3; j++) {
      int tempAMin = (int)(c[j] - (maxRad / span[j]) - 1);
      tempAMin = (tempAMin < 0) ? 0 : tempAMin;
      // amin[j] = (int)(c[j] - (maxRad / span[j]) - 1);
      // amin[j] = (amin[j] < 0) ? 0 : amin[j];
      amin[j] = tempAMin;
      int tempAMax = (int)(c[j] + (maxRad / span[j]) + 1);
      tempAMax = (tempAMax > gdim[j]) ? gdim[j] : tempAMax;
      // amax[j] = (int)(c[j] + (maxRad / span[j]) + 1);
      // amax[j] = (amax[j] > gdim[j]) ? gdim[j] : amax[j];
      amax[j] = tempAMax;
    }
#ifdef USE_MPI
    // (parallel specific)
    bool intersects = false;
    // find out if the atom's bounding cube intersects our portion
    // of the volume
    if (amin[2] >= zoff && amin[2] < (zoff + dim[2])) {
      intersects = true;
      if (amax[2] >= (zoff + dim[2])) {
        amax[2] = zoff + dim[2];
      }
    } else if (amax[2] >= zoff && amax[2] < (zoff + dim[2])) {
      intersects = true;
      if (amin[2] < zoff) {
        amin[2] = zoff;
      }
    } else if (amin[2] < zoff && amax[2] >= (zoff + dim[2])) {
      intersects = true;
      amin[2] = zoff;
      amax[2] = zoff + dim[2];
    }
    // skip to the next atom if this one doesn't intersect our portion
    // of the volume
    if (!intersects) {
      continue;
    }
    // (end parallel specific)
#endif
    // unpack the color
    unsigned char aR, aG, aB;
    float aAlpha;
    // printf("Atom color %u\n", at->m_UserDefinedColor );
    aAlpha = ((float)(at->m_UserDefinedColor & 0xff)) / 255.0;
    aR = (unsigned char)((at->m_UserDefinedColor & 0xff000000) >> 24);
    aG = (unsigned char)((at->m_UserDefinedColor & 0xff0000) >> 16);
    aB = (unsigned char)((at->m_UserDefinedColor & 0xff00) >> 8);
    c[0] = orig[0] + c[0] * span[0];
    c[1] = orig[1] + c[1] * span[1];
    c[2] = orig[2] + c[2] * span[2];
    // begin blurring kernel
    for (l = amin[2]; l < amax[2]; l++) {
      for (k = amin[1]; k < amax[1]; k++) {
        for (j = amin[0]; j < amax[0]; j++) {
          int n = j + k * dim[0] + (l - zoff) * dim[0] * dim[1];
          float pnt[3], density;
          pnt[0] = orig[0] + j * span[0];
          pnt[1] = orig[1] + k * span[1];
          pnt[2] = orig[2] + l * span[2];
          // calculate the electron density / hydrophobicity
          density = evalDensity(/**at,*/ c, radius, pnt, densityType,
                                hydrophobe, charge);
          // if( l==(amax[2]+amin[2])/2 && k==(amax[1]+amin[1])/2 )
          //  printf("%f\n", density );
          if (withColor) {
            // accumulate
            // colormapped blurring: scale the density by alpha
            dens[n] += density * aAlpha;
          } else {
            // no colormap: no scaling
            dens[n] += density;
          }
          // calculate the color
          if (withColor && density > maxDens[n]) {
            // update the maximum density for this voxel
            maxDens[n] = density;
            // assign the color
            red[n] = aR;
            green[n] = aG;
            blue[n] = aB;
          }
        }
      }
    }
  }
}

double getElecDecay(double r) {
  r = fabs(r);
  double decay = 0;
  if (r > 2) {
    decay = 1.0 / r;
  } else {
    decay = 0.5;
  }
  double sigmoidalFunc = 0;
  if (r <= 6) {
    sigmoidalFunc = 4;
  } else if (r >= 8) {
    sigmoidalFunc = 80;
  } else {
    sigmoidalFunc = 38 * r - 224;
  }
  return decay * sigmoidalFunc;
}

float GOABlur::evalDensity(/*const PDBParser::Atom& atom,*/ double ctr[3],
                           float radius, float pnt[3],
                           PDBParser::GroupOfAtoms::FUNCTIONS densityType,
                           double hydrophobicity, double charge) {
  // cout << "b " << mBlobbyness << " " << radius  << endl;
  if (radius < 0.0000001) {
    return 0;
  }
  double expval = 0;
  if (densityType == PDBParser::GroupOfAtoms::ELECTROSTATIC_POTENTIAL) {
    /*double r = sqrt((ctr[0]-pnt[0])*(ctr[0]-pnt[0]) +
            (ctr[1]-pnt[1])*(ctr[1]-pnt[1]) +
            (ctr[2]-pnt[2])*(ctr[2]-pnt[2]));
    return charge*getElecDecay( r );*/
    double blobby = mBlobbyness;
    double r = (ctr[0] - pnt[0]) * (ctr[0] - pnt[0]) +
               (ctr[1] - pnt[1]) * (ctr[1] - pnt[1]) +
               (ctr[2] - pnt[2]) * (ctr[2] - pnt[2]);
    // double r = (ctr[0]-pnt[0])*(ctr[0]-pnt[0]);
    double r0 = radius * radius;
    expval = blobby * r / r0 - blobby;
    return (float)(charge * exp(expval));
  }
  bool useLinear = false;
  if (useLinear) {
    double r = sqrt((ctr[0] - pnt[0]) * (ctr[0] - pnt[0]) +
                    (ctr[1] - pnt[1]) * (ctr[1] - pnt[1]) +
                    (ctr[2] - pnt[2]) * (ctr[2] - pnt[2]));
    expval = mBlobbyness * (r - radius);
  } else {
    double blobby = mBlobbyness;
    //            double blobby = -3.5*log(4)/(radius);   // Wenqi (densityType
    //            = ELECTRON_DENSITY)
    double r = (ctr[0] - pnt[0]) * (ctr[0] - pnt[0]) +
               (ctr[1] - pnt[1]) * (ctr[1] - pnt[1]) +
               (ctr[2] - pnt[2]) * (ctr[2] - pnt[2]);
    // double r = (ctr[0]-pnt[0])*(ctr[0]-pnt[0]);
    double r0 = radius * radius;
    expval = blobby * r / r0 - blobby;
    // cout << "b " << mBlobbyness << " " << radius  << " " << expval << " " <<
    // r << endl;
  }
  if (densityType == PDBParser::GroupOfAtoms::ELECTRON_DENSITY) {
    return (float)(exp(expval));
  } else // hydrophobicity
  {
    return (float)(hydrophobicity * exp(expval));
  }
}
