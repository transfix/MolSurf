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
#ifndef CCV_BLURMAPSDATAMANAGER_H
#define CCV_BLURMAPSDATAMANAGER_H

#ifdef _WIN32
#define LITTLE_ENDIAN
#endif

#include <Utility/utility.h>
#include <Blurmaps/AtomPopulator.h>
#include <Blurmaps/SkinGenerator.h>
#include <Blurmaps/InterfaceAtomDetector.h>
#include <Blurmaps/InterfaceSurface.h>
#include <Blurmaps/F2dGenerator.h>
#include <Blurmaps/CurvaturesGrid.h>
#include <Blurmaps/GOABlur.h>
#include <Blurmaps/MultilevelBlur.h>
#include <Blurmaps/SignDistanceMolSurface.h>
#include <Blurmaps/SkinRegion2.h>
#include <Blurmaps/SkinRegion.h>
#include <Curvatures/SumOfGaussiansCurvature.h>
#include <GeometryFileTypes/GeometryLoader.h>
#include <Geometry/Geometry.h>
#include <GOAFileTypes/GOALoader.h>
#include <PDBParser/Atom.h>
#include <PDBParser/CollectionData.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/MutatedAtomInformation.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>
#include <VolumeFileTypes/VolumeLoader.h>

// arand comment...
// using PDBParser::Atom;
// using PDBParser::GroupOfAtoms;

namespace PDBParser {
class Atom;
};

class Geometry;
class SimpleVolumeData;
class GOABlur;

class BlurMapsDataManager {
public:
  BlurMapsDataManager();
  virtual ~BlurMapsDataManager();

  /*

  static bool flattenGOA(PDBParser::GroupOfAtoms* grp,
                                             vector<PDBParser::Atom*> &
  flatList, PDBParser::CollectionData* collectionData, int
  numberOfTransformations, double* rotations, double* translations,
                                             PDBParser::GroupOfAtoms::RADIUS_TYPE
  radiusType, unsigned int level, bool colorBySubunits);
  */

  static void getBoundingBox(vector<PDBParser::Atom *> &grp, float min[3],
                             float max[3],
                             PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                             int probeRadius, double blobbiness);

  static SimpleVolumeData *
  getVolume(PDBParser::GroupOfAtoms *molecule, const char *volFileName,
            int dim1, int dim2, int dim3,
            PDBParser::GroupOfAtoms::FUNCTIONS densityType, bool writeRawV,
            double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
            const char *cmapFile, int gap,
            PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
            CCVOpenGLMath::Matrix *transformation = 0);

  static SimpleVolumeData *
  getVolume(PDBParser::GroupOfAtoms *molecule, const char *volFileName,
            int dim1, int dim2, int dim3, float bbmin[3], float bbmax[3],
            PDBParser::GroupOfAtoms::FUNCTIONS densityType, bool writeRawV,
            double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
            const char *cmapFile, int gap,
            PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
            CCVOpenGLMath::Matrix *transformation = 0);

  static SimpleVolumeData *getVolumeFixedResolution(
      PDBParser::GroupOfAtoms *molecule, const char *volFileName,
      double resolution, PDBParser::GroupOfAtoms::FUNCTIONS densityType,
      bool writeRawV, double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
      const char *cmapFile, int gap,
      PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
      CCVOpenGLMath::Matrix *transformation = 0);

  static SimpleVolumeData *getMultiLevelBlurredVolume(
      PDBParser::GroupOfAtoms *molecule, const char *volFileName, int dim1,
      int dim2, int dim3, PDBParser::GroupOfAtoms::FUNCTIONS densityType,
      bool writeRawV, double blob, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel,
      const char *cmapFile, int gap,
      PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
      CCVOpenGLMath::Matrix *transformation = 0);

  static bool getCurvatures(const char *proteinFileName, int dim1, int dim2,
                            int dim3, double blob, Geometry *geometry,
                            const char *outputMeanRawSurfaceFileName,
                            const char *outputGaussianRawSurfaceFileName,
                            const char *curvatureFileName,
                            int numberOfGridDivisions, double maxFunctionError,
                            PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                            unsigned int level);
  static bool getCurvaturesFromIsocontourFile(
      const char *proteinFileName, int dim1, int dim2, int dim3, double blob,
      const char *inputRawSurfaceFileName,
      const char *outputMeanRawSurfaceFileName,
      const char *outputGaussianRawSurfaceFileName,
      const char *curvatureFileName, int numberOfGridDivisions,
      double maxFunctionError, PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
      unsigned int level);

  static SimpleVolumeData *
  getSkinRegionVolume(const char *inputFileName, const char *outputFileName,
                      int dim1, int dim2, int dim3, double probeRadius,
                      PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                      int depth);

  static bool
  getInterfaceAtoms(string pqrFileName1, string pqrFileName2,
                    string outputFileName, double interfaceWidth,
                    string atomTypes,
                    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);

  static bool getInterfaceSurface(string surfaceFileName1,
                                  string surfaceFileName2,
                                  double interfaceWidth, bool colored,
                                  double &area1, double &area2);

  static bool generateF2d(string pqrFileName, string xyzFileName,
                          string outputFileName, bool receptor,
                          PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType,
                          double skinAtomRad);

  static bool populateSAS(const char *inputFileName, const char *outputFileName,
                          int dim1, int dim2, int dim3, double probeRadius,
                          double floatingBandWidth,
                          PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);
  static bool
  populateSASUsingMesh(string pqrFileName, string surfaceFileName,
                       string outputFileName, double probeRadius,
                       double floatingBandWidth, double clashDistance,
                       PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType);

  static Geometry *
  getSurfaceFromPDB(PDBParser::GroupOfAtoms *molecule, int dim1, int dim2,
                    int dim3, double isovalue, double blobbiness,
                    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType, int level);

protected:
  static bool outputCurvFiles(double *HandK, Geometry *geometry,
                              const char *outputMeanRawSurfaceFileName,
                              const char *outputGaussianRawSurfaceFileName,
                              const char *curvatureFileName);
  static double getDistanceOfVoxel(int i, int j, int k, int width, int height,
                                   int depth);
};

#endif
