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
#ifndef _MOLSURF_API_H_
#define _MOLSURF_API_H_

#include <Utility/utility.h>
#include <ASpline/ASpline.h>
#include <ASpline/ASplinePrimitives.h>
#include <Blurmaps/BlurMapsDataManager.h>
#include <Blurmaps/MolecularSurface.h>
#include <Blurmaps/SurfaceAtomExtractor.h>
#include <ColorSurfaceByFunction/colorSurfaceByCharge.h>
#include <ComputeNormals/geometry.h>
#include <ComputeNormals/GeometryParser.h>
#include <ComputeNormals/Surface.h>
#include <Decimation/Decimation.h>
#include <GeometryFileTypes/GeometryLoader.h>
#include <Geometry/Geometry.h>
#include <GOAFileTypes/GOALoader.h>
#include <HLevelSet/HLevelSet.h>
#include <LBIE_lib/LBIE_Mesher.h>
#include <libCG/ARM/ARM.h>

#ifdef WITH_CG
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/coarseGrain.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>
#include <libCG/CoarseGrain/PQRParser.h>
#include <libCG/CoarseGrain/surface.h>
#include <MolSurfAPI/pqr.h>
#endif // WITH_CG

#include <libCG/OptimizerLib/geometryParser.h>
#include <PDBParser/GroupOfAtoms.h>
#include <Pocket/Pocket.h>

#ifdef WITH_CGAL
#include <PocketTunnel/pocket_tunnel.h>
using namespace PocketTunnel;
#endif

#include <SignDistanceFunction_v2/DistanceTransform.h>
#include <SignDistanceFunction_v2/Geom3DParser.h>
#include <SimpleVolumeData/SimpleVolumeData.h>
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>
#include <VolumeFileTypes/VolumeLoader.h>

using namespace GEOMETRY;

extern float IsoValueofHLS;

bool area(int argc, char *argv[]);
double area2(char *rawnFile);
bool aSplineQuad(int argc, char *argv[]);
Geometry *aSplineRaw(Geometry *geometry, int resolution);
Geometry *aSplineRawQuad(Geometry *geometry, int numOfPts, int resolution,
                         string output1, string output2, string type);
bool aSplineVolume(int argc, char *argv[]);
bool bornRadius(int argc, char *argv[]);
bool coarseGrainCharge(int argc, char *argv[]);
bool coarseGrain(int argc, char *argv[]);
bool convertToRawnUsingHLS(int argc, char *argv[]);
bool decimate(int argc, char *argv[]);
bool decimatePercentage(int argc, char *argv[]);
bool reOrientate(int argc, char *argv[]);
bool normals(int argc, char *argv[]);
bool optimizeCGCharge(int argc, char *argv[]);
bool optimizeCGGeometry(int argc, char *argv[]);
bool pocket(int argc, char *argv[]);
bool pocketTunnel(int argc, char *argv[]);
bool populateSAS(int argc, char *argv[]);
bool populateSASUsingMesh(int argc, char *argv[]);
bool getInterfaceAtoms(int argc, char *argv[]);
bool getInterfaceSurface(int argc, char *argv[]);
bool colorByCharge(int argc, char *argv[]);
bool generateF2d(int argc, char *argv[]);
bool qualityImprove(int argc, char *argv[]);
bool rotatePDB(int argc, char *argv[]);
bool removeInteriorPockets(int argc, char *argv[]);
bool signedDistanceFunction(int argc, char *argv[]);
bool signedDistanceFunction2(int argc, char *argv[]);
bool surfaceAtoms(int argc, char *argv[]);
bool surfaceFromPDB(int argc, char *argv[]);
bool volume(int argc, char *argv[]);
bool volumeUsingHLS(int argc, char *argv[]);

bool getVolumeWithHydrophobicity(int argc, char *argv[]);
bool getVolumeWithElectron(int argc, char *argv[]);
bool getVolumeWithChargeDensity(int argc, char *argv[]);

bool surfaceUsingLBIE(int argc, char *argv[]);
Geometry *surfaceUsingHLS(PDBParser::GroupOfAtoms *molecule, int size,
                          int iterations, double weight);
Geometry *surfaceUsingAdaptiveGrid(PDBParser::GroupOfAtoms *molecule, int size,
                                   double width);
Geometry *surfaceUsingAdaptiveGrid(PDBParser::GroupOfAtoms *molecule, int size);
Geometry *surfaceUsingGaussianBlur(PDBParser::GroupOfAtoms *molecule,
                                   double isovalue, int size, double blobbiness,
                                   int radiusInt, int level,
                                   bool optimizeRadii);
Geometry *surfaceUsingGaussianBlur2(PDBParser::GroupOfAtoms *molecule,
                                    double isovalue, double resolution,
                                    double blobbiness, int radiusInt,
                                    int level);

PDBParser::GroupOfAtoms *loadGOA(const string &file);
void saveAndDeleteGeometry(const string &file, Geometry *geometry);

typedef struct {
  float p[3];
  float rad;
} weiPoint;

#endif
