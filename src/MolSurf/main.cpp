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
#include <MolSurfAPI/molsurf.h>
// #include <DPG/PG.h>
#include <iostream>
#include <map>

// TODO: cleanup the API's below and use namespaces
//       to avoid future conflicts
#include <OffsetSurface/MakeOffsetSurface.h>
#include <RepairSurfaceMesh/RepairMesh.h>
#include <PartitionMesh/PartitionMesh.h>
#include <CardinalBSpline/CardinalBSplineAPI.h>
#include <LeastRMSD/LeastRMSDInterface.h>
#include <SimpleCuration/SimpleCurationAPI.h>
#include <FastMolSurf/fastMolSurfAPI.h>
#include <Hausdorff/HausdorffAPI.h>

void printUsage();

// aspline raw
void usageASplineRaw() {
  cout << endl
       << "Usage: MolSurf -aSpline -raw <raw input file> <raw output file> "
          "<resolution>"
       << endl
       << "    Generate A-spline model (ASMS) from a triangular mesh." << endl;
}

bool aSplineRaw(int argc, char *argv[]) {
  // First read in the arguments
  if (argc != 6) // If the wrong number of arguments, print usage
  {
    usageASplineRaw();
    exit(0);
  }
  string input = string(argv[3]);
  string output = string(argv[4]);
  int resolution = atoi(argv[5]);
  Geometry *geometry = GeometryLoader().loadFile(input);
  Geometry *newGeo = aSplineRaw(geometry, resolution);

  saveAndDeleteGeometry(output, newGeo);
  return true;
}

// Aspline raw and quadrature
void usageASplineRawQuad() {
  cout << endl
       << "Usage: MolSurf -aSpline -raw+quad <raw input file> <raw output file>"
       << endl
       << "               <raw output w/ quad points> <quadrature file>" << endl
       << "               <resolution> <quadrature rule> <quadrature order>"
       << endl
       << "   Generate A-spline model (ASMS) from a triangular mesh and "
          "generate quadrature points."
       << endl;
}

bool aSplineRawQuad(int argc, char *argv[]) {
  // First read in the arguments
  if (argc != 10) // If the wrong number of arguments, print usage
  {
    usageASplineRawQuad();
    return false;
  }
  string input = string(argv[3]);
  string output0 = string(argv[4]);
  string output1 = string(argv[5]);
  string output2 = string(argv[6]);
  int resolution = atoi(argv[7]);
  string type = string(argv[8]);
  int numOfPts = atoi(argv[9]);
  Geometry *geometry = GeometryLoader().loadFile(input);

  Geometry *newGeo =
      aSplineRawQuad(geometry, numOfPts, resolution, output1, output2, type);

  saveAndDeleteGeometry(output0, geometry);
  return true;
}

// Get surface using gaussian blur
void usageSurfaceUsingGaussianBlur() {
  cout << endl
       << "Usage: MolSurf -getSurfaceUsingGaussianBlur <pdb/pqr input file> "
          "<raw output file> "
       << endl
       << "         <isovalue>   <size> <blobbiness> <radius type> <level> "
          "<optimzeRadii> "
       << endl
       << "    Generate molecular surface using the Gaussian blur algorithm."
       << endl
       << endl;
}

bool surfaceUsingGaussianBlur(int argc, char *argv[]) {
  // First read in the arguments
  if (argc != 10) // If the wrong number of arguments, print usage
  {
    usageSurfaceUsingGaussianBlur();
    exit(0);
  }
  string pdbFile = string(argv[2]);
  string surfaceFile = string(argv[3]);
  double isovalue = atof(argv[4]);
  int size = atoi(argv[5]);
  double blobbiness = (double)atof(argv[6]);
  int radiusInt = atoi(argv[7]);
  int level = atoi(argv[8]);
  bool optimizeRadii = (strcmp(argv[9], "true") == 0);

  // Read a file, call the function, write a file.
  PDBParser::GroupOfAtoms *molecule = loadGOA(pdbFile);
  Geometry *geo = surfaceUsingGaussianBlur(molecule, isovalue, size, blobbiness,
                                           radiusInt, level, optimizeRadii);
  saveAndDeleteGeometry(surfaceFile, geo);

  return true;
}

// Get surface using gaussian blur
void usageSurfaceUsingGaussianBlur2() {
  cout
      << endl
      << "Usage: MolSurf -getSurfaceUsingGaussianBlur2 <pdb/pqr input file> "
         "<raw output file> "
      << endl
      << "         <isovalue>   <resolution> <blobbiness> <radius type> <level>"
      << endl
      << "    Generate molecular surface using the Gaussian blur algorithm."
      << endl
      << endl;
}

bool surfaceUsingGaussianBlur2(int argc, char *argv[]) {
  // First read in the arguments
  if (argc != 9) // If the wrong number of arguments, print usage
  {
    usageSurfaceUsingGaussianBlur2();
    exit(0);
  }
  string pdbFile = string(argv[2]);
  string surfaceFile = string(argv[3]);
  double isovalue = atof(argv[4]);
  double resolution = atof(argv[5]);
  double blobbiness = (double)atof(argv[6]);
  int radiusInt = atoi(argv[7]);
  int level = atoi(argv[8]);

  // Read a file, call the function, write a file.
  PDBParser::GroupOfAtoms *molecule = loadGOA(pdbFile);
  Geometry *geo = surfaceUsingGaussianBlur2(molecule, isovalue, resolution,
                                            blobbiness, radiusInt, level);
  saveAndDeleteGeometry(surfaceFile, geo);

  return true;
}

// Get surface using adaptive grid
void usageSurfaceUsingFastSurf() {
  cout << "Usage: MolSurf -surfaceUsingFastSurf <pdb/pqr input file> <dim0> "
          "<dim1> <dim2> <bool SAS> <output.raw>"
       << endl
       << "    Generate molecular surface using a fast algorithm." << endl;
}

bool surfaceUsingFastSurf(int argc, char *argv[]) {
  if (argc < 8) {
    usageSurfaceUsingFastSurf();
    return false;
  }

  constructFastMolSurface(argc, argv);
  return true;
}

// Get surface using adaptive grid
void usageSurfaceUsingAdaptiveGrid() {

  cout << endl
       << "Usage: MolSurf -surfaceUsingAdaptiveGrid <pdb/pqr input file> <raw "
          "output file> <size> [ <enlarged width> ]"
       << endl
       << "    Generate molecular surface using the adaptive grid algorithm."
       << endl;
}

bool surfaceUsingAdaptiveGrid(int argc, char *argv[]) {
  // First read in the arguments
  if (argc != 5 && argc != 6) // If the wrong number of arguments, print usage
  {
    usageSurfaceUsingAdaptiveGrid();
    return false;
  }
  string pdbFile = string(argv[2]);
  string surfaceFile = string(argv[3]);
  int size = atoi(argv[4]);

  double width = 0.0;
  if (argc == 6) {
    width = atof(argv[5]);
  }

  PDBParser::GroupOfAtoms *molecule = loadGOA(pdbFile);

  Geometry *geometry = surfaceUsingAdaptiveGrid(molecule, size, width);
  delete molecule;
  saveAndDeleteGeometry(surfaceFile, geometry);
  return true;
}

// Get surface using HLS
void usageSurfaceUsingHLS() {
  cout << endl
       << "Usage: MolSurf -surfaceUsingHLS <pdb/pqr input file> <raw output "
          "file> "
       << endl
       << "          <size> [nIterations [weight]]" << endl
       << "    Generate molecular surface using the higher order level-set "
          "algorithm."
       << endl
       << "    nIterations and weight are optional arguments." << endl;
}

bool surfaceUsingHLS(int argc, char *argv[]) {
  if (argc != 5 && argc != 6 && argc != 7) {
    usageSurfaceUsingHLS();
    return false;
  }
  string inputFile = string(argv[2]);
  string outputFile = string(argv[3]);
  int size = atoi(argv[4]);
  int iterations = 3;
  if (argc == 6 || argc == 7) {
    iterations = atoi(argv[5]);
  }
  double weight = 0.0001;
  if (argc == 7) {
    weight = atof(argv[6]);
  }
  PDBParser::GroupOfAtoms *molecule = loadGOA(inputFile);
  Geometry *geometry = surfaceUsingHLS(molecule, size, iterations, weight);
  delete molecule;
  saveAndDeleteGeometry(outputFile, geometry);
  return true;
}

void usageSurfaceUsingCardinalBSpline() {
  cout << endl
       << "Usage: MolSurf -surfaceUsingCardinalBSpline <pdb/pqr input file> "
          "<int size>"
       << endl
       << "		<int order(0-4)> [float isovalue (0-0.5 (default "
          "0.05)) ] <raw output file>"
       << endl
       << " 		Generate molecular surface using cardinal B spline "
          "algorithm "
       << endl;
}

bool surfaceUsingCardinalBSpline(int argc, char *argv[]) {
  if (argc != 6 && argc != 7) {
    usageSurfaceUsingCardinalBSpline();
    return false;
  }

  string outputFile;
  float isovalue = 0.05;
  if (argc == 7) {
    isovalue = atof(argv[5]);
    outputFile = string(argv[6]);
  } else
    outputFile = string(argv[5]);

  string inputFile = string(argv[2]);
  int size = atoi(argv[3]);
  int order = atoi(argv[4]);

  PDBParser::GroupOfAtoms *molecule = loadGOA(inputFile);
  Geometry *geometry =
      surfaceUsingCardinalBSpline(molecule, size, order, isovalue);
  delete molecule;
  saveAndDeleteGeometry(outputFile, geometry);
  return true;
}

void printASplineUsage() {
  cout << endl
       << "Usage: MolSurf -aSpline -[raw|quad|raw+quad|volume] [args]" << endl
       << endl;
}

bool aSpline(int argc, char *argv[]) {

  if (argc < 3) {
    printASplineUsage();
    return false;
  }

  if (strcmp(argv[2], "-raw+quad") == 0) {
    return aSplineRawQuad(argc, argv);
  } else if (strcmp(argv[2], "-raw") == 0) {
    return aSplineRaw(argc, argv);
  } else if (strcmp(argv[2], "-quad") == 0) {
    return aSplineQuad(argc, argv);
  } else if (strcmp(argv[2], "-volume") == 0) {
    return aSplineVolume(argc, argv);
  }

  cout << endl << "Unknown command: " << argv[2] << endl;
  printASplineUsage();
  return false;
}

void printMergeGeometryUsage() {
  cout << endl
       << "Usage: MolSurf -mergeGeometry  <input file1> <input file2 ...> "
          "<outputfle>"
       << endl
       << endl;
}

bool mergeGeometry(int argc, char *argv[]) {
  if (argc < 5) {
    printMergeGeometryUsage();
    return false;
  }

  Geometry *all = GeometryLoader().loadFile(string(argv[2]));

  string inputFile;
  for (unsigned int i = 3; i < (unsigned int)(argc - 1); i++) {
    inputFile = string(argv[i]);
    Geometry *geom = GeometryLoader().loadFile(inputFile);
    if (geom != NULL) {
      all = all->merge(geom);
    } else {
      cout << "wrong" << endl;
      return false;
    }
    delete geom;
  }

  string outFile = string(argv[argc - 1]);

  saveAndDeleteGeometry(outFile, all);
  return true;
}

// using namespace MOLECULE;
#ifdef WITH_CG
void usageMergeProtein() {
  cout << endl
       << "Usage: ./MolSurf -mergeProtein <radius type: vdw or pqr> "
          "<inputFile1> <inputFile2> <outputFile> ";
  cout << endl;
}

bool mergeProtein(int argc, char *argv[]) {
  char inputFileA[256];
  char inputFileB[256];
  char outputFile[256];

  if (argc != 6) {
    usageMergeProtein();
    return false;
  }

  string radiusString = string(argv[2]);

  strcpy(inputFileA, argv[3]);
  strcpy(inputFileB, argv[4]);
  strcpy(outputFile, argv[5]);

  bool withH = 1;
  bool withChain = 1;
  MOLECULE::RADIUS_TYPE radiusType;

  if (radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0)
    radiusType = MOLECULE::VDW_RADIUS;
  else if (radiusString.compare("pqr") == 0 || radiusString.compare("PQR") == 0)
    radiusType = MOLECULE::PQR_RADIUS;
  else
    radiusType = MOLECULE::VDW_RADIUS;

  FILE *fpA = fopen(inputFileA, "r");
  FILE *fpB = fopen(inputFileB, "r");
  FILE *fpout = fopen(outputFile, "w");

  MOLECULE::GroupOfAtoms *inputA = new MOLECULE::GroupOfAtoms();
  MOLECULE::GroupOfAtoms *inputB = new MOLECULE::GroupOfAtoms();

  if (radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0) {
    inputA = ParsePDB(fpA, radiusType);
    inputB = ParsePDB(fpB, radiusType);
  } else if (radiusString.compare("pqr") == 0 ||
             radiusString.compare("PQR") == 0) {
    inputA = ParsePQR(fpA, radiusType);
    inputB = ParsePQR(fpB, radiusType);
  } else {
    inputA = ParsePDB(fpA, radiusType);
    inputB = ParsePDB(fpB, radiusType);
  }

  fclose(fpA);
  fclose(fpB);

  append(inputA, inputB);

  if (radiusString.compare("vdw") == 0 || radiusString.compare("VDW") == 0)
    savePDB(inputA, fpout, withH, withChain);
  else if (radiusString.compare("pqr") == 0 || radiusString.compare("PQR") == 0)
    savePQR(inputA, fpout);
  else
    savePDB(inputA, fpout, withH, withChain);

  fclose(fpout);
  return true;
}
void usagePQR2PDB() {
  cout << endl
       << "Usage: ./MolSurf -PQR2PDB -withH -withChain <inputFile.pqr>  "
          "<outputFile.pdb> ";
  cout << endl
       << " e.g.,./MolSurf -PQR2PDB 1 1 <inputFile.pqr>  <outputFile.pdb> with "
          "hydrogen, with chain";
  cout << endl
       << " e.g.,./MolSurf -PQR2PDB 0  0 <inputFile.pqr>  <outputFile.pdb> no "
          "hydrogen, no chain";
  cout << endl;
}

bool pqr2pdb(int argc, char *argv[]) {
  if (argc != 6) {
    usagePQR2PDB();
    return false;
  }

  char inputFile[256];
  char outputFile[256];

  strcpy(inputFile, argv[4]);
  strcpy(outputFile, argv[5]);

  bool withH = atoi(argv[2]);
  bool withChain = atoi(argv[3]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fin = fopen(inputFile, "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();

  input = ParsePQR(fin, radiusType);

  fclose(fin);

  FILE *fout = fopen(outputFile, "w");
  savePDB(input, fout, withH, withChain);
  fclose(fout);
  return true;
}

void usagePQR2XYZR() {
  cout << endl
       << "Usage: ./MolSurf -PQR2XYZR <inputFile.pqr>  <outputFile.xyzr> ";
  cout << endl;
}

bool pqr2xyzr(int argc, char *argv[]) {
  if (argc != 4) {
    usagePQR2XYZR();
    return false;
  }

  char inputFile[256];
  char outputFile[256];

  strcpy(inputFile, argv[2]);
  strcpy(outputFile, argv[3]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fin = fopen(inputFile, "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();

  input = ParsePQR(fin, radiusType);

  fclose(fin);

  fstream fout;
  fout.open(outputFile, ios::out);
  for (int i = 0; i < input->numOfAtoms; i++)
    fout << input->m_Atoms[i]->position[0] << " "
         << input->m_Atoms[i]->position[1] << " "
         << input->m_Atoms[i]->position[2] << " " << input->m_Atoms[i]->radius
         << endl;
  fout.close();
}

void usageSplitPDB() {
  cout << endl
       << "Usage: ./MolSurf -splitPDB2Chains <inputFile.pdb>  <outputFile "
          "(without ext)> ";
  cout << endl;
}

bool splitPDB2Chains(int argc, char *argv[]) {
  char inputFile[256];
  char outputFile[256];

  if (argc != 4) {
    usageSplitPDB();
    return false;
  }

  strcpy(inputFile, argv[2]);
  strcpy(outputFile, argv[3]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::VDW_RADIUS;
  //      if(radiusString.compare("vdw") ==0 || radiusString.compare("VDW") ==
  //      0)
  //              radiusType = MOLECULE::VDW_RADIUS;
  //      else if(radiusString.compare("pqr") == 0 ||
  //      radiusString.compare("PQR") == 0)
  //              radiusType = MOLECULE::PQR_RADIUS;
  //      else
  //             radiusType = MOLECULE::VDW_RADIUS;

  FILE *fpA = fopen(inputFile, "r");

  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();

  input = ParsePDB(fpA, radiusType);

  fclose(fpA);

  char tmpChainID;
  int tmpResSeq;
  if (input->m_Atoms.size() != 0) {
    tmpChainID = input->m_Atoms[0]->chainID;
    tmpResSeq = input->m_Atoms[0]->resSeq;
  }

  strcat(outputFile, "_");
  char tmpFileName[256];
  strcpy(tmpFileName, outputFile);
  int length = strlen(tmpFileName);
  tmpFileName[length] = tmpChainID;
  tmpFileName[length + 1] = '\0';
  strcat(tmpFileName, ".pdb");

  FILE *fpout = fopen(tmpFileName, "w");

  int j = 1;
  for (int i = 0; i < input->m_Atoms.size(); i++) {
    MOLECULE::Atom *at = input->m_Atoms[i];
    if (at->chainID != tmpChainID) {
      fprintf(fpout, "TER %7d      %3s %c%4d\n", j,
              input->m_Atoms[i - 1]->resName, input->m_Atoms[i - 1]->chainID,
              input->m_Atoms[i - 1]->resSeq - tmpResSeq + 1);
      fclose(fpout);
      j = 1;
      tmpChainID = at->chainID;
      tmpResSeq = at->resSeq;
      strcpy(tmpFileName, outputFile);
      tmpFileName[length] = tmpChainID;
      tmpFileName[length + 1] = '\0';
      strcat(tmpFileName, ".pdb");
      FILE *fpout = fopen(tmpFileName, "w");
    }

    fprintf(
        fpout,
        "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n", j,
        at->atomName, at->resName, at->chainID, at->resSeq - tmpResSeq + 1,
        at->position[0], at->position[1], at->position[2], at->occupancy,
        at->tempFactor, at->atomID);
    j++;
  }

  fprintf(fpout, "TER %7d      %3s %c%4d\n", j, input->m_Atoms.back()->resName,
          input->m_Atoms.back()->chainID,
          input->m_Atoms.back()->resSeq - tmpResSeq + 1);

  fclose(fpout);
  return true;
}

using namespace DPG;

void usageColorMeshByAtomType() {
  cout << endl
       << "Usage: ./MolSurf -colorMeshByAtomType <inputFile.pdb> "
          "<inputMesh.raw>  <outputMesh.rawc> ";
  cout << endl;
}

bool ColorMeshByAtomType(int argc, char *argv[]) {

  if (argc != 5) {
    usageColorMeshByAtomType();
    return false;
  }

  char inputPDB[256];
  char inputMesh[256];
  char outputMesh[256];

  strcpy(inputPDB, argv[2]);
  strcpy(inputMesh, argv[3]);
  strcpy(outputMesh, argv[4]);

  FILE *fpA = fopen(inputPDB, "r");

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::VDW_RADIUS;
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();

  input = ParsePDB(fpA, radiusType);
  fclose(fpA);

  PG *pg = new PG(10.0, 1000.0, 3.0);

  vector<DPG::WeightedPoint *> pts;

  //      std::map<DPG::WeightedPoint, int> pdbmap;

  //      int pt_index = 0;
  for (int i = 0; i < input->m_Atoms.size(); i++) {
    DPG::WeightedPoint *pt = new DPG::WeightedPoint(
        input->m_Atoms[i]->position[0], input->m_Atoms[i]->position[1],
        input->m_Atoms[i]->position[2], (double)i);
    //              const DPG::WeighPoint pp = *pt;

    //              (pt->getX(), pt->getY(), pt->getZ());

    //              std::cout<<"pp.x: "<< pp.getX() <<" " << pp.getY() << " "
    //              <<pp.getZ() << endl;
    pts.push_back(pt);
    pg->addPoint(pts[i]);
    //     if(pdbmap.find(pp) == pdbmap.end())
    //     {
    //                      pdbmap[pp] = pt_index;
    //                      pt_index++;
    //              }
  }

  Geometry *geometry = GeometryLoader().loadFile(inputMesh);

  float color[][3] = {
      {0.3, 0.3, 0.3}, // C
      {1, 0, 0},       // O
      {0, 0, 1},       // N
      {0, 1, 0},       // others
  };

  vector<DPG::Point *> results;

  char atomType;

  int ColorID;

  geometry->AllocateTriVertColors();

  for (int i = 0; i < geometry->m_NumTriVerts; i++) {

    DPG::Point *p = new DPG::Point(geometry->m_TriVerts[3 * i + 0],
                                   geometry->m_TriVerts[3 * i + 1],
                                   geometry->m_TriVerts[3 * i + 2]);
    geometry->m_TriVertColorsTransparent[4 * i + 0] = 0.0;
    geometry->m_TriVertColorsTransparent[4 * i + 1] = 0.0;
    geometry->m_TriVertColorsTransparent[4 * i + 2] = 0.0;
    geometry->m_TriVertColorsTransparent[4 * i + 3] = 0.0;

    results = pg->range(p, 2.0);
    vector<int> indices;

    if (results.size() == 0)
      results = pg->range(p, 3.0);

    for (int j = 0; j < results.size(); j++) {
      int index = (int)(((WeightedPoint *)results[j])->getWeight());
      indices.push_back(index);
    }
    results.clear();

    for (int j = 0; j < indices.size(); j++) {
      int ascii;
      for (int t = 0; t < 6; t++) {
        if (indices[j] > input->m_Atoms.size())
          printf(" indices %d\n", indices[j]);
        ascii = (int)(input->m_Atoms[indices[j]]->atomName[t]);
        if (ascii <= 90 && ascii >= 65)
          break;
      }
      if (ascii == 67)
        ColorID = 0;
      else if (ascii == 79)
        ColorID = 1;
      else if (ascii == 78)
        ColorID = 2;
      else
        ColorID = 3;

      geometry->m_TriVertColorsTransparent[4 * i + 0] +=
          color[ColorID][0] / (float)indices.size();
      geometry->m_TriVertColorsTransparent[4 * i + 1] +=
          color[ColorID][1] / (float)indices.size();
      geometry->m_TriVertColorsTransparent[4 * i + 2] +=
          color[ColorID][2] / (float)indices.size();
    }
    indices.clear();
  }

  if (!GeometryLoader().saveFile(outputMesh, geometry)) {
    delete geometry;
    return false;
  }
  return true;
}

void usageAtomCentroid() {
  cout << endl << "Usage: ./MolSurf -centroid <inputPQR>" << endl;
}

bool atomCentroid(int argc, char *argv[]) {
  if (argc != 3) {
    usageAtomCentroid();
  }

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;
  FILE *fp = fopen(argv[2], "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fp, radiusType);
  fclose(fp);

  double centroid[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < input->m_Atoms.size(); i++)
    for (int j = 0; j < 3; j++)
      centroid[j] += input->m_Atoms[i]->position[j];

  for (int j = 0; j < 3; j++)
    centroid[j] /= (double)input->m_Atoms.size();

  cout << "Centroid: " << centroid[0] << " " << centroid[1] << " "
       << centroid[2] << endl;
}

void usageAtomsBoundingbox() {
  cout << endl
       << "Usage: ./MolSurf -bbox <pqrlist?(0 or 1)> <inputPQR/PQRlistfile>"
       << endl;
}

bool atomsBoundingbox(int argc, char *argv[]) {
  if (argc != 4) {
    usageAtomsBoundingbox();
  }

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;
  bool hasListfile = (atoi(argv[2]) == 1);
  if (hasListfile)
    fprintf(stderr, "has list file: %s\n", argv[3]);
  else
    fprintf(stderr, "single file: %s\n", argv[3]);

  double min[3] = {100000000.0, 100000000.0, 100000000.0};
  double max[3] = {-100000000.0, -100000000.0, -100000000.0};
  if (hasListfile) {
    FILE *lfp = fopen(argv[3], "r");
    char line[1025];
    while (NULL != fgets(line, 1024, lfp)) {
      char pqrfile[1025];
      sscanf(line, "%s", pqrfile);
      FILE *fp = fopen(pqrfile, "r");
      if (!fp) {
        fprintf(stderr, "%s file open fail\n", pqrfile);
        continue;
      }
      MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
      input = ParsePQR(fp, radiusType);
      fclose(fp);

      for (int i = 0; i < input->m_Atoms.size(); i++)
        for (int j = 0; j < 3; j++) {
          if (min[j] > input->m_Atoms[i]->position[j])
            min[j] = input->m_Atoms[i]->position[j];
          if (max[j] < input->m_Atoms[i]->position[j])
            max[j] = input->m_Atoms[i]->position[j];
        }
      delete input;
    }
    fclose(lfp);
  } else {
    FILE *fp = fopen(argv[3], "r");
    MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
    input = ParsePQR(fp, radiusType);
    fclose(fp);

    for (int i = 0; i < input->m_Atoms.size(); i++)
      for (int j = 0; j < 3; j++) {
        if (min[j] > input->m_Atoms[i]->position[j])
          min[j] = input->m_Atoms[i]->position[j];
        if (max[j] < input->m_Atoms[i]->position[j])
          max[j] = input->m_Atoms[i]->position[j];
      }
    delete input;
  }
  cout << "Bbox min: " << min[0] << " " << min[1] << " " << min[2] << endl;
  cout << "Bbox max: " << max[0] << " " << max[1] << " " << max[2] << endl;
}

void usageAtom2SurfaceDistance() {
  cout << endl
       << "Usage: ./Molsurf -atom2surfaceDistance <inputPQR> <inputMesh> "
          "<outputCDR>";
  cout << endl
       << "Warning. outputCDR is (P)ositio replaced by (C)olor, (Q)charge "
          "replace by (D)istance	to the surface and R invarant";
  cout << endl;
}

bool Atom2SurfaceDistance(int argc, char *argv[]) {
  if (argc != 5) {
    usageAtom2SurfaceDistance();
    return false;
  }

  char inputPQR[256];
  char inputMesh[256];
  char outputCDR[256];

  strcpy(inputPQR, argv[2]);
  strcpy(inputMesh, argv[3]);
  strcpy(outputCDR, argv[4]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fout = fopen(outputCDR, "w");
  FILE *fpA = fopen(inputPQR, "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fpA, radiusType);
  fclose(fpA);

  Geometry *geometry = GeometryLoader().loadFile(inputMesh);
  //      double color[5][3]={
  //              {255, 0, 0},  //C
  //              {0, 255, 0},  // O
  //              {0, 0, 255},  // N
  //              {255, 255, 0},  // H
  //              {0, 255, 255},  // others
  //      };

  PG *pg = new PG(10.0, 1000.0, 3.0);

  //      vector <DPG::Point* > meshPoints;
  for (int i = 0; i < geometry->m_NumTriVerts; i++) {
    DPG::Point *p = new DPG::Point(geometry->m_TriVerts[3 * i + 0],
                                   geometry->m_TriVerts[3 * i + 1],
                                   geometry->m_TriVerts[3 * i + 2]);
    //              meshPoints.push_back(p);
    pg->addPoint(p);
  }

  vector<DPG::Point *> results;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    DPG::Point *pt = new DPG::Point(input->m_Atoms[i]->position[0],
                                    input->m_Atoms[i]->position[1],
                                    input->m_Atoms[i]->position[2]);
    float threshold0 = 0.5;
    float threshold = threshold0;
    int j = 1;
    do {
      results = pg->range(pt, threshold);
      j++;
      threshold += threshold0 * j;
    } while (results.size() == 0);

    float dist = pt->distance(results[0]);
    float temp;
    for (int k = 1; k < results.size(); k++) {
      temp = pt->distance(results[k]);
      if (temp < dist)
        dist = temp;
    }
    input->m_Atoms[i]->setCharge(dist);
    int ascii;
    float r, g, b;
    for (int t = 0; t < 6; t++) {
      ascii = (int)(input->m_Atoms[i]->atomName[t]);
      if (ascii <= 90 && ascii >= 65)
        break;
    }
    switch (ascii) {
    case 67:
      input->m_Atoms[i]->setOccupancy(0);
      break;
    case 79:
      input->m_Atoms[i]->setOccupancy(1);
      break;
    case 78:
      input->m_Atoms[i]->setOccupancy(2);
      break;
    case 72:
      input->m_Atoms[i]->setOccupancy(3);
      break;
    default:
      input->m_Atoms[i]->setOccupancy(4);
      break;
    }
  }
  saveCDR(input, fout);
  fclose(fout);

  return true;
}

void readUntilNewline(ifstream &fin) {

  char c = 'a';
  while (c != '\n') {
    fin.get(c);
  }
}

void usageTinkerPDBXYZ2PQR() {
  cout << endl
       << "Usage: ./Molsurf -tinkerPDBXYZ2PQR <inputPDB> <inputXYZ> "
          "<inputBornRadiiFile> <outputPQR>";
  cout << endl << "inputPDB is the output of Tinker xyzpdb";
  cout << endl;
}

bool TinkerPDBXYZ2PQR(int argc, char *argv[]) {
  if (argc != 6) {
    usageTinkerPDBXYZ2PQR();
    return false;
  }

  char inputPDB[256];
  char inputXYZ[256];
  char outputPQR[256];

  strcpy(inputPDB, argv[2]);
  strcpy(inputXYZ, argv[3]);
  string inputBorn = argv[4];
  strcpy(outputPQR, argv[5]);

  FILE *fout = fopen(outputPQR, "w");
  FILE *fpA = fopen(inputPDB, "r");

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::VDW_RADIUS;

  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePDB(fpA, radiusType);
  fclose(fpA);

  ifstream fpXYZ(inputXYZ);
  ifstream fpBorn(inputBorn.c_str());

  string tmp;

  int atmNum;
  fpXYZ >> atmNum;

  while (getline(fpBorn, tmp)) {
    if (strncmp(tmp.c_str(), " Born", 5) == 0)
      break;
  }
  readUntilNewline(fpBorn);

  //      fscanf(fpXYZ, "%d\n", atmNum);
  if (atmNum != input->m_Atoms.size()) {
    cout << "The number of atoms in XYZ file is different from PDB file, exit!"
         << endl;
    return false;
  }
  float weightPoints[atmNum][4];
  int atnInd;
  char atnName[10];
  vector<int> Index;
  vector<int>::iterator iit;

  vector<string> text;
  for (int i = 0; i < atmNum; i++) {
    if (i % 5 == 0) {
      text.clear();
      //                      cout<<"i= " << i << endl;
      getline(fpBorn, tmp);
      //                      cout<<tmp<<endl;

      istringstream iss(tmp);
      string s;
      while (iss >> s)
        text.push_back(s);
      //    cout<<text[5] << endl;
      if (!(atoi(text[0].c_str()) == i + 1)) {
        cout << "Error reading the Born file, exit!" << endl;
        return false;
      }
    }

    fpXYZ >> atnInd >> atnName >> weightPoints[i][0] >> weightPoints[i][1] >>
        weightPoints[i][2];
    weightPoints[i][3] = atof(text[2 * (i % 5) + 1].c_str());
    readUntilNewline(fpXYZ);
    //              cout<<  weightPoints[i][0] << " " <<  weightPoints[i][1] <<
    //              " " <<  weightPoints[i][2] << " " <<  weightPoints[i][3] <<
    //              endl;
    Index.push_back(i);
  }

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    input->m_Atoms[i]->setCharge(1.0);
    for (iit = Index.begin(); iit < Index.end(); iit++) {
      //    cout<< *iit <<" " << weightPoints[*iit][0] << weightPoints[*iit][1]
      //    << weightPoints[*iit][2]<<" " <<weightPoints[*iit][3]  << endl;
      //    cout<<input->m_Atoms[i]->position[0]<<" "
      //    <<input->m_Atoms[i]->position[1] <<" "
      //    <<input->m_Atoms[i]->position[2]  << endl;
      if (fabs(weightPoints[*iit][0] - input->m_Atoms[i]->position[0]) <
              0.001 &&
          fabs(weightPoints[*iit][1] - input->m_Atoms[i]->position[1]) <
              0.001 &&
          fabs(weightPoints[*iit][2] - input->m_Atoms[i]->position[2]) <
              0.001) {
        input->m_Atoms[i]->setDefinedRadius(weightPoints[*iit][3]);
        Index.erase(iit);

        //            cout<<"find it." << endl;
        break;
      }
    }
    //              cout<<"atom i has "<< " " << input->m_Atoms[i]->radius <<" "
    //              << input->m_Atoms[i]->charge << endl;
  }

  savePQR(input, fout);
  return true;
}

void usageAmberAtom2surfaceDistPQR() {
  cout << endl
       << "Usage: ./Molsurf -amberAtom2surfaceDistPQR <inputPQR>  <inputMesh> "
          "<inputAmberOutput> <outputCDR>";
  cout << endl
       << "Info. outputCDR is (Q)charge replace by (D)istance	to the surface "
          "and (R) replaced by  Born radius from Amber output.";
  cout << endl << "Warning: Not a standard PQR output.";
  cout << endl;
}

bool AmberAtom2surfaceDistPQR(int argc, char *argv[]) {
  if (argc != 6) {
    usageAmberAtom2surfaceDistPQR();
    return false;
  }

  char inputPQR[256];
  char inputMesh[256];
  char outputCDR[256];

  strcpy(inputPQR, argv[2]);
  strcpy(inputMesh, argv[3]);
  string inputBorn = argv[4];
  strcpy(outputCDR, argv[5]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fout = fopen(outputCDR, "w");
  FILE *fpA = fopen(inputPQR, "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fpA, radiusType);
  fclose(fpA);

  Geometry *geometry = GeometryLoader().loadFile(inputMesh);
  PG *pg = new PG(10.0, 1000.0, 3.0);

  //      vector <DPG::Point* > meshPoints;
  for (int i = 0; i < geometry->m_NumTriVerts; i++) {
    DPG::Point *p = new DPG::Point(geometry->m_TriVerts[3 * i + 0],
                                   geometry->m_TriVerts[3 * i + 1],
                                   geometry->m_TriVerts[3 * i + 2]);
    //              meshPoints.push_back(p);
    pg->addPoint(p);
  }
  vector<DPG::Point *> results;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    DPG::Point *pt = new DPG::Point(input->m_Atoms[i]->position[0],
                                    input->m_Atoms[i]->position[1],
                                    input->m_Atoms[i]->position[2]);
    float threshold0 = 0.5;
    float threshold = threshold0;
    int j = 1;
    do {
      results = pg->range(pt, threshold);
      j++;
      threshold += threshold0 * j;
    } while (results.size() == 0);

    float dist = pt->distance(results[0]);
    float temp;
    for (int k = 1; k < results.size(); k++) {
      temp = pt->distance(results[k]);
      if (temp < dist)
        dist = temp;
    }
    input->m_Atoms[i]->setCharge(dist);
    int ascii;
    for (int t = 0; t < 6; t++) {
      ascii = (int)(input->m_Atoms[i]->atomName[t]);
      if (ascii <= 90 && ascii >= 65)
        break;
    }
    switch (ascii) {
    case 67:
      input->m_Atoms[i]->setOccupancy(0);
      break;
    case 79:
      input->m_Atoms[i]->setOccupancy(1);
      break;
    case 78:
      input->m_Atoms[i]->setOccupancy(2);
      break;
    case 72:
      input->m_Atoms[i]->setOccupancy(3);
      break;
    default:
      input->m_Atoms[i]->setOccupancy(4);
      break;
    }
  }

  ifstream fpBorn(inputBorn.c_str());
  string tmp;
  while (getline(fpBorn, tmp)) {
    if (strncmp(tmp.c_str(), "STATISTICS OF EFFECTIVE BORN RADII", 34) == 0)
      break;
  }
  readUntilNewline(fpBorn);

  vector<string> text;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    getline(fpBorn, tmp);
    istringstream iss(tmp);
    string s;
    while (iss >> s)
      text.push_back(s);

    if ((i < 9999) && !(atoi(text[0].c_str()) == i + 1)) {
      cout << "Error reading the Born file, exit! i= " << i << " ID is  "
           << text[0] << endl;
      return false;
    }
    input->m_Atoms[i]->setDefinedRadius(atof(text[1].c_str()));
    text.clear();
  }

  saveCDR(input, fout);
  fclose(fout);

  return true;
}

void usageamberMDExtractInfo() {
  cout << endl
       << "Usage: ./Molsurf -amberMDExtractInfo <inputPQR>  <input [force "
          "(.for)], [acceleration (.acc)], [velocity (.vel)]> <output root>";
  cout << endl
       << "Or  ./Molsurf -amberMDExtractInfo <inputPQR> <input velocity "
          "(.vel)> <input [force (.for)], [acceleration (.acc)], [velocity "
          "(.vel)]> <float dt> <output root>";
  cout << endl << "Info. input needs extension, output without extension.";
  cout << endl;
}

bool AmberMDExtractInfo(int argc, char *argv[]) {
  if (argc < 5) {
    usageamberMDExtractInfo();
    return false;
  }

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fpin = fopen(argv[2], "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fpin, radiusType);
  fclose(fpin);

  // FILE * fout = fopen(argv[4], "w");

  ifstream fin(argv[3]);

  string filename(argv[argc - 1]);
  string tmp;

  double acc[3];
  double force[3];
  double vel[3];
  double KineticEnergy = 0.0;

  size_t found;
  size_t founda;
  size_t foundb;

  string inputfile(argv[3]);

  string Force(".for");
  string Veloc(".vel");
  string Accel(".acc");

  if (argc == 5) {
    found = inputfile.find(Force);
    if (found != string::npos) {
      tmp = filename + ".acc";
      ofstream fout(tmp.c_str());
      for (int i = 0; i < input->m_Atoms.size(); i++) {
        fin >> force[0] >> force[1] >> force[2];
        for (int j = 0; j < 3; j++)
          acc[j] = force[j] / input->m_Atoms[i]->mass;
        fout << acc[0] << " " << acc[1] << " " << acc[2] << endl;
      }
      return 1;
    }

    found = inputfile.find(Veloc);
    if (found != string::npos) {
      tmp = filename + ".kin";
      ofstream fout(tmp.c_str());
      for (int i = 0; i < input->m_Atoms.size(); i++) {
        fin >> vel[0] >> vel[1] >> vel[2];
        double tempK = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
        KineticEnergy += tempK * input->m_Atoms[i]->mass * 0.5;
      }
      fout << KineticEnergy << "  kcal/mol." << endl;
      return 1;
    }

    found = inputfile.find(Accel);
    if (found != string::npos) {
      tmp = filename + ".for";
      ofstream fout(tmp.c_str());
      for (int i = 0; i < input->m_Atoms.size(); i++) {
        fin >> acc[0] >> acc[1] >> acc[2];
        for (int j = 0; j < 3; j++)
          force[j] = acc[j] * input->m_Atoms[i]->mass;
        fout << force[0] << " " << force[1] << " " << force[2] << endl;
      }
      return 1;
    }
  } else if (argc == 7) {
    string inputfileb(argv[4]);
    float dt = atof(argv[5]);
    float velb[3];
    float accb[3];
    float forceb[3];

    string tmpb;
    string tmpc;
    string tmpd;

    double kinea = 0.0;
    double kineb = 0.0;

    ifstream finb(argv[4]);

    founda = inputfile.find(Veloc);
    if (founda != string::npos) {
      foundb = inputfileb.find(Veloc);
      if (foundb != string::npos) {
        tmp = filename + ".acc";
        ofstream fout(tmp.c_str());
        tmpb = filename + ".for";
        ofstream foutb(tmpb.c_str());
        tmpc = filename + ".kin";
        ofstream foutc(tmpc.c_str());

        for (int i = 0; i < input->m_Atoms.size(); i++) {
          fin >> vel[0] >> vel[1] >> vel[2];
          finb >> velb[0] >> velb[1] >> velb[2];
          kinea += 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) *
                   input->m_Atoms[i]->mass;
          kineb += 0.5 *
                   (velb[0] * velb[0] + velb[1] * velb[1] + velb[2] * velb[2]) *
                   input->m_Atoms[i]->mass;

          for (int j = 0; j < 3; j++) {
            accb[j] = (velb[j] - vel[j]) / dt;
            forceb[j] = accb[j] * input->m_Atoms[i]->mass;
          }
          fout << accb[0] << " " << accb[1] << " " << accb[2] << endl;
          foutb << forceb[0] << " " << forceb[1] << " " << forceb[2] << endl;
        }
        foutc << kinea << "  " << kineb << "    kcal/mol" << endl;
        return 1;
      }

      foundb = inputfileb.find(Accel);
      if (foundb != string::npos) {
        tmpb = filename + ".for";
        ofstream foutb(tmpb.c_str());
        tmpc = filename + ".kin";
        ofstream foutc(tmpc.c_str());
        tmpd = filename + ".vel";
        ofstream foutd(tmpd.c_str());

        for (int i = 0; i < input->m_Atoms.size(); i++) {
          fin >> vel[0] >> vel[1] >> vel[2];
          finb >> accb[0] >> accb[1] >> accb[2];
          kinea += 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) *
                   input->m_Atoms[i]->mass;
          for (int j = 0; j < 3; j++) {
            velb[j] = vel[j] + accb[j] * dt;
            forceb[j] = input->m_Atoms[i]->mass * accb[j];
          }
          kineb += 0.5 *
                   (velb[0] * velb[0] + velb[1] * velb[1] + velb[2] * velb[2]) *
                   input->m_Atoms[i]->mass;
          foutd << velb[0] << " " << velb[1] << " " << velb[2] << endl;
          foutb << forceb[0] << " " << forceb[1] << " " << forceb[2] << endl;
        }
        foutc << kinea << "  " << kineb << "    kcal/mol" << endl;
        return 1;
      }

      foundb = inputfileb.find(Force);
      if (foundb != string::npos) {
        tmp = filename + ".acc";
        ofstream fout(tmp.c_str());
        tmpc = filename + ".kin";
        ofstream foutc(tmpc.c_str());
        tmpd = filename + ".vel";
        ofstream foutd(tmpd.c_str());

        for (int i = 0; i < input->m_Atoms.size(); i++) {
          fin >> vel[0] >> vel[1] >> vel[2];
          finb >> forceb[0] >> forceb[1] >> forceb[2];
          kinea += 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) *
                   input->m_Atoms[i]->mass;
          for (int j = 0; j < 3; j++) {
            accb[j] = forceb[j] / input->m_Atoms[i]->mass;
            velb[j] = vel[j] + accb[j] * dt;
          }
          kineb += 0.5 *
                   (velb[0] * velb[0] + velb[1] * velb[1] + velb[2] * velb[2]) *
                   input->m_Atoms[i]->mass;
          foutd << velb[0] << " " << velb[1] << " " << velb[2] << endl;
          fout << accb[0] << " " << accb[1] << " " << accb[2] << endl;
        }
        foutc << kinea << "  " << kineb << "    kcal/mol" << endl;
        return 1;
      }
    }

  } else {
    usageamberMDExtractInfo();
    return 0;
  }
}

void usagecvcAtom2surfaceDistPQR() {
  cout << endl
       << "Usage: ./Molsurf -cvcAtom2surfaceDistPQR <inputPQR>  <inputMesh> "
          "<inputCVCBornOutput> <outputCDR>";
  cout << endl
       << "Info. outputCDR is (Q)charge replace by (D)istance	to the surface "
          "and (R) replaced by Born radius from CVC Born radius output.";
  cout << endl << "Warning: Not a standard PQR output.";
  cout << endl;
}

bool CVCAtom2surfaceDistPQR(int argc, char *argv[]) {
  if (argc != 6) {
    usagecvcAtom2surfaceDistPQR();
    return false;
  }

  char inputPQR[256];
  char inputMesh[256];
  char inputBorn[256];
  char outputCDR[256];

  strcpy(inputPQR, argv[2]);
  strcpy(inputMesh, argv[3]);
  strcpy(inputBorn, argv[4]);
  strcpy(outputCDR, argv[5]);

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  FILE *fout = fopen(outputCDR, "w");
  FILE *fpA = fopen(inputPQR, "r");
  FILE *fpBorn = fopen(inputBorn, "r");
  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fpA, radiusType);
  fclose(fpA);

  Geometry *geometry = GeometryLoader().loadFile(inputMesh);
  PG *pg = new PG(10.0, 1000.0, 3.0);

  //      vector <DPG::Point* > meshPoints;
  for (int i = 0; i < geometry->m_NumTriVerts; i++) {
    DPG::Point *p = new DPG::Point(geometry->m_TriVerts[3 * i + 0],
                                   geometry->m_TriVerts[3 * i + 1],
                                   geometry->m_TriVerts[3 * i + 2]);
    //              meshPoints.push_back(p);
    pg->addPoint(p);
  }
  vector<DPG::Point *> results;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    DPG::Point *pt = new DPG::Point(input->m_Atoms[i]->position[0],
                                    input->m_Atoms[i]->position[1],
                                    input->m_Atoms[i]->position[2]);
    float threshold0 = 0.5;
    float threshold = threshold0;
    int j = 1;
    do {
      results = pg->range(pt, threshold);
      j++;
      threshold += threshold0 * j;
    } while (results.size() == 0);

    float dist = pt->distance(results[0]);
    float temp;
    for (int k = 1; k < results.size(); k++) {
      temp = pt->distance(results[k]);
      if (temp < dist)
        dist = temp;
    }
    input->m_Atoms[i]->setCharge(dist);
    int ascii;
    for (int t = 0; t < 6; t++) {
      ascii = (int)(input->m_Atoms[i]->atomName[t]);
      if (ascii <= 90 && ascii >= 65)
        break;
    }
    switch (ascii) {
    case 67:
      input->m_Atoms[i]->setOccupancy(0);
      break;
    case 79:
      input->m_Atoms[i]->setOccupancy(1);
      break;
    case 78:
      input->m_Atoms[i]->setOccupancy(2);
      break;
    case 72:
      input->m_Atoms[i]->setOccupancy(3);
      break;
    default:
      input->m_Atoms[i]->setOccupancy(4);
      break;
    }
  }

  int atmNum;
  fscanf(fpBorn, "%d\n", &atmNum);
  if (atmNum != input->m_Atoms.size()) {
    cout << "Atom numbers in PQR file and Born radii file are not consistent!"
         << endl;
    return false;
  }

  float bornRadius;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    fscanf(fpBorn, "%f\n", &bornRadius);
    input->m_Atoms[i]->setDefinedRadius(bornRadius);
  }

  saveCDR(input, fout);
  fclose(fout);

  return true;
}

void usageChangeCharge4PQR() {
  cout << "MolSurf -changeCharge4PQR <inputPQRA> <inputPQRB> <outputPQRC> "
       << endl;
  cout << "A's charge column will be replaced by charge column of B" << endl;
}

bool ChangeCharge4PQR(int argc, char *argv[]) {
  if (argc != 5) {
    usageChangeCharge4PQR();
    return false;
  }

  char inputPQRA[256];
  char outputPQR[256];

  strcpy(inputPQRA, argv[2]);
  string inputPQRB = argv[3];
  strcpy(outputPQR, argv[4]);

  FILE *fout = fopen(outputPQR, "w");
  FILE *fpA = fopen(inputPQRA, "r");

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParseOptimalPQR(fpA, radiusType);
  fclose(fpA);

  ifstream fpPQRB(inputPQRB.c_str());

  string tmp;

  vector<weiPoint> WeightPoints;
  vector<weiPoint>::iterator wit;
  weiPoint wp;

  vector<string> text;

  while (getline(fpPQRB, tmp)) {
    if (strncmp(tmp.c_str(), "ATOM", 4) == 0) {
      istringstream iss(tmp);
      string s;
      while (iss >> s)
        text.push_back(s);

      wp.p[0] = atof(text[5].c_str());
      wp.p[1] = atof(text[6].c_str());
      wp.p[2] = atof(text[7].c_str());
      wp.rad = atof(text[8].c_str());
      WeightPoints.push_back(wp);

      text.clear();
    }
  }

  int ascii;
  float r, g, b;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    //             cout<< "atom i: " <<i<<" "<<
    //             input->m_Atoms[i]->position[0]<<" "
    //             <<input->m_Atoms[i]->position[1] <<" "
    //             <<input->m_Atoms[i]->position[2] <<" " <<
    //             input->m_Atoms[i]->radius << endl;

    for (wit = WeightPoints.begin(); wit < WeightPoints.end(); wit++) {
      //             cout<< (*wit).p[0] << " " <<  (*wit).p[1] << " " <<
      //             (*wit).p[2] << " " <<  (*wit).rad << endl;

      if (fabs((*wit).p[0] - input->m_Atoms[i]->position[0]) < 0.001 &&
          fabs((*wit).p[1] - input->m_Atoms[i]->position[1]) < 0.001 &&
          fabs((*wit).p[2] - input->m_Atoms[i]->position[2]) < 0.001) {
        input->m_Atoms[i]->setCharge((*wit).rad);
        WeightPoints.erase(wit);
        break;
      }
      if (wit == WeightPoints.end()) {
        cout << "Warning: no match for atom: " << i << "using 1.0 as default"
             << endl;
        input->m_Atoms[i]->setCharge(1.0);
      }
    }

    for (int t = 0; t < 6; t++) {
      ascii = (int)(input->m_Atoms[i]->atomName[t]);
      if (ascii <= 90 && ascii >= 65)
        break;
    }
    switch (ascii) {
    case 67:
      input->m_Atoms[i]->setOccupancy(0);
      break;
    case 79:
      input->m_Atoms[i]->setOccupancy(1);
      break;
    case 78:
      input->m_Atoms[i]->setOccupancy(2);
      break;
    case 72:
      input->m_Atoms[i]->setOccupancy(3);
      break;
    default:
      input->m_Atoms[i]->setOccupancy(4);
      break;
    }
  }
  saveCDR(input, fout);
  fclose(fout);

  return true;
}

void usageSetVDWParametersFromAmber() {
  cout << "MolSurf -setVDWParametersFromAmber <inputPQR> <input aminolib> "
          "<input amino nt lib> <input amino ct lib> <input param file>  "
       << endl;
}

typedef struct _amberVDWPara {
  char name[6];
  double r;
  double e;
} amberVDWPara;

bool SetVDWParametersFromAmber(int argc, char *argv[]) {
  if (argc < 7) {
    usageSetVDWParametersFromAmber();
    return false;
  }

  char inputPQR[256];

  strcpy(inputPQR, argv[2]);

  string aminolib = argv[3];
  string aminontlib = argv[4];
  string aminoctlib = argv[5];

  string parafile = argv[6];

  FILE *ofp = NULL;
  if (argc == 8) {
    ofp = fopen(argv[7], "w");
    if (!ofp) {
      fprintf(stderr, "%s file open fail\n", argv[7]);
      return false;
    }
  }

  ifstream fp(aminolib.c_str());
  ifstream fpnt(aminontlib.c_str());
  ifstream fpct(aminoctlib.c_str());
  ifstream fparam(parafile.c_str());

  FILE *fpPQR = fopen(inputPQR, "r");

  MOLECULE::RADIUS_TYPE radiusType = MOLECULE::PQR_RADIUS;

  MOLECULE::GroupOfAtoms *input = new MOLECULE::GroupOfAtoms();
  input = ParsePQR(fpPQR, radiusType);
  fclose(fpPQR);

  char chain0 = input->m_Atoms[0]->chainID;
  int resIndex0 = input->m_Atoms[0]->getResIndex();

  char currentChain;
  int currentResInd;

  int resNum = 1;
  vector<int> chainLength;
  for (int i = 0; i < input->m_Atoms.size(); i++) {
    currentChain = input->m_Atoms[i]->chainID;
    if (currentChain == chain0) {
      currentResInd = input->m_Atoms[i]->getResIndex();
      if (strcmp(input->m_Atoms[i]->getResName(), "HIS") == 0)
        input->m_Atoms[i]->setResName("HIE");
      if (currentResInd > resIndex0) {
        resNum++;
        resIndex0 = currentResInd;
      }
    } else {
      chainLength.push_back(resNum);
      resIndex0 = input->m_Atoms[i]->getResIndex();
      chain0 = currentChain;
      resNum = 1;
    }
  }
  chainLength.push_back(resNum);

  for (int i = 0; i < chainLength.size(); i++)
    cout << " size : " << chainLength[i] << endl;

  chain0 = input->m_Atoms[0]->chainID;
  resIndex0 = input->m_Atoms[0]->getResIndex();

  int resTmpIndex = input->m_Atoms[0]->getResIndex();

  string tmp;
  vector<string> text;

  char tmpAtomName[15];
  char tmpAtomType[15];

  char currentAtomName[10];
  char currentResName[5];
  int chainIndex = 0;

  map<string, string> myMap;

  bool tag = 0;

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    currentChain = input->m_Atoms[i]->chainID;

    if (currentChain != chain0) {
      chain0 = currentChain;
      currentResInd = input->m_Atoms[i]->getResIndex();
      resIndex0 = input->m_Atoms[i]->getResIndex();
      chainIndex++;
      tag = 0;
    }

    {
      currentResInd = input->m_Atoms[i]->getResIndex();
      if (resTmpIndex != currentResInd) {
        resTmpIndex = currentResInd;
        myMap.clear();
        tag = 0;
      }

      if (currentResInd == resIndex0) {
        strcpy(currentResName, input->m_Atoms[i]->getResName());

        if (tag == 0) {
          while (getline(fpnt, tmp)) {
            if (strncmp(tmp.c_str() + 8, currentResName, 3) == 0)
              break;
          }
          while (getline(fpnt, tmp)) {
            istringstream iss(tmp);
            string s;
            while (iss >> s)
              text.push_back(s);
            if (text.size() == 8) {
              //  It also works.
              //                                                      text[0].erase(text[0].end()-1);
              //                                                      text[0].erase(text[0].begin());
              //                                                      text[1].erase(text[1].end()-1);
              //                                                      text[1].erase(text[1].begin());

              text[0].erase(remove(text[0].begin(), text[0].end(), '"'),
                            text[0].end());
              text[1].erase(remove(text[1].begin(), text[1].end(), '"'),
                            text[1].end());

              myMap[text[0]] = text[1];
            } else {
              tag = 1;
              fpnt.seekg(0, ios::beg);
              break;
            }
            text.clear();
          }
        }
        strcpy(currentAtomName, input->m_Atoms[i]->atomName);
        {
          string ttt = string(currentAtomName);
          ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
          strcpy(input->m_Atoms[i]->atomType, myMap[ttt].c_str());
        }
      } else if (currentResInd == resIndex0 + chainLength[chainIndex] - 1) {
        strcpy(currentResName, input->m_Atoms[i]->getResName());
        if (tag == 0) {
          while (getline(fpct, tmp)) {
            if (strncmp(tmp.c_str() + 8, currentResName, 3) == 0)
              break;
          }
          while (getline(fpct, tmp)) {
            istringstream iss(tmp);
            string s;
            while (iss >> s)
              text.push_back(s);
            if (text.size() == 8) {
              text[0].erase(remove(text[0].begin(), text[0].end(), '"'),
                            text[0].end());
              text[1].erase(remove(text[1].begin(), text[1].end(), '"'),
                            text[1].end());

              myMap[text[0]] = text[1];
            } else {
              tag = 1;
              fpct.seekg(0, ios::beg);
              break;
            }
            text.clear();
          }
        }
        strcpy(currentAtomName, input->m_Atoms[i]->atomName);
        {
          string ttt = string(currentAtomName);
          ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
          strcpy(input->m_Atoms[i]->atomType, myMap[ttt].c_str());
        }

      } else {
        strcpy(currentResName, input->m_Atoms[i]->getResName());
        if (tag == 0) {
          while (getline(fp, tmp)) {
            if (strncmp(tmp.c_str() + 7, currentResName, 3) == 0)
              break;
          }
          while (getline(fp, tmp)) {
            istringstream iss(tmp);
            string s;
            while (iss >> s)
              text.push_back(s);
            if (text.size() == 8) {
              text[0].erase(remove(text[0].begin(), text[0].end(), '"'),
                            text[0].end());
              text[1].erase(remove(text[1].begin(), text[1].end(), '"'),
                            text[1].end());

              myMap[text[0]] = text[1];
            } else {
              tag = 1;
              fp.seekg(0, ios::beg);
              break;
            }
            text.clear();
          }
        }
        strcpy(currentAtomName, input->m_Atoms[i]->atomName);
        {
          string ttt = string(currentAtomName);
          ttt.erase(remove(ttt.begin(), ttt.end(), ' '), ttt.end());
          strcpy(input->m_Atoms[i]->atomType, myMap[ttt].c_str());
        }
      }

      text.clear();
    }

    //              cout<<"atom id: atom name, atomtype,  res id, res name: "<<
    //              i <<" " << input->m_Atoms[i]->atomName << " "
    //              <<" " << input->m_Atoms[i]->atomType << " "<<
    //              input->m_Atoms[i]->getResIndex() << " " <<
    //              input->m_Atoms[i]->getResName() << endl;
  }

  while (getline(fparam, tmp)) {
    if (strncmp(tmp.c_str(), "MOD4      RE", 12) == 0)
      break;
  }
  std::vector<amberVDWPara> params;
  while (getline(fparam, tmp)) {
    istringstream iss(tmp);
    string s;
    while (iss >> s)
      text.push_back(s);
    if (text.size() >= 3) {
      amberVDWPara param;
      strcpy(param.name, text[0].c_str());
      param.r = atof(text[1].c_str());
      param.e = atof(text[2].c_str());
      params.push_back(param);
      text.clear();
    } else {
      text.clear();
      break;
    }
  }

  for (int i = 0; i < input->m_Atoms.size(); i++) {
    char *tmpAtomType = input->m_Atoms[i]->atomType;
    if (tmpAtomType[0] == 'C') {
      if (strcmp(tmpAtomType, "CT") && strcmp(tmpAtomType, "C") &&
          strcmp(tmpAtomType, "C0") && strcmp(tmpAtomType, "Cs") &&
          strcmp(tmpAtomType, "Cl"))
        strcpy(tmpAtomType, "C*");
    }
    if (tmpAtomType[0] == 'N')
      tmpAtomType[1] = '\0';

    for (int j = 0; j < params.size(); j++) {
      if (!strcmp(params[j].name, tmpAtomType)) {
        input->m_Atoms[i]->ljradius = params[j].r;
        input->m_Atoms[i]->ljwellth = params[j].e;
        break;
      }
    }
    //    cout<<"atom id: atom name, atomtype,   LJ r, LJ e: "<< i<< " "<<
    //    input->m_Atoms[i]->atomName << " "<<  input->m_Atoms[i]->atomType <<"
    //    " << input->m_Atoms[i]->ljradius << " "<< input->m_Atoms[i]->ljwellth
    //    << endl;
  }

  if (ofp) {
    fprintf(ofp, "%d\n", input->m_Atoms.size());
    for (int i = 0; i < input->m_Atoms.size(); i++) {
      fprintf(ofp, "%lf %lf\n", input->m_Atoms[i]->ljradius,
              input->m_Atoms[i]->ljwellth);
      // fprintf( stderr, "%lf %lf\n", input->m_Atoms[i]->ljradius,
      // input->m_Atoms[i]->ljwellth );
    }
    fclose(ofp);
    fprintf(stderr, "%s file dumped\n", argv[7]);
  }
}

#endif

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printUsage();
  }

  else if (strcmp(argv[1], "-area") == 0) {
    area(argc, argv);
  }

  else if (strcmp(argv[1], "-aSpline") == 0) {
    aSpline(argc, argv);
  }

  else if (strcmp(argv[1], "-bornRadius") == 0) {
    bornRadius(argc, argv);
  }

  else if (strcmp(argv[1], "-coarseGrain") == 0) {
    coarseGrain(argc, argv);
  }

  else if (strcmp(argv[1], "-coarseGrainCharge") == 0) {
    coarseGrainCharge(argc, argv);
  }

  else if (strcmp(argv[1], "-convertToRawnUsingHLS") == 0) {
    convertToRawnUsingHLS(argc, argv);
  }

  else if (strcmp(argv[1], "-decimate") == 0) {
    decimate(argc, argv);
  }

  else if (strcmp(argv[1], "-decimatePercentage") == 0) {
    decimatePercentage(argc, argv);
  }

  else if (strcmp(argv[1], "-hausdorff") == 0) {
    hausdorff(argc, argv);
  }

  else if (strcmp(argv[1], "-normals") == 0) {
    normals(argc, argv);
  } else if (strcmp(argv[1], "-reoriente") == 0) {
    reOrientate(argc, argv);
  }

  else if (strcmp(argv[1], "-optimizeCGCharge") == 0) {
    optimizeCGCharge(argc, argv);
  }

  else if (strcmp(argv[1], "-optimizeCGGeometry") == 0) {
    optimizeCGGeometry(argc, argv);
  }

  else if (strcmp(argv[1], "-pocket") == 0) {
    pocket(argc, argv);
  }

  else if (strcmp(argv[1], "-pocketTunnel") == 0) {
    pocketTunnel(argc, argv);
  }

  else if (strcmp(argv[1], "-populateSAS") == 0) {
    populateSAS(argc, argv);
  } else if (strcmp(argv[1], "-populateSASUsingMesh") == 0) {
    populateSASUsingMesh(argc, argv);
  } else if (strcmp(argv[1], "-getInterfaceAtoms") == 0) {
    getInterfaceAtoms(argc, argv);
  } else if (strcmp(argv[1], "-getInterfaceSurface") == 0) {
    getInterfaceSurface(argc, argv);
  } else if (strcmp(argv[1], "-generateF2d") == 0) {
    generateF2d(argc, argv);
  } else if (strcmp(argv[1], "-offsetSurface") == 0) {
    offsetSurface(argc, argv);
  } else if (strcmp(argv[1], "-colorByCharge") == 0) {
    colorByCharge(argc, argv);
  } else if (strcmp(argv[1], "-leastRMSD") == 0) {
    getLeastRMSD(argc, argv);
  } else if (strcmp(argv[1], "-qualityImprove") == 0) {
    qualityImprove(argc, argv);
  } else if (strcmp(argv[1], "-removeInteriorPockets") == 0) {
    removeInteriorPockets(argc, argv);
  } else if (strcmp(argv[1], "-partitionMesh") == 0) {
    partitionMesh(argc, argv);
  } else if (strcmp(argv[1], "-repairSurface") == 0) {
    repairSurface(argc, argv);
  } else if (strcmp(argv[1], "-rotatePDB") == 0) {
    rotatePDB(argc, argv);
  } else if (strcmp(argv[1], "-signedDistanceFunction") == 0) {
    signedDistanceFunction(argc, argv);
  } else if (strcmp(argv[1], "-signedDistanceFunction2") == 0) {
    signedDistanceFunction2(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceAtoms") == 0) {
    surfaceAtoms(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceFromPDB") == 0) {
    surfaceFromPDB(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceUsingAdaptiveGrid") == 0) {
    surfaceUsingAdaptiveGrid(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceUsingFastSurf") == 0) {
    surfaceUsingFastSurf(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceUsingGaussianBlur") == 0) {
    surfaceUsingGaussianBlur(argc, argv);
  } else if (strcmp(argv[1], "-surfaceUsingGaussianBlur2") == 0) {
    surfaceUsingGaussianBlur2(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceUsingHLS") == 0) {
    surfaceUsingHLS(argc, argv);
  } else if (strcmp(argv[1], "-surfaceUsingCardinalBSpline") == 0) {
    surfaceUsingCardinalBSpline(argc, argv);
  }

  else if (strcmp(argv[1], "-surfaceUsingLBIE") == 0) {
    surfaceUsingLBIE(argc, argv);
  }

  else if (strcmp(argv[1], "-volume") == 0) {
    volume(argc, argv);
  }

  else if (strcmp(argv[1], "-volumeUsingHLS") == 0) {
    volumeUsingHLS(argc, argv);
  } else if (strcmp(argv[1], "-getVolumeWithHydrophobicity") == 0) {
    getVolumeWithHydrophobicity(argc, argv);
  } else if (strcmp(argv[1], "-getVolumeWithElectron") == 0) {
    getVolumeWithElectron(argc, argv);
  } else if (strcmp(argv[1], "-getVolumeWithChargeDensity") == 0) {
    getVolumeWithChargeDensity(argc, argv);
  } else if (strcmp(argv[1], "-mergeGeometry") == 0) {
    mergeGeometry(argc, argv);
  }
#ifdef WITH_CG
  else if (strcmp(argv[1], "-mergeProtein") == 0) {
    mergeProtein(argc, argv);
  } else if (strcmp(argv[1], "-splitPDB2Chains") == 0) {
    splitPDB2Chains(argc, argv);
  } else if (strcmp(argv[1], "-PQR2PDB") == 0) {
    pqr2pdb(argc, argv);
  } else if (strcmp(argv[1], "-PQR2XYZR") == 0) {
    pqr2xyzr(argc, argv);
  }

  else if (strcmp(argv[1], "-colorMeshByAtomType") == 0) {
    ColorMeshByAtomType(argc, argv);
  } else if (strcmp(argv[1], "-atom2surfaceDistance") == 0) {
    Atom2SurfaceDistance(argc, argv);
  } else if (strcmp(argv[1], "-tinkerPDBXYZ2PQR") == 0) {
    TinkerPDBXYZ2PQR(argc, argv);
  } else if (strcmp(argv[1], "-changeCharge4PQR") == 0) {
    ChangeCharge4PQR(argc, argv);
  } else if (strcmp(argv[1], "-amberAtom2surfaceDistPQR") == 0) {
    AmberAtom2surfaceDistPQR(argc, argv);
  } else if (strcmp(argv[1], "-vdwvolume") == 0) {
    // VDWVolume( argc, argv );
  } else if (strcmp(argv[1], "-cvcAtom2surfaceDistPQR") == 0) {
    CVCAtom2surfaceDistPQR(argc, argv);
  } else if (strcmp(argv[1], "-setVDWParametersFromAmber") == 0) {
    SetVDWParametersFromAmber(argc, argv);
  } else if (strcmp(argv[1], "-amberMDExtractInfo") == 0) {
    AmberMDExtractInfo(argc, argv);
  }
#endif
  else if (strcmp(argv[1], "-h") == 0) {
    printUsage();
  }
  /*      else if(strcmp(argv[1], "-centroid") == 0)
          {
                  atomCentroid(argc, argv);
          }
          else if(strcmp(argv[1], "-atomBbox") == 0)
          {
                  atomsBoundingbox(argc, argv);
          }*/
  else {
    printUsage();
    return -1;
  }
  return 0;
}

void printUsage() {

  cout << endl << "Usage:  MolSurf -<command> [arguments]" << endl << endl;

  cout << "      MolSurf -<command> " << endl;
  cout << "        produces a description and argument listing" << endl
       << "        for that command." << endl;

  cout << endl;
  cout << "List of commands:" << endl;
  cout << "      surfaceUsingAdaptiveGrid, surfaceUsingGaussianBlur" << endl;
  cout << "      surfaceUsingHSL, surfaceUsingLBIE, surfaceUsingFastSurf "
       << endl;
  cout << "      area, volume, qualityImprove, removeInteriorPockets, "
          "normals,reoriente, aSpline"
       << endl;
  cout << "      offsetSurface, mergeGeometry" << endl;
  cout << "      volumeUsingHLS, getVolumeWithElectron, "
          "getVolumeWithHydrophobicity, getVoluemWithChargeDensity"
       << endl;
  cout << "if WITH_CG, more commands,  mergeProtein, splitPDB2Chains, "
          "PQR2PDB,PQR2XYZR, colorMeshByAtomType, atom2surfaceDistance"
       << endl;
  cout << " tinkerPDBXYZ2PQR,changeCharge4PQR, amberAtom2surfaceDistPQR, "
          "cvcAtom2surfaceDistPQR, setVDWParametersFromAmber "
       << endl;
  cout << " amberMDExtractInfo" << endl;
  cout << " centroid" << endl;
  cout << " atomBbox" << endl;
  cout << endl;

  /*    cout    << endl << "Usage:  MolSurf -h | --help"
                << endl << "    Display this summary."
                << endl ;
  */
  /*      usageArea();
          usageASplineQuad();
          usageASplineRaw();
          usageASplineRawQuad();
          usageASplineVolume();
          usageBornRadius();
          usageCoarseGrain();
          usageCoarseGrainCharge();
          usageConvertToRawnUsingHLS();
          usageDecimate();
          usageDecimatePercentage();
          usageNormals();
          usageOptimizeCGCharge();
          usageOptimizeCGGeometry();
          usagePocket();
          usagePocketTunnel();
          usagePopulateSAS();
          usageQualityImprove();
          usageSignedDistanceFunction();
          usageSurfaceAtoms();
          usageSurfaceFromPDB();
          usageSurfaceUsingAdaptiveGrid();
          usageSurfaceUsingGaussianBlur();
          usageSurfaceUsingHLS();
          usageSurfaceUsingLBIE();
          usageVolume();
          usageVolumeUsingHLS();
  */
}
