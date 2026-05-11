#include <PDBParser/Atom.h>
#include <DPG/PG.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <PDBParser/ElementInformation.h>

class SurfaceHydrophobicity {
private:
  bool dataReady;
  vector<PDBParser::Atom *> atms;
  vector<double> hydrophobicityContributionPerAtom;
  vector<DPG::WeightedPoint *> quadPoints;

  double xlatePG;
  DPG::PG *pgPoints;

  int *hydroIndex;
  int nHydroIndex;

  int strcmp_nospace(const char *s1, const char *s2);
  void preprocessElementInformationTable(int **index, int *nIndex);
  float getHydrophobicity(const char *atomName, const char *residueName,
                          int *index, int nIndex, bool useInterfacePropensity,
                          bool perResidueHydrophobicity);

  void computeXlateForPG();
  void buildPG();
  bool readPQR(string pqrFileName);
  bool readQuad(string quadFileName);

public:
  SurfaceHydrophobicity(string pqrFileName, string quadFileName);
  ~SurfaceHydrophobicity();

  double getSurfaceHydrophobicity();
};
