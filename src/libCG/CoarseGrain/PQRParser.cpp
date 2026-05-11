#include <Utility/utility.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/PQRParser.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>

using namespace MOLECULE;

PQRParser::PQRParser() {}

GroupOfAtoms *PQRParser::ParserPDBtoGOA(FILE *fp,
                                        MOLECULE::RADIUS_TYPE radiusType) {
  char atomname[10];
  int serial;
  char name[10];
  char altLoc;
  char resName[10];
  char chainID;
  int resSeq;
  char iCode;
  double p[3];
  float occupancy = 1.0;
  float tempFactor = 20.0;
  double q;
  double r;
  char segID[10];
  char element[10];
  char charge[10];
  char line[82];
  int numofres = 0;
  int resSeq0 = -1;
  int resAdd = 0;
  char iCode0 = ' ';
  char chainID0 = ' ';
  GroupOfAtoms *molecule = new GroupOfAtoms();
  int atomIndex = 0;

  while (fgets(line, 81, fp) != 0) {
    if ((strncmp(line, "ATOM", 4) == 0) && (strncmp(line, "REMARK", 6) != 0)) {
      getString(line, atomname, 0, 6);
      getInt(line, &serial, 6, 5);
      getString(line, name, 12, 4);
      getChar(line, &altLoc, 16);
      getString(line, resName, 17, 3);
      getChar(line, &chainID, 21);
      getInt(line, &resSeq, 22, 4);
      getChar(line, &iCode, 26);
      getDouble(line, &p[0], 30, 8);
      getDouble(line, &p[1], 38, 8);
      getDouble(line, &p[2], 46, 8);
      getFloat(line, &occupancy, 54, 6);
      getFloat(line, &tempFactor, 60, 6);
      getString(line, segID, 72, 4);
      getString(line, element, 76, 2);
      getString(line, charge, 78, 2);
      if (chainID != chainID0) {
        resAdd = 0;
        chainID0 = chainID;
      }
      if (iCode != iCode0) {
        iCode0 = iCode;
        resAdd++;
      }
      resSeq = resSeq + resAdd / 2 + resAdd % 2;
      Atom *a = new Atom();
      a->getData(atomIndex, name, resName, chainID, resSeq, p, q, r, radiusType,
                 occupancy, tempFactor); // translate resName to resID
      molecule->addAtom(a);
      atomIndex++;
      // count number of residues
      if (resSeq != resSeq0) {
        numofres++;
        resSeq0 = resSeq;
      }
    }
  }
  molecule->numOfAtoms = atomIndex;
  molecule->numOfRes = numofres;
  return molecule;
}

GroupOfAtoms *PQRParser::ParserPQRtoGOA(FILE *fp,
                                        MOLECULE::RADIUS_TYPE radiusType) {
  char atomname[10];
  int serial;
  char name[10];
  char resName[10];
  char chainID = ' ';
  int resSeq;
  double p[3];
  double q;
  double r;
  char line[82];
  int numofres = 0;
  int resSeq0 = -1;
  GroupOfAtoms *molecule = new GroupOfAtoms();
  int atomIndex = 0;
  char chainID0 = 'A';

  while (fgets(line, 81, fp) != 0) {
    if ((strncmp(line, "ATOM", 4) == 0) && (strncmp(line, "REMARK", 6) != 0)) {
      getString(line, atomname, 0, 6);
      getInt(line, &serial, 6, 5);
      getString(line, name, 12, 4);
      getString(line, resName, 17, 3);
      getChar(line, &chainID, 21);
      getInt(line, &resSeq, 22, 4);
      getDouble(line, &p[0], 30, 8);
      getDouble(line, &p[1], 38, 8);
      getDouble(line, &p[2], 46, 8);
      getDouble(line, &q, 54, 8);
      getDouble(line, &r, 63, 8);
      if (chainID == ' ') {
        chainID = chainID0;
        if (resSeq < resSeq0) {
          chainID++;
          chainID0 = chainID;
        }
      }
      Atom *a = new Atom();
      a->getData(atomIndex, name, resName, chainID, resSeq, p, q, r,
                 radiusType); // translate resName to resID
      molecule->addAtom(a);
      atomIndex++;
      // count number of residues
      if (resSeq != resSeq0) {
        numofres++;
        resSeq0 = resSeq;
      }
    }
  }
  molecule->numOfAtoms = atomIndex;
  molecule->numOfRes = numofres;
  return molecule;
}

GroupOfAtoms *
PQRParser::ParserOptimalPQRtoGOA(FILE *fp, MOLECULE::RADIUS_TYPE radiusType) {
  char atomname[10];
  int serial;
  char name[10];
  char resName[10];
  char chainID;
  int resSeq;
  double p[3];
  double q;
  double r;
  char line[82];
  int numofres = 0;
  int resSeq0 = -1;
  GroupOfAtoms *molecule = new GroupOfAtoms();
  int atomIndex = 0;
  while (fgets(line, 81, fp) != 0) {
    if ((strncmp(line, "ATOM", 4) == 0) && (strncmp(line, "REMARK", 6) != 0)) {
      getString(line, atomname, 0, 6);
      getInt(line, &serial, 6, 5);
      getString(line, name, 12, 4);
      getString(line, resName, 17, 3);
      getChar(line, &chainID, 21);
      getInt(line, &resSeq, 22, 4);
      getDouble(line, &p[0], 30, 8);
      getDouble(line, &p[1], 38, 8);
      getDouble(line, &p[2], 46, 8);
      getDouble(line, &q, 54, 9);
      getDouble(line, &r, 63, 8);
      Atom *a = new Atom();
      a->getData(atomIndex, name, resName, chainID, resSeq, p, q, r,
                 radiusType); // translate resName to resID
      molecule->addAtom(a);
      atomIndex++;
      if (resSeq != resSeq0) {
        numofres++;
        resSeq0 = resSeq;
      }
    }
  }
  molecule->numOfAtoms = atomIndex;
  molecule->numOfRes = numofres;
  return molecule;
}

void PQRParser::getString(char *line, char *str, int pos, int len) {
  char buffer[80];
  assert(len < 80);
  strcpy(str, "");
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  int i;
  for (i = 0; i < len; i++) {
    str[i] = buffer[i];
  }
  str[len] = '\0';
}

void PQRParser::getChar(char *line, char *str, int pos) {
  char buffer[80];
  *str = ' ';
  strncpy(buffer, &line[pos], 1);
  buffer[1] = '\0';
  sscanf(buffer, "%c", str);
}

void PQRParser::getInt(char *line, int *i, int pos, int len) {
  char buffer[80];
  *i = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%d", i);
}

void PQRParser::getDouble(char *line, double *d, int pos, int len) {
  char buffer[80];
  *d = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%lf", d);
}

void PQRParser::getFloat(char *line, float *d, int pos, int len) {
  char buffer[80];
  *d = -1;
  strncpy(buffer, &line[pos], len);
  buffer[len] = '\0';
  sscanf(buffer, "%f", d);
}
