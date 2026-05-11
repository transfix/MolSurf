#ifndef PQRPARSER_H
#define PQRPARSER_H

#include "atom2.h"
#include <stdio.h>

namespace MOLECULE {
class Atom;
class GroupOfAtoms;
class PQRParser {
public:
  PQRParser();
  GroupOfAtoms *ParserPQRtoGOA(FILE *fp, MOLECULE::RADIUS_TYPE radiusType);
  GroupOfAtoms *ParserPDBtoGOA(FILE *fp, MOLECULE::RADIUS_TYPE radiusType);
  GroupOfAtoms *ParserOptimalPQRtoGOA(FILE *fp,
                                      MOLECULE::RADIUS_TYPE radiusType);

private:
  void getString(char *line, char *str, int pos, int len);
  void getChar(char *line, char *str, int pos);
  void getInt(char *line, int *i, int pos, int len);
  void getDouble(char *line, double *d, int pos, int len);
  void getFloat(char *line, float *d, int pos, int len);
};
} // namespace MOLECULE

#endif // PQRPARSER_H
