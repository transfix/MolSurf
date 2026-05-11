#ifndef ___ATOM_H
#define ___ATOM_H

#include <libCG/CoarseGrain/surface.h>

namespace MOLECULE {
typedef struct densityFunc {
  double para_c;
  double coeffOfr[4];
} density;

enum RADIUS_TYPE { VDW_RADIUS, PQR_RADIUS };

class Atom {
public:
  Atom();
  Atom(Atom *);
  virtual ~Atom();
  void setAtomID(char *atom_name);
  void setVDWradius();
  void setMass();
  void setOccupancy(float);
  void setTempFactor(float);
  void setAtomName(char *atom_name);
  void setResName(char *res_name);
  void setChainID(char chain_id);
  void setIndexInPQR(int i);
  void setResIndex(int res_seq);
  void setPosition(double *p);
  void setCharge(double q);
  void setDefinedRadius(double r);
  char *getAtomType();
  float getLJradius();
  float getLJwellth();
  int getResIndex();
  char *getResName();

  void getData(int i, char *atom_name, char *res_name, char chainID,
               int res_seq, double *p, double q, double r,
               RADIUS_TYPE radiusType, float occupancy, float tempFactor);
  void getData(int i, char *atom_name, char *res_name, char chainID,
               int res_seq, double *p, double q, double r,
               RADIUS_TYPE radiusType);

  void getData(double *p, double q, double r, char *res_name);
  void getBornRadius(GEOMETRY::Surface *surface);

  int indexInPQR; // start from 0
  char atomID;    // name

  char atomName[10];
  char resName[5];
  char chainID;
  int resSeq;

  double position[3];
  double radius;
  double charge;
  double Born_radius;
  float occupancy;
  float tempFactor;
  char atomType[10];
  double ljradius;
  double ljwellth;

  double mass;
  int *child, numOfChildren;

  density *rho;

private:
  void getDensityInfo();
  double DotProduct(double *, double *);
};
}; // namespace MOLECULE

#endif //___ATOM_H
