#include <Utility/utility.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/surface.h>
#include <libCG/CoarseGrain/geometry.h>

using namespace MOLECULE;
using namespace GEOMETRY;

Atom::Atom() {
  indexInPQR = -1;
  //     indexInBackbone = -1;
  //     indexInResidue = -1;
  //     residueID = -1;
  //     residueSeq = -1;
  chainID = 0;
  resSeq = 0;
  position[0] = 0;
  position[1] = 0;
  position[2] = 0;
  radius = 0;
  charge = 0;
  atomID = NULL;
  mass = 0;
  child = 0;
  occupancy = 1.00;
  tempFactor = 20.00;
  ljradius = 0.0;
  ljwellth = 0.0;
  //     overlapNeighbors = new int[5000];
}

Atom::Atom(Atom *a) {
  indexInPQR = a->indexInPQR;
  //     indexInBackbone = a->indexInBackbone;
  //     indexInResidue = a->indexInResidue;
  //     residueID = a->residueID;
  //     residueSeq = a->residueSeq;
  strcpy(atomName, a->atomName);
  strcpy(resName, a->resName);
  resSeq = a->resSeq;
  position[0] = a->position[0];
  position[1] = a->position[1];
  position[2] = a->position[2];
  radius = a->radius;
  charge = a->charge;
  atomID = a->atomID;
  chainID = a->chainID;
  rho = a->rho;
  mass = a->mass;
  child = a->child;
  occupancy = a->occupancy;
  tempFactor = a->tempFactor;
  ljradius = a->ljradius;
  ljwellth = a->ljwellth;
  strcpy(atomType, a->atomType);

  //     overlapNeighbors = a->overlapNeighbors;
}

Atom::~Atom() {
  if (rho)
    delete rho;
  if (child)
    delete[] child;
}

void Atom::setAtomID(char *atom_name) {
  if ((strstr(atom_name, "H")) != NULL) {
    atomID = 'H';
  } else if ((strstr(atom_name, "C")) != NULL) {
    atomID = 'C';
  } else if ((strstr(atom_name, "N")) != NULL) {
    atomID = 'N';
  } else if ((strstr(atom_name, "O")) != NULL) {
    atomID = 'O';
  } else if ((strstr(atom_name, "S")) != NULL) {
    atomID = 'S';
  } else if ((strstr(atom_name, "P")) != NULL) {
    atomID = 'P';
  } else {
    cout << indexInPQR << " Unknown Atom! " << endl;
    exit(0);
  }
}

void Atom::setVDWradius() {
  switch (atomID) {
  case 'N':
    radius = 1.55;
    break;
  case 'C':
    radius = 1.7;
    break;
  case 'O':
    radius = 1.4;
    break;
  case 'H':
    radius = 1.2;
    break;
  case 'S':
    radius = 1.85;
    break;
  case 'P':
    radius = 1.9;
    break;
  default:
    break;
  }
}

void Atom::setMass() {
  switch (atomID) {
  case 'N':
    mass = 14.01;
    break;
  case 'C':
    mass = 12.01;
    break;
  case 'O':
    mass = 16.0;
    break;
  case 'H':
    mass = 1.008;
    break;
  case 'S':
    mass = 32.06;
    break;
  case 'P':
    mass = 30.97;
    break;
  default:
    break;
  }
}

void Atom::setAtomName(char *atom_name) { strcpy(atomName, atom_name); }

void Atom::setResName(char *res_name) { strcpy(resName, res_name); }

void Atom::setChainID(char chain_id) { chainID = chain_id; }

void Atom::setIndexInPQR(int i) { indexInPQR = i; }

void Atom::setResIndex(int res_seq) { resSeq = res_seq; }

void Atom::setPosition(double *p) {
  position[0] = p[0];
  position[1] = p[1];
  position[2] = p[2];
}

void Atom::setCharge(double q) { charge = q; }

void Atom::setDefinedRadius(double r) { radius = r; }

void Atom::setOccupancy(float occu) { occupancy = occu; }

void Atom::setTempFactor(float temp) { tempFactor = temp; }

char *Atom::getAtomType() { return atomType; }
char *Atom::getResName() { return resName; }
int Atom::getResIndex() { return resSeq; }
float Atom::getLJradius() { return ljradius; }
float Atom::getLJwellth() { return ljwellth; }

void Atom::getData(int i, char *atom_name, char *res_name, char chain_id,
                   int res_seq, double *p, double q, double r,
                   MOLECULE::RADIUS_TYPE radiusType, float occupancy,
                   float tempFactor) {
  setIndexInPQR(i);
  setAtomName(atom_name);
  setAtomID(atom_name);
  setResName(res_name);
  setChainID(chain_id);
  setResIndex(res_seq);
  setPosition(p);
  setCharge(q);
  if (radiusType == VDW_RADIUS) {
    setVDWradius();
  }
  if (radiusType == PQR_RADIUS) {
    setDefinedRadius(r);
  }
  setMass();
  setOccupancy(occupancy);
  setTempFactor(tempFactor);
  getDensityInfo();
}

void Atom::getData(int i, char *atom_name, char *res_name, char chain_id,
                   int res_seq, double *p, double q, double r,
                   MOLECULE::RADIUS_TYPE radiusType) {
  setIndexInPQR(i);
  setAtomName(atom_name);
  setAtomID(atom_name);
  setResName(res_name);
  setChainID(chain_id);
  setResIndex(res_seq);
  setPosition(p);
  setCharge(q);
  if (radiusType == VDW_RADIUS) {
    setVDWradius();
  }
  if (radiusType == PQR_RADIUS) {
    setDefinedRadius(r);
  }
  setMass();
  getDensityInfo();
}

void Atom::getData(double *p, double q, double r, char *res_name) {
  setPosition(p);
  setCharge(q);
  setDefinedRadius(r);
  setResName(res_name);
}

void Atom::getDensityInfo() {
  rho = new density;
  const double w = 1.4;
  double x = w / radius;
  rho->para_c = 10.0 / (2.0 * x * x * x + 9.0 * x * x + 15.0 * x + 10.0);
  rho->coeffOfr[0] = 1.0 - 3.0 / (x * x) - 2.0 / (x * x * x);
  rho->coeffOfr[1] = 6.0 / (w * x) * (1.0 + 1.0 / x);
  rho->coeffOfr[2] = -3.0 / (w * w) * (1.0 + 2.0 / x);
  rho->coeffOfr[3] = 2.0 / (w * w * w);
  rho->coeffOfr[0] *= rho->para_c;
  rho->coeffOfr[1] *= rho->para_c;
  rho->coeffOfr[2] *= rho->para_c;
  rho->coeffOfr[3] *= rho->para_c;
}

void Atom::getBornRadius(Surface *surface) {
  const double pi = 3.1415926;
  double xx[3], n[3], nn, w, sum1, sum2;
  double t, b, b4, b7;
  double A4, A7;
  int i, j, k;
  sum1 = 0.0;
  sum2 = 0.0;
  /* Fast Summation */
  // GBMV model Wenqi added on 2/28/07
  const double S = 0.9114, C0 = 0.2966, C1 = 1.0369, D = -0.0637;
  /* Wenqi modified on 4/23 to avoid NAN */
  (sum1 > 0.0) ? A4 = sum1 / (4.0 * pi) : A4 = 0.0;
  //     printf("A4=%f\n", A4);
  (sum2 > 0.0) ? A7 = sum2 / (16.0 * pi) : A7 = 0.0;
  A7 = sqrt(A7);
  A7 = sqrt(A7);
  //     printf("A7=%f\n", A7);
  Born_radius = 1.0 / ((1.0 - 1.0 / sqrt(2.0)) * A4 + A7); // this is GBSW model
  //    Born_radius = S / ( C0*A4 + C1*A7 ) + D; // GBMV
  if (Born_radius < 0.0 || Born_radius == 0.0) {
    printf("%d %f\n", indexInPQR, Born_radius);
  }
}

double Atom::DotProduct(double *v1, double *v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
