
#include "Octree.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

#define VDW
// #define GAUSS

using namespace std;

vector<double> mvx;
vector<double> mvy;
vector<double> mvz;
vector<int> tri0;
vector<int> tri1;
vector<int> tri2;

double getTime(void) {
#ifdef _WIN32
  time_t ltime;
  _timeb tstruct;
  time(&ltime);
  _ftime(&tstruct);
  return (double)(ltime + 1e-3 * (tstruct.millitm));
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)(t.tv_sec + 1e-6 * t.tv_usec);
#endif
}

//
// Octree surface
//
// TODO: Gaussian surface
//       Solvent Exposed Surface
//
//       Fix raw output (eliminate duplicate vertices)
//       Print normals?
//       Fix pqr reader...
//       ...
//

#include "MarchingCubes.cpp"
#include "ImprovedMarchingCubes.cpp"

Octree *OTNode::tree = NULL;

// old
/*
double f(double * xyz) {
  return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
}
*/

double Octree::fvdw(OTNode *nd) {

  // double * xyz = nd->xyz;
  // return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] -400;

  double mydist = -50.0;
  bool inside = false;
  int size = nd->atoms.size();

  for (int iA = 0; iA < size; iA++) {
    double distsq = (nd->xyz[0] - atoms[nd->atoms[iA]]->xyz[0]) *
                        (nd->xyz[0] - atoms[nd->atoms[iA]]->xyz[0]) +
                    (nd->xyz[1] - atoms[nd->atoms[iA]]->xyz[1]) *
                        (nd->xyz[1] - atoms[nd->atoms[iA]]->xyz[1]) +
                    (nd->xyz[2] - atoms[nd->atoms[iA]]->xyz[2]) *
                        (nd->xyz[2] - atoms[nd->atoms[iA]]->xyz[2]);

    /*
                    double val = sqrt(distsq) - atoms[nd->atoms[iA]].r - 1.4;
                    if(val<mydist)
                            mydist = val;
    */

    double val = atoms[nd->atoms[iA]]->r + 1.4 - sqrt(distsq);
    if (val > mydist)
      mydist = val;
  }
  return mydist;
}

double Octree::fgauss(OTNode *nd) {

  // double * xyz = nd->xyz;
  // return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] -400;

  double val = 0.0;
  //  double blobby = -2.3;

  for (int iA = 0; iA < nd->atoms.size(); iA++) {
    double distsq = (nd->xyz[0] - atoms[nd->atoms[iA]]->xyz[0]) *
                        (nd->xyz[0] - atoms[nd->atoms[iA]]->xyz[0]) +
                    (nd->xyz[1] - atoms[nd->atoms[iA]]->xyz[1]) *
                        (nd->xyz[1] - atoms[nd->atoms[iA]]->xyz[1]) +
                    (nd->xyz[2] - atoms[nd->atoms[iA]]->xyz[2]) *
                        (nd->xyz[2] - atoms[nd->atoms[iA]]->xyz[2]);
    double expval =
        blobby * distsq / (atoms[nd->atoms[iA]]->r * atoms[nd->atoms[iA]]->r) -
        blobby;
    val += (exp(expval));
  }

  return val;
}

void Octree::fvdw_update_add(OTNode *nd, Atom *a) {
  double mydist = 10.0;
  double dx = nd->xyz[0] - a->xyz[0];
  double dy = nd->xyz[1] - a->xyz[1];
  double dz = nd->xyz[2] - a->xyz[2];

  double distsq = dx * dx + dy * dy + dz * dz;

  double val = sqrt(distsq) - a->r;
  if (val < nd->v)
    nd->v = val;
}

void Octree::fgauss_update_add(OTNode *nd, Atom *a) {
  double val = 0.0;
  //  double blobby = -2.3;

  double dx = nd->xyz[0] - a->xyz[0];
  double dy = nd->xyz[1] - a->xyz[1];
  double dz = nd->xyz[2] - a->xyz[2];
  double distsq = dx * dx + dy * dy + dz * dz;

  double expval = blobby * distsq / (a->r * a->r) - blobby;

  nd->v += (exp(expval));
}

void Octree::fgauss_update_remove(OTNode *nd, Atom *a) {
  double val = 0.0;
  //  double blobby = -2.3;

  double dx = nd->xyz[0] - a->xyz[0];
  double dy = nd->xyz[1] - a->xyz[1];
  double dz = nd->xyz[2] - a->xyz[2];
  double distsq = dx * dx + dy * dy + dz * dz;

  double expval = blobby * distsq / (a->r * a->r) - blobby;

  nd->v -= (exp(expval));
}

OTNode::OTNode() {

  c = NULL;
  flag = -1;
  splitDisabled = false;
  sdfComputed = false;
  markedForLocalUpdate = false;
  v = -100000.0;
}

OTNode::OTNode(int l1, double x, double y, double z, double w1) {
  c = NULL;
  l = l1;
  w = w1;
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;

  std::vector<int> atoms;

  flag = -1;
  splitDisabled = false;
  interfaceFlag = false;
  sdfComputed = false;
  markedForLocalUpdate = false;
  v = -100000.0;
}

void OTNode::SetData(double v1) { v = v1; }

/*
void OTNode::Merge() {
  //checking if merging is valid. This cell must have children which are all
leaves if(c==NULL) return; for (int i=1; i<8; i++) if(c[i].c!=NULL) return;

  //checking if this cell need to be merged
  bool allAbove = c[0]->v > contourValue;
  for (int i=1; i<8; i++) {
    if(allAbove && c[i]->v < contourValue )
      return; // no need to merge
    else if(!allAbove && c[i]->v > contourValue)
      return;
  }

  //merging and deleting the children
  atoms.clear();
  for (int i=0; i<8; i++) {
    int numAtoms = c[i]->atoms.size();
    for (int iA=0; iA<numAtoms; iA++)
      atoms.push_back(c[i]->atoms[iA]);
    pg->removePoint(c[i]);
    delete c[i];
  }
  delete c;
//  pg->addPoint(c[i]);

  //merging the parent if necessary
  Merge(parent);
}
*/

void OTNode::Split() {

  if (c != NULL) {
    cout << "WARNING: Split() called on an already split node." << endl;
    return;
  }

  c = new OTNode *[8];

  c[0] = new OTNode(l + 1, xyz[0] - w / 4.0, xyz[1] - w / 4.0, xyz[2] - w / 4.0,
                    w / 2.0);
  c[1] = new OTNode(l + 1, xyz[0] - w / 4.0, xyz[1] - w / 4.0, xyz[2] + w / 4.0,
                    w / 2.0);
  c[2] = new OTNode(l + 1, xyz[0] - w / 4.0, xyz[1] + w / 4.0, xyz[2] - w / 4.0,
                    w / 2.0);
  c[3] = new OTNode(l + 1, xyz[0] - w / 4.0, xyz[1] + w / 4.0, xyz[2] + w / 4.0,
                    w / 2.0);
  c[4] = new OTNode(l + 1, xyz[0] + w / 4.0, xyz[1] - w / 4.0, xyz[2] - w / 4.0,
                    w / 2.0);
  c[5] = new OTNode(l + 1, xyz[0] + w / 4.0, xyz[1] - w / 4.0, xyz[2] + w / 4.0,
                    w / 2.0);
  c[6] = new OTNode(l + 1, xyz[0] + w / 4.0, xyz[1] + w / 4.0, xyz[2] - w / 4.0,
                    w / 2.0);
  c[7] = new OTNode(l + 1, xyz[0] + w / 4.0, xyz[1] + w / 4.0, xyz[2] + w / 4.0,
                    w / 2.0);

  for (int iA = 0; iA < atoms.size(); iA++) {
    for (int iC = 0; iC < 8; iC++) {
      if (tree->DoesOverlap(c[iC], atoms[iA])) {
        c[iC]->atoms.push_back(atoms[iA]);
        c[iC]->c = NULL;
      }
      c[iC]->parent = this; // needed during merge
      //      pg->addPoint(c[iC]);    //needed to quickly find affected cells
      //      during update
    }
  }

  // clear the atoms list for the master cell.. I don't think I need it any more
  // atoms.clear();
}

void OTNode::Merge(std::vector<Atom *> &allAtoms) {

  if (c == NULL) {
    cout << "WARNING: Merge() called on a leaf." << endl;
    return;
  }
  for (int i = 1; i < 8; i++)
    if (c[i]->c != NULL)
      return;

  atoms.clear();

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < c[i]->atoms.size(); j++) {
      int index = c[i]->atoms[j];

      if (!allAtoms[index]->copyDisabled) {
        atoms.push_back(index);
        allAtoms[index]->copyDisabled = true;
      }
    }
    c[i]->atoms.clear();
    delete c[i];
  }
  delete c;

  for (int i = 0; i < atoms.size(); i++)
    allAtoms[atoms[i]]->copyDisabled = false;
}

void OTNode::getInterfaceArea(OTNode *otherNode, double *trans,
                              double interfaceWidth, double &areaInner,
                              double &areaOuter, double isoValue) {
  double sumRad =
      (w + otherNode->w) * 0.866; // sum of radii of circumscribed spheres
  if (l > 3)
    sumRad += interfaceWidth;
  double sumRad2 = sumRad * sumRad;

  double oldPos[4] = {otherNode->xyz[0], otherNode->xyz[1], otherNode->xyz[2],
                      1.0};
  double newPos[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      newPos[i] += oldPos[j] * trans[i * 4 + j];

  double dx = xyz[0] - newPos[0];
  double dy = xyz[1] - newPos[1];
  double dz = xyz[2] - newPos[2];
  double d2 = dx * dx + dy * dy + dz * dz;

  if (d2 < sumRad2) // close enough
  {
    bool thisIsALeaf = (c == NULL);
    bool otherIsALeaf = (otherNode->c == NULL);

    if (thisIsALeaf && otherIsALeaf) // base case
    {
      if (flag == 1 && otherNode->flag == 1) // both contains boundary
      //                      if(l == 6 && otherNode->l == 6) // both contains
      //                      boundary
      {
        double interfaceWidthSq = interfaceWidth * interfaceWidth;

        if (interfaceWidthSq > d2) {
          if (v > isoValue) // rough approx
            areaInner += w * w;
          else
            areaOuter += w * w;

          if (interfaceFlag ==
              false) // marking for contouring based area computation
          {
            interfaceFlag = true;
          }
        }
      }
    } else if (thisIsALeaf) // break down other
    {
      for (int j = 0; j < 8; j++)
        getInterfaceArea(otherNode->c[j], trans, interfaceWidth, areaInner,
                         areaOuter, isoValue);
    } else if (otherIsALeaf) // break down this
    {
      for (int i = 0; i < 8; i++)
        c[i]->getInterfaceArea(otherNode, trans, interfaceWidth, areaInner,
                               areaOuter, isoValue);
    } else // break down both
    {
      for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
          c[i]->getInterfaceArea(otherNode->c[j], trans, interfaceWidth,
                                 areaInner, areaOuter, isoValue);
    }
  }

  return;
}

void OTNode::Print() {
  cout << "Node: (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ") , " << w
       << endl;
}

void Cell::Print() {
  cout << "Cell: "; //("<< xyz[0] << "," << xyz[1] << "," << xyz[2] << ")";
  for (int i = 0; i < 8; i++) {
    cout << ", " << values[i];
  }
  cout << endl;
}

void Atom::Print() {
  cout << "Atom: (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")";
  cout << " " << r << endl;
}

// Octree::Octree(int n, double w, double x, double y, double z) {
Octree::Octree(int n, string fname) {
  paircount = 0;
  triplecount = 0;
  tripleprob = 0;
  pairprob = 0;
  numoverwrite = 0;
  numSDFCalled = 0;
  numPairFound = 0;
  numOldPairSearched = 0;
  numOldPairUsed = 0;
  numPairComputed = 0;
  numTripleFound = 0;
  numOldTripleSearched = 0;
  numOldTripleUsed = 0;
  numTripleComputed = 0;

  ReadPQR(fname);

  basesize = n;

  double scale = 1.5;
  size *= scale;
  size += 12;

  basewidth = size / n;

  t = new OTNode[n * n * n];

  t[0].tree = this;

  for (int i = 0; i < 3; i++) {
    base[i] = center[i] - size / 2.0;
  }

  cout << "Basewidth: " << basewidth << endl;
  cout << "Base: " << base[0] << " " << base[1] << " " << base[2] << endl;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        OTNode *nd = GetNode(i, j, k);
        nd->l = 0;
        nd->w = basewidth;
        nd->xyz[0] = base[0] + basewidth / 2.0 + i * basewidth;
        nd->xyz[1] = base[1] + basewidth / 2.0 + j * basewidth;
        nd->xyz[2] = base[2] + basewidth / 2.0 + k * basewidth;
        nd->c = NULL;
        nd->splitDisabled = false;
        //      pg->addPoint(nd);
      }
    }
  }
}

Octree::Octree(int n, string fname, string fname2) {
  paircount = 0;
  triplecount = 0;
  tripleprob = 0;
  pairprob = 0;
  numoverwrite = 0;
  numSDFCalled = 0;
  numPairFound = 0;
  numOldPairSearched = 0;
  numOldPairUsed = 0;
  numPairComputed = 0;
  numTripleFound = 0;
  numOldTripleSearched = 0;
  numOldTripleUsed = 0;
  numTripleComputed = 0;
  ReadPQR(fname, fname2);

  basesize = n;

  double scale = 1.5;
  size *= scale;
  size += 12;

  basewidth = size / n;

  t = new OTNode[n * n * n];

  t[0].tree = this;

  for (int i = 0; i < 3; i++) {
    base[i] = center[i] - size / 2.0;
  }

  cout << "Basewidth: " << basewidth << endl;
  cout << "Base: " << base[0] << " " << base[1] << " " << base[2] << endl;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        OTNode *nd = GetNode(i, j, k);
        nd->l = 0;
        nd->w = basewidth;
        nd->xyz[0] = base[0] + basewidth / 2.0 + i * basewidth;
        nd->xyz[1] = base[1] + basewidth / 2.0 + j * basewidth;
        nd->xyz[2] = base[2] + basewidth / 2.0 + k * basewidth;
        nd->c = NULL;
        nd->splitDisabled = false;
        //      pg->addPoint(nd);
      }
    }
  }
}

OTNode *Octree::BaseLocate(double *xyz) {
  int i = (xyz[0] - base[0]) / basewidth;
  int j = (xyz[1] - base[1]) / basewidth;
  int k = (xyz[2] - base[2]) / basewidth;

  if (i < 0 || j < 0 || k < 0 || i >= basesize || j >= basesize ||
      k >= basesize) {
    return NULL;
  }

  return GetNode(i, j, k);
}

OTNode *Octree::Locate(double *xyz) { return Locate(xyz, 1000); }

OTNode *Octree::Locate(double *xyz, int depth) {
  OTNode *current = BaseLocate(xyz);

  if (current == NULL) {
    //              cout<<"base locate returned null"<<endl;
    //              cout<<"["<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<"]"<<endl;
    return NULL;
  }

  int curD = 0;
  bool done = false;
  while (!done && curD < depth) {

    curD++;
    // current->Print();
    if (current->c == NULL) {
      done = true;
    } else {

      if (xyz[0] < current->xyz[0] && xyz[1] < current->xyz[1] &&
          xyz[2] < current->xyz[2]) {
        current = current->c[0];
      } else if (xyz[0] < current->xyz[0] && xyz[1] < current->xyz[1] &&
                 xyz[2] >= current->xyz[2]) {
        current = current->c[1];
      } else if (xyz[0] < current->xyz[0] && xyz[1] >= current->xyz[1] &&
                 xyz[2] < current->xyz[2]) {
        current = current->c[2];
      } else if (xyz[0] < current->xyz[0] && xyz[1] >= current->xyz[1] &&
                 xyz[2] >= current->xyz[2]) {
        current = current->c[3];
      } else if (xyz[0] >= current->xyz[0] && xyz[1] < current->xyz[1] &&
                 xyz[2] < current->xyz[2]) {
        current = current->c[4];
      } else if (xyz[0] >= current->xyz[0] && xyz[1] < current->xyz[1] &&
                 xyz[2] >= current->xyz[2]) {
        current = current->c[5];
      } else if (xyz[0] >= current->xyz[0] && xyz[1] >= current->xyz[1] &&
                 xyz[2] < current->xyz[2]) {
        current = current->c[6];
      } else if (xyz[0] >= current->xyz[0] && xyz[1] >= current->xyz[1] &&
                 xyz[2] >= current->xyz[2]) {
        current = current->c[7];
      }
    }
  }
  return current;
}

void Octree::SetInitialData() {

  // loop over the atoms
  for (int iA = 0; iA < atoms.size(); iA++) {
    // put them in the appropriate lists...
    // FIXME: this is currently inefficient
    for (int i = 0; i < basesize; i++) {
      for (int j = 0; j < basesize; j++) {
        for (int k = 0; k < basesize; k++) {
          OTNode *nd = GetNode(i, j, k);
          // FIXME: this is just a VDW hack...
          if (DoesOverlap(nd, iA)) {
            nd->atoms.push_back(iA);
          }
        }
      }
    }
  }

  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
#ifdef VDW
        computeSDF(nd, 1.4);
        nd->sdfComputed = true;
        // nd->SetData(fvdw(nd));
#endif
#ifdef GAUSS
        nd->SetData(fgauss(nd));
#endif
        newNodes.push_back(nd);
        // cout << fvdw(nd) << ", ";
        // nd->Print();

      } // cout << endl;
    } // cout << endl;
  }
}

void Octree::Neighbors(OTNode *nd, OTNode **nbs) {
  double width = nd->w;
  int depth = nd->l;
  double xyz[3];

  xyz[0] = nd->xyz[0] - width;
  xyz[1] = nd->xyz[1];
  xyz[2] = nd->xyz[2];
  nbs[0] = Locate(xyz, depth);

  xyz[0] = nd->xyz[0] + width;
  xyz[1] = nd->xyz[1];
  xyz[2] = nd->xyz[2];
  nbs[1] = Locate(xyz, depth);

  xyz[0] = nd->xyz[0];
  xyz[1] = nd->xyz[1] - width;
  xyz[2] = nd->xyz[2];
  nbs[2] = Locate(xyz, depth);

  xyz[0] = nd->xyz[0];
  xyz[1] = nd->xyz[1] + width;
  xyz[2] = nd->xyz[2];
  nbs[3] = Locate(xyz, depth);

  xyz[0] = nd->xyz[0];
  xyz[1] = nd->xyz[1];
  xyz[2] = nd->xyz[2] - width;
  nbs[4] = Locate(xyz, depth);

  xyz[0] = nd->xyz[0];
  xyz[1] = nd->xyz[1];
  xyz[2] = nd->xyz[2] + width;
  nbs[5] = Locate(xyz, depth);
}

void Octree::MarkExtendedNeighbors(OTNode *nd) {
  double width = nd->w;
  int depth = nd->l;
  double xyz[3];

  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        xyz[0] = nd->xyz[0] + i * width;
        xyz[1] = nd->xyz[1] + j * width;
        xyz[2] = nd->xyz[2] + k * width;
        OTNode *nd1 = Locate(xyz, depth);
        if (nd1 != NULL)
          nd1->flag = 1;
      }
    }
  }
}

void Octree::MarkExtendedNeighbors2(OTNode *nd) {
  double width = nd->w;
  int depth = nd->l;
  double xyz[3];

  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        xyz[0] = nd->xyz[0] + i * width;
        xyz[1] = nd->xyz[1] + j * width;
        xyz[2] = nd->xyz[2] + k * width;
        OTNode *nd1 = Locate(xyz, depth);
        if (nd1 != NULL) {
          if (nd1->flag == 1 && nd->flag == 1)
            continue;
          else if ((nd1->v - contourValue) * (nd->v - contourValue) < 0.0) {
            if (nd1->flag != 1) {
              nd1->flag = 1;
              markedNodes.push_back(nd1);
            }
            if (nd->flag != 1) {
              nd->flag = 1;
              markedNodes.push_back(nd);
            }
          }
        }
      }
    }
  }
}

void Octree::GetExtendedNeighbors(OTNode *nd, OTNode **nbs) {
  double width = nd->w;
  int depth = nd->l;
  double xyz[3];
  int count = 0;
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        xyz[0] = nd->xyz[0] + i * width;
        xyz[1] = nd->xyz[1] + j * width;
        xyz[2] = nd->xyz[2] + k * width;
        nbs[count++] = Locate(xyz, depth);
      }
    }
  }
}

void Octree::GetAffectedCells(
    Atom *a,
    vector<OTNode *> &intersectedCells) // to be used instead of PG-range query
{
  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);

        if (DoesOverlap(nd, a)) {
          //                                      cout<<"found base"<<endl;
          LocateAffectedCells(nd, a, intersectedCells);
        }
      }
    }
  }
}

void Octree::LocateAffectedCells(OTNode *nd, Atom *a,
                                 vector<OTNode *> &intersectedCells) {
  if (nd == NULL) {
    cout << "caught a null" << endl;
    return;
  }

  if (!DoesOverlap(nd, a))
    return;
  else {
    if (nd->c == NULL) {
      intersectedCells.push_back(nd);
      return;
    } else {
      for (int i = 0; i < 8; i++) {
        if (nd->c[i] != NULL)
          LocateAffectedCells(nd->c[i], a, intersectedCells);
      }
    }
  }
}

void Octree::AddAtom(Atom *a) {
  vector<OTNode *> intersectedCells;
  //      intersectedCells = pg->range(new Point(a->xyz[0],a->xyz[1],a->xyz[2]),
  //      a->r);

  GetAffectedCells(a, intersectedCells);

  int numCells = intersectedCells.size();

  vector<OTNode *> cellsToBeUpdated;

  for (int i = 0; i < numCells; i++) {
    OTNode *nd = (OTNode *)intersectedCells[i];
    //              nd->atoms.push_back(a);

    double prevVal = nd->v;
    fgauss_update_add(nd, a);
    double newVal = nd->v;

    if ((prevVal > contourValue && newVal < contourValue) ||
        (prevVal > contourValue && newVal < contourValue)) {
      nd->markedForLocalUpdate = true;
      cellsToBeUpdated.push_back(nd);
    }
  }

  LocalRefinement(cellsToBeUpdated);
}

void Octree::AddAtoms(vector<Atom *> atms) {
  vector<OTNode *> cellsToBeUpdated;

  int numAtoms = atms.size();

  for (int j = 0; j < numAtoms; j++) {
    vector<OTNode *> intersectedCells;
    //     intersectedCells = pg->range(new
    //     Point(a->xyz[0],a->xyz[1],a->xyz[2]), a->r);

    GetAffectedCells(atms[j], intersectedCells);

    int numCells = intersectedCells.size();

    for (int i = 0; i < numCells; i++) {
      OTNode *nd = (OTNode *)intersectedCells[i];
      //             nd->atoms.push_back(a);

      double prevVal = nd->v;
      fgauss_update_add(nd, atms[j]);
      double newVal = nd->v;

      if ((prevVal > contourValue && newVal < contourValue) ||
          (prevVal > contourValue && newVal < contourValue)) {
        if (nd->markedForLocalUpdate == false) {
          cellsToBeUpdated.push_back(nd);
        }
        nd->markedForLocalUpdate = true;
      }
    }
  }

  LocalRefinement(cellsToBeUpdated);
}

void Octree::RemoveAtom(Atom *a) {
  vector<OTNode *> intersectedCells;
  //      intersectedCells = pg->range(new Point(a->xyz[0],a->xyz[1],a->xyz[2]),
  //      a->r);

  GetAffectedCells(a, intersectedCells);

  int numCells = intersectedCells.size();

  //      cout<<"NumCells: "<<numCells<<endl;

  vector<OTNode *> cellsToBeUpdated;

  for (int i = 0; i < numCells; i++) {
    OTNode *nd = (OTNode *)intersectedCells[i];

    double prevVal = nd->v;
    fgauss_update_remove(nd, a);
    double newVal = nd->v;

    if ((prevVal > contourValue && newVal < contourValue) ||
        (prevVal > contourValue && newVal < contourValue)) {
      cellsToBeUpdated.push_back(nd);
      nd->markedForLocalUpdate = true;
    }
  }

  //      cout<<"marked cells for refinement"<<endl;

  LocalRefinement(cellsToBeUpdated);
}

void Octree::RemoveAtoms(vector<Atom *> atms) {
  vector<OTNode *> cellsToBeUpdated;

  int numAtoms = atms.size();

  for (int j = 0; j < numAtoms; j++) {
    vector<OTNode *> intersectedCells;
    //     intersectedCells = pg->range(new
    //     Point(a->xyz[0],a->xyz[1],a->xyz[2]), a->r);

    GetAffectedCells(atms[j], intersectedCells);

    int numCells = intersectedCells.size();

    //     cout<<"NumCells: "<<numCells<<endl;

    for (int i = 0; i < numCells; i++) {
      OTNode *nd = (OTNode *)intersectedCells[i];

      double prevVal = nd->v;
      fgauss_update_remove(nd, atms[j]);
      double newVal = nd->v;

      if ((prevVal > contourValue && newVal < contourValue) ||
          (prevVal > contourValue && newVal < contourValue)) {
        if (nd->markedForLocalUpdate == false) {
          cellsToBeUpdated.push_back(nd);
        }
        nd->markedForLocalUpdate = true;
      }
    }
  }

  //      cout<<"marked cells for refinement"<<endl;

  LocalRefinement(cellsToBeUpdated);
}

void Octree::LocalRefinement(vector<OTNode *> cellsToBeUpdated) {
  vector<OTNode *> cellsToBeSplit;
  vector<OTNode *> cellsToBeMerged;

  int iteration = 0;

  // Keep splitting cells if necessary. Mark the cells to be merged.
  // If a cell is marked for merge, it cannot be considered for split.
  // So all the merge are only one level deep. No iteration needed.
  // All merge are done at the end.

  while (1) {
    int numCellsToBeUpdated = cellsToBeUpdated.size();
    for (int i = 0; i < numCellsToBeUpdated; i++) {
      OTNode *nd = (OTNode *)cellsToBeUpdated[i];

      if (mergeRequired(nd)) {
        cellsToBeMerged.push_back(nd);
        nd->splitDisabled = true;
      }

      else if (splitRequired(nd)) {
        OTNode *nbs[27];
        GetExtendedNeighbors(nd, nbs);
        for (int j = 0; j < 27; j++) {
          nbs[j]->flag = 1;
          cellsToBeSplit.push_back(nbs[j]);
        }
      }
    }

    //              cout<<"Iteration "<<iteration<<endl;
    //              cout<<"Marked for split "<<cellsToBeSplit.size()<<endl;
    //              cout<<"Marked for merge "<<cellsToBeMerged.size()<<endl;

    int numCellsToBeSplit = cellsToBeSplit.size();
    for (int i = 0; i < numCellsToBeSplit; i++)
      SplitMarkedCells(cellsToBeSplit[i]);
    cellsToBeSplit.clear();

    int newCells =
        uncomputedCells
            .size(); // These are newly created cells due to split.
                     // We compute the values in these cells and consider these
                     // as our new list of cellsToBeUpdated. If this list is
                     // empty, then we break
    if (newCells > 0) {
      cellsToBeUpdated.clear();
      for (int i = 0; i < newCells; i++) {
        cellsToBeUpdated.push_back(uncomputedCells[i]);
      }

      ComputeNewCellValues();
    } else
      break;
  }

  int numCellsToBeMerged = cellsToBeMerged.size();
  for (int i = 0; i < numCellsToBeMerged; i++) {
    cellsToBeMerged[i]->Merge(atoms);
    cellsToBeMerged[i]->splitDisabled = false;
  }
  cellsToBeMerged.clear();
}

bool Octree::mergeRequired(OTNode *nd) {
  if (nd->c == NULL) // a leaf cannot be merged
    return false;

  for (int j = 0; j < 6; j++) {
    OTNode *child = nd->c[j];
    if (child != NULL)
      if (splitRequired(child))
        return false;
  }

  return true;
}

bool Octree::splitRequired(OTNode *nd) {
  if (nd->c != NULL) // already split
  {
    return false;
  }

  if (nd->splitDisabled) // cell is currently marked for merge.
  {
    return false;
  }

  OTNode *nbs[6];
  Neighbors(nd, nbs);
  if ((nd->v > contourValue &&
       ((nbs[0] != NULL && nbs[0]->v < contourValue) ||
        (nbs[1] != NULL && nbs[1]->v < contourValue) ||
        (nbs[2] != NULL && nbs[2]->v < contourValue) ||
        (nbs[3] != NULL && nbs[3]->v < contourValue) ||
        (nbs[4] != NULL && nbs[4]->v < contourValue) ||
        (nbs[5] != NULL && nbs[5]->v < contourValue))) ||
      (nd->v < contourValue &&
       ((nbs[0] != NULL && nbs[0]->v > contourValue) ||
        (nbs[1] != NULL && nbs[1]->v > contourValue) ||
        (nbs[2] != NULL && nbs[2]->v > contourValue) ||
        (nbs[3] != NULL && nbs[3]->v > contourValue) ||
        (nbs[4] != NULL && nbs[4]->v > contourValue) ||
        (nbs[5] != NULL && nbs[5]->v > contourValue)))) {

    return true;
  }
  return false;
}

void Octree::MarkIsoContour(OTNode *nd) {
  if (nd->c == NULL) {
    if (splitRequired(nd))
      MarkExtendedNeighbors(nd);
    // MarkExtendedNeighbors2(nd);
  } else {
    for (int i = 0; i < 8; i++) {
      MarkIsoContour(nd->c[i]);
    }
  }
}

void Octree::MarkIsoContour() {

  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        MarkIsoContour(nd);
      }
    }
  }

  /*
    markedNodes.clear();
    int size = newNodes.size();
    for(int i=0;i<size;i++)
          MarkIsoContour(newNodes[i]);
    newNodes.clear();
  */
}

void Octree::SplitMarkedCells(OTNode *nd) {
  if (nd->c == NULL) {
    if (nd->flag != -1) {
      //      if (nd->l < targetLevel || nd->w > targetRes) {
      if (nd->l < targetLevel) {
        nd->Split();
        for (int i = 0; i < 8; i++) {
          uncomputedCells.push_back(nd->c[i]);
        }
      }
    }
  } else {
    for (int i = 0; i < 8; i++) {
      SplitMarkedCells(nd->c[i]);
    }
  }
  nd->flag = -1;
}

void Octree::SplitMarkedCells() {

  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        SplitMarkedCells(nd);
      }
    }
  }

  /*
    int size = markedNodes.size();
    for(int i=0;i<size;i++)
          SplitMarkedCells(markedNodes[i]);
    markedNodes.clear();
   */
}

bool Octree::ComputeNewCellValues() {
  int size = uncomputedCells.size();

  cout << "Computing " << size << " new values" << endl;

  for (int i = 0; i < size; i++) {
    OTNode *nd = uncomputedCells[i];

#ifdef VDW
    if (nd->l > (targetLevel - 3)) {
      computeSDF(nd, 1.4);
      nd->sdfComputed = true;
    } else
      nd->SetData(fvdw(nd));
#endif
#ifdef GAUSS
    nd->SetData(fgauss(nd));
#endif
    newNodes.push_back(nd);
  }

  uncomputedCells.clear();

  if (size > 0) {
    return true;
  }
  return false;
}

void Octree::GetFinalCell(OTNode *nd, OTNode **nbs) {

  double width = nd->w;
  int depth = nd->l;
  double xyz[3];

  int count = 0;

  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++) {
      for (int k = 0; k <= 1; k++) {
        xyz[0] = nd->xyz[0] + i * width;
        xyz[1] = nd->xyz[1] + j * width;
        xyz[2] = nd->xyz[2] + k * width;

        OTNode *nd1 = Locate(xyz, depth);

        nbs[count] = nd1;
        count++;
      }
    }
  }
}

void Octree::GetFinalCellsAndDestroyTree(OTNode *nd) {
  //  if (nd->l == targetLevel || nd->w < targetRes) {
  if (nd->l == targetLevel) {
    // check and see if this corresponds to a valid cell...

    // get the 8 neighbors to make a cell...
    OTNode *nbs[8];
    GetFinalCell(nd, nbs);

    bool sameSize = true;
    bool containsContour = false;

    for (int i = 1; i < 8; i++) {
      if ((nbs[0]->v < contourValue && nbs[i]->v > contourValue) ||
          (nbs[0]->v > contourValue && nbs[i]->v < contourValue)) {
        containsContour = true;
      }
      if (nbs[i]->l != nbs[0]->l) {
        sameSize = false;
      }
    }

    if (sameSize && containsContour) {

      /* // old style
      tmp.xyz[0] = nd->xyz[0];
      tmp.xyz[1] = nd->xyz[1];
      tmp.xyz[2] = nd->xyz[2];
      */

      Cell tmp;
      for (int i = 0; i < 8; i++) {
        if (nbs[i]->flag < 0) {
          cout << "ERROR: flag is " << nbs[i]->flag << endl;
        }
      }

      // arand: neighbors function and the marching cubes codes use
      //        different orderings...

      tmp.vertices[0] = nbs[0]->flag;
      tmp.vertices[4] = nbs[1]->flag;
      tmp.vertices[2] = nbs[2]->flag;
      tmp.vertices[6] = nbs[3]->flag;
      tmp.vertices[1] = nbs[4]->flag;
      tmp.vertices[5] = nbs[5]->flag;
      tmp.vertices[3] = nbs[6]->flag;
      tmp.vertices[7] = nbs[7]->flag;

      // gauss or sdf from SAS

      tmp.values[0] = nbs[0]->v;
      tmp.values[4] = nbs[1]->v;
      tmp.values[2] = nbs[2]->v;
      tmp.values[6] = nbs[3]->v;
      tmp.values[1] = nbs[4]->v;
      tmp.values[5] = nbs[5]->v;
      tmp.values[3] = nbs[6]->v;
      tmp.values[7] = nbs[7]->v;

      // sdf from SAS
      /*
            tmp.values[0] = nbs[0]->sdf;
            tmp.values[4] = nbs[1]->sdf;
            tmp.values[2] = nbs[2]->sdf;
            tmp.values[6] = nbs[3]->sdf;
            tmp.values[1] = nbs[4]->sdf;
            tmp.values[5] = nbs[5]->sdf;
            tmp.values[3] = nbs[6]->sdf;
            tmp.values[7] = nbs[7]->sdf;
      */

      //      cout<<"Cell values ["<< tmp.values[0] << ", "<< tmp.values[1] <<
      //      ", "<< tmp.values[2] << ", "<< tmp.values[3] << ", "<<
      //      tmp.values[4] << ", "<< tmp.values[5] << ", "<< tmp.values[6] <<
      //      ", "<< tmp.values[7] << "]"<< endl;

      // for the regular marching cubes..
      /*
      tmp.vertices[0] = nbs[0]->flag;
      tmp.vertices[4] = nbs[1]->flag;
      tmp.vertices[3] = nbs[2]->flag;
      tmp.vertices[7] = nbs[3]->flag;
      tmp.vertices[1] = nbs[4]->flag;
      tmp.vertices[5] = nbs[5]->flag;
      tmp.vertices[2] = nbs[6]->flag;
      tmp.vertices[6] = nbs[7]->flag;

      tmp.values[0] = nbs[0]->v;
      tmp.values[4] = nbs[1]->v;
      tmp.values[3] = nbs[2]->v;
      tmp.values[7] = nbs[3]->v;
      tmp.values[1] = nbs[4]->v;
      tmp.values[5] = nbs[5]->v;
      tmp.values[2] = nbs[6]->v;
      tmp.values[6] = nbs[7]->v;
      */
      /*
      tmp.vertices[0] = nbs[0]->flag;
      tmp.vertices[1] = nbs[1]->flag;
      tmp.vertices[2] = nbs[2]->flag;
      tmp.vertices[3] = nbs[3]->flag;
      tmp.vertices[4] = nbs[4]->flag;
      tmp.vertices[5] = nbs[5]->flag;
      tmp.vertices[6] = nbs[6]->flag;
      tmp.vertices[7] = nbs[7]->flag;

      tmp.values[0] = nbs[0]->v;
      tmp.values[1] = nbs[1]->v;
      tmp.values[2] = nbs[2]->v;
      tmp.values[3] = nbs[3]->v;
      tmp.values[4] = nbs[4]->v;
      tmp.values[5] = nbs[5]->v;
      tmp.values[6] = nbs[6]->v;
      tmp.values[7] = nbs[7]->v;
      */

      tmp.span = nd->w;
      nd->dual = &tmp;
      finalCells.push_back(tmp);
      // tmp.Print();
    }
  }

  if (nd->c != NULL) {
    // recurse
    for (int i = 0; i < 8; i++) {
      GetFinalCellsAndDestroyTree(nd->c[i]);
    }

    // delete children
    /*    for (int i=0; i<8; i++) {
          delete nd->c[i];
        }*/
  }
}

void Octree::GetFinalCellsAndDestroyTree() {
  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        GetFinalCellsAndDestroyTree(nd);
      }
    }
  }
  cout << "Number of cells for isocontouring: " << finalCells.size() << endl;
}

void Octree::ComputeFinalValues() {
  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        ComputeFinalValues(nd);
      }
    }
  }
  cout << "Number of final vertices: " << finalVertices.size() << endl;
}

void Octree::ComputeFinalValues(OTNode *nd) {
  if (nd->l == targetLevel) {
    if (nd->sdfComputed == false) {
      computeSDF(nd, 1.4);
      nd->sdfComputed = true;
    }
  }
  if (nd->c != NULL) {
    // recurse
    for (int i = 0; i < 8; i++) {
      ComputeFinalValues(nd->c[i]);
    }
  }
}

void Octree::LabelFinalVertices() {
  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        LabelFinalVertices(nd);
      }
    }
  }
  cout << "Number of final vertices: " << finalVertices.size() << endl;
}

void Octree::LabelFinalVertices(OTNode *nd) {
  //  if (nd->l == targetLevel || nd->w < targetRes) {

  if (nd->l == targetLevel) {
    // check and see if this corresponds to a valid cell...

    // get the 8 neighbors to make a cell...
    OTNode *nbs[8];
    GetFinalCell(nd, nbs);

    bool sameSize = true;
    bool containsContour = true;

    for (int i = 1; i < 8; i++) {
      if ((nbs[0]->v < contourValue && nbs[i]->v > contourValue) ||
          (nbs[0]->v > contourValue && nbs[i]->v < contourValue)) {
        containsContour = true;
      }
      if (nbs[i]->l != nbs[0]->l) {
        sameSize = false;
      }
    }

    if (sameSize && containsContour) {
      for (int i = 0; i < 8; i++) {
        if (nbs[i]->flag < 0) {
          nbs[i]->flag = finalVertices.size();
          Vertex tmp;
          tmp.xyz[0] = nbs[i]->xyz[0];
          tmp.xyz[1] = nbs[i]->xyz[1];
          tmp.xyz[2] = nbs[i]->xyz[2];

#ifdef VDW
          // sdf computation.
          if (nbs[i]->sdfComputed == false) {
            computeSDF(nbs[i], 1.4);
            nbs[i]->sdfComputed = true;
          }

          // nbs[i]->SetData(fvdw(nbs[i]));
#endif
#ifdef GAUSS
          nbs[i]->SetData(fgauss(nbs[i]));
#endif
          finalVertices.push_back(tmp);
        }
      }
    }
  }

  if (nd->c != NULL) {
    // recurse
    for (int i = 0; i < 8; i++) {
      LabelFinalVertices(nd->c[i]);
    }
  }
}

double getPointToCircleDistance(
    Atom *a1, Atom *a2, double x, double y, double z,
    double solventRad) // circle specified by two intersecting spheres
{
  double cx1 = a1->xyz[0];
  double cy1 = a1->xyz[1];
  double cz1 = a1->xyz[2];
  double R1 = a1->r + solventRad;

  double cx2 = a2->xyz[0];
  double cy2 = a2->xyz[1];
  double cz2 = a2->xyz[2];
  double R2 = a2->r + solventRad;

  double dx = cx1 - cx2;
  double dy = cy1 - cy2;
  double dz = cz1 - cz2;

  double dx1 = x - cx1;
  double dy1 = y - cy1;
  double dz1 = z - cz1;

  double dx2 = x - cx2;
  double dy2 = y - cy2;
  double dz2 = z - cz2;

  double distSq = (dx * dx + dy * dy + dz * dz);
  double dist = sqrt(distSq);

  double R1Sq = R1 * R1;
  double R2Sq = R2 * R2;

  // equation of plane containing the circle
  double A = dx / dist;
  double B = dy / dist;
  double C = dz / dist;
  double D = -((cx1 * cx1 - cx2 * cx2) + (cy1 * cy1 - cy2 * cy2) +
               (cz1 * cz1 - cz2 * cz2) - (R1Sq - R2Sq)) /
             (2 * dist);

  // center and radius of circle
  double distanceFromCenter1ToCircle = (R1Sq - R2Sq + distSq) / (2 * dist);
  double distanceFromCenter2ToCircle = dist - distanceFromCenter1ToCircle;
  double circleCenterX =
      (cx1 * distanceFromCenter2ToCircle + cx2 * distanceFromCenter1ToCircle) /
      dist;
  double circleCenterY =
      (cy1 * distanceFromCenter2ToCircle + cy2 * distanceFromCenter1ToCircle) /
      dist;
  double circleCenterZ =
      (cz1 * distanceFromCenter2ToCircle + cz2 * distanceFromCenter1ToCircle) /
      dist;
  double circleRad =
      sqrt(R1Sq - distanceFromCenter1ToCircle * distanceFromCenter1ToCircle);

  // projecting the point onto the plane
  double pointPlaneDist = A * x + B * y + C * z + D;
  double projectedX = x - A * pointPlaneDist;
  double projectedY = y - B * pointPlaneDist;
  double projectedZ = z - C * pointPlaneDist;

  // distance from the projected point to the center along the plane
  double onPlaneDirFromCCX = projectedX - circleCenterX;
  double onPlaneDirFromCCY = projectedY - circleCenterY;
  double onPlaneDirFromCCZ = projectedZ - circleCenterZ;
  double onPlaneDistanceToCircleCenter =
      sqrt(onPlaneDirFromCCX * onPlaneDirFromCCX +
           onPlaneDirFromCCY * onPlaneDirFromCCY +
           onPlaneDirFromCCZ * onPlaneDirFromCCZ);

  // extrapolation/interpolation to get the closest point on the circle from the
  // projected point
  double ratio = circleRad / onPlaneDistanceToCircleCenter;
  double closestPointOnCircleX = circleCenterX + ratio * onPlaneDirFromCCX;
  double closestPointOnCircleY = circleCenterY + ratio * onPlaneDirFromCCY;
  double closestPointOnCircleZ = circleCenterZ + ratio * onPlaneDirFromCCZ;

  // distance from the point on the circle to the given point
  double dPtX = x - closestPointOnCircleX;
  double dPtY = y - closestPointOnCircleY;
  double dPtZ = z - closestPointOnCircleZ;

  return sqrt(dPtX * dPtX + dPtY * dPtY + dPtZ * dPtZ);
}

bool isPointInsideTet(
    double P1[], double P2[], double P3[], double P4[],
    double P[]) // returns true if P is inside the tet formed by P1,P2,P3,P4
{
  double vecP1P2[3] = {P2[0] - P1[0], P2[1] - P1[1], P2[2] - P1[2]};
  double vecP1P3[3] = {P3[0] - P1[0], P3[1] - P1[1], P3[2] - P1[2]};
  double vecP1P4[3] = {P4[0] - P1[0], P4[1] - P1[1], P4[2] - P1[2]};
  double vecP2P4[3] = {P4[0] - P2[0], P4[1] - P2[1], P4[2] - P2[2]};
  double vecP2P3[3] = {P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2]};

  // plane1 P1 P2 P3. normal is n1 = vecP1P2 X vecP1P3. plane equation: n1.x -
  // n1.P1 = 0. sign of a point P: n1.P - n1.P1
  double n1[3];
  n1[0] = vecP1P2[1] * vecP1P3[2] - vecP1P2[2] * vecP1P3[1];
  n1[1] = vecP1P2[2] * vecP1P3[0] - vecP1P2[0] * vecP1P3[2];
  n1[2] = vecP1P2[0] * vecP1P3[1] - vecP1P2[1] * vecP1P3[0];

  double n1P1 =
      n1[0] * P1[0] + n1[1] * P1[1] + n1[2] * P1[2]; // a point on the plane
  double n1P4 =
      n1[0] * P4[0] + n1[1] * P4[1] + n1[2] * P4[2]; // other corner of tet
  double n1P = n1[0] * P[0] + n1[1] * P[1] + n1[2] * P[2]; // given point
  if (!(((n1P4 >= n1P1) && (n1P >= n1P1)) ||
        ((n1P4 <= n1P1) && (n1P <= n1P1)))) // not same sign
    return false;

  // plane2 P1 P2 P4. normal is n2 = vecP1P2 X vecP1P4
  double n2[3];
  n2[0] = vecP1P2[1] * vecP1P4[2] - vecP1P2[2] * vecP1P4[1];
  n2[1] = vecP1P2[2] * vecP1P4[0] - vecP1P2[0] * vecP1P4[2];
  n2[2] = vecP1P2[0] * vecP1P4[1] - vecP1P2[1] * vecP1P4[0];

  double n2P1 = n2[0] * P1[0] + n2[1] * P1[1] + n2[2] * P1[2];
  double n2P3 = n2[0] * P3[0] + n2[1] * P3[1] + n2[2] * P3[2];
  double n2P = n2[0] * P[0] + n2[1] * P[1] + n2[2] * P[2];
  if (!(((n2P3 >= n2P1) && (n2P >= n2P1)) ||
        ((n2P3 <= n2P1) && (n2P <= n2P1)))) // not same sign
    return false;

  // plane3 P1 P3 P4. normal is n3 = vecP1P3 X vecP1P4
  double n3[3];
  n3[0] = vecP1P3[1] * vecP1P4[2] - vecP1P3[2] * vecP1P4[1];
  n3[1] = vecP1P3[2] * vecP1P4[0] - vecP1P3[0] * vecP1P4[2];
  n3[2] = vecP1P3[0] * vecP1P4[1] - vecP1P3[1] * vecP1P4[0];

  double n3P1 = n3[0] * P1[0] + n3[1] * P1[1] + n3[2] * P1[2];
  double n3P2 = n3[0] * P2[0] + n3[1] * P2[1] + n3[2] * P2[2];
  double n3P = n3[0] * P[0] + n3[1] * P[1] + n3[2] * P[2];
  if (!(((n3P2 >= n3P1) && (n3P >= n3P1)) ||
        ((n3P2 <= n3P1) && (n3P <= n3P1)))) // not same sign
    return false;

  // plane4 P2 P3 P4. normal is n4 = vecP2P3 X vecP2P4
  double n4[3];
  n4[0] = vecP2P3[1] * vecP2P4[2] - vecP2P3[2] * vecP2P4[1];
  n4[1] = vecP2P3[2] * vecP2P4[0] - vecP2P3[0] * vecP2P4[2];
  n4[2] = vecP2P3[0] * vecP2P4[1] - vecP2P3[1] * vecP2P4[0];

  double n4P2 = n4[0] * P2[0] + n4[1] * P2[1] + n4[2] * P2[2];
  double n4P1 = n4[0] * P1[0] + n4[1] * P1[1] + n4[2] * P1[2];
  double n4P = n4[0] * P[0] + n4[1] * P[1] + n4[2] * P[2];
  if (!(((n4P1 >= n4P2) && (n4P >= n4P2)) ||
        ((n4P1 <= n4P2) && (n4P <= n4P2)))) // not same sign
    return false;

  return true;
}

bool Atom::isTripleIntersect(Atom *a1, Atom *a2, Atom *a3, double solventRad) {
  double distFromCenter3ToCircle12 = getPointToCircleDistance(
      a1, a2, a3->xyz[0], a3->xyz[1], a3->xyz[2], solventRad);
  return distFromCenter3ToCircle12 < (a3->r + solventRad);
}

bool Octree::isBuried(double *P) {
  OTNode *nd = Locate(P);

  if (nd == NULL) {
    int a;
    cin >> a;
    return false;
  }
  //      cout<<nd->atoms.size();

  return nd->isBuried(atoms, P, 1.4);
}

bool OTNode::isBuried(std::vector<Atom *> &atms, double *P, double solventRad) {
  int size = atoms.size();
  //      cout<<"checking if P is buried"<<endl;

  for (int i = 0; i < size; i++) {
    int index = atoms[i];
    //              cout<<index<<" ";
    Atom *a = atms[index];

    double dist = sqrt((a->xyz[0] - P[0]) * (a->xyz[0] - P[0]) +
                       (a->xyz[1] - P[1]) * (a->xyz[1] - P[1]) +
                       (a->xyz[2] - P[2]) * (a->xyz[2] - P[2]));
    //              cout<<dist<<endl;

    if (dist < a->r + solventRad - 0.001) {
      //                      cout<<"returning true"<<endl;
      return true;
    }
  }
  //      cout<<"returning false"<<endl;
  return false;
}

void Octree::computeSDF(OTNode *nd, double solventRad) {
  numSDFCalled++;
  int dum;
  int size = (nd->atoms).size();

  double sdf = 50.0;
  double dist;

  double singleSDF = 50.0;
  double pairSDF = 50.0;
  double tripleSDF = 50.0;

  bool insideVDW = false;
  bool inside = false;
  bool inCone = false;
  bool inTet = false;

  bool *insideAtom = new bool[size];
  bool **pairIntersect = new bool *[size];

  for (int i = 0; i < size; i++) {
    pairIntersect[i] = new bool[size];

    for (int j = 0; j < size; j++)
      pairIntersect[i][j] = false;

    pairIntersect[i][i] = true;
  }

  //      cout<<"gridpoint = ["<<nd->xyz[0]<<", "<<nd->xyz[1]<<",
  //      "<<nd->xyz[2]<<"]"<<endl; cout<<"number of atoms: "<<size<<endl;

  // Checking single atom cases. Finding closest distance to the SAS.

  for (int i = 0; i < size; i++) {
    Atom *a = atoms[(nd->atoms)[i]];

    //              cout<<a->xyz[0]<<" "<<a->xyz[1]<<" "<<a->xyz[2]<<"
    //              "<<a->r<<endl;

    double dx = nd->xyz[0] - a->xyz[0];
    double dy = nd->xyz[1] - a->xyz[1];
    double dz = nd->xyz[2] - a->xyz[2];

    double dist2 = sqrt(dx * dx + dy * dy + dz * dz);
    dist = a->r + solventRad - dist2;

    if (dist > 0) // ignore points outside SAS.
    {
      if (dist < singleSDF) {
        double ratio = (a->r + solventRad) / dist2;
        double closestPointOnSphere[3] = {a->xyz[0] + dx * ratio,
                                          a->xyz[1] + dy * ratio,
                                          a->xyz[2] + dz * ratio};

        //                              cout<<"closest point on sphere =
        //                              ["<<closestPointOnSphere[0]<<",
        //                              "<<closestPointOnSphere[1]<<",
        //                              "<<closestPointOnSphere[2]<<"].
        //                              Distance: "<<dist<<endl;

        if (!isBuried(closestPointOnSphere))
          singleSDF = dist;
      }
      inside = true;
      insideAtom[i] = true;

      if (dist > solventRad) // inside vdw
      {
        insideVDW = true;
      }
    }
  }
  //      cout << "Computed single atom SDF: "<<singleSDF<<endl;

  // Checking pairwise intersection cases. Finding closest distance to the SAS.

  if (inside) {
    for (int i = 0; i < size; i++) {
      for (int j = i + 1; j < size; j++) {
        int atomA = (nd->atoms)[i];
        int atomB = (nd->atoms)[j];

        bool pairSDFFound = false;

        //                              cout<<"checking for pairwise
        //                              intersection between "<<atomA<<" and
        //                              "<<atomB<<endl;

        if (atoms[atomA]->intersect(atoms[atomB], solventRad)) {
          numPairFound++;
          pairIntersect[i][j] = true;
          pairIntersect[j][i] = true;
          double tempsdf;
          bool previousPairFound = false;
          bool pairIntersectionReady = false;
          PairIntersection *PI = new PairIntersection();

          int numPairIntersections = atoms[atomA]->pairIntersections.size();
          //                                      cout<<"previous pairwise
          //                                      intersection count:
          //                                      "<<numPairIntersections<<endl;
          for (int npi = 0; npi < numPairIntersections; npi++) {
            numOldPairSearched++;
            PI = atoms[atomA]->pairIntersections[npi];

            if ((PI->atomA == atomA && PI->atomB == atomB) ||
                (PI->atomA == atomB &&
                 PI->atomB == atomA)) // already identified the intersection
            {
              //                                                      cout<<"previous
              //                                                      pairwise
              //                                                      intersection
              //                                                      between
              //                                                      "<<atomA<<"
              //                                                      and
              //                                                      "<<atomB<<"
              //                                                      found"<<endl;
              numOldPairUsed++;
              previousPairFound = true;
              pairIntersectionReady = true;
              break;
            }
          }

          if (!previousPairFound) {
            PI = new PairIntersection();
            if (getPairwiseIntersection(atomA, atomB, nd->xyz[0], nd->xyz[1],
                                        nd->xyz[2], solventRad, PI)) {
              numPairComputed++;
              //                                                      cout<<"new
              //                                                      pairwise
              //                                                      intersection
              //                                                      between
              //                                                      "<<atomA<<"
              //                                                      and
              //                                                      "<<atomB<<"
              //                                                      computed"<<endl;
              pairIntersectionReady = true;
            }
          }

          if (pairIntersectionReady) {
            if (getSDFForIntersectingSpheres(PI, nd->xyz[0], nd->xyz[1],
                                             nd->xyz[2], solventRad, tempsdf)) {
              //                                                      cout<<"new
              //                                                      pairwise
              //                                                      sdf
              //                                                      between
              //                                                      "<<atomA<<"
              //                                                      and
              //                                                      "<<atomB<<"
              //                                                      computed"<<endl;
              if (tempsdf > atoms[atomA]->r + solventRad &&
                  tempsdf > atoms[atomB]->r + solventRad) {
                cout << "In cone, but distance is too large " << tempsdf
                     << endl;
                pairprob++;
              }

              else if (tempsdf < pairSDF) {
                pairSDF = tempsdf;
                inCone = true;
              }

              paircount++;
            }
          }

        } // end if double intersect

        //                              cout<<"Finished checking for pairwise
        //                              intersection between "<<atomA<<" and
        //                              "<<atomB<<endl;

      } // end for j

    } // end for i

  } // end if inside

  //      cout<<"Pairwise SDF computed: "<<pairSDF<<endl;

  // Checking triple intersection cases. Finding closest distance to the SAS.

  if (inside) {
    for (int i = 0; i < size; i++) {
      for (int j = i + 1; j < size; j++) {
        for (int k = j + 1; k < size; k++) {
          int atomA = (nd->atoms)[i];
          int atomB = (nd->atoms)[j];
          int atomC = (nd->atoms)[k];

          //                                      cout<<"checking for triple
          //                                      intersection between
          //                                      "<<atomA<<", "<<atomB<<",
          //                                      "<<atomC<<endl;

          double tempsdf;
          bool previousTripleFound = false;
          bool tripleIntersectionReady = false;
          TripleIntersection *TI = new TripleIntersection();
          int numTripleIntersections = atoms[atomA]->tripleIntersections.size();

          //                                      cout<<"previous intersections:
          //                                      "<<numTripleIntersections<<endl;
          //                                      cout<<"pairIntersect["<<i<<"]["<<j<<"]
          //                                      =
          //                                      ";cout<<pairIntersect[i][j]<<endl;
          //                                      cout<<"pairIntersect["<<j<<"]["<<k<<"]
          //                                      =
          //                                      ";cout<<pairIntersect[j][k]<<endl;
          //                                      cout<<"pairIntersect["<<i<<"]["<<k<<"]
          //                                      =
          //                                      ";cout<<pairIntersect[i][k]<<endl;

          if (pairIntersect[i][j] && pairIntersect[j][k] &&
              pairIntersect[i][k]) // all pairs intersect. so possible triple
                                   // intersection
          {
            numTripleFound++;

            //                                              cout<<"previous
            //                                              intersections:
            //                                              "<<numTripleIntersections<<endl;
            for (int nti = 0; nti < numTripleIntersections; nti++) {
              numOldTripleSearched++;
              TI = atoms[atomA]->tripleIntersections[nti];

              if ((TI->atomA == atomA || TI->atomB == atomA ||
                   TI->atomC == atomA) &&
                  (TI->atomA == atomB || TI->atomB == atomB ||
                   TI->atomC == atomB) &&
                  (TI->atomA == atomC || TI->atomB == atomC ||
                   TI->atomC == atomC)) // already identified the intersection
              {
                numOldTripleUsed++;
                //                                                              cout<<"previous triple intersection between "<<atomA<<", "<<atomB<<", "<<atomC<<" found"<<endl;
                previousTripleFound = true;
                tripleIntersectionReady = true;
                break;
              }
            }

            if (!previousTripleFound) {
              TI = new TripleIntersection();
              if (atoms[atomA]->isTripleIntersect(atoms[atomA], atoms[atomB],
                                                  atoms[atomC], solventRad)) {
                if (computeTripleIntersection(atomA, atomB, atomC, nd->xyz[0],
                                              nd->xyz[1], nd->xyz[2],
                                              solventRad, TI)) {
                  numTripleComputed++;
                  //                                                                      cout<<"new triple intersection between "<<atomA<<", "<<atomB<<", "<<atomC<<" computed"<<endl;
                  tripleIntersectionReady = true;
                }
              }
            }

            if (tripleIntersectionReady) {
              if (getSDFForTripleIntersect(TI, nd->xyz[0], nd->xyz[1],
                                           nd->xyz[2], tempsdf)) {
                if (tempsdf > atoms[atomA]->r + solventRad &&
                    tempsdf > atoms[atomB]->r + solventRad &&
                    tempsdf > atoms[atomC]->r + solventRad) {
                  tripleprob++;
                  // cout<<"In tet, but distance is too large "<<tempsdf<<endl;
                }

                else if (tripleSDF > tempsdf) {
                  tripleSDF = tempsdf;
                  inTet = true;
                }

                triplecount++;
              }
            }
          }

          //                                      cout<<"finished checking for
          //                                      triple intersection between
          //                                      "<<atomA<<", "<<atomB<<",
          //                                      "<<atomC<<endl;

        } // end for k

      } // end for j

    } // end for i

  } // end if inside

  //      cout<<"Triple SDF computed: "<<tripleSDF<<endl;

  /*
          if(inTet)
                  sdf = tripleSDF;
          else if(inCone)
                  sdf = pairSDF;
          else
                  sdf = singleSDF;
  */

  // cout << "Node " << sdf << " " << singleSDF << " " << pairSDF << " " <<
  // tripleSDF << endl;

  if (inside && sdf > singleSDF)
    sdf = singleSDF;
  if (inCone && sdf > pairSDF)
    sdf = pairSDF;
  if (inTet && sdf > tripleSDF)
    sdf = tripleSDF;
  // if (singleSDF < 49.0 && pairSDF < 49.0) {
  //   cout << "???" << endl;
  // }

  if (!inside) {
    sdf = 0.0;
  }

  /*
  if(!insideVDW && sdf > 49.0) {
    cout << "arand " << sdf << " " << inside << " " << inCone << " " << inTet <<
  " "
         << singleSDF << " " << pairSDF << " " << tripleSDF << endl;

    for (int i=0; i<size; i++) {
      cout << "Atom " << i << " "<< insideAtom[i] << endl;
    }

      sdf = 1.41;
  }
  */
  /*
          if(sdf>3 && inside)
          {
                  cout << "Node1 " << sdf << " " << singleSDF << " " << pairSDF
     << " " << tripleSDF << endl;

                  cout<<sdf<<endl;
          }
  */

  nd->v = sdf;
  nd->sdf = sdf;

  delete[] insideAtom;
  for (int i = 0; i < size; i++)
    delete[] pairIntersect[i];
  delete[] pairIntersect;
}

bool Octree::computeTripleIntersection(int atomA, int atomB, int atomC,
                                       double x, double y, double z,
                                       double solventRad,
                                       TripleIntersection *TI) {
  Atom *a1 = atoms[atomA];
  Atom *a2 = atoms[atomB];
  Atom *a3 = atoms[atomC];

  //      cout<<"Inside triple intersection sdf"<<endl;
  // enlarged radii
  double R1 = a1->r + solventRad;
  double R2 = a2->r + solventRad;
  double R3 = a3->r + solventRad;

  // computing distances
  double d12[3], d32[3];
  for (int i = 0; i < 3; i++) {
    d12[i] = a1->xyz[i] - a2->xyz[i];
    d32[i] = a3->xyz[i] - a2->xyz[i];
  }

  double dist12, dist32;
  dist12 = sqrt(d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2]);
  dist32 = sqrt(d32[0] * d32[0] + d32[1] * d32[1] + d32[2] * d32[2]);

  // computing the planes. Two planes are enough
  double plane12[4], plane32[4];

  plane12[3] = 0;
  plane32[3] = 0;
  for (int i = 0; i < 3; i++) {
    plane12[i] = d12[i] / dist12;
    plane32[i] = d32[i] / dist32;

    plane12[3] += a2->xyz[i] * a2->xyz[i] - a1->xyz[i] * a1->xyz[i];
    plane32[3] += a2->xyz[i] * a2->xyz[i] - a3->xyz[i] * a3->xyz[i];
  }

  plane12[3] += R1 * R1 - R2 * R2;
  plane32[3] += R3 * R3 - R2 * R2;

  plane12[3] /= 2 * dist12;
  plane32[3] /= 2 * dist32;

  //      cout<<"sphere 1: ["<<a1.xyz[0]<<", "<<a1.xyz[1]<<", "<<a1.xyz[2]<<",
  //      "<<R1<<"] "<<endl; cout<<"sphere 2: ["<<a2.xyz[0]<<", "<<a2.xyz[1]<<",
  //      "<<a2.xyz[2]<<", "<<R2<<"] "<<endl; cout<<"sphere 3: ["<<a3.xyz[0]<<",
  //      "<<a3.xyz[1]<<", "<<a3.xyz[2]<<", "<<R3<<"] "<<endl;

  //      cout<<"distance between sphere 1 and 2: "<<dist12<<" "<<a1.dist(a2)<<"
  //      "<<a2.dist(a1)<<endl; cout<<"distance between sphere 3 and 2:
  //      "<<dist32<<" "<<a3.dist(a2)<<" "<<a2.dist(a3)<<endl;

  //      cout<<"plane 1: ["<<plane12[0]<<", "<<plane12[1]<<", "<<plane12[2]<<",
  //      "<<plane12[3]<<"] "<<endl; cout<<"plane 2: ["<<plane32[0]<<",
  //      "<<plane32[1]<<", "<<plane32[2]<<", "<<plane32[3]<<"] "<<endl;

  // intersecting the planes to get a line
  double cross[3];
  cross[0] = plane12[1] * plane32[2] - plane12[2] * plane32[1];
  cross[1] = plane12[2] * plane32[0] - plane12[0] * plane32[2];
  cross[2] = plane12[0] * plane32[1] - plane12[1] * plane32[0];

  double dot[3];
  dot[0] = plane12[0] * plane32[0];
  dot[1] = plane12[1] * plane32[1];
  dot[2] = plane12[2] * plane32[2];

  //      cout<<"cross: ["<<cross[0]<<", "<<cross[1]<<", "<<cross[2]<<"]
  //      "<<endl; cout<<"dot: ["<<dot[0]<<", "<<dot[1]<<", "<<dot[2]<<"]
  //      "<<endl;

  double pointOnPlane[3];

  if (abs(cross[0]) >= abs(cross[1]) &&
      abs(cross[0]) >= abs(cross[2])) // set x coordinate to 0
  {
    // cout << "AAA" << endl;
    pointOnPlane[0] = 0.0;
    pointOnPlane[1] =
        (plane32[3] * plane12[2] - plane12[3] * plane32[2]) / cross[0];
    pointOnPlane[2] =
        (plane12[3] * plane32[1] - plane32[3] * plane12[1]) / cross[0];
  } else if (abs(cross[1]) >= abs(cross[0]) &&
             abs(cross[1]) >= abs(cross[2])) // set y coordinate to 0
  {
    // cout << "BBB" << endl;
    pointOnPlane[0] =
        (plane12[3] * plane32[2] - plane32[3] * plane12[2]) / cross[1];
    pointOnPlane[1] = 0.0;
    pointOnPlane[2] =
        (plane32[3] * plane12[0] - plane12[3] * plane32[0]) / cross[1];
  } else if (abs(cross[2]) >= abs(cross[1]) &&
             abs(cross[2]) >= abs(cross[0])) // set z coordinate to 0
  {
    // cout << "CCC" << endl;
    pointOnPlane[0] =
        (plane32[3] * plane12[1] - plane12[3] * plane32[1]) / cross[2];
    pointOnPlane[1] =
        (plane12[3] * plane32[0] - plane32[3] * plane12[0]) / cross[2];
    pointOnPlane[2] = 0.0;
  }

  // The parametric equation of the line is pointOnPlane + t*cross, i.e x(t) =
  // pointOnPlane[0] + t*cross[0] ... plug this into the equation of a sphere,
  // and solve for t t^2 (i^2 + j^2 + k^2) + 2t( (x0-x1)i + (y0-y1)j + (z0-z1)k
  // ) + ( (x0-x1)^2 + (y0-y1)^2 + (z0-z1)^2 - r1^2) = 0 where, cross = [i,j,k];
  // pointOnPlane = [x0,y0,z0]; a1.xyz = [x1,y1,z1]; R1 = r1;

  //      cout<<"point on plane : ["<<pointOnPlane[0]<<", "<<pointOnPlane[1]<<",
  //      "<<pointOnPlane[2]<<"] "<<endl;

  double dx = pointOnPlane[0] - a1->xyz[0];
  double dy = pointOnPlane[1] - a1->xyz[1];
  double dz = pointOnPlane[2] - a1->xyz[2];

  double a = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
  double b = 2.0 * (dx * cross[0] + dy * cross[1] + dz * cross[2]);
  double c = dx * dx + dy * dy + dz * dz - R1 * R1;

  //      cout<<"[dx, dy, dz]: "<<dx<<" "<<dy<<" "<<dz<<endl;
  //      cout<<"quadratic equation: "<<a<<"t^2 + "<<b<<"t +"<<c<<" = 0 "<<endl;

  double det = b * b - 4.0 * a * c;
  if (det < 0) {
    //              cout<<"Error: negative det"<<endl;
    //              cin>>det;
    return false;
  } else
    det = sqrt(det);

  double t1 = (-b + det) / (2.0 * a);
  double t2 = (-b - det) / (2.0 * a);

  double point1[3] = {pointOnPlane[0] + t1 * cross[0],
                      pointOnPlane[1] + t1 * cross[1],
                      pointOnPlane[2] + t1 * cross[2]};
  double point2[3] = {pointOnPlane[0] + t2 * cross[0],
                      pointOnPlane[1] + t2 * cross[1],
                      pointOnPlane[2] + t2 * cross[2]};

  TI->atomA = atomA;
  TI->atomB = atomB;
  TI->atomC = atomC;
  TI->x1 = point1[0];
  TI->y1 = point1[1];
  TI->z1 = point1[2];
  TI->x2 = point2[0];
  TI->y2 = point2[1];
  TI->z2 = point2[2];
  TI->isBuried1 = isBuried(point1);
  TI->isBuried2 = isBuried(point2);
  atoms[atomA]->tripleIntersections.push_back(TI);
  atoms[atomB]->tripleIntersections.push_back(TI);
  atoms[atomC]->tripleIntersections.push_back(TI);

  return true;
}

bool Octree::getSDFForTripleIntersect(TripleIntersection *TI, double x,
                                      double y, double z, double &sdf) {
  double dist1 = 10.0;
  double dist2 = 10.0;
  bool sdfcomputed = false;

  // cout<<"given point: ["<<givenPoint[0]<<", "<<givenPoint[1]<<",
  // "<<givenPoint[2]<<"] "<<endl;

  if (!TI->isBuried1) {
    // cout << "Buried 1" << endl;
    {
      double dx1 = TI->x1 - x;
      double dy1 = TI->y1 - y;
      double dz1 = TI->z1 - z;
      dist1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
      // cout<<"point on tet 1: ["<<point1[0]<<", "<<point1[1]<<",
      // "<<point1[2]<<"]. Distance: "<<dist1<<endl;
      sdfcomputed = true;
    }
  }
  if (!TI->isBuried2) {

    // cout << "Buried 2" << endl;
    {
      double dx2 = TI->x2 - x;
      double dy2 = TI->y2 - y;
      double dz2 = TI->z2 - z;
      dist2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
      // cout<<"point 2: ["<<point2[0]<<", "<<point2[1]<<", "<<point2[2]<<"].
      // Distance: "<<dist2<<endl;
      sdfcomputed = true;
    }
  }

  // cout<<"dist1: "<<dist1<<", dist2: "<<dist2<<endl;

  sdf = dist1 < dist2 ? dist1 : dist2;

  return sdfcomputed;
}

bool Octree::getPairwiseIntersection(
    int atomA, int atomB, double x, double y, double z, double solventRad,
    PairIntersection *PI) // returns true if pairwise intersection of the atoms
                          // are computed successfully.
{
  double cx1 = atoms[atomA]->xyz[0];
  double cy1 = atoms[atomA]->xyz[1];
  double cz1 = atoms[atomA]->xyz[2];
  double R1 = atoms[atomA]->r + solventRad;

  double cx2 = atoms[atomB]->xyz[0];
  double cy2 = atoms[atomB]->xyz[1];
  double cz2 = atoms[atomB]->xyz[2];
  double R2 = atoms[atomB]->r + solventRad;

  double dx = cx1 - cx2;
  double dy = cy1 - cy2;
  double dz = cz1 - cz2;

  double distSq = (dx * dx + dy * dy + dz * dz);
  double dist = sqrt(distSq);

  double R1Sq = R1 * R1;
  double R2Sq = R2 * R2;

  double distanceFromCenter1ToCircle =
      (R1Sq - R2Sq + distSq) /
      (2 *
       dist); // solution of R1^2 - d1^2 = R2^2 - d2^2 = r^2 and d1+d2 = dist
  double distanceFromCenter2ToCircle = dist - distanceFromCenter1ToCircle;

  double cosTheta1 =
      distanceFromCenter1ToCircle / R1; // half-angles subtended from the
                                        // centers to the circle of intersection
  double cosTheta2 = distanceFromCenter2ToCircle / R2;

  double vecC2C1x = dx / dist; // normalized vector pointing from center of atom
                               // 2 to center of atom 1
  double vecC2C1y = dy / dist;
  double vecC2C1z = dz / dist;

  double A = vecC2C1x; // equation of plane containing the circle
  double B = vecC2C1y;
  double C = vecC2C1z;
  double D = ((cx2 * cx2 - cx1 * cx1) + (cy2 * cy2 - cy1 * cy1) +
              (cz2 * cz2 - cz1 * cz1) + (R1Sq - R2Sq)) /
             (2 * dist);

  double circleCenterX =
      (cx1 * distanceFromCenter2ToCircle + cx2 * distanceFromCenter1ToCircle) /
      dist; // center and radius of circle
  double circleCenterY =
      (cy1 * distanceFromCenter2ToCircle + cy2 * distanceFromCenter1ToCircle) /
      dist;
  double circleCenterZ =
      (cz1 * distanceFromCenter2ToCircle + cz2 * distanceFromCenter1ToCircle) /
      dist;
  double circleRad =
      sqrt(R1Sq - distanceFromCenter1ToCircle * distanceFromCenter1ToCircle);

  PI->atomA = atomA;
  PI->atomB = atomB;
  PI->cx = circleCenterX;
  PI->cy = circleCenterY;
  PI->cz = circleCenterZ;
  PI->r = circleRad;
  PI->A = A;
  PI->B = B;
  PI->C = C;
  PI->D = D;
  PI->thetaA = cosTheta1;
  PI->thetaB = cosTheta2;
  atoms[atomA]->pairIntersections.push_back(PI);
  atoms[atomB]->pairIntersections.push_back(PI);

  return true;
}

bool Octree::getSDFForIntersectingSpheres(
    PairIntersection *PI, double x, double y, double z, double solventRad,
    double &sdf) // returns true if reentrant sdf computed (point falls inside
                 // cone).
{
  double cx1 = atoms[PI->atomA]->xyz[0];
  double cy1 = atoms[PI->atomA]->xyz[1];
  double cz1 = atoms[PI->atomA]->xyz[2];
  double R1 = atoms[PI->atomA]->r + solventRad;

  double cx2 = atoms[PI->atomB]->xyz[0];
  double cy2 = atoms[PI->atomB]->xyz[1];
  double cz2 = atoms[PI->atomB]->xyz[2];
  double R2 = atoms[PI->atomB]->r + solventRad;

  double dx1 = cx1 - x;
  double dy1 = cy1 - y;
  double dz1 = cz1 - z;

  double dx2 = cx2 - x;
  double dy2 = cy2 - y;
  double dz2 = cz2 - z;

  double dist1Sq = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  double dist1 = sqrt(dist1Sq);

  double vecPtC1x = dx1 / dist1; // normalized vector pointing from the given
                                 // Point to center of atom 1
  double vecPtC1y = dy1 / dist1;
  double vecPtC1z = dz1 / dist1;

  double dist2Sq = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  double dist2 = sqrt(dist2Sq);

  double vecPtC2x = dx2 / dist2; // normalized vector pointing from the given
                                 // point to center of atom 2
  double vecPtC2y = dy2 / dist2;
  double vecPtC2z = dz2 / dist2;

  double cosThetaPt1 =
      vecPtC1x * PI->A + vecPtC1y * PI->B +
      vecPtC1z *
          PI->C; // angle between the axis and the line connecting Pt and C1
  double cosThetaPt2 =
      -(vecPtC2x * PI->A + vecPtC2y * PI->B +
        vecPtC2z *
            PI->C); // angle between the axis and the line connecting Pt and C2

  if (cosThetaPt1 > PI->thetaA &&
      cosThetaPt2 >
          PI->thetaB) // inside both cones. SDF is distance from Pt to circle
  {
    double pointPlaneDist = PI->A * x + PI->B * y + PI->C * z + PI->D;
    double projectedX = x - PI->A * pointPlaneDist;
    double projectedY = y - PI->B * pointPlaneDist;
    double projectedZ = z - PI->C * pointPlaneDist;

    double onPlaneDirFromCCX = projectedX - PI->cx;
    double onPlaneDirFromCCY = projectedY - PI->cy;
    double onPlaneDirFromCCZ = projectedZ - PI->cz;

    double onPlaneDistanceToCircleCenter =
        sqrt(onPlaneDirFromCCX * onPlaneDirFromCCX +
             onPlaneDirFromCCY * onPlaneDirFromCCY +
             onPlaneDirFromCCZ * onPlaneDirFromCCZ);
    double ratio =
        PI->r /
        onPlaneDistanceToCircleCenter; // possible numerical problem when the
                                       // projection is close to the center.
                                       // have to fix this.

    double closestPointOnCircleX = PI->cx + ratio * onPlaneDirFromCCX;
    double closestPointOnCircleY = PI->cy + ratio * onPlaneDirFromCCY;
    double closestPointOnCircleZ = PI->cz + ratio * onPlaneDirFromCCZ;

    double closestPointOnCircle[3] = {
        closestPointOnCircleX, closestPointOnCircleY, closestPointOnCircleZ};

    // cout<<"Point: ["<<x<<", "<<y<<", "<<z<<"]"<<endl;
    if (!isBuried(closestPointOnCircle)) {
      double dPtX = x - closestPointOnCircleX;
      double dPtY = y - closestPointOnCircleY;
      double dPtZ = z - closestPointOnCircleZ;

      sdf = sqrt(dPtX * dPtX + dPtY * dPtY + dPtZ * dPtZ);

      return true; // reentrant case found
    }
  }

  return false;
}

void Octree::ReadPQR(string fname, string fname2) {

  double xyzmin[3];
  double xyzmax[3];

  for (int i = 0; i < 3; i++) {
    xyzmin[3] = 1000000000.0;
    xyzmax[3] = -1000000000.0;
  }

  ifstream input(fname.c_str());

  string atomname, resname, chain;
  int atomindex, resindex;
  float q;

  string keyword;
  while (input >> keyword) {
    if (strcmp(keyword.c_str(), "ATOM") == 0) {

      Atom *a = new Atom();
      input >> atomindex >> atomname >> resname >> chain >> resindex >>
          a->xyz[0] >> a->xyz[1] >> a->xyz[2] >> q >> a->r;

      // FIXME: this is a hack...
      // a.r += 1.4;

      //      cout << a->xyz[0] << " " << a->xyz[1] << " " << a->xyz[2] << " "
      //      << a->r << endl;

      if (a->r == 0.0)
        a->r = 1.4;

      atoms.push_back(a);
      // a.Print();
      for (int i = 0; i < 3; i++) {
        if (xyzmin[i] > a->xyz[i]) {
          xyzmin[i] = a->xyz[i];
        }
        if (xyzmax[i] < a->xyz[i]) {
          xyzmax[i] = a->xyz[i];
        }
      }
    }
  }

  ifstream input2(fname2.c_str());

  double x, y, z, r;
  int id = 0;
  while (input2 >> keyword) {
    if (strcmp(keyword.c_str(), "ATOM") == 0) {

      input2 >> atomindex >> atomname >> resname >> resindex >> x >> y >> z >>
          q >> r;

      atoms[id]->r = r;

      cout << atoms[id]->xyz[0] << " " << atoms[id]->xyz[1] << " "
           << atoms[id]->xyz[2] << " " << atoms[id]->r << endl;
      id++;
    }
  }

  size = 0.0;
  for (int i = 0; i < 3; i++) {
    center[i] = (xyzmin[i] + xyzmax[i]) * .5;
    if (xyzmax[i] - xyzmin[i] > size) {
      size = xyzmax[i] - xyzmin[i];
    }
  }

  cout << "Number of atoms " << atoms.size() << endl;
  cout << "Size: " << size << endl;
  cout << "Center: " << center[0] << " " << center[1] << " " << center[2]
       << endl;
}

void Octree::ReadPQR(string fname) {

  double xyzmin[3];
  double xyzmax[3];

  for (int i = 0; i < 3; i++) {
    xyzmin[3] = 1000000000.0;
    xyzmax[3] = -1000000000.0;
  }

  ifstream input(fname.c_str());

  string atomname, resname, chain;
  int atomindex, resindex;
  float q;

  string keyword;
  while (input >> keyword) {
    if (strcmp(keyword.c_str(), "ATOM") == 0) {

      Atom *a = new Atom();
      input >> atomindex >> atomname >> resname >> resindex >> a->xyz[0] >>
          a->xyz[1] >> a->xyz[2] >> q >> a->r;

      // FIXME: this is a hack...
      // a.r += 1.4;

      // cout << a.xyz[0] << " " << a.xyz[1] << " " << a.xyz[2] << " " << a.r <<
      // endl;

      atoms.push_back(a);
      // a.Print();
      for (int i = 0; i < 3; i++) {
        if (xyzmin[i] > a->xyz[i]) {
          xyzmin[i] = a->xyz[i];
        }
        if (xyzmax[i] < a->xyz[i]) {
          xyzmax[i] = a->xyz[i];
        }
      }
    }
  }

  size = 0.0;
  for (int i = 0; i < 3; i++) {
    center[i] = (xyzmin[i] + xyzmax[i]) * .5;
    if (xyzmax[i] - xyzmin[i] > size) {
      size = xyzmax[i] - xyzmin[i];
    }
  }

  cout << "Number of atoms " << atoms.size() << endl;
  cout << "Size: " << size << endl;
  cout << "Center: " << center[0] << " " << center[1] << " " << center[2]
       << endl;
}

bool Octree::DoesOverlap(OTNode *nd, int iA) {
  // double epsilon = 0.1;
  double epsilon = 1.4;
  double distsq =
      (nd->xyz[0] - atoms[iA]->xyz[0]) * (nd->xyz[0] - atoms[iA]->xyz[0]) +
      (nd->xyz[1] - atoms[iA]->xyz[1]) * (nd->xyz[1] - atoms[iA]->xyz[1]) +
      (nd->xyz[2] - atoms[iA]->xyz[2]) * (nd->xyz[2] - atoms[iA]->xyz[2]);
  // FIXME: only need .70 not 1.41?
  // if (distsq < (atoms[iA].r + 1.41*nd->w)*(atoms[iA].r + 1.41*nd->w)) {
  //  return true;
  //}

  // for Gaussian?
  // if (distsq < (1.52*atoms[iA].r + 0.71*nd->w)*(1.52*atoms[iA].r +
  // 0.71*nd->w)) {

  // for SAS
  if (distsq < (atoms[iA]->r + 1.4 + 0.87 * nd->w + epsilon) *
                   (atoms[iA]->r + 1.4 + 0.87 * nd->w + epsilon)) {
    return true;
  }
  return false;
}

bool Octree::DoesOverlap(OTNode *nd, Atom *a) {
  // double epsilon = 0.1;
  double epsilon = 1.4;

  double distsq = (nd->xyz[0] - a->xyz[0]) * (nd->xyz[0] - a->xyz[0]) +
                  (nd->xyz[1] - a->xyz[1]) * (nd->xyz[1] - a->xyz[1]) +
                  (nd->xyz[2] - a->xyz[2]) * (nd->xyz[2] - a->xyz[2]);
  // FIXME: only need .70 not 1.41?
  // if (distsq < (atoms[iA].r + 1.41*nd->w)*(atoms[iA].r + 1.41*nd->w)) {
  //  return true;
  //}

  // for SAS
  if (distsq < (a->r + 1.4 + 0.87 * nd->w + epsilon) *
                   (a->r + 1.4 + 0.87 * nd->w + epsilon)) {
    return true;
  }
  return false;
}

void Octree::MarchingCubes() {
  cout << finalCells.size() << endl;
  for (int i = 0; i < finalCells.size(); i++) {
    // vmarchcube(&finalCells[i], contourValue);
    // int result = ImprovedMC(&finalCells[i],
    // finalVertices[finalCells[i].vertices[0]].xyz, contourValue);
    int result = 0;
    if (result == 0) {
      // do some swapping first...

      double tmp;
      tmp = finalCells[i].values[2];
      finalCells[i].values[2] = finalCells[i].values[3];
      finalCells[i].values[3] = tmp;
      tmp = finalCells[i].values[6];
      finalCells[i].values[6] = finalCells[i].values[7];
      finalCells[i].values[7] = tmp;

      int tmp1;
      tmp1 = finalCells[i].vertices[2];
      finalCells[i].vertices[2] = finalCells[i].vertices[3];
      finalCells[i].vertices[3] = tmp1;
      tmp1 = finalCells[i].vertices[6];
      finalCells[i].vertices[6] = finalCells[i].vertices[7];
      finalCells[i].vertices[7] = tmp1;

      vmarchcube(&finalCells[i], contourValue);
    }
  }
}

void Octree::PrintMesh(char *fname) {

  ofstream fout(fname);
  fout << mvx.size() << " " << tri0.size() << endl;
  for (int i = 0; i < mvx.size(); i++) {
    fout << mvx[i] << " " << mvy[i] << " " << mvz[i] << endl;
  }
  for (int i = 0; i < tri0.size(); i++) {
    fout << tri0[i] << " " << tri1[i] << " " << tri2[i] << endl;
  }

  /*
  ofstream fout(fname);
  fout << vertices.size()/3 << " " << triangles.size()/3 << endl;
  for (int i=0; i<vertices.size(); i+=3) {
    fout << vertices[i] << " " << vertices[i+1] << " " << vertices[i+2] << endl;
  }
  for (int i=0; i<triangles.size(); i+=3) {
    fout << triangles[i] << " " << triangles[i+1] << " " << triangles[i+2] <<
  endl;
  }
  */
}

double getArea(double x1, double y1, double z1, double x2, double y2, double z2,
               double x3, double y3, double z3) {
  double dx12 = x1 - x2;
  double dx23 = x2 - x3;
  double dx31 = x3 - x1;

  double dy12 = y1 - y2;
  double dy23 = y2 - y3;
  double dy31 = y3 - y1;

  double dz12 = z1 - z2;
  double dz23 = z2 - z3;
  double dz31 = z3 - z1;

  double A2 = dx12 * dx12 + dy12 * dy12 + dy12 * dy12;
  double B2 = dx23 * dx23 + dy23 * dy23 + dy23 * dy23;
  double C2 = dx31 * dx31 + dy31 * dy31 + dy31 * dy31;

  //      return 0.25 * sqrt( 2 * (A2*B2 + B2*C2 + C2*A2) - (A2*A2 + B2*B2 +
  //      C2*C2) );

  double A = sqrt(A2);
  double B = sqrt(B2);
  double C = sqrt(C2);

  double S = (A + B + C) / 2.0;

  return sqrt(S * (S - A) * (S - B) * (S - C));
}

double Octree::areaMesh() {

  double area = 0.0;
  cout << "Number of vertices " << mvx.size() << endl;
  cout << "Number of triangles " << tri0.size() << endl;

  for (int i = 0; i < tri0.size(); i++) {
    int v1 = tri0[i];
    int v2 = tri1[i];
    int v3 = tri2[i];

    area += getArea(mvx[v1], mvy[v1], mvz[v1], mvx[v2], mvy[v2], mvz[v2],
                    mvx[v3], mvy[v3], mvz[v3]);
  }

  return area;
}

void Octree::ClearIntermediateData() {

  finalCells.clear();
  uncomputedCells.clear();
  atoms.clear();
  finalVertices.clear();
  edgeToVert.clear();

  std::vector<double> vertices;
  std::vector<int> triangles;
}

double Octree::getInterfaceAreaUnderRoot(OTNode *root, Octree *otherTree,
                                         double *trans, double interfaceWidth) {
  double areaInner = 0.0;
  double areaOuter = 0.0;

  int otherBaseSize = otherTree->basesize;

  for (int i = 0; i < otherBaseSize; i++) {
    for (int j = 0; j < otherBaseSize; j++) {
      for (int k = 0; k < otherBaseSize; k++) {
        OTNode *otherRoot = otherTree->GetNode(i, j, k);
        double areaI = 0.0;
        double areaO = 0.0;
        root->getInterfaceArea(otherRoot, trans, interfaceWidth, areaI, areaO,
                               contourValue);
        areaInner += areaI;
        areaOuter += areaO;
      }
    }
  }

  return (areaInner + areaOuter) / 2.0;
}

double Octree::getInterfaceArea(Octree *otherTree, double *trans,
                                double interfaceWidth) {
  markedNodes.clear();
  double area = 0.0;

  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *aRoot = GetNode(i, j, k);
        area +=
            getInterfaceAreaUnderRoot(aRoot, otherTree, trans, interfaceWidth);
      }
    }
  }

  return area;
}

void Octree::clear() {
  int size = atoms.size();
  for (int i = 0; i < size; i++) {
    atoms[i]->pairIntersections.clear();
    atoms[i]->tripleIntersections.clear();
  }
  atoms.clear();
}

int main(int argc, char **argv) {

  if (argc < 8) {
    cout << "Usage: Octree <input pqr file> <output raw file> <initial "
            "gridsize> <tree depth> <target resolution> <isovalue> <blobbyness>"
         << endl;
  }

  double startT = getTime();

  //  Octree T(atoi(argv[3]),argv[1],argv[8]);
  Octree T(atoi(argv[3]), argv[1]);

  T.targetLevel = atoi(argv[4]);
  T.targetRes = atof(argv[5]);
  T.contourValue = atof(argv[6]);
  T.blobby = atof(argv[7]);

  cout << "Target resolution " << T.targetRes << endl;
  cout << "Target level " << T.targetLevel << endl;

  T.SetInitialData();

  //  for (int i=0; i<T.targetLevel+3; i++) {
  bool notFinished = true;
  int i = 0;
  while (notFinished) {
    cout << "Iteration " << ++i << endl;
    T.MarkIsoContour();
    T.SplitMarkedCells();
    notFinished = T.ComputeNewCellValues();
  }

#ifdef GAUSS
//  T.MarkIsoContour(); // boundary cells are marked
#endif

  double endT = getTime();
  cout << "Gridding time " << endT - startT << endl;

#ifdef VDW
  startT = getTime();
  T.ComputeFinalValues();
  endT = getTime();
  cout << "SDF computation time " << endT - startT << endl;

  cout << "Number of times single atom SDF is overwritten by pair/triplewise "
          "SDF: "
       << T.numoverwrite << endl;
  cout << T.triplecount << " " << T.paircount << endl;
  cout << T.tripleprob << " " << T.pairprob << endl;
#endif

  startT = getTime();

  cout << "A" << endl;
  T.LabelFinalVertices();
  cout << "B" << endl;
  T.GetFinalCellsAndDestroyTree();
  cout << "C" << endl;
  T.MarchingCubes();
  cout << "D" << endl;
  T.ClearIntermediateData();
  T.PrintMesh(argv[2]);
  endT = getTime();

  cout << T.numSDFCalled << endl;
  cout << T.numPairFound << endl;
  cout << T.numOldPairSearched << endl;
  cout << T.numOldPairUsed << endl;
  cout << T.numPairComputed << endl;
  cout << T.numTripleFound << endl;
  cout << T.numOldTripleSearched << endl;
  cout << T.numOldTripleUsed << endl;
  cout << T.numTripleComputed << endl;

  T.clear();

  cout << "Meshing Time " << endT - startT << endl;
}

void Octree::constructOctree(int targLev, double isoVal, double bl) {
  targetLevel = targLev;
  contourValue = isoVal;
  blobby = bl;
  SetInitialData();

  bool notFinished = true;
  int i = 0;
  while (notFinished) {
    cout << "Iteration " << ++i << endl;
    MarkIsoContour();
    SplitMarkedCells();
    notFinished = ComputeNewCellValues();
  }
  MarkIsoContour(); // boundary cells are marked
}

double Octree::computeMarkedArea(OTNode *nd) {
  if (nd->l == targetLevel && nd->interfaceFlag == true)

  {
    OTNode *nbs[8];
    GetFinalCell(nd, nbs);

    bool sameSize = true;
    bool containsContour = false;

    for (int i = 1; i < 8; i++) {
      if ((nbs[0]->v < contourValue && nbs[i]->v > contourValue) ||
          (nbs[0]->v > contourValue && nbs[i]->v < contourValue)) {
        containsContour = true;
      }
      if (nbs[i]->l != nbs[0]->l) {
        sameSize = false;
      }
    }

    int marchingCubeLayout[8] = {0, 4, 2, 6, 1, 5, 3, 7};
    Cell c;

    if (sameSize && containsContour) {
      for (int i = 0; i < 8; i++) {
        //                              if (nbs[i]->flag < 0)
        {
          nbs[i]->flag = finalVertices.size();
          Vertex tmp;
          tmp.xyz[0] = nbs[i]->xyz[0];
          tmp.xyz[1] = nbs[i]->xyz[1];
          tmp.xyz[2] = nbs[i]->xyz[2];
          finalVertices.push_back(tmp);

          int index = marchingCubeLayout[i];
          c.vertices[index] = nbs[i]->flag;
          c.values[index] = nbs[i]->v;
          //                                      cout<<c.values[index]<<" ";
        }
      }
      c.span = nd->w;
      //                      cout<<nd->w<<" "<<endl;
      finalCells.push_back(c);
      nd->dual = &c;
    }
  }
  if (nd->c != NULL) {
    for (int i = 0; i < 8; i++) {
      computeMarkedArea(nd->c[i]);
    }
  }
}

double Octree::computeMarkedArea() {
  for (int i = 0; i < basesize; i++) {
    for (int j = 0; j < basesize; j++) {
      for (int k = 0; k < basesize; k++) {
        OTNode *nd = GetNode(i, j, k);
        computeMarkedArea(nd);
      }
    }
  }
  cout << "Number of final dual cells: " << finalCells.size() << endl;

  MarchingCubes();

  return areaMesh();
}

/* Use for comparing interface area computation */
/*
int main(int argc, char ** argv) {

  if (argc < 9) {
    cout << "Usage: Octree <input pqr file1> <input pqr file2> <output interface
raw> <initial gridsize> <tree depth> <isovalue> <blobbyness> <interfaceWidth>"
<< endl;
  }

  double startT = getTime();
  Octree* T1 = new Octree(atoi(argv[4]), argv[1]);
  T1->constructOctree(atoi(argv[5]), atof(argv[6]), atof(argv[7]));
  double endT = getTime();
  double time = endT-startT;
  cout<<"Time for constructing octree 1 = "<<time<<endl;

  startT = getTime();
  Octree* T2 = new Octree(atoi(argv[4]), argv[2]);
  T2->constructOctree(atoi(argv[5]), atof(argv[6]), atof(argv[7]));
  endT = getTime();
  time = endT-startT;
  cout<<"Time for constructing octree 2 = "<<time<<endl;

  startT = getTime();
  double trans[12] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
0.0}; double area = T1->getInterfaceArea(T2, trans, atof(argv[8])); endT =
getTime(); time = endT-startT; cout<<"Time for approximating interface area =
"<<time<<endl; cout<<"Interface area = "<<area<<endl;

  startT = getTime();
  area = T1->computeMarkedArea();
  endT = getTime();
  time = endT-startT;
  cout<<"Time for computing interface area = "<<time<<endl;
  cout<<"Interface area = "<<area<<endl;

  cout<<argv[3]<<endl;
  T1->PrintMesh(argv[3]);
}*/

/*
int main(int argc, char ** argv) {

  if (argc < 8) {
    cout << "Usage: Octree <input pqr file> <output raw file> <initial gridsize>
<tree depth> <target resolution> <isovalue> <blobbyness>" << endl;
  }

  double startT = getTime();

  Octree T(atoi(argv[3]),argv[1]);

  T.targetLevel = atoi(argv[4]);
  T.targetRes = atof(argv[5]);
  T.contourValue = atof(argv[6]);
  T.blobby = atof(argv[7]);

  cout<<"Target resolution "<<T.targetRes<<endl;
  cout<<"Target level "<<T.targetLevel<<endl;

  T.SetInitialData();

  bool notFinished = true;
  int i=0;
  while (notFinished) {
    cout << "Iteration " << ++i << endl;
    T.MarkIsoContour();
    T.SplitMarkedCells();
    notFinished = T.ComputeNewCellValues();
  }
//  T.MarkIsoContour(); // boundary cells are marked

  double endT = getTime();
  cout<<"Gridding time "<<endT-startT<<endl;

  startT = getTime();
  cout << "A" << endl;
  T.LabelFinalVertices();
  cout << "B" << endl;
  T.GetFinalCellsAndDestroyTree();
  cout << "C" << endl;
  T.MarchingCubes();
  T.ClearIntermediateData();
  T.PrintMesh(argv[2]);
  endT = getTime();

  cout<<"Meshing Time "<<endT-startT<<endl;

}
*/

/*
int main(int argc, char ** argv) {

  if (argc < 9) {
    cout << "Usage: Octree <input pqr file> <output raw file> <initial gridsize>
<tree depth> <target resolution> <isovalue> <blobbyness> <atomList>" << endl;
  }

  double startT = getTime();

  Octree T(atoi(argv[3]),argv[1]);

  T.targetLevel = atoi(argv[4]);
  T.targetRes = atof(argv[5]);
  T.contourValue = atof(argv[6]);
  T.blobby = atof(argv[7]);

  cout<<"Target resolution "<<T.targetRes<<endl;
  cout<<"Target level "<<T.targetLevel<<endl;

  T.SetInitialData();

  //  for (int i=0; i<T.targetLevel+3; i++) {
  bool notFinished = true;
  int i=0;
  while (notFinished) {
    cout << "Iteration " << ++i << endl;
    T.MarkIsoContour();
    T.SplitMarkedCells();
    notFinished = T.ComputeNewCellValues();
  }
//  T.MarkIsoContour(); // boundary cells are marked

  double endT = getTime();
  cout<<"Gridding time "<<endT-startT<<endl;

  startT = getTime();
  T.LabelFinalVertices();
  T.GetFinalCellsAndDestroyTree();
  T.MarchingCubes();
  T.ClearIntermediateData();
  T.PrintMesh(argv[2]);
  endT = getTime();

  cout<<"Meshing Time "<<endT-startT<<endl;


        FILE *fp = fopen(argv[8],"rt");

        if(fp == NULL)
        {
                cout<<"Could not open interface atom list"<<endl;
                return -1;
        }

        vector<int> interfaceAtomIds;
        int a;
        while(fscanf(fp,"%d",&a)==1)
                interfaceAtomIds.push_back(a);

        int numInterfaceAtoms = interfaceAtomIds.size();
        int numAtoms = atoms.size();

        for(int k=400; k>=4; k/=2)
        {
                int numAtomsToMove = numAtoms/k;

                startT = getTime();

                for(int j=1; j<numAtomsToMove; j++)
                {
                        int randId = rand()%numInterfaceAtoms;
                        int randMovementX = rand()%2 - rand()%2;
                        int randMovementY = rand()%2 - rand()%2;
                        int randMovementZ = rand()%2 - rand()%2;

                        Atom originalAtom = atoms[randId];
                        Atom movedAtom;
                        movedAtom.xyz[0] = originalAtom.xyz[0] +
(double)randMovementX*0.5; movedAtom.xyz[1] = originalAtom.xyz[1] +
(double)randMovementY*0.5; movedAtom.xyz[2] = originalAtom.xyz[2] +
(double)randMovementZ*0.5; movedAtom.r = originalAtom.r;

                        T.RemoveAtom(&originalAtom);
                        T.AddAtom(&movedAtom);

                        atoms[randId].xyz[0] = movedAtom.xyz[0];
                        atoms[randId].xyz[1] = movedAtom.xyz[1];
                        atoms[randId].xyz[2] = movedAtom.xyz[2];
                }

                endT = getTime();

                cout<<"Time to move "<<numAtomsToMove<<" atoms: "<<endT -
startT<<endl;
        }
}
*/
