/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <DynamicAdaptiveGrid/DynamicAdaptiveGrid.h>
#include <PDBParser/Atom.h>
#include <time.h>
#include <iostream>
using namespace std;

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

int main(int argc, char *argv[]) {
  /*      DynamicAdaptiveGrid::DynamicAdaptiveGrids *dag = new
  DynamicAdaptiveGrid::DynamicAdaptiveGrids(0.7,3.0,1.4,-10,-10,-10,10,10,10);
          cout<<"Created DAG"<<endl;

          PDBParser::Atom* a1 = new PDBParser::Atom();
          a1->setPosition(0,0,0);
          a1->setRadius(1.0);
          dag->addAtom(a1);
          cout<<"Added 1 Atom"<<endl;

  //      dag->printSurface("test");

          PDBParser::Atom* a2 = new PDBParser::Atom();
          a2->setPosition(0,2.0,0);
          a2->setRadius(2.0);
          dag->addAtom(a2);
          cout<<"Added 2 Atoms"<<endl;

  //      dag->printSurface("test2");

          PDBParser::Atom* a3 = new PDBParser::Atom();
          a3->setPosition(2.0,0,0);
          a3->setRadius(1.5);
          dag->addAtom(a3);
          cout<<"Added 3 Atoms"<<endl;

          dag->printSurface("test3");

  //      dag->removeAtom(a3);
  //      cout<<"Removed 1 Atom"<<endl;
  //      dag->printSurface("test4.xyzr");*/

  double xyzmin[3];
  double xyzmax[3];

  for (int i = 0; i < 3; i++) {
    xyzmin[i] = 1000.0;
    xyzmax[i] = -1000.0;
  }

  string origPQRName(argv[1]);
  string origRawName(argv[2]);

  cout << origPQRName << endl;
  cout << origPQRName << endl;
  //      string movementFileName(argv[3]);
  //      string movedRawName(argv[4]);

  ifstream input(origPQRName.c_str());

  string atomname, resname;
  int atomindex, resindex;
  float q, x, y, z, r;
  vector<PDBParser::Atom *> atoms;

  string keyword;
  while (input >> keyword) {
    if (strcmp(keyword.c_str(), "ATOM") == 0) {
      PDBParser::Atom *a = new PDBParser::Atom();
      input >> atomindex >> atomname >> resname >> resindex >> x >> y >> z >>
          q >> r;

      a->setPosition(x, y, z);
      a->setRadius(r);
      atoms.push_back(a);

      if (x < xyzmin[0])
        xyzmin[0] = x;
      if (y < xyzmin[1])
        xyzmin[1] = y;
      if (z < xyzmin[2])
        xyzmin[2] = z;
      if (x > xyzmax[0])
        xyzmax[0] = x;
      if (y > xyzmax[1])
        xyzmax[1] = y;
      if (z > xyzmax[2])
        xyzmax[2] = z;

      cout << x << " " << y << " " << z << endl;
      if (x < -100 || y < -100 || z < -100)
        cin >> x;
    }
  }

  //      cout<<xyzmin[0]<<" "<<xyzmin[1]<<" "<<xyzmin[2]<<" "<<xyzmax[0]<<"
  //      "<<xyzmax[1]<<" "<<xyzmax[2]<<endl;

  int numAtoms = atoms.size();
  double startT = getTime();

  DynamicAdaptiveGrid::DynamicAdaptiveGrids *dag2 =
      new DynamicAdaptiveGrid::DynamicAdaptiveGrids(
          0.7, 5.0, atoi(argv[3]), 1.4, 1.4, xyzmin[0] - 4, xyzmin[1] - 4,
          xyzmin[2] - 4, xyzmax[0] + 4, xyzmax[1] + 4, xyzmax[2] + 4);
  dag2->addAtoms(atoms);

  double endT = getTime();

  //      cout<<"Time for inserting " << numAtoms <<" atoms ";
  cout << numAtoms << " " << endT - startT << " " << endl;

  dag2->printSurface(origRawName);

  /*      ifstream input2(movementFileName.c_str());
          double dist;


          vector<PDBParser::Atom*> atoms1;
          vector<PDBParser::Atom*> atoms2;

          int numMoved = 0;
          while(input2>>atomindex)
          {
                  PDBParser::Atom* a = new PDBParser::Atom();

                  input2 >> x >> y >> z >> dist;
  //              cout<<atoms[atomindex-1]->getX()<<"
  "<<atoms[atomindex-1]->getY()<<" "<<atoms[atomindex-1]->getZ()<<endl;
                  atoms1.push_back(atoms[atomindex-1]);

  //              cout<<x<<" "<<y<<" "<<z<<endl;
                  a->setPosition(x,y,z);
                  a->setRadius(atoms[atomindex-1]->getRadius());

                  atoms2.push_back(a);
                  numMoved++;
          }
          startT = getTime();
          dag2->removeAtoms(atoms1);
  //      dag2->printSurface(movedRawName);

          dag2->addAtoms(atoms2);
          endT = getTime();
  //      cout<<"Time for moving " << numAtoms <<" atoms ";
          cout<<numMoved<<" "<< endT - startT<<endl;

          dag2->printSurface(movedRawName);       */

  return 0;
}
