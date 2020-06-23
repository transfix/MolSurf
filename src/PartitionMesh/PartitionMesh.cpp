/*
  Copyright 2011 The University of Texas at Austin

	Authors: Alex Rand <arand@ices.utexas.edu>
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
#include <Utility/utility.h>
#include <PartitionMesh/PartitionMesh.h>

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>

using namespace CVCUtility;
using namespace std;

void usagePartitionMesh()
{
	cout	<< endl << "Usage: MolSurf -partitionMesh <raw input file> <# of partitions>"
		<< endl << "    Partition mesh for parallel computation."
		<< endl ;
}

bool partitionMesh(int argc, char ** argv) {
  if(argc != 4) {
    usagePartitionMesh();
    return false;
  }

  int nP = atoi(argv[3]);

  int nV,nT;
  vector<double>x;
  vector<double>y;
  vector<double>z;
  vector<int> t1;
  vector<int> t2;
  vector<int> t3;  

  vector<double> n1;
  vector<double> n2;
  vector<double> n3;

  bool haveNormal = endsWith(argv[2],".rawn") || endsWith(argv[2],".rawnc");

  if (endsWith(argv[2],".rawnc")) {
    cout << "WARNING: .rawnc not yet supported..." << endl;
  }


  // read the surface
  ifstream fin(argv[2]);
  fin >> nV >> nT;
  for (int i=0; i<nV; i++) {
    double xx,yy,zz,nx,ny,nz;
    fin >> xx >> yy >>zz;
    x.push_back(xx);
    y.push_back(yy);
    z.push_back(zz);
    if (haveNormal) {
      fin >> nx >> ny >> nz;
      n1.push_back(nx);
      n2.push_back(ny);
      n3.push_back(nz);
    }
  }

  for (int i=0; i<nT; i++) {
    int tt1, tt2, tt3;
    fin >> tt1 >> tt2 >> tt3;
    t1.push_back(tt1);
    t2.push_back(tt2);
    t3.push_back(tt3);
  }


  // do more stuff...
  int * tstart = new int[nP];
  int * tend = new int[nP];

  cout << "nT: " << nT << endl;
  cout << "nP: " << nP << endl;
  for (int i=0; i<nP; i++) {
    if (i > 0)
      tstart[i] = tend[i-1]+1;
    else
      tstart[i] = 0;

    if (i<nP-1) {
      tend[i] = tstart[i] + nT/nP;
      if (nT%nP <= i)
	tend[i]--;
    } else {
      tend[i] = nT-1;
    }
  }

  cout << "Intervals:";
  for (int i=0; i<nP; i++) {
    cout << " " << tend[i] - tstart[i] +1;
  }cout << endl;


  vector<int> partId;
  for (int i=0; i<x.size(); i++) {
    partId.push_back(-1);
  }
  int * old2new = new int[nV];
  int * new2old = new int[nV];

  vector<set<int> > ghostverts;
  ghostverts.resize(nP);

  int * vstart = new int[nP];
  int * vend = new int[nP];


  // assign vertices
  int count = 0;
  int * nVPar;
  nVPar = new int[nP];
  cout << "Vertex allocations:";
  for (int iP=0; iP<nP; iP++) {
    nVPar[iP] = 0;
    vstart[iP] = count;
    for (int iT=tstart[iP]; iT<=tend[iP]; iT++) {
      if (partId[t1[iT]] == -1) {
	nVPar[iP]++;
	partId[t1[iT]]=iP;
	old2new[t1[iT]] = count;
	new2old[count] = t1[iT];
	count++;
      } else if (partId[t1[iT]] != iP) {
	ghostverts[iP].insert(t1[iT]);
      }
      if (partId[t2[iT]] == -1) {
	nVPar[iP]++;
	partId[t2[iT]]=iP;
	old2new[t2[iT]] = count;
	new2old[count] = t2[iT];
	count++;
      } else if (partId[t2[iT]] != iP) {
	ghostverts[iP].insert(t2[iT]);
      }
      if (partId[t3[iT]] == -1) {
	nVPar[iP]++;
	partId[t3[iT]]=iP;
	old2new[t3[iT]] = count;
	new2old[count] = t3[iT];
	count++;
      } else if (partId[t3[iT]] != iP) {
	ghostverts[iP].insert(t3[iT]);
      }     
    }
    vend[iP] = count-1;

    cout << " : " << nVPar[iP] << "," << ghostverts[iP].size();
  } cout << endl;
  
  int sum = 0;
  for (int iP=0; iP<nP; iP++) {
    sum += nVPar[iP];
  }
  cout << "nV: " << nV << " " << sum << endl;


  // write out the results...
  for (int iP=0; iP<nP; iP++) {
    char fname[256];
    if (haveNormal)
      sprintf(fname,"%s.%d.rawn",argv[2],iP);
    else 
      sprintf(fname,"%s.%d.raw",argv[2],iP);
    ofstream fout(fname);
    
    fout << nVPar[iP] + ghostverts[iP].size() << " " 
	 << tend[iP] - tstart[iP] +1 <<endl;

    // print the vertices
    for (int i=vstart[iP]; i<= vend[iP]; i++) {
      fout << x[new2old[i]] << " "
	   << y[new2old[i]] << " "
	   << z[new2old[i]];

      if (haveNormal) {
	fout << " " << n1[new2old[i]] 
	     << " " << n2[new2old[i]] 
	     << " " << n3[new2old[i]];	
      }
     
      fout << endl;

    }

    map<int,int> gv2new;
    int count = nVPar[iP];
    // print the ghost verts
    set<int>::iterator it;
    for (it=ghostverts[iP].begin(); it!=ghostverts[iP].end(); it++) {
      fout << x[*it] << " " << y[*it] << " " << z[*it];
      if (haveNormal) {
	fout << " " << n1[*it] << " " << n2[*it] << " " << n3[*it];	
      }
      fout << endl;
      gv2new[*it] = count;
      count++;
    }

    for (int i=tstart[iP]; i<=tend[iP]; i++) {
      if (t1[i]>= vstart[iP] && t1[i] <= vend[iP]) {
	fout << old2new[t1[i]]-vstart[iP];
      } else {
	fout << gv2new[t1[i]];
      }
      fout << " ";
      if (t2[i]>= vstart[iP] && t2[i] <= vend[iP]) {
	fout << old2new[t2[i]]-vstart[iP];
      } else {
	fout << gv2new[t2[i]];
      }
      fout << " ";
      if (t3[i]>= vstart[iP] && t3[i] <= vend[iP]) {
	fout << old2new[t3[i]]-vstart[iP];
      } else {
	fout << gv2new[t3[i]];
      }
      fout << endl;


    }

    // print global indices for the ghost vertices...
    fout << ghostverts[iP].size() << endl;
    for (it=ghostverts[iP].begin(); it!=ghostverts[iP].end(); it++) {
      fout << old2new[*it] << endl;
      
    }

    fout.close();
  }

  return true;
} 
