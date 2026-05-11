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

#include <InterfaceStats/InterfaceFootprintOverlap.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char **argv) {
  if (argc < 9) {
    printf("Usage: testInterfaceOverlap <receptor.rawn> <ligand1.rawn> "
           "<ligand2.rawn> <generate highlighted surfaces? [0=no, 1=yes]> "
           "<interface width> <xformFile> <startIndex> <endIndex>");
    return -1;
  }

  FILE *xformFile = fopen(argv[6], "rt");
  int startIndex = atoi(argv[7]);
  int endIndex = atoi(argv[8]);
  int genSurfMode = atoi(argv[4]);
  double interfaceWidth = atof(argv[5]);

  InterfaceFootprintOverlap *ifo =
      new InterfaceFootprintOverlap(std::string(argv[1]), std::string(argv[2]),
                                    std::string(argv[3]), interfaceWidth, true);

  Matrix mtx;

  if (xformFile == NULL) {
    printf("Could not open xform file\n");
    return -2;
  }

  int numXForm;
  fscanf(xformFile, "%d", &numXForm);

  if (startIndex >= numXForm || startIndex < 0 || endIndex >= numXForm ||
      endIndex < 0 || endIndex < startIndex) {
    printf("Xform Index out of bounds\n");
    return -3;
  }

  double dummy;
  for (int i = 0; i <= endIndex; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        double mtx_jk;
        fscanf(xformFile, "%lf", &mtx_jk);
        // printf("%lf ", mtx_jk);
        mtx.set(j, k, mtx_jk);
      }
    }

    fscanf(xformFile, "%lf", &dummy);
    fscanf(xformFile, "%lf", &dummy);

    if (i >= startIndex) {
      if (ifo->computeInterfaces(mtx)) {
        printf("Area of the interface of receptor with ligand 1: %lf\n",
               ifo->getAreaOfReceptorInterfaceWithStaticLigand());
        printf("Area of the interface of receptor with ligand 2: %lf\n",
               ifo->getAreaOfReceptorInterfaceWithMovingLigand());
        printf("Area of the interface of ligand 1 with receptor: %lf\n",
               ifo->getAreaOfStaticLigandInterface());
        printf("Area of the interface of ligand 2 with receptor: %lf\n",
               ifo->getAreaOfMovingLigandInterface());
        printf("Area of the overlap: %lf\n", ifo->getAreaOfOverlap());

        if (genSurfMode == 1) {
          ifo->prepareStaticInterfaceSurfaceFiles();
          ifo->prepareMovingInterfaceSurfaceFiles(i + 1);
        }
      } else {
        printf("ERROR.\n ERROR.ERROR.\n ERROR.ERROR.ERROR.\n");
      }
    }
  }

  fclose(xformFile);
}
