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
#include <InterfaceStats/ComputeInterfaceStats.h>
#include <UsefulMath/Matrix.h>

int main(int argc, char **argv) {
  if (argc < 6) {
    printf("Usage: testInterfaceStats <mol1.pqr> <mol2.pqr> <mol1.quad> "
           "<mol2.quad> <xforms.txt>  \n\n (output will be printed to standard "
           "out. use redirection if you prefer)\n\n");
    return -1;
  }

  FILE *xformFile = fopen(argv[5], "rt");
  double score, rmsd;

  if (xformFile != NULL) {
    int numXForm;
    fscanf(xformFile, "%d", &numXForm);

    printf("Trying to read %d xfroms\n\n", numXForm);

    for (int i = 0; i < numXForm; i++) {
      vector<Matrix> transformations;
      Matrix mtx;

      printf("Reading matrix %d\n", i);

      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          double mtx_jk;
          fscanf(xformFile, "%lf", &mtx_jk);
          printf("%lf ", mtx_jk);
          mtx.set(j, k, mtx_jk);
        }
        printf("\n");
      }

      fscanf(xformFile, "%lf", &score);
      fscanf(xformFile, "%lf", &rmsd);

      transformations.push_back(mtx);

      ComputeInterfaceStats *compStats = new ComputeInterfaceStats(
          argv[1], argv[2], argv[3], argv[4], transformations);
    }

    fclose(xformFile);
    return 0;
  }
  return -2;
}
