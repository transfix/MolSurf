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
#include <stdio.h>
#include <DPG/PG.h>
#include <string.h>

using namespace DPG;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Input text file not specified!\n");
    return 1;
  }

  char *Molecule;
  char *Update;
  char *Range;
  char *Output;

  char s[2000];
  char key[500], val[500];
  FILE *fp;

  fp = fopen(argv[1], "r");

  if (fp == NULL) {
    printf("Failed to open parameter file %s!\n", argv[1]);
    return 1;
  }

  Molecule = NULL;
  Update = NULL;
  Range = NULL;
  Output = NULL;

  while (fgets(s, 1999, fp) != NULL) {
    if (sscanf(s, "%s %s", key, val) != 2)
      continue;

    if (!strcasecmp(key, "moleculeFile"))
      Molecule = strdup(val);
    else if (!strcasecmp(key, "updateFile"))
      Update = strdup(val);
    else if (!strcasecmp(key, "queryFile"))
      Range = strdup(val);
    else if (!strcasecmp(key, "outFile"))
      Output = strdup(val);
  }

  if (fp)
    fclose(fp);

  if (Molecule == NULL) {
    printf("Missing molecule file name for the static molecule!\n");
    return 1;
  }

  /* Creating DPG */

  PG *pg = new PG(10.0, 1000.0, 3.0);
  vector<Point *> pts;
  int numPoints = 0;
  double x, y, z, r;

  /* inserting atoms into DPG */

  fp = fopen(Molecule, "r");
  if (fp == NULL) {
    printf("Molecule File not found\n");
    return false;
  }
  while (!feof(fp)) {
    fscanf(fp, "%lf", &x);
    fscanf(fp, "%lf", &y);
    fscanf(fp, "%lf", &z);
    fscanf(fp, "%lf", &r);

    Ball *p = new Ball(x, y, z, r);

    pts.push_back(p);

    pg->addPoint(p);
    numPoints++;
  }
  if (fp)
    fclose(fp);

  /* updating atoms */

  if (Update == NULL) {
    printf("No updates specified\n");
  } else {
    fp = fopen(Update, "r");

    if (fp == NULL) {
      printf("Update File not found\n");
      return 1;
    }

    double x1, y1, z1, r1;
    int type;

    while (!feof(fp)) {
      fscanf(fp, "%d", &type);

      if (type == 1) {
        fscanf(fp, "%lf", &x);
        fscanf(fp, "%lf", &y);
        fscanf(fp, "%lf", &z);
        fscanf(fp, "%lf", &r);

        Ball *p = new Ball(x, y, z, r);

        pg->removePoint(p);
      } else if (type == 2) {
        fscanf(fp, "%lf", &x);
        fscanf(fp, "%lf", &y);
        fscanf(fp, "%lf", &z);
        fscanf(fp, "%lf", &r);

        Ball *p = new Ball(x, y, z, r);

        fscanf(fp, "%lf", &x1);
        fscanf(fp, "%lf", &y1);
        fscanf(fp, "%lf", &z1);
        fscanf(fp, "%lf", &r1);

        Ball *p1 = new Ball(x1, y1, z1, r1);

        pg->move(p, p1);
      }
    }

    if (fp)
      fclose(fp);
  }

  /* Performing queries */

  if (Range == NULL) {
    printf("No queries specified\n");
  } else {
    if (Output == NULL) {
      printf("Output file not specified. Using default filename.\n");
      Output = strdup("testOutput.txt");
    }

    fp = fopen(Range, "r");

    if (fp == NULL) {
      fclose(fp);
      printf("Query File not found\n");
      return 1;
    }

    FILE *fp1 = fopen(Output, "w");

    if (fp1 == NULL) {
      printf("Output File not found\n");
      if (fp)
        fclose(fp);
      return 1;
    }

    vector<Point *> results;

    while (!feof(fp)) {
      results.clear();

      fscanf(fp, "%lf", &x);
      fscanf(fp, "%lf", &y);
      fscanf(fp, "%lf", &z);
      fscanf(fp, "%lf", &r);

      Ball *p = new Ball(x, y, z, r);

      pg->range(p, r, results);

      int size = results.size();

      if (size == 0) {
        fprintf(
            fp1,
            "No points found within distance %lf of point %lf, %lf, %lf\n\n", r,
            x, y, z);
      } else {
        fprintf(fp1,
                "%d points found within distance %lf of point %lf, %lf, %lf\n",
                size, r, x, y, z);

        for (int i = 0; i < size; i++) {
          Ball *tempRes = (Ball *)results[i];
          fprintf(fp1, "%lf %lf %lf %lf\n", tempRes->getX(), tempRes->getY(),
                  tempRes->getZ(), tempRes->getRadius());
        }
      }
    }

    if (fp)
      fclose(fp);
    if (fp1)
      fclose(fp1);
  }

  return 0;
}
