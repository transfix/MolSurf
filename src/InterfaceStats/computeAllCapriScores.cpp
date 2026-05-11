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
#include <InterfaceStats/CapriCriteria.h>
#include <UsefulMath/Matrix.h>

int main(int argc, char *argv[]) {
  if (argc < 6) {
    cout << "Usage: " << endl;
    cout << "computeAllCapriScores <receptor pqr> <correct ligand pqr> "
            "<docking output file> <capri score file> <capri summary file> "
            "[optional:poseindex]"
         << endl;
    cout << "Docking output file must start with the number of docking "
            "predictions (let's call it n) in the file, followed by n rows. "
            "Each row reports the 4x4 rigid body tranformation matrix of a "
            "prediction. The rows are ordered in increasing order of rank, "
            "i.e. the first row is supposed to contain the best result."
         << endl;
    cout << "The capri score file will contain the irmsd, lrmsd, fnat, fnonnat "
            "and the final capri score for each row of the docking output file"
         << endl;
    cout << "The capri summary file will contain the number of high, medium, "
            "acceptable and incorrect results for different ranges of ranks"
         << endl;
    cout << "If the poseIndex (starts at 1 (top rank)) is given, then only "
            "that particular pose will be evaluated"
         << endl;

    return 1;
  }

  string receptor(argv[1]);
  string ligand(argv[2]);
  string output(argv[3]);
  string score(argv[4]);
  string summary(argv[5]);

  int index = 0;
  bool indexGiven = argc == 7;
  if (indexGiven)
    index = atoi(argv[6]) - 1;

  cout << "Given receptor file: " << receptor << endl;
  cout << "Given ligand file: " << ligand << endl;
  cout << "Given output file: " << output << endl;
  cout << "Score file: " << score << endl;
  cout << "Summary file: " << summary << endl;

  int verdictsInRange[8][4];
  char verdictNames[4][15] = {"Incorrect", "Acceptable", "Medium", "High"};

  FILE *fp = fopen(argv[3], "rt");

  if (!fp) {
    cout << "Could not open output file" << endl;
    return -1;
  }

  FILE *scoreFile = fopen(argv[4], "wt");
  if (!scoreFile) {
    cout << "Could not create score file" << endl;
    return -2;
  }

  fprintf(scoreFile, "%d \t %lf \t %lf \t %lf \t %lf \t %s\n", index, 0.0, 0.0,
          10.0, 10.0, verdictNames[0]);
  fclose(scoreFile);

  FILE *summaryFile = fopen(argv[5], "wt");
  if (!summaryFile) {
    cout << "Could not create summary file" << endl;
    fclose(scoreFile);
    return -3;
  }

  int numRes;
  fscanf(fp, "%d", &numRes);

  vector<Matrix> transformations;

  for (int i = 0; i < numRes; i++) {
    Matrix mtx;
    double temp;

    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        fscanf(fp, "%lf", &temp);
        mtx.set(j, k, temp);
      }
    }
    fscanf(fp, "%lf", &temp);
    fscanf(fp, "%lf", &temp);

    transformations.push_back(mtx);
  }

  fclose(fp);

  scoreFile = fopen(argv[4], "wt");
  if (!scoreFile) {
    cout << "Could not create score file" << endl;
    return -2;
  }

  CapriCriteria *capri = new CapriCriteria(receptor, ligand);

  int verdict;
  double fnat, fnonnat, irmsd, lrmsd;
  int rangeEndsAt[8] = {1, 10, 50, 100, 500, 1000, 2000, numRes};

  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 4; j++)
      verdictsInRange[i][j] = 0;

  //      fprintf(scoreFile, "Rank \t F_nat \t F_nonnat \t I_RMSD \t L_RMSD \t
  //      score\n");

  if (indexGiven && index >= numRes)
    fprintf(scoreFile, "%d \t %lf \t %lf \t %lf \t %lf \t %s\n", index, 0.0,
            0.0, 10.0, 10.0, verdictNames[0]);

  else {
    for (int i = 0; i < numRes; i++) {
      //             cout<<"\n\nCOMPUTING CAPRI SCORE FOR PREDICTION "<<i<<endl;

      if (indexGiven) {
        if (i == index)
          capri->computeCapriScore(transformations[i]);
        else
          continue;
      } else
        capri->computeCapriScore(transformations[i]);

      verdict = capri->getVerdict();
      fnat = capri->getFnat();
      fnonnat = capri->getFnonnat();
      irmsd = capri->getIRMSD();
      lrmsd = capri->getLRMSD();

      fprintf(scoreFile, "%d \t %lf \t %lf \t %lf \t %lf \t %s\n", i, fnat,
              fnonnat, irmsd, lrmsd, verdictNames[verdict]);

      for (int j = 0; j < 8; j++) {
        if (i < rangeEndsAt[j])
          verdictsInRange[j][verdict]++;
      }
    }
  }

  //      fprintf(summaryFile, "Range \t\t High \t Medium \t Acceptable \t
  //      Incorrect\n");

  for (int i = 0; i < 8; i++) {
    //              fprintf(summaryFile, "[1, %d] \t %d \t %d \t %d \t %d\n",
    //              rangeEndsAt[i], verdictsInRange[i][3],
    //              verdictsInRange[i][2], verdictsInRange[i][1],
    //              verdictsInRange[i][0]);
    fprintf(summaryFile, "%d \t %d \t %d \t %d\n", verdictsInRange[i][3],
            verdictsInRange[i][2], verdictsInRange[i][1],
            verdictsInRange[i][0]);
  }

  fclose(scoreFile);
  fclose(summaryFile);

  return 0;
}
