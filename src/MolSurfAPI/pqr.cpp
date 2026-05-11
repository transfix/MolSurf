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

#ifdef WITH_CG

#include <MolSurfAPI/pqr.h>
#include <string>
#include <iostream>
#include <algorithm>
MOLECULE::GroupOfAtoms *ParsePQR(FILE *fp, MOLECULE::RADIUS_TYPE radiusType) {
  MOLECULE::PQRParser *parser = new MOLECULE::PQRParser();
  MOLECULE::GroupOfAtoms *molecule = parser->ParserPQRtoGOA(fp, radiusType);
  delete parser;
  return molecule;
}

MOLECULE::GroupOfAtoms *ParsePDB(FILE *fp, MOLECULE::RADIUS_TYPE radiusType) {
  MOLECULE::PQRParser *parser = new MOLECULE::PQRParser();
  MOLECULE::GroupOfAtoms *molecule = parser->ParserPDBtoGOA(fp, radiusType);
  delete parser;
  return molecule;
}

MOLECULE::GroupOfAtoms *ParseOptimalPQR(FILE *fp,
                                        MOLECULE::RADIUS_TYPE radiusType) {
  MOLECULE::PQRParser *parser = new MOLECULE::PQRParser();
  MOLECULE::GroupOfAtoms *molecule =
      parser->ParserOptimalPQRtoGOA(fp, radiusType);
  delete parser;
  return molecule;
}

void append(MOLECULE::GroupOfAtoms *goaA, MOLECULE::GroupOfAtoms *goaB) {
  int resSeq0 = -1;
  for (int i = 0; i < goaB->numOfAtoms; i++) {
    MOLECULE::Atom *a = new MOLECULE::Atom();
    a = goaB->m_Atoms[i];
    goaA->addAtom(a);
    goaA->numOfAtoms++;
    if (a->resSeq != resSeq0) {
      goaA->numOfRes++;
      resSeq0 = a->resSeq;
    }
  }
}

void savePQR(MOLECULE::GroupOfAtoms *goa, FILE *fp) {
  for (int i = 0; i < goa->numOfAtoms; i++)
    fprintf(fp, "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f%8.4f%7.4f\n", i + 1,
            goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
            goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
            goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
            goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
            goa->m_Atoms[i]->radius);
}

void saveCDR(MOLECULE::GroupOfAtoms *goa, FILE *fp) {
  FILE *fpH = fopen("OutH.txt", "w");
  FILE *fpO = fopen("OutO.txt", "w");
  FILE *fpC = fopen("OutC.txt", "w");
  FILE *fpN = fopen("OutN.txt", "w");

  char temp[10];
  for (int i = 0; i < goa->numOfAtoms; i++) {
    switch ((int)goa->m_Atoms[i]->occupancy) {
    case 0:
      strcpy(temp, "0xff0000");
      fprintf(fpC,
              "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f %8.4f %8.4f %8s\n",
              i + 1, goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
              goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
              goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
              goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
              goa->m_Atoms[i]->radius, temp);

      break;
    case 1:
      strcpy(temp, "0x00ff00");
      fprintf(fpO,
              "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f %8.4f %8.4f %8s\n",
              i + 1, goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
              goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
              goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
              goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
              goa->m_Atoms[i]->radius, temp);

      break;
    case 2:
      strcpy(temp, "0x0000ff");
      fprintf(fpN,
              "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f %8.4f %8.4f %8s\n",
              i + 1, goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
              goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
              goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
              goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
              goa->m_Atoms[i]->radius, temp);

      break;
    case 3:
      strcpy(temp, "0x000000");
      fprintf(fpH,
              "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f %8.4f %8.4f %8s\n",
              i + 1, goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
              goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
              goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
              goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
              goa->m_Atoms[i]->radius, temp);

      break;
    case 4:
      strcpy(temp, "0x00ff00");
      break;
    }

    fprintf(fp, "ATOM%7d %-3s %3s %c%4d    %8.3f%8.3f%8.3f %8.4f %8.4f %8s\n",
            i + 1, goa->m_Atoms[i]->atomName, goa->m_Atoms[i]->resName,
            goa->m_Atoms[i]->chainID, goa->m_Atoms[i]->resSeq,
            goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
            goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->charge,
            goa->m_Atoms[i]->radius, temp);
  }
  fclose(fpC);
  fclose(fpH);
  fclose(fpO);
  fclose(fpN);
}

void savePDB(MOLECULE::GroupOfAtoms *goa, FILE *fp, bool withH,
             bool withChain) {
  char setChain = 'A';
  char setChain0 = 'A';
  char atomType;
  int resIndex = 0;
  int resTemp = 0;
  int flag = 0;

  int j = 1;
  for (int i = 0; i < goa->numOfAtoms; i++) {
    for (int t = 0; t < 10; t++) {
      int ascii = (int)(goa->m_Atoms[i]->atomName[t]);
      if (ascii <= 90 && ascii >= 65) {
        atomType = goa->m_Atoms[i]->atomName[t];
        break;
      }
    }
    string atomname = string(goa->m_Atoms[i]->atomName);
    atomname.erase(std::remove(atomname.begin(), atomname.end(), ' '),
                   atomname.end());

    if (goa->m_Atoms[i]->chainID != ' ') {
      if (setChain != goa->m_Atoms[i]->chainID) {
        flag = 0;
        if (goa->m_Atoms[i]->resSeq > resIndex) {
          flag = 1;
          resIndex = goa->m_Atoms[i]->resSeq - 1;
        }
        setChain = goa->m_Atoms[i]->chainID;
      }
      flag = 0; // problem here.
      if (flag == 1)
        resTemp = goa->m_Atoms[i]->resSeq - resIndex;
      else
        resTemp = goa->m_Atoms[i]->resSeq;
      //             cout<<"resTemp: "<< goa->m_Atoms[i]->resSeq<< " " <<
      //             resTemp<< endl;

      if (withChain == 0) {
        setChain = ' ';
        goa->m_Atoms[i]->chainID = ' ';
      }
      if (strstr(goa->m_Atoms[i]->atomName, "OXT") == NULL) {

        if (withH == 0 && atomType == 'H')
          continue;
        else {
          fprintf(fp,
                  "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f         "
                  "  %c\n",
                  j,
                  //   goa->m_Atoms[i]->atomName,
                  atomname.c_str(), goa->m_Atoms[i]->resName,
                  goa->m_Atoms[i]->chainID, resTemp,
                  goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
                  goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->occupancy,
                  goa->m_Atoms[i]->tempFactor, atomType);
          j++;
        }
      } else {
        if (withChain == 0) {
          setChain = ' ';
          goa->m_Atoms[i]->chainID = ' ';
        }

        while (i < goa->numOfAtoms - 1 &&
               (goa->m_Atoms[i]->resSeq == goa->m_Atoms[i + 1]->resSeq)) {

          if (withH == 0 && atomType == 'H') {
            i++;
            for (int t = 0; t < 10; t++) {
              int ascii = (int)(goa->m_Atoms[i]->atomName[t]);
              if (ascii <= 90 && ascii >= 65) {
                atomType = goa->m_Atoms[i]->atomName[t];
                break;
              }
            }
            continue;
          } else {

            string atomname1 = string(goa->m_Atoms[i]->atomName);
            atomname1.erase(
                std::remove(atomname1.begin(), atomname1.end(), ' '),
                atomname1.end());

            fprintf(
                fp,
                "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           "
                "%c\n",
                j,
                //                                      goa->m_Atoms[i]->atomName,
                atomname1.c_str(), goa->m_Atoms[i]->resName, setChain, resTemp,
                goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
                goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->occupancy,
                goa->m_Atoms[i]->tempFactor, atomType);
            i++;
            j++;
            for (int t = 0; t < 10; t++) {
              int ascii = (int)(goa->m_Atoms[i]->atomName[t]);
              if (ascii <= 90 && ascii >= 65) {
                atomType = goa->m_Atoms[i]->atomName[t];
                break;
              }
            }
          }
        }
        if (withChain == 0) {
          setChain = ' ';
          goa->m_Atoms[i]->chainID = ' ';
        }

        /*                    fprintf(fp, "ATOM%7d %-4s%3s %c%4d
      %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n", j,
      //                              goa->m_Atoms[i]->atomName,
                                  atomname.c_str(),
                                      goa->m_Atoms[i]->resName,
                                      goa->m_Atoms[i]->chainID,
                                      resTemp,
                                      goa->m_Atoms[i]->position[0],
                                      goa->m_Atoms[i]->position[1],
                                      goa->m_Atoms[i]->position[2],
                                      goa->m_Atoms[i]->occupancy,
                                      goa->m_Atoms[i]->tempFactor,
                                      atomType);
                                      */
        fprintf(fp, "TER%8d      %3s %c%4d            \n", j,
                goa->m_Atoms[i]->resName, goa->m_Atoms[i]->chainID,
                //                              goa->m_Atoms[i]->resSeq);
                resTemp);
        j += 1;
        //                     resIndex = goa->m_Atoms[i]->resSeq;
      }
    } else {

      if (setChain0 != setChain) {

        flag = 0;
        if (goa->m_Atoms[i]->resSeq > resIndex) {
          flag = 1;
          resIndex = goa->m_Atoms[i]->resSeq - 1;
        }
        setChain0 = setChain;
      }
      if (flag == 1)
        resTemp = goa->m_Atoms[i]->resSeq - resIndex;
      else
        resTemp = goa->m_Atoms[i]->resSeq;

      if (withChain == 0) {
        setChain = ' ';
        goa->m_Atoms[i]->chainID = ' ';
      }

      if (strstr(goa->m_Atoms[i]->atomName, "OXT") == NULL) {
        if (withH == 0 && atomType == 'H')
          continue;
        else {

          fprintf(fp,
                  "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f         "
                  "  %c\n",
                  j,
                  //                              goa->m_Atoms[i]->atomName,
                  atomname.c_str(), goa->m_Atoms[i]->resName, setChain, resTemp,
                  goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
                  goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->occupancy,
                  goa->m_Atoms[i]->tempFactor, atomType);
          j++;
        }
      } else {

        fprintf(fp,
                "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           "
                "%c\n",
                j,
                //                              goa->m_Atoms[i]->atomName,
                atomname.c_str(), goa->m_Atoms[i]->resName, setChain, resTemp,
                goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
                goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->occupancy,
                goa->m_Atoms[i]->tempFactor, atomType);
        while (i < goa->numOfAtoms - 1 &&
               (goa->m_Atoms[i]->resSeq == goa->m_Atoms[i + 1]->resSeq)) {

          if (withH == 0 && atomType == 'H') {
            i++;
            for (int t = 0; t < 10; t++) {
              int ascii = (int)(goa->m_Atoms[i]->atomName[t]);
              if (ascii <= 90 && ascii >= 65) {
                atomType = goa->m_Atoms[i]->atomName[t];
                break;
              }
            }
            continue;
          } else {

            string atomname2 = string(goa->m_Atoms[i]->atomName);
            atomname2.erase(
                std::remove(atomname2.begin(), atomname2.end(), ' '),
                atomname2.end());

            fprintf(
                fp,
                "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           "
                "%c\n",
                j,
                //                                      goa->m_Atoms[i]->atomName,
                atomname2.c_str(), goa->m_Atoms[i]->resName, setChain, resTemp,
                goa->m_Atoms[i]->position[0], goa->m_Atoms[i]->position[1],
                goa->m_Atoms[i]->position[2], goa->m_Atoms[i]->occupancy,
                goa->m_Atoms[i]->tempFactor, atomType);
            i++;
            j++;

            for (int t = 0; t < 10; t++) {
              int ascii = (int)(goa->m_Atoms[i]->atomName[t]);
              if (ascii <= 90 && ascii >= 65) {
                atomType = goa->m_Atoms[i]->atomName[t];
                break;
              }
            }
          }
        }
        if (withChain == 0) {
          setChain = ' ';
          goa->m_Atoms[i]->chainID = ' ';
        }

        fprintf(fp, "TER%8d      %3s %c%4d            \n", j,
                goa->m_Atoms[i]->resName, setChain,
                (goa->m_Atoms[i]->resSeq - resIndex));

        j += 1;
        setChain += 1;
        //                  resIndex = goa->m_Atoms[i]->resSeq;
      }
    }
  }
}

#endif // WITH_CG
