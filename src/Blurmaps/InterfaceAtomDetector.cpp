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

#include <Blurmaps/InterfaceAtomDetector.h>
#include <PDBParser/FlattenGOA.h>
#include <DPG/PG.h>
#include <string>
#include <sstream>

InterfaceAtomDetector::InterfaceAtomDetector(
    PDBParser::GroupOfAtoms *molecule1, PDBParser::GroupOfAtoms *molecule2,
    string ofn, double iw, string at,
    PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType) {
  PDBParser::CollectionData *collectionData1 = 0;
  if (molecule1->type == PDBParser::COLLECTION_TYPE) {
    collectionData1 = molecule1->m_CollectionData;
  }
  FlattenGOA(molecule1, m_AtomList1, collectionData1, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);

  PDBParser::CollectionData *collectionData2 = 0;
  if (molecule2->type == PDBParser::COLLECTION_TYPE) {
    collectionData2 = molecule2->m_CollectionData;
  }
  FlattenGOA(molecule2, m_AtomList2, collectionData2, 0, 0, 0, radiusType,
             PDBParser::ATOM_TYPE, false);

  interfaceWidth = iw;
  outputFileName = ofn;
  atomTypes = at;
  selectedIds.clear();
}

InterfaceAtomDetector::~InterfaceAtomDetector() {
  m_AtomList1.clear();
  m_AtomList2.clear();
  selectedIds.clear();
}

bool InterfaceAtomDetector::detect() {
  FILE *outputFile, *xyzrFile;

  string xyzrFileName;
  xyzrFileName.assign(outputFileName, 0, outputFileName.length() - 3);
  xyzrFileName += "xyzr";

  outputFile = fopen(outputFileName.c_str(), "wt");

  if (!outputFile) {
    cout << "could not open output file" << endl;
    return false;
  }

  else {
    xyzrFile = fopen(xyzrFileName.c_str(), "wt");

    if (!xyzrFile) {
      cout << "could not create xyzr file" << endl;
      fclose(outputFile);
      return false;
    }

    // Deciding which types of atoms to print

    bool printAll = false;
    bool printBackBone = false;
    bool printC = false;
    bool printCA = false;
    bool printN = false;
    bool printO = false;

    string C(" C ");
    string CA(" CA ");

    if (atomTypes.compare("all") == 0) {
      cout << "Printing all" << endl;
      printAll = true;
    } else if (atomTypes.compare("backbone") == 0) {
      cout << "Printing backbone only" << endl;
      printBackBone = true;
    } else {
      int startPos = 0;
      int endPos = atomTypes.find(",");

      while (1) {
        //                              cout<<startPos<<" "<<endPos<<endl;

        int length;

        if (endPos == string::npos)
          length = atomTypes.length() - startPos;
        else
          length = endPos - startPos;

        string sub = atomTypes.substr(startPos, length);

        //                              cout<<sub<<endl;

        if (sub.compare("C") == 0) {
          cout << "Printing C only" << endl;
          printC = true;
        } else if (sub.compare("CA") == 0) {
          cout << "Printing CA only" << endl;
          printCA = true;
        }

        else if (sub.compare("O") == 0) {
          cout << "Printing O only" << endl;
          printO = true;
        }

        else if (sub.compare("N") == 0) {
          cout << "Printing N only" << endl;
          printN = true;
        }

        if (endPos == string::npos)
          break;

        startPos = endPos + 1;
        endPos = atomTypes.find(",", startPos);
      }
    }

    // inserting receptor atoms into dpg

    int numAtoms1 = m_AtomList1.size();
    DPG::PG *pg = new DPG::PG(10.0, 1000.0, 3.0);

    for (int i = 0; i < numAtoms1; i++) {
      DPG::Point *p = new DPG::Point(m_AtomList1[i]->m_Position[0],
                                     m_AtomList1[i]->m_Position[1],
                                     m_AtomList1[i]->m_Position[2]);
      pg->addPoint(p);
    }

    // Identifying interface atoms of ligand and printing to the output

    vector<DPG::Point *> results;

    int numAtoms2 = m_AtomList2.size();

    int yesCount = 0;
    int noCount = 0;
    int CCount = 0;
    int CACount = 0;
    int NCount = 0;
    int OCount = 0;

    for (int i = 0; i < numAtoms2; i++) {
      PDBParser::Atom *atom = m_AtomList2[i];

      bool isC = false;
      bool isCA = false;
      bool isO = false;
      bool isN = false;

      string name(atom->getName());

      if (name.compare(" N  ") == 0) {
        isN = true;
        NCount++;
      } else if (name.compare(" O  ") == 0) {
        isO = true;
        OCount++;
      } else if (name.compare(" C  ") == 0) {
        isC = true;
        CCount++;
      } else if (name.compare(" CA ") == 0) {
        isCA = true;
        CACount++;
      }

      bool doPrint = false;

      if (printAll)
        doPrint = true;
      else {
        if (printBackBone && (isC || isCA || isN))
          doPrint = true;

        else if (printC && isC)
          doPrint = true;
        else if (printCA && isCA)
          doPrint = true;
        else if (printO && isO)
          doPrint = true;
        else if (printN && isN)
          doPrint = true;
      }

      if (!doPrint)
        continue;

      DPG::Point *p = new DPG::Point(m_AtomList2[i]->m_Position[0],
                                     m_AtomList2[i]->m_Position[1],
                                     m_AtomList2[i]->m_Position[2]);

      results.clear();
      results = pg->range(p, interfaceWidth);

      if (results.size() > 0) {
        yesCount++;
        selectedIds.push_back(i);
      }
    }

    int interfaceSize = selectedIds.size();

    fprintf(outputFile, "%d\n", numAtoms2);
    fprintf(outputFile, "%d\n", interfaceSize);

    set<string> ResidueIndex;
    map<char, int> ChainIDNamemap;
    int chainIndex = 0;

    set<string>::iterator it;
    vector<string> ResidueIds[20];
    pair<set<string>::iterator, bool> ret;

    int stat[20];

    // 0 -- non-polar and neutral
    // 1 -- polar and neutral
    // 2 -- acidic and polar
    // 3 -- basic and polar

    int charges[4];

    for (int i = 0; i < 20; i++)
      stat[i] = 0;
    for (int i = 0; i < 4; i++)
      charges[i] = 0;

    ofstream myfile;
    myfile.open(strcat((char *)outputFileName.c_str(), ".cmap"));

    myfile << "# cmap file " << endl;
    myfile << "# \"moleID\" \"chainID\" \"chainName\" \"ssID\" \"residueID\" "
              "\"residueName\" \"atomName\" \"atomID\" \"r\" \"g\" \"b\" "
              "\"alpha\" "
           << endl;

    for (int i = 0; i < interfaceSize; i++) {
      int index = selectedIds[i];
      PDBParser::Atom *atom = m_AtomList2[index];
      int temp = atom->getResidueNum();
      char cha = atom->getChain();

      stringstream out;

      if (cha != ' ') {
        out << cha << temp;
        if (ChainIDNamemap.find(cha) == ChainIDNamemap.end()) {
          ChainIDNamemap[cha] = chainIndex;
          chainIndex++;
        }
      } else
        out << temp;
      string chainres = out.str();
      ret = ResidueIndex.insert(chainres);

      if (ret.second == true) {
        string s = atom->getResidueName();
        myfile << " \"0\" \"" << ChainIDNamemap[cha] << "\" \"" << cha
               << "\" \"0\" \"" << temp << "\" \"" << s << "\""
               << " \"\" \"\" \"0\" \"0\" \"1\" \"1.0\" " << endl;
        //                       myfile<<" \"0\" \"" <<ChainIDNamemap[cha] <<"\"
        //                       \"" << cha <<"\" \"0\" \"" <<temp << "\" \""
        //                       <<s <<"\""
        //                     <<" \"\" \"\" \"1\" \"0\" \"0\" \"1.0\" "<<endl;

        /*
                                                if(strncmp(s.c_str(),"GLY",3)==0)
                                                {
                                                        stat[0]++;
                                                        charges[0]++;
                                                        cout<<"chainres: " <<
           chainres <<endl; ResidueIds[0].push_back(chainres); myfile<<" \"0\"
           \"" <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp <<
           "\" \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; continue;
                                                }

                                                if(strncmp(s.c_str(),"ALA",3)==0)
                                                {
                                                        stat[1]++;
                                                        charges[0]++;
                                                        ResidueIds[1].push_back(chainres);
                                                   myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                                        <<" \"\" \"\" \"0.5\"
           \"0.5\" \"0.5\" \"1.0\" "<<endl; continue;
                                                }
                                                if(strncmp(s.c_str(),"SER",3)==0)
                                                {
                                                        stat[2]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[2].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"THR",3)==0)
                                                {
                                                        stat[3]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[3].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"CYS",3)==0)
                                                {
                                                        stat[4]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[4].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"VAL",3)==0)
                                                {
                                                        stat[5]++;
                                                        charges[0]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[5].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"ILE",3)==0)
                                                {
                                                        stat[6]++;
                                                        charges[0]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[6].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"LEU",3)==0)
                                                {
                                                        stat[7]++;
                                                        charges[0]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[7].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"PRO",3)==0)
                                                {
                                                        stat[8]++;
                                                        charges[0]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[8].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"PHE",3)==0)
                                                {
                                                        stat[9]++;
                                                        charges[0]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[9].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"TYR",3)==0)
                                                {
                                                        stat[10]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[10].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"MET",3)==0)
                                                {
                                                        stat[11]++;
                                                        ResidueIds[11].push_back(chainres);
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.5\" \"0.5\" \"0.5\"
           \"1.0\" "<<endl; charges[0]++; continue;
                                                }
                                                if(strncmp(s.c_str(),"TRP",3)==0)
                                                {
                                                        stat[12]++;
                                                        ResidueIds[12].push_back(chainres);
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; continue;
                                                }
                                                if(strncmp(s.c_str(),"ASN",3)==0)
                                                {
                                                        stat[13]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[13].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"GLN",3)==0)
                                                {
                                                        stat[14]++;
                                                        charges[1]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"1\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[14].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"HIS",3)==0)
                                                {
                                                        stat[15]++;
                                                        charges[3]++;
                                                        ResidueIds[15].push_back(chainres);
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.25\" \"0.25\" \"0.65\"
           \"1.0\" "<<endl; continue;
                                                }
                                                if(strncmp(s.c_str(),"ASP",3)==0)
                                                {
                                                        stat[16]++;
                                                        charges[2]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"0\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[16].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"GLU",3)==0)
                                                {
                                                        stat[17]++;
                                                        charges[2]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"1\" \"0\" \"0.5\"
           \"1.0\" "<<endl; ResidueIds[17].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"LYS",3)==0)
                                                {
                                                        stat[18]++;
                                                        charges[3]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.25\" \"0.25\" \"0.65\"
           \"1.0\" "<<endl; ResidueIds[18].push_back(chainres); continue;
                                                }
                                                if(strncmp(s.c_str(),"ARG",3)==0)
                                                {
                                                        stat[19]++;
                                                        charges[3]++;
                                                    myfile<<" \"0\" \""
           <<ChainIDNamemap[cha] <<"\" \"" << cha <<"\" \"0\" \"" <<temp << "\"
           \"" <<s <<"\""
                                        <<" \"\" \"\" \"0.25\" \"0.25\" \"0.65\"
           \"1.0\" "<<endl; ResidueIds[19].push_back(chainres); continue;
                                                }
                                                */
      }
      // cout << "chain: "<< cha <<"  "<< (atom->getResidueNum()) << " " << s <<
      // endl;

      fprintf(outputFile, "%d   %5.3lf   %5.3lf   %5.3lf\n", index + 1,
              atom->m_Position[0], atom->m_Position[1], atom->m_Position[2]);
      fprintf(xyzrFile, "%5.3lf   %5.3lf   %5.3lf   %5.3lf\n",
              atom->m_Position[0], atom->m_Position[1], atom->m_Position[2],
              atom->getRadius());
    }

    for (int i = 0; i < numAtoms2; i++) {
      PDBParser::Atom *atom = m_AtomList2[i];
      int temp = atom->getResidueNum();
      char cha = atom->getChain();

      stringstream out;

      if (cha != ' ')
        out << cha << temp;
      else
        out << temp;
      string chainres = out.str();
      if (ResidueIndex.find(chainres) == ResidueIndex.end()) {
        ResidueIndex.insert(chainres);
        myfile << " \"0\" \"" << ChainIDNamemap[cha] << "\" \"" << cha
               << "\" \"0\" \"" << temp << "\" \"" << atom->getResidueName()
               << "\"" << " \"\" \"\" \"1\" \"0.65\" \"0.85\" \"1.0\" " << endl;
        //                    myfile<<" \"0\" \"" <<ChainIDNamemap[cha] <<"\"
        //                    \"" << cha <<"\" \"0\" \"" <<temp << "\" \"" <<
        //                    atom->getResidueName() <<"\""
        //            <<" \"\" \"\" \"0\" \"1\" \"1\" \"1.0\" "<<endl;
      }
    }
    myfile.close();

    cout << "Total: " << ResidueIndex.size() << " residues." << endl;

    cout << "Neutral nonpolar: " << charges[0] << ", of which: " << endl;
    cout << "  ALA: " << stat[1] << ", they are: ";
    for (int i = 0; i < ResidueIds[1].size(); i++)
      cout << ResidueIds[1][i] << " ";
    cout << endl;
    cout << "  GLY: " << stat[0] << ", they are: ";
    for (int i = 0; i < ResidueIds[0].size(); i++)
      cout << ResidueIds[0][i] << " ";
    cout << endl;
    cout << "  ILE: " << stat[6] << ", they are: ";
    for (int i = 0; i < ResidueIds[6].size(); i++)
      cout << ResidueIds[6][i] << " ";
    cout << endl;
    cout << "  LEU: " << stat[7] << ", they are: ";
    for (int i = 0; i < ResidueIds[7].size(); i++)
      cout << ResidueIds[7][i] << " ";
    cout << endl;
    cout << "  MET: " << stat[11] << ", they are: ";
    for (int i = 0; i < ResidueIds[11].size(); i++)
      cout << ResidueIds[11][i] << " ";
    cout << endl;
    cout << "  PHE: " << stat[9] << ", they are: ";
    for (int i = 0; i < ResidueIds[9].size(); i++)
      cout << ResidueIds[9][i] << " ";
    cout << endl;
    cout << "  PRO: " << stat[8] << ", they are: ";
    for (int i = 0; i < ResidueIds[8].size(); i++)
      cout << ResidueIds[8][i] << " ";
    cout << endl;
    cout << "  VAL: " << stat[5] << ", they are: ";
    for (int i = 0; i < ResidueIds[5].size(); i++)
      cout << ResidueIds[5][i] << " ";
    cout << endl;

    cout << "Neutral polar: " << charges[1] << ", of which: " << endl;
    cout << "  ASN: " << stat[13] << ", they are: ";
    for (int i = 0; i < ResidueIds[13].size(); i++)
      cout << ResidueIds[13][i] << " ";
    cout << endl;
    cout << "  CYS: " << stat[4] << ", they are: ";
    for (int i = 0; i < ResidueIds[4].size(); i++)
      cout << ResidueIds[4][i] << " ";
    cout << endl;
    cout << "  GLN: " << stat[14] << ", they are: ";
    for (int i = 0; i < ResidueIds[14].size(); i++)
      cout << ResidueIds[14][i] << " ";
    cout << endl;
    cout << "  SER: " << stat[2] << ", they are: ";
    for (int i = 0; i < ResidueIds[2].size(); i++)
      cout << ResidueIds[2][i] << " ";
    cout << endl;
    cout << "  TYR: " << stat[10] << ", they are: ";
    for (int i = 0; i < ResidueIds[10].size(); i++)
      cout << ResidueIds[10][i] << " ";
    cout << endl;
    cout << "  THR: " << stat[3] << ", they are: ";
    for (int i = 0; i < ResidueIds[3].size(); i++)
      cout << ResidueIds[3][i] << " ";
    cout << endl;
    cout << "  TRP: " << stat[12] << ", they are: ";
    for (int i = 0; i < ResidueIds[12].size(); i++)
      cout << ResidueIds[12][i] << " ";
    cout << endl;

    cout << "Acidic polar: " << charges[2] << ", of which: " << endl;
    cout << "  ASP: " << stat[16] << ", they are: ";
    for (int i = 0; i < ResidueIds[16].size(); i++)
      cout << ResidueIds[16][i] << " ";
    cout << endl;
    cout << "  GLU: " << stat[17] << ", they are: ";
    for (int i = 0; i < ResidueIds[17].size(); i++)
      cout << ResidueIds[17][i] << " ";
    cout << endl;

    cout << "Basic polar: " << charges[3] << ", of which: " << endl;
    cout << "  ARG: " << stat[19] << ", they are: ";
    for (int i = 0; i < ResidueIds[19].size(); i++)
      cout << ResidueIds[19][i] << " ";
    cout << endl;
    cout << "  HIS: " << stat[15] << ", they are: ";
    for (int i = 0; i < ResidueIds[15].size(); i++)
      cout << ResidueIds[15][i] << " ";
    cout << endl;
    cout << "  LYS: " << stat[18] << ", they are: ";
    for (int i = 0; i < ResidueIds[18].size(); i++)
      cout << ResidueIds[18][i] << " ";
    cout << endl;

    cout << "Printed " << yesCount << " interface atoms out of " << numAtoms2
         << " total atoms" << endl;

    //              cout<<"Oxygens = "<<OCount<<endl;
    //              cout<<"Carbons = "<<CCount<<endl;
    //              cout<<"Alpha Carbons = "<<CACount<<endl;
    //              cout<<"Nitrogens = "<<NCount<<endl;

    fclose(outputFile);
    fclose(xyzrFile);
    return true;
  }
}
