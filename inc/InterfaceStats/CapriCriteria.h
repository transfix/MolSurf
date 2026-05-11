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
/*
Given the correct docking configuration and a candidate docking result, this
class can compute the following evalutaion metric used in CAPRI-


1. f_nat = rescont_nat_pred / rescont_nat_correct
2. f_nonnat = rescont_nonnat_pred / rescont_nat_correct
        where,
        rescont_nat_pred = number of correct residue-residue contacts in the
candidate soln rescont_nonnat_pred = number of incorrect residue-residue
contacts in the candidate soln rescont_nat_correct = number of residue-residue
contacts in the correct soln note that, correct contacts are the ones found in
the correct solution.
3. I_RMSD. If IA is a the set of backbone atoms on the ligand which are on the
interface of the correct solution, then I_RMSD is the RMSD of the same set of
atoms in their predicted positions.
4. L_RMSD. If LA is a the set of all backbone atoms on the ligand, then L_RMSD
is the RMSD of the same set of atoms in their predicted positions.

Note that, F2Dock previously reported only I_RMSD.
*/

#ifndef CAPRI_CRITERIA_H
#define CAPRI_CRITERIA_H

#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>
#include <InterfaceStats/InterfaceStats.h>
#include <InterfaceStats/AtomsOnInterface.h>
#include <InterfaceStats/ComputeInterfaceStats.h>
#include <Blurmaps/InterfaceAtomDetector.h>
#include <PDBParser/Atom.h>
#include <PDBParser/FlattenGOA.h>
#include <DPG/PG.h>

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

class CapriCriteria {
private:
  string receptorFileName;
  string correctLigandFileName;
  string predictedLigandFileName;

  bool state; // false means error occured during earlier stages
  bool mode; // true = predicted ligand is pretransformed, false = prediction is
             // given as a xform
  Matrix mtx;
  Matrix identity;

  int numNat;
  int numCommon;
  int numPredMinusNat;

  double fnat;
  double fnonnat;
  double irmsd;
  double lrmsd;

  int verdict;

  PDBParser::GroupOfAtoms *receptor;
  PDBParser::GroupOfAtoms *correctLigand;
  PDBParser::GroupOfAtoms *predictedLigand;

  vector<PDBParser::Atom *> receptorAtoms;
  vector<PDBParser::Atom *> correctLigandAtoms;
  vector<PDBParser::Atom *> predictedLigandAtoms;

  PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;

  AtomsOnInterface *atomInterfaceCorrect;
  AtomsOnInterface *atomInterfacePredicted;

  vector<int> correctResiduePairs;
  vector<int> predictedResiduePairs;
  vector<int> interfaceAtoms;

  void init();

  bool getCorrectResiduePairs();
  bool getPredictedResiduePairs();
  bool computeIRMSD();
  bool computeLRMSD();
  bool computeNats();
  bool readAtomsFromFiles();

  bool isBackBoneAtom(PDBParser::Atom *atom);

public:
  CapriCriteria(string receptorFileName, string correctLigandFileName,
                string predictedLigandFileName);
  CapriCriteria(string receptorFileName, string correctLigandFileName,
                Matrix xform);
  CapriCriteria(string receptorFileName, string correctLigandFileName);
  ~CapriCriteria();

  bool computeCapriScore();
  bool computeCapriScore(Matrix xform);

  int getVerdict() { return verdict; }
  double getFnat() { return fnat; }
  double getFnonnat() { return fnonnat; }
  double getIRMSD() { return irmsd; }
  double getLRMSD() { return lrmsd; }
};

#endif
