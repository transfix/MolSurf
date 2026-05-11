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
#include <PDBParser/FlattenGOA.h>

//
// arand, 11-19-2010: file created
// this code has been moved from the BlurMapsDataManager
// so that the GOALoader doesn't have "upward" dependencies
//
namespace PDBParser {

bool FlattenGOA(GroupOfAtoms *grp, vector<Atom *> &flatList,
                CollectionData *collectionData, int numberOfTransformations,
                double *rotations, double *translations,
                GroupOfAtoms::RADIUS_TYPE radiusType, unsigned int level,
                bool colorBySubunits) {

  if (grp->type == CHAIN_TYPE) {
    char cID =
        (char)grp->resID; // resID is overloaded in GOA to be chainID for chains
    if (collectionData) {
      if (collectionData->isTransformationPresent(cID)) {
        int nTrans = collectionData->getNumberOfTransformations(cID);
        if (nTrans > 0) {
          numberOfTransformations = nTrans;
          rotations = new double[numberOfTransformations * 9];
          translations = new double[numberOfTransformations * 3];
          if (!collectionData->getRotations(cID, rotations) ||
              !collectionData->getTranslations(cID, translations)) {
            numberOfTransformations = 0;
            delete[] rotations;
            rotations = 0;
            delete[] translations;
            translations = 0;
          }
        }
      }
    }
  }

  if (grp->type == COLLECTION_TYPE) {
    assert(!((grp->m_SubGroups).empty()));
    FlattenGOA(grp->m_SubGroups.at(0), flatList, collectionData,
               numberOfTransformations, rotations, translations, radiusType,
               level, colorBySubunits);
  } else {
    for (int i = 0; i < grp->m_SubGroups.size(); i++) {
      FlattenGOA(grp->m_SubGroups.at(i), flatList, collectionData,
                 numberOfTransformations, rotations, translations, radiusType,
                 level, colorBySubunits);
    }
  }

  int m = grp->m_Atoms.size();
  for (int i = 0; i < m; i++) {

    Atom *at = grp->m_Atoms.at(i);
    float x = at->m_Position[0];
    float y = at->m_Position[1];
    float z = at->m_Position[2];
    if (numberOfTransformations > 0 && rotations && translations) {
      for (int j = 0; j < numberOfTransformations; j++) {
        double rot[9];
        rot[0] = rotations[j * 9 + 0];
        rot[1] = rotations[j * 9 + 1];
        rot[2] = rotations[j * 9 + 2];
        rot[3] = rotations[j * 9 + 3];
        rot[4] = rotations[j * 9 + 4];
        rot[5] = rotations[j * 9 + 5];
        rot[6] = rotations[j * 9 + 6];
        rot[7] = rotations[j * 9 + 7];
        rot[8] = rotations[j * 9 + 8];
        float x1 = x * rot[0] + y * rot[1] + z * rot[2];
        float y1 = x * rot[3] + y * rot[4] + z * rot[5];
        float z1 = x * rot[6] + y * rot[7] + z * rot[8];
        double trans[3];
        trans[0] = translations[j * 3 + 0];
        trans[1] = translations[j * 3 + 1];
        trans[2] = translations[j * 3 + 2];
        x1 += trans[0];
        y1 += trans[1];
        z1 += trans[2];
        Atom *dupAtom = new Atom(*at);
        dupAtom->m_Position[0] = x1;
        dupAtom->m_Position[1] = y1;
        dupAtom->m_Position[2] = z1;
        {
          // make only one subunit colored!
          if (colorBySubunits) {
            srand(j * 7);
            double t1 = rand();
            double t2 = rand();
            unsigned int r, g, b, a;
            double red = rand() / ((double)(RAND_MAX));
            double green = rand() / ((double)(RAND_MAX));
            double blue = rand() / ((double)(RAND_MAX));
            double alpha = 1.0;
            r = (unsigned int)(red * 255.0);
            g = (unsigned int)(green * 255.0);
            b = (unsigned int)(blue * 255.0);
            a = (unsigned int)(alpha * 255.0);
            dupAtom->m_UserDefinedColor = ((r << 24) & 0xff000000) |
                                          ((g << 16) & 0xff0000) |
                                          ((b << 8) & 0xff00) | (a & 0xff);
          }
        }
        flatList.push_back(dupAtom);
      }
    } else {
      Atom *dupAtom = new Atom(*at);
      // arand: constructor does this?
      // dupAtom->m_Position[0] = at->m_Position[0];
      //  dupAtom->m_Position[1] = at->m_Position[1];
      //  dupAtom->m_Position[2] = at->m_Position[2];
      flatList.push_back(dupAtom);
    }
  }
  return true;
}

} // namespace PDBParser
