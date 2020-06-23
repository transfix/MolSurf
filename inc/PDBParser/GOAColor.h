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
#ifndef GOA_COLOR_H
#define GOA_COLOR_H

#include <Utility/utility.h>
#include <PDBParser/GroupOfAtoms.h>
#include <PDBParser/Atom.h>
#include <PDBParser/ResidueInformation.h>
#include <ColorManager/PastelColorGenerator.h>

namespace ColorManager
{
	class ColorGenerator;
}

namespace PDBParser
{
	class GroupOfAtoms;
	class Atom;
	const int NumTokens = 12;
	class GOAColor
	{
		public:
			GOAColor();
			~GOAColor();
			void parseColormapFile(char*);
			void parseColormapString(vector<std::string> strings);
			void printColormap() const;
			void ApplyColormap(PDBParser::GroupOfAtoms*, int transIndex = 0) const;
			void ApplyColorsByGOALevel(PDBParser::GroupOfAtoms* molecule, PDBParser::GroupOfAtoms::GOA_TYPE colorLevel);
			void ApplyElementColors(PDBParser::GroupOfAtoms*) const;
			void ApplyResidueColors(PDBParser::GroupOfAtoms*) const;
			void ApplyChainColors(PDBParser::GroupOfAtoms* grp, ColorManager::ColorGenerator* colorGenerator) const;
			bool UsesDefaultColoring(int transIndex) const;
			struct ColormapEntry
			{
				int transformSection;
				int molecule;
				int chainIndex;
				char chainID;
				int ssIndex;
				int resIndex;
				char resName[4];
				char atomName[5];
				int atomIndex;
				unsigned int lineNum;
				double red;
				double green;
				double blue;
				double alpha;
			};

			void ApplyColormapEntryToGOA(const ColormapEntry&,PDBParser::GroupOfAtoms*) const;
			void ApplyColormapEntryToGOA_old(const ColormapEntry&,PDBParser::GroupOfAtoms*) const;
			void ApplyColorToGOA(double,double,double,double, PDBParser::GroupOfAtoms*) const;
			void ApplyColorToAtom(float,float,float,float, PDBParser::Atom*) const;

		private:
			PDBParser::Atom* getAnAtom(PDBParser::GroupOfAtoms* grp) const;
			void printEntry(const ColormapEntry&) const;
			char* eatWhitespace(char*);
			char** parseEntryLine(char*);
			int parseSectionMarker(char*);
			bool genEntryFromTokens(char**, ColormapEntry*);
			void freeTokens(char**);
			enum TokenFields
			{
				MolIndex=0, ChainIndex,ChainID, SSIndex, ResidueIndex,ResidueName,
				AtomName,AtomIndex, Red,Green,Blue,Alpha
			};
			vector<ColormapEntry> m_EntList;
	};
};

#endif
