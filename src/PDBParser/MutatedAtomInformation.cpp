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
#include <PDBParser/MutatedAtomInformation.h>

using PDBParser::MutatedAtomInformation;

MutatedAtomInformation::MutatedAtomInformation(float radius, const char* name, unsigned char residueID, float charge)
{
	m_Radius = radius;
	if(name)
	{
		strcpy(m_Name, name);
	}
	else
	{
		m_Name[0] = '\0';
	}
	m_ResidueId = residueID;
	m_Charge = charge;
	resetColor();
}

MutatedAtomInformation::MutatedAtomInformation(const MutatedAtomInformation& m)
{
	m_Radius = m.m_Radius;
	strcpy(m_Name, m.m_Name);
	m_ResidueId = m.m_ResidueId;
	m_Charge = m.m_Charge;
	resetColor();
}

MutatedAtomInformation::~MutatedAtomInformation()
{
}

unsigned char MutatedAtomInformation::getResidueId()
{
	return m_ResidueId;
}

void MutatedAtomInformation::resetColor()
{
	// search through element table and set color.
	// arand - commented to remove compiler warnings
	//if( m_ResidueId >= 0 )
	//{
	parserPDBtoGOA* p = new parserPDBtoGOA();
	int elementIndex = p->mapPQRElementAndResidue(m_Name, residueTable[m_ResidueId].residueName);
	delete p;
	if(elementIndex >=0 && elementIndex < PDBParser::MAX_BIOCHEM_ELEMENTS)
	{
		color[0] = PDBParser::elementTable[elementIndex].red;
		color[1] = PDBParser::elementTable[elementIndex].green;
		color[2] = PDBParser::elementTable[elementIndex].blue;
		return;
	}
	// arand commented
	//}
	// if not found, try this
	if(m_Name[1] == 'C')
	{
		color[0] = 0.3f;
		color[1] = 0.3f;
		color[2] = 0.3f;
		return;
	}
	if(m_Name[1] == 'O')
	{
		color[0] = 1.0f;
		color[1] = 0.0f;
		color[2] = 0.0f;
		return;
	}
	if(m_Name[1] == 'N')
	{
		color[0] = 0.0f;
		color[1] = 0.0f;
		color[2] = 1.0f;
		return;
	}
	if(m_Name[1] == 'P')
	{
		color[0] = 0.0f;
		color[1] = 1.0f;
		color[2] = 0.0f;
		return;
	}
	if(m_Name[1] == 'S')
	{
		color[0] = 1.0f;
		color[1] = 1.0f;
		color[2] = 0.0f;
		return;
	}
	if(m_Name[1] == 'H')
	{
		color[0] = 0.8f;
		color[1] = 0.8f;
		color[2] = 0.8f;
		return;
	}
	color[0] = PDBParser::defaultInformation.red;
	color[1] = PDBParser::defaultInformation.blue;
	color[2] = PDBParser::defaultInformation.green;
}

void MutatedAtomInformation::setResidueId(unsigned char residueId)
{
	m_ResidueId = residueId;
	resetColor();
}

// first look up using name and residue name and see if we get a match!
bool MutatedAtomInformation::getColor(float* red, float* green, float* blue)
{
	*red = color[0];
	*green  = color[1];
	*blue = color[2];
	return true;
}
