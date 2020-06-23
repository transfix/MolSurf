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
#if !defined(AFX_RESIDUEDATA_H__380D8ECF_6A09_4DB9_9D91_6C93ACAED005__INCLUDED_)
#define AFX_RESIDUEDATA_H__380D8ECF_6A09_4DB9_9D91_6C93ACAED005__INCLUDED_

#include <Utility/utility.h>

namespace PDBParser
{
	class ResidueData
	{
		public:
			ResidueData();
			ResidueData(const ResidueData& residueData);
			virtual ~ResidueData();
			// Torsion angles for nucleic acids (i.e. RNA)
			float m_alpha, m_beta, m_gamma, m_delta, m_epsilon, m_zeta, m_nu0, m_nu1, m_nu2, m_nu3, m_nu4, m_chi;
			// Torsion angles for amino acids
			float m_Phi, m_Psi;
			// Reduced torsion angles for nucleic acids
			float m_Eta, m_Theta;
	};
};

#endif
