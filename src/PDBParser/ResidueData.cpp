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
#include <PDBParser/ResidueData.h>

using PDBParser::ResidueData;

ResidueData::ResidueData()
{
	m_alpha = 0;
	m_beta = 0;
	m_gamma = 0;
	m_delta = 0;
	m_epsilon = 0;
	m_zeta = 0;
	m_nu0 = 0;
	m_nu1 = 0;
	m_nu2 = 0;
	m_nu3 = 0;
	m_nu4 = 0;
	m_chi = 0;
	m_Phi = 0;
	m_Psi = 0;
	m_Eta = 0;
	m_Theta = 0;
}

ResidueData::ResidueData(const ResidueData& residueData)
{
	m_alpha = residueData.m_alpha;
	m_beta = residueData.m_beta;
	m_gamma = residueData.m_gamma;
	m_delta = residueData.m_delta;
	m_epsilon = residueData.m_epsilon;
	m_zeta = residueData.m_zeta;
	m_nu0 = residueData.m_nu0;
	m_nu1 = residueData.m_nu1;
	m_nu2 = residueData.m_nu2;
	m_nu3 = residueData.m_nu3;
	m_nu4 = residueData.m_nu4;
	m_chi = residueData.m_chi;
	m_Phi = residueData.m_Phi;
	m_Psi = residueData.m_Psi;
	m_Eta = residueData.m_Eta;
	m_Theta = residueData.m_Theta;
}

ResidueData::~ResidueData()
{
}
