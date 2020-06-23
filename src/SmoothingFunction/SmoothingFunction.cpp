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
#include <SmoothingFunction/SmoothingFunction.h>

SmoothingFunction::SmoothingFunction(double alpha1, int m1, int n1, int N1)
{
	alpha = alpha1;
	m = m1;
	n = n1;
	N = N1;
	m_Phi = 0;
	m_PhiHat = 0;
	m_PhiLength = N1*(m1+1);
	m_PhiHatLength = n1;
	m_Initialized = false;
}

SmoothingFunction::~SmoothingFunction()
{
	if(m_Phi)
	{
		delete [] m_Phi;
		m_Phi = 0;
	}
	if(m_PhiHat)
	{
		delete [] m_PhiHat;
		m_PhiHat = 0;
	}
}
