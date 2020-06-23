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
#if !defined(AFX_SMOOTHINGFUNCTION_H__29D14D8B_2356_4D85_86BF_AAD558BB98E5__INCLUDED_)
#define AFX_SMOOTHINGFUNCTION_H__29D14D8B_2356_4D85_86BF_AAD558BB98E5__INCLUDED_

#include <Utility/utility.h>

class SmoothingFunction
{
	public:
		SmoothingFunction(double alpha, int m, int n, int N);
		virtual ~SmoothingFunction();
		inline double getPhi(double index);
		inline double getPhiHat(int index);
		virtual inline double getPhiAtRealPos(double pos) = 0;
		virtual inline double getPhiHatAtRealPos(double pos) = 0;

	protected:
		virtual bool precompute() = 0;
		double* m_Phi;
		double* m_PhiHat;
		int m_PhiLength;
		int m_PhiHatLength;
		double alpha;
		int m;
		int n;
		int N; // used to compare with large FFT.
		bool m_Initialized;
};

// this index is a real number between 0 and m+1.
// When multiplied by N, this should be an integer.
// This is useful to compare with a FFT method.
double SmoothingFunction::getPhi(double index)
{
	if(!m_Initialized)
	{
		precompute();
	}
	if(!m_Phi)
	{
		return 0;
	}
	if((index<-(m+1)) || (index >= m+1))
	{
		return 0;
	}
	return m_Phi[(int) fabs((double)(index*N))];
}

double SmoothingFunction::getPhiHat(int index)
{
	if(!m_Initialized)
	{
		precompute();
	}
	if(!m_PhiHat)
	{
		return 0;
	}
	if((index<0) || (index >= m_PhiHatLength))
	{
		return 0;
	}
	return m_PhiHat[index];
}

#endif
