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
#if !defined(AFX_GAUSSIAN_H__C3E0B6A2_14F1_4872_9F56_D70142C89D07__INCLUDED_)
#define AFX_GAUSSIAN_H__C3E0B6A2_14F1_4872_9F56_D70142C89D07__INCLUDED_

#include <Utility/utility.h>
#include <SmoothingFunction/SmoothingFunction.h>

class Gaussian
 : public SmoothingFunction
{
	public:
		Gaussian(double alpha, int m, int n, int N);
		virtual ~Gaussian();
		inline double getPhiAtRealPos(double pos);
		inline double getPhiHatAtRealPos(double pos);
	protected:
		bool precompute();
};

double Gaussian::getPhiAtRealPos(double pos)
{
	return 0;
}

double Gaussian::getPhiHatAtRealPos(double pos)
{
	return 0;
}

#endif
