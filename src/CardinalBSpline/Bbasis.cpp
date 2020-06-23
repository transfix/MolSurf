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
#include <CardinalBSpline/Bbasis.h>

float CardinalBBasis(int order, float x, float r)
{
	float absx = fabs(x);
	float tem = absx/r;
	if(order ==0)
	{
		if(absx < 0.5f*r)
		{
			return 1.0f;
		}
		else
		{
			return 0.0f;
		}
	}
	if(order == 1)
	{
		if(absx < r)
		{
			return (1.0f - tem);
		}
		else
		{
			return 0.0f;
		}
	}
	if(order == 2)
	{
		if(absx < 0.5f*r)
		{
			return (0.5f*(1.5f - 2.0f*tem*tem));
		}
		if((absx >= 0.5f*r) && (absx< 1.5f*r))
		{
			return (0.5f*(1.5f-tem)*(1.5f-tem));
		}
		else
		{
			return 0.0f;
		}
	}
	if(order == 3)
	{
		if(absx < r)
		{
			return (2.0f/3.0f -tem*tem + 0.5f*tem*tem*tem);
		}
		if((absx >= r) && (absx< 2.0f*r))
		{
			return (1.0f/6.0f *(2.0f-tem)*(2.0f-tem)*(2.0f-tem));
		}
		else
		{
			return 0.0f;
		}
	}
	if(order == 4)
	{
		if(x < -2.5*r)
		{
			return 0.0f;
		}
		if((x >= -2.5*r) && (x < -1.5*r))
		{
			return ((2.0*x/r+5.0)*(2.0*x/r+5.0)*(2.0*x/r+5.0)*(2.0*x/r+5.0)/384.0);
		}
		if((x >= -1.5*r) && (x < -0.5*r))
		{
			return (55.0/96.0 - 5.0 *x/(24.0*r) - 5.0*tem*tem/4.0 - 5.0*tem*tem*x/(6.0*r) -tem*tem*tem*tem/6.0);
		}
		if((x >= -0.5*r) && (x < 0.5*r))
		{
			return (115.0/192.0 - 5.0*tem * tem /8.0 + 0.25*  tem*tem *tem *tem);
		}
		if((x >= 0.5*r) && (x < 1.5*r))
		{
			return (55.0/96.0 + 5.0*x/(24.0*r) - 5.0*tem*tem/4.0 + 5.0*tem*tem*x/(6.0*r) - tem*tem*tem*tem/6.0);
		}
		if((x >= 1.5*r) && (x < 2.5*r))
		{
			return (625.0/384.0 - 125.0 * x/(48.0*r)+25.0*tem*tem/16.0 -5.0*tem*tem*x/(12.0*r) + tem*tem*tem*tem/24.0);
		}
		if(x >= 2.5*r)
		{
			return 0.0f;
		}
	}
}
