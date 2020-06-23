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
#ifndef MISC_H
#define MISC_H

#include <cmath>
#include <cstdio>
#include <iomanip>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

template <typename T> T Distance_Two_points3D(T* P1, T* P2)
{
	T result;
	result = sqrt((P1[0] - P2[0])*(P1[0] - P2[0]) +
				  (P1[1] - P2[1])*(P1[1] - P2[1]) +
				  (P1[2] - P2[2])*(P1[2] - P2[2]));
	return(result);
}

//Template of Max_Of_Two and Min_Of_Two which return the biggest and smallest of three
template  <class T> T Max_Of_Two(T a, T b)
{
	if(a>=b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

template  <class T> T  Min_Of_Two(T a, T b)
{
	if(a<=b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

//Template of Max_Of_Three, Min_Of_Three, Mid_Of_Three  which return the biggest, smallest and middle of three
template <class T> T Max_Of_Three(T a, T b, T c)
{
	if((a <= b) && (c <= b))
	{
		return b;
	}
	if((a <= c) && (b <= c))
	{
		return c;
	}
	if((b <= a) && (c <= a))
	{
		return a;
	}
}

template <class T> T Min_Of_Three(T a, T b, T c)
{
	if((a >=b) && (c >= b))
	{
		return b;
	}
	if((a >=c) && (b >= c))
	{
		return c;
	}
	if((b >=a) && (c >= a))
	{
		return a;
	}
}

template <class T> T Mid_Of_Three(T a, T b, T c)
{
	if((a >=b) && (b >= c))
	{
		return b;
	}
	if((a >=c) && (c >= b))
	{
		return c;
	}
	if((b >=a) && (a >= c))
	{
		return a;
	}
}

#endif
