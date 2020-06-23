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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float   InitialAntiCausalCoefficient(float*, int, float);
float	InitialCausalCoefficient(float*, int, float, float);
void	ConvertToInterpolationCoefficients(float*, int, float*, int ,float);
void    TransImg2Spline(float*, float*, int, int, int);
double  BS_Fun(double);
double  BS_GraFun(double);
void	GradientAtPoint(float*,float , float , float , int, int , int, float*);
