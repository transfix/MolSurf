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
#include <Utility/utility.h>
#include <Decimation/Decimation.h>
#include <Decimation/arandDecimation.h>

// note: control is no longer used since arandDecimate has replaced XuDecimate
void Decimation::Decimate(const char* in_file, const char* out_file, int niter, float control, float vari, float inangle)
{
	int showpoly;           // show poly or not (doesn't matter)
	int grouping;           // ?? some flag
	int retriways;          // if (retriways  == 1) MinimalAngleConnection
	int savenorm;            // save file flag
	//Geoframe input = Geoframe(in_file);
	Geoframe input(in_file);
	showpoly = 0;
	grouping = 1;
	retriways = 1;
	savenorm = 1;
	float controllor = control;
	float variation = vari;
	float innerangle = inangle;
	ArandDecimation ad;
	ad.arandDecimate(input, niter, 0.0, inangle, vari);
	ad.printRaw(out_file);
}

// note: control is no longer used since arandDecimate has replaced XuDecimate
void Decimation::DecimatePercentage(const char* in_file, const char* out_file, float percentage, float control, float vari, float inangle)
{
	int showpoly;           // show poly or not (doesn't matter)
	int grouping;           // ?? some flag
	int retriways;          // if (retriways  == 1) MinimalAngleConnection
	int savenorm;            // save file flag
	//Geoframe input = Geoframe(in_file);
	Geoframe input(in_file);
	showpoly = 0;
	grouping = 1;
	retriways = 1;
	savenorm = 1;
	float controllor = control;
	float variation = vari;
	float innerangle = inangle;
	ArandDecimation ad;
	ad.arandDecimate(input, 40, percentage, inangle, vari);
	ad.printRaw(out_file);
}
