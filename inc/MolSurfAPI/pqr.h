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
// Tools for working with PQR files

#ifndef _PQR_H_
#define _PQR_H_

#include <MolSurfAPI/molsurf.h>
#include <GOAFileTypes/GOALoader.h>

MOLECULE::GroupOfAtoms* ParsePQR(FILE* fp, MOLECULE::RADIUS_TYPE radiusType);
MOLECULE::GroupOfAtoms* ParsePDB(FILE* fp, MOLECULE::RADIUS_TYPE radiusType);
MOLECULE::GroupOfAtoms* ParseOptimalPQR(FILE* fp, MOLECULE::RADIUS_TYPE radiusType);
void savePQR(MOLECULE::GroupOfAtoms* goa, FILE* fp);
void savePDB(MOLECULE::GroupOfAtoms* goa, FILE* fp, bool withH, bool withChain);
void saveCDR(MOLECULE::GroupOfAtoms* goa, FILE* fp);
void append(MOLECULE::GroupOfAtoms* goaA, MOLECULE::GroupOfAtoms* goaB);

#endif
