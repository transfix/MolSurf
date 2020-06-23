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
#ifndef TCOCONE_H
#define TCOCONE_H

#include <Utility/utility.h>
#include <PocketTunnel/datastruct_pt.h>
#include <PocketTunnel/util.h>
#include <PocketTunnel/robust_cc.h>

void compute_poles(Triangulation& triang);
void mark_flat_vertices(Triangulation& triang, double ratio, double cocone_phi, double flat_phi);
void tcocone(const double DEFAULT_ANGLE, const double DEFAULT_SHARP, const double DEFAULT_FLAT, const double DEFAULT_RATIO, Triangulation& triang);

#endif
