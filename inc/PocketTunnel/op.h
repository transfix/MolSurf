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
#ifndef _OP_H
#define _OP_H

#include <Utility/utility.h>
#include <PocketTunnel/datastruct_pt.h>
#include <Geometry/Geometry.h>

void write_wt(const Triangulation& triang, const char* file_prefix);
void convert_pocket_tunnel_to_rawc_geometry(Geometry** PTV, const Triangulation& triang, map<int, cell_cluster> &cluster_set, const vector<int> &sorted_cluster_index_vector, const int& opc, const int& otc);
void convert_curated_surf_to_raw_geometry(Geometry** curated_surf, const Triangulation& triang, map<int, cell_cluster> &cluster_set, const vector<int> &sorted_cluster_index_vector, const int& opc, const int& otc);

#endif
