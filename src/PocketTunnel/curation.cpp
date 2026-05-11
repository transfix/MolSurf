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
#include <PocketTunnel/curation.h>

using std::endl;

void curate(Triangulation &triang, map<int, cell_cluster> &cluster_set,
            const vector<int> &sorted_cluster_index_vector,
            const int output_tunnel_count, const int output_pocket_count) {
  int tun_void_cnt = 0;
  for (int i = 0; i < (int)sorted_cluster_index_vector.size(); i++) {
    if (i >= (int)sorted_cluster_index_vector.size()) {
      break;
    }
    int cl_id = sorted_cluster_index_vector[i];
    if (cluster_set[cl_id].mouth_cnt < 1) {
      tun_void_cnt++;
    } else if (cluster_set[cl_id].mouth_cnt > 1) {
      tun_void_cnt++;
    } else {
      continue;
    }
    if (tun_void_cnt <= output_tunnel_count) {
      cerr << "Tunnel/Void number " << tun_void_cnt << " is not to be curated."
           << endl;
      continue;
    }
    // currently we do not curate more than 100 tunnel/void.
    if (tun_void_cnt > 100) {
      break;
    }
    // curate this tunnel/void.
    for (FCI cit = triang.finite_cells_begin();
         cit != triang.finite_cells_end(); cit++) {
      if (cluster_set[cit->id].find() != cl_id) {
        continue;
      }
      cit->outside = false;
    }
  }
  int pocket_cnt = 0;
  for (int i = 0; i < (int)sorted_cluster_index_vector.size(); i++) {
    if (i >= (int)sorted_cluster_index_vector.size()) {
      break;
    }
    int cl_id = sorted_cluster_index_vector[i];
    if (cluster_set[cl_id].mouth_cnt < 1) {
      continue;
    } else if (cluster_set[cl_id].mouth_cnt > 1) {
      continue;
    } else {
      pocket_cnt++;
    }
    if (pocket_cnt <= output_pocket_count) {
      cerr << "Pocket number " << pocket_cnt << " is not to be curated."
           << endl;
      continue;
    }
    // currently we do not curate more than 100 tunnel/void.
    if (pocket_cnt > 100) {
      break;
    }
    // curate this pocket.
    for (FCI cit = triang.finite_cells_begin();
         cit != triang.finite_cells_end(); cit++) {
      if (cluster_set[cit->id].find() != cl_id) {
        continue;
      }
      cit->outside = false;
    }
  }
  // remove bubbles near the surface.
  for (FCI cit = triang.finite_cells_begin(); cit != triang.finite_cells_end();
       cit++)
    if (!cluster_set[cit->id].outside) {
      cit->outside = false;
    }
}
