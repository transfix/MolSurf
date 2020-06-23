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
#include <PocketTunnel/handle.h>

void detect_handle(Triangulation& triang, map<int, cell_cluster>& cluster_set)
{
	for(FFI fit = triang.finite_facets_begin(); fit != triang.finite_facets_end(); fit ++)
	{
		Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
		int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
		c[0]->f_visited[id[0]] = false;
		c[1]->f_visited[id[1]] = false;
	}
	// A facet belongs to a mouth, if
	// both its adjacents cells are outside and
	// one of them is in cluster while the other is not.
	for(FFI fit = triang.finite_facets_begin(); fit != triang.finite_facets_end(); fit ++)
	{
		Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
		int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
		if(c[0]->f_visited[id[0]])
		{
			continue;
		}
		// if any adjacent cell is inside, continue.
		if(!c[0]->outside || !c[1]->outside)
		{
			continue;
		}
		// if both the adjacent cells are in valid cluster, continue.
		// or none of them is in valid cluster, continue.
		if(cluster_set[c[0]->id].in_cluster ==
				cluster_set[c[1]->id].in_cluster)
		{
			continue;
		}
		// the following will be true.
		CGAL_assertion(cluster_set[c[0]->id].find() != cluster_set[c[1]->id].find());
		CGAL_assertion(!cluster_set[c[0]->id].in_cluster ||
					   !cluster_set[c[1]->id].in_cluster);
		// if it is a surface facet continue.
		CGAL_assertion(cluster_set[c[0]->id].in_cluster !=
					   cluster_set[c[1]->id].in_cluster);
		// the cluster id is of the owner of the cell which is in cluster.
		int cl_id = cluster_set[c[0]->id].in_cluster?
					cluster_set[c[0]->id].find():
					cluster_set[c[1]->id].find();
		vector<Facet> mouth;
		mouth.push_back((*fit));
		c[0]->f_visited[id[0]] = true;
		c[1]->f_visited[id[1]] = true;
		c[0]->set_mouth(id[0], true);
		c[1]->set_mouth(id[1], true);
		while(! mouth.empty())
		{
			Facet f = mouth.back();
			mouth.pop_back();
			// circulate around each edge to collect a face that belongs to
			// the same mouth.
			Cell_handle _c = f.first;
			int _id = f.second;
			for(int i = 1; i < 4; i ++)
			{
				int u = (_id+i)%4, v = (_id+i+1)%4;
				if(i==3)
				{
					v = (_id+1)%4;
				}
				Facet_circulator fcirc = triang.incident_facets(Edge(_c,u,v));
				Facet_circulator begin = fcirc;
				do
				{
					Cell_handle cur_c[2] = {(*fcirc).first, (*fcirc).first->neighbor((*fcirc).second)
										   };
					int cur_id[2] = {(*fcirc).second, cur_c[1]->index(cur_c[0])};
					// check if this facet belongs to the same mouth.
					if(cur_c[0]->f_visited[cur_id[0]])
					{
						fcirc++;
						continue;
					}
					// if any adjacent cell is inside, continue.
					if(!cur_c[0]->outside || !cur_c[1]->outside)
					{
						fcirc++;
						continue;
					}
					// if both the adjacent cells are in valid cluster
					// or none of them is in valid cluster, continue.
					if(cluster_set[cur_c[0]->id].in_cluster ==
							cluster_set[cur_c[1]->id].in_cluster)
					{
						fcirc++;
						continue;
					}
					int cur_cl_id = cluster_set[cur_c[0]->id].in_cluster?
									cluster_set[cur_c[0]->id].find():
									cluster_set[cur_c[1]->id].find();
					if(cur_cl_id != cl_id)
					{
						fcirc++;
						continue;
					}
					// push this facet to the mouth.
					mouth.push_back((*fcirc));
					cur_c[0]->f_visited[cur_id[0]] = true;
					cur_c[1]->f_visited[cur_id[1]] = true;
					cur_c[0]->set_mouth(cur_id[0], true);
					cur_c[1]->set_mouth(cur_id[1], true);
					fcirc++;
				}
				while(fcirc != begin);
			}
		}
		cluster_set[cl_id].mouth_cnt++;
	}
}
