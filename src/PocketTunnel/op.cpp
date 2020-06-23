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
#include <PocketTunnel/op.h>

// Write out the boundary between inside and outside tetrehedra as surface.
void write_wt(const Triangulation& triang, const char* file_prefix)
{
	char filename[100];
	strcat(strcpy(filename, file_prefix), ".surf");
	ofstream fout;
	fout.open(filename);
	if(! fout)
	{
		cerr << "Can not open " << filename << " for writing. " << endl;
		exit(1);
	}
	// Count number of facets on the surface.
	int num_facets = 0;
	for(FFI fit = triang.finite_facets_begin();
			fit != triang.finite_facets_end(); ++fit)
		if((*fit).first->cocone_flag((*fit).second))
		{
			num_facets ++;
		}
	fout <<"OFF" << endl;
	fout << triang.number_of_vertices() << " " << num_facets << " 0" << endl;
	// Write the vertices.
	for(FVI vit = triang.finite_vertices_begin(); vit != triang.finite_vertices_end(); ++vit)
	{
		fout << vit->point() << endl;
	}
	// Write the facets.
	for(FFI fit = triang.finite_facets_begin();
			fit != triang.finite_facets_end(); ++fit)
	{
		Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
		int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
		if(! c[0]->cocone_flag(id[0]))
		{
			continue;
		}
		CGAL_assertion(c[0]->bdy[id[0]] && c[1]->bdy[id[1]]);
		CGAL_assertion(c[0]->outside != c[1]->outside);
		Vertex_handle vh[3] = { c[0]->vertex((id[0]+1)%4), c[0]->vertex((id[0]+2)%4), c[0]->vertex((id[0]+3)%4)
							  };
		if(! c[0]->outside)
			if(CGAL::is_negative(Tetrahedron(vh[0]->point(), vh[1]->point(), vh[2]->point(), c[0]->vertex(id[0])->point()).volume()))
			{
				fout << "3\t" << vh[0]->id << " " << vh[1]->id << " " << vh[2]->id << " ";
			}
			else
			{
				fout << "3\t" << vh[1]->id << " " << vh[0]->id << " " << vh[2]->id << " ";
			}
		else if(CGAL::is_negative(Tetrahedron(vh[0]->point(), vh[1]->point(), vh[2]->point(), c[1]->vertex(id[1])->point()).volume()))
		{
			fout << "3\t" << vh[0]->id << " " << vh[1]->id << " " << vh[2]->id << " ";
		}
		else
		{
			fout << "3\t" << vh[1]->id << " " << vh[0]->id << " " << vh[2]->id << " ";
		}
		fout << "1 1 1 0.3" << endl;
	}
	fout.close();
}

void convert_pocket_tunnel_to_rawc_geometry(Geometry** PTV, const Triangulation& triang, map<int, cell_cluster> &cluster_set, const vector<int> &sorted_cluster_index_vector, const int& opc, const int& otc)
{
	// Create a single rawc geometry from the handles.
	vector<int> facet_ids;
	vector<float> vertex_colors;
	vertex_colors.resize(3*triang.number_of_vertices(), 1);
	int facet_count = 0;
	// Load the geometry of the tunnels and voids.
	int tun_void_cnt = 0;
	int i = -1;
	for(; tun_void_cnt < otc;)
	{
		i++;
		if(i >= (int)sorted_cluster_index_vector.size())
		{
			cerr << endl << "The number of tunnel/voids are less than " << otc << endl;
			break;
		}
		if(tun_void_cnt >= 100)
		{
			cerr << "More than 100 tunnel/voids will not be output." << endl;
			break;
		}
		int cl_id = sorted_cluster_index_vector[i];
		if(cluster_set[cl_id].mouth_cnt < 1)
		{
			cerr << "Cluster " << cl_id << " is a void." << endl;
		}
		else if(cluster_set[cl_id].mouth_cnt > 1)
			cerr << "Cluster " << cl_id << " is a tunnel with "
				 << cluster_set[cl_id].mouth_cnt << " mouths." << endl;
		else
		{
			continue;
		}
		tun_void_cnt++;
		for(FFI fit = triang.finite_facets_begin();
				fit != triang.finite_facets_end(); fit ++)
		{
			Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
			int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
			if(cluster_set[c[0]->id].find() ==
					cluster_set[c[1]->id].find())
			{
				continue;
			}
			if(cluster_set[c[0]->id].find() != cl_id &&
					cluster_set[c[1]->id].find() != cl_id)
			{
				continue;
			}
			// this facet is part of this tunnel/void.
			// point the face away from the cell which is in this cluster.
			Cell_handle temp_c = c[0];
			int temp_id = c[0]->index(c[1]);
			if(cluster_set[c[0]->id].find() != cl_id)
			{
				temp_c = c[1];
				temp_id = c[1]->index(c[0]);
			}
			int vid[3] = {temp_c->vertex((temp_id+1)%4)->id,
						  temp_c->vertex((temp_id+2)%4)->id,
						  temp_c->vertex((temp_id+3)%4)->id
						 };
			Tetrahedron tet(temp_c->vertex((temp_id+1)%4)->point(),
							temp_c->vertex((temp_id+2)%4)->point(),
							temp_c->vertex((temp_id+3)%4)->point(),
							temp_c->vertex(temp_id)->point());
			if(CGAL::to_double(tet.volume()) > 0)
			{
				facet_ids.push_back(vid[0]);
				facet_ids.push_back(vid[2]);
				facet_ids.push_back(vid[1]);
				facet_count++;
			}
			else
			{
				facet_ids.push_back(vid[0]);
				facet_ids.push_back(vid[1]);
				facet_ids.push_back(vid[2]);
				facet_count++;
			}
			// color.
			CGAL_assertion(cluster_set[cl_id].mouth_cnt != 1); // it is not a pocket.
			float r,g,b;
			if(cluster_set[cl_id].mouth_cnt > 1)  // tunnel
			{
				if(cluster_set[c[0]->id].outside != cluster_set[c[1]->id].outside)
				{
					r = 1;
					g = 1;
					b = 0; // yellow
				}
				else
				{
					r = 1;
					g = 0;
					b = 0; // mouth - red
				}
			}
			else
			{
				CGAL_assertion(cluster_set[cl_id].mouth_cnt == 0); // void
				if(c[0]->cocone_flag(id[0]))
				{
					r = 1;
					g = 0;
					b = 1; // purple
				}
				else
				{
					r = 0;
					g = 0;
					b = 1; // blue
				}
			}
			// set the color in vertex_colors.
			for(int j = 0; j < 3; j ++)
			{
				if(vertex_colors[3*vid[j]] == 1 &&
						vertex_colors[3*vid[j]+1] == 0 &&
						vertex_colors[3*vid[j]+2] == 0)
				{
					continue;
				}
				vertex_colors[3*vid[j]+0] = r;
				vertex_colors[3*vid[j]+1] = g;
				vertex_colors[3*vid[j]+2] = b;
			}
		}
	}
	// Now create the geometry for the pockets.
	int pocket_cnt = 0;
	i = -1;
	for(; pocket_cnt < opc;)
	{
		i++;
		if(i >= (int)sorted_cluster_index_vector.size())
		{
			cerr << endl << "No more pockets." << endl;
			break;
		}
		if(pocket_cnt >= 100)
		{
			cerr << "More than 100 pockets will not be output." << endl;
			break;
		}
		int cl_id = sorted_cluster_index_vector[i];
		if(cluster_set[cl_id].mouth_cnt < 1)
		{
			continue;
		}
		else if(cluster_set[cl_id].mouth_cnt > 1)
		{
			continue;
		}
		else
		{
			cerr << "Cluster " << cl_id << " is a pocket." << endl;
		}
		pocket_cnt++;
		for(FFI fit = triang.finite_facets_begin();
				fit != triang.finite_facets_end(); fit ++)
		{
			Cell_handle c[2] = {(*fit).first, (*fit).first->neighbor((*fit).second)};
			int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
			if(cluster_set[c[0]->id].find() ==
					cluster_set[c[1]->id].find())
			{
				continue;
			}
			if(cluster_set[c[0]->id].find() != cl_id &&
					cluster_set[c[1]->id].find() != cl_id)
			{
				continue;
			}
			Cell_handle temp_c = c[0];
			int temp_id = c[0]->index(c[1]);
			if(cluster_set[c[0]->id].find() != cl_id)
			{
				temp_c = c[1];
				temp_id = c[1]->index(c[0]);
			}
			int vid[3] = {temp_c->vertex((temp_id+1)%4)->id,
						  temp_c->vertex((temp_id+2)%4)->id,
						  temp_c->vertex((temp_id+3)%4)->id
						 };
			Tetrahedron tet(temp_c->vertex((temp_id+1)%4)->point(),
							temp_c->vertex((temp_id+2)%4)->point(),
							temp_c->vertex((temp_id+3)%4)->point(),
							temp_c->vertex(temp_id)->point());
			if(CGAL::to_double(tet.volume()) > 0)
			{
				facet_ids.push_back(vid[0]);
				facet_ids.push_back(vid[2]);
				facet_ids.push_back(vid[1]);
				facet_count++;
			}
			else
			{
				facet_ids.push_back(vid[0]);
				facet_ids.push_back(vid[1]);
				facet_ids.push_back(vid[2]);
				facet_count++;
			}
			// color.
			CGAL_assertion(cluster_set[cl_id].mouth_cnt == 1); // it is a pocket.
			float r,g,b;
			if(cluster_set[c[0]->id].outside !=
					cluster_set[c[1]->id].outside)
			{
				r = 0;
				g = 1;
				b = 0; // surface - green
			}
			else
			{
				r = 1;
				g = 0;
				b = 1; // mouth - purple
			}
			// set the color in vertex_colors.
			for(int j = 0; j < 3; j ++)
			{
				if(vertex_colors[3*vid[j]] == 1 &&
						vertex_colors[3*vid[j]+1] == 0 &&
						vertex_colors[3*vid[j]+2] == 1)
				{
					continue;
				}
				vertex_colors[3*vid[j]] = r;
				vertex_colors[3*vid[j]+1] = g;
				vertex_colors[3*vid[j]+2] = b;
			}
		}
	}
	// Now we have the vertex_ids, facet_ids and vertex_colors all set.
	// Load these in Geometry.
	// (*PTV) = new Geometry();
	(*PTV)->AllocateTris(triang.number_of_vertices(), facet_count);
	(*PTV)->AllocateTriVertColors();
	for(FVI vit = triang.finite_vertices_begin();
			vit != triang.finite_vertices_end(); vit ++)
	{
		(*PTV)->m_TriVerts[3*vit->id] = CGAL::to_double(vit->point().x());
		(*PTV)->m_TriVerts[3*vit->id+1] = CGAL::to_double(vit->point().y());
		(*PTV)->m_TriVerts[3*vit->id+2] = CGAL::to_double(vit->point().z());
		(*PTV)->m_TriVertColorsTransparent[4*vit->id] = vertex_colors[3*vit->id];
		(*PTV)->m_TriVertColorsTransparent[4*vit->id+1] = vertex_colors[3*vit->id+1];
		(*PTV)->m_TriVertColorsTransparent[4*vit->id+2] = vertex_colors[3*vit->id+2];
		(*PTV)->m_TriVertColorsTransparent[4*vit->id+3] = 1.0f;
	}
	for(i = 0; i < facet_count; i++)
	{

		(*PTV)->m_Tris[3*i] = facet_ids[3*i];
		(*PTV)->m_Tris[3*i+1] = facet_ids[3*i+1];
		(*PTV)->m_Tris[3*i+2] = facet_ids[3*i+2];
	}
/* for debug	
	ofstream fout;
	fout.open("debug_output/temp_PTV.coff");
	fout << "COFF" << endl;
	fout << (*PTV)->m_NumTriVerts << " " << (*PTV)->m_NumTris << " 0" << endl;
	for(i = 0; i < (*PTV)->m_NumTriVerts; i ++)
	{
		fout << (*PTV)->m_TriVerts[3*i] << " "
			 << (*PTV)->m_TriVerts[3*i+1] << " "
			 << (*PTV)->m_TriVerts[3*i+2] << " "
			 << (*PTV)->m_TriVertColorsTransparent[4*i] << " "
			 << (*PTV)->m_TriVertColorsTransparent[4*i+1] << " "
			 << (*PTV)->m_TriVertColorsTransparent[4*i+2] << " 1" << endl;
	}
	for(i = 0; i < (*PTV)->m_NumTris; i ++)
	{
		fout << "3\t";
		fout << (*PTV)->m_Tris[3*i] << " "
			 << (*PTV)->m_Tris[3*i+1] << " "
			 << (*PTV)->m_Tris[3*i+2] << endl;
	}
*/

  // original code of curate_tr
   tun_void_cnt = 0;
   for(int i = 0; i < (int)sorted_cluster_index_vector.size(); i++)
   {
      if(i >= (int)sorted_cluster_index_vector.size())
         break;
      int cl_id = sorted_cluster_index_vector[i];
      if( cluster_set[cl_id].mouth_cnt < 1 )
         tun_void_cnt++;
      else if( cluster_set[cl_id].mouth_cnt > 1 )
         tun_void_cnt++;
      else
         continue;
      if( tun_void_cnt <= otc )
      {
         cerr << "Tunnel/Void number " << tun_void_cnt << " is not to be curated." << endl;
         continue;
      }
      // currently we do not curate more than 100 tunnel/void.
      if( tun_void_cnt > 100 ) break;
      // curate this tunnel/void.
      for(FCI cit = triang.finite_cells_begin();
         cit != triang.finite_cells_end(); cit ++)
      {
	 if(cluster_set[cit->id].find() != cl_id ) continue;
         cit->outside = false;
         // tag this cell.
//         cit->c_tag = true;
      }
   }

   pocket_cnt = 0;
   for(int i = 0; i < (int)sorted_cluster_index_vector.size(); i++)
   {
      if(i >= (int)sorted_cluster_index_vector.size())
         break;
      int cl_id = sorted_cluster_index_vector[i];
      if( cluster_set[cl_id].mouth_cnt < 1 )
         continue;
      else if( cluster_set[cl_id].mouth_cnt > 1 )
         continue;
      else
         pocket_cnt++;

      if( pocket_cnt <= opc )
      {
         cerr << "Pocket number " << pocket_cnt << " is not to be curated." << endl;
         continue;
      }
      // currently we do not curate more than 100 tunnel/void.
      if( pocket_cnt > 100 ) break;

      // curate this pocket.
      for(FCI cit = triang.finite_cells_begin();
         cit != triang.finite_cells_end(); cit ++)
      {
	 if(cluster_set[cit->id].find() != cl_id ) continue;
         cit->outside = false;
         // tag this cell.
  //       cit->c_tag = true;
      }
   }



}

void convert_curated_surf_to_raw_geometry(Geometry** curated_surf, const Triangulation& triang, map<int, cell_cluster> &cluster_set, const vector<int> &sorted_cluster_index_vector, const int& opc, const int& otc)
{
	
	vector<int> facet_ids;
	int facet_count = 0;
	for(FFI fit = triang.finite_facets_begin(); fit != triang.finite_facets_end(); fit ++)
	{
		Cell_handle c[2] = {fit->first, fit->first->neighbor(fit->second)};
		int id[2] = {c[0]->index(c[1]), c[1]->index(c[0])};
		if(c[0]->outside != c[1]->outside)
		{
			facet_count++;
			Vertex_handle vh[3] = { c[0]->vertex((id[0]+1)%4), c[0]->vertex((id[0]+2)%4), c[0]->vertex((id[0]+3)%4)};
			if(! c[0]->outside)
			{
				if(CGAL::is_negative(Tetrahedron(vh[0]->point(), vh[1]->point(), vh[2]->point(), c[0]->vertex(id[0])->point()).volume()))
				{
					facet_ids.push_back(vh[0]->id);
					facet_ids.push_back(vh[1]->id);
					facet_ids.push_back(vh[2]->id);
				}
				else
				{
					facet_ids.push_back(vh[1]->id);
					facet_ids.push_back(vh[0]->id);
					facet_ids.push_back(vh[2]->id);
				}
			}
			else
			{
				if(CGAL::is_negative(Tetrahedron(vh[0]->point(), vh[1]->point(), vh[2]->point(), c[1]->vertex(id[1])->point()).volume()))
				{
					facet_ids.push_back(vh[0]->id);
					facet_ids.push_back(vh[1]->id);
					facet_ids.push_back(vh[2]->id);
				}
				else
				{
					facet_ids.push_back(vh[1]->id);
					facet_ids.push_back(vh[0]->id);
					facet_ids.push_back(vh[2]->id);
				}
			}
		}
		else continue;
	}
	cout<<"facet_count: "<< facet_count << endl;
	// convert to raw geometry.
	(*curated_surf)->AllocateTris(triang.number_of_vertices(), facet_count);
//	(*curated_surf)->AllocateTriVertColors();
	for(FVI vit = triang.finite_vertices_begin();
			vit != triang.finite_vertices_end(); vit ++)
	{
		(*curated_surf)->m_TriVerts[3*vit->id] = CGAL::to_double(vit->point().x());
		(*curated_surf)->m_TriVerts[3*vit->id+1] = CGAL::to_double(vit->point().y());
		(*curated_surf)->m_TriVerts[3*vit->id+2] = CGAL::to_double(vit->point().z());
	}
	for(int i = 0; i < facet_count; i++)
	{
		(*curated_surf)->m_Tris[3*i] = facet_ids[3*i];
		(*curated_surf)->m_Tris[3*i+1] = facet_ids[3*i+1];
		(*curated_surf)->m_Tris[3*i+2] = facet_ids[3*i+2];
	}

/* //for debug
    ofstream fout;
	fout.open("debug_output/temp_curated.off");
	fout << "OFF" << endl;
	fout << (*curated_surf)->m_NumTriVerts << " " << (*curated_surf)->m_NumTris << " 0" << endl;
	for(int i = 0; i < (*curated_surf)->m_NumTriVerts; i ++)
		fout << (*curated_surf)->m_TriVerts[3*i] << " "
			 << (*curated_surf)->m_TriVerts[3*i+1] << " "
			 << (*curated_surf)->m_TriVerts[3*i+2] << endl;
	for(int i = 0; i < (*curated_surf)->m_NumTris; i ++)
		fout << (*curated_surf)->m_Tris[3*i] << " "
			 << (*curated_surf)->m_Tris[3*i+1] << " "
			 << (*curated_surf)->m_Tris[3*i+2] << endl;  
*/
	return;
}
