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
#ifndef POCKET_TUNNEL_UTIL_H
#define POCKET_TUNNEL_UTIL_H

#include <Utility/utility.h>
#include <PocketTunnel/datastruct_pt.h>
#include <PocketTunnel/robust_cc.h>

bool check_del_vor_property(Triangulation &triang);
bool does_allow_ball_R(const Facet &f, const double &R);
bool find_0_volume_tetrahedron(Triangulation &triang);
bool identify_cospherical_neighbor(Triangulation &triang);
bool is_inf_VF(const Triangulation &triang, const Cell_handle &c, const int uid,
               const int vid);
bool is_obtuse(const Point &p0, const Point &p1, const Point &p2);
double cell_volume(const Cell_handle &c);
double cosine(const Vector &v, const Vector &w);
double length_of_seg(const Segment &s);
int edge_index(const int facet_index, const int first_vertex_index,
               const int second_vertex_index);
int find_third_vertex_index(const Facet &f, Vertex_handle v, Vertex_handle w);
void cluster_cospherical_tetrahedra(Triangulation &triang);
void find_flow_direction(Triangulation &triang);
void identify_sink_and_saddle(Triangulation &triang);
void normalize(Vector &v);
void vertex_indices(const int facet_index, const int edge_index,
                    int &first_vertex, int &second_vertex);

#endif
