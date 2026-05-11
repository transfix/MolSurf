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
#include <PocketTunnel/hfn_util.h>

bool is_maxima(const Cell_handle &c) {
  return (!is_outflow(Facet(c, 0)) && !is_outflow(Facet(c, 1)) &&
          !is_outflow(Facet(c, 2)) && !is_outflow(Facet(c, 3)));
}

bool is_outflow(const Facet &f) {
  Cell_handle c = f.first;
  int id = f.second;
  Point VV = c->voronoi();
  Point p[4] = {c->vertex((id + 1) % 4)->point(),
                c->vertex((id + 2) % 4)->point(),
                c->vertex((id + 3) % 4)->point(), c->vertex(id)->point()};
  return (CGAL::to_double((Tetrahedron(p[0], p[1], p[2], p[3]).volume()) *
                          (Tetrahedron(p[0], p[1], p[2], VV).volume())) < 0);
}

bool is_transversal_flow(const Facet &f) {
  Cell_handle c = f.first;
  int id = f.second;
  Point p[3];
  for (int i = 0; i < 3; i++) {
    p[i] = c->vertex((id + i + 1) % 4)->point();
  }
  for (int i = 0; i < 3; i++) {
    Vector v0 = p[(i + 1) % 3] - p[i];
    Vector v1 = p[(i + 2) % 3] - p[i];
    if (CGAL::to_double(v0 * v1) < 0) {
      return true;
    }
  }
  return false;
}

bool find_acceptor(const Cell_handle &c, const int &id, int &uid, int &vid,
                   int &wid) {
  if (!is_transversal_flow(Facet(c, id))) {
    return false;
  }
  Point p[3] = {c->vertex((id + 1) % 4)->point(),
                c->vertex((id + 2) % 4)->point(),
                c->vertex((id + 3) % 4)->point()};
  for (int i = 0; i < 3; i++)
    if (is_obtuse(p[(i + 1) % 3], p[(i + 2) % 3], p[i])) {
      wid = (id + i + 1) % 4;
    }
  vertex_indices(id, wid, uid, vid);
  return true;
}

bool is_i2_saddle(const Facet &f) {
  Cell_handle c[2];
  int id[2];
  c[0] = f.first;
  id[0] = f.second;
  c[1] = c[0]->neighbor(id[0]);
  id[1] = c[1]->index(c[0]);
  Point p[3];
  p[0] = c[0]->vertex((id[0] + 1) % 4)->point();
  p[1] = c[0]->vertex((id[0] + 2) % 4)->point();
  p[2] = c[0]->vertex((id[0] + 3) % 4)->point();
  Tetrahedron t[2];
  t[0] = Tetrahedron(p[0], p[1], p[2], c[0]->voronoi());
  t[1] = Tetrahedron(p[0], p[1], p[2], c[1]->voronoi());
  if (CGAL::to_double(t[0].volume() * t[1].volume()) >= 0) {
    return false;
  }
  for (int i = 0; i < 3; i++) {
    Vector v1 = c[0]->vertex((id[0] + (i + 1) % 3 + 1) % 4)->point() -
                c[0]->vertex((id[0] + i + 1) % 4)->point();
    Vector v2 = c[0]->vertex((id[0] + (i + 2) % 3 + 1) % 4)->point() -
                c[0]->vertex((id[0] + i + 1) % 4)->point();
    if (cosine(v1, v2) < 0) {
      return false;
    }
  }
  return true;
}

bool is_i1_saddle(const Edge &e, const Triangulation &triang) {
  // create the VF from e and triang.
  vector<Point> VF;
  Facet_circulator fcirc = triang.incident_facets(e);
  Facet_circulator begin = fcirc;
  do {
    if (triang.is_infinite((*fcirc).first)) {
      return false;
    }
    VF.push_back((*fcirc).first->voronoi());
    fcirc++;
  } while (fcirc != begin);
  return does_intersect_convex_polygon_segment_3_in_3d(
      VF, Segment(e.first->vertex(e.second)->point(),
                  e.first->vertex(e.third)->point()));
}

bool is_acceptor_for_any_VE(const Triangulation &triang, const Edge &e) {
  Cell_handle cell = e.first;
  int uid = e.second, vid = e.third;
  Facet_circulator fcirc = triang.incident_facets(e);
  Facet_circulator begin = fcirc;
  do {
    Cell_handle c[2];
    int id[2];
    c[0] = (*fcirc).first;
    id[0] = (*fcirc).second;
    c[1] = c[0]->neighbor(id[0]);
    id[1] = c[1]->index(c[0]);
    // VE is between c[0]->voronoi() and c[1]->voronoi().
    // dual is f(c[0], id[0]) = f(c[1], id[1]).
    if (!is_transversal_flow(Facet(c[0], id[0]))) {
      CGAL_assertion(!is_transversal_flow(Facet(c[1], id[1])));
      fcirc++;
      continue;
    }
    // find if this VF is acceptor for this VE.
    int cur_uid = -1, cur_vid = -1, cur_wid = -1;
    if (find_acceptor(c[0], id[0], cur_uid, cur_vid, cur_wid)) {
      if ((c[0]->vertex(cur_uid)->id == cell->vertex(uid)->id &&
           c[0]->vertex(cur_vid)->id == cell->vertex(vid)->id) ||
          (c[0]->vertex(cur_uid)->id == cell->vertex(vid)->id &&
           c[0]->vertex(cur_vid)->id == cell->vertex(uid)->id)) {
        return true;
      }
    }
    fcirc++;
  } while (fcirc != begin);
  return false;
}
