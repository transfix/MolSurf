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
#ifndef INTERSECT_H
#define INTERSECT_H

#include <PocketTunnel/datastruct_pt.h>
#include <PocketTunnel/util.h>

bool does_intersect_ray3_seg3_in_plane(const cgal_Ray &r, const Segment &s);
Point intersect_ray3_seg3(const cgal_Ray &r, const Segment &s,
                          bool &is_correct_intersection);
bool does_intersect_convex_polygon_segment_3_in_3d(
    const vector<Point> &conv_poly, const Segment &s);

#endif
