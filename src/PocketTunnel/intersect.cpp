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
#include <PocketTunnel/intersect.h>

// Compute the intersection of one line segment and one ray<->.
// store the result in the hit_p point.
bool does_intersect_ray3_seg3_in_plane(const cgal_Ray& r, const Segment& s)
{
	bool is_correct_intersection = false;
	// steps are
	// transform start_p and driver to this frame
	// find the intersection between trans(p1 - p2) and
	//                     trans(driver -> start_p)
	// step 1 : set up a local frame with p1 as origin
	Point org = s.point(0);
	Vector xaxis, yaxis, zaxis = CGAL::NULL_VECTOR;
	//   normalized([p1 -> p2]) as x axis.
	xaxis = s.point(1) - s.point(0);
	normalize(xaxis);
	//   normalized(vec(p1,p2) X vec(p1, start_p)) as z axis
	zaxis = CGAL::cross_product(xaxis, r.source() - org);
	normalize(zaxis);
	yaxis = CGAL::cross_product(zaxis, xaxis);
	normalize(yaxis);
	// convert the points into this 2d frame.
	// condition imposed :
	//    1. segment is aligned with xaxis.
	CGAL::Point_2<Rep> _s0 = CGAL::Point_2<Rep>(0, 0);
	CGAL::Point_2<Rep> _s1 = CGAL::Point_2<Rep>(length_of_seg(s), 0);
	CGAL::Segment_2<Rep> _s = CGAL::Segment_2<Rep>(_s0, _s1);
	// condition imposed :
	//    2. ray is transformed in 2d.
	Point r0 = r.source();
	CGAL::Point_2<Rep> _r0 = CGAL::Point_2<Rep>(CGAL::to_double((r0 - org)*xaxis), CGAL::to_double((r0 - org)*yaxis));
	Point r1 = r0 + r.to_vector();
	CGAL::Point_2<Rep> _r1 = CGAL::Point_2<Rep>(CGAL::to_double((r1 - org)*xaxis), CGAL::to_double((r1 - org)*yaxis));
	CGAL::Ray_2<Rep> _r = CGAL::Ray_2<Rep>(_r0, _r1);
	// find the intersection between transformed ray and transformed segment.
	CGAL::Object result = CGAL::intersection(_r, _s);
	CGAL::Point_2<Rep> _p;
	// the computation is correct if the intersection is a point.
	is_correct_intersection = CGAL::assign(_p, result);
	return  is_correct_intersection;
}

// Compute the intersection of one line segment and one ray<->.
// store the result in the hit_p point.
Point intersect_ray3_seg3(const cgal_Ray& r, const Segment& s, bool& is_correct_intersection)
{
	// steps are
	// transform start_p and driver to this frame
	// find the intersection between trans(p1 - p2) and trans(driver -> start_p)
	// step 1 : set up a local frame with p1 as origin
	Point org = s.point(0);
	Vector xaxis, yaxis, zaxis = CGAL::NULL_VECTOR;
	//   normalized([p1 -> p2]) as x axis.
	xaxis = s.point(1) - s.point(0);
	normalize(xaxis);
	//   normalized(vec(p1,p2) X vec(p1, start_p)) as z axis
	zaxis = CGAL::cross_product(xaxis, r.source() - org);
	normalize(zaxis);
	yaxis = CGAL::cross_product(zaxis, xaxis);
	normalize(yaxis);
	// convert the points into this 2d frame.
	// condition imposed :
	//    1. segment is aligned with xaxis.
	CGAL::Point_2<Rep> _s0 = CGAL::Point_2<Rep>(0, 0);
	CGAL::Point_2<Rep> _s1 = CGAL::Point_2<Rep>(length_of_seg(s), 0);
	CGAL::Segment_2<Rep> _s = CGAL::Segment_2<Rep>(_s0, _s1);
	// condition imposed :
	//    2. ray is transformed in 2d.
	Point r0 = r.source();
	CGAL::Point_2<Rep> _r0 = CGAL::Point_2<Rep>(CGAL::to_double((r0 - org)*xaxis), CGAL::to_double((r0 - org)*yaxis));
	Point r1 = r0 + r.to_vector();
	CGAL::Point_2<Rep> _r1 = CGAL::Point_2<Rep>(CGAL::to_double((r1 - org)*xaxis), CGAL::to_double((r1 - org)*yaxis));
	CGAL::Ray_2<Rep> _r = CGAL::Ray_2<Rep>(_r0, _r1);
	// find the intersection between transformed ray and transformed segment.
	CGAL::Object result = CGAL::intersection(_r, _s);
	CGAL::Point_2<Rep> _p;
	// the computation is correct if the intersection is a point.
	is_correct_intersection = CGAL::assign(_p, result);
	if(! is_correct_intersection)
	{
		cerr << "x";
		return CGAL::ORIGIN;
	}
	else
	{
		// convert _p into 3d.
		CGAL_assertion(_s.has_on(_p));
		double t = sqrt(CGAL::to_double((_p - _s0)*(_p - _s0)))/length_of_seg(s);
		return s.point(0) + t*(s.point(1) - s.point(0));
	}
}

//a flag that will be set if the call to circumsphere fails
static bool cgal_failed;

//the handler that will flag the operation as having failed
//oh, what the heck, let's make it static while we're at it
static void failure_func(const char* type, const char* exp, const char* file, int line, const char* expl)
{
	//bad CGAL!
	cgal_failed = true;
	return;
}


bool does_intersect_convex_polygon_segment_3_in_3d(const vector<Point>& conv_poly, const Segment& s)
{
	// dissect the polygon into a set of triangles with common apex at the centroid.
	Vector centroid = CGAL::NULL_VECTOR;
	vector<Triangle_3> t;
	for(int i = 0; i < (int)conv_poly.size(); i ++)
	{
		centroid = centroid + (conv_poly[i] - CGAL::ORIGIN);
	}
	centroid = (1./(int)conv_poly.size()) * centroid;
	for(int i = 0; i < (int)conv_poly.size(); i ++)
		t.push_back(Triangle_3(CGAL::ORIGIN + centroid, 						   conv_poly[i], conv_poly[(i+1)%((int)conv_poly.size())]));
	//CGAL::Failure_function old_ff = CGAL::set_error_handler(failure_func);
	//CGAL::Failure_behaviour old_fb = CGAL::set_error_behaviour(CGAL::CONTINUE);
	// Now for every triangle in t check if the Delaunay edge e intersects it.
	// use CGAL routine.
	bool is_degenerate_poly = true;
	for(int i = 0; i < (int)t.size(); i ++)
	{
		cgal_failed = false;
		bool intersect_result = CGAL::do_intersect(t[i], s);
		if(cgal_failed)
		{
			continue;
		}
		is_degenerate_poly = false;
		if(intersect_result)
		{
			//CGAL::set_error_handler(old_ff);
			//CGAL::set_error_behaviour(old_fb);
			is_degenerate_poly = false;
			return true;
		}
	}
	//CGAL::set_error_handler(old_ff);
	//CGAL::set_error_behaviour(old_fb);
	/*if(is_degenerate_poly)
	{
		cerr << "+";
	}*/
	return false;
}
