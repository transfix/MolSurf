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
#ifndef __DATASTRUCT_PT_H__
#define __DATASTRUCT_PT_H__

#include <Utility/utility.h>
#include <CGAL/basic.h>// for compiler discrepancies
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>// implements Traits.
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/MP_Float.h>
#include <CGAL/number_utils.h>
#include <CGAL/Point_3.h>
#include <CGAL/Quotient.h>
#include <CGAL/Timer.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>// check
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/utility.h>

#define __OUTSIDE__
#define __INSIDE__

// Define a new Vertex_base with additional attributes.
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
class Tn_vertex : public Vb
{
		typedef Vb Base;
	public:
		typedef typename Vb::Point Point;
		typedef typename Point::R Rep;
		typedef typename Vb::Vertex_handle Vertex_handle;
		typedef typename Vb::Cell_handle Cell_handle;
		typedef CGAL::Vector_3<Rep> Vector;

		template < typename TDS2 >
		struct Rebind_TDS
		{
			typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
			typedef Tn_vertex<Gt,Vb2> Other;
		};

		Tn_vertex()
		{
			init();
		}
		Tn_vertex(const Point& p) : Base(p)
		{
			init();
		}
		Tn_vertex(const Point& p, void* cell)
			: Base(p, cell)
		{
			init();
		}
		Tn_vertex(void* cell) : Base(cell)
		{
			init();
		}

		void set_normal(const Vector& v)
		{
			v_normal = v;
		}
		Vector normal() const
		{
			return v_normal;
		}

		void set_pole(const Point& p)
		{
			v_pole = p;
		}
		Point pole() const
		{
			return v_pole;
		}

		void set_width(const double& w)
		{
			v_width = w;
		}
		double width() const
		{
			return v_width;
		}

		void set_diameter(const double& d)
		{
			v_diameter = d;
		}
		double diameter() const
		{
			return v_diameter;
		}

		void set_flat(const bool b)
		{
			v_flat = b;
		}
		bool is_flat() const
		{
			return v_flat;
		}

		void set_convex_hull(const bool b)
		{
			v_chull = b;
		}
		bool on_convex_hull() const
		{
			return v_chull;
		}

		void set_isolated(const bool b)
		{
			v_isolated = b;
		}
		bool is_isolated() const
		{
			return v_isolated;
		}

		void set_umbrella(const bool b)
		{
			v_umbrella = b;
		}
		bool has_umbrella() const
		{
			return v_umbrella;
		}

		// for robust cocone
		void set_on_smooth_surface(const bool b)
		{
			v_on_smooth_surface = b;
		}
		bool on_smooth_surface() const
		{
			return v_on_smooth_surface;
		}

		int id;
		bool tag;// Some tag needed during computations.
		bool visited;
		bool bad;
		bool bad_neighbor;
		vector<int > bad_v;

		// for robust cocone
		vector<double> nnd_vector;

	private:
		inline void init()
		{
			v_normal = CGAL::NULL_VECTOR;
			v_pole = this->point();
			v_width = 0.0;
			v_diameter = 0.0;
			v_flat = true;
			v_chull = false;
			v_umbrella = false;
			v_isolated = true;
			// for robust cocone
			v_on_smooth_surface = false;
			nnd_vector.clear();
			nnd_vector.resize(3, HUGE_VAL);
		}

		Vector v_normal;// Estimated normal at the vertex.
		Point v_pole;// The negative pole.
		double v_width;// Radius of the cocone.
		double v_diameter;// Diamter of the the incident umbrella.
		bool v_flat;// Indicates whether the vertex is flat.
		bool v_chull;// Indicates wheter the vertex lies on the boundary of the convex hull.
		bool v_umbrella;// Indicates if v has a (non sharp) umbrella.
		bool v_isolated;// Indicates whether the vertex is isolated.
		bool v_on_smooth_surface;// robust cocone

};




// Define a new Cell_base storing associated Voronoi vertices
// and cocone flags.
 
template < class Gt, class Cb = CGAL::Triangulation_cell_base_3<Gt> >
class Tn_cell : public Cb
{
		typedef Cb Base;
	public:
		typedef typename Cb::Point Point;
		typedef typename Cb::Vertex_handle Vertex_handle;
		typedef typename Cb::Cell_handle Cell_handle;

		Tn_cell()
		{
			init();
		}
		Tn_cell(Vertex_handle v0,
				Vertex_handle v1,
				Vertex_handle v2,
				Vertex_handle v3)
			: Base(v0,v1,v2,v3)
		{
			init();
		}
		Tn_cell(Vertex_handle v0,
				Vertex_handle v1,
				Vertex_handle v2,
				Vertex_handle v3,
				Cell_handle c0,
				Cell_handle c1,
				Cell_handle c2,
				Cell_handle c3)
			: Base(v0,v1,v2,v3,c0,c1,c2,c3)
		{
			init();
		}

		template < typename TDS2 >
		struct Rebind_TDS
		{
			typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
			typedef Tn_cell<Gt,Cb2> Other;
		};

		void set_voronoi(const Point& p)
		{
			c_voronoi = p;
		}
		Point voronoi() const
		{
			return c_voronoi;
		}

		void set_dirty(const bool& b)
		{
			c_dirty = b;
		}
		bool dirty() const
		{
			return c_dirty;
		}

		void set_cell_radius(const double r)
		{
			c_cell_radius = r;
		}
		double cell_radius() const
		{
			return c_cell_radius;
		}

		void set_cocone_flag(int i, bool new_value)
		{
			CGAL_precondition(0 <= i && i < 4);
			c_cocone_flag[i] = new_value;
		}

		bool cocone_flag(int i) const
		{
			CGAL_precondition(0 <= i && i < 4);
			return c_cocone_flag[i];
		}

		void set_removable(int i, bool new_value)
		{
			CGAL_precondition(0 <= i && i < 4);
			c_removable[i] = new_value;
		}

		bool removable(int i) const
		{
			CGAL_precondition(0 <= i && i < 4);
			return c_removable[i];
		}

		void set_mouth(const int& id, const bool& b)
		{
			c_mouth[id] = b;
		}
		bool mouth(int id) const
		{
			return c_mouth[id];
		}

		// Functions for Medial Axis
		void set_VF_on_medax(const int i, const int j, const bool val)
		{
			CGAL_precondition(0 <= i && i <= 3);
			CGAL_precondition(0 <= j && j <= 3);
			CGAL_precondition(i != j);
			c_VF_on_medax[i][j] = val;
			c_VF_on_medax[j][i] = val;
		}

		bool VF_on_medax(int i, int j) const
		{
			CGAL_precondition(0 <= i && i <= 3);
			CGAL_precondition(0 <= j && j <= 3);
			CGAL_precondition(i != j);
			return c_VF_on_medax[i][j];
		}

		void set_VV_on_medax(const bool& b)
		{
			c_VV_on_medax = b;
		}
		bool VV_on_medax() const
		{
			return c_VV_on_medax;
		}

		void set_bb(const bool flag)
		{
			c_bb = flag;
		}
		bool bb() const
		{
			return c_bb;
		}

		void set_deep_int(const int i, const bool flag)
		{
			c_deep_int[i] = flag;
		}
		bool deep_int(int i) const
		{
			return c_deep_int[i];
		}

		bool tag[4];// Some tag needed in computations.
		bool visited;
		bool f_visited[4];
		int id;
		bool flag;

		bool transp;
		bool outside;

		bool bdy[4];
		bool opaque[4];
		int umbrella_member[4][4];

		inline void init()
		{
			c_voronoi = CGAL::ORIGIN;
			c_dirty = false;
			c_cell_radius = 0.0;
			c_bb = false;
			c_VV_on_medax = false;
			for(int i = 0; i < 4; i ++)
			{
				c_cocone_flag[i] = false;
				c_removable[i] = true;
				c_deep_int[i] = false;
				f_visited[i] = false;
				c_mouth[i] = false;
				for(int j = 0; j < 4; j ++)
				{
					c_VF_on_medax[i][j] = false;
				}
			}
		}
	private:
		Point c_voronoi;
		bool c_dirty;
		double c_cell_radius;
		bool c_VV_on_medax;
		bool c_VF_on_medax[4][4];
		bool c_cocone_flag[4];
		bool c_removable[4];
		bool c_deep_int[4];
		bool c_bb;
		bool c_mouth[4];
};

//Craig: these typedefs really pollute the namespace, can they be encapuslated somehow?
struct K : CGAL::Exact_predicates_inexact_constructions_kernel	{};
typedef Tn_vertex< K >	 					myVertex;
typedef Tn_cell< K >	 					myCell;
typedef CGAL::Triangulation_data_structure_3<myVertex, myCell>	TDS;
typedef CGAL::Delaunay_triangulation_3<K,TDS>	 		Triangulation;
typedef Triangulation::Point 					Point;
typedef Point::R						Rep;
typedef CGAL::Tetrahedron_3<Rep>				Tetrahedron;
typedef CGAL::Triangle_3<Rep> 					Triangle_3;
typedef CGAL::Vector_3<Rep>					Vector;
typedef CGAL::Ray_3<Rep>					cgal_Ray;
typedef CGAL::Segment_3<Rep>					Segment;
typedef Triangulation::Vertex_handle				Vertex_handle;
typedef Triangulation::Edge					Edge;
typedef Triangulation::Facet					Facet;
typedef Triangulation::Cell_handle				Cell_handle;
typedef Triangulation::Finite_vertices_iterator			FVI;
typedef Triangulation::Finite_edges_iterator			FEI;
typedef Triangulation::Finite_facets_iterator			FFI;
typedef Triangulation::Finite_cells_iterator		 	FCI;
typedef Triangulation::All_vertices_iterator		 	AVI;
typedef Triangulation::All_edges_iterator		 	AEI;
typedef Triangulation::All_facets_iterator		 	AFI;
typedef Triangulation::All_cells_iterator		 	ACI;
typedef Triangulation::Facet_circulator				Facet_circulator;
typedef Triangulation::Cell_circulator				Cell_circulator;

// SEGMENT :
// Clusters are to be maintained in a UNION data structure.
 struct cell_cluster
{
	cell_cluster()
		: rep(0), nxt(0), tail(0), in_cluster(false), outside(true), mouth_cnt(0) {}

	cell_cluster(int i)
		: rep(i), nxt(0), tail(i), in_cluster(false), outside(true), mouth_cnt(0) {}

	int find() const
	{
		return rep;
	}

	bool is_in_cluster() const
	{
		return in_cluster;
	}

	// Fields for traversal and ownership.
	int rep;
	cell_cluster* nxt;
	int tail;

	// Fields for attributes.
	bool in_cluster;
	bool outside;
	int mouth_cnt;
	vector<bool> is_big_mouth_vector;
};


#endif
