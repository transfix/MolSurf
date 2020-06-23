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
#include <PocketTunnel/pocket_tunnel.h>

namespace PocketTunnel
{
	// Constants
	// cocone
	const double DEFAULT_ANGLE = M_PI / 8.0;      // Half of the co-cone angle.
	const double DEFAULT_SHARP = 2 * M_PI / 3.0;  // Angle of sharp edges.
	const double DEFAULT_RATIO = 1.2 * 1.2;       // Squared thinness factor.
	const double DEFAULT_FLAT  = M_PI / 3.0;      // Angle for flatness Test
	// robust cocone 
	const double DEFAULT_BIGBALL_RATIO  = 4.*4.;  // parameter to choose big balls.
	const double DEFAULT_THETA_IF_d  = 5.0;       // parameter for infinite-finite deep intersection.
	const double DEFAULT_THETA_FF_d  = 10.0;      // parameter for finite-finite deep intersection.
	const double DEFAULT_THETA_II_d  = 30.0;      // parameter for infinite-infinite deep intersection.
	int pocket_tunnel_fromsurf(Geometry* molsurf,  // input surface data.
				Geometry** PTV, Geometry** curated_surf, // outputs.
				const int& opc, const int& otc) // parameters.
	{
		map<int, cell_cluster> cluster_set;
		// robust cocone parameters.
		double bb_ratio = DEFAULT_BIGBALL_RATIO;
		double theta_ff = M_PI/180.0*DEFAULT_THETA_FF_d;
		double theta_if = M_PI/180.0*DEFAULT_THETA_IF_d;
		vector<Point> pts_list;
		for(int i = 0; i < molsurf->m_NumTriVerts*3;)
		{
			float x = molsurf->m_TriVerts[i++],
				  y = molsurf->m_TriVerts[i++],
				  z = molsurf->m_TriVerts[i++];
			pts_list.push_back(Point(x,y,z));
		}
		CGAL::Timer timer;
		timer.start();
		//cerr << "Delaunay ";
		Triangulation triang;
		triang.insert(pts_list.begin(), pts_list.end());
		//cerr << "done." << endl;
		//cerr << "Time: " << timer.time() << endl;
		timer.reset();
		// Initialization of all the required fields
		//cerr << "Initialization ";
		initialize(triang);
		//cerr << ".";
		// compute voronoi vertex
		compute_voronoi_vertex_and_cell_radius(triang);
		//cerr << ". done." << endl;
		//cerr << "Time: " << timer.time() << endl;
		timer.reset();
		// Surface Reconstruction using Tight Cocone
		//cerr << "Surface Reconstruction ";
		tcocone(DEFAULT_ANGLE, DEFAULT_SHARP, DEFAULT_FLAT, DEFAULT_RATIO, triang);
		//cerr << " done." << endl;
		//cerr << "Time: " << timer.time() << endl;
		timer.reset();
		//cerr << "Computing S_MAX ";
		double mr = 1.3; // not used in this routine.
		vector<int> sorted_smax_index_vector = compute_smax(triang, cluster_set, mr);
		//cerr << " done." << endl;
		// detect pocket, tunnel, void.
		//cerr << "Computing Pocket-Tunnels ";
		detect_handle(triang, cluster_set);
		//cerr << " done." << endl;
		//cerr << "Time: " << timer.time() << endl;
		timer.reset();

		// convert the handles into rawc geometries to be viewed by TexMol.
		convert_pocket_tunnel_to_rawc_geometry(PTV, triang, cluster_set, sorted_smax_index_vector, opc, otc);
		convert_curated_surf_to_raw_geometry(curated_surf, triang, cluster_set,sorted_smax_index_vector, opc, otc);
		return 0;
	}
}
