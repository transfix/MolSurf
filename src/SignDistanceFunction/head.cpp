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
#include <SignDistanceFunction/common.h>

namespace SDFLibrary
{
	double MAX_DIST;
	int size;
	int total_points, total_triangles, all_verts_touched;
	double minx, miny, minz, maxx, maxy, maxz;
	double TOLERANCE;
	int octree_depth;
	int flipNormals;
	SDFLibrary::triangle* surface;
	SDFLibrary::myVert* vertices;
	SDFLibrary::myPoint* normals;
	SDFLibrary::cell** * sdf;
	SDFLibrary::voxel* values;
	double* distances;
	char* ifname;
	bool* bverts;
	int* queues;
	double minext[3];
	double maxext[3];
	double span[3];
};

void SDFLibrary::init_all_vars()
{
	SDFLibrary::TOLERANCE = 1e-5;
	SDFLibrary::size =64;
	SDFLibrary::flipNormals =0;
	SDFLibrary::ifname= NULL;
	SDFLibrary::surface = NULL;
	SDFLibrary::vertices = NULL;
	SDFLibrary::normals = NULL;
	SDFLibrary::distances = NULL;
	SDFLibrary::sdf = NULL;
	SDFLibrary::values = NULL;
	SDFLibrary::bverts = NULL;
	SDFLibrary::queues = NULL;
	SDFLibrary::minext[0] = SDFLibrary::minext[1] = SDFLibrary::minext[2] = 10000.00;
	SDFLibrary::maxext[0] = SDFLibrary::maxext[1] = SDFLibrary::maxext[2] = -10000.00;
	SDFLibrary::span[0] = SDFLibrary::span[1] = SDFLibrary::span[2] = 1.0;
	SDFLibrary::total_points = SDFLibrary::total_triangles = SDFLibrary::all_verts_touched= 0;
}
