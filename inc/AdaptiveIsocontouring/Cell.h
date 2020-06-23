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
#ifndef _CELL_H_
#define _CELL_H_

#include <Utility/utility.h>
#include <AdaptiveIsocontouring/comp.h>

namespace AdaptiveIsocontouring {

class Cell
{
	private:
		vector<float*> particles;
		Cell* children[8];
		Cell* parent;
		float* func_vals[8];
		float corners[8][3];
		bool hasContour;
		bool hasChildren;
		float mesh_vertex[3];
		// for particle decimation
		float decimate_err;
		float dec_particle[3];
		unsigned int mesh_vert_index;
		int child_index; // stores the index of this child 0-7
		int level;
		bool edges_tagged[12]; //all the edges of cells that has, tag if edge has been checked

	public:
		Cell(float data[8][3], Cell* parent, int ci, int level);
		Cell();
		~Cell();
		float max_X();
		float min_X();
		float max_Y();
		float min_Y();
		float max_Z();
		float min_Z();
		static void setParticlesForCell(Cell* c, vector<float*> in);
		static void findMinMax(float* min, float* max, vector<float*> &in);
		static void findMinMax(float* min, float* max, vector<double*> &in);
};

}

#endif
