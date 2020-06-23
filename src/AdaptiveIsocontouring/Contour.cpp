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
#include <AdaptiveIsocontouring/Contour.h>

using namespace AdaptiveIsocontouring;

Contour::Contour()
{
	bi_vert_out = NULL;
	bi_vert_in = NULL;
	vd_fine = vd_coarse = NULL;
}

Contour::~Contour()
{
	if(bi_vert_out!=NULL)
	{
		delete [] bi_vert_out;
	}
	if(bi_vert_in!=NULL)
	{
		delete [] bi_vert_in;
	}
}

void Contour::initialize(int fine_level_, int coarse_level_, float* vd_coarse_, float* min, float* max, LittleGeoframe* geo_)
{
	fine_level = fine_level_;
	coarse_level = coarse_level_;
	vd_coarse = vd_coarse_;
	geo = geo_;
	int coarse_size = (1 << coarse_level) + 1;
	copyVect(min, in_min, 3);
	copyVect(max, in_max, 3);
	for(int j=0; j<3; j++)
	{
		dim[j] = (in_max[j] - in_min[j]) / (coarse_size -1);
	}
}

void Contour::initialize(int fine_level, int coarse_level, float* vol, float* min, float* max, LittleGeoframe& geo)
{
	initialize(fine_level, coarse_level, vol, min, max, &geo);
}

void Contour::initialize(int fine_level_, int coarse_level_, double* vol, LittleGeoframe& geo_)
{
	fine_level = fine_level_;
	coarse_level = coarse_level_;
	geo = &geo_;
	int coarse_size = (1 << coarse_level) + 1;
	if(vd_coarse != NULL)
	{
		delete [] vd_coarse;
	}
	vd_coarse = new float[coarse_size * coarse_size * coarse_size];
	for(int j=0; j<coarse_size ; j++)
	{
		vd_coarse[j] = vol[j];
	}
}


void Contour::set_isovalue_outer(float iso)
{
	iso_val_outer = iso;
}

void Contour::set_isovalue_inner(float iso)
{
	iso_val_inner = iso;
}

void Contour::put_bishoulder()
{
	geo->Clear2();
	int c_size = (1<< coarse_level);
	int c_osize = c_size + 1;
	int plus_x = 1;
	int plus_y = c_osize ;
	int plus_z = c_osize * c_osize;
	float all_mid[3] = {in_min[0] + dim[0] / 2, in_min[1] + dim[1] / 2, in_min[2] + dim[2] / 2 };
	float tmp[8];
	float bi_vert_tmp[3];
	for(int j=0; j<c_size * c_size * c_size ; j++)
	{
		Vtx holder;
		int c_x = j % c_size;
		int c_y = (j / c_size) % c_size;
		int c_z = (j / c_size) / c_size;
		int vd_ind = c_x + c_y * c_osize + c_z * c_osize * c_osize;
		tmp[0] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind + plus_x;   // 100
		tmp[1] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind + plus_z;   // 101
		tmp[2] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind - plus_x;   // 001
		tmp[3] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind - plus_z + plus_y;   // 010
		tmp[4] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind + plus_x;   // 110
		tmp[5] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind + plus_z;   // 111
		tmp[6] = (float) vd_coarse[vd_ind];
		vd_ind = vd_ind - plus_x;   // 011
		tmp[7] = (float) vd_coarse[vd_ind];
		bool all_less = true;
		bool all_greater = true;
		for(int j=0; j<8; j++)
		{
			if(tmp[j] > iso_val_outer)
			{
				all_less = false;
			}
			if(tmp[j] < iso_val_outer)
			{
				all_greater = false;
			}
		}
		if((all_less || all_greater))
		{
			bi_vert_tmp[0] = in_min[0] + ((float) c_x +.5) * dim[0] ;
			bi_vert_tmp[1] = in_min[1] + ((float) c_y +.5) * dim[1] ;
			bi_vert_tmp[2] = in_min[2] + ((float) c_z +.5) * dim[2] ;
		}
		else
		{
			GetBishoulder(tmp, holder, iso_val_outer);
			bi_vert_tmp[0] =  in_min[0] + ((float) c_x +  holder.x) * dim[0];
			bi_vert_tmp[1] =  in_min[1] + ((float) c_y + holder.y) * dim[1];
			bi_vert_tmp[2] =  in_min[2] + ((float) c_z +  holder.z) * dim[2];
		}
		float c[4] = {.5, .9, .8, 1.0};
		float norm[3] = {0, 0, 0};
		geo->addVert(bi_vert_tmp, norm, c);
	}
}

bool check_iso_range(float iso_val, float val1, float val2)
{
	return (iso_val <val2 && iso_val > val1);
}

void Contour::polygonize()
{
	// try all X-edges, note that iso-value has to be defined
	cout<<"called Contour::polygonize, geo_size: "<<geo->getNumVerts()<<" isovalue: "<<iso_val_outer<<endl;
	int c_size = 1<<coarse_level;
	int c_osize = c_size +1;
	int ind1, ind2;
	unsigned int quad[4];
	for(int y=1; y<c_osize-1; y++)
	{
		for(int z=1; z<c_osize-1; z++)
		{
			for(int x=0; x<c_size; x++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + 1 + y * c_osize + z * c_osize * c_osize;
				if(check_iso_range(iso_val_outer, vd_coarse[ind1], vd_coarse[ind2]))
				{
					// if this is a sign change edge, polgonize;
					quad[0] = x + (y-1) * c_size + (z-1) * c_size * c_size;
					quad[1] = x + (y-1) * c_size + (z) * c_size * c_size;
					quad[2] = x + (y) * c_size + (z) * c_size * c_size;
					quad[3] = x + (y) * c_size + (z-1) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_outer, vd_coarse[ind2], vd_coarse[ind1]))
				{
					quad[0] = x + (y-1) * c_size + (z-1) * c_size * c_size;
					quad[1] = x + (y) * c_size + (z-1) * c_size * c_size;
					quad[2] = x + (y) * c_size + (z) * c_size * c_size;
					quad[3] = x + (y-1) * c_size + (z) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	for(int x=1; x<c_osize-1; x++)
	{
		for(int z=1; z<c_osize-1; z++)
		{
			for(int y=0; y<c_size; y++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + (y+1) * c_osize + z * c_osize * c_osize;
				if(check_iso_range(iso_val_outer, vd_coarse[ind1], vd_coarse[ind2]))
				{
					quad[0] = x-1 + y * c_size + (z-1) * c_size * c_size;
					quad[1] = (x) + y * c_size + (z-1) * c_size * c_size;
					quad[2] = x + y * c_size + z * c_size * c_size;
					quad[3] = (x-1) + y * c_size + (z) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_outer, vd_coarse[ind2], vd_coarse[ind1]))
				{
					quad[0] = x-1 + y * c_size + (z-1) * c_size * c_size;
					quad[1] = (x-1) + y * c_size + (z) * c_size * c_size;
					quad[2] = x + y * c_size + z * c_size * c_size;
					quad[3] = x + y * c_size + (z-1) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	for(int y=1; y<c_osize-1; y++)
	{
		for(int x=1; x<c_osize-1; x++)
		{
			for(int z=0; z<c_size; z++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + y * c_osize + (z+1) * c_osize * c_osize;
				if(check_iso_range(iso_val_outer, vd_coarse[ind1], vd_coarse[ind2]))
				{
					quad[0] = x-1 + (y-1) * c_size + z * c_size * c_size;
					quad[1] = (x-1) + (y) * c_size + z * c_size * c_size;
					quad[2] = x + y * c_size + z * c_size * c_size;
					quad[3] = x + (y-1) * c_size + z * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_outer, vd_coarse[ind2], vd_coarse[ind1]))
				{
					quad[0] = x-1 + (y-1) * c_size + z * c_size * c_size;
					quad[1] = (x) + (y-1) * c_size + z * c_size * c_size;
					quad[2] = x + y * c_size + z * c_size * c_size;
					quad[3] = (x-1) + (y) * c_size + z * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	num_outer_quads = geo->getNumQuads();
}


void Contour::update_quad_render()
{
	int c_size = 1<< coarse_level;
	bool render_quad = true;
	for(int j=0; j<geo->getNumQuads(); j++)
	{
		render_quad = true;
		for(int k=0; k<4; k++)
		{
			int vert = geo->getQuad(j,k);
			int loc[3] = { vert % c_size, vert / c_size % c_size, vert / c_size / c_size};
			if(loc[0] > box_rel_min[0]+1 && loc[0] < box_rel_max[0]-1 &&
					loc[1] > box_rel_min[1]+1 && loc[1] < box_rel_max[1]-1 &&
					loc[2] > box_rel_min[2]+1 && loc[2] < box_rel_max[2]-1)
			{
				geo->setQuadRender(j, false);
				render_quad = false;
				break;
			}
		}
		if(render_quad)
		{
			geo->setQuadRender(j, true);
		}
	}
}

void Contour::initialize_redbox(int resol, float* b_min, float* b_max, float* sum_in)
{
	vd_fine = sum_in;
	box_resol = resol; // assume that it is uniform
	copyVect(b_min, box_min, 3);
	copyVect(b_max, box_max, 3);
	int c_size = 1 << coarse_level;
	for(int j=0; j<3; j++)
	{
		box_dim[j] = (box_max[j] - box_min[j]) / (box_resol);
		box_rel_min[j] = (int) floor((double)(b_min[j] - in_min[j]) / dim[j]);  // the coarse index where the red_box begins
		box_rel_max[j] = (int) floor((double)(b_max[j] - in_min[j]) / dim[j]);   // the coarse index where the red_box ends
		if(box_rel_min[j] < 0)
		{
			box_rel_min[j] =0;
		}
		else if(box_rel_min[j] >= c_size)
		{
			box_rel_min[j] = c_size -1;
		}
		if(box_rel_max[j] < 0)
		{
			box_rel_max[j] = 0;
		}
		else if(box_rel_max[j] >= c_size)
		{
			box_rel_max[j] = c_size -1;
		}
	}
	vd_fine = sum_in;
}

void Contour::polygonize_redbox()
{
	// try all X-edges, note that iso-value has to be defined
	int c_size = box_resol - 1;
	int c_osize = c_size +1;
	int index_offset = 1<<coarse_level;
	index_offset = (index_offset * index_offset * index_offset);
	geo->setNumQuads(num_outer_quads);
	int ind1, ind2;
	unsigned int quad[4];
	for(int y=1; y<c_osize-1; y++)
	{
		for(int z=1; z<c_osize-1; z++)
		{
			for(int x=0; x<c_size; x++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + 1 + y * c_osize + z * c_osize * c_osize;
				for(int j=0; j< 4; j++)
				{
					quad[j] = index_offset;
				}
				if(check_iso_range(iso_val_inner, vd_fine[ind1], vd_fine[ind2]))
				{
					// if this is a sign change edge, polgonize;
					quad[0] += x + (y-1) * c_size + (z-1) * c_size * c_size;
					quad[1] += x + (y-1) * c_size + (z) * c_size * c_size;
					quad[2] += x + (y) * c_size + (z) * c_size * c_size;
					quad[3] += x + (y) * c_size + (z-1) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_inner, vd_fine[ind2], vd_fine[ind1]))
				{
					quad[0] += x + (y-1) * c_size + (z-1) * c_size * c_size;
					quad[1] += x + (y) * c_size + (z-1) * c_size * c_size;
					quad[2] += x + (y) * c_size + (z) * c_size * c_size;
					quad[3] += x + (y-1) * c_size + (z) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	for(int x=1; x<c_osize-1; x++)
	{
		for(int z=1; z<c_osize-1; z++)
		{
			for(int y=0; y<c_size; y++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + (y+1) * c_osize + z * c_osize * c_osize;
				for(int j=0; j< 4; j++)
				{
					quad[j] = index_offset;
				}
				if(check_iso_range(iso_val_inner, vd_fine[ind1], vd_fine[ind2]))
				{
					quad[0] += x-1 + y * c_size + (z-1) * c_size * c_size;
					quad[1] += (x) + y * c_size + (z-1) * c_size * c_size;
					quad[2] += x + y * c_size + z * c_size * c_size;
					quad[3] += (x-1) + y * c_size + (z) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_inner, vd_fine[ind2], vd_fine[ind1]))
				{
					quad[0] += x-1 + y * c_size + (z-1) * c_size * c_size;
					quad[1] += (x-1) + y * c_size + (z) * c_size * c_size;
					quad[2] += x + y * c_size + z * c_size * c_size;
					quad[3] += x + y * c_size + (z-1) * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	for(int y=1; y<c_osize-1; y++)
	{
		for(int x=1; x<c_osize-1; x++)
		{
			for(int z=0; z<c_size; z++)
			{
				ind1 = x + y * c_osize + z * c_osize * c_osize;
				ind2 = x + y * c_osize + (z+1) * c_osize * c_osize;
				for(int j=0; j< 4; j++)
				{
					quad[j] = index_offset;
				}
				if(check_iso_range(iso_val_inner, vd_fine[ind1], vd_fine[ind2]))
				{
					quad[3] += x-1 + (y-1) * c_size + z * c_size * c_size;
					quad[2] += (x-1) + (y) * c_size + z * c_size * c_size;
					quad[1] += x + y * c_size + z * c_size * c_size;
					quad[0] += x + (y-1) * c_size + z * c_size * c_size;
					geo->addQuad(quad, 4);
				}
				else if(check_iso_range(iso_val_inner, vd_fine[ind2], vd_fine[ind1]))
				{
					quad[0] += x-1 + (y-1) * c_size + z * c_size * c_size;
					quad[1] += (x) + (y-1) * c_size + z * c_size * c_size;
					quad[2] += x + y * c_size + z * c_size * c_size;
					quad[3] += (x-1) + (y) * c_size + z * c_size * c_size;
					geo->addQuad(quad, 4);
				}
			}
		}
	}
	geo->calculateQuadNormals();
	update_quad_render();
}

void Contour::put_bishoulder_redbox()
{
	int c_size = box_resol - 1;
	int c_osize = c_size + 1;
	int plus_x = 1;
	int plus_y = c_osize ;
	int plus_z = c_osize * c_osize;
	float all_mid[3] = {in_min[0] + dim[0] / 2, in_min[1] + dim[1] / 2, in_min[2] + dim[2] / 2};
	float tmp[8];
	float bi_vert_tmp[3];
	int sizetmp = 1<<coarse_level ;
	geo->setNumVerts(sizetmp * sizetmp * sizetmp);
	for(int j=0; j<c_size * c_size * c_size ; j++)
	{
		Vtx holder;
		int c_x = j % c_size;
		int c_y = (j / c_size) % c_size;
		int c_z = (j / c_size) / c_size;
		int vd_ind = c_x + c_y * c_osize + c_z * c_osize * c_osize;
		tmp[0] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind + plus_x;   // 100
		tmp[1] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind + plus_z;   // 101
		tmp[2] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind - plus_x;   // 001
		tmp[3] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind - plus_z + plus_y;   // 010
		tmp[4] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind + plus_x;   // 110
		tmp[5] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind + plus_z;   // 111
		tmp[6] = (float) vd_fine[vd_ind];
		vd_ind = vd_ind - plus_x;   // 011
		tmp[7] = (float) vd_fine[vd_ind];
		bool all_less = true;
		bool all_greater = true;
		{
			for(int k=0; k<8; k++)
			{
				if(tmp[k] > iso_val_inner)
				{
					all_less = false;
				}
				if(tmp[k] < iso_val_inner)
				{
					all_greater = false;
				}
			}
		}
		if((all_less || all_greater))
		{
			bi_vert_tmp[0] = box_min[0] + ((float) c_x + .5) * box_dim[0];
			bi_vert_tmp[1] = box_min[1] + ((float) c_y + .5) * box_dim[1];
			bi_vert_tmp[2] = box_min[2] + ((float) c_z + .5) * box_dim[2];
		}
		else
		{
			GetBishoulder(tmp, holder, iso_val_inner);
			bi_vert_tmp[0] =  box_min[0] + ((float) c_x +  holder.x) * box_dim[0];
			bi_vert_tmp[1] =  box_min[1] + ((float) c_y + holder.y) * box_dim[1];
			bi_vert_tmp[2] =  box_min[2] + ((float) c_z +  holder.z) * box_dim[2];
		}
		float c[4] = {.9, .5, .2, .6};
		float norm[3] = {0, 0, 0};
		geo->addVert(bi_vert_tmp, norm, c);
	}
}
