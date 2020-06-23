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
#include <AdaptiveIsocontouring/comp.h>

#define X_SLICE 0
#define Y_SLICE 1
#define Z_SLICE 2

#define min2(x,y) (x<y) ? (x) : (y)
#define max2(x,y) (x>y) ? (x) : (y)

Vtx cube_pos[8] =
{
	{0,0,0} , {1,0,0} , {1,0,1} , {0,0,1} ,
	{0,1,0} , {1,1,0} , {1,1,1} , {0,1,1}
};
int face_edge[4][2] = { {0,1},{1,3},{3,2},{2,0} };
int computeFaceShoulder(float* face_val, Vtx& face_shoulder,
int slice_direction, float displacement,float isovalue);
void getSliceValue(float* slice_val , float* val , int slice_direction ,
float displacement);

int face_vtx[6][4] =
{
	{0,1,4,5},
	{0,1,3,2},
	{1,5,2,6},
	{4,5,7,6},
	{0,4,3,7},
	{3,2,7,6}
};

float dist2(Vtx a, Vtx b)
{
	float x,y,z;
	x= a.x-b.x;
	y= a.y-b.y;
	return sqrt(x*x + y*y);
}

// val[8] order
// val[0] : (0,0,0)
// val[1] : (1,0,0)
// val[2] : (1,0,1)
// val[3] : (0,0,1)
// val[4] : (0,1,0)
// val[5] : (1,1,0)
// val[6] : (1,1,1)
// val[7] : (0,1,1)
void GetBishoulder(float* val , Vtx& bishoulder , float isovalue)
{
	float face_val[4] , slice_val[4];
	Vtx intersection[4];
	int f;
	int edge_idx=0;
	Vtx face_shoulder1;
	float min,max;
	for(f = 0 ; f < 6 ; f++)
	{
		for(int j = 0 ; j < 4 ; j++)
		{
			face_val[j] = val[face_vtx[f][j]];
		}
		for(int e = 0 ; e < 4 ; e++)
		{
			min=min2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
			max=max2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
			if(min < isovalue && max >= isovalue)
			{
				edge_idx++;
			}
		}
		if(edge_idx > 0)
		{
			break;
		}
	}
	float slope;
	switch(f)
	{
		case 0 :
			computeFaceShoulder(face_val,face_shoulder1,Z_SLICE,0,isovalue);
			getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
			computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
			break;
		case 1 :
			computeFaceShoulder(face_val,face_shoulder1,Y_SLICE,0,isovalue);
			getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
			computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
			break;
		case 2 :
			computeFaceShoulder(face_val,face_shoulder1,X_SLICE,1,isovalue);
			getSliceValue(slice_val,val,Y_SLICE,face_shoulder1.y);
			computeFaceShoulder(slice_val,bishoulder,Y_SLICE,face_shoulder1.y,isovalue);
			break;
		case 3 :
			computeFaceShoulder(face_val,face_shoulder1,Y_SLICE,1,isovalue);
			getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
			computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
			break;
		case 4 :
			computeFaceShoulder(face_val,face_shoulder1,X_SLICE,0,isovalue);
			getSliceValue(slice_val,val,Y_SLICE,face_shoulder1.y);
			computeFaceShoulder(slice_val,bishoulder,Y_SLICE,face_shoulder1.y,isovalue);
			break;
		case 5 :
			computeFaceShoulder(face_val,face_shoulder1,Z_SLICE,1,isovalue);
			getSliceValue(slice_val,val,X_SLICE,face_shoulder1.x);
			computeFaceShoulder(slice_val,bishoulder,X_SLICE,face_shoulder1.x,isovalue);
			break;
	}
}

int computeFaceShoulder(float* face_val, Vtx& face_shoulder, int slice_direction, float displacement,float isovalue)
{
	int e , edge_idx;
	float min,max,ratio;
	Vtx face_pos[4]= {{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
	Vtx intersection[4];
	edge_idx=0;
	for(e = 0 ; e < 4 ; e++)
	{
		min=min2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
		max=max2(face_val[face_edge[e][0]] , face_val[face_edge[e][1]]);
		if(min < isovalue && max >= isovalue)
		{
			ratio = (isovalue - min) / (max - min);
			if(face_val[face_edge[e][0]] < face_val[face_edge[e][1]])
			{
				intersection[edge_idx].x = face_pos[face_edge[e][0]].x * (1 - ratio) + face_pos[face_edge[e][1]].x * ratio;
				intersection[edge_idx].y = face_pos[face_edge[e][0]].y * (1 - ratio) + face_pos[face_edge[e][1]].y * ratio;
			}
			else
			{
				intersection[edge_idx].x = face_pos[face_edge[e][1]].x * (1 - ratio) + face_pos[face_edge[e][0]].x * ratio;
				intersection[edge_idx].y = face_pos[face_edge[e][1]].y * (1 - ratio) + face_pos[face_edge[e][0]].y * ratio;
			}
			edge_idx++;
		}
	}
	if(edge_idx!=2)
	{
		face_shoulder.x=0.5;
		face_shoulder.y=0.5;
		face_shoulder.z=0.5;
		return false;
	}
	double p,q,r,s;
	float tface_val[4];
	tface_val[0]=face_val[0]-isovalue;
	tface_val[1]=face_val[1]-isovalue;
	tface_val[2]=face_val[2]-isovalue;
	tface_val[3]=face_val[3]-isovalue;
	double slope_sign;
	p=tface_val[0];
	q=-tface_val[0]+tface_val[1];
	r=-tface_val[0]+tface_val[2];
	s=tface_val[0]-tface_val[1]-tface_val[2]+tface_val[3];
	double slope;
	Vtx middle_point,face_candidate1,face_candidate2,face_candidate;
	slope=(intersection[1].y-intersection[0].y)/(intersection[1].x-intersection[0].x);
	if(slope>0)
	{
		slope_sign=1;
	}
	else
	{
		slope_sign=-1;
	}
	middle_point.x = (intersection[0].x + intersection[1].x)/2.;
	middle_point.y = (intersection[0].y + intersection[1].y)/2.;
	if((s<0.01 && s>-0.01) ||
	(slope<0.01 && slope>-0.01) ||
	(slope<-1000 || slope>1000))
	{
		face_candidate=middle_point;
	}
	else
	{
		face_candidate1.x=(float)((-r+sqrt(-(q*r-s*p)/slope))/s);
		face_candidate1.y=(float)((-q-slope_sign*sqrt(-(q*r-s*p)*slope))/s);
		face_candidate2.x=(float)((-r-sqrt(-(q*r-s*p)/slope))/s);
		face_candidate2.y=(float)((-q+slope_sign*sqrt(-(q*r-s*p)*slope))/s);
		if(face_candidate1.x>=0 && face_candidate1.x<=1 && face_candidate1.y>=0 && face_candidate1.y<=1)
		{
			face_candidate=face_candidate1;
		}
		else if(face_candidate2.x>=0 && face_candidate2.x<=1 && face_candidate2.y>=0 && face_candidate2.y<=1)
		{
			face_candidate=face_candidate2;
		}
		else
		{
			if(dist2(middle_point,face_candidate1)<dist2(middle_point,face_candidate2))
			{
				face_candidate=face_candidate1;
			}
			else
			{
				face_candidate=face_candidate2;
			}
		}
	}
	if(face_candidate.x<0)
	{
		face_candidate.x=0;
		face_candidate=middle_point;
	}
	if(face_candidate.x>1)
	{
		face_candidate.x=1;
		face_candidate=middle_point;
	}
	if(face_candidate.y<0)
	{
		face_candidate.y=0;
		face_candidate=middle_point;
	}
	if(face_candidate.y>1)
	{
		face_candidate.y=1;
		face_candidate=middle_point;
	}
	switch(slice_direction)
	{
		case X_SLICE :
			face_shoulder.x=displacement;
			face_shoulder.y=face_candidate.x;
			face_shoulder.z=face_candidate.y;
			break;
		case Y_SLICE :
			face_shoulder.x=face_candidate.x;
			face_shoulder.y=displacement;
			face_shoulder.z=face_candidate.y;
			break;
		case Z_SLICE :
			face_shoulder.x=face_candidate.x;
			face_shoulder.y=face_candidate.y;
			face_shoulder.z=displacement;
			break;
	}
	return true;
}

void getSliceValue(float* slice_val , float* val , int slice_direction , float displacement)
{
	int x_edge[4][2] = { {0,1},{4,5},{3,2},{7,6} };
	int y_edge[4][2] = { {0,4},{1,5},{3,7},{2,6} };
	int z_edge[4][2] = { {0,3},{1,2},{4,7},{5,6} };
	int i;
	switch(slice_direction)
	{
		case X_SLICE :
			for(i=0; i<4; i++)
			{
				slice_val[i] = val[x_edge[i][0]]*(1-displacement) + val[x_edge[i][1]]*displacement;
			}
			break;
		case Y_SLICE :
			for(i=0; i<4; i++)
			{
				slice_val[i] = val[y_edge[i][0]]*(1-displacement) + val[y_edge[i][1]]*displacement;
			}
			break;
		case Z_SLICE :
			for(i=0; i<4; i++)
			{
				slice_val[i] = val[z_edge[i][0]]*(1-displacement) + val[z_edge[i][1]]*displacement;
			}
			break;
	}
}

// input : normalized_bishoulder - the location of the computed bishoulder point
//         which is assumed to be in a cell (0,0,0)-(1,1,1)
//         x,y,z - minimum coordinate of a cell which has the bishoulder point
// output: bishoulder - actual coordinate of the bishoulder point
void Norm2Read_Bishoulder(Vtx normalized_bishoulder,float* bishoulder,float x,float y,float z,float cell_size)
{
	bishoulder[0] = x + cell_size*normalized_bishoulder.x;
	bishoulder[1] = y + cell_size*normalized_bishoulder.y;
	bishoulder[2] = z + cell_size*normalized_bishoulder.z;
}

void Norm2Read_Bishoulder2(Vtx normalized_bishoulder,float* bishoulder, float xyz[3] , float cell_size[3])
{
	bishoulder[0] = xyz[0] + cell_size[0]*normalized_bishoulder.x;
	bishoulder[1] = xyz[1] + cell_size[1]*normalized_bishoulder.y;
	bishoulder[2] = xyz[2] + cell_size[2]*normalized_bishoulder.z;
}

bool is_ambiguous_case(float* val,float iso_val,int Case)
{
	CellQueue vtx_queue;
	int vtx_connectivity[8][3]= {{1,3,4} , {0,2,5} , {1,3,6} , {0,2,7} , {0,5,7} , {1,4,6} , {2,5,7} , {3,4,6}};
	int vtx_idx=-1;
	int code[8] = {0,0,0,0,0,0,0,0} ;
	for(int i = 0 ; i < 8 ; i++)
	{
		if(Case==0)
		{
			if(val[i] >= iso_val)
			{
				code[i] = 1;
				vtx_idx=i;
			}
		}
		else
		{
			if(val[i] < iso_val)
			{
				code[i] = 1;
				vtx_idx=i;
			}
		}
	}
	if(vtx_idx==-1)
	{
		return false;
	}
	else
	{
		vtx_queue.Add(vtx_idx);
		code[vtx_idx]=0;
		while(vtx_queue.Get(vtx_idx)>=0)
		{
			for(int i = 0 ; i < 3 ; i++)
			{
				if(code[vtx_connectivity[vtx_idx][i]]==1)
				{
					vtx_queue.Add(vtx_connectivity[vtx_idx][i]);
					code[vtx_connectivity[vtx_idx][i]]=0;
				}
			}
		}
	}
	for(int i = 0 ; i < 8 ; i++)
	{
		if(code[i]==1)
		{
			return true;
		}
	}
	return false;
}

// Decides if an isosurface in the cell specified by its value is ambiguous.
// if a cell has an ambiguous cell, you must refine the cell until isosurfaces in the cells are not ambiguous
bool is_ambiguous(float* val,float iso_val)
{
	bool flag_val_case1,flag_val_case2;
	flag_val_case1=is_ambiguous_case(val,iso_val,0);
	flag_val_case2=is_ambiguous_case(val,iso_val,1);
	if(flag_val_case1==true || flag_val_case2==true)
	{
		return true;
	}
	else
	{
		return false;
	}
}
