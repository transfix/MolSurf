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
#ifndef CCV_SDF_HEAD_H
#define CCV_SDF_HEAD_H

#include <Utility/utility.h>

namespace SDFLibrary
{

#define MAX_TRIS_PER_VERT 100

	typedef struct _Pt_
	{
		double x;
		double y;
		double z;
		char isNull;

	}  myPoint;

	typedef struct _Vt_
	{
		double x;
		double y;
		double z;
		char isNull;

		int tris[MAX_TRIS_PER_VERT]; //not more than MAX_TRIS_PER_VERT triangles can share a vertex.
		int trisUsed; //elements used in the above array.

	}  myVert;


	typedef struct _tri_
	{
		int v1;
		int v2;
		int v3;

		int type; // default = -1; done =1.	wrong =3;
	} triangle;

	typedef struct listnodedef
	{
		int index;	//index of the triangle
		struct listnodedef* next;
	} listnode;

	typedef struct nodedef
	{
		char useful;	//  0 - no triangles in it	; 1 - there are triangles in it
		char type;		//	0 - interior node		; 1 - leaf node, containing triangles
		long int no;
		listnode* tindex;
	} cell;

	typedef struct
	{
		double ox;
		double oy;
		double oz;

		double dx;
		double dy;
		double dz;
	} ray;

	typedef struct _voxel_
	{
		float value;
		signed char signe;  //-1 = inside,		1 = outside
		bool processed;		// 1 = propagated distance FROM here. 0 = not
		int closestV;		//the closest triangle on the surface
	} voxel;


}; //namespace SDFLibrary

#endif
