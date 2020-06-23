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
// the head file of fatsfs 

typedef struct Face1
{
	int         Index[3];             /* the indexs of a triangle-increase order*/
	int         Orien;                /* 1--orientaioned, othwise no            */
	int         AdjTri[3];            /* neighbor triangles index               */
	float       Center[3];            /* the center of the triangle             */
	float       Normal[3];            /* face normal                            */
	int         whichnml[3];          /* Index[0] normal index                  */
	int         Tou;                  /* the begin index of the inside points   */
	int         Wei;                  /* the end index+1 of the inside points   */
	int         State;                /* 1--removed; 2-- new, 0--old;           */
} *Face;

typedef struct Neighbor1
{
	int*         NB;
} *Neighbor;

typedef struct VertixNormal1
{
	float   Point[3];
	float   Normal[3];
	float*   Hessian;                 /* the hesssian at the vertices           */
	int     NumFaces;
	float*   SharpNormal;             /* the normal for sharp vertex           */
} *VertixNormal;

typedef struct Triangle1
{
	float       Vertex1[3];
	float       Vertex2[3];
	float       Vertex3[3];
	float       Normal1[3];
	float       Normal2[3];
	float       Normal3[3];
} *Triangle;
