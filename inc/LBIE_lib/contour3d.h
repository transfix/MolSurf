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
#ifndef __CONTOUR3D_H__
#define __CONTOUR3D_H__

class Contour3d
{
	private:
		int nvert, ntri, nquad;
		int vsize,tsize;
		float(*vert)[3];
		unsigned int (*tri)[3];
		unsigned int (*quad)[4];


	public :
		Contour3d();
		~Contour3d();

		void Clear();
		int getNTri(void)
		{
			return ntri;
		}
		int getNQuad(void)
		{
			return nquad;
		}
		int getNVert(void)
		{
			return nvert;
		}
		int AddTri(unsigned int v1,unsigned int v2,unsigned int v3);
		int AddVert(float v_pos[3], float norm[3]);
		int center_vtx(int v1,int v2,int v3);
		void display();
};
#endif //__CONTOUR3D_H__
