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
#if !defined(AFX_INOUTOCTREE_H__E846DFC4_E4AB_469C_8A3F_AD4E4D54962C__INCLUDED_)
#define AFX_INOUTOCTREE_H__E846DFC4_E4AB_469C_8A3F_AD4E4D54962C__INCLUDED_

#include <Utility/utility.h>

class InOutOctree {
public:
  InOutOctree(int depth, float *min, float *max, float *orig, float *span);
  virtual ~InOutOctree();

  void updateOctree(float cx, float cy, float cz, float r, int iIndex,
                    int jIndex, int kIndex);
  double getClosestDistance(int iIndex, int jIndex, int kIndex);

  int depth;
  float m_Min[3], m_Max[3];
  float m_Orig[3], m_Span[3];
  unsigned int m_Dim[3];
};

#endif
