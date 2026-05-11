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
#if !defined(                                                                  \
    AFX_MARCHINGCUBESBUFFERS_H__2F1981CF_3DB7_402D_B078_AB961D4B6189__INCLUDED_)
#define AFX_MARCHINGCUBESBUFFERS_H__2F1981CF_3DB7_402D_B078_AB961D4B6189__INCLUDED_

namespace CCVContouring {
class MarchingCubesBuffers {
public:
  MarchingCubesBuffers();
  virtual ~MarchingCubesBuffers();
  bool allocateEdgeBuffers(unsigned int width, unsigned int height);
  void destroyEdgeBuffers();
  void swapEdgeBuffers();
  // cached edges to avoid recomputation
  unsigned int *m_EdgeCaches[5];
  unsigned int *m_VertClassifications[2];

protected:
  void setDefaults();
  bool forceAllocateEdgeBuffers(unsigned int width, unsigned int height);
  unsigned int m_AmountAllocated;
};
}; // namespace CCVContouring

#endif
