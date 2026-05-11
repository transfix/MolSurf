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
#include <Summation/UniformOutputGrid.h>

using namespace CCVSummationModule;

UniformOutputGrid::UniformOutputGrid(float *output, float *origin, float *span,
                                     unsigned int *dimensions) {
  m_Output = output;
  m_Origin = origin;
  m_Span = span;
  m_Dimensions = dimensions;
}

UniformOutputGrid::~UniformOutputGrid() {}

bool UniformOutputGrid::initialized() {
  if (!m_Output || !m_Origin || !m_Span || !m_Dimensions) {
    return false;
  }
  return true;
}
