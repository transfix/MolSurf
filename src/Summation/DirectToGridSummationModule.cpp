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
#include <Summation/DirectToGridSummationModule.h>

using namespace CCVSummationModule;

DirectToGridSummationModule::DirectToGridSummationModule(
    double *centers, double *radii, double *weights, int numberOfCenters,
    Kernel *kernel, OutputGrid *outputGrid)
    : SummationModule(centers, radii, weights, numberOfCenters, kernel) {
  m_UniformOutputGrid = (UniformOutputGrid *)outputGrid;
}

DirectToGridSummationModule::~DirectToGridSummationModule() {}

bool DirectToGridSummationModule::sum() {
  if (!m_Centers || !m_Radii || !m_Weights || !m_Kernel ||
      !m_UniformOutputGrid) {
    return false;
  }
  if (!m_UniformOutputGrid->initialized()) {
    return false;
  }
  int minExtent[3];
  int maxExtent[3];
  minExtent[0] = minExtent[1] = minExtent[2] = 0;
  maxExtent[0] = m_UniformOutputGrid->m_Dimensions[0];
  maxExtent[1] = m_UniformOutputGrid->m_Dimensions[1];
  maxExtent[2] = m_UniformOutputGrid->m_Dimensions[2];
  {
    int i;
    int size = m_UniformOutputGrid->m_Dimensions[0] *
               m_UniformOutputGrid->m_Dimensions[1] *
               m_UniformOutputGrid->m_Dimensions[2];
    for (i = 0; i < size; i++) {
      m_UniformOutputGrid->m_Output[i] = 0;
    }
  }
  int c;
  for (c = 0; c < m_NumberOfCenters; c++) {
    if (m_Kernel->isDecayingKernel()) {
      double lengthOfDecayingKernel = -1;
      if (!m_Kernel->getLength(&lengthOfDecayingKernel, m_Radii[c])) {
        return false;
      }
      if (lengthOfDecayingKernel <= 0) {
        return false;
      }
      double minx = m_Centers[c * 3 + 0] - lengthOfDecayingKernel;
      if (minx < m_UniformOutputGrid->m_Origin[0]) {
        minx = m_UniformOutputGrid->m_Origin[0];
      }
      double miny = m_Centers[c * 3 + 1] - lengthOfDecayingKernel;
      if (miny < m_UniformOutputGrid->m_Origin[1]) {
        miny = m_UniformOutputGrid->m_Origin[1];
      }
      double minz = m_Centers[c * 3 + 2] - lengthOfDecayingKernel;
      if (minz < m_UniformOutputGrid->m_Origin[2]) {
        minz = m_UniformOutputGrid->m_Origin[2];
      }
      minExtent[0] = (int)((minx - m_UniformOutputGrid->m_Origin[0]) /
                           m_UniformOutputGrid->m_Span[0]);
      minExtent[1] = (int)((miny - m_UniformOutputGrid->m_Origin[1]) /
                           m_UniformOutputGrid->m_Span[1]);
      minExtent[2] = (int)((minz - m_UniformOutputGrid->m_Origin[2]) /
                           m_UniformOutputGrid->m_Span[2]);
      maxExtent[0] =
          (int)(ceil((m_Centers[c * 3 + 0] + lengthOfDecayingKernel) /
                     m_UniformOutputGrid->m_Span[0]));
      if (maxExtent[0] >= (int)(m_UniformOutputGrid->m_Dimensions[0])) {
        maxExtent[0] = m_UniformOutputGrid->m_Dimensions[0] - 1;
      }
      maxExtent[1] =
          (int)(ceil((m_Centers[c * 3 + 1] + lengthOfDecayingKernel) /
                     m_UniformOutputGrid->m_Span[1]));
      if (maxExtent[1] >= (int)(m_UniformOutputGrid->m_Dimensions[1])) {
        maxExtent[1] = m_UniformOutputGrid->m_Dimensions[1] - 1;
      }
      maxExtent[2] =
          (int)(ceil((m_Centers[c * 3 + 2] + lengthOfDecayingKernel) /
                     m_UniformOutputGrid->m_Span[2]));
      if (maxExtent[2] >= (int)(m_UniformOutputGrid->m_Dimensions[2])) {
        maxExtent[2] = m_UniformOutputGrid->m_Dimensions[2] - 1;
      }
    }
    int i, j, k;
    for (i = minExtent[0]; i <= maxExtent[0]; i++) {
      for (j = minExtent[1]; j <= maxExtent[1]; j++) {
        for (k = minExtent[2]; k <= maxExtent[2]; k++) {
          int gridPoint = CCVOpenGLMath::TrilinearGrid::xyz2vtx(
              i, j, k, m_UniformOutputGrid->m_Dimensions);
          m_UniformOutputGrid->m_Output[gridPoint] +=
              (float)(sumAtPoint(i, j, k, c));
        }
      }
    }
    printf("Finished %3.1lf%%\r", c / ((float)(m_NumberOfCenters)) * 100);
  }
  return true;
}

double DirectToGridSummationModule::sumAtPoint(int iIndex, int jIndex,
                                               int kIndex, int center) {
  double x = m_UniformOutputGrid->m_Origin[0] +
             iIndex * m_UniformOutputGrid->m_Span[0];
  double y = m_UniformOutputGrid->m_Origin[1] +
             jIndex * m_UniformOutputGrid->m_Span[1];
  double z = m_UniformOutputGrid->m_Origin[2] +
             kIndex * m_UniformOutputGrid->m_Span[2];
  return m_Weights[center] *
         m_Kernel->getFunctionAt(x, y, z, m_Centers[center * 3 + 0],
                                 m_Centers[center * 3 + 1],
                                 m_Centers[center * 3 + 2], m_Radii[center]);
}
