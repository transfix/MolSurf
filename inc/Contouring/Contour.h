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
#if !defined(AFX_CONTOUR_H__25453561_B42B_46A3_B5AB_15FA0A287BEF__INCLUDED_)
#define AFX_CONTOUR_H__25453561_B42B_46A3_B5AB_15FA0A287BEF__INCLUDED_

#include <Contouring/ContourExtractor.h>
#include <Contouring/ContourGeometry.h>

class Geometry;

namespace CCVContouring
{
	class Contour
	{
		public:
			Contour();
			virtual ~Contour();
			void renderContour(const ContourExtractor& contourExtractor);
			void extract(const ContourExtractor& contourExtractor);
			void addToGeometry(const ContourExtractor& contourExtractor, Geometry* geometry, const CCVOpenGLMath::Matrix& matrix, int& nextVert, int& nextTri);
			int getNumVerts(const ContourExtractor& contourExtractor);
			int getNumTris(const ContourExtractor& contourExtractor);
			void setID(int id);
			int getID();
			void setIsovalue(float isovalue);
			void setSingleColor(float R, float G, float B, bool clobber);
			void resetContour();
			void setWireframeMode(bool state);
			bool useColors();
		protected:
			void setDefaults();
			int m_ID;
			float m_Isovalue;
			float m_R,m_G,m_B;
			bool m_ContourReady;
			ContourGeometry m_ContourGeometry;
	};
};

#endif
