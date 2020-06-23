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
#if !defined(AFX_MULTICONTOUR_H__9DC23294_BACA_494B_AF88_C41B94D54D4F__INCLUDED_)
#define AFX_MULTICONTOUR_H__9DC23294_BACA_494B_AF88_C41B94D54D4F__INCLUDED_

#include <Contouring/ContourExtractor.h>
#include <Contouring/Contour.h>
#include <Contouring/RGBAExtractor.h>
#include <Contouring/SingleExtractor.h>
#include <Geometry/Geometry.h>
#include <UsefulMath/Matrix.h>

struct iPoly;
class Geometry;

namespace CCVContouring
{
	class MultiContour
	{
		public:
			MultiContour();
			virtual ~MultiContour();
			void setData(unsigned char* data,
						 unsigned int width, unsigned int height, unsigned int depth,
						 double aspectX, double aspectY, double aspectZ,
						 double subMinX, double subMinY, double subMinZ,
						 double subMaxX, double subMaxY, double subMaxZ,
						 double minX, double minY, double minZ,
						 double maxX, double maxY, double maxZ);
			void setData(unsigned char* data, unsigned char* red, unsigned char* green, unsigned char* blue,
						 unsigned int width, unsigned int height, unsigned int depth,
						 double aspectX, double aspectY, double aspectZ,
						 double subMinX, double subMinY, double subMinZ,
						 double subMaxX, double subMaxY, double subMaxZ,
						 double minX, double minY, double minZ,
						 double maxX, double maxY, double maxZ);
			void addContour(int ID, float isovalue, float R, float G, float B);
			void removeContour(int ID);
			void removeAll();
			void setIsovalue(int ID, float isovalue);
			void setColor(int ID, float R, float G, float B);
			virtual void setWireframeMode(bool state);
			void renderContours() const;
			void forceExtraction() const;
			virtual bool render(double a, double b, double c, double d);
			Geometry* getGeometry() const ;
			int getNumVerts() const ;
			int getNumTris() const ;
		private:
			class MultiContourNode
			{
				public:
					MultiContourNode(int ID, float isovalue, float R, float G, float B, bool wire, MultiContourNode* next) : m_Next(next)
					{
						m_Contour.setID(ID);
						m_Contour.setIsovalue(isovalue);
						m_Contour.setWireframeMode(wire);
						m_Contour.setSingleColor(R,G,B, true);
					};
					virtual ~MultiContourNode()
					{
						delete m_Next;
					};
					MultiContourNode* m_Next;
					Contour m_Contour;
			};
		protected:
			void setDefaults();
			void resetContours();
			SingleExtractor m_SingleExtractor;
			RGBAExtractor m_RGBAExtractor;
			ContourExtractor* m_ContourExtractor;
			MultiContourNode* m_Head;
			bool m_DataLoaded;
			bool m_WireframeRender;
			// matrix used to scale and center the contour into the right place
			CCVOpenGLMath::Matrix* m_Matrix;
			// matrix used to align the data to the correct space for saving
			CCVOpenGLMath::Matrix* m_SaveMatrix;
	};

};

#endif
