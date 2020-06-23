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
#if !defined(AFX_EXTENTS_H__541859B1_96FC_453D_8524_1D6789C4F864__INCLUDED_)
#define AFX_EXTENTS_H__541859B1_96FC_453D_8524_1D6789C4F864__INCLUDED_

#include <UsefulMath/Vector.h>

using namespace CCVOpenGLMath;

class Extents
{
	public:
		Extents();
		Extents(
			double xMin, double xMax,
			double yMin, double yMax,
			double zMin, double zMax
		);
		virtual ~Extents();

		void setExtents(
			double xMin, double xMax,
			double yMin, double yMax,
			double zMin, double zMax
		);
		Vector getOrigin() const;
		void setOrigin(Vector vector, const Extents& boundaryExtents);
		void move(const Vector& vector);
		bool withinCube(const Vector& vector) const;
		void clampTo(const Extents& boundaryExtents);
		double getXMin() const;
		double getYMin() const;
		double getZMin() const;
		double getXMax() const;
		double getYMax() const;
		double getZMax() const;
		void setXMin(double xMin);
		void setYMin(double yMin);
		void setZMin(double zMin);
		void setXMax(double xMax);
		void setYMax(double yMax);
		void setZMax(double zMax);

	protected:
		double m_XMin, m_XMax;
		double m_YMin, m_YMax;
		double m_ZMin, m_ZMax;
};

#endif
