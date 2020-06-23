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
#if !defined(AFX_RAY_H__8760CFFE_A656_460F_A893_4C81FF806264__INCLUDED_)
#define AFX_RAY_H__8760CFFE_A656_460F_A893_4C81FF806264__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Vector.h>

namespace CCVOpenGLMath
{
	class Ray
	{
		public:
			Ray();
			Ray(const Vector& origin, const Vector& dir);
			virtual ~Ray();
			void print();
			Vector getPointOnRay(float t) const;
			float nearestTOnXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			float nearestTOnYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			float nearestTOnZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			Vector nearestPointOnXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			Vector nearestPointOnYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			Vector nearestPointOnZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			float distanceToXAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			float distanceToYAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			float distanceToZAxis(Vector Origin = Vector(0.0, 0.0, 0.0, 1.0)) const;
			bool intersectSphere(Vector center, float radius, Vector* point1, Vector* point2, float* distance1, float* distance2);
			Vector m_Origin;
			Vector m_Dir;
	};
};

#endif
