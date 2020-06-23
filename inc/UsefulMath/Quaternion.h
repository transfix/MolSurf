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
#if !defined(AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_)
#define AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Tuple.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Ray.h>
#include <UsefulMath/Matrix.h>

namespace CCVOpenGLMath
{
	class Vector;
	class Matrix;
	class Ray;

	class Quaternion : public Tuple
	{
		public:
			Quaternion();
			virtual ~Quaternion();
			Quaternion(const Quaternion& copy);
			Quaternion& operator=(const Quaternion& copy);
			Quaternion(float w, float x, float y, float z);
			Quaternion& set(float w, float x, float y, float z);
			Quaternion& set(float* array);
			Quaternion& set(const Quaternion& copy);
			Quaternion operator*(const Quaternion& quat) const;
			Quaternion operator*(float scalar) const;
			Quaternion& operator*=(float scalar);
			Quaternion operator/(float scalar) const;
			Quaternion& operator/=(float scalar);
			Quaternion& preMultiply(const Quaternion& quat);
			Quaternion& postMultiply(const Quaternion& quat);
			Quaternion& rotate(float angle, float x, float y, float z);
			Quaternion& normalize();
			Quaternion conjugate() const;
			Quaternion inverse() const;
			float norm() const;
			Vector applyRotation(const Vector& vec) const;
			Ray applyRotation(const Ray& ray) const;
			Matrix buildMatrix() const;
			Quaternion power(double scalar);

			Quaternion slerp(const Quaternion& rhs, const float t)
			{
				Quaternion lhs = *this;
				Quaternion rv = lhs * ((lhs.inverse() * rhs).power(t));
				return rv;
			}

			static Quaternion rotation(float angle, float x, float y, float z);
			static Quaternion rotation(float angle, const Vector& axis);
		protected:
			explicit Quaternion(const Vector& vec);
	};
};

#endif
