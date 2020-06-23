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
// Implements a 3D point class and corresponding Vector

#ifndef __Point_Three_H__
#define __Point_Three_H__

#include <Utility/utility.h>

using std::istream;
using std::ostream;

class Vector_3;

class Point_3
{
		friend class Vector_3;
		int _ptr;
	public:
		double xyz[3];
		Point_3()
		{
			xyz[0] = 0;
			xyz[1] = 0;
			xyz[2] = 0;
			_ptr = 0;
		} Point_3(int a, int b, int c)
		{
			xyz[0] = a;
			xyz[1] = b;
			xyz[2] = c;
			_ptr = 0;
		}
		Point_3(double a, double b, double c)
		{
			xyz[0] = a;
			xyz[1] = b;
			xyz[2] = c;
			_ptr = 0;
		}
		Point_3(double a, double b)
		{
			xyz[0] = a;
			xyz[1] = b;
			xyz[2] = 0;
		}
		~Point_3()
		{
		};
		friend istream& operator>>(istream&, Point_3&);
		friend ostream& operator<<(ostream&, Point_3);
		friend bool operator < (const Point_3& p1, const Point_3& p2);
		Point_3& operator =(const Point_3& p);
		int operator==(Point_3);
		int operator!=(Point_3);
		double& operator[](int);
		double operator[](int) const;
		double at(int index) const;
		void push_back(double);
		inline int size() const
		{
			return 3;
		}
		Vector_3 operator-(Point_3);
		Point_3 operator+(Vector_3);
		Point_3 operator-(Vector_3);
		Point_3& operator+=(Vector_3);
		Point_3& operator-=(Vector_3);
		// Scalar Multiplication
		friend Point_3 operator*(int, Point_3);
		friend Point_3 operator*(double, Point_3);
		friend Point_3 operator*(Point_3, int);
		friend Point_3 operator*(Point_3, double);
		// Scalar Division
		friend Point_3 operator/(Point_3, int);
		friend Point_3 operator/(Point_3, double);
		// Point_3 Addition
		friend Point_3 operator+(Point_3, Point_3);
		// Point_3 Relations
		friend double d(Point_3, Point_3);       // Distance
		friend double d2(Point_3, Point_3);      // Distance Squared
};

class Vector_3 : public Point_3
{
	public:
		Vector_3() : Point_3()
		{
		};
		Vector_3(int a, int b, int c) : Point_3(a, b, c)
		{
		};
		Vector_3(double a, double b, double c) : Point_3(a, b, c)
		{
		};

		Vector_3(double a, double b): Point_3(a,b)
		{
		};
		~Vector_3()
		{
			//this->print();
		};
		friend ostream& operator<<(ostream&, Vector_3);
		Vector_3 operator-();
		// Scalar Multiplication
		friend Vector_3 operator*(int, Vector_3);
		friend Vector_3 operator*(double, Vector_3);
		friend Vector_3 operator*(Vector_3, int);
		friend Vector_3 operator*(Vector_3, double);
		// Scalar Division
		friend Vector_3 operator/(Vector_3, int);
		friend Vector_3 operator/(Vector_3, double);
		// Vector_3 Arithmetic Operations
		Vector_3 operator+(Vector_3);
		Vector_3 operator-(Vector_3);
		double operator*(Vector_3);    // dot product
		double operator*(Point_3);     // dot product
		Vector_3 operator^(Vector_3);    // cross product
		Vector_3& operator*=(double);
		Vector_3& operator/=(double);
		Vector_3& operator+=(Vector_3);
		Vector_3& operator-=(Vector_3);
		Vector_3& operator^=(Vector_3);
		double len()                 // vector length
		{
			return sqrt(xyz[0] *xyz[0] + xyz[1] *xyz[1] + xyz[2] *xyz[2]);
		}
		double len2()                // vector length squared
		{
			return (xyz[0] *xyz[0] + xyz[1] *xyz[1] + xyz[2] *xyz[2]);
		}
		void normalize();
		void print()
		{
			cout << "<" << xyz[0] <<"," << xyz[1] << "," << xyz[2] << ">" << endl;
		}
};

#endif
