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
#if !defined(AFX_LINEARALGEBRA_H__AD644F94_D742_4858_A926_AC7E6964B1E2__INCLUDED_)
#define AFX_LINEARALGEBRA_H__AD644F94_D742_4858_A926_AC7E6964B1E2__INCLUDED_

#include <Utility/utility.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Matrix.h>

namespace CCVOpenGLMath
{

	class LinearAlgebra
	{
		public:
			LinearAlgebra();
			virtual ~LinearAlgebra();
			static bool getCylinderFit(int n, double* x, double* y, double* z, CCVOpenGLMath::Vector* p1, CCVOpenGLMath::Vector* p2, double* radius);
			// fit a line y = mx + c minimizing the least squares norm.
			static bool leastSquares(int n, double* x, double* y, double* m, double* c, double* radius);
			static bool mean(double* x, int n, double* mean);
			static bool summation(double* x, int n, double* sum);
			static bool sumOfSquares(double* x, int n, double* sumSquare);
			static bool dotProduct(double* x, double* y, int n, double* dotProd);
			static bool correlate(double* x, double* y, int n, double* correlationCoefficient);
			static bool selectivelyCorrelate(double* x, int rangeToCorrelate, double* y, int n, double* correlationCoefficient);
			static bool discretize(double* x, int n, double posVal, double negVal);
			// solve ax=b
			static bool solveSystem(double a11, double a12, double a13,
					double a21, double a22, double a23,
					double a31, double a32, double a33,
					double b1,  double b2,  double b3,
					double* x,  double* y,  double* z);
			static bool solve2x2System(double a11, double a12, double b1,
					double a21, double a22, double b2,
					double* x, double* y);
			static bool solveEigenSystem(double c11, double c12, double c13,
					double c21, double c22, double c23,
					double c31, double c32, double c33,
					double* k1Vec, double* k2Vec, double g1, double g2);
			static bool solveDependentEquations(double c11, double c12, double c13,
					double c21, double c22, double c23,
					double c31, double c32, double c33,
					double* vec);

			static void tridiagonalHouseholder( double mat[ 3 ][ 3 ], double val[ 3 ], double val2[ 3 ] ); 
			static void tridiagonalQL( double mat[ 3 ][ 3 ], double val[ 3 ], double val2[ 3 ] );
			static void eigenSymmetric3x3( const Matrix &mat, std::vector< Vector > &vec, std::vector< float > &val );


	};
};

#endif
