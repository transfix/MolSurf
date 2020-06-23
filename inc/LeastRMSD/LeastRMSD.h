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
#include "UsefulMath/LinearAlgebra.h"
#include "UsefulMath/Vector.h"
#include "UsefulMath/Matrix.h"
#include <math.h>

#if !defined(LEAST_RMSD_H_INCLUDED_)
#define LEAST_RMSD_H_INCLUDED_

namespace CCVOpenGLMath {

class Vector;
class Matrix;
class LinearAlgebra;

	class LeastRMSD
	{
	public:
		LeastRMSD( );
		virtual ~LeastRMSD( );


                static void getCentroid( Vector & ctr, const std::vector< Vector > & pts );
                static float getRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2 );
                static void moveToOrigin( std::vector< Vector > & pts1, std::vector< Vector > & pts2 );
                static void getRotationMatrix( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2, Matrix & R );
                static void computeEigenValuesAndVectors( Matrix & R, std::vector< Vector > & evec, std::vector< float > & eval );
                static float computeLeastRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2, const Matrix & R, 
											  const std::vector< Vector > & evec, const std::vector< float > &eval );
                static float getLeastRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2 );
        };

};

#endif // !defined(LEAST_RMSD_H_INCLUDED_)
