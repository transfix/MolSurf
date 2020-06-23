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
#include "LeastRMSD/LeastRMSD.h"
#include "UsefulMath/LinearAlgebra.h"
#include "UsefulMath/Vector.h"
#include "UsefulMath/Matrix.h"
#include <math.h>

using CCVOpenGLMath::Vector;
using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::LinearAlgebra;
using CCVOpenGLMath::LeastRMSD;

//////////////////////////////////////////////////////////////////////
// Constructor and Destructor
//////////////////////////////////////////////////////////////////////

LeastRMSD::LeastRMSD( )
{

}

LeastRMSD::~LeastRMSD( )
{

}


/*
  compute centroid
*/
void LeastRMSD::getCentroid( Vector & ctr, const std::vector< Vector > & pts )
{
  ctr[ 0 ] = ctr[ 1 ] = ctr[ 2 ] = 0;
  
  int n = pts.size( );

  for ( int i = 0; i < n; i++ )
    {
      Vector v = pts.at( i );
      ctr[ 0 ] += v[ 0 ];
      ctr[ 1 ] += v[ 1 ];
      ctr[ 2 ] += v[ 2 ];
    }

  ctr[ 0 ] /= n;
  ctr[ 1 ] /= n;
  ctr[ 2 ] /= n;
}


/*
  calculate rmsd between two aligned structures (trivial)
  Params: pts1 - first structure
          pts2 - second structure
  Returns: rmsd
*/
float LeastRMSD::getRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2 )
{
  float rmsd2 = 0;

  int n = pts1.size( );

  for ( int i = 0; i < n; i++ )
    {
     Vector v1 = pts1.at( i );
     Vector v2 = pts2.at( i );
          
     for ( int j = 0; j < 3; j++ ) 
       {
        float a = v1[ j ] - v2[ j ];
        rmsd2 += a * a;
       }
    }
    
  return sqrt( rmsd2 / n );
}


void LeastRMSD::moveToOrigin( std::vector< Vector > & pts1, std::vector< Vector > & pts2 )
{
  int n = pts1.size( );
  Vector c1, c2;
  
  getCentroid( c1, pts1 );
  getCentroid( c2, pts2 );
  
  for ( int i = 0; i < n; i++ )
    {
     pts1[ i ] -= c1;
     pts2[ i ] -= c2;     
    }
}



void LeastRMSD::getRotationMatrix( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2, Matrix & R )
{
  R.set( 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 1 );

  int n = pts1.size( );

  for ( int k = 0; k < n; k++ )
    {
     Vector v1 = pts1.at( k );
     Vector v2 = pts2.at( k );
          
     for ( int i=0; i < 3; i++ )
       for ( int j = 0; j < 3; j++ ) 
         {
          float f = R.get( j, i ) + v1[ i ] * v2[ j ];
          R.set( j, i, f );
         } 
    }
}


void LeastRMSD::computeEigenValuesAndVectors( Matrix & R, std::vector< Vector > & evec, std::vector< float > & eval )
{
  Matrix Rt, RtR;
  
  Rt = R.transpose( );
  RtR = R.preMultiplication( Rt );
  
  LinearAlgebra::eigenSymmetric3x3( RtR, evec, eval );
  
  evec[ 2 ] = evec[ 0 ].cross( evec[ 1 ] );
  
/*  for ( int i = 0; i < 3; i++ )
     printf( "Eigenvector %d = < %f, %f, %f >\n", i, evec[ i ][ 0 ], evec[ i ][ 1 ], evec[ i ][ 2 ] );
     
  printf( "Eigenvalues = < %f, %f, %f >\n", eval[ 0 ], eval[ 1 ], eval[ 2 ] );   */
}



float LeastRMSD::computeLeastRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2, const Matrix & R, 
                             const std::vector< Vector > & evec, const std::vector< float > &eval )
{
  float E0 = 0;
  int n = pts1.size( );

  for ( int i = 0; i < n; i++ )
     for ( int j = 0; j < 3; j++ )
        E0 += ( pts1[ i ][ j ] * pts1[ i ][ j ] + pts2[ i ][ j ] * pts2[ i ][ j ] ) / 2.0;
        
//  printf( "n = %d\n", n );      
//  printf( "E0 = %f\n", E0 );        
        
  std:vector< Vector > bvec( 3 ); 

  for ( int i = 0; i < 3; i++ )
    {
     bvec[ i ] = R * evec[ i ];
     bvec[ i ].normalize( );
    }

  bvec[ 2 ] = bvec[ 0 ].cross( bvec[ 1 ] );
  
/*  for ( int i = 0; i < 3; i++ )
     printf( "bVector %d = < %f, %f, %f >\n", i, bvec[ i ][ 0 ], bvec[ i ][ 1 ], bvec[ i ][ 2 ] );
*/  
  
  Matrix U( 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 1 );
            
  for ( int k = 0; k < 3; k++ )
    for ( int i = 0; i < 3; i++ )  
      for ( int j = 0; j < 3; j++ )
        {
         float f = U.get( j, i ) + bvec[ k ][ i ] * evec[ k ][ j ];
         U.set( j, i, f );
        }    

  float s3 = 0;
  
  if ( bvec[ 2 ].dot( R * evec[ 2 ] ) < 0 ) s3 = -1;
  else s3 = 1;
  
  float rmsd2 = 2 * ( E0 - sqrt( eval[ 0 ] ) - sqrt( eval[ 1 ] ) - s3 * sqrt( eval[ 2 ] ) ) / n ;
  
//  printf( "Least RMSD^2 = %f\n", rmsd2 );
  
  if ( rmsd2 < 0 ) rmsd2 = 0;
  
  return sqrt( rmsd2 );
}


float LeastRMSD::getLeastRMSD( const std::vector< Vector > & pts1, const std::vector< Vector > & pts2 )
{
  if ( pts1.size( ) != pts2.size( ) ) 
    {
     std::cerr << "getLeastRMSD: Both lists must have the same number of points!\n";    
     return -1;
    }

  std::vector< Vector > spts1( pts1 );
  std::vector< Vector > spts2( pts2 );  
  Matrix R;
  std::vector< Vector > evec( 3 );
  std::vector< float > eval( 3 );  
  float rmsd;

  moveToOrigin( spts1, spts2 );    
  getRotationMatrix( spts1, spts2, R );
  computeEigenValuesAndVectors( R, evec, eval );
  rmsd = computeLeastRMSD( spts1, spts2, R, evec, eval );
  
  return rmsd;
}

