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
// this class is based on the theory in
// "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds"
// Mark Meyer, Mathieu Desbrun, Peter Schroder, Alan H Barr.

#ifndef CCV_TEXMOL_GEOMETRY_MESH_DERIVATIVES_H
#define CCV_TEXMOL_GEOMETRY_MESH_DERIVATIVES_H

#include <vector>

using std::vector;

class MeshDerivatives
{
	public:
		MeshDerivatives(
			unsigned int numTris,
			unsigned int numTriVerts,
			unsigned int* tris,
			float* triVerts,
			float* m_TriVertNormals,
			float* triMeanCurv,
			float* triGaussianCurv,
			float* k1,
			float* k2);
		virtual ~MeshDerivatives();
		bool computeDerivatives();
	protected:
		bool computeFaceSet();
		bool computeDerivatives(unsigned int vertIndex);
		double computeArea(unsigned int vertIndex);
		double computeThetaSum(unsigned int vertIndex);
		bool findVerts(unsigned int triIndex, unsigned int vertIndex, unsigned int* vleft, unsigned int* vright);
		double areaOfInfinitesimalTriangle(int triangleIdx, unsigned int vertIndex, unsigned int vleft, unsigned int vright);
		double getEdgeLength(unsigned int v1, unsigned int v2);
		double voronoiArea(unsigned int vertIndex, unsigned int vleft, unsigned int vright, double angleVLeft, double angleVRight);
		double getEdgeLengthSQ(unsigned int v1, unsigned int v2);
		bool isObtuseTriangle(double angle1, double angle2, double angle3);
		double getTriangleArea(unsigned int vertIndex, unsigned int vleft, unsigned int vright);
		bool computeAngles();
		void getNeighbors(unsigned int vertIndex, vector<unsigned int> *neighboringVertices);
		bool getLaplaceBeltramiOperator(double* Kx, double* Ky, double* Kz, unsigned int vertIndex, double A_mixed);
		void unionVertex(unsigned int v, vector<unsigned int> *neighboringVertices);
		void getNeighborsAngles(unsigned int x_i, unsigned int x_j, double* alpha, double* beta, int* count);
		unsigned int m_NumTris;
		unsigned int m_NumTriVerts;
		unsigned int* m_Tris;
		float* m_TriVerts;
		float* m_TriVertNormals;
		float* m_TriMeanCurv;
		float* m_TriGaussianCurv;
		float* m_K1;
		float* m_K2;
		vector<int>* m_FaceSet;
		double* m_Angles;
};

#endif
