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
#include <ASpline/ASpline.h>

using namespace ASPLINE;

void ASpline::Triangle_Bezier_Patch(int resolution)
{
	int count =0;
	double b300[3], b210[3], b120[3], b030[3], b021[3], b012[3], b003[3], b102[3], b201[3];
	double b111_1[9],b111_2[9],b111_3[9];
	for(int i = 0; i < numbtris; i++)
	{
		double p1[3], p2[3], p3[3], n1[3], n2[3], n3[3];
		for(int j = 0; j<3; j++)
		{
			p1[j] = verticies[faces[i]->indexInRAW[0]]->point[j];
			p2[j] = verticies[faces[i]->indexInRAW[1]]->point[j];
			p3[j] = verticies[faces[i]->indexInRAW[2]]->point[j];
			n1[j] = verticies[faces[i]->indexInRAW[0]]->normal[j];
			n2[j] = verticies[faces[i]->indexInRAW[1]]->normal[j];
			n3[j] = verticies[faces[i]->indexInRAW[2]]->normal[j];
		}

		/* Compute Bezier coefficients    */
		Cubic_Bezier_Coeffs(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201);
		/* Compute three b_{111}          */
		Simple_Bezier_B111(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b201,b021,b111_3);
		Simple_Bezier_B111(p2,p3,p1,n2,n3,n1,b030,b021,b012,b003,b120,b102,b111_1);
		Simple_Bezier_B111(p3,p1,p2,n3,n1,n2,b003,b102,b201,b300,b012,b210,b111_2);
		/* Compute the inner nodes				*/
		Triangle_Patch(i, p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3, resolution);
	}
}

void ASpline::Triangle_Bezier_Patch_With_Quadrature(const char* type, int numOfPts, int resolution)
{
	int count =0;
	Quadrature* my_quadr = new Quadrature(type, numOfPts);
	double b300[3], b210[3], b120[3], b030[3], b021[3], b012[3], b003[3], b102[3], b201[3];
	double b111_1[9],b111_2[9],b111_3[9];
	for(int i = 0; i < numbtris; i++)
	{
		double p1[3], p2[3], p3[3], n1[3], n2[3], n3[3];
		int j;
		for(j = 0; j<3; j++)
		{
			p1[j] = verticies[faces[i]->indexInRAW[0]]->point[j];
			p2[j] = verticies[faces[i]->indexInRAW[1]]->point[j];
			p3[j] = verticies[faces[i]->indexInRAW[2]]->point[j];
			n1[j] = verticies[faces[i]->indexInRAW[0]]->normal[j];
			n2[j] = verticies[faces[i]->indexInRAW[1]]->normal[j];
			n3[j] = verticies[faces[i]->indexInRAW[2]]->normal[j];
		}
		/* Compute Bezier coefficients    */
		Cubic_Bezier_Coeffs(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201);
		/* Compute three b_{111}          */
		Simple_Bezier_B111(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b201,b021,b111_3);
		Simple_Bezier_B111(p2,p3,p1,n2,n3,n1,b030,b021,b012,b003,b120,b102,b111_1);
		Simple_Bezier_B111(p3,p1,p2,n3,n1,n2,b003,b102,b201,b300,b012,b210,b111_2);
		/* Compute the inner nodes				*/
		Triangle_Patch_With_Quadrature(my_quadr,i,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3, resolution);
	}
}


void ASpline::Triangle_Bezier_Patch_Nodes(const char* type, int numOfPts)
{
	double b300[3], b210[3], b120[3], b030[3], b021[3], b012[3], b003[3], b102[3], b201[3];
	double b111_1[9],b111_2[9],b111_3[9];
	int i;
	Quadrature* my_quadr = new Quadrature(type, numOfPts);
	for(i = 0; i < numbtris; i++)
	{
		double p1[3], p2[3], p3[3], n1[3], n2[3], n3[3];
		int j;
		for(j = 0; j<3; j++)
		{
			p1[j] = verticies[faces[i]->indexInRAW[0]]->point[j];
			p2[j] = verticies[faces[i]->indexInRAW[1]]->point[j];
			p3[j] = verticies[faces[i]->indexInRAW[2]]->point[j];
			n1[j] = verticies[faces[i]->indexInRAW[0]]->normal[j];
			n2[j] = verticies[faces[i]->indexInRAW[1]]->normal[j];
			n3[j] = verticies[faces[i]->indexInRAW[2]]->normal[j];
		}
		// Compute Bezier coefficients
		Cubic_Bezier_Coeffs(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201);
		// Compute three b_{111}
		Simple_Bezier_B111(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b201,b021,b111_3);
		Simple_Bezier_B111(p2,p3,p1,n2,n3,n1,b030,b021,b012,b003,b120,b102,b111_1);
		Simple_Bezier_B111(p3,p1,p2,n3,n1,n2,b003,b102,b201,b300,b012,b210,b111_2);
		// Compute Quadrature nodes and weights
		Triangle_Patch_Nodes(my_quadr, i,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3);
		// Compute Gaussian nodes and weights over linear piece
		//        Gaussian_Triangle(i,p1,p2,p3,n1,n2,n3);
		// Compute Gaussian nodes and weights over cubic patch
		//	Gaussian_Triangle_Patch(i,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3);
	}
}


void ASpline::Triangle_Patch_Volume()
{
	double b300[3], b210[3], b120[3], b030[3], b021[3], b012[3], b003[3], b102[3], b201[3];
	double b111_1[9],b111_2[9],b111_3[9];
	int i;
	Quadrature* my_quadr = new Quadrature((char*)"Gaussian", 3);  // 3-pt Gaussian
	volume = 0.0;
	for(i = 0; i < numbtris; i++)
	{
		double p1[3], p2[3], p3[3], n1[3], n2[3], n3[3];
		int j;
		for(j = 0; j<3; j++)
		{
			p1[j] = verticies[faces[i]->indexInRAW[0]]->point[j];
			p2[j] = verticies[faces[i]->indexInRAW[1]]->point[j];
			p3[j] = verticies[faces[i]->indexInRAW[2]]->point[j];
			n1[j] = verticies[faces[i]->indexInRAW[0]]->normal[j];
			n2[j] = verticies[faces[i]->indexInRAW[1]]->normal[j];
			n3[j] = verticies[faces[i]->indexInRAW[2]]->normal[j];
		}
		// Compute Bezier coefficients
		Cubic_Bezier_Coeffs(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201);
		// Compute three b_{111}
		Simple_Bezier_B111(p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b201,b021,b111_3);
		Simple_Bezier_B111(p2,p3,p1,n2,n3,n1,b030,b021,b012,b003,b120,b102,b111_1);
		Simple_Bezier_B111(p3,p1,p2,n3,n1,n2,b003,b102,b201,b300,b012,b210,b111_2);
		// Compute Quadrature nodes and weights
		volume += Patch_Volume(my_quadr, i,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3);
		// Compute Gaussian nodes and weights over linear piece
		//        Gaussian_Triangle(i,p1,p2,p3,n1,n2,n3);
		// Compute Gaussian nodes and weights over cubic patch
		//	Gaussian_Triangle_Patch(i,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102, b201,b111_1,b111_2,b111_3);
	}
}
void ASpline::Cubic_Bezier_Coeffs(double* p1, double* p2, double* p3, double* n1,double* n2,double* n3,
								  double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201)
{
	int    j;
	double  v21[3],  v32[3],  v13[3], n21[3], n32[3], n13[3];
	double  onethird = 1.0/3.0;
	for(j = 0; j < 3; j++)
	{
		v21[j] = p2[j] - p1[j];
		v32[j] = p3[j] - p2[j];
		v13[j] = p1[j] - p3[j];
		n21[j] = n2[j] - n1[j];
		n32[j] = n3[j] - n2[j];
		n13[j] = n1[j] - n3[j];
	}
	b300[0] = 0.0;
	b300[1] = DotProduct(n1,n1);
	b300[2] = 0.0;
	b030[0] = 0.0;
	b030[1] = DotProduct(n2,n2);
	b030[2] = 0.0;
	b003[0] = 0.0;
	b003[1] = DotProduct(n3,n3);
	b003[2] = 0.0;
	b210[0] = b300[0] + onethird*DotProduct(n1, v21);
	b210[1] = b300[1] + onethird*DotProduct(n1, n21);
	b210[2] = b300[2];
	b120[0] = b030[0] - onethird*DotProduct(n2, v21);
	b120[1] = b030[1] - onethird*DotProduct(n2, n21);
	b120[2] = b030[2];
	b021[0] = b030[0] + onethird*DotProduct(n2, v32);
	b021[1] = b030[1] + onethird*DotProduct(n2, n32);
	b021[2] = b030[2];
	b012[0] = b003[0] - onethird*DotProduct(n3, v32);
	b012[1] = b003[1] - onethird*DotProduct(n3, n32);
	b012[2] = b003[2];
	b102[0] = b003[0] + onethird*DotProduct(n3, v13);
	b102[1] = b003[1] + onethird*DotProduct(n3, n13);
	b102[2] = b003[2];
	b201[0] = b300[0] - onethird*DotProduct(n1, v13);
	b201[1] = b300[1] - onethird*DotProduct(n1, n13);
	b201[2] = b300[2];
}


void ASpline::Simple_Bezier_B111(double* p1, double* p2, double* p3, double* n1, double* n2, double* n3,
								 double* b300, double* b210, double* b120, double* b030, double* b201, double* b021, double* b111)
{
	int   i;
	double v21[3],v13[3],v23[3],
		   n21[3],n13[3],n23[3],an12[3], m1[3],m2[3];
	double Fb1[3], Fb2[3], Flmd[2], w[3],w1[3],w2[3],w3[4],left[5];
	double det[3],                        grad[6];
	/*
	for (i = 0; i < 3; i++) {
	   grad[i] = 0.5*(g1[i] + g2[i]);
	}
	*/
	for(i = 0; i < 3; i++)
	{
		v21[i] = p2[i] - p1[i];
		v13[i] = p1[i] - p3[i];
		v23[i] = p2[i] - p3[i];
		n21[i] = n2[i] - n1[i];
		n13[i] = n1[i] - n3[i];
		n23[i] = n2[i] - n3[i];
		an12[i] = 0.5*(n1[i] + n2[i]);
	}
	/* compute d3                */
	CrossProduct(v21, an12, m1);
	CrossProduct(n21, an12, m2);
	/* compute ||d3||^2          */
	b111[5] = DotProduct(m1,m1);
	b111[6] = 2*DotProduct(m1,m2);
	b111[7] = DotProduct(m2,m2);
	/* compute det(T)            */
	det[0] = Determinant(v13,v23,an12);
	det[1] = Determinant(n13,v23,an12)+
			 Determinant(v13,n23,an12);
	det[2] = Determinant(n13,n23,an12);
	/* compute partial derivative*/
	for(i = 0; i < 3; i++)
	{
		Fb1[i] = 0.75*(b300[i] - b201[i] + 2*b210[i] + b120[i] - b021[i]);
		Fb2[i] = 0.75*(b210[i] - b201[i] + 2*b120[i] + b030[i] - b021[i]);
	}
	Flmd[0] = (b300[1] + 3*(b210[1] + b120[1]) + b030[1])/8.0;  // constant term
	Flmd[1] = (b300[2] + 3*(b210[2] + b120[2]) + b030[2])/4.0;  // 1st degree term
	/* Wenqi changed on 2/8
	           // compute directional deriv
	    w[0] = DotProduct(m1,an12);
	    w[1] = DotProduct(m2,an12);
	    w[2] = 0.0;

	           // times det(T) in both sides
	    left[0] = w[0]*det[0];
	    left[1] = w[0]*det[1] + w[1]*det[0];
	    left[2] = w[0]*det[2] + w[1]*det[1] + w[2]*det[0];
	    left[3] = w[1]*det[2] + w[2]*det[1];
	    left[4] = w[2]*det[2];
	*/
	/* Wenqi added on 2/8  Begin */
	left[0] = 0.0;
	left[1] = 0.0;
	left[2] = 0.0;
	left[3] = 0.0;
	left[4] = 0.0;
	/*                     End  */
	/* The above change did on 2/8 does not affect the result because d3 dot an12 = 0 */
	/*compute d3^TT^{-1} Gradient*/
	CrossProduct(v23, an12, w);            // w,w1=BxC
	CrossProduct(n23, an12, w1);
	/* Coefficient of Fb1        */
	w2[0] = DotProduct(m1,w);
	w2[1] = DotProduct(m1,w1) + DotProduct(m2,w);
	w2[2] = DotProduct(m2,w1);
	/* times Fb1                 */
	left[0] = left[0] -(w2[0]*Fb1[0]);
	left[1] = left[1] -(w2[0]*Fb1[1] + w2[1]*Fb1[0]);
	left[2] = left[2] -(w2[0]*Fb1[2] + w2[1]*Fb1[1] + w2[2]*Fb1[0]);
	left[3] = left[3] -(w2[1]*Fb1[2] + w2[2]*Fb1[1]);
	left[4] = left[4] -(w2[2]*Fb1[2]);        // it turns out that Fb1[2]=0
	/* Coefficient of Fb2        */
	CrossProduct(an12,v13, w);
	CrossProduct(an12,n13, w1);
	w2[0] = DotProduct(m1,w);
	w2[1] = DotProduct(m1,w1) + DotProduct(m2,w);
	w2[2] = DotProduct(m2,w1);
	/* times Fb2                 */
	left[0] = left[0] -(w2[0]*Fb2[0]);
	left[1] = left[1] -(w2[0]*Fb2[1] + w2[1]*Fb2[0]);
	left[2] = left[2] -(w2[0]*Fb2[2] + w2[1]*Fb2[1] + w2[2]*Fb2[0]);
	left[3] = left[3] -(w2[1]*Fb2[2] + w2[2]*Fb2[1]);
	left[4] = left[4] -(w2[2]*Fb2[2]);
	/* Coefficient of Flmd       */
	CrossProduct(v13,v23, w);
	CrossProduct(v13,n23, w1);
	CrossProduct(n13,v23, w2);
	for(i =0; i < 3; i++)
	{
		w1[i] = w1[i] + w2[i];
	}
	CrossProduct(n13,n23, w2);
	w3[0] = DotProduct(m1,w);
	w3[1] = DotProduct(m1,w1) + DotProduct(m2,w);
	w3[2] = DotProduct(m1,w2) + DotProduct(m2,w1);
	w3[3] = DotProduct(m2,w2);
	/* times Flmd                */
	left[0] = left[0] -(w3[0]*Flmd[0]);
	left[1] = left[1] -(w3[0]*Flmd[1] + w3[1]*Flmd[0]);
	left[2] = left[2] -(w3[1]*Flmd[1] + w3[2]*Flmd[0]);
	left[3] = left[3] -(w3[2]*Flmd[1] + w3[3]*Flmd[0]);
	left[4] = left[4] -(w3[3]*Flmd[1]);        // it turns out that Flmd[1]=0
	b111[0] = -2.0/3.0*left[0];
	b111[1] = -2.0/3.0*left[1];
	b111[2] = -2.0/3.0*left[2];
	b111[3] = -2.0/3.0*left[3];
	b111[4] = -2.0/3.0*left[4];   // so b111[4]=0
}


// arand: 4/22/2010. new function for computing refined a-spline triangulation and quadrature together
void ASpline::Triangle_Patch_With_Quadrature(Quadrature* my_quadr,int t, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3, double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, int resolution)
{
	// do the subdivision...
	Triangle_Patch(t,p1,p2,p3,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,resolution);
	double b1, b2, b3;
	int i,j,k,l;
	double scale = 1.0/resolution;
	double lmd;
	double p1_lmd[3], p2_lmd[3], p3_lmd[3];
	double n[3];
	double x_b1, y_b1, z_b1, x_b2, y_b2, z_b2, E, F, G;
	// get quadrature for each subfacet
	//
	//           ----
	//          |\   |
	// lower    | \  | upper triangle
	// triangle |  \ |
	//          |   \|
	//           ----
	//       (i,j)
	//
	for(j = 0; j < resolution; j++)
	{
		b2 = (double)j/resolution;
		for(i = 0; i < resolution-1-j; i++)
		{
			b1 = (double)i/resolution;
			// b1,b2,b3 is the barycentric coords for the lower left corner...
			for(k = 0; k < my_quadr->numOfPts_perTri; k++)
			{
				// add quadrature points for two triangles
				double b1lower = b1 + my_quadr->xi1[k]*scale;
				double b2lower = b2 + my_quadr->xi2[k]*scale;
				double b3lower = 1.0 - b1lower - b2lower;
				lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1lower,b2lower,b3lower);
				for(l = 0; l < 3; l++)
				{
					p1_lmd[l] = p1[l] + lmd*n1[l];
					p2_lmd[l] = p2[l] + lmd*n2[l];
					p3_lmd[l] = p3[l] + lmd*n3[l];
				}
				BB_Gradient(p1_lmd,p2_lmd,p3_lmd,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1lower,b2lower,b3lower,lmd,n);
				Vertex* node = new Vertex();
				for(l = 0; l < 3; l++)
				{
					node->point[l] = b1lower*p1_lmd[l] + b2lower*p2_lmd[l] + b3lower*p3_lmd[l];
					node->normal[l] = n[l];
				}
				x_b1 = p1_lmd[0]-p3_lmd[0];
				y_b1 = p1_lmd[1]-p3_lmd[1];
				z_b1 = p1_lmd[2]-p3_lmd[2];
				x_b2 = p2_lmd[0]-p3_lmd[0];
				y_b2 = p2_lmd[1]-p3_lmd[1];
				z_b2 = p2_lmd[2]-p3_lmd[2];
				E = x_b1*x_b1 + y_b1*y_b1 + z_b1*z_b1;
				F = x_b1*x_b2 + y_b1*y_b2 + z_b1*z_b2;
				G = x_b2*x_b2 + y_b2*y_b2 + z_b2*z_b2;
				node->weight = my_quadr->w[k] * sqrt(E*G-F*F)*scale*scale/2.0; // arand: added scale
				faces[t]->node.push_back(node);
				double b1upper = b1 + scale - my_quadr->xi1[k]*scale;
				double b2upper = b2 + scale - my_quadr->xi2[k]*scale;
				double b3upper = 1.0 - b1upper - b2upper;
				lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1upper,b2upper,b3upper);
				for(l = 0; l < 3; l++)
				{
					p1_lmd[l] = p1[l] + lmd*n1[l];
					p2_lmd[l] = p2[l] + lmd*n2[l];
					p3_lmd[l] = p3[l] + lmd*n3[l];
				}
				BB_Gradient(p1_lmd,p2_lmd,p3_lmd,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1upper,b2upper,b3upper,lmd,n);
				node = new Vertex();
				for(l = 0; l < 3; l++)
				{
					node->point[l] = b1upper*p1_lmd[l] + b2upper*p2_lmd[l] + b3upper*p3_lmd[l];
					node->normal[l] = n[l];
				}
				x_b1 = p1_lmd[0]-p3_lmd[0];
				y_b1 = p1_lmd[1]-p3_lmd[1];
				z_b1 = p1_lmd[2]-p3_lmd[2];
				x_b2 = p2_lmd[0]-p3_lmd[0];
				y_b2 = p2_lmd[1]-p3_lmd[1];
				z_b2 = p2_lmd[2]-p3_lmd[2];
				E = x_b1*x_b1 + y_b1*y_b1 + z_b1*z_b1;
				F = x_b1*x_b2 + y_b1*y_b2 + z_b1*z_b2;
				G = x_b2*x_b2 + y_b2*y_b2 + z_b2*z_b2;
				node->weight = my_quadr->w[k] * sqrt(E*G-F*F)*scale*scale/2.0; // arand: added scale
				faces[t]->node.push_back(node);
			}
		}
		// add quadrature points for last triangle
		b1 = ((double)resolution - 1.0 - (double)j)/resolution;
		for(k = 0; k < my_quadr->numOfPts_perTri; k++)
		{
			// add quadrature points for two triangles
			double b1lower = b1 + my_quadr->xi1[k]*scale;
			double b2lower = b2 + my_quadr->xi2[k]*scale;
			double b3lower = 1.0 - b1lower - b2lower;
			lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1lower,b2lower,b3lower);
			for(l = 0; l < 3; l++)
			{
				p1_lmd[l] = p1[l] + lmd*n1[l];
				p2_lmd[l] = p2[l] + lmd*n2[l];
				p3_lmd[l] = p3[l] + lmd*n3[l];
			}
			BB_Gradient(p1_lmd,p2_lmd,p3_lmd,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1lower,b2lower,b3lower,lmd,n);
			Vertex* node = new Vertex();
			for(l = 0; l < 3; l++)
			{
				node->point[l] = b1lower*p1_lmd[l] + b2lower*p2_lmd[l] + b3lower*p3_lmd[l];
				node->normal[l] = n[l];
			}
			x_b1 = p1_lmd[0]-p3_lmd[0];
			y_b1 = p1_lmd[1]-p3_lmd[1];
			z_b1 = p1_lmd[2]-p3_lmd[2];
			x_b2 = p2_lmd[0]-p3_lmd[0];
			y_b2 = p2_lmd[1]-p3_lmd[1];
			z_b2 = p2_lmd[2]-p3_lmd[2];
			E = x_b1*x_b1 + y_b1*y_b1 + z_b1*z_b1;
			F = x_b1*x_b2 + y_b1*y_b2 + z_b1*z_b2;
			G = x_b2*x_b2 + y_b2*y_b2 + z_b2*z_b2;
			node->weight = my_quadr->w[k] * sqrt(E*G-F*F)*scale*scale/2.0; // arand: added scale
			faces[t]->node.push_back(node);
		}
	}
}

// compute refined triangulation on the A-Spline
void ASpline::Triangle_Patch(int t, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3, double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, int resolution)
{
	double lmd;
	double b1, b2, b3;
	double p1_lmd[3], p2_lmd[3], p3_lmd[3];
	double n[3];
	int I[3];
	int i, j, k;
	int* p_index = new int[(resolution+1)*(resolution+2)/2]; // new vertices
	int numOfVts_now = verticies.size(); // original vertices
	int count1 = 0;
	int count2 = 0;
	for(j = 0; j < resolution+1; j++)
	{
		b2 = (double)j/resolution;
		for(i = 0; i < resolution+1-j; i++)
		{
			b1 = (double)i/resolution;
			b3 = 1.0 - b1 - b2;
			//cout << "Bary " << b1 << " " << b2 << " " << b3 << endl;
			if(b1 != 1.0 && b2 != 1.0 && b3 != 1.0)
			{
				// compute lambda value in A-spline representation given bezier coefficients
				// and the barycentric coordinates
				lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1,b2,b3);
				for(k = 0; k < 3; k++)
				{
					p1_lmd[k] = p1[k] + lmd*n1[k];
					p2_lmd[k] = p2[k] + lmd*n2[k];
					p3_lmd[k] = p3[k] + lmd*n3[k];
				}
				// get the normal to the surface at the new point
				BB_Gradient(p1_lmd,p2_lmd,p3_lmd,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1,b2,b3,lmd,n);
				Vertex* p = new Vertex();
				// compute the new vertex location and assign its normal
				for(k = 0; k < 3; k++)
				{
					p->point[k] = b1*p1_lmd[k] + b2*p2_lmd[k] + b3*p3_lmd[k];
					p->normal[k] = n[k];
					p->color[k] = verticies[faces[t]->indexInRAW[0]]->color[k];
				}
				//cout << endl;
				addPoint(p);
				p_index[count1] = count2 + numOfVts_now;
				count2++;
			}
			if(b1 == 1.0)
			{
				p_index[count1] = faces[t]->indexInRAW[0];
			}
			if(b2 == 1.0)
			{
				p_index[count1] = faces[t]->indexInRAW[1];
			}
			if(b3 == 1.0)
			{
				p_index[count1] = faces[t]->indexInRAW[2];
			}
			count1++;
		}
	}
	for(j = 0; j < resolution; j++)
	{
		for(i = 0; i < resolution-1-j; i++)
		{
			Face* f_low = new Face();
			I[0] = p_index[computeTriIndex(i, j, resolution)];
			I[1] = p_index[computeTriIndex(i+1, j, resolution)];
			I[2] = p_index[computeTriIndex(i, j+1, resolution)];
			for(k = 0; k < 3; k++)
			{
				f_low->indexInRAW[k] = I[k];
			}
			faces[t]->subFacets.push_back(f_low);
			Face* f_up = new Face();
			I[0] = p_index[computeTriIndex(i+1, j, resolution)];
			I[1] = p_index[computeTriIndex(i+1, j+1, resolution)];
			I[2] = p_index[computeTriIndex(i, j+1, resolution)];
			for(k = 0; k < 3; k++)
			{
				f_up->indexInRAW[k] = I[k];
			}
			faces[t]->subFacets.push_back(f_up);
		}
		Face* f_end = new Face();
		I[0] = p_index[computeTriIndex(resolution-1-j, j, resolution)];
		I[1] = p_index[computeTriIndex(resolution-j, j, resolution)];
		I[2] = p_index[computeTriIndex(resolution-1-j, j+1, resolution)];
		for(k = 0; k < 3; k++)
		{
			f_end->indexInRAW[k] = I[k];
		}
		faces[t]->subFacets.push_back(f_end);
	}
}


void ASpline::Triangle_Patch_Nodes(Quadrature* my_quadr, int i, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3,
								   double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003,
								   double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3)
{
	double lmd;
	double p1_lmd[3], p2_lmd[3], p3_lmd[3];
	double n[3];
	double x_b1, y_b1, z_b1, x_b2, y_b2, z_b2, E, F, G;
	double b1, b2, b3;
	int j;
	int k;
	for(k = 0; k < my_quadr->numOfPts_perTri; k++)
	{
		b1 = my_quadr->xi1[k];
		b2 = my_quadr->xi2[k];
		b3 = my_quadr->xi3[k];
		lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1,b2,b3);
		for(j = 0; j < 3; j++)
		{
			p1_lmd[j] = p1[j] + lmd*n1[j];
			p2_lmd[j] = p2[j] + lmd*n2[j];
			p3_lmd[j] = p3[j] + lmd*n3[j];
		}
		BB_Gradient(p1_lmd,p2_lmd,p3_lmd,n1,n2,n3,b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1,b2,b3,lmd,n);
		Vertex* node = new Vertex();
		for(j = 0; j < 3; j++)
		{
			node->point[j] = b1*p1_lmd[j] + b2*p2_lmd[j] + b3*p3_lmd[j];
			node->normal[j] = n[j];
		}
		x_b1 = p1_lmd[0]-p3_lmd[0];
		y_b1 = p1_lmd[1]-p3_lmd[1];
		z_b1 = p1_lmd[2]-p3_lmd[2];
		x_b2 = p2_lmd[0]-p3_lmd[0];
		y_b2 = p2_lmd[1]-p3_lmd[1];
		z_b2 = p2_lmd[2]-p3_lmd[2];
		E = x_b1*x_b1 + y_b1*y_b1 + z_b1*z_b1; // see Wenqi Zhao's thesis to make sense of this...
		F = x_b1*x_b2 + y_b1*y_b2 + z_b1*z_b2;
		G = x_b2*x_b2 + y_b2*y_b2 + z_b2*z_b2;
		node->weight = my_quadr->w[k] * sqrt(E*G-F*F)/2.0; // added /2.0 on 12-6-2011
		faces[i]->node.push_back(node);
	}
}

double ASpline::Patch_Volume(Quadrature* my_quadr, int i, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3,
							 double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003,
							 double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3)
{
	double lmd;
	double p1_lmd[3], p2_lmd[3], p3_lmd[3];
	double n[3];
	double x_b1, y_b1, z_b1, x_b2, y_b2, z_b2, E, F, G;
	double b1, b2, b3;
	int j;
	double p_volume = 0.0;
	double tri_area = Triangle_area(p1, p2, p3);
	int k;
	for(k = 0; k < my_quadr->numOfPts_perTri; k++)
	{
		b1 = my_quadr->xi1[k];
		b2 = my_quadr->xi2[k];
		b3 = my_quadr->xi3[k];
		lmd = BB2lmd(b300,b210,b120,b030,b021,b012,b003,b102,b201,b111_1,b111_2,b111_3,b1,b2,b3);
		p_volume += tri_area * my_quadr->w[k] * lmd;
	}
	return p_volume;
}

double ASpline::BB2lmd(double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003,
					   double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, double b1, double b2, double b3)
{
	double lmd;
	int i;
	double bb300, bb210, bb120, bb030, bb201, bb111, bb021, bb102, bb012, bb003;
	double bb200, bb110, bb020, bb101, bb011, bb002;
	double w, w1, w2, w3;
	double pol[3], p31[4], p32[4], p33[4], q21[3], q22[3], q23[3];
	for(i = 0; i < 3; i++)
	{
		q21[i] = b111_1[5+i];
		q22[i] = b111_2[5+i];
		q23[i] = b111_3[5+i];
	}
	/* quadratic Bezier base          */
	bb200 = b1*b1;
	bb020 = b2*b2;
	bb002 = b3*b3;
	bb110 = b1*b2;
	w3 = bb110;
	bb110 = bb110 + bb110;
	bb101 = b1*b3;
	w2 = bb101;
	bb101 = bb101 + bb101;
	bb011 = b2*b3;
	w1 = bb011;
	bb011 = bb011 + bb011;
	/* cubic Bezier base              */
	bb300 = bb200*b1;
	bb210 = 3*bb200*b2;
	bb120 = 3*b1*bb020;
	bb030 = bb020*b2;
	bb201 = 3*bb200*b3;
	bb102 = 3*b1*bb002;
	bb021 = 3*bb020*b3;
	bb012 = 3*b2*bb002;
	bb003 = bb002*b3;
	bb111 = 3*bb110*b3;
	for(i = 0; i < 3; i++)
	{
		pol[i] = bb300*b300[i] + bb210*b210[i] + bb120*b120[i] +
				 bb030*b030[i] + bb201*b201[i] + bb102*b102[i] +
				 bb021*b021[i] + bb012*b012[i] + bb003*b003[i];
	}
	w1 = w1*w1;
	w2 = w2*w2;
	w3 = w3*w3;
	w = w1 + w2 + w3;
	for(i = 0; i < 4; i++)
	{
		p31[i] = bb111*w1*b111_1[i]/w;
		p32[i] = bb111*w2*b111_2[i]/w;
		p33[i] = bb111*w3*b111_3[i]/w;
	}
	lmd = Bezier_LambdaPolyZeroC1(pol,p31,p32,p33,q21,q22,q23);
	return lmd;
}


double ASpline::Bezier_LambdaPolyZeroC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2, double* D3)
{
	double DF1[2], DN1[7], DN2[7], DN3[7], W1[5], W2[5], W3[5], WW1[7], WW2[7], WW3[7];
	double lmd, lmd1, w1, f, df;
	int i,k;
	PolyDerivative(F1, 1, DF1);
	PolyDerivative(N1, 3, W1);
	PolyDerivative(N2, 3, W2);
	PolyDerivative(N3, 3, W3);
	PolyProduct(W1, 2, D1, 2, DN1);
	PolyProduct(W2, 2, D2, 2, DN2);
	PolyProduct(W3, 2, D3, 2, DN3);
	PolyDerivative(D1, 2, W1);
	PolyDerivative(D2, 2, W2);
	PolyDerivative(D3, 2, W3);
	PolyProduct(W1, 1, N1, 3, WW1);
	PolyProduct(W2, 1, N2, 3, WW2);
	PolyProduct(W3, 1, N3, 3, WW3);
	for(i = 0; i <= 4; i++)
	{
		DN1[i] = DN1[i] - WW1[i];
		DN2[i] = DN2[i] - WW2[i];
		DN3[i] = DN3[i] - WW3[i];
	}
	/* compute initial value of root */
	lmd1 = -F1[0]/F1[1];
	w1 = lmd1;
	/* Newton iteration, here should combine besection */
	lmd = 10000.0;
	k = 0;
	while(fabs(lmd1 - lmd) > 1.0e-5)    // 1.0e-15
	{
		lmd = lmd1;
		f = Bezier_ThreeLambdaPolyValueC1(F1,N1,N2,N3,D1,D2,D3,lmd);
		df = Bezier_ThreeLambdaDervPolyValueC1(DF1,DN1,DN2,DN3,D1,D2,D3,lmd);
		lmd1 = lmd - f/df;
		k = k + 1;
		if(df == 0.0)
		{
			printf("lmd,lmd1 = %lf,  %lf\n", lmd,lmd1);
		}
		if(k > 10)
		{
			printf("k > 10 in Newton iteration\n");
			return w1;
		}
	}
	return lmd1;
}


double ASpline::Bezier_ThreeLambdaPolyValueC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2,double* D3, double lmd)
{
	double w, f1,f2,f3,F;
	w = HornerEvaluate(F1,1,lmd);
	f1 = HornerEvaluate(N1,3,lmd)/HornerEvaluate(D1,2,lmd);
	f2 = HornerEvaluate(N2,3,lmd)/HornerEvaluate(D2,2,lmd);
	f3 = HornerEvaluate(N3,3,lmd)/HornerEvaluate(D3,2,lmd);
	F = w + f1 + f2 + f3;;
	return F;
}


void ASpline::BB_Gradient(double* p1_lmd, double* p2_lmd, double* p3_lmd, double* n1, double* n2, double* n3,
						  double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003,
						  double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3,
						  double b1, double b2, double b3, double lmd, double* n)
{
	double bb200, bb020, bb002, bb110, bb101, bb011, w1, w2, w3, w;
	double bb300, bb030, bb003, bb210, bb201, bb120, bb021, bb102, bb012, bb111;
	double pol[3], p31[4], p32[4], p33[4], q21[3], q22[3], q23[3];
	double dpol[2], DN1[7], DN2[7], DN3[7], W1[5], W2[5], W3[5], WW1[7], WW2[7], WW3[7];
	double A[3], B[3], C[3], det_T, T1[3], T2[3], T3[3];
	double df_b1, df_b2, df_lmd;
	int i;
	for(i = 0; i < 3; i++)
	{
		q21[i] = b111_1[5+i];
		q22[i] = b111_2[5+i];
		q23[i] = b111_3[5+i];
	}
	/* quadratic Bezier base          */
	bb200 = b1*b1;
	bb020 = b2*b2;
	bb002 = b3*b3;
	bb110 = b1*b2;
	w3 = bb110;
	bb110 = bb110 + bb110;
	bb101 = b1*b3;
	w2 = bb101;
	bb101 = bb101 + bb101;
	bb011 = b2*b3;
	w1 = bb011;
	bb011 = bb011 + bb011;
	w1 = w1*w1;
	w2 = w2*w2;
	w3 = w3*w3;
	w = w1 + w2 + w3;
	/* DF_b1  */
	for(i = 0; i < 3; i++)
	{
		pol[i] = 3*(bb200*b300[i] + bb110*b210[i] + bb020*b120[i] -
					bb020*b021[i] - bb011*b012[i] - bb002*b003[i] +
					(bb002-bb101)*b102[i] + (bb101-bb200)*b201[i]
				   );
	}
	for(i = 0; i < 4; i++)
	{
		p31[i] = 3*(bb011-bb110)*w1*b111_1[i]/w;
		p32[i] = 3*(bb011-bb110)*w2*b111_2[i]/w;
		p33[i] = 3*(bb011-bb110)*w3*b111_3[i]/w;
	}
	df_b1 = Bezier_ThreeLambdaPolyValueC1(pol,p31,p32,p33,q21,q22,q23,lmd);
	/* DF_b2  */
	for(i = 0; i < 3; i++)
	{
		pol[i] = 3*(bb200*b210[i] + bb110*b120[i] + bb020*b030[i] +
					(bb011-bb020)*b021[i] + (bb002-bb011)*b012[i] -
					bb002*b003[i] - bb101*b102[i] - bb200*b201[i]
				   );
	}
	for(i = 0; i < 4; i++)
	{
		p31[i] = 3*(bb101-bb110)*w1*b111_1[i]/w;
		p32[i] = 3*(bb101-bb110)*w2*b111_2[i]/w;
		p33[i] = 3*(bb101-bb110)*w3*b111_3[i]/w;
	}
	df_b2 = Bezier_ThreeLambdaPolyValueC1(pol,p31,p32,p33,q21,q22,q23,lmd);
	/* DF_lmd   */
	/* cubic Bezier base              */
	bb300 = bb200*b1;
	bb210 = 3*bb200*b2;
	bb120 = 3*b1*bb020;
	bb030 = bb020*b2;
	bb201 = 3*bb200*b3;
	bb102 = 3*b1*bb002;
	bb021 = 3*bb020*b3;
	bb012 = 3*b2*bb002;
	bb003 = bb002*b3;
	bb111 = 3*bb110*b3;
	for(i = 0; i < 3; i++)
	{
		pol[i] = bb300*b300[i] + bb210*b210[i] + bb120*b120[i] +
				 bb030*b030[i] + bb201*b201[i] + bb102*b102[i] +
				 bb021*b021[i] + bb012*b012[i] + bb003*b003[i];
	}
	for(i = 0; i < 4; i++)
	{
		p31[i] = bb111*w1*b111_1[i]/w;
		p32[i] = bb111*w2*b111_2[i]/w;
		p33[i] = bb111*w3*b111_3[i]/w;
	}
	PolyDerivative(pol, 1, dpol);
	PolyDerivative(p31, 3, W1);
	PolyDerivative(p32, 3, W2);
	PolyDerivative(p33, 3, W3);
	PolyProduct(W1, 2, q21, 2, DN1);
	PolyProduct(W2, 2, q22, 2, DN2);
	PolyProduct(W3, 2, q23, 2, DN3);
	PolyDerivative(q21, 2, W1);
	PolyDerivative(q22, 2, W2);
	PolyDerivative(q23, 2, W3);
	PolyProduct(W1, 1, p31, 3, WW1);
	PolyProduct(W2, 1, p32, 3, WW2);
	PolyProduct(W3, 1, p33, 3, WW3);
	for(i = 0; i <= 4; i++)
	{
		DN1[i] = DN1[i] - WW1[i];
		DN2[i] = DN2[i] - WW2[i];
		DN3[i] = DN3[i] - WW3[i];
	}
	df_lmd = Bezier_ThreeLambdaDervPolyValueC1(dpol,DN1,DN2,DN3,q21,q22,q23,lmd);
	/* Compute T^{-1}  */
	for(i = 0; i< 3; i++)
	{
		A[i] = p1_lmd[i] - p3_lmd[i];
		B[i] = p2_lmd[i] - p3_lmd[i];
		C[i] = b1*n1[i] + b2*n2[i] + b3*n3[i];
	}
	det_T = Determinant(A,B,C);
	CrossProduct(B,C,T1);
	CrossProduct(C,A,T2);
	CrossProduct(A,B,T3);
	/* Compute Gradient    */
	for(i = 0; i< 3; i++)
	{
		n[i] = 1.0/det_T * (df_b1*T1[i] + df_b2*T2[i] + df_lmd*T3[i]);
	}
}


double ASpline::Bezier_ThreeLambdaDervPolyValueC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2,double* D3, double lmd)
{
	double w1,w2,w3,F;
	w1 = HornerEvaluate(D1,2,lmd);
	w2 = HornerEvaluate(D2,2,lmd);
	w3 = HornerEvaluate(D3,2,lmd);
	F = F1[0] + HornerEvaluate(N1,4,lmd)/(w1*w1) +
		HornerEvaluate(N2,4,lmd)/(w2*w2) +
		HornerEvaluate(N3,4,lmd)/(w3*w3);
	return F;
}


double ASpline::DotProduct(double* v1, double* v2)
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void ASpline::CrossProduct(double* v1, double* v2, double* v)
{
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double ASpline::Determinant(double* v1, double* v2, double* v3)
{
	return v1[0]*v2[1]*v3[2] + v2[0]*v3[1]*v1[2] + v3[0]*v1[1]*v2[2] -
		   v1[2]*v2[1]*v3[0] - v2[2]*v3[1]*v1[0] - v3[2]*v1[1]*v2[0];
}

double ASpline::Triangle_area(double* v1, double* v2, double* v3)
{
	/* Heron's formula  */
	double a,b,c,s;
	a = sqrt((v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) + (v1[2]-v2[2])*(v1[2]-v2[2]));
	b = sqrt((v2[0]-v3[0])*(v2[0]-v3[0]) + (v2[1]-v3[1])*(v2[1]-v3[1]) + (v2[2]-v3[2])*(v2[2]-v3[2]));
	c = sqrt((v1[0]-v3[0])*(v1[0]-v3[0]) + (v1[1]-v3[1])*(v1[1]-v3[1]) + (v1[2]-v3[2])*(v1[2]-v3[2]));
	s = (a+b+c)/2.0;
	return sqrt(s*(s-a)*(s-b)*(s-c));
}

void  ASpline::PolyDerivative(double* p, int n, double* dp)
{
	int   i;
	for(i = 1; i <= n; i++)
	{
		dp[i-1] = i*p[i];
	}
}

void  ASpline::PolyProduct(double* p1, int m, double* p2, int n, double* p3)
{
	int   k, i;
	double sum;
	for(k = 0; k <= m+n; k++)
	{
		sum = 0.0;
		for(i = 0; i <= k; i++)
		{
			if(i <= m && k - i <= n)
			{
				sum = sum + p1[i]*p2[k-i];
			}
		}
		p3[k] = sum;
	}
}

double ASpline::HornerEvaluate(double* p, int n, double t)
{
	double result;
	int i;
	result = p[n];
	for(i = 0; i < n; i++)
	{
		result = result*t + p[n-i-1];
	}
	return result;
}

int ASpline::computeTriIndex(int i, int j, int resolution)
{
	int index = (2*resolution + 3 - j) * j / 2 + i;
	return index;
}
