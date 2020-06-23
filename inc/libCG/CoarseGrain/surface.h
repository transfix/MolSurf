#ifndef CGSURFACE_H
#define CGSURFACE_H

#include <vector>

namespace GEOMETRY
{
	class Vertex;
	class Face;
	class Surface
	{
		public:
			Surface();
			~Surface();
			void addPoint(Vertex* point);
			void addFacet(Face* facet);
			void Triangle_Bezier_Patch();
			std::vector<Vertex*> m_Vertices;
			std::vector<Face*> m_Faces;
			int numbpts, numbtris, nodes_per_face;
			double min[3], max[3];
			double center[3];
		private:
			void Cubic_Bezier_Coeffs(double* p1, double* p2, double* p3, double* n1,double* n2,double* n3, double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201);
			void Simple_Bezier_B111(double* p1, double* p2, double* p3, double* n1, double* n2, double* n3, double* b300, double* b210, double* b120, double* b030, double* b201, double* b021, double* b111);
			void Gaussian_Triangle(int i, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3, int NUM);
			void Gaussian_Triangle_Patch(int i, double* p1, double* p2, double* p3, double* n1, double* n2, double* n3, double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, int NUM);
			double BB2lmd(double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, double b1, double b2, double b3);
			double Bezier_LambdaPolyZeroC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2, double* D3);
			double Bezier_ThreeLambdaPolyValueC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2,double* D3, double lmd);
			double Bezier_ThreeLambdaDervPolyValueC1(double* F1, double* N1, double* N2, double* N3, double* D1, double* D2,double* D3, double lmd);
			void BB_Gradient(double* p1_lmd, double* p2_lmd, double* p3_lmd, double* n1, double* n2, double* n3, double* b300, double* b210, double* b120, double* b030, double* b021, double* b012, double* b003, double* b102, double* b201, double* b111_1, double* b111_2, double* b111_3, double b1, double b2, double b3, double lmd, double* n);
			double DotProduct(double* v1, double* v2);
			void CrossProduct(double* v1, double* v2, double* v);
			double Determinant(double* v1, double* v2, double* v3);
			double Triangle_area(double* v1, double* v2, double* v3);
			void PolyDerivative(double* p, int n, double* dp);
			void PolyProduct(double* p1, int m, double* p2, int n, double* p3);
			double HornerEvaluate(double* p, int n, double t);
			void BoundingBox();
	};
}

#endif //CGSURFACE_H
