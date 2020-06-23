#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace GEOMETRY
{
	class Vertex
	{
		public:
			Vertex()
			{
			}
			double   Point[3];
			double   Normal[3];
			int     NumFaces;
	};

	class Face
	{
		public:
			Face()
			{
			}
			int Index[3];      /* the indexs of a triangle-increase order*/
			Vertex node[12];    /* Gaussian integration nodes             */        // Wenqi 2/9: 4->12
			double weight[12];   /* Gaussian weights                       */       // Wenqi 2/9: 4->12
			int Orien;         /* 1--orientaioned, othwise no            */
			int AdjTri[3];     /* neighbor triangles index               */
			double Center[3];   /* the center of the triangle             */
			double Normal[3];   /* face normal                            */
			int whichnml[3];   /* Index[0] normal index                  */
			int Tou;           /* the begin index of the insided points  */
			int Wei;           /* the end index+1 of the insided points  */
			int State;         /* 1--removed; 2-- new, 0--old;           */
	};
}
#endif //GEOMETRY_H

