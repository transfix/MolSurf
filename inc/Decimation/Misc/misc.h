
#ifndef MISC_H
#define MISC_H

#include <string.h>
#include <Decimation/Bezier/BezierMath.h>
#include <Decimation/Linear/display.h>
#if 0
#include <X11/StringDefs.h>
#endif

/*
#ifdef __cplusplus
extern "C" {
#endif
*/

void CrossProduct(double[], double[], double[], double[]),
    CrossProduct1(double[], double[], double[]),
    CrossProduct2(float[], float[], float[]),
    CrossProduct2D(double[], double[], double[]),
    CrossProductF(float[], float[], float[], float[]),
    ProjectVectorToPlane(float[], float[], float[], float[]),
    SortingInteger(int[], int[], int), SortingDouble(double[], int[], int),
    /*MaxOfOneVarFunc(float (*)(), float,float,int,float,float *,float *),
     */
    RefineConvex(float (*)(float), float, float, float, float, float, float,
                 float, float *, float *),
    NormalizationOf3Vector(float[]), DNormalizationOf3Vector(double[3]),
    BernPolyFormalBartial2D(BernPoly, BernPoly, int, int),
    XuDrawPolygon(int, Vertex *, Normal *),
    DrawLevelSphere(float[], float, int),
    /*DrawTube(float [],float [], float, int), */
    /*DrawCone(float [],float [], float, float, int), */
    TriangleNormal(float[], float[], float[], float[]),
    RotateMatrix_z(double, double, double, double[3][3]),
    RotateMatrix(double[], double[], double[], double[3][3]),
    gaussinverse(double *, int, double, int *),
    gaussinverse3(double *, double, int *),
    MatrixTimesMatrix(double[], double[], double[], int, int, int),
    MinimalLeastSquare6(double[], double[], double[]),
    CompLevelSphere(float, int), CompLevelOctahedron(float, int);

extern int
/*ReadHowmanyComponents(String),  /* can not compile  with String*/
ReadHowmanyComponents(char *),
    dsvdc_(double *, int *, int *, int *, double *, double *, double *, int *,
           double *, int *, double *, int *, int *),
    LinearSystemOrder2(float, float, float, float, float, float, float *,
                       float *),
    SearchingInSequence(int[], int, int),
    CircumCenter(float[], float[], float[], float[], float[]);

float ThreePointInterpMax(float, float, float, float, float, float),
    Determinant(float[], float[], float[]), DotProduct(float[], float[]),
    LongDotProduct(float[], float[], int),
    DistanceOfTwoLines(float[], float[], float[], float[], float *, float *),
    SizeOfTriangle(float[], float[], float[]),
    AreaOfTriangle(float[], float[], float[]),
    Angle_Of_Three_Points(float[], float[], float[]),
    Distance_Two_points(float[], float[]),
    Distance_Two_points3D(float[], float[]);

double InnerProduct(double[], double[]),
    DDeterminant(double[], double[], double[]);

/*
#ifdef __cplusplus
}
#endif
*/

#endif
