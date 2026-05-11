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
#if !defined(AFX_GEOMETRY_H__A99E4F7E_DEED_48DB_9F49_0DBB21EC4211__INCLUDED_)
#define AFX_GEOMETRY_H__A99E4F7E_DEED_48DB_9F49_0DBB21EC4211__INCLUDED_

#include <Utility/utility.h>
#include <Utility/glUtility.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/AreaVolumes.h>
#include <Geometry/MeshDerivatives.h>
#include <Geometry/Texture2D.h>

class Texture2D;
class MyExtensions;

class Geometry {
public:
  Geometry();
  Geometry(const Geometry &copy);
  Geometry &operator=(const Geometry &copy);
  virtual ~Geometry();
  void ZeroMembers();
  void InitializeColors();
  void CopyPoints(const Geometry &copy);
  void CopyPointColors(const Geometry &copy);
  void CopyLines(const Geometry &copy);
  void CopyLineColors(const Geometry &copy);
  void CopyQuads(const Geometry &copy);
  void CopyQuadTexCoords(const Geometry &copy);
  void CopyTris(const Geometry &copy);
  void CopyTriTexCoords(const Geometry &copy);
  void CopyTriVertColors(const Geometry &copy);
  void CopyColors(const Geometry &copy);
  void ClearGeometry();
  void ClearPoints();
  void ClearPointColors();
  void ClearLines();
  void ClearLineColors();
  void ClearQuads();
  void ClearQuadTexCoords();
  void ClearTris();
  void ClearTriTexCoords();
  void ClearTriTexCoords2D();
  void ClearTriVertColors();
  void AllocatePoints(unsigned int NumPoints);
  void AllocatePointColors();
  void AllocateLines(unsigned int NumLineVerts, unsigned int NumLines);
  void AllocateLineColors();
  void AllocateQuads(unsigned int NumQuadVerts, unsigned int NumQuads);
  void AllocateQuadTexCoords();
  void AllocateTris(unsigned int NumTriVerts, unsigned int NumTris);
  void AllocateTriTexCoords();
  void AllocateTriTexCoords2D();
  void AllocateTriVertColors();
  void allocateTriangleDerivates();
  void CalculateQuadSmoothNormals();
  void CalculateTriSmoothNormals();
  void CalculateQuadFlatNormals();
  void CalculateTriFlatNormals();
  void SetTriNormalsReady();
  void SetQuadNormalsReady();
  void computeDerivatives();
  void computeTriangleMeshDerivatives();
  bool printTriangleDerivatives(FILE *fp);
  void setTriangleColorsByNormals();
  void setTriangleColorsByMeanCurv();
  void setTriangleColorsByGaussianCurv();
  void setTriangleColors(int renderingMode);
  bool printDerivatives(const char *filename);
  void setColors(int renderingMode);
  void SetDiffusedColor(float r, float g, float b);
  void SetSpecularColor(float r, float g, float b);
  void SetAmbientColor(float r, float g, float b);
  void SetShininess(float s);
  void SetWireframeWidth(float wireframeWidth);
  void setWireframeColor(float r, float g, float b);
  void setPointSize(float pointSize);
  void setPointColor(float r, float g, float b);
  bool useUniquePointColors();
  void SetLineWidth(float lineWidth);
  void setLineColor(float r, float g, float b);
  bool useUniqueLineColors();
  bool useUniqueTriangleColors();
  bool useUniqueWireframeColors();
  bool use2DTriangleTexture();
  bool set2DTexture(int width, int height, unsigned char *bits,
                    int numBytesPerEntry);
  bool enable2DTriangleTexture(bool enable);
  void GetReadyToDrawWire();
  void GetReadyToDrawSmooth();
  void GetReadyToDrawFlat();
  void getExtents(double *minx, double *miny, double *minz, double *maxx,
                  double *maxy, double *maxz);
  bool separateComponents(vector<Geometry *> *components);
  Geometry *merge(Geometry *geometry);
  void print();
  static bool getPointInTriangle(float x1, float y1, float z1, float x2,
                                 float y2, float z2, float x3, float y3,
                                 float z3, float *xinterp, float *yinterp,
                                 float *zinterp);

  bool m_bExtentsReady;
  bool m_bQuadFlatNormalsReady;
  bool m_bQuadSmoothNormalsReady;
  bool m_bTriFlatNormalsReady;
  bool m_bTriSmoothNormalsReady;
  bool m_RenderLines;
  bool m_RenderPoints;
  bool m_RenderQuads;
  bool m_RenderTriangles;
  bool m_RenderWireframe;
  bool m_TriangleDerivativesValid;
  bool m_Use2DTriangleTexture;
  bool m_UseLineColors;
  bool m_UsePointColors;
  bool m_UseTriangleColors;
  bool m_UseWireframeColors;
  float m_AmbientColor[4];
  float m_Center[3];
  float m_DiffuseColor[4];
  float *m_K1;
  float *m_K2;
  float *m_LineColors;
  float *m_LineVerts;
  float m_LineWidth;
  float m_Max[3];
  float m_Min[3];
  float *m_PointColorsArray;
  float *m_Points;
  float m_PointSize;
  float *m_QuadFlatNormals;
  float *m_QuadFlatTexCoords;
  float *m_QuadFlatVerts;
  float *m_QuadVertNormals;
  float *m_QuadVerts;
  float *m_QuadVertTexCoords;
  float m_Shininess;
  float m_SpecularColor[4];
  float *m_TriFlatNormals;
  float *m_TriFlatTexCoords;
  float *m_TriFlatVerts;
  float *m_TriGaussianCurv;
  float *m_TriMeanCurv;
  float *m_TriVertColorsTransparent;
  float *m_TriVertNormals;
  float *m_TriVerts;
  float *m_TriVertTexCoords;
  float *m_TriVertTexCoords2D;
  float m_UniqueLineColors[3];
  float m_UniquePointColors[3];
  float m_UniqueWireframeColors[3];
  float m_WireframeWidth;
  Texture2D *m_Texture2D;
  unsigned int *m_Lines;
  unsigned int m_NumLines;
  unsigned int m_NumLineVerts;
  unsigned int m_NumPoints;
  unsigned int m_NumQuads;
  unsigned int m_NumQuadVerts;
  unsigned int m_NumTris;
  unsigned int m_NumTriVerts;
  unsigned int *m_Quads;
  unsigned int *m_Tris;

  void CalculateExtents();
  static double getHausdorffDistance(float p11[3], float p12[3], float p13[3],
                                     float p21[3], float p22[3], float p23[3]);
  /* VBO support */
  void bindColorBuffer(bool colorsChanged);
  void bindNormalBuffer();
  void bindVertexBuffer();
  void unbindBuffer();
  /* depth sorting support */
  void sortGeometryTris();
  CCVOpenGLMath::Vector triangleIndexToCenter(int triangleIndex);
  GLuint *getDepthOrderedIndices() { return m_DepthOrderedIndices; }
  double getArea(float *point, double radius); // moved from SurfaceData
  double getArea();                            // moved from SurfaceData
  double getVolume();                          // moved from SurfaceData

protected:
  void mergeComponents(int *componentIndex, int v1, int v2, int n);
  bool getNextComponentSize(int *numVerts, int *numTris,
                            int *nextComponentIndex, int *componentIndex);
  bool addComponent(Geometry *component, int nextComponentIndex,
                    int *componentIndex);
  void resetUsedComponentIndex(int nextComponentIndex, int *componentIndex);

private:
  void CalculateTriangleNormal(float *norm, unsigned int v0, unsigned int v1,
                               unsigned int v2);
  void CalculateQuadNormal(float *norm, unsigned int v0, unsigned int v1,
                           unsigned int v2);
  /* VBO support */
  void allocateBufferObjects();
  void prepareBuffersTris();
  void deallocateBufferObjects();
  bool m_UsingBuffers;
  bool m_BuffersInitialized;
  GLuint m_VertexObject;
  GLuint m_NormalObject;
  GLuint m_ColorObject;
  MyExtensions *m_Extensions;
  /* Sorted triangles (alpha) support */
  bool m_SortInitialized;
  /* triangle index for depth-ordered rendering when transparent */
  std::vector<int> m_DepthOrderedTriangles;
  GLuint *m_DepthOrderedIndices;
  double
  distVtxToPoint(int vertex,
                 float *point); // Private because it shouldn't even BE here
  void CrossProduct(double *v1, double *v2,
                    double *v); // Private because it shouldn't even BE here
  double DotProduct(double *v1,
                    double *v2); // Private because it shouldn't even BE here
};

#endif
