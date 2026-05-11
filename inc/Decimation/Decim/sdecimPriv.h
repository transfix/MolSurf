#ifndef SDECIMPRIV_H
#define SDECIMPRIV_H

#include <Decimation/Decim/sdecim.h>
#include <sys/types.h>

typedef struct {
  float *below, /* errors for each face, for one variable */
      *above;
} SDecimError, *SDecimErrorP;

typedef struct {
  SDecimErrorP err; /* errors for each variable */
} SDecimErrors, *SDecimErrorsP;

typedef struct SDecimMesh {
  VTmeshP mesh;       /* the unstructured result mesh */
  SDecimErrorsP *err; /* error values for faces, over each time */

  /* pointers into data for current time */
  VTsurf3DmdataP tdata; /* current time data */
  SDecimErrorsP curerr;
  int nvert, nface, usedfaces;
  int *vtof;

  /* loop structure */
  int loopsize; /* allocated size of arrays    */
  int nv;       /* number of verts in loop     */
  int *verts;   /* indices of vertices in loop */
  int *face;    /* indices of faces in loop    */
  int *adj;     /* adj faces around loop       */
  int *vprev, *vnext;
  float *preverr;
  SDecimErrorP lerr;   /* errors for variables        */
  SDecimError looperr; /* error for a single face     */
  SDecimError lasterr;
  float loopnorm[3];

  int nedges;
  int (*edgeadj)[2];
  int (*edgeij)[2];
  float *edgecost;
  float (*edgenorm)[3];
  SDecimErrorP edgeerr;

  int nf;
  SDecimErrorP ferr;
  int *newfaces;
  int (*newadj)[3];
  int (*newf)[3];

  /* deleted indices in mesh */
  int dvert;
  int dtri;
} SDecimMesh;

int SDecimMeshDecimateVert(SDecimMeshP, int, float *, int, float, float);
int SDecimMeshRemoveSimpleVert(SDecimMeshP, int, float *, int, float, float);
int SDecimMeshRemoveEdgeVert(SDecimMeshP, int, float *, int, float, float);
void SDecimMeshFreeVert(SDecimMeshP, int);
void SDecimMeshFreeFace(SDecimMeshP, int);
void SDecimMeshGetPlaneEq(SDecimMeshP, int, int, float[3], double[4]);
double SDecimMeshPtPlane(SDecimMeshP, int, double[4]);

#endif
