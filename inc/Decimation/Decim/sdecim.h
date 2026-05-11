#ifndef SDECIM_H
#define SDECIM_H

#include <sys/types.h>
#include <Decimation/Decim/mesh.h>
typedef struct SDecimMesh *SDecimMeshP;

SDecimMeshP SDecimMeshCreate(VTmeshP);
void SDecimMeshClear(SDecimMeshP, int);
void SDecimMeshBegin(SDecimMeshP, int, int, int);
void SDecimMeshEnd(SDecimMeshP, int);
int SDecimMeshAddVert(SDecimMeshP, double, double, double, float *);
int SDecimMeshAddFace(SDecimMeshP, int[3], int[3], float *, float *);
int SDecimMeshSetFace(SDecimMeshP, int, int[3], int[3], float *, float *);
void SDecimMeshGetFaces(SDecimMeshP, int *, int);
void SDecimMeshDecimateVertIndex(SDecimMeshP, int *, int, int, float *, float,
                                 float);
void SDecimMeshDecimate(SDecimMeshP, float *, int, float, float);

#endif
