#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <Decimation/Decim/sdecimPriv.h>

#define MAX2(x, y) ((x) > (y) ? (x) : (y))
#define SQR(x) ((x) * (x))

#define DEBUG_EDGENo

int SDecimMeshRemoveEdgeVert(SDecimMeshP dmp, int v, float *err, int boundangle,
                             float aspect, float dihedral) {
  double val1, val2, val, neww, newerr;
  double splane[4], n[3], p1[3], *p2;
  int i;

  if (dmp->nv < 3)
    return (0);

  /*
   * take the plane perpendicular to the "new" edge, through the
   * edge point to be deleted
   */
  splane[0] = dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt[0] -
              dmp->tdata->verts[dmp->verts[0]].pt[0];
  splane[1] = dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt[1] -
              dmp->tdata->verts[dmp->verts[0]].pt[1];
  splane[2] = dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt[2] -
              dmp->tdata->verts[dmp->verts[0]].pt[2];
  splane[3] = -splane[0] * dmp->tdata->verts[v].pt[0] -
              splane[1] * dmp->tdata->verts[v].pt[1] -
              splane[2] * dmp->tdata->verts[v].pt[2];

  /* use this plane to interpolate between old mesh and new */
  val1 = fabs(SDecimMeshPtPlane(dmp, dmp->verts[0], splane));
  val2 = fabs(SDecimMeshPtPlane(dmp, dmp->verts[dmp->nv - 1], splane));
  val = val1 / (val1 + val2);

  /* interpolate each variable to obtain added error */
  for (i = 0; i < VTmeshNvars(dmp->mesh); i++) {
    dmp->lerr[i].above[dmp->nv - 1] =
        MAX2(dmp->lerr[i].above[0], dmp->lerr[i].above[dmp->nv - 2]);
    dmp->lerr[i].below[dmp->nv - 1] =
        MAX2(dmp->lerr[i].below[0], dmp->lerr[i].below[dmp->nv - 2]);

    switch (VTmeshVarType(dmp->mesh, i)) {
    case VT_UCHAR:
      neww = (1.0 - val) * dmp->tdata->data[i]->data.ucdata[dmp->verts[0]] +
             val * dmp->tdata->data[i]->data.ucdata[dmp->verts[dmp->nv - 1]];
      newerr = dmp->tdata->data[i]->data.ucdata[v] - neww;
      break;
    case VT_SHORT:
      neww = (1.0 - val) * dmp->tdata->data[i]->data.sdata[dmp->verts[0]] +
             val * dmp->tdata->data[i]->data.sdata[dmp->verts[dmp->nv - 1]];
      newerr = dmp->tdata->data[i]->data.sdata[v] - neww;
      break;
    case VT_FLOAT:
      neww = (1.0 - val) * dmp->tdata->data[i]->data.fdata[dmp->verts[0]] +
             val * dmp->tdata->data[i]->data.fdata[dmp->verts[dmp->nv - 1]];
      newerr = dmp->tdata->data[i]->data.fdata[v] - neww;
      break;
    }
    if (newerr < 0) {
      dmp->lerr[i].below[dmp->nv - 1] -= newerr;
      if (dmp->lerr[i].below[dmp->nv - 1] > err[i]) {
        /*
        printf("error below edge (%d) is %lf\n", v,
        dmp->lerr[i].below[dmp->nv-1]);
        */
        return (0);
      }
    } else {
      dmp->lerr[i].above[dmp->nv - 1] += newerr;
      if (dmp->lerr[i].above[dmp->nv - 1] > err[i]) {
        /*
        printf("error above edge (%d) is %lf\n", v,
        dmp->lerr[i].above[dmp->nv-1]);
        */
        return (0);
      }
    }
  }

  /* compute propagated geometric error */
  dmp->lerr[VTmeshNvars(dmp->mesh)].above[dmp->nv - 1] =
      MAX2(dmp->lerr[VTmeshNvars(dmp->mesh)].above[0],
           dmp->lerr[VTmeshNvars(dmp->mesh)].above[dmp->nv - 2]);
  dmp->lerr[VTmeshNvars(dmp->mesh)].below[dmp->nv - 1] =
      MAX2(dmp->lerr[VTmeshNvars(dmp->mesh)].below[0],
           dmp->lerr[VTmeshNvars(dmp->mesh)].below[dmp->nv - 2]);

  /* compute added geometric error */
  for (i = 0; i < 3; i++) {
    p1[i] = (1.0 - val) * dmp->tdata->verts[dmp->verts[0]].pt[i] +
            val * dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt[i];
  }
  p2 = dmp->tdata->verts[v].pt;
  newerr = sqrt(SQR(p1[0] - p2[0]) + SQR(p1[1] - p2[1]) + SQR(p1[2] - p2[2]));
  if (splane[0] * p1[0] + splane[1] * p1[1] + splane[2] * p1[2] + splane[3] <
      0.0) {
    dmp->lerr[VTmeshNvars(dmp->mesh)].below[dmp->nv - 1] += newerr;
    if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[dmp->nv - 1] >
        err[VTmeshNvars(dmp->mesh)])
      return (0);
  } else {
    dmp->lerr[VTmeshNvars(dmp->mesh)].above[dmp->nv - 1] += newerr;
    if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[dmp->nv - 1] >
        err[VTmeshNvars(dmp->mesh)])
      return (0);
  }

  /* boundaries of edge pass error constraints, triangulate as
     for simple vertices */
  dmp->face[dmp->nv - 1] = -1;
  dmp->adj[dmp->nv - 1] = -1;
#ifdef DEBUG_EDGE
  printf("Removing Edge:\n");
  for (i = 0; i < dmp->nv; i++) {
    printf(" %3d  %3d  %3d\n", dmp->verts[i], dmp->face[i], dmp->adj[i]);
    printf(" %lf %lf\n", dmp->lerr[0].above[i], dmp->lerr[0].below[i]);
  }
#endif
  return (
      SDecimMeshRemoveSimpleVert(dmp, v, err, boundangle, aspect, dihedral));
#ifdef DEBUG_EDGE
  sleep(0);
#endif
}
