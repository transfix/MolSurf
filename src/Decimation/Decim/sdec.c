#include <stdio.h>
#include <unistd.h>
#include <math.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <memory.h>
#include <string.h>
#include <sys/time.h>
#include <Decimation/Decim/sdecimPriv.h>

#define SIMPLE 0
#define EDGE 1
#define COMPLEX 4

#define ERR_HACK

#define DEBUG_LOOP_No

#define SQR(x) ((x) * (x))

/*-------------------------------------------------------------
 * SDecimMeshFaceVertIndex - find index of a vertex on a face,
 *                           or return -1 if vertex not on face
 *-------------------------------------------------------------
 */
int SDecimMeshFaceVertIndex(SDecimMeshP dmp, int f, int v) {
  if (dmp->tdata->tri[f].vert[0] == v)
    return (0);
  if (dmp->tdata->tri[f].vert[1] == v)
    return (1);
  if (dmp->tdata->tri[f].vert[2] == v)
    return (2);
  return (-1);
}

void SDecimMeshResizeLoop(SDecimMeshP dmp) {
  int i, k;

  dmp->loopsize *= 2;
  dmp->verts = (int *)realloc(dmp->verts, sizeof(int) * dmp->loopsize);
  dmp->adj = (int *)realloc(dmp->adj, sizeof(int) * dmp->loopsize);
  dmp->vnext = (int *)realloc(dmp->vnext, sizeof(int) * dmp->loopsize);
  dmp->vprev = (int *)realloc(dmp->vprev, sizeof(int) * dmp->loopsize);
  dmp->face = (int *)realloc(dmp->face, sizeof(int) * dmp->loopsize);
  for (k = 0; k < VTmeshNvars(dmp->mesh) + 1; k++) {
    dmp->lerr[k].above =
        (float *)realloc(dmp->lerr[k].above, sizeof(float) * dmp->loopsize);
    dmp->lerr[k].below =
        (float *)realloc(dmp->lerr[k].below, sizeof(float) * dmp->loopsize);
  }

  for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
    dmp->edgeerr[i].above =
        (float *)realloc(dmp->edgeerr[i].above, sizeof(float) * dmp->loopsize);
    dmp->edgeerr[i].below =
        (float *)realloc(dmp->edgeerr[i].below, sizeof(float) * dmp->loopsize);
  }
  dmp->edgecost =
      (float *)realloc(dmp->edgecost, sizeof(float) * dmp->loopsize);
  dmp->edgeadj =
      (int(*)[2])realloc(dmp->edgeadj, sizeof(int[2]) * dmp->loopsize);
  dmp->edgeij = (int(*)[2])realloc(dmp->edgeij, sizeof(int[2]) * dmp->loopsize);
  dmp->edgenorm =
      (float(*)[3])realloc(dmp->edgenorm, sizeof(float[3]) * dmp->loopsize);

  dmp->newfaces = (int *)realloc(dmp->newfaces, sizeof(int) * dmp->loopsize);
  dmp->newadj = (int(*)[3])realloc(dmp->newadj, sizeof(int[3]) * dmp->loopsize);
  dmp->newf = (int(*)[3])realloc(dmp->newf, sizeof(int[3]) * dmp->loopsize);
  for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
    dmp->ferr[i].below =
        (float *)realloc(dmp->ferr[i].below, sizeof(float) * dmp->loopsize);
    dmp->ferr[i].above =
        (float *)realloc(dmp->ferr[i].above, sizeof(float) * dmp->loopsize);
  }
}

/*-------------------------------------------------------------
 * SDecimMeshGetLoop - get the loop around the given vertex
 *-------------------------------------------------------------
 */
int SDecimMeshGetLoop(SDecimMeshP dmp, int v) {
  int i, j, k, f, nextf, firstface, nv, n;
  double y1, y2, x1, x2;
  double *(p[4]);
  double len;

  nv = 0;

  /* find a face with which to begin */
  f = dmp->vtof[v];
  i = SDecimMeshFaceVertIndex(dmp, f, v);

  if (f < 0) {
    fprintf(stderr, "SDecimMeshGetLoop: couldn't find a starting place\n");
#ifdef DEBUG_LOOP
    sleep(0);
#endif
    return (COMPLEX);
  }

  firstface = f;
  nextf = dmp->tdata->tri[f].adjtri[i];

  /* loop around counterclockwise to find edge (if there is one) */
  j = 0;
  while (nextf != -1 && nextf != firstface) {
    i = SDecimMeshFaceVertIndex(dmp, nextf, v);
    f = nextf;
    nextf = dmp->tdata->tri[f].adjtri[i];
    j++;
    if (j > 200) {
      fprintf(stderr, "detected infinite loop, abort vertex\n");
      return (COMPLEX);
    }
  }

  if (nextf == -1)
    firstface = f;

    /* 'i' is the index of vertex v on face f */
#ifdef DEBUG_LOOP
  printf("starting face is %d\n", f);
#endif

  /* copy the first two vertices into the loop */
  if (i == 2)
    j = 0;
  else
    j = i + 1;
  dmp->verts[nv] = dmp->tdata->tri[f].vert[j];
#ifdef DEBUG_LOOP
  printf("vert %d\n", dmp->verts[nv]);
#endif
  dmp->adj[nv] = dmp->tdata->tri[f].adjtri[j];
  dmp->face[nv] = f;
  for (k = 0; k < VTmeshNvars(dmp->mesh) + 1; k++) {
    dmp->lerr[k].above[nv] = dmp->curerr->err[k].above[f];
    dmp->lerr[k].below[nv] = dmp->curerr->err[k].below[f];
  }
  nv++;

  if (j == 2)
    j = 0;
  else
    j++;
  dmp->verts[nv] = dmp->tdata->tri[f].vert[j];
#ifdef DEBUG_LOOP
  printf("vert %d\n", dmp->verts[nv]);
#endif
  nv++;

  firstface = f;
  f = dmp->tdata->tri[f].adjtri[j];

  while (f != -1 && f != firstface) {
#ifdef DEBUG_LOOP
    printf("moving to face %d (first = %d)\n", f, firstface);
#endif
    if (nv >= dmp->loopsize - 8)
      SDecimMeshResizeLoop(dmp);

    i = SDecimMeshFaceVertIndex(dmp, f, dmp->verts[nv - 1]);
    if (i == -1) {
      fprintf(stderr, "SDecimMeshGetLoop: next face not valid\n");
      sleep(0);
      return (COMPLEX);
    }
    dmp->face[nv - 1] = f;
    dmp->adj[nv - 1] = dmp->tdata->tri[f].adjtri[i];
    for (k = 0; k < VTmeshNvars(dmp->mesh) + 1; k++) {
      dmp->lerr[k].above[nv - 1] = dmp->curerr->err[k].above[f];
      dmp->lerr[k].below[nv - 1] = dmp->curerr->err[k].below[f];
    }
    if (i == 2)
      j = 0;
    else
      j = i + 1;

    dmp->verts[nv] = dmp->tdata->tri[f].vert[j];

    /* compute the edge normal for the last 3 verts */
    p[0] = dmp->tdata->verts[dmp->verts[nv - 2]].pt;
    p[1] = dmp->tdata->verts[dmp->verts[nv - 1]].pt;
    p[2] = dmp->tdata->verts[dmp->verts[nv - 0]].pt;
    p[3] = dmp->tdata->verts[v].pt;

    memset(dmp->edgenorm[nv - 1], '\0', sizeof(float) * 3);
    for (k = 3, n = 0; n < 4; k = n, n++) {
      dmp->edgenorm[nv - 1][0] += (p[k][1] + p[n][1]) * (p[k][2] - p[n][2]);
      dmp->edgenorm[nv - 1][1] += (p[k][2] + p[n][2]) * (p[k][0] - p[n][0]);
      dmp->edgenorm[nv - 1][2] += (p[k][0] + p[n][0]) * (p[k][1] - p[n][1]);
    }
    len = sqrt(SQR(dmp->edgenorm[nv - 1][0]) + SQR(dmp->edgenorm[nv - 1][1]) +
               SQR(dmp->edgenorm[nv - 1][2]));
    dmp->edgenorm[nv - 1][0] /= len;
    dmp->edgenorm[nv - 1][1] /= len;
    dmp->edgenorm[nv - 1][2] /= len;

#ifdef DEBUG_LOOP
    printf("vert %d\n", dmp->verts[nv]);
#endif
    nv++;
    f = dmp->tdata->tri[f].adjtri[j];

    if (nv > 200) {
      fprintf(stderr, "Hmm. nv is 200\n");
      sleep(0);
      return (COMPLEX);
    }
  }

  /* compute the normal for the loop */
  memset(dmp->loopnorm, '\0', sizeof(float) * 3);
  for (k = nv - 1, n = 0; n < nv; k = n, n++) {
    p[0] = dmp->tdata->verts[dmp->verts[k]].pt;
    p[1] = dmp->tdata->verts[dmp->verts[n]].pt;
    dmp->loopnorm[0] += (p[0][1] + p[1][1]) * (p[0][2] - p[1][2]);
    dmp->loopnorm[1] += (p[0][2] + p[1][2]) * (p[0][0] - p[1][0]);
    dmp->loopnorm[2] += (p[0][0] + p[1][0]) * (p[0][1] - p[1][1]);
  }
  len = sqrt(SQR(dmp->loopnorm[0]) + SQR(dmp->loopnorm[1]) +
             SQR(dmp->loopnorm[2]));
  dmp->loopnorm[0] /= len;
  dmp->loopnorm[1] /= len;
  dmp->loopnorm[2] /= len;

  /* see if we've completed a full loop */
  if (f == firstface) {
    dmp->nv = nv - 1;

    /* compute the edge normal for the last 3 verts */
    p[0] = dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt;
    p[1] = dmp->tdata->verts[dmp->verts[0]].pt;
    p[2] = dmp->tdata->verts[dmp->verts[1]].pt;
    p[3] = dmp->tdata->verts[v].pt;

    memset(dmp->edgenorm[0], '\0', sizeof(float) * 3);
    for (k = 3, n = 0; n < 4; k = n, n++) {
      dmp->edgenorm[0][0] += (p[k][1] + p[n][1]) * (p[k][2] - p[n][2]);
      dmp->edgenorm[0][1] += (p[k][2] + p[n][2]) * (p[k][0] - p[n][0]);
      dmp->edgenorm[0][2] += (p[k][0] + p[n][0]) * (p[k][1] - p[n][1]);
    }
    len = sqrt(SQR(dmp->edgenorm[0][0]) + SQR(dmp->edgenorm[0][1]) +
               SQR(dmp->edgenorm[0][2]));
    dmp->edgenorm[0][0] /= len;
    dmp->edgenorm[0][1] /= len;
    dmp->edgenorm[0][2] /= len;

    return (SIMPLE);
  }

  /* if not, this is a boundary face... return edge */
  dmp->nv = nv;

  p[0] = dmp->tdata->verts[dmp->verts[0]].pt;
  p[1] = dmp->tdata->verts[dmp->verts[1]].pt;
  p[2] = dmp->tdata->verts[v].pt;
  memset(dmp->edgenorm[0], '\0', sizeof(float) * 3);
  for (k = 2, n = 0; n < 3; k = n, n++) {
    dmp->edgenorm[0][0] += (p[k][1] + p[n][1]) * (p[k][2] - p[n][2]);
    dmp->edgenorm[0][1] += (p[k][2] + p[n][2]) * (p[k][0] - p[n][0]);
    dmp->edgenorm[0][2] += (p[k][0] + p[n][0]) * (p[k][1] - p[n][1]);
  }
  len = sqrt(SQR(dmp->edgenorm[0][0]) + SQR(dmp->edgenorm[0][1]) +
             SQR(dmp->edgenorm[0][2]));
  dmp->edgenorm[0][0] /= len;
  dmp->edgenorm[0][1] /= len;
  dmp->edgenorm[0][2] /= len;

  p[0] = dmp->tdata->verts[dmp->verts[dmp->nv - 2]].pt;
  p[1] = dmp->tdata->verts[dmp->verts[dmp->nv - 1]].pt;
  p[2] = dmp->tdata->verts[v].pt;
  memset(dmp->edgenorm[dmp->nv - 1], '\0', sizeof(float) * 3);
  for (k = 2, n = 0; n < 3; k = n, n++) {
    dmp->edgenorm[dmp->nv - 1][0] += (p[k][1] + p[n][1]) * (p[k][2] - p[n][2]);
    dmp->edgenorm[dmp->nv - 1][1] += (p[k][2] + p[n][2]) * (p[k][0] - p[n][0]);
    dmp->edgenorm[dmp->nv - 1][2] += (p[k][0] + p[n][0]) * (p[k][1] - p[n][1]);
  }
  len = sqrt(SQR(dmp->edgenorm[dmp->nv - 1][0]) +
             SQR(dmp->edgenorm[dmp->nv - 1][1]) +
             SQR(dmp->edgenorm[dmp->nv - 1][2]));
  dmp->edgenorm[dmp->nv - 1][0] /= len;
  dmp->edgenorm[dmp->nv - 1][1] /= len;
  dmp->edgenorm[dmp->nv - 1][2] /= len;

  return (EDGE);
}

/*-------------------------------------------------------------
 * SDecimMeshDecimateVert - attempt to decimate a mesh at a
 *                          given vertex, given error limits
 *-------------------------------------------------------------
 */
int SDecimMeshDecimateVert(SDecimMeshP dmp, int v, float *err, int boundangle,
                           float aspect, float dihedral) {
  int type;
  static int ncall = 0;
  static int simp = 0;
  static int edge = 0;
  static int corner = 0;

  /*
  printf("%d calls\n", ++ncall);
  */
  type = SDecimMeshGetLoop(dmp, v);
  switch (type) {
  case SIMPLE:
    /*
    simp++;
    printf("%d simple (%d %d)\n", simp, edge, corner);
    */
    return (
        SDecimMeshRemoveSimpleVert(dmp, v, err, boundangle, aspect, dihedral));
    break;
  case EDGE:
    /*
             edge++;
    printf("%d edges (%d %d)\n", edge, simp, corner);
    */
    return (
        SDecimMeshRemoveEdgeVert(dmp, v, err, boundangle, aspect, dihedral));
    break;
  case COMPLEX:
    /*
             corner++;
    printf("%d corners (%d %d)\n", corner, simp, edge);
    */
    break;
  default:
    printf("unknown vert type!!\n");
    sleep(0);
  }
  return (0);
}

void SDecimMeshDecimate(SDecimMeshP dmp, float *err, int niter, float aspect,
                        float dihedral) {
  int i, v, j, p, tmp;
  int *vindex;
  int *deleted, nd;
  struct timeval tp;
  struct timezone tzp;
  float *curerr;

  vindex = (int *)malloc(sizeof(int) * dmp->tdata->nvert);
  curerr = (float *)malloc(sizeof(float) * (VTmeshNvars(dmp->mesh) + 1));
  deleted = (int *)malloc(sizeof(int) * dmp->tdata->nvert);

#ifdef ERR_HACK
  /* if one variable, named initerr, initialize geometry error based on it */
  if (VTmeshNvars(dmp->mesh) == 1 &&
      !strcmp(VTmeshVarName(dmp->mesh, 0), "initerr")) {
    int i, v, vert;
    float init;
    for (i = 0; i < dmp->tdata->ntri; i++) {
      for (v = 0; v < 3; v++) {
        vert = dmp->tdata->tri[i].vert[v];
        init = dmp->tdata->data[0]->data.fdata[vert];
        if (dmp->curerr->err[1].above[i] < init * err[1])
          dmp->curerr->err[1].above[i] = init * err[1];
        if (dmp->curerr->err[1].below[i] < init * err[1])
          dmp->curerr->err[1].below[i] = init * err[1];
      }
    }
  }
#endif

  /* compute a random perturbation of the current vertices */
  gettimeofday(&tp, &tzp);
  srand48(tp.tv_usec);
  for (v = 0; v < dmp->tdata->nvert; v++)
    vindex[v] = v;
  for (v = dmp->tdata->nvert - 1; v > 0; v--) {
    j = lrand48() % (v + 1);
    tmp = vindex[v];
    vindex[v] = vindex[j];
    vindex[j] = tmp;
  }

  memset(deleted, 0, sizeof(int) * dmp->tdata->nvert);
  for (p = 0; p < niter; p++) {
    /* compute the current error value to use */
    for (v = 0; v < VTmeshNvars(dmp->mesh) + 1; v++)
      curerr[v] = SQR((p + 1.0) / niter) * err[v];

    /* decimate the verts */
    for (v = 0; v < dmp->tdata->nvert; v++)
      if (!deleted[vindex[v]])
        SDecimMeshDecimateVert(dmp, vindex[v], curerr, 1, aspect, dihedral);

    /* compress the deleted parts to start over again.. */
    if (p != niter - 1) {
      nd = 0;
      memset(deleted, 0, sizeof(int) * dmp->tdata->nvert);
      for (i = dmp->dvert; i != -1; nd++, i = dmp->tdata->verts[i].pt[0])
        deleted[i] = 1;
      printf("%d are deleted\n", nd);
      /*
               SDecimMeshCompress(dmp);
      */
    }
  }

  free(deleted);
  free(vindex);
  free(curerr);
}
