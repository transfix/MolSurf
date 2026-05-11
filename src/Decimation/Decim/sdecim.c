#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <string.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <Decimation/Decim/mesh.h>
#include <Decimation/Decim/sdecimPriv.h>

#define PROP_ERR
static int ncorner, nintedge, nboundary, nsimple, ncomplex;

#define max3(a, b, c)                                                          \
  ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))
#define max2(a, b) ((a) > (b) ? (b) : (a))

#define CORNER 0
#define INT_EDGE 1
#define BOUNDARY 2
#define SIMPLE 3
#define COMPLEX 4

/*-------------------------------------------------------------
 * SDecimMeshCreate - initialize a decimation mesh
 *-------------------------------------------------------------
 */
SDecimMeshP SDecimMeshCreate(VTmeshP mesh) {
  SDecimMeshP dmp = (SDecimMeshP)malloc(sizeof(SDecimMesh));
  int i;

  /* initialize to some default values */
  dmp->mesh = mesh;

  /* currently not decimating anything */
  dmp->tdata = NULL;
  dmp->curerr = NULL;

  /* allocate the errors array, but don't put anything in yet */
  dmp->err = (SDecimErrorsP *)malloc(sizeof(SDecimErrorsP) * VTmeshNtime(mesh));
  memset(dmp->err, '\0', sizeof(SDecimErrorsP) * VTmeshNtime(mesh));

  dmp->loopsize = 80;
  dmp->nv = 0;
  dmp->verts = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->face = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->adj = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->vprev = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->vnext = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->lerr =
      (SDecimErrorP)malloc(sizeof(SDecimError) * (VTmeshNvars(mesh) + 1));
  for (i = 0; i < VTmeshNvars(mesh) + 1; i++) {
    dmp->lerr[i].above = (float *)malloc(sizeof(float) * dmp->loopsize);
    dmp->lerr[i].below = (float *)malloc(sizeof(float) * dmp->loopsize);
  }
  dmp->preverr = (float *)malloc(sizeof(float) * (VTmeshNvars(mesh) + 1));
  dmp->lasterr.below = (float *)malloc(sizeof(float) * (VTmeshNvars(mesh) + 1));
  dmp->lasterr.above = (float *)malloc(sizeof(float) * (VTmeshNvars(mesh) + 1));

  dmp->looperr.below = (float *)malloc(sizeof(float) * (VTmeshNvars(mesh) + 1));
  dmp->looperr.above = (float *)malloc(sizeof(float) * (VTmeshNvars(mesh) + 1));

  dmp->nedges = 0;
  dmp->edgeerr =
      (SDecimErrorP)malloc(sizeof(SDecimError) * (VTmeshNvars(mesh) + 1));
  for (i = 0; i < VTmeshNvars(mesh) + 1; i++) {
    dmp->edgeerr[i].above = (float *)malloc(sizeof(float) * dmp->loopsize);
    dmp->edgeerr[i].below = (float *)malloc(sizeof(float) * dmp->loopsize);
  }
  dmp->edgecost = (float *)malloc(sizeof(float) * dmp->loopsize);
  dmp->edgeadj = (int(*)[2])malloc(sizeof(int[2]) * dmp->loopsize);
  dmp->edgeij = (int(*)[2])malloc(sizeof(int[2]) * dmp->loopsize);
  dmp->edgenorm = (float(*)[3])malloc(sizeof(float[3]) * dmp->loopsize);

  dmp->newfaces = (int *)malloc(sizeof(int) * dmp->loopsize);
  dmp->newadj = (int(*)[3])malloc(sizeof(int[3]) * dmp->loopsize);
  dmp->newf = (int(*)[3])malloc(sizeof(int[3]) * dmp->loopsize);
  dmp->ferr =
      (SDecimErrorP)malloc(sizeof(SDecimError) * (VTmeshNvars(mesh) + 1));
  for (i = 0; i < VTmeshNvars(mesh) + 1; i++) {
    dmp->ferr[i].below = (float *)malloc(sizeof(float) * dmp->loopsize);
    dmp->ferr[i].above = (float *)malloc(sizeof(float) * dmp->loopsize);
  }

  dmp->dvert = -1;
  dmp->dtri = -1;
  return (dmp);
}

/*-------------------------------------------------------------
 * SDecimMeshClear - clear all results for a decimated mesh
 *-------------------------------------------------------------
 */
void SDecimMeshClear(SDecimMeshP dmp, int time) {
  int i;

  if (dmp->err[time] != NULL) {
    if (dmp->err[time]->err != NULL) {
      for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
        free(dmp->err[time]->err[i].below);
        free(dmp->err[time]->err[i].above);
      }
      free(dmp->err[time]->err);
    }
    free(dmp->err[time]);
    dmp->err[time] = NULL;
  }
  dmp->curerr = NULL;
}

/*-------------------------------------------------------------
 * SDecimMeshBegin - prepare to decimate mesh for a given time
 *-------------------------------------------------------------
 */
void SDecimMeshBegin(SDecimMeshP dmp, int time, int vapprox, int fapprox) {
  int i;

  /* create the error structure for this timestep */
  if (dmp->err[time] == NULL) {
    dmp->curerr = dmp->err[time] = (SDecimErrorsP)malloc(sizeof(SDecimErrors));
    dmp->curerr->err = (SDecimErrorP)malloc(sizeof(SDecimError) *
                                            (VTmeshNvars(dmp->mesh) + 1));
    for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
      dmp->curerr->err[i].below = (float *)malloc(sizeof(float) * fapprox);
      dmp->curerr->err[i].above = (float *)malloc(sizeof(float) * fapprox);
    }
  }

  /* create the mesh data for this timestep */
  dmp->tdata = dmp->mesh->mesh.surf3d->data[time] =
      (VTsurf3DmdataP)malloc(sizeof(VTsurf3Dmdata));
  dmp->tdata->nvert = 0;
  dmp->tdata->ntri = 0;
  dmp->tdata->verts = (VTsurf3DvdataP)malloc(sizeof(VTsurf3Dvdata) * vapprox);
  dmp->vtof = (int *)malloc(sizeof(int) * vapprox);
  dmp->tdata->tri = (VTsurf3DtdataP)malloc(sizeof(VTsurf3Dtdata) * fapprox);
  dmp->tdata->data = (VTsurf3DdataP *)malloc(sizeof(VTsurf3DdataP) *
                                             (VTmeshNvars(dmp->mesh) + 1));
  for (i = 0; i < VTmeshNvars(dmp->mesh); i++) {
    dmp->tdata->data[i] = (VTsurf3DdataP)malloc(sizeof(VTsurf3Ddata));
    switch (VTmeshVarType(dmp->mesh, i)) {
    case VT_UCHAR:
      dmp->tdata->data[i]->data.ucdata =
          (unsigned char *)malloc(sizeof(unsigned char) * vapprox);
      break;
    case VT_CHAR:
      dmp->tdata->data[i]->data.cdata = (char *)malloc(sizeof(char) * vapprox);
      break;
    case VT_SHORT:
      dmp->tdata->data[i]->data.sdata =
          (short *)malloc(sizeof(short) * vapprox);
      break;
    case VT_LONG:
      dmp->tdata->data[i]->data.ldata = (long *)malloc(sizeof(long) * vapprox);
      break;
    case VT_FLOAT:
      dmp->tdata->data[i]->data.fdata =
          (float *)malloc(sizeof(float) * vapprox);
      break;
    default:
      fprintf(stderr, "Can't handle data type %d in decimation\n",
              VTmeshVarType(dmp->mesh, i));
    }
  }
  dmp->nvert = vapprox;
  dmp->nface = fapprox;
  dmp->usedfaces = 0;
}

void SDecimMeshCompress(SDecimMeshP dmp) {
  int i, j, deleted, copyto;
  int *vmap, *tmap;

  vmap = (int *)malloc(sizeof(int) * dmp->tdata->nvert);
  memset(vmap, 0, sizeof(int) * dmp->tdata->nvert);

  /* mark the deleted vertices */
  for (i = dmp->dvert; i != -1; i = (int)dmp->tdata->verts[i].pt[0])
    vmap[i] = -1;

  deleted = 0;
  copyto = 0;
  for (i = 0; i < dmp->tdata->nvert; i++) {
    if (vmap[i] == -1)
      deleted++;
    else {
      vmap[i] = copyto;
      memcpy(&dmp->tdata->verts[copyto], &dmp->tdata->verts[i],
             sizeof(VTsurf3Dvdata));
      copyto++;
    }
  }
  for (i = 0; i < VTmeshNvars(dmp->mesh); i++) {
    switch (VTmeshVarType(dmp->mesh, i)) {
    case VT_UCHAR:
      for (j = 0; j < dmp->tdata->nvert; j++) {
        if (vmap[j] != -1) {
          dmp->tdata->data[i]->data.ucdata[vmap[j]] =
              dmp->tdata->data[i]->data.ucdata[j];
        }
      }
      break;
    case VT_SHORT:
      for (j = 0; j < dmp->tdata->nvert; j++) {
        if (vmap[j] != -1) {
          dmp->tdata->data[i]->data.sdata[vmap[j]] =
              dmp->tdata->data[i]->data.sdata[j];
        }
      }
      break;
    case VT_FLOAT:
      for (j = 0; j < dmp->tdata->nvert; j++) {
        if (vmap[j] != -1) {
          dmp->tdata->data[i]->data.fdata[vmap[j]] =
              dmp->tdata->data[i]->data.fdata[j];
        }
      }
      break;
    }
  }
  printf("%d verts, %d deleted, %d total\n", dmp->tdata->nvert, deleted,
         dmp->tdata->nvert - deleted);
  dmp->dvert = -1;
  /* put all deleted vertices on the new deleted list */
  if (deleted == 0)
    dmp->dvert = -1;
  else {
    dmp->dvert = dmp->tdata->nvert - deleted;
    for (i = dmp->tdata->nvert - deleted; i < dmp->tdata->nvert - 1; i++)
      dmp->tdata->verts[i].pt[0] = (float)(i + 1);
    dmp->tdata->verts[dmp->tdata->nvert - 1].pt[0] = -1;
  }
  dmp->tdata->nvert -= deleted;

  tmap = (int *)malloc(sizeof(int) * dmp->tdata->ntri);
  memset(tmap, '\0', sizeof(int) * dmp->tdata->ntri);

  /* mark the deleted triangles */
  for (i = dmp->dtri; i != -1; i = dmp->tdata->tri[i].vert[1])
    tmap[i] = -1;

  deleted = 0;
  copyto = 0;
  i = 0;

  for (i = 0; i < dmp->tdata->ntri; i++) {
    /*  WHY THE -2???  shouldn't have this case here..
          if (tmap[i] == -1 || dmp->tdata->tri[i].vert[0] == -2) {
             tmap[i] = -1;
             deleted++;
          }
    */
    if (tmap[i] == -1)
      deleted++;
    else {
      tmap[i] = copyto;
      /* map the verts of this tri to the compressed mesh */
      dmp->tdata->tri[i].vert[0] = vmap[dmp->tdata->tri[i].vert[0]];
      dmp->tdata->tri[i].vert[1] = vmap[dmp->tdata->tri[i].vert[1]];
      dmp->tdata->tri[i].vert[2] = vmap[dmp->tdata->tri[i].vert[2]];
      if (dmp->tdata->tri[i].vert[0] == -1 ||
          dmp->tdata->tri[i].vert[1] == -1 ||
          dmp->tdata->tri[i].vert[2] == -1) {
        printf("tri references deleted vert\n");
        tmap[i] = -1;
        deleted++;
        continue;
      }
      /*
      if (dmp->tdata->tri[i].vert[0] > 50000) {
         fprintf(stderr, "warning: vert 0 maps to %d\n",
                 dmp->tdata->tri[i].vert[0]);
         sleep(0);
      }
      if (dmp->tdata->tri[i].vert[1] > 50000) {
         fprintf(stderr, "warning: vert 1 maps to %d\n",
                 dmp->tdata->tri[i].vert[1]);
         sleep(0);
      }
      if (dmp->tdata->tri[i].vert[2] > 50000) {
         fprintf(stderr, "warning: vert 2 maps to %d\n",
                 dmp->tdata->tri[i].vert[2]);
         sleep(0);
      }
      */
      /* and move this triangle back */
      memcpy(&dmp->tdata->tri[copyto], &dmp->tdata->tri[i],
             sizeof(VTsurf3Dtdata));
#define COMPRESS_ALL
#ifdef COMPRESS_ALL
      for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
        dmp->curerr->err[j].below[copyto] = dmp->curerr->err[j].below[i];
        dmp->curerr->err[j].above[copyto] = dmp->curerr->err[j].above[i];
      }
      dmp->vtof[dmp->tdata->tri[copyto].vert[0]] = copyto;
      dmp->vtof[dmp->tdata->tri[copyto].vert[1]] = copyto;
      dmp->vtof[dmp->tdata->tri[copyto].vert[2]] = copyto;
#endif
      copyto++;
    }
  }
  printf("%d tri, %d deleted, %d total\n", dmp->tdata->ntri, deleted,
         dmp->tdata->ntri - deleted);
  dmp->dtri = -1;
  dmp->usedfaces = dmp->tdata->ntri - deleted;
#if 1
  if (deleted == 0)
    dmp->dtri = -1;
  else {
    dmp->dtri = dmp->tdata->ntri - deleted;
    for (i = dmp->tdata->ntri - deleted; i < dmp->tdata->ntri - 1; i++) {
      dmp->tdata->tri[i].vert[0] = -1;
      dmp->tdata->tri[i].vert[1] = (float)(i + 1);
    }
    dmp->tdata->tri[dmp->tdata->ntri - 1].vert[0] = -1;
    dmp->tdata->tri[dmp->tdata->ntri - 1].vert[1] = -1;
  }
#endif
  dmp->tdata->ntri -= deleted;

  /* finally change the triangle adjacencies to remain consistent */
  for (i = 0; i < dmp->tdata->ntri; i++) {
    if (dmp->tdata->tri[i].adjtri[0] != -1)
      dmp->tdata->tri[i].adjtri[0] = tmap[dmp->tdata->tri[i].adjtri[0]];
    if (dmp->tdata->tri[i].adjtri[1] != -1)
      dmp->tdata->tri[i].adjtri[1] = tmap[dmp->tdata->tri[i].adjtri[1]];
    if (dmp->tdata->tri[i].adjtri[2] != -1)
      dmp->tdata->tri[i].adjtri[2] = tmap[dmp->tdata->tri[i].adjtri[2]];
  }

  free(vmap);
  free(tmap);
}

void SDecimMeshEnd(SDecimMeshP dmp, int time) {
  SDecimMeshCompress(dmp);

  /* unlink the lists of deleted vertices */
  dmp->dvert = -1;
  dmp->dtri = -1;

  dmp->curerr = NULL;
  dmp->nvert = 0;
  dmp->nface = 0;
  dmp->tdata = NULL;
}

/*-------------------------------------------------------------
 * SDecimMeshAddVert - add a vertex to a decimated mesh
 *-------------------------------------------------------------
 */
int SDecimMeshAddVert(SDecimMeshP dmp, double x, double y, double z,
                      float *data) {
  int n;
  int i;

  if (dmp->dvert != -1) {
    /* unlink a free vertex from the linked list */
    n = dmp->dvert;
    dmp->dvert = (int)dmp->tdata->verts[n].pt[0];
  } else {
    /* reallocate if we need more verts */
    if ((n = dmp->tdata->nvert++) >= dmp->nvert) {
      dmp->nvert *= 2;
      dmp->tdata->verts = (VTsurf3DvdataP)realloc(
          dmp->tdata->verts, sizeof(VTsurf3Dvdata) * dmp->nvert);
      dmp->vtof = (int *)realloc(dmp->vtof, sizeof(int) * dmp->nvert);
      for (i = 0; i < VTmeshNvars(dmp->mesh); i++) {
        switch (VTmeshVarType(dmp->mesh, i)) {
        case VT_UCHAR:
          dmp->tdata->data[i]->data.ucdata =
              (unsigned char *)realloc(dmp->tdata->data[i]->data.ucdata,
                                       sizeof(unsigned char) * dmp->nvert);
          break;
        case VT_CHAR:
          dmp->tdata->data[i]->data.cdata = (char *)realloc(
              dmp->tdata->data[i]->data.cdata, sizeof(char) * dmp->nvert);
          break;
        case VT_SHORT:
          dmp->tdata->data[i]->data.sdata = (short *)realloc(
              dmp->tdata->data[i]->data.sdata, sizeof(short) * dmp->nvert);
          break;
        case VT_LONG:
          dmp->tdata->data[i]->data.ldata = (long *)realloc(
              dmp->tdata->data[i]->data.ldata, sizeof(long) * dmp->nvert);
          break;
        case VT_FLOAT:
          dmp->tdata->data[i]->data.fdata = (float *)realloc(
              dmp->tdata->data[i]->data.fdata, sizeof(float) * dmp->nvert);
          break;
        default:
          fprintf(stderr, "Can't handle data type %d in decimation\n",
                  VTmeshVarType(dmp->mesh, i));
        }
      }
    }
  }

  dmp->tdata->verts[n].pt[0] = x;
  dmp->tdata->verts[n].pt[1] = y;
  dmp->tdata->verts[n].pt[2] = z;
  dmp->vtof[n] = -1;

  for (i = 0; i < VTmeshNvars(dmp->mesh); i++) {
    switch (VTmeshVarType(dmp->mesh, i)) {
    case VT_UCHAR:
      dmp->tdata->data[i]->data.ucdata[n] = data[i];
      break;
    case VT_CHAR:
      dmp->tdata->data[i]->data.cdata[n] = data[i];
      break;
    case VT_SHORT:
      dmp->tdata->data[i]->data.sdata[n] = data[i];
      break;
    case VT_LONG:
      dmp->tdata->data[i]->data.ldata[n] = data[i];
      break;
    case VT_FLOAT:
      dmp->tdata->data[i]->data.fdata[n] = data[i];
      break;
    default:
      fprintf(stderr, "Can't handle data type %d in decimation\n",
              VTmeshVarType(dmp->mesh, i));
    }
  }

  return (n);
}

/*-------------------------------------------------------------
 * SDecimMeshFreeVert - add a vertex to the internal free list
 *                      it is assumed that no triangles in the
 *                      mesh reference the vertex
 *-------------------------------------------------------------
 */
void SDecimMeshFreeVert(SDecimMeshP dmp, int v) {
  dmp->tdata->verts[v].pt[0] = (double)dmp->dvert;
  dmp->dvert = v;
}

/*-------------------------------------------------------------
 * SDecimMeshSetFace  - set the vertices, adjacent triangles
 *                      and error values for a face
 *-------------------------------------------------------------
 */
int SDecimMeshSetFace(SDecimMeshP dmp, int f, int v[3], int adj[3],
                      float *above, float *below) {
  int i;

  if (f >= dmp->tdata->ntri)
    dmp->tdata->ntri = f + 1;

#ifdef DEBUG
  if (v[0] == v[1] || v[0] == v[2] || v[1] == v[2]) {
    fprintf(stderr, "adding bad face!\n");
    sleep(0);
  }
  if (adj[0] < -1 || adj[1] < -1 || adj[2] < -1) {
    fprintf(stderr, "adding neg adj!\n");
    sleep(0);
  }
  if ((adj[0] == adj[1] && adj[0] != -1) ||
      (adj[0] == adj[2] && adj[0] != -1) ||
      (adj[1] == adj[2] && adj[1] != -1)) {
    fprintf(stderr, "adding bad adj!\n");
    sleep(0);
  }
#endif
  dmp->vtof[v[0]] = f;
  dmp->vtof[v[1]] = f;
  dmp->vtof[v[2]] = f;
  memcpy(dmp->tdata->tri[f].vert, v, sizeof(int) * 3);
  memcpy(dmp->tdata->tri[f].adjtri, adj, sizeof(int) * 3);

#ifdef PROP_ERR
  if (above == NULL) {
    for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
      dmp->curerr->err[i].above[f] = 0.0;
      dmp->curerr->err[i].below[f] = 0.0;
    }
  } else {
    for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
      dmp->curerr->err[i].above[f] = above[i];
      dmp->curerr->err[i].below[f] = below[i];
    }
  }
#else
  for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
    dmp->curerr->err[i].above[f] = 0.0;
    dmp->curerr->err[i].below[f] = 0.0;
  }
#endif
}

/*-------------------------------------------------------------
 * SDecimMeshAddFace  - add a face
 *-------------------------------------------------------------
 */
int SDecimMeshAddFace(SDecimMeshP dmp, int v[3], int adj[3], float *above,
                      float *below) {
  int n;
  int i;

  if (dmp->dtri != -1) {
    n = dmp->dtri;
    dmp->dtri = dmp->tdata->tri[n].vert[1];
  } else {
    if ((n = dmp->usedfaces++) >= dmp->nface) {
      dmp->nface *= 2;
      dmp->tdata->tri = (VTsurf3DtdataP)realloc(
          dmp->tdata->tri, sizeof(VTsurf3Dtdata) * dmp->nface);
      for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
        dmp->curerr->err[i].below = (float *)realloc(
            dmp->curerr->err[i].below, sizeof(float) * dmp->nface);
        dmp->curerr->err[i].above = (float *)realloc(
            dmp->curerr->err[i].above, sizeof(float) * dmp->nface);
      }
    }
  }

  SDecimMeshSetFace(dmp, n, v, adj, above, below);

  return (n);
}

/*-------------------------------------------------------------
 * SDecimMeshFreeFace  - add a face to the internal free list
 *                      it is assumed that the face is no longer
 *                      used
 *-------------------------------------------------------------
 */
void SDecimMeshFreeFace(SDecimMeshP dmp, int f) {
  dmp->tdata->tri[f].vert[0] = -1;
  dmp->tdata->tri[f].vert[1] = dmp->dtri;
  dmp->dtri = f;
}

/*-------------------------------------------------------------
 * SDecimMeshGetFaces  - reserve 'nf' faces for use and return
 *                       the indices for them.  This is provided
 *                       so the user can reserve enough faces to
 *                       construct the adjacency information needed
 *                       to add faces to the mesh.  After this
 *                       function is used, SDecimMeshSetFace should
 *                       be called to fill in the information for
 *                       each face allocated.
 *-------------------------------------------------------------
 */
void SDecimMeshGetFaces(SDecimMeshP dmp, int *faces, int nf) {
  int i, j, remain;

  i = 0;

#ifdef DEBUG
  printf("Getting faces: ");
#endif
  /* take as many as possible from the free list of triangles */
  while (i < nf && dmp->dtri != -1) {
#ifdef DEBUG
    printf("%d ", dmp->dtri);
#endif
    faces[i++] = dmp->dtri;
    dmp->dtri = dmp->tdata->tri[dmp->dtri].vert[1];
  }

  remain = nf - i;
#ifdef DEBUG_
  printf("Getting faces: %d from free list, %d new\n", i, remain);
  if (remain)
    printf("%d->%d\n", dmp->usedfaces, dmp->usedfaces + nf - i - 1);
  else
    putchar('\n');
#endif
  for (j = 0; j < remain; j++)
    faces[i++] = dmp->usedfaces + j;

  if (remain > 0) {
    dmp->usedfaces += remain;

    if (dmp->usedfaces >= dmp->nface) {
      dmp->nface *= 2;
      dmp->tdata->tri = (VTsurf3DtdataP)realloc(
          dmp->tdata->tri, sizeof(VTsurf3Dtdata) * dmp->nface);
      for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
        dmp->curerr->err[i].below = (float *)realloc(
            dmp->curerr->err[i].below, sizeof(float) * dmp->nface);
        dmp->curerr->err[i].above = (float *)realloc(
            dmp->curerr->err[i].above, sizeof(float) * dmp->nface);
      }
    }
  }

  /* mark the faces as allocated, but not filled */
  for (i = 0; i < nf; i++)
    dmp->tdata->tri[faces[i]].vert[0] = -2;
}

/*-------------------------------------------------------------
 * SDecimMeshDecimateVertIndex  - Decimate the vertices in the
 *                       array 'v', from index 'i' to 'j'.
 *-------------------------------------------------------------
 */
void SDecimMeshDecimateVertIndex(SDecimMeshP dmp, int *v, int i, int j,
                                 float *err, float aspect, float dihedral) {
  int start, skip, k, lastdone;

  /*
  printf("decimating from %d to %d\n", v[i], v[j]);
  */
  start = i;
  skip = 2;
  lastdone = 0;
  while (start < j) {
    for (k = start; k <= j; k += skip)
      SDecimMeshDecimateVert(dmp, v[k], err, 1.0, aspect, dihedral);
    if (k == j + skip)
      lastdone = 1;
    start += (skip >> 1);
    skip <<= 1;
  }
  if (!lastdone)
    SDecimMeshDecimateVert(dmp, v[j], err, 1.0, aspect, dihedral);
}
