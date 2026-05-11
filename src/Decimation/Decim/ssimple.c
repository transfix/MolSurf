#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <Decimation/Decim/sdecimPriv.h>

#define DEBUG_SUMMARY
#define DEBUG_REJECTNo
#define DEBUG_LOOPNo
#define DEBUG_SWAPNo
#define DEBUG_TRINo
#define DEBUG_GEOMNo
#define SWAP

#ifdef DEBUG_SUMMARY
static int NATT = 0, NDEL = 0, NRJTRI = 0, NRJERR = 0, NRJFTR = 0, NRJTRI2 = 0,
           NRJTERR = 0;
#endif

#ifdef DEBUG_REJECT
static int NRJ = 0;
#endif

#define SWAP_BY_ANGLES

#define DIHEDRAL_CHECK
/* feature angle of 120 deg */
/*
#define DIHEDRAL       .3
*/
/* feature angle of 150 deg */
/*
#define DIHEDRAL       .5
*/
/* feature angle of 60 deg */
/*
#define DIHEDRAL       -.3
*/
/* #define DIHEDRAL       0.54 */
/*
#define DIHEDRAL       -0.7
*/

#define MINANG_CRITERIA
/*
#define MINANG 0.98
*/

#define MAX3(x, y, z)                                                          \
  ((x) > (y) ? ((x) > (z) ? (x) : (z)) : ((y) > (z) ? (y) : (z)))
#define SQR(x) ((x) * (x))

void SDecimMeshGetPlaneEq(SDecimMeshP dmp, int v1, int v2, float norm[3],
                          double eq[4]) {
  double n[3];

  n[0] = dmp->tdata->verts[v2].pt[0] - dmp->tdata->verts[v1].pt[0];
  n[1] = dmp->tdata->verts[v2].pt[1] - dmp->tdata->verts[v1].pt[1];
  n[2] = dmp->tdata->verts[v2].pt[2] - dmp->tdata->verts[v1].pt[2];

  eq[0] = n[1] * norm[2] - n[2] * norm[1];
  eq[1] = n[2] * norm[0] - n[0] * norm[2];
  eq[2] = n[0] * norm[1] - n[1] * norm[0];
  eq[3] = -eq[0] * dmp->tdata->verts[v1].pt[0] -
          eq[1] * dmp->tdata->verts[v1].pt[1] -
          eq[2] * dmp->tdata->verts[v1].pt[2];
}

double SDecimMeshPtPlane(SDecimMeshP dmp, int v, double eq[4]) {
  return (eq[0] * dmp->tdata->verts[v].pt[0] +
          eq[1] * dmp->tdata->verts[v].pt[1] +
          eq[2] * dmp->tdata->verts[v].pt[2] + eq[3]);
}

int SDecimMeshValidPlane(SDecimMeshP dmp, int v1, int v2, double eq[4]) {
  double val, val1, val2;
  int i;

  if (v1 > v2) {
    i = v1;
    v1 = v2;
    v2 = i;
  }

  /* first, see if all pts in (v1,v2) lie on same side */
  val1 = SDecimMeshPtPlane(dmp, dmp->verts[i = dmp->vnext[v1]], eq);

  for (i = dmp->vnext[i]; i != v2; i = dmp->vnext[i]) {
    val = SDecimMeshPtPlane(dmp, dmp->verts[i], eq);
    if (val1 * val <= 0.0)
      return (0);
  }

  /* next, check that side 1 and side 2 differ */
  val2 = SDecimMeshPtPlane(dmp, dmp->verts[i = dmp->vnext[v2]], eq);

  if (val1 * val2 >= 0.0)
    return (0);

  /* finally, check that all of side 2 is on same side */
  for (i = dmp->vnext[i]; i != v1; i = dmp->vnext[i]) {
    val = SDecimMeshPtPlane(dmp, dmp->verts[i], eq);
    if (val2 * val <= 0.0)
      return (0);
  }

  /* line is valid to split loop */
  return (1);
}

void SDecimMeshCompHoleErr(SDecimMeshP dmp, int v, int i, int j, int n,
                           double eq[4], float *err) {
  double vval, kval, ival, jval, kvval, ijval, kvijerr;
  double rplane[4];
  double p1[3], p2[3];
  float *above, *below;
  int k, l, var, prev;

  for (var = 0; var < VTmeshNvars(dmp->mesh) + 1; var++) {
    dmp->preverr[var] = 0.0;
    dmp->edgeerr[var].below[n] = 0.0;
    dmp->edgeerr[var].above[n] = 0.0;
  }

  dmp->edgecost[n] = 0.0;

  vval = SDecimMeshPtPlane(dmp, v, eq);

#ifdef DEBUG_LOOP
  printf("checking error on (%d->%d)\n", i, j);
#endif
  if (vval > 0.0) {
#ifdef DEBUG_LOOP
    printf("vval is positive\n");
#endif
    /* vertex lies to side of hole 1 */
    for (prev = i, k = (i + 1) % dmp->nv; k != j;
         prev = k, k = (k == dmp->nv - 1 ? 0 : k + 1)) {
#ifdef DEBUG_LOOP
      printf("testing error to vert %d\n", k);
#endif
      kval = fabs(SDecimMeshPtPlane(dmp, dmp->verts[k], eq));
      SDecimMeshGetPlaneEq(dmp, dmp->verts[k], v, dmp->edgenorm[k], rplane);
      ival = fabs(SDecimMeshPtPlane(dmp, dmp->verts[i], rplane));
      jval = fabs(SDecimMeshPtPlane(dmp, dmp->verts[j], rplane));
#ifdef DEBUG_LOOP
      printf("vval, kval, ival, jval = %lf, %lf, %lf, %lf\n", vval, kval, ival,
             jval);
#endif
      for (var = 0; var < VTmeshNvars(dmp->mesh); var++) {
        above = dmp->edgeerr[var].above;
        below = dmp->edgeerr[var].below;
        switch (VTmeshVarType(dmp->mesh, var)) {
        case VT_UCHAR:
          kvval =
              (1.0 - kval / (kval + vval)) *
                  dmp->tdata->data[var]->data.ucdata[dmp->verts[k]] +
              (kval / (kval + vval)) * dmp->tdata->data[var]->data.ucdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.ucdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.ucdata[dmp->verts[j]];
          break;
        case VT_SHORT:
          kvval = (1.0 - kval / (kval + vval)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[k]] +
                  (kval / (kval + vval)) * dmp->tdata->data[var]->data.sdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[j]];
          break;
        case VT_FLOAT:
          kvval = (1.0 - kval / (kval + vval)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[k]] +
                  (kval / (kval + vval)) * dmp->tdata->data[var]->data.fdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[j]];
          break;
        }
        kvijerr = kvval - ijval;
#ifdef DEBUG_LOOP
        printf("   var %d: kvijerr = %lf\n", var, kvijerr);
#endif
        if (kvijerr < dmp->preverr[var]) {
          if (dmp->lerr[var].below[prev] - kvijerr > below[n])
            below[n] = dmp->lerr[var].below[prev] - kvijerr;
          if (dmp->lerr[var].above[prev] + dmp->preverr[var] > above[n])
            above[n] = dmp->lerr[var].above[prev] + dmp->preverr[var];
        } else {
          if (dmp->lerr[var].below[prev] - dmp->preverr[var] > below[n])
            below[n] = dmp->lerr[var].below[prev] - dmp->preverr[var];
          if (dmp->lerr[var].above[prev] + kvijerr > above[n])
            above[n] = dmp->lerr[var].above[prev] + kvijerr;
        }
        dmp->preverr[var] = kvijerr;
        if (above[n] > err[var] || below[n] > err[var]) {
#ifdef DEBUG_GEOM
          printf("0: failed due to var cost\n");
#endif
          dmp->edgecost[n] = 2.0;
          return;
        }
        if (above[n] > below[n]) {
          if (above[n] / err[var] > dmp->edgecost[n])
            dmp->edgecost[n] = above[n] / err[var];
        } else {
          if (below[n] / err[var] > dmp->edgecost[n])
            dmp->edgecost[n] = below[n] / err[var];
        }
      }

      /* compute geometric error at this point */
      above = dmp->edgeerr[VTmeshNvars(dmp->mesh)].above;
      below = dmp->edgeerr[VTmeshNvars(dmp->mesh)].below;
      for (l = 0; l < 3; l++) {
        p1[l] = (1.0 - kval / (kval + vval)) *
                    dmp->tdata->verts[dmp->verts[k]].pt[l] +
                (kval / (kval + vval)) * dmp->tdata->verts[v].pt[l];
        p2[l] = (1.0 - ival / (jval + ival)) *
                    dmp->tdata->verts[dmp->verts[i]].pt[l] +
                (ival / (jval + ival)) * dmp->tdata->verts[dmp->verts[j]].pt[l];
      }
#ifdef DEBUG_LOOP
      printf("p1 = %lf %lf %lf\n", p1[0], p1[1], p1[2]);
      printf("p2 = %lf %lf %lf\n", p2[0], p2[1], p2[2]);
#endif
      kvijerr =
          sqrt(SQR(p1[0] - p2[0]) + SQR(p1[1] - p2[1]) + SQR(p1[2] - p2[2]));
#ifdef DEBUG_LOOP
      printf("GEOM kvijerr = %lf\n", kvijerr);
#endif
      if (dmp->edgenorm[k][0] * p2[0] + dmp->edgenorm[k][1] * p2[1] +
              dmp->edgenorm[k][2] * p2[2] - dmp->edgenorm[k][0] * p1[0] -
              dmp->edgenorm[k][1] * p1[1] - dmp->edgenorm[k][2] * p1[2] >
          0.0) {
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] + kvijerr > below[n])
          below[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] + kvijerr;
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] +
                dmp->preverr[VTmeshNvars(dmp->mesh)] >
            above[n])
          above[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] +
                     dmp->preverr[VTmeshNvars(dmp->mesh)];
      } else {
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] +
                dmp->preverr[VTmeshNvars(dmp->mesh)] >
            below[n])
          below[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] +
                     dmp->preverr[VTmeshNvars(dmp->mesh)];
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] + kvijerr > above[n])
          above[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] + kvijerr;
      }
      dmp->preverr[VTmeshNvars(dmp->mesh)] = kvijerr;
      if (above[n] > err[VTmeshNvars(dmp->mesh)] ||
          below[n] > err[VTmeshNvars(dmp->mesh)]) {
#ifdef DEBUG_GEOM
        printf("1: failed due to GEOM cost p1=(%lf %lf %lf) p2=(%lf %lf %lf)\n",
               p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
#endif
        dmp->edgecost[n] = 2.0;
        return;
      }
      if (above[n] > below[n]) {
        if (above[n] / err[VTmeshNvars(dmp->mesh)] > dmp->edgecost[n])
          dmp->edgecost[n] = above[n] / err[VTmeshNvars(dmp->mesh)];
      } else {
        if (below[n] / err[VTmeshNvars(dmp->mesh)] > dmp->edgecost[n])
          dmp->edgecost[n] = below[n] / err[VTmeshNvars(dmp->mesh)];
      }
    }
    for (var = 0; var < VTmeshNvars(dmp->mesh); var++) {
      above = dmp->edgeerr[var].above;
      below = dmp->edgeerr[var].below;
      if (dmp->preverr[var] < 0.0) {
        if (dmp->lerr[var].below[prev] - dmp->preverr[var] > below[n])
          below[n] = dmp->lerr[var].below[prev] - dmp->preverr[var];
        if (dmp->lerr[var].above[prev] + 0.0 > above[n])
          above[n] = dmp->lerr[var].above[prev] + 0.0;
      } else {
        if (dmp->lerr[var].below[prev] - 0.0 > below[n])
          below[n] = dmp->lerr[var].below[prev] - 0.0;
        if (dmp->lerr[var].above[prev] + dmp->preverr[var] > above[n])
          above[n] = dmp->lerr[var].above[prev] + dmp->preverr[var];
      }
      if (above[n] > err[var] || below[n] > err[var]) {
#ifdef DEBUG_GEOM
        printf("2: failed due to var cost\n");
#endif
        dmp->edgecost[n] = 2.0;
        return;
      }
      if (above[n] > below[n]) {
        if (above[n] / err[var] > dmp->edgecost[n])
          dmp->edgecost[n] = above[n] / err[var];
      } else {
        if (below[n] / err[var] > dmp->edgecost[n])
          dmp->edgecost[n] = below[n] / err[var];
      }
    }

    /* compute geometric error at this point */
    return;
  } else {
    vval = -vval;
#ifdef DEBUG_LOOP
    printf("vval is negative\n");
#endif
    /* vertex lies to side of hole 2 */
    for (prev = j, k = (j + 1) % dmp->nv; k != i;
         prev = k, k = (k == dmp->nv - 1 ? 0 : k + 1)) {
#ifdef DEBUG_LOOP
      printf("testing error to vert %d\n", k);
#endif
      kval = fabs(SDecimMeshPtPlane(dmp, dmp->verts[k], eq));
      SDecimMeshGetPlaneEq(dmp, dmp->verts[k], v, dmp->edgenorm[k], rplane);
      ival = fabs(SDecimMeshPtPlane(dmp, dmp->verts[i], rplane));
      jval = fabs(SDecimMeshPtPlane(dmp, dmp->verts[j], rplane));
#ifdef DEBUG_LOOP
      printf("vval, kval, ival, jval = %lf, %lf, %lf, %lf\n", vval, kval, ival,
             jval);
#endif
      for (var = 0; var < VTmeshNvars(dmp->mesh); var++) {
        above = dmp->edgeerr[var].above;
        below = dmp->edgeerr[var].below;
        switch (VTmeshVarType(dmp->mesh, var)) {
        case VT_UCHAR:
          kvval =
              (1.0 - kval / (kval + vval)) *
                  dmp->tdata->data[var]->data.ucdata[dmp->verts[k]] +
              (kval / (kval + vval)) * dmp->tdata->data[var]->data.ucdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.ucdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.ucdata[dmp->verts[j]];
          break;
        case VT_SHORT:
          kvval = (1.0 - kval / (kval + vval)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[k]] +
                  (kval / (kval + vval)) * dmp->tdata->data[var]->data.sdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.sdata[dmp->verts[j]];
          break;
        case VT_FLOAT:
          kvval = (1.0 - kval / (kval + vval)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[k]] +
                  (kval / (kval + vval)) * dmp->tdata->data[var]->data.fdata[v];
          ijval = (1.0 - ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[i]] +
                  (ival / (jval + ival)) *
                      dmp->tdata->data[var]->data.fdata[dmp->verts[j]];
          break;
        }
        kvijerr = kvval - ijval;
#ifdef DEBUG_LOOP
        printf("   var %d: kvijerr = %lf\n", var, kvijerr);
#endif
        if (kvijerr < dmp->preverr[var]) {
          if (dmp->lerr[var].below[prev] - kvijerr > below[n])
            below[n] = dmp->lerr[var].below[prev] - kvijerr;
          if (dmp->lerr[var].above[prev] + dmp->preverr[var] > above[n])
            above[n] = dmp->lerr[var].above[prev] + dmp->preverr[var];
        } else {
          if (dmp->lerr[var].below[prev] - dmp->preverr[var] > below[n])
            below[n] = dmp->lerr[var].below[prev] - dmp->preverr[var];
          if (dmp->lerr[var].above[prev] + kvijerr > above[n])
            above[n] = dmp->lerr[var].above[prev] + kvijerr;
        }
        dmp->preverr[var] = kvijerr;
        if (above[n] > err[var] || below[n] > err[var]) {
          dmp->edgecost[n] = 2.0;
#ifdef DEBUG_GEOM
          printf("3: failed due to var cost\n");
#endif
          return;
        }
        if (above[n] > below[n]) {
          if (above[n] / err[var] > dmp->edgecost[n])
            dmp->edgecost[n] = above[n] / err[var];
        } else {
          if (below[n] / err[var] > dmp->edgecost[n])
            dmp->edgecost[n] = below[n] / err[var];
        }
      }
      /* compute geometric error at this point */
      above = dmp->edgeerr[VTmeshNvars(dmp->mesh)].above;
      below = dmp->edgeerr[VTmeshNvars(dmp->mesh)].below;
      for (l = 0; l < 3; l++) {
        p1[l] = (1.0 - kval / (kval + vval)) *
                    dmp->tdata->verts[dmp->verts[k]].pt[l] +
                (kval / (kval + vval)) * dmp->tdata->verts[v].pt[l];
        p2[l] = (1.0 - ival / (jval + ival)) *
                    dmp->tdata->verts[dmp->verts[i]].pt[l] +
                (ival / (jval + ival)) * dmp->tdata->verts[dmp->verts[j]].pt[l];
      }
      kvijerr =
          sqrt(SQR(p1[0] - p2[0]) + SQR(p1[1] - p2[1]) + SQR(p1[2] - p2[2]));
#ifdef DEBUG_LOOP
      printf("GEOM kvijerr = %lf\n", kvijerr);
#endif
      if (dmp->edgenorm[k][0] * p2[0] + dmp->edgenorm[k][1] * p2[1] +
              dmp->edgenorm[k][2] * p2[2] - dmp->edgenorm[k][0] * p1[0] -
              dmp->edgenorm[k][1] * p1[1] - dmp->edgenorm[k][2] * p1[2] >
          0.0) {
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] + kvijerr > below[n])
          below[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] + kvijerr;
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] +
                dmp->preverr[VTmeshNvars(dmp->mesh)] >
            above[n])
          above[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] +
                     dmp->preverr[VTmeshNvars(dmp->mesh)];
      } else {
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] +
                dmp->preverr[VTmeshNvars(dmp->mesh)] >
            below[n])
          below[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].below[prev] +
                     dmp->preverr[VTmeshNvars(dmp->mesh)];
        if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] + kvijerr > above[n])
          above[n] = dmp->lerr[VTmeshNvars(dmp->mesh)].above[prev] + kvijerr;
      }
      dmp->preverr[VTmeshNvars(dmp->mesh)] = kvijerr;
      if (above[n] > err[VTmeshNvars(dmp->mesh)] ||
          below[n] > err[VTmeshNvars(dmp->mesh)]) {
        dmp->edgecost[n] = 2.0;
#ifdef DEBUG_GEOM
        printf("4: failed due to GEOM cost p1=(%lf %lf %lf) p2=(%lf %lf %lf)\n",
               p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
#endif
        return;
      }
      if (above[n] > below[n]) {
        if (above[n] / err[VTmeshNvars(dmp->mesh)] > dmp->edgecost[n])
          dmp->edgecost[n] = above[n] / err[VTmeshNvars(dmp->mesh)];
      } else {
        if (below[n] / err[VTmeshNvars(dmp->mesh)] > dmp->edgecost[n])
          dmp->edgecost[n] = below[n] / err[VTmeshNvars(dmp->mesh)];
      }
    }
    for (var = 0; var < VTmeshNvars(dmp->mesh); var++) {
      above = dmp->edgeerr[var].above;
      below = dmp->edgeerr[var].below;
      if (dmp->preverr[var] < 0.0) {
        if (dmp->lerr[var].below[prev] - dmp->preverr[var] > below[n])
          below[n] = dmp->lerr[var].below[prev] - dmp->preverr[var];
        if (dmp->lerr[var].above[prev] + 0.0 > above[n])
          above[n] = dmp->lerr[var].above[prev] + 0.0;
      } else {
        if (dmp->lerr[var].below[prev] - 0.0 > below[n])
          below[n] = dmp->lerr[var].below[prev] - 0.0;
        if (dmp->lerr[var].above[prev] + dmp->preverr[var] > above[n])
          above[n] = dmp->lerr[var].above[prev] + dmp->preverr[var];
      }
      if (above[n] > err[var] || below[n] > err[var]) {
        dmp->edgecost[n] = 2.0;
#ifdef DEBUG_GEOM
        printf("5: failed due to var cost\n");
#endif
        return;
      }
      if (above[n] > below[n]) {
        if (above[n] / err[var] > dmp->edgecost[n])
          dmp->edgecost[n] = above[n] / err[var];
      } else {
        if (below[n] / err[var] > dmp->edgecost[n])
          dmp->edgecost[n] = below[n] / err[var];
      }
    }
    /* compute geometric error at this point */
    return;
  }
}

int SDecimMeshCompTriErr(SDecimMeshP dmp, int v, int tri, float *err) {
  int i, var, v1, v2, v3;
  double u1[3], u2[3], u3[3], cross[3];
  double area1, area2, area3, sum;
  double newval, newerr, newpt[3];
  double *p;

  v1 = dmp->verts[dmp->newf[tri][0]];
  v2 = dmp->verts[dmp->newf[tri][1]];
  v3 = dmp->verts[dmp->newf[tri][2]];
  for (i = 0; i < 3; i++) {
    u1[i] = dmp->tdata->verts[v1].pt[i] - dmp->tdata->verts[v].pt[i];
    u2[i] = dmp->tdata->verts[v2].pt[i] - dmp->tdata->verts[v].pt[i];
    u3[i] = dmp->tdata->verts[v3].pt[i] - dmp->tdata->verts[v].pt[i];
  }
  cross[0] = u1[1] * u2[2] - u1[2] * u2[1];
  cross[1] = u1[2] * u2[0] - u1[0] * u2[2];
  cross[2] = u1[0] * u2[1] - u1[1] * u2[0];
  area1 = sqrt(SQR(cross[0]) + SQR(cross[1]) + SQR(cross[2]));

  cross[0] = u2[1] * u3[2] - u2[2] * u3[1];
  cross[1] = u2[2] * u3[0] - u2[0] * u3[2];
  cross[2] = u2[0] * u3[1] - u2[1] * u3[0];
  area2 = sqrt(SQR(cross[0]) + SQR(cross[1]) + SQR(cross[2]));

  cross[0] = u3[1] * u1[2] - u3[2] * u1[1];
  cross[1] = u3[2] * u1[0] - u3[0] * u1[2];
  cross[2] = u3[0] * u1[1] - u3[1] * u1[0];
  area3 = sqrt(SQR(cross[0]) + SQR(cross[1]) + SQR(cross[2]));

  sum = area1 + area2 + area3;

  for (var = 0; var < VTmeshNvars(dmp->mesh); var++) {
    /* new face has at least the errors of original faces in loop */
    for (i = 0; i < dmp->nv; i++) {
      if (dmp->lerr[var].above[i] > dmp->ferr[var].above[tri])
        dmp->ferr[var].above[tri] = dmp->lerr[var].above[i];
      if (dmp->lerr[var].below[i] > dmp->ferr[var].below[tri])
        dmp->ferr[var].below[tri] = dmp->lerr[var].below[i];
    }
    switch (VTmeshVarType(dmp->mesh, var)) {
    case VT_UCHAR:
      newval = area1 / sum * dmp->tdata->data[var]->data.ucdata[v3] +
               area2 / sum * dmp->tdata->data[var]->data.ucdata[v1] +
               area3 / sum * dmp->tdata->data[var]->data.ucdata[v2];
      newerr = dmp->tdata->data[var]->data.ucdata[v] - newval;
      break;
    case VT_SHORT:
      newval = area1 / sum * dmp->tdata->data[var]->data.sdata[v3] +
               area2 / sum * dmp->tdata->data[var]->data.sdata[v1] +
               area3 / sum * dmp->tdata->data[var]->data.sdata[v2];
      newerr = dmp->tdata->data[var]->data.sdata[v] - newval;
      break;
    case VT_FLOAT:
      newval = area1 / sum * dmp->tdata->data[var]->data.fdata[v3] +
               area2 / sum * dmp->tdata->data[var]->data.fdata[v1] +
               area3 / sum * dmp->tdata->data[var]->data.fdata[v2];
      newerr = dmp->tdata->data[var]->data.fdata[v] - newval;
      break;
    }
#ifdef DEBUG_LOOP
    printf("VERT var %d err=%lf\n", var, newerr);
#endif
    if (newerr < 0.0) {
      dmp->ferr[var].above[tri] -= newerr;
      if (dmp->ferr[var].above[tri] > err[var]) {
#ifdef DEBUG_GEOM
        printf("6: failed due to var cost\n");
#endif
#ifdef DEBUG_LOOP
        printf("VAR %d exceeded\n", var);
#endif
        return (0);
      }
    } else {
      dmp->ferr[var].below[tri] += newerr;
      if (dmp->ferr[var].below[tri] > err[var]) {
#ifdef DEBUG_GEOM
        printf("7: failed due to var cost\n");
#endif
#ifdef DEBUG_LOOP
        printf("VAR %d exceeded\n", var);
#endif
        return (0);
      }
    }
  }

  /* compute the propagated geometric error to this vertex */
  for (i = 0; i < dmp->nv; i++) {
    if (dmp->lerr[VTmeshNvars(dmp->mesh)].above[i] >
        dmp->ferr[VTmeshNvars(dmp->mesh)].above[tri])
      dmp->ferr[VTmeshNvars(dmp->mesh)].above[tri] =
          dmp->lerr[VTmeshNvars(dmp->mesh)].above[i];
    if (dmp->lerr[VTmeshNvars(dmp->mesh)].below[i] >
        dmp->ferr[VTmeshNvars(dmp->mesh)].below[tri])
      dmp->ferr[VTmeshNvars(dmp->mesh)].below[tri] =
          dmp->lerr[VTmeshNvars(dmp->mesh)].below[i];
  }

  /* now evaluate the geometric error at the point */
  for (i = 0; i < 3; i++) {
    newpt[i] = area1 / sum * dmp->tdata->verts[v3].pt[i] +
               area2 / sum * dmp->tdata->verts[v1].pt[i] +
               area3 / sum * dmp->tdata->verts[v2].pt[i];
  }

  p = dmp->tdata->verts[v].pt;
  newerr =
      sqrt(SQR(newpt[0] - p[0]) + SQR(newpt[1] - p[1]) + SQR(newpt[2] - p[2]));
  if (dmp->loopnorm[0] * newpt[0] + dmp->loopnorm[1] * newpt[1] +
          dmp->loopnorm[2] * newpt[2] - dmp->loopnorm[0] * p[0] -
          dmp->loopnorm[1] * p[1] - dmp->loopnorm[2] * p[2] <
      0.0) {
    dmp->ferr[VTmeshNvars(dmp->mesh)].above[tri] += newerr;
    if (dmp->ferr[VTmeshNvars(dmp->mesh)].above[tri] >
        err[VTmeshNvars(dmp->mesh)]) {
#ifdef DEBUG_GEOM
      printf("8: failed due to GEOM cost p=(%lf %lf %lf) newpt=(%lf %lf %lf)\n",
             p[0], p[1], p[2], newpt[0], newpt[1], newpt[2]);
#endif
      return (0);
    }
  } else {
    dmp->ferr[VTmeshNvars(dmp->mesh)].below[tri] += newerr;
    if (dmp->ferr[VTmeshNvars(dmp->mesh)].below[tri] >
        err[VTmeshNvars(dmp->mesh)]) {
#ifdef DEBUG_GEOM
      printf("9: failed due to GEOM cost p=(%lf %lf %lf) newpt=(%lf %lf %lf)\n",
             p[0], p[1], p[2], newpt[0], newpt[1], newpt[2]);
#endif
      return (0);
    }
  }
#ifdef DEBUG_LOOP
  printf("VERT newerr = %lf\n", newerr);
#endif
  return (1);
}

double SDecimMeshMinTriAng(SDecimMeshP dmp, int v1, int v2, int v3) {
  int j;
  double u[3], v[3], dot1, dot2, dot3, len;

  for (j = 0; j < 3; j++) {
    u[j] = dmp->tdata->verts[v2].pt[j] - dmp->tdata->verts[v1].pt[j];
    v[j] = dmp->tdata->verts[v3].pt[j] - dmp->tdata->verts[v1].pt[j];
  }
  dot1 = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  len = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) *
        sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  dot1 /= len;

  for (j = 0; j < 3; j++) {
    u[j] = dmp->tdata->verts[v1].pt[j] - dmp->tdata->verts[v2].pt[j];
    v[j] = dmp->tdata->verts[v3].pt[j] - dmp->tdata->verts[v2].pt[j];
  }
  dot2 = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  len = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) *
        sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  dot2 /= len;

  for (j = 0; j < 3; j++) {
    u[j] = dmp->tdata->verts[v1].pt[j] - dmp->tdata->verts[v3].pt[j];
    v[j] = dmp->tdata->verts[v2].pt[j] - dmp->tdata->verts[v3].pt[j];
  }
  dot3 = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  len = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) *
        sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  dot3 /= len;

  if (dot1 > dot2 && dot1 > dot3)
    return (dot1);
  else if (dot2 > dot3)
    return (dot2);
  else
    return (dot3);
}

double SDecimDihedral(SDecimMeshP dmp, int v1, int v2, int v3, int v4) {
  int i;
  double u[3], v[3], norm1[3], norm2[3], len;

  /* compute normal for v1,v2,v3 */
  for (i = 0; i < 3; i++) {
    u[i] = dmp->tdata->verts[v2].pt[i] - dmp->tdata->verts[v1].pt[i];
    v[i] = dmp->tdata->verts[v3].pt[i] - dmp->tdata->verts[v1].pt[i];
  }
  norm1[0] = u[1] * v[2] - u[2] * v[1];
  norm1[1] = u[2] * v[0] - u[0] * v[2];
  norm1[2] = u[0] * v[1] - u[1] * v[0];
  len = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1] + norm1[2] * norm1[2]);
  norm1[0] /= len;
  norm1[1] /= len;
  norm1[2] /= len;

  /* compute normal for v1,v3,v4 */
  for (i = 0; i < 3; i++) {
    u[i] = dmp->tdata->verts[v3].pt[i] - dmp->tdata->verts[v1].pt[i];
    v[i] = dmp->tdata->verts[v4].pt[i] - dmp->tdata->verts[v1].pt[i];
  }
  norm2[0] = u[1] * v[2] - u[2] * v[1];
  norm2[1] = u[2] * v[0] - u[0] * v[2];
  norm2[2] = u[0] * v[1] - u[1] * v[0];
  len = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1] + norm2[2] * norm2[2]);
  norm2[0] /= len;
  norm2[1] /= len;
  norm2[2] /= len;

  return (norm1[0] * norm2[0] + norm1[1] * norm2[1] + norm1[2] * norm2[2]);
}

int SDecimMeshFillHole(SDecimMeshP dmp, int v, float *err, int boundangle,
                       float aspect, float dihedral) {
  int niter, nedges, edges[255];
  int i, j, k, e, f, nexttry, next, var;
  int changed;
  double splitplane[4];
  double minloopang, ang;
  static int startv = 0;
  int f1, f2, p1, p2, p3, p4;
  int newf1[3], newf2[3], newadj1[3], newadj2[3];
  int vtri;
  double trierr;
  double min1, min2, min3, min4, min11, min21;

  if (dmp->nv == 3) {
#ifdef DEBUG_LOOP
    printf("filling hole of 3 verts\n");
#endif
    /* compute propagated error from previous faces */
    for (var = 0; var < VTmeshNvars(dmp->mesh) + 1; var++) {
      dmp->ferr[var].below[0] =
          MAX3(dmp->lerr[var].below[0], dmp->lerr[var].below[1],
               dmp->lerr[var].below[2]);
      dmp->ferr[var].above[0] =
          MAX3(dmp->lerr[var].above[0], dmp->lerr[var].above[1],
               dmp->lerr[var].above[2]);
      if (dmp->ferr[var].below[0] > err[var] ||
          dmp->ferr[var].above[0] > err[var]) {
#ifdef DEBUG_SUMMARY
        NRJERR++;
#endif
        return (0);
      }
    }

    /* check the added error by projection */
    dmp->nf = 1;
    dmp->newf[0][0] = 0;
    dmp->newf[0][1] = 1;
    dmp->newf[0][2] = 2;
    if (!SDecimMeshCompTriErr(dmp, v, 0, err)) {
#ifdef DEBUG_SUMMARY
      NRJERR++;
#endif
      return (0);
    }

    /* get index for new face */
    SDecimMeshGetFaces(dmp, dmp->newfaces, 1);

    /* face is ok to add... */
    for (i = 0; i < dmp->nv; i++) {
      f = dmp->adj[i];
      if (f != -1) {
        if (dmp->tdata->tri[f].adjtri[0] == dmp->face[i])
          dmp->tdata->tri[f].adjtri[0] = dmp->newfaces[0];
        else if (dmp->tdata->tri[f].adjtri[1] == dmp->face[i])
          dmp->tdata->tri[f].adjtri[1] = dmp->newfaces[0];
        else if (dmp->tdata->tri[f].adjtri[2] == dmp->face[i])
          dmp->tdata->tri[f].adjtri[2] = dmp->newfaces[0];
      }
    }
    for (var = 0; var < VTmeshNvars(dmp->mesh) + 1; var++) {
      dmp->looperr.below[var] = dmp->ferr[var].below[0];
      dmp->looperr.above[var] = dmp->ferr[var].above[0];
    }
    SDecimMeshSetFace(dmp, dmp->newfaces[0], dmp->verts, dmp->adj,
                      dmp->looperr.below, dmp->looperr.above);
    return (1);
  }

  /* to start triangulation, we have no edges */
  dmp->nedges = 0;

  /* set up the indexing for shrinking the loop */
  dmp->vprev[0] = dmp->nv - 1;
  dmp->vnext[0] = 1;
  for (i = 1; i < dmp->nv - 1; i++) {
    dmp->vprev[i] = i - 1;
    dmp->vnext[i] = i + 1;
  }
  dmp->vprev[dmp->nv - 1] = dmp->nv - 2;
  dmp->vnext[i] = 0;

  if ((nexttry = startv++) >= dmp->nv)
    nexttry = startv = 0;

#ifdef DEBUG_LOOP
  printf("Loop: size %d\n  ", dmp->nv);
  for (i = 0; i < dmp->nv; i++)
    printf("   %d (%lf, %lf, %lf %lf)\n", dmp->verts[i],
           dmp->tdata->verts[dmp->verts[i]].pt[0],
           dmp->tdata->verts[dmp->verts[i]].pt[1],
           dmp->tdata->verts[dmp->verts[i]].pt[2],
           dmp->tdata->data[0]->data.fdata[dmp->verts[i]]);
  printf("\n");
  for (i = 0; i < dmp->nv; i++)
    printf(" %d", dmp->adj[i]);
  printf("\n");
  for (i = 0; i < dmp->nv; i++)
    printf(" %d", dmp->face[i]);
  printf("\n");
  for (i = 0; i < dmp->nv; i++)
    printf(" %lf %lf", dmp->lerr[0].below[i], dmp->lerr[0].above[i]);
  printf("\n");
#endif

  dmp->nf = 0;

  /* check for an edge boundary in loop */
#ifdef DIHEDRAL_CHECK
  nedges = 0;
  for (i = 0; i < dmp->nv; i++) {
    /* check dihedral angle between v and verts[i] */
    if (SDecimDihedral(dmp, dmp->verts[i], dmp->verts[(i + 1) % dmp->nv], v,
                       dmp->verts[i == 0 ? dmp->nv - 1 : i - 1]) < dihedral) {
      edges[nedges] = i;
      nedges++;
    }
  }

  if (nedges > 2) {
    printf("%d edges in loop, rejected\n", nedges);
    printf("vert: %lf %lf %lf\n", dmp->tdata->verts[v].pt[0],
           dmp->tdata->verts[v].pt[1], dmp->tdata->verts[v].pt[2]);
    for (i = 0; i < dmp->nv; i++)
      printf("   %d (%lf, %lf, %lf %lf)\n", dmp->verts[i],
             dmp->tdata->verts[dmp->verts[i]].pt[0],
             dmp->tdata->verts[dmp->verts[i]].pt[1],
             dmp->tdata->verts[dmp->verts[i]].pt[2],
             dmp->tdata->data[0]->data.fdata[dmp->verts[i]]);
#ifdef DEBUG_SUMMARY
    NRJFTR++;
#endif
    return (0);
  }
#endif

  /* perform an initial triangulation */
  vtri = -1;
  for (dmp->nedges = 0; dmp->nedges < dmp->nv - 3; dmp->nedges++) {
    /* choose another edge to chop off */
    for (i = 0; i < dmp->nv - dmp->nedges; i++) {
      next = dmp->vnext[nexttry];
#ifdef DEBUG_LOOP
      printf("testing vertex %d (%d, %d)\n", nexttry, dmp->vprev[nexttry],
             dmp->vnext[nexttry]);
#endif
#ifdef DIHEDRAL_CHECK
      if (nedges != 2 || (((dmp->vprev[nexttry] >= edges[1] ||
                            dmp->vprev[nexttry] <= edges[0]) &&
                           (dmp->vnext[nexttry] >= edges[1] ||
                            dmp->vnext[nexttry] <= edges[0])) ||
                          ((dmp->vprev[nexttry] >= edges[0] &&
                            dmp->vprev[nexttry] <= edges[1]) &&
                           (dmp->vnext[nexttry] >= edges[0] &&
                            dmp->vnext[nexttry] <= edges[1])))) {
#endif
        SDecimMeshGetPlaneEq(dmp, dmp->verts[dmp->vprev[nexttry]],
                             dmp->verts[dmp->vnext[nexttry]], dmp->loopnorm,
                             splitplane);
        if (SDecimMeshValidPlane(dmp, dmp->vprev[nexttry], dmp->vnext[nexttry],
                                 splitplane)) {
#ifdef DEBUG_LOOP
          printf("chopping vertex %d\n", nexttry);
#endif
          /* save information for face created */
          dmp->newf[dmp->nf][0] = dmp->vprev[nexttry];
          dmp->newf[dmp->nf][1] = nexttry;
          dmp->newf[dmp->nf][2] = dmp->vnext[nexttry];
          dmp->newadj[dmp->nf][0] = dmp->adj[dmp->vprev[nexttry]];
          dmp->newadj[dmp->nf][1] = dmp->adj[nexttry];
          dmp->newadj[dmp->nf][2] = -2; /* don't know adjacency yet */

          dmp->adj[dmp->vprev[nexttry]] = -dmp->nf - 3;

          /* save information for edge created */
          dmp->edgeadj[dmp->nedges][0] = dmp->nf;
          dmp->edgeadj[dmp->nedges][1] = -1;
          dmp->edgeij[dmp->nedges][0] = dmp->vprev[nexttry];
          dmp->edgeij[dmp->nedges][1] = dmp->vnext[nexttry];

          dmp->vprev[dmp->vnext[nexttry]] = dmp->vprev[nexttry];
          dmp->vnext[dmp->vprev[nexttry]] = dmp->vnext[nexttry];
          dmp->vprev[nexttry] = dmp->vnext[nexttry] = -1;

          nexttry = next;

          /* check to see if this triangle contains the deleted point */
          if (vtri == -1 && SDecimMeshPtPlane(dmp, v, splitplane) < 0.0)
            vtri = dmp->nf;

          dmp->nf++;

          break;
        }
#ifdef DIHEDRAL_CHECK
      }
#endif

      nexttry = next;
    }

    /* did we lop off another vertex? */
    if (i == dmp->nv - dmp->nedges) {
#ifdef DIHEDRAL_CHECK
#ifdef DEBUG_TRI
      if (nedges == 2)
        printf("no initial triangulation for edge hole\n");
      printf("can't find first triangulation (startv = %d)!\n", startv);
      for (i = 0; i < dmp->nv; i++) {
        printf("   %lf %lf %lf\n", dmp->tdata->verts[dmp->verts[i]].pt[0],
               dmp->tdata->verts[dmp->verts[i]].pt[1],
               dmp->tdata->verts[dmp->verts[i]].pt[2]);
      }
      for (i = 0; i < dmp->nv - dmp->nedges; i++) {
        printf(" %d ", nexttry);
        nexttry = dmp->vnext[nexttry];
      }
      printf("\n");
#endif
#endif
#ifdef DEBUG_TRI
      printf("can't find first triangulation (startv = %d)!\n", startv);
      for (i = 0; i < dmp->nv; i++) {
        printf("   %lf %lf %lf\n", dmp->tdata->verts[dmp->verts[i]].pt[0],
               dmp->tdata->verts[dmp->verts[i]].pt[1],
               dmp->tdata->verts[dmp->verts[i]].pt[2]);
      }
      for (i = 0; i < dmp->nv - dmp->nedges; i++) {
        printf(" %d ", nexttry);
        nexttry = dmp->vnext[nexttry];
      }
#endif
#ifdef DEBUG_SUMMARY
#ifdef DIHEDRAL_CHECK
      if (nedges == 2)
        NRJTRI2++;
      else
        NRJTRI++;
#else
      NRJTRI++;
#endif
#endif
      return (0);
    }
  }

  /* The last split forms one more triangle.. */
  dmp->newf[dmp->nf][0] = dmp->newf[dmp->nf - 1][0];
  dmp->newf[dmp->nf][1] = dmp->newf[dmp->nf - 1][2];
  dmp->newf[dmp->nf][2] = dmp->vnext[dmp->newf[dmp->nf - 1][2]];
  dmp->newadj[dmp->nf][0] = dmp->adj[dmp->newf[dmp->nf][0]];
  dmp->newadj[dmp->nf][1] = dmp->adj[dmp->newf[dmp->nf][1]];
  dmp->newadj[dmp->nf][2] = dmp->adj[dmp->newf[dmp->nf][2]];

  /* if we hadn't cut off the original vertex yet, it must be
     contained in this final triangle */
  if (vtri == -1)
    vtri = dmp->nf;

  dmp->nf++;

#ifdef DEBUG_TRI
  printf("deleted vert is in tri %d\n", vtri);
#endif

#ifdef DEBUG_LOOP
  printf("faces formed:\n");
  for (i = 0; i < dmp->nf; i++)
    printf("(%d %d %d) (%d %d %d)\n", dmp->newf[i][0], dmp->newf[i][1],
           dmp->newf[i][2], dmp->newadj[i][0], dmp->newadj[i][1],
           dmp->newadj[i][2]);

  printf("edges formed:\n");
  for (i = 0; i < dmp->nedges; i++)
    printf("(%d %d)\n", dmp->edgeadj[i][0], dmp->edgeadj[i][1]);
#endif

  for (i = 0; i < dmp->nf; i++) {
    for (j = 0; j < 3; j++) {
      if (dmp->newadj[i][j] <= -3) {
        f = -dmp->newadj[i][j] - 3;
        if (dmp->newadj[f][2] == -2)
          dmp->newadj[f][2] = -i - 3;
        for (e = 0; e < dmp->nedges; e++) {
          if (dmp->edgeadj[e][0] == f)
            dmp->edgeadj[e][1] = i;
        }
      }
    }
  }

#ifdef DEBUG_LOOP
  printf("faces formed:\n");
  for (i = 0; i < dmp->nf; i++)
    printf("(%d %d %d) (%d %d %d)\n", dmp->newf[i][0], dmp->newf[i][1],
           dmp->newf[i][2], dmp->newadj[i][0], dmp->newadj[i][1],
           dmp->newadj[i][2]);

  printf("edges formed:\n");
  for (i = 0; i < dmp->nedges; i++)
    printf("(%d %d)\n", dmp->edgeadj[i][0], dmp->edgeadj[i][1]);
#endif

    /*
     * compute the costs for the edges in initial triangulation
     */
#ifdef DEBUG_LOOP
  printf("edge costs:\n");
#endif
  for (i = 0; i < dmp->nedges; i++) {
    SDecimMeshGetPlaneEq(dmp, dmp->verts[dmp->edgeij[i][0]],
                         dmp->verts[dmp->edgeij[i][1]], dmp->loopnorm,
                         splitplane);

    SDecimMeshCompHoleErr(dmp, v, dmp->edgeij[i][0], dmp->edgeij[i][1], i,
                          splitplane, err);
#ifdef DEBUG_LOOP
    printf("%lf (%lf %lf)", dmp->edgecost[i], dmp->edgeerr[0].below[i],
           dmp->edgeerr[0].above[i]);
#endif
  }
#ifdef DEBUG_LOOP
  printf("END\n");
  printf("\n");
#endif

  minloopang = -1e100;
  for (i = 0; i < dmp->nv; i++) {
    ang = SDecimMeshMinTriAng(dmp, dmp->verts[i], v,
                              dmp->verts[(i + 1) % dmp->nv]);
    if (ang > minloopang)
      minloopang = ang;
  }

#ifdef SWAP
  /*
   * flip the edges to try to lower the error
   */
  niter = 0;
  do {
    changed = 0;
    for (i = 0; i < dmp->nedges; i++) {
#ifdef DIHEDRAL_CHECK
      if (nedges == 2 &&
          ((dmp->edgeij[i][0] == edges[0] && dmp->edgeij[i][1] == edges[1]) ||
           (dmp->edgeij[i][0] == edges[1] && dmp->edgeij[i][1] == edges[0])))
        continue;
#endif
      f1 = dmp->edgeadj[i][0];
      f2 = dmp->edgeadj[i][1];
      for (p1 = 0; p1 < 3; p1++)
        if (dmp->newf[f1][p1] == dmp->edgeij[i][1])
          break;
      if (p1 == 3) {
        fprintf(stderr, "Ack! EdgeVert1 not found on face!\n");
        sleep(0);
        return (0);
      }
      for (p2 = 0; p2 < 3; p2++)
        if (dmp->newf[f2][p2] == dmp->edgeij[i][0])
          break;
      if (p2 == 3) {
        fprintf(stderr, "Ack! EdgeVert2 not found on face!\n");
        sleep(0);
        return (0);
      }

      /* found the of the faces which meet.. */
      p3 = (p1 == 0 ? 2 : p1 - 1);
      p4 = (p2 == 0 ? 2 : p2 - 1);

      /* get the error for the flipped edge */
      SDecimMeshGetPlaneEq(dmp, dmp->verts[dmp->newf[f1][p3]],
                           dmp->verts[dmp->newf[f2][p4]], dmp->loopnorm,
                           splitplane);
      /* is the flip valid? */
      if (SDecimMeshPtPlane(dmp, dmp->verts[dmp->newf[f1][p1]], splitplane) *
              SDecimMeshPtPlane(dmp, dmp->verts[dmp->newf[f2][p2]],
                                splitplane) >=
          0.0)
        continue;
      SDecimMeshCompHoleErr(dmp, v, dmp->newf[f1][p3], dmp->newf[f2][p4],
                            dmp->nedges, splitplane, err);

#ifdef SWAP_BY_ANGLES
      min1 = SDecimMeshMinTriAng(dmp, dmp->verts[dmp->newf[f1][0]],
                                 dmp->verts[dmp->newf[f1][1]],
                                 dmp->verts[dmp->newf[f1][2]]);
      min2 = SDecimMeshMinTriAng(dmp, dmp->verts[dmp->newf[f2][0]],
                                 dmp->verts[dmp->newf[f2][1]],
                                 dmp->verts[dmp->newf[f2][2]]);
      min11 = (min1 > min2) ? min1 : min2;
      min3 = SDecimMeshMinTriAng(dmp, dmp->verts[dmp->newf[f1][p1]],
                                 dmp->verts[dmp->newf[f1][p3]],
                                 dmp->verts[dmp->newf[f2][p4]]);
      min4 = SDecimMeshMinTriAng(dmp, dmp->verts[dmp->newf[f2][p2]],
                                 dmp->verts[dmp->newf[f2][p4]],
                                 dmp->verts[dmp->newf[f1][p3]]);
      min21 = (min3 > min4) ? min3 : min4;
      /* use angle plus cost as measure */
      if ((!boundangle && dmp->edgecost[dmp->nedges] < dmp->edgecost[i]) ||
          (min21 + dmp->edgecost[dmp->nedges] < min11 + dmp->edgecost[i])) {
/* cost must improve, angles must stay bounded */
/*
         if ((!boundangle || ((min21 < min11) || min21 < aspect)) &&
               dmp->edgecost[dmp->nedges] < dmp->edgecost[i]) {
*/
/* angles must improve, cost must stay bounded */
/*
         if ((!boundangle || (min21 < min11)) &&
              ((dmp->edgecost[dmp->nedges] < dmp->edgecost[i]) ||
                dmp->edgecost[dmp->nedges] < 1.0)) {
*/
#else
      if (dmp->edgecost[dmp->nedges] < dmp->edgecost[i]) {
#endif
        /* perform a swap */
        changed = 1;

#ifdef DEBUG_SWAP
        printf("flipping edge %d (%d,%d) (%lf->%lf)\n", i, dmp->edgeij[i][0],
               dmp->edgeij[i][1], dmp->edgecost[i], dmp->edgecost[dmp->nedges]);
        printf("p1, p2, p3, p4 = %d %d %d %d\n", dmp->newf[f1][p1],
               dmp->newf[f2][p2], dmp->newf[f1][p3], dmp->newf[f2][p4]);
#endif

        /* first update the edgeadj for other edges */
        /* two adjacencies will change */
        for (e = 0; e < dmp->nedges; e++) {
          if (e != i) {
            if (dmp->edgeij[e][0] == dmp->newf[f1][p1] &&
                dmp->edgeij[e][1] == dmp->newf[f1][p3]) {
              if (dmp->edgeadj[e][0] != f1) {
                fprintf(stderr, "Hmm.. adjacency not right\n");
                sleep(0);
                return (0);
              }
              /* update the edge adjacency and face adjacency for
                 adjacent face */
              dmp->edgeadj[e][0] = f2;
              f = dmp->edgeadj[e][1];
              for (j = 0; j < 3; j++) {
                if (dmp->newadj[f][j] == -f1 - 3) {
                  dmp->newadj[f][j] = -f2 - 3;
                  break;
                }
              }
              if (j == 3) {
                fprintf(stderr, "Hmm.. fadjacency not right\n");
                sleep(0);
                return (0);
              }
            } else if (dmp->edgeij[e][1] == dmp->newf[f1][p1] &&
                       dmp->edgeij[e][0] == dmp->newf[f1][p3]) {
              if (dmp->edgeadj[e][1] != f1) {
                fprintf(stderr, "Hmm.. adjacency not right\n");
                sleep(0);
                return (0);
              }
              dmp->edgeadj[e][1] = f2;
              f = dmp->edgeadj[e][0];
              for (j = 0; j < 3; j++) {
                if (dmp->newadj[f][j] == -f1 - 3) {
                  dmp->newadj[f][j] = -f2 - 3;
                  break;
                }
              }
              if (j == 3) {
                fprintf(stderr, "Hmm.. fadjacency not right\n");
                sleep(0);
                return (0);
              }
            }
            if (dmp->edgeij[e][0] == dmp->newf[f2][p2] &&
                dmp->edgeij[e][1] == dmp->newf[f2][p4]) {
              if (dmp->edgeadj[e][0] != f2) {
                fprintf(stderr, "Hmm.. adjacency not right\n");
                sleep(0);
                return (0);
              }
              /* update the edge adjacency and face adjacency for
                 adjacent face */
              dmp->edgeadj[e][0] = f1;
              f = dmp->edgeadj[e][1];
              for (j = 0; j < 3; j++) {
                if (dmp->newadj[f][j] == -f2 - 3) {
                  dmp->newadj[f][j] = -f1 - 3;
                  break;
                }
              }
              if (j == 3) {
                fprintf(stderr, "Hmm.. fadjacency not right\n");
                sleep(0);
                return (0);
              }
            } else if (dmp->edgeij[e][1] == dmp->newf[f2][p2] &&
                       dmp->edgeij[e][0] == dmp->newf[f2][p4]) {
              if (dmp->edgeadj[e][1] != f2) {
                fprintf(stderr, "Hmm.. adjacency not right\n");
                sleep(0);
                return (0);
              }
              dmp->edgeadj[e][1] = f1;
              f = dmp->edgeadj[e][0];
              for (j = 0; j < 3; j++) {
                if (dmp->newadj[f][j] == -f2 - 3) {
                  dmp->newadj[f][j] = -f1 - 3;
                  break;
                }
              }
              if (j == 3) {
                fprintf(stderr, "Hmm.. fadjacency not right\n");
                sleep(0);
                return (0);
              }
            }
          }
        }

        /* now change the internal parts */
        /* edgeadj[i][0,1] stay the same */
        dmp->edgeij[i][0] = dmp->newf[f2][p4];
        dmp->edgeij[i][1] = dmp->newf[f1][p3];
#ifdef DEBUG_SWAP
        printf("new edgeij: %d %d\n", dmp->edgeij[i][0], dmp->edgeij[i][1]);
#endif
        for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
          dmp->edgeerr[j].below[i] = dmp->edgeerr[j].below[dmp->nedges];
          dmp->edgeerr[j].above[i] = dmp->edgeerr[j].above[dmp->nedges];
        }
        dmp->edgecost[i] = dmp->edgecost[dmp->nedges];

        /* and change the face info */
        newf1[0] = dmp->newf[f1][p3];
        newf1[1] = dmp->newf[f2][p4];
        newf1[2] = dmp->newf[f2][p2];
        newadj1[0] = -f2 - 3;
        newadj1[1] = dmp->newadj[f2][p4];
        newadj1[2] = dmp->newadj[f1][p3 == 0 ? 2 : p3 - 1];

        newf2[0] = dmp->newf[f1][p3];
        newf2[1] = dmp->newf[f1][p1];
        newf2[2] = dmp->newf[f2][p4];
        newadj2[0] = dmp->newadj[f1][p3];
        newadj2[1] = dmp->newadj[f2][p4 == 0 ? 2 : p4 - 1];
        newadj2[2] = -f1 - 3;

        memcpy(dmp->newf[f1], newf1, sizeof(int[3]));
        memcpy(dmp->newadj[f1], newadj1, sizeof(int[3]));
        memcpy(dmp->newf[f2], newf2, sizeof(int[3]));
        memcpy(dmp->newadj[f2], newadj2, sizeof(int[3]));

#ifdef DEBUG_SWAP
        printf("new t1: %d %d %d\n", dmp->newf[f1][0], dmp->newf[f1][1],
               dmp->newf[f1][2]);
        printf("new t2: %d %d %d\n", dmp->newf[f2][0], dmp->newf[f2][1],
               dmp->newf[f2][2]);
#endif
        /*
         * if deleted vertex was in one of these triangles, we must
         * determine which new triangle contains it
         */
        if (vtri == f1 || vtri == f2) {
          /* check to see which triangle contains the deleted point */
          if (SDecimMeshPtPlane(dmp, v, splitplane) < 0.0)
            vtri = f2;
          else
            vtri = f1;
        }
      }
    }
    if (++niter > 50) {
      printf("swapping endlessly!\n");
      ;
      break;
    }
  } while (changed);
#ifdef DEBUG_SWAP
  printf("***DONE FLIPPING\n");
#endif
#endif

  /* do all edges meet the error criteria? */
  for (i = 0; i < dmp->nedges; i++) {
    if (dmp->edgecost[i] > 1.0) {
#ifdef DEBUG_SUMMARY
      NRJERR++;
#endif
      return (0);
    }
  }

  /* zero errors on faces */
  for (i = 0; i < VTmeshNvars(dmp->mesh) + 1; i++) {
    memset(dmp->ferr[i].below, '\0', sizeof(float) * dmp->nf);
    memset(dmp->ferr[i].above, '\0', sizeof(float) * dmp->nf);
  }

  /*
   * compute the cost for the new triangle containing the deleted
   * vertex
   */
  if (!SDecimMeshCompTriErr(dmp, v, vtri, err)) {
#ifdef DEBUG_SUMMARY
    NRJTERR++;
#endif
    return (0);
  }

  /* get errors from interior edges */
  for (i = 0; i < dmp->nedges; i++) {
    f1 = dmp->edgeadj[i][0];
    f2 = dmp->edgeadj[i][1];
    for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
      if (dmp->edgeerr[j].below[i] > dmp->ferr[j].below[f1])
        dmp->ferr[j].below[f1] = dmp->edgeerr[j].below[i];
      if (dmp->edgeerr[j].above[i] > dmp->ferr[j].above[f1])
        dmp->ferr[j].above[f1] = dmp->edgeerr[j].above[i];

      if (dmp->edgeerr[j].below[i] > dmp->ferr[j].below[f2])
        dmp->ferr[j].below[f2] = dmp->edgeerr[j].below[i];
      if (dmp->edgeerr[j].above[i] > dmp->ferr[j].above[f2])
        dmp->ferr[j].above[f2] = dmp->edgeerr[j].above[i];
    }
  }

  /* and error from outer edges */
  for (i = 0; i < dmp->nf; i++) {
    if (dmp->newf[i][1] == (dmp->newf[i][0] + 1) % dmp->nv) {
      for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
        if (dmp->lerr[j].below[dmp->newf[i][0]] > dmp->ferr[j].below[i])
          dmp->ferr[j].below[i] = dmp->lerr[j].below[dmp->newf[i][0]];
        if (dmp->lerr[j].above[dmp->newf[i][0]] > dmp->ferr[j].above[i])
          dmp->ferr[j].above[i] = dmp->lerr[j].above[dmp->newf[i][0]];
      }
    }
    if (dmp->newf[i][2] == (dmp->newf[i][1] + 1) % dmp->nv) {
      for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
        if (dmp->lerr[j].below[dmp->newf[i][1]] > dmp->ferr[j].below[i])
          dmp->ferr[j].below[i] = dmp->lerr[j].below[dmp->newf[i][1]];
        if (dmp->lerr[j].above[dmp->newf[i][1]] > dmp->ferr[j].above[i])
          dmp->ferr[j].above[i] = dmp->lerr[j].above[dmp->newf[i][1]];
      }
    }
    if (dmp->newf[i][0] == (dmp->newf[i][2] + 1) % dmp->nv) {
      for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
        if (dmp->lerr[j].below[dmp->newf[i][2]] > dmp->ferr[j].below[i])
          dmp->ferr[j].below[i] = dmp->lerr[j].below[dmp->newf[i][2]];
        if (dmp->lerr[j].above[dmp->newf[i][2]] > dmp->ferr[j].above[i])
          dmp->ferr[j].above[i] = dmp->lerr[j].above[dmp->newf[i][2]];
      }
    }
  }

  for (i = 0; i < dmp->nf; i++) {
    dmp->newf[i][0] = dmp->verts[dmp->newf[i][0]];
    dmp->newf[i][1] = dmp->verts[dmp->newf[i][1]];
    dmp->newf[i][2] = dmp->verts[dmp->newf[i][2]];
  }

#ifdef MINANG_CRITERIA
  if (boundangle)
    for (i = 0; i < dmp->nf; i++) {
      double v1[3], v2[3], dot, len;
      for (j = 0; j < 3; j++) {
        v1[j] = dmp->tdata->verts[dmp->newf[i][1]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][0]].pt[j];
        v2[j] = dmp->tdata->verts[dmp->newf[i][2]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][0]].pt[j];
      }
      dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      len = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) *
            sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
      dot /= len;
      if (fabs(dot) > 1.0)
        printf("dot is %lf\n", dot);
      if (dot > aspect && dot > minloopang) {
#ifdef DEBUG_REJECT
        printf("rejected %d %lf\n", ++NRJ, dot);
#endif
#ifdef DEBUG_SUMMARY
        NRJFTR++;
#endif
        return (0);
      }

      for (j = 0; j < 3; j++) {
        v1[j] = dmp->tdata->verts[dmp->newf[i][0]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][1]].pt[j];
        v2[j] = dmp->tdata->verts[dmp->newf[i][2]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][1]].pt[j];
      }
      dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      len = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) *
            sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
      dot /= len;
      if (fabs(dot) > 1.0)
        printf("dot is %lf\n", dot);
      if (dot > aspect && dot > minloopang) {
#ifdef DEBUG_REJECT
        printf("rejected %d %lf\n", ++NRJ, dot);
#endif
#ifdef DEBUG_SUMMARY
        NRJFTR++;
#endif
        return (0);
      }

      for (j = 0; j < 3; j++) {
        v1[j] = dmp->tdata->verts[dmp->newf[i][0]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][2]].pt[j];
        v2[j] = dmp->tdata->verts[dmp->newf[i][1]].pt[j] -
                dmp->tdata->verts[dmp->newf[i][2]].pt[j];
      }
      dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      len = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]) *
            sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
      dot /= len;
      if (fabs(dot) > 1.0)
        printf("dot is %lf\n", dot);
      if (dot > aspect && dot > minloopang) {
#ifdef DEBUG_REJECT
        printf("rejected %d %lf\n", ++NRJ, dot);
#endif
#ifdef DEBUG_SUMMARY
        NRJFTR++;
#endif
        return (0);
      }
    }
#endif

  /* reserve indices for new faces */
  SDecimMeshGetFaces(dmp, dmp->newfaces, dmp->nv - 2);

  /* now change the vertices and adjacencies to actually be correct */
  for (i = 0; i < dmp->nf; i++) {
    if (dmp->newadj[i][0] <= -3)
      dmp->newadj[i][0] = dmp->newfaces[-dmp->newadj[i][0] - 3];
    else if (dmp->newadj[i][0] >= 0) {
      f = dmp->newadj[i][0];
      if (dmp->tdata->tri[f].vert[0] == dmp->newf[i][1])
        dmp->tdata->tri[f].adjtri[0] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[1] == dmp->newf[i][1])
        dmp->tdata->tri[f].adjtri[1] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[2] == dmp->newf[i][1])
        dmp->tdata->tri[f].adjtri[2] = dmp->newfaces[i];
    }

    if (dmp->newadj[i][1] <= -3)
      dmp->newadj[i][1] = dmp->newfaces[-dmp->newadj[i][1] - 3];
    else if (dmp->newadj[i][1] >= 0) {
      f = dmp->newadj[i][1];
      if (dmp->tdata->tri[f].vert[0] == dmp->newf[i][2])
        dmp->tdata->tri[f].adjtri[0] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[1] == dmp->newf[i][2])
        dmp->tdata->tri[f].adjtri[1] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[2] == dmp->newf[i][2])
        dmp->tdata->tri[f].adjtri[2] = dmp->newfaces[i];
    }

    if (dmp->newadj[i][2] <= -3)
      dmp->newadj[i][2] = dmp->newfaces[-dmp->newadj[i][2] - 3];
    else if (dmp->newadj[i][2] >= 0) {
      f = dmp->newadj[i][2];
      if (dmp->tdata->tri[f].vert[0] == dmp->newf[i][0])
        dmp->tdata->tri[f].adjtri[0] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[1] == dmp->newf[i][0])
        dmp->tdata->tri[f].adjtri[1] = dmp->newfaces[i];
      else if (dmp->tdata->tri[f].vert[2] == dmp->newf[i][0])
        dmp->tdata->tri[f].adjtri[2] = dmp->newfaces[i];
    }
  }

#ifdef DEBUG_LOOP
  printf("faces formed:\n");
  for (i = 0; i < dmp->nf; i++)
    printf("(%d %d %d) (%d %d %d) (%lf %lf)\n", dmp->newf[i][0],
           dmp->newf[i][1], dmp->newf[i][2], dmp->newadj[i][0],
           dmp->newadj[i][1], dmp->newadj[i][2], dmp->ferr[0].below[i],
           dmp->ferr[0].above[i]);
#endif

  for (i = 0; i < dmp->nf; i++) {
    for (j = 0; j < VTmeshNvars(dmp->mesh) + 1; j++) {
      dmp->looperr.below[j] = dmp->ferr[j].below[i];
      dmp->looperr.above[j] = dmp->ferr[j].above[i];
    }
    SDecimMeshSetFace(dmp, dmp->newfaces[i], dmp->newf[i], dmp->newadj[i],
                      dmp->looperr.below, dmp->looperr.above);
  }

  return (1);
}

int SDecimMeshRemoveSimpleVert(SDecimMeshP dmp, int v, float *err,
                               int boundangle, float aspect, float dihedral) {
  int i;
  static int nok = 0;

  /*
  if (nok > 0)
     return(0);
  */

#ifdef DEBUG_SUMMARY
  if (NATT % 100 == 0) {
    printf(
        "att %d del, perf %d (%d TRI) (%d TRI2) (%d ERR) (%d TERR) (%d FTR)\n",
        NATT, NDEL, NRJTRI, NRJTRI2, NRJERR, NRJTERR, NRJFTR);
  }
  NATT++;
#endif

  if (dmp->nv < 3) {
    fprintf(stderr, "can't remove simple vert w/ %d adj\n", dmp->nv);
    return (0);
  }

  /* see if the hole can be filled */
  if (!SDecimMeshFillHole(dmp, v, err, boundangle, aspect, dihedral)) {
    return (0);
  }
#ifdef DEBUG_SUMMARY
  NDEL++;
#endif
  /*
  sleep(0);
  if ((++nok) % 1000 == 0)
     printf("%d verts deleted\n", nok);
  */

  /* free the vertex */
  SDecimMeshFreeVert(dmp, v);

  /* free the old triangles */
  for (i = 0; i < dmp->nv - 1; i++) {
#ifdef DEBUG_LOOP
    printf("freeing face (%d) %d\n", i, dmp->face[i]);
#endif
    SDecimMeshFreeFace(dmp, dmp->face[i]);
  }
  if (dmp->face[dmp->nv - 1] != -1)
    SDecimMeshFreeFace(dmp, dmp->face[dmp->nv - 1]);
#ifdef DEBUG_LOOP
  printf("done\n");

#endif

  return (1);
}
