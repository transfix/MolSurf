#include <stdio.h>
#include <math.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include <stdlib.h>
#include <Decimation/DecimXu/xudecimate.h>
#include <Decimation/Linear/display.h>
#include <Decimation/Misc/misc.h>
#include <Decimation/geometry.h>

extern VertixNormal *pointnormals; /* original points and normal */
extern Face *facetetras;           /* original triangulation     */
extern int *facevalence, *edgevalence, numbpts, numbtris, newntri, steps,
    hvnormal;
extern Neighbor *faceneighbor, *edgeneighbor;

int *insidepts, /* grouping point set in tris */
    *vertDel,   /* -1-- removed; 0--not       */
    newnpts;
Face *facecopy;           /* decimated triangulation    */
VertixNormal *vertexcopy; /* decimated vertices         */

void LambdaPolyForVoume();
static int ReTriListAtVertex();
static void SecondLevelFaceList(), ModifyDataList(), MinimalAngleConnection(),
    MinimalDistanceConnection(), MinimalPathConnection_in(),
    MinimalPathConnection();
static float MiniAngleAroundVertex();
extern int Last_oid, animate_index, CurrentWin;

#define Model(i, n) ((i) < (n) ? (i) : (i - n))

/*-----------------------------------------------------------------------------
 * CallXuDecimate - Call Xu's decimation code
 *---------------------------------------------------------------------------*/
void CallXuDecimate(
    int niter,        // number of iteration
    float controllor, // ? (*PI radian measure)
    float variation,  // ? (*PI radian measure)
    float innerangle, // ? (*PI radian measure)
    int showpoly,     // show poly or not (doesn't matter)
    int grouping,     // ?? some flag
    int retriways,    // if (retriways  == 1) MinimalAngleConnection
                      // if (retriways  == 2) MinimalDistanceConnection
                      // if (retriways  == 3) MinimalPathConnection;
    int savenorm,     // save file flag
    const char *fname) {
  int i, j, k, i0, i1, i2, ii, jj, kk, v, t, goodtri;
  int trilist[300], info, ntri, deleat, facelist[200], nface, jjj,
      whichnllist[100], aroundvertex[100], numbver;
  int *tolptsindex, numbeachtri[100], tolnumb, running, maxtol;
  int *vertMap, *insideptsnew, msize;
  float p[3], p1[3], p2[3], p3[3], pp[3], n1[3], n2[3], n3[3], avernormal[3];
  float px[3], py[3], pz[3], q[4], bb[3], lambda, length, innercontrol;
  float v1[3], v2[3], v3[3], control, det, miniangle, angle;
  // GrObject      *decimate;
  Vertex *vertices;
  Normal *normals;
  FILE *fppts;

  int arandcount, arandtmp;
  float arandflt;

  arandcount = 0;

  /* copy the face data             */
  if (steps < 3) {
    facecopy = (Face *)malloc(numbtris * sizeof(Face));
    insidepts = (int *)malloc(numbpts * sizeof(int));
  }
  printf("step:%d numtris:%d numbpts:%d\n", steps, numbtris, numbpts);
  for (i = 0; i < numbtris; i++) {
    if (steps < 3) {
      facecopy[i] = (Face)malloc(sizeof(struct Face1));
    }
    facecopy[i]->Index[0] = facetetras[i]->Index[0];
    facecopy[i]->Index[1] = facetetras[i]->Index[1];
    facecopy[i]->Index[2] = facetetras[i]->Index[2];
    facecopy[i]->Normal[0] = facetetras[i]->Normal[0];
    facecopy[i]->Normal[1] = facetetras[i]->Normal[1];
    facecopy[i]->Normal[2] = facetetras[i]->Normal[2];

    // printf("Triangle %d %d %d\n", facecopy[i]->Index[0],
    // facecopy[i]->Index[1],facecopy[i]->Index[2]);

    // arand: sort indices...
    if (facecopy[i]->Index[0] > facecopy[i]->Index[1]) {
      arandtmp = facecopy[i]->Index[1];
      facecopy[i]->Index[1] = facecopy[i]->Index[0];
      facecopy[i]->Index[0] = arandtmp;
    }
    if (facecopy[i]->Index[1] > facecopy[i]->Index[2]) {
      arandtmp = facecopy[i]->Index[2];
      facecopy[i]->Index[2] = facecopy[i]->Index[1];
      facecopy[i]->Index[1] = arandtmp;
    }
    if (facecopy[i]->Index[0] > facecopy[i]->Index[1]) {
      arandtmp = facecopy[i]->Index[1];
      facecopy[i]->Index[1] = facecopy[i]->Index[0];
      facecopy[i]->Index[0] = arandtmp;
    }

    // printf("Triangle %d %d %d\n", facecopy[i]->Index[0],
    // facecopy[i]->Index[1],facecopy[i]->Index[2]);

    facecopy[i]->whichnml[0] = facetetras[i]->whichnml[0];
    facecopy[i]->whichnml[1] = facetetras[i]->whichnml[1];
    facecopy[i]->whichnml[2] = facetetras[i]->whichnml[2];
    facecopy[i]->Tou = 0;
    facecopy[i]->Wei = 0;
    facecopy[i]->State = 0;
  }

  msize = 3 * numbpts;
  vertDel = (int *)malloc(numbpts * sizeof(int));
  if (grouping == 1) {
    insideptsnew = (int *)malloc(msize * sizeof(int));
  }
  for (i = 0; i < numbpts; i++) {
    vertDel[i] = 0;
  }
  steps = 3;
  maxtol = 100;
  tolptsindex = (int *)malloc(maxtol * sizeof(int));
  det = (1.0 - controllor) / niter;
  printf("det:%f\n", det);
  /* Outer interation of decimation */
  for (t = 0; t < niter; t++) {
    printf("Begin Out Loop %d\n", t);
    running = 0;
    control = 1.0 - (t + 1) * det;
    if (t / 2 + t / 2 == t) {
      innercontrol = innerangle;
    } else {
      innercontrol = 2 * innerangle * innerangle - 1.0;
    }
    angle = 1.1;

    for (v = 0; v < numbpts; v++) {
      /*
        if ( (vertDel[v] > -1)  && (pointnormals[v]->SharpNormal == NULL) ) {
      */
      if ((vertDel[v] > -1) && (pointnormals[v]->SharpNormal == NULL) &&
          (facevalence[v] == edgevalence[v])) {
        if (t / 2 + t / 2 < t) {
          angle = MiniAngleAroundVertex(v);
        }
        if (angle > innercontrol) {

          /* produce triangle lists at v    */
          info = ReTriListAtVertex(v, trilist, &ntri, &miniangle, t, retriways,
                                   aroundvertex, &numbver);
          // printf("info:%d innercontrol:%f miniangle:%f
          // ntri:%d\n",info,innercontrol,miniangle, ntri);
          if ((info == 1) && (innercontrol > miniangle)) {

            /* if  every new triangle is good */
            // printf("new triangle is good\n");
            goodtri = 1;
            for (i = 0; i < ntri; i++) {
              ii = trilist[3 * i];
              jj = trilist[3 * i + 1];
              kk = trilist[3 * i + 2];
              for (j = 0; j < 3; j++) {
                n1[j] = pointnormals[ii]->Normal[j];
                n2[j] = pointnormals[jj]->Normal[j];
                n3[j] = pointnormals[kk]->Normal[j];
                v1[j] = pointnormals[jj]->Point[j] - pointnormals[ii]->Point[j];
                v2[j] = pointnormals[kk]->Point[j] - pointnormals[ii]->Point[j];
              }
              CrossProduct2(v1, v2, v3);
              length = sqrt(DotProduct(v3, v3));
              if ((fabs(DotProduct(n1, v3)) < control * length) ||
                  (fabs(DotProduct(n2, v3)) < control * length) ||
                  (fabs(DotProduct(n3, v3)) < control * length)) {
                goodtri = 0;
                // printf("%f %f %f %f\n", fabs(DotProduct(n1, v3)),
                // fabs(DotProduct(n2, v3)), fabs(DotProduct(n3, v3)),
                // control*length); printf("%f %f %f\n", n1[0], n1[1], n1[2]);

                break;
              }
            } /* end i loop                     */
            // printf("goodtri %d\n", goodtri);

            /* find the volume p located in   */
            deleat = 0;
            if (goodtri == 1) {
              // printf("find volume p located in\n");
              p[0] = pointnormals[v]->Point[0];
              p[1] = pointnormals[v]->Point[1];
              p[2] = pointnormals[v]->Point[2];

              /* count how many point are there */
              /* for malloc memery only         */
              if (grouping == 1) {
                tolnumb = 1; /* here 1 means p      */
                for (i = 0; i < facevalence[v]; i++) {
                  j = faceneighbor[v]->NB[i];
                  tolnumb = tolnumb + facecopy[j]->Wei - facecopy[j]->Tou;
                }
                if (tolnumb > maxtol) {
                  printf("tolnumb = %d > maxtol = %d\n", tolnumb, maxtol);
                  maxtol = tolnumb + 100;
                  tolptsindex =
                      (int *)realloc(tolptsindex, maxtol * sizeof(int));
                  if (tolptsindex == NULL) {
                    printf("Could not realloc\n");
                    return;
                  }
                }
              }

              tolnumb = 0;
              for (i = 0; i < ntri; i++) {
                // printf("tri %d\n",i);
                numbeachtri[i] = 0;
                ii = trilist[3 * i];
                jj = trilist[3 * i + 1];
                kk = trilist[3 * i + 2];
                for (j = 0; j < 3; j++) {
                  p1[j] = pointnormals[ii]->Point[j];
                  p2[j] = pointnormals[jj]->Point[j];
                  p3[j] = pointnormals[kk]->Point[j];
                  n1[j] = pointnormals[ii]->Normal[j];
                  n2[j] = pointnormals[jj]->Normal[j];
                  n3[j] = pointnormals[kk]->Normal[j];
                }
                LambdaPolyForVoume(p1, p2, p3, n1, n2, n3, px, py, pz, q);
                if (GlobalToLocal(p1, p2, p3, n1, n2, n3, p, px, py, pz, q, bb,
                                  &lambda) > 0) {

                  avernormal[0] = bb[0] * n1[0] + bb[1] * n2[0] + bb[2] * n3[0];
                  avernormal[1] = bb[0] * n1[1] + bb[1] * n2[1] + bb[2] * n3[1];
                  avernormal[2] = bb[0] * n1[2] + bb[1] * n2[2] + bb[2] * n3[2];
                  length = sqrt(DotProduct(avernormal, avernormal));
                  length =
                      DotProduct(avernormal, pointnormals[v]->Normal) / length;
                  if (length > variation) {
                    deleat = 1;
                    tolptsindex[tolnumb] = v;
                    tolnumb = tolnumb + 1;
                    numbeachtri[i] = 1;
                  } else {
                    deleat = 0;
                    goto L1;
                  }
                } /* end if                 */

                if (grouping == 1) {
                  for (j = 0; j < facevalence[v]; j++) {
                    k = faceneighbor[v]->NB[j];
                    i1 = facecopy[k]->Tou;
                    i2 = facecopy[k]->Wei;
                    if (facecopy[k]->State == 0) {
                      for (k = i1; k < i2; k++) {
                        i0 = insidepts[k];
                        pp[0] = pointnormals[i0]->Point[0];
                        pp[1] = pointnormals[i0]->Point[1];
                        pp[2] = pointnormals[i0]->Point[2];
                        if (GlobalToLocal(p1, p2, p3, n1, n2, n3, pp, px, py,
                                          pz, q, bb, &lambda) > 0) {

                          avernormal[0] =
                              bb[0] * n1[0] + bb[1] * n2[0] + bb[2] * n3[0];
                          avernormal[1] =
                              bb[0] * n1[1] + bb[1] * n2[1] + bb[2] * n3[1];
                          avernormal[2] =
                              bb[0] * n1[2] + bb[1] * n2[2] + bb[2] * n3[2];
                          length = sqrt(DotProduct(avernormal, avernormal));
                          length =
                              DotProduct(avernormal, pointnormals[v]->Normal) /
                              length;

                          if (length > variation) {
                            deleat = 1;
                            tolptsindex[tolnumb] = i0;
                            tolnumb = tolnumb + 1;
                            numbeachtri[i] = numbeachtri[i] + 1;
                          } else {
                            deleat = 0;
                            goto L1;
                          }
                        }
                      } /* k loop */
                    } /* end if */
                    if (facecopy[faceneighbor[v]->NB[j]]->State == 2) {
                      for (k = i1; k < i2; k++) {
                        i0 = insideptsnew[k];
                        pp[0] = pointnormals[i0]->Point[0];
                        pp[1] = pointnormals[i0]->Point[1];
                        pp[2] = pointnormals[i0]->Point[2];
                        if (GlobalToLocal(p1, p2, p3, n1, n2, n3, pp, px, py,
                                          pz, q, bb, &lambda) > 0) {

                          avernormal[0] =
                              bb[0] * n1[0] + bb[1] * n2[0] + bb[2] * n3[0];
                          avernormal[1] =
                              bb[0] * n1[1] + bb[1] * n2[1] + bb[2] * n3[1];
                          avernormal[2] =
                              bb[0] * n1[2] + bb[1] * n2[2] + bb[2] * n3[2];
                          length = sqrt(DotProduct(avernormal, avernormal));
                          length =
                              DotProduct(avernormal, pointnormals[v]->Normal) /
                              length;

                          if (length > variation) {
                            deleat = 1;
                            tolptsindex[tolnumb] = i0;
                            tolnumb = tolnumb + 1;
                            numbeachtri[i] = numbeachtri[i] + 1;
                          } else {
                            deleat = 0;
                            goto L1;
                          }
                        }
                      } /* end k loop            */
                    } /* end if                */
                  } /* end j loop            */
                } /* end if grouping == 1  */

              } /* end i loop                     */
            } /* end if (goodtri == 1 )         */

            /* deleat point                   */
          L1:;
            if ((deleat == 1) && (goodtri == 1)) {

              /* Adjust face normal, right-hand */
              /* record the which normal used   */
              for (i = 0; i < numbver; i++) {
                if (pointnormals[aroundvertex[i]]->SharpNormal != NULL) {
                  for (j = 0; j < facevalence[v]; j++) {
                    jjj = faceneighbor[v]->NB[j];
                    if (aroundvertex[i] == facecopy[jjj]->Index[0]) {
                      whichnllist[i] = facecopy[jjj]->whichnml[0];
                      break;
                    }
                    if (aroundvertex[i] == facecopy[jjj]->Index[1]) {
                      whichnllist[i] = facecopy[jjj]->whichnml[1];
                      break;
                    }
                    if (aroundvertex[i] == facecopy[jjj]->Index[2]) {
                      whichnllist[i] = facecopy[jjj]->whichnml[2];
                      break;
                    }
                  }
                }
              }

              /* Adjust data = -1 means deleated*/
              vertDel[v] = -1;

              SecondLevelFaceList(v, facelist, &nface);

              // arand - seems ok to here

              /* Adjust face normal, right-hand */
              for (i = 0; i < ntri; i++) {
                ii = trilist[3 * i];
                jj = trilist[3 * i + 1];
                kk = trilist[3 * i + 2];
                for (j = 0; j < 3; j++) {
                  v1[j] =
                      pointnormals[jj]->Point[j] - pointnormals[ii]->Point[j];
                  v2[j] =
                      pointnormals[kk]->Point[j] - pointnormals[ii]->Point[j];
                }
                CrossProduct2(v1, v2, n2);
                length = 1.0 / sqrt(DotProduct(n2, n2));
                jj = faceneighbor[v]->NB[i];
                facecopy[jj]->Normal[0] = n2[0] * length;
                facecopy[jj]->Normal[1] = n2[1] * length;
                facecopy[jj]->Normal[2] = n2[2] * length;
              }

              /* order  triangle index          */
              for (i = 0; i < ntri; i++) {
                ii = trilist[3 * i];
                jj = trilist[3 * i + 1];
                kk = trilist[3 * i + 2];
                if (ii < kk && kk < jj) {
                  trilist[3 * i + 1] = kk;
                  trilist[3 * i + 2] = jj;
                }
                if (jj < ii && ii < kk) {
                  trilist[3 * i] = jj;
                  trilist[3 * i + 1] = ii;
                }
                if (jj < kk && kk < ii) {
                  trilist[3 * i] = jj;
                  trilist[3 * i + 1] = kk;
                  trilist[3 * i + 2] = ii;
                }
                if (kk < ii && ii < jj) {
                  trilist[3 * i] = kk;
                  trilist[3 * i + 1] = ii;
                  trilist[3 * i + 2] = jj;
                }
                if (kk < jj && jj < ii) {
                  trilist[3 * i] = kk;
                  trilist[3 * i + 1] = jj;
                  trilist[3 * i + 2] = ii;
                }
              }

              /* Adjust  face triangle index    */
              tolnumb = 0;

              for (j = 0; j < ntri; j++) {
                jj = faceneighbor[v]->NB[j];
                // printf("Adding %d : %d %d %d\n",
                // jj,trilist[3*j],trilist[3*j+1],trilist[3*j+2]);

                facecopy[jj]->Index[0] = trilist[3 * j];
                facecopy[jj]->Index[1] = trilist[3 * j + 1];
                facecopy[jj]->Index[2] = trilist[3 * j + 2];
                facecopy[jj]->whichnml[0] = -1;
                facecopy[jj]->whichnml[1] = -1;
                facecopy[jj]->whichnml[2] = -1;

                /* give the normal index          */
                for (i = 0; i < numbver; i++) {
                  if (pointnormals[aroundvertex[i]]->SharpNormal != NULL) {
                    if (trilist[3 * j] == aroundvertex[i])
                      facecopy[jj]->whichnml[0] = whichnllist[i];
                    if (trilist[3 * j + 1] == aroundvertex[i])
                      facecopy[jj]->whichnml[1] = whichnllist[i];
                    if (trilist[3 * j + 2] == aroundvertex[i])
                      facecopy[jj]->whichnml[2] = whichnllist[i];
                  }
                }
                if (grouping == 1) {
                  i0 = facecopy[jj]->Tou;          /* old tou  */
                  i1 = facecopy[jj]->Wei;          /* old wei, need to change */
                  if (i1 - i0 >= numbeachtri[j]) { /* no need to malloc, use
                                                      original space , keep the
                                                      State unchanged   */
                    facecopy[jj]->Wei = i0 + numbeachtri[j]; /* new wei   */
                    if (facecopy[jj]->State == 0) {
                      for (k = 0; k < numbeachtri[j]; k++) {
                        insidepts[i0 + k] = tolptsindex[tolnumb];
                        tolnumb = tolnumb + 1;
                      }
                    }
                    if (facecopy[jj]->State == 2) {
                      for (k = 0; k < numbeachtri[j]; k++) {
                        insideptsnew[i0 + k] = tolptsindex[tolnumb];
                        tolnumb = tolnumb + 1;
                      }
                    }
                  } else { /* malloc new space */
                    facecopy[jj]->State = 2;
                    facecopy[jj]->Tou = running;
                    facecopy[jj]->Wei = running + numbeachtri[j];

                    for (k = 0; k < numbeachtri[j]; k++) {
                      if (running >= msize) {
                        msize = msize + numbpts;
                        insideptsnew =
                            (int *)realloc(insideptsnew, msize * sizeof(int));
                        if (insideptsnew == NULL) {
                          printf("Could not realloc\n");
                          return;
                        }
                      }

                      insideptsnew[running] = tolptsindex[tolnumb];
                      running = running + 1;
                      tolnumb = tolnumb + 1;
                    }
                  }
                } /* end if (grouping == 1) */
              } /* end j loop             */

              for (j = ntri; j < facevalence[v]; j++) {
                jj = faceneighbor[v]->NB[j];
                facecopy[jj]->State = 1; /* means deleated */
              }

              /* Adjust  data faceneighbor      */
              ModifyDataList(v, facelist, nface, trilist, ntri);

              // arandcount++;
              // if (arandcount >=2)
              // goto L2;

            } /* end if ( (deleat == 1) && ...  */
          } /* end if (info == 1)             */
        } /* end if (angle > miniangle)     */
      } /* end if vertDel[v] > -1         */
    } /* end v loop                     */

    /* complet the insideptsnew       */
    if (grouping == 1) {
      for (i = 0; i < numbtris; i++) {
        if (facecopy[i]->State == 0) {
          i0 = facecopy[i]->Tou;
          i1 = facecopy[i]->Wei;
          facecopy[i]->Tou = running;
          facecopy[i]->Wei = running + (i1 - i0);
          for (j = i0; j < i1; j++) {
            if (running >= msize) {
              msize = msize + numbpts;
              insideptsnew = (int *)realloc(insideptsnew, msize * sizeof(int));
              if (insideptsnew == NULL) {
                printf("Could not realloc\n");
                return;
              }
            }
            insideptsnew[running] = insidepts[j];
            running = running + 1;
          }
        }
        if (facecopy[i]->State == 2) {
          facecopy[i]->State = 0;
        }
      } /* end i loop                     */
      printf("running = %d \n", running);

      /* move the new to old            */
      running = 0;
      for (i = 0; i < numbtris; i++) {
        if (facecopy[i]->State == 0) {
          i0 = facecopy[i]->Tou;
          i1 = facecopy[i]->Wei;
          facecopy[i]->Tou = running;
          facecopy[i]->Wei = running + (i1 - i0);
          for (j = i0; j < i1; j++) {
            insidepts[running] = insideptsnew[j];
            running = running + 1;
          }
        }
      } /* end i loop                     */
    } /* end (grouping == 1)            */
    printf("grouping points = %d \n", running);

  } /* end t loop                     */

L2: // for debug purposes only...

  /* index map of the two meshes    */
  vertMap = (int *)malloc(numbpts * sizeof(int));
  newnpts = 0;
  for (i = 0; i < numbpts; i++) {
    if (vertDel[i] > -1) {
      vertMap[i] = newnpts;
      newnpts = newnpts + 1;
    }
  }

  /* move vertex to vertexcopy         */
  /*
    vertexcopy = (VertixNormal *) malloc(newnpts*sizeof(VertixNormal));
    newnpts = 0;
    for (i = 0; i < numbpts; i++) {
    if (vertDel[i] > -1) {
    vertexcopy[newnpts] = (VertixNormal)malloc(sizeof(struct VertixNormal));
    vertexcopy[newnpts]->Point[0] = pointnormals[i]->Point[0];
    vertexcopy[newnpts]->Point[1] = pointnormals[i]->Point[1];
    vertexcopy[newnpts]->Point[2] = pointnormals[i]->Point[2];

    vertexcopy[newnpts]->Normal[0] = pointnormals[i]->Normal[0];
    vertexcopy[newnpts]->Normal[1] = pointnormals[i]->Normal[1];
    vertexcopy[newnpts]->Normal[2] = pointnormals[i]->Normal[2];
    newnpts = newnpts + 1;
    }
    }
  */

  /* move face data forward          */
  newntri = 0;
  for (i = 0; i < numbtris; i++) {
    if (facecopy[i]->State != 1) {
      facecopy[newntri]->Index[0] = facecopy[i]->Index[0];
      facecopy[newntri]->Index[1] = facecopy[i]->Index[1];
      facecopy[newntri]->Index[2] = facecopy[i]->Index[2];
      facecopy[newntri]->Normal[0] = facecopy[i]->Normal[0];
      facecopy[newntri]->Normal[1] = facecopy[i]->Normal[1];
      facecopy[newntri]->Normal[2] = facecopy[i]->Normal[2];
      facecopy[newntri]->Tou = facecopy[i]->Tou;
      facecopy[newntri]->Wei = facecopy[i]->Wei;
      facecopy[newntri]->State = 0;
      facecopy[newntri]->whichnml[0] = facecopy[i]->whichnml[0];
      facecopy[newntri]->whichnml[1] = facecopy[i]->whichnml[1];
      facecopy[newntri]->whichnml[2] = facecopy[i]->whichnml[2];
      newntri = newntri + 1;
    }
  }
  /* save the decimated trianglation*/
  // fppts = fopen("output.rawn","w");
  fppts = fopen(fname, "w");
  fprintf(fppts, "%d %d\n", newnpts, newntri);
  // printf("in xuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu\n");

  /* write vertices normals  first  */
  if (savenorm == 1) {
    for (i = 0; i < numbpts; i++) {
      if (vertDel[i] > -1) {
        fprintf(fppts, "%lf %lf %lf\n", pointnormals[i]->Point[0],
                pointnormals[i]->Point[1], pointnormals[i]->Point[2]);

        // fprintf(fppts, "%lf %lf %lf %lf %lf %lf\n",
        // pointnormals[i]->Point[0], pointnormals[i]->Point[1],
        // pointnormals[i]->Point[2],
        // pointnormals[i]->Normal[0],pointnormals[i]->Normal[1],pointnormals[i]->Normal[2]);
      }
    }
  }
  if (savenorm == 0) {
    for (i = 0; i < numbpts; i++) {
      if (vertDel[i] > -1) {
        fprintf(fppts, "%lf %lf %lf\n",

                pointnormals[i]->Point[0], pointnormals[i]->Point[1],
                pointnormals[i]->Point[2]);
        /* for brain
           pointnormals[i]->Point[2], pointnormals[i]->Point[0],
           -pointnormals[i]->Point[1]);
        */
      }
    }
  }

  /* write triangle indices second  */
  for (i = 0; i < newntri; i++) {
    /*
      fprintf(fppts, "%d %d %d\n", vertMap[facecopy[i]->Index[0]],
      vertMap[facecopy[i]->Index[1]],
      vertMap[facecopy[i]->Index[2]]);
    */

    // check orientation?
    float v1[3], v2[3], v3[3];
    v1[0] = pointnormals[facecopy[i]->Index[1]]->Point[0] -
            pointnormals[facecopy[i]->Index[0]]->Point[0];
    v1[1] = pointnormals[facecopy[i]->Index[1]]->Point[1] -
            pointnormals[facecopy[i]->Index[0]]->Point[1];
    v1[2] = pointnormals[facecopy[i]->Index[1]]->Point[2] -
            pointnormals[facecopy[i]->Index[0]]->Point[2];
    v2[0] = pointnormals[facecopy[i]->Index[2]]->Point[0] -
            pointnormals[facecopy[i]->Index[0]]->Point[0];
    v2[1] = pointnormals[facecopy[i]->Index[2]]->Point[1] -
            pointnormals[facecopy[i]->Index[0]]->Point[1];
    v2[2] = pointnormals[facecopy[i]->Index[2]]->Point[2] -
            pointnormals[facecopy[i]->Index[0]]->Point[2];
    CrossProduct2(v1, v2, v3);

    // printf("CP %f %f %f %f %f %f %f %f %f\n", v1[0], v1[1], v1[2], v2[0],
    // v2[1], v2[2], v3[0], v3[1], v3[2]); printf("Nor %f %f %f\n",
    // facecopy[i]->Normal[0],facecopy[i]->Normal[1],facecopy[i]->Normal[2]);
    ii = vertMap[facecopy[i]->Index[0]];
    jj = vertMap[facecopy[i]->Index[1]];
    kk = vertMap[facecopy[i]->Index[2]];
    /*
      if (ii > 9968) ii = ii - 1;
      if (jj > 9968) jj = jj - 1;
      if (kk > 9968) kk = kk - 1;
    */

    if (DotProduct(v3, facecopy[i]->Normal) > 0) {
      fprintf(fppts, "%d %d %d\n", ii, kk, jj);
    } else {
      fprintf(fppts, "%d %d %d\n", ii, jj, kk);
    }
  }

  /*chnage the triangle index to new*/
  /*
    for (i = 0; i < newntri; i++) {
    facecopy[i]->Index[0] = vertMap[facecopy[i]->Index[0]];
    facecopy[i]->Index[1] = vertMap[facecopy[i]->Index[1]];
    facecopy[i]->Index[2] = vertMap[facecopy[i]->Index[2]];
    }
  */

  fclose(fppts);
  printf("Decimated Mesh points =%d, triangles =%d origin_vertex =%d,  "
         "origin_tri =%d\n",
         newnpts, newntri, numbpts, numbtris);

  /* the following is for produce cut data, from the decimate output tmp*/
  /* write to tmp2                                                       */
  /*
    if (showpoly == 1) {
    int  jige,tol;
    Triangle tri;
    float p1[3],p2[3],p3[3],n1[3],n2[3],n3[3];
    fppts = fopen("../../data/Triangles/tmp2","w");

    tol = 0;
    fprintf(fppts,"%10d\n", tol);
    for (i = 0; i < newntri; i++) {
    ii = facecopy[i]->Index[0];
    jj = facecopy[i]->Index[1];
    kk = facecopy[i]->Index[2];

    for (j = 0; j < 3; j++) {
    p1[j] = pointnormals[ii]->Point[j];
    p2[j] = pointnormals[jj]->Point[j];
    p3[j] = pointnormals[kk]->Point[j];
    n1[j] = pointnormals[ii]->Normal[j];
    n2[j] = pointnormals[jj]->Normal[j];
    n3[j] = pointnormals[kk]->Normal[j];
    }
    tri = (Triangle) malloc(sizeof(struct Triangle));

    jige = CutTriangleByTwoPlane(p1,p2,p3,n1,n2,n3,0.5,1.2,tri);
    if (jige == 1) {
    fprintf(fppts, "%lf %lf %lf %lf %lf %lf\n",
    tri->Vertex1[0], tri->Vertex1[1],tri->Vertex1[2],
    tri->Normal1[0], tri->Normal1[1],tri->Normal1[2]);

    fprintf(fppts, "%lf %lf %lf %lf %lf %lf\n",
    tri->Vertex2[0], tri->Vertex2[1],tri->Vertex2[2],
    tri->Normal2[0], tri->Normal2[1],tri->Normal2[2]);

    fprintf(fppts, "%lf %lf %lf %lf %lf %lf\n",
    tri->Vertex3[0], tri->Vertex3[1],tri->Vertex3[2],
    tri->Normal3[0], tri->Normal3[1],tri->Normal3[2]);
    tol = tol + 1;
    }
    }
    rewind(fppts);
    fprintf(fppts,"%10d\n", tol);
    printf("%d\n", tol);
    fclose(fppts);
    TransferFile_raw();
    TransferToFat();
    }
  */

  /* draw the decimated mesh        */
  /* free memories                  */
  free(tolptsindex);
  free(vertMap);
  if (grouping == 1) {
    free(insideptsnew);
  }
} /* end the routine                */

/*-----------------------------------------------------------------------------
 * ModifyDataList - modify the data lists
 *---------------------------------------------------------------------------*/
static void ModifyDataList(v, facelist, nface, trilist, ntri) int v, facelist[],
    nface, trilist[], ntri;
{
  int vala, i, j, k, l, m, ii, jj, kk, ifequal, iface, cur;
  Neighbor faceneigh, edgeneigh;
  int *facenb, *edgenb;

  // printf("ModifyDataList\n");
  // for (i=0; i<ntri; i++) {
  //   printf("Tri %d %d %d \n", trilist[3*i],trilist[3*i+1],trilist[3*i+2]);
  // }

  vala = edgevalence[v];
  for (i = 0; i < vala; i++) {
    j = edgeneighbor[v]->NB[i]; /*  vertex index */
    facevalence[j] = 0;
    edgevalence[j] = 0;
    free(faceneighbor[j]->NB);
    free(edgeneighbor[j]->NB);
  }
  /* add the facevalence for second level tri */
  for (i = 0; i < nface; i++) {
    j = facelist[i]; /* face index    */
    ii = facecopy[j]->Index[0];
    jj = facecopy[j]->Index[1];
    kk = facecopy[j]->Index[2];
    for (k = 0; k < vala; k++) {
      j = edgeneighbor[v]->NB[k]; /* vertex index  */
      if ((ii == j) || (jj == j) || (kk == j))
        facevalence[j] = facevalence[j] + 1;
    }
  }
  /* add the facevalence for first level tri  */
  for (i = 0; i < ntri; i++) {
    ii = trilist[3 * i];
    jj = trilist[3 * i + 1];
    kk = trilist[3 * i + 2];
    for (k = 0; k < vala; k++) {
      j = edgeneighbor[v]->NB[k]; /* vertex index  */
      if ((ii == j) || (jj == j) || (kk == j))
        facevalence[j] = facevalence[j] + 1;
    }
  }

  /* maccloc neighbor tri index                   */
  for (i = 0; i < vala; i++) {
    j = edgeneighbor[v]->NB[i]; /* vertex index  */
    faceneigh = (Neighbor)malloc(sizeof(Neighbor));
    if (faceneigh == NULL) {
      printf("malloc fail in ModifyDataList 1\n");
      exit(0);
    }
    facenb = (int *)malloc(facevalence[j] * sizeof(int));
    if (facenb == NULL) {
      printf("malloc fail in ModifyDataList 2\n");
      exit(0);
    }
    facenb[facevalence[j] - 1] = -1;
    faceneigh->NB = facenb;
    faceneighbor[j] = faceneigh;

    /* if surface has boundary, the
       facevalence[j] is incoorect                  */
    edgeneigh = (Neighbor)malloc(sizeof(Neighbor));
    if (edgeneigh == NULL) {
      printf("malloc fail in ModifyDataList 3\n");
      exit(0);
    }
    edgenb = (int *)malloc((facevalence[j] + 2) * sizeof(int));
    if (edgenb == NULL) {
      printf("malloc fail in ModifyDataList 4\n");
      exit(0);
    }
    edgeneigh->NB = edgenb;
    edgeneighbor[j] = edgeneigh;
  }

  /* compute faceneighbor index of second level   */
  for (i = 0; i < nface; i++) {
    j = facelist[i];            /* face index    */
    ii = facecopy[j]->Index[0]; /* vertex index  */
    jj = facecopy[j]->Index[1];
    kk = facecopy[j]->Index[2];
    for (k = 0; k < vala; k++) {
      j = edgeneighbor[v]->NB[k]; /* vertex index  */
      if ((ii == j) || (jj == j) || (kk == j)) {
        l = faceneighbor[j]->NB[facevalence[j] - 1] + 1;
        faceneighbor[j]->NB[facevalence[j] - 1] = l;
        faceneighbor[j]->NB[l] = facelist[i];
      }
    }
  }
  /* add the facenebour  for first level tri     */
  for (i = 0; i < ntri; i++) {
    ii = trilist[3 * i];
    jj = trilist[3 * i + 1];
    kk = trilist[3 * i + 2];
    for (k = 0; k < vala; k++) {
      j = edgeneighbor[v]->NB[k]; /* vertex index  */
      if ((ii == j) || (jj == j) || (kk == j)) {
        l = faceneighbor[j]->NB[facevalence[j] - 1] + 1;
        faceneighbor[j]->NB[facevalence[j] - 1] = l;
        faceneighbor[j]->NB[l] = faceneighbor[v]->NB[i];
      }
    }
  }

  /* compute neighbor vertex index                 */
  for (k = 0; k < vala; k++) {
    i = edgeneighbor[v]->NB[k];
    cur = 0;
    for (j = 0; j < facevalence[i]; j++) {
      iface = faceneighbor[i]->NB[j];

      for (l = 0; l < 3; l++) {
        ii = facecopy[iface]->Index[l];
        if (ii != i) {
          ifequal = 0;
          for (m = 0; m < cur; m++) {
            if (ii == edgeneighbor[i]->NB[m])
              ifequal = 1;
          }
          if (ifequal == 0 && cur <= facevalence[i] + 1) {
            edgeneighbor[i]->NB[cur] = ii;
            cur = cur + 1;
          }
        }
      }
    }
    edgevalence[i] = cur;
  }
}

/*-----------------------------------------------------------------------------
 * SecondLevelFaceList - Compute second level face list around vertex v
 *---------------------------------------------------------------------------*/
static void SecondLevelFaceList(l, facelist, nface) int l, facelist[], *nface;
{
  int i, j, k, jj, kk, ll, curentvalence, curvalence2, numb, ifequal;

  // printf("Enter SecondLevel\n");

  curentvalence = edgevalence[l];
  curvalence2 = facevalence[l];
  // printf("V %d %d\n", curentvalence,curvalence2);

  numb = 0;
  for (i = 0; i < curentvalence; i++) {
    j = edgeneighbor[l]->NB[i]; /* vertex index  */

    // printf("Immediate Vertex Neighbor: %d %d %d\n", curentvalence, i,j);

    /* level 2 data, limit to 200 points         */
    kk = facevalence[j];
    for (k = 0; k < kk; k++) {
      jj = faceneighbor[j]->NB[k]; /* face index    */

      ifequal = 0;
      for (ll = 0; ll < curentvalence; ll++) {
        if (faceneighbor[l]->NB[ll] == jj) {
          ifequal = 1;
          break;
        }
      }
      if (ifequal == 0) {
        for (ll = 0; ll < numb; ll++) {
          if (facelist[ll] == jj) {
            ifequal = 1;
            break;
          }
        }
      }
      //  printf("Indirect Face Neighbor: %d %d %d, ie %d\n", kk, k,jj,ifequal);

      if (ifequal == 0 && numb < 200) {
        // printf("Save face %d, %d\n", jj,numb);
        facelist[numb] = jj;
        numb = numb + 1;
      }
    } /* end k loop      */
  } /* end i loop      */
  *nface = numb;
  // printf("SLFL %d %d : ", l,numb);
  // for (i=0; i<numb; i++) {
  //   printf("%d ", facelist[i]);
  //  }
  // printf("\n");
}

/*-----------------------------------------------------------------------------
 * ReTriListAtVertex - Compute the re-triangle list at a vertex v
 *---------------------------------------------------------------------------*/
static int ReTriListAtVertex(l, trilist, ntri, smallangle, tt, retriways,
                             around, numbver)
int l, trilist[], *ntri, tt, retriways, around[], *numbver;
float *smallangle;
{
  int i, j, k, ii, jj, kk, iii, jjj, kkk, i0, j0, k0, i1, i2, curentvalence;
  int aroundvertex[100], tou, wei, newtou, newwei, curnum, w;
  float v1[3], v2[3], v3[3], angle;
  float n1[3], sum;

  curentvalence = facevalence[l];
  if (curentvalence < 2) {
    printf("facevalence[%d] is less than 2\n", l);
    return (0);
  }
  if (curentvalence > 100) {
    printf("facevalence[%d] is larger than 100\n", l);
    return (0);
  }
  /* Compute around vertex lists    */
  /* initialization                 */
  i0 = faceneighbor[l]->NB[0]; /* current triangle */
  ii = facecopy[i0]->Index[0];
  jj = facecopy[i0]->Index[1];
  kk = facecopy[i0]->Index[2];
  // printf("i0 %d, ii %d, jj %d, kk %d, l %d\n",i0,ii,jj,kk,l);

  if (ii == l) {
    aroundvertex[0] = jj;
    aroundvertex[1] = kk;
    tou = jj;
    wei = kk;
  }
  if (jj == l) {
    aroundvertex[0] = ii;
    aroundvertex[1] = kk;
    tou = ii;
    wei = kk;
  }
  if (kk == l) {
    aroundvertex[0] = ii;
    aroundvertex[1] = jj;
    tou = ii;
    wei = jj;
  }
  curnum = 2;

  /* begin loop                     */
  i = 0;
  while (wei != tou && curnum <= curentvalence + 1 && i < curentvalence) {
    i = i + 1;
    for (j = 0; j < curentvalence; j++) {
      j0 = faceneighbor[l]->NB[j];
      if (i0 != j0) {
        iii = facecopy[j0]->Index[0];
        jjj = facecopy[j0]->Index[1];
        kkk = facecopy[j0]->Index[2];
        w = -1;
        if (wei == iii) {
          w = iii;
          if (jjj == l)
            newwei = kkk;
          if (kkk == l)
            newwei = jjj;
        }
        if (wei == jjj) {
          w = jjj;
          if (iii == l)
            newwei = kkk;
          if (kkk == l)
            newwei = iii;
        }
        if (wei == kkk) {
          w = kkk;
          if (iii == l)
            newwei = jjj;
          if (jjj == l)
            newwei = iii;
        }

        if (w >= 0) {
          wei = newwei;
          if (newwei != tou) {
            i0 = j0;
            aroundvertex[curnum] = wei;
            curnum = curnum + 1;
          }
        }
        // printf("i0 %d, iii %d, jjj %d, kkk %d, l %d, j0
        // %d\n",i0,iii,jjj,kkk,l,j0);
      }
    } /* end j loop  */
  } /* end while   */

  if (curnum < curentvalence) {
    /* search in other direction      */
    i = 0;
    i0 = faceneighbor[l]->NB[0];
    k0 = 0;
    while ((curnum <= curentvalence + 1) && (i < curentvalence)) {
      i = i + 1;
      for (j = 0; j < curentvalence; j++) {
        j0 = faceneighbor[l]->NB[j];
        if (i0 != j0) {
          iii = facecopy[j0]->Index[0];
          jjj = facecopy[j0]->Index[1];
          kkk = facecopy[j0]->Index[2];
          w = -1;
          if (tou == iii) {
            w = iii;
            if (jjj == l)
              newtou = kkk;
            if (kkk == l)
              newtou = jjj;
          }
          if (tou == jjj) {
            w = jjj;
            if (iii == l)
              newtou = kkk;
            if (kkk == l)
              newtou = iii;
          }
          if (tou == kkk) {
            w = kkk;
            if (iii == l)
              newtou = jjj;
            if (jjj == l)
              newtou = iii;
          }

          if (w >= 0) {
            tou = newtou;
            i0 = j0;
            aroundvertex[curnum] = tou;
            curnum = curnum + 1;
            k0 = k0 + 1;
          }
          // printf("2* i0 %d, iii %d, jjj %d, kkk %d, l %d, j0
          // %d\n",i0,iii,jjj,kkk,l,j0);
        }
      } /* end j loop                     */
    } /* end while                      */

    printf("            At vertex %d, edgenumber, facenumber = %d,%d \n", l,
           curnum, curentvalence);

    /* reorder the last part data     */
    for (i = 0; i < k0 / 2; i++) {
      i1 = aroundvertex[curnum - k0 + i];
      i2 = aroundvertex[curnum - 1 - i];
      aroundvertex[curnum - k0 + i] = i2;
      aroundvertex[curnum - 1 - i] = i1;
    }
  } /* end if curnum < valence        */

  *numbver = curnum;
  for (i = 0; i < curnum; i++)
    around[i] = aroundvertex[i];

  if (curnum > curentvalence + 1)
    return (0);
  if (curnum <= 2)
    return (0);

  /* check it is left or right hand */
  for (iii = 0; iii < 3; iii++) {
    v1[iii] =
        pointnormals[aroundvertex[0]]->Point[iii] - pointnormals[l]->Point[iii];
    v2[iii] =
        pointnormals[aroundvertex[1]]->Point[iii] - pointnormals[l]->Point[iii];
  }
  CrossProduct2(v1, v2, n1);
  sum = DotProduct(n1, facecopy[faceneighbor[l]->NB[0]]->Normal);
  if (sum < 0.0) {

    /* reverse the aroundvertex list  */
    k = curnum / 2;
    for (i = 0; i < k; i++) {
      i1 = aroundvertex[i];
      i2 = aroundvertex[curnum - 1 - i];
      aroundvertex[i] = i2;
      aroundvertex[curnum - 1 - i] = i1;
    }
  }
  k = curnum;
  // printf("%d AroundVertex ",l);
  // for (i = 0; i < k; i++ ) {
  //   printf("%d ", aroundvertex[i]);
  // }
  // printf("\n");

  /* form triangles lists           */
  *ntri = 0;
  if (retriways == 1)
    MinimalAngleConnection(aroundvertex, curnum, ntri, trilist);
  if (retriways == 2)
    MinimalDistanceConnection(aroundvertex, ntri, trilist, curnum);
  if (retriways == 3)
    MinimalPathConnection(l, aroundvertex, ntri, trilist, curnum);

  // for (i = 0; i < *ntri; i++) {
  //   printf("NT %d %d %d\n", trilist[3*i], trilist[3*i+1], trilist[3*i+2]);
  // }

  /* compute the minimal angle      */
  *smallangle = -1.0;
  for (i = 0; i < *ntri; i++) {
    ii = trilist[3 * i];
    jj = trilist[3 * i + 1];
    kk = trilist[3 * i + 2];

    for (j = 0; j < 3; j++) {
      v1[j] = pointnormals[jj]->Point[j] - pointnormals[ii]->Point[j];
      v2[j] = pointnormals[kk]->Point[j] - pointnormals[ii]->Point[j];
      v3[j] = pointnormals[kk]->Point[j] - pointnormals[jj]->Point[j];
    }
    angle = DotProduct(v1, v2) / sqrt(DotProduct(v1, v1) * DotProduct(v2, v2));
    if (*smallangle < angle)
      *smallangle = angle;

    angle = -DotProduct(v1, v3) / sqrt(DotProduct(v1, v1) * DotProduct(v3, v3));
    if (*smallangle < angle)
      *smallangle = angle;

    angle = DotProduct(v2, v3) / sqrt(DotProduct(v2, v2) * DotProduct(v3, v3));
    if (*smallangle < angle)
      *smallangle = angle;
  }
  return (1);
}

/*-----------------------------------------------------------------------------
 * MinimalPathConnection_in -- Form the triangulation by connect ponts that
 *                              has minimal geodesic distance
 *---------------------------------------------------------------------------*/
static void MinimalPathConnection_in(aroundvertex, ntri, trilist, number,
                                     minidis, index,
                                     tolnumber) int aroundvertex[],
    *ntri, trilist[], number, index[], tolnumber;
float minidis[];
{
  int i, j, i0, i1, i2, anglesdir[100], number1, number2, minii1, minii2, truei,
      truej;
  int aroundvertex1[100], aroundvertex2[100], index1[100], index2[100], nn2,
      nni2;
  float mindis, dis, v1[3], v2[3], v3[3], n2[3];

  if (number == 3) {
    trilist[3 * *ntri] = aroundvertex[0];
    trilist[3 * *ntri + 1] = aroundvertex[1];
    trilist[3 * *ntri + 2] = aroundvertex[2];
    *ntri = *ntri + 1;
    return;
  }
  /* Compute angles at each vertex  */
  for (i = 0; i < number; i++) {
    anglesdir[i] = 0;
    i1 = i - 1;
    if (i == 0)
      i1 = number - 1;
    i2 = i + 1;
    if (i == number - 1)
      i2 = 0;

    v3[0] = pointnormals[aroundvertex[i]]->Normal[0];
    v3[1] = pointnormals[aroundvertex[i]]->Normal[1];
    v3[2] = pointnormals[aroundvertex[i]]->Normal[2];

    v1[0] = pointnormals[aroundvertex[i1]]->Point[0] -
            pointnormals[aroundvertex[i]]->Point[0];
    v1[1] = pointnormals[aroundvertex[i1]]->Point[1] -
            pointnormals[aroundvertex[i]]->Point[1];
    v1[2] = pointnormals[aroundvertex[i1]]->Point[2] -
            pointnormals[aroundvertex[i]]->Point[2];

    v2[0] = pointnormals[aroundvertex[i2]]->Point[0] -
            pointnormals[aroundvertex[i]]->Point[0];
    v2[1] = pointnormals[aroundvertex[i2]]->Point[1] -
            pointnormals[aroundvertex[i]]->Point[1];
    v2[2] = pointnormals[aroundvertex[i2]]->Point[2] -
            pointnormals[aroundvertex[i]]->Point[2];

    CrossProduct2(v2, v1, n2);
    if (DotProduct(n2, v3) < 0.0)
      anglesdir[i] = -1;
  }

  /* Compute minimal distance       */
  mindis = 10000.0;
  nn2 = (tolnumber - 1) * (tolnumber - 2) / 2;
  for (i = 0; i < number - 2; i++) {
    if (anglesdir[i] == -1) {
      truei = index[i];
      nni2 = (tolnumber - 1 - truei) * (tolnumber - 2 - truei) / 2;
      for (j = i + 2; (j < number) && (j - i <= number - 2); j++) {
        i0 = Model(j, number);
        if (anglesdir[i0] == -1) {
          truej = index[i0];
          dis = minidis[nn2 - nni2 + truej - truei - 2];
          if (dis < mindis) {
            mindis = dis;
            minii1 = i;
            minii2 = i0;
            number1 = j - i + 1;
          }
        }
      }
    }
  }

  /* form two aroundvertex lists    */
  if (mindis < 9000.0) {
    number2 = number - number1 + 2;
    for (i = 0; i < number1; i++) {
      i1 = Model(i + minii1, number);
      aroundvertex1[i] = aroundvertex[i1];
      index1[i] = index[i1];
    }
    for (i = 0; i < number2; i++) {
      i1 = Model(i + minii2, number);
      aroundvertex2[i] = aroundvertex[i1];
      index2[i] = index[i1];
    }

    /* recursive calling              */
    MinimalPathConnection_in(aroundvertex1, ntri, trilist, number1, minidis,
                             index1, tolnumber);
    MinimalPathConnection_in(aroundvertex2, ntri, trilist, number2, minidis,
                             index2, tolnumber);
    return;
  }

  for (i = 0; i < number - 2; i++) {
    if (anglesdir[i] == -1) {
      truei = index[i];
      nni2 = (tolnumber - 1 - truei) * (tolnumber - 2 - truei) / 2;
      for (j = i + 2; (j < number) && (j - i <= number - 2); j++) {
        i0 = Model(j, number);
        truej = index[i0];
        dis = minidis[nn2 - nni2 + truej - truei - 2];
        if (dis < mindis) {
          mindis = dis;
          minii1 = i;
          minii2 = i0;
          number1 = j - i + 1;
        }
      }
    }
  }

  /* form two aroundvertex lists    */
  if (mindis < 9000.0) {
    number2 = number - number1 + 2;
    for (i = 0; i < number1; i++) {
      i1 = Model(i + minii1, number);
      aroundvertex1[i] = aroundvertex[i1];
      index1[i] = index[i1];
    }
    for (i = 0; i < number2; i++) {
      i1 = Model(i + minii2, number);
      aroundvertex2[i] = aroundvertex[i1];
      index2[i] = index[i1];
    }

    /* recursive calling              */
    MinimalPathConnection_in(aroundvertex1, ntri, trilist, number1, minidis,
                             index1, tolnumber);
    MinimalPathConnection_in(aroundvertex2, ntri, trilist, number2, minidis,
                             index2, tolnumber);
    return;
  }

  MinimalAngleConnection(aroundvertex, number, ntri, trilist);
}

/*-----------------------------------------------------------------------------
 * MinimalPathConnection -- Form the triangulation by connect ponts that
 *                              has minimal geodesic distance
 *---------------------------------------------------------------------------*/
static void MinimalPathConnection(v, aroundvertex, ntri, trilist, number) int v,
    aroundvertex[], *ntri, trilist[], number;
{
  int i, j, k, l, i0, i1, i2, numb, nn2, nni2, cur, index[100], over;
  float p0[3], p1[3], p2[3], p3[3], p4[3], p5[3], p[3], minidis[4851],
      unfoldpts[303], unfoldright[303];
  float work[303], normal[3];
  float s, t, u, w, pp1[3], p21[3], p24[3], p45[3], p40[3], sum;
  double n1[3], n2[3], n3[3], matrix[3][3];

  /* check if there nonconvex point */
  over = 0;
  for (i = 0; i < number; i++) {
    i1 = i - 1;
    if (i == 0)
      i1 = number - 1;
    i2 = i + 1;
    if (i == number - 1)
      i2 = 0;

    p3[0] = pointnormals[aroundvertex[i]]->Normal[0];
    p3[1] = pointnormals[aroundvertex[i]]->Normal[1];
    p3[2] = pointnormals[aroundvertex[i]]->Normal[2];

    p1[0] = pointnormals[aroundvertex[i1]]->Point[0] -
            pointnormals[aroundvertex[i]]->Point[0];
    p1[1] = pointnormals[aroundvertex[i1]]->Point[1] -
            pointnormals[aroundvertex[i]]->Point[1];
    p1[2] = pointnormals[aroundvertex[i1]]->Point[2] -
            pointnormals[aroundvertex[i]]->Point[2];

    p2[0] = pointnormals[aroundvertex[i2]]->Point[0] -
            pointnormals[aroundvertex[i]]->Point[0];
    p2[1] = pointnormals[aroundvertex[i2]]->Point[1] -
            pointnormals[aroundvertex[i]]->Point[1];
    p2[2] = pointnormals[aroundvertex[i2]]->Point[2] -
            pointnormals[aroundvertex[i]]->Point[2];

    CrossProduct2(p2, p1, p);
    if (DotProduct(p, p3) < 0.0) {
      over = 1;
      break;
    }
  }

  if (over == 0) {
    MinimalPathConnection_in(aroundvertex, ntri, trilist, number, minidis,
                             index, number);
    return;
  }

  nn2 = (number - 1) * (number - 2) / 2;
  /* compute the unfold point list  */
  for (i = 0; i < 3; i++) {
    p[i] = pointnormals[v]->Point[i];
    p1[i] = pointnormals[aroundvertex[0]]->Point[i] - p[i];
    p2[i] = pointnormals[aroundvertex[1]]->Point[i] - p[i];
  }
  CrossProduct2(p1, p2, normal);

  /* 1. move all the point to unfold*/
  for (i = 0; i <= number; i++) {
    i0 = 3 * i;
    unfoldpts[i0] = pointnormals[aroundvertex[i % number]]->Point[0];
    unfoldpts[i0 + 1] = pointnormals[aroundvertex[i % number]]->Point[1];
    unfoldpts[i0 + 2] = pointnormals[aroundvertex[i % number]]->Point[2];
  }
  /* 2. begin loop on i             */
  for (i = 2; i <= number; i++) {
    i0 = i * 3;
    i1 = i0 - 3;
    i2 = i1 - 3;
    for (j = 0; j < 3; j++) {
      p0[j] = unfoldpts[i2 + j];
      p1[j] = unfoldpts[i1 + j];
      p2[j] = unfoldpts[i0 + j];
      p21[j] = p2[j] - p1[j];
      pp1[j] = p[j] - p1[j];
    }
    w = DotProduct(pp1, pp1);
    t = DotProduct(p21, pp1) / w;

    p4[0] = p1[0] + t * pp1[0];
    p4[1] = p1[1] + t * pp1[1];
    p4[2] = p1[2] + t * pp1[2];

    p40[0] = p4[0] - p0[0];
    p40[1] = p4[1] - p0[1];
    p40[2] = p4[2] - p0[2];

    /* 3. compute p5                  */
    s = DotProduct(p40, pp1) / w + t;
    for (j = 0; j < 3; j++) {
      p5[j] = 0.5 * (p0[j] + p1[j] + s * pp1[j]);
      p24[j] = p2[j] - p4[j];
      p45[j] = p4[j] - p5[j];
    }
    /* 4. normalize the vactors       */
    s = sqrt(DotProduct(p24, p24));
    t = sqrt(DotProduct(p45, p45));
    u = sqrt(w);
    for (j = 0; j < 3; j++) {
      n1[j] = p24[j] / s;
      n2[j] = p45[j] / t;
      n3[j] = pp1[j] / u;
    }

    /* 5. compute rotate matrix       */
    RotateMatrix(n1, n2, n3, matrix);

    /* 6. rotate all the points       */
    for (k = i; k <= number; k++) {
      i0 = k + k + k;

      p2[0] = unfoldpts[i0];
      p2[1] = unfoldpts[i0 + 1];
      p2[2] = unfoldpts[i0 + 2];

      p21[0] = p2[0] - p1[0];
      p21[1] = p2[1] - p1[1];
      p21[2] = p2[2] - p1[2];

      t = DotProduct(p21, pp1) / w;

      p4[0] = p1[0] + t * pp1[0];
      p4[1] = p1[1] + t * pp1[1];
      p4[2] = p1[2] + t * pp1[2];

      p24[0] = p2[0] - p4[0];
      p24[1] = p2[1] - p4[1];
      p24[2] = p2[2] - p4[2];

      unfoldpts[i0] = p4[0] + matrix[0][0] * p24[0] + matrix[0][1] * p24[1] +
                      matrix[0][2] * p24[2];
      unfoldpts[i0 + 1] = p4[1] + matrix[1][0] * p24[0] +
                          matrix[1][1] * p24[1] + matrix[1][2] * p24[2];
      unfoldpts[i0 + 2] = p4[2] + matrix[2][0] * p24[0] +
                          matrix[2][1] * p24[1] + matrix[2][2] * p24[2];
    } /* end k loop                     */
  } /* end i loop                     */

  /* Compute minimal distance       */
  /* C1. left round                 */
  for (i = 0; i < number - 2; i++) {
    nni2 = (number - 1 - i) * (number - 2 - i) / 2;
    i0 = i + i + i;
    i1 = i0 + 3;
    for (j = 0; j < 3; j++) {
      work[j] = unfoldpts[i0 + j];
      work[j + 3] = unfoldpts[i1 + j];
      p3[j] = p[j] - work[j];
    }
    u = sqrt(DotProduct(p3, p3));
    cur = 2;
    over = 0;
    for (j = i + 2; j < number; j++) {
      i0 = j + j + j;
      i1 = cur + cur + cur;

      work[i1] = unfoldpts[i0];
      work[i1 + 1] = unfoldpts[i0 + 1];
      work[i1 + 2] = unfoldpts[i0 + 2];
      cur = cur + 1;
      /* C2. convexity test             */
      numb = 1;
      while ((numb > 0) && (cur > 2) && (over == 0)) {
        numb = 0;
        for (k = 1; k < cur - 1; k++) {
          i0 = k + k + k;
          i1 = i0 - 3;
          i2 = i0 + 3;

          p0[0] = work[i1] - work[i0];
          p0[1] = work[i1 + 1] - work[i0 + 1];
          p0[2] = work[i1 + 2] - work[i0 + 2];

          p2[0] = work[i2] - work[i0];
          p2[1] = work[i2 + 1] - work[i0 + 1];
          p2[2] = work[i2 + 2] - work[i0 + 2];

          CrossProduct2(p2, p0, p1);
          if (DotProduct(p1, normal) > 0.0) {

            /* C3. remove the point k         */
            numb = 1;
            for (l = k + 1; l < cur; l++) {
              i0 = l + l + l;
              i1 = i0 - 3;
              work[i1] = work[i0];
              work[i1 + 1] = work[i0 + 1];
              work[i1 + 2] = work[i0 + 2];
            }
            cur = cur - 1;
            break;
          }
        } /* end k loop                     */
      } /* end while(numb > 0 && cur > 2) */

      /* C4.compute minimal dis for(i,j)*/
      i1 = 3 * (cur - 1);
      p1[0] = work[i1] - work[0];
      p1[1] = work[i1 + 1] - work[1];
      p1[2] = work[i1 + 2] - work[2];
      CrossProduct2(p1, p3, p0);
      if (DotProduct(p0, normal) < 0.0) {
        over = 1;
        p1[0] = work[i1] - p[0];
        p1[1] = work[i1 + 1] - p[1];
        p1[2] = work[i1 + 2] - p[2];
        sum = u + sqrt(DotProduct(p1, p1));
      } else {
        sum = 0.0;
        for (l = 0; l < cur - 1; l++) {
          i0 = l + l + l;
          i1 = i0 + 3;
          p1[0] = work[i1] - work[i0];
          p1[1] = work[i1 + 1] - work[i0 + 1];
          p1[2] = work[i1 + 2] - work[i0 + 2];
          sum = sum + sqrt(DotProduct(p1, p1));
        }
      }
      minidis[nn2 - nni2 + j - i - 2] = sum;
    } /* end j loop                     */
  } /* end i loop                     */

  /* Rotate the points on place     */
  /* R1. form the ratate matrix     */
  p1[0] = unfoldpts[0] - p[0];
  p1[1] = unfoldpts[1] - p[1];
  p1[2] = unfoldpts[2] - p[2];
  u = sqrt(DotProduct(p1, p1));
  n1[0] = p1[0] / u;
  n1[1] = p1[1] / u;
  n1[2] = p1[2] / u;

  i0 = 3 * number;
  p2[0] = unfoldpts[i0] - p[0];
  p2[1] = unfoldpts[i0 + 1] - p[1];
  p2[2] = unfoldpts[i0 + 2] - p[2];
  u = sqrt(DotProduct(p2, p2));
  n2[0] = p2[0] / u;
  n2[1] = p2[1] / u;
  n2[2] = p2[2] / u;
  w = sqrt(DotProduct(normal, normal));
  n3[0] = normal[0] / w;
  n3[1] = normal[1] / w;
  n3[2] = normal[2] / w;

  RotateMatrix(n2, n1, n3, matrix);

  /* R2. rotate all the points      */
  for (i = 0; i <= number; i++) {
    i0 = 3 * i;
    i1 = i0 + 1;
    i2 = i0 + 2;
    unfoldright[i0] = matrix[0][0] * (unfoldpts[i0] - p[0]) +
                      matrix[0][1] * (unfoldpts[i1] - p[1]) +
                      matrix[0][2] * (unfoldpts[i2] - p[2]) + p[0];

    unfoldright[i1] = matrix[1][0] * (unfoldpts[i0] - p[0]) +
                      matrix[1][1] * (unfoldpts[i1] - p[1]) +
                      matrix[1][2] * (unfoldpts[i2] - p[2]) + p[1];

    unfoldright[i2] = matrix[2][0] * (unfoldpts[i0] - p[0]) +
                      matrix[2][1] * (unfoldpts[i1] - p[1]) +
                      matrix[2][2] * (unfoldpts[i2] - p[2]) + p[2];
  }
  /* Right  round distance          */
  for (i = 0; i < number - 2; i++) {
    nni2 = (number - 1 - i) * (number - 2 - i) / 2;

    /* R1. read point i, i-1, ..., 0  */
    for (j = i; j >= 0; j--) {
      i0 = 3 * j;
      i1 = 3 * (i - j);
      work[i1] = unfoldpts[i0];
      work[i1 + 1] = unfoldpts[i0 + 1];
      work[i1 + 2] = unfoldpts[i0 + 2];
    }
    p3[0] = p[0] - work[0];
    p3[1] = p[1] - work[1];
    p3[2] = p[2] - work[2];
    u = sqrt(DotProduct(p3, p3));
    cur = i + 1;
    over = 0;
    for (j = number - 1; j >= i + 2; j--) {
      i0 = j + j + j;
      i1 = cur + cur + cur;

      work[i1] = unfoldright[i0];
      work[i1 + 1] = unfoldright[i0 + 1];
      work[i1 + 2] = unfoldright[i0 + 2];

      cur = cur + 1;
      /* R2. convexity test             */
      numb = 1;
      while ((numb > 0) && (cur > 2) && (over == 0)) {
        numb = 0;
        for (k = 1; k < cur - 1; k++) {
          i0 = k + k + k;
          i1 = i0 - 3;
          i2 = i0 + 3;

          p0[0] = work[i1] - work[i0];
          p0[1] = work[i1 + 1] - work[i0 + 1];
          p0[2] = work[i1 + 2] - work[i0 + 2];

          p2[0] = work[i2] - work[i0];
          p2[1] = work[i2 + 1] - work[i0 + 1];
          p2[2] = work[i2 + 2] - work[i0 + 2];

          CrossProduct2(p0, p2, p1);
          if (DotProduct(p1, normal) > 0.0) {

            /* R3. remove the point k         */
            numb = 1;
            for (l = k + 1; l < cur; l++) {
              i0 = l + l + l;
              i1 = i0 - 3;
              work[i1] = work[i0];
              work[i1 + 1] = work[i0 + 1];
              work[i1 + 2] = work[i0 + 2];
            }
            cur = cur - 1;
            break;
          }
        } /* end k loop                     */
      } /* end while(numb > 0 && cur > 2) */

      /* R4.compute minimal dis for(i,j)*/
      i1 = 3 * (cur - 1);
      p1[0] = work[i1] - work[0];
      p1[1] = work[i1 + 1] - work[1];
      p1[2] = work[i1 + 2] - work[2];
      CrossProduct2(p3, p1, p0);
      if (DotProduct(p0, normal) < 0.0) {
        over = 1;
        p1[0] = work[i1] - p[0];
        p1[1] = work[i1 + 1] - p[1];
        p1[2] = work[i1 + 2] - p[2];
        sum = u + sqrt(DotProduct(p1, p1));
      } else {
        sum = 0.0;
        for (l = 0; l < cur - 1; l++) {
          i0 = l + l + l;
          i1 = i0 + 3;
          p1[0] = work[i1] - work[i0];
          p1[1] = work[i1 + 1] - work[i0 + 1];
          p1[2] = work[i1 + 2] - work[i0 + 2];
          sum = sum + sqrt(DotProduct(p1, p1));
        }
      }
      if (sum < minidis[nn2 - nni2 + j - i - 2])
        minidis[nn2 - nni2 + j - i - 2] = sum;
    } /* end j loop                     */
  } /* end i loop                     */

  /*
  for (i = 0; i < 6; i++)
      printf("dis = %f\n", minidis[i]);
  */

  for (i = 0; i < number; i++)
    index[i] = i;
  MinimalPathConnection_in(aroundvertex, ntri, trilist, number, minidis, index,
                           number);
}

/*-----------------------------------------------------------------------------
 * MinimalDistanceConnection -- Form the triangulation by connect ponts that
 *                              has minimal distance
 *---------------------------------------------------------------------------*/
static void MinimalDistanceConnection(aroundvertex, ntri, trilist,
                                      number) int aroundvertex[],
    *ntri, trilist[], number;
{
  int i, j, k, i0, i1, i2, anglesdir[100], number1, number2, minii1, minii2;
  int aroundvertex1[100], aroundvertex2[100];
  float mindis, dis, v1[3], v2[3], v3[3], n2[3];

  if (number == 3) {
    trilist[3 * *ntri] = aroundvertex[0];
    trilist[3 * *ntri + 1] = aroundvertex[1];
    trilist[3 * *ntri + 2] = aroundvertex[2];
    *ntri = *ntri + 1;
    return;
  }

  /* Compute angles at each vertex  */
  for (i = 0; i < number; i++) {
    anglesdir[i] = 0;
    i1 = i - 1;
    if (i == 0)
      i1 = number - 1;
    i2 = i + 1;
    if (i == number - 1)
      i2 = 0;
    for (j = 0; j < 3; j++) {
      v3[j] = pointnormals[aroundvertex[i]]->Normal[j];
      v1[j] = pointnormals[aroundvertex[i1]]->Point[j] -
              pointnormals[aroundvertex[i]]->Point[j];
      v2[j] = pointnormals[aroundvertex[i2]]->Point[j] -
              pointnormals[aroundvertex[i]]->Point[j];
    }
    CrossProduct2(v2, v1, n2);
    if (DotProduct(n2, v3) < 0.0)
      anglesdir[i] = -1;
  }

  /* Compute minimal distance       */
  mindis = 10000.0;
  for (i = 0; i < number - 2; i++) {
    if (anglesdir[i] == -1) {
      for (j = i + 2; (j < number) && (j - i <= number - 2); j++) {
        i0 = Model(j, number);
        if (anglesdir[i0] == -1) {
          for (k = 0; k < 3; k++) {
            v1[k] = pointnormals[aroundvertex[i]]->Point[k] -
                    pointnormals[aroundvertex[i0]]->Point[k];
          }
          dis = DotProduct(v1, v1);
          if (dis < mindis) {
            mindis = dis;
            minii1 = i;
            minii2 = i0;
            number1 = j - i + 1;
          }
        }
      }
    }
  }

  /* form two aroundvertex lists    */
  if (mindis < 9000.0) {
    number2 = number - number1 + 2;
    for (i = 0; i < number1; i++) {
      i1 = Model(i + minii1, number);
      aroundvertex1[i] = aroundvertex[i1];
    }
    for (i = 0; i < number2; i++) {
      i1 = Model(i + minii2, number);
      aroundvertex2[i] = aroundvertex[i1];
    }

    /* recursive calling              */
    MinimalDistanceConnection(aroundvertex1, ntri, trilist, number1);
    MinimalDistanceConnection(aroundvertex2, ntri, trilist, number2);
    return;
  }

  for (i = 0; i < number - 2; i++) {
    if (anglesdir[i] == -1) {
      for (j = i + 2; (j < number) && (j - i <= number - 2); j++) {
        i0 = Model(j, number);
        for (k = 0; k < 3; k++) {
          v1[k] = pointnormals[aroundvertex[i]]->Point[k] -
                  pointnormals[aroundvertex[i0]]->Point[k];
        }
        dis = DotProduct(v1, v1);
        if (dis < mindis) {
          mindis = dis;
          minii1 = i;
          minii2 = i0;
          number1 = j - i + 1;
        }
      }
    }
  }
  /* form two aroundvertex lists    */
  if (mindis < 9000.0) {
    number2 = number - number1 + 2;
    for (i = 0; i < number1; i++) {
      i1 = Model(i + minii1, number);
      aroundvertex1[i] = aroundvertex[i1];
    }
    for (i = 0; i < number2; i++) {
      i1 = Model(i + minii2, number);
      aroundvertex2[i] = aroundvertex[i1];
    }

    /* recursive calling              */
    MinimalDistanceConnection(aroundvertex1, ntri, trilist, number1);
    MinimalDistanceConnection(aroundvertex2, ntri, trilist, number2);
    return;
  }

  MinimalAngleConnection(aroundvertex, number, ntri, trilist);
  /*
  for (i = 0; i < number; i++) {
      for (j = i + 2; j < i + number - 1; j++) {
          i0 = Model(j, number);
          for (k = 0; k < 3; k++) {
              v1[k] = pointnormals[aroundvertex[i ]]->Point[k];
                    - pointnormals[aroundvertex[i0]]->Point[k];
          }
          dis = DotProduct(v1,v1);
          if (dis < mindis) {
             mindis = dis;
             minii1 = i;
             minii2 = i0;
             number1 = j - i + 1;
          }
      }
  }
  */

  /* form two aroundvertex lists    */
  /*
  number2 = number - number1 + 2;
  for (i = 0; i < number1; i++) {
      i1 = Model(i + minii1, number);
      aroundvertex1[i] = aroundvertex[i1];
  }
  for (i = 0; i < number2; i++) {
      i1 = Model(i + minii2, number);
      aroundvertex2[i] = aroundvertex[i1];
  }
  */
  /* recursive calling              */
  /*
  MinimalDistanceConnection(aroundvertex1,ntri,trilist,number1);
  MinimalDistanceConnection(aroundvertex2,ntri,trilist,number2);
  return;
  */
}

/*-----------------------------------------------------------------------------
 * MinimalAngleConnection -- Form the triangulation by connect ponts that
 *                           has minimal angle
 *---------------------------------------------------------------------------*/
static void MinimalAngleConnection(aroundvertex, curnum, ntri,
                                   trilist) int aroundvertex[],
    curnum, *ntri, trilist[];
{
  int i, j, maxi, minii;
  float miniangle, maxangle, angle, v1[3], v2[3], v3[3], n2[3];

  while (curnum > 3) {
    if (curnum == 4) {
      /* find the one with maximal angle*/
      maxangle = 0.0;
      miniangle = 6.2831853;
      for (i = 0; i < curnum; i++) {
        for (j = 0; j < 3; j++) {
          v3[j] = pointnormals[aroundvertex[i]]->Normal[j];
          if (i == 0) {
            v1[j] = pointnormals[aroundvertex[curnum - 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          } else {
            v1[j] = pointnormals[aroundvertex[i - 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          }
          if (i == curnum - 1) {
            v2[j] = pointnormals[aroundvertex[0]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          } else {
            v2[j] = pointnormals[aroundvertex[i + 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          }
        }
        angle = acos(DotProduct(v1, v2) /
                     sqrt(DotProduct(v1, v1) * DotProduct(v2, v2)));
        CrossProduct2(v2, v1, n2);
        if (DotProduct(n2, v3) < 0.0) {
          angle = 6.2831853 - angle;
        }
        if (maxangle < angle) {
          maxangle = angle;
          maxi = i;
        }
        if (miniangle > angle) {
          miniangle = angle;
          minii = i;
        }

      } /* end i loop                     */
      if (maxangle >= 3.141) {
        if (maxi > 0)
          minii = maxi - 1;
        if (maxi == 0)
          minii = maxi + 1;
      }
    } /* end if curnum == 4             */

    /* find the one with minimal angle*/

    if (curnum > 4) {

      /*
            if (curnum > 3) {
      */
      miniangle = 6.2831853;
      for (i = 0; i < curnum; i++) {
        for (j = 0; j < 3; j++) {
          v3[j] = pointnormals[aroundvertex[i]]->Normal[j];
          if (i == 0) {
            v1[j] = pointnormals[aroundvertex[curnum - 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          } else {
            v1[j] = pointnormals[aroundvertex[i - 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          }
          if (i == curnum - 1) {
            v2[j] = pointnormals[aroundvertex[0]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          } else {
            v2[j] = pointnormals[aroundvertex[i + 1]]->Point[j] -
                    pointnormals[aroundvertex[i]]->Point[j];
          }
        }
        angle = acos(DotProduct(v1, v2) /
                     sqrt(DotProduct(v1, v1) * DotProduct(v2, v2)));
        CrossProduct2(v2, v1, n2);
        if (DotProduct(n2, v3) < 0.0) {
          angle = 6.2831853 - angle;
        }
        if (miniangle > angle) {
          miniangle = angle;
          minii = i;
        }
      } /* end i loop                     */
    } /* end if curnum > 4              */

    /* form triangle at minimal angle */

    if (minii == 0) {
      trilist[3 * *ntri] = aroundvertex[curnum - 1];
      trilist[3 * *ntri + 1] = aroundvertex[minii];
      trilist[3 * *ntri + 2] = aroundvertex[minii + 1];
    }
    if (minii == curnum - 1) {
      trilist[3 * *ntri] = aroundvertex[minii - 1];
      trilist[3 * *ntri + 1] = aroundvertex[minii];
      trilist[3 * *ntri + 2] = aroundvertex[0];
    }
    if (minii > 0 && minii < curnum - 1) {
      trilist[3 * *ntri] = aroundvertex[minii - 1];
      trilist[3 * *ntri + 1] = aroundvertex[minii];
      trilist[3 * *ntri + 2] = aroundvertex[minii + 1];
    }
    *ntri = *ntri + 1;

    /* adjust the around list         */
    for (i = minii + 1; i < curnum; i++) {
      aroundvertex[i - 1] = aroundvertex[i];
    }
    curnum = curnum - 1;
  } /* end while                      */

  /* the final triangle             */
  for (j = 0; j < 3; j++) {
    v3[j] = pointnormals[aroundvertex[0]]->Normal[j];
    v1[j] = pointnormals[aroundvertex[2]]->Point[j] -
            pointnormals[aroundvertex[0]]->Point[j];
    v2[j] = pointnormals[aroundvertex[1]]->Point[j] -
            pointnormals[aroundvertex[0]]->Point[j];
  }
  angle =
      acos(DotProduct(v1, v2) / sqrt(DotProduct(v1, v1) * DotProduct(v2, v2)));
  CrossProduct2(v2, v1, n2);
  /*
     if (DotProduct(n2,v3) > 0.0) {
  */

  trilist[3 * *ntri] = aroundvertex[0];
  trilist[3 * *ntri + 1] = aroundvertex[1];
  trilist[3 * *ntri + 2] = aroundvertex[2];
  *ntri = *ntri + 1;
  /*
     }
  */
}

/* -----------------------------------------------------------------------------*/
/* MiniAngleAroundVertex -- Compute the mnimal angle for all the  triangles */
/*                          around a vertex */
/* -----------------------------------------------------------------------------*/
static float MiniAngleAroundVertex(v)
int v;
{
  int i, j, k, ii, jj, kk;
  float v1[3], v2[3], v3[3], smallangle, angle;

  smallangle = -1.0;
  for (i = 0; i < facevalence[v]; i++) {
    k = faceneighbor[v]->NB[i];
    ii = facecopy[k]->Index[0];
    jj = facecopy[k]->Index[1];
    kk = facecopy[k]->Index[2];
    for (j = 0; j < 3; j++) {
      v1[j] = pointnormals[jj]->Point[j] - pointnormals[ii]->Point[j];
      v2[j] = pointnormals[kk]->Point[j] - pointnormals[ii]->Point[j];
      v3[j] = pointnormals[kk]->Point[j] - pointnormals[jj]->Point[j];
    }
    angle = DotProduct(v1, v2) / sqrt(DotProduct(v1, v1) * DotProduct(v2, v2));
    if (smallangle < angle)
      smallangle = angle;

    angle = -DotProduct(v1, v3) / sqrt(DotProduct(v1, v1) * DotProduct(v3, v3));
    if (smallangle < angle)
      smallangle = angle;

    angle = DotProduct(v2, v3) / sqrt(DotProduct(v2, v2) * DotProduct(v3, v3));
    if (smallangle < angle)
      smallangle = angle;
  }
  return (smallangle);
}

/* -----------------------------------------------------------------------------*/
/* LambdaPolyForVoume -- Compute the polynomials that computing b1,b2,b3,lambda
 */
/*                       for given point in the volume D_{ijk} */
/* -----------------------------------------------------------------------------*/
void LambdaPolyForVoume(v1, v2, v3, n1, n2, n3, px, py, pz, q) float v1[], v2[],
    v3[], n1[], n2[], n3[], px[], py[], pz[], q[];
{
  float v13[3], v23[3], n13[3], n23[3];
  int i;

  for (i = 0; i < 3; i++) {
    v13[i] = v1[i] - v3[i];
    v23[i] = v2[i] - v3[i];
    n13[i] = n1[i] - n3[i];
    n23[i] = n2[i] - n3[i];
  }
  px[0] = v13[1] * v23[2] - v13[2] * v23[1];
  px[1] = v13[1] * n23[2] + v23[2] * n13[1] - v13[2] * n23[1] - v23[1] * n13[2];
  px[2] = n13[1] * n23[2] - n13[2] * n23[1];

  py[0] = -(v13[0] * v23[2] - v13[2] * v23[0]);
  py[1] =
      -(v13[0] * n23[2] + v23[2] * n13[0] - v13[2] * n23[0] - v23[0] * n13[2]);
  py[2] = -(n13[0] * n23[2] - n13[2] * n23[0]);

  pz[0] = v13[0] * v23[1] - v13[1] * v23[0];
  pz[1] = v13[0] * n23[1] + v23[1] * n13[0] - v13[1] * n23[0] - v23[0] * n13[1];
  pz[2] = n13[0] * n23[1] - n13[1] * n23[0];

  q[0] = Determinant(v1, v2, v3);
  q[1] = Determinant(n1, v2, v3) + Determinant(v1, n2, v3) +
         Determinant(v1, v2, n3);
  q[2] = Determinant(v1, n2, n3) + Determinant(n1, v2, n3) +
         Determinant(n1, n2, v3);
  q[3] = Determinant(n1, n2, n3);
}

/* -----------------------------------------------------------------------------*/
/* GlobalToLocal -- Compute the local corrdinate b1,b2,b3,lambda for given
 * point*/
/*              in R^3, it return 0 if the point does not in the volume D_{ijk}
 */
/* -----------------------------------------------------------------------------*/
int GlobalToLocal(v1, v2, v3, n1, n2, n3, p, px, py, pz, q, bb, lambda)
float v1[], v2[], v3[], n1[], n2[], n3[];
float p[], px[], py[], pz[], q[], bb[], *lambda;
{
  double poly[4], root[3], matrix[9], minroot;
  int i, j, k, number, info;

  for (i = 0; i < 4; i++) {
    poly[i] = -q[i];
  }
  for (i = 0; i < 3; i++) {
    poly[i] = poly[i] + p[0] * px[i] + p[1] * py[i] + p[2] * pz[i];
  }
  FindRoots3(poly[3], poly[2], poly[1], poly[0], root, &number);

  minroot = 10000.0;
  for (i = 0; i < number; i++) {
    if (fabs(minroot) > fabs(root[i]))
      minroot = root[i];
  }
  /* find b1,b2,b3    */
  for (j = 0; j < 3; j++) {
    matrix[3 * j] = v1[j] + minroot * n1[j];
    matrix[3 * j + 1] = v2[j] + minroot * n2[j];
    matrix[3 * j + 2] = v3[j] + minroot * n3[j];
  }
  gaussinverse3(matrix, 0.0000000000001, &info);
  if (info != 0) {
    for (j = 0; j < 3; j++) {
      bb[j] = 0;
      for (k = 0; k < 3; k++) {
        bb[j] = bb[j] + matrix[j * 3 + k] * p[k];
      }
    }
    if (bb[0] >= 0.00000 && bb[1] >= 0.00000 && bb[2] >= 0.00000) {
      *lambda = minroot;
      return (1);
    }
  }
  return (0);
}

void generateNormal(int item_no) {
  ProduceNormalS(pointnormals, facetetras, numbtris);

  facevalence = (int *)malloc(numbpts * sizeof(int));
  edgevalence = (int *)malloc(numbpts * sizeof(int));
  faceneighbor = (Neighbor *)malloc(numbpts * sizeof(Neighbor));
  edgeneighbor = (Neighbor *)malloc(numbpts * sizeof(Neighbor));

  Valence_Neighbor(facetetras, numbpts, numbtris, faceneighbor, edgeneighbor,
                   facevalence, edgevalence);
  AdjacentTriangles(facetetras, numbtris, facevalence, faceneighbor);
}
