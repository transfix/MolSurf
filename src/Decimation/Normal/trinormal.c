#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <Decimation/Bezier/BezierMath.h>
#include <Decimation/Misc/misc.h>
#include <Decimation/geometry.h>

#include <Decimation/Normal/normal.h>

static void RemoveOneVertexWork();
static void ReviseAllLists();
static void Face_Gradient();

/*-----------------------------------------------------------------------------
 * NewOreintationS-- Orientate all the face normals
 *---------------------------------------------------------------------------*/
void NewOreintationS(pointnormals, numbpts, facetetras, numbtris, facevalence,
                     edgevalence, faceneighbor,
                     edgeneighbor) VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
{
  int i, j, k, ii, jj, kk, ll, iii, jjj, kkk, lll;
  int u, v, w, x, count1;
  float pu[3], pv[3], pw[3], px[3], normal[3], normal1[3], sum;
  int numwell, numwell1, nworking, curentvertex, curentvalence, i0, j0;
  int *workvalence, *workindex, *inworkindex, *newold, nworking1;
  char buf[50];
  /*
  Normal   *normals;
  Vertex   *vertices;
  GrObject *tripatch;
  */

  workvalence = (int *)malloc(numbpts * sizeof(int));
  workindex = (int *)malloc(numbpts * sizeof(int));
  inworkindex = (int *)malloc(numbpts * sizeof(int));
  newold = (int *)malloc(numbpts * sizeof(int));

  /* set to zero                 */
  for (i = 0; i < numbpts; i++) {
    workvalence[i] = 0;
    workindex[i] = 0;
    inworkindex[i] = 0;
    newold[i] = 0;
  }

  /* initilization                */
  for (j = 0; j < numbtris; j++) {
    facetetras[j]->Orien = 0;
  }

  lll = 0;
  for (ll = 0; ll < numbtris; ll++) /* these two line are added for
                                       handling multiple sheet surface */
    if (facetetras[ll]->Orien == 0) {
      facetetras[ll]->Orien = 1;
      lll = lll + 1;

      /* Patching  */
      /*
      tripatch = StartGrObjectDescription();
      DrawColor(cagd_colors[(current_color+multi)%NUM_COLORS]);
      DrawShade((current_shade + multi - 1)%NUM_SHADES);

         ii = facetetras[ll]->Index[0];
         jj = facetetras[ll]->Index[1];
         kk = facetetras[ll]->Index[2];

         vertices = (Vertex *) malloc (3 * sizeof (Vertex));
         normals = (Normal *) malloc (3 * sizeof (Normal));
         for (j = 0; j < 3; j++) {
            vertices[0][j] = pointnormals[ii]->Point[j];
            normals[0][j] = pointnormals[ii]->Normal[j];

            vertices[1][j] = pointnormals[jj]->Point[j];
            normals[1][j] = pointnormals[jj]->Normal[j];

            vertices[2][j] = pointnormals[kk]->Point[j];
            normals[2][j] = pointnormals[kk]->Normal[j];
         }
         DrawPolygon(3, vertices, normals);


      */

      nworking = 0;
      for (i = 0; i < 3; i++) {
        ii = facetetras[ll]->Index[i];
        if (facevalence[ii] > 1) {
          workindex[nworking] = ii;
          workvalence[nworking] = 1;
          inworkindex[ii] = nworking;
          newold[ii] = 1; /* ==1 -- old    */
          nworking = nworking + 1;
        }
      }

      printf("ll = %d\n", ll);
      nworking1 = nworking;
      count1 = 0;
      /* out  loop                    */
      while (nworking > 0) {
        curentvertex = workindex[nworking - 1]; /* the last one  */
        curentvalence = facevalence[curentvertex];

        /* how many has done            */
        numwell = 0;
        numwell1 = -1;
        for (i = 0; i < curentvalence; i++) {
          if (facetetras[faceneighbor[curentvertex]->NB[i]]->Orien == 1)
            numwell = numwell + 1;
        }

        /* around one vertex loop       */
        while (numwell < curentvalence && numwell1 < numwell) { /* else done */
          numwell1 = numwell;
          for (i = 0; i < curentvalence; i++) {
            i0 = faceneighbor[curentvertex]->NB[i];
            for (j = 0; j < curentvalence; j++) {
              j0 = faceneighbor[curentvertex]->NB[j];
              if ((facetetras[i0]->Orien == 1) && (i0 != j0) &&
                  (facetetras[j0]->Orien == 0)) {
                ii = facetetras[i0]->Index[0];
                jj = facetetras[i0]->Index[1];
                kk = facetetras[i0]->Index[2];

                iii = facetetras[j0]->Index[0];
                jjj = facetetras[j0]->Index[1];
                kkk = facetetras[j0]->Index[2];
                w = -1;
                if (ii == iii && jj == jjj) {
                  u = ii;
                  v = jj;
                  w = kkk;
                  x = kk;
                }
                if (ii == jjj && jj == kkk) {
                  u = ii;
                  v = jj;
                  w = iii;
                  x = kk;
                }
                if (ii == iii && jj == kkk) {
                  u = ii;
                  v = jj;
                  w = jjj;
                  x = kk;
                }

                if (jj == iii && kk == jjj) {
                  u = jj;
                  v = kk;
                  w = kkk;
                  x = ii;
                }
                if (jj == jjj && kk == kkk) {
                  u = jj;
                  v = kk;
                  w = iii;
                  x = ii;
                }
                if (jj == iii && kk == kkk) {
                  u = jj;
                  v = kk;
                  w = jjj;
                  x = ii;
                }

                if (ii == iii && kk == jjj) {
                  u = ii;
                  v = kk;
                  w = kkk;
                  x = jj;
                }
                if (ii == jjj && kk == kkk) {
                  u = ii;
                  v = kk;
                  w = iii;
                  x = jj;
                }
                if (ii == iii && kk == kkk) {
                  u = ii;
                  v = kk;
                  w = jjj;
                  x = jj;
                }
                if (w >= 0) {
                  facetetras[j0]->Orien = 1;
                  numwell = numwell + 1;
                  for (k = 0; k < 3; k++) {
                    pu[k] = pointnormals[u]->Point[k];
                    pv[k] = pointnormals[v]->Point[k];
                    pw[k] = pointnormals[w]->Point[k];
                    px[k] = pointnormals[x]->Point[k];
                  }
                  TriangleNormal(pv, pu, pw, normal);
                  TriangleNormal(pu, pv, px, normal1);
                  sum = 0.0;
                  for (k = 0; k < 3; k++) {
                    sum = sum + normal1[k] * facetetras[i0]->Normal[k];
                  }
                  if (sum < 0.0) { /* change direction        */
                    for (k = 0; k < 3; k++) {
                      normal[k] = -normal[k];
                    }
                  }
                  for (k = 0; k < 3; k++) { /* sign the correct normal */
                    facetetras[j0]->Normal[k] = normal[k];
                  }

                  /* revise all the lists       */
                  ReviseAllLists(&nworking, iii, jjj, kkk, workvalence,
                                 workindex, inworkindex, newold, facevalence);
                }
              }
            } /* end j loop    */
          } /* end i loop    */
        } /* end while     */
        if (numwell < curentvalence) {
          sprintf(buf, "The suarface is not manifold at vertex %d",
                  curentvertex);
          // AppendOutputLine(buf);
          RemoveOneVertexWork(&nworking, curentvertex, workvalence, workindex,
                              inworkindex);
        }
        if (nworking1 == nworking) {
          count1 = count1 + 1;
          if (count1 > 2000) {
            printf(" count > 2000, curentvalence = %d,%d\n", curentvalence, ll);
            goto L1;
          }
        } else {
          nworking1 = nworking;
          count1 = 0;
        }
      }
    L1:;

      /*
      EndGrObjectDescription();
      UIAddObject(tripatch, gtImplicitPlaneCurve, 0, NULL, auto_select);
      */

    } /* newly added */
  printf("This surface has %d components2222\n", lll);

  free(workvalence);
  free(workindex);
  free(inworkindex);
  free(newold);
}

/*-----------------------------------------------------------------------------
 * NunmberCompnents-- Orientate all the face normals
 *---------------------------------------------------------------------------*/
int NunmberCompnents(pointnormals, numbpts, facetetras, numbtris, facevalence,
                     edgevalence, faceneighbor, edgeneighbor, rearangeorder,
                     grouptou)
VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
int *rearangeorder, grouptou[];
{
  int i, j, ii, jj, kk, ll, iii, jjj, kkk, lll;
  int u, v, w, x;
  int numwell, numwell1, nworking, curentvertex, curentvalence, i0, j0, count;
  int *workvalence, *workindex, *inworkindex, *newold;
  char buf[50];

  workvalence = (int *)malloc(numbpts * sizeof(int));
  workindex = (int *)malloc(numbpts * sizeof(int));
  inworkindex = (int *)malloc(numbpts * sizeof(int));
  newold = (int *)malloc(numbpts * sizeof(int));

  /* set to zero                 */
  for (i = 0; i < numbpts; i++) {
    workvalence[i] = 0;
    workindex[i] = 0;
    inworkindex[i] = 0;
    newold[i] = 0;
  }

  /* initilization                */
  for (j = 0; j < numbtris; j++) {
    facetetras[j]->Orien = 0;
  }

  lll = 0;
  count = 0;
  for (ll = 0; ll < numbtris; ll++) /* these two line are added for
                                       handling multiple sheet surface */
    if (facetetras[ll]->Orien == 0) {
      facetetras[ll]->Orien = 1;
      rearangeorder[count] = ll;
      count = count + 1;
      grouptou[lll] = ll;
      lll = lll + 1;

      nworking = 0;
      for (i = 0; i < 3; i++) {
        ii = facetetras[ll]->Index[i];
        if (facevalence[ii] > 1) {
          workindex[nworking] = ii;
          workvalence[nworking] = 1;
          inworkindex[ii] = nworking;
          newold[ii] = 1; /* ==1 -- old    */
          nworking = nworking + 1;
        }
      }

      /* out  loop                    */
      while (nworking > 0) {
        curentvertex = workindex[nworking - 1]; /* the last one  */
        curentvalence = facevalence[curentvertex];

        /* how many has done            */
        numwell = 0;
        numwell1 = -1;
        for (i = 0; i < curentvalence; i++) {
          if (facetetras[faceneighbor[curentvertex]->NB[i]]->Orien == 1)
            numwell = numwell + 1;
        }

        /* around one vertex loop       */
        while (numwell < curentvalence && numwell1 < numwell) { /* else done */
          numwell1 = numwell;
          for (i = 0; i < curentvalence; i++) {
            i0 = faceneighbor[curentvertex]->NB[i];
            for (j = 0; j < curentvalence; j++) {
              j0 = faceneighbor[curentvertex]->NB[j];
              if ((facetetras[i0]->Orien == 1) && (i0 != j0) &&
                  (facetetras[j0]->Orien == 0)) {
                ii = facetetras[i0]->Index[0];
                jj = facetetras[i0]->Index[1];
                kk = facetetras[i0]->Index[2];

                iii = facetetras[j0]->Index[0];
                jjj = facetetras[j0]->Index[1];
                kkk = facetetras[j0]->Index[2];
                w = -1;
                if (ii == iii && jj == jjj) {
                  u = ii;
                  v = jj;
                  w = kkk;
                  x = kk;
                }
                if (ii == jjj && jj == kkk) {
                  u = ii;
                  v = jj;
                  w = iii;
                  x = kk;
                }
                if (ii == iii && jj == kkk) {
                  u = ii;
                  v = jj;
                  w = jjj;
                  x = kk;
                }

                if (jj == iii && kk == jjj) {
                  u = jj;
                  v = kk;
                  w = kkk;
                  x = ii;
                }
                if (jj == jjj && kk == kkk) {
                  u = jj;
                  v = kk;
                  w = iii;
                  x = ii;
                }
                if (jj == iii && kk == kkk) {
                  u = jj;
                  v = kk;
                  w = jjj;
                  x = ii;
                }

                if (ii == iii && kk == jjj) {
                  u = ii;
                  v = kk;
                  w = kkk;
                  x = jj;
                }
                if (ii == jjj && kk == kkk) {
                  u = ii;
                  v = kk;
                  w = iii;
                  x = jj;
                }
                if (ii == iii && kk == kkk) {
                  u = ii;
                  v = kk;
                  w = jjj;
                  x = jj;
                }
                if (w >= 0) {
                  facetetras[j0]->Orien = 1;

                  rearangeorder[count] = j0;
                  count = count + 1;

                  numwell = numwell + 1;
                  /* revise all the lists       */
                  ReviseAllLists(&nworking, iii, jjj, kkk, workvalence,
                                 workindex, inworkindex, newold, facevalence);
                }
              }
            } /* end j loop    */
          } /* end i loop    */
        } /* end while     */
        if (numwell < curentvalence) {
          sprintf(buf, "The suarface is not manifold at vertex %d",
                  curentvertex);
          // AppendOutputLine(buf);
          RemoveOneVertexWork(&nworking, curentvertex, workvalence, workindex,
                              inworkindex);
        }
      }
    } /* newly added */
  printf("This surface has %d components1111\n", lll);
  grouptou[lll] = numbtris;

  free(workvalence);
  free(workindex);
  free(inworkindex);
  free(newold);

  return (lll);
}
/*-----------------------------------------------------------------------------
 * NormalByAveraging-- Compute vertex normal by averaging face normals
 *---------------------------------------------------------------------------*/
void NormalByAveraging(pointnormals, numbpts, facetetras, numbtris, facevalence,
                       edgevalence, faceneighbor,
                       edgeneighbor) VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
{
  int i, j, k, jj;
  float sum, normal[3];
  /* compute normal at vertices */
  /* by averging the face normal*/
  for (i = 0; i < numbpts; i++) {
    if (facevalence[i] > 0) {
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 0.0;
      for (j = 0; j < facevalence[i]; j++) {
        jj = faceneighbor[i]->NB[j];
        for (k = 0; k < 3; k++) {
          normal[k] = normal[k] + facetetras[jj]->Normal[k];
        }
      }
      sum = 0.0;
      for (k = 0; k < 3; k++) {
        sum = sum + normal[k] * normal[k];
      }
      sum = sqrt(sum);
      for (k = 0; k < 3; k++) {
        pointnormals[i]->Normal[k] = normal[k] / sum;
      }
    } else {
      printf("The point %d is a isolated point\n", i);
    }
  }
}

/*-----------------------------------------------------------------------------
 * GradientByAveraging-- Compute vertex gradient by averaging face normals
 *---------------------------------------------------------------------------*/
void GradientByAveraging(pointnormals, numbpts, facetetras, numbtris,
                         facevalence, edgevalence, faceneighbor, edgeneighbor,
                         fungrad) VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
float *fungrad;
{
  int i, j, k, i1, i2, i3;
  float normal[3], p1[3], p2[3], p3[3], grad[3], f1, f2, f3;

  /* compute gradient at vertices */
  /* by averging the face gradient*/
  for (i = 0; i < numbtris; i++) {
    i1 = facetetras[i]->Index[0];
    i2 = facetetras[i]->Index[1];
    i3 = facetetras[i]->Index[2];
    for (j = 0; j < 3; j++) {
      p1[j] = pointnormals[i1]->Point[j];
      p2[j] = pointnormals[i2]->Point[j];
      p3[j] = pointnormals[i3]->Point[j];
      normal[j] = facetetras[i]->Normal[j];
    }
    f1 = fungrad[4 * i1];
    f2 = fungrad[4 * i2];
    f3 = fungrad[4 * i3];
    Face_Gradient(p1, p2, p3, normal, f1, f2, f3, grad);

    fungrad[4 * i1 + 1] = fungrad[4 * i1 + 1] + grad[0];
    fungrad[4 * i1 + 2] = fungrad[4 * i1 + 2] + grad[1];
    fungrad[4 * i1 + 3] = fungrad[4 * i1 + 3] + grad[2];

    fungrad[4 * i2 + 1] = fungrad[4 * i2 + 1] + grad[0];
    fungrad[4 * i2 + 2] = fungrad[4 * i2 + 2] + grad[1];
    fungrad[4 * i2 + 3] = fungrad[4 * i2 + 3] + grad[2];

    fungrad[4 * i3 + 1] = fungrad[4 * i3 + 1] + grad[0];
    fungrad[4 * i3 + 2] = fungrad[4 * i3 + 2] + grad[1];
    fungrad[4 * i3 + 3] = fungrad[4 * i3 + 3] + grad[2];
  }

  /* averaging                       */
  for (i = 0; i < numbpts; i++) {
    if (facevalence[i] > 0) {
      for (k = 1; k <= 3; k++) {
        fungrad[4 * i + k] = fungrad[4 * i + k] / facevalence[i];
      }
    } else {
      printf("The point %d is a isolated point\n", i);
    }
  }
}

void Face_Gradient(p1, p2, p3, n, f1, f2, f3, grad) float p1[], p2[], p3[], n[],
    f1, f2, f3, grad[];
{
  int i, j;
  float left[3];
  double matrix[9];

  left[0] = f2 - f1;
  left[1] = f3 - f1;
  left[2] = 0.0;

  matrix[0] = p2[0] - p1[0];
  matrix[1] = p2[1] - p1[1];
  matrix[2] = p2[2] - p1[2];

  matrix[3] = p3[0] - p1[0];
  matrix[4] = p3[1] - p1[1];
  matrix[5] = p3[2] - p1[2];

  matrix[6] = n[0];
  matrix[7] = n[1];
  matrix[8] = n[2];

  Normalized_GaussInverse(matrix, 3, 0.00000001);

  for (i = 0; i < 3; i++) {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++) {
      grad[i] = grad[i] + matrix[i * 3 + j] * left[j];
    }
  }
}

/*-----------------------------------------------------------------------------
 * NormalByLoop-- Compute vertex normal by  loop's method
 *---------------------------------------------------------------------------*/
void NormalByLoop(pointnormals, numbpts, facetetras, numbtris, facevalence,
                  edgevalence, faceneighbor,
                  edgeneighbor) VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
{
  float points[90], normal[3], sum;
  int l, i, j, k, ii, jj, kk, iii, jjj, kkk, i0, j0, curentvalence, tou, wei,
      newwei, curnum, w;

  /* first by averging               */
  NormalByAveraging(pointnormals, numbpts, facetetras, numbtris, facevalence,
                    edgevalence, faceneighbor, edgeneighbor);

  for (l = 0; l < numbpts; l++) {
    curentvalence = facevalence[l];
    if ((curentvalence > 30) || (curentvalence < 2) ||
        (facevalence[l] != edgevalence[l])) { /* begin if--1      */
      printf("facevalence[%d]  = %d, edgevalence[%d] = %d\n", l, facevalence[l],
             l, edgevalence[l]);
    } else {

      /* initialization                   */
      i0 = faceneighbor[l]->NB[0]; /* current triangle */
      ii = facetetras[i0]->Index[0];
      jj = facetetras[i0]->Index[1];
      kk = facetetras[i0]->Index[2];
      if (ii == l) {
        for (k = 0; k < 3; k++) {
          points[k] = pointnormals[jj]->Point[k];
          points[3 + k] = pointnormals[kk]->Point[k];
        }
        tou = jj;
        wei = kk;
      }
      if (jj == l) {
        for (k = 0; k < 3; k++) {
          points[k] = pointnormals[ii]->Point[k];
          points[3 + k] = pointnormals[kk]->Point[k];
        }
        tou = ii;
        wei = kk;
      }
      if (kk == l) {
        for (k = 0; k < 3; k++) {
          points[k] = pointnormals[ii]->Point[k];
          points[3 + k] = pointnormals[jj]->Point[k];
        }
        tou = ii;
        wei = jj;
      }
      curnum = 2;

      /* begin loop          */
      i = 0;
      while (wei != tou && curnum != curentvalence && i < curentvalence) {
        i = i + 1;
        for (j = 0; j < curentvalence; j++) {
          j0 = faceneighbor[l]->NB[j];
          if (i0 != j0) {
            iii = facetetras[j0]->Index[0];
            jjj = facetetras[j0]->Index[1];
            kkk = facetetras[j0]->Index[2];
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
                for (k = 0; k < 3; k++) {
                  points[curnum * 3 + k] = pointnormals[wei]->Point[k];
                }
                curnum = curnum + 1;
              }
            }
          }
        } /* end j loop  */
      } /* end while   */
      LoopNormal(points, curnum, normal);
      sum = 1.0;
      if (DotProduct(pointnormals[l]->Normal, normal) < 0.0)
        sum = -1.0;
      for (k = 0; k < 3; k++) {
        pointnormals[l]->Normal[k] = sum * normal[k];
      }
    } /* end if-1    */
  } /* end l loop  */
}

/*-----------------------------------------------------------------------------
 * NormalByFitting-- Compute vertex normal by the lease square fitting
 *                   a quadratic
 *---------------------------------------------------------------------------*/
void NormalByFitting(pointnormals, numbpts, facetetras, numbtris, facevalence,
                     edgevalence, faceneighbor, edgeneighbor, errorbd,
                     hessian) VertixNormal *pointnormals;
Face *facetetras;
int numbtris, numbpts, *facevalence, *edgevalence;
Neighbor *faceneighbor, *edgeneighbor;
float errorbd;
int hessian; /* hessian == 1, want Hessian  */
{
  float level0[3], level1[75], level2[150], normal[3], Hessian[6], sum;
  int numb1, numb2, level2list[100], ifequal, ll, info;
  int i, j, k, l, jj, kk, curentvalence;

  /* first by averging               */
  NormalByAveraging(pointnormals, numbpts, facetetras, numbtris, facevalence,
                    edgevalence, faceneighbor, edgeneighbor);

  for (l = 0; l < numbpts; l++) {
    numb2 = 0;

    /* level 0 data                             */
    level0[0] = pointnormals[l]->Point[0];
    level0[1] = pointnormals[l]->Point[1];
    level0[2] = pointnormals[l]->Point[2];

    /* level 1 data, limit to 25 points         */
    curentvalence = edgevalence[l];
    numb1 = curentvalence;
    if (numb1 > 25)
      numb1 = 25;
    for (i = 0; i < numb1; i++) {
      j = edgeneighbor[l]->NB[i];
      level1[3 * i] = pointnormals[j]->Point[0];
      level1[3 * i + 1] = pointnormals[j]->Point[1];
      level1[3 * i + 2] = pointnormals[j]->Point[2];

      /* level 2 data, limit to 50 points         */
      kk = edgevalence[j];
      for (k = 0; k < kk; k++) {
        jj = edgeneighbor[j]->NB[k];
        if (jj != l) {
          ifequal = 0;
          for (ll = 0; ll < curentvalence; ll++) {
            if (edgeneighbor[l]->NB[ll] == jj) {
              ifequal = 1;
              break;
            }
          }
          if (ifequal == 0) {
            for (ll = 0; ll < numb2; ll++) {
              if (level2list[ll] == jj) {
                ifequal = 1;
                break;
              }
            }
          }
          if (ifequal == 0 && numb2 < 50) {
            level2[3 * numb2] = pointnormals[jj]->Point[0];
            level2[3 * numb2 + 1] = pointnormals[jj]->Point[1];
            level2[3 * numb2 + 2] = pointnormals[jj]->Point[2];
            level2list[numb2] = jj;
            numb2 = numb2 + 1;
          }
        }
      } /* end k loop      */
    } /* end i loop      */

    info = LeastSqureFitNormal(level0, level1, numb1, level2, numb2, normal,
                               Hessian, errorbd);
    printf("vertex  %d, Least square Fit. info = %d\n", l, info);

    /* put the normal computed and orien  */
    sum = 1.0;
    if (DotProduct(pointnormals[l]->Normal, normal) < 0.0)
      sum = -1.0;
    for (k = 0; k < 3; k++) {
      pointnormals[l]->Normal[k] = sum * normal[k];
    }
    if (hessian == 1) {
      for (k = 0; k < 6; k++) {
        pointnormals[l]->Hessian[k] = sum * Hessian[k];
      }
    }
  } /* end l loop        */
}

/*-----------------------------------------------------------------------------
 * ReviseAllLists -- Modify all the listd
 *---------------------------------------------------------------------------*/
void ReviseAllLists(nworking, iii, jjj, kkk, workvalence, workindex,
                    inworkindex, newold, valence) int *nworking,
    iii, jjj, kkk, *workvalence, *workindex, *inworkindex, *newold, *valence;
{
  int ijk[3], i, j, k, l;

  ijk[0] = iii;
  ijk[1] = jjj;
  ijk[2] = kkk;
  for (i = 0; i < 3; i++) {
    j = ijk[i];
    if (newold[j] == 0) { /* newly added      */
      workvalence[*nworking] = 1;
      workindex[*nworking] = j;
      inworkindex[j] = *nworking;
      *nworking = *nworking + 1;
      newold[j] = 1;
    } else { /* old vertex      */
      k = inworkindex[j];
      workvalence[k] = workvalence[k] + 1;
    }
  }

  /* moving the lists                         */
  for (i = 0; i < 3; i++) {
    j = ijk[i];
    k = inworkindex[j];
    if (workvalence[k] == valence[j]) { /* need to move */
      for (l = k + 1; l < *nworking; l++) {
        inworkindex[workindex[l]] = inworkindex[workindex[l]] - 1;
        workvalence[l - 1] = workvalence[l];
        workindex[l - 1] = workindex[l];
      }
      workvalence[*nworking - 1] = 0;
      workindex[*nworking - 1] = 0;
      *nworking = *nworking - 1;
    }
  }
}

/*-----------------------------------------------------------------------------
 * RemoveOneVertexWorK -- Remove one vertex(j-th) from the working lists
 *---------------------------------------------------------------------------*/
void RemoveOneVertexWork(nworking, j, workvalence, workindex,
                         inworkindex) int *nworking,
    j, *workvalence, *workindex, *inworkindex;

{
  int k, l;
  k = inworkindex[j];
  for (l = k + 1; l < *nworking; l++) {
    inworkindex[workindex[l]] = inworkindex[workindex[l]] - 1;
    workvalence[l - 1] = workvalence[l];
    workindex[l - 1] = workindex[l];
  }
  workvalence[*nworking - 1] = 0;
  workindex[*nworking - 1] = 0;
  *nworking = *nworking - 1;
}
/*-----------------------------------------------------------------------------
 * ProduceNormalS--compute the normal of all triangles
 *---------------------------------------------------------------------------*/
void ProduceNormalS(pointnormals, facetetras,
                    numbtris) VertixNormal *pointnormals;
Face *facetetras;
int numbtris;
{
  int i, j, ii, jj, kk;
  float p1[3], p2[3], p3[3], normal[3], center[3];

  for (i = 0; i < numbtris; i++) {
    ii = facetetras[i]->Index[0];
    jj = facetetras[i]->Index[1];
    kk = facetetras[i]->Index[2];
    for (j = 0; j < 3; j++) {
      p1[j] = pointnormals[ii]->Point[j];
      p2[j] = pointnormals[jj]->Point[j];
      p3[j] = pointnormals[kk]->Point[j];
      center[j] = (p1[j] + p2[j] + p3[j]) / 3.0;
      facetetras[i]->Center[j] = center[j];
    }
    TriangleNormal(p1, p2, p3, normal);
    for (j = 0; j < 3; j++) {
      facetetras[i]->Normal[j] = normal[j];
    }
  }
}

/*-----------------------------------------------------------------------------
 * OreintationNormalS--compute the normal of all triangle,  suupose vertex
 *                     normals are given
 *---------------------------------------------------------------------------*/
void OreintationNormalS(facetetras, numbtris,
                        pointnormals) VertixNormal *pointnormals;
Face *facetetras;
int numbtris;
{
  int i, j;

  for (i = 0; i < numbtris; i++) {
    if (DotProduct(facetetras[i]->Normal,
                   pointnormals[facetetras[i]->Index[0]]->Normal) +
            DotProduct(facetetras[i]->Normal,
                       pointnormals[facetetras[i]->Index[1]]->Normal) +
            DotProduct(facetetras[i]->Normal,
                       pointnormals[facetetras[i]->Index[2]]->Normal) <
        0.0) {
      for (j = 0; j < 3; j++) {
        facetetras[i]->Normal[j] = -facetetras[i]->Normal[j];
      }
    }
  }
}
/*-----------------------------------------------------------------------------
 * GlobalOreintationS-- Orientate all the normals
 *---------------------------------------------------------------------------*/
void GlobalOreintationS(pointnormals, facetetras,
                        numbtris) VertixNormal *pointnormals;
Face *facetetras;
int numbtris;
{
  int i, j, k, ii, jj, kk, iii, jjj, kkk;
  int u, v, w, x;
  float pu[3], pv[3], pw[3], px[3], normal[3], normal1[3], sum;
  int numwell;

  for (j = 0; j < numbtris; j++) {
    facetetras[j]->Orien = 0;
  }
  facetetras[0]->Orien = 1;
  numwell = 1;

  while (numwell < numbtris) {
    for (i = 0; i < numbtris; i++) {
      for (j = 0; j < numbtris; j++) {
        if ((facetetras[i]->Orien == 1) && (i != j) &&
            (facetetras[j]->Orien == 0)) {
          ii = facetetras[i]->Index[0];
          jj = facetetras[i]->Index[1];
          kk = facetetras[i]->Index[2];

          iii = facetetras[j]->Index[0];
          jjj = facetetras[j]->Index[1];
          kkk = facetetras[j]->Index[2];
          w = -1;
          if (ii == iii && jj == jjj) {
            u = ii;
            v = jj;
            w = kkk;
            x = kk;
          }
          if (ii == jjj && jj == kkk) {
            u = ii;
            v = jj;
            w = iii;
            x = kk;
          }
          if (ii == iii && jj == kkk) {
            u = ii;
            v = jj;
            w = jjj;
            x = kk;
          }

          if (jj == iii && kk == jjj) {
            u = jj;
            v = kk;
            w = kkk;
            x = ii;
          }
          if (jj == jjj && kk == kkk) {
            u = jj;
            v = kk;
            w = iii;
            x = ii;
          }
          if (jj == iii && kk == kkk) {
            u = jj;
            v = kk;
            w = jjj;
            x = ii;
          }

          if (ii == iii && kk == jjj) {
            u = ii;
            v = kk;
            w = kkk;
            x = jj;
          }
          if (ii == jjj && kk == kkk) {
            u = ii;
            v = kk;
            w = iii;
            x = jj;
          }
          if (ii == iii && kk == kkk) {
            u = ii;
            v = kk;
            w = jjj;
            x = jj;
          }
          if (w >= 0) {
            facetetras[j]->Orien = 1;
            numwell = numwell + 1;
            for (k = 0; k < 3; k++) {
              pu[k] = pointnormals[u]->Point[k];
              pv[k] = pointnormals[v]->Point[k];
              pw[k] = pointnormals[w]->Point[k];
              px[k] = pointnormals[x]->Point[k];
            }
            TriangleNormal(pv, pu, pw, normal);
            TriangleNormal(pu, pv, px, normal1);
            sum = 0.0;
            for (k = 0; k < 3; k++) {
              sum = sum + normal1[k] * facetetras[i]->Normal[k];
            }
            if (sum < 0.0) {
              for (k = 0; k < 3; k++) {
                normal[k] = -normal[k];
              }
            }
            for (k = 0; k < 3; k++) {
              facetetras[j]->Normal[k] = normal[k];
            }
          }
        }
      } /* end j loop    */
    } /* end i loop    */
  } /* end while     */
}

/*-----------------------------------------------------------------------------
 * VertexNormalS--compute the normal of all vertices
 *---------------------------------------------------------------------------*/
void VertexNormalS(pointnormals, numbpts, facetetras,
                   numbtris) VertixNormal *pointnormals;
Face *facetetras;
int numbpts, numbtris;
{
  int i, j, k, ii, jj, kk;
  float normal[3], sum;

  for (i = 0; i < numbpts; i++) {
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    for (j = 0; j < numbtris; j++) {
      ii = facetetras[j]->Index[0];
      jj = facetetras[j]->Index[1];
      kk = facetetras[j]->Index[2];
      if (i == ii || i == jj || i == kk) {
        for (k = 0; k < 3; k++) {
          normal[k] = normal[k] + facetetras[j]->Normal[k];
        }
      }
    }
    sum = 0.0;
    for (k = 0; k < 3; k++) {
      sum = sum + normal[k] * normal[k];
    }
    sum = sqrt(sum);
    for (k = 0; k < 3; k++) {
      normal[k] = normal[k] / sum;
      pointnormals[i]->Normal[k] = normal[k];
    }
  }
}

/*-----------------------------------------------------------------------------
 * LoopNormal -- Compute normal by loop's subdivision
 *---------------------------------------------------------------------------*/
void LoopNormal(points, n, normal) float points[], normal[];
int n;
{
  float u[3], v[3], pin, c, d;
  int i, j, k;

  pin = 3.1415926 / n;
  pin = pin + pin;

  for (i = 0; i < 3; i++) {
    u[i] = 0.0;
    v[i] = 0.0;
  }
  for (i = 1; i <= n; i++) {
    j = i + 1;
    if (j > n)
      j = j - n;
    c = cos(i * pin);
    d = cos(j * pin);
    for (k = 0; k < 3; k++) {
      u[k] = u[k] + c * points[(i - 1) * 3 + k];
      v[k] = v[k] + d * points[(i - 1) * 3 + k];
    }
  }
  CrossProduct2(u, v, normal);
  c = DotProduct(normal, normal);
  c = sqrt(c);
  for (k = 0; k < 3; k++) {
    normal[k] = normal[k] / c;
  }
}

/*-----------------------------------------------------------------------------
 * LeastSqureFitNormal -- Compute normal by the least square fitting
 *                        a qudratic polynomail
 *---------------------------------------------------------------------------*/
int LeastSqureFitNormal(level0, level1, numb1, level2, numb2, normal, hessaian,
                        errorbd)
float level0[], level1[], level2[]; /* given fitting point in three levels */
int numb1, numb2;                   /* the numbers of points at level1-2   */
float normal[];                     /* the result of estimated normal      */
float hessaian[];                   /* the Hessian in the order of
                                           xx,xy,xz,yy,yz,zz               */
float errorbd;                      /* fitting error bound                 */
{
  double matrix[760], Q[5776], weight, weight0, x, y, z, x0, y0, z0;
  float error, a, c;
  int i, j, k, info, row;

  if (1 + numb1 + numb2 < 10) {
    printf("The number of fitting points < 10, -- %d\n", 1 + numb1 + numb2);
    return (0);
  }

  x0 = level0[0];
  y0 = level0[1];
  z0 = level0[2];
  weight0 = 1.0;
  error = 1000000.0;

  while (error > errorbd) {
    /* form the coefficient matrix       */
    for (i = 0; i < 760; i++)
      matrix[i] = 0.0;

    /* level 0 data                      */
    matrix[0] = weight0;

    /* level 1 data                      */
    weight = 0.5;
    for (i = 0; i < numb1; i++) {

      k = i + i + i;
      x = level1[k] - x0;
      y = level1[k + 1] - y0;
      z = level1[k + 2] - z0;

      j = i + 1;
      matrix[10 * j] = weight;
      matrix[10 * j + 1] = weight * x * x;
      matrix[10 * j + 2] = weight * x * y;
      matrix[10 * j + 3] = weight * x * z;
      matrix[10 * j + 4] = weight * y * y;
      matrix[10 * j + 5] = weight * y * z;
      matrix[10 * j + 6] = weight * z * z;
      matrix[10 * j + 7] = weight * x;
      matrix[10 * j + 8] = weight * y;
      matrix[10 * j + 9] = weight * z;
    }

    /* level 2 data                      */
    weight = 0.25;
    for (i = 0; i < numb2; i++) {
      k = i + i + i;
      x = level2[k] - x0;
      y = level2[k + 1] - y0;
      z = level2[k + 2] - z0;

      j = i + 1 + numb1;
      matrix[10 * j] = weight;
      matrix[10 * j + 1] = weight * x * x;
      matrix[10 * j + 2] = weight * x * y;
      matrix[10 * j + 3] = weight * x * z;
      matrix[10 * j + 4] = weight * y * y;
      matrix[10 * j + 5] = weight * y * z;
      matrix[10 * j + 6] = weight * z * z;
      matrix[10 * j + 7] = weight * x;
      matrix[10 * j + 8] = weight * y;
      matrix[10 * j + 9] = weight * z;
    }
    /* QR decomposition                  */
    row = 1 + numb1 + numb2;
    info = QRdecomposition(matrix, row, 10, Q);
    if (info == 0) {
      printf("The fitting matrix is column linear dependent\n");
      return (0);
    }

    if (info == 1) {
      double R11[49], R12[21], R22[9];
      double U[9], V[9], S[3], E[3], work[4];
      float x1[7], x2[3], R12x2[7];
      int n, p, ldx, ldu, ldv, job, info;

      /* QR decomposition                  */
      /* Form R11, R12, R22                */
      for (i = 0; i < 7; i++) {
        for (j = 0; j < 7; j++) {
          R11[i * 7 + j] = matrix[i * 10 + j];
        }
        for (j = 0; j < 3; j++) {
          R12[i * 3 + j] = matrix[i * 10 + j + 7];
        }
      }
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {

          /* R22 are stored in clomn to use svd */
          R22[i + j * 3] = matrix[(i + 7) * 10 + j + 7];
        }
      }

      /* SVD decomposition of R22= U\sum V^T*/
      n = 3;   /* the number of rows   */
      p = 3;   /* the number of column */
      ldx = n; /* ldx >= n             */
      ldu = n; /* ldu >= n             */
      ldv = p; /* ldv >= p             */
      job = 11;
      dsvdc_(R22, &ldx, &n, &p, S, E, U, &ldu, V, &ldv, work, &job, &info);

      /* Get x2 = the third column of V     */
      x2[0] = V[6];
      x2[1] = V[7];
      x2[2] = V[8];
      /*
         printf("x2 = %f,%f,%f\n", x2[0],x2[1],x2[2]);
      */

      /* Compute x1 = - R11^{-1} R12 x2     */
      for (i = 0; i < 7; i++) {
        R12x2[i] = 0.0;
        for (j = 0; j < 3; j++) {
          R12x2[i] = R12x2[i] - R12[i * 3 + j] * x2[j];
        }
      }

      /* Solve a triangular equation
                           R11 x1 = R12x2  */
      for (i = 0; i < 7; i++) {
        j = 6 - i;
        x1[j] = 0.0;
        for (k = j + 1; k < 7; k++) {
          x1[j] = x1[j] + R11[j * 7 + k] * x1[k];
        }
        x1[j] = (R12x2[j] - x1[j]) / R11[j * 7 + j];
      }

      /* compute the fitting error          */
      a = x1[1] * x2[0] * x2[0] + x1[2] * x2[0] * x2[1] +
          x1[3] * x2[0] * x2[2] + x1[4] * x2[1] * x2[1] +
          x1[5] * x2[1] * x2[2] + x1[6] * x2[2] * x2[2];
      c = x1[0];
      if (fabs(a) < errorbd / 100.0) {
        error = -c;
      } else {
        c = 1.0 - 4.0 * a * c;
        if (c < 0.0) {
          error = 100000.0;
        } else {
          error = (-1.0 + sqrt(c)) / (a + a);
        }
      }
      error = fabs(error);
      /*
         printf("Fitting error = %f,a,c = %f,%f, weight0 = %f\n", error,a,c,
         weight0);
      */
      normal[0] = x2[0];
      normal[1] = x2[1];
      normal[2] = x2[2];
      hessaian[0] = 2 * x1[1];
      hessaian[1] = x1[2];
      hessaian[2] = x1[3];
      hessaian[3] = 2 * x1[4];
      hessaian[4] = x1[5];
      hessaian[5] = 2 * x1[6];
    }
    weight0 = weight0 + 1.0;
  } /* end while loop                     */
  return (1);
}

/*-----------------------------------------------------------------------------
 * AdjacentTriangles-- Find the adjecent triangles for every triangle
 *---------------------------------------------------------------------------*/
void AdjacentTriangles(facetetras, numbtris, valence,
                       neighbor) Face *facetetras;
int numbtris, *valence;
Neighbor *neighbor;
{
  int i, j, ii, jj, kk, iii, jjj, kkk, j0, curval;

  for (i = 0; i < numbtris; i++) {
    ii = facetetras[i]->Index[0];
    jj = facetetras[i]->Index[1];
    kk = facetetras[i]->Index[2];
    curval = valence[ii];
    for (j = 0; j < curval; j++) {
      j0 = neighbor[ii]->NB[j];
      if (i != j0) {
        iii = facetetras[j0]->Index[0];
        jjj = facetetras[j0]->Index[1];
        kkk = facetetras[j0]->Index[2];

        if (ii == iii && jj == jjj) {
          facetetras[i]->AdjTri[2] = j0;
        }
        if (ii == jjj && jj == kkk) {
          facetetras[i]->AdjTri[2] = j0;
        }
        if (ii == iii && jj == kkk) {
          facetetras[i]->AdjTri[2] = j0;
        }

        if (ii == iii && kk == jjj) {
          facetetras[i]->AdjTri[1] = j0;
        }
        if (ii == jjj && kk == kkk) {
          facetetras[i]->AdjTri[1] = j0;
        }
        if (ii == iii && kk == kkk) {
          facetetras[i]->AdjTri[1] = j0;
        }
      }
    }
    curval = valence[jj];
    for (j = 0; j < curval; j++) {
      j0 = neighbor[jj]->NB[j];
      if (i != j0) {
        iii = facetetras[j0]->Index[0];
        jjj = facetetras[j0]->Index[1];
        kkk = facetetras[j0]->Index[2];

        if (jj == iii && kk == jjj) {
          facetetras[i]->AdjTri[0] = j0;
        }
        if (jj == jjj && kk == kkk) {
          facetetras[i]->AdjTri[0] = j0;
        }
        if (jj == iii && kk == kkk) {
          facetetras[i]->AdjTri[0] = j0;
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * Valence_Neighbor -- Compute adjecent information
 *---------------------------------------------------------------------------*/
void Valence_Neighbor(facetetras, numbpts, numbtris, faceneighbor, edgeneighbor,
                      facevalence, edgevalence) Face *facetetras;
int *facevalence, *edgevalence, numbpts, numbtris;
Neighbor *faceneighbor, *edgeneighbor;

{
  int ii, jj, kk, i, j, k, l, m, cur, iface, ifequal;
  Neighbor faceneigh, edgeneigh;
  int *facenb, *edgenb;

  /* set to zero                                   */
  for (i = 0; i < numbpts; i++) {
    facevalence[i] = 0;
  }
  /* compute the facevalence                      */
  for (i = 0; i < numbtris; i++) {
    ii = facetetras[i]->Index[0];
    jj = facetetras[i]->Index[1];
    kk = facetetras[i]->Index[2];
    facevalence[ii] = facevalence[ii] + 1;
    facevalence[jj] = facevalence[jj] + 1;
    facevalence[kk] = facevalence[kk] + 1;
  }

  // for (i = 0; i < numbpts; i++) {
  //   printf("Valence: %d %d\n", i,facevalence[i]);
  // }

  /* maccloc neighbor tri index                    */
  for (i = 0; i < numbpts; i++) {
    faceneigh = (Neighbor)malloc(sizeof(Neighbor));
    facenb = (int *)malloc(facevalence[i] * sizeof(int));
    facenb[facevalence[i] - 1] = -1;
    faceneigh->NB = facenb;
    faceneighbor[i] = faceneigh;

    /* if surface has boundary, the
       facevalence[i] is incoorect                   */
    edgeneigh = (Neighbor)malloc(sizeof(Neighbor));
    edgenb = (int *)malloc((facevalence[i] + 2) * sizeof(int));
    edgeneigh->NB = edgenb;
    edgeneighbor[i] = edgeneigh;
  }

  /* compute neighbor tri index                    */
  for (i = 0; i < numbtris; i++) {
    for (j = 0; j < 3; j++) {
      k = facetetras[i]->Index[j];
      l = faceneighbor[k]->NB[facevalence[k] - 1] + 1;
      faceneighbor[k]->NB[facevalence[k] - 1] = l;
      faceneighbor[k]->NB[l] = i;
    }
  }

  /* compute neighbor vertex index                 */
  for (i = 0; i < numbpts; i++) {
    cur = 0;
    for (j = 0; j < facevalence[i]; j++) {
      iface = faceneighbor[i]->NB[j];

      for (l = 0; l < 3; l++) {
        ii = facetetras[iface]->Index[l];
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
