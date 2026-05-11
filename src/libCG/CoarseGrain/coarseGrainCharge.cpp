#include <Utility/utility.h>
#include <libCG/CoarseGrain/coarseGrain.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>

using namespace MOLECULE;

void CoarseGrain::CoarseGrainCharge(GroupOfAtoms *groupOfBeads) {
  // for residue charge constraint
  const double h[3] = {1.0, 1.0, 1.0};
  const double cutoff = 10.0;
  int min[3], max[3];
  int Mb, Nb, Nr;
  int i, j, k;
  getBoundingBox(groupOfBeads, cutoff, min, max);
  Nb = groupOfBeads->numOfAtoms;
  Nr = groupOfBeads->numOfRes;
  Mb = Nb - Nr; // number of unknowns
  printf("%d beads, %d residues, %d unknowns\n", Nb, Nr, Mb);
  double *A = new double[Mb * Mb];
  double *b = new double[Mb];
  double *Q = new double[Nr];
  IndexRange *beadRange = new IndexRange[Nb];
  IndexRange *atomRange = new IndexRange[molecule->numOfAtoms];
  getBeadRange(min, max, cutoff, groupOfBeads, beadRange);
  getBeadRange(min, max, cutoff, molecule, atomRange);
  int residueEnd[Nr];
  getResidueEnd(groupOfBeads, residueEnd);
  totalChargeContained(groupOfBeads, Q, Nr);
  setupLinearSystemAnybead1(
      A, b, groupOfBeads, beadRange, atomRange, Q,
      residueEnd); // set up linear system for any-bead model large pqr
  LinearSystemSolver(A, b, Mb);
  assignCharge(groupOfBeads, b, Q, residueEnd);
  check(groupOfBeads, min, max, beadRange, atomRange);
}

void CoarseGrain::check(GroupOfAtoms *groupOfBeads, int *min, int *max,
                        IndexRange *beadRange, IndexRange *atomRange) {
  int i, j, k, l;
  int Nb = groupOfBeads->numOfAtoms;
  int Na = molecule->numOfAtoms;
  double *phi, *phi_cg;
  double x, y, z;
  int r[3];
  int dim1 = max[0] - min[0];
  int dim2 = max[1] - min[1];
  int dim3 = max[2] - min[2];
  // printf("%d %d %d\n", dim1, dim2, dim3);
  phi = new double[dim1 * dim2 * dim3];
  phi_cg = new double[dim1 * dim2 * dim3];
  for (i = 0; i < dim1 * dim2 * dim3; i++) {
    phi[i] = phi_cg[i] = 0.0;
  }
  for (i = 0; i < dim1; i++) {
    r[0] = min[0] + i;
    for (j = 0; j < dim2; j++) {
      r[1] = min[1] + j;
      for (k = 0; k < dim3; k++) {
        r[2] = min[2] + k;
        for (l = 0; l < Nb; l++) {
          if (isInAtomRange(r[0], r[1], r[2], l, beadRange)) {
            x = groupOfBeads->m_Atoms[l]->position[0];
            y = groupOfBeads->m_Atoms[l]->position[1];
            z = groupOfBeads->m_Atoms[l]->position[2];
            phi_cg[i * dim2 * dim3 + j * dim3 + k] +=
                groupOfBeads->m_Atoms[l]->charge /
                sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1]) +
                     (z - r[2]) * (z - r[2]));
          }
          // printf("%f\n", phi_cg[i*dim2*dim3+j*dim3+k]);
        }
        for (l = 0; l < Na; l++) {
          if (isInAtomRange(r[0], r[1], r[2], l, atomRange)) {
            x = molecule->m_Atoms[l]->position[0];
            y = molecule->m_Atoms[l]->position[1];
            z = molecule->m_Atoms[l]->position[2];
            phi[i * dim2 * dim3 + j * dim3 + k] +=
                molecule->m_Atoms[l]->charge /
                sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1]) +
                     (z - r[2]) * (z - r[2]));
          }
        }
      }
    }
  }
  double e = 0.0;
  for (i = 0; i < dim1 * dim2 * dim3; i++) {
    e += (phi[i] - phi_cg[i]) * (phi[i] - phi_cg[i]);
  }
  printf("%f\n", sqrt(e));
}

bool CoarseGrain::isInAtomRange(int x, int y, int z, int l,
                                IndexRange *atomRange) {
  int min[3], max[3];
  for (int i = 0; i < 3; i++) {
    min[i] = atomRange[l].min[i];
    max[i] = atomRange[l].max[i];
  }
  if (x > min[0] && x < max[0] && y > min[1] && y < max[1] && z > min[2] &&
      z < max[2]) {
    return true;
  }
  return false;
}

void CoarseGrain::getBoundingBox(GroupOfAtoms *groupOfBeads, double cutoff,
                                 int *Min, int *Max) {
  double min[3], max[3];
  int i, j;
  for (i = 0; i < 3; i++) {
    min[i] = 9999.0;
    max[i] = -9999.0;
  }
  for (i = 0; i < groupOfBeads->numOfAtoms; i++) {
    double x;
    double p[3];
    double r = groupOfBeads->m_Atoms[i]->radius;
    for (j = 0; j < 3; j++) {
      p[j] = groupOfBeads->m_Atoms[i]->position[j];
      if ((x = p[j] - cutoff) < min[j]) {
        min[j] = x;
      }
      if ((x = p[j] + cutoff) > max[j]) {
        max[j] = x;
      }
    }
  }
  for (i = 0; i < 3; i++) {
    Min[i] = (int)(min[i] - 1.0);
    Max[i] = (int)(max[i] + 1.0);
  }
  printf("BoundingBox: %d %d %d %d %d %d\n", Min[0], Max[0], Min[1], Max[1],
         Min[2], Max[2]);
}

void CoarseGrain::getBeadRange(int *min, int *max, double cutoff,
                               GroupOfAtoms *groupOfBeads,
                               IndexRange *beadRange) {
  int i, j;
  for (i = 0; i < groupOfBeads->numOfAtoms; i++) {
    double p[3];
    double r;
    for (j = 0; j < 3; j++) {
      p[j] = groupOfBeads->m_Atoms[i]->position[j];
    }
    r = groupOfBeads->m_Atoms[i]->radius;
    for (j = 0; j < 3; j++) {
      beadRange[i].min[j] =
          (int)(p[j] - cutoff - 1.0); // - min[j]; for relative 10/20/2008
      if (beadRange[i].min[j] < min[j]) {
        beadRange[i].min[j] = min[j];
      }
      beadRange[i].max[j] =
          (int)(p[j] + cutoff + 1.0); // - min[j]; for relative 10/20/2008
      if (beadRange[i].max[j] > max[j]) {
        beadRange[i].max[j] = max[j];
      }
    }
  }
}

void CoarseGrain::setupLinearSystem(double *A, double *b,
                                    GroupOfAtoms *groupOfBeads,
                                    IndexRange *beadRange,
                                    IndexRange *atomRange) {
  int i, j, k;
  int Mb;
  int Mb_min[3], Mb_max[3];
  double r_Mb, r_iMb, r_ij, r_jMb, aij;
  double y_Mb, y_i, bi;
  int N;
  Mb = groupOfBeads->numOfAtoms;
  N = Mb - 1; // NxN is the size of matrix
  for (k = 0; k < 3; k++) {
    Mb_min[k] = beadRange[Mb - 1].min[k];
    Mb_max[k] = beadRange[Mb - 1].max[k];
  }
  r_Mb = getDistanceSum(Mb - 1, Mb_min, Mb_max, groupOfBeads); // x_k(r_Mb)^2
  y_Mb = getyDistanceSum(Mb - 1, Mb_min, Mb_max, groupOfBeads,
                         atomRange); // y_k * x_k(r_Mb)
  for (i = 0; i < N; i++) {
    int i_min[3], i_max[3];
    for (k = 0; k < 3; k++) {
      i_min[k] = beadRange[i].min[k];
      i_max[k] = beadRange[i].max[k];
    }
    r_iMb =
        getDistanceSum(i, Mb - 1, i_min, i_max, Mb_min, Mb_max, groupOfBeads);
    for (j = i; j < N; j++) {
      int j_min[3], j_max[3];
      for (k = 0; k < 3; k++) {
        j_min[k] = beadRange[j].min[k];
        j_max[k] = beadRange[j].max[k];
      }
      r_ij = getDistanceSum(j, i, j_min, j_max, i_min, i_max, groupOfBeads);
      r_jMb =
          getDistanceSum(j, Mb - 1, j_min, j_max, Mb_min, Mb_max, groupOfBeads);
      aij = r_ij - r_jMb - r_iMb + r_Mb;
      A[i * N + j] = aij;
      /*
                              if (i == j)
                              {       printf("%f %f %f %f\n", r_ij, r_jMb,
         r_iMb, r_Mb);
                              }
      */
    }
    y_i = getyDistanceSum(i, i_min, i_max, groupOfBeads,
                          atomRange); // y_k * x_k(r_i)
    bi = y_i - y_Mb - molecule->total_charge * (r_iMb - r_Mb);
    b[i] = bi;
  }
  /*
          printf("diagonal :");
          for (i = 0; i < Mb-1; i++)
                  printf("%f ", A[i*(Mb-1)+i]);
          printf("\n");
  */
}

void CoarseGrain::setupLinearSystem2bead(double *A, double *b,
                                         GroupOfAtoms *groupOfBeads,
                                         IndexRange *beadRange,
                                         IndexRange *atomRange, double *Q) {
  int i, j, k, l;
  int Mb;
  double r_ij, r_iJ, r_Ij, r_IJ, aij, r_il, r_Il;
  double y_i, y_I, bi;
  int N;
  Mb = groupOfBeads->numOfAtoms;
  if (Mb % 2 != 0) {
    printf("odd number of beads\n");
  }
  N = Mb / 2; // NxN is the size of matrix
  for (i = 0; i < N; i++) {
    int i_min[3], i_max[3], I_min[3], I_max[3];
    int ii = 2 * i;
    int II = ii + 1;
    for (k = 0; k < 3; k++) {
      i_min[k] = beadRange[ii].min[k];
      i_max[k] = beadRange[ii].max[k];
      I_min[k] = beadRange[II].min[k];
      I_max[k] = beadRange[II].max[k];
    }
    for (j = i; j < N; j++) {
      int j_min[3], j_max[3], J_min[3], J_max[3];
      ;
      int jj = 2 * j;
      int JJ = jj + 1;
      for (k = 0; k < 3; k++) {
        j_min[k] = beadRange[jj].min[k];
        j_max[k] = beadRange[jj].max[k];
        J_min[k] = beadRange[JJ].min[k];
        J_max[k] = beadRange[JJ].max[k];
      }
      r_ij = getDistanceSum(jj, ii, j_min, j_max, i_min, i_max, groupOfBeads);
      r_iJ = getDistanceSum(JJ, ii, J_min, J_max, i_min, i_max, groupOfBeads);
      r_Ij = getDistanceSum(jj, II, j_min, j_max, I_min, I_max, groupOfBeads);
      r_IJ = getDistanceSum(JJ, II, J_min, J_max, I_min, I_max, groupOfBeads);
      aij = r_ij - r_iJ - r_Ij + r_IJ;
      A[i * N + j] = aij;
    }
    bi = 0;
    y_i = getyDistanceSum(ii, i_min, i_max, groupOfBeads,
                          atomRange); // y_k * x_k(r_i)
    y_I = getyDistanceSum(II, I_min, I_max, groupOfBeads,
                          atomRange); // y_k * x_k(r_I)
    for (l = 0; l < N; l++) {
      int l_min[3], l_max[3];
      int ll = 2 * l + 1;
      for (k = 0; k < 3; k++) {
        l_min[k] = beadRange[ll].min[k];
        l_max[k] = beadRange[ll].max[k];
      }
      r_il = getDistanceSum(ll, ii, l_min, l_max, i_min, i_max, groupOfBeads);
      r_Il = getDistanceSum(ll, II, l_min, l_max, I_min, I_max, groupOfBeads);
      bi += Q[l] * (r_il - r_Il);
    }
    bi = y_i - y_I - bi;
    b[i] = bi;
  }
}

void CoarseGrain::setupLinearSystemAnybead1(double *A, double *b,
                                            GroupOfAtoms *groupOfBeads,
                                            IndexRange *beadRange,
                                            IndexRange *atomRange, double *Q,
                                            int *residueEnd) {
  int i, j, k, l;
  double r_ij, r_iJ, r_Ij, r_IJ, aij, r_il, r_Il;
  double y_i, y_I, bi;
  int Nb = groupOfBeads->numOfAtoms;
  int Nr = groupOfBeads->numOfRes;
  int N = Nb - Nr; // NxN is the size of matrix
  //     compute the unknown list
  int unknownIndex[N];
  j = 0;
  k = 0;
  for (i = 0; i < Nb; i++) {
    if (i != residueEnd[k]) {
      unknownIndex[j++] = i;
    } else {
      k++;
    }
  }
  // construct the matrix
  int iend = 0;
  for (i = 0; i < N; i++) {
    int i_min[3], i_max[3], I_min[3], I_max[3];
    int ii = unknownIndex[i];
    if (groupOfBeads->m_Atoms[ii]->resSeq !=
        groupOfBeads->m_Atoms[residueEnd[iend]]->resSeq) {
      iend++;
    }
    int II = residueEnd[iend];
    for (k = 0; k < 3; k++) {
      i_min[k] = beadRange[ii].min[k];
      i_max[k] = beadRange[ii].max[k];
      I_min[k] = beadRange[II].min[k];
      I_max[k] = beadRange[II].max[k];
    }
    int jend = 0;
    for (j = i; j < N; j++) {
      int j_min[3], j_max[3], J_min[3], J_max[3];
      ;
      int jj = unknownIndex[j];
      while (groupOfBeads->m_Atoms[jj]->resSeq !=
             groupOfBeads->m_Atoms[residueEnd[jend]]->resSeq) {
        jend++;
      }
      int JJ = residueEnd[jend];
      for (k = 0; k < 3; k++) {
        j_min[k] = beadRange[jj].min[k];
        j_max[k] = beadRange[jj].max[k];
        J_min[k] = beadRange[JJ].min[k];
        J_max[k] = beadRange[JJ].max[k];
      }
      r_ij = getDistanceSum(jj, ii, j_min, j_max, i_min, i_max, groupOfBeads);
      r_iJ = getDistanceSum(JJ, ii, J_min, J_max, i_min, i_max, groupOfBeads);
      r_Ij = getDistanceSum(jj, II, j_min, j_max, I_min, I_max, groupOfBeads);
      r_IJ = getDistanceSum(JJ, II, J_min, J_max, I_min, I_max, groupOfBeads);
      aij = r_ij - r_iJ - r_Ij + r_IJ;
      A[i * N + j] = aij;
    }
    bi = 0;
    y_i = getyDistanceSum(ii, i_min, i_max, groupOfBeads,
                          atomRange); // y_k * x_k(r_i)
    y_I = getyDistanceSum(II, I_min, I_max, groupOfBeads,
                          atomRange); // y_k * x_k(r_I)
    for (l = 0; l < Nr; l++) {
      int l_min[3], l_max[3];
      int ll = residueEnd[l];
      for (k = 0; k < 3; k++) {
        l_min[k] = beadRange[ll].min[k];
        l_max[k] = beadRange[ll].max[k];
      }
      r_il = getDistanceSum(ll, ii, l_min, l_max, i_min, i_max, groupOfBeads);
      r_Il = getDistanceSum(ll, II, l_min, l_max, I_min, I_max, groupOfBeads);
      bi += Q[l] * (r_il - r_Il);
    }
    bi = y_i - y_I - bi;
    b[i] = bi;
  }
}

void CoarseGrain::setupLinearSystemAnybead2(
    double *AtA, double *Atb, GroupOfAtoms *groupOfBeads, IndexRange *beadRange,
    IndexRange *atomRange, double *Q, int *residueEnd, int *min, int *max) {
  int i, j, k, l, t;
  double r[3];
  int beadIndex;
  int Nb = groupOfBeads->numOfAtoms;
  int Nr = groupOfBeads->numOfRes;
  int N = Nb - Nr; // NxN is the size of matrix
  int dim[3];
  dim[0] = max[0] - min[0];
  dim[1] = max[1] - min[1];
  dim[2] = max[2] - min[2];
  //     compute the unknown list
  int unknownIndex[N];
  j = 0;
  k = 0;
  for (i = 0; i < Nb; i++) {
    if (i != residueEnd[k]) {
      unknownIndex[j++] = i;
    } else {
      k++;
    }
  }
  // construct the matrix
  printf("%d\n", dim[0] * dim[1] * dim[2]);
  double A[N][dim[0] * dim[1] * dim[2]];
  double b[dim[0] * dim[1] * dim[2]];
  double temp[Nr];
  for (i = 0; i < N * N; i++) {
    AtA[i] = 0.0;
  }
  for (i = 0; i < N; i++) {
    Atb[i] = 0.0;
  }
  for (i = 0; i < dim[0] * dim[1] * dim[2]; i++) {
    for (j = 0; j < N; j++) {
      A[j][i] = 0.0;
    }
    b[i] = 0.0;
  }
  for (i = 0; i < dim[0]; i++) {
    r[0] = (double)(i + min[0]);
    for (j = 0; j < dim[1]; j++) {
      r[1] = (double)(j + min[1]);
      for (k = 0; k < dim[2]; k++) {
        r[2] = (double)(k + min[2]);
        for (l = 0; l < Nr; l++) {
          beadIndex = residueEnd[l];
          if (isInAtomRange(r[0], r[1], r[2], beadIndex, beadRange)) {
            double x = groupOfBeads->m_Atoms[beadIndex]->position[0];
            double y = groupOfBeads->m_Atoms[beadIndex]->position[1];
            double z = groupOfBeads->m_Atoms[beadIndex]->position[2];
            temp[l] =
                1.0 / sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1]) +
                           (z - r[2]) * (z - r[2]));
          } else {
            temp[l] = 0.0;
          }
        }
        for (t = 0; t < N; t++) {
          beadIndex = unknownIndex[t];
          if (isInAtomRange(r[0], r[1], r[2], beadIndex, beadRange)) {
            double x = groupOfBeads->m_Atoms[beadIndex]->position[0];
            double y = groupOfBeads->m_Atoms[beadIndex]->position[1];
            double z = groupOfBeads->m_Atoms[beadIndex]->position[2];
            A[t][i * dim[1] * dim[2] + j * dim[2] + k] =
                1.0 / sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1]) +
                           (z - r[2]) * (z - r[2]));
          }
          l = getResidueIndex(beadIndex, residueEnd, Nr);
          A[t][i * dim[1] * dim[2] + j * dim[2] + k] -= temp[l];
        }
        for (t = 0; t < molecule->numOfAtoms; t++) {
          if (isInAtomRange(r[0], r[1], r[2], t, atomRange)) {
            double x = molecule->m_Atoms[t]->position[0];
            double y = molecule->m_Atoms[t]->position[1];
            double z = molecule->m_Atoms[t]->position[2];
            b[i * dim[1] * dim[2] + j * dim[2] + k] +=
                molecule->m_Atoms[t]->charge /
                sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1]) +
                     (z - r[2]) * (z - r[2]));
          }
        }
        for (l = 0; l < Nr; l++) {
          b[i * dim[1] * dim[2] + j * dim[2] + k] -= Q[l] * temp[l];
        }
      }
    }
  }
  for (i = 0; i < N; i++) {
    for (j = i; j < N; j++) {
      for (k = 0; k < dim[0] * dim[1] * dim[2]; k++) {
        AtA[i * N + j] += A[i][k] * A[j][k];
      }
    }
    for (k = 0; k < dim[0] * dim[1] * dim[2]; k++) {
      Atb[i] += A[i][k] * b[k];
    }
  }
}

int CoarseGrain::getResidueIndex(int beadIndex, int *residueEnd, int Nr) {
  int i = 0;
  while (residueEnd[i] < beadIndex) {
    i++;
  }
  return i;
}

// r_Mb = getDistanceSum(Mb-1, Mb_min, Mb_max, groupOfBeads);
double CoarseGrain::getDistanceSum(int M, int *Mmin, int *Mmax,
                                   GroupOfAtoms *groupOfBeads) {
  int kx, ky, kz;
  double rx, ry, rz;
  double sum = 0.0;
  for (kz = Mmin[2]; kz <= Mmax[2]; kz++) {
    rz = (groupOfBeads->m_Atoms[M]->position[2] - kz) *
         (groupOfBeads->m_Atoms[M]->position[2] - kz);
    for (ky = Mmin[1]; ky <= Mmax[1]; ky++) {
      ry = (groupOfBeads->m_Atoms[M]->position[1] - ky) *
           (groupOfBeads->m_Atoms[M]->position[1] - ky);
      for (kx = Mmin[0]; kx <= Mmax[0]; kx++) {
        rx = (groupOfBeads->m_Atoms[M]->position[0] - kx) *
             (groupOfBeads->m_Atoms[M]->position[0] - kx);
        double rm = rx + ry + rz;
        if (rm < 1.0e-5) {
          sum += 1.0 / (1.0e-5);
        } else {
          sum += 1.0 / rm;
        }
      }
    }
  }
  return sum;
}

double CoarseGrain::getDistanceSum(int i, int j, int *imin, int *imax,
                                   int *jmin, int *jmax,
                                   GroupOfAtoms *groupOfBeads) {
  int kx, ky, kz;
  double rik[3], rjk[3];
  double sum = 0.0;
  for (kz = imin[2]; kz <= imax[2]; kz++) {
    if (kz < jmin[2] || kz > jmax[2]) {
      continue;
    }
    rik[2] = (groupOfBeads->m_Atoms[i]->position[2] - kz) *
             (groupOfBeads->m_Atoms[i]->position[2] - kz);
    rjk[2] = (groupOfBeads->m_Atoms[j]->position[2] - kz) *
             (groupOfBeads->m_Atoms[j]->position[2] - kz);
    for (ky = imin[1]; ky <= imax[1]; ky++) {
      if (ky < jmin[1] || ky > jmax[1]) {
        continue;
      }
      rik[1] = (groupOfBeads->m_Atoms[i]->position[1] - ky) *
               (groupOfBeads->m_Atoms[i]->position[1] - ky);
      rjk[1] = (groupOfBeads->m_Atoms[j]->position[1] - ky) *
               (groupOfBeads->m_Atoms[j]->position[1] - ky);
      for (kx = imin[0]; kx <= imax[0]; kx++) {
        if (kx < jmin[0] || kx > jmax[0]) {
          continue;
        }
        rik[0] = (groupOfBeads->m_Atoms[i]->position[0] - kx) *
                 (groupOfBeads->m_Atoms[i]->position[0] - kx);
        rjk[0] = (groupOfBeads->m_Atoms[j]->position[0] - kx) *
                 (groupOfBeads->m_Atoms[j]->position[0] - kx);
        double ri =
            sqrt((rik[0] + rik[1] + rik[2]) * (rjk[0] + rjk[1] + rjk[2]));
        if (ri < 1.0e-5) {
          sum += 1.0 / (1.0e-5);
        } else {
          sum += 1.0 / ri;
        }
      }
    }
  }
  return sum;
}

// y_Mb = getyDistanceSum(Mb-1, Mb_min, Mb_max, groupOfBeads, atomRange);
double CoarseGrain::getyDistanceSum(int M, int *Mmin, int *Mmax,
                                    GroupOfAtoms *groupOfBeads,
                                    IndexRange *atomRange) {
  int kx, ky, kz;
  double rx, ry, rz;
  double sum = 0.0;
  for (kz = Mmin[2]; kz <= Mmax[2]; kz++) {
    rz = (groupOfBeads->m_Atoms[M]->position[2] - kz) *
         (groupOfBeads->m_Atoms[M]->position[2] - kz);
    for (ky = Mmin[1]; ky <= Mmax[1]; ky++) {
      ry = (groupOfBeads->m_Atoms[M]->position[1] - ky) *
           (groupOfBeads->m_Atoms[M]->position[1] - ky);
      for (kx = Mmin[0]; kx <= Mmax[0]; kx++) {
        rx = (groupOfBeads->m_Atoms[M]->position[0] - kx) *
             (groupOfBeads->m_Atoms[M]->position[0] - kx);
        double rm = sqrt(rx + ry + rz);
        double yk = getYk(kx, ky, kz, molecule, atomRange);
        if (rm < 1.0e-5) {
          sum += yk / (1.0e-5);
        } else {
          sum += yk / rm;
        }
      }
    }
  }
  return sum;
}

double CoarseGrain::getYk(int kx, int ky, int kz, GroupOfAtoms *goa,
                          IndexRange *atomRange) {
  int i, j;
  int atom_min[3], atom_max[3];
  double p[3], q;
  double sum = 0.0;
  for (i = 0; i < goa->numOfAtoms; i++) {
    for (j = 0; j < 3; j++) {
      atom_min[j] = atomRange[i].min[j];
      atom_max[j] = atomRange[i].max[j];
      p[j] = goa->m_Atoms[i]->position[j];
    }
    q = goa->m_Atoms[i]->charge;
    if (kx < atom_min[0] || kx > atom_max[0] || ky < atom_min[1] ||
        ky > atom_max[1] || kz < atom_min[2] || kz > atom_max[2]) {
      continue;
    }
    double rik = sqrt((p[0] - kx) * (p[0] - kx) + (p[1] - ky) * (p[1] - ky) +
                      (p[2] - kz) * (p[2] - kz));
    if (rik < 1.0e-5) {
      sum += q / (1.0e-5);
    } else {
      sum += q / rik;
    }
  }
  return sum;
}

void CoarseGrain::getResidueEnd(GroupOfAtoms *groupOfBeads, int *residueEnd) {
  int Nb = groupOfBeads->numOfAtoms;
  int Nr = groupOfBeads->numOfRes;
  int res_seq = groupOfBeads->m_Atoms[0]->resSeq;
  int i, j;
  j = 0;
  for (i = 0; i < Nb; i++) {
    if (res_seq != groupOfBeads->m_Atoms[i]->resSeq) {
      residueEnd[j++] = i - 1;
      res_seq = groupOfBeads->m_Atoms[i]->resSeq;
    }
  }
  residueEnd[Nr - 1] = Nb - 1;
  /*
          for ( i = 0; i < Nr; i++)
                  printf("%d\n", residueEnd[i]);
  */
}

void CoarseGrain::totalChargeContained(GroupOfAtoms *groupOfBeads, double *Q,
                                       int N) {
  int i, j, k;
  int resIndex = 0;
  int resSeq;
  for (i = 0; i < N; i++) {
    Q[i] = 0.0;
  }
  resSeq = groupOfBeads->m_Atoms[0]->resSeq;
  for (i = 0; i < groupOfBeads->numOfAtoms; i++) {
    //            printf("%d\n", groupOfBeads->m_Atoms[i]->resSeq);
    if (resSeq == groupOfBeads->m_Atoms[i]->resSeq) {
      for (j = 0; j < groupOfBeads->m_Atoms[i]->numOfChildren; j++) {
        k = groupOfBeads->m_Atoms[i]->child[j];
        Q[resIndex] += molecule->m_Atoms[k]->charge;
      }
    } else {
      resIndex++;
      for (j = 0; j < groupOfBeads->m_Atoms[i]->numOfChildren; j++) {
        k = groupOfBeads->m_Atoms[i]->child[j];
        Q[resIndex] += molecule->m_Atoms[k]->charge;
      }
      resSeq = groupOfBeads->m_Atoms[i]->resSeq;
    }
  }
  if (resIndex != N - 1) {
    printf("N=%d, resIndex=%d\n", N, resIndex);
    exit(0);
  }
}

void CoarseGrain::assignCharge(GroupOfAtoms *groupOfBeads, double *b, double *Q,
                               int *residueEnd) {
  int Nb = groupOfBeads->numOfAtoms;
  int endIndex = 0;
  int ii = 0;
  double qq = 0.0;
  for (int i = 0; i < Nb; i++) {
    if (i != residueEnd[endIndex]) {
      groupOfBeads->m_Atoms[i]->setCharge(b[ii]);
      qq += b[ii];
      ii++;
    } else {
      double charge = Q[endIndex] - qq;
      groupOfBeads->m_Atoms[i]->setCharge(charge);
      endIndex++;
      qq = 0.0;
    }
  }
}
