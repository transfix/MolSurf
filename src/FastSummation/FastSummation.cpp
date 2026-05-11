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
#include <cstdlib>
#include <cstdio>
#include <cstring>
// #include <complex.h>
#include <cmath>
#include <pthread.h>

#include <iostream>
#include "FastSummation/FastSummation.h"

extern "C" {
#include <complex.h>

#include <FastSummation/kernels.h>
#include <FastSummation/fastsum.h>
/* nfft3util.h was removed in NFFT >= 3.3; nfft3.h declares all utility
 * functions used here. */
#include <nfft3.h>
}

using namespace std;

FastSummation::FastSummation() {}

FastSummation::FastSummation(int dim, int N_source, int M_target, int expansion,
                             int cut_off, int p_degree, char *kernel_func,
                             double *kernel_param, double e_I, double e_B) {
  d = dim;
  N = N_source;
  M = M_target;
  n = expansion;
  m = cut_off;
  p = p_degree;

  // printf("%s\n", kernel_func);
  if (strcmp(kernel_func, "gaussian") == 0)
    kernel = gaussian;
  else if (strcmp(kernel_func, "multiquadric") == 0)
    kernel = multiquadric;
  else if (strcmp(kernel_func, "inverse_multiquadric") == 0)
    kernel = inverse_multiquadric;
  else if (strcmp(kernel_func, "logarithm") == 0)
    kernel = logarithm;
  else if (strcmp(kernel_func, "thinplate_spline") == 0)
    kernel = thinplate_spline;
  else if (strcmp(kernel_func, "one_over_square") == 0)
    kernel = one_over_square;
  else if (strcmp(kernel_func, "one_over_modulus") == 0)
    kernel = one_over_modulus;
  else if (strcmp(kernel_func, "one_over_x") == 0)
    kernel = one_over_x;
  else if (strcmp(kernel_func, "inverse_multiquadric3") == 0)
    kernel = inverse_multiquadric3;
  else if (strcmp(kernel_func, "sinc_kernel") == 0)
    kernel = sinc_kernel;
  else if (strcmp(kernel_func, "cosc") == 0)
    kernel = cosc;
  else if (strcmp(kernel_func, "cot") == 0)
    kernel = cot;
  else if (strcmp(kernel_func, "one_over_quartic") == 0)
    kernel = one_over_quartic;
  else if (strcmp(kernel_func, "one_over_seventh") == 0)
    kernel = one_over_seventh;
  else if (strcmp(kernel_func, "one_over_sixth") == 0)
    kernel = one_over_sixth;
  else if (strcmp(kernel_func, "Still_GB") == 0)
    kernel = still_gb;
  else {
    printf("No Kernel Match!\n");
    kernel_func = (char *)"multiquadric";
    kernel = multiquadric;
  }

  c = kernel_param;
  eps_I = e_I;
  eps_B = e_B;
  //      my_fastsum_plan = (fastsum_plan *)malloc(sizeof(fastsum_plan));
}

FastSummation::FastSummation(int dim, int N_source, int M_target, int expansion,
                             int cut_off, int p_degree, char *kernel_func1,
                             char *kernel_func2, double *kernel_param,
                             double e_I, double e_B) {
  d = dim;
  N = N_source;
  M = M_target;
  n = expansion;
  m = cut_off;
  p = p_degree;

  printf("%s and %s\n", kernel_func1, kernel_func2);
  if (strcmp(kernel_func1, "gaussian") == 0)
    kernel = gaussian;
  else if (strcmp(kernel_func1, "multiquadric") == 0)
    kernel = multiquadric;
  else if (strcmp(kernel_func1, "inverse_multiquadric") == 0)
    kernel = inverse_multiquadric;
  else if (strcmp(kernel_func1, "logarithm") == 0)
    kernel = logarithm;
  else if (strcmp(kernel_func1, "thinplate_spline") == 0)
    kernel = thinplate_spline;
  else if (strcmp(kernel_func1, "one_over_square") == 0)
    kernel = one_over_square;
  else if (strcmp(kernel_func1, "one_over_modulus") == 0)
    kernel = one_over_modulus;
  else if (strcmp(kernel_func1, "one_over_x") == 0)
    kernel = one_over_x;
  else if (strcmp(kernel_func1, "inverse_multiquadric3") == 0)
    kernel = inverse_multiquadric3;
  else if (strcmp(kernel_func1, "sinc_kernel") == 0)
    kernel = sinc_kernel;
  else if (strcmp(kernel_func1, "cosc") == 0)
    kernel = cosc;
  else if (strcmp(kernel_func1, "cot") == 0)
    kernel = cot;
  else if (strcmp(kernel_func1, "one_over_quartic") == 0)
    kernel = one_over_quartic;
  else if (strcmp(kernel_func1, "one_over_seventh") == 0)
    kernel = one_over_seventh;
  else if (strcmp(kernel_func1, "one_over_sixth") == 0)
    kernel = one_over_sixth;
  else if (strcmp(kernel_func1, "Still_GB") == 0)
    kernel = still_gb;
  else {
    printf("No Kernel Match!\n");
    kernel_func1 = (char *)"multiquadric";
    kernel = multiquadric;
  }

  if (strcmp(kernel_func2, "gaussian") == 0)
    kernel2 = gaussian;
  else if (strcmp(kernel_func2, "multiquadric") == 0)
    kernel2 = multiquadric;
  else if (strcmp(kernel_func2, "inverse_multiquadric") == 0)
    kernel2 = inverse_multiquadric;
  else if (strcmp(kernel_func2, "logarithm") == 0)
    kernel2 = logarithm;
  else if (strcmp(kernel_func2, "thinplate_spline") == 0)
    kernel2 = thinplate_spline;
  else if (strcmp(kernel_func2, "one_over_square") == 0)
    kernel2 = one_over_square;
  else if (strcmp(kernel_func2, "one_over_modulus") == 0)
    kernel2 = one_over_modulus;
  else if (strcmp(kernel_func2, "one_over_x") == 0)
    kernel2 = one_over_x;
  else if (strcmp(kernel_func2, "inverse_multiquadric3") == 0)
    kernel2 = inverse_multiquadric3;
  else if (strcmp(kernel_func2, "sinc_kernel") == 0)
    kernel2 = sinc_kernel;
  else if (strcmp(kernel_func2, "cosc") == 0)
    kernel2 = cosc;
  else if (strcmp(kernel_func2, "cot") == 0)
    kernel2 = cot;
  else if (strcmp(kernel_func2, "one_over_quartic") == 0)
    kernel2 = one_over_quartic;
  else if (strcmp(kernel_func2, "one_over_seventh") == 0)
    kernel2 = one_over_seventh;
  else if (strcmp(kernel_func2, "one_over_sixth") == 0)
    kernel2 = one_over_sixth;
  else if (strcmp(kernel_func2, "Still_GB") == 0)
    kernel2 = still_gb;
  else {
    printf("No Kernel Match!\n");
    kernel_func2 = (char *)"multiquadric";
    kernel2 = multiquadric;
  }

  c = kernel_param;
  eps_I = e_I;
  eps_B = e_B;
  //      my_fastsum_plan = (fastsum_plan *)malloc(sizeof(fastsum_plan));
}

void FastSummation::fastSumThread(int p, double *alpha, double *sum) {
  for (int k = 0; k < N; k++)
    my_par_fastsum_plan[p].alpha[k] = (double _Complex)(alpha[k] + I * 0.0);

  for (int k = 0; k < N; k++) {

    my_par_fastsum_plan[p].mv1.f[k] = (double _Complex)(alpha[k] + I * 0.0);
    // arand hack: this was wrong elsewhere...
    // my_par_fastsum_plan[ p ].mv1.f[ k ][0] = alpha[ k ];
    // my_par_fastsum_plan[ p ].mv1.f[ k ][1] = 0.0;
  }

  //   double t = getTime( ); //nfft_second( );

  fastsum_trafo(&my_par_fastsum_plan[p]);
  //   t = nfft_second( ) - t;
  //   t = getTime( ) - t;

  //   printf( "thread %d: %f sec\n", p + 1, t );

  for (int k = 0; k < M; k++)
    sum[k] = creal(my_par_fastsum_plan[p].f[k]);
}

void FastSummation::fastSumCorrectionThread(int p, double *alpha, double *sum) {
  for (int k = 0; k < N; k++)
    my_par_corr_fastsum_plan[p].alpha[k] =
        (double _Complex)(alpha[k] + I * 0.0);

  for (int k = 0; k < N; k++) {

    my_par_corr_fastsum_plan[p].mv1.f[k] =
        (double _Complex)(alpha[k] + I * 0.0);
    // arand hack, below is wrong
    // my_par_corr_fastsum_plan[ p ].mv1.f[ k ][0] = alpha[ k ];
    // my_par_corr_fastsum_plan[ p ].mv1.f[ k ][1] = 0.0;
  }

  fastsum_trafo(&my_par_corr_fastsum_plan[p]);
  for (int k = 0; k < M; k++)
    sum[k] = creal(my_par_corr_fastsum_plan[p].f[k]);
}

void FastSummation::threadedFastSumCorrection(
    double *x, double *y, double *alpha, double *alphax, double *alphay,
    double *alphaz, double *sum, double *sumx, double *sumy, double *sumz,
    double *c_sum, double *c_sumx, double *c_sumy, double *c_sumz) {

  double tm;
  fastsum_plan *ths[8];

  for (int q = 0; q < 4; q++) {
    fastsum_init_guru(&my_par_corr_fastsum_plan[q], d, N, M, kernel, c, 0, n, m,
                      p, eps_I, eps_B);
    fastsum_init_guru(&my_par_corr_fastsum_plan[q + 4], d, N, M, kernel2, c, 0,
                      n, m, p, eps_I, eps_B);
  }

  //   tm = nfft_second( );

  BuildTree(d, 0, x, alpha, alphax, alphay, alphaz, N);

  for (int i = 0; i < d * N; i++) {
    my_par_corr_fastsum_plan[0].x[i] = my_par_corr_fastsum_plan[1].x[i] =
        my_par_corr_fastsum_plan[2].x[i] = my_par_corr_fastsum_plan[3].x[i] =
            x[i];
    my_par_corr_fastsum_plan[4].x[i] = my_par_corr_fastsum_plan[5].x[i] =
        my_par_corr_fastsum_plan[6].x[i] = my_par_corr_fastsum_plan[7].x[i] =
            x[i];
  }

  for (int i = 0; i < M * d; i++) {
    my_par_corr_fastsum_plan[0].y[i] = my_par_corr_fastsum_plan[1].y[i] =
        my_par_corr_fastsum_plan[2].y[i] = my_par_corr_fastsum_plan[3].y[i] =
            y[i];
    my_par_corr_fastsum_plan[4].y[i] = my_par_corr_fastsum_plan[5].y[i] =
        my_par_corr_fastsum_plan[6].y[i] = my_par_corr_fastsum_plan[7].y[i] =
            y[i];
  }

  for (int p = 0; p < 8; p++)
    ths[p] = &my_par_corr_fastsum_plan[p];

  if (!(ths[0]->flags & EXACT_NEARFIELD)) {
    if (ths[0]->d == 1) {
      for (int i = -ths[0]->Ad / 2 - 2; i <= ths[0]->Ad / 2 + 2; i++) {
        ths[0]->Add[i + ths[0]->Ad / 2 + 2] =
            ths[1]->Add[i + ths[1]->Ad / 2 + 2] =
                ths[2]->Add[i + ths[2]->Ad / 2 + 2] =
                    ths[3]->Add[i + ths[3]->Ad / 2 + 2] = (double _Complex)(
                        regkern1(ths[0]->k,
                                 ths[0]->eps_I * (double)i / ths[0]->Ad * 2,
                                 ths[0]->p, ths[0]->kernel_param, ths[0]->eps_I,
                                 ths[0]->eps_B) +
                        I * 0.0);
        ths[4]->Add[i + ths[4]->Ad / 2 + 2] =
            ths[5]->Add[i + ths[5]->Ad / 2 + 2] =
                ths[6]->Add[i + ths[6]->Ad / 2 + 2] =
                    ths[7]->Add[i + ths[7]->Ad / 2 + 2] = (double _Complex)(
                        regkern1(ths[0]->k,
                                 ths[0]->eps_I * (double)i / ths[0]->Ad * 2,
                                 ths[0]->p, ths[0]->kernel_param, ths[0]->eps_I,
                                 ths[0]->eps_B) +
                        I * 0.0);
      }
    } else {
      for (int i = 0; i <= ths[0]->Ad + 2; i++) {
        ths[0]->Add[i] = ths[1]->Add[i] = ths[2]->Add[i] = ths[3]->Add[i] =
            (double _Complex)(regkern3(ths[0]->k,
                                       ths[0]->eps_I * (double)i / ths[0]->Ad,
                                       ths[0]->p, ths[0]->kernel_param,
                                       ths[0]->eps_I, ths[0]->eps_B) +
                              I * 0.0);
        ths[4]->Add[i] = ths[5]->Add[i] = ths[6]->Add[i] = ths[7]->Add[i] =
            (double _Complex)(regkern3(ths[0]->k,
                                       ths[0]->eps_I * (double)i / ths[0]->Ad,
                                       ths[0]->p, ths[0]->kernel_param,
                                       ths[0]->eps_I, ths[0]->eps_B) +
                              I * 0.0);
      }
    }
  }

  for (int q = 0; q < 8; q++) {

    for (int i = 0; i < ths[q]->mv1.M_total; i++)
      for (int j = 0; j < ths[q]->mv1.d; j++)
        ths[q]->mv1.x[ths[q]->mv1.d * i + j] =
            -ths[q]->x[ths[q]->mv1.d * i + j];

    if (ths[q]->mv1.flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&(ths[q]->mv1));

    if (ths[q]->mv1.flags & PRE_PSI)
      nfft_precompute_psi(&(ths[q]->mv1));

    if (ths[q]->mv1.flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&(ths[q]->mv1));

    for (int i = 0; i < ths[q]->mv2.M_total; i++)
      for (int j = 0; j < ths[q]->mv2.d; j++)
        ths[q]->mv2.x[ths[q]->mv2.d * i + j] =
            -ths[q]->y[ths[q]->mv2.d * i + j];

    if (ths[q]->mv2.flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&(ths[q]->mv2));

    if (ths[q]->mv2.flags & PRE_PSI)
      nfft_precompute_psi(&(ths[q]->mv2));

    if (ths[q]->mv2.flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&(ths[q]->mv2));
  }

  int n_total = 1;

  for (int i = 0; i < ths[0]->d; i++)
    n_total *= ths[0]->n;

  for (int j = 0; j < n_total; j++) {
    if (ths[0]->d == 1) {
      ths[0]->b[j] = ths[1]->b[j] = ths[2]->b[j] = ths[3]->b[j] =
          (double _Complex)(regkern1(ths[0]->k, (double)j / (ths[0]->n) - 0.5,
                                     ths[0]->p, ths[0]->kernel_param,
                                     ths[0]->eps_I, ths[0]->eps_B) /
                                n_total +
                            I * 0.0);

      ths[0]->b[j] = ths[1]->b[j] = ths[2]->b[j] = ths[3]->b[j] =
          (double _Complex)(regkern1(ths[0]->k, (double)j / (ths[0]->n) - 0.5,
                                     ths[0]->p, ths[0]->kernel_param,
                                     ths[0]->eps_I, ths[0]->eps_B) /
                                n_total +
                            I * 0.0);

      ths[4]->b[j] = ths[5]->b[j] = ths[6]->b[j] = ths[7]->b[j] =
          (double _Complex)(regkern1(ths[0]->k, (double)j / (ths[0]->n) - 0.5,
                                     ths[0]->p, ths[0]->kernel_param,
                                     ths[0]->eps_I, ths[0]->eps_B) /
                                n_total +
                            I * 0.0);
    } else {
      int k = j;

      ths[0]->b[j] = 0.0 + I * 0.0;

      for (int t = 0; t < ths[0]->d; t++) {
        ths[0]->b[j] += (double _Complex)(
            ((double)(k % (ths[0]->n)) / (ths[0]->n) - 0.5) *
                ((double)(k % (ths[0]->n)) / (ths[0]->n) - 0.5) +
            I * 0.0);
        // bad below
        // ths[ 0 ]->b[ j ] +=  (double _Complex ) ( ( ( ( double ) ( k % ( ths[
        // 0 ]->n ) ) / ( ths[ 0 ]->n ) - 0.5 ) * ( ( double ) ( k % ( ths[ 0
        // ]->n ) ) / ( ths[ 0 ]->n ) - 0.5 ) ) + I*0.0);
        k = k / (ths[0]->n);
      }

      ths[0]->b[j] = ths[1]->b[j] = ths[2]->b[j] = ths[3]->b[j] =
          (double _Complex)(
              regkern3(ths[0]->k, sqrt((double)creal(ths[0]->b[j])), ths[0]->p,
                       ths[0]->kernel_param, ths[0]->eps_I, ths[0]->eps_B) /
                  n_total +
              I * 0.0);
      ths[4]->b[j] = ths[5]->b[j] = ths[6]->b[j] = ths[7]->b[j] =
          (double _Complex)(
              regkern3(ths[0]->k, sqrt((double)creal(ths[0]->b[j])), ths[0]->p,
                       ths[0]->kernel_param, ths[0]->eps_I, ths[0]->eps_B) /
                  n_total +
              I * 0.0);
    }
  }

  for (int q = 0; q < 8; q++) {
    nfft_fftshift_complex(ths[q]->b, ths[q]->mv1.d, ths[q]->mv1.N);
    fftw_execute(ths[q]->fft_plan);
    nfft_fftshift_complex(ths[q]->b, ths[q]->mv1.d, ths[q]->mv1.N);
  }

  //   tm = nfft_second( ) - tm;
  //   printf("pre-computation: %f sec\n", tm );

  PARAMS pr[8];
  pthread_t pt[8];

  pr[0].p = 0;
  pr[0].alpha = alpha;
  pr[0].sum = sum;
  pr[0].FS = this;

  pr[1].p = 1;
  pr[1].alpha = alphax;
  pr[1].sum = sumx;
  pr[1].FS = this;

  pr[2].p = 2;
  pr[2].alpha = alphay;
  pr[2].sum = sumy;
  pr[2].FS = this;

  pr[3].p = 3;
  pr[3].alpha = alphaz;
  pr[3].sum = sumz;
  pr[3].FS = this;

  pr[4].p = 4;
  pr[4].alpha = alpha;
  pr[4].sum = c_sum;
  pr[4].FS = this;

  pr[5].p = 5;
  pr[5].alpha = alphax;
  pr[5].sum = c_sumx;
  pr[5].FS = this;

  pr[6].p = 6;
  pr[6].alpha = alphay;
  pr[6].sum = c_sumy;
  pr[6].FS = this;

  pr[7].p = 7;
  pr[7].alpha = alphaz;
  pr[7].sum = c_sumz;
  pr[7].FS = this;

  for (int q = 0; q < 8; q++)
    pthread_create(&pt[q], NULL, startFastSumCorrectionThread, (void *)&pr[q]);

  for (int q = 0; q < 8; q++)
    pthread_join(pt[q], NULL);

  for (int q = 0; q < 8; q++)
    fastsum_finalize(&my_par_corr_fastsum_plan[q]);
}

void FastSummation::threadedFastSum(double *x, double *y, double *alpha,
                                    double *alphax, double *alphay,
                                    double *alphaz, double *sum, double *sumx,
                                    double *sumy, double *sumz) {
  double tm;
  fastsum_plan *ths[4];

  for (int q = 0; q < 4; q++)
    fastsum_init_guru(&my_par_fastsum_plan[q], d, N, M, kernel, c, 0, n, m, p,
                      eps_I, eps_B);

  //   tm = nfft_second( );

  BuildTree(d, 0, x, alpha, alphax, alphay, alphaz, N);

  for (int i = 0; i < d * N; i++)
    my_par_fastsum_plan[0].x[i] = my_par_fastsum_plan[1].x[i] =
        my_par_fastsum_plan[2].x[i] = my_par_fastsum_plan[3].x[i] = x[i];

  for (int i = 0; i < M * d; i++)
    my_par_fastsum_plan[0].y[i] = my_par_fastsum_plan[1].y[i] =
        my_par_fastsum_plan[2].y[i] = my_par_fastsum_plan[3].y[i] = y[i];

  for (int p = 0; p < 4; p++)
    ths[p] = &my_par_fastsum_plan[p];

  if (!(ths[0]->flags & EXACT_NEARFIELD)) {
    if (ths[0]->d == 1)
      for (int i = -ths[0]->Ad / 2 - 2; i <= ths[0]->Ad / 2 + 2; i++)
        ths[0]->Add[i + ths[0]->Ad / 2 + 2] =
            ths[1]->Add[i + ths[1]->Ad / 2 + 2] =
                ths[2]->Add[i + ths[2]->Ad / 2 + 2] =
                    ths[3]->Add[i + ths[3]->Ad / 2 + 2] =
                        regkern1(ths[0]->k,
                                 ths[0]->eps_I * (double)i / ths[0]->Ad * 2,
                                 ths[0]->p, ths[0]->kernel_param, ths[0]->eps_I,
                                 ths[0]->eps_B) +
                        I * 0.0;
    else
      for (int i = 0; i <= ths[0]->Ad + 2; i++)
        ths[0]->Add[i] = ths[1]->Add[i] = ths[2]->Add[i] = ths[3]->Add[i] =
            regkern3(ths[0]->k, ths[0]->eps_I * (double)i / ths[0]->Ad,
                     ths[0]->p, ths[0]->kernel_param, ths[0]->eps_I,
                     ths[0]->eps_B) +
            I * 0.0;
  }

  for (int q = 0; q < 4; q++) {
    for (int i = 0; i < ths[q]->mv1.M_total; i++)
      for (int j = 0; j < ths[q]->mv1.d; j++)
        ths[q]->mv1.x[ths[q]->mv1.d * i + j] =
            -ths[q]->x[ths[q]->mv1.d * i + j];

    if (ths[q]->mv1.flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&(ths[q]->mv1));

    if (ths[q]->mv1.flags & PRE_PSI)
      nfft_precompute_psi(&(ths[q]->mv1));

    if (ths[q]->mv1.flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&(ths[q]->mv1));

    for (int i = 0; i < ths[q]->mv2.M_total; i++)
      for (int j = 0; j < ths[q]->mv2.d; j++)
        ths[q]->mv2.x[ths[q]->mv2.d * i + j] =
            -ths[q]->y[ths[q]->mv2.d * i + j];

    if (ths[q]->mv2.flags & PRE_LIN_PSI)
      nfft_precompute_lin_psi(&(ths[q]->mv2));

    if (ths[q]->mv2.flags & PRE_PSI)
      nfft_precompute_psi(&(ths[q]->mv2));

    if (ths[q]->mv2.flags & PRE_FULL_PSI)
      nfft_precompute_full_psi(&(ths[q]->mv2));
  }

  int n_total = 1;

  for (int i = 0; i < ths[0]->d; i++)
    n_total *= ths[0]->n;

  for (int j = 0; j < n_total; j++) {
    if (ths[0]->d == 1)
      ths[0]->b[j] = ths[1]->b[j] = ths[2]->b[j] = ths[3]->b[j] =
          (double _Complex)(regkern1(ths[0]->k, (double)j / (ths[0]->n) - 0.5,
                                     ths[0]->p, ths[0]->kernel_param,
                                     ths[0]->eps_I, ths[0]->eps_B) /
                                n_total +
                            I * 0.0);
    else {
      int k = j;

      ths[0]->b[j] = 0.0 + I * 0.0;

      for (int t = 0; t < ths[0]->d; t++) {
        ths[0]->b[j] += (((double)(k % (ths[0]->n)) / (ths[0]->n) - 0.5) *
                         ((double)(k % (ths[0]->n)) / (ths[0]->n) - 0.5));
        k = k / (ths[0]->n);
      }

      ths[0]->b[j] = ths[1]->b[j] = ths[2]->b[j] = ths[3]->b[j] =
          (double _Complex)(
              regkern3(ths[0]->k, sqrt((double)creal(ths[0]->b[j])), ths[0]->p,
                       ths[0]->kernel_param, ths[0]->eps_I, ths[0]->eps_B) /
                  n_total +
              I * 0.0);
    }
  }

  for (int q = 0; q < 4; q++) {
    nfft_fftshift_complex(ths[q]->b, ths[q]->mv1.d, ths[q]->mv1.N);
    fftw_execute(ths[q]->fft_plan);
    nfft_fftshift_complex(ths[q]->b, ths[q]->mv1.d, ths[q]->mv1.N);
  }

  //   tm = nfft_second( ) - tm;
  //   printf("pre-computation: %f sec\n", tm );

  PARAMS pr[4];
  pthread_t pt[4];

  pr[0].p = 0;
  pr[0].alpha = alpha;
  pr[0].sum = sum;
  pr[0].FS = this;

  pr[1].p = 1;
  pr[1].alpha = alphax;
  pr[1].sum = sumx;
  pr[1].FS = this;

  pr[2].p = 2;
  pr[2].alpha = alphay;
  pr[2].sum = sumy;
  pr[2].FS = this;

  pr[3].p = 3;
  pr[3].alpha = alphaz;
  pr[3].sum = sumz;
  pr[3].FS = this;

  for (int q = 0; q < 4; q++)
    pthread_create(&pt[q], NULL, startFastSumThread, (void *)&pr[q]);

  for (int q = 0; q < 4; q++)
    pthread_join(pt[q], NULL);

  for (int q = 0; q < 4; q++)
    fastsum_finalize(&my_par_fastsum_plan[q]);
}

void FastSummation::fastSum(double *x, double *y, double *alpha, double *alphax,
                            double *alphay, double *alphaz, double *sum,
                            double *sumx, double *sumy, double *sumz) {
  int j, k, t;
  double time;
  fastsum_plan *ths;

  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, c, 0, n, m, p, eps_I,
                    eps_B);
  /*fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, EXACT_NEARFIELD, n,
   * m, p);*/

  printf("pre-computation:    "); // fflush(NULL);
  // time=nfft_second();
  //        fastsum_precompute(&my_fastsum_plan);

  // sort source knots //
  BuildTree(d, 0, x, alpha, alphax, alphay, alphaz, N);

  for (k = 0; k < d * N; k++) {
    my_fastsum_plan.x[k] = x[k];
  }
  for (k = 0; k < M * d; k++) {
    my_fastsum_plan.y[k] = y[k];
  }

  ths = &my_fastsum_plan;

  // precompute spline values for near field //
  if (!(ths->flags & EXACT_NEARFIELD)) {
    if (ths->d == 1)
      for (k = -ths->Ad / 2 - 2; k <= ths->Ad / 2 + 2; k++)
        ths->Add[k + ths->Ad / 2 + 2] =
            regkern1(ths->k, ths->eps_I * (double)k / ths->Ad * 2, ths->p,
                     ths->kernel_param, ths->eps_I, ths->eps_B) +
            I * 0.0;
    else
      for (k = 0; k <= ths->Ad + 2; k++)
        ths->Add[k] = regkern3(ths->k, ths->eps_I * (double)k / ths->Ad, ths->p,
                               ths->kernel_param, ths->eps_I, ths->eps_B) +
                      I * 0.0;
  }

  // init NFFT plan for transposed transform in first step //
  for (k = 0; k < ths->mv1.M_total; k++) //<==> for(k=0; k<N; k++)?
    for (t = 0; t < ths->mv1.d; t++)
      ths->mv1.x[ths->mv1.d * k + t] =
          -ths->x[ths->mv1.d * k + t]; // note the factor -1 for transposed
                                       // transform instead of adjoint //

  // precompute psi, the entries of the matrix B //
  if (ths->mv1.flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv1));

  if (ths->mv1.flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv1));

  if (ths->mv1.flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv1));
  /*
          // init Fourier coefficients //
          for(k=0; k<ths->mv1.M_total;k++)
                  ths->mv1.f[k] = (double _Complex)alpha[k];
  */
  // init NFFT plan for transform in third step //
  for (j = 0; j < ths->mv2.M_total; j++)
    for (t = 0; t < ths->mv2.d; t++)
      ths->mv2.x[ths->mv2.d * j + t] =
          -ths->y[ths->mv2.d * j + t]; // note the factor -1 for conjugated
                                       // transform instead of standard //

  // precompute psi, the entries of the matrix B //
  if (ths->mv2.flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv2));

  if (ths->mv2.flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv2));

  if (ths->mv2.flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv2));

  // precompute Fourier coefficients of regularised kernel //
  int n_total = 1;
  for (t = 0; t < ths->d; t++)
    n_total *= ths->n;

  for (j = 0; j < n_total; j++) {
    if (ths->d == 1)
      ths->b[j] = (double _Complex)(regkern1(ths->k, (double)j / (ths->n) - 0.5,
                                             ths->p, ths->kernel_param,
                                             ths->eps_I, ths->eps_B) /
                                        n_total +
                                    I * 0.0);
    else {
      k = j;
      ths->b[j] = 0.0 + I * 0.0;

      for (t = 0; t < ths->d; t++) {
        ths->b[j] += (((double)(k % (ths->n)) / (ths->n) - 0.5) *
                      ((double)(k % (ths->n)) / (ths->n) - 0.5));
        k = k / (ths->n);
      }
      ths->b[j] = (double _Complex)(
          regkern3(ths->k, sqrt((double)creal(ths->b[j])), ths->p,
                   ths->kernel_param, ths->eps_I, ths->eps_B) /
              n_total +
          I * 0.0);
    }
  }

  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);
  fftw_execute(ths->fft_plan);
  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);

  // time=nfft_second()-time;
  // printf("%fsec\n",time);

  // first summation //

  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alpha[k] + I * 0.0);
  }
  // init Fourier coefficients //
  for (k = 0; k < N; k++) {
    my_fastsum_plan.mv1.f[k] = alpha[k] + I * 0.0;
  }

  printf("fast-computation:   "); // fflush(NULL);
  // time=nfft_second();
  fastsum_trafo(&my_fastsum_plan);
  // time=nfft_second()-time;
  // printf("%fsec\n",time);

  for (k = 0; k < M; k++) {
    sum[k] = (double)creal(my_fastsum_plan.f[k]);
    //              printf("%f\n", sum[k]);
  }

  // x summation //

  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alphax[k] + I * 0.0);
  }
  // init Fourier coefficients //
  for (k = 0; k < N; k++) {
    my_fastsum_plan.mv1.f[k] = alphax[k] + I * 0.0;
  }

  printf("fast computation2:   "); // fflush(NULL);
  // time=nfft_second();
  fastsum_trafo(&my_fastsum_plan);
  // time=nfft_second()-time;
  // printf("%fsec\n",time);

  for (k = 0; k < M; k++) {
    sumx[k] = (double)creal(my_fastsum_plan.f[k]);
    //              printf("%f\n", sumx[k]);
  }

  // y summation //

  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alphay[k] + I * 0.0);
  }
  // init Fourier coefficients //
  for (k = 0; k < N; k++) {
    my_fastsum_plan.mv1.f[k] = alphay[k] + I * 0.0;
  }

  printf("fast computation3:   "); // fflush(NULL);
  // time=nfft_second();
  fastsum_trafo(&my_fastsum_plan);
  // time=nfft_second()-time;
  // printf("%fsec\n",time);

  for (k = 0; k < M; k++) {
    sumy[k] = (double)creal(my_fastsum_plan.f[k]);
    //              printf("%f\n", sumy[k]);
  }

  // z summation //

  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alphaz[k] + I * 0.0);
  }
  // init Fourier coefficients //
  for (k = 0; k < N; k++) {
    my_fastsum_plan.mv1.f[k] = alphaz[k] + I * 0.0;
  }

  printf("fast computation4:   "); // fflush(NULL);
  // time=nfft_second();
  fastsum_trafo(&my_fastsum_plan);
  // time=nfft_second()-time;
  // printf("%fsec\n",time);

  for (k = 0; k < M; k++) {
    sumz[k] = (double)creal(my_fastsum_plan.f[k]);
    //              printf("%f\n", sumz[k]);
  }

  /*
          error=0.0;
          for (j=0; j<my_fastsum_plan.M_total; j++)
          {
                  if
     (cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j])>error)
                          error=cabs(direct[j]-my_fastsum_plan.f[j])/cabs(direct[j]);
          }
          printf("max relative error: %e\n",error);
  */
  fastsum_finalize(&my_fastsum_plan);
}

void FastSummation::fastSum(double *x, double *y, double *alpha, double *sum) {
  int j, k, t;
  double time;
  fastsum_plan *ths;

  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, c, 0, n, m, p, eps_I,
                    eps_B);

  // sort source knots //
  BuildTree(d, 0, x, alpha, N);

  for (k = 0; k < d * N; k++) {
    my_fastsum_plan.x[k] = x[k];
  }
  for (k = 0; k < M * d; k++) {
    my_fastsum_plan.y[k] = y[k];
  }

  ths = &my_fastsum_plan;

  // precompute spline values for near field //
  if (!(ths->flags & EXACT_NEARFIELD)) {
    if (ths->d == 1)
      for (k = -ths->Ad / 2 - 2; k <= ths->Ad / 2 + 2; k++)
        ths->Add[k + ths->Ad / 2 + 2] =
            regkern1(ths->k, ths->eps_I * (double)k / ths->Ad * 2, ths->p,
                     ths->kernel_param, ths->eps_I, ths->eps_B) +
            I * 0.0;
    else
      for (k = 0; k <= ths->Ad + 2; k++)
        ths->Add[k] = regkern3(ths->k, ths->eps_I * (double)k / ths->Ad, ths->p,
                               ths->kernel_param, ths->eps_I, ths->eps_B) +
                      I * 0.0;
  }

  // init NFFT plan for transposed transform in first step //
  for (k = 0; k < ths->mv1.M_total; k++) //<==> for(k=0; k<N; k++)?
    for (t = 0; t < ths->mv1.d; t++)
      ths->mv1.x[ths->mv1.d * k + t] =
          -ths->x[ths->mv1.d * k + t]; // note the factor -1 for transposed
                                       // transform instead of adjoint //

  // precompute psi, the entries of the matrix B //
  if (ths->mv1.flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv1));

  if (ths->mv1.flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv1));

  if (ths->mv1.flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv1));
  /*
          // init Fourier coefficients //
          for(k=0; k<ths->mv1.M_total;k++)
                  ths->mv1.f[k] = (double _Complex)alpha[k];
  */
  // init NFFT plan for transform in third step //
  for (j = 0; j < ths->mv2.M_total; j++)
    for (t = 0; t < ths->mv2.d; t++)
      ths->mv2.x[ths->mv2.d * j + t] =
          -ths->y[ths->mv2.d * j + t]; // note the factor -1 for conjugated
                                       // transform instead of standard //

  // precompute psi, the entries of the matrix B //
  if (ths->mv2.flags & PRE_LIN_PSI)
    nfft_precompute_lin_psi(&(ths->mv2));

  if (ths->mv2.flags & PRE_PSI)
    nfft_precompute_psi(&(ths->mv2));

  if (ths->mv2.flags & PRE_FULL_PSI)
    nfft_precompute_full_psi(&(ths->mv2));

  // precompute Fourier coefficients of regularised kernel //
  int n_total = 1;
  for (t = 0; t < ths->d; t++)
    n_total *= ths->n;

  for (j = 0; j < n_total; j++) {
    if (ths->d == 1) {
      ths->b[j] = (double _Complex)(regkern1(ths->k, (double)j / (ths->n) - 0.5,
                                             ths->p, ths->kernel_param,
                                             ths->eps_I, ths->eps_B) /
                                        n_total +
                                    I * 0.0);
    } else {
      k = j;
      ths->b[j] = 0.0 + I * 0.0;

      for (t = 0; t < ths->d; t++) {
        ths->b[j] += (((double)(k % (ths->n)) / (ths->n) - 0.5) *
                      ((double)(k % (ths->n)) / (ths->n) - 0.5));
        k = k / (ths->n);
      }
      ths->b[j] = (double _Complex)(
          regkern3(ths->k, sqrt((double)creal(ths->b[j])), ths->p,
                   ths->kernel_param, ths->eps_I, ths->eps_B) /
              n_total +
          I * 0.0);
    }
  }

  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);
  fftw_execute(ths->fft_plan);
  nfft_fftshift_complex(ths->b, ths->mv1.d, ths->mv1.N);

  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alpha[k] + I * 0.0);
  }
  // init Fourier coefficients //
  for (k = 0; k < N; k++) {
    my_fastsum_plan.mv1.f[k] = alpha[k] + I * 0.0;
  }

  fastsum_trafo(&my_fastsum_plan);

  for (k = 0; k < M; k++)
    sum[k] = (double)creal(my_fastsum_plan.f[k]);
  fastsum_finalize(&my_fastsum_plan);
}

void FastSummation::directSum(double *x, double *y, double *alpha,
                              double *sum) {
  int k;
  double time;

  fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, c, 0, n, m, p, eps_I,
                    eps_B);
  /*fastsum_init_guru(&my_fastsum_plan, d, N, M, kernel, &c, EXACT_NEARFIELD, n,
   * m, p);*/
  for (k = 0; k < d * N; k++) {
    my_fastsum_plan.x[k] = x[k];
  }
  for (k = 0; k < N; k++) {
    my_fastsum_plan.alpha[k] = (double _Complex)(alpha[k] + I * 0.0);
  }
  for (k = 0; k < M * d; k++) {
    my_fastsum_plan.y[k] = y[k];
  }

  fastsum_exact(&my_fastsum_plan);

  for (k = 0; k < M; k++) {
    sum[k] = (double)creal(my_fastsum_plan.f[k]);
  }
}

/** recursive sort of source knots dimension by dimension to get tree structure
 */
void FastSummation::BuildTree(int d, int t, double *x, double *alpha,
                              double *alphax, double *alphay, double *alphaz,
                              int N) {
  if (N > 1) {
    int m = N / 2;

    quicksort(d, t, x, alpha, alphax, alphay, alphaz, N);

    BuildTree(d, (t + 1) % d, x, alpha, alphax, alphay, alphaz, m);
    BuildTree(d, (t + 1) % d, x + (m + 1) * d, alpha + (m + 1),
              alphax + (m + 1), alphay + (m + 1), alphaz + (m + 1), N - m - 1);
  }
}

void FastSummation::BuildTree(int d, int t, double *x, double *alpha, int N) {
  if (N > 1) {
    int m = N / 2;

    quicksort(d, t, x, alpha, N);

    BuildTree(d, (t + 1) % d, x, alpha, m);
    BuildTree(d, (t + 1) % d, x + (m + 1) * d, alpha + (m + 1), N - m - 1);
  }
}

/** quicksort algorithm for source knots and associated coefficients */
void FastSummation::quicksort(int d, int t, double *x, double *alpha,
                              double *alphax, double *alphay, double *alphaz,
                              int N) {
  int lpos = 0;
  int rpos = N - 1;
  /*double pivot=x[((N-1)/2)*d+t];*/
  double pivot = x[(N / 2) * d + t];

  int k;
  double temp1;
  double temp2;

  while (lpos <= rpos) {
    while (x[lpos * d + t] < pivot)
      lpos++;
    while (x[rpos * d + t] > pivot)
      rpos--;
    if (lpos <= rpos) {
      for (k = 0; k < d; k++) {
        temp1 = x[lpos * d + k];
        x[lpos * d + k] = x[rpos * d + k];
        x[rpos * d + k] = temp1;
      }
      temp2 = alpha[lpos];
      alpha[lpos] = alpha[rpos];
      alpha[rpos] = temp2;

      temp2 = alphax[lpos];
      alphax[lpos] = alphax[rpos];
      alphax[rpos] = temp2;

      temp2 = alphay[lpos];
      alphay[lpos] = alphay[rpos];
      alphay[rpos] = temp2;

      temp2 = alphaz[lpos];
      alphaz[lpos] = alphaz[rpos];
      alphaz[rpos] = temp2;

      lpos++;
      rpos--;
    }
  }
  if (0 < rpos)
    quicksort(d, t, x, alpha, alphax, alphay, alphaz, rpos + 1);
  if (lpos < N - 1)
    quicksort(d, t, x + lpos * d, alpha + lpos, alphax + lpos, alphay + lpos,
              alphaz + lpos, N - lpos);
}

void FastSummation::quicksort(int d, int t, double *x, double *alpha, int N) {
  int lpos = 0;
  int rpos = N - 1;
  double pivot = x[(N / 2) * d + t];

  int k;
  double temp1;
  double temp2;

  while (lpos <= rpos) {
    while (x[lpos * d + t] < pivot)
      lpos++;
    while (x[rpos * d + t] > pivot)
      rpos--;
    if (lpos <= rpos) {
      for (k = 0; k < d; k++) {
        temp1 = x[lpos * d + k];
        x[lpos * d + k] = x[rpos * d + k];
        x[rpos * d + k] = temp1;
      }
      temp2 = alpha[lpos];
      alpha[lpos] = alpha[rpos];
      alpha[rpos] = temp2;

      lpos++;
      rpos--;
    }
  }
  if (0 < rpos)
    quicksort(d, t, x, alpha, rpos + 1);
  if (lpos < N - 1)
    quicksort(d, t, x + lpos * d, alpha + lpos, N - lpos);
}

/** regularized kernel with K_I arbitrary and K_B periodized
 *  (used in 1D)
 */
double FastSummation::regkern1(double _Complex (*kernel)(double, int,
                                                         const double *),
                               double xx, int p, const double *param, double a,
                               double b) {
  int r;
  double sum = 0.0;

  if (xx < -0.5)
    xx = -0.5;
  if (xx > 0.5)
    xx = 0.5;
  if ((xx >= -0.5 + b && xx <= -a) || (xx >= a && xx <= 0.5 - b)) {
    return (double)creal(kernel(xx, 0, param));
  } else if ((xx > -a) && (xx < a)) {
    for (r = 0; r < p; r++) {
      sum +=
          pow(a, (double)r) *
          ((double)creal(kernel(-a, r, param) * BasisPoly(p - 1, r, xx / a)) +
           (double)creal(kernel(a, r, param) * BasisPoly(p - 1, r, -xx / a) *
                         (r & 1 ? -1 : 1)));
    }
    return sum;
  } else if (xx < -0.5 + b) {
    for (r = 0; r < p; r++) {
      sum += pow(b, (double)r) *
             ((double)creal(kernel(0.5 - b, r, param) *
                            BasisPoly(p - 1, r, (xx + 0.5) / b)) +
              (double)creal(kernel(-0.5 + b, r, param) *
                            BasisPoly(p - 1, r, -(xx + 0.5) / b) *
                            (r & 1 ? -1 : 1)));
    }
    return sum;
  } else if (xx > 0.5 - b) {
    for (r = 0; r < p; r++) {
      sum += pow(b, (double)r) *
             ((double)creal(kernel(0.5 - b, r, param) *
                            BasisPoly(p - 1, r, (xx - 0.5) / b)) +
              (double)creal(kernel(-0.5 + b, r, param) *
                            BasisPoly(p - 1, r, -(xx - 0.5) / b) *
                            (r & 1 ? -1 : 1)));
    }
    return sum;
  }
  return (double)creal(kernel(xx, 0, param));
}

/** regularized kernel for even kernels with K_I even
 *  and K_B mirrored smooth to K(1/2) (used in dD, d>1)
 */
double FastSummation::regkern3(double _Complex (*kernel)(double, int,
                                                         const double *),
                               double xx, int p, const double *param, double a,
                               double b) {
  int r;
  double sum = 0.0;

  xx = fabs(xx);

  if (xx >= 0.5) {
    /*return kern(typ,c,0,0.5);*/
    xx = 0.5;
  }
  /* else */
  if ((a <= xx) && (xx <= 0.5 - b)) {
    return (double)creal(kernel(xx, 0, param));
  } else if (xx < a) {
    for (r = 0; r < p; r++) {
      sum += pow(-a, (double)r) * (double)creal(kernel(a, r, param)) *
             (BasisPoly(p - 1, r, xx / a) + BasisPoly(p - 1, r, -xx / a));
    }
    /*sum=kern(typ,c,0,xx); */
    return sum;
  } else if ((0.5 - b < xx) && (xx <= 0.5)) {
    sum = (double)creal(kernel(0.5, 0, param) *
                        BasisPoly(p - 1, 0, -2.0 * xx / b + (1.0 - b) / b));
    /* sum=regkern2(typ,c,p,a,b, 0.5)*BasisPoly(p-1,0,-2.0*xx/b+(1.0-b)/b); */
    for (r = 0; r < p; r++) {
      sum += pow(b / 2.0, (double)r) *
             (double)creal(kernel(0.5 - b, r, param)) *
             BasisPoly(p - 1, r, 2.0 * xx / b - (1.0 - b) / b);
    }
    return sum;
  }
  printf("hier ");
  return 0.0;
}

/** factorial */
double FastSummation::fak(int n) {
  if (n <= 1)
    return 1.0;
  else
    return (double)n * fak(n - 1);
}

/** binomial coefficient */
double FastSummation::binom(int n, int m) {
  return fak(n) / fak(m) / fak(n - m);
}

/** basis polynomial for regularized kernel */
double FastSummation::BasisPoly(int m, int r, double xx) {
  int k;
  double sum = 0.0;

  for (k = 0; k <= m - r; k++) {
    sum += binom(m + k, k) * pow((xx + 1.0) / 2.0, (double)k);
  }
  return sum * pow((xx + 1.0), (double)r) * pow(1.0 - xx, (double)(m + 1)) /
         (1 << (m + 1)) / fak(r); /* 1<<(m+1) = 2^(m+1) */
}
