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
/*! \file kernels.h
 *  \brief Header file with predefined kernels for the fast summation algorithm.
 */
#ifndef KERNELS_H
#define KERNELS_H

/** 
 * \addtogroup applications_fastsum
 * \{
 */

double _Complex gaussian(double x, int der, const double *param);              /* K(x)=exp(-x^2/c^2) */
double _Complex multiquadric(double x, int der, const double *param);          /* K(x)=sqrt(x^2+c^2) */
double _Complex inverse_multiquadric(double x, int der, const double *param);  /* K(x)=1/sqrt(x^2+c^2) */
double _Complex logarithm(double x, int der, const double *param);             /* K(x)=log |x| */
double _Complex thinplate_spline(double x, int der, const double *param);      /* K(x) = x^2 log |x| */
double _Complex one_over_square(double x, int der, const double *param);       /* K(x) = 1/x^2 */
double _Complex one_over_modulus(double x, int der, const double *param);      /* K(x) = 1/|x| */
double _Complex one_over_x(double x, int der, const double *param);            /* K(x) = 1/x */
double _Complex inverse_multiquadric3(double x, int der, const double *param); /* K(x) = 1/sqrt(x^2+c^2)^3 */
double _Complex sinc_kernel(double x, int der, const double *param);           /* K(x) = sin(cx)/x */
double _Complex cosc(double x, int der, const double *param);                  /* K(x) = cos(cx)/x */
double _Complex cot(double x, int der, const double *param);                   /* K(x) = cot(cx) */
double _Complex one_over_cube(double x, int der, const double *param);                /* K(x) = 1/x^3 */
double _Complex one_over_quartic(double x, int der, const double *param);      /* K(x) = 1/x^4 */
double _Complex one_over_fifth(double x, int der, const double *param);        /* K(x) = 1/|x|^5 */
double _Complex one_over_sixth(double x, int der, const double *param);        /* K(x) = 1/|x|^6 */
double _Complex one_over_seventh(double x, int der, const double *param);      /* K(x) = 1/|x|^7 */
double _Complex truncated_one_over_quartic(double x, int der, const double *param);
double _Complex still_gb(double x, int der, const double *param);	      /* K(x) = 1/sqrt(c1 x*2 + c2 exp(c3 x^2)) */
/* \} */

#endif
/* kernels.h */
