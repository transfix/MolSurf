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
#ifndef FASTSUMMATION_H
#define FASTSUMMATION_H

//#include <complex>

extern "C" {
  //#include <complex.h>
#include <complex.h>
/* nfft3util.h was removed in NFFT >= 3.3; utility prototypes are in nfft3.h. */
#include "nfft3.h"
#include <FastSummation/fastsum.h>
}

typedef double _Complex (*kernel)(double , int , const double *);

/**
 * Constant symbols
 */
#define EXACT_NEARFIELD  (1U<< 0)

class FastSummation {
	private:
	        int d;		/**< number of dimensions             	*/
        	int N;	        /**< number of source knots         	*/
       		int M;	        /**< number of target knots          	*/
	        int n;		/**< expansion degree                	*/
        	int m;		/**< cut-off parameter for the NFFT  	*/
	        int p;		/**< degree of smoothness of regularization */
	        double _Complex (*kernel)(double , int , const double *);
	        double _Complex (*kernel2)(double , int , const double *);
				/**< kernel function				*/
        	double *c;	/**< parameters for kernel functino		*/
	        double eps_I;	/**< inner boundary                  		*/
	        double eps_B;	/**< outer boundary                  		*/
	        
	        typedef struct
                    {
                     int p;
                     double *alpha, *sum;
                     FastSummation *FS;
                    } PARAMS;

		void BuildTree(int d, int t, double *x, double *alpha, double *alphax, double *alphay, double *alphaz, int N);
		void BuildTree(int d, int t, double *x, double *alpha, int N);
		void quicksort(int d, int t, double *x, double *alpha, double *alphax, double *alphay, double *alphaz, int N);
		void quicksort(int d, int t, double *x, double *alpha, int N);
		double regkern1(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b);
		double regkern3(double _Complex (*kernel)(double , int , const double *), double xx, int p, const double *param, double a, double b);
		double fak(int n);
		double binom(int n, int m);
		double BasisPoly(int m, int r, double xx);
                void fastSumThread( int p, double *alpha, double *sum );
		void fastSumCorrectionThread( int p, double *alpha, double *sum );
                
                static void *startFastSumThread( void *v )
                {
                   PARAMS *pr = ( PARAMS * ) v;
                   
                   pr->FS->fastSumThread( pr->p, pr->alpha, pr->sum );   
                }

                static void *startFastSumCorrectionThread( void *v )
                {
                   PARAMS *pr = ( PARAMS * ) v;

                   pr->FS->fastSumCorrectionThread( pr->p, pr->alpha, pr->sum );
                }

	public:
		FastSummation();
		FastSummation(int dim, int N_source, int M_target, int expansion, int cut_off, int p_degree, char *kernel_func, double *kernel_param, double e_I, double e_B);
		FastSummation(int dim, int N_source, int M_target, int expansion, int cut_off, int p_degree, char *kernel_func1, char *kernel_func2, double *kernel_param, double e_I, double e_B);
		void fastSum(double *x, double *y, double *alpha, double *alphax, double *alphay, double *alphaz, double *sum, double *sumx, double *sumy, double *sumz);
		void fastSum(double *x, double *y, double *alpha, double *sum);
                void threadedFastSum(double *x, double *y, double *alpha, double *alphax, double *alphay, double *alphaz, double *sum, double *sumx, double *sumy, double *sumz);		
                void threadedFastSumCorrection(double *x, double *y, double *alpha, double *alphax, double *alphay, double *alphaz, double *sum, double *sumx, double *sumy, double *sumz, double *c_sum, double *c_sumx, double *c_sumy, double *c_sumz);
		void directSum(double *x, double *y, double *alpha, double *sum);

		fastsum_plan my_fastsum_plan;
		fastsum_plan my_par_fastsum_plan[ 4 ];
		fastsum_plan my_par_corr_fastsum_plan[8];
};


#endif //FASTSUMMATION_H
