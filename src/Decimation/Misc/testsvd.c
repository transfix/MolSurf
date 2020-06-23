#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Decimation/Misc/misc.h>

/* Test svd code.    all the matrix are column-wise stored  */
/* ldx, ldu,ldv are leading dimension of x,u,v, leading dimension 
   is the number of rows                                      */
/* The result is x = u s v^T                                */
/* x is destroyed after the calling                         */
void TestSvd()
{
  double w, x[1000000], u[1000000], v[1000000], s[1000], e[1000], work[1000], mid[1000000];
  int    n, p, ldx, ldu, ldv, job, info, i,j,k,m;
  double   maxw, minw;
  
  n = 1000;               /* the number of rows   */
  p = 1000;               /* the number of column */
  ldx = n;                /* ldx >= n             */
  ldu = n;                /* ldu >= n             */
  ldv =p;                 /* ldv >= p             */
  job = 11;
  
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      x[i + j*ldx] = 1.0/(i+j+1.0);
      mid[i*n + j] = 0.0;
    }
    x[i + i*ldx] =  x[i + i*ldx] + 0.0;
  }
  
  dsvdc_ (x, &ldx, &n, &p, s, e, u, &ldu, v, &ldv, work, &job, &info);
  
  printf("test svd info = %d \n", info);
  for (i = 0; i <n; i++) {
    for (j = 0; j <n ; j++) {
      printf("U[%d,%d] = %f\n",i,j, u[i + j*ldu]);
    }
  }
  
  for (i = 0; i <p; i++) {
    for (j = 0; j <p ; j++) {
      printf("V[%d, %d] = %f\n", i,j, v[i + j*ldv]);
    }
  }
  
  
  m = n;
  if (n > p) m = p;
  
  /* recover the matrix */
  for (i = 0; i <m ; i++) {
    printf("S = %f, E = %f\n", s[i], e[i]);
    mid[i*n + i] = s[i];
  }
  
  /* U*S    */
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      w = 0.0; 
      for (k = 0; k < n; k++) {
	w = w + u[i+ldu*k] * mid[k + j*n];
      } 
      x[i + j*n] = w;
    }
  }
  
  maxw = 0.0, minw = 100.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      w = 0.0; 
      for (k = 0; k < p; k++) {
	w = w + x[i+n*k] * v[k*ldv + j];
      }
      
      w = w - 1.0/(i+j+1.0);
      w = w*(i+j+1.0);
      if (maxw < w) maxw = w;
      if (minw > w) minw = w;
      
      printf("matrix[%d,%d] = %f\n", i,j, w - (1.0/(i+j+1.0))); 
    }
  }
  
  printf("max error = %e, %e\n", maxw, minw);
  
  /* U orthgonal xing */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      w = 0.0;
      for (k = 0; k < n; k++) {
	w = w + u[i+ldu*k] * u[k*ldu + j];
      }
      
      printf("U*U^T[%d,%d] = %f\n", i,j, w);
    }
  }
  
  /* V orthgonal xing */
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      w = 0.0;
      for (k = 0; k < p; k++) {
	w = w + v[i+ldv*k] * v[k*ldv + j];
      }
      printf("V*V^T[%d,%d] = %f\n", i,j, w);
    }
  }
} 
