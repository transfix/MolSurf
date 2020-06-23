/* find real roots of a polynamial in a given interval using Sturm method */

/* Define the following to use BRENT's algorithm to compute isoloated roots.
   (it's only faster when high precision is requested (<1e-12)) */
#undef BRENT

/* Don't enable NEWTON, it's incomplete!!! */
#undef NEWTON	/* use Newton's method */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __TURBOC__
  #include <bios.h>
#endif

#include <Decimation/Linear/defs.h>

#include <Decimation/Linear/sturm.h>

#ifdef BRENT
#  include "brent.h"
#endif

/* Lower the precision by 8 bits to avoid some troubles caused by
   round-off errors, especially when deciding if a and b are roots */

#define EPS	(256.0*EPSILON)

static double *sturm_rwk; /* contains the polynomials in the Sturm series */
static int sturm_np;	  /* number of polynomials in the Sturm series */
static int *sturm_iwk;	  /* degrees of the polynomials */
static double *sturm_r;	  /* contains the roots found */
static int sturm_nr;	  /* contains the number of roots found */
static double sturm_tol;  /* tolerance */

static int near sturm1(double *p,int n,double a,double b,double tol,
                       double *rwk,int *iwk);
static void near sturm0(double a,double pa,int ia,double b,double pb,int ib);
static void near sturm_factor(double *p,register int n,double a,double *rwk);
static int near sturm_init(double *p,int n,double tol,double *rwk,int *iwk);
static int near sturm_remainder(int n0,double *p0,int n1,double *p1, double *p2);
static double near poly_eval(register int n,double *p,double x);
static int near sturm_ncis(double x,double *px);

int sturm(double *p,int n,double a,double b,double tol,double *rwk,int *iwk)
/*-------------------------------------------------------------------------
  p  	: contains coefficients of the polynomial as input and the roots
	  found as output.
          p(x)=p[0]*x^n+p[1]*x^(n-1)+...+p[n]
  n	: the order of the polynomial.
  a,b	: search for roots in the interval [a,b].
  rwk	: work vector of length >=((n+1)*(n+2)+4)/2.
  iwk	: work vector of length >=n+1.
          It contains the multiplicity of the roots found.

  Returns number of different roots found (-1 means a null polynomial).
-------------------------------------------------------------------------*/
{
  double pa,pb,*p0;
  register int ia,ib,k,nr=0;

#if 0
  		printf("n=%d:",n);
  		for (ia=0; ia<=n; ia++) printf(" %22.16e",p[ia]);
  		printf("\n");
#endif

  p0=p;
  if (tol<EPS) tol=EPS;
  if (a>b) {pa=a; a=b; b=pa;}

  /* normalize the coefficients */
  pa=fabs(p[0]);
  for (ia=1;ia<=n;ia++) {pb=fabs(p[ia]); if (pa<pb) pa=pb;}
  if (pa==0.0) return(-1);
  for (ia=0;ia<=n;ia++) p[ia]/=pa;

  /* remove leading zeros */
  ia=0;
  while (fabs(p[ia])<=EPS && ia<n) ia++;
  if (ia>=n) goto end;
  if (ia) for (ib=ia;ib<=n;ib++) p[ib-ia]=p[ib];
  n-=ia;

  /* check if a is root of p[x] */
  k=0;
  while (n) {
    sturm_factor(p,n,a,rwk);
    pa=0.0; for (ia=0; ia<n; ia++) pa=max(pa,fabs(rwk[ia]));
    if (fabs(rwk[n]/pa)>EPS) break;
    for (ia=0;ia<n;ia++) p[ia]=rwk[ia];
    k++; n--;
  }
  if (k) {
    for (ia=n+1;ia>0;ia--) p[ia]=p[ia-1];
    p[0]=a; iwk[0]=k; p++; iwk++; nr++;
  }
  if (!n) goto end;
  if (a==b) goto end;

  /* check if b is root of p[x] */
  k=0;
  while (n) {
    sturm_factor(p,n,b,rwk);
    pa=0.0; for (ia=0; ia<n; ia++) pa=max(pa,fabs(rwk[ia]));
    if (fabs(rwk[n]/pa)>EPS) break;
    for (ia=0;ia<n;ia++) p[ia]=rwk[ia];
    k++; n--;
  }
  if (k) {
    for (ia=n+1;ia>0;ia--) p[ia]=p[ia-1];
    p[0]=b; iwk[0]=k; p++; iwk++; nr++;
  }
  if (!n) goto end;

  nr+=sturm1(p,n,a-EPS,b+EPS,tol,rwk,iwk);

end:
  for (ia=0; ia<nr; ia++) {
    if (p0[ia]<a && p0[ia]>=a-EPS) p0[ia]=a;
    else if (p0[ia]>b && p0[ia]<=b+EPS) p0[ia]=b;
  }

#if 0
  		printf("nr=%d:",nr);
  		for (ia=0; ia<nr; ia++) printf(" %22.16e",p0[ia]);
  		printf("\n");
#endif
  return nr;
}

static int near sturm1(double *p,int n,double a,double b,double tol,
                       double *rwk,int *iwk)
{
  double pa,pb;
  register int ia,ib,nr;

  nr=0;

  if (!n) return(0);
  else if (n==1) {		/* linear equation */
    p[0]=-p[1]/p[0]; iwk[0]=1;
check:
    return(((p[0]>=a && p[0]<=b) ? 1:0));
  } else if (n==2) {		/* quadratique equation */
    pa=p[1]*p[1]-4.0*p[0]*p[2];
    if (pa<-EPS) return(0);
    if (pa<=EPS) {		/* two equal roots */
      p[0]=-p[1]/(p[0]+p[0]); iwk[0]=2;
      goto check;
    } else {			/* two different roots */
      pa=sqrt(pa);
      pb=p[0]+p[0];
      p[0]=(-p[1]+pa)/pb; p[1]=(-p[1]-pa)/pb;
      iwk[0]=iwk[1]=1;
      ia=0;
      if (p[0]>=a && p[0]<=b) ia++; else p[0]=p[1];
      if (p[ia]>=a && p[ia]<=b) ia++;
      return(ia);
    }
  }

  if (sturm_init(p,n,tol,rwk,iwk)) return(0);

  /* sturm_np must be greater than or equal to 2 */
  if (sturm_np==2) {
    /* p'(x) is factor of p(x) ==> p(x)=p[0]*( x + (p[n]/p[0])^{1/n} )^n */
    pa=-p[n]/p[0];
    if (pa<0.0) {
      if (pa>=-EPS) {
        pa=0.0;
      } else if (n&1) {
        pa=-exp(log(-pa)/n);
      } else {
        printf("Unknown error in sturm()!\n"); bioskey(0);
        return(-1);
      }
    } else if (pa>0) pa=exp(log(pa)/n);
    p[0]=pa; iwk[0]=n;
    return(pa>=a && pa<=b ? 1:0);
  } else {
    /* find ALL roots of the last polynomial of the Sturm sequence */
    int n0,n1;

    ib=0;
    for (ia=0;ia<sturm_np-1;ia++) ib+=iwk[ia]+1;
    n0=iwk[ia];
    if (n0) {
      int j,k;
      double d;

      for (ia=0;ia<=n0;ia++) {
        rwk[ia]=rwk[ib+ia];
      }
      /* find lower and upper bounds for rwk(x) */
      pa=0.0; pb=0.0;
      for (k=1;k<=n0;k++) {
        d=fabs(rwk[k]);
        pa=max(pa,d); pb+=d;
      }
      d=fabs(rwk[0]);
      pa=1.0+pa/d; pb=max(1.0,pb/d);
      pa=0.0001+min(pa,pb);
      n1=sturm1(rwk,n0,-pa,pa,EPS,rwk+n0+1,iwk);
      /* store roots found with multiplicity increased by 1,
         and eliminate corresponding factors from p(x) */
      for (ib=0;ib<n1;ib++) {
        pa=rwk[ib]; j=iwk[ib];
        for (ia=0;ia<=j;ia++) {
          sturm_factor(p,n,pa,rwk+n1);
          for (k=0;k<n;k++) p[k]=rwk[n1+k];
          n--;
        }
        if (pa>=a && pa<=b) {
          for (ia=n+1;ia>0;ia--) p[ia]=p[ia-1];
          *(p++)=pa; iwk[nr++]=j+1;
        }
      } /* end of "for (ib=0;..." */
      iwk+=nr;
      /* reinit the arrays (this will also restore sturm_rwk, etc.) */
      if (sturm_init(p,n,tol,rwk,iwk)) return(nr);
    }  /* end of "if (n0) ..." */
  }

  /* at this point, p(x) can only have simple roots */
  ia=sturm_ncis(a,&pa), ib=sturm_ncis(b,&pb);
  sturm_nr=0;
  sturm0(a,pa,ia,b,pb,ib);
  for (ia=0;ia<sturm_nr;ia++) iwk[ia]=1;
  return(nr+sturm_nr);
}

static void near sturm_factor(double *p,register int n,double a,double *rwk)
/* divide the polynomial p[x] by (x-a), return result in rwk,
   rwk[n]=remainder. */
{
  register int i;

  rwk[0]=p[0];
  for (i=1;i<=n;i++) rwk[i]=p[i]+a*rwk[i-1];
}

static int near sturm_init(double *p,int n,double tol,double *rwk,int *iwk)
/* returns !0 if p(x) is a constant polynomial */
{
  register int i;
  int p0,p1,p2,n0,n1,n2;

  p0=0;
  /* p0(x)=p(x) */
  for (i=0;i<=n;i++) rwk[i]=p[i];
  iwk[0]=n0=n; p1=n0+1;
  /* p1(x)=p'(x) */
  for (i=0;i<n;i++) rwk[p1+i]=rwk[p0+i]*(n-i);
  iwk[1]=n1=n-1; p2=p1+n1+1;

  sturm_rwk=rwk; sturm_iwk=iwk; sturm_r=p; sturm_tol=tol;
  sturm_np=2;
  /* compute p_{k+2}(x): p_k(x)=q_k(x)*p_{k+1}(x)-p_{k+2}(x), k>=0 */
  do {
    n2=sturm_remainder(n0,rwk+p0,n1,rwk+p1,rwk+p2);
    if (n2<0) break;
    p0=p1; p1=p2; p2+=n2+1;
    n0=n1; n1=n2;
    iwk[sturm_np++]=n2;
  } while (1);

  return(0);
}

static int near sturm_remainder(int n0,double *p0,int n1,double *p1, double *p2)
/* compute p2(x) such that p0(x)=q(x)*p1(x)-p2(x). returns the degree of
   p2(x) (returns -1 if p2(x)==0). The highest-order coefficients of p0
   and p1 are supposed nonzero */
{
  register int i,j,n2;
  double d;

  if (!n0) return(-1);

  for (i=0;i<=n0;i++) p2[i]=p0[i];

  n2=n0; j=0;
  while (n2>=n1) {
    d=p2[j]/p1[0];
    for (i=1;i<=n1;i++) p2[j+i]-=p1[i]*d;
    j++; n2--;
    while (n2>=0 && fabs(p2[j])<=EPS) {
      n2--; j++;
    }
  }
  if (n2>=0) for (i=0;i<=n2;i++) p2[i]=-p2[j+i];
  return(n2);
}

static double near poly_eval(register int n,double *p,register double x)
/* evaluate the value of p(x) */
{
  register int i;
  register double d;

  d=p[0];
  for (i=1;i<=n;i++) d=d*x+p[i];
  return(d);
}

static int near sturm_ncis(double x,double *px)
/* returns the number of changes in sign. *px contains p0(x) */
{
  register int i,s;
  double d0,d1,*p;

  *px=d0=poly_eval(sturm_iwk[0],p=sturm_rwk,x);
  p+=sturm_iwk[0]+1;
  i=1;
  while (d0==0.0 && i<sturm_np) {
    d0=poly_eval(sturm_iwk[i],p,x);
    p+=sturm_iwk[i]+1;
    i++;
  }
  s=0;
  while (i<sturm_np) {
    d1=poly_eval(sturm_iwk[i],p,x);
    p+=sturm_iwk[i]+1;
    i++;
    if (d1==0.0) continue;
    if ((d0>0.0)!=(d1>0.0)) s++;
    d0=d1;
  }
  return(s);
}

#ifndef BRENT	/* =========================================================*/

#ifdef NEWTON
static int near sturm_newton(int n, double *p, double *pprime,
                             double *x, double p0, double a, double b)
/* iterates using Newton's method. returns !0 if converged */
{
  double x0, pp;
  int it=0, maxit=10;

  do {
    if (it>=maxit) return 0;
		printf("it=%d, x=%lf, p(x)=%lf\n", it, *x, p0);
    pp=poly_eval(n-1,pprime,x0=*x);
    if (pp==0.0) {
    		printf("derivative null!\n");
      return 0;
    }

    *x=x0-p0/pp;

#if 0
    if (*x<a || *x>b) return 0;
#endif

    if ( fabs(*x-x0)<=sturm_tol*(fabs(*x)+fabs(x0))) return 1;
    pp=poly_eval(n,p,*x);

#if 1
    if (fabs(pp)>=fabs(p0)) {
    		printf("divergence!\n");
      return 0;
    }
#endif

    p0=pp; it++;
  } while (1);
}
#endif

static void near sturm_bisec(double a,double pa,double b,double pb)
/* compute the root in the interval (a,b) to desired precision using
   a combination of bisection and Newton's method. */
{
  double c;
  register int n;
#ifdef NEWTON
  double *pprime;	/* first derivative of the polynomial */
  double x;
  int i;
#endif

  n=sturm_iwk[0];

#ifdef NEWTON
  pprime=malloc(n*sizeof(double));
  if (pprime!=NULL) {
    for (i=0; i<n; i++) pprime[i]=(i+1)*sturm_rwk[i+1];
  }
#endif

  do {
    c=(a+b)*0.5;
    if ((b-a)<=sturm_tol*(fabs(a)+fabs(b))) {
end:
#if 0
      // arand commenting to remove warnings...
      //for (n=0;n<sturm_nr;n++) if (fabs(c-sturm_r[n])<1e-8???) return;
#endif
      sturm_r[sturm_nr++]=c;
#ifdef NEWTON
      if (pprime!=NULL) free(pprime);
#endif
      return;
    }
    pb=poly_eval(n,sturm_rwk,c);
    if (pb==0.0) goto end;

#ifdef NEWTON
    /* try Newton's method */
    x=c;
    if (sturm_newton(n, sturm_rwk, pprime, &x, pb, a, b)) {
      c=x;
      goto end;
    }
#endif

    if ((pb>0.0)==(pa>0.0)) {
      a=c; pa=pb;
    } else b=c;

  } while (1);
}

#else		/* =========================================================*/

double sturm_func(register double x)
{
  return poly_eval(sturm_iwk[0],sturm_rwk,x);
}

static void near sturm_bisec(double a,double pa,double b,double pb)
{
  sturm_r[sturm_nr++]=zeroin(a, b, pa, pb, sturm_func, sturm_tol);
}

#endif		/* =========================================================*/

static void near sturm0(double a,double pa,int ia,double b,double pb,int ib)
/* locate recursively real roots of p0(x) */
{
  double c,pc;
  register int ic;

  if (ia<=ib) {
    return;
  }
  if (ia==ib+1 && (pa>0.0)!=(pb>0.0)) {
    sturm_bisec(a,pa,b,pb);
  } else {
    c=(a+b)*0.5;
    if ((b-a)<=sturm_tol*(fabs(a)+fabs(b))) {
#if 0
      // arand commenting to remove warnings...
      //for (ic=0;ic<sturm_nr;i++) if (fabs(c-sturm_r[ic])<=1e-8???) return;
#endif
      sturm_r[sturm_nr]=c;
      return;
    }
    ic=sturm_ncis(c,&pc);
    sturm0(a,pa,ia,c,pc,ic);
    sturm0(c,pc,ic,b,pb,ib);
  }
}
