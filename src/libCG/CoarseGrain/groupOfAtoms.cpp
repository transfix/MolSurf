#include <Utility/utility.h>
#include <libCG/CoarseGrain/groupOfAtoms.h>
#include <libCG/CoarseGrain/atom2.h>
#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>
#include <FastSummation/FastSummation.h>

//#define _TRIVIAL_
#define _FASTSUM_
#define _CORRECTION_

using namespace MOLECULE;
using namespace GEOMETRY;

GroupOfAtoms::GroupOfAtoms()
{
}

GroupOfAtoms::GroupOfAtoms(GroupOfAtoms* m)
{
	int i;
	for(i=0; i<m->numOfAtoms; i++)
	{
		addAtom(m->m_Atoms[i]);
	}
	numOfAtoms = m->numOfAtoms;
	numOfRes = m->numOfRes;
}

GroupOfAtoms::~GroupOfAtoms()
{
	std::vector<Atom*>::iterator iter = m_Atoms.begin(), end = m_Atoms.end();
	for(; iter != end; ++iter)
	{
		delete *iter;
	}
	m_Atoms.clear();
}

void GroupOfAtoms::addAtom(Atom* currentAtom)
{
	Atom* a = new Atom(currentAtom);
	m_Atoms.push_back(a);
}

void GroupOfAtoms::getBornRadii(Surface* surface)
{
	int i,k;
	double* x, *y, *coeff, *coeffx, *coeffy, *coeffz;
	double* sum, *sumx, *sumy, *sumz, *c_sum, *c_sumx, *c_sumy, *c_sumz;
	double A4, A7;
	FastSummation* m_fastsum;
	const double pi = 3.1415926;
	int nn = surface->nodes_per_face;
	int N = surface->numbtris*nn;
	int M = numOfAtoms;
	int n = 30;	//100;
	int p = 2;
	int m = 4;
	double param = 2;
	double eps_I = 0.1; //0.01; // for x_4th and x_7th
	double eps_B = 0.0; // 0.01;
	x = (double*)malloc(sizeof(double)*3*N);
	y = (double*)malloc(sizeof(double)*3*M);
	coeff = (double*)malloc(sizeof(double)*N);
	coeffx = (double*)malloc(sizeof(double)*N);
	coeffy = (double*)malloc(sizeof(double)*N);
	coeffz = (double*)malloc(sizeof(double)*N);
	sum = (double*)malloc(sizeof(double)*M);
	sumx = (double*)malloc(sizeof(double)*M);
	sumy = (double*)malloc(sizeof(double)*M);
	sumz = (double*)malloc(sizeof(double)*M);
#if defined( _TRIVIAL_) && defined(_CORRECTION_)
	/*
	        double *t_sum = (double *)malloc(sizeof(double)*M);
	        double *t_sumx = (double *)malloc(sizeof(double)*M);
	        double *t_sumy = (double *)malloc(sizeof(double)*M);
	        double *t_sumz = (double *)malloc(sizeof(double)*M);
	*/
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			x[3*(i*nn+k)] = surface->m_Faces[i]->node[k].Point[0];
			x[3*(i*nn+k)+1] = surface->m_Faces[i]->node[k].Point[1];
			x[3*(i*nn+k)+2] = surface->m_Faces[i]->node[k].Point[2];
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		y[3*i  ] = m_Atoms[i]->position[0];
		y[3*i+1] = m_Atoms[i]->position[1];
		y[3*i+2] = m_Atoms[i]->position[2];
	}
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			coeff[i* nn+k] = dotProduct(surface->m_Faces[i]->node[k].Point, surface->m_Faces[i]->node[k].Normal);
			coeff[i* nn+k] *= surface->m_Faces[i]->weight[k];
			coeffx[i* nn+k] = surface->m_Faces[i]->node[k].Normal[0];
			coeffx[i* nn+k] *= surface->m_Faces[i]->weight[k];
			coeffy[i* nn+k] = surface->m_Faces[i]->node[k].Normal[1];
			coeffy[i* nn+k] *= surface->m_Faces[i]->weight[k];
			coeffz[i* nn+k] = surface->m_Faces[i]->node[k].Normal[2];
			coeffz[i* nn+k] *= surface->m_Faces[i]->weight[k];
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		sum[i] = 0;
		sumx[i] = 0;
		sumy[i] = 0;
		sumz[i] = 0;
		for(k = 0; k < surface->numbtris*nn; k++)
		{
			double denom = (x[3*k  ] - y[3*i  ])*(x[3*k  ] - y[3*i  ]) +
						   (x[3*k+1] - y[3*i+1])*(x[3*k+1] - y[3*i+1]) +
						   (x[3*k+2] - y[3*i+2])*(x[3*k+2] - y[3*i+2]);
			//if (denom > 5*5)	continue;
			denom *= denom;
			sum[i]  += coeff[k] / denom;
			sumx[i] += coeffx[k] / denom;
			sumy[i] += coeffy[k] / denom;
			sumz[i] += coeffz[k] / denom;
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = (sum[i] - m_Atoms[i]->position[0]*sumx[i] - m_Atoms[i]->position[1]*sumy[i] - m_Atoms[i]->position[2]*sumz[i])/(4*pi);
		(A4 > 0.0) ? m_Atoms[i]->Born_radius = A4 : m_Atoms[i]->Born_radius = 0.0;
		if(A4 < 0.0)
		{
			printf("negative\n");
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = m_Atoms[i]->Born_radius;
		(A4 != 0) ? m_Atoms[i]->Born_radius = 1.0/A4 : m_Atoms[i]->Born_radius=100000.0;
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		sum[i] = 0;
		sumx[i] = 0;
		sumy[i] = 0;
		sumz[i] = 0;
		for(k = 0; k < surface->numbtris*nn; k++)
		{
			double denom = (x[3*k  ] - y[3*i  ])*(x[3*k  ] - y[3*i  ]) +
						   (x[3*k+1] - y[3*i+1])*(x[3*k+1] - y[3*i+1]) +
						   (x[3*k+2] - y[3*i+2])*(x[3*k+2] - y[3*i+2]);
			/* Version before 10/18/08
			                        denom = denom*denom*denom;
			                        sum[i]  += coeff[k] / (denom*sqrt(denom));
			                        sumx[i] += coeffx[k] / (denom*sqrt(denom));
			                        sumy[i] += coeffy[k] / (denom*sqrt(denom));
			                        sumz[i] += coeffz[k] / (denom*sqrt(denom));
			*/
			sum[i]  += coeff[k] / (denom*denom*denom*sqrt(denom));
			sumx[i] += coeffx[k] / (denom*denom*denom*sqrt(denom));
			sumy[i] += coeffy[k] / (denom*denom*denom*sqrt(denom));
			sumz[i] += coeffz[k] / (denom*denom*denom*sqrt(denom));
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = 1.0 / m_Atoms[i]->Born_radius;
		A7 = (sum[i]-m_Atoms[i]->position[0]*sumx[i]-m_Atoms[i]->position[1]*sumy[i]-m_Atoms[i]->position[2]*sumz[i])/(16*pi);
		//		printf("%f %f\n", A4, A7);
		if(A7 < 0.0)
		{
			A7 = 0.0;
		}
		//                A7 = sqrt(A7); //A7 = sqrt(A7); // wenqi added on 12/17
		A7 = sqrt(A7);
		A7 = sqrt(A7); // wenqi added on 10/18/08
		if(A4 + A7 > 0.0000001)
		{
			m_Atoms[i]->Born_radius = 1.0 / ((1.0-1.0/sqrt(2.0)) * A4 + A7);    // this is GBSW model
		}
		else
		{
			m_Atoms[i]->Born_radius = 100000.0;
		}
		if(m_Atoms[i]->Born_radius < 0.0)
		{
			printf("negative %d %f\n", i, m_Atoms[i]->Born_radius);
		}
	}
#endif // _TRIVIAL_ and _CORRECTION__
#if defined(_FASTSUM_) && !defined(_CORRECTION_)
	// make x,y within [-1/4, 1/4]
	double abs_max, abs_i;
	abs_max = -1.0;
	for(i = 0; i < 3; i++)
	{
		abs_i = surface->max[i]-surface->min[i];
		if(abs_i > abs_max)
		{
			abs_max = abs_i;
		}
		/*		abs_i = fabs(surface->min[i]);
				if ( abs_i > abs_max)	abs_max = abs_i;
				abs_i = fabs(surface->max[i]);
		                if ( abs_i > abs_max)   abs_max = abs_i;
		*/
	}
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			x[3*(i*nn+k)] = (surface->m_Faces[i]->node[k].Point[0]-surface->center[0]) / (2*abs_max);
			x[3*(i*nn+k)+1] = (surface->m_Faces[i]->node[k].Point[1]-surface->center[1]) / (2*abs_max);
			x[3*(i*nn+k)+2] = (surface->m_Faces[i]->node[k].Point[2]-surface->center[2]) / (2*abs_max);
			if(fabs(x[3*(i*nn+k)]) > 0.25 || fabs(x[3*(i*nn+k)+1]) > 0.25 || fabs(x[3*(i*nn+k)+2]) > 0.25)
			{
				printf("x out of range\n");
			}
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		y[3*i  ] = (m_Atoms[i]->position[0]-surface->center[0]) / (2*abs_max);
		y[3*i+1] = (m_Atoms[i]->position[1]-surface->center[1]) / (2*abs_max);
		y[3*i+2] = (m_Atoms[i]->position[2]-surface->center[2]) / (2*abs_max);
		if(fabs(y[3*i]) > 0.25 | fabs(y[3*i+1]) > 0.25 | fabs(y[3*i+2]) > 0.25)
		{
			printf("y out of range\n");
		}
	}
	// first summation //
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			coeff[i* nn+k] = dotProduct(surface->m_Faces[i]->node[k].Point, surface->m_Faces[i]->node[k].Normal);
			coeff[i* nn+k] *= surface->m_Faces[i]->weight[k] / pow(2*abs_max, 4);
			coeffx[i* nn+k] = surface->m_Faces[i]->node[k].Normal[0];
			coeffx[i* nn+k] *= surface->m_Faces[i]->weight[k] / pow(2*abs_max, 4);
			coeffy[i* nn+k] = surface->m_Faces[i]->node[k].Normal[1];
			coeffy[i* nn+k] *= surface->m_Faces[i]->weight[k] / pow(2*abs_max, 4);
			coeffz[i* nn+k] = surface->m_Faces[i]->node[k].Normal[2];
			coeffz[i* nn+k] *= surface->m_Faces[i]->weight[k] / pow(2*abs_max, 4);
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		sum[i] = 0.0;
		sumx[i] = 0.0;
		sumy[i] = 0.0;
		sumz[i] = 0.0;
	}
	param /= (1*abs_max);
	m_fastsum = new FastSummation(3, N, M, n, m, p, "one_over_quartic", &param, eps_I, eps_B);
	m_fastsum->threadedFastSum(x, y, coeff, coeffx, coeffy, coeffz, sum, sumx, sumy, sumz);
	delete m_fastsum;
	/*
		printf("compare\n");
		for (i = 0; i < 20; i++)
		{
			printf("%f %f %f %f\n", t_sumy[i], sumy[i], t_sumz[i], sumz[i]);
		}
	*/
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = (sum[i] - m_Atoms[i]->position[0]*sumx[i] - m_Atoms[i]->position[1]*sumy[i] - m_Atoms[i]->position[2]*sumz[i])/(4*pi);
		(A4 > 0.0) ? m_Atoms[i]->Born_radius = A4 : m_Atoms[i]->Born_radius = 0.0;
		if(A4 < 0.0)
		{
			printf("negative\n");
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = m_Atoms[i]->Born_radius;
		(A4 != 0) ? m_Atoms[i]->Born_radius = 1.0/A4 : m_Atoms[i]->Born_radius=100000.0;
	}
#endif // _FASTSUM_
	// with correction
#if defined(_FASTSUM_) && defined(_CORRECTION_)
	c_sum = (double*)malloc(sizeof(double)*M);
	c_sumx = (double*)malloc(sizeof(double)*M);
	c_sumy = (double*)malloc(sizeof(double)*M);
	c_sumz = (double*)malloc(sizeof(double)*M);
	/*
	#ifdef _TRIVIAL_
	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
	                        x[3*(i*nn+k)  ] = surface->m_Faces[i]->node[k].Point[0];
	                        x[3*(i*nn+k)+1] = surface->m_Faces[i]->node[k].Point[1];
	                        x[3*(i*nn+k)+2] = surface->m_Faces[i]->node[k].Point[2];
	                }
	        }

	        for ( i = 0; i < numOfAtoms; i++)
	        {
	                y[3*i  ] = m_Atoms[i]->position[0];
	                y[3*i+1] = m_Atoms[i]->position[1];
	                y[3*i+2] = m_Atoms[i]->position[2];
	        }

	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
	                        coeff[i*nn+k] = dotProduct(surface->m_Faces[i]->node[k].Point, surface->m_Faces[i]->node[k].Normal);
	                        coeff[i*nn+k] *= surface->m_Faces[i]->weight[k];

	                        coeffx[i*nn+k] = surface->m_Faces[i]->node[k].Normal[0];
	                        coeffx[i*nn+k] *= surface->m_Faces[i]->weight[k];

	                        coeffy[i*nn+k] = surface->m_Faces[i]->node[k].Normal[1];
	                        coeffy[i*nn+k] *= surface->m_Faces[i]->weight[k];

	                        coeffz[i*nn+k] = surface->m_Faces[i]->node[k].Normal[2];
	                        coeffz[i*nn+k] *= surface->m_Faces[i]->weight[k];
	                }
	        }

	        for (i = 0; i < numOfAtoms; i++)
	        {
	                t_sum[i] = 0;     t_sumx[i] = 0;    t_sumy[i] = 0;    t_sumz[i] = 0;
	                for (k = 0; k < surface->numbtris*nn; k++)
	                {
	                        double denom = (x[3*k  ] - y[3*i  ])*(x[3*k  ] - y[3*i  ]) +
	                                       (x[3*k+1] - y[3*i+1])*(x[3*k+1] - y[3*i+1]) +
	                                       (x[3*k+2] - y[3*i+2])*(x[3*k+2] - y[3*i+2]);
	                        //if (denom > 5*5)      continue;
	                        denom *= denom;
	                        t_sum[i]  += coeff[k] / (denom*denom*denom*sqrt(denom));
	                        t_sumx[i] += coeffx[k] / (denom*denom*denom*sqrt(denom));
	                        t_sumy[i] += coeffy[k] / (denom*denom*denom*sqrt(denom));
	                        t_sumz[i] += coeffz[k] / (denom*denom*denom*sqrt(denom));
	                }
	        }
	#endif // _TRIVIAL_
	*/
	// make x,y within [-1/4, 1/4]
	double abs_max, abs_i;
	abs_max = -1.0;
	for(i = 0; i < 3; i++)
	{
		abs_i = surface->max[i]-surface->min[i];
		if(abs_i > abs_max)
		{
			abs_max = abs_i;
		}
	}
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			x[3*(i*nn+k)] = (surface->m_Faces[i]->node[k].Point[0]-surface->center[0]) / (2*abs_max);
			x[3*(i*nn+k)+1] = (surface->m_Faces[i]->node[k].Point[1]-surface->center[1]) / (2*abs_max);
			x[3*(i*nn+k)+2] = (surface->m_Faces[i]->node[k].Point[2]-surface->center[2]) / (2*abs_max);
			if(fabs(x[3*(i*nn+k)]) > 0.25 || fabs(x[3*(i*nn+k)+1]) > 0.25 || fabs(x[3*(i*nn+k)+2]) > 0.25)
			{
				printf("x out of range\n");
			}
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		y[3*i  ] = (m_Atoms[i]->position[0]-surface->center[0]) / (2*abs_max);
		y[3*i+1] = (m_Atoms[i]->position[1]-surface->center[1]) / (2*abs_max);
		y[3*i+2] = (m_Atoms[i]->position[2]-surface->center[2]) / (2*abs_max);
		if(fabs(y[3*i]) > 0.25 | fabs(y[3*i+1]) > 0.25 | fabs(y[3*i+2]) > 0.25)
		{
			printf("y out of range\n");
		}
	}
	// first summation //
	for(i = 0; i < surface->numbtris; i++)
	{
		for(k = 0; k < nn; k++)
		{
			coeff[i* nn+k] = dotProduct(surface->m_Faces[i]->node[k].Point, surface->m_Faces[i]->node[k].Normal);
			coeff[i* nn+k] *= surface->m_Faces[i]->weight[k]; // / pow(2*abs_max, 4);
			coeffx[i* nn+k] = surface->m_Faces[i]->node[k].Normal[0];
			coeffx[i* nn+k] *= surface->m_Faces[i]->weight[k]; // / pow(2*abs_max, 4);
			coeffy[i* nn+k] = surface->m_Faces[i]->node[k].Normal[1];
			coeffy[i* nn+k] *= surface->m_Faces[i]->weight[k]; // / pow(2*abs_max, 4);
			coeffz[i* nn+k] = surface->m_Faces[i]->node[k].Normal[2];
			coeffz[i* nn+k] *= surface->m_Faces[i]->weight[k]; // / pow(2*abs_max, 4);
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		sum[i] = 0.0;
		sumx[i] = 0.0;
		sumy[i] = 0.0;
		sumz[i] = 0.0;
		c_sum[i] = 0.0;
		c_sumx[i] = 0.0;
		c_sumy[i] = 0.0;
		c_sumz[i] = 0.0;
	}
	param /= (1*abs_max);
	m_fastsum = new FastSummation(3, N, M, n, m, p, (char*)"one_over_quartic", (char*)"one_over_seventh", &param, eps_I, eps_B);
	m_fastsum->threadedFastSumCorrection(x, y, coeff, coeffx, coeffy, coeffz, sum, sumx, sumy, sumz, c_sum, c_sumx, c_sumy, c_sumz);
	delete m_fastsum;
	//	for (i = 0; i < numOfAtoms; i++)
	//		m_Atoms[i]->Born_radius = 1.0;
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = (sum[i] - m_Atoms[i]->position[0]*sumx[i] - m_Atoms[i]->position[1]*sumy[i] - m_Atoms[i]->position[2]*sumz[i])/(4*pi*pow(2*abs_max, 4));
		(A4 > 0.0) ? m_Atoms[i]->Born_radius = 1.0/A4 : m_Atoms[i]->Born_radius = 100000.0;
		if(A4 < 0.0)
		{
			printf("negative\n");
		}
	}
	/*
		for (i = 0; i < numOfAtoms; i++)
		{
			c_sum[i] /= pow(2*abs_max, 3);	c_sumx[i] /= pow(2*abs_max, 3);
			c_sumy[i] /= pow(2*abs_max, 3);	c_sumz[i] /= pow(2*abs_max, 3);
		}
	*/
	for(i = 0; i < numOfAtoms; i++)
	{
		A4 = 1.0 / m_Atoms[i]->Born_radius;
		A7 = (c_sum[i]-m_Atoms[i]->position[0]*c_sumx[i]-m_Atoms[i]->position[1]*c_sumy[i]-m_Atoms[i]->position[2]*c_sumz[i])/(16*pi*pow(2*abs_max, 7));
		//		printf("%f %f\n", A4, sqrt(sqrt(A7)));
		if(A7 < 0.0)
		{
			A7 = 0.0;
		}
		//                A7 = sqrt(A7); //A7 = sqrt(A7); // wenqi added on 12/17
		A7 = sqrt(A7);
		A7 = sqrt(A7); // wenqi modified on 10/18/08
		if(A4 + A7 > 0.0000001)
		{
			m_Atoms[i]->Born_radius = 1.0 / ((1.0-1.0/sqrt(2.0)) * A4 + A7);    // this is GBSW model
		}
		else
		{
			m_Atoms[i]->Born_radius = 100000.0;
		}
		if(m_Atoms[i]->Born_radius < 0.0)
		{
			printf("negative %d %f\n", i, m_Atoms[i]->Born_radius);
		}
	}
#endif // _FASTSUM_ && _CORRECTION_
	/*
		double testsum[20];
		double num, den, xy[3];
		int a;
		for (a = 0; a < 20; a++)
		{
			testsum[a] = 0.0;
			for (i = 0; i < surface->numbtris; i++)
			{
				for (k = 0; k < nn; k++)
				{
					xy[0] = surface->m_Faces[i]->node[k].Point[0] - m_Atoms[a]->position[0];
					xy[1] = surface->m_Faces[i]->node[k].Point[1] - m_Atoms[a]->position[1];
					xy[2] = surface->m_Faces[i]->node[k].Point[2] - m_Atoms[a]->position[2];
					den = dotProduct(xy, xy);
					if (sqrt(den) > 4.0)	continue;
					den = den*den;
					num = dotProduct(xy, surface->m_Faces[i]->node[k].Normal);
					num *= surface->m_Faces[i]->weight[k];
					testsum[a] += num / den;
				}
			}
			printf("%f, %f\n", sum[a]-m_Atoms[a]->position[0]*sumx[a]-m_Atoms[a]->position[1]*sumy[a]-m_Atoms[a]->position[2]*sumz[a], testsum[a]);
		}
	*/
	/*
		// correction
	        p = 4;	// Wenqi added on 12/17
	        m_fastsum = new FastSummation(3, N, M, n, m, p, "one_over_seventh", &param, eps_I, eps_B);
	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
				coeff[i*nn+k] *= 1.0 / (64*abs_max*abs_max*abs_max);
				coeffx[i*nn+k] *= 1.0 / (64*abs_max*abs_max*abs_max);
				coeffy[i*nn+k] *= 1.0 / (64*abs_max*abs_max*abs_max);
				coeffz[i*nn+k] *= 1.0 / (64*abs_max*abs_max*abs_max);
	                }
	        }

		m_fastsum->threadedFastSum(x, y, coeff, coeffx, coeffy, coeffz, sum, sumx, sumy, sumz);

	        for (i = 0; i < numOfAtoms; i++)
	        {
	                A7 = (sum[i]-m_Atoms[i]->position[0]*sumx[i]-m_Atoms[i]->position[1]*sumy[i]-m_Atoms[i]->position[2]*sumz[i])/(16*pi);
			if(A7 < 0.0) A7 = 0.0;
			A7 = sqrt(A7); A7 = sqrt(A7); // wenqi added on 12/17
	        }
		delete m_fastsum;

	        for (i = 0; i < numOfAtoms; i++)
	        {
			A4 = m_Atoms[i]->Born_radius;
			if( A4 + A7 > 0.0000001)
				m_Atoms[i]->Born_radius = 1.0 / ( (1.0-1.0/sqrt(2.0)) * A4 + A7 ); // this is GBSW model
			else
				m_Atoms[i]->Born_radius = 100000.0;
	                if (m_Atoms[i]->Born_radius < 0.0)
	                        printf("%d %f\n", i, m_Atoms[i]->Born_radius);
	        }
	*/
	/*
	///////////////////// Delete Below /////////////////////////////////////////////////////////////////////////
		// x-component of the second summation
	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
	                        coeff[i*nn+k] = surface->m_Faces[i]->node[k].Normal[0];
	                        coeff[i*nn+k] *= surface->m_Faces[i]->weight[k] / pow(4*abs_max, 4);
	                }
	        }
	        m_fastsum->fastSum(x, y, coeff, sum2);

		for (i = 0; i < numOfAtoms; i++)
		{
			sum1[i] -= m_Atoms[i]->position[0]*sum2[i];
		}

	        // y-component of the second summation
	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
	                        coeff[i*nn+k] = surface->m_Faces[i]->node[k].Normal[1];
	                        coeff[i*nn+k] *= surface->m_Faces[i]->weight[k] / pow(4*abs_max, 4);
	                }
	        }
	        m_fastsum->fastSum(x, y, coeff, sum2);

	        for (i = 0; i < numOfAtoms; i++)
	        {
	                sum1[i] -= m_Atoms[i]->position[1]*sum2[i];
	        }

		// z-component of the second summation
	        for (i = 0; i < surface->numbtris; i++)
	        {
	                for (k = 0; k < nn; k++)
	                {
	                        coeff[i*nn+k] = surface->m_Faces[i]->node[k].Normal[2];
	                        coeff[i*nn+k] *= surface->m_Faces[i]->weight[k] / pow(4*abs_max, 4);
	                }
	        }
	        m_fastsum->fastSum(x, y, coeff, sum2);

	        for (i = 0; i < numOfAtoms; i++)
	        {
	                sum1[i] -= m_Atoms[i]->position[2]*sum2[i];
	        }

		for (i = 0; i < numOfAtoms; i++)
		{
			m_Atoms[i]->Born_radius = 4*pi/sum1[i];
			if (m_Atoms[i]->Born_radius < 0.0)
				printf("%d %f\n", i, m_Atoms[i]->Born_radius);
		}
	*/
	/*
	        for (i = 0; i < numOfAtoms; i++)
	        {
	                printf("%d %f\n", i, sum1[i]);
	                m_Atoms[i]->Born_radius = 0.0;
	        }
	*/
	/* direct sum method
	    int i;
	    for (i = 0; i < numOfAtoms; i++)
	    {
	        m_Atoms[i]->getBornRadius(surface);
	//        printf("atom %d: %f\n", i+1, m_Atoms[i]->Born_radius);
	    }
	*/
}


double GroupOfAtoms::dotProduct(double* u, double* v)
{
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}


void GroupOfAtoms::getNeighbors()
{
	int i, j;
	neighbor = new bool[numOfAtoms*numOfAtoms];
	for(i = 0; i < numOfAtoms*numOfAtoms; i++)
	{
		neighbor[i] = 0;
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		for(j = i; j < numOfAtoms; j++)
		{
			if(j == i)
			{
				neighbor[i* numOfAtoms+j] = 1;
				continue;
			}
			if(overlap(i, j))
			{
				neighbor[i* numOfAtoms+j] = 1;
			}
		}
	}
	for(i = 0; i < numOfAtoms; i++)
	{
		for(j = 0; j < i; j++)
		{
			neighbor[i* numOfAtoms+j] = neighbor[j*numOfAtoms+i];
		}
	}
}

bool GroupOfAtoms::overlap(int i, int j)
{
	double xi[3], xj[3];
	double ri, rj;
	double dij;
	const double w = 1.4;
	int k;
	for(k = 0; k < 3; k++)
	{
		xi[k] = m_Atoms[i]->position[k];
		xj[k] = m_Atoms[j]->position[k];
	}
	ri = m_Atoms[i]->radius+w;
	rj = m_Atoms[j]->radius+w;
	dij = (xi[0]-xj[0])*(xi[0]-xj[0]) + (xi[1]-xj[1])*(xi[1]-xj[1]) + (xi[2]-xj[2])*(xi[2]-xj[2]);
	if(dij < (ri+rj)*(ri+rj))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void GroupOfAtoms::assignBRtoChild(GroupOfAtoms* mol)
{
	int i, j;
	for(i = 0; i < numOfAtoms; i++)
	{
		if(m_Atoms[i]->child !=0)
		{
			double BR = m_Atoms[i]->Born_radius;
			for(j = 0; j < m_Atoms[i]->numOfChildren; j++)
			{
				int a = m_Atoms[i]->child[j];
				mol->m_Atoms[a]->Born_radius = BR;
			}
		}
	}
}

