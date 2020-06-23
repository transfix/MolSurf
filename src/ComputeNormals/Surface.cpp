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
#include <ComputeNormals/Surface.h>

Surface::Surface()
{
}

void Surface::addPoint(Vertex* point)
{
	Vertex* p = new Vertex();
	p->Point[0] = point->Point[0];
	p->Point[1] = point->Point[1];
	p->Point[2] = point->Point[2];
	if(normalFlag == true)
	{
		p->Normal[0] = point->Normal[0];
		p->Normal[1] = point->Normal[1];
		p->Normal[2] = point->Normal[2];
	}
	if(colorFlag == true)
	{
		p->Color[0] = point->Color[0];
		p->Color[1] = point->Color[1];
		p->Color[2] = point->Color[2];
	}
	p->weight = point->weight;
	m_Vertices.push_back(p);
}

void Surface::addFacet(Face* facet)
{
	Face* f = new Face();
	f->Index[0] = facet->Index[0];
	f->Index[1] = facet->Index[1];
	f->Index[2] = facet->Index[2];
	f->IndexInRAW[0] = facet->IndexInRAW[0];
	f->IndexInRAW[1] = facet->IndexInRAW[1];
	f->IndexInRAW[2] = facet->IndexInRAW[2];
	f->AdjTri[0] = facet->AdjTri[0];
	f->AdjTri[1] = facet->AdjTri[1];
	f->AdjTri[2] = facet->AdjTri[2];
	f->whichnml[0] = facet->whichnml[0];
	f->whichnml[1] = facet->whichnml[1];
	f->whichnml[2] = facet->whichnml[2];
	m_Faces.push_back(f);
}

void Surface::BoundingBox()
{
	int i, j;
	min[0] = 100000.0;
	min[1] = 100000.0;
	min[2] = 100000.0;
	max[0] = -100000.0;
	max[1] = -100000.0;
	max[2] = -100000.0;
	for(i = 0; i < numbtris; i++)
	{
		double p[3];
		for(j = 0; j < 3; j++)
		{
			p[j] = m_Vertices[m_Faces[i]->Index[0]]->Point[j];
			if(p[j] < min[j])
			{
				min[j] = p[j];
			}
			if(p[j] > max[j])
			{
				max[j] = p[j];
			}
		}
	}
}

double Surface::DotProduct(double* v1, double* v2)
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void Surface::CrossProduct(double* v1, double* v2, double* v)
{
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

// Compute adjacent information
void Surface::Valence_Neighbor()
{
	int ii, jj, kk, i, j, k, l, m, cur, iface;
	// set to zero
	for(i = 0; i < numbpts; i++)
	{
		m_Vertices[i]->facevalence = 0;
	}
	// compute the facevalence
	for(i = 0; i < numbtris; i++)
	{
		ii = m_Faces[i]->IndexInRAW[0];
		jj = m_Faces[i]->IndexInRAW[1];
		kk = m_Faces[i]->IndexInRAW[2];
		m_Vertices[ii]->facevalence += 1;
		m_Vertices[jj]->facevalence += 1;
		m_Vertices[kk]->facevalence += 1;
	}
	// compute neighbor tri index
	for(i = 0; i < numbtris; i++)
	{
		for(j = 0; j < 3; j++)
		{
			k = m_Faces[i]->IndexInRAW[j];
			m_Vertices[k]->faceneighbor.push_back(i);
		}
	}
	// compute neighbor vertex index
	for(i = 0; i < numbpts; i++)
	{
		for(j = 0; j < m_Vertices[i]->facevalence; j++)
		{
			iface = m_Vertices[i]->faceneighbor[j];
			for(l = 0; l < 3; l++)
			{
				ii = m_Faces[iface]->IndexInRAW[l];
				if(ii != i)
				{
					for(m = 0; m < m_Vertices[i]->edgeneighbor.size(); m++)
					{
						if(ii == m_Vertices[i]->edgeneighbor[m])
						{
							break;
						}
					}
					if(m >= m_Vertices[i]->edgeneighbor.size())
					{
						m_Vertices[i]->edgeneighbor.push_back(ii);
					}
				}
			}
		}
		m_Vertices[i]->edgevalence = m_Vertices[i]->edgeneighbor.size();
	}
}

void Surface::Produce_Normal(bool item_no)
{
	ProduceFaceNormal();
	Valence_Neighbor();
	if(item_no)
	{
		NormalByLoop();
	}
	else
	{
		NormalByAveraging();
	}

}

// TriangleNormal--compute the normal of a triangle, by right-handed rule
void Surface::TriangleNormal(double* p1, double* p2, double* p3, double* normal)
{
	double       x1,y1,z1,x2,y2,z2;
	x1 = p1[0] - p3[0];
	y1 = p1[1] - p3[1];
	z1 = p1[2] - p3[2];
	x2 = p2[0] - p3[0];
	y2 = p2[1] - p3[1];
	z2 = p2[2] - p3[2];
	normal[0] = y1*z2 - y2*z1;
	normal[1] = x2*z1 - x1*z2;
	normal[2] = x1*y2 - x2*y1;
	x1 = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
	x1 = (double)sqrt(x1);
	if(x1 != 0.0)
	{
		normal[0] = normal[0]/x1;
		normal[1] = normal[1]/x1;
		normal[2] = normal[2]/x1;
	}
	else
	{
		printf("normal = %f , %f, %f\n", normal[0], normal[1], normal[2]);
	}
}

void Surface::ProduceFaceNormal()
{
	int    i,j,ii,jj,kk;
	double  p1[3], p2[3], p3[3], normal[3], center[3];
	for(i = 0; i < numbtris; i++)
	{
		ii = m_Faces[i]->IndexInRAW[0];
		jj = m_Faces[i]->IndexInRAW[1];
		kk = m_Faces[i]->IndexInRAW[2];
		for(j = 0; j < 3; j++)
		{
			p1[j] = m_Vertices[ii]->Point[j];
			p2[j] = m_Vertices[jj]->Point[j];
			p3[j] = m_Vertices[kk]->Point[j];
			center[j] = (double)((p1[j] + p2[j] + p3[j])/3.0);
			m_Faces[i]->Center[j] = center[j];
		}
		TriangleNormal(p1,p2,p3,normal);
		for(j = 0; j < 3; j++)
		{
			m_Faces[i]->Normal[j] = normal[j];
		}
	}
}

// LoopNormal -- Compute normal by loop's subdivision
void  Surface::LoopNormal(double* points, int n, double* normal)
{
	double u[3],v[3],pin,c,d;
	int   i, j, k;
	pin = (double)(3.1415926/n);
	pin = pin + pin;
	for(i = 0; i < 3; i++)
	{
		u[i] = 0.0;
		v[i] = 0.0;
	}
	for(i = 1; i <= n; i++)
	{
		j = i + 1;
		if(j > n)
		{
			j = j - n;
		}
		c = (double)cos(i*pin);
		d = (double)cos(j*pin);
		for(k = 0; k < 3; k++)
		{
			u[k] = u[k] + c*points[(i-1)*3+k];
			v[k] = v[k] + d*points[(i-1)*3+k];
		}
	}
	CrossProduct(u,v,normal);
	c = DotProduct(normal,normal);
	c = (double)sqrt(c);
	for(k = 0; k < 3; k++)
	{
		normal[k] = normal[k]/c;
	}
}

//NormalByLoop - Compute vertex normal by loop's method
void Surface::NormalByLoop()
{
	double points[90], normal[3], sum;
	int   l, i, j, k, ii, jj, kk, iii, jjj, kkk, i0, j0,
	  currentvalence, tou, wei, newwei, curnum, w;
	// first by averging
	NormalByAveraging();
	for(l = 0; l < numbpts; l++)
	{
		currentvalence = m_Vertices[l]->facevalence;
		if((currentvalence > 30) ||
				(currentvalence < 2) ||
				(m_Vertices[l]->facevalence != m_Vertices[l]->edgevalence))   // begin if--1
		{
			printf("facevalence[%d]  = %d, edgevalence[%d] = %d\n",
				   l, m_Vertices[l]->facevalence, l, m_Vertices[l]->edgevalence);
		}
		else
		{
			// initialization
			i0 = m_Vertices[l]->faceneighbor[0]; // current triangle
			ii = m_Faces[i0]->IndexInRAW[0];
			jj = m_Faces[i0]->IndexInRAW[1];
			kk = m_Faces[i0]->IndexInRAW[2];
			if(ii == l)
			{
				for(k = 0; k < 3; k++)
				{
					points[k] = m_Vertices[jj]->Point[k];
					points[3+k] = m_Vertices[kk]->Point[k];
				}
				tou = jj;
				wei = kk;
			}
			if(jj == l)
			{
				for(k = 0; k < 3; k++)
				{
					points[k] = m_Vertices[ii]->Point[k];
					points[3+k] = m_Vertices[kk]->Point[k];
				}
				tou = ii;
				wei = kk;
			}
			if(kk == l)
			{
				for(k = 0; k < 3; k++)
				{
					points[k] = m_Vertices[ii]->Point[k];
					points[3+k] = m_Vertices[jj]->Point[k];
				}
				tou = ii;
				wei = jj;
			}
			curnum = 2;
			i = 0;
			while(wei != tou && curnum != currentvalence && i < currentvalence)
			{
				i = i + 1;
				for(j = 0; j < currentvalence; j++)
				{
					j0 = m_Vertices[l]->faceneighbor[j];
					if(i0 != j0)
					{
						iii = m_Faces[j0]->IndexInRAW[0];
						jjj = m_Faces[j0]->IndexInRAW[1];
						kkk = m_Faces[j0]->IndexInRAW[2];
						w = -1;
						if(wei == iii)
						{
							w = iii;
							if(jjj == l)
							{
								newwei = kkk;
							}
							if(kkk == l)
							{
								newwei = jjj;
							}
						}
						if(wei == jjj)
						{
							w = jjj;
							if(iii == l)
							{
								newwei = kkk;
							}
							if(kkk == l)
							{
								newwei = iii;
							}
						}
						if(wei == kkk)
						{
							w = kkk;
							if(iii == l)
							{
								newwei = jjj;
							}
							if(jjj == l)
							{
								newwei = iii;
							}
						}
						if(w >= 0)
						{
							wei = newwei;
							if(newwei != tou)
							{
								i0 = j0;
								for(k = 0; k < 3; k++)
								{
									points[curnum*3+k] = m_Vertices[wei]->Point[k];
								}
								curnum = curnum + 1;
							}
						}
					}
				}
			}
			LoopNormal(points,curnum,normal);
			sum = 1.0;
			if(DotProduct(m_Vertices[l]->Normal, normal) < 0.0)
			{
				sum = -1.0;
			}
			for(k = 0; k < 3; k++)
			{
				m_Vertices[l]->Normal[k] = sum*normal[k];
			}
		}
	}
}

// NormalByAveraging-- Compute vertex normal by averaging face normals
void Surface::NormalByAveraging()
{
	int	i, j, k, jj;
	float	sum, normal[3];
	/* compute normal at vertices */
	/* by averging the face normal*/
	for(i = 0; i < numbpts; i++)
	{
		if(m_Vertices[i]->facevalence > 0)
		{
			normal[0] = 0.0;
			normal[1] = 0.0;
			normal[2] = 0.0;
			for(j = 0; j < m_Vertices[i]->facevalence; j++)
			{
				jj = m_Vertices[i]->faceneighbor[j];
				for(k = 0; k < 3; k++)
				{
					normal[k] = normal[k] + m_Faces[jj]->Normal[k];
				}
			}
			sum = 0.0;
			for(k = 0; k < 3; k++)
			{
				sum = sum + normal[k]*normal[k];
			}
			sum = (float)sqrt(sum);
			for(k = 0; k < 3; k++)
			{
				m_Vertices[i]->Normal[k] = normal[k]/sum;
			}
		}
		else
		{
			printf("The point %d is a isolated point\n", i);
		}
	}
}
