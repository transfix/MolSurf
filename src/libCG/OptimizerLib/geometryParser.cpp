#include <libCG/CoarseGrain/geometry.h>
#include <libCG/CoarseGrain/surface.h>
#include <libCG/OptimizerLib/geometryParser.h>
#include <Utility/utility.h>

using namespace CVCUtility;

using namespace GEOMETRY;

GeometryParser::GeometryParser()
{
}

int GeometryParser::ReadHowmanyComponents(FILE* fp)
{
	char  str[256];
	double  f0, f1, f2, f3, f4, f5, f6,f7,f8,f9,f10,f11;
	int    nvert, nface,component;
	fgetsSafely(str, 256, fp);
	sscanf(str, "%d %d", &nvert, &nface);
	fgetsSafely(str, 256, fp);
	component = sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		&f0, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11);
	return(component);
}

Surface* GeometryParser::SurfaceTriangulation(FILE* fp)
{
	double maxx, maxy, maxz, minx, miny, minz;
	Surface* surface = new Surface();
	int           i,ii,jj,kk,i1,i2,i3,minindex,nscan;
	double        x,y,z,f,nx,ny,nz,red,green,blue;
	double         w;
	double         maxf, minf;
	howmany = ReadHowmanyComponents(fp);
	if(howmany != 3 && howmany != 4 && howmany != 6 && howmany !=9)
	{
		printf("howmany = %d\n", howmany);
		printf("Wrong data type in geometry\n");
		exit(0);
	}
	if(howmany == 3)
	{
		normalFlag = 0;
	}
	if(howmany == 4)
	{
		normalFlag = 0;    /* have function value on vertex */
	}
	if(howmany == 6)
	{
		normalFlag = 1;
	}
	if(howmany == 9)
	{
		normalFlag = 1;    /* have rgb value on vertex */
	}
	maxx = -100000.0;
	maxy = -100000.0;
	maxz = -100000.0;
	minx = 100000.0;
	miny = 100000.0;
	minz = 100000.0;
	minindex = 1000;
	rewind(fp);
	if(fscanf(fp,"%d %d\n",&surface->numbpts, &surface->numbtris) != 2) error("Couldn't read the number of triangles and vertices!");
	if(howmany == 4)
	{
	}
	for(i = 0 ; i < surface->numbpts; i++)
	{
		if(howmany == 3)
		{
			nscan = fscanf(fp,"%lf %lf %lf\n", &x, &y, &z);
			if(nscan != 3)
			{
				printf("Wrong data type in geometry\n");
				exit(0);
			}
		}
		if(howmany == 4)
		{
			nscan = fscanf(fp,"%lf %lf %lf %lf\n", &x, &y, &z, &f);
			if(nscan != 4)
			{
				printf("Wrong data type in geometry\n");
				exit(0);
			}
		}
		if(normalFlag == 1)
		{
			if(howmany == 6)
				nscan = fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",
							   &x, &y, &z, &nx, &ny, &nz);
			if(howmany == 9)
				nscan = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
							   &x, &y, &z, &nx, &ny, &nz, &red, &green, &blue);
			if(nscan != 6 && nscan != 9)
			{
				printf("Wrong data type in geometry at %d\n", i);
				exit(0);
			}
		}
		if(x < minx)
		{
			minx = x;
		}
		if(y < miny)
		{
			miny = y;
		}
		if(z < minz)
		{
			minz = z;
		}
		if(x > maxx)
		{
			maxx = x;
		}
		if(y > maxy)
		{
			maxy = y;
		}
		if(z > maxz)
		{
			maxz = z;
		}
		Vertex* point = new Vertex();
		point->Point[0] = x;
		point->Point[1] = y;
		point->Point[2] = z;
		if(normalFlag == 1)
		{
			w = sqrt(nx*nx + ny*ny + nz*nz);
			point->Normal[0] = nx/w;
			point->Normal[1] = ny/w;
			point->Normal[2] = nz/w;
		}
		else
		{
			point->Normal[0] = 0.0;
			point->Normal[1] = 0.0;
			point->Normal[2] = 0.0;
		}
		surface->addPoint(point);
	} /* end point (i) loop */
	for(i = 0 ; i < surface->numbtris; i++)
	{
		if(fscanf(fp,"%d %d %d\n", &ii, &jj, &kk) != 3) error("Error reading traingle!");
		if(ii < jj && jj < kk)
		{
			i1 = ii;
			i2 = jj;
			i3 = kk;
		}
		if(ii < kk && kk < jj)
		{
			i1 = ii;
			i2 = kk;
			i3 = jj;
		}
		if(jj < kk && kk < ii)
		{
			i1 = jj;
			i2 = kk;
			i3 = ii;
		}
		if(jj < ii && ii < kk)
		{
			i1 = jj;
			i2 = ii;
			i3 = kk;
		}
		if(kk < jj && jj < ii)
		{
			i1 = kk;
			i2 = jj;
			i3 = ii;
		}
		if(kk < ii && ii < jj)
		{
			i1 = kk;    // result in i1 <= i2 <= i3
			i2 = ii;
			i3 = jj;
		}
		Face* facet = new Face();
		facet->Index[0] = i1;
		facet->Index[1] = i2;
		facet->Index[2] = i3;
		facet->AdjTri[0] = -1;
		facet->AdjTri[1] = -1;
		facet->AdjTri[2] = -1;
		facet->whichnml[0] = -1;
		facet->whichnml[1] = -1;
		facet->whichnml[2] = -1;
		surface->addFacet(facet);
		if(i1 < minindex)
		{
			minindex = i1;
		}
	}   /* end facet (i) loop */
	for(i = 0 ; i < surface->numbtris; i++)
	{
		surface->m_Faces[i]->Index[0] = surface->m_Faces[i]->Index[0] - minindex;
		surface->m_Faces[i]->Index[1] = surface->m_Faces[i]->Index[1] - minindex;
		surface->m_Faces[i]->Index[2] = surface->m_Faces[i]->Index[2] - minindex;
	}
	//printf("minindex = %d,  numbtris = %d,  numbpts = %d\n", minindex,surface->numbtris,surface->numbpts);
	if(howmany == 4)  /* with function value */
	{
	}
	const double layer = 1.0;
	surface->min[0] = minx-layer;
	surface->min[1] = miny-layer;
	surface->min[2] = minz-layer;
	surface->max[0] = maxx+layer;
	surface->max[1] = maxy+layer;
	surface->max[2] = maxz+layer;
	for(i = 0; i < 3; i++)
	{
		surface->center[i] = (surface->min[i] + surface->max[i])*0.5;
	}
	return surface;
}
