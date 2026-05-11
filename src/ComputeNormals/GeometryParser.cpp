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
#include <ComputeNormals/GeometryParser.h>

#include <iostream>
using namespace std;
using namespace CVCUtility;
using namespace COMPUTENORMALS;

GeometryParser::GeometryParser() {}

int GeometryParser::ReadHowmanyComponents(FILE *fp) {
  char str[256];
  double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;
  int nvert, nface, component;
  // fGetString(str, 256, fp);
  fgetsSafely(str, 256, fp);
  sscanf(str, "%d %d", &nvert, &nface);
  fgetsSafely(str, 256, fp);
  component =
      sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &f0, &f1,
             &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11);
  return (component);
}

Surface *GeometryParser::SurfaceTriangulation(FILE *fp) {
  double maxx, maxy, maxz, minx, miny, minz;
  Surface *surface = new Surface();
  int i, ii, jj, kk, i1, i2, i3, minindex, nscan;
  double x, y, z, f, nx, ny, nz, red, green, blue;
  double w;
  double maxf, minf;
  howmany = ReadHowmanyComponents(fp);
  if (howmany != 3 && howmany != 4 && howmany != 6 && howmany != 9) {
    printf("howmany = %d\n", howmany);
    printf("Wrong data type in geometry\n");
    exit(0);
  }
  if (howmany == 3) {
    normalFlag = 0;
    colorFlag = 0;
  }
  if (howmany == 4) {
    normalFlag = 0; // have function value on vertex
    colorFlag = 0;
  }
  if (howmany == 6) {
    // arand: assume colors not normals...
    // normalFlag = 1;
    //   colorFlag = 0;
    normalFlag = 0;
    colorFlag = 1;
  }
  if (howmany == 9) {
    normalFlag = 1; // have rgb value on vertex
    colorFlag = 1;
  }
  (normalFlag == 0) ? surface->normalFlag = false : surface->normalFlag = true;
  (colorFlag == 0) ? surface->colorFlag = false : surface->colorFlag = true;

  maxx = -100000.0;
  maxy = -100000.0;
  maxz = -100000.0;
  minx = 100000.0;
  miny = 100000.0;
  minz = 100000.0;
  minindex = 1000;
  rewind(fp);
  if (fscanf(fp, "%d %d\n", &surface->numbpts, &surface->numbtris) != 2)
    error("Could not read number of vertices and triangles!");
  for (i = 0; i < surface->numbpts; i++) {
    if (howmany == 3) {
      nscan = fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);
      if (nscan != 3) {
        printf("Wrong data type in geometry\n");
        exit(0);
      }
    }
    if (howmany == 4) {
      nscan = fscanf(fp, "%lf %lf %lf %lf\n", &x, &y, &z, &f);
      if (nscan != 4) {
        printf("Wrong data type in geometry\n");
        exit(0);
      }
    }
    if (normalFlag == 1) {
      if (howmany == 6)
        nscan =
            fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &nx, &ny, &nz);
      if (howmany == 9)
        nscan = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x, &y, &z,
                       &nx, &ny, &nz, &red, &green, &blue);
      if (nscan != 6 && nscan != 9) {
        printf("Wrong data type in geometry at %d\n", i);
        exit(0);
      }
    } else if (colorFlag == 1) {
      if (howmany == 6)
        nscan = fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &red,
                       &green, &blue);
    }
    if (x < minx) {
      minx = x;
    }
    if (y < miny) {
      miny = y;
    }
    if (z < minz) {
      minz = z;
    }
    if (x > maxx) {
      maxx = x;
    }
    if (y > maxy) {
      maxy = y;
    }
    if (z > maxz) {
      maxz = z;
    }
    Vertex *point = new Vertex();
    point->Point[0] = x;
    point->Point[1] = y;
    point->Point[2] = z;
    if (normalFlag == 1) {
      w = 1.0;
      point->Normal[0] = nx / w;
      point->Normal[1] = ny / w;
      point->Normal[2] = nz / w;
    } else {
      point->Normal[0] = 0.0;
      point->Normal[1] = 0.0;
      point->Normal[2] = 0.0;
    }
    if (colorFlag == 1) {
      point->Color[0] = red;
      point->Color[1] = green;
      point->Color[2] = blue;
    }
    point->weight = 1.0; // default node weight
    surface->addPoint(point);
  }
  for (i = 0; i < surface->numbtris; i++) {
    if (fscanf(fp, "%d %d %d\n", &ii, &jj, &kk) != 3)
      error("Could not read triangle!");
    if (ii < jj && jj < kk) {
      i1 = ii;
      i2 = jj;
      i3 = kk;
    }
    if (ii < kk && kk < jj) {
      i1 = ii;
      i2 = kk;
      i3 = jj;
    }
    if (jj < kk && kk < ii) {
      i1 = jj;
      i2 = kk;
      i3 = ii;
    }
    if (jj < ii && ii < kk) {
      i1 = jj;
      i2 = ii;
      i3 = kk;
    }
    if (kk < jj && jj < ii) {
      i1 = kk;
      i2 = jj;
      i3 = ii;
    }
    if (kk < ii && ii < jj) {
      i1 = kk; // result in i1 <= i2 <= i3
      i2 = ii;
      i3 = jj;
    }
    Face *facet = new Face();
    facet->Index[0] = i1;
    facet->Index[1] = i2;
    facet->Index[2] = i3;
    facet->IndexInRAW[0] = ii;
    facet->IndexInRAW[1] = jj;
    facet->IndexInRAW[2] = kk;
    facet->AdjTri[0] = -1;
    facet->AdjTri[1] = -1;
    facet->AdjTri[2] = -1;
    facet->whichnml[0] = -1;
    facet->whichnml[1] = -1;
    facet->whichnml[2] = -1;
    surface->addFacet(facet);
    if (i1 < minindex) {
      minindex = i1;
    }
  }
  for (i = 0; i < surface->numbtris; i++) {
    surface->m_Faces[i]->Index[0] = surface->m_Faces[i]->Index[0] - minindex;
    surface->m_Faces[i]->Index[1] = surface->m_Faces[i]->Index[1] - minindex;
    surface->m_Faces[i]->Index[2] = surface->m_Faces[i]->Index[2] - minindex;
    surface->m_Faces[i]->IndexInRAW[0] =
        surface->m_Faces[i]->IndexInRAW[0] - minindex;
    surface->m_Faces[i]->IndexInRAW[1] =
        surface->m_Faces[i]->IndexInRAW[1] - minindex;
    surface->m_Faces[i]->IndexInRAW[2] =
        surface->m_Faces[i]->IndexInRAW[2] - minindex;
  }
  // printf("minindex = %d,  numbtris = %d,  numbpts = %d\n",
  // minindex,surface->numbtris,surface->numbpts);
  return surface;
}
