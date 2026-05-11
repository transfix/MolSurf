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
#include <ASpline/ASpline.h>

using namespace CVCUtility;
using namespace ASPLINE;

int ReadHowmanyComponents(FILE *fp) {
  char str[256];
  fgetsSafely(str, 256, fp);
  fgetsSafely(str, 256, fp);
  double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11;
  int component =
      sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &f0, &f1,
             &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11);
  return (component);
}

ASpline::ASpline(FILE *fp) // from .raw
{
  // duplicates code in LBIE_lib/geoframe::read_raw()
  int howmany = ReadHowmanyComponents(fp);
  if (howmany != 3 && howmany != 4 && howmany != 6 && howmany != 9) {
    error("Wrong number of components in the geometry file!\n");
  }
  if (howmany == 3 | howmany == 4) {
    normalFlag = false; /* have function value on vertex */
    colorFlag = false;
  } else if (howmany == 6) {
    normalFlag = true;
    colorFlag = false;
  } else if (howmany == 9) {
    normalFlag = true; /* have rgb value on vertex */
    colorFlag = true;
  }
  // Bounding box dimensions
  max[0] = max[1] = max[2] = -boundingBoxSize;
  min[0] = min[1] = min[2] = boundingBoxSize;
  int minIndex = 1000;
  // Get the number of points and triangles from the file
  rewind(fp);
  if (fscanf(fp, "%d %d\n", &numbpts, &numbtris) != 2)
    error("Could not read the number of vertices and triangles!");
  // Iterate over points, adding Vertex*s to our Aspline
  for (int i = 0; i < numbpts; i++) {
    // Read in data
    double x, y, z, f, nx, ny, nz, red, green, blue = 0;
    if (howmany == 3) {
      int nscan = fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);
      if (nscan != 3) {
        error("Wrong data type in geometry! Expected 3 data points.\n");
      }
    }
    if (howmany == 4) {
      int nscan = fscanf(fp, "%lf %lf %lf %lf\n", &x, &y, &z, &f);
      if (nscan != 4) {
        error("Wrong data type in geometry! Expected 4 data points.\n");
      }
    }
    if (howmany == 6) {
      int nscan =
          fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &nx, &ny, &nz);
      if (nscan != 6) {
        error("Wrong data type in geometry! Expected 6 data points.\n");
      }
    }
    if (howmany == 9) {
      int nscan = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x, &y,
                         &z, &nx, &ny, &nz, &red, &green, &blue);
      if (nscan != 9) {
        error("Wrong data type in geometry! Expected 9 data points.\n");
      }
    }
    // Refine bounding box
    if (x < min[0]) {
      min[0] = x;
    }
    if (y < min[1]) {
      min[1] = y;
    }
    if (z < min[2]) {
      min[2] = z;
    }
    if (x > max[0]) {
      max[0] = x;
    }
    if (y > max[1]) {
      max[1] = y;
    }
    if (z > max[2]) {
      max[2] = z;
    }

    // arand: normalize the normals
    double len = sqrt(nx * nx + ny * ny + nz * nz);
    if (len > 0.0001) {
      nx /= len;
      ny /= len;
      nz /= len;
    }

    // Add our vertex
    addPoint(new Vertex(x, y, z, nx, ny, nz, red, green, blue));
  }
  // Iterate over triangles, adding Face*s to our Aspline
  for (int i = 0; i < numbtris; i++) {
    int index0, index1, index2;
    if (fscanf(fp, "%d %d %d\n", &index0, &index1, &index2) != 3)
      error("Could not read triangle!");
    Face *facet = new Face(index0, index1, index2);
    addFacet(facet);
    // index[0] is the least index of the three
    if (facet->index[0] < minIndex) {
      minIndex = facet->index[0];
    }
  }
  // ?
  for (int i = 0; i < numbtris; i++) {
    faces[i]->index[0] = faces[i]->index[0] - minIndex;
    faces[i]->index[1] = faces[i]->index[1] - minIndex;
    faces[i]->index[2] = faces[i]->index[2] - minIndex;
  }
  cout << "minIndex = " << minIndex << " numbtris = " << numbtris
       << " numbpts = " << numbpts << "\n";
}

// Coerce back into a Geometry
// In ASpline since Geoframe has more public members
// There are four classes which store the same data in different ways. We only
// need one. The four classes are Surface, Geometry, Geoframe, and ASpline
Geometry *ASpline::getGeometry() {

  //
  // WARNING: this only works AFTER constructing the a-spline surface
  //          i.e. after the subfacets are created...
  //
  //

  Geometry *geometry = new Geometry();

  float r = 0;
  float g = 0;
  float b = 0;

  // geometry->AllocatePoints(verticies.size()); // unneded, line below handles
  // this...
  geometry->AllocateTris(verticies.size(),
                         faces.size() * faces[0]->subFacets.size());
  geometry->AllocateTriVertColors();

  for (int i = 0; i < verticies.size(); i++) {
    // geometry->m_Points[i*3+0] = verticies[i]->point[0];
    // geometry->m_Points[i*3+1] = verticies[i]->point[1];
    // geometry->m_Points[i*3+2] = verticies[i]->point[2];

    geometry->m_TriVerts[i * 3 + 0] = verticies[i]->point[0];
    geometry->m_TriVerts[i * 3 + 1] = verticies[i]->point[1];
    geometry->m_TriVerts[i * 3 + 2] = verticies[i]->point[2];
    geometry->m_TriVertNormals[i * 3 + 0] = verticies[i]->normal[0];
    geometry->m_TriVertNormals[i * 3 + 1] = verticies[i]->normal[1];
    geometry->m_TriVertNormals[i * 3 + 2] = verticies[i]->normal[2];
    // Ignoring alpha channel
    geometry->m_TriVertColorsTransparent[i * 4 + 0] = verticies[i]->color[0];
    geometry->m_TriVertColorsTransparent[i * 4 + 1] = verticies[i]->color[1];
    geometry->m_TriVertColorsTransparent[i * 4 + 2] = verticies[i]->color[2];

    /* // eventually duplicate this code
       // which is more robust to different case (including vertex normals or
    colors... if(aspline->normalFlag == true && aspline->colorFlag == true)
    {
            fprintf(fp, "%f %f %f %f %f %f %f %f %f\n",
    P[0],P[1],P[2],N[0],N[1],N[2],C[0],C[1],C[2]);
    }
    else if(aspline->normalFlag == true && aspline->colorFlag == false)
    {
            fprintf(fp, "%f %f %f %f %f %f\n", P[0],P[1],P[2],N[0],N[1],N[2]);
    }
    else if(aspline->normalFlag == false && aspline->colorFlag == true)
    {
            fprintf(fp, "%f %f %f %f %f %f\n", P[0],P[1],P[2],C[0],C[1],C[2]);
    }
    else if(aspline->normalFlag == false && aspline->colorFlag == false)
    {
            fprintf(fp, "%f %f %f\n", P[0],P[1],P[2]);
    }
    */
  }

  int count = 0;
  for (int i = 0; i < faces.size(); i++) {
    for (int t = 0; t < faces[i]->subFacets.size(); t++) {
      int index[3];
      for (int j = 0; j < 3; j++) {
        index[j] = faces[i]->subFacets[t]->indexInRAW[j];
      }
      geometry->m_Tris[count * 3 + 0] = index[0];
      geometry->m_Tris[count * 3 + 1] = index[1];
      geometry->m_Tris[count * 3 + 2] = index[2];
      count++;
    }
  }

  /* // old incorrect version...
  for(int i=0; i<numbtris; i++)
  {
          geometry->m_Tris[i*3+0] = faces[i]->index[0];
          geometry->m_Tris[i*3+1] = faces[i]->index[1];
          geometry->m_Tris[i*3+2] = faces[i]->index[2];
  }
  */

  return geometry;
}

// Coerce back into a Geoframe
// In ASpline since Geoframe has more public members
Geoframe *ASpline::getGeoframe() {
  Geoframe *geoframe = new Geoframe();
  geoframe->meshType = TRI;

  float r = 0;
  float g = 0;
  float b = 0;

  geoframe->numVerts = numbpts;
  geoframe->numTris = numbtris;
  geoframe->verts = (float(*)[3])malloc(sizeof(float[3]) * geoframe->numVerts);
  geoframe->normals = (float(*)[3])malloc(
      sizeof(float[3]) * minimum(geoframe->numVerts, geoframe->numTris));
  geoframe->triangles = (uint(*)[3])malloc(sizeof(uint[3]) * geoframe->numTris);
  geoframe->bound_sign =
      (uint(*))realloc(geoframe->bound_sign, sizeof(uint) * geoframe->vsize);
  geoframe->bound_tri =
      (uint(*))realloc(geoframe->bound_tri, sizeof(uint) * geoframe->numTris);

  for (int i = 0; i < geoframe->numVerts; i++) {
    geoframe->verts[i][0] = verticies[i]->point[0];
    geoframe->verts[i][1] = verticies[i]->point[1];
    geoframe->verts[i][2] = verticies[i]->point[2];
    geoframe->normals[i][0] = verticies[i]->normal[0];
    geoframe->normals[i][1] = verticies[i]->normal[1];
    geoframe->normals[i][2] = verticies[i]->normal[2];
    geoframe->colors[i][0] = verticies[i]->color[0];
    geoframe->colors[i][1] = verticies[i]->color[1];
    geoframe->colors[i][2] = verticies[i]->color[2];
  }
  for (int i = 0; i < geoframe->numTris; i++) {
    geoframe->bound_tri[i] = 1;
    geoframe->triangles[i][0] = faces[i]->index[0];
    geoframe->triangles[i][1] = faces[i]->index[1];
    geoframe->triangles[i][2] = faces[i]->index[2];
  }

  return geoframe;
}

// From a .raw that has already been opened
ASpline::ASpline(Geoframe *geoframe) {
  // If there is no color in the file, Geoframe creates it
  colorFlag = true;
  // Normals only if the geoframe was of a type that held them.
  normalFlag =
      (geoframe->getFileType() == RAWN || geoframe->getFileType() == RAWNC);
  // Bounding box dimensions
  max[0] = max[1] = max[2] = -boundingBoxSize;
  min[0] = min[1] = min[2] = boundingBoxSize;
  int minIndex = 1000;
  numbpts = geoframe->getNumVerts();
  numbtris = geoframe->getNumTris();
  // Iterate over points, adding Vertex*s to our Aspline
  for (int i = 0; i < numbpts; i++) {
    // Read in data
    double x, y, z, f, nx, ny, nz, red, green, blue = 0;
    x = geoframe->verts[i][0];
    y = geoframe->verts[i][1];
    z = geoframe->verts[i][2];
    red = geoframe->colors[i][0];
    green = geoframe->colors[i][1];
    blue = geoframe->colors[i][2];
    if (normalFlag) {
      nx = geoframe->normals[i][0];
      ny = geoframe->normals[i][1];
      nz = geoframe->normals[i][2];
    }

    // Refine bounding box
    if (x < min[0]) {
      min[0] = x;
    }
    if (y < min[1]) {
      min[1] = y;
    }
    if (z < min[2]) {
      min[2] = z;
    }
    if (x > max[0]) {
      max[0] = x;
    }
    if (y > max[1]) {
      max[1] = y;
    }
    if (z > max[2]) {
      max[2] = z;
    }
    // Add our vertex
    addPoint(new Vertex(x, y, z, nx, ny, nz, red, green, blue));
  }
  // Iterate over triangles, adding Face*s to our Aspline
  for (int i = 0; i < numbtris; i++) {
    int index0, index1, index2;
    index0 = geoframe->triangles[i][0];
    index1 = geoframe->triangles[i][1];
    index2 = geoframe->triangles[i][2];
    Face *facet = new Face(index0, index1, index2);
    addFacet(facet);
    // index[0] is the least index of the three
    if (facet->index[0] < minIndex) {
      minIndex = facet->index[0];
    }
  }
  // ?
  for (int i = 0; i < numbtris; i++) {
    faces[i]->index[0] = faces[i]->index[0] - minIndex;
    faces[i]->index[1] = faces[i]->index[1] - minIndex;
    faces[i]->index[2] = faces[i]->index[2] - minIndex;
  }
  cout << "minIndex = " << minIndex << " numbtris = " << numbtris
       << " numbpts = " << numbpts << "\n";
}

// From a .raw that has already been opened
ASpline::ASpline(Geometry *geometry) {
  // In a Geometry, if an array is set to zero, there is none of that type of
  // data
  colorFlag = geometry->m_TriVertColorsTransparent;
  normalFlag = geometry->m_TriVertNormals;
  // Bounding box dimensions
  max[0] = max[1] = max[2] = -boundingBoxSize;
  min[0] = min[1] = min[2] = boundingBoxSize;
  int minIndex = 1000;
  // numbpts = geometry->m_NumPoints;
  numbpts = geometry->m_NumTriVerts; // arand change to prevent crash
  numbtris = geometry->m_NumTris;
  // Iterate over points, adding Vertex*s to our Aspline
  for (int i = 0; i < numbpts; i++) {
    // Read in data
    double x, y, z, f, nx, ny, nz, red, green, blue = 0;
    // x   = geometry->m_Points[i*3+0];
    // y   = geometry->m_Points[i*3+1];
    // z   = geometry->m_Points[i*3+2];
    x = geometry->m_TriVerts[i * 3 + 0];
    y = geometry->m_TriVerts[i * 3 + 1];
    z = geometry->m_TriVerts[i * 3 + 2];

    if (geometry->m_UseTriangleColors) {
      red = geometry->m_TriVertColorsTransparent[i * 4 + 0];
      green = geometry->m_TriVertColorsTransparent[i * 4 + 1];
      blue = geometry->m_TriVertColorsTransparent[i * 4 + 2];
    }
    if (normalFlag) {
      nx = geometry->m_TriVertNormals[i * 3 + 0];
      ny = geometry->m_TriVertNormals[i * 3 + 1];
      nz = geometry->m_TriVertNormals[i * 3 + 2];

      // arand: normalize the normals...
      double len = sqrt(nx * nx + ny * ny + nz * nz);
      if (len > 0.00001) {
        nx /= len;
        ny /= len;
        nz /= len;
      }
    }

    // Refine bounding box
    if (x < min[0]) {
      min[0] = x;
    }
    if (y < min[1]) {
      min[1] = y;
    }
    if (z < min[2]) {
      min[2] = z;
    }
    if (x > max[0]) {
      max[0] = x;
    }
    if (y > max[1]) {
      max[1] = y;
    }
    if (z > max[2]) {
      max[2] = z;
    }
    // Add our vertex
    addPoint(new Vertex(x, y, z, nx, ny, nz, red, green, blue));
  }
  // Iterate over triangles, adding Face*s to our Aspline
  for (int i = 0; i < numbtris; i++) {
    int index0, index1, index2;
    index0 = geometry->m_Tris[i * 3 + 0];
    index1 = geometry->m_Tris[i * 3 + 1];
    index2 = geometry->m_Tris[i * 3 + 2];
    Face *facet = new Face(index0, index1, index2);
    addFacet(facet);
    // index[0] is the least index of the three
    if (facet->index[0] < minIndex) {
      minIndex = facet->index[0];
    }
  }
  // ?
  for (int i = 0; i < numbtris; i++) {
    faces[i]->index[0] = faces[i]->index[0] - minIndex;
    faces[i]->index[1] = faces[i]->index[1] - minIndex;
    faces[i]->index[2] = faces[i]->index[2] - minIndex;
  }
  // cout << "minIndex = " << minIndex << " numbtris = " << numbtris << "
  // numbpts = " << numbpts << "\n";
}

ASpline::ASpline() {}

void ASpline::addPoint(Vertex *v_) {
  Vertex *v = new Vertex();
  v = v_;
  verticies.push_back(v);
}

void ASpline::addFacet(Face *f_) {
  Face *f = new Face();
  f = f_;
  faces.push_back(f);
}

void ASpline::BoundingBox() {
  min[0] = 100000.0;
  min[1] = 100000.0;
  min[2] = 100000.0;
  max[0] = -100000.0;
  max[1] = -100000.0;
  max[2] = -100000.0;
  for (int i = 0; i < numbtris; i++) {
    double p[3];
    for (int j = 0; j < 3; j++) {
      p[j] = verticies[faces[i]->index[0]]->point[j];
      if (p[j] < min[j]) {
        min[j] = p[j];
      }
      if (p[j] > max[j]) {
        max[j] = p[j];
      }
    }
  }
}

// Craig: moved the following 3 functions here from /MolSurf/rawn.cpp

void ASpline::saveNewRawncFile(FILE *fp, int resolution) {
  fprintf(fp, "%lu %lu\n", verticies.size(),
          faces.size() * resolution * resolution);
  // vertices loop
  for (int i = 0; i < verticies.size(); i++) {
    double P[3], N[3], C[3];
    for (int j = 0; j < 3; j++) {
      P[j] = verticies[i]->point[j];
      N[j] = verticies[i]->normal[j];
    }
    C[0] = verticies[i]->color[0];
    C[1] = verticies[i]->color[1];
    C[2] = verticies[i]->color[2];
    if (normalFlag == true && colorFlag == true) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2],
              N[0], N[1], N[2], C[0], C[1], C[2]);
    } else if (normalFlag == true && colorFlag == false) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], N[0], N[1],
              N[2]);
    } else if (normalFlag == false && colorFlag == true) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], C[0], C[1],
              C[2]);
    } else if (normalFlag == false && colorFlag == false) {
      fprintf(fp, "%lf %lf %lf\n", P[0], P[1], P[2]);
    }
  }
  // face index loop
  for (int i = 0; i < faces.size(); i++) {
    for (int t = 0; t < faces[i]->subFacets.size(); t++) {
      int index[3];
      for (int j = 0; j < 3; j++) {
        index[j] = faces[i]->subFacets[t]->indexInRAW[j];
      }
      fprintf(fp, "%d %d %d\n", index[0], index[1], index[2]);
    }
  }
}

void ASpline::printRawnWithQuadrature(const char *fname, int resolution,
                                      int numQuadPts) {
  FILE *fp = fopen(fname, "w");
  if (!fp) {
    printf("could not open the file for write\n");
    exit(0);
  }
  fprintf(fp, "%lu %lu\n",
          verticies.size() +
              faces.size() * resolution * resolution * numQuadPts,
          faces.size() * resolution * resolution);
  // vertices loop
  for (int i = 0; i < verticies.size(); i++) {
    double P[3], N[3], C[3];
    for (int j = 0; j < 3; j++) {
      P[j] = verticies[i]->point[j];
      N[j] = verticies[i]->normal[j];
    }
    C[0] = verticies[i]->color[0];
    C[1] = verticies[i]->color[1];
    C[2] = verticies[i]->color[2];
    if (normalFlag == true && colorFlag == true) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2],
              N[0], N[1], N[2], C[0], C[1], C[2]);
    } else if (normalFlag == true && colorFlag == false) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], N[0], N[1],
              N[2]);
    } else if (normalFlag == false && colorFlag == true) {
      fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], C[0], C[1],
              C[2]);
    } else if (normalFlag == false && colorFlag == false) {
      fprintf(fp, "%lf %lf %lf\n", P[0], P[1], P[2]);
    }
  }
  // quad points... not part of any triangles
  for (int i = 0; i < faces.size(); i++) {
    for (int j = 0; j < faces[i]->node.size(); j++) {
      double P[3], N[3], C[3];
      for (int k = 0; k < 3; k++) {
        P[k] = faces[i]->node[j]->point[k];
        N[k] = faces[i]->node[j]->normal[k];
      }
      C[0] = 1.0;
      C[1] = 1.0;
      C[2] = 1.0;
      if (normalFlag == true && colorFlag == true) {
        fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2],
                N[0], N[1], N[2], C[0], C[1], C[2]);
      } else if (normalFlag == true && colorFlag == false) {
        fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], N[0], N[1],
                N[2]);
      } else if (normalFlag == false && colorFlag == true) {
        fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], C[0], C[1],
                C[2]);
      } else if (normalFlag == false && colorFlag == false) {
        fprintf(fp, "%lf %lf %lf\n", P[0], P[1], P[2]);
      }
    }
  }
  // face index loop
  for (int i = 0; i < faces.size(); i++) {
    for (int t = 0; t < faces[i]->subFacets.size(); t++) {
      int index[3];
      for (int j = 0; j < 3; j++) {
        index[j] = faces[i]->subFacets[t]->indexInRAW[j];
      }
      fprintf(fp, "%d %d %d\n", index[0], index[1], index[2]);
    }
  }
}

void ASpline::printQuadratureResults(const char *fname) {
  FILE *fp = fopen(fname, "w");
  if (!fp) {
    printf("could not open the file for write\n");
    exit(0);
  }
  fprintf(fp, "%lu\n", faces.size() * faces[0]->node.size());
  for (int i = 0; i < faces.size(); i++) {
    for (int j = 0; j < faces[i]->node.size(); j++) {
      double P[3], N[3];
      for (int k = 0; k < 3; k++) {
        P[k] = faces[i]->node[j]->point[k];
        N[k] = faces[i]->node[j]->normal[k];
      }
      fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", P[0], P[1], P[2], N[0], N[1],
              N[2], faces[i]->node[j]->weight);
    }
  }
}
