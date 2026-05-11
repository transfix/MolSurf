
#include <OctreeSurfaceNew/Octree.h>
#include "MarchingCubes.cpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

using namespace std;

Octree::Octree(std::string pqrFN, std::string rawFN, int mindpth, int maxdpth,
               float blob, float iso) {
  minDepth = mindpth;
  maxDepth = maxdpth;
  blobbyness = blob;
  isoValue = iso;

  if (readPQR(pqrFN)) {
    cout << "PQR read" << endl;

    cout << "x range: [" << minbbox[0] << " " << maxbbox[0] << " ]" << endl;
    cout << "y range: [" << minbbox[1] << " " << maxbbox[1] << " ]" << endl;
    cout << "z range: [" << minbbox[2] << " " << maxbbox[2] << " ]" << endl;

    if (initDataStructures()) {
      cout << "initialized" << endl;

      if (buildTree()) // computes values and splits the cells recursively.
                       // after termination the activeNodes contains all leaves
                       // at depth = maxDepth
      {
        cout << "tree built" << endl;
        cout << "Total number of cells: " << allNodes.size() << endl;
        cout << "Total number of vertices: " << allVertices.size() << endl;
        cout << "Total number of leafs: " << activeNodes.size() << endl;

        MarchingCubes();

        cout << "mesh ready" << endl;
        cout << "Total number of triangles: " << allMeshSimplices.size()
             << endl;
        cout << "Total number of mesh vertices: " << allMeshVertices.size()
             << endl;

        PrintMesh(rawFN);

        cout << "mesh printed" << endl;
      }
    }
  }

  cleanUp();

  cout << "cleanup done" << endl;
}

void Octree::cleanUp() {
  allVertices.clear();
  allNodes.clear();
  allMeshVertices.clear();
  allMeshSimplices.clear();
  while (!activeNodes.empty())
    activeNodes.pop();
  allAtoms.clear();
  delete pg;
}

bool Octree::initDataStructures() {
  int dpgxlate = minbbox[0];
  if (minbbox[1] < dpgxlate)
    dpgxlate = minbbox[1];
  if (minbbox[2] < dpgxlate)
    dpgxlate = minbbox[2];

  if (dpgxlate < 0)
    dpgxlate = -dpgxlate;

  dpgxlate += size;

  cout << "dpg xlate " << dpgxlate << endl;

  pg = new DPG::PG(5, dpgxlate, 1.75);
  for (int i = 0; i < allAtoms.size(); i++)
    pg->addPoint(allAtoms[i]);

  OTVertex *v1 = new OTVertex(minbbox[0], minbbox[1], minbbox[2]);
  computeValue(v1);
  OTVertex *v2 = new OTVertex(minbbox[0], minbbox[1], maxbbox[2]);
  computeValue(v2);
  OTVertex *v3 = new OTVertex(minbbox[0], maxbbox[1], minbbox[2]);
  computeValue(v3);
  OTVertex *v4 = new OTVertex(minbbox[0], maxbbox[1], maxbbox[2]);
  computeValue(v4);
  OTVertex *v5 = new OTVertex(maxbbox[0], minbbox[1], minbbox[2]);
  computeValue(v5);
  OTVertex *v6 = new OTVertex(maxbbox[0], minbbox[1], maxbbox[2]);
  computeValue(v6);
  OTVertex *v7 = new OTVertex(maxbbox[0], maxbbox[1], minbbox[2]);
  computeValue(v7);
  OTVertex *v8 = new OTVertex(maxbbox[0], maxbbox[1], maxbbox[2]);
  computeValue(v8);

  allVertices.push_back(v1);
  allVertices.push_back(v2);
  allVertices.push_back(v3);
  allVertices.push_back(v4);
  allVertices.push_back(v5);
  allVertices.push_back(v6);
  allVertices.push_back(v7);
  allVertices.push_back(v8);

  int vs[8];
  for (int i = 0; i < 8; i++)
    vs[i] = i;

  root = 0;
  allNodes.clear();
  allNodes.push_back(new OTCell(vs, (minbbox[0] + maxbbox[0]) / 2.0,
                                (minbbox[1] + maxbbox[1]) / 2.0,
                                (minbbox[2] + maxbbox[2]) / 2.0, size, 1));
  activeNodes.push(root);
}

void Octree::split(OTCell *cell) {
  if (!cell->isLeaf) {
    //              cout<<"Already split"<<endl;
    return;
  }

  /* coordinate values */

  float xs[3] = {allVertices[cell->vertices[0]]->xyz[0],
                 (allVertices[cell->vertices[0]]->xyz[0] +
                  allVertices[cell->vertices[7]]->xyz[0]) /
                     2.0,
                 allVertices[cell->vertices[7]]->xyz[0]};
  float ys[3] = {allVertices[cell->vertices[0]]->xyz[1],
                 (allVertices[cell->vertices[0]]->xyz[1] +
                  allVertices[cell->vertices[7]]->xyz[1]) /
                     2.0,
                 allVertices[cell->vertices[7]]->xyz[1]};
  float zs[3] = {allVertices[cell->vertices[0]]->xyz[2],
                 (allVertices[cell->vertices[0]]->xyz[2] +
                  allVertices[cell->vertices[7]]->xyz[2]) /
                     2.0,
                 allVertices[cell->vertices[7]]->xyz[2]};

  /* creating 27 new vertices (8 of which are copied from old ones) */

  int verts[3][3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if ((i | j | k) & 1) {
          verts[i][j][k] = allVertices.size();
          OTVertex *v = new OTVertex(xs[i], ys[j], zs[k]);
          computeValue(v);
          allVertices.push_back(v);

          //                                      cout<<"vertex "<<i<<" "<<j<<"
          //                                      "<<k<<" : [ "<<xs[i]<<"
          //                                      "<<ys[j]<<" "<<zs[k]<<" ]"<<"
          //                                      with id
          //                                      "<<verts[i][j][k]<<endl;
        }
      }
    }
  }

  verts[0][0][0] = cell->vertices[0];
  verts[0][0][2] = cell->vertices[1];
  verts[0][2][0] = cell->vertices[2];
  verts[0][2][2] = cell->vertices[3];
  verts[2][0][0] = cell->vertices[4];
  verts[2][0][2] = cell->vertices[5];
  verts[2][2][0] = cell->vertices[6];
  verts[2][2][2] = cell->vertices[7];

  /* Creating the children */

  for (int ii = 0; ii < 2; ii++) {
    for (int jj = 0; jj < 2; jj++) {
      for (int kk = 0; kk < 2; kk++) {
        int vs[8];

        //                              cout<<"child "<<ii<<" "<<jj<<" "<<kk<<"
        //                              : [ ";

        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
              vs[i * 4 + j * 2 + k] = verts[i + ii][j + jj][k + kk];
              //                                                      cout<<verts[i+ii][j+jj][k+kk]<<
              //                                                      " ";
            }
          }
        }

        //                              cout<<" ]"<<endl;

        cell->children[ii * 4 + jj * 2 + kk] = allNodes.size();
        allNodes.push_back(new OTCell(
            vs, (xs[ii] + xs[ii + 1]) / 2.0, (ys[jj] + ys[jj + 1]) / 2.0,
            (zs[kk] + zs[kk + 1]) / 2.0, cell->span / 2.0, cell->depth + 1));
      }
    }
  }

  for (int i = 0; i < 8; i++) {
    int childCellIndex = cell->children[i];

    if (childCellIndex != -1)
      activeNodes.push(childCellIndex);
  }

  cell->isLeaf = false;
}

void Octree::splitIfRequried(OTCell *cell) {
  if (!cell->isLeaf)
    return;

  if (cell->depth < minDepth) {
    cell->splitRequired = true;
  } else {
    bool insideFound = false;
    bool outsideFound = false;

    for (int i = 0; i < 8; i++) {
      OTVertex *v = allVertices[cell->vertices[i]];

      if (!v->valueReady)
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             << endl;

      //            cout<<v->value<<" ";
      if (v->value >= isoValue)
        insideFound = true;
      if (v->value <= isoValue)
        outsideFound = true;
    }

    //    cout<<(insideFound && outsideFound)<<endl;

    if (insideFound && outsideFound) {
      cell->splitRequired = true;
      splitNeighbors(cell);
      //    split(cell);
    }
  }

  if (cell->splitRequired)
    split(cell);
}

OTCell *Octree::locate(float x, float y, float z, int depth) {
  OTCell *result = allNodes[root];

  //      cout<<"Trying to locate "<<x<<"  "<<y<<"  "<<z<<endl;

  while (result->depth < depth && !result->isLeaf) {
    float cx = result->xyz[0];
    float cy = result->xyz[1];
    float cz = result->xyz[2];

    int childIndex = 0;

    if (cx < x)
      childIndex += 4;
    if (cy < y)
      childIndex += 2;
    if (cz < z)
      childIndex += 1;

    //              cout<<"New cell at depth "<<depth<<endl;
    //              cout<<"Position "<<cx<<"  "<<cy<<"  "<<cz<<endl;
    //              cout<<"Moving to child index "<<childIndex<<endl;
    //              cout<<"with node id
    //              "<<result->children[childIndex]<<endl<<endl;

    result = allNodes[result->children[childIndex]];
  }

  //      cout<<"Located node "<<result->xyz[0]<<"  "<<result->xyz[1]<<"
  //      "<<result->xyz[2]<<endl<<endl<<endl;

  return result;
}

void Octree::splitNeighbors(OTCell *c) {
  float cx = c->xyz[0];
  float cy = c->xyz[1];
  float cz = c->xyz[2];
  float sp = c->span;
  int dpth = c->depth;

  //      cout<<"Marking neighbors of cell at position " << cx << " " << cy << "
  //      " << cz << "   with span " << sp << "  and depth " << dpth << endl;

  // 6 face neighbors

  float change[3] = {-sp, 0, sp};

  OTCell *cell;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        cell = locate(cx + change[i], cy + change[j], cz + change[k], dpth);
        split(cell);
      }
    }
  }
}

bool Octree::buildTree() {
  int numLeaves, i;

  while (!activeNodes.empty()) {
    OTCell *cell = allNodes[activeNodes.front()];

    //              cout<<"Id of current cell: "<<activeNodes.front()<<endl;
    //              cout<<"Depth of current cell: "<<cell->depth<<endl;

    if (cell->depth == maxDepth)
      break;

    splitIfRequried(cell);

    activeNodes.pop();
  }

  return true;
}

void Octree::MarchingCubes() {
  allMeshVertices.clear();
  allMeshSimplices.clear();

  while (!activeNodes.empty()) {
    vmarchcube(allNodes[activeNodes.front()]);
    activeNodes.pop();
  }
}

void Octree::vmarchcube(OTCell *cell) {
  int iCorner, iVertex, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
  float fOffset;
  iPoint asEdgeVertex[12];

  // Find which vertices are inside of the surface and which are outside

  iFlagIndex = 0;

  int marchingCubeLayout[8] = {0, 4, 6, 2, 1, 5, 7, 3};
  int marchingCubeIndex;

  for (iVertex = 0; iVertex < 8; iVertex++) {
    marchingCubeIndex = marchingCubeLayout[iVertex];
    if (allVertices[cell->vertices[marchingCubeIndex]]->value <= isoValue)
      iFlagIndex |= 1 << iVertex;
  }

  // Find which edges are intersected by the surface
  iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

  // If the cube is entirely inside or outside of the surface, then there will
  // be no intersections
  if (iEdgeFlags == 0)
    return;

  // Find the point of intersection of the surface with each edge
  // Then find the normal to the surface at those points

  int vertind[12];

  for (iEdge = 0; iEdge < 12; iEdge++) {
    vertind[iEdge] = -1;

    if (iEdgeFlags & (1 << iEdge)) {
      int endpoint1index = marchingCubeLayout[a2iEdgeConnection[iEdge][0]];
      int endpoint2index = marchingCubeLayout[a2iEdgeConnection[iEdge][1]];

      fOffset = fGetOffset(allVertices[cell->vertices[endpoint1index]]->value,
                           allVertices[cell->vertices[endpoint2index]]->value,
                           isoValue);

      float x =
          (1 - fOffset) * allVertices[cell->vertices[endpoint1index]]->xyz[0] +
          fOffset * allVertices[cell->vertices[endpoint2index]]->xyz[0];
      float y =
          (1 - fOffset) * allVertices[cell->vertices[endpoint1index]]->xyz[1] +
          fOffset * allVertices[cell->vertices[endpoint2index]]->xyz[1];
      float z =
          (1 - fOffset) * allVertices[cell->vertices[endpoint1index]]->xyz[2] +
          fOffset * allVertices[cell->vertices[endpoint2index]]->xyz[2];

      vertind[iEdge] = allMeshVertices.size();
      MeshVertex mv(x, y, z);
      allMeshVertices.push_back(mv);
    }
  }

  // Add the triangles that were found.  There can be up to five per cube
  for (iTriangle = 0; iTriangle < 5; iTriangle++) {
    if (a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle] < 0)
      break;

    int iVertex1 =
        vertind[a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle]];
    int iVertex2 =
        vertind[a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle + 1]];
    int iVertex3 =
        vertind[a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle + 2]];

    MeshSimplex ms(iVertex1, iVertex2, iVertex3);

    allMeshSimplices.push_back(ms);
  }
}

bool Octree::computeValue(OTVertex *v) {
  float x = v->xyz[0], y = v->xyz[1], z = v->xyz[2];
  DPG::Point *p = new DPG::Point(x, y, z);

  vector<DPG::Point *> nearAtoms;

  pg->range(p, 4.6, nearAtoms);

  if (nearAtoms.size() < 0)
    return false;

  float val = 0.0;

  for (int i = 0; i < nearAtoms.size(); i++) {
    float ax = nearAtoms[i]->getX(), ay = nearAtoms[i]->getY(),
          az = nearAtoms[i]->getZ(),
          r = ((DPG::Ball *)nearAtoms[i])->getRadius();
    float distsq =
        (x - ax) * (x - ax) + (y - ay) * (y - ay) + (z - az) * (z - az);
    float expval = (blobbyness * distsq) / (r * r) - blobbyness;
    val += (exp(expval));
  }

  v->value = val;
  v->valueReady = true;

  //      cout<<"value = "<<val<<endl;

  return true;
}

void Octree::PrintMesh(string fname) {
  ofstream fout(fname.c_str());

  fout << allMeshVertices.size() << " " << allMeshSimplices.size() << endl;

  for (int i = 0; i < allMeshVertices.size(); i++) {
    fout << allMeshVertices[i].xyz[0] << " " << allMeshVertices[i].xyz[1] << " "
         << allMeshVertices[i].xyz[2] << endl;
  }
  for (int i = 0; i < allMeshSimplices.size(); i++) {
    fout << allMeshSimplices[i].verts[0] << " " << allMeshSimplices[i].verts[1]
         << " " << allMeshSimplices[i].verts[2] << endl;
  }
}

bool Octree::readPQR(string fname) {
  float minX = 100000, minY = 100000, minZ = 100000;
  float maxX = -100000, maxY = -100000, maxZ = -100000;

  ifstream input(fname.c_str());
  string atomname, resname, chain;
  int atomindex, resindex;
  float x, y, z, r;
  float q;
  string keyword;

  while (input >> keyword) {
    if (strcmp(keyword.c_str(), "ATOM") == 0) {
      input >> atomindex >> atomname >> resname >> resindex >> x >> y >> z >>
          q >> r;

      if (x < minX)
        minX = x;
      if (x > maxX)
        maxX = x;

      if (y < minY)
        minY = y;
      if (y > maxY)
        maxY = y;

      if (z < minZ)
        minZ = z;
      if (z > maxZ)
        maxZ = z;

      DPG::Ball *a = new DPG::Ball(x, y, z, r);
      allAtoms.push_back(a);
    }
  }

  float dimX = maxX - minX;
  float dimY = maxY - minY;
  float dimZ = maxZ - minZ;

  size = dimX;
  if (size < dimY)
    size = dimY;
  if (size < dimZ)
    size = dimZ;

  size += 10.0;

  minbbox[0] = minX - 5.0;
  minbbox[1] = minY - 5.0;
  minbbox[2] = minZ - 5.0;

  maxbbox[0] = minX + size;
  maxbbox[1] = minY + size;
  maxbbox[2] = minZ + size;

  cout << "Number of atoms " << allAtoms.size() << endl;

  //      cout << "x range: [" << minX << " " << maxX << " ]" << endl;
  //      cout << "y range: [" << minY << " " << maxY << " ]" << endl;
  //      cout << "z range: [" << minZ << " " << maxZ << " ]" << endl;

  //      cout << "Size: " << size << endl;

  return true;
}
