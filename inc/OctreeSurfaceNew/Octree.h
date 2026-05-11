#include <vector>
#include <queue>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <DPG/PG.h>

using namespace std;

class Octree;

//
// future: remove almost all the information from
//         each cell since I can compute the location
//         and size during the tree traversal,
//         that should halve to total memory usage
//         for a minimal time cost

class OTVertex {
public:
  OTVertex(float a, float b, float c) {
    xyz[0] = a;
    xyz[1] = b;
    xyz[2] = c;
    value = 0.0;
    valueReady = false;
  }
  ~OTVertex();

  float xyz[3];
  float value;
  bool valueReady;
};

class OTCell {
public:
  float xyz[3];
  int vertices[8];
  float span;
  int depth;
  int children[8];
  bool isLeaf;
  bool splitRequired;

  OTCell() {
    isLeaf = true;
    splitRequired = false;
    for (int i = 0; i < 8; i++) {
      vertices[i] = -1;
      children[i] = -1;
      span = 0.0;
    }
  }
  OTCell(int *vs, float a, float b, float c, float sp, int d) {
    isLeaf = true;
    splitRequired = false;
    xyz[0] = a;
    xyz[1] = b;
    xyz[2] = c;
    span = sp;
    depth = d;
    for (int i = 0; i < 8; i++) {
      vertices[i] = vs[i];
      children[i] = -1;
    }
  }
};

class MeshVertex {
public:
  MeshVertex(float x, float y, float z) {
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  float xyz[3];
};

class MeshSimplex {
public:
  MeshSimplex(float a, float b, float c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
  int verts[3];
};

class Octree {
private:
  vector<OTVertex *> allVertices;
  vector<OTCell *> allNodes;

  vector<MeshVertex> allMeshVertices;
  vector<MeshSimplex> allMeshSimplices;

  int root;
  queue<int> activeNodes;

  vector<DPG::Ball *> allAtoms;

  DPG::PG *pg;

  float minbbox[3];
  float maxbbox[3];
  float size;

  //              int depth;
  int minDepth;
  int maxDepth;

  float isoValue;
  float blobbyness;

  void cleanUp();

  bool readPQR(string fname);

  bool initDataStructures();
  bool computeValue(OTVertex *v);
  bool buildTree();
  OTCell *locate(float x, float y, float z, int depth);
  void splitNeighbors(OTCell *cell);
  void splitIfRequried(OTCell *cell);
  void split(OTCell *cell);

  void MarchingCubes();
  void vmarchcube(OTCell *cell);
  void PrintMesh(string fname);

public:
  Octree(std::string pqrFN, std::string rawFN, int minDepth, int maxDepth,
         float blobbyness, float isovalue);
};
