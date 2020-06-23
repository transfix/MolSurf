#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cmath>
using namespace std;

class Octree;

//
// future: remove almost all the information from
//         each cell since I can compute the location 
//         and size during the tree traversal,
//         that should halve to total memory usage
//         for a minimal time cost

class Cell {
 public:
  int vertices[8];
  double values[8];
  double span;
  void Print();
};

class Vertex {
 public:
  double xyz[3];
};


class PairIntersection
{
	public:
	PairIntersection() {}
	int atomA;
	int atomB;
	double cx,cy,cz, r;	//center of circle of intersection
	double A, B, C, D;	//plane containing the circle
	double thetaA;		//angle subtended by the cone from center of atom A onto the circle
	double thetaB;		//angle subtended by the cone from center of atom B onto the circle
};

class TripleIntersection
{
	public:
	TripleIntersection() {}
	int atomA;
	int atomB;
	int atomC;
	double x1, y1, z1;	//point 1
	double x2, y2, z2;	//point 2
	bool isBuried1;
	bool isBuried2;
};

class Atom {
 public:
  double xyz[3];
  double r;  
  void Print();
  bool copyDisabled;

  Atom(){copyDisabled = false; pairIntersections.clear(); tripleIntersections.clear();}

  double dist(Atom* a) { return sqrt ((xyz[0]-a->xyz[0])*(xyz[0]-a->xyz[0]) + (xyz[1]-a->xyz[1])*(xyz[1]-a->xyz[1]) + (xyz[2]-a->xyz[2])*(xyz[2]-a->xyz[2])); }
  double pointDist(double x, double y, double z, double R) {return ((xyz[0]-x)*(xyz[0]-x) + (xyz[1]-y)*(xyz[1]-y) + (xyz[2]-z)*(xyz[2]-z) - (r+R)*(r+R));}
  bool intersect(Atom* a, double solventRad) { return dist(a) < (r + a->r + 2*solventRad);}
  bool isTripleIntersect(Atom* a1, Atom* a2, Atom* a3, double solventRad);

  vector<PairIntersection *> pairIntersections;
  vector<TripleIntersection *> tripleIntersections;
};


class OTNode{
 public:
  OTNode ** c;
  OTNode *parent;
  int l;         // cell depth
  double xyz[3]; // cell center
  double w;      // cell width
  double v;      // data
  double sdf;

  int flag;
  bool interfaceFlag;
  bool splitDisabled;
  bool sdfComputed;
  bool markedForLocalUpdate;

  OTNode ();
  OTNode (int l1, double x, double y, double z, double w1);
  void getInterfaceArea(OTNode *otherNode, double* trans, double interfaceWidth, double& areaInner, double& areaOuter, double isoValue);

  void Split();
  void Merge(std::vector<Atom*>& atoms);
  void Print();

  void SetData(double v1);

  bool isBuried(vector<Atom*>& atoms, double *point, double solventRad);

  std::vector<int> atoms;

  static Octree * tree;
  Cell* dual;

 private:

};




typedef std::pair<int,int> Edge;

class Octree {
 public:
  //  Octree(int n, double w, double x, double y, double z);

  Octree(int n, std::string fname);
  Octree(int n, std::string fname,  std::string fname2);
  void constructOctree(int targLev, double isoVal, double bl);

  inline OTNode * GetNode(int i, int j, int k);
  OTNode * BaseLocate(double * xyz);
  OTNode * Locate(double * xyz);
  OTNode * Locate(double * xyz, int depth);

  double getInterfaceArea(Octree *otherTree, double* trans, double interfaceWidth);
  double getInterfaceAreaUnderRoot(OTNode* root, Octree *otherTree, double* trans, double interfaceWidth);

  void GetAffectedCells(Atom *a, std::vector<OTNode*>& intersectedCells);
  void LocateAffectedCells(OTNode * nd, Atom *a, std::vector<OTNode*>& intersectedCells);

  void AddAtom(int id) {AddAtom(atoms[id]);}
  void RemoveAtom(int id){RemoveAtom(atoms[id]);}
  void AddAtom(Atom *a);
  void RemoveAtom(Atom *a);
  void AddAtoms(vector<Atom*> atms);
  void RemoveAtoms(vector<Atom*> atms);
  void LocalRefinement(std::vector<OTNode *> cellsToBeUpdated);

  bool mergeRequired(OTNode * nd);
  bool splitRequired(OTNode * nd);

  void Neighbors(OTNode * nd, OTNode ** nbs);// nbs must be an array of length 6
  void GetExtendedNeighbors(OTNode*, OTNode**);

  bool isBuried(double *point);

  void computeSDF(OTNode *nd, double solventRad);

  bool getPairwiseIntersection(int atomA, int atomB, double x, double y, double z, double solventRad, PairIntersection *PI);
  bool getSDFForIntersectingSpheres(PairIntersection* PI, double x, double y, double z, double solventRad,  double &sdf);

  bool computeTripleIntersection(int atomA,  int atomB, int atomC, double x, double y, double z, double solventRad, TripleIntersection * TI);
  bool getSDFForTripleIntersect(TripleIntersection* TI, double x, double y, double z, double &sdf);

  void MarkIsoContour();
  void MarkIsoContour(OTNode * nd);
  void MarkExtendedNeighbors(OTNode * nd);
  void MarkExtendedNeighbors2(OTNode * nd);

  void SplitMarkedCells();
  void SplitMarkedCells(OTNode * nd);

  bool ComputeNewCellValues();

  void ComputeFinalValues();
  void ComputeFinalValues(OTNode * nd);

  void LabelFinalVertices();
  void LabelFinalVertices(OTNode * nd);

  void GetFinalCellsAndDestroyTree();
  void GetFinalCellsAndDestroyTree(OTNode * nd);
  void GetFinalCell(OTNode * nd,OTNode ** nbs);

  void MarchingCubes();
  void vmarchcube(Cell* cell, float  iso);

  void ClearIntermediateData();
  void PrintMesh(char * fname);

  int basesize;
  double basewidth;
  OTNode * t;
  double base[3];

  void SetInitialData();

  double contourValue;
  double targetLevel;
  double targetRes;
  double blobby;

  void ReadPQR(std::string fname);
  void ReadPQR(std::string fname,std::string fname2);



  bool DoesOverlap(OTNode * nd, int iA); 
  bool DoesOverlap(OTNode * nd, Atom* a);


  double fvdw(OTNode * nd);
  double fgauss(OTNode * nd);
  void fvdw_update_add(OTNode * nd, Atom* a);
  void fgauss_update_add(OTNode * nd, Atom* a);
  void fgauss_update_remove(OTNode * nd, Atom* a);


  double computeMarkedArea();
  double computeMarkedArea(OTNode* nd);
  double areaMesh();
  

  std::vector<Cell> finalCells;
  int paircount;
  int paircount1;
  int triplecount;
  int triplecount1;
  int tripleprob;
  int pairprob;
  int numoverwrite;

  std::vector<Atom *> atoms;
  std::vector<OTNode *> uncomputedCells;
  std::vector<OTNode *> newNodes;
  std::vector<OTNode *> markedNodes;
 
  std::vector<Vertex> finalVertices;

  std::vector<double> vertices;
  std::vector<int> triangles;
  
  std::map<Edge,int> edgeToVert;

  int numSDFCalled;
  int numPairFound;
  int numOldPairSearched;
  int numOldPairUsed;
  int numPairComputed;
  int numTripleFound;
  long numOldTripleSearched;
  int numOldTripleUsed;
  int numTripleComputed;

  void clear();

 private:
 
  double center[3];
  double size;

};

inline OTNode * Octree::GetNode(int i, int j, int k) {
  return &t[basesize*basesize*i + basesize*j + k];

}
