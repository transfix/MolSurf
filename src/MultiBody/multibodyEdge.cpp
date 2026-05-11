#include <iostream>
#include <vector>
#include <MultiBody/multibodyEdge.h>
#include <MultiBody/weights.h>

/** List of private variables
int edgeId;
int source;
int destination;
std::set< SimpleEdge, compareSimpleEdge > simpleEdges;
bool isPartOfAssembly;
std::vector<AssemblyEdge> assemblyEdges;
**/

MultiBodyEdge::MultiBodyEdge(int src, int dest, int id,
                             std::string xformFileName, ScoreWeight *sw,
                             ScoreComparator *ec, bool v) {
  verbose = v;
  edgeId = id;
  source = src;
  destination = dest;
  numSimpleEdges = 0;
  currentSimpleEdge = 0;
  isPartOfAssembly = false;
  edgeDataReady = false;
  edgeWeights = sw;
  edgeComp = ec;

  if (verbose)
    std::cout << "reading xforms and scores from " << xformFileName << "\n";

  if (readXformFile(xformFileName)) {
    std::cout << "Finished creating edge from " << src << " to " << dest
              << "\n";

    if (verbose) {
      simpleEdges[0]->getTransformation().print();
      simpleEdges[0]->getScores()->print();
    }

    edgeDataReady = true;
  }
}

MultiBodyEdge::~MultiBodyEdge() {
  assemblyEdges.clear();
  simpleEdges.clear();
  sortedEdgeIndices.clear();
}

bool MultiBodyEdge::readXformFile(string xformFileName) {
  char str[300];
  int numXforms;
  int numScores;

  FILE *fp = fopen(xformFileName.c_str(), "rt");

  if (!fp) {
    std::cout << "Could not open xform file " << xformFileName << "\n";
    return false;
  }

  if (fscanf(fp, "%d", &numXforms) != 1) {
    std::cout << "Error reading xform file. Failed to read numXforms\n";
    return false;
  }

  if (fscanf(fp, "%d", &numScores) != 1) {
    std::cout << "Error reading xform file. Failed to read numScores\n";
    return false;
  }

  double sc;
  double t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31,
      t32, t33;

  for (int i = 0; i < numXforms; i++) {
    if (fscanf(
            fp,
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &t00, &t01, &t02, &t03, &t10, &t11, &t12, &t13, &t20, &t21, &t22,
            &t23, &t30, &t31, &t32, &t33) != 16) {
      std::cout << "Error reading xform file. Failed to read transformation "
                << i << "\n";
      return false;
    }
    CCVOpenGLMath::Matrix T(t00, t01, t02, t03, t10, t11, t12, t13, t20, t21,
                            t22, t23, t30, t31, t32, t33);

    vector<double> scores;

    for (int j = 0; j < numScores; j++) {
      if (fscanf(fp, "%lf", &sc) != 1) {
        std::cout << "Error reading xform file. Failed to read score " << i
                  << " " << j << "\n";
        return false;
      }
      scores.push_back(sc);
    }

    addSimpleEdge(T, scores);
  }

  fclose(fp);

  return true;
}

void MultiBodyEdge::addSimpleEdge(CCVOpenGLMath::Matrix mtx,
                                  vector<double> &scores) {
  SimpleEdge *simpedge = new SimpleEdge(numSimpleEdges, source, destination,
                                        edgeWeights, scores, mtx, verbose);
  simpleEdges.push_back(simpedge);
  numSimpleEdges++;
}

void MultiBodyEdge::makePartOfAssembly(int assemblyId, int simpleEdgeIndex) {
  AssemblyEdge ae(assemblyId, simpleEdgeIndex);
  assemblyEdges.push_back(ae);
}

int MultiBodyEdge::getAssemblyEdgeIndex(int assemblyId) const {
  int numAssemblyEdges = assemblyEdges.size();

  for (int i = 0; i < numAssemblyEdges; i++) {
    if (assemblyEdges[i].getAssemblyId() == assemblyId) {
      return assemblyEdges[i].getSelectedSimpleEdgeIndex();
    }
  }
  return -1;
}

bool MultiBodyEdge::getTransformation(int simpleEdgeIndex,
                                      CCVOpenGLMath::Matrix &T) const {
  int numSimpleEdges = simpleEdges.size();

  if (simpleEdgeIndex < 0 || simpleEdgeIndex >= numSimpleEdges)
    return false;

  T = simpleEdges[simpleEdgeIndex]->getTransformation();
  return true;
}

bool MultiBodyEdge::getScore(int simpleEdgeIndex, int scoreIndex,
                             double &score) const {
  int numSimpleEdges = simpleEdges.size();

  if (simpleEdgeIndex < 0 || simpleEdgeIndex >= numSimpleEdges)
    return false;

  return simpleEdges[simpleEdgeIndex]->getScore(scoreIndex, score);
}

void MultiBodyEdge::sortSimpleEdges() {
  std::vector<SimpleEdge *> edgestemp = simpleEdges;

  std::sort(edgestemp.begin(), edgestemp.end(),
            compareSimpleEdges(edgeComp, false));

  sortedEdgeIndices.clear();
  for (int i = 0; i < simpleEdges.size(); i++)
    sortedEdgeIndices.push_back(edgestemp[i]->getId());

  if (verbose) {
    for (int i = 0; i < 10; i++) {
      int index = sortedEdgeIndices[i];
      std::cout << "Rank " << i << " Edge = " << index << std::endl;
      simpleEdges[index]->getScores()->print();
    }
  }
}
