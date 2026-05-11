#ifndef MULTI_BODY_EDGE_H
#define MULTI_BODY_EDGE_H

#include <iostream>
#include <set>
#include <vector>
#include <UsefulMath/Matrix.h>
#include <MultiBody/scoreComparator.h>

class SimpleEdge {
private:
  bool verbose;
  int simpleEdgeId;
  int source;
  int destination;
  Score *edgeScore;
  CCVOpenGLMath::Matrix transformation;

public:
  SimpleEdge(int id, int src, int dest, ScoreWeight *sw, vector<double> &scores,
             CCVOpenGLMath::Matrix xform, bool v = false) {
    verbose = v;
    simpleEdgeId = id;
    source = src;
    destination = dest;
    transformation = xform;
    edgeScore = new Score(sw);

    int numScores = scores.size();

    if (sw == NULL || scores.size() == 0)
      std::cout << "Error here\n";

    for (int i = 0; i < numScores; i++)
      edgeScore->addScore(scores[i]);
  }

  ~SimpleEdge() { delete edgeScore; }

  int getId() const { return simpleEdgeId; }

  int getSource() const { return source; }

  int getDestination() const { return destination; }

  const Score *getScores() const { return edgeScore; }

  bool getScore(int scoreIndex, double &score) {
    return edgeScore->getScore(scoreIndex, score);
  }

  CCVOpenGLMath::Matrix getTransformation() const { return transformation; }

  bool updateScore(int scoreIndex, double score) {
    return edgeScore->updateScore(scoreIndex, score);
  }
};

class compareSimpleEdges {
  ScoreComparator *edgeComp;
  bool verbose;

public:
  compareSimpleEdges(ScoreComparator *ec, bool v = false) {
    verbose = v;
    edgeComp = ec;
  }

  int compare(SimpleEdge *lhs, SimpleEdge *rhs) const {
    return edgeComp->compare(lhs->getScores(), rhs->getScores());
  }

  bool operator()(const SimpleEdge *lhs, const SimpleEdge *rhs) const {
    if (verbose)
      std::cout << "Comparing edges " << lhs->getId() << " and " << rhs->getId()
                << std::endl;
    return edgeComp->compare(lhs->getScores(), rhs->getScores()) >= 0;
  }
};

class AssemblyEdge {
private:
  int assemblyId;
  int selectedSimpleEdgeIndex;

public:
  AssemblyEdge(int aid, int ssei) {
    assemblyId = aid;
    selectedSimpleEdgeIndex = ssei;
  }
  ~AssemblyEdge() {}

  int getAssemblyId() const { return assemblyId; }
  int getSelectedSimpleEdgeIndex() const { return selectedSimpleEdgeIndex; }

  void setAssemblyId(int ai) { assemblyId = ai; }
  void setSelectedSimpleEdgeIndex(int ssei) { selectedSimpleEdgeIndex = ssei; }
};

/*
class compareSimpleEdge
{
        bool reverse;

        public:
                compareSimpleEdge(const bool& revparam=false)
                {
                        reverse=revparam;
                }

                bool operator() (const SimpleEdge& lhs, const SimpleEdge& rhs)
const
                {
                        double l = lhs.getTotalScore();
                        double r = rhs.getTotalScore();
                        if (reverse) return (l>=r);
                        else return (l<r);
                }
};
*/

class MultiBodyEdge {
private:
  bool verbose;
  int edgeId;
  int source;
  int destination;

  ScoreWeight *edgeWeights;
  ScoreComparator *edgeComp;

  int numSimpleEdges;
  std::vector<SimpleEdge *> simpleEdges;
  std::vector<int> sortedEdgeIndices;
  int currentSimpleEdge;

  bool isPartOfAssembly;
  std::vector<AssemblyEdge> assemblyEdges;

  bool edgeDataReady;

  bool readXformFile(string xformFileName);

public:
  MultiBodyEdge(int nodeIndex1, int nodeIndex2, int id, string xformFileName,
                ScoreWeight *sw, ScoreComparator *edgeComp,
                bool verbose = false);
  ~MultiBodyEdge();

  int getId() const { return edgeId; }

  bool isReady() const { return edgeDataReady; }

  void sortSimpleEdges();

  void addSimpleEdge(CCVOpenGLMath::Matrix mtx, vector<double> &);

  int getSourceNodeIndex() const { return source; }
  int getDestinationNodeIndex() const { return destination; }

  void makePartOfAssembly(int assemblyId, int simpleEdgeIndex);
  bool isPartOfAnyAssembly() const { return isPartOfAssembly; }
  int getAssemblyEdgeIndex(int assemblyId) const;

  int getNumEdges() const { return simpleEdges.size(); }
  bool getTransformation(int simpleEdgeIndex, CCVOpenGLMath::Matrix &T) const;
  bool getScore(int simpleEdgeIndex, int scoreIndex, double &score) const;

  SimpleEdge *getBestSimpleEdge() const {
    return simpleEdges[sortedEdgeIndices[0]];
  }
  SimpleEdge *getNextBestSimpleEdge() {
    return simpleEdges[sortedEdgeIndices[currentSimpleEdge++]];
  } // increments currentSimpleEdge

  int getIndexOfNextBestSimpleEdge() const {
    return sortedEdgeIndices[currentSimpleEdge];
  }
  SimpleEdge *getSimpleEdge(int index) const {
    if (index >= 0 && index < numSimpleEdges)
      return simpleEdges[index];
    return NULL;
  }
  void moveToNextBestSimpleEdge() { currentSimpleEdge++; }
};

class compareMultiEdges {
  compareSimpleEdges *compSimpEdges;
  bool verbose;

public:
  compareMultiEdges(compareSimpleEdges *cse, bool v = false) {
    compSimpEdges = cse;
    verbose = v;
  }

  bool operator()(const MultiBodyEdge *lhs, const MultiBodyEdge *rhs) const {
    return compSimpEdges->compare(lhs->getBestSimpleEdge(),
                                  rhs->getBestSimpleEdge()) >= 0;
  }
};

#endif
