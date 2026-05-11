#include <iostream>
#include <vector>
#include <algorithm>

#include <MultiBody/multibodyGraph.h>

MultiBodyGraph::MultiBodyGraph(bool v) {
  verbose = v;
  numNodes = 0;
  numEdges = 0;
  comparatorsReady = false;
}

MultiBodyGraph::MultiBodyGraph(ScoreWeight *nw, ScoreWeight *ew,
                               ScoreComparator *nc, ScoreComparator *ec,
                               bool v) {
  verbose = v;
  numNodes = 0;
  numEdges = 0;
  nodeWeights = nw;
  edgeWeights = ew;
  nodeComp = nc;
  edgeComp = ec;
  comparatorsReady = true;
}

MultiBodyGraph::~MultiBodyGraph() {
  nodes.clear();
  edges.clear();
}

int MultiBodyGraph::addNode(std::string nodeFileName) {
  MultiBodyNode *node =
      new MultiBodyNode(numNodes, nodeFileName, nodeWeights, verbose);
  nodes.push_back(node);
  return ++numNodes;
}

int MultiBodyGraph::addMultiEdge(int src, int tgt, string xformFileName) {
  MultiBodyEdge *edge = new MultiBodyEdge(src, tgt, numEdges, xformFileName,
                                          edgeWeights, edgeComp, verbose);
  edges.push_back(edge);
  nodes[src]->addEdge(numEdges);
  return ++numEdges;
}

MultiBodyNode *MultiBodyGraph::getNode(int nodeIndex) {
  if (nodeIndex < 0 || nodeIndex >= numNodes)
    return NULL;

  return nodes[nodeIndex];
}

MultiBodyEdge *MultiBodyGraph::getEdge(int edgeIndex) {
  if (edgeIndex < 0 || edgeIndex >= numEdges)
    return NULL;

  return edges[edgeIndex];
}

void MultiBodyGraph::sortNodes() {
  std::vector<MultiBodyNode *> nodestemp = nodes;

  std::sort(nodestemp.begin(), nodestemp.end(),
            compareNodes(nodeComp, verbose));

  sortedNodeIndices.clear();
  for (int i = 0; i < nodes.size(); i++)
    sortedNodeIndices.push_back(nodestemp[i]->getId());

  if (verbose) {
    for (int i = 0; i < nodes.size(); i++) {
      int index = sortedNodeIndices[i];
      std::cout << "Rank " << i << " Node = " << index << std::endl;
      nodes[index]->getScores()->print();
    }
  }
}

void MultiBodyGraph::sortMultiEdges() {
  if (verbose)
    std::cout << "Sorting " << numEdges << " edges\n";

  for (int i = 0; i < numEdges; i++)
    edges[i]->sortSimpleEdges();

  std::vector<MultiBodyEdge *> edgestemp = edges;

  std::sort(edgestemp.begin(), edgestemp.end(),
            compareMultiEdges(new compareSimpleEdges(edgeComp, false), false));

  sortedEdgeIndices.clear();
  for (int i = 0; i < edges.size(); i++)
    sortedEdgeIndices.push_back(edgestemp[i]->getId());

  char chr[8] = "ABCDEFG";

  if (verbose) {
    for (int i = 0; i < edges.size(); i++) {
      int index = sortedEdgeIndices[i];
      std::cout << "Rank " << i << " MultiEdge = " << index << ".   "
                << chr[edges[index]->getSourceNodeIndex()] << " - "
                << chr[edges[index]->getDestinationNodeIndex()]
                << " with best simple edge index "
                << edges[index]->getIndexOfNextBestSimpleEdge() << std::endl;
      // edges[index]->getBestSimpleEdge()->getScores()->print();
    }
  }
}
