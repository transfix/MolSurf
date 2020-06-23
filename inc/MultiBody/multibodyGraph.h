#ifndef MULTI_BODY_GRAPH_H
#define MULTI_BODY_GRAPH_H

#include <iostream>
#include <vector>

#include<MultiBody/multibodyNode.h>
#include<MultiBody/multibodyEdge.h>


class MultiBodyGraph
{
	private:
		bool verbose;
		bool comparatorsReady;
		ScoreWeight* nodeWeights; 
		ScoreWeight* edgeWeights; 
		ScoreComparator* nodeComp; 
		ScoreComparator* edgeComp;

		std::vector<MultiBodyNode*> nodes;
		std::vector<MultiBodyEdge*> edges;
		int numNodes;
		int numEdges;

		std::vector<int> sortedNodeIndices;
		std::vector<int> sortedEdgeIndices;

	public:
		MultiBodyGraph(bool verbose = false);
		MultiBodyGraph(ScoreWeight* nodeWeights, ScoreWeight* edgeWeights, ScoreComparator* nodeComp, ScoreComparator* edgeComp, bool verbose = false);
		~MultiBodyGraph();

		int getNumNodes() const;
		int getNumEdges() const;

		int addNode(string nodeFileName);
		int addMultiEdge(int src, int tgt, string xformFileName);

		MultiBodyNode* getNode(int nodeIndex);
		MultiBodyEdge* getEdge(int edgeIndex);

		void sortNodes();
		void sortMultiEdges();
};

class compareNodes
{
	ScoreComparator* nodeComp; 
	bool verbose;

	public:
		compareNodes(ScoreComparator* nc, bool v = false)
		{
			verbose = v;
			nodeComp = nc;
		}

		bool operator() (const MultiBodyNode* lhs, const MultiBodyNode* rhs) const
		{
			if(verbose) std::cout<<"Comparing nodes "<<lhs->getId()<< " and "<<rhs->getId()<<std::endl;
			return nodeComp->compare(lhs->getScores(), rhs->getScores()) >= 0;
		}
};

#endif
