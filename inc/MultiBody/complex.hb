#ifndef MULTI_COMPLEX_H
#define MULTI_COMPLEX_H

#include <iostream>
#include <vector>

#include<MultiBody/multibodyGraph.h>

class Pose
{
	
}

class Complex
{
	private:
		double score;
		std::vector<SimpleEdge*> selectedEdges;

		void updateScore(int selectedEdgeIndex);	// updates the score based on the interactions introduced 
								// by the edge corresponding to selectedEdgeIndex with 
								// all edges coressponding to 0 ... selectedEdgeIndex-1

		void computeScore();

	public:
		Complex();
		~Complex();

		bool makeComplex(std::vector<SimpleEdge*> selectedEdges);
		bool addToComplex(SimpleEdge* selectedEdge);

		void reinitialize();

		double getScore();
};

#endif
