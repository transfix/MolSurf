#ifndef MULTI_CHAIN_XFORM_H
#define MULTI_CHAIN_XFORM_H

#include <iostream>
#include <vector>

#include<MultiBody/multibodyNode.h>
#include<MultiBody/multibodyEdge.h>


class ChainTransform
{
	private:
		std::vector<SimpleEdge*> simpleEdgeChain;

	public:
		ChainTransform();
		~ChainTransform();

		bool addSimpleEdge(SimpleEdge* se);
		void reset();
		CCVOpenGLMath::Matrix getFinalTransform();
};

#endif
