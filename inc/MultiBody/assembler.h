#ifndef MULTI_ASSEMBLER_H
#define MULTI_ASSEMBLER_H

#include <iostream>
#include <vector>

#include<MultiBody/multibodyGraph.h>

//#include<MultiBody/complex.h>
//#include<MultiBody/chainTransform.h>


class Assembler
{
	private:
		bool verbose;
		MultiBodyGraph* graph;

		bool readScoreComparatorSettings(string scoreSettingsFileName, int& mode, ScoreWeight& sw);
 
	public:
		Assembler(bool verbose = false);
		~Assembler();

		bool makeGraph(string nodeDataFileName, string edgeDataFileName, string nodeScoreSettingsFileName, string edgeScoreSettingsFileName);
		bool assemble();
};

#endif
