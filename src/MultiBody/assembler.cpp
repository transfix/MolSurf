
#include<MultiBody/assembler.h>


Assembler::Assembler(bool v)
{
	verbose = v;
	graph = new MultiBodyGraph(verbose);
}


Assembler::~Assembler()
{
	delete graph;
}


bool Assembler::makeGraph(string nodeDataFileName, string edgeDataFileName, string nodeScoreSettingsFileName, string edgeScoreSettingsFileName)
{
	int numNodes;
	int numMultiEdges;
	int source, target;
	char str[300];

	FILE* fp = fopen(nodeScoreSettingsFileName.c_str(), "rt");
	if(!fp)
	{
		std::cout<<"Could not open nodeScoreSettings file\n";
		return false;
	}

	int mode;
	if(fscanf(fp, "%d\n", &mode) != 1)
	{
		std::cout<<"Error reading nodeScoreSettings file. Failed to read mode\n";
		return false;
	}

	ScoreWeight* nodeWeights = new ScoreWeight();

	if( mode == WEIGHTEDSUM || mode == MAJORITY || mode == LEXICOGRAPHIC )
	{
		int numWeights;
		if(fscanf(fp, "%d\n", &numWeights) != 1)
		{
			std::cout<<"Error reading nodeScoreSettings file. Failed to read numWeights\n";
			return false;
		}

		for(int i=0; i<numWeights; i++)
		{
			double w;
			if(fscanf(fp, "%lf", &w) != 1)
			{
				std::cout<<"Error reading nodeScoreSettings file. Failed to read weight\n";
				return false;
			}
			nodeWeights->addWeight(w);
		}
	}

	ScoreComparator* nodeComp = new ScoreComparator(mode, false);

	fclose(fp);

	std::cout<<"Finished reading node weights file\n";
	if(verbose) nodeWeights->print();



	fp = fopen(edgeScoreSettingsFileName.c_str(), "rt");
	if(!fp)
	{
		std::cout<<"Could not open edgeScoreSettings file\n";
		return false;
	}

	if(fscanf(fp, "%d\n", &mode) != 1)
	{
		std::cout<<"Error reading edgeScoreSettings file. Failed to read mode\n";
		return false;
	}

	ScoreWeight* edgeWeights = new ScoreWeight();

	if( mode == WEIGHTEDSUM || mode == MAJORITY || mode == LEXICOGRAPHIC )
	{
		int numWeights;
		if(fscanf(fp, "%d\n", &numWeights) != 1)
		{
			std::cout<<"Error reading edgeScoreSettings file. Failed to read numWeights\n";
			return false;
		}

		for(int i=0; i<numWeights; i++)
		{
			double w;
			if(fscanf(fp, "%lf", &w) != 1)
			{
				std::cout<<"Error reading edgeScoreSettings file. Failed to read weight\n";
				return false;
			}
			edgeWeights->addWeight(w);
		}
	}

	ScoreComparator* edgeComp = new ScoreComparator(mode, false);

	fclose(fp);

	std::cout<<"Finished reading edge weights file\n";
	if(verbose) edgeWeights->print();
	
	graph = new MultiBodyGraph(nodeWeights, edgeWeights, nodeComp, edgeComp, verbose);

	fp = fopen(nodeDataFileName.c_str(), "rt");

	if(!fp)
	{
		std::cout<<"Could not open nodedata file\n";
		return false;
	}

	if(fscanf(fp, "%d\n", &numNodes) != 1)
	{
		std::cout<<"Error reading nodedata file. Failed to read numNodes\n";
		return false;
	}

	
	for(int i=0; i<numNodes; i++)
	{
		if ( fgets (str , 300 , fp) != NULL )
		{
			std::string nstr(str);
			graph->addNode(nstr.erase(nstr.size()-1));
		}
		else
		{
			std::cout<<"Error reading nodedata file. Failed to read nfo file name "<<i<<"\n";
			return false;
		}
	}

	fclose(fp);

	std::cout<<"Finished reading node data file\n";



	fp = fopen(edgeDataFileName.c_str(), "rt");

	if(!fp)
	{
		std::cout<<"Could not open edgedata file\n";
		return false;
	}

	if(fscanf(fp, "%d", &numMultiEdges) != 1)
	{
		std::cout<<"Error reading edgedata file. Failed to read numMultiEdges\n";
		return false;
	}

	for(int i=0; i<numMultiEdges; i++)
	{
		if(fscanf(fp, "%d %d ", &source, &target) != 2)
		{
			std::cout<<"Error reading edgedata file. Failed to read source and target\n";
			return false;
		}

		if(source < 0 || source >=numNodes || target < 0 || target>= numNodes)
		{
			std::cout<<"Error in edgedata file. Source or target out of bounds\n";
			return false;
		}

		if ( fgets (str , 300 , fp) != NULL )
		{
			std::string nstr(str);
			graph->addMultiEdge(source, target, nstr.erase(nstr.size()-1));
		}
		else
		{
			std::cout<<"Error reading edgedata file. Failed to read nfo file name "<<i<<"\n";
			return false;
		}
	}

	fclose(fp);

	std::cout<<"Finished reading edge data file\n";

	return true;
}


bool Assembler::assemble()
{
	graph->sortNodes();
	graph->sortMultiEdges();
}
