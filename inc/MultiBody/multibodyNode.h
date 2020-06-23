#ifndef MULTI_BODY_NODE_H
#define MULTI_BODY_NODE_H

#include <vector>
#include <string>
#include<MultiBody/scoreComparator.h>

class MultiBodyNode
{
	private:
		bool verbose;
		int nodeId;
		bool nodeDataReady;

		std::string pointDataFileName;
		std::string surfaceDataFileName;
		std::string quadDataFileName;
		std::string volumeDataFileName;

		Score* nodeScore;

		std::vector<int> multiEdgeIndices;


		bool openAllFiles();
		bool readInfoFile(std::string nodeFileName);

		
	public:
		MultiBodyNode(int id, std::string nodeFileName, ScoreWeight* sw, bool verbose = false);
		~MultiBodyNode();

		int getId() const {return nodeId;}

		std::string getPointDataFileName() {return pointDataFileName;}
		std::string getSurfaceDataFileName() {return surfaceDataFileName;}
		std::string getQuadDataFileName() {return quadDataFileName;}
		std::string getVolumeDataFileName() {return volumeDataFileName;}

		int addEdge(int multiEdgeIndex);
		bool updateScore(int scoreIndex, double score);
		bool isReady() {return nodeDataReady;}

		int getMultiEdgeIndex(int indexNum) const;
		bool getScore(int scoreIndex, double& score) const;
		const Score* getScores() const;
};

#endif
