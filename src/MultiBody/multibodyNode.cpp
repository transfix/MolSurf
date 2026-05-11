#include <iostream>
#include <vector>

#include <MultiBody/multibodyNode.h>
#include <cstdio>
#include <cstring>

/** List of private variables
        int nodeId
        vector<int> multiEdgeIndices
        std::vector<double> scores
**/

MultiBodyNode::MultiBodyNode(int id, std::string nfn, ScoreWeight *sw, bool v) {
  verbose = v;
  nodeId = id;
  nodeDataReady = false;

  nodeScore = new Score(sw);

  if (readInfoFile(nfn)) {
    if (openAllFiles()) {
      std::cout << "Finished creating node " << id << "\n";

      if (verbose) {
        double s;
        if (nodeScore->getScore(0, s))
          std::cout << s << std::endl;
      }

      nodeDataReady = true;
    }
  }
}

MultiBodyNode::~MultiBodyNode() {
  delete nodeScore;
  multiEdgeIndices.clear();
}

bool MultiBodyNode::readInfoFile(std::string nodeFileName) {
  char str[300];

  FILE *fp = fopen(nodeFileName.c_str(), "rt");

  if (!fp) {
    std::cout << "Could not open node info file " << nodeFileName << "\n";
    return false;
  }

  if (fgets(str, 300, fp) != NULL) {
    std::string nstr(str);
    pointDataFileName = nstr.erase(nstr.size() - 1);
  } else {
    std::cout
        << "Error reading node info file. Failed to read point file name from "
        << nodeFileName << "\n";
    return false;
  }

  if (fgets(str, 300, fp) != NULL) {
    std::string nstr(str);
    surfaceDataFileName = nstr.erase(nstr.size() - 1);
  } else {
    std::cout
        << "Error reading node info file. Failed to surface file name from "
        << nodeFileName << "\n";
    return false;
  }

  if (fgets(str, 300, fp) != NULL) {
    std::string nstr(str);
    quadDataFileName = nstr.erase(nstr.size() - 1);
  } else {
    std::cout
        << "Error reading node info file. Failed to quad point file name from "
        << nodeFileName << "\n";
    return false;
  }

  if (fgets(str, 300, fp) != NULL) {
    std::string nstr(str);
    volumeDataFileName = nstr.erase(nstr.size() - 1);
  } else {
    std::cout
        << "Error reading node info file. Failed to volume file name from "
        << nodeFileName << "\n";
    return false;
  }

  int numScores;
  if (fscanf(fp, "%d", &numScores) != 1) {
    std::cout << "Error reading node info file. Failed to read numScores\n";
    return false;
  }

  if (numScores != nodeScore->getWeights()->getNumWeights()) {
    std::cout << "Error reading node info file. Too few scores\n";
    return false;
  }

  double sc;
  for (int i = 0; i < numScores; i++) {
    if (fscanf(fp, "%lf", &sc) != 1) {
      std::cout << "Error reading node info file. Failed to read score " << i
                << "\n";
      return false;
    }
    nodeScore->addScore(sc);
  }

  fclose(fp);

  return true;
}

bool MultiBodyNode::openAllFiles() {
  bool oneFileFailed = false;
  bool oneFileOpened = false;

  FILE *fp;

  if (strcmp("none", pointDataFileName.c_str())) {
    fp = fopen(pointDataFileName.c_str(), "rt");

    if (fp)
      oneFileOpened = true;
    else
      oneFileFailed = true;

    if (fp)
      fclose(fp);
    else
      std::cout << "Cound not open " << pointDataFileName << "\n";
  }

  if (strcmp("none", surfaceDataFileName.c_str())) {
    fp = fopen(surfaceDataFileName.c_str(), "rt");

    if (fp)
      oneFileOpened = true;
    else
      oneFileFailed = true;

    if (fp)
      fclose(fp);
    else
      std::cout << "Cound not open " << surfaceDataFileName << "\n";
  }

  if (strcmp("none", quadDataFileName.c_str())) {
    fp = fopen(quadDataFileName.c_str(), "rt");

    if (fp)
      oneFileOpened = true;
    else
      oneFileFailed = true;

    if (fp)
      fclose(fp);
    else
      std::cout << "Cound not open " << quadDataFileName << "\n";
  }

  if (strcmp("none", volumeDataFileName.c_str())) {
    fp = fopen(volumeDataFileName.c_str(), "rt");

    if (fp)
      oneFileOpened = true;
    else
      oneFileFailed = true;

    if (fp)
      fclose(fp);
    else
      std::cout << "Cound not open " << volumeDataFileName << "\n";
  }

  return oneFileOpened & (!oneFileFailed);
}

int MultiBodyNode::addEdge(int multiEdgeIndex) {
  multiEdgeIndices.push_back(multiEdgeIndex);
}

bool MultiBodyNode::updateScore(int scoreIndex, double score) {
  return nodeScore->updateScore(scoreIndex, score);
}

int MultiBodyNode::getMultiEdgeIndex(int indexNum) const {
  if (indexNum >= 0 && indexNum < multiEdgeIndices.size())
    return multiEdgeIndices[indexNum];
  return -1;
}

bool MultiBodyNode::getScore(int scoreIndex, double &score) const {
  return nodeScore->getScore(scoreIndex, score);
}

const Score *MultiBodyNode::getScores() const { return nodeScore; }
