#ifndef MULTI_BODY_COMPARATOR_H
#define MULTI_BODY_COMPARATOR_H

#include <iostream>
#include <set>
#include <vector>

enum 
{ 	
	WEIGHTEDSUM,
	MAJORITY, 
	LEXICOGRAPHIC, 
	GROUPWEIGHTEDSUM, 
	HIERARCHICAL
};

class ScoreWeight
{
	private:
		std::vector<double> weights;
	public:
		ScoreWeight() {}

		~ScoreWeight() {weights.clear();}

		void addWeight(double w) {weights.push_back(w);}	

		const int getNumWeights() const {return weights.size();}

		void getWeights(std::vector<double>& ws) const { ws = weights;}

		bool getWeight(int index, double& w) const
		{
			if(index>=0 && index<weights.size()) 
			{
				w = weights[index];
				return true;
			}
			w = 0.0;
			return false;
		}

		bool updateWeight(int index, double w) 
		{
			if(index>=0 && index<weights.size()) 
			{
				weights[index] = w;
				return true;
			}
			return false;
		}

		bool removeWeight()
		{
			if(weights.size() > 0)
			{
				weights.pop_back();
				return true;
			}
			return false;
		}

		bool removeWeight(int index) 
		{
			if(index>=0 && index<weights.size()) 
			{
				weights.erase(weights.begin()+index);
				return true;
			}
			return false;
		}

		void print() const
		{
			std::cout<<"Number of weights: "<<weights.size()<<"\n";
			
			for(int i=0; i<weights.size(); i++)
				std::cout<<weights[i]<<" ";
			std::cout<<"\n\n";
		}
};

/*
class HierarchicalWeight
{
	private:
		ScoreWeight* weightsForThisLayer;
		std::vector<HierarchicalWeight*> weightsForNextLayer;
		bool isLeaf;

	public:
		HierarchicalWeight(ScoreWeight* w) {weightsForThisLayer = w; isLeaf = true;}

		~HierarchicalWeight() {weightsForNextLayer.clear();}

		void setWeightsForThisLayer(ScoreWeight* w) { weightsForThisLayer = w; }

		void addWeightsForNextLayer(HierarchicalWeight* w) { weightsForNextLayer.push_back(w); isLeaf = false;}

		bool getWeightsForNextLayer(int index, HierarchicalWeight* hw)
		{
			if(index>=0 && index<weightsForNextLayer.size()) 
			{
				hw = weightsForNextLayer[index];
				return true;
			}
			hw = NULL;
			return false;	
		}

		ScoreWeight* getWeightForThisLayer() { return weightsForThisLayer; }
}
*/

class Score
{
	private:
		std::vector<double> scores;
		ScoreWeight* weights;

	public:
		Score(ScoreWeight* w) {weights = w;}

		~Score() {scores.clear();}

		void setWeights(ScoreWeight* w) {weights = w;}

		const ScoreWeight* getWeights() const {return weights;}

		void addScore(double s) {scores.push_back(s);}	

		const int getNumScores() const {return scores.size();}

		void getScores(std::vector<double>& ss) const { ss = scores;}

		bool getScore(int index, double& s) const
		{
			if(index>=0 && index<scores.size()) 
			{
				s = scores[index];
				return true;
			}
			s = 0.0;
			return false;
		}

		bool updateScore(int index, double s) 
		{
			if(index>=0 && index<scores.size()) 
			{
				scores[index] = s;
				return true;
			}
			return false;
		}

		bool removeScore()
		{
			if(scores.size() > 0)
			{
				scores.pop_back();
				return true;
			}
			return false;
		}

		bool removeScore(int index) 
		{
			if(index>=0 && index<scores.size()) 
			{
				scores.erase(scores.begin()+index);
				return true;
			}
			return false;
		}

		void print() const
		{
			std::cout<<"Number of scores: "<<scores.size()<<"\n";
			
			for(int i=0; i<scores.size(); i++)
				std::cout<<scores[i]<<" ";
			std::cout<<"\n\n";
		}
};


class ScoreComparator
{
	private:
		int mode;	
		
		bool verbose;

		int compareUsingWeightedSum(const Score* S1, const Score* S2) const;

		int compareUsingMajority(const Score* S1, const Score* S2) const;

		int compareUsingLexicographicOrder(const Score* S1, const Score* S2) const;

		int compareUsingGroupWeightedSum(const Score* S1, const Score* S2) const;

	public:
		ScoreComparator(int m, bool verbose = false);
		~ScoreComparator() {}

		int compare(const Score* S1, const Score* S2) const; 	// 0 = equal, -1 = S1 < S2, +1 = S1 > S2

		bool operator() (const Score* lhs, const Score* rhs) const	// lets us use the comparator in sorted data structures like sets, priority queues etc.
		{
			return compare(lhs, rhs) >= 0;
		}
};


#endif
