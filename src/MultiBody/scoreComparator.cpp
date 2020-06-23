
#include <MultiBody/scoreComparator.h>


ScoreComparator::ScoreComparator(int m, bool v)
{
	verbose = v;
	mode = m;
}

int ScoreComparator::compareUsingWeightedSum(const Score* S1, const Score* S2) const
{
	if(verbose) std::cout<<"Using weighted sum"<< std::endl;

	double totalscore1 = 0.0;
	double totalscore2 = 0.0;
	double s1, s2, w1, w2;

	const ScoreWeight* W1 = S1->getWeights();
	const ScoreWeight* W2 = S2->getWeights();

	int numScores = W1->getNumWeights();

	for(int i=0; i<numScores; i++)
	{
		if( ! W1->getWeight(i, w1) ) return 0;
		if( ! W2->getWeight(i, w2) ) return 0;
		if( ! S1->getScore(i, s1) ) return 0;
		if( ! S2->getScore(i, s2) ) return 0;

		if(verbose) std::cout << w1 << " " << s1 << " " << w2 << " " << s2 << "\n";

		totalscore1 += w1*s1;
		totalscore2 += w2*s2;
	}

	if(verbose) std::cout << totalscore1 << " " << totalscore2 << "\n ";

	if(totalscore1 > totalscore2) return 1;
	if(totalscore1 < totalscore2) return -1;

	return 0;
}


int ScoreComparator::compareUsingMajority(const Score* S1, const Score* S2) const
{
	if(verbose) std::cout<<"Using majority"<< std::endl;

	int vote1 = 0;
	int vote2 = 0;
	double s1, s2, w1, w2;

	const ScoreWeight* W1 = S1->getWeights();
	const ScoreWeight* W2 = S2->getWeights();

	int numScores = W1->getNumWeights();

	for(int i=0; i<numScores; i++)
	{
		if( ! W1->getWeight(i, w1) ) return 0;
		if( ! W2->getWeight(i, w2) ) return 0;
		if( ! S1->getScore(i, s1) ) return 0;
		if( ! S2->getScore(i, s2) ) return 0;

		if(verbose) std::cout << w1 << " " << s1 << " " << w2 << " " << s2 << "\n";

		if( w1*s1 > w2*s2) vote1++;
		else if( w1*s1 < w2*s2) vote2++;
	}

	if(verbose) std::cout << vote1 << " " << vote2 << "\n ";

	if(vote1 > vote2) return 1;
	if(vote1 < vote2) return -1;

	return 0;
}

int ScoreComparator::compareUsingLexicographicOrder(const Score* S1, const Score* S2) const
{
	if(verbose) std::cout<<"Using lexicography"<< std::endl;

	double s1, s2, w1, w2;

	const ScoreWeight* W1 = S1->getWeights();
	const ScoreWeight* W2 = S2->getWeights();

	int numScores = W1->getNumWeights();

	for(int i=0; i<numScores; i++)
	{
		if( ! W1->getWeight(i, w1) ) return 0;
		if( ! W2->getWeight(i, w2) ) return 0;
		if( ! S1->getScore(i, s1) ) return 0;
		if( ! S2->getScore(i, s2) ) return 0;

		if(verbose) std::cout << w1 << " " << s1 << " " << w2 << " " << s2 << "\n";

		if( w1*s1 > w2*s2) return 1;
		else if( w1*s1 < w2*s2) return -1;
	}

	return 0;
}

int ScoreComparator::compareUsingGroupWeightedSum(const Score* S1, const Score* S2) const 
{
	if(verbose) std::cout<<"Using grouped weighted sum"<< std::endl;

	return 0;
}

int ScoreComparator::compare(const Score* S1, const Score* S2) const
{
	if(verbose)
	{
		std::cout <<"Comparing "<<"\n"; 
		S1->print();
		std::cout << "with "<<"\n"; 
		S2->print();
		std::cout<<"Mode "<< mode<<std::endl;
	}

	int numScores = S1->getWeights()->getNumWeights();

	if(numScores != S2->getWeights()->getNumWeights())
		return 0;

	if(numScores != S1->getNumScores())
		return 0;

	if(numScores != S2->getNumScores())
		return 0;

	switch(mode)
	{
		case WEIGHTEDSUM:
			return compareUsingWeightedSum(S1, S2);
		case MAJORITY:
			return compareUsingMajority(S1, S2);
		case LEXICOGRAPHIC:
			return compareUsingLexicographicOrder(S1, S2);
		case GROUPWEIGHTEDSUM:
			return compareUsingGroupWeightedSum(S1, S2);
		default:
			return 0;
	}
}
