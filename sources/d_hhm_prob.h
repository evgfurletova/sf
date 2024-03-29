#pragma once

#include <vector>
#include <string>

//Let P(q,a) be a probability function where q - state, a - symbol
//�(q,a) - transition function
class D_HHM_Prob {
public:
	static int AssociatedCondFlag(int q1,std::string word, int q2);	//1, if there exists path prom q1 to q2 with a mark word
																	// else - 0 
	static int AssociatedFlag(std::string word, int q);				//1,if there exists path from initial state to q with a mark word		
	static std::vector<int> CalcAssociatedStates(std::string word);
	static double TermCondProb(int q1,std::string word,int q2);		//calculates Prob_q1(word|q2)
	static double TermProb(std::string word, int q);				//calculates Prob_q(word)
	static double TransitionProb(int q1, int q2);					//calculates probability of transition from q1 to q2
	static std::vector<int> ConsistStates(int q);						//gets all states q' such that exisis a symbol a for that q'=�(q,a)
	static void Debug(void);
	D_HHM_Prob(void);
	~D_HHM_Prob(void);
};
