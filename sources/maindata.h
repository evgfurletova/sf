#pragma once

#ifndef MAINDATA_H
#define MAINDATA_H

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <deque>
using namespace std;

#define		FileLineMax		10000
#define		AlpMax			100
#define		KeyWordMax		500000
#define		WordLenMax		2000
#define		MaxOccur		1000
#define		NMaxStates		2000

class MainData
{
public:
	// 1. Input files
	static std::string AlpFileName;
	static std::string PatternFileName;
	static std::string PSSMFileName;
	static std::string FootPrintFileName;
	static std::string ConsAlpFileName;
	static std::string OutName;

	// 2. Output files
	static ofstream ResLOG;						// output files
	static ofstream ResROG;
	static ofstream ResF;
	static ofstream ResWords;
	

	// 3. Input Parameters
	// 3.1. Alphabet and model of letters distribution 
	static int			 AlpSize;					// size of the alphabet
	static char		     AlpMas[AlpMax];			//alphabet
	static int		 	 order;						//Order of Markovian model (0 in the case of Bernoulli model)
	static double		 BernProb[AlpMax];			//Array of Bernulli probabilities for letters

	// 3.2 Parameters for the pattern
	//3.2.1 common parameters for all type of input parameters respect to the pattern
	static int	mode;						//type of input parameters
	static int NWords;						// number of words in pattern
	static int WordLen;						// length of pattern
	static int CrDistribFlag;
	//3.2.2. parameters for the pattern described by PSSM
	static double**	PssmMas;		// Pssm matrix
	static double Thr;						// Threhold for Pssm
	//3.2.3 parameters for the pattern described by a motif, number of replacements and number of constant positions
	static int* motif;							//motif
	static int Nreplace;					//number of replacements in the pattern
	static int*	ConstPositions;		//list of constant positions		
	//3.2.4 parameters of the pattern described by a consensus
	static int* consensus;			//consensus
	static vector< vector<char> > ConsAlp;	//alphabet for consensus
	//3.2.5 parameters for random pattern
	static double* RandPatProbs;			//probabilities to generate a random pattern	
	// 3.3 Information for calculating of probabilities
	static int TLen;						//length of random text 
	static int NOccur;						// Maximal number of patterns occurrences in the text
	
	// 4. Variables
	static int MAX;
	// 4.1. Variables for graphs creating
	
	// 4.2	Variables for probabilities computation
	static double	ProbRes;				//Sum_q(Prob_q(R(TLen,NOccur,HH)))
	static double	Pvalue;					//B(TLen,NOccur)
	// 4.2.1. for a Bernoulli model

	//Let Q be a set of states 
	// 4.2.2. for a Markov model
	static double**	MarkovProbs;	//Matrix ()of symbol probabilities for different states;
											// For Markov model the states are in alphabetical order;
	// 4.2.3. for a non-deterministic HHM
	static double*** ND_HHMProbs;			//Mairix (Q*Q*Alp->�(q,a,q')) of probabilities  �(q,a,q') being in a state q to generate a symbol a and traverse to a state q�.
	static vector<int>** ND_HHMTrans;		//Matrix (Q*Alp->Q)of sets of states such that �(q,a,q')>0
	// 4.2.4. for a deterministic HHM
	static double** D_HHMProbs;				//Mairix (Q*Alp->P(q,a)) of probabilities  P(q,a) being in a state q to generate a symbol a
	static int** D_HHMTrans;				//Mairix (Q*Alp->�(q,a)) state q=�(q,a) in which system passes from a state q generating a symbol a



public:

	static	int		AToi(char Let);	//gives number of a symbol Let in AlpMas
	static	char    IToa(int pos);	//gives symbol Let having number pos in AlpMas
	static  int		in_Alp(char *str);
	static  int	    GenRanWords(void);
	static  double  CountThr(char *word);
	static	void    SetScorMas(double* SMas);
	static  int     GenPssmWords(int *word, int i, double score, double* SMas);
	static  void    MotifVariations(int num, int pos, int* word);
	static  int     Pos_In_Cons_Alp(char Let);
	static  int     ConsVariations(int pos, int* word);
	static  int		GetInput(void);	//sets input parameters
	static	int		ComLineParse(int argc, char **argv); //parsing of the command line
	static	void	CrOutFiles(void); //prints information to the out files
	static	void    ErrorDetect(int Error);
};

//MainData Data;
extern std::string Words;
#endif
