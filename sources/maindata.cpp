#include "maindata.h"
#include "ac_trie.h"
#include "ov_graf.h"
#include "bprob_graf.h"
#include "h_m_ovgraf.h"
#include "mmodel_prob.h"
#include "m_trtree.h"
#include "minov_graf.h"
#include "hhm_ovgraf.h"
#include "m_ovgraf.h"
#include <cstring>

std::string MainData::AlpFileName;
std::string MainData::PatternFileName;
std::string MainData::PSSMFileName;
std::string MainData::FootPrintFileName;
std::string MainData::ConsAlpFileName;

	// 2. Output files
ofstream MainData::ResLOG;						
ofstream MainData::ResROG;
ofstream MainData::ResF;
ofstream MainData::ResWords;
std::string MainData::OutName;

	// 3. Input Parameters
	// 3.1. Alphabet and model of letters distribution 
int MainData::AlpSize = 0;					
char    MainData::AlpMas[AlpMax];			
int	    MainData::order = 0;						
double	MainData::BernProb[AlpMax];			
	// 3.2 Parameters for pattern
		//3.2.1 common parameters for all type of input parameters respect to the pattern
int	MainData::mode = 0;	
int MainData::NWords = 0;						
int MainData::WordLen = 0;	
int MainData::CrDistribFlag = 0;
		//3.2.2. parameters for the pattern described by PSSM
double** MainData::PssmMas;		
double MainData::Thr = 0;						
		//3.2.3 parameters for the pattern described by a motif, number of replacements and number of constant positions
int* MainData::motif;				
int MainData::Nreplace = 0;					
int* MainData::ConstPositions = NULL;	
		//3.2.4 parameters of the pattern described by a consensus
int* MainData::consensus;			
vector< vector<char> > MainData::ConsAlp;
		//3.2.5 parameters for random pattern
double* MainData::RandPatProbs = NULL;

	// 3.3 Information for calculating of probabilities
int MainData::TLen	   = 0;						
int MainData::NOccur  = 0;					
	
	// 4. Variables
int MainData::MAX	   = 0;
						
	// 4.2	Variables for probabilities computation
double	MainData::ProbRes  = 0;
double	MainData::Pvalue   = 0;
// 4.2.1. for Markov model
double**	MainData::MarkovProbs = NULL;
// 4.2.3. for a non-deterministic HHM
double*** MainData::ND_HHMProbs = NULL; 
vector<int>** MainData::ND_HHMTrans = NULL;
// 4.2.4. for a deterministic HHM
double** MainData::D_HHMProbs = NULL;
int** MainData::D_HHMTrans = NULL;
std::string Words;	
//sets default parameters
void Default(void){
	MainData::order  = 0;
	MainData::AlpSize = 4;
	MainData::AlpMas[0] = 'A';
	MainData::AlpMas[1] = 'C';
	MainData::AlpMas[2] = 'G';
	MainData::AlpMas[3] = 'T';
	MainData::TLen   = 100;
	MainData::NOccur = 1;
	MainData::OutName = "Res.txt";
}

//checks is the line empty or not
int EmptyLine(char* line){
	size_t s = strlen(line);
	size_t i;
	for(i = 0; i < s; i++){
		if(line[i]== '!'){
			return 1;
		}
		if(line[i] != ' ')
			return 0;
	}
	return 1;
}

// finds first non empty line
void FNonEmptyLine(ifstream *ff,char* &line){
	while((ff->getline(line,FileLineMax))&& (EmptyLine(line) == 1)){		
	};
	if(EmptyLine(line) == 1){
        line[0] = '\0';
	}

	return;
}

//chacks is the line a figure or not
int Figure(char* line){
	int i;
	int s = (int)strlen(line);
	int ind = 0;
	int k = 0;
	if(line[0]=='-'){
			k++;
	}
	for(i = k; i < s; i++){
		int c = line[i];
		if((c < 48)||(c > 57)){
			if(line[i] == '.'){
				ind++;
				if(ind > 1)
					return 1;
			}
			else
				return 1;

		} 
	}
	return 0;
}


//Get alphabet
/////////For HHM/////////////////////
//reads information about HHM from the alphabet file
int Get_ND_HHM(ifstream *ff){
	int i,j,k;
	char* line = new char[FileLineMax];	
	char *pch;
	char delim[] =", ; \t";
	FNonEmptyLine(ff,line);
	H_M_OVGraf::NumAllStates = atoi(line);

	MainData::ND_HHMProbs = new double**[H_M_OVGraf::NumAllStates];
	MainData::ND_HHMTrans = new vector<int>*[H_M_OVGraf::NumAllStates];
	
	for(i = 0; i < H_M_OVGraf::NumAllStates; i++){
		MainData::ND_HHMProbs[i] = new double*[H_M_OVGraf::NumAllStates];
		MainData::ND_HHMTrans[i] = new vector<int>[MainData::AlpSize];
		for(j = 0; j < H_M_OVGraf::NumAllStates; j++){
			MainData::ND_HHMProbs[i][j] = new double[MainData::AlpSize];
			for(k = 0; k < MainData::AlpSize; k++){
				MainData::ND_HHMProbs[i][j][k] = 0;
			}
		}
	}	
	MainData::CrDistribFlag = 1;
	FNonEmptyLine(ff,line);
	FNonEmptyLine(ff,line);
	i = 0;
	while((strlen(line) != 0)&&(i < H_M_OVGraf::NumAllStates)){
		j = 0;
		while((strlen(line) != 0)&&(j < H_M_OVGraf::NumAllStates)){
			if(line[0] == '-'){
				delete[] line;
				return 7;
			}
			pch = strtok(line,delim);
			k = 0;
			while ((pch != NULL)&&(k < MainData::AlpSize))
			{
				if(Figure(pch) == 0){
					MainData::ND_HHMProbs[i][j][k] = atof(pch);
					if(MainData::ND_HHMProbs[i][j][k] != 0){
						MainData::ND_HHMTrans[i][k].push_back(j);
					}
				}
				else{
					delete[] line;
					return 4;
				}
				pch = strtok(NULL,delim);
				k++;
			}
		
			if(k != MainData::AlpSize){
				delete[] line;
				return 6;
			}
			FNonEmptyLine(ff,line);
			j++;
		}
		FNonEmptyLine(ff,line);
		i++;
	}
	if(i != H_M_OVGraf::NumAllStates){
		delete line;
		return 6;
	}
	delete[] line;
	return 0;
}

int Get_D_HHM(ifstream *ff){
	int i,j;
	char* line = new char[FileLineMax];	
	char *pch;
	char delim[] =", ; \t";
	FNonEmptyLine(ff,line);
	H_M_OVGraf::NumAllStates = atoi(line);

	MainData::D_HHMProbs = new double*[H_M_OVGraf::NumAllStates];
	MainData::D_HHMTrans = new int*[H_M_OVGraf::NumAllStates];
	
	for(i = 0; i < H_M_OVGraf::NumAllStates; i++){
		MainData::D_HHMProbs[i] = new double[MainData::AlpSize];
		MainData::D_HHMTrans[i] = new int[MainData::AlpSize];
		for(j = 0; j < MainData::AlpSize; j++){
			MainData::D_HHMProbs[i][j] = 0;
			MainData::D_HHMTrans[i][j] = 0;	
		}
	}	
	MainData::CrDistribFlag = 1;
	FNonEmptyLine(ff,line);
	FNonEmptyLine(ff,line);
	i = 0;
	while((line[0] != '-')&&(i < H_M_OVGraf::NumAllStates)){
		j = 0;
		pch = strtok(line,delim);
		while ((pch != NULL)&&(j < MainData::AlpSize))
		{
			if(Figure(pch) == 0){
				MainData::D_HHMProbs[i][j] = atof(pch);
			}
			else{
				delete[] line;
				return 4;
			}
			pch = strtok(NULL,delim);
			j++;
		}
		
		if(j != MainData::AlpSize){
			delete[] line;
			return 6;
		}
		FNonEmptyLine(ff,line);
		i++;
	}

	if(i != H_M_OVGraf::NumAllStates){
		delete[] line;
		return 6;
	}

	FNonEmptyLine(ff,line);
	i = 0;
	while((strlen(line) != 0)&&(i < H_M_OVGraf::NumAllStates)){
		j = 0;
		pch = strtok(line,delim);
		while ((pch != NULL)&&(j < MainData::AlpSize))
		{
			if(Figure(pch) == 0){
				MainData::D_HHMTrans[i][j] = atoi(pch);
			}
				else{
					delete[] line;
					return 4;
				}
				pch = strtok(NULL,delim);
				j++;
			}
		
			if(j != MainData::AlpSize){
				delete[] line;
				return 6;
			}
			FNonEmptyLine(ff,line);
			i++;
		}
	if(i != H_M_OVGraf::NumAllStates){
		delete[] line;
		return 8;
	}


	delete[] line;
	return 0;
}

//////////////////////////////////////
//sets parameters from alphabet file

int GetAlp(void){

	ifstream ff(MainData::AlpFileName.c_str());

	if(!ff){
		return 2;
	}	
	char* line = new char[FileLineMax];
	char* pch;
	char delim[] =", ; \t";
	FNonEmptyLine(&ff,line);
	pch = strtok(line,delim);
	if(Figure(pch) == 0){
		MainData::order = atoi(line);
	}
	else{
		delete[] line;
		return 3;
	}

	FNonEmptyLine(&ff,line);
	pch = strtok(line,delim);
	int i = 0;
	while (pch != NULL)
	{
		MainData::AlpMas[i]= pch[0];
		pch = strtok(NULL,delim);
		i++;
	}
	MainData::AlpSize = i;

	double norm = 0;
	if(MainData::order == 0){
		FNonEmptyLine(&ff,line);
		pch = strtok(line,delim);
		i = 0;
		while ((pch != NULL)&&(i < MainData::AlpSize))
		{
			if(Figure(pch) == 0){
				MainData::BernProb[i] = atof(pch);
				norm += MainData::BernProb[i];
			}
			else{
				delete[] line;
				return 4;
			}
			pch = strtok(NULL,delim);
			i++;
		}
		
		if((i != MainData::AlpSize)||(norm <= 0)){
			delete[] line;
			return 5;
		}
		for(i = 0; i < MainData::AlpSize; i++){
			MainData::BernProb[i] = MainData::BernProb[i]/norm;
		}
	}

	if(MainData::order > 0){
		FNonEmptyLine(&ff,line);
		int s = 1;
		int j;
		for(i = 0; i < MainData::order; i++){
			s = s*MainData::AlpSize;
		}
		MainData::MarkovProbs = new double*[MainData::AlpSize];
		for(i = 0; i < MainData::AlpSize; i++){
			MainData::MarkovProbs[i] = new double[s];
		}
		MainData::CrDistribFlag = 1;
		i = 0;
		while((strlen(line) != 0)&&(i < s)){
			pch = strtok(line,delim);
			j = 0;
			while ((pch != NULL)&&(j < MainData::AlpSize))
			{
				if(Figure(pch) == 0){
					MainData::MarkovProbs[j][i] = atof(pch);
					//norm += MainData::MarkovProbs[j][i];
				}
				else{
					delete[] line;
					return 4;
				}
				pch = strtok(NULL,delim);
				j++;
			}
		
			if(j != MainData::AlpSize){
				delete[] line;
				return 6;
			}
			FNonEmptyLine(&ff,line);
			i++;
		}
		if(i != s){
			delete[] line;
			return 6;
		}
		
	}
	if(MainData::order == -2){
		int Error = Get_ND_HHM(&ff);
		if(Error > 0){
			return Error;
		}
	}
	if(MainData::order == -1){
		int Error = Get_D_HHM(&ff);
		if(Error > 0){
			return Error;
		}
	}

	ff.close();
	delete[] line;
  return 0;

}


//gives number of a symbol Let in AlpMas
int MainData::AToi(char Let){

	int i;
	for(i = 0; i < AlpSize; i++){
		if(Let == AlpMas[i]){
			return i;
		}
	}
	return -1;
}

//gives symbol Let having number pos in AlpMas
char MainData::IToa(int pos){
	if(pos == -1){
		return '\0';
	}
	return AlpMas[pos];
}

//chacks is str a string under given alphabet
int MainData::in_Alp(char *str){
	size_t i;
	for(i = 0; i < strlen(str); i++){
		int k = MainData::AToi(str[i]);
		if((k < 0)|| (k > MainData::AlpSize)){
			return 0;
		}
	}
	return 1;
}

//It gets list of words and adds to trie

int GetKeyWords(void){

int i=0;
int k,j;
int Error = 0;
char *line = new char[WordLenMax];  
char *pch;
int* DigitLine;
char delim[] =", ; \t";
ifstream ff(MainData::PatternFileName.c_str());

	if (!ff){
		delete[] line;
		return 11;
	}
		
	FNonEmptyLine(&ff,line);
	if(strlen(line) !=  0){
		pch = strtok(line,delim);
		MainData::WordLen = (int)strlen(pch);
		DigitLine = new int[MainData::WordLen + 1];
		DigitLine[MainData::WordLen] = -1;
		for(j = 0; j< MainData::WordLen; j++){
			k = MainData::AToi(pch[j]);
			if(k != -1){
				DigitLine[j] = k;
			}
			else{
				delete[] line;
				delete[] DigitLine;
				return 12;
			}
		}
		i++;	
		AC_Trie::gTrie->InsertWord(DigitLine, i);
	}

	FNonEmptyLine(&ff,line);
	while (strlen(line) != 0){
		pch = strtok(line,delim);
		if(strlen(pch) == MainData::WordLen){	
			for(j = 0; j< MainData::WordLen; j++){
				k = MainData::AToi(pch[j]);
				if(k != -1){
					DigitLine[j] = k;
				}
				else{
					delete[] line;
					delete[] DigitLine;
					return 12;
				}
			}
			i++;
			AC_Trie::gTrie->InsertWord(DigitLine,i);
			FNonEmptyLine(&ff,line);
		}
		else{
			delete[] line;
			delete[] DigitLine;
			return 13;
		}
	}


	MainData::NWords = i;
	ff.close();
	delete[] line;
	delete[] DigitLine;
	return 0;

}


/////////////////////////////////////
//It generates a random word
void RanWord(int* &rword){ 
	
	int i;
	double r;
 
	for (i = 0; i< MainData::WordLen; i++){
	
		r = (double) rand ()/RAND_MAX;
		double p = 0;
		double p1 = MainData::RandPatProbs[0];
		int flag = 0;
		int k = 0;
		while( flag == 0){
			if((p <= r)&(r <= p1)){
				rword[i] = k;
				flag = 1;
			}else{
				p += MainData::RandPatProbs[k];
				p1 += MainData::RandPatProbs[k+1];
				k++;
			}
		}

	}
	return;
}




//It gives 1 if the word is in the trie else -0
 int In_Trie(int *word){
	
	 int k;
	 int i;
	 AC_Trie *T = AC_Trie::gTrie;

	 for(i = 0; i < MainData::WordLen; i++){
	 
		 k = word[i]; 
		if(T->Childs[k] == NULL)
			return 0;
		
		T = T->Childs[k];
	
	  }
	  
	  return 1;
  } 

//It generates list of random words and addes them to the trie
 int MainData::GenRanWords(void){

	int i = 0;
	int k;
	int *rword = new int[MainData::WordLen + 1];
	rword[MainData::WordLen] = -1;
	time_t time1 = time(NULL);
	int time = int(time1);
	srand(time);
	rand(); 
	
	while (i < MainData::NWords){

		RanWord(rword); 
		if(In_Trie(rword)==0){	
			i++;
			k = AC_Trie::gTrie->InsertWord(rword, i);
			if(k == 1){
				delete[] rword;
				return 1;
			}
		}
	}
	delete[] rword;
	return 0; 
}

/////////////////////////////////////////////////
//It gives 1 if the word requred pssm
int PssmWord(int *Word){

double Weight = 0;
int flag;
int i;
	
for (i=0; i < MainData::WordLen; i++){
		
		int j = Word[i];
		Weight = Weight + MainData::PssmMas[i][j];
	}
	
if(Weight < MainData::Thr)
		flag = 0;
	else
		flag = 1;
	
	return	flag;
}



//It gets pssm from the file
/*
int GetPssm(int* error){

	int i = 0;
	int j;
	float value = 0;

	ifstream ff(MainData::PSSMFileName.c_str());
	
	if (!ff){
		return -1;
		*(error) = 14;
	}
	
	char *line = new char[FileLineMax];
	FNonEmptyLine(&ff,line);
	vector<double*> PSSM;
	while (ff.eof() == 0 ){
		double*  row = new double[MainData::AlpSize]; 
		PSSM.push_back(row);
		for (j = 0; j < MainData::AlpSize; j++){
			ff>>value;
			PSSM[i][j]= value;
		}
		i++;
	}
	int len = i;
	MainData::PssmMas = new double*[len];
	for(i = 0; i < len; i++){
		MainData::PssmMas[i] = new double[MainData::AlpSize];
		for(j = 0; j < MainData::AlpSize; j++){
			MainData::PssmMas[i][j] = PSSM[i][j];
		}
		delete[] PSSM[i];
	}
	PSSM.clear();
	delete[] line;
	ff.close();
	return len;
}
*/

// It prints pssm
void PrintPssm(ofstream *ff){
	int i,j;

	*(ff)<<"\t   A\t\t  C\t\t  G\t\t  T\n";
	for( i = 0; i < MainData::WordLen; i++){
		*(ff)<<i+1<<'\t';	
		for(j = 0; j < MainData::AlpSize; j++){
			*(ff)<<setw(8)<<MainData::PssmMas[i][j]<<'\t';
		}
		*(ff)<<'\n';
	}
}


int GetPssm(int* error){

	int i = 0;
	int j;
	double value = 0; 
	char *line = new char[FileLineMax];
	char *pch;
	char delim[] =", ; \t";
	ifstream ff(MainData::PSSMFileName.c_str());
	
	if (!ff){
		error[0] = 14;
		return -1;
	}
	FNonEmptyLine(&ff,line);
	FNonEmptyLine(&ff,line);
	vector<double*> PSSM;
	while (strlen(line) != 0){
		double*  row = new double[MainData::AlpSize]; 
		PSSM.push_back(row);	
		pch = strtok(line,delim);
		if(Figure(pch) == 1){
			error[0] = 30;
			return -1;
		}
		value = atof(pch);
		PSSM[i][0] = value;
		for (j = 1; j < MainData::AlpSize; j++){
			pch = strtok(NULL,delim);
			if(Figure(pch) == 1){
				error[0] = 30;
				return -1;
			}
			value = atof(pch);
			PSSM[i][j]= value;
		}
		pch = NULL;
		i++;
		FNonEmptyLine(&ff,line);
	}
	int len = i;
	MainData::PssmMas = new double*[len];
	for(i = 0; i < len; i++){
		MainData::PssmMas[i] = new double[MainData::AlpSize];
		for(j = 0; j < MainData::AlpSize; j++){
			MainData::PssmMas[i][j] = PSSM[i][j];
		}
		delete[] PSSM[i];
	}
	PSSM.clear();
	delete[] line;
	ff.close();
	return len;
}





void MainData::SetScorMas(double* SMas){
	int i, j;
	for( i = 0; i<MainData::WordLen; i++){
	
		SMas[i] = -100;
		
		for( j = 0; j < MainData::AlpSize; j++){
			if(MainData::PssmMas[i][j] > SMas[i]){
				SMas[i] = MainData::PssmMas[i][j];
			}
		}
	}

	i = MainData::WordLen - 1;

	while(i > 0){
		i--;
		SMas[i] = SMas[i] + SMas[i+1]; 
	}

	return;
}



//It generates pssm words and adds them to the trie
/*
int MainData::GenPssmWords(int *word, int i, double score, double* SMas){

	double score1 = score;
	int j;
	for(j = 0; j < MainData::AlpSize; j++){

		score = score1;
		
		if(i == (MainData::WordLen - 1)){
			
			word[i] = j;
			
			if (PssmWord(word) == 1){
				MainData::NWords++;
			    AC_Trie::gTrie->InsertWord(word, MainData::NWords);
			}
		} 
		else{
		
			double ms = score + MainData::PssmMas[i][j]+SMas[i+1];
	
			if((ms > MainData::Thr)||(ms == MainData::Thr)){
				word[i] = j;
				score = score + MainData::PssmMas[i][j];
				MainData::GenPssmWords(word, i+1, score, SMas);
			}
		}
	}
	return 0;
}
*/

int MainData::GenPssmWords(int *word, int i, double score, double* SMas){

	double score1 = score;
	int j;
	for(j = 0; j < MainData::AlpSize; j++){

		score = score1;
		
		if(i == (MainData::WordLen - 1)){
			word[i] = j;
			score = score + MainData::PssmMas[i][j];
			if((score > MainData::Thr)||(score == MainData::Thr)){
				MainData::NWords++;
				AC_Trie::gTrie->InsertWord(word, MainData::NWords);
			}
		} 
		else{
		
			double ms = score + MainData::PssmMas[i][j]+SMas[i+1];
	
			if((ms > MainData::Thr)||(ms == MainData::Thr)){
				word[i] = j;
				score = score + MainData::PssmMas[i][j];
				MainData::GenPssmWords(word, i+1, score, SMas);
			}
		}
	}
	return 0;
}


//It calculates cut-off for a footprint

double MainData::CountThr(char *word){
	
	double t = -100;
	double Weight = 0;
	int i;
	int j,k;	
	size_t s1 = strlen(word);
	int s = (int)s1;
	if(s < MainData::WordLen)
		return -100;


	for(i = 0; i < s - MainData::WordLen + 1; i++){
	
		Weight = 0;
		for(k = 0; k < MainData::WordLen; k++){
			
			j = MainData::AToi(word[i+k]);
			Weight = Weight + MainData::PssmMas[k][j];
		}
		
		if(Weight > t)
			t = Weight;
	}
	return t;
}


//It count cut-off for a list of footprints
int GetThr(void){

int i=1;
double t;
char *line = new char[WordLenMax];
ifstream ff(MainData::FootPrintFileName.c_str());	

	if (!ff){
		delete[] line;
		return 15;
	}
	

	ff.getline(line,WordLenMax);

	if((line[0]!=0)&(line[0]!='!'))

		t = MainData::CountThr(line);
	
	MainData::Thr = t;

	while (ff.getline(line,WordLenMax))
		
		if((line[0]!=0)&(line[0]!='!')){
	
			t = MainData::CountThr(line);
			 if((t < MainData::Thr) & (t!= -100))
				 MainData::Thr = t;
			 i++;
		}
		
	delete[] line;
	ff.close();
	return 0;

}


///////////////////////////////
//for the given motif and given number of replacements gives possible variations of words

void MainData::MotifVariations(int num, int pos, int* word){
	int i;	
	if((num > 0)&&(pos < WordLen)){
		if((ConstPositions == NULL)||(ConstPositions[pos] == 0)){
			int Let = word[pos];
			for( i = 0; i < MainData::AlpSize; i++){
				if(Let != i){
					word[pos] = i;
					NWords++;
					AC_Trie::gTrie->InsertWord(word, MainData::NWords);
					MotifVariations(num - 1, pos + 1, word);
					word[pos] = Let;
				}else{
					MotifVariations(num, pos + 1, word);
				}
			}
		}
		else{
			MotifVariations(num, pos + 1, word);
		}
	}
	return;
} 
//////////////////////////////// 
//get alphabet of consensus
int GetConsensusAlp(void){
int i=0;
char *line = new char[WordLenMax];  
char *pch;
char delim[] =", ; \t { } =";
ifstream ff(MainData::ConsAlpFileName.c_str());

	for(i = 0; i < MainData::AlpSize; i++){
		vector<char> vec;
		vec.push_back(MainData::AlpMas[i]);
		vec.push_back(MainData::AlpMas[i]);
		MainData::ConsAlp.push_back(vec);
	}

	if (!ff){
		delete[] line;
		return 1;
	}
	FNonEmptyLine(&ff,line);
	char Let;
	if(strlen(line) !=  0){
		pch = strtok(line,delim);
		Let = pch[0];
		std::vector<char> vec;
		vec.push_back(Let);
		pch = strtok(NULL,delim);
		while(pch != NULL){
			if(MainData::in_Alp(pch) == 1){
				Let = pch[0];
				vec.push_back(Let);
			}
			else{
				delete[] line;
				return 17;
			}
			pch = strtok(NULL,delim);
		}
		MainData::ConsAlp.push_back(vec);
	}

	FNonEmptyLine(&ff,line);
	while (strlen(line) != 0){
		pch = strtok(line,delim);
		Let = pch[0];
		std::vector<char> vec;
		vec.push_back(Let);
		pch = strtok(NULL,delim);
		while(pch != NULL){
			if(MainData::in_Alp(pch) == 1){
				Let = pch[0];
				vec.push_back(Let);
			}
			else{
				delete[] line;
				return 17;
			}
			pch = strtok(NULL,delim);
		}
		MainData::ConsAlp.push_back(vec);
		FNonEmptyLine(&ff,line);
	}
	delete[] line;
	return 0;
};

//for  letter Let gives positions in the consensus alphabet 
int MainData::Pos_In_Cons_Alp(char Let){
	int i;
	int s = (int)MainData::ConsAlp.size();
	for(i = 0; i < s; i++){
		if(Let == MainData::ConsAlp[i][0]){
			return i;
		}
	}
	return -1;
}

//gives all words satisfying to the consensus
int MainData::ConsVariations(int pos, int* word){
	if(pos == MainData::WordLen){
		MainData::NWords++;
		AC_Trie::gTrie->InsertWord(word, MainData::NWords);
		return 0;
	}
	int pos1 = consensus[pos];
	int s = (int)MainData::ConsAlp[pos1].size();
	int i;
	for(i = 1; i < s; i++){
		char Let1 = MainData::ConsAlp[pos1][i];
		word[pos] = MainData::AToi(Let1);
		ConsVariations(pos + 1,word);
	}
	return 0;
}

////////////////////////////

int AlpFlag = 0;
std::string ForConsens;
//sets input parameters
int MainData::GetInput(void){
	int i;
	if(AlpFlag == 1){
		i = GetAlp();
		if(i > 0){
		 return i;
		}
	}
	AC_Trie::gTrie = new AC_Trie;
	if (mode ==0){
		i = GetKeyWords();
		if(i > 0){
			return i;
		}
	}

	if(mode == 1){
		if((WordLen > 0)&&(NWords > 0)){
			if(WordLen < 31){
				int pow = MModel_Prob::NumPower(AlpSize,WordLen);
				if(NWords > pow){
					return 29;
				}
			}

			GenRanWords();
			/*std::string FileName = OutName;
			int s = (int)FileName.length();
			if(FileName.at(s - 4) == '.'){
				FileName.erase(s-4,4);
			}
			FileName += "_Words.txt";
			ResWords.open(FileName.c_str());
			AC_Trie::gTrie->PrintTrie(&ResWords);
			*/
		}
	}
	
	if((mode == 2)||(mode == 3)){
		int* error = new int[1];
		i = GetPssm(error);
		if(i == -1){
			int er = error[0];
			delete[] error;
			return er;
		}
		delete[] error;
		WordLen = i;

		if(mode == 3){
			i = GetThr();
			if(i == -1)
				return 100;
		}
		
		int *word = new int[WordLen +1];
		word[WordLen] = -1;
		double* SMas = new double[WordLen + 1];
		SetScorMas(SMas);
		SMas[WordLen] = 0;
		GenPssmWords(word, 0, 0,SMas);
		delete[] word;
		for(i = 0; i < WordLen; i++){
			delete[] PssmMas[i];
		}
		delete[] PssmMas;
		PssmMas = NULL;
		delete[] SMas;
		
	}
	if(mode == 4){
		NWords ++;
		AC_Trie::gTrie->InsertWord(motif, NWords);
		int* word = new int[MainData::WordLen + 1];
		word[MainData::WordLen]= -1;
		for(i = 0; i < WordLen; i++){
			word[i] = motif[i];
		}
		MotifVariations(Nreplace, 0, word);
		delete[] MainData::motif;
		delete[] word;
		delete[] MainData::ConstPositions;
	}
	if(MainData::mode == 5){
		i = GetConsensusAlp();
		if(i > 0){
			return i;
		}
		consensus = new int[WordLen + 1];
		consensus[WordLen] = -1;
		int j;
		for(j = 0; j < WordLen; j++){
			int Let = Pos_In_Cons_Alp(ForConsens.at(j));
			if( Let != -1){
				consensus[j] = Let;
			}
			else{
				delete[] consensus;
				return 18;
			}
		}
				
		int* word = new int[WordLen + 1];
		word[WordLen] = -1;
		i = ConsVariations(0, word);
		if(i > 0){
			return i;
		}
		delete[] word;
		delete[] MainData::consensus;
	}

	if(NWords == 0){
		return 23;
	}

	if(WordLen == 0){
		return 24;
	}
	return 0;
}

///////////////////////////////////////////////////////
int NParam = 1;

//Parsing of the command line
int ParseLine(int i, char **argv, int argc){

	if(argv[i][1]=='a'){
		if((NParam <= argc - 2) && (argv[i+1][0]!='-')){
			NParam += 2;
			if((NParam <= argc - 3)&&(argv[i+2][0] != '-')){
				return 1;
			}
			MainData::AlpFileName = argv[i+1];
			AlpFlag = 1;	
		}
		else{
			return 1;
		}
	}

	if(argv[i][1]=='m'){
		if(NParam <= argc - 3){
			if(atoi(argv[i+1])==0){
				if(argv[i+2][0]!='-'){
					if((NParam <= argc - 4)&&(argv[i+3][0] != '-')){
						return 9;
					}
					NParam += 3;
					MainData::mode = 0;
					MainData::PatternFileName = argv[i+2];
				} 
				else{
					return 9;
				}
			}
			if((atoi(argv[i+1])==1)&&(NParam <= argc - 4)){
				if((argv[i+2][0]!='-')&(argv[i+3][0]!='-')){
					NParam += 4;
					MainData::mode = 1;
					if((Figure(argv[i+2]) == 1)||(Figure(argv[i+3]) == 1)){
						return 10;
					}
					MainData::WordLen = atoi(argv[i+2]);
					MainData::NWords = atoi(argv[i+3]);
					int j;
					MainData::RandPatProbs = new double[MainData::AlpSize];
					if((NParam <= argc - 5)&&(argv[i + 4][0] != '-')){
						for(j = 0; j < MainData::AlpSize; j++){
							if((Figure(argv[i+j+4])==0)&&(NParam < argc)){
								MainData::RandPatProbs[j]= (double)atof(argv[i+j+4]);
								NParam++;
							}
							else{
								return 10;
							}
						}
					}
					else{
						for(j = 0; j < MainData::AlpSize; j++){
							MainData::RandPatProbs[j] = (double)1/MainData::AlpSize;
						}
					}
				}
				else{
					return 9;
				}
			}
			if((atoi(argv[i+1])==2)&&(NParam <= argc - 4)){
				if((argv[i+2][0]!='-')&&(Figure(argv[i+3])==0)){	
					if((NParam <= argc - 5)&&(argv[i+4][0] != '-')){
						return 9;
					}
					NParam += 4; 			
					MainData::mode = 2;
					MainData::PSSMFileName = argv[i+2];
					MainData::Thr = atof(argv[i+3]);
				} 
				else{
					return 9;
				}
			}
			if((atoi(argv[i+1])== 3)&&(NParam <= argc - 4)){
				if((argv[i+2][0]!='-')&(argv[i+3][0]!='-')){
					if((NParam <= argc - 5)&&(argv[i+4][0] != '-')){
						return 9;
					}
					NParam += 4;
					MainData::mode = 3;
					MainData::PSSMFileName = argv[i+2];
					MainData::FootPrintFileName = argv[i+3];
				} else{
					return 9;
				}
			}
			if((atoi(argv[i+1])==4)&&(NParam <= argc - 4)){
				int j;
				MainData::mode = 4;
				if((argv[i+2][0] != '-')&&(argv[i+3][0] != '-')){
					NParam += 4;
					MainData::WordLen = (int)strlen(argv[i+2]);
					MainData::motif = new int[MainData::WordLen + 1];
					MainData::motif[MainData::WordLen] = -1;
					for(j = 0; j < MainData::WordLen; j++){
						int Let = MainData::AToi(argv[i+2][j]);
						if( Let != -1){
							MainData::motif[j] = Let;
						}
						else{
							delete[] MainData::motif;
							return 12;
						}
					}
					if(Figure(argv[i+3]) == 1){
						return 10;
					}
					MainData::Nreplace = atoi(argv[i+3]);
				}
				else{
					return 9;
				}
				int l = 4;
				int pos;
				if(argv[i+4][0] != '-'){
					MainData::ConstPositions = new int[MainData::WordLen];
					for(j = 0; j < MainData::WordLen; j++){
					 	MainData::ConstPositions[j] = 0;
					}
				}
				
				while((NParam < argc)&&(argv[i+l][0] != '-')){
					if(Figure(argv[i+1]) == 1){
						return 10;
					}
					pos = atoi(argv[i+l]);
					if((pos >= 0)&&(pos < MainData::WordLen)){
						MainData::ConstPositions[pos] = 1;
						NParam++;
						l++;
					}
					else{
						delete[] MainData::ConstPositions;
						delete[] MainData::motif;
						return 10;
					}
				}
			}
			if((atoi(argv[i+1])== 5)&&(NParam <= argc - 4)){
				if((argv[i+2][0]!='-')&(argv[i+3][0]!='-')){
					if((NParam <= argc - 5)&&(argv[i+4][0] != '-')){
						return 9;
					}
					NParam += 4;
					MainData::mode = 5;
					MainData::WordLen = (int)strlen(argv[i+2]);
					ForConsens = argv[i+2];
					MainData::ConsAlpFileName = argv[i+3];
				} else{
					return 9;
				}
			}

		}else{
			return 9;
		}
	}

	if(argv[i][1]=='o'){
		if(NParam <= argc - 2){
			NParam += 2;
			MainData::OutName = argv[i+1];
		}else{
			return 21;
		}
	}

	if(argv[i][1]== 'p'){
		if(NParam <= argc - 3){
			if((argv[i+1][0]!='-')&&(argv[i+2][0]!='-')){
				if((NParam <= argc - 4)&&(argv[i+3][0] != '-')){
					return 19;
				}
				NParam += 3;
				if((Figure(argv[i+1]) == 1)||(Figure(argv[i+2]) == 1)){
					return 20;
				}
				MainData::NOccur = atoi(argv[i+1]);
				MainData::TLen = atoi(argv[i+2]);
			} 
			else{
				return 20;
			}
		} else{
			return 19;
		}
	}

	if(argv[i][1] == 'b'){
		if(NParam <= argc - 5){
			if((argv[i+1][0]!='-')&(argv[i+2][0]!='-')&(argv[i+3][0]!='-')&(argv[i+4][0]!='-')){
				
				NParam += 5;
				double sum = atof(argv[i+1]) + atof(argv[i+2]) + atof(argv[i+3]) + atof(argv[i+4]);

				MainData::BernProb[0] = atof(argv[i+1])/sum;
				MainData::BernProb[1] = atof(argv[i+2])/sum;
				MainData::BernProb[2] = atof(argv[i+3])/sum;
				MainData::BernProb[3] = atof(argv[i+4])/sum;
			}
			else{
				return 22;
			} 
		}
		else{
			return 22;
		}
	}	
	return 0;
}

int ParseParam(void){

return 0;
}


int MainData::ComLineParse(int argc, char **argv){
	int i;
	int flag1 = 0;
	int flag2 = 0;
	int code;

	Default();
	
	for(i = 1; i< argc; i++){
		if(argv[i][0] == '-'){
			code = ParseLine(i, argv, argc);
			if(code > 0)
				return code;
			if(argv[i][1]=='m')
				flag1 = 1;
			if(argv[i][1]=='p')
				flag2 = 1;
		}
	}

	if(flag1==0){
		return 9;
	}
	if(flag2==0){
		return 19;
	}
	
	if(ParseParam() == 1)
		return 1;
	return 0;
}
/////////////////////////////////////////////////////


////////////////////////////////////////////////

//It prints main data


void PrintMain(ofstream *ff){
		

	*(ff)<<"\tProgram SuffixPrefix, v.1.2.\n";
	*(ff)<<"1. Input prameters\n";
	*(ff)<<"1.1. Alphabet and probabilities distribution of letters:\n";
	int i;
	*(ff)<<"Size of the alphabet: "<<MainData::AlpSize<<'\n';
	*(ff)<<"Alphabet: "<<'\n';
	for(i = 0; i < MainData::AlpSize; i++){
		*(ff)<<MainData::AlpMas[i]<<'\t';
	}
	if(MainData::order == 0){
		*(ff)<<"\n\n Bernoulli distribution of letters:\n";

		for(i = 0; i < MainData::AlpSize; i++){
			*(ff)<<MainData::BernProb[i]<<'\t';
		}
		*(ff)<<"\n\n";
	}
	if(MainData::order > 0){
		*(ff)<<"\n\nOrder of Markovian model: "<<MainData::order<<"\n";
		*(ff)<<"Distribution of letters (in the prefix order) \n";
		int s = 1;
		int j;
		for(i = 0; i <MainData:: order; i++){
			s = s*MainData::AlpSize;
		}
		for(i = 0; i < s; i++){
			for(j = 0; j < MainData::AlpSize; j++){
				*(ff)<<MainData::MarkovProbs[j][i]<<'\t';
			}
			*(ff)<<'\n';
		}
	}
	if(MainData::order == -1){
		*(ff)<<"\n\n Probability model: Deterministic HHM\n";
	}
	if(MainData::order == -2){
		*(ff)<<"\n\n Probability model: HHM\n";
	}

	*(ff)<<"\n1.2. Pattern description: \n";

	if(MainData::mode==0){
	  *(ff)<<"Input mode: List of words\n";
	}
	if(MainData::mode==1){
	
      *(ff)<<"Input mode: Random words\n";
	}
	if(MainData::mode>1){
	
      *(ff)<<"Input mode: Given PSSM\n";
	  *(ff)<<"Cut-off: "<<MainData::Thr;	 
	  *(ff)<<"\nPSSM:  \n";
	  //PrintPssm(ff);
	  *(ff)<<"\n\n";
	}
	
	*(ff)<<"Word length: "<<MainData::WordLen<<'\n';
	*(ff)<<"Number of words: "<<MainData::NWords<<"\n\n";

	*(ff)<<"1.3. Parameters for probabilities calculation: \n";
	*(ff)<<"Size of the text: "<<MainData::TLen<<'\n';
	*(ff)<<"Number of occurences: "<<MainData::NOccur<<"\n\n";

	*(ff)<<"2. Results \n";
	*(ff)<<"2.1. Number of overlap classes: "<<OV_Graf::NClasses<<"\n\n";
	*(ff)<<"2.2. Number of nodes in graphs: "<<OV_Graf::NumOVNodes - MainData::NWords + OV_Graf::NClasses<<'\n';
	if(OV_Graf::NumOVNodes != 1){
		*(ff)<<"Number of  Internal nodes in LOG/ROG: "<<OV_Graf::NumOVNodes - MainData::NWords<<'\n';
	}else{
		*(ff)<<"Number of  Internal nodes in LOG/ROG: 0 \n";
	}
	*(ff)<<"Number of  nodes in AC trie: "<<AC_Trie::NumACNodes<<"\n\n";
	//ff->setf(ios::fixed,ios::floatfield); 
	ff->precision(15);
	*(ff)<<"2.3. Probabilities:\n";
	*(ff)<<"Probability r("<<MainData::TLen<<','<<MainData::NOccur<<",H):"<<MainData::ProbRes<<'\n';
	*(ff)<<"pValue: "<<MainData::Pvalue;
	//PrintProbs();
	return;

}     

//prints information to the out files
void MainData::CrOutFiles(void){
//	int i;
	size_t len1 = OutName.length();
	int len = (int)len1; 
	std::string str;
	str = OutName;

	str +=  "_res.txt";
	ResF.open(str.c_str());
	PrintMain(&ResF);

	str = OutName;
	str += "LOG.xml";
	ResLOG.open(str.c_str());
	
	str = OutName;
	str += "ROG.xml";
	ResROG.open(str.c_str());

	if(MainData::order == 0){
		BProb_Graf::gBOVG->PrintLOG(0, &ResLOG);	
		BProb_Graf::gBOVG->PrintROG(0, &ResROG);
	}
	str = OutName;
	str += "Words.txt";
	ResWords.open(str.c_str());
	//AC_Trie::gTrie->PrintTrie(&ResWords);	
	ResWords.close();

	ResF.close();
	ResLOG.close();
	ResROG.close();

}

///////////////////////Detector of Errors/////////////////////////////////////////
void MainData::ErrorDetect(int Error){
	if(Error == 1){
		cerr<<"Error1: Wrong number of parameters associated with the alphabet"<< '\n';
		return;
	}	
	if(Error == 2){
		cerr<<"Error2: Error of opening of file with alphabet description"<<'\n';
		return;
	}	
	if(Error == 3){
		cerr<<"Error3: Invalid value of model order. This value must be integer"<<'\n';
		return;
	}
	if(Error == 4){
		cerr<<"	Error4: Incorrect distribution of letters. Frequences of letters must be real numbers"<<'\n';
		return;
	}
	if(Error == 5){
		cerr<<"Error5: Invalid distribution. Number of frequences must be equal to the size of alphabet"<<'\n';
		return;
	}
	if(Error == 6){
		cerr<<"Error6: Invalid distribution. Incorrect size of matrix with frequences of letters"<< '\n';
		return;
	}
	if(Error == 7){
		cerr<<"Error7: Incorrect alphabet distribution"<< '\n';
		return;
	}
	if(Error == 8){
		cerr<<"Error8: Invalid distribution. Incorrect size of matrix with transition states"<< '\n';
		return;
	}
	if(Error == 9){
		cerr<<"Error9: Wrong number of parameters associated with the pattern"<<'\n';
		return;
	}
	if(Error == 10){
		cerr<<"Error10: Incorrect parameters associated with the pattern"<<'\n';
		return;
	}
	if(Error == 11){
		cerr<<"Error11: Error of opening of file with pattern description"<<'\n';
		return;
	}
	if(Error == 12){
		cerr<<"Error12: Incorrect words in the pattern. Letters of words must be in the alphabet"<<'\n';
		return;
	}
	if(Error == 13){
		cerr<<"Error13: Incorrect length of the words in the pattern"<<'\n';
		return;
	}
	if(Error == 14){
		cerr<<"Error14: Error of opening of file with PSSM"<<'\n';
		return;
	}
	if(Error == 15){
		cerr<<"Error15: Error of opening of file with footprints"<<'\n';
		return;
	}
	if(Error == 16){
		cerr<<"Error16: Error of opening of file with consensus alphabet description"<<'\n';
		return;
	}
	if(Error == 17){
		cerr<<"Error17: Incorrect description of consensus alphabet"<<'\n';
		return;
	}
	if(Error == 18){
		cerr<<"Error18: Incorrect letters of the consensus"<<'\n';
		return;
	}
	if(Error == 19){
		cerr<<"Error19: Wrong number of parameters associated with text length or number of occurences"<<'\n';
		return;
	}
	if(Error == 20){
		cerr<<"Error20: Incorrect parameters associated with text length or number of occurences"<<'\n';
		return;
	}
	if(Error == 21){
		cerr<<"Error21: Wrong number of parameters associated with the name of output file"<<'\n';
		return;
	}
	if(Error == 22){
		cerr<<"Error22: Wrong number of probabilities associated with Bernoulli probabilities"<<'\n';
		return;
	}
	if(Error == 23){
		cerr<<"Error23: Incorrect number of words in the pattern"<<'\n';
		return;
	}
	if(Error == 24){
		cerr<<"Error24: Incorrect length of the pattern"<<'\n';
		return;
	}
	if(Error == 25){
		cerr<<"Error25: Incorrect size of the alphabet"<<'\n';
		return;
	}
	if(Error == 26){
		cerr<<"Error26: Incorrect order of the probability model"<<'\n';
		return;
	}
	if(Error == 27){
		cerr<<"Error27: Incorrect number of states of the HHM"<<'\n';
		return;
	}
	if(Error == 28){
		cerr<<"Error28: Incorrect size of the alphabet"<<'\n';
		return;
	}
	if(Error == 29){
		cerr<<"Error29: The number of words must be smaller then  alphabet size in power of the pattern length"<<'\n';
		return;
	}
	if(Error == 30){
		cerr<<"Error30: Incorrect PSSM description"<<'\n';
	}
	return;
}

