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


//////////////the function for the server////////////////////////
int func_set_input_data( int order, int mode, int TLen, int NOccur)
{

	if(order < -2){
		return 3;
	}
	
	if((mode < 0)||(mode > 5)){
		return 10;
	}
	if((TLen < 0)||(NOccur < 0)){
		return 8;
	}
	MainData::order = order;
	MainData::mode = mode;
	MainData::TLen = TLen;
	MainData::NOccur = NOccur;
	return 0;
}


int func_analis_alp_bern_data(int AlpSize, char* AlpMas, int order, double* BernProb)
{
	int i;

	if(AlpSize <= 0){
		return 25;
	}
	
	MainData::AlpSize = AlpSize;
	for(i = 0; i < AlpSize; i++){
		MainData::AlpMas[i] = AlpMas[i];
	}
	
	double sum = 0;
	for(i = 0; i < AlpSize; i++){
		MainData::BernProb[i] = BernProb[i];
		sum = sum + BernProb[i];
	}
	if(sum != 1){
		return 4;
	}
	return 0;
}

int func_analis_alp_mark_data(int AlpSize, char* AlpMas, double**	MarkovProbs)
{
	int i,j;

	if(AlpSize <= 0){
		return 25;
	}
	MainData::AlpSize = AlpSize;
	for(i = 0; i < AlpSize; i++){
		MainData::AlpMas[i] = AlpMas[i];
	}

	int power = MModel_Prob::NumPower(AlpSize,MainData::order);
	MainData::MarkovProbs = new double*[MainData::AlpSize];
	for(i = 0; i < AlpSize; i++){
		MainData::MarkovProbs[i] = new double[power];
	}
	for(i = 0; i < AlpSize; i++){
		for(j = 0; j < power; j++){
			MainData::MarkovProbs[i][j] = MarkovProbs[i][j];
		}
	}
	MainData::CrDistribFlag = 1;
	return 0;
}

int func_analis_alp_dhhm_data(int AlpSize, char* AlpMas, int NumAllStates,
							  double** D_HHMProbs, int** D_HHMTrans)
{
	int i,j;

	if(AlpSize <= 0){
		return 25;
	}
	MainData::AlpSize = AlpSize;
	for(i = 0; i < AlpSize; i++){
		MainData::AlpMas[i] = AlpMas[i];
	}
	
	if(NumAllStates < 0){
		return 27;
	}
	H_M_OVGraf::NumAllStates = NumAllStates;
	MainData::D_HHMProbs = new double*[NumAllStates];
	MainData::D_HHMTrans = new int*[NumAllStates];
	
	for(i = 0; i < NumAllStates; i++){
		MainData::D_HHMProbs[i] = new double[AlpSize];
		MainData::D_HHMTrans[i] = new int[AlpSize];
		for(j = 0; j < AlpSize; j++){
			MainData::D_HHMProbs[i][j] = D_HHMProbs[i][j];
			MainData::D_HHMTrans[i][j] = D_HHMTrans[i][j];	
		}
	}	
	MainData::CrDistribFlag = 1;

	return 0;
}


int func_analis_alp_hhm_data(int AlpSize, char* AlpMas,
							 int NumAllStates, double*** ND_HHMProbs)
{
	int i,j;

	if(AlpSize <= 0){
		return 25;
	}
	MainData::AlpSize = AlpSize;
	for(i = 0; i < AlpSize; i++){
		MainData::AlpMas[i] = AlpMas[i];
	}
	
	if(NumAllStates < 0){
		return 27;
	}
	H_M_OVGraf::NumAllStates = NumAllStates;
	MainData::ND_HHMProbs = new double**[NumAllStates];
	MainData::ND_HHMTrans = new vector<int>*[NumAllStates];
	
	for(i = 0; i < NumAllStates; i++){
		MainData::ND_HHMProbs[i] = new double*[NumAllStates];
			MainData::ND_HHMTrans[i] = new vector<int>[AlpSize];			
			for(j = 0; j < NumAllStates; j++){
				MainData::ND_HHMProbs[i][j] = new double[AlpSize];
				int k;
				for(k = 0; k < AlpSize; k++){
					MainData::ND_HHMProbs[i][j][k] = ND_HHMProbs[i][j][k];
					if(ND_HHMProbs[i][j][k] != 0){
						MainData::ND_HHMTrans[i][k].push_back(j);
					}
				}
			}
		}
	MainData::CrDistribFlag = 1;

	return 0;
}

int func_analis_pattern_data_0(int NWords, char **WordsList)
{
	int i,j,k;
	AC_Trie::gTrie = new AC_Trie;
	if(MainData::mode == 0){
		MainData::NWords  = NWords;
		if(MainData::NWords == 0)
			return 23;
		MainData::WordLen = (int)strlen(WordsList[0]);
		int* DigitLine = new int[MainData::WordLen + 1];
		DigitLine[MainData::WordLen] = -1;
	
		for(i = 1; i <= MainData::NWords; i++){
			if(strlen(WordsList[i]) == MainData::WordLen){					
				for(j = 0; j< MainData::WordLen; j++){
					k = MainData::AToi(WordsList[i][j]);
					if(k != -1){
						DigitLine[j] = k;
					}
					else{
						delete[] DigitLine;
						return 12;
					}
				}
				AC_Trie::gTrie->InsertWord(DigitLine,i);
			}
			else{
				delete[] DigitLine;
				return 13;
			}
		}
		delete[] DigitLine;
	}
	return 0;
}

int func_analis_pattern_data_1(int NWords, int WordLen, double* RandPatProbs)
{
	int i;

	AC_Trie::gTrie = new AC_Trie;
	
	if(MainData::mode == 1){
		if(WordLen <= 0){
			return 24;
		}
		if(NWords <= 0){
			return 23;
		}
		MainData::WordLen = WordLen;
		MainData::NWords = NWords;
		if(RandPatProbs == NULL){
			for(i = 0; i < MainData::AlpSize; i++){
				MainData::RandPatProbs[i] = (double)1/MainData::AlpSize;
			}
		}
		else{
			double sum = 0;
			for(i = 0; i < MainData::AlpSize; i++){	
				MainData::RandPatProbs[i] = RandPatProbs[i];
				sum = sum + RandPatProbs[i];
			}
			if(sum != 1){
				return 10;
			}
		}

		MainData::GenRanWords();
	}
	return 0;
}

int func_analis_pattern_data_2_3(int WordLen, int NFootPrints, char **FootPrints, double** PssmMas, double Thr)
{
	int i;
	AC_Trie::gTrie = new AC_Trie;
	
	if((MainData::mode == 2)||(MainData::mode == 3)){
		if(WordLen <= 0){
			return 24;
		}
		MainData::WordLen = WordLen;
		MainData::PssmMas = new double*[MainData::WordLen];
		for(i = 0; i < MainData::WordLen; i++){
			MainData::PssmMas[i] = new double[MainData::AlpSize];
			int j;
			for(j = 0; j < MainData::AlpSize; j++){
				MainData::PssmMas[i][j] = PssmMas[i][j];
			}
		}
		if(MainData::mode == 2){
			MainData::Thr = Thr;
		}
		else{
			double t;
			MainData::Thr = 100;
			for(i = 0; i < NFootPrints; i++){
				t = MainData::CountThr(FootPrints[i]);
				 if((t < MainData::Thr) & (t!= -100))
					 MainData::Thr = t;
			}
		}
		int *word = new int[WordLen +1];
		word[WordLen] = -1;
		double* SMas = new double[MainData::WordLen];
		MainData::SetScorMas(SMas);
		MainData::GenPssmWords(word, 0, 0, SMas);
		delete[] word;
		for(i = 0; i < WordLen; i++){
			delete[] MainData::PssmMas[i];
		}
		delete[] MainData::PssmMas;
		PssmMas = NULL;
		delete[] SMas;
		SMas = NULL;
	}
	return 0;
}

int func_analis_pattern_data_4(char* motif, int Nreplace, int NConstPositions, 
									  int *ConstPositions)
{
	int i;
	
	AC_Trie::gTrie = new AC_Trie;
	
	if(MainData::mode == 4){
		if((strlen(motif) <= 0)||(Nreplace < 0)){
			return 10;
		}
		MainData::WordLen = (int)strlen(motif);
		MainData::motif = new int[MainData::WordLen + 1];
		MainData::motif[MainData::WordLen] = -1;
		for(i = 0; i < MainData::WordLen; i++){
			int Let = MainData::AToi(motif[i]);
			if( Let != -1){
				MainData::motif[i] = Let;
			}
			else{
				delete[] MainData::motif;
				return 12;
			}
		}
		
		MainData::Nreplace = Nreplace;
		MainData::NWords ++;
		AC_Trie::gTrie->InsertWord(MainData::motif, MainData::NWords);

		int s = NConstPositions;
		if(s > 0){
			MainData::ConstPositions = new int[MainData::WordLen];
			for(i = 0; i < MainData::WordLen; i++){
				MainData::ConstPositions[i] = 0;
			}
			for(i = 0; i < s; i++){
				int pos = ConstPositions[i];
				if((0 <= pos)&&(pos < MainData::WordLen)){
					MainData::ConstPositions[pos] = 1;
				}
				else{
					delete[] MainData::motif;
					delete[] MainData::ConstPositions;
					return 10;
				}
			}
		}
		int* word = new int[MainData::WordLen + 1];
		word[MainData::WordLen]= -1;
		for(i = 0; i < MainData::WordLen; i++){
			word[i] = motif[i];
		}
		MainData::MotifVariations(Nreplace, 0, word);
		delete[] MainData::motif;
		delete[] word;
		delete[] MainData::ConstPositions;
	}
	return 0;
}

int func_analis_pattern_data_5(char *consensus, int NSymbols, char **ConsAlp)
{
	int i;
	int Error;
	AC_Trie::gTrie = new AC_Trie;
	
	if(MainData::mode == 5){
		for(i = 0; i < MainData::AlpSize; i++){
			vector<char> vec;
			vec.push_back(MainData::AlpMas[i]);
			vec.push_back(MainData::AlpMas[i]);
			MainData::ConsAlp.push_back(vec);
		}

		int s,j;
		for(i = 0; i < NSymbols; i++){
			s = (int)strlen(ConsAlp[i]);
			vector<char> vec;
			for(j = 0; j < s; j++){
				vec.push_back(ConsAlp[i][j]);
			}
			MainData::ConsAlp.push_back(vec);
		}	
		MainData::WordLen = (int)strlen(consensus);
		MainData::consensus = new int[MainData::WordLen + 1];
		MainData::consensus[MainData::WordLen] = -1;
		for(j = 0; j < MainData::WordLen; j++){
			int Let = MainData::Pos_In_Cons_Alp(consensus[j]);
			if( Let != -1){
				MainData::consensus[j] = Let;
			}
			else{
				delete[] MainData::consensus;
				return 18;
			}
		}
		int* word = new int[MainData::WordLen + 1];
		word[MainData::WordLen] = '\0';
		Error = MainData::ConsVariations(0, word);
		if(Error > 0){
			return Error;
		}
		delete[] word;
		delete[] MainData::consensus;
	}
	return 0;
}


int func_main(double* pvalue){

	int i;
	int ExitFlag = 0;
	int Error = 0;

	if((MainData::NOccur == 0)&&(ExitFlag == 0)){
		*(pvalue) = 1;
		ExitFlag = 1;
	}
	
	if((MainData::TLen < MainData::WordLen)&&(ExitFlag == 0)){
		*(pvalue) = 1;
		ExitFlag = 1;
	}
	if((MainData::TLen == MainData::WordLen)&&(MainData::NOccur > 1)&&(ExitFlag == 0)){
		*(pvalue) = 0;
		ExitFlag = 1;
	}

	if((MainData::AlpSize == 1)&&(ExitFlag == 0)){
		if(MainData::TLen >= MainData::WordLen + MainData::NOccur - 1){
			*(pvalue) = 1;
			ExitFlag = 1;
		}
		else{
			*(pvalue) = 0;
			ExitFlag = 1;
		}
	}

    if((MainData::order > 0)&&(ExitFlag == 0)){
		if(MainData::order > MainData::WordLen){
			ExitFlag =  1;
		}
        M_TrTree::ForLeaves = new int[MainData::NWords];
    }

	if(ExitFlag == 0){
		AC_Trie::CreateTrie();
		OV_Graf::gOVG = new OV_Graf(0);
		//OV_Graf::CreateGraf();
		MinOV_Graf::gMinOVG = new MinOV_Graf();
		MinOV_Graf::gMinOVG->Copy(1,OV_Graf::gOVG);
		delete OV_Graf::gOVG;
		MinOV_Graf::MinimizeGraf();
    
	    if(MainData::order < 0){
			HHM_OVGraf::gHHM = new HHM_OVGraf();
	        HHM_OVGraf::gHHM->Copy(4,MinOV_Graf::gMinOVG);  
	        HHM_OVGraf::gHHM->Preprocessing();
	        delete AC_Trie::gTrie;
			delete MinOV_Graf::gMinOVG;
			 HHM_OVGraf::gHHM->ProbCalc();
			delete HHM_OVGraf::gHHM;
		}   
		if(MainData::order == 0){
			BProb_Graf::gBOVG = new BProb_Graf();
			BProb_Graf::gBOVG->Copy(2,MinOV_Graf::gMinOVG);
			BProb_Graf::PredStep();
			delete AC_Trie::gTrie;
			delete MinOV_Graf::gMinOVG;
			BProb_Graf::ProbCalc();
			delete BProb_Graf::gBOVG;
			delete[] BProb_Graf::PHCl;
			delete[] BProb_Graf::BSumProb;
		}

	    if(MainData::order > 0){
			M_TrTree::gMTr = new M_TrTree();
			M_OVGraf::gMOVG = new M_OVGraf();
	        M_OVGraf::gMOVG->Copy(3,MinOV_Graf::gMinOVG);
	        M_OVGraf::gMOVG->Preprocessing();
			delete AC_Trie::gTrie;
	        delete MinOV_Graf::gMinOVG;
		    M_OVGraf::gMOVG->ProbCalc();
	        delete M_OVGraf::gMOVG;
			delete M_TrTree::gMTr;
	    }
	}
	if(MainData::mode == 1){
		delete[] MainData::RandPatProbs;
	}
	if((MainData::order > 0)&&(MainData::CrDistribFlag == 1)){
		for(i = 0; i < MainData::AlpSize; i++){
			delete[] MainData::MarkovProbs[i];
		}
		delete[] MainData::MarkovProbs;
	}
	if(ExitFlag == 1){
		int j;
		if((MainData::order == -2)&&(MainData::CrDistribFlag == 1)){
			for(i = 0; i < H_M_OVGraf::NumAllStates; i++){
				for(j = 0; j < H_M_OVGraf::NumAllStates; j++){
					delete[] MainData::ND_HHMProbs[i][j];
				}
				delete[] MainData::ND_HHMProbs[i];
				MainData::ND_HHMTrans[i]->clear();
				delete[] MainData::ND_HHMTrans[i];
			}
			delete[] MainData::ND_HHMProbs;
			delete[] MainData::ND_HHMTrans;
		}
		if((MainData::order == -1)&&(MainData::CrDistribFlag == 1)){
			for(i = 0; i < H_M_OVGraf::NumAllStates; i++){	
				delete[] MainData::D_HHMProbs[i];
				delete[] MainData::D_HHMTrans[i];
				delete[] H_M_OVGraf::TransStepProbMatrix[i];
			}
			delete[] MainData::D_HHMProbs;
			delete[] MainData::D_HHMTrans;
		}
		if(AC_Trie::gTrie != NULL){
			delete AC_Trie::gTrie;
		}
		return Error;
	}
	pvalue = &MainData::Pvalue;
	return 0;
}

