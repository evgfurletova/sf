//#include "stdafx.h"
#include "maindata.h"
#include "ac_trie.h"
#include "ov_graf.h"
#include "minov_graf.h"
#include "bprob_graf.h"
#include "mmodel_prob.h"
#include "m_ovgraf.h" 
#include "m_trtree.h"
#include "nd_hhm_prob.h"
#include "d_hhm_prob.h"
#include "hhm_ovgraf.h"
#include "sf_main.h"

int main(int argc, char* argv[])
{
	int i;
	time_t time0 = time(NULL);
	int ExitFlag = 0;
    int Error = MainData::ComLineParse(argc, argv);
    if(Error > 0)
        ExitFlag = 1;

	if(ExitFlag == 0){
		Error = MainData::GetInput();
		if(Error > 0)
			ExitFlag = 1;
	}
	
	ofstream out;
	out.open(MainData::OutName.c_str(), ios::app);

	if((MainData::NOccur == 0)&&(ExitFlag == 0)){
		out<<1<<'\n';
		out.close();
		ExitFlag = 1;
	}
	
	if((MainData::TLen < MainData::WordLen)&&(ExitFlag == 0)){
		out<<0<<'\n';
		out.close();
		ExitFlag = 1;
	}
	if((MainData::TLen == MainData::WordLen)&&(MainData::NOccur > 1)&&(ExitFlag == 0)){
		out<<0<<'\n';
		out.close();
		ExitFlag = 1;
	}

	if((MainData::AlpSize == 1)&&(ExitFlag == 0)){
		if(MainData::TLen >= MainData::WordLen + MainData::NOccur - 1){
			out<<1<<'\n';
			out.close();
			ExitFlag = 1;
		}
		else{
			out<<0<<'\n';
			out.close();
			ExitFlag = 1;
		}
	}

    if((MainData::order > 0)&&(ExitFlag == 0)){
		if(MainData::order > MainData::WordLen){
			cout<<"The pattern length has to be equal or bigger then order of the Markovian model";
			ExitFlag =  1;
		}
        M_TrTree::ForLeaves = new int[MainData::NWords];
    }

///////// 
	if(ExitFlag == 0){
			std::string FileName = MainData::OutName;
			int s = (int)FileName.length();
			if(FileName.at(s - 4) == '.'){
				FileName.erase(s-4,4);
			}
			FileName += "_Words.txt";
			if(MainData::mode == 1){
			    MainData::ResWords.open(FileName.c_str());
			    AC_Trie::gTrie->PrintTrie(&MainData::ResWords);
			    MainData::ResWords.close();
			}
	}
/////////
    

	if(ExitFlag == 0){
		AC_Trie::CreateTrie();
		
		/*OV_Graf::gOVG = new OV_Graf(0);
		OV_Graf::CreateGraf();
		OV_Graf::gOVG->PrintLOG(0,&MainData::ResLOG);	
		MinOV_Graf::gMinOVG = new MinOV_Graf();
		MinOV_Graf::gMinOVG->Copy(1,OV_Graf::gOVG);
		delete OV_Graf::gOVG;
		MinOV_Graf::MinimizeGraf();
*/

	    if(MainData::order < 0){
			HHM_OVGraf::gHHM = new HHM_OVGraf();
			OV_Graf::CreateGraf(HHM_OVGraf::gHHM);
	        HHM_OVGraf::gHHM->Preprocessing();
	        delete AC_Trie::gTrie;
			delete MinOV_Graf::gMinOVG;
			 HHM_OVGraf::gHHM->ProbCalc();
			delete HHM_OVGraf::gHHM;
		}   
		if(MainData::order == 0){
			BProb_Graf::gBOVG = new BProb_Graf(0);
			//BProb_Graf::gBOVG->Copy(2,MinOV_Graf::gMinOVG);
			//BProb_Graf::PredStep();
			OV_Graf::CreateGraf(BProb_Graf::gBOVG);
			
			delete AC_Trie::gTrie;
			//delete MinOV_Graf::gMinOVG;
			BProb_Graf::ProbCalc();
			delete BProb_Graf::gBOVG;
			delete[] BProb_Graf::PHCl;
			delete[] BProb_Graf::BSumProb;
		}

	    if(MainData::order > 0){
			M_TrTree::gMTr = new M_TrTree();
			M_OVGraf::gMOVG = new M_OVGraf();
			OV_Graf::CreateGraf(M_OVGraf::gMOVG);
	        M_OVGraf::gMOVG->Preprocessing();
			delete AC_Trie::gTrie;
		    M_OVGraf::gMOVG->ProbCalc();
	        delete M_OVGraf::gMOVG;
			delete M_TrTree::gMTr;
	    }
		out<<MainData::Pvalue;
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
			}
			delete[] MainData::D_HHMProbs;
			delete[] MainData::D_HHMTrans;
		}
		if(AC_Trie::gTrie != NULL){
			delete AC_Trie::gTrie;
		}
		MainData::ErrorDetect(Error);
		return 1;
	}
	time0 = time(NULL) - time0;
	out<<'\t'<<time0<<'\n';
	out.close();
	return 0;
}
