#ifdef __cplusplus
extern "C" int func_set_input_data( int order, int mode, int TLen, int NOccur);
extern "C" int func_analis_alp_bern_data(int AlpSize, char* AlpMas, int order, double* BernProb);
extern "C" int func_analis_alp_mark_data(int AlpSize, char* AlpMas, double** MarkovProbs);
extern "C" int func_analis_alp_dhhm_data(int AlpSize, char* AlpMas, int NumAllStates,
							  double** D_HHMProbs, int** D_HHMTrans);
extern "C" int func_analis_alp_hhm_data(int AlpSize, char* AlpMas,
							 int NumAllStates, double*** ND_HHMProbs);
extern "C" int func_analis_pattern_data_0(int NWords, char **WordsList);
extern "C" int func_analis_pattern_data_1(int NWords, int WordLen, double* RandPatProbs);
extern "C" int func_analis_pattern_data_2_3(int WordLen, int NFootPrints, char **FootPrints, 
										double** PssmMas, double Thr);
extern "C" int func_analis_pattern_data_4(char* motif, int Nreplace, int NConstPositions, 
									  int *ConstPositions);
extern "C" int func_analis_pattern_data_5(char *consensus, int NSymbols, char **ConsAlp);
extern "C" int func_main(double* pvalue);
#else
extern int func_set_input_data( int order, int mode, int TLen, int NOccur);
extern int func_analis_alp_bern_data(int AlpSize, char* AlpMas, int order, double* BernProb);
extern int func_analis_alp_mark_data(int AlpSize, char* AlpMas, double** MarkovProbs);
extern int func_analis_alp_dhhm_data(int AlpSize, char* AlpMas, int NumAllStates,
							  double** D_HHMProbs, int** D_HHMTrans);
extern int func_analis_alp_hhm_data(int AlpSize, char* AlpMas,
							 int NumAllStates, double*** ND_HHMProbs);
extern int func_analis_pattern_data_0(int NWords, char **WordsList);
extern int func_analis_pattern_data_1(int NWords, int WordLen, double* RandPatProbs);
extern int func_analis_pattern_data_2_3(int WordLen, int NFootPrints, char **FootPrints, 
										double** PssmMas, double Thr);
extern int func_analis_pattern_data_4(char* motif, int Nreplace, int NConstPositions, 
									  int *ConstPositions);
extern int func_analis_pattern_data_5(char *consensus, int NSymbols, char **ConsAlp);
extern int func_main(double* pvalue);
#endif

