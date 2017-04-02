#ifndef DEF_ZeeTiming
#define DEF_ZeeTiming

#include "RazorAnalyzer.h"

class ZeeTiming: public RazorAnalyzer {
    public: 
        uint start_run_tmp;
        uint end_run_tmp;	
        vector <float> *IC_time_all;
        vector <int> *detID_all;

        ZeeTiming(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
	float getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID);
};

#endif
