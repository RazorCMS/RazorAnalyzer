#ifndef DEF_ZeeTiming
#define DEF_ZeeTiming

#include "RazorAnalyzer.h"

class ZeeTiming: public RazorAnalyzer {
    public: 
        ZeeTiming(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
