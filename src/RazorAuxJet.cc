#include "RazorAnalyzer.h"

//B-tag working points for 53X Summer13 studies at 8 TeV (with 22Jan2013 ReReco Data)
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
bool RazorAnalyzer::isCSVL(int i){
    return jetCSV[i] > 0.244;
}

bool RazorAnalyzer::isCSVM(int i){
    return jetCSV[i] > 0.679;
}

bool RazorAnalyzer::isCSVT(int i){
    return jetCSV[i] > 0.898;
}
