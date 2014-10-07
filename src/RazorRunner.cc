//C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

#include <RazorAnalyzer.h>

using namespace std;

int main(int argc, char* argv[]){

    //get input files and analysis type from command line
    if(argc < 3){
        cerr << "Usage: RazorRunner <input list> <analysis type>" << endl;
        cerr << "Analyses available: " << endl << "dummy   --   do nothing useful" << endl;
        return -1;
    }
    string inputFileName(argv[1]);
    string analysisType(argv[2]);

    //build the TChain
    TChain *theChain = new TChain("ntuples/outputTree");
    string curFileName;
    ifstream *inputFile = new ifstream(inputFileName.c_str());
    while(getline(*inputFile, curFileName)){
        theChain->Add(curFileName.c_str());
        std::cout << "chaining " << curFileName << std::endl;
    }

    RazorAnalyzer analyzer(theChain);
    
    //enable desired branches and execute analysis
    if(analysisType == "dummy"){
        cout << "Executing dummy analysis..." << endl;
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.DummyAnalysis();
    }
    else{ //analysis not found
        cerr << "Error: the given analysis type is not defined in RazorTestAnalysis.cc!" << endl;
    }

    //this is a joke; feel free to remove it.
    cerr << "rm: cannot remove `thisiswhyitcrashed*': No such file or directory" << endl;

    return 0;
}
