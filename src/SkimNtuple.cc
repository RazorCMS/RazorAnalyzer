#include <fstream>
#include <sstream>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeFormula.h"
#include "SimpleTable.h"
#include "TKey.h"

using namespace std;


//get list of files to open, add normalization branch to the tree in each file
int main(int argc, char* argv[]) {

    //parse input list to get names of ROOT files
    if(argc < 5){
        cerr << "usage SkimNtuple inputList.txt <outputDirectory> <outputfileLabel> <Skim Cut String>" << endl;
        return -1;
    }
    string inputList(argv[1]);
    
    string outputDir = argv[2];
    string outputfileLabel = argv[3];
    string SkimCutString = argv[4];

    ifstream filein(inputList.c_str());
    string curFilename;
    vector<string> inputLines;
    while(getline(filein, curFilename)){
        if(curFilename.at(0) != '#') inputLines.push_back(curFilename); //'#' denotes a comment
        else cout << "(Skipping commented line in input)" << endl;
    }

    //open each ROOT file and add the normalization branch
    for(auto& line : inputLines){
        //parse input -- input lines should be in the form datasetName fileName
        istringstream buf(line);
        istream_iterator<std::string> beg(buf), end;
        vector<std::string> inputs(beg, end);
        
        string fileName = inputs[0];

        //create output file
	string outputfilename = Form("%s/%s_%s.root", outputDir.c_str(), 
				     fileName.substr(fileName.find_last_of("/")+1, fileName.find_last_of(".")).c_str(),
				     outputfileLabel.c_str());
	cout << "Output file: " << outputfilename << "\n";
        TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");
	
        //loop over all TTrees in the file and add the weight branch to each of them
        TFile inputFile(fileName.c_str(), "READ");
        inputFile.cd();
        inputFile.Purge(); //purge unwanted TTree cycles in file
        TIter nextkey(inputFile.GetListOfKeys());
        TKey *key;
        while((key = (TKey*)nextkey())){
            string className = key->GetClassName();
            cout << "Getting key from file.  Class type: " << className << endl;
            if(className.compare("TTree") != 0){
                cout << "Skipping key (not a TTree)" << endl;
                continue;
            }

            TTree *inputTree = (TTree*)key->ReadObj();
            cout << "Processing tree " << inputTree->GetName() << endl;

            //create new normalized tree
            outputFile->cd();
            TTree *outputTree = inputTree->CloneTree(0);  
            cout << "Events in the ntuple: " << inputTree->GetEntries() << endl;

	    TTreeFormula *formula = new TTreeFormula("SkimCutString", SkimCutString.c_str(), inputTree);
	    int EventsPassed = 0;

            //store the weights
            for (int n=0;n<inputTree->GetEntries();n++) { 
	      if (n%1000000==0) cout << "Processed Event " << n << "\n";
                inputTree->GetEntry(n);

		bool passSkim = false;
		
		passSkim = formula->EvalInstance();

		if (passSkim) {
		  EventsPassed++;
		  outputTree->Fill(); 
		}
            }

	    delete formula;
	    cout << "Skim Efficiency : " << EventsPassed << " / " << inputTree->GetEntries() 
		 << " = " << float(EventsPassed ) / float(inputTree->GetEntries()) 
		 << " \n";

            //save
            outputTree->Write();
            inputFile.cd();
        }
        inputFile.Close();
        cout << "Closing output file." << endl;

        outputFile->Close();
        delete outputFile;
    }
}
