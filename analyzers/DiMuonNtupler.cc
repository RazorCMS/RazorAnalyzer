#include "DiMuonNtupler.h"
#include "DiMuonTree.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void DiMuonNtupler::Analyze(bool isData, int Option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  cout << "Getting JEC parameters from " << pathname << endl;

  std::vector<JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "DiMuonNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    DiMuonTree* mumuTree = new DiMuonTree;
    mumuTree->CreateTree();
    mumuTree->tree_->SetAutoFlush(0);

    std::cout << "Run With Option = " << Option << "\n";

    UInt_t NDiMuonsFilled = 0;

    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    cout << "nentries = " << nentries << "\n";
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {
        //***********************
        //Fill Jet
        //***********************
        //NJetsFilled++;
        mumuTree->tree_->Fill();
      }
    }
    //end of event loop

    //cout << "Filled Total of " << NJetsFilled << " Jets\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}
