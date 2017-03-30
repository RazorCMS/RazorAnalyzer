#include "ZeeTiming.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns

void ZeeTiming::Analyze(bool isData, int option, string outFileName, string label)
{

  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  bool doPhotonScaleCorrection = true;

  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;
  
  if ( outFileName.empty() ) {
    std::cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << std::endl;
    outFileName = "HggRazor.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );

  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("ZeeTiming", "Info on selected razor inclusive events");

  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float mass;
  float t1, t2;
  float t1_seed, t2_seed;
  float t1raw_seed, t2raw_seed;
  float ele1E, ele1Pt, ele1Eta, ele1Phi;
  float ele2E, ele2Pt, ele2Eta, ele2Phi;
  int NPU;
  //int nPV;
  unsigned int run, lumi, event;


  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("nPV", &nPV, "nPV/i");
  outputTree->Branch("mass", &mass, "mass/F");
  outputTree->Branch("t1", &t1, "t1/F");
  outputTree->Branch("t2", &t2, "t2/F");
  outputTree->Branch("t1_seed", &t1_seed, "t1_seed/F");
  outputTree->Branch("t2_seed", &t2_seed, "t2_seed/F");
  outputTree->Branch("t1raw_seed", &t1raw_seed, "t1raw_seed/F");
  outputTree->Branch("t2raw_seed", &t2raw_seed, "t2raw_seed/F");
  outputTree->Branch("ele1E", &ele1E, "ele1E/F");
  outputTree->Branch("ele1Pt", &ele1Pt, "ele1Pt/F");
  outputTree->Branch("ele1Eta", &ele1Eta, "ele1Eta/F");
  outputTree->Branch("ele1Phi", &ele1Phi, "ele1Phi/F");
  outputTree->Branch("ele2E", &ele2E, "ele2E/F");
  outputTree->Branch("ele2Pt", &ele2Pt, "ele2Pt/F");
  outputTree->Branch("ele2Eta", &ele2Eta, "ele2Eta/F");
  outputTree->Branch("ele2Phi", &ele2Phi, "ele2Phi/F");

  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);


  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //begin event
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //fill normalization histogram
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;


    //get NPU
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	NPU = nPUmean[i];
      }   
    }
    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    
    double pvX = 0;


    int nEle = 0;
    TLorentzVector ele1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector ele2 = makeTLorentzVector(0,0,0,0);
    double ele1_time = 0;
    double ele2_time = 0;
    double ele1_seedtime = 0;
    double ele2_seedtime = 0;
    double ele1_seedtimeraw = 0;
    double ele2_seedtimeraw = 0;
    for(int i = 0; i < nElectrons; i++){
      // if(elePt[i] < 35) continue;
      // if(fabs(eleEta[i]) > 2.5) continue;
      // if(!(isEGammaPOGTightElectron(i))) continue;

      nEle++;
      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << "\n";
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
      double time = 0;
    
      
      uint seedhitIndex =  (*ele_SeedRechitIndex)[i];
      double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];
      double timeSeedHit = rawSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;;

      // cout << ele_NEcalRechitID[i] << "\n";
      for (int k=0; k<(*ele_EcalRechitIndex)[i].size(); ++k) {
      	
	uint rechitIndex = (*ele_EcalRechitIndex)[i][k];
		  
	double rawT = (*ecalRechit_T)[rechitIndex];

	//correct for TOF
	double corrT = rawT + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
	 	
	// cout << "\n";
      }
            
      if (thisElectron.Pt() > ele1.Pt()) {
	ele1 = thisElectron;
	ele1_time = time;
	ele1_seedtime = timeSeedHit;
	ele1_seedtimeraw = rawSeedHitTime;
      } else if (thisElectron.Pt() > ele2.Pt()) {
	ele2 = thisElectron;
	ele2_time = time;
	ele2_seedtime = timeSeedHit; 
 	ele2_seedtimeraw = rawSeedHitTime;
     }	
    }
    
    if (nEle >= 2) {
      ele1E = ele1.E();
      ele1Pt = ele1.Pt();
      ele1Eta = ele1.Eta();
      ele1Phi = ele1.Phi();
      ele2E = ele2.E();
      ele2Pt = ele2.Pt();
      ele2Eta = ele2.Eta();
      ele2Phi = ele2.Phi();

      mass = (ele1+ele2).M();
      t1 = ele1_time;
      t2 = ele2_time;
      t1_seed = ele1_seedtime;
      t2_seed = ele2_seedtime;
      t1raw_seed = ele1_seedtimeraw;
      t2raw_seed = ele2_seedtimeraw;
       //cout << "ele2: " << ele2.Pt() << " " << ele2_seedtime << "\n";
    }

    //Fill Event
    if (mass > 60 && mass < 120) {
      outputTree->Fill();
    }

  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  outputTree->Write();
  NEvents->Write();

  outFile->Close();
}
