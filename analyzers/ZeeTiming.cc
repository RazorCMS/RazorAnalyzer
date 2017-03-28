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
  int NPU;
  int nPV;
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
    double ele1_time;
    double ele2_time;
    double ele1_seedtime;
    double ele2_seedtime;
    for(int i = 0; i < nElectrons; i++){
      // if(elePt[i] < 35) continue;
      // if(fabs(eleEta[i]) > 2.5) continue;
      // if(!(isEGammaPOGTightElectron(i))) continue;

      nEle++;
      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << "\n";
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
      double time = 0;
      double timeSeedHit = 0;
      
      double maxHitEnergy = 0;
      
      // cout << ele_NEcalRechitID[i] << "\n";
      for (int k=0; k<ele_NEcalRechitID[i]; ++k) {
      	
	uint tmpID = ele_EcalRechitID[i][k];
	// cout << "ele " << k << " hit: " << tmpID << " ";
	bool found = false;
	uint foundIndex = 0;
	for (int l=0; l < nEcalRechits; ++l) {
	  if (ecalRechit_ID[l] == tmpID) {
	    found = true;
	    foundIndex = l;
	    break;
	  }
	}	
	if (found) {
	  // cout << ecalRechit_E[foundIndex] << " " << ecalRechit_T[foundIndex] << " ";
	  //find the max energy hit as seed...for now
	  if (ecalRechit_E[foundIndex] > maxHitEnergy) {
	    maxHitEnergy = ecalRechit_E[foundIndex];
	    
	    double rawT = ecalRechit_T[foundIndex];
	    //correct for TOF
	    timeSeedHit = rawT + (std::sqrt(pow(ecalRechit_X[foundIndex],2)+pow(ecalRechit_Y[foundIndex],2)+pow(ecalRechit_Z[foundIndex],2))-std::sqrt(pow(ecalRechit_X[foundIndex]-pvX,2)+pow(ecalRechit_Y[foundIndex]-pvY,2)+pow(ecalRechit_Z[foundIndex]-pvZ,2)))/SPEED_OF_LIGHT;
	    timeSeedHit = rawT;
	  }
	}
	// cout << "\n";
      }
 
      if (thisElectron.Pt() > ele1.Pt()) {
	ele1 = thisElectron;
	ele1_time = time;
	ele1_seedtime = timeSeedHit;
      } else if (thisElectron.Pt() > ele2.Pt()) {
	ele2 = thisElectron;
	ele2_time = time;
	ele2_seedtime = timeSeedHit; 
      }	
    }
    
    if (nEle >= 2) {
      mass = (ele1+ele2).M();
      t1 = ele1_time;
      t2 = ele2_time;
      t1_seed = ele1_seedtime;
      t2_seed = ele2_seedtime;
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
