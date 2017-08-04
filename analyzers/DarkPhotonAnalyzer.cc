//LOCAL INCLUDES
#include "DarkPhotonAnalyzer.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 
//ROOT INCLUDES
#include <TH1F.h>
#include <TMath.h>
#include <TH2D.h>
#include "TRandom3.h"

using namespace std;

const int NUM_PDF_WEIGHTS = 60;

void PhotonNtupler::Analyze(bool isData, int Option, string outputFilename, string label)
{
  cout << "Initializing..." << endl;

  //--------------------------------
  //Initialize helper
  //--------------------------------
  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;
  bool isFastsimSMS = false;
  RazorHelper *helper = 0;
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, isFastsimSMS);
  else if (analysisTag == "Razor2016_80X") helper = new RazorHelper("Razor2016_80X", isData, isFastsimSMS);
  else helper = new RazorHelper(analysisTag, isData, isFastsimSMS);


  //--------------------------------
  //Initialize Output
  //--------------------------------
  string outfilename = outputFilename;
  if (outfilename == "") outfilename = "/afs/cern.ch/user/j/jbamber/scratch1/CMSSW_7_4_15/src/ROOT_outputs/DarkPhoton.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");

  cout << "Run With Option = " << Option << "\n";

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);
    
  //--------------
  //tree variables
  //--------------
  unsigned int run, lumi, event;
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float triggerEffWeight;
  float triggerEffSFWeight;
  float photonEffSF;

  int NPU;
  float MET, METPhi;
  float PhotonPt, PhotonEta, PhotonPhi;
  float MT;
  float j1_Eta;
  float j2_Eta;
  // working variables:
  float Copy_Jet_PT[900];

  TTree *outputTree = new TTree("HggRazor", "Info on selected razor inclusive events");

  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("triggerEffWeight", &triggerEffWeight, "triggerEffWeight/F");
  outputTree->Branch("triggerEffSFWeight", &triggerEffSFWeight, "triggerEffSFWeight/F");
  outputTree->Branch("photonEffSF", &photonEffSF, "photonEffSF/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("MET", &MET, "MET/F");
  outputTree->Branch("METPhi", &METPhi, "METPhi/F");
  outputTree->Branch("PhotonPt", &PhotonPt, "PhotonPt/F");
  outputTree->Branch("PhotonEta", &PhotonEta, "PhotonEta/F");
  outputTree->Branch("PhotonPhi", &PhotonPhi, "PhotonPhi/F");
  outputTree->Branch("MT", &MT, "MT/F");
  outputTree->Branch("j1_Eta", &j1_Eta, "j1_Eta/F");
  outputTree->Branch("j2_Eta", &j2_Eta, "j2_Eta/F");

  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double leadingPhotonPt = -999;
    double leadingPhotonEta = -999;
    double leadingPhotonPhi = -999;

    // define variables for finding the highest PT photon
    int i_best=0;
    
    //Select highest pT photon
    for(int i = 0; i < nPhotons; i++){

      //require medium photon selection
      if (!isMediumPhoton(i)) continue;

      //find the highest pT photon here
      //*******
      if (phoPt[i] > phoPt[i_best]) {
      	  i_best = i;
      }
      //*******
    }

    //Fill MET, MetPhi, MT here
    MET = metPt;
    METPhi = metPhi;
    PhotonPt = phoPt[i_best];
    PhotonEta = phoEta[i_best];
    PhotonPhi = phoPhi[i_best];
    
    MT = sqrt(2*MET*PhotonPt*(1 - cos(PhotonPhi - METPhi)));
    
    //Fill j1, j2 eta values here
    
    // ## Find two jets with highest pT: j1, j2 with j1_pT > j2_pT
	// make array copy of Jet
	for(int i=0; i<900; i++) {
		Copy_Jet_PT[i] = jetPt[i];	
	}
	// find indices
	int j1_index = TMath::LocMax(900,jetPt);	// find j1 index
	Copy_Jet_PT[j1_index] = 0;				// set j1 value in Copy to zero
	int j2_index = TMath::LocMax(900,Copy_Jet_PT);	// find j2 index
	// assign variables
	j1_Eta = jetEta[j1_index];
	j2_Eta = jetEta[j2_index];
    
    //fill normalization histogram    
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;
    SumWeights->Fill(1.0, weight);    
    run = runNum;
    lumi = lumiNum; 
    event = eventNum; 

    //------------------
    //Pileup reweighting
    //------------------
    pileupWeight = 1.0;
    if( !isData ) {
      //Get number of PU interactions
      for (int i = 0; i < nBunchXing; i++) {
	if (BunchXing[i] == 0) {
	  NPU = nPUmean[i];
	}
      }
      pileupWeight = helper->getPileupWeight(NPU);
      pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
      pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;	
    }
    
    /////////////////////////////////
    //Scale and PDF variations
    /////////////////////////////////
    if( !isData ) {
      if ( (*scaleWeights).size() >= 9 ) 
	{
	  // sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
	  // sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
	  // sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
	  // sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
	  // sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
	  // sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;
	    
	    
	  SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
	  SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
	  SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
	  SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
	  SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
	  SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
	}
	
      // sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
      for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) {
      // 	  sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
      }
    }
    
    // Put Trigger stuff here ?
      
    outputTree->Fill();

  } // loop over events 
    
  cout << "Writing output trees..." << endl;
  outFile->Write();
  outFile->Close();
    
}
