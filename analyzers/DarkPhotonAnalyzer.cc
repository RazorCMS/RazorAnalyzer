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

void DarkPhotonAnalyzer::Analyze(bool isData, int Option, string outputFilename, string label)
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
  float l1_PT;
  float l1_Eta;
  float l1_Phi;
  float l2_PT;
  float l2_Eta;
  float l2_Phi;
  int Iso_lepton;
  // working variables:
  float Copy_Jet_PT[900];
  float Copy_Electron_PT[700];
  float Copy_Muon_PT[700];
  
  TTree *outputTree = new TTree("VBFTree", "Info on selected razor inclusive events");

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
  outputTree->Branch("l1_PT", &l1_PT, "l1_PT/F");
  outputTree->Branch("l1_Eta", &l1_Eta, "l1_Eta/F");
  outputTree->Branch("l1_Phi", &l1_Phi, "l1_Phi/F");
  outputTree->Branch("l2_PT", &l2_PT, "l2_PT/F");
  outputTree->Branch("l2_Eta", &l2_Eta, "l2_Eta/F");
  outputTree->Branch("l2_Phi", &l2_Phi, "l2_Phi/F");
  outputTree->Branch("Iso_lepton", &Iso_lepton, "Iso_lepton/i");	// is there an isolated lepton? (defined according to cuts standards) 1=YES, 0=NO

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
	
	//Fill & select lepton PT, Eta and Phi for the two highest PT leptons here
	// copy arrays
	for(int i=0; i<700; i++) {
		Copy_Electron_PT[i] = elePt[i];	
	}
	for(int i=0; i<700; i++) {
		Copy_Muon_PT[i] = muonPt[i];	
	}
	// find indices
	int e1_index = TMath::LocMax(700,elePt);	// find e1 index
	Copy_Electron_PT[e1_index] = 0;				// set e1 value in Copy to zero
	int e2_index = TMath::LocMax(700,Copy_Electron_PT);	// find e2 index
	//
	int m1_index = TMath::LocMax(700,muonPt);	// find m1 index
	Copy_Muon_PT[m1_index] = 0;				// set m1 value in Copy to zero
	int m2_index = TMath::LocMax(700,Copy_Muon_PT);	// find m2 index
	//	 	
	if ( (muonPt[m2_index]<20) && (elePt[e2_index]>20) ) {
		// select the two electrons
		l1_PT=elePt[e1_index];
		l1_Eta=eleEta[e1_index];
		l1_Phi=elePhi[e1_index];
		l2_PT=elePt[e2_index];
		l2_Eta=eleEta[e2_index];
		l2_Phi=elePhi[e2_index];
	} else if ( (muonPt[m2_index]>20) && (elePt[e2_index]<20) ) {
		// select the two muons
		l1_PT=muonPt[m1_index];
		l1_Eta=muonEta[m1_index];
		l1_Phi=muonPhi[m1_index];
		l2_PT=muonPt[m2_index];
		l2_Eta=muonEta[m2_index];
		l2_Phi=muonPhi[m2_index];
	} else if ( (muonPt[m2_index]<20) && (elePt[e2_index]<20) ) {
		continue;
	} else if ( (muonPt[m2_index]>20) && (elePt[e2_index]>20) ) {				// not sure how to deal with both muons and electrons being possible candidate leptons
		// select the two muons
		l1_PT=muonPt[m1_index];
		l1_Eta=muonEta[m1_index];
		l1_Phi=muonPhi[m1_index];
		l2_PT=muonPt[m2_index];
		l2_Eta=muonEta[m2_index];
		l2_Phi=muonPhi[m2_index];	// default to selecting the muons at the moment
	}
    
	// Isolated lepton criterion
	Iso_lepton = 0;
	Float_t Phi_l;
	Float_t Eta_l;
	Float_t DR;
	Float_t Pi = TMath::Pi();
	for (Int_t j=0; j<700; j++) {	
		Phi_l = TMath::Abs(elePhi[j] - PhotonPhi);
		if (Phi_l > Pi) {
			Phi_l = 2*Pi - Phi_l;
		}
		Eta_l = TMath::Abs(eleEta[j] - PhotonEta);
		DR = sqrt(pow(Phi_l,2) + pow((Eta_l - PhotonEta),2));
		if (DR > 0.3) {									// Isolated electron criterion
			Iso_lepton = 1;
		}
	}
	for (Int_t j=0; j<700; j++) {	
		Phi_l = TMath::Abs(muonPhi[j] - PhotonPhi);
		Eta_l = TMath::Abs(muonEta[j] - PhotonEta);
		DR = sqrt(pow(Phi_l,2) + pow(Eta_l,2));
		if (DR > 0.3) {									// Isolated muon criterion
			Iso_lepton = 1;
		}
	}
	
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
