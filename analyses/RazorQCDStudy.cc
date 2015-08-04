
#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
#include "ControlSampleEvents.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};
 
void RazorAnalyzer::RazorQCDStudy( string outputfilename, int option, bool isData, bool isRunOne)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    cout << "IsData = " << isData << "\n";

    TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    bool printSyncDebug = false;
    std::vector<JetCorrectorParameters> correctionParameters;

    if (isRunOne) {
      if (isData) {
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Winter14_V8_DATA_L1FastJet_AK5PF.txt"));
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Winter14_V8_DATA_L2Relative_AK5PF.txt"));
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Winter14_V8_DATA_L3Absolute_AK5PF.txt")); 
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Winter14_V8_DATA_L2L3Residual_AK5PF.txt")); 
      } else {
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Summer13_V4_MC_L1FastJet_AK5PF.txt"));
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Summer13_V4_MC_L2Relative_AK5PF.txt"));
	correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Summer13_V4_MC_L3Absolute_AK5PF.txt")); 
      }
    } else {
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L2Relative_AK4PFchs.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt"));    
    }
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
    JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/JetResolutionInputAK5PF.txt");
    SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);


    //*************************************************************************
    //Set up Output File
    //*************************************************************************
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorControlRegions.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    TTree *outTree = new TTree("QCDTree", "Info on selected razor inclusive events");
     //tree variables
    Float_t                 weight;
    UInt_t                  run;
    UInt_t                  lumi;
    UInt_t                  event;
    UInt_t                  NPU_0;
    UInt_t                  NPU_Minus1;
    UInt_t                  NPU_Plus1;
    UInt_t                  NPV;
    Float_t                 Rho;
    Bool_t                  HLTDecision[100];
    Float_t                 MR;
    Float_t                 Rsq;
    Float_t                 minDPhi;
    Float_t                 minDPhiN;
    Float_t                 dPhiRazor;
    Float_t                 MET;
    Float_t                 HT;
    UInt_t                  NJets40;
    UInt_t                  NJets80;
    UInt_t                  NBJetsLoose;
    UInt_t                  NBJetsMedium;
    UInt_t                  NBJetsTight;
    Float_t                 genJetMR;
    Float_t                 genJetRsq;
    Float_t                 genJetDPhiRazor;    
    Float_t                 genJetHT;
    int NJets;
    float JetE[99];
    float JetPt[99];
    float JetEta[99];
    float JetPhi[99];
    float GenJetE[99];
    float GenJetPt[99];
    float GenJetEta[99];
    float GenJetPhi[99];
    bool  JetIDTight[99];
    float PtNeutrinoClosestToJet[99]; 
    float DRNeutrinoClosestToJet[99];
    float JetEMFraction[99];
    float JetPartonFlavor[99];
    float JetPileupID[99];
    
    //book the branches that go in all types of trees
    outTree->Branch("weight",&weight,"weight/F");
    outTree->Branch("run",&run,"run/i");
    outTree->Branch("lumi",&lumi,"lumi/i");
    outTree->Branch("event",&event,"event/i");
    outTree->Branch("NPU_0",&NPU_0,"NPU_0/i");
    outTree->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
    outTree->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
    outTree->Branch("NPV",&NPV,"NPV/i");
    outTree->Branch("Rho",&Rho,"Rho/F");
    outTree->Branch("HLTDecision",&HLTDecision,"HLTDecision[100]/O");
    outTree->Branch("MR",&MR,"MR/F");
    outTree->Branch("Rsq",&Rsq,"Rsq/F");
    outTree->Branch("minDPhi",&minDPhi,"minDPhi/F"); 
    outTree->Branch("minDPhiN",&minDPhiN,"minDPhiN/F"); 
    outTree->Branch("dPhiRazor",&dPhiRazor,"dPhiRazor/F");
    outTree->Branch("MET",&MET,"MET/F");
    outTree->Branch("HT",&HT,"HT/F");	  
    outTree->Branch("NJets40",&NJets40,"NJets40/i");
    outTree->Branch("NJets80",&NJets80,"NJets80/i");
    outTree->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
    outTree->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
    outTree->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
    outTree->Branch("genJetMR",&genJetMR,"genJetMR/F");
    outTree->Branch("genJetRsq",&genJetRsq,"genJetRsq/F");
    outTree->Branch("genJetDPhiRazor",&genJetDPhiRazor,"genJetDPhiRazor/F");
    outTree->Branch("genJetHT",&genJetHT,"genJetHT/F");	  
    outTree->Branch("genMET",&genMetPt,"genMET/F");	         
    outTree->Branch("genMETPhi",&genMetPhi,"genMETPhi/F");	         
    outTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
    outTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter,"Flag_CSCTightHaloFilter/O");
    outTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter,"Flag_hcalLaserEventFilter/O");
    outTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    outTree->Branch("Flag_goodVertices", &Flag_goodVertices,"Flag_goodVertices/O");
    outTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter,"Flag_trackingFailureFilter/O");
    outTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
    outTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
    outTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters,"Flag_trkPOGFilters/O");
    outTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X,"Flag_trkPOG_manystripclus53X/O");
    outTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X,"Flag_trkPOG_toomanystripclus53X/O");
    outTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters,"Flag_trkPOG_logErrorTooManyClusters/O");
    outTree->Branch("Flag_METFilters", &Flag_METFilters,"Flag_METFilters/O");	
    outTree->Branch("nJets", &NJets,"nJets/I");
    outTree->Branch("JetE", JetE,"JetE[nJets]/F");
    outTree->Branch("JetPt", JetPt,"JetPt[nJets]/F");
    outTree->Branch("JetEta", JetEta,"JetEta[nJets]/F");
    outTree->Branch("JetPhi", JetPhi,"JetPhi[nJets]/F");
    outTree->Branch("GenJetE", GenJetE,"GenJetE[nJets]/F");
    outTree->Branch("GenJetPt", GenJetPt,"GenJetPt[nJets]/F");
    outTree->Branch("GenJetEta", GenJetEta,"GenJetEta[nJets]/F");
    outTree->Branch("GenJetPhi", GenJetPhi,"GenJetPhi[nJets]/F");
    outTree->Branch("JetIDTight", JetIDTight,"JetIDTight[nJets]/O");
    outTree->Branch("PtNeutrinoClosestToJet", PtNeutrinoClosestToJet,"PtNeutrinoClosestToJet[nJets]/F");
    outTree->Branch("DRNeutrinoClosestToJet", DRNeutrinoClosestToJet,"DRNeutrinoClosestToJet[nJets]/F");
    outTree->Branch("JetEMFraction", JetEMFraction,"JetEMFraction[nJets]/F");
    outTree->Branch("JetPartonFlavor", JetPartonFlavor,"JetPartonFlavor[nJets]/F");
    outTree->Branch("JetPileupID", JetPileupID,"JetPileupID[nJets]/F");


    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  

    //*************************************************************************
    //Look over Input File Events
    //*************************************************************************
    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    //for (Long64_t jentry=0; jentry<1000;jentry++) {
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      printSyncDebug = false;

      //fill normalization histogram
      NEvents->Fill(1.0);
     
	//event info
      weight = 1.0;
      run = runNum;
      lumi = lumiNum;
      event = eventNum;
      
      //get NPU
      for (int i=0; i < nBunchXing; ++i) {
	if (BunchXing[i] == 0) {
	  NPU_0 = nPUmean[i];
	}
	if (BunchXing[i] == -1) {
	  NPU_Minus1 = nPUmean[i];
	}
	if (BunchXing[i] == 1) {
	  NPU_Plus1 = nPUmean[i];
	}	  
      }
      NPV = nPV;
 
      //************************************************************************
      //Find all Relevent Jets	
      //************************************************************************
      vector<TLorentzVector> GoodJets;
      int numJetsAbove80GeV = 0;
      int numJetsAbove40GeV = 0;
      double MetX_Type1Corr = 0;
      double MetY_Type1Corr = 0;
      int nBJetsLoose20GeV = 0;
      int nBJetsMedium20GeV = 0;
      int nBJetsTight20GeV = 0;
      minDPhi = 9999;
      minDPhiN = 9999;

      if (printSyncDebug) cout << "NJets: " << nJets << "\n";
      NJets = nJets;
      for(int i = 0; i < nJets; i++){

	//*******************************************************
	//apply jet iD
	//*******************************************************
	//if (!jetPassIDTight[i]) continue;

	//*******************************************************
	//Correct Jet Energy Scale and Resolution
	//*******************************************************
	double tmpRho = fixedGridRhoFastjetAll;
	if (isRunOne) tmpRho = fixedGridRhoAll;
	double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i], 
					       JetCorrector);   
	Rho = tmpRho;

	double jetEnergySmearFactor = 1.0;
	if (!isData) {
	  std::vector<float> fJetEta, fJetPtNPU;
	  fJetEta.push_back(jetEta[i]);  
	  fJetPtNPU.push_back(jetPt[i]*JEC); 
	  fJetPtNPU.push_back(NPU_0); 
	  if (printSyncDebug) {
	    cout << "Jet: " << jetPt[i] << " " << jetEta[i] << " " << jetPhi[i] << "\n";
	    cout << "Jet Resolution : " << jetPt[i]*JEC << " " << jetEta[i] << " " << jetPhi[i] << " : " 
		 << JetResolutionCalculator->resolution(fJetEta,fJetPtNPU) << "\n";
	  }
	  jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], NPU_0, JetResolutionCalculator, random);
	}
	if (printSyncDebug) {
	  cout << "Jet Smearing Factor " << jetEnergySmearFactor << "\n";
	}

	TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
	TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

	//*******************************
	//Add to Type1 Met Correction
	//*******************************
	if (jetPt[i]*JEC > 20) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - UnCorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - UnCorrJet.Py()  );
	  if (printSyncDebug) cout << "Met Type1 Corr: " << thisJet.Px() - UnCorrJet.Px() << " " << thisJet.Py() - UnCorrJet.Py() << "\n";
	}

	//*******************************************************
	//Jet Variab;es
	//*******************************************************
	JetE[i] = jetE[i]*JEC*jetEnergySmearFactor;
	JetPt[i] = jetPt[i]*JEC*jetEnergySmearFactor;
	JetEta[i] = jetEta[i];
	JetPhi[i] = jetPhi[i];
	JetIDTight[i] = jetPassIDTight[i];
	//JetEMFraction[i] = jetEMFraction[i];
	JetPartonFlavor[i] = jetPartonFlavor[i];
	JetPileupID[i] = jetPileupId[i];

	//Match To GenJet
	//cout << "Jet: " << jetPt[i]*JEC*jetEnergySmearFactor << " " << jetEta[i] << " " << jetPhi[i] << " : " << jetPt[i]*JEC << " " << jetPt[i] << " \n";
	GenJetE[i] = -999;
	GenJetPt[i] = -999;
	GenJetEta[i] = -999;
	GenJetPhi[i] = -999;
	double minDRToGenJet = 9999;	
	for(int j = 0; j < nGenJets; j++){
	  double DR = deltaR( genJetEta[j], genJetPhi[j], jetEta[i], jetPhi[i]);
	  if (DR > 0.4) continue;
	  //cout << "genJet: " << genJetPt[j] << " " << genJetEta[j] << " " << genJetPhi[j] << " : " << DR << "\n";
	  if (DR < minDRToGenJet) {
	    minDRToGenJet = DR;
	    GenJetE[i] = genJetE[j];
	    GenJetPt[i] = genJetPt[j];
	    GenJetEta[i] = genJetEta[j];
	    GenJetPhi[i] = genJetPhi[j];
	  }     
	}

	
	// if ( jetPt[i]*JEC*jetEnergySmearFactor > 50 && GenJetPt[i] > 0 && (jetPt[i]*JEC*jetEnergySmearFactor - GenJetPt[i])/GenJetPt[i] > 1) {
	//   cout << "DEBUG\n";
	//   cout << "GenJet: " << GenJetPt[i] << " " << GenJetEta[i] << " " << GenJetPhi[i] << " : " << (jetPt[i]*JEC*jetEnergySmearFactor - GenJetPt[i])/GenJetPt[i] << " " << minDRToGenJet << "\n";

	//   cout << "Gen Particles\n";
	//   for(int j = 0; j < nGenParticle; j++){
	//     if (deltaR(gParticleEta[j] , gParticlePhi[j], jetEta[i], jetPhi[i]) < 0.4) {
	//       cout << "particle: " << gParticleId[j] << " " << gParticleStatus[j] << " : " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << "\n";
	//     }
	//   }
	//   cout << "\n\n";

	// }

	//Check for neutrinos near jet
	PtNeutrinoClosestToJet[i] = -999;
	double minDRNeutrino = 9999;
	for(int j = 0; j < nGenParticle; j++){
	  if (!(abs(gParticleId[j]) == 12 || abs(gParticleId[j]) == 14 || abs(gParticleId[j]) == 16)) continue;
	  double DR = deltaR( jetEta[i], jetPhi[i], gParticleEta[j], gParticlePhi[j] );	  
	  if ( DR < minDRNeutrino ) {
	    minDRNeutrino = DR;
	    PtNeutrinoClosestToJet[i] = gParticlePt[j];
	  }
	}
	DRNeutrinoClosestToJet[i] = minDRNeutrino;
	

	if (jetPt[i]*JEC > 20 && isCSVL(i)) nBJetsLoose20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVM(i)) nBJetsMedium20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVT(i)) nBJetsTight20GeV++;
	
	if (jetPt[i]*JEC*jetEnergySmearFactor > 40 && fabs(jetEta[i]) < 3) {
	  numJetsAbove40GeV++;
	  if(jetPt[i]*JEC*jetEnergySmearFactor > 80) numJetsAbove80GeV++;	  
	  GoodJets.push_back(thisJet);	  
	}

      } //loop over jets


      ///sort good jets
      TLorentzVector tmpjet;
      for (int i=0;i<int(GoodJets.size());i++) {
	for (int j=0;j<int(GoodJets.size()-1);j++) {
	  if (GoodJets[j+1].Pt() > GoodJets[j].Pt()) {
	    tmpjet = GoodJets[j]; //swap elements
	    GoodJets[j] = GoodJets[j+1];
	    GoodJets[j+1] = tmpjet;
	  }
	}
      }
  
      //compute minDPhi variables
      double JER = 0.1; //average jet energy resolution
      for (int i=0;i<int(GoodJets.size());i++) {
	if (i>2) continue; //only consider first 3 jets for minDPhi 	  
	double dPhi = fmin( fabs(deltaPhi(metPhi,GoodJets[i].Phi())) , fabs( 3.1415926 - fabs(deltaPhi(metPhi,GoodJets[i].Phi()))));
	if (dPhi < minDPhi) minDPhi = dPhi;
	double deltaT = 0;
	for(auto& jet2 : GoodJets) {
	  deltaT += pow(JER*jet2.Pt()*sin(fabs(deltaPhi(GoodJets[i].Phi(),jet2.Phi()))),2);
	}
	double dPhiN = dPhi / atan2( sqrt(deltaT) , metPt);
	if (dPhiN < minDPhiN) minDPhiN = dPhiN;
      }


      //Compute the razor variables using the selected jets and possibly leptons
      vector<TLorentzVector> GoodPFObjects;
      for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);

      double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
      double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;

      TLorentzVector PFMET; PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));
      TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

      if (printSyncDebug) {
	cout << "UnCorrectedMET: " << PFMETUnCorr.Pt() << " " << PFMETUnCorr.Phi() << "\n";
	cout << "Corrected PFMET: " << PFMET.Pt() << " " << PFMET.Phi() << " | X,Y Correction :  " << MetX_Type1Corr << " " << MetY_Type1Corr << "\n";
      }



      //*************************************************************************
      //Make hemispheres
      //*************************************************************************	
      MR = 0;
      Rsq = 0;
	
      //only compute razor variables if we have 2 jets above 80 GeV
      TLorentzVector hemisphereVector0;
      TLorentzVector hemisphereVector1;
      if (GoodPFObjects.size() >= 2 && GoodJets.size() < 20
	  //&& numJetsAbove80GeV >= 2 //Si: I think we don't need this requirement at this point
	  ) {
	vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	MR = computeMR(hemispheres[0], hemispheres[1]); 
	Rsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
	hemisphereVector0 = hemispheres[0];
	hemisphereVector1 = hemispheres[1];
      }
         
      MET = PFMET.Pt();
      NJets40 = numJetsAbove40GeV;
      NJets80 = numJetsAbove80GeV;
      NBJetsLoose = nBJetsLoose20GeV;
      NBJetsMedium = nBJetsMedium20GeV;
      NBJetsTight = nBJetsTight20GeV;
      
      HT = 0;
      for(auto& pfobj : GoodPFObjects) HT += pfobj.Pt();


  
      //*************************************************************************
      //Make GenJet vector for hemispheres
      //*************************************************************************
      double GenMetX = genMetPt*cos(genMetPhi);
      double GenMetY = genMetPt*sin(genMetPhi);
      TLorentzVector GenMET; GenMET.SetPxPyPzE(GenMetX, GenMetY, 0, sqrt(GenMetX*GenMetX + GenMetY*GenMetY));

      genJetMR = 0;
      genJetRsq = 0;
      genJetDPhiRazor = 0;
      genJetHT = 0;
      vector<TLorentzVector> GenJetObjects;
      for(int j = 0; j < nGenJets; j++){
	if (genJetPt[j] > 40 && fabs(genJetEta[j]) < 3) {
	  TLorentzVector thisGenJet = makeTLorentzVector(genJetPt[j], genJetEta[j], genJetPhi[j], genJetE[j]);
	  GenJetObjects.push_back(thisGenJet);
	  genJetHT += genJetPt[j];
	}
      }
      if (GenJetObjects.size() >= 2 ) {
	vector<TLorentzVector> tmpHemispheres = getHemispheres(GenJetObjects);
	genJetMR = computeMR(tmpHemispheres[0], tmpHemispheres[1]); 
	genJetRsq = computeRsq(tmpHemispheres[0], tmpHemispheres[1], GenMET);
	genJetDPhiRazor = deltaPhi(tmpHemispheres[0].Phi(),tmpHemispheres[1].Phi());
      }
      
     
      //*************************************************************************
      //save HLT Decisions
      //*************************************************************************
      for(int k=0; k<100; ++k) {
	HLTDecision[k] = HLTDecision[k];
      }
      
      //*************************************************************************
      //Skimming
      //*************************************************************************
      bool passSkim = true;
      if (option == 1) {
      	if (!(MR>300 && Rsq>0.05)) {
      	  passSkim = false;
      	}
      }

      //*************************************************************************
      //Fill Tree
      //*************************************************************************
       if (passSkim) {
	outTree->Fill();
      }
      

      //*************************************************************
      //DEBUG
      //*************************************************************
      if (MR > 4000) printSyncDebug = true;

      if (printSyncDebug) {
	cout << "\n****************************************************************\n";
	cout << "Debug Event : " << runNum << " " << lumiNum << " " << eventNum << "\n";
	cout << "PFObjs: " << GoodPFObjects.size() << "\n";
	cout << "MR,Rsq: " << MR << " " << Rsq << " " << dPhiRazor << "\n";
	cout << "HT,MET: " << HT << " " << MET << "\n";
	cout << "\n";

	cout << "Objects: " << "\n";
	for (int i=0; i < int(GoodPFObjects.size()) ; ++i) {
	  cout << "Obj " << i << " : " << GoodPFObjects[i].Pt() << " " << GoodPFObjects[i].Eta() << " " << GoodPFObjects[i].Phi() << " \n";
	}
	cout << "\n";

	cout << "hemispheres\n";
	cout << hemisphereVector0.Pt() << " " << hemisphereVector0.Eta() << " " << hemisphereVector0.Phi() << " " << hemisphereVector0.M() << " \n";
	cout << hemisphereVector1.Pt() << " " << hemisphereVector1.Eta() << " " << hemisphereVector1.Phi() << " " << hemisphereVector1.M() << " \n";
	cout << "\n";

	cout << "NJets: " << nJets << "\n";
	for(int i = 0; i < nJets; i++){

	  //*******************************************************
	  //Correct Jet Energy Scale and Resolution
	  //*******************************************************
	  double tmpRho = fixedGridRhoFastjetAll;
	  if (isRunOne) tmpRho = fixedGridRhoAll;
	  double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						 tmpRho, jetJetArea[i], 
						 JetCorrector);   
	  Rho = tmpRho;

	  double jetEnergySmearFactor = 1.0;
	  if (!isData) {
	    std::vector<float> fJetEta, fJetPtNPU;
	    fJetEta.push_back(jetEta[i]);  
	    fJetPtNPU.push_back(jetPt[i]*JEC); 
	    fJetPtNPU.push_back(NPU_0); 
		  
	    jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], NPU_0, JetResolutionCalculator, random);
	  }
	
	  //Match To GenJet
	  cout << "Jet: " << jetPt[i]*JEC*jetEnergySmearFactor << " " << jetEta[i] << " " << jetPhi[i] << " : " << jetPt[i]*JEC << " " << jetPt[i] << " \n";

	  double minDRToGenJet = 9999;	
	  for(int j = 0; j < nGenJets; j++){
	    double DR = deltaR( genJetEta[j], genJetPhi[j], jetEta[i], jetPhi[i]);
	    if (DR > 0.4) continue;
	    cout << "found genJet: " << genJetPt[j] << " " << genJetEta[j] << " " << genJetPhi[j] << " : " << DR << "\n";
	    if (DR < minDRToGenJet) {
	      minDRToGenJet = DR;
	      GenJetE[i] = genJetE[j];
	      GenJetPt[i] = genJetPt[j];
	      GenJetEta[i] = genJetEta[j];
	      GenJetPhi[i] = genJetPhi[j];
	    }     
	  }
	
	  cout << "Closest GenJet: " << GenJetPt[i] << " " << GenJetEta[i] << " " << GenJetPhi[i] << " : " << (jetPt[i]*JEC*jetEnergySmearFactor - GenJetPt[i])/GenJetPt[i] << " " << minDRToGenJet << "\n";
	
	  // cout << "Gen Particles\n";
	  // for(int j = 0; j < nGenParticle; j++){
	  //   if (deltaR(gParticleEta[j] , gParticlePhi[j], jetEta[i], jetPhi[i]) < 0.4) {
	  //     cout << "particle: " << gParticleId[j] << " " << gParticleStatus[j] << " : " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << "\n";
	  //   }
	  // }
	  // cout << "\n\n";
	} //loop over jets

      }





    }//end of event loop
    

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

