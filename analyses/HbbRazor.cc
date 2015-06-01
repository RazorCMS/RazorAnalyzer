#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

void RazorAnalyzer::HbbRazor(string outFileName, bool combineTrees, bool isData, bool isRunOne)
{
  //initialization: create one TTree for each analysis box 
  cout << "Initializing..." << endl;
  bool printdebug = false;

  //Pileup Weights
  TFile *pileupWeightFile = 0;
  TH1D *pileupWeightHist = 0;
  if (isRunOne) {
    pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
    pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run1");
    assert(pileupWeightHist);
  }

  //Lepton Efficiency Correction Factors
  TH2D *eleLooseEffSFHist = 0;
  TH2D *eleTightEffSFHist = 0;
  if (isRunOne) {
    TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/Run1/ElectronSelection_Run2012ReReco_53X.root","READ");
    eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
    assert(eleLooseEffSFHist);
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("sfTIGHT");
    assert(eleTightEffSFHist);
  }

  if (outFileName.empty()){
    cout << "HbbRazor: Output filename not specified!" << endl << "Using default output name HbbRazor.root" << endl;
    outFileName = "HbbRazor.root";
  }
  TFile outFile(outFileName.c_str(), "RECREATE");
    
  //one tree to hold all events
  TTree *razorTree = new TTree("HbbRazor", "Info on selected razor inclusive events");
    
  //initialize jet energy corrections
  TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333
  std::vector<JetCorrectorParameters> correctionParameters;
  //get correct directory for JEC files (different for lxplus and t3-higgs)
  struct stat sb;
  string dir;
  if(stat("/afs/cern.ch/work/s/sixie/public", &sb) == 0 && S_ISDIR(sb.st_mode)){ //check if Si's directory exists
    if (isRunOne) {
      dir = "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data";
    } else {
      dir = "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data";
    }
    cout << "Getting JEC parameters from " << dir << endl;
  }
  else{ //we are on t3-higgs (for running locally on your laptop we need a separate solution)
    dir = Form("%s/src/RazorAnalyzer/data/", getenv("CMSSW_BASE"));
    cout << "Getting JEC parameters from " << dir << endl;
  }

  if (isRunOne) {
    if (isData) {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L1FastJet_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2Relative_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L3Absolute_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2L3Residual_AK5PF.txt", dir.c_str())));
    } else {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L1FastJet_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L2Relative_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L3Absolute_AK5PF.txt", dir.c_str())));
    }
  } else {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt", dir.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L2Relative_AK4PFchs.txt", dir.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt", dir.c_str())));
  }
  
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
  JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",dir.c_str()));
  SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  //tree variables
  int nSelectedJets, nBTaggedJets;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  int nVetoMuons, nVetoElectrons, nLooseTaus;
  float dPhiRazor;
  float theMR;
  float theRsq;  
  float met;
  float HT;
  float weight = 1.0;
  float pileupWeight = 1.0;
  float lepEffCorrFactor = 1.0;
  //float lepTrigCorrFactor = 1.0;
  float btagCorrFactor = 1.0;
  float b1pt = 0;
  float b1eta = 0;
  float b1phi = 0;
  float b2pt = 0;
  float b2eta = 0;
  float b2phi = 0;
  float mbb = 0;
  float ptbb = 0;
  float etabb = 0;      

  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("mbb", &mbb, "mbb/F");
    razorTree->Branch("ptbb", &ptbb, "ptbb/F");
    razorTree->Branch("etabb", &etabb, "etabb/F");
    razorTree->Branch("b1pt", &b1pt, "b1pt/F");
    razorTree->Branch("b1eta", &b1eta, "b1eta/F");
    razorTree->Branch("b1phi", &b1phi, "b1phi/F");
    razorTree->Branch("b2pt", &b2pt, "b2pt/F");
    razorTree->Branch("b2eta", &b2eta, "b2eta/F");
    razorTree->Branch("b2phi", &b2phi, "b2phi/F");


    //MET Filters
    // razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
    // razorTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
    // razorTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
    // razorTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    // razorTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
    // razorTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
    // razorTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
    // razorTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
    // razorTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
    // razorTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
    // razorTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
    // razorTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
    // razorTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
    //razorTree->Branch("HLTDecision", &HLTDecision, "HLTDecision[100]/O");

    if (!isData) {    
      razorTree->Branch("weight", &weight, "weight/F");
      //   razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
      //   razorTree->Branch("lepEffCorrFactor", &lepEffCorrFactor, "lepEffCorrFactor/F");
      //   razorTree->Branch("lepTrigCorrFactor", &lepTrigCorrFactor, "lepTrigCorrFactor/F");
      //   razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
    } else {
      razorTree->Branch("run", &runNum, "run/i");
      razorTree->Branch("lumi", &lumiNum, "lumi/i");
    }
  }

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
    printdebug = false;

    //fill normalization histogram
    NEvents->Fill(1.0);

    //reset tree variables
    nSelectedJets = 0;
    nBTaggedJets = 0;
    nVetoMuons = 0;
    nLooseMuons = 0;
    nTightMuons = 0;
    nVetoElectrons = 0;
    nLooseElectrons = 0;
    nTightElectrons = 0;
    nLooseTaus = 0;
    nTightTaus = 0;
    theMR = -1;
    theRsq = -1;
    weight = 1.0;

    //*****************************************
    //TODO: triggers!
    //*****************************************
    bool passedDileptonTrigger = false;
    bool passedSingleLeptonTrigger = false;
    bool passedLeptonicTrigger = false;
    bool passedHadronicTrigger= false;
    if (isRunOne) {
      if (HLTDecision[46] || HLTDecision[47] ||HLTDecision[48] ||HLTDecision[49] ||HLTDecision[50]) passedHadronicTrigger = true;
      if (HLTDecision[3] || HLTDecision[4] || HLTDecision[6] || HLTDecision[7]|| HLTDecision[12]) passedDileptonTrigger = true;
      if (isData) {
	if (HLTDecision[0] || HLTDecision[1] ||HLTDecision[8] ||HLTDecision[9]) passedSingleLeptonTrigger = true;
      } else {
	if (HLTDecision[0] || HLTDecision[1] || HLTDecision[9]) passedSingleLeptonTrigger = true;
      }
      passedLeptonicTrigger = passedSingleLeptonTrigger || passedDileptonTrigger;
      if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
    } else {
      passedDileptonTrigger = true;
      passedSingleLeptonTrigger = true;
      passedLeptonicTrigger = passedSingleLeptonTrigger || passedDileptonTrigger;
      passedHadronicTrigger = true;
    }
        
    //*****************************************
    //Get Pileup Information
    //*****************************************
    double NPU = 0;
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	NPU = nPUmean[i];
      }
      if (BunchXing[i] == -1) {
	NPU = nPUmean[i];
      }
      if (BunchXing[i] == 1) {
	NPU = nPUmean[i];
      }	  
    }

    if (isRunOne) {
      pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }


    //*****************************************
    //Select Leptons
    //*****************************************
    lepEffCorrFactor  = 1.0;

    vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
    for(int i = 0; i < nMuons; i++){

      if(muonPt[i] < 5) continue;
      if(abs(muonEta[i]) > 2.4) continue;

      //Calculate MC->Data Scale Factors
      if (RazorAnalyzer::matchesGenMuon(muonEta[i],muonPhi[i])) {	
	//apply muon efficiency scale factors
      }

      if(isVetoMuon(i)) nVetoMuons++;
      if(isLooseMuon(i) && muonPt[i] >= 10 ) nLooseMuons++;
      if(isTightMuon(i) && muonPt[i] >= 10) nTightMuons++;

      if(!isVetoMuon(i)) continue;  
      TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
      GoodLeptons.push_back(thisMuon);           
    }
    for(int i = 0; i < nElectrons; i++){
      if(elePt[i] < 5) continue;
      if(fabs(eleEta[i]) > 2.5) continue;

      if (isRunOne) {
	if (RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i])) {
	  double effLoose = getElectronEfficiencyRunOne("loose",elePt[i],eleEta[i]);	  
	  double effLooseSF = eleLooseEffSFHist->GetBinContent( eleLooseEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
								eleLooseEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
	  double tmpLooseSF = 1.0;
	  if( (!isRunOne && isLooseElectron(i)) || (isRunOne && isRunOneLooseElectron(i)) ) {
	    tmpLooseSF = effLooseSF;
	  } else {
	    tmpLooseSF = ( 1/effLoose - effLooseSF) / ( 1/effLoose - 1);
	  }

	  if (tmpLooseSF != tmpLooseSF) cout << tmpLooseSF << " " << effLoose << " " << effLooseSF << "\n";

	  double effTight = getElectronEfficiencyRunOne("tight",elePt[i],eleEta[i]);	  
	  double effTightSF = eleTightEffSFHist->GetBinContent( eleTightEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
								eleTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
	  double tmpTightSF = 1.0;
	  if( (!isRunOne && isTightElectron(i)) || (isRunOne && isRunOneTightElectron(i)) ) {
	    tmpTightSF = effTightSF;
	  } else {
	    tmpTightSF = ( 1/effTight - effTightSF) / ( 1/effTight - 1);
	  }
	  if (tmpTightSF != tmpTightSF) cout << tmpTightSF << " " << effTight << " " << effTightSF << "\n";

	  lepEffCorrFactor *= tmpLooseSF;
	  lepEffCorrFactor *= tmpTightSF;
	}
      }

      if(isMVANonTrigVetoElectron(i)) nVetoElectrons++;
      if( ( (!isRunOne && isLooseElectron(i)) || (isRunOne && isRunOneLooseElectron(i)) ) && elePt[i] > 10 ) nLooseElectrons++;
      if( ( (!isRunOne && isTightElectron(i)) || (isRunOne && isRunOneTightElectron(i)) ) && elePt[i] > 10 ) nTightElectrons++;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : GoodLeptons){
	if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
      }
      if(overlap) continue;

      if(!isMVANonTrigVetoElectron(i)) continue; 
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
      GoodLeptons.push_back(thisElectron);            
    }

    //******************************
    //Only Do Taus for Run2
    //******************************
    if (!isRunOne) {
      for(int i = 0; i < nTaus; i++){	 
	if (tauPt[i] < 20) continue;
	if (fabs(tauEta[i]) > 2.4) continue;
	
	if(isLooseTau(i)) nLooseTaus++;
	if(isTightTau(i)) nTightTaus++;
	
	//remove overlaps
	bool overlap = false;
	for(auto& lep : GoodLeptons){
	  if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
	}
	if(overlap) continue;
	
	if (!isLooseTau(i)) continue;
	TLorentzVector thisTau; thisTau.SetPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
	GoodLeptons.push_back(thisTau);  
      }
    }
        

    vector<TLorentzVector> GoodJets;
    int numJetsAbove80GeV = 0;

    // //***********************************************
    // //use genjets instead , for debugging
    // //***********************************************
    // for(int j = 0; j < nGenJets; j++){

    //   if(genJetPt[j] < 40) continue;
    //   if(fabs(genJetEta[j]) > 3.0) continue;

    //   //exclude selected muons and electrons from the jet collection
    //   double deltaR = -1;
    //   TLorentzVector thisJet = makeTLorentzVector(genJetPt[j], genJetEta[j], genJetPhi[j], genJetE[j]);
    //   for(auto& lep : GoodLeptons){
    // 	double thisDR = thisJet.DeltaR(lep);
    // 	if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
    //   }
    //   if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

    //   if(genJetPt[j] > 80) numJetsAbove80GeV++;
    //   GoodJets.push_back(thisJet);
    //   nSelectedJets++;

    //   bool isBJet = false;
    //   for(int i = 0; i < nJets; i++){
    // 	double tmpDR = RazorAnalyzer::deltaR( genJetEta[j], genJetPhi[j], jetEta[i], jetPhi[i] );
    // 	if ( tmpDR < 0.4 && abs(jetPartonFlavor[i]) == 5) isBJet = true;
    //   }
    //   if(isBJet) nBTaggedJets++;
     
    // }


    //initialize B-Tagging Correction Factor
    btagCorrFactor = 1.0;

    //***********************************************
    //Variables for Type1 Met Correction
    //***********************************************
    double MetX_Type1Corr = 0;
    double MetY_Type1Corr = 0;


    vector<TLorentzVector> BJetCandidates; //candidates for H->bb

    //***********************************************
    //Select Jets
    //***********************************************
    for(int i = 0; i < nJets; i++){


      //Don't exclude leptons from jets in this analysis
      // //*****************************************************************
      // //exclude selected muons and electrons from the jet collection
      // //*****************************************************************
      // double deltaR = -1;
      // for(auto& lep : GoodLeptons){
      // 	double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());  
      // 	if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      // }
      // if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
      

      //*****************************************************************
      //apply Jet ID
      //*****************************************************************
      if (isRunOne) {
	if (!jetPassIDTight[i]) continue;
      }

      //*****************************************************************
      //Apply Jet Energy and Resolution Corrections
      //*****************************************************************
      double tmpRho = fixedGridRhoFastjetAll;
      if (isRunOne) tmpRho = fixedGridRhoAll;
      double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
    					     tmpRho, jetJetArea[i], 
    					     JetCorrector);   

      double jetEnergySmearFactor = 1.0;
      if (!isData) {
	if (isRunOne) {
	  jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random);
	}
      }
      
      TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
      TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);      
      double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
      //double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;

      // cout << "Jet " << i << " : " << jetCorrPt << " " << jetEta[i] << " " << jetPhi[i] << " : " <<jetPartonFlavor[i] 
      // 	   << "\n";
      
      //*******************************
      //B-Tagging Correction Factor
      //*******************************
      if (isRunOne) {
	if (abs(jetPartonFlavor[i]) == 5 &&jetCorrPt > 20) {
	  double tmpBTagCorrFactor = 1.0;
	
	  double tmpCorrFactor = 0.938887 + 0.00017124 * jetCorrPt + (-2.76366e-07) * jetCorrPt * jetCorrPt ;
	  double MCEff = 1.0;
	  if (jetCorrPt < 50) MCEff = 0.65;
	  else if (jetCorrPt < 80) MCEff = 0.70;
	  else if (jetCorrPt < 120) MCEff = 0.73;
	  else if (jetCorrPt < 210) MCEff = 0.73;
	  else MCEff = 0.66;				 
	
	  //if pass CSV Medium
	  if((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) {
	    tmpBTagCorrFactor = tmpCorrFactor;
	  } else {
	    tmpBTagCorrFactor = ( 1/MCEff - tmpCorrFactor) / ( 1/MCEff - 1);
	  }

	  btagCorrFactor *= tmpBTagCorrFactor;
	}
      }

      //*******************************
      //Add to Type1 Met Correction
      //*******************************
      if (isRunOne) {
	if (jetPt[i]*JEC*jetEnergySmearFactor > 20) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - UnCorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - UnCorrJet.Py()  );
	}
      }

      //*******************************************************
      //apply  Pileup Jet ID
      //*******************************************************
      if (isRunOne) {
	int level = 2; //loose jet ID
	if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;
      }
      
      //*******************************************************
      //apply Jet cuts
      //*******************************************************
      if(jetCorrPt < 40) continue;
      if(fabs(jetEta[i]) > 2.5) continue;
            
      if(jetCorrPt > 80) numJetsAbove80GeV++;
      GoodJets.push_back(thisJet);
      nSelectedJets++;

      if((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))){ 
    	nBTaggedJets++;
	TLorentzVector thisBJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);       
	BJetCandidates.push_back(thisBJet);
      }
    }


    //*************************************************************
    //find Higgs->bb candidates
    //*************************************************************
    // for (int i=0; i<BJetCandidates.size(); ++i) {
    //   cout << "BJet " << i << " : " << BJetCandidates[i].Pt() << " " << BJetCandidates[i].Eta()  << " " << BJetCandidates[i].Phi() 
    // 	   << "\n";
    // }

    b1pt = 0;
    b1eta = 0;
    b1phi = 0;
    b2pt = 0;
    b2eta = 0;
    b2phi = 0;
    mbb = 0;
    ptbb = 0;
    etabb = 0;
    bool foundHiggsCandidate = false;
    pair<TLorentzVector,TLorentzVector> HiggsToBBCandidate;
    double highestHiggsPt = 0;
    for (int i=0; i<int(BJetCandidates.size()); ++i) {
      for (int j=i+1; j< int (BJetCandidates.size()); ++j) {
	//cout << "BJet Pair " << i << " " << j << " : " << (BJetCandidates[i]+BJetCandidates[j]).Pt() << "\n";
	if ( (BJetCandidates[i]+BJetCandidates[j]).Pt() > highestHiggsPt) {
	  //cout << "highest pt\n";
	  foundHiggsCandidate = true;
	  highestHiggsPt = (BJetCandidates[i]+BJetCandidates[j]).Pt();
	  HiggsToBBCandidate.first = BJetCandidates[i];
	  HiggsToBBCandidate.second = BJetCandidates[j];
	}
      }
    }

    TLorentzVector HiggsCandidate;
    if (foundHiggsCandidate) {
      HiggsCandidate = HiggsToBBCandidate.first + HiggsToBBCandidate.second;
      mbb = HiggsCandidate.M();
      ptbb = HiggsCandidate.Pt();
      etabb = HiggsCandidate.Eta();
      b1pt = HiggsToBBCandidate.first.Pt();
      b1eta = HiggsToBBCandidate.first.Eta();
      b1phi = HiggsToBBCandidate.first.Phi();
      b2pt = HiggsToBBCandidate.second.Pt();
      b2eta = HiggsToBBCandidate.second.Eta();
      b2phi = HiggsToBBCandidate.second.Phi();
    }

    //if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

    //Compute the razor variables using the selected jets and possibly leptons
    vector<TLorentzVector> GoodPFObjects;
    if (foundHiggsCandidate) {
      GoodPFObjects.push_back(HiggsCandidate);
    }
    for(auto& jet : GoodJets) {
      if (foundHiggsCandidate) {
       	if ( (jet.Pt() == HiggsToBBCandidate.first.Pt() && jet.Eta() == HiggsToBBCandidate.first.Eta() && jet.Phi() == HiggsToBBCandidate.first.Phi())
	     || 
	     (jet.Pt() == HiggsToBBCandidate.second.Pt() && jet.Eta() == HiggsToBBCandidate.second.Eta() && jet.Phi() == HiggsToBBCandidate.second.Phi())
	     ) {
	  continue;
	}
      }
      GoodPFObjects.push_back(jet);
    }


    //*************************************************************
    //Apply Type1 Met Correction
    //*************************************************************
    double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
    double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;
    TLorentzVector PFMET; 
    if (isRunOne) {
      PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));      
    } else {
      PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));      
    }
    TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

    HT = 0;
    for(auto& obj : GoodPFObjects) HT += obj.Pt();


    vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
    theMR = computeMR(hemispheres[0], hemispheres[1]); 
    theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
    dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
    met = metPt;


    //**********************************************************************
    //Apply ECAL Dead Cells Filter
    //**********************************************************************
    if (isRunOne) {
      if (isData) {
	if (!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) {
	  cout << "Fail noise filter\n";
	  continue;
	}
      }
    } else {
      if (Flag_EcalDeadCellTriggerPrimitiveFilter == false) continue;
    }

    //**********************************************************************
    //Compute correction factor weight
    //**********************************************************************
    if (isRunOne) {
      weight *= pileupWeight;
      weight *= lepEffCorrFactor;
      weight *= btagCorrFactor;    
    }


    //**********************************************************************
    //Categorize Events into Razor Boxes 
    //**********************************************************************

    if (foundHiggsCandidate 
	&& mbb > 0 && theMR > 0 && theRsq > 0.0
	) {
      razorTree->Fill();
    }


    //******************************
    //Print Debug
    //******************************
    if (printdebug) {
      cout << "\nNew Event\n";
      for(int j = 0; j < nGenParticle; j++){
	cout << "GenParticle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";
      }
    }



  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  if(combineTrees) razorTree->Write();
  NEvents->Write();

  outFile.Close();
}
