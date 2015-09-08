
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
 
void RazorAnalyzer::RazorControlRegions( string outputfilename, int option, bool isData)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    cout << "IsData = " << isData << "\n";

    // TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    bool printSyncDebug = false;
    std::vector<JetCorrectorParameters> correctionParameters;
    char* cmsswPath;
    cmsswPath = getenv("CMSSW_BASE");
    string pathname;
    if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
    cout << "Getting JEC parameters from " << pathname << endl;

    if (isData) {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
      //correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
    } else {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    }
    
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
    JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
    SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);


    //*************************************************************************
    //Set up Output File
    //*************************************************************************
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorControlRegions.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ControlSampleEvents *events = new ControlSampleEvents;
    
    //**********************************************
    //Options
    //tens and ones digits refer to tree type.  
    //1: Single-Lepton
    //2: Single-Lepton Add To MET
    //3: Dilepton
    //4: Dilepton Add To MET
    //5: Photon Add To MET
    //6: Zero Lepton
    //7: Single Veto-Lepton
    //8: Loose Lepton + Veto-Lepton
    //9: Single Tau
    //11: Single-Lepton Reduced
    //12: Single-Lepton Add To MET Reduced
    //13: Dilepton Reduced
    //14: Dilepton Add To MET Reduced
    //15: Photon Reduced
    //17: Single Veto-Lepton Reduced
    //18: Loose Lepton + Veto-Lepton Reduced
    //19: Single Tau Reduced
    //hundreds digit refers to lepton skim option
    //0: no skim
    //1: 1 lepton skim
    //2: 2 lepton skim
    //5: photon skim
    //thousand digit refers to razor skim option
    //0: no skim
    //1: razor skim MR > 300 Rsq > 0.15
    //2: razor skim MR_NoW > 300 && Rsq_NoW > 0.15
    //3: razor skim MR_NoZ > 300 && Rsq_NoZ > 0.15
    //5: razor skim MR_NoPho > 300 && Rsq_NoPho > 0.15
    //*********************************************
    int razorSkimOption = floor(float(option) / 1000);
    int leptonSkimOption = floor( float(option - razorSkimOption*1000) / 100);
    int treeTypeOption = option - razorSkimOption*10000 - leptonSkimOption*100;

    if (treeTypeOption == 1)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (treeTypeOption == 11)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Reduced);
    else if (treeTypeOption == 2)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLeptonAdd2MET_Full);
    else if (treeTypeOption == 3)
      events->CreateTree(ControlSampleEvents::kTreeType_Dilepton_Full);
    else if (treeTypeOption == 4)
      events->CreateTree(ControlSampleEvents::kTreeType_DileptonAdd2MET_Full);
    else if (treeTypeOption == 5)
      events->CreateTree(ControlSampleEvents::kTreeType_Photon_Full);
    else if (treeTypeOption == 6)
      events->CreateTree(ControlSampleEvents::kTreeType_ZeroLepton_Full);
    else if (treeTypeOption == 7)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (treeTypeOption == 8)
      events->CreateTree(ControlSampleEvents::kTreeType_Dilepton_Full);
    else {
      events->CreateTree(ControlSampleEvents::kTreeType_Default);
    }
    events->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    int NEventProcessed = 0;
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  

    //*************************************************************************
    //Look over Input File Events
    //*************************************************************************
    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
      //initialize all varabiles
      events->InitVariables();

      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      printSyncDebug = false;
      if (printSyncDebug) {
	cout << "\n****************************************************************\n";
	cout << "Debug Event : " << runNum << " " << lumiNum << " " << eventNum << "\n";
      }


      //fill normalization histogram
      NEventProcessed++;
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      
      //event info
      events->option = option;
      events->weight = genWeight;
      events->run = runNum;
      events->lumi = lumiNum;
      events->event = eventNum;
      events->processID = 0;
 
      //get NPU
      for (int i=0; i < nBunchXing; ++i) {
	if (BunchXing[i] == 0) {
	  events->NPU_0 = nPUmean[i];
	}
	if (BunchXing[i] == -1) {
	  events->NPU_Minus1 = nPUmean[i];
	}
	if (BunchXing[i] == 1) {
	  events->NPU_Plus1 = nPUmean[i];
	}	  
      }
      events->NPV = nPV;
 

        
      //******************************************
      //Find Generated leptons
      //******************************************
      vector<int> genLeptonIndex;

      for(int j = 0; j < nGenParticle; j++){
	//look for electrons
	if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	    && abs(gParticleEta[j]) < 2.5 && gParticlePt[j] > 5
	    ) {
	  if ( abs(gParticleMotherId[j]) == 24 
	       || abs(gParticleMotherId[j]) == 23 
	       || ( (abs(gParticleMotherId[j]) == 15 || abs(gParticleMotherId[j]) == 13 || abs(gParticleMotherId[j]) == 11 )
		    && gParticleMotherIndex[j] >= 0 
		    && (abs(gParticleMotherId[gParticleMotherIndex[j]]) == 24 || 
			abs(gParticleMotherId[gParticleMotherIndex[j]]) == 23)
		    )
	       )  {	      
	    genLeptonIndex.push_back(j);	      
	  }	    
	}

	//look for muons
	if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	    && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 5
	    ) {
	  if ( abs(gParticleMotherId[j]) == 24 
	       || abs(gParticleMotherId[j]) == 23 
	       || ( (abs(gParticleMotherId[j]) == 15 || abs(gParticleMotherId[j]) == 13 || abs(gParticleMotherId[j]) == 11 )
		    && gParticleMotherIndex[j] >= 0 
		    && (abs(gParticleMotherId[gParticleMotherIndex[j]]) == 24 || 
			abs(gParticleMotherId[gParticleMotherIndex[j]]) == 23)
		    )
	       ) {	     
	    genLeptonIndex.push_back(j);
	  }	    
	}

	//look for hadronic taus
	if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
	    && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
	    && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 20
	    ) {
	  bool isLeptonicTau = false;
	  for(int k = 0; k < nGenParticle; k++){
	    if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == j) {
	      isLeptonicTau = true;
	      break;
	    }

	    //handle weird case
	    //the status 2 tau has a status 23 tau as mother. the e or mu have the status 23 tau as its mother
	    if ( gParticleMotherIndex[j] >= 0 && gParticleId[gParticleMotherIndex[j]] == gParticleId[j] 
		 && gParticleStatus[gParticleMotherIndex[j]] == 23
		 && (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13)
		 && gParticleMotherIndex[k] == gParticleMotherIndex[j]) {
	      isLeptonicTau = true;
	      break;
	    }

	  }
	  if (!isLeptonicTau) {
	    genLeptonIndex.push_back(j);
	  }	   
	}

      } //loop over gen particles
	

	//sort gen leptons by pt
      int tempIndex = -1;
      for(uint i = 0; i < genLeptonIndex.size() ; i++) {
	for (uint j=0; j < genLeptonIndex.size()-1; j++) {
	  if (gParticlePt[genLeptonIndex[j+1]] > gParticlePt[genLeptonIndex[j]]) { 
	    tempIndex = genLeptonIndex[j];             // swap elements
	    genLeptonIndex[j] = genLeptonIndex[j+1];
	    genLeptonIndex[j+1] = tempIndex;
	  }
	}
      }

      events->genlep1.SetPtEtaPhiM(0,0,0,0);
      events->genlep2.SetPtEtaPhiM(0,0,0,0);
      events->genlep1Type=0;
      events->genlep2Type=0;
      for (uint i=0;i<genLeptonIndex.size();i++) {
	if(i==0) {
	  double mass = 0.000511;
	  if (abs(gParticleId[genLeptonIndex[i]]) == 13) mass = 0.1057;
	  if (abs(gParticleId[genLeptonIndex[i]]) == 15) mass = 1.777;	    
	  events->genlep1.SetPtEtaPhiM(gParticlePt[genLeptonIndex[i]],gParticleEta[genLeptonIndex[i]],gParticlePhi[genLeptonIndex[i]],mass);
	  events->genlep1Type = gParticleId[genLeptonIndex[i]];
	}
	if(i==1) {
	  double mass = 0.000511;
	  if (abs(gParticleId[genLeptonIndex[i]]) == 13) mass = 0.1057;
	  if (abs(gParticleId[genLeptonIndex[i]]) == 15) mass = 1.777;	    
	  events->genlep2.SetPtEtaPhiM(gParticlePt[genLeptonIndex[i]],gParticleEta[genLeptonIndex[i]],gParticlePhi[genLeptonIndex[i]],mass);
	  events->genlep2Type = gParticleId[genLeptonIndex[i]];
	}				      
      }




      //*************************************************************************
      //Find Reconstructed Leptons
      //*************************************************************************
	
      vector<int> VetoLeptonIndex; 
      vector<int> VetoLeptonType;
      vector<int> VetoLeptonPt;
      vector<int> LooseLeptonIndex; 
      vector<int> LooseLeptonType;
      vector<int> LooseLeptonPt;
      vector<int> TightLeptonIndex; 
      vector<int> TightLeptonType;
      vector<int> TightLeptonPt;
      vector<TLorentzVector> GoodLeptons;//leptons used to compute hemispheres
      vector<int> GoodLeptonType;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsTight;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsLoose;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsVeto;//leptons used to compute hemispheres



      for(int i = 0; i < nMuons; i++){

	if(muonPt[i] < 5) continue;
	if(fabs(muonEta[i]) > 2.4) continue;

	//don't count duplicate muons 
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), muonEta[i],muonPhi[i]) < 0.1) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if(isTightMuon(i) && muonPt[i] >= 10) {
	  TightLeptonType.push_back(13 * -1 * muonCharge[i]);
	  TightLeptonIndex.push_back(i);
	  TightLeptonPt.push_back(muonPt[i]);
	}
	else if(isLooseMuon(i) && muonPt[i] >= 10 ) {
	  LooseLeptonType.push_back(13 * -1 * muonCharge[i]);
	  LooseLeptonIndex.push_back(i);
	  LooseLeptonPt.push_back(muonPt[i]);
	}
	else if(isVetoMuon(i)) {
	  VetoLeptonType.push_back(13 * -1 * muonCharge[i]);
	  VetoLeptonIndex.push_back(i);
	  VetoLeptonPt.push_back(muonPt[i]);
	}

	if (printSyncDebug) cout << "muon " << i << " " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : Tight = " << isTightMuon(i) << " Loose = " << isLooseMuon(i) << " Veto = " << isVetoMuon(i) << " \n";
                        	   
	TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	thisMuon.SetPtEtaPhiM( muonPt[i], muonEta[i], muonPhi[i], 0.1057);

	//*******************************************************
	//For Single Lepton Options, use only tight leptons
	//For Dilepton Options, use only loose leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12) {
	  if (isTightMuon(i)) isGoodLepton = true;
	}
	if (treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14) {
	  if (isLooseMuon(i)) isGoodLepton = true;
	}
	if (treeTypeOption == 7 || treeTypeOption == 8 || treeTypeOption == 17 || treeTypeOption == 18
	    || treeTypeOption == 9 || treeTypeOption == 19 ) {
	  if (isVetoMuon(i) || isLooseMuon(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisMuon);
	  GoodLeptonType.push_back(13 * -1 * muonCharge[i]);
	  GoodLeptonIsTight.push_back( isTightMuon(i) );
	  GoodLeptonIsLoose.push_back( isLooseMuon(i) );
	  GoodLeptonIsVeto.push_back( isVetoMuon(i) );
	}

      } // loop over muons


      for(int i = 0; i < nElectrons; i++){

	if(elePt[i] < 5) continue;
	if(fabs(eleEta[i]) > 2.5) continue;

	//don't count electrons that were already selected as muons
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), eleEta[i],elePhi[i]) < 0.1) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if( isTightElectron(i) && elePt[i] > 10 ) {
	  TightLeptonType.push_back(11 * -1 * eleCharge[i]);
	  TightLeptonIndex.push_back(i);
	  TightLeptonPt.push_back(elePt[i]);
	}
	else if( isLooseElectron(i) && elePt[i] > 10 ) {
	  LooseLeptonType.push_back(11 * -1 * eleCharge[i]);
	  LooseLeptonIndex.push_back(i);
	  LooseLeptonPt.push_back(elePt[i]);
	}
	else if(isVetoElectron(i)) {
	  VetoLeptonType.push_back(11 * -1 * eleCharge[i]);
	  VetoLeptonIndex.push_back(i);
	  VetoLeptonPt.push_back(elePt[i]);
	}
            
	if (printSyncDebug) cout << "ele " << i << " " << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " : Tight = " << isTightElectron(i) << " Loose = " << isLooseElectron(i) << " Veto = " << isVetoElectron(i) << " \n";

	if(!isVetoElectron(i)) continue; 

	TLorentzVector thisElectron; // = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
	thisElectron.SetPtEtaPhiM( elePt[i], eleEta[i], elePhi[i], 0.000511);

	//*******************************************************
	//For Single Lepton Options, use only tight leptons
	//For Dilepton Options, use only loose leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12) {
	  if (isTightElectron(i)) isGoodLepton = true;
	}
	if (treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14) {
	  if (isLooseElectron(i)) isGoodLepton = true;
	}
	if (treeTypeOption == 7 || treeTypeOption == 8 || treeTypeOption == 17 || treeTypeOption == 18
	    || treeTypeOption == 9 || treeTypeOption == 19 ) {
	  if (isVetoElectron(i) || isLooseElectron(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisElectron);        
	  GoodLeptonType.push_back(11 * -1 * muonCharge[i]);
	  GoodLeptonIsTight.push_back( isTightElectron(i) );
	  GoodLeptonIsLoose.push_back( isLooseElectron(i) );
	  GoodLeptonIsVeto.push_back( isVetoElectron(i) );
	}
      }

      for(int i = 0; i < nTaus; i++){
	if (tauPt[i] < 20) continue;
	if (fabs(tauEta[i]) > 2.4) continue;

	//don't count taus that were already selected as muons or electrons
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), tauEta[i],tauPhi[i]) < 0.1 ) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if(isTightTau(i)){
	  TightLeptonType.push_back(15);
	  TightLeptonIndex.push_back(i);
	  TightLeptonPt.push_back(tauPt[i]);
	}
	else if(isLooseTau(i)){
	  LooseLeptonType.push_back(15);
	  LooseLeptonIndex.push_back(i);
	  LooseLeptonPt.push_back(tauPt[i]);
	}

	TLorentzVector thisTau;
	thisTau.SetPtEtaPhiM( tauPt[i], tauEta[i], tauPhi[i], 1.777);

	//*******************************************************
	//For Single Lepton Options, use only tight leptons
	//For Dilepton Options, use only loose leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 9 || treeTypeOption == 19 ) {
	  if (isLooseTau(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisTau);        
	  GoodLeptonType.push_back(15);
	  GoodLeptonIsTight.push_back( isTightTau(i) );
	  GoodLeptonIsLoose.push_back( isLooseTau(i) );
	  GoodLeptonIsVeto.push_back( isLooseTau(i) );
	}	  
      }

      //************************************************************************
      //Sort the collections by pT
      //************************************************************************
      for(uint i = 0; i < GoodLeptons.size() ; i++) {
	for (uint j=0; j < GoodLeptons.size()-1; j++) {
	  if (GoodLeptons[j+1].Pt() > GoodLeptons[j].Pt()) { 

	    // swap elements
	    TLorentzVector tmpV = GoodLeptons[j]; 
	    int tmpType = GoodLeptonType[j];
	    bool tmpIsTight = GoodLeptonIsTight[j];
	    bool tmpIsLoose = GoodLeptonIsLoose[j];
	    bool tmpIsVeto = GoodLeptonIsVeto[j];

	    GoodLeptons[j] = GoodLeptons[j+1];
	    GoodLeptonType[j] = GoodLeptonType[j+1];
	    GoodLeptonIsTight[j] = GoodLeptonIsTight[j+1];
	    GoodLeptonIsLoose[j] = GoodLeptonIsLoose[j+1];
	    GoodLeptonIsVeto[j] = GoodLeptonIsVeto[j+1];

	    GoodLeptons[j+1] = tmpV;
	    GoodLeptonType[j+1] = tmpType;
	    GoodLeptonIsTight[j+1] = tmpIsTight;
	    GoodLeptonIsLoose[j+1] = tmpIsLoose;
	    GoodLeptonIsVeto[j+1] = tmpIsVeto;	  
	  }
	}
      }

   

      //************************************************************************
      //Fill Lepton Information using Good Leptons collection
      //************************************************************************
      events->lep1.SetPtEtaPhiM(0,0,0,0);
      events->lep2.SetPtEtaPhiM(0,0,0,0);
      events->lep1Type = 0;
      events->lep2Type = 0;
      events->lep1MatchedGenLepIndex = -1;
      events->lep2MatchedGenLepIndex = -1;
      events->lep1PassVeto = false;
      events->lep1PassLoose = false;
      events->lep1PassTight = false;
      events->lep2PassVeto = false;
      events->lep2PassLoose = false;
      events->lep2PassTight = false;

      bool lep1Found = false;
      bool lep2Found = false;
      
      //************************************************************************
      //For single lepton and dilepton Control Regions
      //************************************************************************
      if ( treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12 
	   || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
	   || treeTypeOption == 7 || treeTypeOption == 17
	   || treeTypeOption == 9 || treeTypeOption == 19
	   ) {
	for (uint i=0; i<GoodLeptons.size(); i++) {

	  //For tau control region, consider only tau leptons
	  if (treeTypeOption == 9 || treeTypeOption == 19) {
	    if ( abs(GoodLeptonType[i]) != 15) continue;
	  }

	  if (!lep1Found) {
	    events->lep1Type = GoodLeptonType[i];
	    events->lep1.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  
	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(GoodLeptonType[i]) == abs(events->genlep1Type) &&
		 deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		 ) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = GoodLeptonIsTight[i];
	    events->lep1PassLoose = GoodLeptonIsLoose[i];
	    events->lep1PassVeto = GoodLeptonIsVeto[i];
	    lep1Found = true;
	  } else {
	    if (!lep2Found) {
	      events->lep2Type = GoodLeptonType[i];
	      events->lep2.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  	  	   
	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(GoodLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = GoodLeptonIsTight[i];
	      events->lep2PassLoose = GoodLeptonIsLoose[i];
	      events->lep2PassVeto = GoodLeptonIsVeto[i];
	      lep2Found = true;
	    }
	  }
	} //loop over good leptons	

	
      } // if using 1-lepton or 2-lepton control region option
      
      
      //************************************************************************
      //For Loose Lepton + Veto lepton control region option
      //First select a Loose Lepton, assign it to lep1
      //Next select a Veto Lepton
       //************************************************************************
      if ( treeTypeOption == 8 || treeTypeOption == 18 ) {
	
	//look for a loose lepton first
	for (uint i=0; i<GoodLeptons.size(); i++) {
	  if (!GoodLeptonIsLoose[i]) continue;
	  if (!lep1Found) {
	    events->lep1Type = GoodLeptonType[i];
	    events->lep1.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  
	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(GoodLeptonType[i]) == abs(events->genlep1Type) &&
		 deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		 ) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = GoodLeptonIsTight[i];
	    events->lep1PassLoose = GoodLeptonIsLoose[i];
	    events->lep1PassVeto = GoodLeptonIsVeto[i];
	    lep1Found = true;
	    break;
	  }
	}

	//next look for a veto lepton
	if (lep1Found) {
	  for (uint i=0; i<GoodLeptons.size(); i++) {
	    if (!GoodLeptonIsVeto[i]) continue; 

	    //skip the object already saved as lep1
	    if ( deltaR( GoodLeptons[i].Eta(),GoodLeptons[i].Phi(), events->lep1.Eta(), events->lep1.Phi()) < 0.01) continue;

	    if (!lep2Found) {
	      events->lep2Type = GoodLeptonType[i];
	      events->lep2.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  	  	   
	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(GoodLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = GoodLeptonIsTight[i];
	      events->lep2PassLoose = GoodLeptonIsLoose[i];
	      events->lep2PassVeto = GoodLeptonIsVeto[i];
	      lep2Found = true;
	    }
	  }
	}
		
      } //if using loose lepton + veto lepton control region option

      //************************************************************************
      //save this for the 1-lepton mini ntuples
      //************************************************************************
      if (lep1Found) {
	events->lep1Pt = events->lep1.Pt();
	events->lep1Eta = events->lep1.Eta();
      } else {
	events->lep1Pt = -99;
	events->lep1Eta = -99;
      }
      
      //************************************************************************
      //save for dilepton ntuples
      //************************************************************************
      if (lep1Found && lep2Found) {
	events->mll = (events->lep1 + events->lep2).M();
      }
      
      if (printSyncDebug) {
	cout << "\n\n";
	cout << "lep1: " << events->lep1Type << " | " << events->lep1.Pt() << " " << events->lep1.Eta() << " " << events->lep1.Phi() << " | Tight = " << events->lep1PassTight << " Loose = " << events->lep1PassLoose << " Veto = " << events->lep1PassVeto << "\n";
	cout << "lep2: " << events->lep2Type << " | " << events->lep2.Pt() << " " << events->lep2.Eta() << " " << events->lep2.Phi() << " | Tight = " << events->lep2PassTight << " Loose = " << events->lep2PassLoose << " Veto = " << events->lep2PassVeto << "\n";	 
      }
    

      //****************************************************//
      //             Select photons                         //
      //****************************************************//
      vector<TLorentzVector> GoodPhotons;
      int nPhotonsAbove40GeV = 0;

      if (treeTypeOption == 5 || treeTypeOption == 15) {
	for(int i = 0; i < nPhotons; i++){


	  if(phoPt[i] < 10) continue;
	  if(fabs(phoEta[i]) > 2.5) continue;
	  if(!isTightPhoton(i)) continue;

	  if(phoPt[i] > 40) nPhotonsAbove40GeV++;
	  TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], phoE[i]);
	  GoodPhotons.push_back(thisPhoton);
	}
      }     
      events->nSelectedPhotons = nPhotonsAbove40GeV;
      
      //Sort Photon Collection
      sort(GoodPhotons.begin(), GoodPhotons.end(), greater_than_pt());
    
      //****************************************************//
      // Set Photon HLT Bits
      //****************************************************//
      for(int i = 62; i <= 67; i++){
	if(HLTDecision[i] == 1) events->HLT_Photon = true;
      }
	
      if(isData && events->HLT_Photon){
	//save the trigger bits
	if(HLTDecision[67] == 1){
	  events->HLT_Photon165 = true;
	}
	if(HLTDecision[66] == 1){
	  events->HLT_Photon120 = true;
	}
	if(HLTDecision[65] == 1){
	  events->HLT_Photon90 = true;
	}
	if(HLTDecision[64] == 1){
	  events->HLT_Photon75 = true;
	}
	if(HLTDecision[63]== 1){
	  events->HLT_Photon50 = true;
	}
	if(HLTDecision[62]== 1){
	  events->HLT_Photon36 = true;
	}
      }

    


      //************************************************************************
      //Select Jets
      //************************************************************************
      bool bjet1Found = false;
      bool bjet2Found = false;
      vector<TLorentzVector> GoodJets;
      int numJetsAbove80GeV = 0;
      int numJetsAbove40GeV = 0;
      double MetX_Type1Corr = 0;
      double MetY_Type1Corr = 0;
      int nBJetsLoose20GeV = 0;
      int nBJetsMedium20GeV = 0;
      int nBJetsTight20GeV = 0;
      events->bjet1.SetPtEtaPhiM(0,0,0,0);
      events->bjet2.SetPtEtaPhiM(0,0,0,0);
      events->bjet1PassLoose = false;
      events->bjet1PassMedium = false;
      events->bjet1PassTight = false;
      events->bjet2PassLoose = false;
      events->bjet2PassMedium = false;
      events->bjet2PassTight = false;       
      events->minDPhi = 9999;
      events->minDPhiN = 9999;

      float mhx = 0., mhy = 0., mhx_nohf = 0., mhy_nohf = 0.;

      for(int i = 0; i < nJets; i++){

	//exclude Good leptons from the jet collection
	//NB: Currently the code excludes ALL identified "Good" leptons from the jet collection
	//even if the event is interpreted as a single lepton or two lepton event. This may cause
	//inconsistencies in the future.
	double dR = -1;
	for(auto& lep : GoodLeptons){
	  double thisDR = deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());
	  if(dR < 0 || thisDR < dR) dR = thisDR;
	}
	if(dR > 0 && dR < 0.4) continue; //jet matches a selected lepton
	  
	  
	if (printSyncDebug)  {
	  cout << "jet " << i << " : " << jetPt[i] << " " << jetEta[i] << " " << jetPhi[i] 
	       << " : rho = " << fixedGridRhoAll << " area = " << jetJetArea[i] << " "
	       << " | " 
	       << "correctedPt = " << jetPt[i]*JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
									 fixedGridRhoAll, jetJetArea[i], 
									 JetCorrector) << " "
	       << " | passID = " << jetPassIDTight[i] << " passPUJetID = " << bool((jetPileupIdFlag[i] & (1 << 2)) != 0) 
	       << " | csv = " << jetCSV[i] << " passCSVL = " << isCSVL(i) << " passCSVM = " << isCSVM(i) << " " << "\n";
	}

	//*******************************************************
	//apply jet iD
	//*******************************************************
	if (!jetPassIDTight[i]) continue;

	//*******************************************************
	//Correct Jet Energy Scale and Resolution
	//*******************************************************
	double tmpRho = fixedGridRhoFastjetAll;

	double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i], 
					       JetCorrector);   

	double jetEnergySmearFactor = 1.0;
	if (!isData) {
	  std::vector<float> fJetEta, fJetPtNPU;
	  fJetEta.push_back(jetEta[i]);  
	  fJetPtNPU.push_back(jetPt[i]*JEC); 
	  fJetPtNPU.push_back(events->NPU_0); 
	  if (printSyncDebug) {
	    cout << "Jet Resolution : " << jetPt[i]*JEC << " " << jetEta[i] << " " << jetPhi[i] << " : " 
		 << JetResolutionCalculator->resolution(fJetEta,fJetPtNPU) << "\n";
	  }
	  //jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], events->NPU_0, JetResolutionCalculator, random);
	}
	if (printSyncDebug) {
	  cout << "Jet Smearing Factor " << jetEnergySmearFactor << "\n";
	}

	TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
	TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

	//*******************************
	//Add to Type1 Met Correction
	//*******************************
	if (jetPt[i]*JEC > 20 && fabs(jetEta[i]) < 3.0) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - UnCorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - UnCorrJet.Py()  );
	  if (printSyncDebug) cout << "Met Type1 Corr: " << thisJet.Px() - UnCorrJet.Px() << " " << thisJet.Py() - UnCorrJet.Py() << "\n";
	}


	//**********************************************************************************************
	//apply  Pileup Jet ID, don't do this yet for Run2 because we don't know the cut to make yet
	//**********************************************************************************************
	// int level = 2; //loose jet ID
	//if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;


	if (abs(jetPartonFlavor[i]) == 5) {
	  if (!bjet1Found) {
	    bjet1Found = true;
	    events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	    events->bjet1PassLoose = false;
	    events->bjet1PassMedium = false;
	    events->bjet1PassTight = false;
	    if(isCSVL(i)) events->bjet1PassLoose = true;	      
	    if(isCSVM(i)) events->bjet1PassMedium = true;
	    if(isCSVT(i)) events->bjet1PassTight = true;	      
	  } else if ( thisJet.Pt() > events->bjet1.Pt() ) {
	    events->bjet2.SetPtEtaPhiM(events->bjet1.Pt(), events->bjet1.Eta(), events->bjet1.Phi(), events->bjet1.M());
	    events->bjet2PassLoose = events->bjet1PassLoose;
	    events->bjet2PassMedium = events->bjet1PassMedium;
	    events->bjet2PassTight = events->bjet1PassTight;

	    events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	    events->bjet1PassLoose = false;
	    events->bjet1PassMedium = false;
	    events->bjet1PassTight = false;
	    if(isCSVL(i)) events->bjet1PassLoose = true;	      
	    if(isCSVM(i)) events->bjet1PassMedium = true;
	    if(isCSVT(i)) events->bjet1PassTight = true;	      
	  } else {
	    if (!bjet2Found || thisJet.Pt() > events->bjet2.Pt() ) {
	      bjet2Found = true;
	      events->bjet2.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	      events->bjet2PassLoose = false;
	      events->bjet2PassMedium = false;
	      events->bjet2PassTight = false;
	      if(isCSVL(i)) events->bjet2PassLoose = true;	      
	      if(isCSVM(i)) events->bjet2PassMedium = true;
	      if(isCSVT(i)) events->bjet2PassTight = true;	   	    
	    }
	  }
	} //if it's a bjet

	if (jetPt[i]*JEC > 20 && isCSVL(i) ) nBJetsLoose20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVM(i) ) nBJetsMedium20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVT(i) ) nBJetsTight20GeV++;

	if(jetPt[i]*JEC*jetEnergySmearFactor < 40) continue;

	// calculate MHT
	if(jetPt[i]*JEC*jetEnergySmearFactor > 20)
	  {
	    mhx -= jetPt[i]*JEC*jetEnergySmearFactor*cos(jetPhi[i]);
	    mhy -= jetPt[i]*JEC*jetEnergySmearFactor*sin(jetPhi[i]);
	    
	    if (fabs(jetEta[i]) < 3.0)
	      {
		mhx_nohf -= jetPt[i]*JEC*jetEnergySmearFactor*cos(jetPhi[i]);
		mhy_nohf -= jetPt[i]*JEC*jetEnergySmearFactor*sin(jetPhi[i]);
	      }
	  }
	
	if(fabs(jetEta[i]) > 3.0) continue;
	 
	numJetsAbove40GeV++;
	if(jetPt[i]*JEC*jetEnergySmearFactor > 80) numJetsAbove80GeV++;	  
	GoodJets.push_back(thisJet);	  
	  
      } //loop over jets
    
      events->MHT = sqrt(mhx*mhx + mhy*mhy);
      events->MHTnoHF = sqrt(mhx_nohf*mhx_nohf + mhy_nohf*mhy_nohf);

      //sort good jets
      sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());
   
      //Fill two leading jets
      events->jet1.SetPtEtaPhiM(0,0,0,0);
      events->jet2.SetPtEtaPhiM(0,0,0,0);
      events->jet1PassCSVLoose = false;
      events->jet1PassCSVMedium = false;
      events->jet1PassCSVTight = false;
      events->jet2PassCSVLoose = false;
      events->jet2PassCSVMedium = false;
      events->jet2PassCSVTight = false;       
      for (int i=0;i<int(GoodJets.size());i++) {
	if (i==0) {
	  events->jet1.SetPtEtaPhiM(GoodJets[i].Pt(), GoodJets[i].Eta(),GoodJets[i].Phi(),GoodJets[i].M());
	  if(isCSVL(i)) events->jet1PassCSVLoose = true;	      
	  if(isCSVM(i)) events->jet1PassCSVMedium = true;
	  if(isCSVT(i)) events->jet1PassCSVTight = true;	
	}
	if (i==1) {
	  events->jet2.SetPtEtaPhiM(GoodJets[i].Pt(), GoodJets[i].Eta(),GoodJets[i].Phi(),GoodJets[i].M());
	  if(isCSVL(i)) events->jet2PassCSVLoose = true;	      
	  if(isCSVM(i)) events->jet2PassCSVMedium = true;
	  if(isCSVT(i)) events->jet2PassCSVTight = true; 	   	
	}
      }
    
      //***************************************************************
      //compute minDPhi variables
      //***************************************************************
    
      double JER = 0.1; //average jet energy resolution
    
      for (int i=0;i<int(GoodJets.size());i++) {
	if (i>2) continue; //only consider first 3 jets for minDPhi 	  
	double dPhi = fmin( fabs(deltaPhi(metPhi,GoodJets[i].Phi())) , fabs( 3.1415926 - fabs(deltaPhi(metPhi,GoodJets[i].Phi()))));
	if (dPhi < events->minDPhi) events->minDPhi = dPhi;
	double deltaT = 0;
	for(auto& jet2 : GoodJets) {
	  deltaT += pow(JER*jet2.Pt()*sin(fabs(deltaPhi(GoodJets[i].Phi(),jet2.Phi()))),2);
	}      
	double dPhiN = dPhi / atan2( sqrt(deltaT) , metPt);
	if (dPhiN < events->minDPhiN) events->minDPhiN = dPhiN;
      }

    
      //***************************************************************************
      //Compute the razor variables using the selected jets and good leptons
      //***************************************************************************
      vector<TLorentzVector> GoodPFObjects;
      for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
      for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);

      double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
      double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;

      double PFMetnoHFX = metNoHFPt*cos(metNoHFPhi) + MetX_Type1Corr;
      double PFMetnoHFY = metNoHFPt*sin(metNoHFPhi) + MetY_Type1Corr;
    
      TLorentzVector PFMET; PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));
      TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
      TLorentzVector PFMETType0Plus1 = makeTLorentzVectorPtEtaPhiM(metType0Plus1Pt, 0, metType0Plus1Phi, 0);
      // TLorentzVector MyMET = PFMETType1;

      TLorentzVector PFMETnoHFType1;
      PFMETnoHFType1.SetPxPyPzE(PFMetnoHFX, PFMetnoHFY, 0, sqrt(PFMetnoHFX*PFMetnoHFX + PFMetnoHFY*PFMetnoHFY));
      
      TLorentzVector MyMET = PFMETnoHFType1;
	
      if (printSyncDebug) {
	cout << "UnCorrectedMET: " << PFMETUnCorr.Pt() << " " << PFMETUnCorr.Phi() << "\n";
	cout << "Corrected PFMET: " << PFMET.Pt() << " " << PFMET.Phi() << " | X,Y Correction :  " << MetX_Type1Corr << " " << MetY_Type1Corr << "\n";
      }

      events->MR = 0;
      events->Rsq = 0;
	
      //only compute razor variables if we have 2 jets
      //Use Type 1 PFMet for MET
      if (GoodPFObjects.size() >= 2 
	  && GoodJets.size() < 20 //this is a protection again crazy events
	  ) {
	vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	events->MR = computeMR(hemispheres[0], hemispheres[1]); 
	events->Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
	events->RsqnoHF = computeRsq(hemispheres[0], hemispheres[1], PFMETnoHFType1);
      }
	

      if (printSyncDebug)  {
	cout << "MR = " << events->MR << " Rsq = " << events->Rsq << " | "
	     << " Mll = " << (events->lep1 + events->lep2).M() << " | " 
	     << " NJets80 = " << numJetsAbove80GeV << " NJets40 = " << numJetsAbove40GeV << " GoodPFObjects.size() = " << GoodPFObjects.size() << " "
	     << " MET = " << MyMET.Pt() << " MetPhi = " << MyMET.Phi() << " nBTagsMedium = " << nBJetsMedium20GeV << "\n";
      }

      events->MET = MyMET.Pt();
      events->METnoHF = PFMETnoHFType1.Pt();
      events->NJets40 = numJetsAbove40GeV;
      events->NJets80 = numJetsAbove80GeV;
      events->NBJetsLoose = nBJetsLoose20GeV;
      events->NBJetsMedium = nBJetsMedium20GeV;
      events->NBJetsTight = nBJetsTight20GeV;

      events->HT = 0;
      for(auto& pfobj : GoodPFObjects) events->HT += pfobj.Pt();
	
      //compute M_T for lep1 and MET
      events->lep1MT = sqrt(events->lep1.M2() + 2*MyMET.Pt()*events->lep1.Pt()*(1 - cos(deltaPhi(MyMET.Phi(),events->lep1.Phi()))));
      events->lep1MTnoHF = sqrt(events->lep1.M2() + 2*PFMETnoHFType1.Pt()*events->lep1.Pt()*(1 - cos(deltaPhi(PFMETnoHFType1.Phi(),events->lep1.Phi()))));
	
      //save HLT Decisions
      for(int k=0; k<150; ++k) {
	events->HLTDecision[k] = HLTDecision[k];	
      }
    	
      //MET Filter
      events->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
      events->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
      events->Flag_hcalLaserEventFilter = Flag_hcalLaserEventFilter;
      events->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      events->Flag_goodVertices = Flag_goodVertices;
      events->Flag_trackingFailureFilter = Flag_trackingFailureFilter;
      events->Flag_eeBadScFilter = Flag_eeBadScFilter;
      events->Flag_ecalLaserCorrFilter = Flag_ecalLaserCorrFilter;
      events->Flag_trkPOGFilters = Flag_trkPOGFilters;
      events->Flag_trkPOG_manystripclus53X = Flag_trkPOG_manystripclus53X;
      events->Flag_trkPOG_toomanystripclus53X = true;
      events->Flag_trkPOG_logErrorTooManyClusters = Flag_trkPOG_logErrorTooManyClusters;
      events->Flag_METFilters = Flag_METFilters;	


      //***********************************************************************//
      //    Compute razor vars adding photon to the MET
      //***********************************************************************//
	
      if(GoodPhotons.size()>0){

	events->pho1.SetPtEtaPhiM(GoodPhotons[0].Pt(),GoodPhotons[0].Eta(),GoodPhotons[0].Phi(),GoodPhotons[0].M());

	//compute MET with leading photon added
	TLorentzVector m1 = GoodPhotons[0];
	TLorentzVector m2 = MyMET;
	TLorentzVector photonPlusMet_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), 0.0);

	events->MET_NoPho = photonPlusMet_perp.Pt();
	events->METPhi_NoPho = photonPlusMet_perp.Phi();

	//remove leading photon from collection of selected jets
	vector<TLorentzVector> GoodJetsNoLeadPhoton = GoodJets;
	int subtractedIndex = SubtractParticleFromCollection(GoodPhotons[0], GoodJetsNoLeadPhoton);
	if(subtractedIndex >= 0){
	  if(GoodJetsNoLeadPhoton[subtractedIndex].Pt() < 40){ //erase this jet
	    GoodJetsNoLeadPhoton.erase(GoodJetsNoLeadPhoton.begin()+subtractedIndex);
	  }
	}
	    
	//count the number of jets above 80 GeV now
	int numJets80_noPho = 0.;
	for(auto& jet : GoodJetsNoLeadPhoton){
	  if(jet.Pt() > 80) numJets80_noPho++;
	}
	events->NJets80_NoPho = numJets80_noPho;
	    
	//count jets and compute HT
	events->NJets_NoPho = GoodJetsNoLeadPhoton.size();
	float ht_noPho = 0.;
	for(auto& pf : GoodJetsNoLeadPhoton) ht_noPho += pf.Pt();
	events->HT_NoPho = ht_noPho;
	    
	if(GoodJetsNoLeadPhoton.size() >= 2 && GoodJetsNoLeadPhoton.size() <20){
	  //remake the hemispheres using the new jet collection
	  vector<TLorentzVector> hemispheresNoLeadPhoton = getHemispheres(GoodJetsNoLeadPhoton);
	  TLorentzVector PFMET_NOPHO = makeTLorentzVectorPtEtaPhiM(events->MET_NoPho, 0, events->METPhi_NoPho, 0);
	  events->MR_NoPho = computeMR(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1]); 
	  events->Rsq_NoPho = computeRsq(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1], PFMET_NOPHO);
	  events->dPhiRazor_NoPho = fabs(hemispheresNoLeadPhoton[0].DeltaPhi(hemispheresNoLeadPhoton[1]));
	}

      }
    
      //***********************************************************************//
      //    Compute razor vars adding good leptons to the MET
      //***********************************************************************//
      //remove selected muons from collection of selected jets and add them to the MET
      vector<TLorentzVector> GoodJetsNoLeptons = GoodJets;
      TLorentzVector TotalLepVec;
      int NLeptonsAdded = 0;
      for(auto& lep : GoodLeptons){
	NLeptonsAdded++;
	TotalLepVec = TotalLepVec + lep; //add this muon's momentum to the sum
	int subtractedIndex = SubtractParticleFromCollection(lep, GoodJetsNoLeptons);
	if(subtractedIndex >= 0){
	  if(GoodJetsNoLeptons[subtractedIndex].Pt() < 40){ //erase this jet
	    GoodJetsNoLeptons.erase(GoodJetsNoLeptons.begin()+subtractedIndex);
	  }
	}
	
	//only count 1 lepton for 1-lepton option. 
	if (treeTypeOption == 2) {
	  if (NLeptonsAdded >= 1) break;
	}

	//only count 2 lepton for 2-lepton option. 
	if (treeTypeOption == 4) {
	  if (NLeptonsAdded >= 2) break;
	}

     }

      //make the MET vector with the muons (or gen muons) added
      TLorentzVector LepPlusMet_perp = makeTLorentzVector((TotalLepVec + MyMET).Pt(), 0., (TotalLepVec + MyMET).Phi(), 0.);

      //count jets and compute HT
      int njets80NoLep = 0;
      float ht_NoLep   = 0.;
      for(auto& jet : GoodJetsNoLeptons){
	ht_NoLep += jet.Pt();
	if(jet.Pt() > 80) njets80NoLep++;
      }


      //**********************************************************
      //1-lepton Add to MET Variables
      //**********************************************************
      if (treeTypeOption == 2) {
	events->MET_NoW = LepPlusMet_perp.Pt();
	events->METPhi_NoW = LepPlusMet_perp.Phi();
	events->NJets_NoW = GoodJetsNoLeptons.size();
	events->HT_NoW = ht_NoLep;
	events->NJets80_NoW = njets80NoLep;

	//compute razor variables 
	if(events->NJets_NoW > 1 && GoodJets.size()<20){
	  vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodJetsNoLeptons);
	  events->Rsq_NoW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], LepPlusMet_perp);
	  events->MR_NoW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 
	  events->dPhiRazor_NoW = fabs(hemispheresNoW[0].DeltaPhi(hemispheresNoW[1])); 
	}

	//more W kinematics
	if (GoodLeptons.size() > 0) {
	  events->recoWpt = (GoodLeptons[0] + MyMET).Pt();
	  events->recoWphi = (GoodLeptons[0] + MyMET).Phi();
	}
      }


      //**********************************************************
      //2-lepton Add to MET Variables
      //**********************************************************
      if (treeTypeOption == 4) {
	events->MET_NoZ = LepPlusMet_perp.Pt();
	events->METPhi_NoZ = LepPlusMet_perp.Phi();
	events->NJets_NoZ = GoodJetsNoLeptons.size();
	events->HT_NoZ = ht_NoLep;
	events->NJets80_NoZ = njets80NoLep;
    events->METPhi = MyMET.Phi();

    // Compute the recoil for Z
    TVector2 vZPt((TotalLepVec.Pt())*cos(TotalLepVec.Phi()),(TotalLepVec.Pt())*sin(TotalLepVec.Phi())); // Z bosont pT
    TVector2 vMet((MyMET.Pt())*cos(MyMET.Phi()), (MyMET.Pt())*sin(MyMET.Phi())); // MET vector
    TVector2 vU = -1.0*(vZPt+vMet); // recoil vector
    events->u1 = ((TotalLepVec.Px())*(vU.Px()) + (TotalLepVec.Py())*(vU.Py()))/(TotalLepVec.Pt());  // u1 = (pT . u)/|pT|
    events->u2 = ((TotalLepVec.Px())*(vU.Py()) - (TotalLepVec.Py())*(vU.Px()))/(TotalLepVec.Pt());  // u2 = (pT x u)/|pT|

	//Z kinematics
	if(NLeptonsAdded >= 2){
	  events->recoZpt = TotalLepVec.Pt();
	  events->recoZmass = TotalLepVec.M();
	} else {
	  events->recoZpt = -1;
	  events->recoZmass = -1;
	}

	//compute razor variables 
	if(events->NJets_NoZ > 1 && GoodJets.size()<20)
	  {
	    vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodJetsNoLeptons);
	    events->Rsq_NoZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], LepPlusMet_perp);
	    events->MR_NoZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
	    events->dPhiRazor_NoZ = fabs(hemispheresNoZ[0].DeltaPhi(hemispheresNoZ[1])); 
	  }
      }
    


      //****************************************************************************
      //Event Skimming
      //****************************************************************************
      bool passSkim = true;

      // Dilepton skim
      if (leptonSkimOption == 2) {
	if (!(
	      (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	      && (abs(events->lep2Type) == 11  || abs(events->lep2Type) == 13 )
	      && events->lep1PassLoose && events->lep2PassLoose
	      && events->lep1.Pt() > 20 && events->lep2.Pt() > 20
	      )
	    ) passSkim = false;
      }

      //single lepton skim
      if (leptonSkimOption == 1) {
	if (!( 
	      (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	      && events->lep1PassLoose
	      && events->lep1.Pt() > 30
	       )
	    ) passSkim = false;
      }

      //Photon Skim
      if(leptonSkimOption == 5)
	{
	  if (!(
		GoodPhotons.size() > 0
		&& GoodPhotons[0].Pt() > 50	
		)
	      ) passSkim = false;	
	}
      
      //razor skim
      if (razorSkimOption == 1) {
	if ( !( 
	       (events->MR > 300 && events->Rsq > 0.1) || GoodJets.size() >= 20
		)
	     ) passSkim = false;
      }

      //1-lepton add to MET razor skim
      if(razorSkimOption == 2)
	{
	  if (!(
		events->MR_NoW > 300 && events->Rsq_NoW > 0.15
		)
	      ) passSkim = false;		
	}
      

      //Dilepton add to MET razor skim
      if(razorSkimOption == 3) 
	{
	  if (!(
		events->MR_NoZ > 300 && events->Rsq_NoZ > 0.15
		)
	      ) passSkim = false;
	}
      
      //Photon Skim
      if(razorSkimOption == 5)
	{
	  if (!(
		events->MR_NoPho < 300
		&& events->Rsq_NoPho > 0.15
		)
	      ) passSkim = false;	
	}

      //fill event 
      if (passSkim) {
	if (printSyncDebug) cout<<"Filling the tree... " <<endl;
	events->tree_->Fill();
      }

    }//end of event loop


    cout << "Filled Total of " << NEventProcessed << " Events\n";
    cout << "Sum of event weights = " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

