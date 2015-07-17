
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
 
void RazorAnalyzer::RazorControlRegions( string outputfilename, int option, bool isData, bool isRunOne)
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
    ControlSampleEvents *events = new ControlSampleEvents;
    
    if (option == 1)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (option == 101)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (option == 201)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (option == 2)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLeptonAdd2MET_Full);
    else if (option == 102)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLeptonAdd2MET_Full);
    else if (option == 3)
      events->CreateTree(ControlSampleEvents::kTreeType_Dilepton_Full);
    else if (option == 103)
      events->CreateTree(ControlSampleEvents::kTreeType_Dilepton_Full);
    else if (option == 4)
      events->CreateTree(ControlSampleEvents::kTreeType_DileptonAdd2MET_Full);
    else if (option == 104)
      events->CreateTree(ControlSampleEvents::kTreeType_DileptonAdd2MET_Full);
    else if (option == 5)
      events->CreateTree(ControlSampleEvents::kTreeType_Photon_Full);
    else if (option == 105)
      events->CreateTree(ControlSampleEvents::kTreeType_Photon_Full);
    else if (option == 6)
      events->CreateTree(ControlSampleEvents::kTreeType_ZeroLepton_Full);
    else {
      events->CreateTree(ControlSampleEvents::kTreeType_Default);
    }
    events->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  

    //*************************************************************************
    //Look over Input File Events
    //*************************************************************************
    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

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
      NEvents->Fill(1.0);
     
	//event info
	events->weight = 1.0;
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
                        	   
	    if(!isVetoMuon(i)) continue;  
	    TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
            GoodLeptons.push_back(thisMuon);
        }


	for(int i = 0; i < nElectrons; i++){

            if(elePt[i] < 5) continue;
            if(fabs(eleEta[i]) > 2.5) continue;

	    //don't count electrons that were already selected as muons
	    bool alreadySelected = false;
	    for (uint j=0; j<GoodLeptons.size(); j++) {
	      if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), eleEta[i],elePhi[i]) < 0.1) alreadySelected = true;
	    }
	    if (alreadySelected) continue;

	    if( ( (!isRunOne && isTightElectron(i)) || (isRunOne && isRunOneTightElectron(i)) )
		&& elePt[i] > 10 ) {
	      TightLeptonType.push_back(11 * -1 * eleCharge[i]);
	      TightLeptonIndex.push_back(i);
	      TightLeptonPt.push_back(elePt[i]);
	    }
	    else if( ( (!isRunOne && isLooseElectron(i)) || (isRunOne && isRunOneLooseElectron(i)) )
		     && elePt[i] > 10 ) {
	      LooseLeptonType.push_back(11 * -1 * eleCharge[i]);
	      LooseLeptonIndex.push_back(i);
	      LooseLeptonPt.push_back(elePt[i]);
	    }
	    else if(isMVANonTrigVetoElectron(i)) {
	      VetoLeptonType.push_back(11 * -1 * eleCharge[i]);
	      VetoLeptonIndex.push_back(i);
	      VetoLeptonPt.push_back(elePt[i]);
	    }
            
	    if (printSyncDebug) cout << "ele " << i << " " << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " : Tight = " << isTightElectron(i) << " Loose = " << isLooseElectron(i) << " Veto = " << isMVANonTrigVetoElectron(i) << " \n";

	    if(!isMVANonTrigVetoElectron(i)) continue; 
	    TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            GoodLeptons.push_back(thisElectron);        
        }

	for(int i = 0; i < nTaus; i++){

	  if (isRunOne) continue;

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
	  
        }


	//Sort Leptons
	tempIndex = -1;
	int tempType = -1;
	int tempPt = -1;
	for(uint i = 0; i < TightLeptonIndex.size() ; i++) {
	  for (uint j=0; j < TightLeptonIndex.size()-1; j++) {
	    if (TightLeptonPt[j] > TightLeptonPt[j]) { 
	      tempIndex = TightLeptonIndex[j];             // swap elements
	      tempType = TightLeptonType[j];
	      tempPt = TightLeptonPt[j];
	      TightLeptonIndex[j] = TightLeptonIndex[j+1];
	      TightLeptonType[j] = TightLeptonType[j+1];
	      TightLeptonPt[j] = TightLeptonPt[j+1]; 
	      TightLeptonIndex[j+1] = tempIndex;
	      TightLeptonType[j+1] = tempType;
	      TightLeptonPt[j+1] = tempPt;
	    }
	  }
	}
	tempIndex = -1;
	tempType = -1;
	tempPt = -1;
	for(uint i = 0; i < LooseLeptonIndex.size() ; i++) {
	  for (uint j=0; j < LooseLeptonIndex.size()-1; j++) {
	    if (LooseLeptonPt[j] > LooseLeptonPt[j]) { 
	      tempIndex = LooseLeptonIndex[j];             // swap elements
	      tempType = LooseLeptonType[j];
	      tempPt = LooseLeptonPt[j];
	      LooseLeptonIndex[j] = LooseLeptonIndex[j+1];
	      LooseLeptonType[j] = LooseLeptonType[j+1];
	      LooseLeptonPt[j] = LooseLeptonPt[j+1]; 
	      LooseLeptonIndex[j+1] = tempIndex;
	      LooseLeptonType[j+1] = tempType;
	      LooseLeptonPt[j+1] = tempPt;
	    }
	  }
	}
	tempIndex = -1;
	tempType = -1;
	tempPt = -1;
	for(uint i = 0; i < VetoLeptonIndex.size() ; i++) {
	  for (uint j=0; j < VetoLeptonIndex.size()-1; j++) {
	    if (VetoLeptonPt[j] > VetoLeptonPt[j]) { 
	      tempIndex = VetoLeptonIndex[j];             // swap elements
	      tempType = VetoLeptonType[j];
	      tempPt = VetoLeptonPt[j];
	      VetoLeptonIndex[j] = VetoLeptonIndex[j+1];
	      VetoLeptonType[j] = VetoLeptonType[j+1];
	      VetoLeptonPt[j] = VetoLeptonPt[j+1]; 
	      VetoLeptonIndex[j+1] = tempIndex;
	      VetoLeptonType[j+1] = tempType;
	      VetoLeptonPt[j+1] = tempPt;
	    }
	  }
	}
	
	
	//************************************************************************
	//Fill Lepton Information, in order of tight, loose, veto
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
	for (uint i=0; i<TightLeptonType.size(); i++) {
	  if (!lep1Found) {
	    events->lep1Type = TightLeptonType[i];
	    double mass = 0.000511;
	    double tmpEta = 0;
	    double tmpPhi = 0;
	    if (abs(TightLeptonType[i]) == 11) {
	      events->lep1.SetPtEtaPhiM(elePt[TightLeptonIndex[i]], eleEta[TightLeptonIndex[i]], elePhi[TightLeptonIndex[i]],mass);	      
	      tmpEta = eleEta[TightLeptonIndex[i]];
	      tmpPhi = elePhi[TightLeptonIndex[i]];
	    } else if (abs(TightLeptonType[i]) == 13) {
	      mass = 0.1057;
	      events->lep1.SetPtEtaPhiM(muonPt[TightLeptonIndex[i]], muonEta[TightLeptonIndex[i]], muonPhi[TightLeptonIndex[i]],mass);
	      tmpEta = muonEta[TightLeptonIndex[i]];
	      tmpPhi = muonPhi[TightLeptonIndex[i]];
	    } else if (abs(TightLeptonType[i]) == 15) {
	      mass = 1.777;
	      events->lep1.SetPtEtaPhiM(tauPt[TightLeptonIndex[i]], tauEta[TightLeptonIndex[i]], tauPhi[TightLeptonIndex[i]],mass);
	      tmpEta = tauEta[TightLeptonIndex[i]];
	      tmpPhi = tauPhi[TightLeptonIndex[i]];
	    }

	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(TightLeptonType[i]) == abs(events->genlep1Type) &&
		 deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(TightLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = true;
	    events->lep1PassLoose = true;
	    events->lep1PassVeto = true;
	    lep1Found = true;
	  } else {
	    if (!lep2Found) {
	      events->lep2Type = TightLeptonType[i];
	      double mass = 0.000511;
	      double tmpEta = 0;
	      double tmpPhi = 0;
	      if (abs(TightLeptonType[i]) == 11) {
		events->lep2.SetPtEtaPhiM(elePt[TightLeptonIndex[i]], eleEta[TightLeptonIndex[i]], elePhi[TightLeptonIndex[i]],mass);	      
		tmpEta = eleEta[TightLeptonIndex[i]];
		tmpPhi = elePhi[TightLeptonIndex[i]];
	      } else if (abs(TightLeptonType[i]) == 13) {
		mass = 0.1057;
		events->lep2.SetPtEtaPhiM(muonPt[TightLeptonIndex[i]], muonEta[TightLeptonIndex[i]], muonPhi[TightLeptonIndex[i]],mass);
		tmpEta = muonEta[TightLeptonIndex[i]];
		tmpPhi = muonPhi[TightLeptonIndex[i]];
	      } else if (abs(TightLeptonType[i]) == 15) {		
		mass = 1.777;
		events->lep2.SetPtEtaPhiM(tauPt[TightLeptonIndex[i]], tauEta[TightLeptonIndex[i]], tauPhi[TightLeptonIndex[i]],mass);
		tmpEta = tauEta[TightLeptonIndex[i]];
		tmpPhi = tauPhi[TightLeptonIndex[i]];
	      }

	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(TightLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(TightLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = true;
	      events->lep2PassLoose = true;
	      events->lep2PassVeto = true;
	      lep2Found = true;
	    }
	  }
	} //loop over Tight leptons	


	for (uint i=0; i<LooseLeptonType.size(); i++) {
	  if (!lep1Found) {
	    events->lep1Type = LooseLeptonType[i];
	    double mass = 0.000511;
	    double tmpEta = 0;
	    double tmpPhi = 0;
	    if (abs(LooseLeptonType[i]) == 11) {
	      events->lep1.SetPtEtaPhiM(elePt[LooseLeptonIndex[i]], eleEta[LooseLeptonIndex[i]], elePhi[LooseLeptonIndex[i]],mass);	      
	      tmpEta = eleEta[LooseLeptonIndex[i]];
	      tmpPhi = elePhi[LooseLeptonIndex[i]];
	    } else if (abs(LooseLeptonType[i]) == 13) {
	      mass = 0.1057;
	      events->lep1.SetPtEtaPhiM(muonPt[LooseLeptonIndex[i]], muonEta[LooseLeptonIndex[i]], muonPhi[LooseLeptonIndex[i]],mass);
	      tmpEta = muonEta[LooseLeptonIndex[i]];
	      tmpPhi = muonPhi[LooseLeptonIndex[i]];
	    } else if (abs(LooseLeptonType[i]) == 15) {
	      mass = 1.777;
	      events->lep1.SetPtEtaPhiM(tauPt[LooseLeptonIndex[i]], tauEta[LooseLeptonIndex[i]], tauPhi[LooseLeptonIndex[i]],mass);
	      tmpEta = tauEta[LooseLeptonIndex[i]];
	      tmpPhi = tauPhi[LooseLeptonIndex[i]];
	    }

	    events->lep1MatchedGenLepIndex = -1;
	    if (abs(LooseLeptonType[i]) == abs(events->genlep1Type)
		&& deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(LooseLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = false;
	    events->lep1PassLoose = true;
	    events->lep1PassVeto = true;
	    lep1Found = true;
	  } else {
	    if (!lep2Found) {
	      events->lep2Type = LooseLeptonType[i];
	      double mass = 0.000511;
	      double tmpEta = 0;
	      double tmpPhi = 0;
	      if (abs(LooseLeptonType[i]) == 11) {
		events->lep2.SetPtEtaPhiM(elePt[LooseLeptonIndex[i]], eleEta[LooseLeptonIndex[i]], elePhi[LooseLeptonIndex[i]],mass);	      
		tmpEta = eleEta[LooseLeptonIndex[i]];
		tmpPhi = elePhi[LooseLeptonIndex[i]];
	      } else if (abs(LooseLeptonType[i]) == 13) {
		mass = 0.1057;
		events->lep2.SetPtEtaPhiM(muonPt[LooseLeptonIndex[i]], muonEta[LooseLeptonIndex[i]], muonPhi[LooseLeptonIndex[i]],mass);
		tmpEta = muonEta[LooseLeptonIndex[i]];
		tmpPhi = muonPhi[LooseLeptonIndex[i]];
	      } else if (abs(LooseLeptonType[i]) == 15) {
		mass = 1.777;
		events->lep2.SetPtEtaPhiM(tauPt[LooseLeptonIndex[i]], tauEta[LooseLeptonIndex[i]], tauPhi[LooseLeptonIndex[i]],mass);
		tmpEta = tauEta[LooseLeptonIndex[i]];
		tmpPhi = tauPhi[LooseLeptonIndex[i]];
	      }

	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(LooseLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(LooseLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = false;
	      events->lep2PassLoose = true;
	      events->lep2PassVeto = true;
	      lep2Found = true;
	    }
	  }
	} //loop over Loose leptons

	for (uint i=0; i<VetoLeptonType.size(); i++) {
	  if (!lep1Found) {
	    events->lep1Type = VetoLeptonType[i];
	    double mass = 0.000511;
	    double tmpEta = 0;
	    double tmpPhi = 0;
	    if (abs(VetoLeptonType[i]) == 11) {
	      events->lep1.SetPtEtaPhiM(elePt[VetoLeptonIndex[i]], eleEta[VetoLeptonIndex[i]], elePhi[VetoLeptonIndex[i]],mass);	      
	      tmpEta = eleEta[VetoLeptonIndex[i]];
	      tmpPhi = elePhi[VetoLeptonIndex[i]];
	    } else if (abs(VetoLeptonType[i]) == 13) {
	      mass = 0.1057;
	      events->lep1.SetPtEtaPhiM(muonPt[VetoLeptonIndex[i]], muonEta[VetoLeptonIndex[i]], muonPhi[VetoLeptonIndex[i]],mass);
	      tmpEta = muonEta[VetoLeptonIndex[i]];
	      tmpPhi = muonPhi[VetoLeptonIndex[i]];
	    } else if (abs(VetoLeptonType[i]) == 15) {
	      mass = 1.777;
	      events->lep1.SetPtEtaPhiM(tauPt[VetoLeptonIndex[i]], tauEta[VetoLeptonIndex[i]], tauPhi[VetoLeptonIndex[i]],mass);
	      tmpEta = tauEta[VetoLeptonIndex[i]];
	      tmpPhi = tauPhi[VetoLeptonIndex[i]];
	    }

	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(VetoLeptonType[i]) == abs(events->genlep1Type)
		 && deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(VetoLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = false;
	    events->lep1PassLoose = false;
	    events->lep1PassVeto = true;
	    lep1Found = true;
	  } else {
	    if (!lep2Found) {
	      events->lep2Type = VetoLeptonType[i];
	      double mass = 0.000511;
	      double tmpEta = 0;
	      double tmpPhi = 0;
	      if (abs(VetoLeptonType[i]) == 11) {
		events->lep2.SetPtEtaPhiM(elePt[VetoLeptonIndex[i]], eleEta[VetoLeptonIndex[i]], elePhi[VetoLeptonIndex[i]],mass);	      
		tmpEta = eleEta[VetoLeptonIndex[i]];
		tmpPhi = elePhi[VetoLeptonIndex[i]];
	      } else if (abs(VetoLeptonType[i]) == 13) {
		mass = 0.1057;
		events->lep2.SetPtEtaPhiM(muonPt[VetoLeptonIndex[i]], muonEta[VetoLeptonIndex[i]], muonPhi[VetoLeptonIndex[i]],mass);
		tmpEta = muonEta[VetoLeptonIndex[i]];
		tmpPhi = muonPhi[VetoLeptonIndex[i]];
	      } else if (abs(VetoLeptonType[i]) == 15) {
		mass = 1.777;
		events->lep2.SetPtEtaPhiM(tauPt[VetoLeptonIndex[i]], tauEta[VetoLeptonIndex[i]], tauPhi[VetoLeptonIndex[i]],mass);
		tmpEta = tauEta[VetoLeptonIndex[i]];
		tmpPhi = tauPhi[VetoLeptonIndex[i]];
	      }

	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(VetoLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(VetoLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(tmpEta,tmpPhi,events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = false;
	      events->lep2PassLoose = false;
	      events->lep2PassVeto = true;
	      lep2Found = true;
	    }
	  }
	} //loop over Veto leptons

	//save this for the 1-lepton mini ntuples
	if (lep1Found) {
	  events->lep1Pt = events->lep1.Pt();
	  events->lep1Eta = events->lep1.Eta();
	} else {
	  events->lep1Pt = -99;
	  events->lep1Eta = -99;
	}

	if (printSyncDebug) {
	  cout << "\n\n";
	  cout << "lep1: " << events->lep1Type << " | " << events->lep1.Pt() << " " << events->lep1.Eta() << " " << events->lep1.Phi() << " | Tight = " << events->lep1PassTight << " Loose = " << events->lep1PassLoose << " Veto = " << events->lep1PassVeto << "\n";
	  cout << "lep2: " << events->lep2Type << " | " << events->lep2.Pt() << " " << events->lep2.Eta() << " " << events->lep2.Phi() << " | Tight = " << events->lep2PassTight << " Loose = " << events->lep2PassLoose << " Veto = " << events->lep2PassVeto << "\n";	 
	}



	//************************************************************************
	//Find all Relevent Jets	
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

       for(int i = 0; i < nJets; i++){

	  //exclude selected muons and electrons from the jet collection
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
		 << " | csv = " << jetCSV[i] << " passCSVL = " << isOldCSVL(i) << " passCSVM = " << isOldCSVM(i) << " " << "\n";
	  }



	  //*******************************************************
	  //apply jet iD
	  //*******************************************************
	  if (!jetPassIDTight[i]) continue;

	  //*******************************************************
	  //Correct Jet Energy Scale and Resolution
	  //*******************************************************
	  double tmpRho = fixedGridRhoFastjetAll;
	  if (isRunOne) tmpRho = fixedGridRhoAll;
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
	    jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], events->NPU_0, JetResolutionCalculator, random);
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
	  //apply  Pileup Jet ID
	  //*******************************************************
	  int level = 2; //loose jet ID
	  //if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;


	  if (abs(jetPartonFlavor[i]) == 5) {
	    if (!bjet1Found) {
	      bjet1Found = true;
	      events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	      events->bjet1PassLoose = false;
	      events->bjet1PassMedium = false;
	      events->bjet1PassTight = false;
	      if((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) events->bjet1PassLoose = true;	      
	      if((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) events->bjet1PassMedium = true;
	      if((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) events->bjet1PassTight = true;	      
	    } else if ( thisJet.Pt() > events->bjet1.Pt() ) {
	      events->bjet2.SetPtEtaPhiM(events->bjet1.Pt(), events->bjet1.Eta(), events->bjet1.Phi(), events->bjet1.M());
	      events->bjet2PassLoose = events->bjet1PassLoose;
	      events->bjet2PassMedium = events->bjet1PassMedium;
	      events->bjet2PassTight = events->bjet1PassTight;

	      events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	      events->bjet1PassLoose = false;
	      events->bjet1PassMedium = false;
	      events->bjet1PassTight = false;
	      if((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) events->bjet1PassLoose = true;	      
	      if((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) events->bjet1PassMedium = true;
	      if((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) events->bjet1PassTight = true;	      
	    } else {
	      if (!bjet2Found || thisJet.Pt() > events->bjet2.Pt() ) {
		bjet2Found = true;
		events->bjet2.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
		events->bjet2PassLoose = false;
		events->bjet2PassMedium = false;
		events->bjet2PassTight = false;
		if((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) events->bjet2PassLoose = true;	      
		if((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) events->bjet2PassMedium = true;
		if((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) events->bjet2PassTight = true;	   	    
	      }
	    }
	  } //if it's a bjet

	  if (jetPt[i]*JEC > 20 && ((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) ) nBJetsLoose20GeV++;
	  if (jetPt[i]*JEC > 20 && ((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) ) nBJetsMedium20GeV++;
	  if (jetPt[i]*JEC > 20 && ((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) ) nBJetsTight20GeV++;

	  if(jetPt[i]*JEC*jetEnergySmearFactor < 40) continue;
	  if(fabs(jetEta[i]) > 3.0) continue;
	 

	  numJetsAbove40GeV++;
	  if(jetPt[i]*JEC*jetEnergySmearFactor > 80) numJetsAbove80GeV++;	  
	  GoodJets.push_back(thisJet);	  
	  
        } //loop over jets


	//sort good jets
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
	    if(((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) ) events->jet1PassCSVLoose = true;	      
	    if(((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) ) events->jet1PassCSVMedium = true;
	    if(((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) ) events->jet1PassCSVTight = true;	
	  }
	  if (i==1) {
	    events->jet2.SetPtEtaPhiM(GoodJets[i].Pt(), GoodJets[i].Eta(),GoodJets[i].Phi(),GoodJets[i].M());
	    if(((!isRunOne && isCSVL(i)) || (isRunOne && isOldCSVL(i))) ) events->jet2PassCSVLoose = true;	      
	    if(((!isRunOne && isCSVM(i)) || (isRunOne && isOldCSVM(i))) ) events->jet2PassCSVMedium = true;
	    if(((!isRunOne && isCSVT(i)) || (isRunOne && isOldCSVT(i))) ) events->jet2PassCSVTight = true; 	   	
	  }
	}

	//compute minDPhi variables
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


        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);

	double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
	double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;

	TLorentzVector PFMET; PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));
        TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

	if (printSyncDebug) {
	  cout << "UnCorrectedMET: " << PFMETUnCorr.Pt() << " " << PFMETUnCorr.Phi() << "\n";
	  cout << "Corrected PFMET: " << PFMET.Pt() << " " << PFMET.Phi() << " | X,Y Correction :  " << MetX_Type1Corr << " " << MetY_Type1Corr << "\n";
	}

	events->MR = 0;
	events->Rsq = 0;
	

	//only compute razor variables if we have 2 jets above 80 GeV
	if (GoodPFObjects.size() >= 2 && GoodJets.size() < 20
	    //&& numJetsAbove80GeV >= 2 //Si: I think we don't need this requirement at this point
	    ) {
	  vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	  events->MR = computeMR(hemispheres[0], hemispheres[1]); 
	  events->Rsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	}
	

	if (printSyncDebug)  {
	  cout << "MR = " << events->MR << " Rsq = " << events->Rsq << " | "
	       << " Mll = " << (events->lep1 + events->lep2).M() << " | " 
	       << " NJets80 = " << numJetsAbove80GeV << " NJets40 = " << numJetsAbove40GeV << " GoodPFObjects.size() = " << GoodPFObjects.size() << " "
	       << " MET = " << metPt << " MetPhi = " << metPhi << " nBTagsMedium = " << nBJetsMedium20GeV << "\n";
	}

	events->MET = PFMET.Pt();
	events->NJets40 = numJetsAbove40GeV;
	events->NJets80 = numJetsAbove80GeV;
	events->NBJetsLoose = nBJetsLoose20GeV;
	events->NBJetsMedium = nBJetsMedium20GeV;
	events->NBJetsTight = nBJetsTight20GeV;

	events->HT = 0;
	for(auto& pfobj : GoodPFObjects) events->HT += pfobj.Pt();
	
	//compute M_T for lep1 and MET
	events->lep1MT = sqrt(events->lep1.M2() + 2*PFMET.Pt()*events->lep1.Pt()*(1 - cos(deltaPhi(PFMET.Phi(),events->lep1.Phi()))));
	
	//save HLT Decisions
	for(int k=0; k<100; ++k) {
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


	///////////////////////////////
	////// Photon Ntuple Part /////
	///////////////////////////////
        for(int i = 29; i <= 34; i++){
	  if(HLTDecision[i] == 1) events->HLT_Photon = true;
        }
	
        if(isData && events->HLT_Photon){
            //save the trigger bits
            if(HLTDecision[34] == 1){
                events->HLT_Photon160 = true;
            }
            if(HLTDecision[33] == 1){
                events->HLT_Photon150 = true;
            }
            if(HLTDecision[32] == 1){
                events->HLT_Photon135 = true;
            }
            if(HLTDecision[31] == 1){
                events->HLT_Photon90 = true;
            }
            if(HLTDecision[30] == 1){
                events->HLT_Photon75 = true;
            }
            if(HLTDecision[29] == 1){
                events->HLT_Photon50 = true;
            }
        }

        //****************************************************//
        //             Select photons                         //
        //****************************************************//
        vector<TLorentzVector> GoodPhotons;
        int nPhotonsAbove40GeV = 0;

        for(int i = 0; i < nPhotons; i++){


	  if(phoPt[i] < 10) continue;
            if(fabs(phoEta[i]) > 2.5) continue;

            if(isRunOne){
                // if(!isTightRunOnePhoton(i)) continue;
	      if(!isMediumRunOnePhoton(i)) continue;
            }
            else{
                if(!isTightPhoton(i)) continue;
            }

            if(phoPt[i] > 40) nPhotonsAbove40GeV++;
            TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], pho_RegressionE[i]);
            GoodPhotons.push_back(thisPhoton);
        }

	events->nSelectedPhotons = nPhotonsAbove40GeV;

	//****************************************************//
        //    Compute razor vars for DY, W, Gamma samples     //
        //****************************************************//
        //photons
	events->pho1.SetPtEtaPhiM(0,0,0,0);
	
        if(GoodPhotons.size()>0){
            sort(GoodPhotons.begin(), GoodPhotons.end(), greater_than_pt());

            //compute MET with leading photon added
            TLorentzVector m1 = GoodPhotons[0];
            TLorentzVector m2 = PFMET;
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

	    events->pho1 = GoodPhotons[0];
        }
	
	//****************************************************//
        //               Select muons                         //
        //****************************************************//
        vector<TLorentzVector> GoodMuons; 
        vector<TLorentzVector> GoodMuonsTight;
	int nvetomuons = 0;
	int nloosemuons = 0;
	int ntightmuons = 0;

        for(int i = 0; i < nMuons; i++){
	  
	  if(!isLooseMuon(i)) continue;
	  if(muonPt[i] < 10) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
	  TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	  
	  if(isVetoMuon(i)) nvetomuons++;
	  if(isTightMuon(i)){
	    ntightmuons++;
	    GoodMuonsTight.push_back(thisMuon);
	    
	    nloosemuons++;
	    
	    GoodMuons.push_back(thisMuon);
	  }
	}
	events->nVetoMuons  = nvetomuons;
	events->nTightMuons = ntightmuons;
	events->nLooseMuons = nloosemuons;
	
	//remove selected muons from collection of selected jets and add them to the MET
	vector<TLorentzVector> GoodJetsNoMuons = GoodJets;
	TLorentzVector TotalMuonVec;
	for(auto& mu : GoodMuons){
	  TotalMuonVec = TotalMuonVec + mu; //add this muon's momentum to the sum
	  int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoMuons);
	  if(subtractedIndex >= 0){
	    if(GoodJetsNoMuons[subtractedIndex].Pt() < 40){ //erase this jet
	      GoodJetsNoMuons.erase(GoodJetsNoMuons.begin()+subtractedIndex);
	    }
	  }
	}
	
        //remove selected TIGHT muons from collection of selected jets and add them to the MET
        vector<TLorentzVector> GoodJetsNoTightMuons = GoodJets;
        TLorentzVector TotalTightMuonVec;
        for(auto& mu : GoodMuonsTight){
            TotalTightMuonVec = TotalTightMuonVec + mu; //add this muon's momentum to the sum
            int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoTightMuons);
            if(subtractedIndex >= 0){
                if(GoodJetsNoTightMuons[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoTightMuons.erase(GoodJetsNoTightMuons.begin()+subtractedIndex);
                }
            }
        }

        //make the MET vector with the muons (or gen muons) added
        TLorentzVector ZPlusMet_perp = makeTLorentzVector((TotalMuonVec + PFMET).Pt(), 0., (TotalMuonVec + PFMET).Phi(), 0.);
        events->MET_NoZ = ZPlusMet_perp.Pt();
        events->METPhi_NoZ = ZPlusMet_perp.Phi();

        TLorentzVector WPlusMet_perp = makeTLorentzVector((TotalTightMuonVec + PFMET).Pt(), 0., (TotalTightMuonVec + PFMET).Phi(), 0.);
        events->MET_NoW = WPlusMet_perp.Pt();
        events->METPhi_NoW = WPlusMet_perp.Phi(); 

        //count jets and compute HT
        //Z
	int njets80noZ = 0;
	int njets80noW = 0;
	float ht_noZ   = 0.;
	float ht_noW   = 0.;
        events->NJets_NoZ = GoodJetsNoMuons.size();
        for(auto& jet : GoodJetsNoMuons){
	  ht_noZ += jet.Pt();
	  if(jet.Pt() > 80) njets80noZ++;
        }
	events->HT_NoZ = ht_noZ;
	events->NJets80_NoZ = njets80noZ;
	 
        //W
        events->NJets_NoW = GoodJetsNoTightMuons.size();
        for(auto& jet : GoodJetsNoTightMuons){
	  ht_noW += jet.Pt();
	  if(jet.Pt() > 80) njets80noW++; 
        }
	events->HT_NoW = ht_noW;
	events->NJets80_NoW = njets80noW;
 
        //get reco Z information
        if(GoodMuons.size() >= 1){
	  events->recoZpt = TotalMuonVec.Pt();
	  events->recoZmass = TotalMuonVec.M();
        }

        //compute reco Z information and razor variables for DY
        if(events->NJets_NoZ > 1 && GoodJets.size()<20)
        {
            vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodJetsNoMuons);
            events->Rsq_NoZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], ZPlusMet_perp);
            events->MR_NoZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
            events->dPhiRazor_NoZ = fabs(hemispheresNoZ[0].DeltaPhi(hemispheresNoZ[1])); 
        }
        //razor variables using tight muons (for W)
        if(events->NJets_NoW > 1 && GoodJets.size()<20){
            vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodJetsNoTightMuons);
            events->Rsq_NoW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], WPlusMet_perp);
            events->MR_NoW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 
            events->dPhiRazor_NoW = fabs(hemispheresNoW[0].DeltaPhi(hemispheresNoW[1])); 
        }

        //for W, also get the transverse mass of the first tight muon and the MET
        if(GoodMuonsTight.size() > 0) 
        {
            TLorentzVector m1 = GoodMuonsTight[0];
            TLorentzVector m2 = PFMET;
            double deltaPhiLepMet = m1.DeltaPhi(m2);

	    //Only use this lep1MT if we are using the option with leptons added to MET
	    if  (option == 2 || option == 102) {
	      events->lep1MT = sqrt(2*m2.Pt()*m1.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) ); //transverse mass calculation
	    }

            //store reco W information
            events->recoWpt = (m1+m2).Pt();
            events->recoWphi = (m1+m2).Phi();
        }
	
	////////////////////////////////////
	////// End Photon Ntuple Part //////
	////////////////////////////////////
	//skim events
	bool passSkim = false;
	if (option == -1) passSkim = true;


	// Dilepton skim
	if (option == 103) { 
	  if ( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       && (abs(events->lep2Type) == 11  || abs(events->lep2Type) == 13 )
	       && events->lep1PassLoose && events->lep2PassLoose
	       && events->lep1.Pt() > 20 && events->lep2.Pt() > 20) passSkim = true;
	}

	//single tight lepton skim
	if (option == 101) { 
	  if ( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       && events->lep1PassTight
	       && events->lep1.Pt() > 30) passSkim = true;
	}

	//single tight lepton plus razor skim
	if (option == 102) { 
	  if ( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       && events->lep1PassTight
	       && events->lep1.Pt() > 30
	       && ( (events->MR > 300 && events->Rsq > 0.1) || GoodJets.size() >= 20)
	       ) {
	    passSkim = true;
	  }
	}

	//Razor skim
	if (option == 100) { 
	  if ((events->MR > 300 && events->Rsq > 0.1) || GoodJets.size() >= 20) passSkim = true;
	}

	//Single lepton skim for 1L CR with adding lepton to MET
	if(option == 102) // kTreeType_OneLeptonAdd2MET_Reduced
	  {
	    passSkim = true;
            if(events->NJets80 < 2) passSkim = false; //event fails to have two 80 GeV jets
            if(events->MR < 300 && events->MR_NoW < 300) passSkim = false;
            if(events->Rsq < 0.15 && events->Rsq_NoW < 0.15) passSkim = false;
            if(GoodMuons.size() == 0) passSkim = false; //don't save event if no muons or photons
	  }


	//Dilepton skim with adding leptons to MET
	if(option == 104) 
	  {
	    passSkim = true;
            if(events->NJets80 < 2) passSkim = false; //event fails to have two 80 GeV jets
            if(events->MR < 300 && events->MR_NoZ < 300) passSkim = false;
            if(events->Rsq < 0.15 && events->Rsq_NoZ < 0.15) passSkim = false;
            if(GoodMuons.size() < 2 ) passSkim = false; //don't save event if no muons or photons
	  }

	//Photon+Jet Skim
	if(option == 105) // kTreeType_Photon_Reduced
	  {
	    passSkim = true;	    
            if(events->NJets80 < 2) passSkim = false; //event fails to have two 80 GeV jets
            if(events->MR < 300 && events->MR_NoPho < 300) passSkim = false;
            if(events->Rsq < 0.15 && events->Rsq_NoPho < 0.15) passSkim = false;
	  }

	//fill event 
	if (passSkim) {
	  cout<<"Filling the tree... " <<endl;
	  events->tree_->Fill();
	}

    }//end of event loop


    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

