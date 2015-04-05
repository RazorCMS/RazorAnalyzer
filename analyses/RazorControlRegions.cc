
#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
#include "ControlSampleEvents.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

 
void RazorAnalyzer::RazorControlRegions( string outputfilename, int option, bool isRunOne)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    std::vector<JetCorrectorParameters> correctionParameters;

    if (isRunOne) {
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/FT_53_V6_AN1_L1FastJet_AK5PF.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/FT_53_V6_AN1_L2Relative_AK5PF.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/FT_53_V6_AN1_L3Absolute_AK5PF.txt")); 
    } else {
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L2Relative_AK4PFchs.txt"));
      correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt"));    
    }
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorVetoLeptonStudy.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ControlSampleEvents *events = new ControlSampleEvents;
    events->CreateTree();
    events->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
    //begin loop
    if (fChain == 0) return;

    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;

    //for (Long64_t jentry=36834; jentry<fChain->GetEntries();jentry++) {
    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

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

 
        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        
	//******************************************
	//find generated leptons
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
	

	//******************************************
	//sort gen leptons by pt
	//******************************************
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


	// for (int i=0;i<genLeptonIndex.size();i++) cout << "Lepton " << i << " : " << gParticleId[genLeptonIndex[i]] << " | " 
	// 					       << gParticlePt[genLeptonIndex[i]] << " "
	// 					       << gParticleEta[genLeptonIndex[i]] << " "
	// 					       << gParticlePhi[genLeptonIndex[i]] << " " 
	// 					       << " \n";
	// cout << "\n";


	// if (genLeptonIndex.size() != 2) {
	//   cout << "\n";
	//   cout << "\n";
	//   for (int i=0;i<genLeptonIndex.size();i++) cout << "Lepton " << i << " : " << gParticleId[genLeptonIndex[i]] << " | " 
	// 					     << gParticlePt[genLeptonIndex[i]] << " "
	// 					     << gParticleEta[genLeptonIndex[i]] << " "
	// 					     << gParticlePhi[genLeptonIndex[i]] << " " 
	// 					     << " \n";
						  
	//   cout << "\n";

	//   for(int j = 0; j < nGenParticle; j++){
	//     cout << "Particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " | "
	// 	 << gParticlePt[j] << " "
	// 	 << gParticleEta[j] << " "
	// 	 << gParticlePhi[j] << " "
	// 	 << " | " << gParticleMotherId[j] << " , " << gParticleMotherIndex[j] 
	// 	 << "\n";
	//   }
	// } 

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


	//********************************
	//Sort Leptons
	//********************************
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
	
	
	//********************************
	//Fill Leptons
	//********************************
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



	//if ((events->lep1.Pt() > 0 && events->lep1MatchedGenLepIndex < 0) || (events->lep2.Pt() > 0 && events->lep2MatchedGenLepIndex < 0)) {
	// if ((events->lep1.Pt() > 0 && events->lep1MatchedGenLepIndex < 0 && abs(events->lep1Type)==15 && events->lep1PassTight && (events->genlep1Type != 0 || events->genlep2Type != 0))) {
	//   cout << "\n\n";
	//   cout << "lep1: " << events->lep1Type << " | " << events->lep1.Pt() << " " << events->lep1.Eta() << " " << events->lep1.Phi() << " | " << events->lep1PassTight << " " << events->lep1PassLoose << " " << events->lep1PassVeto << "\n";
	//   cout << "lep2: " << events->lep2Type << " | " << events->lep2.Pt() << " " << events->lep2.Eta() << " " << events->lep2.Phi() << " | " << events->lep2PassTight << " " << events->lep2PassLoose << " " << events->lep2PassVeto << "\n";
	//   cout << "genlep1: " << events->genlep1Type << " | " << events->genlep1.Pt() << " " << events->genlep1.Eta() << " " << events->genlep1.Phi() << "\n";
	//   cout << "genlep2: " << events->genlep2Type << " | " << events->genlep2.Pt() << " " << events->genlep2.Eta() << " " << events->genlep2.Phi() << "\n";

	//   for(int j = 0; j < nGenParticle; j++){
	//     cout << "Particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " | "
	// 	 << gParticlePt[j] << " "
	// 	 << gParticleEta[j] << " "
	// 	 << gParticlePhi[j] << " "
	// 	 << " | " << gParticleMotherId[j] << " , " << gParticleMotherIndex[j] 
	// 	 << "\n";	 
	//   }
	//   cout << "\n";
	// }


	bool bjet1Found = false;
	bool bjet2Found = false;
        vector<TLorentzVector> GoodJets;
        int numJetsAbove80GeV = 0;
	int numJetsAbove40GeV = 0;
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
	  
	  // //exclude jet if it doesn't match a selected genjet
	  // bool matchedGenJet = false;
	  // for(uint j = 0; j < nGenJets; j++){
	  //   double thisDR = deltaR(genJetEta[j],genJetPhi[j],jetEta[i],jetPhi[i]);
	  //   if(thisDR < 0.4 && fabs(jetPt[i]-genJetPt[j])/genJetPt[j] < 0.5 ){
	  // 	  matchedGenJet = true;
	  // 	  break;
	  //     }
	  // }
	  // if(!matchedGenJet) continue;
	  
	  //*******************************************************
	  //apply jet iD
	  //*******************************************************
	  int level = 2; //loose jet ID
	  if (!jetPassIDTight[i]) continue;
	  if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;


	  double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						 fixedGridRhoFastjetAll, jetJetArea[i], 
						 JetCorrector);   
	  TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC);

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
		//cout << "jet " << i << " " << jetPartonFlavor[i] << " | " << bjet1Found << " " << bjet2Found << " : " << thisJet.Pt() << "\n";
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


	  if(jetPt[i]*JEC < 40) continue;
	  if(fabs(jetEta[i]) > 3.0) continue;
	  

	  numJetsAbove40GeV++;
	  if(jetPt[i]*JEC > 80) numJetsAbove80GeV++;

	  
	  GoodJets.push_back(thisJet);	  
	  
        } //loop over jets


	//*****************************************************
	//sort good jets
	//*****************************************************
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

	//Make Good Jet Collection, excluding the leading jet
	vector<TLorentzVector> GoodJets_NoLeadJet;
	int leadJetIndex = -1;
	double leadJetPt = 0;
	for (int i=0;i<int(GoodJets.size());i++) {
	  if (GoodJets[i].Pt() > leadJetPt) {
	    leadJetIndex = i;
	    leadJetPt = GoodJets[i].Pt();
	  }
	}
	int numJetsAbove80GeV_NoLeadJet = 0;
	for (int i=0;i<int(GoodJets.size());i++) {
	  if (i != leadJetIndex) {
	    GoodJets_NoLeadJet.push_back(GoodJets[i]);
	    if (GoodJets[i].Pt() > 80) numJetsAbove80GeV_NoLeadJet++;
	  }
	}


        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        vector<TLorentzVector> GoodPFObjects_NoLeadJet;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        for(auto& jet : GoodJets_NoLeadJet) GoodPFObjects_NoLeadJet.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects_NoLeadJet.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
        TLorentzVector PFMET_NoLeadJet = PFMET; if (leadJetIndex >= 0) PFMET_NoLeadJet = PFMET + GoodJets[leadJetIndex];

	events->MR = 0;
	events->Rsq = 0;
	

	//cout << "debug: " << numJetsAbove80GeV << " " << GoodJets.size() << " " << GoodLeptons.size() << " : " << GoodPFObjects.size() << "\n";

	//only compute razor variables if we have 2 jets above 80 GeV
	if (numJetsAbove80GeV >= 2 && GoodJets.size() < 20) {
	  vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	  events->MR = computeMR(hemispheres[0], hemispheres[1]); 
	  events->Rsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	}

	if (numJetsAbove80GeV_NoLeadJet >= 2 && GoodJets.size() < 20) {
	  vector<TLorentzVector> hemispheres_NoLeadJet = getHemispheres(GoodPFObjects_NoLeadJet);
	  events->MR_NoLeadJet = computeMR(hemispheres_NoLeadJet[0], hemispheres_NoLeadJet[1]); 
	  events->Rsq_NoLeadJet = computeRsq(hemispheres_NoLeadJet[0], hemispheres_NoLeadJet[1], PFMET_NoLeadJet);
	}

	events->MET = metPt;
	events->MET_NoLeadJet = PFMET_NoLeadJet.Pt();
	events->NJets40 = numJetsAbove40GeV;
	events->NBJetsLoose = nBJetsLoose20GeV;
	events->NBJetsMedium = nBJetsMedium20GeV;
	events->NBJetsTight = nBJetsTight20GeV;

	events->HT = 0;
	for(auto& pfobj : GoodPFObjects) events->HT += pfobj.Pt();
	
	//compute M_T for lep1 and MET
	events->lep1MT = sqrt(events->lep1.M2() + 2*metPt*events->lep1.Pt()*(1 - cos(deltaPhi(metPhi,events->lep1.Phi()))));
	
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
	

	//skim events
	bool passSkim = false;
	if (option == -1) passSkim = true;
	if (option == 1) {
	  if ( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       && (abs(events->lep2Type) == 11  || abs(events->lep2Type) == 13 )
	      && events->lep1PassLoose && events->lep2PassLoose
	      && events->lep1.Pt() > 20 && events->lep2.Pt() > 20) passSkim = true;
	}
	if (option == 10) {
	  if ((events->MR > 300 && events->Rsq > 0.1) || GoodJets.size() >= 20) passSkim = true;
	}

	//fill event 
	if (passSkim) {
	  events->tree_->Fill();
	}

    }//end of event loop


    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

