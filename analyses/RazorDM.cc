#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorAnalyzer::RazorDM(string outFileName, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    if (outFileName.empty()){
        cout << "RazorDM: Output filename not specified!" << endl << "Using default output name RazorDM.root" << endl;
        outFileName = "RazorDM.root";
    }
    TFile outFile(outFileName.c_str(), "RECREATE");


    //Correct the _uncorr variables' Pt and E                                                                                                                                  
    std::vector<JetCorrectorParameters> correctionParameters;
    //Run1 Corrections
    /*
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L1FastJet_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L2Relative_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L3Absolute_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L2L3Residual_AK5PF.txt"));
    */
    //PHYS14 corrections
    correctionParameters.push_back( JetCorrectorParameters("data/PHYS14_V2_MC_L1FastJet_AK4PF.txt") );
    correctionParameters.push_back( JetCorrectorParameters("data/PHYS14_V2_MC_L2Relative_AK4PF.txt") );
    correctionParameters.push_back( JetCorrectorParameters("data/PHYS14_V2_MC_L3Absolute_AK4PF.txt") );
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector( correctionParameters );

    //one tree to hold all events
    TTree *razorTree = new TTree("RazorDM", "Info on selected razor DM events");
    
    //separate trees for individual boxes
    map<string, TTree*> razorBoxes;
    vector<string> boxNames;
    boxNames.push_back("MuMu");
    boxNames.push_back("EleEle");
    boxNames.push_back("Mu");
    boxNames.push_back("Ele");
    boxNames.push_back("MultiJet");
    boxNames.push_back("TwoBJet");
    boxNames.push_back("OneBJet");
    for(size_t i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nBTaggedJetsL, nBTaggedJetsM, nBTaggedJetsT;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
    float MuonPt[5], MuonEta[5], MuonPhi[5], MuonE[5];
    float ElePt[5], EleEta[5], ElePhi[5], EleE[5];
    float theMR;
    float theRsq;
    float pfMET, pfMETcorr;
    RazorBox box;
    float JetPt_uncorr[30], JetEta_uncorr[30], JetPhi_uncorr[30], JetE_uncorr[30];
    float JetPt[30], JetEta[30], JetPhi[30], JetE[30];
    bool hasMatchingGenJet[30];
    int matchingGenJetIndex[30];
    

    //set branches on big tree
    if(combineTrees){
        razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nBTaggedJetsL", &nBTaggedJetsL, "nBTaggedJetsL/I");
	razorTree->Branch("nBTaggedJetsM", &nBTaggedJetsM, "nBTaggedJetsM/I");
	razorTree->Branch("nBTaggedJetsT", &nBTaggedJetsT, "nBTaggedJetsT/I");
        razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
	razorTree->Branch("pfMET", &pfMET, "pfMET/F");
	razorTree->Branch("pfMETcorr", &pfMETcorr, "pfMETcorr/F");
        razorTree->Branch("box", &box, "box/I");
    }
    //set branches on all trees

    else{ 
        for(auto& box : razorBoxes){
            box.second->Branch("nBTaggedJetsL", &nBTaggedJetsL, "nBTaggedJetsL/I");
	    box.second->Branch("nBTaggedJetsM", &nBTaggedJetsM, "nBTaggedJetsM/I");
	    box.second->Branch("nBTaggedJetsT", &nBTaggedJetsT, "nBTaggedJetsT/I");
            box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
	    box.second->Branch("MuonPt", MuonPt,"MuonPt[nLooseMuons]/F");
	    box.second->Branch("MuonEta", MuonEta,"MuonEta[nLooseMuons]/F");
	    box.second->Branch("MuonPhi", MuonPhi,"MuonPhi[nLooseMuons]/F");
	    box.second->Branch("MuonE", MuonE,"MuonE[nLooseMuons]/F");
            box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
            box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
	    box.second->Branch("ElePt", ElePt,"ElePt[nLooseElectrons]/F");
            box.second->Branch("EleEta", EleEta,"EleEta[nLooseElectrons]/F");
            box.second->Branch("ElePhi", ElePhi,"ElePhi[nLooseElectrons]/F");
            box.second->Branch("EleE", EleE,"EleE[nLooseElectrons]/F");
            box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
            box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
            box.second->Branch("MR", &theMR, "MR/F");
            box.second->Branch("Rsq", &theRsq, "Rsq/F");
	    box.second->Branch("pfMET", &pfMET, "pfMET/F");
	    box.second->Branch("pfMETcorr", &pfMETcorr, "pfMETcorr/F");

	    //ADDED LINES

	    box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
	    box.second->Branch("JetE_uncorr", JetE_uncorr, "JetE_uncorr[nSelectedJets]/F");
	    box.second->Branch("JetPt_uncorr", JetPt_uncorr, "JetPt_uncorr[nSelectedJets]/F");
	    box.second->Branch("JetPhi_uncorr", JetPhi_uncorr, "JetPhi_uncorr[nSelectedJets]/F");
	    box.second->Branch("JetEta_uncorr", JetEta_uncorr, "JetEta_uncorr[nSelectedJets]/F");

	    box.second->Branch("genMetPt", &genMetPt, "genMetPt/F");
            box.second->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
	    box.second->Branch("metPt", &metPt, "metPt/F");
	    box.second->Branch("metPhi", &metPhi, "metPhi/F");
	    box.second->Branch("metType0Pt", &metType0Pt, "metType0Pt/F");
	    box.second->Branch("metType0Phi", &metType0Phi, "metType0Phi/F");
	    box.second->Branch("metType1Pt", &metType1Pt, "metType1Pt/F");
	    box.second->Branch("metType1Phi", &metType1Phi, "metType1Phi/F");
	    box.second->Branch("metType0Plus1Pt", &metType0Plus1Pt, "metType0Plus1Pt/F");
	    box.second->Branch("metType0Plus1Phi", &metType0Plus1Phi, "metType0Plus1Phi/F");

	    box.second->Branch("nGenJets",&nGenJets, "nGenJets/I");
            box.second->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
            box.second->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
            box.second->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
            box.second->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
	    

	    box.second->Branch("JetE", JetE, "JetE[nSelectedJets]/F");
	    box.second->Branch("JetPt", JetPt, "JetPt[nSelectedJets]/F");
	    box.second->Branch("JetEta", JetEta, "JetEta[nSelectedJets]/F");
	    box.second->Branch("JetPhi", JetPhi, "JetPhi[nSelectedJets]/F");
	    
	    box.second->Branch("hasMatchingGenJet", hasMatchingGenJet, "hasMatchingGenJet[nSelectedJets]/O");
	    box.second->Branch("matchingGenJetIndex", matchingGenJetIndex, "matchingGenJetIndex[nSelectedJets]/I");

	    box.second->Branch("HLTDecision", HLTDecision, "HLTDecision[100]/O");
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

        //fill normalization histogram
        NEvents->Fill(1.0);

        //reset tree variables
        nSelectedJets = 0;
        nBTaggedJetsL = 0;
	nBTaggedJetsM = 0;
	nBTaggedJetsT =0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nTightTaus = 0;
        theMR = -1.0;
        theRsq = -1.0;
	for ( int j = 0; j < 30; j++ )
          {
	    JetE_uncorr[j] = 0.0;
	    JetPt_uncorr[j] = 0.0;
	    JetPhi_uncorr[j] = 0.0;
	    JetEta_uncorr[j] = 0.0;
            JetE[j] = 0.0;
	    JetPt[j] = 0.0;
	    JetPhi[j] = 0.0;
	    JetEta[j] = 0.0;
          }
        if(combineTrees) box = NONE;

	//DEBUGGING GENJETS OUTPUT
	/*	if (jentry < 5){
	  cout << "nGenJets: " << nGenJets << endl;
	  for (int j = 0; j < nGenJets; j++){
	    cout << "genJetE: " << genJetE[j] << endl;
	    cout << "genJetPt: " << genJetPt[j] << endl;
	    cout << "genJetEta: " << genJetEta[j] << endl;
	    cout << "genJetPhi: " << genJetPhi[j] << endl;
	  }
	}
	*/



        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        
        vector<TLorentzVector> LooseMu; //Muons use to compute MET
        for(int i = 0; i < nMuons; i++){
	  if(muonPt[i] < 10.0 || fabs(muonEta[i]) > 2.4) continue;  
	  if(isLooseMuon(i)){
	    nLooseMuons++;
	    TLorentzVector thisMuon = 
	      makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	    LooseMu.push_back(thisMuon);
	  }
	  if(isTightMuon(i))nTightMuons++;
        }
	//if(nLooseMuons >= 2)std::cout << "Two Muons" << std::endl;
	vector<TLorentzVector> LooseEle; //Electrons use to compute MET
        for(int i = 0; i < nElectrons; i++){
	  if(elePt[i] < 10.0 || fabs(eleEta[i]) > 2.4) continue;
	  if(isLooseElectron(i)){
	    nLooseElectrons++;
	    TLorentzVector thisElectron = 
	      makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
	    LooseEle.push_back(thisElectron);
	  }
	  if(isTightElectron(i))nTightElectrons++;
        }
	
	//Tau Veto
	//if(nTaus > 0)continue;
	//if(nLooseMuons >= 2)std::cout << "Tau Veto Requirement" << std::endl;
	
	for(int i = 0; i < nTaus; i++){
	  if(!isTightTau(i)) continue; 
	  nTightTaus++;
	}
	if(nTightTaus > 0)continue;
	//if(nLooseMuons >= 2)std::cout << "Tau Veto Requirement" << std::endl;
	
        vector<TLorentzVector> GoodJets;
	vector<TLorentzVector> GoodJets_uncorr;
        int numJetsAbove80GeV = 0;

	
        for(int i = 0; i < nJets; i++){
	  if(jetPt[i] < 30.0 || fabs(jetEta[i]) > 3.0) continue;

	  //ADDED LINES
	  //int level = 2; //3rd bit of jetPileupIdFlag
	  
	  /*
	    In the case of v1p5 ntuples the jetPassIDLoose flag is not implemented,
	    please comment out the continue statement below
	  */
	  if ( !jetPassIDLoose[i] ) continue;
	  // if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;	  
	  
	  
	  //Removing Muons from the jet collection
	  double deltaR = -1.0;
	  TLorentzVector thisJet_uncorr = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);


	  //doing the correction
	  double JEC = JetEnergyCorrectionFactor( jetPt[i], JetEta[i], jetPhi[i], jetE[i],
						  fixedGridRhoAll, jetJetArea[i],
						  JetCorrector );
	  if (jentry == 1) cout << "FixedGridRhoAll: " << fixedGridRhoAll << endl << "jetJetArea[" << i << "]: " << jetJetArea[i] << endl; 

	  TLorentzVector thisJet = makeTLorentzVector(JEC*jetPt[i], jetEta[i], jetPhi[i], JEC*jetE[i]);


	  for(auto& Mu : LooseMu){
	    double thisDR = thisJet_uncorr.DeltaR(Mu);
	    if(deltaR < 0.0 || thisDR < deltaR) deltaR = thisDR;
	  }
	  if(deltaR > 0.0 && deltaR < 0.3) continue; //jet matches a Loose Muon
	  //Removing Electrons form the jet collection
	  deltaR = -1.0;
	  for(auto& Ele : LooseEle){
            double thisDR = thisJet_uncorr.DeltaR(Ele);
            if(deltaR < 0.0 || thisDR < deltaR) deltaR = thisDR;
          }
          if(deltaR > 0.0 && deltaR < 0.3) continue; //jet matches a Loose Electron
          
	  if(thisJet.Pt() > 80.0) numJetsAbove80GeV++;
	  GoodJets_uncorr.push_back(thisJet_uncorr);
	  
	  GoodJets.push_back(thisJet);
	  nSelectedJets++;
	  
	 
	  
	  if(isCSVL(i)) nBTaggedJetsL++;
	  if(isCSVM(i)) nBTaggedJetsM++;
	  if(isCSVT(i)) nBTaggedJetsT++;
        }
	
        if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets
	//if(nLooseMuons >= 2)std::cout << "Two Jet Requirement" << std::endl;
	
	int jIndex = 0;
	for (auto& tmpJet: GoodJets_uncorr)
	  {
	    JetPt_uncorr[jIndex] = tmpJet.Pt();
	    JetE_uncorr[jIndex] = tmpJet.E();
	    JetPhi_uncorr[jIndex] = tmpJet.Phi();
	    JetEta_uncorr[jIndex] = tmpJet.Eta();
	    jIndex=jIndex+1;
	  }  

	//Find whether there is a matching genJet for the reconstructed jet

	
	jIndex = 0;
	for (auto& tmpJet: GoodJets)
	  { 
	      matchingGenJetIndex[jIndex] = findClosestGenJet(JetEta_uncorr[jIndex], JetPhi_uncorr[jIndex]);
	      if (matchingGenJetIndex[jIndex] == -1)
		hasMatchingGenJet[jIndex] = false;
	      else
		hasMatchingGenJet[jIndex] = true;

	      
	      jIndex++;
	  }


	
	jIndex = 0;
	for (auto& tmpJet: GoodJets)
	  {
	    

	    JetPt[jIndex] = tmpJet.Pt();
	    JetE[jIndex] = tmpJet.E();
	    JetEta[jIndex] = tmpJet.Eta();
	    JetPhi[jIndex] = tmpJet.Phi();
	  
	    jIndex++;

	  }
	



        //Compute the razor variables using the selected jets and possibly leptons
	/*
        vector<TLorentzVector> GoodPFObjects;
	//Razor Inclusive
	for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
	*/
        pfMET = metPt;
        vector<TLorentzVector> hemispheres = getHemispheres(GoodJets_uncorr);
        theMR = computeMR(hemispheres[0], hemispheres[1]); 
	
        //theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	//TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);

	//MuMu Box
        if(passedLeptonicTrigger && nLooseMuons > 1 && nLooseElectrons == 0 && nBTaggedJetsL == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
	  int it = 0;
	  for(auto& Mu : LooseMu){
	    PFMET.SetPx(PFMET.Px()+Mu.Px());
	    PFMET.SetPx(PFMET.Px()+Mu.Px());
	    MuonPt[it] = Mu.Pt();
	    MuonEta[it] = Mu.Eta();
	    MuonPhi[it] = Mu.Phi();
	    MuonE[it] = Mu.E();
	    it++;
	  }
	  pfMETcorr = PFMET.Pt();
	  theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = MuMu;
	      razorTree->Fill();
	    }
	    else razorBoxes["MuMu"]->Fill();
	  }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nLooseElectrons > 1 && nLooseMuons == 0 && nBTaggedJetsL == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
	  int it = 0;
          for(auto& Ele : LooseEle){
            PFMET.SetPx(PFMET.Px()+Ele.Px());
            PFMET.SetPx(PFMET.Px()+Ele.Px());
	    ElePt[it] = Ele.Pt();
	    EleEta[it] = Ele.Eta();
	    ElePhi[it] = Ele.Phi();
	    EleE[it] = Ele.E();
          }
          theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  pfMETcorr = PFMET.Pt();
	  if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = EleEle;
	      razorTree->Fill();
	    }
	    else razorBoxes["EleEle"]->Fill();
	  }
        }
	//Mu Box
        else if(passedLeptonicTrigger && nLooseMuons == 1 && nLooseElectrons == 0 && nBTaggedJetsL == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
          int it = 0;
          for(auto& Mu : LooseMu){
            PFMET.SetPx(PFMET.Px()+Mu.Px());
            PFMET.SetPx(PFMET.Px()+Mu.Px());
            MuonPt[it] = Mu.Pt();
            MuonEta[it] = Mu.Eta();
            MuonPhi[it] = Mu.Phi();
            MuonE[it] = Mu.E();
            it++;
          }
          theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  pfMETcorr = PFMET.Pt();
	  if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = MuJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["Mu"]->Fill();
	  }
        }
	//Ele Box
        else if(passedLeptonicTrigger && nLooseElectrons == 1 && nLooseMuons == 0 && nBTaggedJetsL == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
          int it = 0;
          for(auto& Ele : LooseEle){
            PFMET.SetPx(PFMET.Px()+Ele.Px());
            PFMET.SetPx(PFMET.Px()+Ele.Px());
            ElePt[it] = Ele.Pt();
            EleEta[it] = Ele.Eta();
            ElePhi[it] = Ele.Phi();
            EleE[it] = Ele.E();
          }
	  pfMETcorr = PFMET.Pt();
          theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = EleJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["Ele"]->Fill();
	  }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJetsL == 0 && nLooseMuons == 0 && nLooseElectrons == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
	  theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  pfMETcorr = PFMET.Pt();
	  if(passesHadronicRazorBaseline(theMR, theRsq)){  
	    if(combineTrees){
	      box = MultiJet;
	      razorTree->Fill();
	    }
	    else { razorBoxes["MultiJet"]->Fill(); }
	  }
        }
        //TwoBJet Box
        else if(passedHadronicTrigger && nBTaggedJetsT > 1 && nLooseMuons == 0 && nLooseElectrons == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
	  pfMETcorr = PFMET.Pt();
          theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = TwoBJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["TwoBJet"]->Fill();
	  }
        }
        //OneBJet Box
        else if(passedHadronicTrigger && nBTaggedJetsT == 1 && nLooseMuons == 0 && nLooseElectrons == 0){
	  TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);
	  pfMETcorr = PFMET.Pt();
          theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	    if(combineTrees){
	      box = OneBJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["OneBJet"]->Fill();
	  }
        }       

    }//end of event loop
    
    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();
    
    outFile.Close();
}
