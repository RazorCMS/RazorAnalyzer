
#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

// enum RazorVetoLeptonStudy_RazorBox {
//     MuEle, 
//     MuMu,
//     EleEle,
//     MuMultiJet,
//     MuJet,
//     EleMultiJet,
//     EleJet,
//     SoftLeptonMultiJet,
//     MultiJet,
//     TwoBJet,
//     OneBJet,
//     ZeroBJet,
//     NONE
// };

bool RazorVetoLeptonStudy_PassesHadronicRazorBaseline(double MR, double Rsq);
bool RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(double MR, double Rsq);

void RazorAnalyzer::RazorVetoLeptonStudy( string outputfilename, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    std::vector<JetCorrectorParameters> correctionParameters;
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt"));
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L2Relative_AK4PFchs.txt"));
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt"));    
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorVetoLeptonStudy.root";
    TFile outFile(outfilename.c_str(), "RECREATE");
    
    //one tree to hold all events
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");
    
    //separate trees for individual boxes
    map<string, TTree*> razorBoxes;
    vector<string> boxNames;
    boxNames.push_back("MuEle");
    boxNames.push_back("MuMu");
    boxNames.push_back("EleEle");
    boxNames.push_back("MuMultiJet");
    boxNames.push_back("MuJet");
    boxNames.push_back("EleMultiJet");
    boxNames.push_back("EleJet");
    boxNames.push_back("SoftLeptonMultiJet");
    boxNames.push_back("MultiJet");
    boxNames.push_back("TwoBJet");
    boxNames.push_back("OneBJet");
    boxNames.push_back("ZeroBJet");
    for(size_t i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nBTaggedJets;
    int nVetoMuons, nLooseMuons, nTightMuons;
    int nVetoElectrons, nLooseElectrons, nTightElectrons;
    int nLooseTaus, nMediumTaus, nTightTaus;
    int nVetoMVAElectrons;
    int nGenMuons, nGenElectrons, nGenTauMuons, nGenTauElectrons, nGenTaus;
    float theMR;
    float theRsq;
    float met;
    float leadingGenMuonPt, leadingGenElectronPt;
    float leadingGenMuonEta, leadingGenElectronEta;
    float minDRGenLeptonToGenParton;
    int npu;
    // RazorVetoLeptonStudy_RazorBox box;
    RazorBox box;

    //set branches on big tree
    if(combineTrees){
      razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
        razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
        razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nVetoMVAElectrons", &nVetoMVAElectrons, "nVetoMVAElectrons/I");
        razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
        razorTree->Branch("nMediumTaus", &nMediumTaus, "nMediumTaus/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
        razorTree->Branch("nGenElectrons", &nGenElectrons, "nGenElectrons/I");
        razorTree->Branch("nGenTauMuons", &nGenTauMuons, "nGenTauMuons/I");
        razorTree->Branch("nGenTauElectrons", &nGenTauElectrons, "nGenTauElectrons/I");
        razorTree->Branch("nGenTaus", &nGenTaus, "nGenTaus/I");
        razorTree->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
        razorTree->Branch("leadingGenElectronPt", &leadingGenElectronPt, "leadingGenElectronPt/F");
        razorTree->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
        razorTree->Branch("leadingGenElectronEta", &leadingGenElectronEta, "leadingGenElectronEta/F");
        razorTree->Branch("minDRGenLeptonToGenParton", &minDRGenLeptonToGenParton, "minDRGenLeptonToGenParton/F");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("met", &met, "met/F");
        razorTree->Branch("box", &box, "box/I");
        razorTree->Branch("npu", &npu, "npu/I");
    }
    //set branches on all trees
    else{ 
        for(auto& box : razorBoxes){
            box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
            box.second->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
	    box.second->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
            box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
            box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
	    box.second->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
	    box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
            box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
	    box.second->Branch("nVetoMVAElectrons", &nVetoMVAElectrons, "nVetoMVAElectrons/I");
            box.second->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
            box.second->Branch("nMediumTaus", &nMediumTaus, "nMediumTaus/I");
            box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
	    box.second->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
	    box.second->Branch("nGenElectrons", &nGenElectrons, "nGenElectrons/I");
	    box.second->Branch("nGenTauMuons", &nGenTauMuons, "nGenTauMuons/I");
	    box.second->Branch("nGenTauElectrons", &nGenTauElectrons, "nGenTauElectrons/I");
	    box.second->Branch("nGenTaus", &nGenTaus, "nGenTaus/I");
	    box.second->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
	    box.second->Branch("leadingGenElectronPt", &leadingGenElectronPt, "leadingGenElectronPt/F");
	    box.second->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
	    box.second->Branch("leadingGenElectronEta", &leadingGenElectronEta, "leadingGenElectronEta/F");
	    box.second->Branch("minDRGenLeptonToGenParton", &minDRGenLeptonToGenParton, "minDRGenLeptonToGenParton/F");
            box.second->Branch("MR", &theMR, "MR/F");
	    box.second->Branch("met", &met, "met/F");
            box.second->Branch("Rsq", &theRsq, "Rsq/F");
	    box.second->Branch("npu", &npu, "npu/I");
        }
    }

    //begin loop
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //begin event
        if(jentry % 100 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //fill normalization histogram
        NEvents->Fill(1.0);

        //reset tree variables
        nSelectedJets = 0;
        nBTaggedJets = 0;
        nVetoMuons = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nVetoMVAElectrons = 0;
        nVetoElectrons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nLooseTaus = 0;
        nMediumTaus = 0;
        nTightTaus = 0;
	nGenMuons = 0;
	nGenElectrons = 0;
	nGenTauMuons = 0;
	nGenTauElectrons = 0;
	nGenTaus = 0;
        theMR = -1;
	met = metPt;
        theRsq = -1;
	minDRGenLeptonToGenParton = 9999;
        if(combineTrees) box = NONE;
	npu = -1;

	//get NPU
	for (int i=0; i < nBunchXing; ++i) {
	  if (BunchXing[i] == 0) {
	    npu = nPU[i];
	  }
	}

        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        
	//count generated leptons
	leadingGenMuonPt = 0;
	leadingGenElectronPt = 0;
	leadingGenMuonEta = -999;
	leadingGenElectronEta = -999;
	for(int j = 0; j < nGenParticle; j++){
	  if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	      && abs(gParticleEta[j]) < 2.5 && gParticlePt[j] > 5
	      ) {
	    if (  (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23) ) {
	      nGenElectrons++;
	      if (gParticlePt[j] > leadingGenElectronPt) {
		leadingGenElectronPt = gParticlePt[j];
		leadingGenElectronEta = gParticleEta[j];
	      }
	    }
	    if ( abs(gParticleMotherId[j]) == 15) {
	      nGenTauElectrons++;
	    }
	  }
	  if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	      && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 5
	      ) {
	    if ( (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)) {
	      nGenMuons++;
	      if (gParticlePt[j] > leadingGenMuonPt) {
		leadingGenMuonPt = gParticlePt[j];
		leadingGenMuonEta = gParticleEta[j];
	      }
	    }
	    if ( abs(gParticleMotherId[j]) == 15) {
	      nGenTauMuons++;
	    }
	  }
	  if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
	      && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
	      && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 20
	      ) nGenTaus++;
	}
	
	for(int j = 0; j < nGenParticle; j++){
	  float minDRToGenLepton = 9999;
	  int closestLeptonIndex = -1;

	  //only look for outgoing partons
	  if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		  && gParticleStatus[j] == 23)
	       ) continue;
	       		      
	  //look for closest lepton
	  for(int k = 0; k < nGenParticle; k++){
	    if ( 
		(abs(gParticleId[k]) == 11 && gParticleStatus[k] == 1 	      
	      && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
		 && abs(gParticleEta[k]) < 2.5 && gParticlePt[k] > 5)
		||
		(abs(gParticleId[k]) == 13 && gParticleStatus[k] == 1  
	      && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
	      && abs(gParticleEta[k]) < 2.4 && gParticlePt[k] > 5
		 )
		||
		(abs(gParticleId[k]) == 15 && gParticleStatus[k] == 2 
	      && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
	      && abs(gParticleEta[k]) < 2.4 && gParticlePt[k] > 20
		 )
		 ) {	 
	      double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], gParticleEta[k], gParticlePhi[k]);
	      if ( tmpDR < minDRToGenLepton ) {
		minDRToGenLepton = tmpDR;
		closestLeptonIndex = k;
	      }
	    }
	  }

	  // cout << "Parton " << j << " : " << gParticleId[j] << " | " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " | " << gParticleMotherId[j] 
	  //      << " : " << minDRToGenLepton 
	  //      << " | " ;
	  // if (closestLeptonIndex >= 0 && minDRToGenLepton < 0.4) {
	  //   cout << gParticleId[closestLeptonIndex] << " " << gParticlePt[closestLeptonIndex] << " " << gParticleEta[closestLeptonIndex] << " " << gParticlePhi[closestLeptonIndex] << " " << gParticleMotherId[closestLeptonIndex];
	  // }
	  // cout << "\n";

	  if ( minDRToGenLepton < minDRGenLeptonToGenParton) minDRGenLeptonToGenParton = minDRToGenLepton;
	}

	
	//if (nGenTaus+nGenMuons+nGenElectrons == 0) {
	// if (minDRGenLeptonToGenParton < 0.4) {
	//   cout << "\n\nDEBUG\n";
	//   cout << nGenElectrons << " " << nGenMuons << " " << nGenTaus << "\n";
	//   cout << "minDRGenLeptonToGenParton : " << minDRGenLeptonToGenParton << "\n";
	//   for(int j = 0; j < nGenParticle; j++){
	//     cout << "particle " << j <<  " : " << gParticleId[j] << " " << gParticleStatus[j] << " : " <<  gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";
	//   }
	// }

        vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
        for(int i = 0; i < nMuons; i++){

            if(muonPt[i] < 5) continue;
            if(abs(muonEta[i]) > 2.4) continue;

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

	    if(isMVANonTrigVetoElectron(i)) nVetoMVAElectrons++;
            if(isMVANonTrigVetoElectron(i)) {
	      nVetoElectrons++;
	    }
            if(isLooseElectron(i) && elePt[i] > 10 ) nLooseElectrons++;
            if(isTightElectron(i) && elePt[i] > 10 ) nTightElectrons++;

	    if(!isMVANonTrigVetoElectron(i)) continue; 
	    TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            GoodLeptons.push_back(thisElectron);        
        }

        for(int i = 0; i < nTaus; i++){

	  if(isLooseTau(i)){
	    nLooseTaus++;
	  }
	  if(isMediumTau(i)){
	    nMediumTaus++;
	  }
	  if(isTightTau(i)){
	    nTightTaus++;
	  }
	    
        }
        
        vector<TLorentzVector> GoodJets;
        int numJetsAbove80GeV = 0;
        for(int i = 0; i < nJets; i++){

	  double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						 fixedGridRhoFastjetAll, jetJetArea[i], 
						 JetCorrector);   
	  
	  if(jetPt[i]*JEC < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //exclude selected muons and electrons from the jet collection
            double dR = -1;
            TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC);
            for(auto& lep : GoodLeptons){
                double thisDR = thisJet.DeltaR(lep);
                if(dR < 0 || thisDR < dR) dR = thisDR;
            }
            if(dR > 0 && dR < 0.4) continue; //jet matches a selected lepton


	    // //exclude jet if it doesn't match a selected genjet
            // bool matchedGenJet = false;
	    // for(int j = 0; j < nGenJets; j++){
	    //   double thisDR = deltaR(genJetEta[j],genJetPhi[j],jetEta[i],jetPhi[i]);
	    //   if(thisDR < 0.4 && fabs(jetPt[i]-genJetPt[j])/genJetPt[j] < 0.5 ){
	    // 	  matchedGenJet = true;
	    // 	  break;
            //     }
            // }
            // if(!matchedGenJet) continue;

        
            if(jetPt[i]*JEC > 80) numJetsAbove80GeV++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVM(i)){ 
                nBTaggedJets++;
            }
        }
        if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);


	if (GoodPFObjects.size() > 20) {
	  cout << "More than 20 objects : " <<GoodPFObjects.size() << " " << nVetoElectrons << " " << nVetoMuons << "\n"; 

	  if (nVetoMuons > 5) {
	    for(int i = 0; i < nMuons; i++){
	      cout << "Muon " << i << " : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << muonIsLoose[i] << " " << muon_ip3dSignificance[i] << " " << muon_relIso04DBetaCorr[i] << " : " << isVetoMuon(i) << " \n";	      
	    }
	  }

	}

	vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	theMR = computeMR(hemispheres[0], hemispheres[1]); 
	theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);

        //MuEle Box
        if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseMuons > 0 && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuEle;
                    razorTree->Fill();
                }
                else razorBoxes["MuEle"]->Fill();
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMu;
                    razorTree->Fill();
                }
                else razorBoxes["MuMu"]->Fill();
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleEle;
                    razorTree->Fill();
                }
                else razorBoxes["EleEle"]->Fill();
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuMultiJet"]->Fill();
            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuJet"]->Fill();
            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleMultiJet"]->Fill();
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleJet"]->Fill();
            }
        }

	//Soft Lepton + MultiJet Box
        else if(passedHadronicTrigger && nVetoElectrons + nVetoMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorVetoLeptonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){  
                if(combineTrees){
                    box = SoftLeptonMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["SoftLeptonMultiJet"]->Fill();
            }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorVetoLeptonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){  
                if(combineTrees){
                    box = MultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["MultiJet"]->Fill();
            }
        }
        //2BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 1){
            if(RazorVetoLeptonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = TwoBJet;
                    razorTree->Fill();
                }
                else razorBoxes["TwoBJet"]->Fill();
            }
        }
        //1BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0){
            if(RazorVetoLeptonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = OneBJet;
                    razorTree->Fill();
                }
                else razorBoxes["OneBJet"]->Fill();
            }
        }
        //0BJetBox
        else if(passedHadronicTrigger){
            if(RazorVetoLeptonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = ZeroBJet;
                    razorTree->Fill();
                }
                else razorBoxes["ZeroBJet"]->Fill();
            }
        }
    }//end of event loop

    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();

    outFile.Close();
}

bool RazorVetoLeptonStudy_PassesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    //temporarily disable these
    // if(MR < 400 || Rsq < 0.25) passes = false;
    // if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool RazorVetoLeptonStudy_PassesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    // if(MR < 300 || Rsq < 0.15) passes = false;
    // if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}
