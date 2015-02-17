#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

// enum RazorPhotonStudy_RazorBox {
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

bool RazorPhotonStudy_PassesHadronicRazorBaseline(double MR, double Rsq);
bool RazorPhotonStudy_PassesLeptonicRazorBaseline(double MR, double Rsq);

struct greater_than_pt{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
    return p1.Pt() > p2.Pt();
  }
};

void RazorAnalyzer::RazorPhotonStudy( string outputfilename, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorPhotonStudy.root";
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
    //    boxNames.push_back("SoftLeptonMultiJet");
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
    float theMR, MR_noPho, MR_noZ, MR_noW;
    float theRsq, Rsq_noPho, Rsq_noZ, Rsq_noW; 
    float met,genmet,met_noPho, met_noZ, met_noW, genHt, genZmass;
    float leadingGenMuonPt, leadingGenElectronPt, leadingGenPhotonPt;
    float leadingGenMuonEta, leadingGenElectronEta, leadingGenPhotonEta;
    float leadingMuonEta, leadingMuonPhi, leadingMuonPt;
    float j1pt, j2pt, j1eta, j2eta, j1phi, j2phi;
    float metphi, metphi_noZ, metphi_noW, metphi_noPho;
    float genZpt, recoZpt, genZeta, recoZeta, recoWpt, recoWeta;
    float minDRGenLeptonToGenParton;
    float leadingPhotonPt;
    int nSelectedPhotons;    
    // RazorPhotonStudy_RazorBox box;
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
        razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
        razorTree->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
        razorTree->Branch("leadingMuonEta", &leadingMuonEta, "leadingMuonEta/F");
        razorTree->Branch("leadingMuonPt", &leadingMuonPt, "leadingMuonPt/F");
        razorTree->Branch("leadingMuonPhi", &leadingMuonPhi, "leadingMuonPhi/F");
        razorTree->Branch("leadingGenElectronPt", &leadingGenElectronPt, "leadingGenElectronPt/F");
        razorTree->Branch("leadingGenPhotonPt", &leadingGenPhotonPt, "leadingGenPhotonPt/F");
        razorTree->Branch("leadingGenPhotonEta", &leadingGenPhotonEta, "leadingGenPhotonEta/F");
        razorTree->Branch("leadingPhotonPt", &leadingPhotonPt, "leadingPhotonPt/F");
        razorTree->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
        razorTree->Branch("genZpt", &genZpt, "genZpt/F");
        razorTree->Branch("genZeta", &genZeta, "genZeta/F");
        razorTree->Branch("genZmass", &genZmass, "genZmass/F");
        razorTree->Branch("genHt", &genHt, "genHt/F");
        razorTree->Branch("recoZpt", &recoZpt, "recoZpt/F");
        razorTree->Branch("recoZeta", &recoZeta, "recoZeta/F");
        razorTree->Branch("recoWpt", &recoWpt, "recoWpt/F");
        razorTree->Branch("recoWeta", &recoWeta, "recoWeta/F");
        razorTree->Branch("leadingGenElectronEta", &leadingGenElectronEta, "leadingGenElectronEta/F");
        razorTree->Branch("minDRGenLeptonToGenParton", &minDRGenLeptonToGenParton, "minDRGenLeptonToGenParton/F");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("MR_noZ", &MR_noZ, "MR_noZ/F");
        razorTree->Branch("MR_noW", &MR_noW, "MR_noW/F");
        razorTree->Branch("MR_noPho", &MR_noPho, "MR_noPho/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("Rsq_noPho", &Rsq_noPho, "Rsq_noPho/F");
        razorTree->Branch("Rsq_noZ", &Rsq_noZ, "Rsq_noZ/F");
        razorTree->Branch("Rsq_noW", &Rsq_noW, "Rsq_noW/F");
        razorTree->Branch("met", &met, "met/F");
        razorTree->Branch("metphi", &metphi, "metphi/F");
        razorTree->Branch("met_noPho", &met_noPho, "met_noPho/F");
        razorTree->Branch("met_noZ", &met_noZ, "met_noZ/F");
        razorTree->Branch("met_noW", &met_noW, "met_noW/F");
        razorTree->Branch("metphi_noZ", &metphi_noZ, "metphi_noZ/F");
        razorTree->Branch("metphi_noW", &metphi_noW, "metphi_noW/F");
        razorTree->Branch("metphi_noPho", &metphi_noPho, "metphi_noPho/F");
        razorTree->Branch("genmet", &genmet, "genmet/F");
        razorTree->Branch("box", &box, "box/I");
        razorTree->Branch("j1pt", &j1pt, "j1pt/F");
        razorTree->Branch("j2pt", &j2pt, "j2pt/F");
        razorTree->Branch("j1eta", &j1eta, "j1eta/F");
        razorTree->Branch("j2eta", &j2eta, "j2eta/F");
        razorTree->Branch("j1phi", &j1phi, "j1phi/F");
        razorTree->Branch("j2phi", &j2phi, "j2phi/F");
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
            box.second->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
	    box.second->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
	    box.second->Branch("leadingGenElectronPt", &leadingGenElectronPt, "leadingGenElectronPt/F");
	    box.second->Branch("leadingGenPhotonPt", &leadingGenPhotonPt, "leadingGenPhotonPt/F");
	    box.second->Branch("leadingGenPhotonEta", &leadingGenPhotonEta, "leadingGenPhotonEta/F");
	    box.second->Branch("leadingPhotonPt", &leadingPhotonPt, "leadingPhotonPt/F");
	    box.second->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
	    box.second->Branch("genZpt", &genZpt, "genZpt/F");
	    box.second->Branch("genZeta", &genZeta, "genZeta/F");
	    box.second->Branch("genZmass", &genZmass, "genZmass/F");
	    box.second->Branch("genHt", &genHt, "genHt/F");
	    box.second->Branch("recoZpt", &recoZpt, "recoZpt/F");
	    box.second->Branch("recoWpt", &recoWpt, "recoWpt/F");
	    box.second->Branch("recoZeta", &recoZeta, "recoZeta/F");
	    box.second->Branch("recoWeta", &recoWeta, "recoWeta/F");
	    box.second->Branch("leadingGenElectronEta", &leadingGenElectronEta, "leadingGenElectronEta/F");
	    box.second->Branch("minDRGenLeptonToGenParton", &minDRGenLeptonToGenParton, "minDRGenLeptonToGenParton/F");
            box.second->Branch("MR", &theMR, "MR/F");
            box.second->Branch("MR_noZ", &MR_noZ, "MR_noZ/F");
            box.second->Branch("MR_noW", &MR_noW, "MR_noW/F");
            box.second->Branch("MR_noPho", &MR_noPho, "MR_noPho/F");
	    box.second->Branch("met", &met, "met/F");
	    box.second->Branch("metphi", &metphi, "metphi/F");
	    box.second->Branch("genmet", &genmet, "genmet/F");
	    box.second->Branch("met_noPho", &met_noPho, "met_noPho/F");
	    box.second->Branch("met_noZ", &met_noZ, "met_noZ/F");
	    box.second->Branch("metphi_noZ", &metphi_noZ, "metphi_noZ/F");
	    box.second->Branch("met_noW", &met_noW, "met_noW/F");
	    box.second->Branch("metphi_noW", &metphi_noW, "metphi_noW/F");
	    box.second->Branch("metphi_noPho", &metphi_noPho, "metphi_noPho/F");
            box.second->Branch("Rsq", &theRsq, "Rsq/F");
	    box.second->Branch("Rsq_noZ", &Rsq_noZ, "Rsq_noZ/F");
	    box.second->Branch("Rsq_noW", &Rsq_noW, "Rsq_noW/F");
            box.second->Branch("Rsq_noPho", &Rsq_noPho, "Rsq_noPho/F");
	    box.second->Branch("j1pt", &j1pt, "j1pt/F");
	    box.second->Branch("j2pt", &j2pt, "j2pt/F");
	    box.second->Branch("j1eta", &j1eta, "j1eta/F");
	    box.second->Branch("j2eta", &j2eta, "j2eta/F");
	    box.second->Branch("j1phi", &j1phi, "j1phi/F");
	    box.second->Branch("j2phi", &j2phi, "j2phi/F");
        }
    }

    //begin loop
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //begin event
        if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
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
        nSelectedPhotons = 0;
        theMR = -1;
        MR_noZ = -1;
        MR_noW = -1;
        MR_noPho = -1;
	genZpt = -1;
	genZeta = -99;
	genZmass = -1;
	genHt = -1;
	recoZpt = -1;
	recoZeta = -999;
        recoWpt = -1;
        recoWeta = -999;
	met = metPt;
	met_noPho = -1.;
	met_noZ = -1.;
	met_noW = -1.;
	metphi_noZ = -99.;
	metphi_noW = -99.;
	metphi_noPho = -99.;
	metphi = -99.;
	leadingPhotonPt = -1;
	genmet = genMetPt;
        theRsq = -1;
        Rsq_noPho = -1;
        Rsq_noZ = -1;
        Rsq_noW = -1;
	j1pt=-1.;
	j2pt=-1.;
	j1eta=-99.;
	j2eta=-99.;
	j1phi=-99.;
	j2phi=-99.;
	
	minDRGenLeptonToGenParton = 9999;
        if(combineTrees) box = NONE;

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
	  //int closestLeptonIndex = -1;

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
		//closestLeptonIndex = k;
	      }
	    }
	  }

	  if ( minDRToGenLepton < minDRGenLeptonToGenParton) minDRGenLeptonToGenParton = minDRToGenLepton;
	}

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
            if(isMVANonTrigVetoElectron(i)) nVetoElectrons++;
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
            if(jetPt[i] < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //exclude selected muons and electrons from the jet collection
            double deltaR = -1;
            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            for(auto& lep : GoodLeptons){
                double thisDR = thisJet.DeltaR(lep);
                if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
            
            if(jetPt[i] > 80) numJetsAbove80GeV++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVM(i)){ 
                nBTaggedJets++;
            }
        }
	if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

	sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());
	j1pt=GoodJets[0].Pt();
	j2pt=GoodJets[1].Pt();
	j1eta=GoodJets[0].Eta();
	j2eta=GoodJets[1].Eta();
	j1phi=GoodJets[0].Phi();
	j2phi=GoodJets[1].Phi();
	
        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
	metphi = metPhi;

	if (GoodPFObjects.size() > 20) {
	  cout << "More than 20 objects : " <<GoodPFObjects.size() << " " << nVetoElectrons << " " << nVetoMuons << "\n"; 

	  if (nVetoMuons > 5) {
	    for(int i = 0; i < nMuons; i++){
	      cout << "Muon " << i << " : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << muonIsLoose[i] << " " << muon_ip3dSignificance[i] << " " << muon_relIso04DBetaCorr[i] << " : " << isVetoMuon(i) << " \n";	      
	    }
	  }
	}
	
	// compute R and MR
	vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	theMR = computeMR(hemispheres[0], hemispheres[1]); 
	theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	
        //photon selection
	vector<TLorentzVector> GoodPhotons;
	vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
	int nPhotonsAbove40GeV = 0;
	for(int i = 0; i < nPhotons; i++){
	  if(!isMediumPhoton(i)) continue;
	  
	  if(phoPt[i] < 25) continue;
	  if(fabs(phoEta[i]) > 1.479) continue;
	  
	  if(phoPt[i] > 40) nPhotonsAbove40GeV++;
	  TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], pho_RegressionE[i]);
	  GoodPhotons.push_back(thisPhoton);
	  GoodPhotonSigmaE.push_back(pho_RegressionEUncertainty[i]);
	  nSelectedPhotons++;
	}

	// found good reco photon
	if(GoodPhotons.size()>0){
	  sort(GoodPhotons.begin(), GoodPhotons.end(), greater_than_pt());
	  leadingPhotonPt = GoodPhotons.at(0).Pt();

	  // calculate met with removed photon
	  Double_t pho_px = GoodPhotons.at(0).Px();
	  Double_t pho_py = GoodPhotons.at(0).Py();
	  Double_t met_px = metPt*cos(metPhi);
	  Double_t met_py = metPt*sin(metPhi);
	  
	  met_noPho = TMath::Sqrt(pow(met_px + pho_px ,2) + pow(met_py + pho_py,2));
	  metphi_noPho = atan2( met_py+pho_py, met_px+pho_px );
	  
          //remove leading photon from jet collection
          vector<TLorentzVector> GoodPFObjectsNoLeadPhoton;
          for(auto& pf : GoodPFObjects){
              //get DR between jet and leading photon
              double thisDR = pf.DeltaR(GoodPhotons[0]);
              if(thisDR < 0.4) continue;
              
              GoodPFObjectsNoLeadPhoton.push_back(pf);
          }
          //remake the hemispheres using the new jet collection
          vector<TLorentzVector> hemispheresNoLeadPhoton = getHemispheres(GoodPFObjectsNoLeadPhoton);
	  TLorentzVector PFMET_NOPHO = makeTLorentzVectorPtEtaPhiM(met_noPho, 0, metphi_noPho, 0);
          MR_noPho = computeMR(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1]); 
          Rsq_noPho = computeRsq(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1], PFMET_NOPHO);

	  // find gen photon matching to leading reco photon
	  int matchedIndex = -1;
	  float minDR = 9999;
	  
	  for(int j = 0; j < nGenParticle; j++){
	    if (abs(gParticleId[j]) != 22) continue;	     
	    if (! ( (abs(gParticleMotherId[j]) >= 1 && abs(gParticleMotherId[j]) <= 5) || 
		    (abs(gParticleMotherId[j]) == 21) || 
		    (abs(gParticleMotherId[j]) == 2212) ) )  
	      continue;
	    if (abs(gParticleStatus[j]) != 1) continue;	     
	    
	    if ( deltaR( GoodPhotons.at(0).Eta(), GoodPhotons.at(0).Phi(), gParticleEta[j], gParticlePhi[j]) < 0.1
		 && deltaR( GoodPhotons.at(0).Eta(), GoodPhotons.at(0).Phi(), gParticleEta[j], gParticlePhi[j]) < minDR ){
	      matchedIndex = j;
	      minDR = deltaR( GoodPhotons.at(0).Eta(), GoodPhotons.at(0).Phi(), gParticleEta[j], gParticlePhi[j]);
	    }
	  }
	  
	  if (matchedIndex >= 0) {
	    leadingGenPhotonPt = gParticlePt[matchedIndex];
            leadingGenPhotonEta = gParticleEta[matchedIndex];
	  }	  
	}
	// Zll stuff
	vector<TLorentzVector> ZMuons; //leptons used to emulate Z->nn
        vector<TLorentzVector> WMuons; //leptons used to emulate W->ln
	for(int i = 0; i < nMuons; i++){
	  if(!isLooseMuon(i)) continue;  
	  if(muonPt[i] < 10) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
	  
	  TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	  ZMuons.push_back(thisMuon);

          if(!isTightMuon(i)) continue;
          WMuons.push_back(thisMuon);
	}
        //do the same for electrons
        vector<TLorentzVector> ZElectrons;
        vector<TLorentzVector> WElectrons;
        for(int i = 0; i < nElectrons; i++){
            if(!isLooseElectron(i)) continue;
            if(elePt[i] < 10) continue;
            if(abs(eleEta[i]) > 2.5) continue;

            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            ZElectrons.push_back(thisElectron);

            if(!isTightElectron(i)) continue;
            WElectrons.push_back(thisElectron);
        }

	// reco level Z pt
	if(ZElectrons.size() == 0 && ZMuons.size()==2)
	  {
              TLorentzVector m1 = ZMuons[0];
              TLorentzVector m2 = ZMuons[1];
              TLorentzVector theZ_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), (m1+m2).M());
	    
	    Double_t met_px = metPt*cos(metPhi) + theZ_perp.Px();
	    Double_t met_py = metPt*sin(metPhi) + theZ_perp.Py();

	    met_noZ = TMath::Sqrt( pow(met_px, 2) + pow(met_py, 2) );
	    metphi_noZ = atan2( met_py, met_px );
	    
	    vector<TLorentzVector> GoodPFObjectsNoZ;
	    for(auto& jet : GoodJets) GoodPFObjectsNoZ.push_back(jet);

	    vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodPFObjectsNoZ);

	    TLorentzVector PFMET_NOZ = makeTLorentzVectorPtEtaPhiM(met_noZ, 0, metphi_noZ, 0);
	    Rsq_noZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], PFMET_NOZ);
	    MR_noZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
	    
	    recoZpt = (m1+m2).Pt();
	    recoZeta = (m1+m2).Eta();
	  }
        else if(ZElectrons.size()==2 && ZMuons.size() == 0)
        {
            TLorentzVector m1 = ZElectrons[0];
            TLorentzVector m2 = ZElectrons[1];
            TLorentzVector theZ_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), (m1+m2).M());

            Double_t met_px = metPt*cos(metPhi) + theZ_perp.Px();
            Double_t met_py = metPt*sin(metPhi) + theZ_perp.Py();

            met_noZ = TMath::Sqrt( pow(met_px, 2) + pow(met_py, 2) );
            metphi_noZ = atan2( met_py, met_px );

            vector<TLorentzVector> GoodPFObjectsNoZ;
            for(auto& jet : GoodJets) GoodPFObjectsNoZ.push_back(jet);

            vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodPFObjectsNoZ);

            TLorentzVector PFMET_NOZ = makeTLorentzVectorPtEtaPhiM(met_noZ, 0, metphi_noZ, 0);
            Rsq_noZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], PFMET_NOZ);
            MR_noZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 

            recoZpt = (m1+m2).Pt();
            recoZeta = (m1+m2).Eta();
        }
	// reco level W pt
        else if(WMuons.size() == 1 && ZElectrons.size() == 0 && ZMuons.size() == 1) //one tight muon, no extra loose muons, and no loose electrons
        {
            TLorentzVector m1 = WMuons[0];
            TLorentzVector m2 = PFMET;
            double deltaPhiLepMet = m1.DeltaPhi(m2);
            double mTLepMet = sqrt(2*m2.Pt()*m1.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) );
            if(mTLepMet > 30 && mTLepMet < 100){
                TLorentzVector theW_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), (m1+m2).M());

                Double_t met_px = theW_perp.Px();
                Double_t met_py = theW_perp.Py();

                met_noW = TMath::Sqrt( pow(met_px, 2) + pow(met_py, 2) );
                metphi_noW = atan2( met_py, met_px );

                vector<TLorentzVector> GoodPFObjectsNoW;
                for(auto& jet : GoodJets) GoodPFObjectsNoW.push_back(jet);

                vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodPFObjectsNoW);
                TLorentzVector PFMET_NOW = makeTLorentzVectorPtEtaPhiM(met_noW, 0, metphi_noW, 0);
                Rsq_noW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], PFMET_NOW);
                MR_noW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 

                recoWpt = (m1+m2).Pt();
                recoWeta = (m1+m2).Eta();
            }
        }
        else if(WElectrons.size() == 1 && ZElectrons.size()==1 && ZMuons.size() == 0) //one tight ele, no extra loose eles,and no loose muons
        {
            TLorentzVector m1 = WElectrons[0];
            TLorentzVector m2 = PFMET;
            double deltaPhiLepMet = m1.DeltaPhi(m2);
            double mTLepMet = sqrt(2*m2.Pt()*m1.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) );
            if(mTLepMet > 30 && mTLepMet < 100){
                TLorentzVector theW_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), (m1+m2).M());

                Double_t met_px = theW_perp.Px();
                Double_t met_py = theW_perp.Py();

                met_noW = TMath::Sqrt( pow(met_px, 2) + pow(met_py, 2) );
                metphi_noW = atan2( met_py, met_px );

                vector<TLorentzVector> GoodPFObjectsNoW;
                for(auto& jet : GoodJets) GoodPFObjectsNoW.push_back(jet);

                vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodPFObjectsNoW);

                TLorentzVector PFMET_NOW = makeTLorentzVectorPtEtaPhiM(met_noW, 0, metphi_noW, 0);
                Rsq_noW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], PFMET_NOW);
                MR_noW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 

                recoWpt = (m1+m2).Pt();
                recoWeta = (m1+m2).Eta();
            }
        }


	// lepton efficiency
	
	//
	
	
	// gen level Z pt
	for(int j = 0; j < nGenParticle; j++){
	  if (abs(gParticleId[j]) != 23) continue;

	  // cout<<"Z boson "<<jentry<<" "<<gParticleId[j] << " " << gParticleStatus[j] << " : " <<  gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";

	  TLorentzVector Zboson = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]); 
	  
	  if(gParticleStatus[j]==22)
	    {
	      genZpt = gParticlePt[j];
	      genZeta = gParticleEta[j];
	      genZmass = Zboson.M();
	    }
	}

	// calculate gen Ht
	for(int j = 0; j < nGenParticle; j++){
	  
	  if ( (abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) ||  (abs(gParticleId[j]) == 21) )
	    if (gParticleStatus[j] == 23 )
	  genHt += gParticlePt[j];
	}

        //MuEle Box
        if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseMuons > 0 && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuEle;
                    razorTree->Fill();
                }
                else razorBoxes["MuEle"]->Fill();
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMu;
                    razorTree->Fill();
                }
                else razorBoxes["MuMu"]->Fill();
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleEle;
                    razorTree->Fill();
                }
                else razorBoxes["EleEle"]->Fill();
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuMultiJet"]->Fill();
            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuJet"]->Fill();
            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleMultiJet"]->Fill();
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleJet"]->Fill();
            }
        }

	// //Soft Lepton + MultiJet Box
        // else if(passedHadronicTrigger && nVetoElectrons + nVetoMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
        //     if(RazorPhotonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){  
        //         if(combineTrees){
        //             box = SoftLeptonMultiJet;
        //             razorTree->Fill();
        //         }
        //         else razorBoxes["SoftLeptonMultiJet"]->Fill();
        //     }
        // }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
            if(RazorPhotonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){  
                if(combineTrees){
                    box = MultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["MultiJet"]->Fill();
            }
        }
        //2BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 1){
            if(RazorPhotonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = TwoBJet;
                    razorTree->Fill();
                }
                else razorBoxes["TwoBJet"]->Fill();
            }
        }
        //1BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0){
            if(RazorPhotonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = OneBJet;
                    razorTree->Fill();
                }
                else razorBoxes["OneBJet"]->Fill();
            }
        }
        //0BJetBox
        else if(passedHadronicTrigger){
            if(RazorPhotonStudy_PassesHadronicRazorBaseline(theMR, theRsq)){ 
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

bool RazorPhotonStudy_PassesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    if(MR < 0 || Rsq < 0) passes = false;
    return passes;
}

bool RazorPhotonStudy_PassesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    // if(MR < 300 || Rsq < 0.15) passes = false;
    // if(MR < 350 && Rsq < 0.2) passes = false;
    if(MR < 0 || Rsq < 0) passes = false;
    return passes;
}
