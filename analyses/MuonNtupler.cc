#define MuonNtupler_cxx
#include "RazorAnalyzer.h"
#include "MuonTree.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorAnalyzer::MuonNtupler(string outputfilename , int Option)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "MuonNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    MuonTree *muTree = new MuonTree;
    muTree->CreateTree(MuonTree::kMuTreeLight);
    muTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";
    
    UInt_t NMuonsFilled = 0;
 
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
      
      bool printDecay = false;
      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {

	for(int i = 0; i < nMuons; i++){
	  if(muonPt[i] < 5) continue;

	  //***********************
	  //Fill Muon Variables
	  //***********************
   
	  muTree->fWeight = 1;
	  muTree->fRunNumber = runNum;
	  muTree->fLumiSectionNumber = lumiNum;
	  muTree->fEventNumber = eventNum;
	  muTree->fMuEventNumberParity = (eventNum % 2 == 0);
	  muTree->fMuPt = muonPt[i]; 
	  muTree->fMuEta = muonEta[i]; 
	  muTree->fMuPhi = muonPhi[i]; 
	  muTree->fRho = 0; 
	  muTree->fNVertices = nPV; 
	  muTree->fMuD0 = muon_d0[i]; 
	  muTree->fMuIP3d = muon_ip3d[i];
	  muTree->fMuIP3dSig = muon_ip3dSignificance[i];
	  muTree->fMuPFIso04 = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
	  muTree->fMuIsLoose = muonIsLoose[i];
	  muTree->fMuIsTight = muonIsTight[i];
	  muTree->fPassVetoSelection = isVetoMuon(i);
	  muTree->fPassLooseSelection = isLooseMuon(i);
	  muTree->fPassTightSelection = isTightMuon(i);
	  muTree->fPtRel = muon_ptrel[i];
	  muTree->fMiniIso = muon_miniiso[i];

	  //Match to Gen particle
	  int matchedIndex = -1;
	  float minDR = 9999;

	  for(int j = 0; j < nGenParticle; j++){
	    if (abs(gParticleId[j]) != 13) continue;	      
	    if ( deltaR( muonEta[i], muonPhi[i], gParticleEta[j], gParticlePhi[j]) < 0.1
		 && deltaR( muonEta[i], muonPhi[i], gParticleEta[j], gParticlePhi[j]) < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = deltaR( muonEta[i], muonPhi[i], gParticleEta[j], gParticlePhi[j]);
	    }
	  }

	  int matchedID = 0;
	  if (matchedIndex >= 0) {
	    muTree->fMuGenPt = gParticlePt[matchedIndex];
	    muTree->fMuGenEta = gParticleEta[matchedIndex];
	    muTree->fMuGenPhi = gParticlePhi[matchedIndex];
	    if (gParticleMotherId[matchedIndex] > 50 || 
		abs(gParticleMotherId[matchedIndex]) == 15) {
	      matchedID = gParticleMotherId[matchedIndex];
	    } else if (abs(gParticleMotherId[matchedIndex]) == 23 || abs(gParticleMotherId[matchedIndex]) == 24 ) {
	      matchedID = gParticleId[matchedIndex];
	    }
	  }
	  muTree->fPdgId = matchedID;

	  //select only fakes
	  if (Option == 0) {
	    if (!(matchedID == 0 || abs(matchedID) > 50)) continue;
	  }
	  //select only real prompt
	  if (Option == 1) {
	    if (!(abs(matchedID) == 13)) continue;
	  }	   


	  //Find Closest Parton
	  float minDRToParton = 9999;
	  int closestPartonIndex = -1;
	  for(int j = 0; j < nGenParticle; j++){
	      
	    //only look for outgoing partons
	    if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		    && gParticleStatus[j] == 23)
		 ) continue;
	      
	    double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], muonEta[i], muonPhi[i]);
	    if ( tmpDR < minDRToParton ) {
	      closestPartonIndex = j;
	      minDRToParton = tmpDR;
	    }
	  }
	  muTree->fDRToClosestParton = minDRToParton;

	  if (closestPartonIndex > 0 || closestPartonIndex <= 0) {
	  // if (minDRToParton < 0.4) {
	  //   printDecay = true;
	  //   cout << "Muon Close to Parton : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " --> " << gParticleId[closestPartonIndex] << " " << gParticlePt[closestPartonIndex] << " " << gParticleEta[closestPartonIndex] << " " << gParticlePhi[closestPartonIndex] << "\n";
	  // }
	  }

	  //***********************
	  //Fill Muon
	  //***********************
	  NMuonsFilled++;
	  muTree->tree_->Fill();
	}

	if (printDecay) {
	  cout << "Print Decay Table\n";
	  for(int j = 0; j < nGenParticle; j++){
	    cout << j << " " <<  gParticleId[j] << " " <<  gParticleStatus[j] << " | " <<  gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " | " << gParticleMotherId[j] << " " << gParticleMotherIndex[j] << "\n";
	  }
	  cout << "\n\n";
	}


      }


      //********************************************
      //Tree entries based on gen-level objects
      //********************************************
      else if (Option >= 10 && Option < 20) {
	  
	for(int i = 0; i < nGenParticle; i++){
	    
	  //select prompt muons
	  if (Option == 11) {
	    if (abs(gParticleId[i]) != 13) continue;	      
	    if (!(abs(gParticleMotherId[i]) == 23 || 
		  abs(gParticleMotherId[i]) == 24)
		) continue;
	  }

	  if(gParticlePt[i] < 5) continue;

	  //***********************
	  //Fill Muon Variables
	  //***********************
   
	  muTree->fWeight = 1;
	  muTree->fRunNumber = runNum;
	  muTree->fLumiSectionNumber = lumiNum;
	  muTree->fEventNumber = eventNum;
	  muTree->fMuEventNumberParity = (eventNum % 2 == 0);
	  muTree->fMuGenPt = gParticlePt[i];
	  muTree->fMuGenEta = gParticleEta[i];
	  muTree->fMuGenPhi = gParticlePhi[i];
	  muTree->fRho = 0; 
	  muTree->fNVertices = nPV; 
	  muTree->fPdgId = gParticleId[i];

	  //Find Closest Parton
	  float minDRToParton = 9999;
	  for(int j = 0; j < nGenParticle; j++){
	      
	    //only look for outgoing partons
	    if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		    && gParticleStatus[j] == 23)
		 ) continue;
	      
	    double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], gParticleEta[i], gParticlePhi[i]);
	    if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	  }
	  muTree->fDRToClosestParton = minDRToParton;


	  //Find Associated Reco Muon
	  int matchedIndex = -1;
	  float minDR = 9999;
	  for(int j = 0; j < nMuons; j++){
	    if ( deltaR( muonEta[j], muonPhi[j], gParticleEta[i], gParticlePhi[i]) < 0.1
		 && deltaR( muonEta[j], muonPhi[j], gParticleEta[i], gParticlePhi[i]) < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = deltaR( muonEta[j], muonPhi[j], gParticleEta[i], gParticlePhi[i]);
	    }	    
	  }	 

	  if (matchedIndex >= 0) {
	    muTree->fMuPt = muonPt[matchedIndex]; 
	    muTree->fMuEta = muonEta[matchedIndex]; 
	    muTree->fMuPhi = muonPhi[matchedIndex]; 
	    muTree->fMuD0 = muon_d0[matchedIndex]; 
	    muTree->fMuIP3d = muon_ip3d[matchedIndex];
	    muTree->fMuIP3dSig = muon_ip3dSignificance[matchedIndex];
	    muTree->fMuPFIso04 = (muon_chargedIso[matchedIndex] + fmax(0.0,  muon_photonIso[matchedIndex] + muon_neutralHadIso[matchedIndex] - 0.5*muon_pileupIso[matchedIndex])) / muonPt[matchedIndex];
	    muTree->fMuIsLoose = muonIsLoose[matchedIndex];
	    muTree->fMuIsTight = muonIsTight[matchedIndex];
	    muTree->fPassVetoSelection = isVetoMuon(matchedIndex);
	    muTree->fPassLooseSelection = isLooseMuon(matchedIndex);
	    muTree->fPassTightSelection = isTightMuon(matchedIndex);
	    muTree->fPtRel = muon_ptrel[matchedIndex];
	    muTree->fMiniIso = muon_miniiso[matchedIndex];
	  } else {
	    muTree->fMuPt = 0;
	    muTree->fMuEta = 0;
	    muTree->fMuPhi = 0;
	    muTree->fMuD0 = 0;
	    muTree->fMuIP3d = 0;
	    muTree->fMuIP3dSig = 0;
	    muTree->fMuPFIso04 = -1;
	    muTree->fMuIsLoose = false;
	    muTree->fMuIsTight = false;
	    muTree->fPassVetoSelection = false;
	    muTree->fPassLooseSelection = false;
	    muTree->fPassTightSelection = false;
	    muTree->fPtRel = 0;
	    muTree->fMiniIso = 0;
	  }
	  
	  //***********************
	  //Fill Muon
	  //***********************
	  NMuonsFilled++;
	  muTree->tree_->Fill();


	}
	  
      }
	
    }//end of event loop

    cout << "Filled Total of " << NMuonsFilled << " Muons\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}



