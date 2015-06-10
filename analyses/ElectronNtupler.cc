#define ElectronNtupler_cxx
#include "RazorAnalyzer.h"
#include "ElectronTree.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorAnalyzer::ElectronNtupler(string outputfilename , int Option)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "ElectronNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ElectronTree *eleTree = new ElectronTree;
    eleTree->CreateTree(ElectronTree::kEleTreeLight);
    eleTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";
    
    UInt_t NElectronsFilled = 0;
 
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


      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {

        for(int i = 0; i < nElectrons; i++){
            if(elePt[i] < 5) continue;

	    //***********************
	    //Fill Electron Variables
	    //***********************
   
	    eleTree->fWeight = 1;
	    eleTree->fRunNumber = runNum;
	    eleTree->fLumiSectionNumber = lumiNum;
	    eleTree->fEventNumber = eventNum;
	    eleTree->fEleEventNumberParity = (eventNum % 2 == 0);
	    eleTree->fCharge = eleCharge[i] ;
	    eleTree->fElePt = elePt[i]; 
	    eleTree->fEleEta = eleEta[i]; 
	    eleTree->fElePhi = elePhi[i]; 
	    eleTree->fEleSCEta = eleEta_SC[i]; 
	    eleTree->fEleTriggerBit = 0;
	    eleTree->fRho = 0; 
	    eleTree->fNVertices = nPV; 
	    eleTree->fEleD0 = ele_d0[i]; 
	    eleTree->fEleDZ = ele_dZ[i]; 
	    eleTree->fElePassConversion = ele_PassConvVeto[i];
	    eleTree->fEleNMissHits =ele_MissHits[i];
	    eleTree->fEleOneOverEMinusOneOverP = ele_OneOverEminusOneOverP[i];
	    eleTree->fEleDEtaIn = ele_dEta[i];					   
	    eleTree->fEleDPhiIn = ele_dPhi[i];
	    eleTree->fEleSigmaIEtaIEta = eleFull5x5SigmaIetaIeta[i];
	    eleTree->fEleR9 = eleR9[i];
	    eleTree->fEleHoverE = ele_HoverE[i];
	    eleTree->fElePFIso04 = (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i];
	    eleTree->fIDMVATrig = ele_IDMVATrig[i];
	    eleTree->fIDMVANonTrig = ele_IDMVANonTrig[i];
	    eleTree->fPassVetoSelection = isVetoElectron(i);
	    eleTree->fPassLooseSelection = isRunOneLooseElectron(i);
	    eleTree->fPassTightSelection = isRunOneTightElectron(i);
	    eleTree->fPassMVANonTrigVetoSelection = isMVANonTrigVetoElectron(i);
	    eleTree->fPtRel = ele_ptrel[i];
	    eleTree->fMiniIso = ele_miniiso[i];

	    //Match to Gen particle
	    int matchedIndex = -1;
	    float minDR = 9999;

	    for(int j = 0; j < nGenParticle; j++){
	      if (abs(gParticleId[j]) != 11) continue;	      
	      if ( deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) < 0.1
		   && deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) < minDR
		   ) {		
		matchedIndex = j;
		minDR = deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]);
	      }
	    }

	    int matchedID = 0;
	    if (matchedIndex >= 0) {
	      eleTree->fEleGenPt = gParticlePt[matchedIndex];
	      eleTree->fEleGenEta = gParticleEta[matchedIndex];
	      eleTree->fEleGenPhi = gParticlePhi[matchedIndex];
	      if (gParticleMotherId[matchedIndex] > 50 || 
		  abs(gParticleMotherId[matchedIndex]) == 15) {
		matchedID = gParticleMotherId[matchedIndex];
	      } else if (abs(gParticleMotherId[matchedIndex]) == 23 || abs(gParticleMotherId[matchedIndex]) == 24 ) {
		matchedID = gParticleId[matchedIndex];
	      }
	    }
	    eleTree->fPdgId = matchedID;

	    //select only fakes
	    if (Option == 0) {
	      if (!(matchedID == 0 || abs(matchedID) > 50)) continue;
	      //if (abs(matchedID) == 11) continue;
	    }
	    //select only real prompt
	    if (Option == 1) {
	      if (!(abs(matchedID) == 11)) continue;
	    }	   


	    //Find Closest Parton
	    float minDRToParton = 9999;
	    for(int j = 0; j < nGenParticle; j++){
	      
	      //only look for outgoing partons
	      if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		      && gParticleStatus[j] == 23)
		   ) continue;
	      
	      double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], eleEta[i], elePhi[i]);
	      if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	    }
	    eleTree->fDRToClosestParton = minDRToParton;


	    //***********************
	    //Fill Electron
	    //***********************
	    NElectronsFilled++;
	    eleTree->tree_->Fill();	   
        }
      }

      //********************************************
      //Tree entries based on gen-level objects
      //********************************************
      else if (Option >= 10 && Option < 20) {

	for(int i = 0; i < nGenParticle; i++){
	    
	  //select prompt electrons
	  if (Option == 11) {
	    if (abs(gParticleId[i]) != 11) continue;	      
	    if (!(abs(gParticleMotherId[i]) == 23 || 
		  abs(gParticleMotherId[i]) == 24)
		) continue;
	  }

	  if(gParticlePt[i] < 5) continue;

	  //***********************
	  //Fill Electron Variables
	  //***********************
   	  eleTree->fWeight = 1;
	  eleTree->fRunNumber = runNum;
	  eleTree->fLumiSectionNumber = lumiNum;
	  eleTree->fEventNumber = eventNum;
	  eleTree->fEleEventNumberParity = (eventNum % 2 == 0);
	  eleTree->fEleGenPt = gParticlePt[i];
	  eleTree->fEleGenEta = gParticleEta[i];
	  eleTree->fEleGenPhi = gParticlePhi[i];
	  eleTree->fRho = 0; 
	  eleTree->fNVertices = nPV; 
	  eleTree->fPdgId = gParticleId[i];


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
	  eleTree->fDRToClosestParton = minDRToParton;


	  //Find Associated Reco Electron
	  int matchedIndex = -1;
	  float minDR = 9999;
	  for(int j = 0; j < nElectrons; j++){
	    if ( deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]) < 0.1
		 && deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]) < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]);
	    }	    
	  }	 

	  if (matchedIndex >= 0) {
	    eleTree->fCharge = eleCharge[matchedIndex] ;
	    eleTree->fElePt = elePt[matchedIndex]; 
	    eleTree->fEleEta = eleEta[matchedIndex]; 
	    eleTree->fElePhi = elePhi[matchedIndex]; 
	    eleTree->fEleSCEta = eleEta_SC[matchedIndex]; 
	    eleTree->fEleTriggerBit = 0;
	    eleTree->fEleD0 = ele_d0[matchedIndex]; 
	    eleTree->fEleDZ = ele_dZ[matchedIndex]; 
	    eleTree->fElePassConversion = ele_PassConvVeto[matchedIndex];
	    eleTree->fEleNMissHits =ele_MissHits[matchedIndex];
	    eleTree->fEleOneOverEMinusOneOverP = ele_OneOverEminusOneOverP[matchedIndex];
	    eleTree->fEleDEtaIn = ele_dEta[matchedIndex];					   
	    eleTree->fEleDPhiIn = ele_dPhi[matchedIndex];
	    eleTree->fEleSigmaIEtaIEta = eleFull5x5SigmaIetaIeta[matchedIndex];
	    eleTree->fEleR9 = eleR9[matchedIndex];
	    eleTree->fEleHoverE = ele_HoverE[matchedIndex];
	    eleTree->fElePFIso04 = (ele_chargedIso[matchedIndex] + fmax(0.0,  ele_photonIso[matchedIndex] + ele_neutralHadIso[matchedIndex] - 0.5*ele_pileupIso[matchedIndex])) / elePt[matchedIndex];
	    eleTree->fIDMVATrig = ele_IDMVATrig[matchedIndex];
	    eleTree->fIDMVANonTrig = ele_IDMVANonTrig[matchedIndex];
	    eleTree->fPassVetoSelection = isVetoElectron(matchedIndex);
	    eleTree->fPassLooseSelection = isRunOneLooseElectron(matchedIndex);
	    eleTree->fPassTightSelection = isRunOneTightElectron(matchedIndex);
	    eleTree->fPassMVANonTrigVetoSelection = isMVANonTrigVetoElectron(matchedIndex);
	    eleTree->fPtRel = ele_ptrel[matchedIndex];
	    eleTree->fMiniIso = ele_miniiso[matchedIndex];
	  } else {
	    eleTree->fCharge = 0;
	    eleTree->fElePt = 0;
	    eleTree->fEleEta = 0;
	    eleTree->fElePhi = 0;
	    eleTree->fEleSCEta = 0;
	    eleTree->fEleTriggerBit = 0;
	    eleTree->fRho = 0; 
	    eleTree->fEleD0 = 0;
	    eleTree->fEleDZ = 0;
	    eleTree->fElePassConversion = false;
	    eleTree->fEleNMissHits = 0;
	    eleTree->fEleOneOverEMinusOneOverP = 0;
	    eleTree->fEleDEtaIn = 0;
	    eleTree->fEleDPhiIn = 0;
	    eleTree->fEleSigmaIEtaIEta = 0;
	    eleTree->fEleR9 = 0;
	    eleTree->fEleHoverE = 0;
	    eleTree->fElePFIso04 = 0;
	    eleTree->fIDMVATrig = 0;
	    eleTree->fIDMVANonTrig = 0;
	    eleTree->fPassVetoSelection = false;
	    eleTree->fPassLooseSelection = false;
	    eleTree->fPassTightSelection = false;
	    eleTree->fPassMVANonTrigVetoSelection = false;
	    eleTree->fPtRel = 0;
	    eleTree->fMiniIso = 0;
	  }
	  
	  //***********************
	  //Fill Electron
	  //***********************
	  NElectronsFilled++;
	  eleTree->tree_->Fill();
	}

      }


    }//end of event loop

    cout << "Filled Total of " << NElectronsFilled << " Electrons\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}



