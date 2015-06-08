#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "PhotonTree.h"

using namespace std;

void RazorAnalyzer::PhotonNtupler(string outputFilename , int Option)
{
    cout << "Initializing..." << endl;
    string outfilename = outputFilename;
    if (outfilename == "") outfilename = "PhotonNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    PhotonTree *phoTree = new PhotonTree;
    phoTree->CreateTree(PhotonTree::kPhotonTreeLight);
    phoTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";

    UInt_t NPhotonsFilled = 0;


    //begin loop
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

      bool printDecay = false;

      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {

        //fill photon information
        for(int i = 0; i < nPhotons; i++){
	  //try to match photon with a gen photon

	  //cout << "Photon " << i << " : " << phoPt[i] << " " << phoEta[i] << " " << phoPhi[i]<< " : " << phoSigmaIetaIeta[i] << " " << pho_HoverE[i] << " \n";
 
	  //***********************
	  //Fill Photon Variables
	  //***********************
	  phoTree->fWeight = 1;
	  phoTree->fRunNumber = runNum;
	  phoTree->fLumiSectionNumber = lumiNum;
	  phoTree->fEventNumber = eventNum;
	  phoTree->fPhoEventNumberParity = (eventNum % 2 == 0);

	  // cout << "start ? " << bool(!phoTree->fPhoIsLoose && phoSigmaIetaIeta[i] < 0.01 && pho_HoverE[i] < 0.028) << " "
	  //      << bool(phoSigmaIetaIeta[i] < 0.01) << " " << bool(pho_HoverE[i] < 0.028) << " \n";

	  phoTree->fRho = fixedGridRhoFastjetAll;
	  phoTree->fNVertices = nPV;
	  phoTree->fPhoPt = phoPt[i];
	  phoTree->fPhoEta = phoEta[i];
	  phoTree->fPhoPhi = phoPhi[i];
	  phoTree->fPhoSigmaIetaIeta = phoSigmaIetaIeta[i];
	  phoTree->fPhoFull5x5SigmaIetaIeta = phoFull5x5SigmaIetaIeta[i];
	  phoTree->fPhoR9 = phoR9[i];
	  phoTree->fPhoHOverE = pho_HoverE[i];
	  phoTree->fPhoChargedHadronIso = pho_sumChargedHadronPt[i];
	  phoTree->fPhoNeutralHadronIso = pho_sumNeutralHadronEt[i];
	  phoTree->fPhoPhotonIso = pho_sumPhotonEt[i];
	  phoTree->fPhoIsConversion = pho_isConversion[i];
	  phoTree->fPhoPassEleVeto = pho_passEleVeto[i];
	  phoTree->fPhoRegressionE = pho_RegressionE[i];
	  phoTree->fPhoRegressionSigmaE = pho_RegressionEUncertainty[i];
	  phoTree->fPhoIDMVA = pho_IDMVA[i];
	  phoTree->fPhoSCEta = pho_superClusterEta[i];
	  phoTree->fPhoHasPixelSeed = pho_hasPixelSeed[i];


	  phoTree->fPhoIsLoose = isGoodPhotonRun2( i , true, WP::Loose, false );
	  phoTree->fPhoIsMedium = isGoodPhotonRun2( i , true, WP::Medium, false );
	  phoTree->fPhoIsTight = isGoodPhotonRun2( i , true, WP::Tight, false );
	  phoTree->fPhoPassLooseID = isGoodPhotonRun2( i , false, WP::Loose, false );
	  phoTree->fPhoPassMediumID = isGoodPhotonRun2( i , false, WP::Medium, false );
	  phoTree->fPhoPassTightID = isGoodPhotonRun2( i , false, WP::Tight, false );
	  phoTree->fPhoPassLooseIso = photonPassIsoRun2( i, WP::Loose, false);
	  phoTree->fPhoPassMediumIso = photonPassIsoRun2( i, WP::Medium, false);
	  phoTree->fPhoPassTightIso = photonPassIsoRun2( i, WP::Tight, false);

	  bool foundMatch = false;
	  int phoMatchingGenPhotonIndex = -1;
	  double deltaEOverEBest = 999;
	  float minDRToParton = 9999;
	  //int closestPartonIndex = -1;
	  for(int g = 0; g < nGenParticle; g++){

	    //find closest parton to photon
	    if  ( ((abs(gParticleId[g]) >= 1 && abs(gParticleId[g]) <= 5) || abs(gParticleId[g]) == 21)
		  && gParticleStatus[g] == 23
		  ) {
	      double tmpDR = deltaR( gParticleEta[g], gParticlePhi[g], phoEta[i], phoPhi[i]);
	      if ( tmpDR < minDRToParton ) {
		//closestPartonIndex = g;
		minDRToParton = tmpDR;
	      }
	    }


	    //Find Matching Gen Photon
	    if(gParticleStatus[g] != 1) continue;
	    if(gParticleId[g] != 22) continue;
	    if(gParticleE[g] < 1.) continue;

	    if(deltaR(phoEta[i], phoPhi[i], gParticleEta[g], gParticlePhi[g]) > 0.2) continue;
	    float deltaEOverE = fabs(pho_RegressionE[i] - gParticleE[g])/gParticleE[g];
	    if(deltaEOverE > 1.) continue;
	    
	    foundMatch = true;	    
	    if(deltaEOverE < deltaEOverEBest){
	      deltaEOverEBest = deltaEOverE;
	      phoMatchingGenPhotonIndex = g;
	    }
	   
	  }

	  
	  phoTree->fDRToClosestParton = minDRToParton;
	  if (foundMatch) {
	    phoTree->fMotherPdgId = gParticleMotherId[phoMatchingGenPhotonIndex];
	    phoTree->fPhoGenE = gParticleE[phoMatchingGenPhotonIndex];
	    phoTree->fPhoGenPt = gParticlePt[phoMatchingGenPhotonIndex];
	    phoTree->fPhoGenEta = gParticleEta[phoMatchingGenPhotonIndex];
	    phoTree->fPhoGenPhi = gParticlePhi[phoMatchingGenPhotonIndex];
	  } else {
	    phoTree->fMotherPdgId = 0;
	    phoTree->fPhoGenE = 0;
	    phoTree->fPhoGenPt = 0;
	    phoTree->fPhoGenEta = 0;
	    phoTree->fPhoGenPhi = 0;
	  }
        

	  if (printDecay) {
	     for(int g = 0; g < nGenParticle; g++){
	       cout << g << " : " << gParticleId[g] << " "  << gParticleStatus[g] << " : " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " " << gParticleMotherId[g] << "\n";
	     }
	  }

	  bool doFillPhoton = false;
	  bool isPromptPhoton = false;
	  if (abs(phoTree->fMotherPdgId) == 25 || ( abs(phoTree->fMotherPdgId) >= 1 && abs(phoTree->fMotherPdgId) <= 6) || 
	      ( abs(phoTree->fMotherPdgId) >= 11 && abs(phoTree->fMotherPdgId) <= 16)) isPromptPhoton = true;
	  
	  if (Option == 0) {
	    if (!isPromptPhoton) doFillPhoton = true;
	  }
	  if (Option == 1) {
	    if (isPromptPhoton) doFillPhoton = true;
	  }
	  if (doFillPhoton) {

	    // cout << "DEBUG\n";
	    // isGoodPhotonRun2( i , true, WP::Loose, true );

	    // cout << "weird ? " << bool(!phoTree->fPhoIsLoose && phoSigmaIetaIeta[i] < 0.01 && pho_HoverE[i] < 0.028) << " "
	    // 	 << bool(phoSigmaIetaIeta[i] < 0.01) << " " << bool(pho_HoverE[i] < 0.028) << " " << phoTree->fPhoIsLoose << "\n";
	    // if (!phoTree->fPhoIsLoose && phoSigmaIetaIeta[i] < 0.01 && pho_HoverE[i] < 0.028) {
	    //   cout << runNum << " " << eventNum << "\n";
	      
	    //   cout << "check: " << phoTree->fPhoSigmaIetaIeta << " " << phoTree->fPhoHOverE << "\n";
	    // }
	    	    
	    NPhotonsFilled++;	     	    
	    phoTree->tree_->Fill();
	  }
	}
      }
      
      //********************************************
      //Tree entries based on gen-level objects
      //********************************************
      else if (Option >= 10 && Option < 20) {

        //select gen photons
        for(int g = 0; g < nGenParticle; g++){
           if(gParticleStatus[g] != 1) continue;
           if(gParticleId[g] != 22) continue;
           if(gParticlePt[g] < 5.) continue;
           
	   if( abs(gParticleMotherId[g]) == 25 || abs(gParticleMotherId[g]) == 24 || ( abs(gParticleMotherId[g]) >= 1 && abs(gParticleMotherId[g]) <= 5) 
	       || ( abs(gParticleMotherId[g]) >= 11 && abs(gParticleMotherId[g]) <= 16)
	       ) {

	     phoTree->fWeight = 1;
	     phoTree->fRunNumber = runNum;
	     phoTree->fLumiSectionNumber = lumiNum;
	     phoTree->fEventNumber = eventNum;
	     phoTree->fPhoEventNumberParity = (eventNum % 2 == 0);
	     phoTree->fRho = fixedGridRhoFastjetAll;
	     phoTree->fNVertices = nPV;
	     
	     phoTree->fMotherPdgId = gParticleMotherId[g];
	     phoTree->fPhoGenE = gParticleE[g];
	     phoTree->fPhoGenPt = gParticlePt[g];
	     phoTree->fPhoGenEta = gParticleEta[g];
	     phoTree->fPhoGenPhi = gParticlePhi[g];

	     //Find Closest Parton
	     float minDRToParton = 9999;
	     for(int j = 0; j < nGenParticle; j++){	      
	       //only look for outgoing partons
	       if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		       && gParticleStatus[j] == 23)
		    ) continue;	       
	       double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], gParticleEta[g], gParticlePhi[g]);
	       if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	     }
	     phoTree->fDRToClosestParton = minDRToParton;
	     
	     //Find Associated Reco Photon
	     int matchedIndex = -1;
	     double deltaEOverEBest = 999;

	     for(int j = 0; j < nPhotons; j++){
	       float deltaEOverE = fabs(pho_RegressionE[j] - gParticleE[g])/gParticleE[g];
	       if(deltaEOverE > 1.) continue;
	       if ( deltaR( phoEta[j], phoPhi[j], gParticleEta[g], gParticlePhi[g]) < 0.1
		    && deltaEOverE < deltaEOverEBest
		    ) {		
		 matchedIndex = j;
		 deltaEOverEBest = deltaEOverE;
	       }
	     } 
	     
	     if (matchedIndex > -1) {
	       phoTree->fPhoPt = phoPt[matchedIndex];
	       phoTree->fPhoEta = phoEta[matchedIndex];
	       phoTree->fPhoPhi = phoPhi[matchedIndex];
	       phoTree->fPhoSigmaIetaIeta = phoSigmaIetaIeta[matchedIndex];
	       phoTree->fPhoFull5x5SigmaIetaIeta = phoFull5x5SigmaIetaIeta[matchedIndex];
	       phoTree->fPhoR9 = phoR9[matchedIndex];
	       phoTree->fPhoHOverE = pho_HoverE[matchedIndex];
	       phoTree->fPhoChargedHadronIso = pho_sumChargedHadronPt[matchedIndex];
	       phoTree->fPhoNeutralHadronIso = pho_sumNeutralHadronEt[matchedIndex];
	       phoTree->fPhoPhotonIso = pho_sumPhotonEt[matchedIndex];
	       phoTree->fPhoIsConversion = pho_isConversion[matchedIndex];
	       phoTree->fPhoPassEleVeto = pho_passEleVeto[matchedIndex];
	       phoTree->fPhoRegressionE = pho_RegressionE[matchedIndex];
	       phoTree->fPhoRegressionSigmaE = pho_RegressionEUncertainty[matchedIndex];
	       phoTree->fPhoIDMVA = pho_IDMVA[matchedIndex];
	       phoTree->fPhoSCEta = pho_superClusterEta[matchedIndex];
	       phoTree->fPhoHasPixelSeed = pho_hasPixelSeed[matchedIndex];
	       phoTree->fPhoIsLoose = isGoodPhotonRun2( matchedIndex , true, WP::Loose, false );
	       phoTree->fPhoIsMedium = isGoodPhotonRun2( matchedIndex , true, WP::Medium, false );
	       phoTree->fPhoIsTight = isGoodPhotonRun2( matchedIndex , true, WP::Tight, false );
	       phoTree->fPhoPassLooseID = isGoodPhotonRun2( matchedIndex , false, WP::Loose, false );
	       phoTree->fPhoPassMediumID = isGoodPhotonRun2( matchedIndex , false, WP::Medium, false );
	       phoTree->fPhoPassTightID = isGoodPhotonRun2( matchedIndex , false, WP::Tight, false );
	       phoTree->fPhoPassLooseIso = photonPassIsoRun2( matchedIndex, WP::Loose, false);
	       phoTree->fPhoPassMediumIso = photonPassIsoRun2( matchedIndex, WP::Medium, false);
	       phoTree->fPhoPassTightIso = photonPassIsoRun2( matchedIndex, WP::Tight, false);
	     }

	     bool doFillPhoton = false;
	     bool isPromptPhoton = false;
	     if (abs(phoTree->fMotherPdgId) == 25 || ( abs(phoTree->fMotherPdgId) >= 1 && abs(phoTree->fMotherPdgId) <= 6) || 
		 ( abs(phoTree->fMotherPdgId) >= 11 && abs(phoTree->fMotherPdgId) <= 16)) isPromptPhoton = true;
	     
	     if (isPromptPhoton) doFillPhoton = true;
	     
	     if (doFillPhoton) {
	       NPhotonsFilled++;	     	    
	       phoTree->tree_->Fill();
	     }
	   }        
        }
      } 

    } // loop over events 
    
    cout << "Filled Total of " << NPhotonsFilled << " Photons\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
    
}
