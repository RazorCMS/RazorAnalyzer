#include "PhotonNtupler.h"
#include "PhotonTree.h"

using namespace std;

void PhotonNtupler::Analyze(bool isData, int Option, string outputFilename, string label)
{
    cout << "Initializing..." << endl;
    string outfilename = outputFilename;
    if (outfilename == "") outfilename = "PhotonNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    PhotonTree *phoTree = new PhotonTree;
    phoTree->CreateTree(PhotonTree::kPhotonTreeLight);
    phoTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";
    bool use25nsSelection = true;
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
	  phoTree->fpvZ_New = pvZ_New;
	  phoTree->fpvdT_New = pvdT_New;
	  phoTree->fpvdZ_New = pvdZ_New;
	  phoTree->fpvZ_Gen = genVertexZ;
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
	  phoTree->fPhoChargedHadronIso_NewPV_NoTiming = pho_sumChargedHadronPt_NewPV_NoTiming[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing50_TrkVtx = pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing80_TrkVtx = pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing100_TrkVtx = pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing120_TrkVtx = pho_sumChargedHadronPt_NewPV_Timing120_TrkVtx[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing50_TrkPho = pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing80_TrkPho = pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing100_TrkPho = pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[i];
	  phoTree->fPhoChargedHadronIso_NewPV_Timing120_TrkPho = pho_sumChargedHadronPt_NewPV_Timing120_TrkPho[i];
	
	  phoTree->fPhoNeutralHadronIso = pho_sumNeutralHadronEt[i];
	  phoTree->fPhoPhotonIso = pho_sumPhotonEt[i];
	  phoTree->fPhoIsConversion = pho_isConversion[i];
	  phoTree->fPhoPassEleVeto = pho_passEleVeto[i];
	  phoTree->fPhoRegressionE = pho_RegressionE[i];
	  phoTree->fPhoRegressionSigmaE = pho_RegressionEUncertainty[i];
	  phoTree->fPhoIDMVA = pho_IDMVA[i];
	  phoTree->fPhoSCEta = pho_superClusterEta[i];
	  phoTree->fPhoHasPixelSeed = pho_hasPixelSeed[i];


	  phoTree->fPhoIsLoose = isLoosePhoton(i,use25nsSelection);
	  phoTree->fPhoIsMedium = isMediumPhoton(i,use25nsSelection);
	  phoTree->fPhoIsTight = isTightPhoton(i,use25nsSelection);
	  phoTree->fPhoPassLooseID = photonPassLooseID(i,use25nsSelection);
	  phoTree->fPhoPassMediumID = photonPassMediumID(i,use25nsSelection);
	  phoTree->fPhoPassTightID = photonPassTightID(i,use25nsSelection);
	  phoTree->fPhoPassLooseIso = photonPassLooseIso(i,use25nsSelection);
	  phoTree->fPhoPassMediumIso = photonPassMediumIso(i,use25nsSelection);
	  phoTree->fPhoPassTightIso = photonPassTightIso(i,use25nsSelection);
	  phoTree->fPV_ndof = PV_ndof; 

	  bool foundMatch = false;
	  int phoMatchingGenPhotonIndex = 0;
	  double deltaEOverEBest = 999;
	  double deltaRBest = 999;
	  float minDRToParton = 9999;
	  //int closestPartonIndex = -1;

	bool isFake = false;
 	bool closeToPrompt = false;

	for(int g = 0; g < nGenParticle && (!closeToPrompt); g++){

	if(deltaR(phoEta[i], phoPhi[i], gParticleEta[g], gParticlePhi[g]) > 0.3 ) continue;

	bool isGenElectron = (gParticle_isPromptFinalState[g] || gParticleStatus[g]==1) && (abs(gParticleId[g]) == 11);
	bool isGenMuon = (gParticle_isPromptFinalState[g] || gParticleStatus[g]==1 )&& (abs(gParticleId[g]) == 13);

	//bool isGenElectron = gParticle_isPromptFinalState[g] && (abs(gParticleId[g]) == 11);
	//bool isGenMuon = gParticle_isPromptFinalState[g] && (abs(gParticleId[g]) == 13);
	bool isGenTau = gParticle_isPromptDecayed[g] && (abs(gParticleId[g]) == 15);
 	bool isPromptPhoton = gParticle_isPromptFinalState[g] && (gParticleId[g] == 22); //(gParticleId[g] == 22) && (abs(gParticleMotherId[g])==25 || (abs(gParticleMotherId[g])>=1 && abs(gParticleMotherId[g])<=6) || (abs(gParticleMotherId[g])>=11 && abs(gParticleMotherId[g])<=16));

        if(isGenElectron || isGenMuon || isGenTau || isPromptPhoton) closeToPrompt = true;
	}
	
	if(!closeToPrompt)
	{	
	for(int jj = 0; jj < nGenJets && (!isFake); jj++)
	{
	  
	   if(genJetPt[jj]<6.0) continue;
 
	   float deltaR_tmp = deltaR(phoEta[i], phoPhi[i], genJetEta[jj], genJetPhi[jj]);
	  
	   if(deltaR_tmp>0.3) continue;

		bool GenJetCloseToPrompt = false;
		for(int g = 0; g < nGenParticle && (!GenJetCloseToPrompt); g++){
	
		if(deltaR(genJetEta[jj], genJetPhi[jj], gParticleEta[g], gParticlePhi[g]) > 0.3 ) continue;
		
		bool isGenElectron = (gParticle_isPromptFinalState[g] || gParticleStatus[g]==1) && (abs(gParticleId[g]) == 11);
		bool isGenMuon = (gParticle_isPromptFinalState[g] || gParticleStatus[g]==1 )&& (abs(gParticleId[g]) == 13);
		bool isGenTau = gParticle_isPromptDecayed[g] && (abs(gParticleId[g]) == 15);
		bool isPromptPhoton = gParticle_isPromptFinalState[g] && (gParticleId[g] == 22); //(gParticleId[g] == 22) && (abs(gParticleMotherId[g])==25 || (abs(gParticleMotherId[g])>=1 && abs(gParticleMotherId[g])<=6) || (abs(gParticleMotherId[g])>=11 && abs(gParticleMotherId[g])<=16));

        	if(isGenElectron || isGenMuon || isGenTau || isPromptPhoton) GenJetCloseToPrompt = true;
	
		}
		
	   if(GenJetCloseToPrompt) continue;
	   isFake = true;
	}
	}

	if(isFake)
	{

	float min_deltaR_toGen = 9999.;
	phoTree->fNearstHighPtGenPdgId = 999 ;
	phoTree->fNearstHighPtGenMotherId = 999 ;
	for(int g = 0; g < nGenParticle; g++){
	
	if(gParticleStatus[g] != 1) continue;
	if(gParticlePt[g]<10.0) continue;
	float deltaR_tmp = deltaR(phoEta[i], phoPhi[i], gParticleEta[g], gParticlePhi[g]);
	if(deltaR_tmp > 0.3) continue;
        
	if(deltaR_tmp < min_deltaR_toGen ) 
	{
		min_deltaR_toGen = deltaR_tmp; 
		phoTree->fNearstHighPtGenPdgId = gParticleId[g];
		phoTree->fNearstHighPtGenMotherId = gParticleMotherId[g];
	
	}
	}
	
	}

	bool doFillPhoton = false;

	if(Option !=0)
	{
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

 	   float deltaR_tmp = deltaR(phoEta[i], phoPhi[i], gParticleEta[g], gParticlePhi[g]);
	   foundMatch = true;
	   if(deltaR_tmp < deltaRBest)
	   {
		deltaRBest = deltaR_tmp;
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

	  bool isPromptPhoton = false;
 	  bool isHiggsPhoton = false;
	  bool isPi0EtaPhoton = false;

	  if (abs(phoTree->fMotherPdgId) == 25 || ( abs(phoTree->fMotherPdgId) >= 1 && abs(phoTree->fMotherPdgId) <= 6) || 
	      ( abs(phoTree->fMotherPdgId) >= 11 && abs(phoTree->fMotherPdgId) <= 16)) isPromptPhoton = true;
	  if(abs(phoTree->fMotherPdgId) == 25) isHiggsPhoton =  true;
	  if(abs(phoTree->fMotherPdgId) == 111 || abs(phoTree->fMotherPdgId) == 221) isPi0EtaPhoton =  true;
	    	  
	  	  if (Option == 1) {// real photons
	    //if (isPromptPhoton) doFillPhoton = true;
	    if (isHiggsPhoton) doFillPhoton = true;
	  }
	  }
	  
          if (Option == 0) {//fake photons
	    //if (!isPromptPhoton) doFillPhoton = true;
	    if (isFake) doFillPhoton = true;
	  }

	  if (doFillPhoton && !isFakePV && abs(PV_ndof)>0.001 ) {

	    // cout << "DEBUG\n";
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
	       phoTree->fPhoIsLoose = isLoosePhoton(matchedIndex,use25nsSelection);
	       phoTree->fPhoIsMedium = isMediumPhoton(matchedIndex,use25nsSelection);
	       phoTree->fPhoIsTight = isTightPhoton(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassLooseID = photonPassLooseID(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassMediumID = photonPassMediumID(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassTightID = photonPassTightID(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassLooseIso = photonPassLooseIso(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassMediumIso = photonPassMediumIso(matchedIndex,use25nsSelection);
	       phoTree->fPhoPassTightIso = photonPassTightIso(matchedIndex,use25nsSelection);
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
