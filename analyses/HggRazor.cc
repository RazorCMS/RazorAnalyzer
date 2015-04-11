#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

enum HggRazorBox {
    HighPt,
    Hbb,
    Zbb,
    HighRes,
    LowRes
};

void RazorAnalyzer::HggRazor(string outFileName, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    if (outFileName.empty()){
        cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << endl;
        outFileName = "HggRazor.root";
    }
    TFile outFile(outFileName.c_str(), "RECREATE");

    //one tree to hold all events
    TTree *razorTree = new TTree("HggRazor", "Info on selected razor inclusive events");

    //separate trees for individual boxes
    map<string, TTree*> razorBoxes;
    vector<string> boxNames;
    boxNames.push_back("HighPt");
    boxNames.push_back("Hbb");
    boxNames.push_back("Zbb");
    boxNames.push_back("HighRes");
    boxNames.push_back("LowRes");
    for(size_t i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nLooseBTaggedJets, nMediumBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
    float theMR;
    float theRsq;
    int nSelectedPhotons;
    float mGammaGamma, pTGammaGamma;
    float mbbZ, mbbH;
    float sigmaEOverE1, sigmaEOverE2;
    HggRazorBox box;

    //set branches on big tree
    if(combineTrees){
        razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
        razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
        razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
        razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
        razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
        razorTree->Branch("mbbZ", &mbbZ, "mbbZ/F");
        razorTree->Branch("mbbH", &mbbH, "mbbH/F");
        razorTree->Branch("sigmaEOverE1", &sigmaEOverE1, "sigmaEOverE1/F");
        razorTree->Branch("sigmaEOverE2", &sigmaEOverE2, "sigmaEOverE2/F");
        razorTree->Branch("box", &box, "box/I");
    }
    //set branches on all trees
    else{ 
        for(auto& box : razorBoxes){
            box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
            box.second->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
            box.second->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
            box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
            box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
            box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
            box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
            box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
            box.second->Branch("MR", &theMR, "MR/F");
            box.second->Branch("Rsq", &theRsq, "Rsq/F");
            box.second->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
            box.second->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
            box.second->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
            box.second->Branch("mbbZ", &mbbZ, "mbbZ/F");
            box.second->Branch("mbbH", &mbbH, "mbbH/F");
            box.second->Branch("sigmaEOverE1", &sigmaEOverE1, "sigmaEOverE1/F");
            box.second->Branch("sigmaEOverE2", &sigmaEOverE2, "sigmaEOverE2/F");
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
        nLooseBTaggedJets = 0;
        nMediumBTaggedJets = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nTightTaus = 0;
        theMR = -1;
        theRsq = -1;
        nSelectedPhotons = 0;
        mGammaGamma = -1;
        pTGammaGamma = -1;
        mbbZ = 0;
        mbbH = 0;
        sigmaEOverE1 = -1;
        sigmaEOverE2 = -1;
        if(combineTrees) box = LowRes;

        //TODO: triggers!
        bool passedDiphotonTrigger = true;
        if(!passedDiphotonTrigger) continue;

        //muon selection
        for(int i = 0; i < nMuons; i++){
            if(!isLooseMuon(i)) continue;  
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;

            nLooseMuons++;

            if(isTightMuon(i)){ 
                nTightMuons++;
            }
        }
        //electron selection
        for(int i = 0; i < nElectrons; i++){
            if(!isLooseElectron(i)) continue; 
            if(elePt[i] < 10) continue;
            if(abs(eleEta[i]) > 2.5) continue;

            nLooseElectrons++;

            if(isTightElectron(i)){ 
                nTightElectrons++;
            }
        }
        //tau selection
        for(int i = 0; i < nTaus; i++){
            if(!isTightTau(i)) continue; 

            nTightTaus++;
        }
	
        //photon selection
        vector<TLorentzVector> GoodPhotons;
        vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
        vector<bool> GoodPhotonPassesIso; //store whether each photon is isolated
        int nPhotonsAbove40GeV = 0;
        for(int i = 0; i < nPhotons; i++){
	  //ID cuts -- apply isolation after candidate pair selection
	  //if(!isMediumPhotonNoIsoCuts(i)){
	  if ( !isGoodPhotonRun1( i , false ) )
	    {
	      std::cout << "[INFO]: Failed photon ID" << std::endl;
	      continue;
	    }

	  double pho_pt = pho_RegressionE[i]/cosh( pho_superClusterEta[i] );//regressed PT	  
	  //if(phoPt[i] < 25){
	  if ( pho_pt < 25.0 )
	    {
	      continue;
	    }
	  
            if(fabs(pho_superClusterEta[i]) > 2.5){
                //allow photons in the endcap here, but if one of the two leading photons is in the endcap, reject the event
                continue; 
            }

            //photon passes
            if( pho_pt > 40.0 ) nPhotonsAbove40GeV++;
            TLorentzVector thisPhoton = makeTLorentzVector( pho_pt, pho_superClusterEta[i], phoPhi[i], pho_RegressionE[i] );
            GoodPhotons.push_back( thisPhoton );
            GoodPhotonSigmaE.push_back( pho_RegressionEUncertainty[i] );
            //GoodPhotonPassesIso.push_back(photonPassesMediumIsoCuts(i));
	    GoodPhotonPassesIso.push_back( isGoodPhotonRun1( i , true ) );
            nSelectedPhotons++;
        }
        //if there is no photon with pT above 40 GeV, reject the event
        if( nPhotonsAbove40GeV == 0 ){
	  continue;
        }
	
	//find the "best" photon pair
        TLorentzVector HiggsCandidate(0,0,0,0);
        int goodPhoIndex1 = -1;
        int goodPhoIndex2 = -1;
        double bestSumPt = 0;
        for(size_t i = 0; i < GoodPhotons.size(); i++){
	  for(size_t j = i+1; j < GoodPhotons.size(); j++){//I like this logic better, I find it easier to understand
                TLorentzVector pho1 = GoodPhotons[i];
                TLorentzVector pho2 = GoodPhotons[j];
                
                //need one photon in the pair to have pt > 40 GeV
                if( pho1.Pt() < 40.0 && pho2.Pt() < 40.0 ){
		  continue;
                }
                //need diphoton mass between > 100 GeV as in AN (April 1st)
                double diphotonMass = (pho1 + pho2).M();
                //if(diphotonMass < 100 || diphotonMass > 180){
		if( diphotonMass < 100.0 ){
		  continue;
                }
                
                //if the sum of the photon pT's is larger than that of the current Higgs candidate, make this the Higgs candidate
                if( pho1.Pt() + pho2.Pt() > bestSumPt ){
		  bestSumPt = pho1.Pt() + pho2.Pt();
		  HiggsCandidate = pho1 + pho2;
		  goodPhoIndex1 = i;
		  goodPhoIndex2 = j;  
                }
            }
        }   
        //if the best candidate pair has pT < 20 GeV, reject the event
        if( HiggsCandidate.Pt() < 20.0 ){
            continue;
        }

	//if the best candidate pair has a photon in the endcap, reject the event
	//Reject gap photons
        if( fabs(GoodPhotons[goodPhoIndex1].Eta()) > 1.44 || fabs(GoodPhotons[goodPhoIndex2].Eta()) > 1.44 ){
            continue;
        }
	//if the best candidate pair has a non-isolated photon, reject the event
        if( !GoodPhotonPassesIso[goodPhoIndex1] || !GoodPhotonPassesIso[goodPhoIndex2] ){
            continue;
        }
	//record higgs candidate info
        mGammaGamma = HiggsCandidate.M();
        pTGammaGamma = HiggsCandidate.Pt();
        sigmaEOverE1 = GoodPhotonSigmaE[goodPhoIndex1]/GoodPhotons[goodPhoIndex1].E();
        sigmaEOverE2 = GoodPhotonSigmaE[goodPhoIndex2]/GoodPhotons[goodPhoIndex2].E();
		
	//Jets
	vector<TLorentzVector> GoodJets;
        vector<pair<TLorentzVector, bool> > GoodCSVLJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not
        for(int i = 0; i < nJets; i++){
	  //if(jetPt[i] < 40) continue;
	  if( jetPt[i] < 30.0 ) continue;//According to the April 1st 2015 AN
	  if( fabs(jetEta[i]) >= 3.0 ) continue;

            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            //exclude selected photons from the jet collection
            double deltaRJetPhoton = min(thisJet.DeltaR(GoodPhotons[goodPhoIndex1]), thisJet.DeltaR(GoodPhotons[goodPhoIndex2]));
            //if(deltaRJetPhoton < 0.4) continue;
	    if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN

            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVL(i)){
                nLooseBTaggedJets++;
                if(isCSVM(i)){ 
                    nMediumBTaggedJets++;
                    GoodCSVLJets.push_back(make_pair(thisJet, true));
                }
                else{
                    GoodCSVLJets.push_back(make_pair(thisJet, false));
                }
            }
        }
	
        //if there are no good jets, reject the event
        if( nSelectedJets == 0 ){
            continue;
        }

        //Compute the razor variables using the selected jets and the diphoton system
        vector<TLorentzVector> JetsPlusHiggsCandidate;
        for(auto& jet : GoodJets) JetsPlusHiggsCandidate.push_back(jet);
        JetsPlusHiggsCandidate.push_back(HiggsCandidate);

        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

        vector<TLorentzVector> hemispheres = getHemispheres(JetsPlusHiggsCandidate);
        theMR = computeMR(hemispheres[0], hemispheres[1]); 
        theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
        //if MR < 200, reject the event
        if ( theMR < 150.0 ){
            continue;
        }

        //if there are two loose b-tags and one medium b-tag, look for b-bbar resonances
        if( nLooseBTaggedJets > 1 && nMediumBTaggedJets > 0 )
	  {
	    for(int i = 0; i < nLooseBTaggedJets; i++){
	      for(int j = i+1; j < nLooseBTaggedJets; j++){
		//if neither of the b-jets passes CSVM, continue
		if( !GoodCSVLJets[i].second && !GoodCSVLJets[j].second ) continue;
		double mbb = (GoodCSVLJets[i].first + GoodCSVLJets[j].first).M();
		//if mbb is closer to the higgs mass than mbbH, make mbbH = mbb
		if( fabs(mbbH - 125.0) > fabs(mbb - 125.0) ) mbbH = mbb;
		//same for mbbZ
		if( fabs(mbbZ - 91.2) > fabs(mbb - 91.2) ) mbbZ = mbb;
	      }//end second jet loop
	    }//end first jet loop
	  }
	
	if ( sigmaEOverE1 < 0.015 && sigmaEOverE2 < 0.015 ) std::cout << "[INFO]: SigmaEoverE1: " << sigmaEOverE1 << " SigmaEoverE2: " << sigmaEOverE2 << std::endl;
        //HighPt Box
        if ( pTGammaGamma > 110.0 )
	  {
	    if(combineTrees){
	      box = HighPt;
	      razorTree->Fill();
	    }
	    else razorBoxes["HighPt"]->Fill();
	  }
        //Hbb Box
	else if ( mbbH > 110.0 && mbbH < 140.0 )
	  {
	    if(combineTrees){
	      box = Hbb;
	      razorTree->Fill();
	    }
	    else razorBoxes["Hbb"]->Fill();
	  }
        //Zbb Box
        else if( mbbZ > 76.0 && mbbZ < 106.0 )
	  {
	    if(combineTrees){
	      box = Zbb;
	      razorTree->Fill();
	    }
	    else razorBoxes["Zbb"]->Fill();
	  }
        //HighRes Box
        else if( sigmaEOverE1 < 0.015 && sigmaEOverE2 < 0.015 )
	  {
	    if(combineTrees){
	      box = HighRes;
	      razorTree->Fill();
	    }
	    else razorBoxes["HighRes"]->Fill();
	  }
        //LowRes Box
        else
	  {
	    if(combineTrees){
	      box = LowRes;
	      razorTree->Fill();
	    }
	    else razorBoxes["LowRes"]->Fill();
	  }
	
    }//end of event loop
    
    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();
    
    outFile.Close();
}
