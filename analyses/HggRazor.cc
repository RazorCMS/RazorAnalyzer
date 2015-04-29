#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
//#include "RazorAuxPhoton.hh"

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

struct PhotonCandidate
{                                                  
  int   Index;
  TLorentzVector photon;
  float SigmaIetaIeta;                                                                        
  float R9;                                                                                  
  float HoverE;                                                                        
  float sumChargedHadronPt;                                                                
  float sumNeutralHadronEt;                                                     
  float sumPhotonEt;                                            
  float sigmaEOverE;
  bool  _passEleVeto;
  bool  _passIso;
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
    
    //Including Jet Corrections
    std::vector<JetCorrectorParameters> correctionParameters;
    
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L1FastJet_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L2Relative_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L3Absolute_AK5PF.txt"));
    correctionParameters.push_back(JetCorrectorParameters("data/FT53_V10_AN3_L2L3Residual_AK5PF.txt"));
    
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector( correctionParameters );
    
    //one tree to hold all events
    TTree *razorTree = new TTree("HggRazor", "Info on selected razor inclusive events");
    /*
      combine Trees
    */
    combineTrees = true;
      
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
    int n_Jets, nLooseBTaggedJets, nMediumBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
    float theMR;
    float theRsq;
    int nSelectedPhotons;
    float mGammaGamma, pTGammaGamma;
    float mbbZ, mbbH;
    HggRazorBox box;
    int run, event;

    //selected photon variables
    float Pho_E[2], Pho_Pt[2], Pho_Eta[2], Pho_Phi[2], Pho_SigmaIetaIeta[2], Pho_R9[2], Pho_HoverE[2];
    float Pho_sumChargedHadronPt[2], Pho_sumNeutralHadronEt[2], Pho_sumPhotonEt[2], Pho_sigmaEOverE[2];
    bool  Pho_passEleVeto[2], Pho_passIso[2];
    
    //jet information
    float jet_E[10], jet_Pt[10], jet_Eta[10], jet_Phi[10];
    
    //set branches on big tree
    if(combineTrees){
      razorTree->Branch("run", &run, "run/I");
      razorTree->Branch("event", &event, "event/I");
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
      
      razorTree->Branch("pho1E", &Pho_E[0], "pho1E/F");
      razorTree->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
      razorTree->Branch("Pho1Eta", &Pho_Eta[0], "pho1Eta/F");
      razorTree->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
      razorTree->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
      razorTree->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
      razorTree->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
      razorTree->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
      razorTree->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
      razorTree->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
      razorTree->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
      razorTree->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
      razorTree->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");
      
      razorTree->Branch("pho2E", &Pho_E[1], "pho2E/F");
      razorTree->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
      razorTree->Branch("Pho2Eta", &Pho_Eta[1], "pho2Eta/F");
      razorTree->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
      razorTree->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
      razorTree->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
      razorTree->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
      razorTree->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
      razorTree->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
      razorTree->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
      razorTree->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
      razorTree->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
      razorTree->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O)");
      
      razorTree->Branch("mbbZ", &mbbZ, "mbbZ/F");
      razorTree->Branch("mbbH", &mbbH, "mbbH/F");

      razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
      razorTree->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
      razorTree->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
      razorTree->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
      razorTree->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
      
      /*
      razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
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
      */
    }
    //set branches on all trees
    else{ 
        for(auto& box : razorBoxes){
	  box.second->Branch("run", &run, "run/I");
	  box.second->Branch("event", &event, "event/I");
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
	  
	  box.second->Branch("pho1E", &Pho_E[0], "pho1E/F");
	  box.second->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
	  box.second->Branch("Pho1Eta", &Pho_Eta[0], "pho1Eta/F");
	  box.second->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
	  box.second->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
	  box.second->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
	  box.second->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
	  box.second->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
	  box.second->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
	  box.second->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
	  box.second->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
	  box.second->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
	  box.second->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");
	  //box.second->Branch("", , "");
	  
	  box.second->Branch("pho2E", &Pho_E[1], "pho2E/F");
	  box.second->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
          box.second->Branch("Pho2Eta", &Pho_Eta[1], "pho2Eta/F");
          box.second->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
          box.second->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
          box.second->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
          box.second->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
          box.second->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
          box.second->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
          box.second->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
          box.second->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
          box.second->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
          box.second->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O");
	  
	  box.second->Branch("mbbZ", &mbbZ, "mbbZ/F");
	  box.second->Branch("mbbH", &mbbH, "mbbH/F");
	  
	  box.second->Branch("n_Jets", &n_Jets, "n_Jets/I");
	  box.second->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
	  box.second->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
	  box.second->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
	  box.second->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
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
        n_Jets = 0;
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
	run = runNum;
	event = eventNum;

	//selected photons variables
	for ( int i = 0; i < 2; i++ )
	  {
	    Pho_E[i]                  = -99.;
	    Pho_Pt[i]                 = -99.;
	    Pho_Eta[i]                = -99.;
	    Pho_Phi[i]                = -99.;
	    Pho_SigmaIetaIeta[i]      = -99.;
	    Pho_R9[i]                 = -99.;
	    Pho_HoverE[i]             = -99.;
	    Pho_sumChargedHadronPt[i] = -99.;
	    Pho_sumNeutralHadronEt[i] = -99.;
	    Pho_sumPhotonEt[i]        = -99.;
	    Pho_sigmaEOverE[i]        = -99.;
	    Pho_passEleVeto[i]        = false;
	    Pho_passIso[i]            = false;
	  }
	
	//jets
	for ( int i = 0; i < 10; i++ )
	  {
	    jet_E[i]   = -99.;
	    jet_Pt[i]  = -99.;
	    jet_Eta[i] = -99.;
	    jet_Phi[i] = -99.;
	  }
	
	//if ( !( run == 206859 && event == 24345 ) ) continue;
	
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
	std::vector< PhotonCandidate > phoCand;//PhotonCandidate defined in RazorAuxPhoton.hh
	
        int nPhotonsAbove40GeV = 0;
        for(int i = 0; i < nPhotons; i++){
	  //ID cuts -- apply isolation after candidate pair selection
	  //std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] << std::endl;
	  if ( !isGoodPhotonRun1( i , false ) )
	    {
	      //std::cout << "[INFO]: failed ID" << std::endl;
	      continue;
	    }
	  
	  //Defining Corrected Photon momentum
	  //float sf = pho_RegressionE[i]/phoE[i];
	  float pho_pt = phoPt[i];//Try to avoid using cosh to sync with Alex
	  //float pho_pt = phoPt[i]*sf;
	  //float pho_pt2 = pho_RegressionE[i]/cosh(phoEta[i]);
	  //std::cout << "pho# " << i << " phopt1: " << pho_pt << " pho_pt2: " << pho_pt2 << std::endl;
	  TVector3 vec;
	  //vec.SetPtEtaPhi( pho_pt, phoEta[i], phoPhi[i] );
	  vec.SetPtEtaPhi( phoPt[i], phoEta[i], phoPhi[i] );
	  
	  if ( pho_pt < 24.0 )
	    {
	      //std::cout << "[INFO]: failed pt" << std::endl;
	      continue;
	    }
	  
	  if( fabs(pho_superClusterEta[i]) > 2.5 ){
	    //allow photons in the endcap here, but if one of the two leading photons is in the endcap, reject the event
	    //std::cout << "[INFO]: failed eta" << std::endl;
	    continue; 
	  }
	  
	  if ( fabs(pho_superClusterEta[i]) > 1.4442 && fabs(pho_superClusterEta[i]) < 1.566 )
	    {
	      //Removing gap photons
	      //std::cout << "[INFO]: failed gap" << std::endl;
	      continue;
	    }
	  //photon passes
	  if( pho_pt > 32.0 ) nPhotonsAbove40GeV++;
	  //setting up photon 4-momentum with zero mass
	  TLorentzVector thisPhoton;
	  thisPhoton.SetVectM( vec, .0 );
	  //GoodPhotons.push_back( thisPhoton );
	  
	  //Filling Photon Candidate
	  PhotonCandidate tmp_phoCand;
	  tmp_phoCand.Index = i;
	  tmp_phoCand.photon = thisPhoton;
	  tmp_phoCand.SigmaIetaIeta = phoSigmaIetaIeta[i];
	  tmp_phoCand.R9 = phoR9[i];
	  tmp_phoCand.HoverE = pho_HoverE[i];
	  tmp_phoCand.sumChargedHadronPt = pho_sumChargedHadronPt[i];
	  tmp_phoCand.sumNeutralHadronEt = pho_sumNeutralHadronEt[i];
	  tmp_phoCand.sumPhotonEt = pho_sumPhotonEt[i];
	  tmp_phoCand.sigmaEOverE = pho_RegressionEUncertainty[i]/pho_RegressionE[i];
	  tmp_phoCand._passEleVeto = pho_passEleVeto[i];
	  tmp_phoCand._passIso = isGoodPhotonRun1( i , true );
	  phoCand.push_back( tmp_phoCand );
	  
	  nSelectedPhotons++;
	  //std::cout << "adding photon: " << nSelectedPhotons << std::endl;
        }
        //if there is no photon with pT above 40 GeV, reject the event
        if( nPhotonsAbove40GeV == 0 )
	  {
	    continue;
	  }
	
	if ( phoCand.size() < 2 )
	  {
	    continue;
	  }
	
	//find the "best" photon pair, higher Pt!
        TLorentzVector HiggsCandidate(0,0,0,0);
	int goodPhoIndex1 = -1;
        int goodPhoIndex2 = -1;
        double bestSumPt = -99.;
	for(size_t i = 0; i < phoCand.size(); i++){
	  for(size_t j = i+1; j < phoCand.size(); j++){//I like this logic better, I find it easier to understand
                PhotonCandidate pho1 = phoCand[i];
                PhotonCandidate pho2 = phoCand[j];
                
		//std::cout << "pho1Pt: " << pho1.photon.Pt() << " pho2Pt: " << pho2.photon.Pt() << std::endl;
                //need one photon in the pair to have pt > 40 GeV
                if( pho1.photon.Pt() < 40.0 && pho2.photon.Pt() < 40.0 ){
		  continue;
                }
                //need diphoton mass between > 100 GeV as in AN (April 1st)
                double diphotonMass = (pho1.photon + pho2.photon).M();
                if( diphotonMass < 100 ){
		  continue;
                }
                
                //if the sum of the photon pT's is larger than that of the current Higgs candidate, make this the Higgs candidate
                if( pho1.photon.Pt() + pho2.photon.Pt() > bestSumPt ){
		  bestSumPt = pho1.photon.Pt() + pho2.photon.Pt();
		  HiggsCandidate = pho1.photon + pho2.photon;
		  goodPhoIndex1 = pho1.Index;
		  goodPhoIndex2 = pho2.Index;  
                }
            }
        }   
	

	//Filling Selected Photon Information
	TLorentzVector pho_cand_vec[2];
	int _pho_index = 0;
        for ( auto& tmpPho : phoCand )
          {
            if ( !( tmpPho.Index == goodPhoIndex1 || tmpPho.Index == goodPhoIndex2 ) ) continue;
            if( _pho_index > 1 ) std::cerr << "[ERROR]: Photon index larger than 1!" << std::endl;
	    pho_cand_vec[_pho_index]           = tmpPho.photon;
	    Pho_E[_pho_index]                  = tmpPho.photon.E();
            Pho_Pt[_pho_index]                 = tmpPho.photon.Pt();
            Pho_Eta[_pho_index]                = tmpPho.photon.Eta();
            Pho_Phi[_pho_index]                = tmpPho.photon.Phi();
            Pho_SigmaIetaIeta[_pho_index]      = tmpPho.SigmaIetaIeta;
            Pho_R9[_pho_index]                 = tmpPho.R9;
            Pho_HoverE[_pho_index]             = tmpPho.HoverE;
            Pho_sumChargedHadronPt[_pho_index] = tmpPho.sumChargedHadronPt;
            Pho_sumNeutralHadronEt[_pho_index] = tmpPho.sumNeutralHadronEt;
            Pho_sumPhotonEt[_pho_index]        = tmpPho.sumPhotonEt;
            Pho_sigmaEOverE[_pho_index]        = tmpPho.sigmaEOverE;
            Pho_passEleVeto[_pho_index]        = tmpPho._passEleVeto;
            Pho_passIso[_pho_index]            = tmpPho._passIso;
            _pho_index++;
          }
	
	std::cout << "iso1: " << Pho_passIso[0] << " iso2: " << Pho_passIso[1] << std::endl; 
	//if the best candidate pair has pT < 20 GeV, reject the event
	if( HiggsCandidate.Pt() < 20.0 )
	  {
	    continue;
	  }
	
	//if the best candidate pair has a photon in the endcap, reject the event
	if ( fabs( Pho_Pt[0] ) > 1.44 || fabs( Pho_Pt[1] ) > 1.44 )
	  {
	    //allow for now, to sync with alex, probably good idea to keep them to debug
	    //continue;
	  }
       
	//if the best candidate pair has a non-isolated photon, reject the event
        if( !Pho_passIso[0] || !Pho_passIso[1] )
	  {
	    continue;
	  }
	//record higgs candidate info
        mGammaGamma = HiggsCandidate.M();
        pTGammaGamma = HiggsCandidate.Pt();
	

	//Jets
	vector<TLorentzVector> GoodJets;
        vector< pair<TLorentzVector, bool> > GoodCSVLJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not

	//I think I am selecting too many jets!
	//is the jet ID applied correctly
	//???
	//std::cout << "nJets: " << nJets << std::endl;
        for(int i = 0; i < nJets; i++){
	  //Jet Corrections                                                                      
	  double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						 fixedGridRhoAll, jetJetArea[i],
						 JetCorrector );
	  
	  TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
	  //std::cout << i << " pt: " << thisJet.Pt() << " eta: " << thisJet.Eta() << " phi: " << thisJet.Phi() << std::endl; 

	  if( thisJet.Pt() < 30.0 ) continue;//According to the April 1st 2015 AN
	  if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
          int level = 2; //loose jet ID
          if ( !jetPassIDLoose[i] ) continue;
          if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;
	  
	  //exclude selected photons from the jet collection
	  double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
	  if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN
	  
	  GoodJets.push_back(thisJet);
	  n_Jets++;
	  
	  /*
	   Change to isCSVL and isCSVM if you want CISV
	   */
	  if( isOldCSVL(i) ){
	    nLooseBTaggedJets++;
	    if( isOldCSVM(i) ){ 
	      nMediumBTaggedJets++;
	      GoodCSVLJets.push_back(make_pair(thisJet, true));
	    }
	    else{
	      GoodCSVLJets.push_back(make_pair(thisJet, false));
	    }
	  }
        }
	
        //if there are no good jets, reject the event
        if( n_Jets == 0 )
	  {
	    continue;
	  }
	
	int iJet = 0;
	for ( auto tmp_jet : GoodJets )
	  {
	    jet_E[iJet] = tmp_jet.E();
	    jet_Pt[iJet] = tmp_jet.Pt();
	    jet_Eta[iJet] = tmp_jet.Eta();
	    jet_Phi[iJet] = tmp_jet.Phi();
	    iJet++;
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
        if ( theMR < 150.0 )
	  {
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
	
	std::cout << "mbbH: " << mbbH << " mbbZ: " << mbbZ << std::endl;
	//Writing output to tree
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
        else if( Pho_sigmaEOverE[0] < 0.015 && Pho_sigmaEOverE[1] < 0.015 )
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
