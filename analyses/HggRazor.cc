#define RazorAnalyzer_cxx
//LOCAL INCLUDES
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
//ROOT INCLUDES
#include "TH1F.h"

using namespace std;

enum HggRazorBox {
    HighPt,
    Hbb,
    Zbb,
    HighRes,
    LowRes
};

enum HggRazorBoxSimple { HR, LR };

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

struct evt
{
  std::string run;
  std::string event;
};

#define _phodebug  0
#define _debug     1
#define _info      1

void RazorAnalyzer::HggRazor(string outFileName, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
  if ( _info) std::cout << "Initializing..." << std::endl;
  if (outFileName.empty()){
    if ( _info ) std::cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << std::endl;
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
  combineTrees = false;
  bool simpleBoxes = false;
  

  /*
    This is to debug and sync with Alex's events
  */
  std::ifstream ifs ( "HggRazorMissing.txt", std::fstream::in );
  std::string r_run, e_evt;
  std::map< std::string, evt > mymap;
  if ( ifs.is_open() )
    {
      while ( ifs.good() )
	{
	  ifs >> r_run >> e_evt;
	  std::string tmp = r_run + e_evt;
	  evt tmp_evt;
	  tmp_evt.run = r_run;
	  tmp_evt.event = e_evt;
	  if ( mymap.find( tmp ) == mymap.end() )
	    {
	      mymap[tmp] = tmp_evt;
	    }
	}
    }
  else
    {
      std::cout << "[ERROR]: unable to open file" << std::endl;
    }
  
  if( _info ) std::cout << "[INFO]: map size: " << mymap.size() << std::endl;
  
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
  
  map<string, TTree*> razorBoxesSimple;
  vector<string> boxNamesSimple;
  boxNamesSimple.push_back("HighRes");
  boxNamesSimple.push_back("LowRes");
  for(size_t i = 0; i < boxNamesSimple.size(); i++){
    razorBoxesSimple[boxNamesSimple[i]] = new TTree(boxNamesSimple[i].c_str(), boxNamesSimple[i].c_str());
  }
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
  //tree variables
  bool trigger;
  int n_Jets, nLooseBTaggedJets, nMediumBTaggedJets;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  float theMR;
  float theRsq, t1Rsq;
  float MET, t1MET;
  float mHem1, ptHem1, etaHem1, phiHem1, mHem2, ptHem2, etaHem2, phiHem2;
  int nSelectedPhotons;
  float mGammaGamma, pTGammaGamma, etaGammaGamma, phiGammaGamma;
  float mbbZ, mbbH;
  //HggRazorBox box;
  //HggRazorBoxSimple simplebox;
  unsigned int lumi, run, event;
  
  //selected photon variables
  float Pho_E[2], Pho_Pt[2], Pho_Eta[2], Pho_Phi[2], Pho_SigmaIetaIeta[2], Pho_R9[2], Pho_HoverE[2];
  float Pho_sumChargedHadronPt[2], Pho_sumNeutralHadronEt[2], Pho_sumPhotonEt[2], Pho_sigmaEOverE[2];
  bool  Pho_passEleVeto[2], Pho_passIso[2];
  
  //jet information
  float jet_E[10], jet_Pt[10], jet_Eta[10], jet_Phi[10];
  
  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("lumi", &lumi, "lumi/i");
    razorTree->Branch("run", &run, "run/i");
    razorTree->Branch("event", &event, "event/i");
    razorTree->Branch("trigger", &trigger, "trigger/O");
    razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
    razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
    razorTree->Branch("MET", &MET, "MET/F");
    razorTree->Branch("t1MET", &t1MET, "t1MET/F");
    razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
    razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
    razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
    razorTree->Branch("etaGammaGamma", &etaGammaGamma, "etaGammaGamma/F");
    razorTree->Branch("phiGammaGamma", &phiGammaGamma, "phiGammaGamma/F");
    
    razorTree->Branch("pho1E", &Pho_E[0], "pho1E/F");
    razorTree->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
    razorTree->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
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
    razorTree->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
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

    razorTree->Branch("mHem1", &mHem1, "mHem1/F");
    razorTree->Branch("ptHem1", &ptHem1, "ptHem1/F");
    razorTree->Branch("etaHem1", &etaHem1, "etaHem1/F");
    razorTree->Branch("phiHem1", &phiHem1, "phiHem1/F");

    razorTree->Branch("mHem2", &mHem2, "mHem2/F");
    razorTree->Branch("ptHem2", &ptHem2, "ptHem2/F");
    razorTree->Branch("etaHem2", &etaHem2, "etaHem2/F");
    razorTree->Branch("phiHem2", &phiHem2, "phiHem2/F");
  }
  //set branches on all trees
  else if ( !simpleBoxes ){ 
    for(auto& box : razorBoxes){
      box.second->Branch("lumi", &lumi, "lumi/i");
      box.second->Branch("run", &run, "run/i");
      box.second->Branch("event", &event, "event/i");
      box.second->Branch("trigger", &trigger, "trigger/O");
      box.second->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
      box.second->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
      box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      box.second->Branch("MR", &theMR, "MR/F");
      box.second->Branch("Rsq", &theRsq, "Rsq/F");
      box.second->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
      box.second->Branch("MET", &MET, "MET/F");
      box.second->Branch("t1MET", &t1MET, "t1MET/F");
      box.second->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
      box.second->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
      box.second->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
      box.second->Branch("etaGammaGamma", &etaGammaGamma, "etaGammaGamma/F");
      box.second->Branch("phiGammaGamma", &phiGammaGamma, "phiGammaGamma/F");
      
      box.second->Branch("pho1E", &Pho_E[0], "pho1E/F");
      box.second->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
      box.second->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
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
      
      box.second->Branch("pho2E", &Pho_E[1], "pho2E/F");
      box.second->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
      box.second->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
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

      box.second->Branch("mHem1", &mHem1, "mHem1/F");
      box.second->Branch("ptHem1", &ptHem1, "ptHem1/F");
      box.second->Branch("etaHem1", &etaHem1, "etaHem1/F");
      box.second->Branch("phiHem1", &phiHem1, "phiHem1/F");

      box.second->Branch("mHem2", &mHem2, "mHem2/F");
      box.second->Branch("ptHem2", &ptHem2, "ptHem2/F");
      box.second->Branch("etaHem2", &etaHem2, "etaHem2/F");
      box.second->Branch("phiHem2", &phiHem2, "phiHem2/F");
      
    } 
  }
  else
    {
      for ( auto& simplebox : razorBoxesSimple )
	{
	  simplebox.second->Branch("lumi", &lumi, "lumi/i");
	  simplebox.second->Branch("run", &run, "run/i");
	  simplebox.second->Branch("event", &event, "event/i");
	  simplebox.second->Branch("trigger", &trigger, "trigger/O");
	  simplebox.second->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
	  simplebox.second->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
	  simplebox.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
	  simplebox.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
	  simplebox.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
	  simplebox.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
	  simplebox.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
	  simplebox.second->Branch("MR", &theMR, "MR/F");
	  simplebox.second->Branch("Rsq", &theRsq, "Rsq/F");
	  simplebox.second->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
	  simplebox.second->Branch("MET", &MET, "MET/F");
	  simplebox.second->Branch("t1MET", &t1MET, "t1MET/F");
	  simplebox.second->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
	  simplebox.second->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
	  simplebox.second->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
	  simplebox.second->Branch("etaGammaGamma", &etaGammaGamma, "etaGammaGamma/F");
	  simplebox.second->Branch("phiGammaGamma", &phiGammaGamma, "phiGammaGamma/F");

	  simplebox.second->Branch("pho1E", &Pho_E[0], "pho1E/F");
	  simplebox.second->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
	  simplebox.second->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
	  simplebox.second->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
	  simplebox.second->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
	  simplebox.second->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
	  simplebox.second->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
	  simplebox.second->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
	  simplebox.second->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
	  simplebox.second->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
	  simplebox.second->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
	  simplebox.second->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
	  simplebox.second->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");

	  simplebox.second->Branch("pho2E", &Pho_E[1], "pho2E/F");
	  simplebox.second->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
	  simplebox.second->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
	  simplebox.second->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
	  simplebox.second->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
	  simplebox.second->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
	  simplebox.second->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
	  simplebox.second->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
	  simplebox.second->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
	  simplebox.second->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
	  simplebox.second->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
	  simplebox.second->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
	  simplebox.second->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O");
	  
	  simplebox.second->Branch("mbbZ", &mbbZ, "mbbZ/F");
	  simplebox.second->Branch("mbbH", &mbbH, "mbbH/F");

	  simplebox.second->Branch("n_Jets", &n_Jets, "n_Jets/I");
	  simplebox.second->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
	  simplebox.second->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
	  simplebox.second->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
	  simplebox.second->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");

	  simplebox.second->Branch("mHem1", &mHem1, "mHem1/F");
	  simplebox.second->Branch("ptHem1", &ptHem1, "ptHem1/F");
	  simplebox.second->Branch("etaHem1", &etaHem1, "etaHem1/F");
	  simplebox.second->Branch("phiHem1", &phiHem1, "phiHem1/F");
	  
	  simplebox.second->Branch("mHem2", &mHem2, "mHem2/F");
	  simplebox.second->Branch("ptHem2", &ptHem2, "ptHem2/F");
	  simplebox.second->Branch("etaHem2", &etaHem2, "etaHem2/F");
	  simplebox.second->Branch("phiHem2", &phiHem2, "phiHem2/F");
	  
	}
    }
  
  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //begin event
    if( _info && (jentry % 10000 == 0) ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //fill normalization histogram
    NEvents->Fill(1.0);
    
    //reset tree variables
    trigger = false;
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
    mGammaGamma      = -1;
    pTGammaGamma     = -1;
    etaGammaGamma    = -1;
    phiGammaGamma    = -1;
    mbbZ = 0;
    mbbH = 0;
    lumi  = lumiNum;
    run   = runNum;
    event = eventNum;
    
    mHem1   = -1;
    ptHem1  = -1;
    etaHem1 = -1;
    phiHem1 = -1;

    mHem2   = -1;
    ptHem2  = -1;
    etaHem2 = -1;
    phiHem2 = -1;

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
    
    /*
    std::stringstream ss;
    ss << run << event;
    if ( mymap.find( ss.str() ) == mymap.end() )continue;
    //if ( !( run == 206859 && event == 24345 ) ) continue;
    */
    if ( _debug ) std::cout << "=======new event=====" << std::endl;
    if ( _debug ) std::cout << "run == " << run << " && evt == " << event << std::endl;
    
    //if(combineTrees) box = LowRes;
    
    //A p p l y   P h o t o n   T r i g g e r s 
    //-----------------------------------------
    trigger = HLTDecision[25] || HLTDecision[26] || HLTDecision[27] || HLTDecision[28];
    //if ( !trigger ) continue;
    //Save trigger decision in tree instead of rejecting events (this needs to be done at the ntuple level)
    
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
      
      if(isTightElectron(i))
	{ 
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
      if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] << std::endl;
      if ( !isGoodPhotonRun1( i , false, false ) )
	{
	  if ( _phodebug ) std::cout << "[DEBUG]: failed ID" << std::endl;
	  continue;
	}
      
      //Defining Corrected Photon momentum
      //float pho_pt = phoPt[i];//nominal pt
      float pho_pt_corr = pho_RegressionE[i]/cosh(phoEta[i]);//regression corrected pt
      TVector3 vec;
      //vec.SetPtEtaPhi( pho_pt, phoEta[i], phoPhi[i] );
      vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );
      
      if ( phoPt[i] < 24.0 )
      //if ( phoE[i]/cosh( phoEta[i] ) < 24.0 )
	{
	  if ( _phodebug ) std::cout << "[DEBUG]: failed pt" << std::endl;
	  continue;
	}
      
      if( fabs(pho_superClusterEta[i]) > 2.5 ){
	//allow photons in the endcap here, but if one of the two leading photons is in the endcap, reject the event
	if ( _phodebug ) std::cout << "[DEBUG]: failed eta" << std::endl;
	continue; 
      }
      
      if ( fabs(pho_superClusterEta[i]) > 1.4442 && fabs(pho_superClusterEta[i]) < 1.566 )
	{
	  //Removing gap photons
	  if ( _phodebug ) std::cout << "[INFO]: failed gap" << std::endl;
	  continue;
	}
      //photon passes
      if( phoPt[i] > 32.0 ) nPhotonsAbove40GeV++;
      //setting up photon 4-momentum with zero mass
      TLorentzVector thisPhoton;
      thisPhoton.SetVectM( vec, .0 );
      
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
      tmp_phoCand._passIso = isGoodPhotonRun1( i , true, false );
      phoCand.push_back( tmp_phoCand );
      
      nSelectedPhotons++;
    }
    //if there is no photon with pT above 40 GeV, reject the event
    if( nPhotonsAbove40GeV == 0 )
      {
	if ( _debug ) std::cout << "[DEBUG]: no photons above 40 GeV, nphotons: " << phoCand.size() << std::endl;
	continue;
      }
    
    if ( phoCand.size() < 2 )
      {
	if ( _debug ) std::cout << "[INFO]: not enough photon, nphotons: " << phoCand.size() << std::endl;
	for(int i = 0; i < nPhotons; i++)
	  {
	    if ( _debug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] 
				    << " SIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
	    isGoodPhotonRun1( i , false, _debug );
	  }
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
        if ( _debug )
	  {
	    std::cout << "[DEBUG]: pho1-> " << pho1.photon.Pt()
		      << " [DEBUG]: pho2->" << pho2.photon.Pt() 
		      << std::endl;
	  }
	//need one photon in the pair to have pt > 40 GeV
	if ( pho1.photon.Pt() < 40.0 && pho2.photon.Pt() < 40.0 )
	  {
	    if ( _debug ) std::cout << "[DEBUG]: both photons failed PT > 40 GeV" << std::endl; 
	    continue;
	  }
	//need diphoton mass between > 100 GeV as in AN (April 1st)
	double diphotonMass = (pho1.photon + pho2.photon).M();
	if ( _debug )
	  {
	    std::cout << "[DEBUG] Hgg cadidate pT: " << pho1.photon.Pt() + pho2.photon.Pt() << std::endl;
	  }
	
	if( diphotonMass < 100 )
	  {
	    if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 100 GeV: " << std::endl;
	    if ( _debug ) std::cout << "... pho1Pt: " << pho1.photon.Pt()  << " pho2Pt: " << pho2.photon.Pt()  << std::endl;
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
    
    if ( _debug )
      {
	std::cout << "[DEBUG]: best photon pair: " 
		  << "\n-> pho1Pt: " << Pho_Pt[0] 
		  << "\n-> pho2Pt: " << Pho_Pt[1] 
		  << std::endl;
      }

    //if the best candidate pair has pT < 20 GeV, reject the event
    if( HiggsCandidate.Pt() < 20.0 )
      {
	if ( _debug ) std::cout << "[DEBUG]: Higgs Pt < 20 GeV, H pt: " << HiggsCandidate.Pt() << std::endl; 
	continue;
      }
    
    //if the best candidate pair has a photon in the endcap, reject the event
    if ( fabs( Pho_Eta[0] ) > 1.44 || fabs( Pho_Eta[1] ) > 1.44 )
      {
	//allow for now, to sync with alex, probably good idea to keep them to debug
	//continue;
      }
    
    //if the best candidate pair has a non-isolated photon, reject the event
    if( !Pho_passIso[0] || !Pho_passIso[1] )
      {
	if ( _debug ) std::cout << "[DEBUG]: Failed ISO: pho1, pho2: " << Pho_passIso[0] << ", " << Pho_passIso[1] << std::endl;
	if ( _debug ) std::cout << "[DEBUG]: pho1Pt: " << Pho_Pt[0] << " pho2Pt: " << Pho_Pt[1] << std::endl;
	for ( auto& phoC : phoCand )
	  {
	    if ( _debug ) std::cout << "===> phopt: " << phoC.photon.Pt() << " phoEta: " << phoC.photon.Eta() << std::endl;
	    photonPassIsoRun1( phoC.Index, _debug );
	  }
	continue;
      }
    //record higgs candidate info
    mGammaGamma = HiggsCandidate.M();
    pTGammaGamma = HiggsCandidate.Pt();
    etaGammaGamma = HiggsCandidate.Eta();
    phiGammaGamma = HiggsCandidate.Phi();
    
    
    //Jets
    vector<TLorentzVector> GoodJets;
    vector< pair<TLorentzVector, bool> > GoodCSVLJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not
    
    for(int i = 0; i < nJets; i++){
      //Jet Corrections                                                                      
      double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
					      fixedGridRhoAll, jetJetArea[i],
					      JetCorrector );
      
      TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
      
      if( thisJet.Pt() < 30.0 ) continue;//According to the April 1st 2015 AN
      if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
      //int level = 2; //3rd bit of jetPileupIdFlag
      if ( !jetPassIDLoose[i] ) continue;
      //if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;
      
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
	if ( _debug ) std::cout << "[DEBUG]: No Jets Selected" << std::endl;
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
    for( auto& jet : GoodJets ) JetsPlusHiggsCandidate.push_back(jet);
    JetsPlusHiggsCandidate.push_back(HiggsCandidate);
    
    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
    TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType0Plus1Pt, 0, metType0Plus1Phi, 0 );
    
    vector<TLorentzVector> hemispheres = getHemispheres(JetsPlusHiggsCandidate);
    theMR  = computeMR(hemispheres[0], hemispheres[1]); 
    theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
    t1Rsq  = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
    MET = metPt;
    t1MET = metType0Plus1Pt;
    //hem1
    mHem1   = hemispheres[0].M();
    ptHem1  = hemispheres[0].Pt();
    etaHem1 = hemispheres[0].Eta();
    phiHem1 = hemispheres[0].Phi();
    //hem2
    mHem2   = hemispheres[1].M();
    ptHem2  = hemispheres[1].Pt();
    etaHem2 = hemispheres[1].Eta();
    phiHem2 = hemispheres[1].Phi();
    //No MR cut
    
    /*
      if ( theMR < 150.0 )
      {
      if ( _debug ) std::cout << "[INFO]: MR < 150 GeV, MR: " << theMR << std::endl;
      for ( auto& jet : JetsPlusHiggsCandidate )
      {
      if ( _debug ) std::cout << "phoPT: " << pTGammaGamma 
      << " jet pt : " << jet.Pt() << " eta: " << jet.Eta() << " phi: " << jet.Phi() 
      << " h1 pt: " << hemispheres[0].Pt() << " h1 eta: " << hemispheres[0].Eta()
      << " h2 pt: " << hemispheres[1].Pt() << " h2 eta: " << hemispheres[1].Eta() << std::endl;
      }
      continue;
      }
    */

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
    
    if ( _debug ) std::cout << "mbbH: " << mbbH << " mbbZ: " << mbbZ << std::endl;
    //Writing output to tree
    if ( !simpleBoxes )
      {
	//HighPt Box
	if ( pTGammaGamma > 110.0 )
	  {
	    if(combineTrees)
	      {
		//box = HighPt;
		razorTree->Fill();
	      }
	    else razorBoxes["HighPt"]->Fill();
	  }
	//Hbb Box
	else if ( mbbH > 110.0 && mbbH < 140.0 )
	  {
	    if(combineTrees)
	      {
		//box = Hbb;
		razorTree->Fill();
	      }
	    else razorBoxes["Hbb"]->Fill();
	  }
	//Zbb Box
	else if( mbbZ > 76.0 && mbbZ < 106.0 )
	  {
	    if(combineTrees)
	      {
		//box = Zbb;
		razorTree->Fill();
	      }
	    else razorBoxes["Zbb"]->Fill();
	  }
	//HighRes Box
	else if( Pho_sigmaEOverE[0] < 0.015 && Pho_sigmaEOverE[1] < 0.015 )
	  {
	    if(combineTrees)
	      {
		//box = HighRes;
		razorTree->Fill();
	      }
	    else razorBoxes["HighRes"]->Fill();
	  }
	//LowRes Box
	else
	  {
	    if(combineTrees)
	      {
		//box = LowRes;
		razorTree->Fill();
	      }
	    else razorBoxes["LowRes"]->Fill();
	  }
      }
    else
      {
	if( Pho_sigmaEOverE[0] < 0.015 && Pho_sigmaEOverE[1] < 0.015 )
	  {
	    razorBoxesSimple["HighRes"]->Fill();
	  }
	else
	  {
	    razorBoxesSimple["LowRes"]->Fill();
	  }
      }
  }//end of event loop
  
  if ( _info ) std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;
  if ( _info ) std::cout << "[INFO]: Writing output trees..." << std::endl;
  if( !simpleBoxes )
    {
      if( combineTrees )
	{
	  razorTree->Write();
	}
      else
	{
	  for ( auto& box : razorBoxes ) box.second->Write();
	}
    }
  else
    {
      for( auto& box : razorBoxesSimple ) box.second->Write();
    }
  NEvents->Write();
  
  outFile.Close();
}
