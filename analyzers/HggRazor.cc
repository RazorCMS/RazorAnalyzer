//LOCAL INCLUDES
#include "HggRazor.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"

using namespace std;

enum HggRazorBox {
  HighPt = 0,
  Hbb = 1,
  Zbb = 2,
  HighRes = 3,
  LowRes = 4
};

struct PhotonCandidate
{                                                  
  int   Index;
  TLorentzVector photon;
  TLorentzVector photonSC;
  float scEta;
  float scPhi;
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

#define _phodebug 0
#define _debug    0
#define _info     1

const double EB_R = 129.0;
const double EE_Z = 317.0;

const double JET_CUT = 30.;
const int NUM_PDF_WEIGHTS = 60;

//Testing branching and merging
void HggRazor::Analyze(bool isData, int option, string outFileName, string label)
{
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  bool doPhotonScaleCorrection = true;

  string analysisTag = "";
  if (option >= 0 && option <= 9) analysisTag = "2015_76X";
  else if (option >= 10 && option <= 19) analysisTag = "2016_80X";
  else {
    cout << "Error: analysisoption == " << option << " is not supported. Exiting.\n";
    return;
  }

  bool doMRSkim = false;
  bool doRequireIso = false;
  bool doEleVeto = false;

  if (option == 1 || option == 2 || option == 3 || 
      option == 11 || option == 12 || option == 13) doRequireIso = true;
  
  if (option == 2 || option == 3 || 
      option == 12 || option == 13) doMRSkim = true;

  if (option ==3 || option == 13) doEleVeto = true;


  std::cout << "[INFO]: option = " << option << std::endl;
  std::cout << "[INFO]: analysisTag --> " << analysisTag << std::endl;
  std::cout << "[INFO]: doEleVeto --> " << doEleVeto << std::endl;
  std::cout << "[INFO]: doPhotonScaleCorrection --> " << doPhotonScaleCorrection << std::endl;
 
  //initialization: create one TTree for each analysis box 
  if ( _info ) std::cout << "Initializing..." << std::endl;
  
  if ( outFileName.empty() )
    {
      if ( _info ) std::cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << std::endl;
      outFileName = "HggRazor.root";
    }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *razorTree = new TTree("HggRazor", "Info on selected razor inclusive events");
  
  //Get CMSSW Directory
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  if (analysisTag == "2015_76X") helper = new RazorHelper("Razor2015_76X", isData, false);
  else if (analysisTag == "2016_80X") helper = new RazorHelper("Razor2016_80X", isData, false);
  else helper = new RazorHelper("", isData, false);
  
  //--------------------------------
  //Photon Trigger Efficiency
  //--------------------------------
  TFile *photonTriggerEffFile_LeadingLeg = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016_Golden_2p6ifb/PhoHLTLeadingLegEffDenominatorLoose.root");
  TFile *photonTriggerEffFile_TrailingLeg = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016_Golden_2p6ifb/PhoHLTTrailingLegEffDenominatorLoose.root");
  TH2D *photonTriggerEffHist_LeadingLeg = (TH2D*)photonTriggerEffFile_LeadingLeg->Get("hEffEtaPt");
  TH2D *photonTriggerEffHist_TrailingLeg = (TH2D*)photonTriggerEffFile_TrailingLeg->Get("hEffEtaPt");
  assert(photonTriggerEffHist_LeadingLeg);
  assert(photonTriggerEffHist_TrailingLeg);

  //--------------------------------
  //Photon Energy Scale and Resolution Corrections
  //--------------------------------
  std::string photonCorrectionPath = "";
  if ( cmsswPath != NULL ) photonCorrectionPath = string(cmsswPath) + "/src/RazorAnalyzer/data/PhotonCorrections/";

  EnergyScaleCorrection_class *photonCorrector = 0;
  if (analysisTag == "2015_76X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/76X_16DecRereco_2015", photonCorrectionPath.c_str()));
  } else if (analysisTag == "2016_80X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/80X_2016", photonCorrectionPath.c_str()));
  }

  if(!isData) {
    photonCorrector->doScale = false; 
    photonCorrector->doSmearings = true;
  } else {
    photonCorrector->doScale = true; 
    photonCorrector->doSmearings = false;
  }

  //--------------------------------
  //Including Jet Energy Corrections
  //--------------------------------
  FactorizedJetCorrector *JetCorrector = helper->getJetCorrector();

  //---------------
  //btag efficiency
  //---------------
  //Medium
  TH2D *btagMediumEfficiencyHist = 0;
  TH2D *btagMediumCharmEfficiencyHist = 0;
  TH2D *btagMediumLightJetsEfficiencyHist = 0;
  //Loose
  TH2D *btagLooseEfficiencyHist = 0;
  TH2D *btagLooseCharmEfficiencyHist = 0;
  TH2D *btagLooseLightJetsEfficiencyHist = 0;
  if ( !isData )
    {
      //Medium
      TFile *btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/BTagEffFastsimToFullsimCorrectionFactors.root");
      btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("BTagEff_Medium_Fullsim");
      assert(btagMediumEfficiencyHist);
      TFile *btagCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/CharmJetBTagEffFastsimToFullsimCorrectionFactors.root");
      btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("BTagEff_Medium_Fullsim");
      assert(btagMediumCharmEfficiencyHist);
      TFile *btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/LightJetBTagEffFastsimToFullsimCorrectionFactors.root");
      btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("BTagEff_Medium_Fullsim");
      assert(btagMediumLightJetsEfficiencyHist);
      //Loose
      TFile *btagLooseEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVL_Fullsim.root");
      btagLooseEfficiencyHist = (TH2D*)btagLooseEfficiencyFile->Get("Efficiency_PtEta");
      assert(btagLooseEfficiencyHist);
      TFile *btagLooseCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVL_Fullsim.root");
      btagLooseCharmEfficiencyHist = (TH2D*)btagLooseCharmEfficiencyFile->Get("Efficiency_PtEta");
      assert(btagLooseCharmEfficiencyHist);
      TFile *btagLooseLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVL_Fullsim.root");
      btagLooseLightJetsEfficiencyHist = (TH2D*)btagLooseLightJetsEfficiencyFile->Get("Efficiency_PtEta");
      assert(btagLooseLightJetsEfficiencyHist);
    }

  //-----------------------
  //B-tagging scale factors
  //-----------------------
  
  string bTagPathname = "";
  if ( cmsswPath != NULL ) bTagPathname = string(cmsswPath) + "/src/RazorAnalyzer/data/ScaleFactors/";
  else bTagPathname = "data/ScaleFactors/";
  //Fullsim
  BTagCalibration btagcalib("csvv2", Form("%s/CSVv2_76X.csv",bTagPathname.c_str()));
  //Medium WP
  BTagCalibrationReader btagreaderM(&btagcalib,           //calibration instance
				    BTagEntry::OP_MEDIUM, //operating point
				    "mujets",             //measurement type
				    "central");           //systematics type
  BTagCalibrationReader btagreaderM_up(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");   //sys up
  BTagCalibrationReader btagreaderM_do(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down"); //sys down
  BTagCalibrationReader btagreaderMistagM(&btagcalib,            //calibration instance
					  BTagEntry::OP_MEDIUM,  //operating point
					  "comb",                //measurement type
					  "central");            //systematics type
  BTagCalibrationReader btagreaderMistagM_up(&btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    //sys up
  BTagCalibrationReader btagreaderMistagM_do(&btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  //sys down
  //Loose WP
  BTagCalibrationReader btagreaderL(&btagcalib,           //calibration instance
				    BTagEntry::OP_LOOSE,  //operating point
				    "mujets",             //measurement type
				    "central");           //systematics type
  BTagCalibrationReader btagreaderL_up(&btagcalib, BTagEntry::OP_LOOSE, "mujets", "up");  //sys up
  BTagCalibrationReader btagreaderL_do(&btagcalib, BTagEntry::OP_LOOSE, "mujets", "down");  //sys down
  BTagCalibrationReader btagreaderMistagL(&btagcalib,           //calibration instance
					  BTagEntry::OP_LOOSE,  //operating point
					  "comb",               //measurement type
					  "central");           //systematics type
  BTagCalibrationReader btagreaderMistagL_up(&btagcalib, BTagEntry::OP_LOOSE, "comb", "up");    //sys up
  BTagCalibrationReader btagreaderMistagL_do(&btagcalib, BTagEntry::OP_LOOSE, "comb", "down");  //sys down

  //----------
  //pu histo
  //----------
  TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  
  //separate trees for individual boxes
  map<string, TTree*> razorBoxes;
  vector<string> boxNames;
  boxNames.push_back("HighPt");
  boxNames.push_back("Hbb");
  boxNames.push_back("Zbb");
  boxNames.push_back("HighRes");
  boxNames.push_back("LowRes");
  for ( size_t i = 0; i < boxNames.size(); i++)
    {
      razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }
  
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

  //--------------
  //tree variables
  //--------------
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float triggerEffWeight;
  float ISRSystWeightUp, ISRSystWeightDown;
  //For btag scale factor uncertainty
  float btagCorrFactor;
  float sf_btagUp, sf_btagDown;
  float sf_bmistagUp, sf_bmistagDown;
  //For scale variation uncertainties
  float sf_facScaleUp, sf_facScaleDown;
  float sf_renScaleUp, sf_renScaleDown;
  float sf_facRenScaleUp, sf_facRenScaleDown;
  //PDF SF
  std::vector<float> sf_pdf;
  
  int NPU;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  float theMR, theMR_JESUp, theMR_JESDown;
  float theRsq, theRsq_JESUp, theRsq_JESDown, t1Rsq, t1Rsq_JESUp, t1Rsq_JESDown;
  float MET, MET_JESUp, MET_JESDown, t1MET, t1MET_JESUp, t1MET_JESDown;
  float HT;


  int nSelectedPhotons;
  float mGammaGamma, pTGammaGamma, mGammaGammaSC, pTGammaGammaSC, sigmaMoverM;
  float mbbZ, mbbZ_L, mbbH, mbbH_L;
  bool passedDiphotonTrigger;
  HggRazorBox razorbox = LowRes;
  
  unsigned int run, lumi, event;
  
  //selected photon variables
  float Pho_E[2], Pho_Pt[2], Pho_Eta[2], Pho_Phi[2], Pho_SigmaIetaIeta[2], Pho_R9[2], Pho_HoverE[2];
  float PhoSC_E[2], PhoSC_Pt[2], PhoSC_Eta[2], PhoSC_Phi[2];
  float Pho_sumChargedHadronPt[2], Pho_sumNeutralHadronEt[2], Pho_sumPhotonEt[2], Pho_sigmaEOverE[2];
  bool  Pho_passEleVeto[2], Pho_passIso[2];
  int   Pho_motherID[2];

  //jet information
  int n_Jets, nLooseBTaggedJets, nMediumBTaggedJets;
  int n_Jets_JESUp, n_Jets_JESDown; 
  float jet_E[50], jet_Pt[50], jet_Eta[50], jet_Phi[50];
  bool jetIsCSVL[50], jetIsCSVM[50], jetIsCSVT[50];

  //------------------------
  //set branches on big tree
  //------------------------

  razorTree->Branch("weight", &weight, "weight/F");
  razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  razorTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  razorTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  razorTree->Branch("triggerEffWeight", &triggerEffWeight, "triggerEffWeight/F");
  razorTree->Branch("ISRSystWeightUp", &ISRSystWeightUp, "ISRSystWeightUp/F");
  razorTree->Branch("ISRSystWeightDown", &ISRSystWeightDown, "ISRSystWeightDown/F");
      
  razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
  razorTree->Branch("sf_btagUp", &sf_btagUp, "sf_btagUp/F");
  razorTree->Branch("sf_btagDown", &sf_btagDown, "sf_btagDown/F");
  razorTree->Branch("sf_bmistagUp", &sf_bmistagUp, "sf_bmistagUp/F");
  razorTree->Branch("sf_bmistagDown", &sf_bmistagDown, "sf_bmistagDown/F");
      
  razorTree->Branch("sf_facScaleUp", &sf_facScaleUp, "sf_facScaleUp/F");
  razorTree->Branch("sf_facScaleDown", &sf_facScaleDown, "sf_facScaleDown/F");
  razorTree->Branch("sf_renScaleUp", &sf_renScaleUp, "sf_renScaleUp/F");
  razorTree->Branch("sf_renScaleDown", &sf_renScaleDown, "sf_renScaleDown/F");
  razorTree->Branch("sf_facRenScaleUp", &sf_facRenScaleUp, "sf_facRenScaleUp/F");
  razorTree->Branch("sf_facRenScaleDown", &sf_facRenScaleDown, "sf_facRenScaleDown/F");
  razorTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights); //get PDF weights directly from RazorEvents
  razorTree->Branch("sf_pdf", "std::vector<float>",&sf_pdf); //sf PDF
      
  //MET filters
  razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  razorTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  razorTree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
  razorTree->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
  razorTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  razorTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  razorTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  razorTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  razorTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  razorTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  razorTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  razorTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  razorTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  razorTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  razorTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  razorTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
      
  razorTree->Branch("run", &run, "run/i");
  razorTree->Branch("lumi", &lumi, "lumi/i");
  razorTree->Branch("event", &event, "event/i");
  razorTree->Branch("passedDiphotonTrigger", &passedDiphotonTrigger, "passedDiphotonTrigger/O");
  razorTree->Branch("NPU", &NPU, "npu/i");
  razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
  razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
  razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
  razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
  razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
  razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
  razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
  razorTree->Branch("MR", &theMR, "MR/F");
  razorTree->Branch("MR_JESUp", &theMR_JESUp, "MR_JESUp/F");
  razorTree->Branch("MR_JESDown", &theMR_JESDown, "MR_JESDown/F");
  razorTree->Branch("Rsq", &theRsq, "Rsq/F");
  razorTree->Branch("Rsq_JESUp", &theRsq_JESUp, "Rsq_JESUp/F");
  razorTree->Branch("Rsq_JESDown", &theRsq_JESDown, "Rsq_JESDown/F");
  razorTree->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
  razorTree->Branch("t1Rsq_JESUp", &t1Rsq_JESUp, "t1Rsq_JESUp/F");
  razorTree->Branch("t1Rsq_JESDown", &t1Rsq_JESDown, "t1Rsq_JESDown/F");
  razorTree->Branch("MET", &MET, "MET/F");
  razorTree->Branch("MET_JESUp", &MET_JESUp, "MET_JESUp/F");
  razorTree->Branch("MET_JESDown", &MET_JESDown, "MET_JESDown/F");
  razorTree->Branch("t1MET", &t1MET, "t1MET/F");
  razorTree->Branch("t1MET_JESUp", &t1MET_JESUp, "t1MET_JESUp/F");
  razorTree->Branch("t1MET_JESDown", &t1MET_JESDown, "t1MET_JESDown/F");
  razorTree->Branch("HT", &HT, "HT/F");

  razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
  razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
  razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
  razorTree->Branch("mGammaGammaSC", &mGammaGammaSC, "mGammaGammaSC/F");
  razorTree->Branch("pTGammaGammaSC", &pTGammaGammaSC, "pTGammaGammaSC/F");
  razorTree->Branch("sigmaMoverM", &sigmaMoverM, "sigmaMoverM/F");
  razorTree->Branch("box", &razorbox, "box/I");
      
  razorTree->Branch("pho1E", &Pho_E[0], "pho1E/F");
  razorTree->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
  razorTree->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
  razorTree->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
  razorTree->Branch("pho1SC_E", &PhoSC_E[0], "pho1SC_E/F");
  razorTree->Branch("pho1SC_Pt", &PhoSC_Pt[0], "pho1SC_Pt/F");
  razorTree->Branch("pho1SC_Eta", &PhoSC_Eta[0], "pho1SC_Eta/F");
  razorTree->Branch("pho1SC_Phi", &PhoSC_Phi[0], "pho1SC_Phi/F");
  razorTree->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
  razorTree->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
  razorTree->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
  razorTree->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
  razorTree->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
  razorTree->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
  razorTree->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
  razorTree->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
  razorTree->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");
  razorTree->Branch("pho1MotherID", &Pho_motherID[0], "pho1MotherID/I");
      
  razorTree->Branch("pho2E", &Pho_E[1], "pho2E/F");
  razorTree->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
  razorTree->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
  razorTree->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
  razorTree->Branch("pho2SC_E", &PhoSC_E[1], "pho2SC_E/F");
  razorTree->Branch("pho2SC_Pt", &PhoSC_Pt[1], "pho2SC_Pt/F");
  razorTree->Branch("pho2SC_Eta", &PhoSC_Eta[1], "pho2SC_Eta/F");
  razorTree->Branch("pho2SC_Phi", &PhoSC_Phi[1], "pho2SC_Phi/F");
  razorTree->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
  razorTree->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
  razorTree->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
  razorTree->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
  razorTree->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
  razorTree->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
  razorTree->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
  razorTree->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
  razorTree->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O)");
  razorTree->Branch("pho2MotherID", &Pho_motherID[1], "pho2MotherID/I");
      
  razorTree->Branch("mbbZ", &mbbZ, "mbbZ/F");
  razorTree->Branch("mbbH", &mbbH, "mbbH/F");
  razorTree->Branch("mbbZ_L", &mbbZ_L, "mbbZ_L/F");
  razorTree->Branch("mbbH_L", &mbbH_L, "mbbH_L/F");
      
  razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
  razorTree->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
  razorTree->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
  razorTree->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
  razorTree->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
  razorTree->Branch("jetIsCSVL", jetIsCSVL, "jetIsCSVL[n_Jets]/O");
  razorTree->Branch("jetIsCSVM", jetIsCSVM, "jetIsCSVM[n_Jets]/O");
  razorTree->Branch("jetIsCSVT", jetIsCSVT, "jetIsCSVT[n_Jets]/O");
  razorTree->Branch("n_Jets_JESUp", &n_Jets_JESUp, "n_Jets_JESUp/I");
  razorTree->Branch("n_Jets_JESDown", &n_Jets_JESDown, "n_Jets_JESDown/I");
  razorTree->Branch("HLTDecision", HLTDecision, "HLTDecision[300]/O");
      
  //GenParticles
  razorTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  razorTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  razorTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  razorTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  razorTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  razorTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  razorTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  razorTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  razorTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");


  //begin loop
  if ( fChain == 0 ) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for ( Long64_t jentry=0; jentry < nentries; jentry++ )
    {
      //begin event
      if( _info && (jentry % 10000 == 0) ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      Long64_t ientry = LoadTree( jentry );
      if ( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
    
      //fill normalization histogram    
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      weight = genWeight;
      SumWeights->Fill(1.0, weight);
      
      //reset tree variables
      ISRSystWeightUp   = 1.0;
      ISRSystWeightDown = 1.0;
      pileupWeight      = 1.0;
      pileupWeightUp    = 1.0;
      pileupWeightDown  = 1.0;
      triggerEffWeight  = 1.0;
      
      btagCorrFactor    = 1.0;
      sf_btagUp         = 1.0;
      sf_btagDown       = 1.0;
      sf_bmistagUp      = 1.0;
      sf_bmistagDown    = 1.0;
      
      sf_facScaleUp = 1.0;
      sf_facScaleDown = 1.0;
      sf_renScaleUp = 1.0;
      sf_renScaleDown = 1.0;
      sf_facRenScaleUp = 1.0;
      sf_facRenScaleDown = 1.0;
      
      n_Jets = 0;
      n_Jets_JESUp = 0;
      n_Jets_JESDown = 0;
      nLooseBTaggedJets = 0;
      nMediumBTaggedJets = 0;
      nLooseMuons = 0;
      nTightMuons = 0;
      nLooseElectrons = 0;
      nTightElectrons = 0;
      nTightTaus = 0;
      theMR = -666;
      theMR_JESUp   = -666;
      theMR_JESDown = -666;
      theRsq = -666;
      theRsq_JESUp   = -666;
      theRsq_JESDown = -666;
      t1Rsq  = -666;
      t1Rsq_JESUp   = -666;
      t1Rsq_JESDown = -666;
      
      nSelectedPhotons = 0;
      mGammaGamma    = -1;
      pTGammaGamma   = -1;
      mGammaGammaSC  = -1;
      pTGammaGammaSC = -1;
      mbbZ   = 0;
      mbbH   = 0;
      mbbZ_L = 0;
      mbbH_L = 0;
      run = runNum;
      lumi = lumiNum; 
      event = eventNum;
      passedDiphotonTrigger = false;
      
      //selected photons variables
      for ( int i = 0; i < 2; i++ )
	{
	  Pho_E[i]                  = -99.;
	  Pho_Pt[i]                 = -99.;
	  Pho_Eta[i]                = -99.;
	  Pho_Phi[i]                = -99.;
	  PhoSC_E[i]                = -99.;
	  PhoSC_Pt[i]               = -99.;
	  PhoSC_Eta[i]              = -99.;
	  PhoSC_Phi[i]              = -99.;
	  Pho_SigmaIetaIeta[i]      = -99.;
	  Pho_R9[i]                 = -99.;
	  Pho_HoverE[i]             = -99.;
	  Pho_sumChargedHadronPt[i] = -99.;
	  Pho_sumNeutralHadronEt[i] = -99.;
	  Pho_sumPhotonEt[i]        = -99.;
	  Pho_sigmaEOverE[i]        = -99.;
	  Pho_passEleVeto[i]        = false;
	  Pho_passIso[i]            = false;
	  Pho_motherID[i]           = 0;
	}
      
      //jets
      for ( int i = 0; i < 50; i++ )
	{
	  jet_E[i]   = -99.;
	  jet_Pt[i]  = -99.;
	  jet_Eta[i] = -99.;
	  jet_Phi[i] = -99.;
	  jetIsCSVL[i] = 0;  
	  jetIsCSVM[i] = 0;  
	  jetIsCSVT[i] = 0;  
	}


      
      //------------------
      //Pileup reweighting
      //------------------
      pileupWeight = 1.0;
      if( !isData ) {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
	puhisto->Fill(NPU);
	pileupWeight = helper->getPileupWeight(NPU);
	pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
	pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;	
      }
      
      /////////////////////////////////
      //Scale and PDF variations
      /////////////////////////////////

      if ( (*scaleWeights).size() >= 9 ) 
	{
	  sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
	  sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
	  sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
	  sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
	  sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
	  sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;

	  
	  SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
	  SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
	  SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
	  SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
	  SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
	  SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
	}
      
      sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
      for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) 
	{
	  sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	  SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	}
      

      if ( _debug ) std::cout << "============" << std::endl;
      if ( _debug ) std::cout << "run == " << run << " && evt == " << event << std::endl;
      
      razorbox = LowRes;
      
      //TODO: triggers!
      //bool passedDiphotonTrigger = true;      
      // if (dataset == "76XX") {
      // 	passedDiphotonTrigger = ( HLTDecision[65] );
      // } else if (dataset == "80X") {
      // 	passedDiphotonTrigger = ( HLTDecision[83] );
      // }     
      
      //--------------
      //muon selection
      //--------------
      for( int i = 0; i < nMuons; i++ )
	{
	  if(!isLooseMuon(i)) continue;  
	  if(muonPt[i] < 10) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
	  nLooseMuons++;
	  if( isTightMuon(i) ) nTightMuons++;
	}
      //------------------
      //electron selection
      //------------------
      for( int i = 0; i < nElectrons; i++ )
	{
	  if( !isLooseElectron(i) ) continue; 
	  if( elePt[i] < 10 ) continue;
	  if( abs(eleEta[i]) > 2.5 ) continue;
	  nLooseElectrons++;
      	  if( isTightElectron(i) ) nTightElectrons++;
	}
      //-------------
      //tau selection
      //-------------
      for( int i = 0; i < nTaus; i++ )
	{
	  if( !isTightTau(i) ) continue; 
	  nTightTaus++;
	}
      
      //photon selection
      vector<TLorentzVector> GoodPhotons;
      vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
      vector<bool> GoodPhotonPassesIso; //store whether each photon is isolated
      std::vector< PhotonCandidate > phoCand;//PhotonCandidate defined in RazorAuxPhoton.hh
      int nPhotonsAbove40GeV = 0;
      for(int i = 0; i < nPhotons; i++)
	{

	  double scale = photonCorrector->ScaleCorrection(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]));
	  double smear = photonCorrector->getSmearingSigma(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]), 0., 0.); 

	  //ID cuts -- apply isolation after candidate pair selection
	  if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] << std::endl;
	  if ( !photonPassLooseIDWithoutEleVeto(i) ) {
	    if ( _phodebug ) std::cout << "[DEBUG]: failed run2 ID" << std::endl;
	    continue;
	  }
	
	  //**********************************************************
	  //Isolation, electron veto, and Barrel requirements are introduced here 
	  //if we want to use the "regular" selection sequence
	  //**********************************************************
	  if (!(pho_passEleVeto[i] && photonPassLooseIso(i))) continue;	  
	  if (!(fabs(pho_superClusterEta[i]) < 1.4442 )) continue;

	  //Defining Corrected Photon momentum
	  float pho_pt_corr = phoPt[i];
	  if (doPhotonScaleCorrection) {
	    if (isData) {
	      pho_pt_corr = phoPt[i]*scale; 
	      if (_phodebug) std::cout << "[DEBUG] : Photon Energy Scale Corrections: " << phoPt[i] << " * " << scale << " --> " << pho_pt_corr << "\n";
	    } else {
	      pho_pt_corr = phoPt[i]*(1+smear*random.Gaus());
	    }
	  }
	  TVector3 vec;
	  vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );
	
	  if ( phoPt[i] < 20.0 )
	    {
	      if ( _phodebug ) std::cout << "[DEBUG]: failed pt" << std::endl;
	      continue;
	    }
	
	  if( fabs(pho_superClusterEta[i]) > 2.5 )
	    {
	      //allow photons in the endcap here, 
	      //but if one of the two leading photons is in the endcap,
	      //reject the event
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
	  if( phoPt[i] > 40.0 ) nPhotonsAbove40GeV++;
	  //setting up photon 4-momentum with zero mass
	  TLorentzVector thisPhoton;
	  thisPhoton.SetVectM( vec, .0 );

	  //-----------------------------
	  //uncorrected photon 4-momentum
	  //-----------------------------
	  TVector3 vtx( pvX, pvY, pvZ );
	  TVector3 phoPos;
	  if ( fabs( pho_superClusterEta[i] ) < 1.479 )
	    {
	      phoPos.SetXYZ( EB_R*cos( pho_superClusterPhi[i]), EB_R*sin( pho_superClusterPhi[i] ), EB_R*sinh( pho_superClusterEta[i] ) );
	    }
	  else
	    {
	      double R = fabs( EE_Z/sinh( pho_superClusterEta[i] ) );
	    
	      if ( pho_superClusterEta[i] > .0 )
		{
		  phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), EE_Z);
		}
	      else
		{
		  phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), -EE_Z);
		}
	    
	    }
	
	  TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_RegressionE[i] );
	
	  //Filling Photon Candidate
	  PhotonCandidate tmp_phoCand;
	  tmp_phoCand.Index = i;
	  tmp_phoCand.photon = thisPhoton;
	  tmp_phoCand.photonSC = phoSC;
	  tmp_phoCand.scEta = pho_superClusterEta[i];
	  tmp_phoCand.scPhi = pho_superClusterPhi[i];
	  tmp_phoCand.SigmaIetaIeta = phoFull5x5SigmaIetaIeta[i];
	  tmp_phoCand.R9 = phoR9[i];
	  tmp_phoCand.HoverE = pho_HoverE[i];
	  tmp_phoCand.sumChargedHadronPt = pho_pfIsoChargedHadronIso[i];
	  tmp_phoCand.sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[i];
	  tmp_phoCand.sumPhotonEt = pho_pfIsoPhotonIso[i];
	  tmp_phoCand.sigmaEOverE = pho_RegressionEUncertainty[i]/pho_RegressionE[i];
	  tmp_phoCand._passEleVeto = pho_passEleVeto[i];
	  tmp_phoCand._passIso = photonPassLooseIso(i);
	  phoCand.push_back( tmp_phoCand );
	
	  nSelectedPhotons++;
	}
    
      //------------------------------------------------------------
      //if there is no photon with pT above 40 GeV, reject the event
      //------------------------------------------------------------
      if( nPhotonsAbove40GeV == 0 ) {
	if ( _debug ) std::cout << "[DEBUG]: no photons above 40 GeV, nphotons: " 
				<< phoCand.size() << std::endl;
	continue;
      }

      //--------------------------------------
      //Require at least two photon candidates
      //--------------------------------------
      if ( phoCand.size() < 2 ) {
	if ( _debug ) std::cout << "[INFO]: not enough photon, nphotons: " 
				<< phoCand.size() << std::endl;
	for(int i = 0; i < nPhotons; i++) {
	  if ( _debug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] 
				  << " pho_eta: " << phoEta[i] 
				  << " SIetaIeta: " << phoFull5x5SigmaIetaIeta[i] << std::endl;
	}
	continue;
      }
      
      
      if ( _debug ) std::cout << "[DEBUG]: nphotons--> " << phoCand.size() 
			      << " " << nSelectedPhotons << std::endl;
    
      //----------------------------------------
      //find the "best" photon pair, highest Pt!
      //----------------------------------------
      TLorentzVector HiggsCandidate(0,0,0,0);
      TLorentzVector HiggsCandidateSC(0,0,0,0);
      int HiggsPhoIndex1 = -1;
      int HiggsPhoIndex2 = -1;
      double bestSumPt = -99.;
      std::vector< PhotonCandidate > phoSelectedCand;
      PhotonCandidate bestCand[2];
      for ( size_t i = 0; i < phoCand.size(); i++ )
	{
	  for ( size_t j = i+1; j < phoCand.size(); j++ )
	    {
	      PhotonCandidate pho1 = phoCand[i];
	      PhotonCandidate pho2 = phoCand[j];
	      if ( _debug )
		{
		  std::cout << "[DEBUG]: pho1-> " << pho1.photon.Pt()
			    << "\n[DEBUG]: pho2->" << pho2.photon.Pt() 
			    << std::endl;
		}
	      //need one photon in the pair to have pt > 40 GeV
	      if ( pho1.photon.Pt() < 40.0 && pho2.photon.Pt() < 40.0 )
		{
		  if ( _debug ) std::cout << "[DEBUG]: both photons failed PT > 40 GeV" << std::endl; 
		  //continue;
		}
	      //need diphoton mass between > 100 GeV as in AN (April 1st)
	      double diphotonMass = (pho1.photon + pho2.photon).M();
	      if ( _debug )
		{
		  std::cout << "[DEBUG] Diphoton Sum pT: " << pho1.photon.Pt() + pho2.photon.Pt() << std::endl;
		}
	    
	      if( diphotonMass < 50 )
		{
		  if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 100 GeV: mgg-> " << diphotonMass << std::endl;
		  if ( _debug ) std::cout << "... pho1Pt: " << pho1.photon.Pt()  << " pho2Pt: " << pho2.photon.Pt()  << std::endl;
		  continue;
		}
	      //---------------------------------------------
	      //if the sum of the photon pT's is larger than 
	      //that of the current Higgs candidate, 
	      //make this the Higgs candidate
	      //---------------------------------------------
	      if( pho1.photon.Pt() + pho2.photon.Pt() > bestSumPt )
		{
		  bestSumPt = pho1.photon.Pt() + pho2.photon.Pt();
		  HiggsCandidate = pho1.photon + pho2.photon;
		  HiggsCandidateSC = pho1.photonSC + pho2.photonSC;
		  if ( pho1.photon.Pt() >= pho2.photon.Pt() )
		    {
		      if ( _debug ) std::cout << "assign photon candidate, pho1Pt > pho2Pt" << std::endl;
		      bestCand[0] = pho1;
		      bestCand[1] = pho2;
		      HiggsPhoIndex1 = pho1.Index;
		      HiggsPhoIndex2 = pho2.Index;  
		    }
		  else
		    {
		      if ( _debug ) std::cout << "assign photon candidate, pho2Pt > pho1Pt" << std::endl;
		      bestCand[0] = pho2;
		      bestCand[1] = pho1;
		      HiggsPhoIndex1 = pho2.Index;
		      HiggsPhoIndex2 = pho1.Index;
		    }
		}//best pt if
	    }
	}
    
    
      //---------------------------------------
      //just use this container for convenience
      //to parse the data into TTree
      //---------------------------------------
      phoSelectedCand.push_back(bestCand[0]);
      phoSelectedCand.push_back(bestCand[1]);
    
      //-----------------------------------
      //Filling Selected Photon Information
      //-----------------------------------
      TLorentzVector pho_cand_vec[2];
      int _pho_index = 0;
      for ( auto& tmpPho : phoSelectedCand )
	{
	  if ( !( tmpPho.Index == HiggsPhoIndex1 || tmpPho.Index == HiggsPhoIndex2 ) ) continue;
	  if( _pho_index > 1 ) std::cerr << "[ERROR]: Photon index larger than 1!" << std::endl;
	  pho_cand_vec[_pho_index]           = tmpPho.photon;
	  Pho_E[_pho_index]                  = tmpPho.photon.E();
	  Pho_Pt[_pho_index]                 = tmpPho.photon.Pt();
	  Pho_Eta[_pho_index]                = tmpPho.photon.Eta();
	  Pho_Phi[_pho_index]                = tmpPho.photon.Phi();
	  PhoSC_E[_pho_index]                = tmpPho.photonSC.E();
	  PhoSC_Pt[_pho_index]               = tmpPho.photonSC.Pt();
	  PhoSC_Eta[_pho_index]              = tmpPho.photonSC.Eta();
	  PhoSC_Phi[_pho_index]              = tmpPho.photonSC.Phi();
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
    
      //removing events with less than two good photon candidates
      if ( _pho_index < 2 ) continue;
    
      if ( _debug )
	{
	  std::cout << "[DEBUG]: best photon pair: " 
		    << "\n-> pho1Pt: " << Pho_Pt[0] 
		    << "\n-> pho2Pt: " << Pho_Pt[1] 
		    << std::endl;
	}
    

      //------------------------------------------------------------
      //Define Skim Requirements
      //------------------------------------------------------------
      if (doEleVeto) {
	if( !Pho_passEleVeto[0] || !Pho_passEleVeto[1] ) {
	  if ( _debug ) std::cout << "[DEBUG]: Failed EleVeto: pho1, pho2: " << Pho_passEleVeto[0] << ", " << Pho_passEleVeto[1] << std::endl;
	  if ( _debug ) std::cout << "[DEBUG]: pho1Pt: " << Pho_Pt[0] << " pho2Pt: " << Pho_Pt[1] << std::endl;
	  for ( auto& phoC : phoSelectedCand ) {
	    if ( _debug ) std::cout << "===> phopt: " << phoC.photon.Pt() << " phoEta: " << phoC.photon.Eta() << std::endl;
	  }
	  continue;
	}
      }

      //if the best candidate pair has a non-isolated photon, reject the event
      if (doRequireIso) {
	if( !Pho_passIso[0] || !Pho_passIso[1] ) {
	  if ( _debug ) std::cout << "[DEBUG]: Failed ISO: pho1, pho2: " << Pho_passIso[0] << ", " << Pho_passIso[1] << std::endl;
	  if ( _debug ) std::cout << "[DEBUG]: pho1Pt: " << Pho_Pt[0] << " pho2Pt: " << Pho_Pt[1] << std::endl;
	  for ( auto& phoC : phoSelectedCand ) {
	    if ( _debug ) std::cout << "===> phopt: " << phoC.photon.Pt() << " phoEta: " << phoC.photon.Eta() << std::endl;
	  }
	  continue;
	}
      }

      //record higgs candidate info
      mGammaGamma    = HiggsCandidate.M();
      pTGammaGamma   = HiggsCandidate.Pt();
      mGammaGammaSC  = HiggsCandidateSC.M();
      pTGammaGammaSC = HiggsCandidateSC.Pt();
      if ( _debug ) std::cout << "[DEBUG]: mgg-> " << mGammaGamma << " pTgg->" << pTGammaGamma << std::endl;
    
      //******************************************************
      //compute trigger efficiency weights for MC
      //******************************************************
      triggerEffWeight = 1.0;
      double leadPhoPt=0;
      double leadPhoEta=0;
      double trailingPhoPt=0;
      double trailingPhoEta=0;
      if (Pho_Pt[0] > Pho_Pt[1]) {
	leadPhoPt = Pho_Pt[0];
	leadPhoEta = Pho_Eta[0];
	trailingPhoPt = Pho_Pt[1];
	trailingPhoEta= Pho_Eta[1];
      } else {
	leadPhoPt = Pho_Pt[1];
	leadPhoEta = Pho_Eta[1];
	trailingPhoPt = Pho_Pt[0];
	trailingPhoEta= Pho_Eta[0];
      }
      double triggerEffLeadingLeg = 
	photonTriggerEffHist_LeadingLeg->GetBinContent( photonTriggerEffHist_LeadingLeg->GetXaxis()->FindFixBin( fabs(leadPhoEta) ),
							photonTriggerEffHist_LeadingLeg->GetYaxis()->FindFixBin( fmax( fmin(leadPhoPt,99.9), 20.01 ) )
							);
      double triggerEffTrailingLeg = 
	photonTriggerEffHist_TrailingLeg->GetBinContent( photonTriggerEffHist_TrailingLeg->GetXaxis()->FindFixBin( fabs(trailingPhoEta) ),
							 photonTriggerEffHist_TrailingLeg->GetYaxis()->FindFixBin( fmax( fmin(trailingPhoPt,99.9), 20.01 ) )
							 );
      triggerEffWeight = triggerEffLeadingLeg*triggerEffTrailingLeg;

      //***********************************************************
      //get mother ID of photons
      //***********************************************************
      //cout << "Photon1 : " << Pho_Pt[0] << " " << Pho_Eta[0] << " " << Pho_Phi[0] << "\n";
      for(int g = 0; g < nGenParticle; g++){
	if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[0],Pho_Phi[0]) < 0.5) ) continue;
	// status = 22 for Higgs bosons in MadGraph/Pythia8
	//if(gParticleStatus[g] != 1) continue;
	if(gParticleId[g] != 22) continue;
	if(!( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23)) continue;
	Pho_motherID[0] = gParticleMotherId[g];
	//cout << "Nearby GenParticle: " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherId[g] << "\n";
      }

      //cout << "Photon2 : " << Pho_Pt[1] << " " << Pho_Eta[1] << " " << Pho_Phi[1] << "\n";
      for(int g = 0; g < nGenParticle; g++){
	if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[1],Pho_Phi[1]) < 0.5) ) continue;
	// status = 22 for Higgs bosons in MadGraph/Pythia8
	//if(gParticleStatus[g] != 1) continue;
	if(gParticleId[g] != 22) continue;
	if(!( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23)) continue;
	Pho_motherID[1] = gParticleMotherId[g];      
	//cout << "Nearby GenParticle: " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherId[g] << "\n";
      }

      //cout << "\nGenParticles:\n";
       // for(int g = 0; g < nGenParticle; g++){
       //   cout << "GenParticle: " << gParticleId[g] << " " << gParticleStatus[g] << " : " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherId[g] << "\n";
       // }
       // cout << "\n\n";



      //----
      //Jets
      //----
      //Propagate jet uncertainties to MET
      float MetXCorr_JESUp = 0;
      float MetYCorr_JESUp = 0;
      float MetXCorr_JESDown = 0;
      float MetYCorr_JESDown = 0;
    
      vector<TLorentzVector> GoodJets;
      vector<bool> GoodJetsIsCVSL;
      vector<bool> GoodJetsIsCVSM;
      vector<bool> GoodJetsIsCVST;
      vector<TLorentzVector> GoodJetsJESUp;
      vector<TLorentzVector> GoodJetsJESDown;
      vector< pair<TLorentzVector, bool> > GoodCSVLJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not
      for(int i = 0; i < nJets; i++)
	{
	  //Jet Corrections                                                                      
	  double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						  fixedGridRhoAll, jetJetArea[i],
						  JetCorrector );
      
	  TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
	
	  if( thisJet.Pt() < JET_CUT ) continue;//According to the April 1st 2015 AN
	  if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
	  //int level = 2; //loose jet ID
	  if ( !jetPassIDLoose[i] ) continue;
	  //if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;
	
	  //exclude selected photons from the jet collection
	  double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
	  if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN
      
	  GoodJets.push_back(thisJet);
	  GoodJetsIsCVSL.push_back(isCSVL(i));
	  GoodJetsIsCVSM.push_back(isCSVM(i));
	  GoodJetsIsCVST.push_back(isCSVT(i));
	  n_Jets++;
	
	  double jetCorrPt = thisJet.Pt();
	  double jetCorrE  = thisJet.E();
	  if ( !isData )
	    {
	      //****************************************************************************
	      //Apply b-tagging correction factor 
	      //****************************************************************************
	      if ( !isData && abs(jetEta[i]) < 2.4 && jetCorrPt > JET_CUT ) 
		{ 
		  double effMedium = 0;
		  double effLoose  = 0;
		  BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
		  if ( abs(jetPartonFlavor[i]) == 5) 
		    {
		      effMedium = btagMediumEfficiencyHist->GetBinContent( btagMediumEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									   btagMediumEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      effLoose  = btagLooseEfficiencyHist->GetBinContent( btagLooseEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									  btagLooseEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      jetType = BTagEntry::FLAV_B;
		    } 
		  else if ( abs(jetPartonFlavor[i]) == 4) 
		    {
		      effMedium = btagMediumCharmEfficiencyHist->GetBinContent( btagMediumCharmEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
										btagMediumCharmEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      effLoose  = btagLooseCharmEfficiencyHist->GetBinContent( btagLooseCharmEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									       btagLooseCharmEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      jetType = BTagEntry::FLAV_C;
		    } 
		  else 
		    {
		      effMedium = btagMediumLightJetsEfficiencyHist->GetBinContent(
										   btagMediumLightJetsEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0))
										   ,btagMediumLightJetsEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      effLoose = btagLooseLightJetsEfficiencyHist->GetBinContent(btagLooseLightJetsEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0))
										 ,btagLooseLightJetsEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		      jetType = BTagEntry::FLAV_UDSG;
		    }

		  //----------------
		  //get scale factor
		  //----------------
		  double jetSF_Loo = -1;
		  double jetSF_LooUp = -1;
		  double jetSF_LooDown = -1;
		  double jetSF_Med = -1;
		  double jetSF_MedUp = -1;
		  double jetSF_MedDown = -1;  
		  if ( abs(jetPartonFlavor[i]) == 5 || abs(jetPartonFlavor[i]) == 4 )//c,b quarks
		    {
		      if (jetCorrPt < 670.) //670 is the largest pt range listed in the CSV text file
			{
			  //M
			  jetSF_Med     = btagreaderM.eval(jetType, jetEta[i], jetCorrPt); 
			  jetSF_MedUp   = btagreaderM_up.eval(jetType, jetEta[i], jetCorrPt);
			  jetSF_MedDown = btagreaderM_do.eval(jetType, jetEta[i], jetCorrPt);
			  //L
			  jetSF_Loo     = btagreaderL.eval(jetType, jetEta[i], jetCorrPt);
			  jetSF_LooUp   = btagreaderL_up.eval(jetType, jetEta[i], jetCorrPt);
			  jetSF_LooDown = btagreaderL_do.eval(jetType, jetEta[i], jetCorrPt);
			}
		      else 
			{
			  //M
			  jetSF_Med     = btagreaderM.eval(jetType, jetEta[i], 669);
			  jetSF_MedUp   = btagreaderM_up.eval(jetType, jetEta[i], 669);
			  jetSF_MedDown = btagreaderM_do.eval(jetType, jetEta[i], 669);
			  //L
			  jetSF_Loo     = btagreaderL.eval(jetType, jetEta[i], 669);
			  jetSF_LooUp   = btagreaderL_up.eval(jetType, jetEta[i], 669);
			  jetSF_LooDown = btagreaderL_do.eval(jetType, jetEta[i], 669);
			}
		    } 
		  else//rest of the quarks
		    {
		      //M
		      jetSF_Med     = btagreaderMistagM.eval(jetType, jetEta[i], 100);//fix in eta and pt
		      jetSF_MedUp   = btagreaderMistagM_up.eval(jetType, jetEta[i], 100);  
		      jetSF_MedDown = btagreaderMistagM_do.eval(jetType, jetEta[i], 100);
		      //L (to be checked)
		      if ( jetCorrPt < 1000 ) 
			{
			  jetSF_Loo     = btagreaderMistagL.eval(jetType, jetEta[i], jetCorrPt);
			  jetSF_LooUp   = btagreaderMistagL_up.eval(jetType, jetEta[i], jetCorrPt);
			  jetSF_LooDown = btagreaderMistagL_do.eval(jetType, jetEta[i], jetCorrPt);
			}
		      else
			{
			  jetSF_Loo     = btagreaderMistagL.eval(jetType, jetEta[i], 999);
			  jetSF_LooUp   = btagreaderMistagL_up.eval(jetType, jetEta[i], 999);
			  jetSF_LooDown = btagreaderMistagL_do.eval(jetType, jetEta[i], 999);
			}
		    }
		
		  //------------------------
		  //Loose Working point only
		  //------------------------
		  //Apply Scale Factors
		  if ( jetSF_Med <= 0 || jetSF_MedUp <= 0 || jetSF_MedDown <= 0  || jetSF_Loo <= 0 || jetSF_LooUp <= 0 || jetSF_LooDown <= 0 )
		    {
		      std::cout << "Warning: b-tag scale factor is <= 0!" << std::endl;
		      std::cout << jetSF_Med << " " << jetSF_MedUp << " " << jetSF_MedDown << " " << jetSF_Loo << " " << jetSF_LooUp << " " 
				<< jetSF_LooDown << std::endl;
		    }
		  else if ( isCSVL(i) )
		    {
		      btagCorrFactor *= jetSF_Loo;
		      if ( abs( jetPartonFlavor[i] ) == 5 || abs( jetPartonFlavor[i] ) == 4 )
			{
			  sf_btagUp   *= jetSF_LooUp/jetSF_Loo;
			  sf_btagDown *= jetSF_LooDown/jetSF_Loo;
			}
		      else
			{
			  sf_bmistagUp   *= jetSF_LooUp/jetSF_Loo;
			  sf_bmistagDown *= jetSF_LooDown/jetSF_Loo;                                                                                          
			} 
		    }
		  else
		    {
		      //only apply the scale factor on the inefficiency, if the corrected efficiency doesn't go above 100%
		      //only record up/down systematics if the nominal and up and down corrected systematics do not go above 100%
		      double sf = 1.0;
		      if (effLoose*jetSF_Loo < 1.0) sf = (1/effLoose - jetSF_Loo) / (1/effLoose - 1);
		    
		      btagCorrFactor *= sf;
		      if (abs(jetPartonFlavor[i]) == 5 || abs(jetPartonFlavor[i]) == 4) 
			{
			  if (effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooUp < 1.0) 
			    {
			      sf_btagUp *= (1/effLoose - jetSF_LooUp) / (1/effLoose - 1) / sf;
			    }
			  if (effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooDown < 1.0) 
			    {
			      sf_btagDown *= (1/effLoose - jetSF_LooDown) / (1/effLoose - 1) / sf;
			    }
			} 
		      else 
			{
			  if ( effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooUp < 1.0 ) 
			    {
			      sf_bmistagUp *= (1/effLoose - jetSF_LooUp) / (1/effLoose - 1) / sf;
			    } 
			  if ( effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooDown < 1.0)
			    {
			      sf_bmistagDown *= (1/effLoose - jetSF_LooDown) / (1/effLoose - 1) / sf;
			    }
			}
		    
		    }
		
		}//Jetcut
	    }//isData
	
	  if ( !isData )
	    {
	      double unc = helper->getJecUnc( jetCorrPt, jetEta[i] );
	      double jetPtJESUp = jetCorrPt*(1+unc);
	      double jetPtJESDown = jetCorrPt/(1+unc);
	      double jetEJESUp = jetCorrE*(1+unc);
	      double jetEJESDown = jetCorrE/(1+unc);
	      TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
	      TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
	    
	      //Propagate uncertainties to the MET
	      if (jetPtJESUp > 20) 
		{
		  MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
		  MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
		}
	      if (jetPtJESDown > 20) 
		{
		  MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
		  MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
		}
	    
	    
	      if ( jetPtJESUp > JET_CUT )
		{
		  GoodJetsJESUp.push_back(thisJetJESUp);
		  n_Jets_JESUp++;
		}
	      if ( jetPtJESDown > JET_CUT )
		{
		  GoodJetsJESDown.push_back(thisJetJESDown);
		  n_Jets_JESDown++;
		}
	    }
	
	  /*
	    Change to isCSVL and isCSVM if you want CISV
	  */
	  if( isCSVL(i) )
	    {
	      nLooseBTaggedJets++;
	      if( isCSVM(i) )
		{ 
		  nMediumBTaggedJets++;
		  GoodCSVLJets.push_back(make_pair(thisJet, true));
		}
	      else
		{
		  GoodCSVLJets.push_back(make_pair(thisJet, false));
		}
	    }
	}
    
      //if there are no good jets, reject the event
      if( n_Jets == 0 )
	{
	  if ( _debug ) std::cout << "[DEBUG]: No Jets Selected" << std::endl;
	  //continue;
	}

      //std::cout << "njets-->" << n_Jets << std::endl;
    
      for ( int iJet = 0; iJet < int(GoodJets.size()) ; iJet++ ) {
	jet_E[iJet] = GoodJets[iJet].E();
	jet_Pt[iJet] = GoodJets[iJet].Pt();
	jet_Eta[iJet] = GoodJets[iJet].Eta();
	jet_Phi[iJet] = GoodJets[iJet].Phi();
	jetIsCSVL[iJet] = GoodJetsIsCVSL[iJet];
	jetIsCSVM[iJet] = GoodJetsIsCVSM[iJet];
	jetIsCSVT[iJet] = GoodJetsIsCVST[iJet];
	iJet++;
      }
    
      //Compute the razor variables using the selected jets and the diphoton system
      HT = Pho_Pt[0] + Pho_Pt[1]; //HT = sum of photon pT  + jet pT
      vector<TLorentzVector> JetsPlusHiggsCandidate;
      for( auto& jet : GoodJets ) {
	JetsPlusHiggsCandidate.push_back(jet);
	HT += jet.Pt();
      }
      JetsPlusHiggsCandidate.push_back(HiggsCandidate);
    
      TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
      MET = metPt;
      t1MET = metType1Pt;
    
      vector<TLorentzVector> hemispheres = getHemispheres(JetsPlusHiggsCandidate);
      theMR  = computeMR(hemispheres[0], hemispheres[1]); 
      if ( theMR > 0 )
	{
	  theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  t1Rsq  = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
	}
    
      //***********************
      //MR skim
      //***********************
      if (doMRSkim) {
	if (!(theMR > 150 && GoodJets.size() >= 1)) continue;
      }


      //Calculations for JES systematics
      if( !isData )
	{
	  //JES up
	  vector<TLorentzVector> JetsPlusHiggsCandidate_JESUp;
	  for( auto& jet : GoodJetsJESUp ) JetsPlusHiggsCandidate_JESUp.push_back(jet);
	  JetsPlusHiggsCandidate_JESUp.push_back(HiggsCandidate);

	  float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
	  float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
	  float t1PFMetXJESUp = t1PFMET.Px() + MetXCorr_JESUp;
	  float t1PFMetYJESUp = t1PFMET.Py() + MetYCorr_JESUp;
	
	  TLorentzVector PFMET_JESUp(PFMetXJESUp, PFMetYJESUp, 0, sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) )); 
	  TLorentzVector t1PFMET_JESUp(t1PFMetXJESUp, t1PFMetYJESUp, 0, sqrt( pow(t1PFMetXJESUp,2) + pow(t1PFMetYJESUp,2) ));
	  vector<TLorentzVector> hemispheres_JESUp = getHemispheres(JetsPlusHiggsCandidate_JESUp);
	  theMR_JESUp  = computeMR(hemispheres_JESUp[0], hemispheres_JESUp[1]); 
	  theRsq_JESUp = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], PFMET_JESUp);
	  t1Rsq_JESUp  = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], t1PFMET_JESUp);
	  MET_JESUp    = PFMET_JESUp.Pt();
	  t1MET_JESUp  = t1PFMET_JESUp.Pt();

	  //JES down
	  vector<TLorentzVector> JetsPlusHiggsCandidate_JESDown;
	  for( auto& jet : GoodJetsJESDown ) JetsPlusHiggsCandidate_JESDown.push_back(jet);
	  JetsPlusHiggsCandidate_JESDown.push_back(HiggsCandidate);
	
	  float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
	  float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
	  float t1PFMetXJESDown = t1PFMET.Px() + MetXCorr_JESDown;
	  float t1PFMetYJESDown = t1PFMET.Py() + MetYCorr_JESDown;
	
	  TLorentzVector PFMET_JESDown(PFMetXJESDown, PFMetYJESDown, 0, sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) )); 
	  TLorentzVector t1PFMET_JESDown(t1PFMetXJESDown, t1PFMetYJESDown, 0, sqrt( pow(t1PFMetXJESDown,2) + pow(t1PFMetYJESDown,2) ));
	  vector<TLorentzVector> hemispheres_JESDown = getHemispheres(JetsPlusHiggsCandidate_JESDown);
	  theMR_JESDown  = computeMR(hemispheres_JESDown[0], hemispheres_JESDown[1]); 
	  theRsq_JESDown = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], PFMET_JESDown);
	  t1Rsq_JESDown  = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], t1PFMET_JESDown);
	  MET_JESDown    = PFMET_JESDown.Pt();
	  t1MET_JESDown  = t1PFMET_JESDown.Pt();
	}
    
      if ( theMR < 0.0 )
	{
	  if ( _debug ) std::cout << "[INFO]: MR < 150 GeV, MR: " << theMR << std::endl;
	  for ( auto& jet : JetsPlusHiggsCandidate )
	    {
	      if ( _debug ) std::cout << "phoPT: " << pTGammaGamma 
				      << " jet pt : " << jet.Pt() << " eta: " << jet.Eta() << " phi: " << jet.Phi() 
				      << " h1 pt: " << hemispheres[0].Pt() << " h1 eta: " << hemispheres[0].Eta()
				      << " h2 pt: " << hemispheres[1].Pt() << " h2 eta: " << hemispheres[1].Eta() << std::endl;
	    }
	  //continue;
	}
    
      //if there are two loose b-tags and one medium b-tag, look for b-bbar resonances
      if( nLooseBTaggedJets > 1 && nMediumBTaggedJets > 0 )
	{
	  for(int i = 0; i < nLooseBTaggedJets; i++)
	    {
	      for(int j = i+1; j < nLooseBTaggedJets; j++)
		{
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
    
      if( nLooseBTaggedJets >= 2 )//at least two btag jets
	{
	  for(int i = 0; i < nLooseBTaggedJets; i++)
	    {
	      for(int j = i+1; j < nLooseBTaggedJets; j++)
		{
		  double mbb = (GoodCSVLJets[i].first + GoodCSVLJets[j].first).M();
		  //if mbb is closer to the higgs mass than mbbH, make mbbH = mbb
		  if( fabs(mbbH_L - 125.0) > fabs(mbb - 125.0) ) mbbH_L = mbb;
		  //same for mbbZ
		  if( fabs(mbbZ_L - 91.2) > fabs(mbb - 91.2) ) mbbZ_L = mbb;
		}//end second jet loop
	    }//end first jet loop
	}


      //------------------------------------------------
      //I n v a ri a n t   m a s s   r e s o l u t i o n
      //------------------------------------------------
      sigmaMoverM = 0.5*sqrt( Pho_sigmaEOverE[0]*Pho_sigmaEOverE[0] + Pho_sigmaEOverE[1]*Pho_sigmaEOverE[1] );

      //Writing output to tree
      //HighPt Box
      if ( pTGammaGamma > 110.0 )
	{
	  razorbox = HighPt;
	  razorTree->Fill();
	}
      //Hbb Box
      else if ( mbbH > 110.0 && mbbH < 140.0 )
	{
	  razorbox = Hbb;
	  razorTree->Fill();
	}
      //Zbb Box
      else if( mbbZ > 76.0 && mbbZ < 106.0 )
	{
	  razorbox = Zbb;
	  razorTree->Fill();
	}
      //HighRes Box
      else if( Pho_sigmaEOverE[0] < 0.015 && Pho_sigmaEOverE[1] < 0.015 )
	{
	  razorbox = HighRes;
	  razorTree->Fill();
	}
      //LowRes Box
      else
	{
	  razorbox = LowRes;
	  razorTree->Fill();
	}
      //end of event loop
    }
  
  if ( _info ) std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;
  if ( _info ) std::cout << "[INFO]: Writing output trees..." << std::endl;
  
  outFile->cd();
  razorTree->Write();
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  puhisto->Write();
  outFile->Close();


  delete photonCorrector;
  delete helper;

}
