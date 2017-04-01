//LOCAL INCLUDES
#include "HggRazorUpgradeTiming.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <sys/stat.h>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <assert.h>
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"
#include "TMVA/Reader.h"

using namespace std;

enum HggRazorBox {
    HighPt = 0,
    Hbb = 1,
    Zbb = 2,
    HighRes = 3,
    LowRes = 4
};

struct NewPhotonCandidate
{                                                  
  int   Index;
  TLorentzVector photon;
  TLorentzVector photonSC;
  float scE;
  float scPt;
  float scEta;
  float scPhi;
  float scSeedX;
  float scSeedY;
  float scSeedZ;
  float scSeedT;
  float scX;
  float scY;
  float scZ;
  float SigmaIetaIeta;                                                                        
  float R9;                                                                                  
  float HoverE;                                                                        
  float sumChargedHadronPt;                                                                
  float sumNeutralHadronEt;                                                     
  float sumPhotonEt;                                            
  float sigmaEOverE;
  bool  _passEleVeto;
  bool  _passIso;
  int convType;
  float convTrkZ;
  float convTrkClusZ;
  float vtxSumPx[500];
  float vtxSumPy[500];
};

struct evt
{
  std::string run;
  std::string event;
};

//#define _phodebug 0
//#define _debug    0
#define _info     0
#define FeedFirstN_only 0

const double EB_R = 129.0;
const double EE_Z = 317.0;


const double JET_CUT = 30.;
const int NUM_PDF_WEIGHTS = 60;

 float SIGMATRKZ[12] = { 0.0125255,
                              0.0178107,
                              0.0581667,
                              0.152157,
                              0.716301,
                              1.3188,
                              0.38521,
                              0.702755,
                              3.17615,
                              2.23662,
                              1.67937,
                              2.46599
                            };
                            
const float SIGMATRKCLUSZ[12] = { 0.0298574,
                                  0.0935307,
                                  0.180419,
                                  0.577081,
                                  0.414393,
                                  0.756568,
                                  0.494722,
                                  0.892751,
                                  1.06805,
                                  0.62143,
                                  1.21941,
                                  1.56638
                                };

//Testing branching and merging
void HggRazorUpgradeTiming::Analyze(bool isData, bool useTiming, bool usePhoChi2, bool useOddEvent, int option, string outFileName, string label)
{
  std::cout << "in analyzer" << std::endl;
  gROOT->Reset();
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  TRandom3 randomPhotonTime(1111);
  bool combineTrees = true;
  bool doPhotonScaleCorrection = true;

  bool doEleVeto = true;
  if (option == 1 || option == 11) doEleVeto = false;

  std::cout << "[INFO]: option = " << option << std::endl;
  std::cout << "[INFO]: doEleVeto --> " << doEleVeto << std::endl;
  std::cout << "[INFO]: doPhotonScaleCorrection --> " << doPhotonScaleCorrection << std::endl;
 
    //initialization: create one TTree for each analysis box 
  if ( _info) std::cout << "Initializing..." << std::endl;
  cout << "Combine Trees = " << combineTrees << "\n";
  if (outFileName.empty()){
    if ( _info ) std::cout << "HggRazorUpgradeTiming: Output filename not specified!" << endl << "Using default output name HggRazorUpgradeTiming.root" << std::endl;
    outFileName = "HggRazorUpgradeTiming.root";
  }

  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //one tree to hold all events                                                                                           
  TTree *razorTree = new TTree("HggRazor", "Info on selected razor inclusive events");
  
  //Get CMSSW path directory
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");

  //Photon Energy Scale and Resolution Corrections
   std::string photonCorrectionPath = "";
  if ( cmsswPath != NULL ) photonCorrectionPath = string(cmsswPath) + "/src/RazorAnalyzer/data/PhotonCorrections/";
 
  EnergyScaleCorrection_class *photonCorrector = 0;
  if (option == 0 || option == 1) {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/76X_16DecRereco_2015", photonCorrectionPath.c_str()));
  } else if (option == 10 || option == 11) {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/80X_2016", photonCorrectionPath.c_str()));
  }

  if(!isData) {
    photonCorrector->doScale = false; 
    photonCorrector->doSmearings = true;
  } else {
    photonCorrector->doScale = true; 
    photonCorrector->doSmearings = false;
  }
  
  //MVA Vertex Selection
  bool doMVAVertex = true;
  float ptasym = 0.;
  float chi2_pho_vtx_train = 0.;
  float logchi2_pho_vtx_train = 0.;
  float ptbal = 0.;
  float logsumpt2 = 0.;
  Float_t vertexsumpt = 0.;
  float vertexpt = 0.;
  float diphotonpt = 0.;
  float virtualVtxdZ = 0.;
  float pull_conv = 0.;
  float nConv = 0.;

  std::vector<float>   *allTrackdT_Sig;
  std::vector<float>   *allTrackPt_Sig;
  std::vector<float>   *allTrackdT_Bkg;
  std::vector<float>   *allTrackPt_Bkg;
  std::vector<float>   *allTrackdZ_Sig;
  std::vector<float>   *allTrackdZ_Bkg;

  allTrackdT_Sig = 0;
  allTrackPt_Sig = 0;
  allTrackdT_Bkg = 0;
  allTrackPt_Bkg = 0;
  allTrackdZ_Sig = 0;
  allTrackdZ_Bkg = 0;

 
  Int_t isMatchPv[500]; 
  Int_t isMaxbdtPv[500]; 
  Int_t isMaxlogsumpt2Pv[500]; 
  float pho1_vtx_dt_gen=0.0;// = {0.};
  float pho2_vtx_dt_gen=0.0;// = {0.};
  float pho1_vtx_dt_all[500];// = {0.};
  float pho2_vtx_dt_all[500];// = {0.};
  float chi2_pho_vtx_gen = 0.0;// = {0.};
  int   rank_chi2_genVtx = 0;// = {0.};
  int   rank_chi2_min_genVtx = 0;// = {0.};
  float chi2_min_pho_vtx_gen = 0.0;// = {0.};
  float chi2_pho_vtx_all[500];// = {0.};
  Int_t rank_chi2_pho_vtx_all[500];// = {0.};
  Int_t rank_chi2_min_pho_vtx_all[500];// = {0.};
  float chi2_min_pho_vtx_all[500];// = {0.};
  float min_vtxT_all[500];// = {0.};
  float min_vtxT_gen=0.0;// = {0.};
  float ptasym_all[500];// = {0.};
  float ptbal_all[500];// = {0.};
  float logsumpt2_all[500];// = {0.};
  Float_t vertexsumpt_all[500];// = {0.};
  float vertexpt_all[500];// = {0.};
  float pull_conv_all[500];// = {0.};
  float nConv_all[500];// = {0.}; 
  float vtxdZCoordinate_all[500];// = {0.}; 
  float vtxdXCoordinate_all[500];// = {0.}; 
  float vtxdYCoordinate_all[500];// = {0.}; 
  float vtxdTCoordinate_all[500];// = {0.}; 
  float vtxXCoordinate_all[500];// = {0.}; 
  float vtxYCoordinate_all[500];// = {0.}; 
  float vtxZCoordinate_all[500];// = {0.}; 
  float vtxTCoordinate_all[500];// = {0.}; 
  Int_t pvNtrack_all[500];// = {0.};
  Int_t hasLargeTrack[500];// = {0.};
  Int_t hasLargeTrack_dt[500];// = {0.};

  Int_t NPV_HasLargeTrack = 0;
  Int_t NPV_HasLargeTrack_dt = 0;
  Int_t NPV_GEN_HasLargeTrack = 0;
  Int_t NPV_GEN_HasLargeTrack_dt = 0;
  for(int i=0;i<500;i++)
  {
	isMatchPv[i]=0;
	isMaxbdtPv[i]=0;
	isMaxlogsumpt2Pv[i]=0;
	chi2_pho_vtx_all[i]=0.0;
	rank_chi2_pho_vtx_all[i]=1;
	rank_chi2_min_pho_vtx_all[i]=1;
	pho1_vtx_dt_all[i]=0.0;
	pho2_vtx_dt_all[i]=0.0;
	chi2_min_pho_vtx_all[i]=0.0;
	min_vtxT_all[i]=0.0;
	ptasym_all[i]=0.0;
	pvNtrack_all[i]=0;
	hasLargeTrack[i]=0;
	hasLargeTrack_dt[i]=0;
	ptbal_all[i]=0.0;
     	logsumpt2_all[i]=0.0;	
     	vertexsumpt_all[i]=0.0;	
     	vertexpt_all[i]=0.0;	
  	pull_conv_all[i]=0.0;
	nConv_all[i]=0.0;
	vtxdZCoordinate_all[i]=0.0;
	vtxdXCoordinate_all[i]=0.0;
	vtxdYCoordinate_all[i]=0.0;
	vtxdTCoordinate_all[i]=0.0;
	vtxXCoordinate_all[i]=0.0;
	vtxYCoordinate_all[i]=0.0;
	vtxZCoordinate_all[i]=0.0;
	vtxTCoordinate_all[i]=0.0;
  }
  
  TMVA::Reader *vtxmvareader = 0;
  TMVA::Reader *vtxmvareader_chi2 = 0;
  if (doMVAVertex) {
    vtxmvareader = new TMVA::Reader( "!Color:Silent" );
    vtxmvareader->AddVariable("ptasym", &ptasym );
    vtxmvareader->AddVariable("ptbal", &ptbal );
    vtxmvareader->AddVariable("logsumpt2", &logsumpt2 );
    vtxmvareader->AddVariable("limPullToConv", &pull_conv );
    vtxmvareader->AddVariable("nConv", &nConv );
    
    std::string vtxpathname;
    if ( cmsswPath != NULL ) vtxpathname = string(cmsswPath) + "/src/RazorAnalyzer/data/";
    vtxmvareader->BookMVA("BDT",Form("%s/TMVAClassification_BDTVtxId_SL_2016.xml",vtxpathname.c_str()));
 
    vtxmvareader_chi2 = new TMVA::Reader( "!Color:Silent" );
    vtxmvareader_chi2->AddVariable("logchi2_pho_vtx", &logchi2_pho_vtx_train );
    vtxmvareader_chi2->AddVariable("ptasym", &ptasym );
    vtxmvareader_chi2->AddVariable("ptbal", &ptbal );
    vtxmvareader_chi2->AddVariable("logsumpt2", &logsumpt2 );
    vtxmvareader_chi2->AddVariable("limPullToConv", &pull_conv );
    vtxmvareader_chi2->AddVariable("nConv", &nConv );
    if ( cmsswPath != NULL ) vtxpathname = string(cmsswPath) + "/src/RazorAnalyzer/data/";
    vtxmvareader_chi2->BookMVA("BDT",Form("%s/TMVAClassification_BDTVtxId_SL_2016_chi2.xml",vtxpathname.c_str()));

  }

  //Including Jet Corrections
  std::vector<JetCorrectorParameters> correctionParameters;
  //get correct directory for JEC files (different for lxplus and t3-higgs)
  std::string pathname;
  if ( cmsswPath != NULL ) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  std::cout << "Getting JEC parameters from " << pathname << std::endl;
  if ( isData ) 
    {
      std::cout << "[INFO]: getting data JEC" << std::endl;
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
    } 
  else 
    {
      std::cout << "[INFO]: getting MC JEC" << std::endl;
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    }
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector( correctionParameters );
  //------------------------------------------------------------
  //Get JEC uncertainty file and set up JetCorrectionUncertainty
  //------------------------------------------------------------
  string jecUncPath;
  if ( isData ) 
    {
      jecUncPath = pathname+"/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
    }
  else 
    {
      jecUncPath = pathname+"/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
    }
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncPath);
  
  //--------------
  //Pileup Weights
  //--------------
  
  TFile *pileupWeightFile = 0;
  TH1F *pileupWeightHist = 0;
  TH1F *pileupWeightSysUpHist = 0;
  TH1F *pileupWeightSysDownHist = 0;
  if( !isData )
    {
      string pathname;
      if (cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/";
      else 
	{
	  std::cout << "ERROR: CMSSW_BASE not detected. Exiting..." << std::endl;
	  assert(false);
	}
      
      pileupWeightFile = TFile::Open(Form("%s/PileupReweight2015_7_6.root",pathname.c_str()));
      pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
      pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
      pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
      assert(pileupWeightHist);
      assert(pileupWeightSysUpHist);
      assert(pileupWeightSysDownHist);
    }
  
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
      TFile *btagLooseEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVL_Fullsim80X_ICHEP2016.root");
      btagLooseEfficiencyHist = (TH2D*)btagLooseEfficiencyFile->Get("Efficiency_PtEta");
      assert(btagLooseEfficiencyHist);
      TFile *btagLooseCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVL_Fullsim80X_ICHEP2016.root");
      btagLooseCharmEfficiencyHist = (TH2D*)btagLooseCharmEfficiencyFile->Get("Efficiency_PtEta");
      assert(btagLooseCharmEfficiencyHist);
      TFile *btagLooseLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVL_Fullsim80X_ICHEP2016.root");
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
  for(size_t i = 0; i < boxNames.size(); i++){
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
  
  int nSelectedPhotons;
  float mGammaGamma, pTGammaGamma, etaGammaGamma, phiGammaGamma, mGammaGammaSC, pTGammaGammaSC, etaGammaGammaSC, phiGammaGammaSC, sigmaMoverM;
  float mjj, mggj[100], mggjj, mbbZ, mbbZ_L, mbbH, mbbH_L;
  float cosThetaCS;
  bool passedDiphotonTrigger;
  HggRazorBox razorbox = LowRes;
  
  unsigned int run, lumi, event;
  
  //selected photon variables
  float Pho_E[2], Pho_Pt[2], Pho_Eta[2], Pho_Phi[2], Pho_SigmaIetaIeta[2], Pho_R9[2], Pho_HoverE[2], Pho_Time[2];
  float PhoSC_E[2], PhoSC_Pt[2], PhoSC_Eta[2], PhoSC_Phi[2];
  float PhoDefaultSC_E[2], PhoDefaultSC_Pt[2], PhoDefaultSC_Eta[2], PhoDefaultSC_Phi[2];
  float Pho_sumChargedHadronPt[2], Pho_sumNeutralHadronEt[2], Pho_sumPhotonEt[2], Pho_sigmaEOverE[2];
  bool  Pho_passEleVeto[2], Pho_passIso[2];
  int   Pho_motherID[2];
  int   vtxIndex;
  float vtxZCoordinate;
  float vtxdZCoordinate;
  float vtxTCoordinate;
  float vtxdTCoordinate;

  //jet information
  int n_Jets, nLooseBTaggedJets, nMediumBTaggedJets;
  float deltaEtajj;
  int n_Jets_JESUp, n_Jets_JESDown; 
  float jet_E[100], jet_Pt[100], jet_Eta[100], jet_Phi[100];
  //------------------------
  //set branches on big tree
  //------------------------

  /*
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
  //tree variables
  float weight;
  int NPU;
  int n_Jets, nLooseBTaggedJets, nMediumBTaggedJets;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  float theMR;
  float theRsq, t1Rsq;
  float MET, t1MET;
  int nSelectedPhotons;
  float mGammaGamma, pTGammaGamma, mGammaGammaSC, pTGammaGammaSC, sigmaMoverM;
  float mbbZ, mbbH;
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
  float jet_E[100], jet_Pt[100], jet_Eta[100], jet_Phi[100];
  */
  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("weight", &weight, "weight/F");
    razorTree->Branch("rho", &fixedGridRhoAll, "rho/F");
    razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
    razorTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
    razorTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
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
    razorTree->Branch("allTrackdT", "std::vector<float>",&allTrackdT); //get PDF weights directly from RazorEvents
    razorTree->Branch("allTrackPt", "std::vector<float>",&allTrackPt); //get PDF weights directly from RazorEvents
    razorTree->Branch("allTrackdZ", "std::vector<float>",&allTrackdZ); //get PDF weights directly from RazorEvents
    razorTree->Branch("allTrackPvIndex", "std::vector<int>",&allTrackPvIndex); //get PDF weights directly from RazorEvents
    razorTree->Branch("sf_pdf", "std::vector<float>",&sf_pdf); //sf PDF
      
    //MET filters
    razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
    razorTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
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
    razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
    razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
    razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
    razorTree->Branch("etaGammaGamma", &etaGammaGamma, "etaGammaGamma/F");
    razorTree->Branch("phiGammaGamma", &phiGammaGamma, "phiGammaGamma/F");
    razorTree->Branch("mGammaGammaSC", &mGammaGammaSC, "mGammaGammaSC/F");
    razorTree->Branch("pTGammaGammaSC", &pTGammaGammaSC, "pTGammaGammaSC/F");
    razorTree->Branch("etaGammaGammaSC", &etaGammaGammaSC, "etaGammaGammaSC/F");
    razorTree->Branch("phiGammaGammaSC", &phiGammaGammaSC, "phiGammaGammaSC/F");
    razorTree->Branch("sigmaMoverM", &sigmaMoverM, "sigmaMoverM/F");
    razorTree->Branch("mjj", &mjj, "mjj/F");
    razorTree->Branch("deltaEtajj", &deltaEtajj, "deltaEtajj/F");
    razorTree->Branch("mggjj", &mggjj, "mggjj/F");
    razorTree->Branch("cosThetaCS", &cosThetaCS, "cosThetaCS/F");
    razorTree->Branch("box", &razorbox, "box/I");
    
    razorTree->Branch("pho1E", &Pho_E[0], "pho1E/F");
    razorTree->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
    razorTree->Branch("pho1Time", &Pho_Time[0], "pho1Time/F");
    razorTree->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
    razorTree->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
    razorTree->Branch("pho1SC_E", &PhoSC_E[0], "pho1SC_E/F");
    razorTree->Branch("pho1SC_Pt", &PhoSC_Pt[0], "pho1SC_Pt/F");
    razorTree->Branch("pho1SC_Eta", &PhoSC_Eta[0], "pho1SC_Eta/F");
    razorTree->Branch("pho1SC_Phi", &PhoSC_Phi[0], "pho1SC_Phi/F");
    razorTree->Branch("pho1DefaultSC_E", &PhoDefaultSC_E[0], "pho1DefaultSC_E/F");
    razorTree->Branch("pho1DefaultSC_Pt", &PhoDefaultSC_Pt[0], "pho1DefaultSC_Pt/F");
    razorTree->Branch("pho1DefaultSC_Eta", &PhoDefaultSC_Eta[0], "pho1DefaultSC_Eta/F");
    razorTree->Branch("pho1DefaultSC_Phi", &PhoDefaultSC_Phi[0], "pho1DefaultSC_Phi/F");
    
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
    razorTree->Branch("pho2Time", &Pho_Time[1], "pho2Time/F");
    razorTree->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
    razorTree->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
    razorTree->Branch("pho2SC_E", &PhoSC_E[1], "pho2SC_E/F");
    razorTree->Branch("pho2SC_Pt", &PhoSC_Pt[1], "pho2SC_Pt/F");
    razorTree->Branch("pho2SC_Eta", &PhoSC_Eta[1], "pho2SC_Eta/F");
    razorTree->Branch("pho2SC_Phi", &PhoSC_Phi[1], "pho2SC_Phi/F");
    razorTree->Branch("pho2DefaultSC_E", &PhoDefaultSC_E[1], "pho2DefaultSC_E/F");
    razorTree->Branch("pho2DefaultSC_Pt", &PhoDefaultSC_Pt[1], "pho2DefaultSC_Pt/F");
    razorTree->Branch("pho2DefaultSC_Eta", &PhoDefaultSC_Eta[1], "pho2DefaultSC_Eta/F");
    razorTree->Branch("pho2DefaultSC_Phi", &PhoDefaultSC_Phi[1], "pho2DefaultSC_Phi/F");
    
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
    razorTree->Branch("mggj", mggj, "mggj[n_Jets]/F");
    razorTree->Branch("n_Jets_JESUp", &n_Jets_JESUp, "n_Jets_JESUp/I");
    razorTree->Branch("n_Jets_JESDown", &n_Jets_JESDown, "n_Jets_JESDown/I");
    razorTree->Branch("HLTDecision", HLTDecision, "HLTDecision[300]/O");
    razorTree->Branch("vtxIndex", &vtxIndex, "vtxIndex/I");
    razorTree->Branch("vtxZ", &vtxZCoordinate, "vtxZ/F");
    razorTree->Branch("vtxdZ", &vtxdZCoordinate, "vtxdZ/F");
    razorTree->Branch("vtxT", &vtxTCoordinate, "vtxT/F");
    razorTree->Branch("vtxdT", &vtxdTCoordinate, "vtxdT/F");

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
    //MVA Vertex Selection

/*
    razorTree->Branch("ptasym", &ptasym, "ptasym/F");
    razorTree->Branch("ptbal", &ptbal, "ptbal/F");
    razorTree->Branch("logsumpt2", &logsumpt2, "logsumpt2/F");
    razorTree->Branch("pull_conv", &pull_conv, "pull_conv/F");
    razorTree->Branch("nConv", &nConv, "nConv/F");
  */
  
    razorTree->Branch("nPVAll", &nPVAll, "nPVAll/I");

    razorTree->Branch("allTrackdT_Sig", "std::vector<float>",&allTrackdT_Sig); 
    razorTree->Branch("allTrackPt_Sig", "std::vector<float>",&allTrackPt_Sig); 
    razorTree->Branch("allTrackdT_Bkg", "std::vector<float>",&allTrackdT_Bkg); 
    razorTree->Branch("allTrackPt_Bkg", "std::vector<float>",&allTrackPt_Bkg); 
    razorTree->Branch("allTrackdZ_Sig", "std::vector<float>",&allTrackdZ_Sig); 
    razorTree->Branch("allTrackdZ_Bkg", "std::vector<float>",&allTrackdZ_Bkg); 

    razorTree->Branch("isMatchPv", isMatchPv, "isMatchPv[nPVAll]/I");
    razorTree->Branch("isMaxbdtPv", isMaxbdtPv, "isMaxbdtPv[nPVAll]/I");
    razorTree->Branch("isMaxlogsumpt2Pv", isMaxlogsumpt2Pv, "isMaxlogsumpt2Pv[nPVAll]/I");

    razorTree->Branch("chi2_min_pho_vtx", chi2_min_pho_vtx_all, "chi2_min_pho_vtx[nPVAll]/F");
    razorTree->Branch("min_vtxT", min_vtxT_all, "min_vtxT[nPVAll]/F");
    razorTree->Branch("min_vtxT_gen", &min_vtxT_gen, "min_vtxT_gen/F");
    razorTree->Branch("pho1_vtx_dt_gen", &pho1_vtx_dt_gen, "pho1_vtx_dt_gen/F");
    razorTree->Branch("pho2_vtx_dt_gen", &pho2_vtx_dt_gen, "pho2_vtx_dt_gen/F");
    razorTree->Branch("pho1_vtx_dt", pho1_vtx_dt_all, "pho1_vtx_dt[nPVAll]/F");
    razorTree->Branch("pho2_vtx_dt", pho2_vtx_dt_all, "pho2_vtx_dt[nPVAll]/F");
    razorTree->Branch("chi2_pho_vtx_gen", &chi2_pho_vtx_gen, "chi2_pho_vtx_gen/F");
    razorTree->Branch("rank_chi2_genVtx", &rank_chi2_genVtx, "rank_chi2_genVtx/I");
    razorTree->Branch("rank_chi2_min_genVtx", &rank_chi2_min_genVtx, "rank_chi2_min_genVtx/I");
    razorTree->Branch("chi2_min_pho_vtx_gen", &chi2_min_pho_vtx_gen, "chi2_min_pho_vtx_gen/F");
    razorTree->Branch("chi2_pho_vtx", chi2_pho_vtx_all, "chi2_pho_vtx[nPVAll]/F");
    razorTree->Branch("rank_chi2_pho_vtx", rank_chi2_pho_vtx_all, "rank_chi2_pho_vtx[nPVAll]/I");
    razorTree->Branch("rank_chi2_min_pho_vtx", rank_chi2_min_pho_vtx_all, "rank_chi2_min_pho_vtx[nPVAll]/I");
    razorTree->Branch("ptasym", ptasym_all, "ptasym[nPVAll]/F");
    razorTree->Branch("pvNtrack", pvNtrack_all, "pvNtrack[nPVAll]/I");
    razorTree->Branch("hasLargeTrack", hasLargeTrack, "hasLargeTrack[nPVAll]/I");
    razorTree->Branch("hasLargeTrack_dt", hasLargeTrack_dt, "hasLargeTrack_dt[nPVAll]/I");
    razorTree->Branch("NPV_HasLargeTrack", &NPV_HasLargeTrack, "NPV_HasLargeTrack/I");
    razorTree->Branch("NPV_HasLargeTrack_dt", &NPV_HasLargeTrack_dt, "NPV_HasLargeTrack_dt/I");
    razorTree->Branch("NPV_GEN_HasLargeTrack", &NPV_GEN_HasLargeTrack, "NPV_GEN_HasLargeTrack/I");
    razorTree->Branch("NPV_GEN_HasLargeTrack_dt", &NPV_GEN_HasLargeTrack_dt, "NPV_GEN_HasLargeTrack_dt/I");
    razorTree->Branch("ptbal", ptbal_all, "ptbal[nPVAll]/F");
    razorTree->Branch("logsumpt2", logsumpt2_all, "logsumpt2[nPVAll]/F");
    razorTree->Branch("vertexsumpt", vertexsumpt_all, "vertexsumpt[nPVAll]/F");
    razorTree->Branch("vertexpt", vertexpt_all, "vertexpt[nPVAll]/F");
    razorTree->Branch("diphotonpt", &diphotonpt, "diphotonpt/F");
    razorTree->Branch("virtualVtxdZ", &virtualVtxdZ, "virtualVtxdZ/F");
    razorTree->Branch("limPullToConv", pull_conv_all, "limPullToConv[nPVAll]/F");
    razorTree->Branch("nConv", nConv_all, "nConv[nPVAll]/F");
    razorTree->Branch("vtxdZ_all", vtxdZCoordinate_all, "vetdZ_all[nPVAll]/F");
    razorTree->Branch("vtxdX_all", vtxdXCoordinate_all, "vetdX_all[nPVAll]/F");
    razorTree->Branch("vtxdY_all", vtxdYCoordinate_all, "vetdY_all[nPVAll]/F");
    razorTree->Branch("vtxdT_all", vtxdTCoordinate_all, "vetdT_all[nPVAll]/F");
    razorTree->Branch("vtxZ_all", vtxZCoordinate_all, "vetZ_all[nPVAll]/F");
    razorTree->Branch("vtxX_all", vtxXCoordinate_all, "vetX_all[nPVAll]/F");
    razorTree->Branch("vtxY_all", vtxYCoordinate_all, "vetY_all[nPVAll]/F");
    razorTree->Branch("vtxT_all", vtxTCoordinate_all, "vetT_all[nPVAll]/F");


  }
  //set branches on all trees
  else{ 
    for(auto& thisBox : razorBoxes){
      thisBox.second->Branch("weight", &weight, "weight/F");
      thisBox.second->Branch("run", &run, "run/i");
      thisBox.second->Branch("lumi", &lumi, "lumi/i");
      thisBox.second->Branch("event", &event, "event/i");
      thisBox.second->Branch("passedDiphotonTrigger", &passedDiphotonTrigger, "passedDiphotonTrigger/O");
      thisBox.second->Branch("NPU", &NPU, "npu/i");
      thisBox.second->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
      thisBox.second->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
      thisBox.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      thisBox.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      thisBox.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      thisBox.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      thisBox.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      thisBox.second->Branch("MR", &theMR, "MR/F");
      thisBox.second->Branch("Rsq", &theRsq, "Rsq/F");
      thisBox.second->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
      thisBox.second->Branch("MET", &MET, "MET/F");
      thisBox.second->Branch("t1MET", &t1MET, "t1MET/F");
      thisBox.second->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
      thisBox.second->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
      thisBox.second->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
      thisBox.second->Branch("mGammaGammaSC", &mGammaGammaSC, "mGammaGammaSC/F");
      thisBox.second->Branch("pTGammaGammaSC", &pTGammaGammaSC, "pTGammaGammaSC/F");
      thisBox.second->Branch("sigmaMoverM", &sigmaMoverM, "sigmaMoverM/F");
      
      thisBox.second->Branch("pho1E", &Pho_E[0], "pho1E/F");
      thisBox.second->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
      thisBox.second->Branch("Pho1Eta", &Pho_Eta[0], "pho1Eta/F");
      thisBox.second->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
      thisBox.second->Branch("pho1SC_E", &PhoSC_E[0], "pho1SC_E/F");
      thisBox.second->Branch("pho1SC_Pt", &PhoSC_Pt[0], "pho1SC_Pt/F");
      thisBox.second->Branch("pho1SC_Eta", &PhoSC_Eta[0], "pho1SC_Eta/F");
      thisBox.second->Branch("pho1SC_Phi", &PhoSC_Phi[0], "pho1SC_Phi/F");
      thisBox.second->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
      thisBox.second->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
      thisBox.second->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
      thisBox.second->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
      thisBox.second->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
      thisBox.second->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
      thisBox.second->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
      thisBox.second->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
      thisBox.second->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");
      
      thisBox.second->Branch("pho2E", &Pho_E[1], "pho2E/F");
      thisBox.second->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
      thisBox.second->Branch("Pho2Eta", &Pho_Eta[1], "pho2Eta/F");
      thisBox.second->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
      thisBox.second->Branch("pho2SC_E", &PhoSC_E[1], "pho2SC_E/F");
      thisBox.second->Branch("pho2SC_Pt", &PhoSC_Pt[1], "pho2SC_Pt/F");
      thisBox.second->Branch("pho2SC_Eta", &PhoSC_Eta[1], "pho2SC_Eta/F");
      thisBox.second->Branch("pho2SC_Phi", &PhoSC_Phi[1], "pho2SC_Phi/F");
      thisBox.second->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
      thisBox.second->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
      thisBox.second->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
      thisBox.second->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
      thisBox.second->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
      thisBox.second->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
      thisBox.second->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
      thisBox.second->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
      thisBox.second->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O");
      
      thisBox.second->Branch("mbbZ", &mbbZ, "mbbZ/F");
      thisBox.second->Branch("mbbH", &mbbH, "mbbH/F");
      
      thisBox.second->Branch("n_Jets", &n_Jets, "n_Jets/I");
      thisBox.second->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
      thisBox.second->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
      thisBox.second->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
      thisBox.second->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
      thisBox.second->Branch("HLTDecision", HLTDecision, "HLTDecision[300]/O");

      //GenParticles
      thisBox.second->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
      thisBox.second->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
      thisBox.second->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
      thisBox.second->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
      thisBox.second->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
      thisBox.second->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
      thisBox.second->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
      thisBox.second->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
      thisBox.second->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
    }
  }
  
  bool use25nsSelection = true;
  std::cout << "use25nsSelection-->" << use25nsSelection << std::endl;
  //begin loop
  if (fChain == 0) return;
  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  if( _info )  std::cout<<nentries<<std::endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {


  allTrackdT_Sig->clear();
  allTrackPt_Sig->clear();
  allTrackdT_Bkg->clear();
  allTrackPt_Bkg->clear();
  allTrackdZ_Sig->clear();
  allTrackdZ_Bkg->clear();


  NPV_HasLargeTrack = 0;
  NPV_HasLargeTrack_dt = 0;
  NPV_GEN_HasLargeTrack = 0;
  NPV_GEN_HasLargeTrack_dt = 0;

  chi2_pho_vtx_gen = 0.0;
  rank_chi2_genVtx = 0;
  rank_chi2_min_genVtx = 0;
  chi2_min_pho_vtx_gen = 0.0;
  min_vtxT_gen = 0.0;
  pho1_vtx_dt_gen = 0.0;
  pho2_vtx_dt_gen = 0.0;
  for(int i=0;i<500;i++)
  {
	isMatchPv[i]=0;
	isMaxbdtPv[i]=0;
	isMaxlogsumpt2Pv[i]=0;
	chi2_min_pho_vtx_all[i]=0.0;
	min_vtxT_all[i]=0.0;
	chi2_pho_vtx_all[i]=0.0;
	rank_chi2_pho_vtx_all[i]=1;
	rank_chi2_min_pho_vtx_all[i]=1;
	pho1_vtx_dt_all[i]=0.0;
	pho2_vtx_dt_all[i]=0.0;
	ptasym_all[i]=0.0;
	pvNtrack_all[i]=0;
	hasLargeTrack[i]=0;
	hasLargeTrack_dt[i]=0;
	ptbal_all[i]=0.0;
     	logsumpt2_all[i]=0.0;	
     	vertexsumpt_all[i]=0.0;	
     	vertexpt_all[i]=0.0;	
  	pull_conv_all[i]=0.0;
	nConv_all[i]=0.0;
	vtxdZCoordinate_all[i]=0.0;
	vtxdXCoordinate_all[i]=0.0;
	vtxdYCoordinate_all[i]=0.0;
	vtxdTCoordinate_all[i]=0.0;
	vtxZCoordinate_all[i]=0.0;
	vtxXCoordinate_all[i]=0.0;
	vtxYCoordinate_all[i]=0.0;
	vtxTCoordinate_all[i]=0.0;
  }
        diphotonpt = 0.0;
        virtualVtxdZ = 0.0;
 
    //begin event
    if( _info && (jentry % 10 == 0) ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    if(useOddEvent && eventNum%2==0) continue;
 
    if( _info )  std::cout<<"DEBUG_ZZC 001"<<std::endl;

    //fill normalization histogram    
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;
    SumWeights->Fill(1.0, weight);
    
    if( _info )  std::cout<<"DEBUG_ZZC 002"<<std::endl;
    
    //reset tree variables
    ISRSystWeightUp   = 1.0;
    ISRSystWeightDown = 1.0;
    pileupWeight      = 1.0;
    pileupWeightUp    = 1.0;
    pileupWeightDown  = 1.0;
    
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
    mGammaGamma     = -1;
    pTGammaGamma    = -1;
    etaGammaGamma   = -666;
    phiGammaGamma   = -666;
    mGammaGammaSC   = -1;
    pTGammaGammaSC  = -1;
    etaGammaGammaSC = -666;
    phiGammaGammaSC = -666;
    
    mjj    = -1;
    mggjj  = -1;
    mbbZ   = -1;
    mbbH   = -1;
    mbbZ_L = -1;
    mbbH_L = -1;
    deltaEtajj = -1;
    cosThetaCS = -666;
    run = runNum;
    lumi = lumiNum; 
    event = eventNum;
    passedDiphotonTrigger = false;
    vtxIndex = -1;
    vtxZCoordinate = -999;
    vtxdZCoordinate = -999;
    vtxTCoordinate = -999;
    vtxdTCoordinate = -999;

    if( _info )  std::cout<<"DEBUG_ZZC 003"<<std::endl;
    //selected photons variables
    for ( int i = 0; i < 2; i++ )
      {
	Pho_E[i]                  = -99.;
	Pho_Pt[i]                 = -99.;
	Pho_Time[i]                = -99.;
	Pho_Eta[i]                = -99.;
	Pho_Phi[i]                = -99.;
	PhoSC_E[i]                = -99.;
	PhoSC_Pt[i]               = -99.;
	PhoSC_Eta[i]              = -99.;
	PhoSC_Phi[i]              = -99.;
	PhoDefaultSC_E[i]         = -99.;
	PhoDefaultSC_Pt[i]        = -99.;
	PhoDefaultSC_Eta[i]       = -99.;
	PhoDefaultSC_Phi[i]       = -99.;
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
    for ( int i = 0; i < 100; i++ )
      {
	jet_E[i]   = -99.;
	jet_Pt[i]  = -99.;
	jet_Eta[i] = -99.;
	jet_Phi[i] = -99.;
	mggj[i]    = -1;
      }


    if( _info )  std::cout<<"DEBUG_ZZC 004"<<std::endl;
      
    //------------------
    //Pileup reweighting
    //------------------
    //double pileupWeight = 1.0;
    if( !isData )
      {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) 
	  {
	    if (BunchXing[i] == 0) 
	      {
		NPU = nPUmean[i];
	      }
	  }
	puhisto->Fill(NPU);
	//NOTE: reweight with nPV for now
	//pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPV));
	pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
	pileupWeightUp = pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU)) / pileupWeight;
	pileupWeightDown = pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU)) / pileupWeight;    
      }
      
    if( _info )  std::cout<<"DEBUG_ZZC 005"<<std::endl;
    /////////////////////////////////
    //Scale and PDF variations
    /////////////////////////////////
    if( !isData )
      {
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
	else
	  {
	    //std::cerr << "[WARNING]: scaleWeight size less than 9, actual size is: " << (*scaleWeights).size() << std::endl;
	  }
	
	sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
	for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) 
	  {
	    sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	    SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	  }
	
      }
    
    if(combineTrees) razorbox = LowRes;
    
    if( _info )  std::cout<<"DEBUG_ZZC 006"<<std::endl;
    //TODO: triggers!
    // bool passedDiphotonTrigger = true;
    passedDiphotonTrigger = ( HLTDecision[40] );
    //if(!passedDiphotonTrigger) continue;
    
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
      if(!isLooseElectron(i,use25nsSelection)) continue; 
      if(elePt[i] < 10) continue;
      if(abs(eleEta[i]) > 2.5) continue;
      
      nLooseElectrons++;
      
      if(isTightElectron(i,use25nsSelection))
	{ 
	  nTightElectrons++;
	}
    }
    //tau selection
    for(int i = 0; i < nTaus; i++){
      if(!isTightTau(i)) continue; 
      nTightTaus++;
    }


    if( _info )  std::cout<<"DEBUG_ZZC 007"<<std::endl;
    bool  _phodebug = false;
    bool  _debug = false;
    if ( 
	(run == 257969  && event == 862653894) ||
	(run == 258158  && event == 2009129070) ||
	(run == 258443  && event == 79092577 ) ||
	(run == 259721  && event == 69584193 ) ||
	(run == 260536  && event == 246155173 ) ||
	(run == 258136  && event == 72682254 ) ||
	(run == 258443  && event == 179542426 ) ||
	(run == 259822  && event == 191351229 ) ||
	(run == 260427  && event == 86898837 ) 
	 ) {
      _phodebug = true;
      _debug = true;
      std::cout << "Event: " << run << " " << event << "\n";
    }

    if( _info )  std::cout<<"DEBUG_ZZC 008"<<std::endl;
    //photon selection
    vector<TLorentzVector> GoodPhotons;
    vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
    vector<bool> GoodPhotonPassesIso; //store whether each photon is isolated
    std::vector< NewPhotonCandidate > phoCand;//PhotonCandidate defined in RazorAuxPhoton.hh
    
    int nPhotonsAbove40GeV = 0;
    //std::cout << "nphotons: " << nPhotons << std::endl;
    for( int i = 0; i < nPhotons; i++ ) {

      double scale = photonCorrector->ScaleCorrection(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]));
      double smear = photonCorrector->getSmearingSigma(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]), 0., 0.); 
 
      //Require Loose Photon ID
      /*
      if ( _phodebug ) {
	std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] 
		  << " | " << scale << " " << smear << " "
		  << std::endl;
      }
      if ( !photonPassLooseIDWithoutEleVetoExo15004(i) ) {
	if ( _phodebug ) std::cout << "[DEBUG]: failed Exo15004 ID" << std::endl;
	continue;
      }
     */

      //Require electron veto
      if (doEleVeto) {
	if ( !pho_passEleVeto[i] ) {
	  if ( _phodebug ) std::cout << "[DEBUG]: failed Exo15004 Ele Veto" << std::endl;
	  continue;
	}	  
      }
	
      if ( _phodebug ) { 
	std::cout << "[DEBUG]: Exo15004 Isolation" << " "
		  << pho_pfIsoChargedHadronIso[i] << " ";
	double effAreaPhotons = 0.0;
	double eta = pho_superClusterEta[i];
	getPhotonEffAreaExo15004( eta, effAreaPhotons );
	std::cout << pho_pfIsoPhotonIso[i] << " " << fixedGridRhoFastjetAll << " " << effAreaPhotons << " " << pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons << " ";
	double cut = 0;
	if( fabs(pho_superClusterEta[i]) < 1.4442 )
	  {
	    std::cout << (2.75 - 2.5) + 0.0045*phoPt[i] ;
	    cut = (2.75 - 2.5) + 0.0045*phoPt[i];
	  } 
	else if ( abs(pho_superClusterEta[i]) < 2.0 )
	  {
	    std::cout <<  (2.0 - 2.5) + 0.0045*phoPt[i] ;
	    cut = (2.0 - 2.5) + 0.0045*phoPt[i];
	  }
	else
	  {
	    std::cout <<   (2.0 - 2.5) + 0.0045*phoPt[i] ;
	    cut = (2.0 - 2.5) + 0.0045*phoPt[i];
	  }
	    
	std::cout << " "
		  << bool( pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons < cut) << " "
		  << std::endl;
      }
	  

    if( _info )  std::cout<<"DEBUG_ZZC 009"<<std::endl;
      //Require Loose Isolation
      /*
      if ( !photonPassLooseIsoExo15004(i) )
	{
	  if ( _phodebug ) { 
	    std::cout << "[DEBUG]: failed Exo15004 Isolation" << " "
		      << pho_pfIsoChargedHadronIso[i] << " ";
	    double effAreaPhotons = 0.0;
	    double eta = pho_superClusterEta[i];
	    getPhotonEffAreaExo15004( eta, effAreaPhotons );
	    std::cout << pho_pfIsoPhotonIso[i] << " " << fixedGridRhoAll << " " << effAreaPhotons << " " << pho_pfIsoPhotonIso[i] - fixedGridRhoAll*effAreaPhotons << " ";
	    if( fabs(pho_superClusterEta[i]) < 1.4442 )
	      {
		std::cout << (2.75 - 2.5) + 0.0045*phoPt[i] ;
	      } 
	    else if ( abs(pho_superClusterEta[i]) < 2.0 )
	      {
		std::cout <<  (2.0 - 2.5) + 0.0030*phoPt[i] ;
	      }
	    else
	      {
		std::cout <<   (2.0 - 2.5) + 0.0030*phoPt[i] ;
	      }

	    std::cout << " "
		      << std::endl;
	  }
	  continue;
	}
	*/
      //Defining Corrected Photon momentum
      float pho_pt_corr = phoPt[i];//nominal pt
      if (doPhotonScaleCorrection) {
	if (isData) {
	  pho_pt_corr = phoPt[i]*scale; 
	  if (_phodebug) std::cout << "[DEBUG] : Photon Energy Scale Corrections: " << phoPt[i] << " * " << scale << " --> " << pho_pt_corr << "\n";
	} else {
	  pho_pt_corr = phoPt[i]*(1+smear*random.Gaus());
	}
      }

     

    if( _info )  std::cout<<"DEBUG_ZZC 010"<<std::endl;
      TVector3 vec;
      vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );
	
      if ( phoPt[i] < 20.0 )
	//if ( phoE[i]/cosh( phoEta[i] ) < 24.0 )
	{
	  if ( _phodebug ) std::cout << "[DEBUG]: failed pt" << std::endl;
	  //continue;
	}
	
      if( fabs(pho_superClusterEta[i]) > 2.5 )
	{
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
      if( phoPt[i] > 40.0 ) nPhotonsAbove40GeV++;
      //setting up photon 4-momentum with zero mass
      TLorentzVector thisPhoton;
      thisPhoton.SetVectM( vec, .0 );

    if( _info )  std::cout<<"DEBUG_ZZC 011"<<std::endl;
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
	
      //TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_superClusterEnergy[i] );
      TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_RegressionE[i] ); 
	
      //std::cout << "phoSC_Pt: " << phoSC.Pt() << " phoCorrPt: " << thisPhoton.Pt() << std::endl;
      //std::cout << "phoSC_Eta: " << phoSC.Eta() << " originalSC_Eta: " << pho_superClusterEta[i] << std::endl;
      //Filling Photon Candidate
      NewPhotonCandidate tmp_phoCand;
      tmp_phoCand.Index = i;
      tmp_phoCand.photon = thisPhoton;
      tmp_phoCand.photonSC = phoSC;
      //std::cout << pho_superClusterEta[i] << std::endl;
      //std::cout << pho_superClusterEnergy[i] << std::endl;
      tmp_phoCand.scE   = pho_superClusterEnergy[i];
      //std::cout << pho_superClusterEnergy[i] << std::endl;
      tmp_phoCand.scPt  = pho_superClusterEnergy[i]/cosh( pho_superClusterEta[i] );
      tmp_phoCand.scEta = pho_superClusterEta[i];
      tmp_phoCand.scPhi = pho_superClusterPhi[i];
      tmp_phoCand.scSeedX = pho_superClusterSeedX[i];
      tmp_phoCand.scSeedY = pho_superClusterSeedY[i];
      tmp_phoCand.scSeedZ = pho_superClusterSeedZ[i];
      tmp_phoCand.scSeedT = pho_superClusterSeedT[i];
      tmp_phoCand.scX = pho_superClusterX[i];
      tmp_phoCand.scY = pho_superClusterY[i];
      tmp_phoCand.scZ = pho_superClusterZ[i];
      tmp_phoCand.SigmaIetaIeta = phoSigmaIetaIeta[i];
      tmp_phoCand.R9 = phoR9[i];
      tmp_phoCand.HoverE = pho_HoverE[i];
      /*tmp_phoCand.sumChargedHadronPt = pho_sumChargedHadronPt[i];
	tmp_phoCand.sumNeutralHadronEt = pho_sumNeutralHadronEt[i];
	tmp_phoCand.sumPhotonEt = pho_sumPhotonEt[i];
      */
      tmp_phoCand.sumChargedHadronPt = pho_pfIsoChargedHadronIso[i];
      tmp_phoCand.sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[i];
      tmp_phoCand.sumPhotonEt  = pho_pfIsoPhotonIso[i];
      tmp_phoCand.sigmaEOverE  = pho_RegressionEUncertainty[i]/pho_RegressionE[i];
      tmp_phoCand._passEleVeto = pho_passEleVeto[i];
      tmp_phoCand._passIso = photonPassLooseIsoExo15004(i);
      tmp_phoCand.convType = pho_convType[i];
      tmp_phoCand.convTrkZ = pho_convTrkZ[i];
      tmp_phoCand.convTrkClusZ = pho_convTrkClusZ[i];
      for (int ipv=0; ipv<nPVAll; ++ipv) {
        tmp_phoCand.vtxSumPx[ipv] = pho_vtxSumPx[i][ipv];
        tmp_phoCand.vtxSumPy[ipv] = pho_vtxSumPy[i][ipv];
      }
      //std::cout << "phoCand: " << phoCand.size() << std::endl;
      phoCand.push_back( tmp_phoCand );
      //std::cout << "phoCand: " << phoCand.size() << std::endl;
      nSelectedPhotons++;
    }
    
    //std::cout << "[DEBUG]: nphotons--> " << phoCand.size() << " " << nSelectedPhotons << std::endl;
    
    //if there is no photon with pT above 40 GeV, reject the event
    if( nPhotonsAbove40GeV == 0 )
      {
	if ( _debug ) std::cout << "[DEBUG]: no photons above 40 GeV, nphotons: " << phoCand.size() << std::endl;
	//continue;//offline
      }
    
    if ( phoCand.size() < 2 )
      {
	if ( _debug ) std::cout << "[INFO]: not enough photon, nphotons: " << phoCand.size() << std::endl;
	for(int i = 0; i < nPhotons; i++)
	  {
	    if ( _debug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] 
				    << " SIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
	  }
	continue;
      }
    

    if( _info )  std::cout<<"DEBUG_ZZC 012"<<std::endl;
    if ( _debug ) std::cout << "[DEBUG]: nphotons--> " << phoCand.size() << " " << nSelectedPhotons << std::endl;
    //find the "best" photon pair, higher Pt!
    TLorentzVector HiggsCandidate(0,0,0,0);
    TLorentzVector HiggsCandidateSC(0,0,0,0);
    int HiggsPhoIndex1 = -1;
    int HiggsPhoIndex2 = -1;
    double bestSumPt = -99.;
    std::vector< NewPhotonCandidate > phoSelectedCand;
    NewPhotonCandidate bestCand[2];
    
    for( size_t i = 0; i < phoCand.size(); i++ )
      {
	for( size_t j = i+1; j < phoCand.size(); j++ )//I like this logic better, I find it easier to understand
	  {
	    NewPhotonCandidate pho1 = phoCand[i];
	    NewPhotonCandidate pho2 = phoCand[j];
            
            if (doMVAVertex) {

              float pho1E = pho1.photon.E();
              float pho2E = pho2.photon.E();
              
              float maxbdtval = -99.;
              float maxlogsumpt2 = -99.;
              int ipvmax = 0;
              int ipvmaxlogsumpt2 = 0;
                            
              for (int ipv=0; ipv<nPVAll; ++ipv) {
                float vtxX = pvAllX[ipv];
                float vtxY = pvAllY[ipv];
                float vtxZ = pvAllZ[ipv];
                float vtxT = pvAllT[ipv];
                
                TVector3 pho1dir(pho1.scX-vtxX,pho1.scY-vtxY,pho1.scZ-vtxZ);
                TVector3 pho2dir(pho2.scX-vtxX,pho2.scY-vtxY,pho2.scZ-vtxZ);
                TVector3 diphomom = pho1E*pho1dir.Unit() + pho2E*pho2dir.Unit();
                TVector3 diphoPt(diphomom.x(), diphomom.y(), 0.);
                if(!useTiming)
		{ 
                TVector3 vtxSumPt(pvAllSumPx[ipv]-pho1.vtxSumPx[ipv]-pho2.vtxSumPx[ipv], pvAllSumPy[ipv]-pho1.vtxSumPy[ipv]-pho2.vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		
		}
		logsumpt2 = pvAllLogSumPtSq[ipv];
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());

		}
	  	else	
		{ 
                TVector3 vtxSumPt(pvAllSumPx_dt[ipv]-pho1.vtxSumPx[ipv]-pho2.vtxSumPx[ipv], pvAllSumPy_dt[ipv]-pho1.vtxSumPy[ipv]-pho2.vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		}


                logsumpt2 = pvAllLogSumPtSq_dt[ipv];
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
                }
	
                bool hasconv1 = pho1.convType>=0;
                bool hasconv2 = pho2.convType>=0;
                
                float convz1 = -99.;
                float convsz1 = -99.;
                float convz2 = -99.;
                float convsz2 = -99.;
                                
                nConv = 0.;
                float convz = -99.;
                float convsz = -99.;
                if (hasconv1) {
                  convz1 = SIGMATRKZ[pho1.convType] < SIGMATRKCLUSZ[pho1.convType] ? pho1.convTrkZ : pho1.convTrkClusZ;
                  convsz1 = std::min(SIGMATRKZ[pho1.convType],SIGMATRKCLUSZ[pho1.convType]);                  
                  
                  convz = convz1;
                  convsz = convsz1;
                  nConv = 1;
                }
                if (hasconv2) {
                  convz2 = SIGMATRKZ[pho2.convType] < SIGMATRKCLUSZ[pho2.convType] ? pho2.convTrkZ : pho2.convTrkClusZ;
                  convsz2 = std::min(SIGMATRKZ[pho2.convType],SIGMATRKCLUSZ[pho2.convType]);
                  
                  convz = convz2;
                  convsz = convsz2;
                  nConv = 1;
                }
                if (hasconv1 && hasconv2) {
                  double w1 = 1./convsz1/convsz1;
                  double w2 = 1./convsz2/convsz2;
                  convz = (w1*convz1 + w2*convz2)/(w1+w2);
                  convsz = sqrt(1./(w1+w2));
                  nConv = 2;
                }
                
                if (hasconv1 || hasconv2) {
                  pull_conv = std::abs(convz-vtxZ)/convsz;
                  pull_conv = std::min(10.f,pull_conv); 
                }
                else {
                  pull_conv = 10.;
                }
              
		if(usePhoChi2)
		{
			chi2_pho_vtx_train = 0.0;
			logchi2_pho_vtx_train = 0.0;
		TVector3 pho1_000(pho1.scSeedX, pho1.scSeedY,pho1.scSeedZ);
		TVector3 pho2_000(pho2.scSeedX,pho2.scSeedY,pho2.scSeedZ);
		TVector3 pho1vtx(pho1.scSeedX-vtxX,pho1.scSeedY-vtxY,pho1.scSeedZ-vtxZ);
                TVector3 pho2vtx(pho2.scSeedX-vtxX,pho2.scSeedY-vtxY,pho2.scSeedZ-vtxZ);
		double CmToNs = 0.1/2.99792458;
		double globalOffset = 0.0111;//11.1ps global offset from the result...
		float pho1Time_m_NoSmear = pho1.scSeedT + CmToNs*pho1_000.Mag() - globalOffset;	
		float pho2Time_m_NoSmear = pho2.scSeedT + CmToNs*pho2_000.Mag() - globalOffset;	
		double smear_phoT = 0.03;//30 ps smearing  
		float pho1Time_m = randomPhotonTime.Gaus(pho1Time_m_NoSmear, smear_phoT);
		float pho2Time_m = randomPhotonTime.Gaus(pho2Time_m_NoSmear, smear_phoT);
			if(useTiming)
			{
		//assuming vertex timing known
		float pho1Time_e = vtxT + CmToNs*pho1vtx.Mag();
		float pho2Time_e = vtxT + CmToNs*pho2vtx.Mag();
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e)/(smear_phoT), 2.0);	
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
		}
			else
			{
		double vtxT_min_chi2 = -0.5*(CmToNs*pho1vtx.Mag()+CmToNs*pho2vtx.Mag()-pho1Time_m-pho2Time_m);
		double pho1Time_e_tmp = vtxT_min_chi2 + CmToNs*pho1vtx.Mag();
		double pho2Time_e_tmp = vtxT_min_chi2 + CmToNs*pho2vtx.Mag();
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);	
	
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
			}	
		}  
                float bdtval = vtxmvareader->EvaluateMVA("BDT");
		if(usePhoChi2) bdtval = vtxmvareader_chi2->EvaluateMVA("BDT");
		if (_debug) std::cout << "Vertex: " << ipv << " " << pvAllZ[ipv] << " : " 
				      << ptasym << " " << ptbal << " " << logsumpt2 << " " << pull_conv << " " << nConv << " : "
				      << bdtval << "\n";
                              
                if (bdtval > maxbdtval) {
                  maxbdtval = bdtval;
                  ipvmax = ipv;
                }
		if(logsumpt2 > maxlogsumpt2)
		{
		  maxlogsumpt2 = logsumpt2;
		  ipvmaxlogsumpt2 = ipv;
		}
                
              }

	     
	       
              int ipvmatch = 0;
	      
              double minDist = 999999.0;
                            
              for (int ipv=0; ipv<nPVAll; ++ipv) {
                float vtxX = pvAllX[ipv];
                float vtxY = pvAllY[ipv];
                float vtxZ = pvAllZ[ipv];
                float vtxT = pvAllT[ipv];
	  	if(!useTiming)
		{	
		double dist_this = sqrt(pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0));
		if(dist_this < minDist)
		{
		    ipvmatch = ipv;
		    minDist = dist_this;
		}
		}
		else
		{
		  double chi2_zt = (pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0))/(0.1*0.1)+ pow(vtxT - genVertexT, 2.0)/(0.06*0.06);
		if(chi2_zt < minDist)
		{
		     ipvmatch = ipv;
		     minDist = chi2_zt;
		}	
		}
               }
 
	      vtxIndex = ipvmax;
	      
	       
	      if (_debug) std::cout << "maxbdt Vertex: " << ipvmax << "\n\n";
	      if (_debug) std::cout << "gen-match Vertex: " << ipvmatch << "   dz: "<<sqrt(pow(pvAllZ[ipvmatch] - genVertexZ, 2.0)+pow(pvAllY[ipvmatch] - genVertexY, 2.0)+pow(pvAllX[ipvmatch] - genVertexX, 2.0))<<" ("<<(pow(pvAllZ[ipvmatch] - genVertexZ, 2.0)+pow(pvAllY[ipvmatch] - genVertexY, 2.0)+pow(pvAllX[ipvmatch] - genVertexX, 2.0))/(0.1*0.1)<<")  dT: "<<pvAllT[ipvmatch] - genVertexT<<"  ("<<pow(pvAllT[ipvmatch] - genVertexT,2.0)/(0.06*0.06)<<")  chi2_zt: "<< minDist<<"\n\n";
	      if (_debug) std::cout << "maxlogsumpt2 Vertex: " << ipvmaxlogsumpt2 << "\n\n";
	      
	      if (_debug) std::cout << "Vertex Selected: " << vtxIndex << "\n\n";

	      vtxZCoordinate =  pvAllZ[vtxIndex];
	      vtxdZCoordinate =  sqrt(pow(pvAllZ[vtxIndex] - genVertexZ, 2.0) + pow(pvAllX[vtxIndex] - genVertexX, 2.0) + pow(pvAllY[vtxIndex] - genVertexY, 2.0) );
	      if(useTiming)
		{
		vtxTCoordinate =  pvAllT[vtxIndex];
	      	vtxdTCoordinate =  pvAllT[vtxIndex] - genVertexT;
		} 
              float vtxX = pvAllX[vtxIndex];
              float vtxY = pvAllY[vtxIndex];
              float vtxZ = pvAllZ[vtxIndex];
              
              TVector3 pho1dir(pho1.scX-vtxX,pho1.scY-vtxY,pho1.scZ-vtxZ);
              pho1dir = pho1dir.Unit();
              pho1.photon.SetPxPyPzE(pho1dir.x()*pho1E,pho1dir.y()*pho1E,pho1dir.z()*pho1E,pho1E);
              
              TVector3 pho2dir(pho2.scX-vtxX,pho2.scY-vtxY,pho2.scZ-vtxZ);
              pho2dir = pho2dir.Unit();
              pho2.photon.SetPxPyPzE(pho2dir.x()*pho2E,pho2dir.y()*pho2E,pho2dir.z()*pho2E,pho2E);
              
            } //end if do MVA Vertex
            
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
		//continue;
	      }
	    //need diphoton mass between > 100 GeV as in AN (April 1st)
	    double diphotonMass = (pho1.photon + pho2.photon).M();
	    if ( _debug )
	      {
		std::cout << "[DEBUG] Diphoton Sum pT: " << pho1.photon.Pt() + pho2.photon.Pt() << std::endl;
	      }
	    
	    if( diphotonMass < 100 )
	      {
		if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 100 GeV: mgg->" << diphotonMass << std::endl;
		if ( _debug ) std::cout << "... pho1Pt: " << pho1.photon.Pt()  << " pho2Pt: " << pho2.photon.Pt()  << std::endl;
		//continue;//apply offline
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
	    
	  }//second photon loop
      }//first photon loop
  

    //apply pt/m cut and mgg cut
    
    if(bestCand[0].photon.Pt()/HiggsCandidate.M() < (1.0/3.0)) 
	{
	if ( _debug ) std::cout << "rejected: pt1/mgg = "<< bestCand[0].photon.Pt()/HiggsCandidate.M()<<std::endl;
	continue;
	}
    if(bestCand[1].photon.Pt()/HiggsCandidate.M() < (1.0/4.0)) 
	{
	if ( _debug ) std::cout << "rejected: pt2/mgg = "<< bestCand[1].photon.Pt()/HiggsCandidate.M()<<std::endl;
	continue;
	}
    if(HiggsCandidate.M()<100.0) 
	{
	if ( _debug ) std::cout << "rejected: mgg = "<< HiggsCandidate.M()<<std::endl;
	continue;
	}
    //reject fake photons
    bool matchGenPhoton1 = false;
    bool matchGenPhoton2 = false;

    for(int g=0;g<nGenParticle;g++)
	{
    	if (deltaR(gParticleEta[g] , gParticlePhi[g], bestCand[0].photon.Eta(),bestCand[0].photon.Phi()) < 0.1) matchGenPhoton1 = true;
    	if (deltaR(gParticleEta[g] , gParticlePhi[g], bestCand[1].photon.Eta(),bestCand[1].photon.Phi()) < 0.1) matchGenPhoton2 = true;
	}
    if(!(matchGenPhoton1&&matchGenPhoton2))
	{
	if ( _debug ) std::cout<<"rejected: matchGenPhoton1 = "<<matchGenPhoton1<<"  matchGenPhoton2 = "<<matchGenPhoton2<<std::endl;
	continue;
	}
    //save MVA vertex

            if (doMVAVertex) {

              float pho1E = bestCand[0].photon.E();
              float pho2E = bestCand[1].photon.E();
             
  	      float maxbdtval = -99.;
              float maxlogsumpt2 = -99.;
              int ipvmax = 0;
              int ipvmatch = 0;
              int ipvmaxlogsumpt2 = 0;
	      
              double minDist = 999999.0;
                           
 		 
              for (int ipv=0; ipv<nPVAll; ++ipv) {
                float vtxX = pvAllX[ipv];
                float vtxY = pvAllY[ipv];
                float vtxZ = pvAllZ[ipv];
                float vtxT = pvAllT[ipv];
	  	if(!useTiming)
		{	
		double dist_this = sqrt(pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0));
		if(dist_this < minDist)
		{
		    ipvmatch = ipv;
		    minDist = dist_this;
		}
		}
		else
		{
		  double chi2_zt = (pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0))/(0.1*0.1)+ pow(vtxT - genVertexT, 2.0)/(0.06*0.06);
		if(chi2_zt < minDist)
		{
		     ipvmatch = ipv;
		     minDist = chi2_zt;
		}	
		}
                


                TVector3 pho1dir(bestCand[0].scX-vtxX,bestCand[0].scY-vtxY,bestCand[0].scZ-vtxZ);
                TVector3 pho2dir(bestCand[1].scX-vtxX,bestCand[1].scY-vtxY,bestCand[1].scZ-vtxZ);
                TVector3 diphomom = pho1E*pho1dir.Unit() + pho2E*pho2dir.Unit();
                TVector3 diphoPt(diphomom.x(), diphomom.y(), 0.);
                
                if(!useTiming)
		{ 
                TVector3 vtxSumPt(pvAllSumPx[ipv]-bestCand[0].vtxSumPx[ipv]-bestCand[1].vtxSumPx[ipv], pvAllSumPy[ipv]-bestCand[0].vtxSumPy[ipv]-bestCand[1].vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		}
		logsumpt2 = pvAllLogSumPtSq[ipv];
		vertexsumpt = pvAllSumPt[ipv];
		vertexpt = vtxSumPt.Mag();
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
                }
		
		else
		{
		TVector3 vtxSumPt(pvAllSumPx_dt[ipv]-bestCand[0].vtxSumPx[ipv]-bestCand[1].vtxSumPx[ipv], pvAllSumPy_dt[ipv]-bestCand[0].vtxSumPy[ipv]-bestCand[1].vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		}
		logsumpt2 = pvAllLogSumPtSq_dt[ipv];
		vertexsumpt = pvAllSumPt_dt[ipv];
		//cout<<"vertexsumpt "<<ipv<<"  "<<pvAllSumPt[ipv]<<endl;	
		vertexpt = vtxSumPt.Mag();
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
		}

		diphotonpt = diphoPt.Mag();

                bool hasconv1 = bestCand[0].convType>=0;
                bool hasconv2 = bestCand[1].convType>=0;
                
                float convz1 = -99.;
                float convsz1 = -99.;
                float convz2 = -99.;
                float convsz2 = -99.;
                                
                nConv = 0.;
                float convz = -99.;
                float convsz = -99.;
                if (hasconv1) {
                  convz1 = SIGMATRKZ[bestCand[0].convType] < SIGMATRKCLUSZ[bestCand[0].convType] ? bestCand[0].convTrkZ : bestCand[0].convTrkClusZ;
                  convsz1 = std::min(SIGMATRKZ[bestCand[0].convType],SIGMATRKCLUSZ[bestCand[0].convType]);                  
                  
                  convz = convz1;
                  convsz = convsz1;
                  nConv = 1;
                }
                if (hasconv2) {
                  convz2 = SIGMATRKZ[bestCand[1].convType] < SIGMATRKCLUSZ[bestCand[1].convType] ? bestCand[1].convTrkZ : bestCand[1].convTrkClusZ;
                  convsz2 = std::min(SIGMATRKZ[bestCand[1].convType],SIGMATRKCLUSZ[bestCand[1].convType]);
                  
                  convz = convz2;
                  convsz = convsz2;
                  nConv = 1;
                }
                if (hasconv1 && hasconv2) {
                  double w1 = 1./convsz1/convsz1;
                  double w2 = 1./convsz2/convsz2;
                  convz = (w1*convz1 + w2*convz2)/(w1+w2);
                  convsz = sqrt(1./(w1+w2));
                  nConv = 2;
                }
                
                if (hasconv1 || hasconv2) {
                  pull_conv = std::abs(convz-vtxZ)/convsz;
                  pull_conv = std::min(10.f,pull_conv); 
                }
                else {
                  pull_conv = 10.;
                }
       
	     
                // chi2 of photon-vertex timing
                chi2_pho_vtx_train = 0.0;
                logchi2_pho_vtx_train = 0.0;

		TVector3 pho1vtx(bestCand[0].scSeedX-vtxX,bestCand[0].scSeedY-vtxY,bestCand[0].scSeedZ-vtxZ);
                TVector3 pho2vtx(bestCand[1].scSeedX-vtxX,bestCand[1].scSeedY-vtxY,bestCand[1].scSeedZ-vtxZ);
              	
		TVector3 pho1_000(bestCand[0].scSeedX,bestCand[0].scSeedY,bestCand[0].scSeedZ);
		TVector3 pho2_000(bestCand[1].scSeedX,bestCand[1].scSeedY,bestCand[1].scSeedZ);
	
		double CmToNs = 0.1/2.99792458;

		double globalOffset = 0.0111;//11.54ps global offset from the result...

		float pho1Time_m_NoSmear = bestCand[0].scSeedT + CmToNs*pho1_000.Mag() - globalOffset;	
		float pho2Time_m_NoSmear = bestCand[1].scSeedT + CmToNs*pho2_000.Mag() - globalOffset;	
		
		double smear_phoT = 0.03;//30 ps smearing  
			
		float pho1Time_m = randomPhotonTime.Gaus(pho1Time_m_NoSmear, smear_phoT);
		float pho2Time_m = randomPhotonTime.Gaus(pho2Time_m_NoSmear, smear_phoT);
		//assuming vertex timing known
                {
		float pho1Time_e = vtxT + CmToNs*pho1vtx.Mag();
		float pho2Time_e = vtxT + CmToNs*pho2vtx.Mag();
		if(ipv<500)
		{
			pho1_vtx_dt_all[ipv] = pho1Time_m_NoSmear - pho1Time_e; 
			pho2_vtx_dt_all[ipv] = pho2Time_m_NoSmear - pho2Time_e; 

			chi2_pho_vtx_all[ipv] = pow((pho1Time_m - pho1Time_e)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e)/(smear_phoT), 2.0);	
		}
		
		if(useTiming)
		{
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e)/(smear_phoT), 2.0);
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
		}
		}
		{
		//assuming vertex timing not known
		double vtxT_min_chi2 = -0.5*(CmToNs*pho1vtx.Mag()+CmToNs*pho2vtx.Mag()-pho1Time_m-pho2Time_m);
		double pho1Time_e_tmp = vtxT_min_chi2 + CmToNs*pho1vtx.Mag();
		double pho2Time_e_tmp = vtxT_min_chi2 + CmToNs*pho2vtx.Mag();
	
		if(ipv<500)
		{
	        chi2_min_pho_vtx_all[ipv] = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);	
		}
		if(!useTiming)
		{
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
		}

		}

                float bdtval = vtxmvareader->EvaluateMVA("BDT");
		if(usePhoChi2) bdtval = vtxmvareader_chi2->EvaluateMVA("BDT");

		if (_debug) std::cout << "Vertex: " << ipv << " " << pvAllZ[ipv] << " : " 
				      << chi2_pho_vtx_train<<"   "<<ptasym << " " << ptbal << " " << logsumpt2 << " " << pull_conv << " " << nConv << " : "
				      << bdtval << "\n";
                              
                if (bdtval > maxbdtval) {
                  maxbdtval = bdtval;
                  ipvmax = ipv;
                }

		if(logsumpt2 > maxlogsumpt2)
                {
                  maxlogsumpt2 = logsumpt2;
                  ipvmaxlogsumpt2 = ipv;
                }

             	if(ipv<500)
		{ 
		ptasym_all[ipv]=ptasym;
		ptbal_all[ipv]=ptbal;
		logsumpt2_all[ipv]=logsumpt2;
		vertexsumpt_all[ipv]=vertexsumpt;
		vertexpt_all[ipv]=vertexpt;
		pull_conv_all[ipv]=pull_conv;
		nConv_all[ipv]=nConv;	  
		vtxdZCoordinate_all[ipv]=pvAllZ[ipv] - genVertexZ;	  
		vtxdXCoordinate_all[ipv]=pvAllX[ipv] - genVertexX;	  
		vtxdYCoordinate_all[ipv]=pvAllY[ipv] - genVertexY;	  
		vtxdTCoordinate_all[ipv]=pvAllT[ipv] - genVertexT;	  
		vtxZCoordinate_all[ipv]=pvAllZ[ipv];	  
		vtxXCoordinate_all[ipv]=pvAllX[ipv];	  
		vtxYCoordinate_all[ipv]=pvAllY[ipv];	  
		vtxTCoordinate_all[ipv]=pvAllT[ipv];	  
		if(useTiming) 
		{
		pvNtrack_all[ipv] = pvNtrack_orig_dt[ipv];
		//cout<<"pvNtrack "<<ipv<<"  "<< pvNtrack_orig_dt[ipv]<<endl;
		}
		else 
		{
		pvNtrack_all[ipv] = pvNtrack_orig[ipv];
		//cout<<"pvNtrack "<<ipv<<"  "<< pvNtrack_orig[ipv]<<endl;
		}
		}
              



	
		if(ipv==0)
		{
			TVector3 pho1vtx_gen(bestCand[0].scSeedX-genVertexX,bestCand[0].scSeedY-genVertexY,bestCand[0].scSeedZ-genVertexZ);
                	TVector3 pho2vtx_gen(bestCand[1].scSeedX-genVertexX,bestCand[1].scSeedY-genVertexY,bestCand[1].scSeedZ-genVertexZ);
              	
			float pho1Time_e_gen = genVertexT + CmToNs*pho1vtx_gen.Mag();
                	float pho2Time_e_gen = genVertexT + CmToNs*pho2vtx_gen.Mag();
			chi2_pho_vtx_gen = pow((pho1Time_m - pho1Time_e_gen)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_gen)/(smear_phoT), 2.0);
		

  			pho1_vtx_dt_gen = pho1Time_m_NoSmear - pho1Time_e_gen;
  			pho2_vtx_dt_gen = pho2Time_m_NoSmear - pho2Time_e_gen;
 	
			//min_chi2_tmp = std::numeric_limits<double>::max();
                	//vtxT_tmp = -2.0;
			min_vtxT_gen = -0.5*(CmToNs*pho1vtx_gen.Mag()+CmToNs*pho2vtx_gen.Mag()-pho1Time_m-pho2Time_m);
                	double pho1Time_e_tmp = min_vtxT_gen + CmToNs*pho1vtx_gen.Mag();
                	double pho2Time_e_tmp = min_vtxT_gen + CmToNs*pho2vtx_gen.Mag();
			chi2_min_pho_vtx_gen = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);
				
			/*
			for(int iStep=0;iStep<8000;iStep++)
			{	
			float pho1Time_e_tmp = vtxT_tmp + CmToNs*pho1vtx_gen.Mag();
			float pho2Time_e_tmp = vtxT_tmp + CmToNs*pho2vtx_gen.Mag();

			double chi2_pho_vtx_this = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);	
			if(chi2_pho_vtx_this < min_chi2_tmp) 
			{
			min_chi2_tmp = chi2_pho_vtx_this;	
			min_vtxT_gen = vtxT_tmp;	
			}
			
				vtxT_tmp = vtxT_tmp + 0.0005;//step width is 1 ps

			}
			chi2_min_pho_vtx_gen = min_chi2_tmp;
			*/
		}	
   		if(ipv==0)
		{
		//std::cout<<"entering z scan to find virtual photon..."<<std::endl;
	        //assuming vertex location and time not known, find the virtual vertex...
		//scan: vertex z: [-15,15]cm, with 0.1mm step, 3000 steps
		double min_deltaT_tmp = std::numeric_limits<double>::max();
		double virtual_vtx_x_this = beamSpotX;//-0.005;
		double virtual_vtx_y_this = beamSpotY;//-0.005;
		double virtual_vtx_z_this = -15.0;
		double virtual_vtx_z_selected = 0.0;
			
		for(int i_z=0;i_z<3000;i_z++){
			TVector3 pho1vtx_this(bestCand[0].scSeedX-virtual_vtx_x_this,bestCand[0].scSeedY-virtual_vtx_y_this,bestCand[0].scSeedZ-virtual_vtx_z_this);
			TVector3 pho2vtx_this(bestCand[1].scSeedX-virtual_vtx_x_this,bestCand[1].scSeedY-virtual_vtx_y_this,bestCand[1].scSeedZ-virtual_vtx_z_this);
			double deltaT_pho_e = CmToNs*pho1vtx_this.Mag() - CmToNs*pho2vtx_this.Mag();
			double deltaT_pho_m = pho1Time_m - pho2Time_m;	
			double deltaT_this = std::abs(deltaT_pho_e - deltaT_pho_m);	
			if(deltaT_this < min_deltaT_tmp)
			{
			  min_deltaT_tmp = deltaT_this;	
			  virtual_vtx_z_selected = virtual_vtx_z_this;	
			}
			virtual_vtx_z_this += 0.01;
		}
		
		//std::cout<<"found the virtual vertex(x, y, z, deltaT):"<<std::endl;
		//std::cout<<virtual_vtx_x_this<<"   "<<virtual_vtx_y_this<<"   "<<virtual_vtx_z_selected<<"  "<<min_deltaT_tmp<<std::endl;
		 
	        virtualVtxdZ =  virtual_vtx_z_selected - genVertexZ;//sqrt(pow(virtual_vtx_z_selected - genVertexZ, 2.0) + pow(virtual_vtx_x_this - genVertexX, 2.0) + pow(virtual_vtx_y_this - genVertexY, 2.0) );
		//std::cout<<"dz (virtual, gen): "<<virtualVtxdZ<<std::endl;
		
		}


         	}

	     //vertex chi2 rank
	     rank_chi2_genVtx = 1;
	     rank_chi2_min_genVtx = 1;
	     for(int ipv=0;ipv<nPVAll;ipv++)
	     {
		if(ipv<500)
		{
			if(chi2_pho_vtx_all[ipv]<chi2_pho_vtx_gen) rank_chi2_genVtx += 1;
			if(chi2_min_pho_vtx_all[ipv]<chi2_min_pho_vtx_gen) rank_chi2_min_genVtx += 1;
			for(int jpv=0;jpv<nPVAll;jpv++)
			{
				if(chi2_pho_vtx_all[jpv]<chi2_pho_vtx_all[ipv]) rank_chi2_pho_vtx_all[ipv] += 1;
				if(chi2_min_pho_vtx_all[jpv]<chi2_min_pho_vtx_all[ipv]) rank_chi2_min_pho_vtx_all[ipv] += 1;
			}
		}
	     }

	
	        if (_debug) std::cout << "gen-match Vertex: " << ipvmatch << "   dz: "<<sqrt(pow(pvAllZ[ipvmatch] - genVertexZ, 2.0)+pow(pvAllY[ipvmatch] - genVertexY, 2.0)+pow(pvAllX[ipvmatch] - genVertexX, 2.0))<<" ("<<(pow(pvAllZ[ipvmatch] - genVertexZ, 2.0)+pow(pvAllY[ipvmatch] - genVertexY, 2.0)+pow(pvAllX[ipvmatch] - genVertexX, 2.0))/(0.004*0.004*0.01*minDist)<<")  dT: "<<pvAllT[ipvmatch] - genVertexT<<"  ("<<pow(pvAllT[ipvmatch] - genVertexT,2.0)/(0.01*0.01*0.01*minDist)<<")  chi2_zt: "<< minDist<<"\n\n";


 		vtxIndex = ipvmax;

		if (_debug) std::cout << "maxbdt Vertex: " << ipvmax << "\n\n";
                if (_debug) std::cout << "gen-match Vertex: " << ipvmatch << "\n\n";
                if (_debug) std::cout << "maxlogsumpt2 Vertex: " << ipvmaxlogsumpt2 << "\n\n";

              	if (_debug) std::cout << "Vertex Selected: " << vtxIndex << "\n\n";

              	vtxZCoordinate =  pvAllZ[vtxIndex];
              	//vtxdZCoordinate =  pvAllZ[vtxIndex] - genVertexZ;
	        vtxdZCoordinate =  sqrt(pow(pvAllZ[vtxIndex] - genVertexZ, 2.0) + pow(pvAllX[vtxIndex] - genVertexX, 2.0) + pow(pvAllY[vtxIndex] - genVertexY, 2.0) );
              	if(useTiming)
                {
                vtxTCoordinate =  pvAllT[vtxIndex];
                vtxdTCoordinate =  pvAllT[vtxIndex] - genVertexT;
                }
	     
	      isMatchPv[ipvmatch]=1; 
	      isMaxbdtPv[ipvmax]=1; 
	      isMaxlogsumpt2Pv[ipvmaxlogsumpt2]=1; 
            

  
            } //end if do MVA Vertex
           
	//feed only the first N vertices (from chi2 rank) to the vertex selection algorithm
	
            if (doMVAVertex && FeedFirstN_only) {

              float pho1E = bestCand[0].photon.E();
              float pho2E = bestCand[1].photon.E();
             
  	      float maxbdtval = -99.;
              float maxlogsumpt2 = -99.;
              int ipvmax = 0;
              int ipvmatch = 0;
              int ipvmaxlogsumpt2 = 0;
	      
              double minDist = 999999.0;
 		 
              for (int ipv=0; ipv<nPVAll; ++ipv) {
		if(useTiming && rank_chi2_pho_vtx_all[ipv] >11 && abs(Pho_Eta[0]-Pho_Eta[1])>0.8) continue;
		if(useTiming && rank_chi2_pho_vtx_all[ipv] >18 && abs(Pho_Eta[0]-Pho_Eta[1])<0.8) continue;
		if(!useTiming && rank_chi2_pho_vtx_all[ipv] >41 && abs(Pho_Eta[0]-Pho_Eta[1])>0.8) continue;
		if(!useTiming && rank_chi2_pho_vtx_all[ipv] >43 && abs(Pho_Eta[0]-Pho_Eta[1])<0.8) continue;

                float vtxX = pvAllX[ipv];
                float vtxY = pvAllY[ipv];
                float vtxZ = pvAllZ[ipv];
                float vtxT = pvAllT[ipv];
	  	if(!useTiming)
		{	
		double dist_this = sqrt(pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0));
		if(dist_this < minDist)
		{
		    ipvmatch = ipv;
		    minDist = dist_this;
		}
		}
		else
		{
		  double chi2_zt = (pow(vtxZ - genVertexZ, 2.0)+pow(vtxY - genVertexY, 2.0)+pow(vtxX - genVertexX, 2.0))/(0.1*0.1)+ pow(vtxT - genVertexT, 2.0)/(0.06*0.06);
		if(chi2_zt < minDist)
		{
		     ipvmatch = ipv;
		     minDist = chi2_zt;
		}	
		}
                


                TVector3 pho1dir(bestCand[0].scX-vtxX,bestCand[0].scY-vtxY,bestCand[0].scZ-vtxZ);
                TVector3 pho2dir(bestCand[1].scX-vtxX,bestCand[1].scY-vtxY,bestCand[1].scZ-vtxZ);
                TVector3 diphomom = pho1E*pho1dir.Unit() + pho2E*pho2dir.Unit();
                TVector3 diphoPt(diphomom.x(), diphomom.y(), 0.);
                
                if(!useTiming)
		{ 
                TVector3 vtxSumPt(pvAllSumPx[ipv]-bestCand[0].vtxSumPx[ipv]-bestCand[1].vtxSumPx[ipv], pvAllSumPy[ipv]-bestCand[0].vtxSumPy[ipv]-bestCand[1].vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		}
		logsumpt2 = pvAllLogSumPtSq[ipv];
		vertexsumpt = pvAllSumPt[ipv];
		vertexpt = vtxSumPt.Mag();
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
                }
		
		else
		{
		TVector3 vtxSumPt(pvAllSumPx_dt[ipv]-bestCand[0].vtxSumPx[ipv]-bestCand[1].vtxSumPx[ipv], pvAllSumPy_dt[ipv]-bestCand[0].vtxSumPy[ipv]-bestCand[1].vtxSumPy[ipv],0.);
                

		if (_debug) {
		  std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			    << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			    << "\n";
		}
		logsumpt2 = pvAllLogSumPtSq_dt[ipv];
		vertexsumpt = pvAllSumPt_dt[ipv];
		//cout<<"vertexsumpt "<<ipv<<"  "<<pvAllSumPt[ipv]<<endl;	
		vertexpt = vtxSumPt.Mag();
                ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
                ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
		}

		diphotonpt = diphoPt.Mag();

                bool hasconv1 = bestCand[0].convType>=0;
                bool hasconv2 = bestCand[1].convType>=0;
                
                float convz1 = -99.;
                float convsz1 = -99.;
                float convz2 = -99.;
                float convsz2 = -99.;
                                
                nConv = 0.;
                float convz = -99.;
                float convsz = -99.;
                if (hasconv1) {
                  convz1 = SIGMATRKZ[bestCand[0].convType] < SIGMATRKCLUSZ[bestCand[0].convType] ? bestCand[0].convTrkZ : bestCand[0].convTrkClusZ;
                  convsz1 = std::min(SIGMATRKZ[bestCand[0].convType],SIGMATRKCLUSZ[bestCand[0].convType]);                  
                  
                  convz = convz1;
                  convsz = convsz1;
                  nConv = 1;
                }
                if (hasconv2) {
                  convz2 = SIGMATRKZ[bestCand[1].convType] < SIGMATRKCLUSZ[bestCand[1].convType] ? bestCand[1].convTrkZ : bestCand[1].convTrkClusZ;
                  convsz2 = std::min(SIGMATRKZ[bestCand[1].convType],SIGMATRKCLUSZ[bestCand[1].convType]);
                  
                  convz = convz2;
                  convsz = convsz2;
                  nConv = 1;
                }
                if (hasconv1 && hasconv2) {
                  double w1 = 1./convsz1/convsz1;
                  double w2 = 1./convsz2/convsz2;
                  convz = (w1*convz1 + w2*convz2)/(w1+w2);
                  convsz = sqrt(1./(w1+w2));
                  nConv = 2;
                }
                
                if (hasconv1 || hasconv2) {
                  pull_conv = std::abs(convz-vtxZ)/convsz;
                  pull_conv = std::min(10.f,pull_conv); 
                }
                else {
                  pull_conv = 10.;
                }
       
	     
                // chi2 of photon-vertex timing
                chi2_pho_vtx_train = 0.0;
                logchi2_pho_vtx_train = 0.0;

		TVector3 pho1vtx(bestCand[0].scSeedX-vtxX,bestCand[0].scSeedY-vtxY,bestCand[0].scSeedZ-vtxZ);
                TVector3 pho2vtx(bestCand[1].scSeedX-vtxX,bestCand[1].scSeedY-vtxY,bestCand[1].scSeedZ-vtxZ);
              	
		TVector3 pho1_000(bestCand[0].scSeedX,bestCand[0].scSeedY,bestCand[0].scSeedZ);
		TVector3 pho2_000(bestCand[1].scSeedX,bestCand[1].scSeedY,bestCand[1].scSeedZ);
	
		double CmToNs = 0.1/2.99792458;

		double globalOffset = 0.0111;//11.54ps global offset from the result...

		float pho1Time_m_NoSmear = bestCand[0].scSeedT + CmToNs*pho1_000.Mag() - globalOffset;	
		float pho2Time_m_NoSmear = bestCand[1].scSeedT + CmToNs*pho2_000.Mag() - globalOffset;	
		
		double smear_phoT = 0.03;//30 ps smearing  
			
		float pho1Time_m = randomPhotonTime.Gaus(pho1Time_m_NoSmear, smear_phoT);
		float pho2Time_m = randomPhotonTime.Gaus(pho2Time_m_NoSmear, smear_phoT);
		//assuming vertex timing known
                {
		float pho1Time_e = vtxT + CmToNs*pho1vtx.Mag();
		float pho2Time_e = vtxT + CmToNs*pho2vtx.Mag();
		if(ipv<500)
		{
			pho1_vtx_dt_all[ipv] = pho1Time_m_NoSmear - pho1Time_e; 
			pho2_vtx_dt_all[ipv] = pho2Time_m_NoSmear - pho2Time_e; 

			chi2_pho_vtx_all[ipv] = pow((pho1Time_m - pho1Time_e)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e)/(smear_phoT), 2.0);	
		}
		
		if(useTiming)
		{
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e)/(smear_phoT), 2.0);
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
		}
		}
		{
		//assuming vertex timing not known
		double vtxT_min_chi2 = -0.5*(CmToNs*pho1vtx.Mag()+CmToNs*pho2vtx.Mag()-pho1Time_m-pho2Time_m);
		double pho1Time_e_tmp = vtxT_min_chi2 + CmToNs*pho1vtx.Mag();
		double pho2Time_e_tmp = vtxT_min_chi2 + CmToNs*pho2vtx.Mag();
	
		if(ipv<500)
		{
	        chi2_min_pho_vtx_all[ipv] = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);	
		}
		if(!useTiming)
		{
		chi2_pho_vtx_train = pow((pho1Time_m - pho1Time_e_tmp)/(smear_phoT), 2.0) + pow((pho2Time_m - pho2Time_e_tmp)/(smear_phoT), 2.0);
		if(chi2_pho_vtx_train>0.0) logchi2_pho_vtx_train = log(chi2_pho_vtx_train);
		}

		}

                float bdtval = vtxmvareader->EvaluateMVA("BDT");
		if(usePhoChi2) bdtval = vtxmvareader_chi2->EvaluateMVA("BDT");

		if (_debug) std::cout << "Vertex: " << ipv << " " << pvAllZ[ipv] << " : " 
				      << chi2_pho_vtx_train<<"   "<<ptasym << " " << ptbal << " " << logsumpt2 << " " << pull_conv << " " << nConv << " : "
				      << bdtval << "\n";
                              
                if (bdtval > maxbdtval) {
                  maxbdtval = bdtval;
                  ipvmax = ipv;
                }

		if(logsumpt2 > maxlogsumpt2)
                {
                  maxlogsumpt2 = logsumpt2;
                  ipvmaxlogsumpt2 = ipv;
                }

         	}

 		vtxIndex = ipvmax;

		if (_debug) std::cout << "maxbdt Vertex: " << ipvmax << "\n\n";
                if (_debug) std::cout << "gen-match Vertex: " << ipvmatch << "\n\n";
                if (_debug) std::cout << "maxlogsumpt2 Vertex: " << ipvmaxlogsumpt2 << "\n\n";

              	if (_debug) std::cout << "Vertex Selected: " << vtxIndex << "\n\n";

              	vtxZCoordinate =  pvAllZ[vtxIndex];
              	//vtxdZCoordinate =  pvAllZ[vtxIndex] - genVertexZ;
	        vtxdZCoordinate =  sqrt(pow(pvAllZ[vtxIndex] - genVertexZ, 2.0) + pow(pvAllX[vtxIndex] - genVertexX, 2.0) + pow(pvAllY[vtxIndex] - genVertexY, 2.0) );
              	if(useTiming)
                {
                vtxTCoordinate =  pvAllT[vtxIndex];
                vtxdTCoordinate =  pvAllT[vtxIndex] - genVertexT;
                }
	     
	      isMatchPv[ipvmatch]=1; 
	      isMaxbdtPv[ipvmax]=1; 
	      isMaxlogsumpt2Pv[ipvmaxlogsumpt2]=1; 
  
            } //end if do MVA Vertex
    

 
    for(int i=0;i<allTrackdT->size();i++)
	{
		if(isMatchPv[(*allTrackPvIndex)[i]] == 0)
		{
			allTrackdT_Bkg->push_back((*allTrackdT)[i]);
			allTrackPt_Bkg->push_back((*allTrackPt)[i]);
			allTrackdZ_Bkg->push_back((*allTrackdZ)[i]);
		}
		else
		{
			allTrackdT_Sig->push_back((*allTrackdT)[i]);
			allTrackPt_Sig->push_back((*allTrackPt)[i]);
			allTrackdZ_Sig->push_back((*allTrackdZ)[i]);

		}	

		if((*allTrackPt)[i] > 10.0 )
		{

			hasLargeTrack[(*allTrackPvIndex)[i]] = 1;	
			if(std::abs((*allTrackdT)[i]) < 0.06) hasLargeTrack_dt[(*allTrackPvIndex)[i]] = 1;
		}

	}

     for(int i=0;i<nPVAll;i++)
	{
	   if(hasLargeTrack[i] == 1 && isMatchPv[i] ==0 ) NPV_HasLargeTrack += 1;
	   if(hasLargeTrack_dt[i] == 1 && isMatchPv[i] ==0 ) NPV_HasLargeTrack_dt += 1;
	   if(hasLargeTrack[i] == 1 && isMatchPv[i] ==1 ) NPV_GEN_HasLargeTrack += 1;
	   if(hasLargeTrack_dt[i] == 1 && isMatchPv[i] ==1 ) NPV_GEN_HasLargeTrack_dt += 1;
	}
 
    //---------------------------------------
    //just use this container for convenience
    //to parse the data into TTree
    //---------------------------------------
    phoSelectedCand.push_back(bestCand[0]);
    phoSelectedCand.push_back(bestCand[1]);
    
    if( _info )  std::cout<<"DEBUG_ZZC 013"<<std::endl;
    //Filling Selected Photon Information
    TLorentzVector pho_cand_vec[2];
    int _pho_index = 0;
    for ( auto& tmpPho : phoSelectedCand )
      {
	if ( !( tmpPho.Index == HiggsPhoIndex1 || tmpPho.Index == HiggsPhoIndex2 ) ) continue;
	if( _pho_index > 1 ) std::cerr << "[ERROR]: Photon index larger than 1!" << std::endl;
	pho_cand_vec[_pho_index]           = tmpPho.photon;
	Pho_E[_pho_index]                  = tmpPho.photon.E();
	Pho_Pt[_pho_index]                 = tmpPho.photon.Pt();
	Pho_Time[_pho_index]                = tmpPho.scSeedT;
	Pho_Eta[_pho_index]                = tmpPho.photon.Eta();
	Pho_Phi[_pho_index]                = tmpPho.photon.Phi();
	PhoSC_E[_pho_index]                = tmpPho.photonSC.E();
	PhoSC_Pt[_pho_index]               = tmpPho.photonSC.Pt();
	PhoSC_Eta[_pho_index]              = tmpPho.photonSC.Eta();
	PhoSC_Phi[_pho_index]              = tmpPho.photonSC.Phi();
	PhoDefaultSC_E[_pho_index]         = tmpPho.scE;
	PhoDefaultSC_Pt[_pho_index]        = tmpPho.scPt;
	PhoDefaultSC_Eta[_pho_index]       = tmpPho.scEta;
	PhoDefaultSC_Phi[_pho_index]       = tmpPho.scPhi;
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

    //if the best candidate pair has pT < 20 GeV, reject the event
    if( HiggsCandidate.Pt() < 20.0 )
      {
	if ( _debug ) std::cout << "[DEBUG]: Higgs Pt < 20 GeV, H pt: " << HiggsCandidate.Pt() << std::endl; 
	//continue;//apply offline
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
	for ( auto& phoC : phoSelectedCand )
	  {
	    if ( _debug ) std::cout << "===> phopt: " << phoC.photon.Pt() << " phoEta: " << phoC.photon.Eta() << std::endl;
	  }
	//continue;//apply offline
      }
    
    //record higgs candidate info
    mGammaGamma     = HiggsCandidate.M();
    pTGammaGamma    = HiggsCandidate.Pt();
    etaGammaGamma   = HiggsCandidate.Eta();
    phiGammaGamma   = HiggsCandidate.Phi();
    mGammaGammaSC   = HiggsCandidateSC.M();
    pTGammaGammaSC  = HiggsCandidateSC.Pt();
    etaGammaGammaSC = HiggsCandidateSC.Eta();
    phiGammaGammaSC = HiggsCandidateSC.Phi();

    //cosTheta*CS
    float pho1Z = Pho_Pt[0]*sinh( Pho_Eta[0] );
    float pho2Z = Pho_Pt[1]*sinh( Pho_Eta[1] );
    
    float pplus      = Pho_E[0] + pho1Z;
    float pminus     = Pho_E[0] - pho1Z;
    float antipplus  = Pho_E[1] + pho2Z;
    float antipminus = Pho_E[1] - pho2Z;
    
    cosThetaCS = signbit( HiggsCandidate.Pz() )*(pplus*antipminus-pminus*antipplus)/( mGammaGamma*sqrt( mGammaGamma*mGammaGamma + pTGammaGamma*pTGammaGamma ) );
    if ( _debug ) std::cout << "[DEBUG]: mgg-> " << mGammaGamma << " pTgg->" << pTGammaGamma << std::endl;
    

    if( _info )  std::cout<<"DEBUG_ZZC 014"<<std::endl;
    //***********************************************************
    //get mother ID of photons
    //***********************************************************
    // cout << "Photon1 : " << Pho_Pt[0] << " " << Pho_Eta[0] << " " << Pho_Phi[0] << "\n";
    for(int g = 0; g < nGenParticle; g++){
      if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[0],Pho_Phi[0]) < 0.1) ) continue;
      if(gParticleStatus[g] != 1) continue;
      if(gParticleId[g] != 22) continue;
      Pho_motherID[0] = gParticleMotherId[g];
      //cout << "Nearby GenParticle: " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherId[g] << "\n";
    }

    // cout << "Photon2 : " << Pho_Pt[1] << " " << Pho_Eta[1] << " " << Pho_Phi[1] << "\n";
    for(int g = 0; g < nGenParticle; g++){
      if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[1],Pho_Phi[1]) < 0.1) ) continue;
      if(gParticleStatus[g] != 1) continue;
      if(gParticleId[g] != 22) continue;
      Pho_motherID[1] = gParticleMotherId[g];      
      //cout << "Nearby GenParticle: " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherId[g] << "\n";
    }

    // cout << "\nGenParticles:\n";
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
    vector<TLorentzVector> GoodJetsJESUp;
    vector<TLorentzVector> GoodJetsJESDown;
    vector< pair<TLorentzVector, bool> > GoodCSVLJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not
    
    for( int i = 0; i < nJets; i++ )
      {
      if( _info )  std::cout<<"DEBUG_ZZC 014.01 #Jets "<<i<<"   "<<nJets<<std::endl;
	//Jet Corrections                                                                      
	double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						fixedGridRhoAll, jetJetArea[i],
						JetCorrector );
	
	TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
	
	if( thisJet.Pt() < 30.0 ) continue;//According to the April 1st 2015 AN
	if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
	//int level = 2; //loose jet ID
	if ( !jetPassIDLoose[i] ) continue;
	//if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;
	
	//exclude selected photons from the jet collection
	double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
	if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN
	
	GoodJets.push_back(thisJet);
	n_Jets++;
	
	//-------------------------------
	//b - t a g   c o r r e c t i o n
	//-------------------------------
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
	

	//-----------------------------
	//J E S   U n c e r t a i n t y
	//-----------------------------
	if ( !isData )
	  {
	    jecUnc->setJetEta(jetEta[i]);
	    jecUnc->setJetPt(jetCorrPt);
	    double unc = jecUnc->getUncertainty(true);
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
	  }//isData
	
    if( _info )  std::cout<<"DEBUG_ZZC 015"<<std::endl;
	//------------------------------
	// b - t a g   s e l e c t i o n
	//------------------------------
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
	  }//isCVSL
      }
    
    //if there are no good jets, reject the event
    if( n_Jets == 0 )
      {
	if ( _debug ) std::cout << "[DEBUG]: No Jets Selected" << std::endl;
	//continue;//apply offline
      }

    
    //sort jets by PT
    auto sortTLV = []( TLorentzVector a, TLorentzVector b ){ return a.Pt() > b.Pt() ?  true : false; };
    std::sort( GoodJets.begin(), GoodJets.end(), sortTLV );
    
    if( _info )  std::cout<<"DEBUG_ZZC 016 GoodJets.size() = "<<GoodJets.size()<<std::endl;
    //dijet mass and gg+jjj mass
    if ( GoodJets.size() >= 2 )
      {
	mjj = (GoodJets[0] + GoodJets[1]).M();
	mggjj = (HiggsCandidate + GoodJets[0] + GoodJets[1]).M();
	deltaEtajj = GoodJets[0].Eta() - GoodJets[1].Eta();
      }

    int iJet = 0;
    for ( auto tmp_jet : GoodJets )
      {
	jet_E[iJet] = tmp_jet.E();
	jet_Pt[iJet] = tmp_jet.Pt();
	jet_Eta[iJet] = tmp_jet.Eta();
	jet_Phi[iJet] = tmp_jet.Phi();
	mggj[iJet]    = (HiggsCandidate + tmp_jet).M();
	iJet++;
      }
    
    if( _info )  std::cout<<"DEBUG_ZZC 016.01 GoodJets.size() = "<<GoodJets.size()<<std::endl;
    //Compute the razor variables using the selected jets and the diphoton system
    vector<TLorentzVector> JetsPlusHiggsCandidate;
    for( auto& jet : GoodJets ) JetsPlusHiggsCandidate.push_back(jet);
    JetsPlusHiggsCandidate.push_back(HiggsCandidate);
    
    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
    TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
    MET = metPt;
    t1MET = metType1Pt;
    
    if( _info )  std::cout<<"DEBUG_ZZC 016.02 JetsPlusHiggsCandidate.size() = "<<JetsPlusHiggsCandidate.size()<<std::endl;
    vector<TLorentzVector> hemispheres = getHemispheres(JetsPlusHiggsCandidate);
    if( _info )  std::cout<<"DEBUG_ZZC 016.03 GoodJets.size() = "<<GoodJets.size()<<std::endl;
    theMR  = computeMR(hemispheres[0], hemispheres[1]); 
    if( _info )  std::cout<<"DEBUG_ZZC 016.04 GoodJets.size() = "<<GoodJets.size()<<std::endl;
    if ( theMR > 0 )
      {
    if( _info )  std::cout<<"DEBUG_ZZC 016.05 GoodJets.size() = "<<GoodJets.size()<<std::endl;
	theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	t1Rsq  = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
      }
    
    if( _info )  std::cout<<"DEBUG_ZZC 017"<<std::endl;
    //-----------------------------------
    //Propagate JES to MET and MR and Rsq
    //-----------------------------------
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
      }//isData
    
    if( _info )  std::cout<<"DEBUG_ZZC 018"<<std::endl;
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
	//continue;//inclusive now, apply razor cuts offline
      }
    
    //if there are two loose b-tags and one medium b-tag, look for b-bbar resonances
    //--------------------------------------------
    //mbbH and mbbZ from 1Medium and 1Loose b-tags
    //-------------------------------------------- 
    if( nLooseBTaggedJets > 1 && nMediumBTaggedJets > 0 )
      {
	for( int i = 0; i < nLooseBTaggedJets; i++ )
	  {
	    for( int j = i+1; j < nLooseBTaggedJets; j++ )
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
    
    //-----------------------------------
    //mbbH and mbbZ from two Loose b-tags
    //-----------------------------------
    if( nLooseBTaggedJets >= 2 )//at least two btag jets
      {
        for( int i = 0; i < nLooseBTaggedJets; i++ )
          {
            for( int j = i+1; j < nLooseBTaggedJets; j++ )
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
	if(combineTrees){
	  if ( _debug ) std::cout << "[DEBUG]: combineTrees" << std::endl;
	  razorbox = HighPt;
	  if ( _debug ) std::cout << "[DEBUG]: combineTrees 1" << std::endl;
	  razorTree->Fill();
	  if ( _debug ) std::cout << "[DEBUG]: combineTrees 2" << std::endl;
	}
	else razorBoxes["HighPt"]->Fill();
      }
    //Hbb Box
    else if ( mbbH > 110.0 && mbbH < 140.0 )
      {
	if(combineTrees){
	  razorbox = Hbb;
	  razorTree->Fill();
	}
	else razorBoxes["Hbb"]->Fill();
      }
    //Zbb Box
    else if( mbbZ > 76.0 && mbbZ < 106.0 )
      {
	if(combineTrees){
	  razorbox = Zbb;
	  razorTree->Fill();
	}
	else razorBoxes["Zbb"]->Fill();
      }
    //HighRes Box
    else if( Pho_sigmaEOverE[0] < 0.015 && Pho_sigmaEOverE[1] < 0.015 )
      {
	if(combineTrees){
	  razorbox = HighRes;
	  razorTree->Fill();
	}
	else razorBoxes["HighRes"]->Fill();
      }
    //LowRes Box
    else
      {
	if(combineTrees){
	  razorbox = LowRes;
	  razorTree->Fill();
	}
	else razorBoxes["LowRes"]->Fill();
      }
    if ( _debug ) std::cout << "[DEBUG]: combineTrees 3?" << std::endl;
    //end of event loop
  }
  
  if ( _info ) std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;
  if ( _info ) std::cout << "[INFO]: Writing output trees..." << std::endl;

  //----------------------------
  //Writing output and finishing
  //----------------------------
  outFile->cd();
  if(combineTrees) razorTree->Write();
  else for(auto& box : razorBoxes) box.second->Write();
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  puhisto->Write();
  outFile->Close();

  delete photonCorrector;

}
