#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>
#include <sstream>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

const int NUM_PDF_WEIGHTS = 60;

void RazorAnalyzer::FullRazorInclusive(string outFileName, bool isData, bool isFastsimSMS)
{
    /////////////////////////////////
    //Basic setup
    /////////////////////////////////

    cout << "Initializing..." << endl;
    TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    //Output file
    if (outFileName.empty()){
        cout << "FullRazorInclusive: Output filename not specified!" << endl << "Using default output name FullRazorInclusive.root" << endl;
        outFileName = "FullRazorInclusive.root";
    }
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

    //Output tree
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

    //For signal samples, create one output file and tree per signal mass point
    map<pair<int,int>, TFile*> smsFiles;
    map<pair<int,int>, TTree*> smsTrees;
    map<pair<int,int>, TH1F*> smsNEvents;
    map<pair<int,int>, TH1F*> smsSumWeights;
    map<pair<int,int>, TH1F*> smsSumScaleWeights;
    //map<pair<int,int>, TH1F*> smsSumPdfWeights;

    //Histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
    TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
    TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
    //TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

    char* cmsswPath;
    cmsswPath = getenv("CMSSW_BASE");
    if (cmsswPath == NULL) {
        cout << "Warning: CMSSW_BASE not detected. Exiting..." << endl;
        return;
    }

    /////////////////////////////////
    //Pileup Weights
    /////////////////////////////////

    TFile *pileupWeightFile = 0;
    TH1D *pileupWeightHist = 0;
    if(!isData){
        pileupWeightFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors//PileupWeights/NVtxReweight_ZToMuMu_2015D_1264ipb.root");
        pileupWeightHist = (TH1D*)pileupWeightFile->Get("NVtxReweight");
        assert(pileupWeightHist);
    }

    /////////////////////////////////
    //Lepton and b-tag Efficiency 
    /////////////////////////////////

    TH2D *eleTightEfficiencyHist = 0;
    TH2D *muTightEfficiencyHist = 0;
    TH2D *eleVetoEfficiencyHist = 0;
    TH2D *muVetoEfficiencyHist = 0;
    TH2D *tauLooseEfficiencyHist = 0;
    TH2D *btagMediumEfficiencyHist = 0;
    TH2D *eleTightEffFastsimSFHist = 0;
    TH2D *muTightEffFastsimSFHist = 0;
    TH2D *btagMediumEffFastsimSFHist = 0;
    if(!isData){
        TFile *eleEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.root");
        eleTightEfficiencyHist = (TH2D*)eleEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
        eleVetoEfficiencyHist = (TH2D*)eleEfficiencyFile->Get("ElectronEff_Veto_Fullsim");
        assert(eleTightEfficiencyHist);
        eleTightEffFastsimSFHist = (TH2D*)eleEfficiencyFile->Get("ElectronTight_FastsimScaleFactor");
        assert(eleTightEffFastsimSFHist);
        TFile *muEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.root");
        muTightEfficiencyHist = (TH2D*)muEfficiencyFile->Get("MuonEff_Tight_Fullsim");
        muVetoEfficiencyHist = (TH2D*)muEfficiencyFile->Get("MuonEff_Veto_Fullsim");
        assert(muTightEfficiencyHist);
        muTightEffFastsimSFHist = (TH2D*)muEfficiencyFile->Get("MuonTight_FastsimScaleFactor");
        assert(muTightEffFastsimSFHist);
        TFile *tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
        tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");
        assert(tauLooseEfficiencyHist);
        TFile *btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/BTagEffFastsimToFullsimCorrectionFactors.root");
        btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("BTagEff_Medium_Fullsim");
        assert(btagMediumEfficiencyHist);
        btagMediumEffFastsimSFHist = (TH2D*)btagEfficiencyFile->Get("BTagMedium_FastsimScaleFactor");
    }

    /////////////////////////////////
    //Lepton Efficiency Correction Factors
    /////////////////////////////////

    TH2D *eleTightEffSFHist = 0;
    TH2D *muTightEffSFHist = 0;
    TH2D *eleVetoEffSFHist = 0;
    TH2D *muVetoEffSFHist = 0;
    if(!isData){
        TFile *eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root");
        eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
        assert(eleTightEffSFHist);
        TFile *muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root"); 
        muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
        assert(muTightEffSFHist);
        TFile *vetoEleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_VetoElectronSelectionEffDenominatorReco_2015D_Golden.root");
        eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
        assert(eleVetoEffSFHist);
        TFile *vetoMuEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_VetoMuonSelectionEffDenominatorReco_2015D_Golden.root"); 
        muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
        assert(muVetoEffSFHist);
    }

    /////////////////////////////////
    //Trigger Efficiency Correction Factors
    /////////////////////////////////

    TH2D *eleTrigSFHist = 0;
    TH2D *muTrigSFHist = 0;
    if(!isData){
        TFile *eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root");
        eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
        assert(eleTrigSFHist);
        TFile *muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root"); 
        muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");
        assert(muTrigSFHist);
    }

    TH2D *eleTrigEffFromFullsimHist = 0;
    TH2D *muTrigEffFromFullsimHist = 0;
    if(isFastsimSMS){
        TFile *eleTrigEffFromFullsimFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleElectronTriggerEfficiencyFromFullsim.root");
        eleTrigEffFromFullsimHist = (TH2D*)eleTrigEffFromFullsimFile->Get("hEffEtaPt");
        assert(eleTrigEffFromFullsimHist);
        TFile *muTrigEffFromFullsimFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleMuonTriggerEfficiencyFromFullsim.root"); 
        muTrigEffFromFullsimHist = (TH2D*)muTrigEffFromFullsimFile->Get("hEffEtaPt");
        assert(muTrigEffFromFullsimHist);
    }
    
    /////////////////////////////////
    //Jet Energy Corrections
    /////////////////////////////////

    //Get directory for JEC files 
    string pathname;
    if (cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
    else {
        cout << "Warning: CMSSW_BASE not detected.  Looking for JEC parameters in data/JEC";
        pathname = "data/JEC/";
    }

    //Get JEC files
    cout << "Getting JEC parameters from " << pathname << endl;
    std::vector<JetCorrectorParameters> correctionParameters;
    if (isData) {
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
    } 
    else {
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    }
    //Set up JEC machinery 
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
    JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
    SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);

    //Get JEC uncertainty file and set up JetCorrectionUncertainty
    string jecUncPath;
    if (isData) {
        jecUncPath = pathname+"/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt";
    }
    else {
        jecUncPath = pathname+"/Summer15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncPath);

    /////////////////////////////////
    //B-tagging scale factors
    /////////////////////////////////

    string bTagPathname = "";
    if (cmsswPath != NULL) bTagPathname = string(cmsswPath) + "/src/RazorAnalyzer/data/ScaleFactors/";
    else bTagPathname = "data/ScaleFactors/";
    BTagCalibration btagcalib("csvv2", Form("%s/CSVv2.csv",bTagPathname.c_str()));
    BTagCalibrationReader btagreader(&btagcalib,               // calibration instance
            BTagEntry::OP_MEDIUM,  // operating point
            "mujets",               // measurement type
            "central");           // systematics type
    BTagCalibrationReader btagreader_up(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    BTagCalibrationReader btagreader_do(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down

    /////////////////////////////////
    //Tree Initialization
    /////////////////////////////////

    //Basic tree variables
    int nVtx, nSelectedJets, nBTaggedJets, nJets80; 
    float dPhiRazor, MR, Rsq, mT, mTLoose, 
          leadingJetPt, subleadingJetPt, leadingTightMuPt, leadingTightElePt;
    float weight = 1.0;
    RazorBox box;
    //For lepton efficiency scale factor uncertainty
    float sf_muonEffUp, sf_muonEffDown;
    float sf_eleEffUp, sf_eleEffDown;
    float sf_muonTrigUp, sf_muonTrigDown;
    float sf_eleTrigUp, sf_eleTrigDown;
    float sf_tauEffUp, sf_tauEffDown;
    float sf_vetoMuonEffUp, sf_vetoMuonEffDown;
    float sf_vetoEleEffUp, sf_vetoEleEffDown;
    //For btag scale factor uncertainty
    float sf_btagUp, sf_btagDown;
    //For scale variation uncertainties
    float sf_facScaleUp, sf_facScaleDown;
    float sf_renScaleUp, sf_renScaleDown;
    float sf_facRenScaleUp, sf_facRenScaleDown;
    //For Fastsim scale factor uncertainties
    float sf_muonEffFastsimSFUp, sf_muonEffFastsimSFDown;
    float sf_eleEffFastsimSFUp, sf_eleEffFastsimSFDown;
    float sf_btagFastsimSFUp, sf_btagFastsimSFDown;
    //For jet uncertainties
    float MR_JESUp, Rsq_JESUp, dPhiRazor_JESUp, leadingJetPt_JESUp, subleadingJetPt_JESUp; 
    float MR_JESDown, Rsq_JESDown, dPhiRazor_JESDown, leadingJetPt_JESDown, subleadingJetPt_JESDown;
    float MR_MESUp, Rsq_MESUp, dPhiRazor_MESUp, leadingTightMuPt_MESUp; 
    float MR_MESDown, Rsq_MESDown, dPhiRazor_MESDown, leadingTightMuPt_MESDown;
    float MR_EESUp, Rsq_EESUp, dPhiRazor_EESUp, leadingTightElePt_EESUp; 
    float MR_EESDown, Rsq_EESDown, dPhiRazor_EESDown, leadingTightElePt_EESDown;
    float MR_JERUp, Rsq_JERUp, dPhiRazor_JERUp, leadingJetPt_JERUp, subleadingJetPt_JERUp;
    float MR_JERDown, Rsq_JERDown, dPhiRazor_JERDown, leadingJetPt_JERDown, subleadingJetPt_JERDown;
    int nSelectedJets_JESUp, nSelectedJets_JESDown, nSelectedJets_JERUp, nSelectedJets_JERDown;
    int nBTaggedJets_JESUp, nBTaggedJets_JESDown, nBTaggedJets_JERUp, nBTaggedJets_JERDown;
    int nSelectedJets_MESUp, nSelectedJets_MESDown, nSelectedJets_EESUp, nSelectedJets_EESDown;
    int nBTaggedJets_MESUp, nBTaggedJets_MESDown, nBTaggedJets_EESUp, nBTaggedJets_EESDown;
    int nJets80_JESUp, nJets80_JESDown, nJets80_JERUp, nJets80_JERDown;
    int nJets80_MESUp, nJets80_MESDown, nJets80_EESUp, nJets80_EESDown;
    RazorBox box_JESUp, box_JESDown, box_JERUp, box_JERDown;
    RazorBox box_MESUp, box_MESDown, box_EESUp, box_EESDown;
    float mT_JESUp, mT_JESDown, mT_JERUp, mT_JERDown;
    float mTLoose_JESUp, mTLoose_JESDown, mTLoose_JERUp, mTLoose_JERDown;
    float mT_MESUp, mT_MESDown, mT_EESUp, mT_EESDown;
    float mTLoose_MESUp, mTLoose_MESDown, mTLoose_EESUp, mTLoose_EESDown;
    float leadingJetPt_MESUp, leadingJetPt_MESDown, leadingJetPt_EESUp, leadingJetPt_EESDown;
    float subleadingJetPt_MESUp, subleadingJetPt_MESDown, subleadingJetPt_EESUp, subleadingJetPt_EESDown;
    int nVetoMuons, nTightMuons, nVetoElectrons, nTightElectrons, nLooseTaus;
    float met, HT;
    //SMS parameters 
    int mGluino, mLSP;

    //Set branches
    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nJets80", &nJets80, "nJets80/I");
    razorTree->Branch("MR", &MR, "MR/F");
    razorTree->Branch("Rsq", &Rsq, "Rsq/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("mT", &mT, "mT/F");
    razorTree->Branch("mTLoose", &mTLoose, "mTLoose/F");//for LooseLepton boxes
    razorTree->Branch("leadingTightMuPt", &leadingTightMuPt, "leadingTightMuPt/F");
    razorTree->Branch("leadingTightElePt", &leadingTightElePt, "leadingTightElePt/F");
    razorTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
    razorTree->Branch("subleadingJetPt", &subleadingJetPt, "subleadingJetPt/F");
    razorTree->Branch("box", &box, "box/I");
    razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
    razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("HLTDecision", &HLTDecision, "HLTDecision[150]/O");
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

    if (!isData) {    
        razorTree->Branch("genWeight", &genWeight, "genWeight/F");
        razorTree->Branch("weight", &weight, "weight/F");
        razorTree->Branch("sf_muonEffUp", &sf_muonEffUp, "sf_muonEffUp/F");
        razorTree->Branch("sf_muonEffDown", &sf_muonEffDown, "sf_muonEffDown/F");
        razorTree->Branch("sf_vetoMuonEffUp", &sf_vetoMuonEffUp, "sf_vetoMuonEffUp/F");
        razorTree->Branch("sf_vetoMuonEffDown", &sf_vetoMuonEffDown, "sf_vetoMuonEffDown/F");
        razorTree->Branch("sf_eleEffUp", &sf_eleEffUp, "sf_eleEffUp/F");
        razorTree->Branch("sf_eleEffDown", &sf_eleEffDown, "sf_eleEffDown/F");
        razorTree->Branch("sf_vetoEleEffUp", &sf_vetoEleEffUp, "sf_vetoEleEffUp/F");
        razorTree->Branch("sf_vetoEleEffDown", &sf_vetoEleEffDown, "sf_vetoEleEffDown/F");
        razorTree->Branch("sf_tauEffUp", &sf_tauEffUp, "sf_tauEffUp/F");
        razorTree->Branch("sf_tauEffDown", &sf_tauEffDown, "sf_tauEffDown/F");
        razorTree->Branch("sf_muonTrigUp", &sf_muonTrigUp, "sf_muonTrigUp/F");
        razorTree->Branch("sf_muonTrigDown", &sf_muonTrigDown, "sf_muonTrigDown/F");
        razorTree->Branch("sf_eleTrigUp", &sf_eleTrigUp, "sf_eleTrigUp/F");
        razorTree->Branch("sf_eleTrigDown", &sf_eleTrigDown, "sf_eleTrigDown/F");
        razorTree->Branch("sf_btagUp", &sf_btagUp, "sf_btagUp/F");
        razorTree->Branch("sf_btagDown", &sf_btagDown, "sf_btagDown/F");
        razorTree->Branch("sf_muonEffFastsimSFUp", &sf_muonEffFastsimSFUp, "sf_muonEffFastsimSFUp/F");
        razorTree->Branch("sf_muonEffFastsimSFDown", &sf_muonEffFastsimSFDown, "sf_muonEffFastsimSFDown/F");
        razorTree->Branch("sf_eleEffFastsimSFUp", &sf_eleEffFastsimSFUp, "sf_eleEffFastsimSFUp/F");
        razorTree->Branch("sf_eleEffFastsimSFDown", &sf_eleEffFastsimSFDown, "sf_eleEffFastsimSFDown/F");
        razorTree->Branch("sf_btagFastsimSFUp", &sf_btagFastsimSFUp, "sf_btagFastsimSFUp/F");
        razorTree->Branch("sf_btagFastsimSFDown", &sf_btagFastsimSFDown, "sf_btagFastsimSFDown/F");
        razorTree->Branch("sf_facScaleUp", &sf_facScaleUp, "sf_facScaleUp/F");
        razorTree->Branch("sf_facScaleDown", &sf_facScaleDown, "sf_facScaleDown/F");
        razorTree->Branch("sf_renScaleUp", &sf_renScaleUp, "sf_renScaleUp/F");
        razorTree->Branch("sf_renScaleDown", &sf_renScaleDown, "sf_renScaleDown/F");
        razorTree->Branch("sf_facRenScaleUp", &sf_facRenScaleUp, "sf_facRenScaleUp/F");
        razorTree->Branch("sf_facRenScaleDown", &sf_facRenScaleDown, "sf_facRenScaleDown/F");
        //razorTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights); //get PDF weights directly from RazorEvents
        razorTree->Branch("MR_JESUp", &MR_JESUp, "MR_JESUp/F");
        razorTree->Branch("Rsq_JESUp", &Rsq_JESUp, "Rsq_JESUp/F");
        razorTree->Branch("dPhiRazor_JESUp", &dPhiRazor_JESUp, "dPhiRazor_JESUp/F");
        razorTree->Branch("leadingJetPt_JESUp", &leadingJetPt_JESUp, "leadingJetPt_JESUp/F");
        razorTree->Branch("subleadingJetPt_JESUp", &subleadingJetPt_JESUp, "subleadingJetPt_JESUp/F");
        razorTree->Branch("nSelectedJets_JESUp", &nSelectedJets_JESUp, "nSelectedJets_JESUp/I");
        razorTree->Branch("nBTaggedJets_JESUp", &nBTaggedJets_JESUp, "nBTaggedJets_JESUp/I");
        razorTree->Branch("nJets80_JESUp", &nJets80_JESUp, "nJets80_JESUp/I");
        razorTree->Branch("mT_JESUp", &mT_JESUp, "mT_JESUp/F");
        razorTree->Branch("mTLoose_JESUp", &mTLoose_JESUp, "mTLoose_JESUp/F");
        razorTree->Branch("box_JESUp", &box_JESUp, "box_JESUp/I");
        razorTree->Branch("MR_JESDown", &MR_JESDown, "MR_JESDown/F");
        razorTree->Branch("Rsq_JESDown", &Rsq_JESDown, "Rsq_JESDown/F");
        razorTree->Branch("dPhiRazor_JESDown", &dPhiRazor_JESDown, "dPhiRazor_JESDown/F");
        razorTree->Branch("leadingJetPt_JESDown", &leadingJetPt_JESDown, "leadingJetPt_JESDown/F");
        razorTree->Branch("subleadingJetPt_JESDown", &subleadingJetPt_JESDown, "subleadingJetPt_JESDown/F");
        razorTree->Branch("nSelectedJets_JESDown", &nSelectedJets_JESDown, "nSelectedJets_JESDown/I");
        razorTree->Branch("nBTaggedJets_JESDown", &nBTaggedJets_JESDown, "nBTaggedJets_JESDown/I");
        razorTree->Branch("nJets80_JESDown", &nJets80_JESDown, "nJets80_JESDown/I");
        razorTree->Branch("mT_JESDown", &mT_JESDown, "mT_JESDown/F");
        razorTree->Branch("mTLoose_JESDown", &mTLoose_JESDown, "mTLoose_JESDown/F");
        razorTree->Branch("box_JESDown", &box_JESDown, "box_JESDown/I");
        razorTree->Branch("MR_EESUp", &MR_EESUp, "MR_EESUp/F");
        razorTree->Branch("Rsq_EESUp", &Rsq_EESUp, "Rsq_EESUp/F");
        razorTree->Branch("dPhiRazor_EESUp", &dPhiRazor_EESUp, "dPhiRazor_EESUp/F");
        razorTree->Branch("leadingTightElePt_EESUp", &leadingTightElePt_EESUp, "leadingTightElePt_EESUp/F");
        razorTree->Branch("nSelectedJets_EESUp", &nSelectedJets_EESUp, "nSelectedJets_EESUp/I");
        razorTree->Branch("nBTaggedJets_EESUp", &nBTaggedJets_EESUp, "nBTaggedJets_EESUp/I");
        razorTree->Branch("nJets80_EESUp", &nJets80_EESUp, "nJets80_EESUp/I");
        razorTree->Branch("mT_EESUp", &mT_EESUp, "mT_EESUp/F");
        razorTree->Branch("mTLoose_EESUp", &mTLoose_EESUp, "mTLoose_EESUp/F");
        razorTree->Branch("leadingJetPt_EESUp", &leadingJetPt_EESUp, "leadingJetPt_EESUp/F");
        razorTree->Branch("subleadingJetPt_EESUp", &subleadingJetPt_EESUp, "subleadingJetPt_EESUp/F");
        razorTree->Branch("box_EESUp", &box_EESUp, "box_EESUp/I");
        razorTree->Branch("MR_EESDown", &MR_EESDown, "MR_EESDown/F");
        razorTree->Branch("Rsq_EESDown", &Rsq_EESDown, "Rsq_EESDown/F");
        razorTree->Branch("dPhiRazor_EESDown", &dPhiRazor_EESDown, "dPhiRazor_EESDown/F");
        razorTree->Branch("leadingTightElePt_EESDown", &leadingTightElePt_EESDown, "leadingTightElePt_EESDown/F");
        razorTree->Branch("nSelectedJets_EESDown", &nSelectedJets_EESDown, "nSelectedJets_EESDown/I");
        razorTree->Branch("nBTaggedJets_EESDown", &nBTaggedJets_EESDown, "nBTaggedJets_EESDown/I");
        razorTree->Branch("nJets80_EESDown", &nJets80_EESDown, "nJets80_EESDown/I");
        razorTree->Branch("mT_EESDown", &mT_EESDown, "mT_EESDown/F");
        razorTree->Branch("mTLoose_EESDown", &mTLoose_EESDown, "mTLoose_EESDown/F");
        razorTree->Branch("leadingJetPt_EESDown", &leadingJetPt_EESDown, "leadingJetPt_EESDown/F");
        razorTree->Branch("subleadingJetPt_EESDown", &subleadingJetPt_EESDown, "subleadingJetPt_EESDown/F");
        razorTree->Branch("box_EESDown", &box_EESDown, "box_EESDown/I");
        razorTree->Branch("MR_MESUp", &MR_MESUp, "MR_MESUp/F");
        razorTree->Branch("Rsq_MESUp", &Rsq_MESUp, "Rsq_MESUp/F");
        razorTree->Branch("dPhiRazor_MESUp", &dPhiRazor_MESUp, "dPhiRazor_MESUp/F");
        razorTree->Branch("leadingTightMuPt_MESUp", &leadingTightMuPt_MESUp, "leadingTightMuPt_MESUp/F");
        razorTree->Branch("nSelectedJets_MESUp", &nSelectedJets_MESUp, "nSelectedJets_MESUp/I");
        razorTree->Branch("nBTaggedJets_MESUp", &nBTaggedJets_MESUp, "nBTaggedJets_MESUp/I");
        razorTree->Branch("nJets80_MESUp", &nJets80_MESUp, "nJets80_MESUp/I");
        razorTree->Branch("mT_MESUp", &mT_MESUp, "mT_MESUp/F");
        razorTree->Branch("mTLoose_MESUp", &mTLoose_MESUp, "mTLoose_MESUp/F");
        razorTree->Branch("leadingJetPt_MESUp", &leadingJetPt_MESUp, "leadingJetPt_MESUp/F");
        razorTree->Branch("subleadingJetPt_MESUp", &subleadingJetPt_MESUp, "subleadingJetPt_MESUp/F");
        razorTree->Branch("box_MESUp", &box_MESUp, "box_MESUp/I");
        razorTree->Branch("MR_MESDown", &MR_MESDown, "MR_MESDown/F");
        razorTree->Branch("Rsq_MESDown", &Rsq_MESDown, "Rsq_MESDown/F");
        razorTree->Branch("dPhiRazor_MESDown", &dPhiRazor_MESDown, "dPhiRazor_MESDown/F");
        razorTree->Branch("leadingTightMuPt_MESDown", &leadingTightMuPt_MESDown, "leadingTightMuPt_MESDown/F");
        razorTree->Branch("nSelectedJets_MESDown", &nSelectedJets_MESDown, "nSelectedJets_MESDown/I");
        razorTree->Branch("nBTaggedJets_MESDown", &nBTaggedJets_MESDown, "nBTaggedJets_MESDown/I");
        razorTree->Branch("nJets80_MESDown", &nJets80_MESDown, "nJets80_MESDown/I");
        razorTree->Branch("mT_MESDown", &mT_MESDown, "mT_MESDown/F");
        razorTree->Branch("mTLoose_MESDown", &mTLoose_MESDown, "mTLoose_MESDown/F");
        razorTree->Branch("leadingJetPt_MESDown", &leadingJetPt_MESDown, "leadingJetPt_MESDown/F");
        razorTree->Branch("subleadingJetPt_MESDown", &subleadingJetPt_MESDown, "subleadingJetPt_MESDown/F");
        razorTree->Branch("box_MESDown", &box_MESDown, "box_MESDown/I");
        razorTree->Branch("MR_JERUp", &MR_JERUp, "MR_JERUp/F");
        razorTree->Branch("Rsq_JERUp", &Rsq_JERUp, "Rsq_JERUp/F");
        razorTree->Branch("dPhiRazor_JERUp", &dPhiRazor_JERUp, "dPhiRazor_JERUp/F");
        razorTree->Branch("leadingJetPt_JERUp", &leadingJetPt_JERUp, "leadingJetPt_JERUp/F");
        razorTree->Branch("subleadingJetPt_JERUp", &subleadingJetPt_JERUp, "subleadingJetPt_JERUp/F");
        razorTree->Branch("nSelectedJets_JERUp", &nSelectedJets_JERUp, "nSelectedJets_JERUp/I");
        razorTree->Branch("nBTaggedJets_JERUp", &nBTaggedJets_JERUp, "nBTaggedJets_JERUp/I");
        razorTree->Branch("nJets80_JERUp", &nJets80_JERUp, "nJets80_JERUp/I");
        razorTree->Branch("mT_JERUp", &mT_JERUp, "mT_JERUp/F");
        razorTree->Branch("mTLoose_JERUp", &mTLoose_JERUp, "mTLoose_JERUp/F");
        razorTree->Branch("box_JERUp", &box_JERUp, "box_JERUp/I");
        razorTree->Branch("MR_JERDown", &MR_JERDown, "MR_JERDown/F");
        razorTree->Branch("Rsq_JERDown", &Rsq_JERDown, "Rsq_JERDown/F");
        razorTree->Branch("dPhiRazor_JERDown", &dPhiRazor_JERDown, "dPhiRazor_JERDown/F");
        razorTree->Branch("leadingJetPt_JERDown", &leadingJetPt_JERDown, "leadingJetPt_JERDown/F");
        razorTree->Branch("subleadingJetPt_JERDown", &subleadingJetPt_JERDown, "subleadingJetPt_JERDown/F");
        razorTree->Branch("nSelectedJets_JERDown", &nSelectedJets_JERDown, "nSelectedJets_JERDown/I");
        razorTree->Branch("nBTaggedJets_JERDown", &nBTaggedJets_JERDown, "nBTaggedJets_JERDown/I");
        razorTree->Branch("nJets80_JERDown", &nJets80_JERDown, "nJets80_JERDown/I");
        razorTree->Branch("mT_JERDown", &mT_JERDown, "mT_JERDown/F");
        razorTree->Branch("mTLoose_JERDown", &mTLoose_JERDown, "mTLoose_JERDown/F");
        razorTree->Branch("box_JERDown", &box_JERDown, "box_JERDown/I");
        if(isFastsimSMS){
            razorTree->Branch("mGluino", &mGluino, "mGluino/I");
            razorTree->Branch("mLSP", &mLSP, "mLSP/I");
        }
    } 
    else {
        razorTree->Branch("run", &runNum, "run/i");
        razorTree->Branch("lumi", &lumiNum, "lumi/i");
        razorTree->Branch("event", &eventNum, "event/i");
    }

    /////////////////////////////////
    //Event loop
    /////////////////////////////////

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        /////////////////////////////////
        //Begin event
        /////////////////////////////////

        //Initialize
        if (jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);

        //Reset tree variables
        nVtx = nPV;
        nSelectedJets = 0;
        nJets80 = 0;
        nBTaggedJets = 0;
        MR = -1;
        Rsq = -1;
        dPhiRazor = -9;
        mT = -1;
        mTLoose = -1;
        leadingJetPt = -1;
        subleadingJetPt = -1;
        leadingTightMuPt = -1;
        leadingTightElePt = -1;
        box = NONE;
        weight = genWeight;
        nVetoMuons = 0;
        nTightMuons = 0;
        nVetoElectrons = 0;
        nTightElectrons = 0;
        nLooseTaus = 0;
        if(!isData){
            sf_muonEffUp = 1.0;
            sf_muonEffDown = 1.0;
            sf_vetoMuonEffUp = 1.0;
            sf_vetoMuonEffDown = 1.0;
            sf_eleEffUp = 1.0;
            sf_eleEffDown = 1.0;
            sf_vetoEleEffUp = 1.0;
            sf_vetoEleEffDown = 1.0;
            sf_muonTrigUp = 1.0;
            sf_muonTrigDown = 1.0;
            sf_eleTrigUp = 1.0;
            sf_eleTrigDown = 1.0;
            sf_tauEffUp = 1.0;
            sf_tauEffDown = 1.0;
            sf_btagUp = 1.0;
            sf_btagDown = 1.0;
            sf_facScaleUp = 1.0;
            sf_facScaleDown = 1.0;
            sf_renScaleUp = 1.0;
            sf_renScaleDown = 1.0;
            sf_facRenScaleUp = 1.0;
            sf_facRenScaleDown = 1.0;
            sf_muonEffFastsimSFUp = 1.0;
            sf_muonEffFastsimSFDown = 1.0;
            sf_eleEffFastsimSFUp = 1.0;
            sf_eleEffFastsimSFDown = 1.0;
            sf_btagFastsimSFUp = 1.0;
            sf_btagFastsimSFDown = 1.0;
            MR_JESUp = -1;
            Rsq_JESUp = -1;
            dPhiRazor_JESUp = -9;
            leadingJetPt_JESUp = -1;
            subleadingJetPt_JESUp = -1;
            nSelectedJets_JESUp = 0;
            nBTaggedJets_JESUp = 0;
            nJets80_JESUp = 0;
            mT_JESUp = -1;
            mTLoose_JESUp = -1;
            box_JESUp = NONE;
            MR_JESDown = -1;
            Rsq_JESDown = -1;
            dPhiRazor_JESDown = -9;
            leadingJetPt_JESDown = -1;
            subleadingJetPt_JESDown = -1;
            nSelectedJets_JESDown = 0;
            nBTaggedJets_JESDown = 0;
            nJets80_JESDown = 0;
            mT_JESDown = -1;
            mTLoose_JESDown = -1;
            box_JESDown = NONE;
            MR_JERUp = -1;
            Rsq_JERUp = -1;
            dPhiRazor_JERUp = -9;
            leadingJetPt_JERUp = -1;
            subleadingJetPt_JERUp = -1;
            nSelectedJets_JERUp = 0;
            nBTaggedJets_JERUp = 0;
            nJets80_JERUp = 0;
            mT_JERUp = -1;
            mTLoose_JERUp = -1;
            box_JERUp = NONE;
            MR_JERDown = -1;
            Rsq_JERDown = -1;
            dPhiRazor_JERDown = -9;
            leadingJetPt_JERDown = -1;
            subleadingJetPt_JERDown = -1;
            nSelectedJets_JERDown = 0;
            nBTaggedJets_JERDown = 0;
            nJets80_JERDown = 0;
            mT_JERDown = -1;
            mTLoose_JERDown = -1;
            box_JERDown = NONE;
            MR_MESUp = -1;
            Rsq_MESUp = -1;
            dPhiRazor_MESUp = -9;
            leadingTightMuPt_MESUp = -1;
            nSelectedJets_MESUp = 0;
            nBTaggedJets_MESUp = 0;
            nJets80_MESUp = 0;
            mT_MESUp = -1;
            mTLoose_MESUp = -1;
            leadingJetPt_MESUp = -1;
            subleadingJetPt_MESUp = -1;
            box_MESUp = NONE;
            MR_MESDown = -1;
            Rsq_MESDown = -1;
            dPhiRazor_MESDown = -9;
            leadingTightMuPt_MESDown = -1;
            nSelectedJets_MESDown = 0;
            nBTaggedJets_MESDown = 0;
            nJets80_MESDown = 0;
            mT_MESDown = -1;
            mTLoose_MESDown = -1;
            leadingJetPt_MESDown = -1;
            subleadingJetPt_MESDown = -1;
            box_MESDown = NONE;
            MR_EESUp = -1;
            Rsq_EESUp = -1;
            dPhiRazor_EESUp = -9;
            leadingTightElePt_EESUp = -1;
            nSelectedJets_EESUp = 0;
            nBTaggedJets_EESUp = 0;
            nJets80_EESUp = 0;
            mT_EESUp = -1;
            mTLoose_EESUp = -1;
            leadingJetPt_EESUp = -1;
            subleadingJetPt_EESUp = -1;
            box_EESUp = NONE;
            MR_EESDown = -1;
            Rsq_EESDown = -1;
            dPhiRazor_EESDown = -9;
            leadingTightElePt_EESDown = -1;
            nSelectedJets_EESDown = 0;
            nBTaggedJets_EESDown = 0;
            nJets80_EESDown = 0;
            mT_EESDown = -1;
            mTLoose_EESDown = -1;
            leadingJetPt_EESDown = -1;
            subleadingJetPt_EESDown = -1;
            box_EESDown = NONE;
            if(isFastsimSMS){
                mGluino = -1;
                mLSP = -1;
            }
        }

        //Reset TLorentzVector collections
        vector<TLorentzVector> GoodJets; //jets used to compute hemispheres
        vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
        //For systematics
        vector<TLorentzVector> GoodJetsJESUp;
        vector<TLorentzVector> GoodJetsJESDown;
        vector<TLorentzVector> GoodJetsJERUp;
        vector<TLorentzVector> GoodJetsJERDown;
        vector<TLorentzVector> GoodJetsMESUp;
        vector<TLorentzVector> GoodJetsMESDown;
        vector<TLorentzVector> GoodJetsEESUp;
        vector<TLorentzVector> GoodJetsEESDown;
        vector<TLorentzVector> GoodLeptonsMESUp;
        vector<TLorentzVector> GoodLeptonsMESDown;
        vector<TLorentzVector> GoodLeptonsEESUp;
        vector<TLorentzVector> GoodLeptonsEESDown;

        /////////////////////////////////
        //Trigger
        /////////////////////////////////

        bool passedDileptonTrigger = false;
        bool passedSingleLeptonTrigger = false;
        bool passedHadronicTrigger= false;

        if (isData) {
            passedDileptonTrigger = bool( HLTDecision[41] || HLTDecision[43] 
                    || HLTDecision[30] || HLTDecision[31] 
                    || HLTDecision[47] || HLTDecision[48] || HLTDecision[49] || HLTDecision[50] );
            passedSingleLeptonTrigger = bool(HLTDecision[2] || HLTDecision[7] || HLTDecision[12] || HLTDecision[11] || HLTDecision[15]
                    || HLTDecision[22] || HLTDecision[23] || HLTDecision[24] || HLTDecision[25] || 
                    HLTDecision[26] || HLTDecision[27] ||
                    HLTDecision[28] || HLTDecision[29]);      
            passedHadronicTrigger = bool(HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
                    || HLTDecision[137] || HLTDecision[138] || HLTDecision[139] 
                    || HLTDecision[140] || HLTDecision[141] || HLTDecision[142] 
                    || HLTDecision[143] || HLTDecision[144]);
        } else {
            passedDileptonTrigger = bool(HLTDecision[41] || HLTDecision[43] 
                    || HLTDecision[30] || HLTDecision[31] 
                    || HLTDecision[47] || HLTDecision[48] || HLTDecision[49] || HLTDecision[50] );
            passedSingleLeptonTrigger = bool( HLTDecision[2] || HLTDecision[7] || HLTDecision[12] 
                    || HLTDecision[11] || HLTDecision[15] 
                    || HLTDecision[18] || HLTDecision[19] || HLTDecision[20] 
                    || HLTDecision[21] || HLTDecision[28] || HLTDecision[29]);
            passedHadronicTrigger = bool(HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
                    || HLTDecision[137] || HLTDecision[138] || HLTDecision[139] 
                    || HLTDecision[140] || HLTDecision[141] || HLTDecision[142] 
                    || HLTDecision[143] || HLTDecision[144]);    
        }

        //ignore trigger for Fastsim
        if(isFastsimSMS){
            passedDileptonTrigger = true;
            passedSingleLeptonTrigger = true;
            passedHadronicTrigger = true;
        }

        /////////////////////////////////
        //Pileup reweighting
        /////////////////////////////////

        //double NPU = 0;
        double pileupWeight = 1.0;
        if(!isData){
            //Get number of PU interactions
            for (int i = 0; i < nBunchXing; i++) {
                if (BunchXing[i] == 0) {
                    //NPU = nPUmean[i];
                }
            }
            //NOTE: reweight with nPV for now
            pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPV));
            //pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
        }

        /////////////////////////////////
        //Muon selection
        /////////////////////////////////

        int nVetoMuons_MESUp = 0;
        int nVetoMuons_MESDown = 0;
        int nTightMuons_MESUp = 0;
        int nTightMuons_MESDown = 0;
        TLorentzVector leadingTightMu, leadingTightMuUp, leadingTightMuDown;
        //Scale factor
        float muonEffCorrFactor = 1.0;
        float muonTrigCorrFactor = 1.0;
        float vetoMuonEffCorrFactor = 1.0;
        //MET corrections
        float MetXCorr_MESUp = 0;
        float MetYCorr_MESUp = 0;
        float MetXCorr_MESDown = 0;
        float MetYCorr_MESDown = 0;
        //Cut parameters
        const float MUON_VETO_CUT = 5;
        const float MUON_TIGHT_CUT = 20;
        //Loop muons
        for (int i = 0; i < nMuons; i++){

            //TLorentzVector for this muon
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

            //acceptance cut
            if (abs(muonEta[i]) > 2.4) continue;

            //lepton scale uncertainty
            if (!isData) {
                float muonPtUp = muonPt[i];
                float muonPtDown = muonPt[i];
                float muonEUp = muonE[i];
                float muonEDown = muonE[i];

                //get up/down muon momenta
                float sfUp = 1.0;
                float sfDown = 1.0;
                if (muonPt[i] < 100) {
                    sfUp = 1.002;
                    sfDown = 0.998;
                }
                else {
                    sfUp = 1.05;
                    sfDown = 0.95;
                }
                muonPtUp *= sfUp;
                muonPtDown *= sfDown;
                muonEUp *= sfUp;
                muonEDown *= sfDown;
                TLorentzVector thisMuonUp = makeTLorentzVector(muonPtUp, muonEta[i], muonPhi[i], muonEUp); 
                TLorentzVector thisMuonDown = makeTLorentzVector(muonPtDown, muonEta[i], muonPhi[i], muonEDown); 

                //Veto selection
                if (isVetoMuon(i)) {
                    if (muonPtUp > MUON_VETO_CUT) {
                        nVetoMuons_MESUp++;
                        GoodLeptonsMESUp.push_back(thisMuonUp); 
                        MetXCorr_MESUp += -1 * (thisMuon.Px() - thisMuonUp.Px());
                        MetYCorr_MESUp += -1 * (thisMuon.Py() - thisMuonUp.Py());
                    }
                    if (muonPtDown > MUON_VETO_CUT) {
                        nVetoMuons_MESDown++;
                        GoodLeptonsMESDown.push_back(thisMuonDown);
                        MetXCorr_MESDown += -1 * (thisMuon.Px() - thisMuonDown.Px());
                        MetYCorr_MESDown += -1 * (thisMuon.Py() - thisMuonDown.Py());
                    }
                }
                //Tight selection
                if (isTightMuon(i)) {
                    if (muonPtUp >= MUON_TIGHT_CUT) {
                        nTightMuons_MESUp++;
                    }
                    if (muonPtDown >= MUON_TIGHT_CUT) {
                        nTightMuons_MESDown++;
                    }
                    if (muonPtUp > leadingTightMuPt_MESUp){
                        leadingTightMuUp = thisMuonUp;
                        leadingTightMuPt_MESUp = muonPtUp;
                    }
                    if (muonPtDown > leadingTightMuPt_MESDown){
                        leadingTightMuDown = thisMuonDown;
                        leadingTightMuPt_MESDown = muonPtDown;
                    }
                }
            }

            //baseline pt cut
            if (muonPt[i] < MUON_VETO_CUT) continue;

            //tight lepton efficiency scale factor
            if (!isData && RazorAnalyzer::matchesGenMuon(muonEta[i], muonPhi[i]) && passedSingleLeptonTrigger && muonPt[i] > MUON_TIGHT_CUT) {
                double effTight = muTightEfficiencyHist->GetBinContent(
                        muTightEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEfficiencyHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effTight_FastsimSFUp = effTight;
                double effTight_FastsimSFDown = effTight;
                if (isFastsimSMS) { //correct efficiency for Fastsim
                    double sf = muTightEffFastsimSFHist->GetBinContent(
                        muTightEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                    double sfErr = muTightEffFastsimSFHist->GetBinError(
                        muTightEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                    effTight *= sf; 
                    effTight_FastsimSFUp *= (sf + sfErr);
                    effTight_FastsimSFDown *= (sf - sfErr);
                    //edge case: efficiency goes above 1: panic and reset
                    if (effTight_FastsimSFUp >= 1.0) {
                        effTight_FastsimSFUp = effTight;
                    }
                }
                double effTightSF = muTightEffSFHist->GetBinContent( 
                        muTightEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEffSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effTightSFErr = muTightEffSFHist->GetBinError( 
                        muTightEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEffSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effTightSFUp = effTightSF + effTightSFErr;
                double effTightSFDown = effTightSF - effTightSFErr;
                double tmpTightSF = 1.0;
                double tmpTightSFUp = 1.0;
                double tmpTightSFDown = 1.0;
                double tmpTightSF_FastsimSFUp = 1.0;
                double tmpTightSF_FastsimSFDown = 1.0;

                if (isTightMuon(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                    tmpTightSF_FastsimSFUp = effTightSF;
                    tmpTightSF_FastsimSFDown = effTightSF;
                } 
                else {
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                    tmpTightSF_FastsimSFUp = (1/effTight_FastsimSFUp - effTightSF) / (1/effTight_FastsimSFUp - 1);
                    tmpTightSF_FastsimSFDown = (1/effTight_FastsimSFDown - effTightSF) / (1/effTight_FastsimSFDown - 1);
                }
                muonEffCorrFactor *= tmpTightSF;
                sf_muonEffUp *= tmpTightSFUp/tmpTightSF;
                sf_muonEffDown *= tmpTightSFDown/tmpTightSF;
                sf_muonEffFastsimSFUp *= tmpTightSF_FastsimSFUp/tmpTightSF;
                sf_muonEffFastsimSFDown *= tmpTightSF_FastsimSFDown/tmpTightSF;
            }

            //veto lepton efficiency scale factor
            if (!isData && RazorAnalyzer::matchesGenMuon(muonEta[i], muonPhi[i]) && passedHadronicTrigger && muonPt[i] > 20) { //NOTE: do not use these scale factors below 20 GeV for now
                double effVeto = muVetoEfficiencyHist->GetBinContent(
                        muVetoEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muVetoEfficiencyHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effVetoSF = muVetoEffSFHist->GetBinContent( 
                        muVetoEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muVetoEffSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effVetoSFErr = muVetoEffSFHist->GetBinError( 
                        muVetoEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muVetoEffSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double effVetoSFUp = effVetoSF + effVetoSFErr;
                double effVetoSFDown = effVetoSF - effVetoSFErr;
                double tmpVetoSF = 1.0;
                double tmpVetoSFUp = 1.0;
                double tmpVetoSFDown = 1.0;

                if (isVetoMuon(i) && muonPt[i]) {
                    tmpVetoSF = effVetoSF;
                    tmpVetoSFUp = effVetoSFUp;
                    tmpVetoSFDown = effVetoSFDown;
                } 
                else {
                    tmpVetoSF = (1/effVeto - effVetoSF) / (1/effVeto - 1);
                    tmpVetoSFUp = (1/effVeto - effVetoSFUp) / (1/effVeto - 1);
                    tmpVetoSFDown = (1/effVeto - effVetoSFDown) / (1/effVeto - 1);
                }
                vetoMuonEffCorrFactor *= tmpVetoSF;
                sf_vetoMuonEffUp *= tmpVetoSFUp/tmpVetoSF;
                sf_vetoMuonEffDown *= tmpVetoSFDown/tmpVetoSF;
            }

            //Trigger scale factor
            if(!isData && muonPt[i] >= MUON_TIGHT_CUT){
                double trigSF = muTrigSFHist->GetBinContent( 
                        muTrigSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTrigSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double trigSFErr = muTrigSFHist->GetBinError( 
                        muTrigSFHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTrigSFHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
                double trigSFUp = trigSF + trigSFErr;
                double trigSFDown = trigSF - trigSFErr;
                if (passedSingleLeptonTrigger && isTightMuon(i)){
                    muonTrigCorrFactor *= trigSF;
                    sf_muonTrigUp *= trigSFUp/trigSF;
                    sf_muonTrigDown *= trigSFDown/trigSF;
                }
		if (isFastsimSMS) {
		  if (passedSingleLeptonTrigger && isTightMuon(i)) {
		    double singleMuonTriggerEfficiencyFromFullsim = 
		      muTrigEffFromFullsimHist->GetBinContent( muTrigEffFromFullsimHist->GetXaxis()->FindFixBin(fabs(muonEta[i])),
							       muTrigEffFromFullsimHist->GetYaxis()->FindFixBin(fmax(fmin(muonPt[i],999.9),15.01))); 
		    muonTrigCorrFactor *= singleMuonTriggerEfficiencyFromFullsim;
		  }		  
		}
            }

            //Veto selection
            if (isVetoMuon(i)){
                nVetoMuons++;
                GoodLeptons.push_back(thisMuon); 
            }
            //Tight selection
            if (isTightMuon(i) && muonPt[i] >= MUON_TIGHT_CUT){
                nTightMuons++;
                if (muonPt[i] > leadingTightMuPt){
                    leadingTightMu = thisMuon;
                    leadingTightMuPt = muonPt[i];
                }
            }
        }

        /////////////////////////////////
        //Electron selection
        /////////////////////////////////
        int nVetoElectrons_EESUp = 0;
        int nVetoElectrons_EESDown = 0;
        int nTightElectrons_EESUp = 0;
        int nTightElectrons_EESDown = 0;
        TLorentzVector leadingTightEle, leadingTightEleUp, leadingTightEleDown; //used for mT calculation
        float eleEffCorrFactor = 1.0;
        float vetoEleEffCorrFactor = 1.0;
        float eleTrigCorrFactor = 1.0;
        //MET correction
        float MetXCorr_EESUp = 0;
        float MetYCorr_EESUp = 0;
        float MetXCorr_EESDown = 0;
        float MetYCorr_EESDown = 0;
        //Cut parameters
        const float ELE_VETO_CUT = 5;
        const float ELE_TIGHT_CUT = 25;
        //Loop electrons
        for (int i = 0; i < nElectrons; i++){

            //TLorentzVector for this electron
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

            //acceptance cut
            if (fabs(eleEta[i]) > 2.5) continue;

            //lepton scale uncertainty
            if (!isData) {
                float elePtUp = elePt[i];
                float elePtDown = elePt[i];
                float eleEUp = eleE[i];
                float eleEDown = eleE[i];

                //get up/down ele momenta
                float sfUp = 1.0;
                float sfDown = 1.0;
                if (fabs(eleEta[i]) < 1.5) {
                    sfUp = 1.006;
                    sfDown = 0.994;
                }
                else {
                    sfUp = 1.015;
                    sfDown = 0.985;
                }
                elePtUp *= sfUp;
                elePtDown *= sfDown;
                eleEUp *= sfUp;
                eleEDown *= sfDown;
                TLorentzVector thisElectronUp = makeTLorentzVector(elePtUp, eleEta[i], elePhi[i], eleEUp); 
                TLorentzVector thisElectronDown = makeTLorentzVector(elePtDown, eleEta[i], elePhi[i], eleEDown); 

                //Veto selection
                if (isVetoElectron(i)) {
                    if (elePtUp > ELE_VETO_CUT) {
                        nVetoElectrons_EESUp++;
                        GoodLeptonsEESUp.push_back(thisElectronUp); 
                        MetXCorr_EESUp += -1 * (thisElectron.Px() - thisElectronUp.Px());
                        MetYCorr_EESUp += -1 * (thisElectron.Py() - thisElectronUp.Py());
                    }
                    if (elePtDown > ELE_VETO_CUT) {
                        nVetoElectrons_EESDown++;
                        GoodLeptonsEESDown.push_back(thisElectronDown);
                        MetXCorr_EESDown += -1 * (thisElectron.Px() - thisElectronDown.Px());
                        MetYCorr_EESDown += -1 * (thisElectron.Py() - thisElectronDown.Py());
                    }
                }
                //Tight selection
                if (isTightElectron(i)) {
                    if (elePtUp >= ELE_TIGHT_CUT) {
                        nTightElectrons_EESUp++;
                    }
                    if (elePtDown >= ELE_TIGHT_CUT) {
                        nTightElectrons_EESDown++;
                    }
                    if (elePtUp > leadingTightElePt_EESUp){
                        leadingTightEleUp = thisElectronUp;
                        leadingTightElePt_EESUp = elePtUp;
                    }
                    if (elePtDown > leadingTightElePt_EESDown){
                        leadingTightEleDown = thisElectronDown;
                        leadingTightElePt_EESDown = elePtDown;
                    }
                }
            }


            //baseline pt cut
            if (elePt[i] < ELE_VETO_CUT) continue;

            //Calculate MC->Data scale factors
            if (!isData && RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i]) && passedSingleLeptonTrigger && elePt[i] > ELE_TIGHT_CUT) {
                //Tight scale factor
                double effTight = eleTightEfficiencyHist->GetBinContent(
                        eleTightEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)),
                        eleTightEfficiencyHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effTight_FastsimSFUp = effTight;
                double effTight_FastsimSFDown = effTight;
                if (isFastsimSMS) { //correct efficiency for Fastsim
                    double sf = eleTightEffFastsimSFHist->GetBinContent(
                        eleTightEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)),
                        eleTightEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                    double sfErr = eleTightEffFastsimSFHist->GetBinError(
                        eleTightEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)),
                        eleTightEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                    effTight *= sf;
                    effTight_FastsimSFUp *= (sf + sfErr);
                    effTight_FastsimSFDown *= (sf - sfErr);
                    //edge case: efficiency goes above 1: panic and reset
                    if (effTight_FastsimSFUp >= 1.0) {
                        effTight_FastsimSFUp = effTight;
                    }
                }
                double effTightSF = eleTightEffSFHist->GetBinContent( 
                        eleTightEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleTightEffSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effTightSFErr = eleTightEffSFHist->GetBinError( 
                        eleTightEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleTightEffSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effTightSFUp = effTightSF + effTightSFErr;
                double effTightSFDown = effTightSF - effTightSFErr;
                double tmpTightSF = 1.0;
                double tmpTightSFUp = 1.0;
                double tmpTightSFDown = 1.0;
                double tmpTightSF_FastsimSFUp = 1.0;
                double tmpTightSF_FastsimSFDown = 1.0;

                if (isTightElectron(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                    tmpTightSF_FastsimSFUp = effTightSF;
                    tmpTightSF_FastsimSFDown = effTightSF;
                } 
                else { 
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                    tmpTightSF_FastsimSFUp = (1/effTight_FastsimSFUp - effTightSF) / (1/effTight_FastsimSFUp - 1);
                    tmpTightSF_FastsimSFDown = (1/effTight_FastsimSFDown - effTightSF) / (1/effTight_FastsimSFDown - 1);
                }
                eleEffCorrFactor *= tmpTightSF;
                sf_eleEffUp *= tmpTightSFUp/tmpTightSF;
                sf_eleEffDown *= tmpTightSFDown/tmpTightSF;
                sf_eleEffFastsimSFUp *= tmpTightSF_FastsimSFUp/tmpTightSF;
                sf_eleEffFastsimSFDown *= tmpTightSF_FastsimSFDown/tmpTightSF;
            }

            //Veto scale factor
            if (!isData && RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i]) && passedHadronicTrigger && elePt[i] > 20) { //NOTE: only use scale factors for electrons above 20 GeV for now
                double effVeto = eleVetoEfficiencyHist->GetBinContent(
                        eleVetoEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)),
                        eleVetoEfficiencyHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effVetoSF = eleVetoEffSFHist->GetBinContent( 
                        eleVetoEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleVetoEffSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effVetoSFErr = eleVetoEffSFHist->GetBinError( 
                        eleVetoEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleVetoEffSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double effVetoSFUp = effVetoSF + effVetoSFErr;
                double effVetoSFDown = effVetoSF - effVetoSFErr;

                double tmpVetoSF = 1.0;
                double tmpVetoSFUp = 1.0;
                double tmpVetoSFDown = 1.0;

                if (isVetoElectron(i)) {
                    tmpVetoSF = effVetoSF;
                    tmpVetoSFUp = effVetoSFUp;
                    tmpVetoSFDown = effVetoSFDown;
                } 
                else { 
                    tmpVetoSF = (1/effVeto - effVetoSF) / (1/effVeto - 1);
                    tmpVetoSFUp = (1/effVeto - effVetoSFUp) / (1/effVeto - 1);
                    tmpVetoSFDown = (1/effVeto - effVetoSFDown) / (1/effVeto - 1);
                }
                vetoEleEffCorrFactor *= tmpVetoSF;
                sf_vetoEleEffUp *= tmpVetoSFUp/tmpVetoSF;
                sf_vetoEleEffDown *= tmpVetoSFDown/tmpVetoSF;
            }

            //Trigger scale factor
            if(!isData && elePt[i] > ELE_TIGHT_CUT){
                double trigSF = eleTrigSFHist->GetBinContent( 
                        eleTrigSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleTrigSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double trigSFErr = eleTrigSFHist->GetBinError( 
                        eleTrigSFHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)), 
                        eleTrigSFHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
                double trigSFUp = trigSF + trigSFErr;
                double trigSFDown = trigSF - trigSFErr;
                if (passedSingleLeptonTrigger && isTightElectron(i)){
                    eleTrigCorrFactor *= trigSF;
                    sf_eleTrigUp *= trigSFUp/trigSF;
                    sf_eleTrigDown *= trigSFDown/trigSF;
                }
		
		if (isFastsimSMS) {
		  if (passedSingleLeptonTrigger && isTightElectron(i) && elePt[i] > ELE_TIGHT_CUT){
		    double singleElectronTriggerEfficiencyFromFullsim = 
		      eleTrigEffFromFullsimHist->GetBinContent( eleTrigEffFromFullsimHist->GetXaxis()->FindFixBin(fabs(eleEta[i])), 
								eleTrigEffFromFullsimHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],999.9),25.01))); 
		    eleTrigCorrFactor *= singleElectronTriggerEfficiencyFromFullsim;		    
		  }
		}
            }

            //Remove overlaps
            bool overlap = false;
            for (auto& lep : GoodLeptons){
                if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;

            //Veto selection
            if (isVetoElectron(i)){
                nVetoElectrons++;
                GoodLeptons.push_back(thisElectron);            
            }
            if (isTightElectron(i) && elePt[i] > ELE_TIGHT_CUT){
                nTightElectrons++;
                if (elePt[i] > leadingTightElePt){
                    leadingTightEle = thisElectron;
                    leadingTightElePt = elePt[i];
                }
            }
        }

        /////////////////////////////////
        //Tau selection
        /////////////////////////////////

        const float TAU_LOOSE_CUT = 20;
        float tauEffCorrFactor = 1.0;
        //Loop taus
        for (int i = 0; i < nTaus; i++){	 
            //Baseline cuts
            if (tauPt[i] < TAU_LOOSE_CUT) continue;
            if (fabs(tauEta[i]) > 2.4) continue;

            //Calculate MC->Data scale factors
            /*if (!isData && RazorAnalyzer::matchesGenTau(tauEta[i],tauPhi[i]) && passedHadronicTrigger) {
                //Loose scale factor
                double effLoose = tauLooseEfficiencyHist->GetBinContent(
                        tauLooseEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(tauPt[i],199.9),10.01)),
                        tauLooseEfficiencyHist->GetYaxis()->FindFixBin(fabs(tauEta[i]))); 
                double effLooseSF = tauLooseEffSFHist->GetBinContent( 
                        tauLooseEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(tauPt[i],199.9),10.01)), 
                        tauLooseEffSFHist->GetYaxis()->FindFixBin(fabs(tauEta[i]))); 
                double effLooseSFErr = tauLooseEffSFHist->GetBinError( 
                        tauLooseEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(tauPt[i],199.9),10.01)), 
                        tauLooseEffSFHist->GetYaxis()->FindFixBin(fabs(tauEta[i]))); 
                double effLooseSFUp = effLooseSF + effLooseSFErr;
                double effLooseSFDown = effLooseSF - effLooseSFErr;
                double tmpLooseSF = 1.0;
                double tmpLooseSFUp = 1.0;
                double tmpLooseSFDown = 1.0;

                if (isLooseTau(i)) {
                    tmpLooseSF = effLooseSF;
                    tmpLooseSFUp = effLooseSFUp;
                    tmpLooseSFDown = effLooseSFDown;
                } 
                else { 
                    tmpLooseSF = (1/effLoose - effLooseSF) / (1/effLoose - 1);
                    tmpLooseSFUp = (1/effLoose - effLooseSFUp) / (1/effLoose - 1);
                    tmpLooseSFDown = (1/effLoose - effLooseSFDown) / (1/effLoose - 1);
                }
                tauEffCorrFactor *= tmpLooseSF;
                sf_tauEffUp *= tmpLooseSFUp/tmpLooseSF;
                sf_tauEffDown *= tmpLooseSFDown/tmpLooseSF;
            }*/

            //remove overlaps
            bool overlap = false;
            for (auto& lep : GoodLeptons){
                if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;

            //Loose selection
            if (isLooseTau(i)){
                nLooseTaus++;
                TLorentzVector thisTau = makeTLorentzVectorPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
                GoodLeptons.push_back(thisTau);  
            }
        }

        /////////////////////////////////
        //Jet selection
        /////////////////////////////////

        //Type 1 MET correction 
        double MetX_Type1Corr = 0;
        double MetY_Type1Corr = 0;
        //BTag scale factor
        float btagCorrFactor = 1.0;
        //Propagate jet uncertainties to MET
        float MetXCorr_JESUp = 0;
        float MetYCorr_JESUp = 0;
        float MetXCorr_JESDown = 0;
        float MetYCorr_JESDown = 0;
        float MetXCorr_JERUp = 0;
        float MetYCorr_JERUp = 0;
        float MetXCorr_JERDown = 0;
        float MetYCorr_JERDown = 0;
        //Hadronic trigger efficiency scale factor
        float hadronicTrigCorrFactor = 1.0; //flat trigger scale factor
        if (isFastsimSMS) {
            hadronicTrigCorrFactor *= 0.975;
        }
        //Loop jets
        for (int i = 0; i < nJets; i++){

            //Apply Jet ID
            if (!jetPassIDTight[i]) continue;

            //Apply pileup jet ID 
            //UNDER CONSTRUCTION (No working point yet for Run2)
            //int level = 2; //loose jet ID
            //if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;

            //Get jet energy correction
            double tmpRho = fixedGridRhoFastjetAll;
            double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], JetCorrector);   

            //Get jet energy resolution correction, with up/down variants
            double jetEnergySmearFactor = 1.0;
            double jetEnergySmearFactorUp = 1.0;
            double jetEnergySmearFactorDown = 1.0;
            //UNDER CONSTRUCTION
            //if (!isData) {
            //    jetEnergySmearFactor = JetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random);
            //    jetEnergySmearFactorUp = UpDownJetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, jetPt[i]*JEC*jetEnergySmearFactor, "up");
            //    jetEnergySmearFactorDown = UpDownJetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, jetPt[i]*JEC*jetEnergySmearFactor, "down");
            //}

            //Get L1-only jet energy correction
            double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], 
                    JetCorrector, 0);   

            //TLorentzVector for this jet
            double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
            double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;
            TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);
            TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);

            //do up/down lepton scale uncertainties
            double deltaR = -1;
            if (!isData) {
                for (auto& lep : GoodLeptonsMESUp) {
                    double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                    if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
                }
                if (deltaR <= 0 || deltaR > 0.4) { 
                    if (jetCorrPt > 15 && jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                        MetXCorr_MESUp += -1 * (thisJet.Px() - L1CorrJet.Px());
                        MetYCorr_MESUp += -1 * (thisJet.Py() - L1CorrJet.Py());
                    }
                    if (jetCorrPt > 40 && fabs(jetEta[i]) < 3.0) {
                        GoodJetsMESUp.push_back(thisJet);
                        nSelectedJets_MESUp++;
                        if (isCSVM(i)) nBTaggedJets_MESUp++;
                        if (jetCorrPt > 80) nJets80_MESUp++;
                    }
                }

                deltaR = -1;
                for (auto& lep : GoodLeptonsMESDown) {
                    double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                    if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
                }
                if (deltaR <= 0 || deltaR > 0.4){ 
                    if (jetCorrPt > 15 && jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                        MetXCorr_MESDown += -1 * (thisJet.Px() - L1CorrJet.Px());
                        MetYCorr_MESDown += -1 * (thisJet.Py() - L1CorrJet.Py());
                    }
                    if (jetCorrPt > 40 && fabs(jetEta[i]) < 3.0) {
                        GoodJetsMESDown.push_back(thisJet);
                        nSelectedJets_MESDown++;
                        if (isCSVM(i)) nBTaggedJets_MESDown++;
                        if (jetCorrPt > 80) nJets80_MESDown++;
                    }
                }

                deltaR = -1;
                for (auto& lep : GoodLeptonsEESUp) {
                    double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                    if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
                }
                if (deltaR <= 0 || deltaR > 0.4) { 
                    if (jetCorrPt > 15 && jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                        MetXCorr_EESUp += -1 * (thisJet.Px() - L1CorrJet.Px());
                        MetYCorr_EESUp += -1 * (thisJet.Py() - L1CorrJet.Py());
                    }
                    if (jetCorrPt > 40 && fabs(jetEta[i]) < 3.0) {
                        GoodJetsEESUp.push_back(thisJet);
                        nSelectedJets_EESUp++;
                        if (isCSVM(i)) nBTaggedJets_EESUp++;
                        if (jetCorrPt > 80) nJets80_EESUp++;
                    }
                }

                deltaR = -1;
                for (auto& lep : GoodLeptonsEESDown) {
                    double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                    if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
                }
                if (deltaR <= 0 || deltaR > 0.4) { 
                    if (jetCorrPt > 15 && jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                        MetXCorr_EESDown += -1 * (thisJet.Px() - L1CorrJet.Px());
                        MetYCorr_EESDown += -1 * (thisJet.Py() - L1CorrJet.Py());
                    }
                    if (jetCorrPt > 40 && fabs(jetEta[i]) < 3.0) {
                        GoodJetsEESDown.push_back(thisJet);
                        nSelectedJets_EESDown++;
                        if (isCSVM(i)) nBTaggedJets_EESDown++;
                        if (jetCorrPt > 80) nJets80_EESDown++;
                    }
                }
            }

            //Remove overlaps
            deltaR = -1;
            for (auto& lep : GoodLeptons){
                double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if (deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

            //Propagate L1 JEC to type1 MET
            if (jetCorrPt > 15 && jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                MetX_Type1Corr += -1 * (thisJet.Px() - L1CorrJet.Px());
                MetY_Type1Corr += -1 * (thisJet.Py() - L1CorrJet.Py());
            }

            //Apply b-tagging correction factor 
            if (!isData && abs(jetPartonFlavor[i]) == 5 && abs(jetEta[i]) < 2.4 && jetCorrPt > 40) { //NOTE: b-tags only go on 40 GeV jets (this may change)
                double effMedium = btagMediumEfficiencyHist->GetBinContent(
                        btagMediumEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
                        btagMediumEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
                double effMedium_FastsimSFUp = effMedium;
                double effMedium_FastsimSFDown = effMedium;
                if (isFastsimSMS) { //correct efficiency for Fastsim
                    double sf = btagMediumEffFastsimSFHist->GetBinContent(
                        btagMediumEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.01)),
                        btagMediumEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(jetEta[i]))); 
                    double sfErr = btagMediumEffFastsimSFHist->GetBinError(
                        btagMediumEffFastsimSFHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.01)),
                        btagMediumEffFastsimSFHist->GetYaxis()->FindFixBin(fabs(jetEta[i]))); 
                    effMedium *= sf; 
                    effMedium_FastsimSFUp *= (sf + sfErr);
                    effMedium_FastsimSFDown *= (sf - sfErr);
                    //edge case: efficiency goes above 1: panic and reset
                    if (effMedium_FastsimSFUp >= 1.0) {
                        effMedium_FastsimSFUp = effMedium;
                    }
                }
                //get scale factor
                double jet_scalefactor = -1;
                double jet_scalefactorUp = -1;
                double jet_scalefactorDown = -1;
                if (jetCorrPt < 670.) { //670 is the largest pt range listed in the CSV text file
                    jet_scalefactor = btagreader.eval(BTagEntry::FLAV_B, jetEta[i], jetCorrPt); 
                    jet_scalefactorUp = btagreader_up.eval(BTagEntry::FLAV_B, jetEta[i], jetCorrPt);
                    jet_scalefactorDown = btagreader_do.eval(BTagEntry::FLAV_B, jetEta[i], jetCorrPt);
                }
                else {
                    jet_scalefactor = btagreader.eval(BTagEntry::FLAV_B, jetEta[i], 669);
                }
                //apply scale factor
                if (jet_scalefactor <= 0 || jet_scalefactorUp <= 0 || jet_scalefactorDown <= 0){
                    // cout << "Warning: b-tag scale factor is <= 0!" << endl;
                }
                else if (isCSVM(i)){
                    btagCorrFactor *= jet_scalefactor;
                    sf_btagUp *= jet_scalefactorUp/jet_scalefactor;
                    sf_btagDown *= jet_scalefactorDown/jet_scalefactor;
                }
                else {
                    double sf = (1/effMedium - jet_scalefactor) / (1/effMedium - 1);
                    btagCorrFactor *= sf;
                    sf_btagUp *= (1/effMedium - jet_scalefactorUp) / (1/effMedium - 1) / sf;
                    sf_btagDown *= (1/effMedium - jet_scalefactorDown) / (1/effMedium - 1) / sf;
                    sf_btagFastsimSFUp *= (1/effMedium_FastsimSFUp - jet_scalefactor) / (1/effMedium_FastsimSFUp - 1) / sf;
                    sf_btagFastsimSFDown *= (1/effMedium_FastsimSFDown - jet_scalefactor) / (1/effMedium_FastsimSFDown - 1) / sf;
                }
            } 

            //Cut on jet eta
            if (fabs(jetEta[i]) > 3.0) continue;

            //Get uncertainty on JEC and JER
            if(!isData){
                jecUnc->setJetEta(jetEta[i]);
                jecUnc->setJetPt(jetCorrPt);
                double unc = jecUnc->getUncertainty(true);
                double jetPtJESUp = jetCorrPt*(1+unc);
                double jetPtJESDown = jetCorrPt/(1+unc);
                double jetPtJERUp = jetPt[i]*JEC*jetEnergySmearFactorUp;
                double jetPtJERDown = jetPt[i]*JEC*jetEnergySmearFactorDown;
                double jetEJESUp = jetCorrE*(1+unc);
                double jetEJESDown = jetCorrE/(1+unc);
                double jetEJERUp = jetE[i]*JEC*jetEnergySmearFactorUp;
                double jetEJERDown = jetE[i]*JEC*jetEnergySmearFactorDown;
                TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
                TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
                TLorentzVector thisJetJERUp = makeTLorentzVector(jetPtJERUp, jetEta[i], jetPhi[i], jetEJERUp);
                TLorentzVector thisJetJERDown = makeTLorentzVector(jetPtJERDown, jetEta[i], jetPhi[i], jetEJERDown);

                //Propagate uncertainties to the MET
                if (jetPtJESUp > 20) {
                    MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
                    MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
                }
                if (jetPtJESDown > 20) {
                    MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
                    MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
                }
                if (jetPtJERUp > 20) {
                    MetXCorr_JERUp += -1 * (thisJetJERUp.Px() - thisJet.Px());
                    MetYCorr_JERUp += -1 * (thisJetJERUp.Py() - thisJet.Py());
                }
                if (jetPtJERDown > 20) {
                    MetXCorr_JERDown += -1 * (thisJetJERDown.Px() - thisJet.Px());
                    MetYCorr_JERDown += -1 * (thisJetJERDown.Py() - thisJet.Py());
                }
                //Record jets that pass the cut
                if(jetPtJESUp > 40){
                    GoodJetsJESUp.push_back(thisJetJESUp);
                    nSelectedJets_JESUp++;
                    if (isCSVM(i)) nBTaggedJets_JESUp++;
                    if (thisJetJESUp.Pt() > 80) nJets80_JESUp++;
                }
                if(jetPtJESDown > 40){
                    GoodJetsJESDown.push_back(thisJetJESDown);
                    nSelectedJets_JESDown++;
                    if (isCSVM(i)) nBTaggedJets_JESDown++;
                    if (thisJetJESDown.Pt() > 80) nJets80_JESDown++;
                }
                if(jetPtJERUp > 40){
                    GoodJetsJERUp.push_back(thisJetJERUp);
                    nSelectedJets_JERUp++;
                    if (isCSVM(i)) nBTaggedJets_JERUp++;
                    if (thisJetJERUp.Pt() > 80) nJets80_JERUp++;
                }
                if(jetPtJERDown > 40){
                    GoodJetsJERDown.push_back(thisJetJERDown);
                    nSelectedJets_JERDown++;
                    if (isCSVM(i)) nBTaggedJets_JERDown++;
                    if (thisJetJERDown.Pt() > 80) nJets80_JERDown++;
                }
            }

            //Cut on jet pt
            if (jetCorrPt < 40) continue;

            //Record this jet
            GoodJets.push_back(thisJet);
            nSelectedJets++;
            if (isCSVM(i)){ 
                nBTaggedJets++;
            }

            //Count 80 GeV jets
            if (jetCorrPt > 80) nJets80++;
        }

        //Get leading and subleading jet pt
        for (auto &jet : GoodJets){
            if (jet.Pt() > leadingJetPt){
                subleadingJetPt = leadingJetPt;
                leadingJetPt = jet.Pt();
            }
            else if (jet.Pt() > subleadingJetPt){
                subleadingJetPt = jet.Pt();
            }
        }
        //Get leading and subleading jet pt for JES/JER/MES/EES up/down
        if (!isData){
            for (auto &jet : GoodJetsJESUp){
                if (jet.Pt() > leadingJetPt_JESUp){
                    subleadingJetPt_JESUp = leadingJetPt_JESUp;
                    leadingJetPt_JESUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JESUp){
                    subleadingJetPt_JESUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJESDown){
                if (jet.Pt() > leadingJetPt_JESDown){
                    subleadingJetPt_JESDown = leadingJetPt_JESDown;
                    leadingJetPt_JESDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JESDown){
                    subleadingJetPt_JESDown = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJERUp){
                if (jet.Pt() > leadingJetPt_JERUp){
                    subleadingJetPt_JERUp = leadingJetPt_JERUp;
                    leadingJetPt_JERUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JERUp){
                    subleadingJetPt_JERUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJERDown){
                if (jet.Pt() > leadingJetPt_JERDown){
                    subleadingJetPt_JERDown = leadingJetPt_JERDown;
                    leadingJetPt_JERDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JERDown){
                    subleadingJetPt_JERDown = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsMESUp){
                if (jet.Pt() > leadingJetPt_MESUp){
                    subleadingJetPt_MESUp = leadingJetPt_MESUp;
                    leadingJetPt_MESUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_MESUp){
                    subleadingJetPt_MESUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsMESDown){
                if (jet.Pt() > leadingJetPt_MESDown){
                    subleadingJetPt_MESDown = leadingJetPt_MESDown;
                    leadingJetPt_MESDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_MESDown){
                    subleadingJetPt_MESDown = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsEESUp){
                if (jet.Pt() > leadingJetPt_EESUp){
                    subleadingJetPt_EESUp = leadingJetPt_EESUp;
                    leadingJetPt_EESUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_EESUp){
                    subleadingJetPt_EESUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsEESDown){
                if (jet.Pt() > leadingJetPt_EESDown){
                    subleadingJetPt_EESDown = leadingJetPt_EESDown;
                    leadingJetPt_EESDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_EESDown){
                    subleadingJetPt_EESDown = jet.Pt();
                }
            }
        }

        /////////////////////////////////
        //Compute razor variables and mT
        /////////////////////////////////

        //Combine jet and lepton collections
        for (auto& lep : GoodLeptons) {
            GoodJets.push_back(lep);
            if(!isData){
                GoodJetsJESUp.push_back(lep);
                GoodJetsJESDown.push_back(lep);
                GoodJetsJERUp.push_back(lep);
                GoodJetsJERDown.push_back(lep);
            }
        }
        if (!isData) {
            for (auto& lep : GoodLeptonsMESUp) {
                GoodJetsMESUp.push_back(lep);
            }
            for (auto& lep : GoodLeptonsMESDown) {
                GoodJetsMESDown.push_back(lep);
            }
            for (auto& lep : GoodLeptonsEESUp) {
                GoodJetsEESUp.push_back(lep);
            }
            for (auto& lep : GoodLeptonsEESDown) {
                GoodJetsEESDown.push_back(lep);
            }
        }

        //Get HT
        HT = 0;
        for (auto& jet : GoodJets) HT += jet.Pt();

        //Get MET
        double PFMetCustomType1CorrectedX = metPt*cos(metPhi) + MetX_Type1Corr;
        double PFMetCustomType1CorrectedY = metPt*sin(metPhi) + MetY_Type1Corr;
        TLorentzVector PFMETCustomType1Corrected; 
        PFMETCustomType1Corrected.SetPxPyPzE(PFMetCustomType1CorrectedX, PFMetCustomType1CorrectedY, 0, 
                sqrt( pow(PFMetCustomType1CorrectedX,2) + pow(PFMetCustomType1CorrectedY,2)));  
        TLorentzVector MyMET = PFMETCustomType1Corrected; //This is the MET that will be used below.
        //TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
        //TLorentzVector MyMET = PFMETType1; //This is the MET that will be used below.
        met = MyMET.Pt();

        //Compute razor variables and dPhiRazor
        vector<TLorentzVector> hemispheres = getHemispheres(GoodJets);
        MR = computeMR(hemispheres[0], hemispheres[1]); 
        Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
        dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());

        //Propagate up/down jet uncertainties to MET and recompute razor variables
        if(!isData){
            //JES up
            float PFMetXJESUp = MyMET.Px() + MetXCorr_JESUp;
            float PFMetYJESUp = MyMET.Py() + MetYCorr_JESUp;
            TLorentzVector PFMET_JESUp(PFMetXJESUp, PFMetYJESUp, 0, sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) )); 
            vector<TLorentzVector> hemispheres_JESUp = getHemispheres(GoodJetsJESUp);
            MR_JESUp = computeMR(hemispheres_JESUp[0], hemispheres_JESUp[1]); 
            Rsq_JESUp = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], PFMET_JESUp);
            dPhiRazor_JESUp = deltaPhi(hemispheres_JESUp[0].Phi(),hemispheres_JESUp[1].Phi());

            //JES down
            float PFMetXJESDown = MyMET.Px() + MetXCorr_JESDown;
            float PFMetYJESDown = MyMET.Py() + MetYCorr_JESDown;
            TLorentzVector PFMET_JESDown(PFMetXJESDown, PFMetYJESDown, 0, sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) )); 
            vector<TLorentzVector> hemispheres_JESDown = getHemispheres(GoodJetsJESDown);
            MR_JESDown = computeMR(hemispheres_JESDown[0], hemispheres_JESDown[1]); 
            Rsq_JESDown = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], PFMET_JESDown);
            dPhiRazor_JESDown = deltaPhi(hemispheres_JESDown[0].Phi(),hemispheres_JESDown[1].Phi());

            //JER up
            float PFMetXJERUp = MyMET.Px() + MetXCorr_JERUp;
            float PFMetYJERUp = MyMET.Py() + MetYCorr_JERUp;
            TLorentzVector PFMET_JERUp(PFMetXJERUp, PFMetYJERUp, 0, sqrt( pow(PFMetXJERUp,2) + pow(PFMetYJERUp,2) )); 
            vector<TLorentzVector> hemispheres_JERUp = getHemispheres(GoodJetsJERUp);
            MR_JERUp = computeMR(hemispheres_JERUp[0], hemispheres_JERUp[1]); 
            Rsq_JERUp = computeRsq(hemispheres_JERUp[0], hemispheres_JERUp[1], PFMET_JERUp);
            dPhiRazor_JERUp = deltaPhi(hemispheres_JERUp[0].Phi(),hemispheres_JERUp[1].Phi());

            //JER down
            float PFMetXJERDown = MyMET.Px() + MetXCorr_JERDown;
            float PFMetYJERDown = MyMET.Py() + MetYCorr_JERDown;
            TLorentzVector PFMET_JERDown(PFMetXJERDown, PFMetYJERDown, 0, sqrt( pow(PFMetXJERDown,2) + pow(PFMetYJERDown,2) )); 
            vector<TLorentzVector> hemispheres_JERDown = getHemispheres(GoodJetsJERDown);
            MR_JERDown = computeMR(hemispheres_JERDown[0], hemispheres_JERDown[1]); 
            Rsq_JERDown = computeRsq(hemispheres_JERDown[0], hemispheres_JERDown[1], PFMET_JERDown);
            dPhiRazor_JERDown = deltaPhi(hemispheres_JERDown[0].Phi(),hemispheres_JERDown[1].Phi());

            //MES up
            float PFMetXMESUp = MyMET.Px() + MetXCorr_MESUp;
            float PFMetYMESUp = MyMET.Py() + MetYCorr_MESUp;
            TLorentzVector PFMET_MESUp(PFMetXMESUp, PFMetYMESUp, 0, sqrt( pow(PFMetXMESUp,2) + pow(PFMetYMESUp,2) )); 
            vector<TLorentzVector> hemispheres_MESUp = getHemispheres(GoodJetsMESUp);
            MR_MESUp = computeMR(hemispheres_MESUp[0], hemispheres_MESUp[1]); 
            Rsq_MESUp = computeRsq(hemispheres_MESUp[0], hemispheres_MESUp[1], PFMET_MESUp);
            dPhiRazor_MESUp = deltaPhi(hemispheres_MESUp[0].Phi(),hemispheres_MESUp[1].Phi());

            //MES down
            float PFMetXMESDown = MyMET.Px() + MetXCorr_MESDown;
            float PFMetYMESDown = MyMET.Py() + MetYCorr_MESDown;
            TLorentzVector PFMET_MESDown(PFMetXMESDown, PFMetYMESDown, 0, sqrt( pow(PFMetXMESDown,2) + pow(PFMetYMESDown,2) )); 
            vector<TLorentzVector> hemispheres_MESDown = getHemispheres(GoodJetsMESDown);
            MR_MESDown = computeMR(hemispheres_MESDown[0], hemispheres_MESDown[1]); 
            Rsq_MESDown = computeRsq(hemispheres_MESDown[0], hemispheres_MESDown[1], PFMET_MESDown);
            dPhiRazor_MESDown = deltaPhi(hemispheres_MESDown[0].Phi(),hemispheres_MESDown[1].Phi());

            //EES up
            float PFMetXEESUp = MyMET.Px() + MetXCorr_EESUp;
            float PFMetYEESUp = MyMET.Py() + MetYCorr_EESUp;
            TLorentzVector PFMET_EESUp(PFMetXEESUp, PFMetYEESUp, 0, sqrt( pow(PFMetXEESUp,2) + pow(PFMetYEESUp,2) )); 
            vector<TLorentzVector> hemispheres_EESUp = getHemispheres(GoodJetsEESUp);
            MR_EESUp = computeMR(hemispheres_EESUp[0], hemispheres_EESUp[1]); 
            Rsq_EESUp = computeRsq(hemispheres_EESUp[0], hemispheres_EESUp[1], PFMET_EESUp);
            dPhiRazor_EESUp = deltaPhi(hemispheres_EESUp[0].Phi(),hemispheres_EESUp[1].Phi());

            //EES down
            float PFMetXEESDown = MyMET.Px() + MetXCorr_EESDown;
            float PFMetYEESDown = MyMET.Py() + MetYCorr_EESDown;
            TLorentzVector PFMET_EESDown(PFMetXEESDown, PFMetYEESDown, 0, sqrt( pow(PFMetXEESDown,2) + pow(PFMetYEESDown,2) )); 
            vector<TLorentzVector> hemispheres_EESDown = getHemispheres(GoodJetsEESDown);
            MR_EESDown = computeMR(hemispheres_EESDown[0], hemispheres_EESDown[1]); 
            Rsq_EESDown = computeRsq(hemispheres_EESDown[0], hemispheres_EESDown[1], PFMET_EESDown);
            dPhiRazor_EESDown = deltaPhi(hemispheres_EESDown[0].Phi(),hemispheres_EESDown[1].Phi());

            //compute various mT's
            if(nTightMuons + nTightElectrons > 0){
                TLorentzVector leadingLepton;
                if (leadingTightMuPt > leadingTightElePt) leadingLepton = leadingTightMu;
                else leadingLepton = leadingTightEle;

                float deltaPhiLepMet_JESUp = leadingLepton.DeltaPhi(PFMET_JESUp);
                mT_JESUp = sqrt(2*leadingLepton.Pt()*PFMET_JESUp.Pt()*(1.0 - cos(deltaPhiLepMet_JESUp))); 

                float deltaPhiLepMet_JESDown = leadingLepton.DeltaPhi(PFMET_JESDown);
                mT_JESDown = sqrt(2*leadingLepton.Pt()*PFMET_JESDown.Pt()*(1.0 - cos(deltaPhiLepMet_JESDown))); 

                float deltaPhiLepMet_JERUp = leadingLepton.DeltaPhi(PFMET_JERUp);
                mT_JERUp = sqrt(2*leadingLepton.Pt()*PFMET_JERUp.Pt()*(1.0 - cos(deltaPhiLepMet_JERUp))); 
                
                float deltaPhiLepMet_JERDown = leadingLepton.DeltaPhi(PFMET_JERDown);
                mT_JERDown = sqrt(2*leadingLepton.Pt()*PFMET_JERDown.Pt()*(1.0 - cos(deltaPhiLepMet_JERDown))); 
            }
            if (nTightMuons_MESUp + nTightElectrons > 0) {
                TLorentzVector leadingLepton_MESUp;
                if (leadingTightMuPt_MESUp > leadingTightElePt) leadingLepton_MESUp = leadingTightMuUp;
                else leadingLepton_MESUp = leadingTightEle;

                float deltaPhiLepMet_MESUp = leadingLepton_MESUp.DeltaPhi(PFMET_MESUp);
                mT_MESUp = sqrt(2*leadingLepton_MESUp.Pt()*PFMET_MESUp.Pt()*(1.0 - cos(deltaPhiLepMet_MESUp))); 
            }
            if (nTightMuons_MESDown + nTightElectrons > 0) {
                TLorentzVector leadingLepton_MESDown;
                if (leadingTightMuPt_MESDown > leadingTightElePt) leadingLepton_MESDown = leadingTightMuDown;
                else leadingLepton_MESDown = leadingTightEle;

                float deltaPhiLepMet_MESDown = leadingLepton_MESDown.DeltaPhi(PFMET_MESDown);
                mT_MESDown = sqrt(2*leadingLepton_MESDown.Pt()*PFMET_MESDown.Pt()*(1.0 - cos(deltaPhiLepMet_MESDown))); 
            }
            if (nTightMuons + nTightElectrons_EESUp > 0) {
                TLorentzVector leadingLepton_EESUp;
                if (leadingTightMuPt > leadingTightElePt_EESUp) leadingLepton_EESUp = leadingTightMu;
                else leadingLepton_EESUp = leadingTightEleUp;

                float deltaPhiLepMet_EESUp = leadingLepton_EESUp.DeltaPhi(PFMET_EESUp);
                mT_EESUp = sqrt(2*leadingLepton_EESUp.Pt()*PFMET_EESUp.Pt()*(1.0 - cos(deltaPhiLepMet_EESUp))); 
            }
            if (nTightMuons + nTightElectrons_EESDown > 0) {
                TLorentzVector leadingLepton_EESDown;
                if (leadingTightMuPt > leadingTightElePt_EESDown) leadingLepton_EESDown = leadingTightMu;
                else leadingLepton_EESDown = leadingTightEleDown;

                float deltaPhiLepMet_EESDown = leadingLepton_EESDown.DeltaPhi(PFMET_EESDown);
                mT_EESDown = sqrt(2*leadingLepton_EESDown.Pt()*PFMET_EESDown.Pt()*(1.0 - cos(deltaPhiLepMet_EESDown))); 
            }
            if (GoodLeptons.size() > 0){
                //get the highest pt lepton
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < GoodLeptons.size(); i++){
                    if (GoodLeptons[i].Pt() > maxLepPt){
                        maxLepPt = GoodLeptons[i].Pt();
                        maxLepIndex = i;
                    }
                }
                //compute MT with highest pt lepton
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet_JESUp = GoodLeptons[maxLepIndex].DeltaPhi(PFMET_JESUp);
                    mTLoose_JESUp = sqrt(2*GoodLeptons[maxLepIndex].Pt()*PFMET_JESUp.Pt()*(1.0 - cos(deltaPhiLepMet_JESUp)));

                    float deltaPhiLepMet_JESDown = GoodLeptons[maxLepIndex].DeltaPhi(PFMET_JESDown);
                    mTLoose_JESDown = sqrt(2*GoodLeptons[maxLepIndex].Pt()*PFMET_JESDown.Pt()*(1.0 - cos(deltaPhiLepMet_JESDown)));

                    float deltaPhiLepMet_JERUp = GoodLeptons[maxLepIndex].DeltaPhi(PFMET_JERUp);
                    mTLoose_JERUp = sqrt(2*GoodLeptons[maxLepIndex].Pt()*PFMET_JERUp.Pt()*(1.0 - cos(deltaPhiLepMet_JERUp)));
                    
                    float deltaPhiLepMet_JERDown = GoodLeptons[maxLepIndex].DeltaPhi(PFMET_JERDown);
                    mTLoose_JERDown = sqrt(2*GoodLeptons[maxLepIndex].Pt()*PFMET_JERDown.Pt()*(1.0 - cos(deltaPhiLepMet_JERDown)));
                }
            }
            if (GoodLeptonsMESUp.size() > 0) {
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < GoodLeptonsMESUp.size(); i++){
                    if (GoodLeptonsMESUp[i].Pt() > maxLepPt){
                        maxLepPt = GoodLeptonsMESUp[i].Pt();
                        maxLepIndex = i;
                    }
                }
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet_MESUp = GoodLeptonsMESUp[maxLepIndex].DeltaPhi(PFMET_MESUp);
                    mTLoose_MESUp = sqrt(2*GoodLeptonsMESUp[maxLepIndex].Pt()*PFMET_MESUp.Pt()*(1.0 - cos(deltaPhiLepMet_MESUp)));
                }
            }
            if (GoodLeptonsMESDown.size() > 0) {
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < GoodLeptonsMESDown.size(); i++){
                    if (GoodLeptonsMESDown[i].Pt() > maxLepPt){
                        maxLepPt = GoodLeptonsMESDown[i].Pt();
                        maxLepIndex = i;
                    }
                }
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet_MESDown = GoodLeptonsMESDown[maxLepIndex].DeltaPhi(PFMET_MESDown);
                    mTLoose_MESDown = sqrt(2*GoodLeptonsMESDown[maxLepIndex].Pt()*PFMET_MESDown.Pt()*(1.0 - cos(deltaPhiLepMet_MESDown)));
                }
            }
            if (GoodLeptonsEESUp.size() > 0) {
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < GoodLeptonsEESUp.size(); i++){
                    if (GoodLeptonsEESUp[i].Pt() > maxLepPt){
                        maxLepPt = GoodLeptonsEESUp[i].Pt();
                        maxLepIndex = i;
                    }
                }
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet_EESUp = GoodLeptonsEESUp[maxLepIndex].DeltaPhi(PFMET_EESUp);
                    mTLoose_EESUp = sqrt(2*GoodLeptonsEESUp[maxLepIndex].Pt()*PFMET_EESUp.Pt()*(1.0 - cos(deltaPhiLepMet_EESUp)));
                }
            }
            if (GoodLeptonsEESDown.size() > 0) {
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < GoodLeptonsEESDown.size(); i++){
                    if (GoodLeptonsEESDown[i].Pt() > maxLepPt){
                        maxLepPt = GoodLeptonsEESDown[i].Pt();
                        maxLepIndex = i;
                    }
                }
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet_EESDown = GoodLeptonsEESDown[maxLepIndex].DeltaPhi(PFMET_EESDown);
                    mTLoose_EESDown = sqrt(2*GoodLeptonsEESDown[maxLepIndex].Pt()*PFMET_EESDown.Pt()*(1.0 - cos(deltaPhiLepMet_EESDown)));
                }
            }
        }

        //Compute transverse mass 
        if (nTightMuons + nTightElectrons > 0){
            TLorentzVector leadingLepton;
            if (leadingTightMuPt > leadingTightElePt) leadingLepton = leadingTightMu;
            else leadingLepton = leadingTightEle;

            float deltaPhiLepMet = leadingLepton.DeltaPhi(MyMET);
            mT = sqrt(2*leadingLepton.Pt()*MyMET.Pt()*(1.0 - cos(deltaPhiLepMet))); 
        }
        //Transverse mass with leading lepton, regardless of quality
        if (GoodLeptons.size() > 0){
            //get the highest pt lepton
            float maxLepPt = -1;
            int maxLepIndex = -1;
            for (uint i = 0; i < GoodLeptons.size(); i++){
                if (GoodLeptons[i].Pt() > maxLepPt){
                    maxLepPt = GoodLeptons[i].Pt();
                    maxLepIndex = i;
                }
            }
            if (maxLepIndex >= 0){
                float deltaPhiLepMet = GoodLeptons[maxLepIndex].DeltaPhi(MyMET);
                mTLoose = sqrt(2*GoodLeptons[maxLepIndex].Pt()*MyMET.Pt()*(1.0 - cos(deltaPhiLepMet)));
            }
        }

        /////////////////////////////////
        //Categorize into boxes
        /////////////////////////////////

        //Get correct box under up/down JES/JER/MES/EES systematic
        if(!isData){
            //MESUp
            if (passedDileptonTrigger && nTightElectrons > 0 && nTightMuons_MESUp > 0){
                box_MESUp = MuEle;
            }
            else if (passedDileptonTrigger && nTightMuons_MESUp > 1){
                box_MESUp = MuMu;
            }
            else if (passedDileptonTrigger && nTightElectrons > 1){
                box_MESUp = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons_MESUp > 0){
                if (nSelectedJets_MESUp > 5) box_MESUp = MuSixJet;
                else if (nSelectedJets_MESUp > 3) box_MESUp = MuFourJet;
                else box_MESUp = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons > 0){
                if (nSelectedJets_MESUp > 5) box_MESUp = EleSixJet;
                else if (nSelectedJets_MESUp > 3) box_MESUp = EleFourJet;
                else box_MESUp = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons_MESUp > 0 && nJets80_MESUp >= 2){
                if (nSelectedJets_MESUp > 5) box_MESUp = LooseLeptonSixJet;
                else if (nSelectedJets_MESUp > 3) box_MESUp = LooseLeptonFourJet;
                else box_MESUp = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger && nJets80_MESUp >= 2){
                if (nSelectedJets_MESUp > 5) box_MESUp = SixJet;
                else if (nSelectedJets_MESUp > 3) box_MESUp = FourJet;
                else box_MESUp = DiJet;
            }
            //MESDown
            if (passedDileptonTrigger && nTightElectrons > 0 && nTightMuons_MESDown > 0){
                box_MESDown = MuEle;
            }
            else if (passedDileptonTrigger && nTightMuons_MESDown > 1){
                box_MESDown = MuMu;
            }
            else if (passedDileptonTrigger && nTightElectrons > 1){
                box_MESDown = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons_MESDown > 0){
                if (nSelectedJets_MESDown > 5) box_MESDown = MuSixJet;
                else if (nSelectedJets_MESDown > 3) box_MESDown = MuFourJet;
                else box_MESDown = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons > 0){
                if (nSelectedJets_MESDown > 5) box_MESDown = EleSixJet;
                else if (nSelectedJets_MESDown > 3) box_MESDown = EleFourJet;
                else box_MESDown = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons_MESDown > 0 && nJets80_MESDown >= 2){
                if (nSelectedJets_MESDown > 5) box_MESDown = LooseLeptonSixJet;
                else if (nSelectedJets_MESDown > 3) box_MESDown = LooseLeptonFourJet;
                else box_MESDown = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger && nJets80_MESDown >= 2){
                if (nSelectedJets_MESDown > 5) box_MESDown = SixJet;
                else if (nSelectedJets_MESDown > 3) box_MESDown = FourJet;
                else box_MESDown = DiJet;
            }
            //EESUp
            if (passedDileptonTrigger && nTightElectrons_EESUp > 0 && nTightMuons > 0){
                box_EESUp = MuEle;
            }
            else if (passedDileptonTrigger && nTightMuons > 1){
                box_EESUp = MuMu;
            }
            else if (passedDileptonTrigger && nTightElectrons_EESUp > 1){
                box_EESUp = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons > 0){
                if (nSelectedJets_EESUp > 5) box_EESUp = MuSixJet;
                else if (nSelectedJets_EESUp > 3) box_EESUp = MuFourJet;
                else box_EESUp = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons_EESUp > 0){
                if (nSelectedJets_EESUp > 5) box_EESUp = EleSixJet;
                else if (nSelectedJets_EESUp > 3) box_EESUp = EleFourJet;
                else box_EESUp = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons_EESUp + nVetoMuons > 0 && nJets80_EESUp >= 2){
                if (nSelectedJets_EESUp > 5) box_EESUp = LooseLeptonSixJet;
                else if (nSelectedJets_EESUp > 3) box_EESUp = LooseLeptonFourJet;
                else box_EESUp = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger && nJets80_EESUp >= 2){
                if (nSelectedJets_EESUp > 5) box_EESUp = SixJet;
                else if (nSelectedJets_EESUp > 3) box_EESUp = FourJet;
                else box_EESUp = DiJet;
            }
            //EESDown
            if (passedDileptonTrigger && nTightElectrons_EESDown > 0 && nTightMuons > 0){
                box_EESDown = MuEle;
            }
            else if (passedDileptonTrigger && nTightMuons > 1){
                box_EESDown = MuMu;
            }
            else if (passedDileptonTrigger && nTightElectrons_EESDown > 1){
                box_EESDown = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons > 0){
                if (nSelectedJets_EESDown > 5) box_EESDown = MuSixJet;
                else if (nSelectedJets_EESDown > 3) box_EESDown = MuFourJet;
                else box_EESDown = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons_EESDown > 0){
                if (nSelectedJets_EESDown > 5) box_EESDown = EleSixJet;
                else if (nSelectedJets_EESDown > 3) box_EESDown = EleFourJet;
                else box_EESDown = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons_EESDown + nVetoMuons > 0 && nJets80_EESDown >= 2){
                if (nSelectedJets_EESDown > 5) box_EESDown = LooseLeptonSixJet;
                else if (nSelectedJets_EESDown > 3) box_EESDown = LooseLeptonFourJet;
                else box_EESDown = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger && nJets80_EESDown >= 2){
                if (nSelectedJets_EESDown > 5) box_EESDown = SixJet;
                else if (nSelectedJets_EESDown > 3) box_EESDown = FourJet;
                else box_EESDown = DiJet;
            }

            //JES/JER
            if(passedDileptonTrigger && nTightElectrons > 0 && nTightMuons > 0){
                box_JESUp = MuEle;
                box_JESDown = MuEle;
                box_JERUp = MuEle;
                box_JERDown = MuEle;
            }
            else if(passedDileptonTrigger && nTightMuons > 1){
                box_JESUp = MuMu;
                box_JESDown = MuMu;
                box_JERUp = MuMu;
                box_JERDown = MuMu;
            }
            else if(passedDileptonTrigger && nTightElectrons>1){
                box_JESUp = EleEle;
                box_JESDown = EleEle;
                box_JERUp = EleEle;
                box_JERDown = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons > 0){
                if (nSelectedJets_JESUp > 5) box_JESUp = MuSixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = MuFourJet;
                else box_JESUp = MuJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = MuSixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = MuFourJet;
                else box_JESDown = MuJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = MuSixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = MuFourJet;
                else box_JERUp = MuJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = MuSixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = MuFourJet;
                else box_JERDown = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons > 0){
                if (nSelectedJets_JESUp > 5) box_JESUp = EleSixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = EleFourJet;
                else box_JESUp = EleJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = EleSixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = EleFourJet;
                else box_JESDown = EleJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = EleSixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = EleFourJet;
                else box_JERUp = EleJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = EleSixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = EleFourJet;
                else box_JERDown = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0){
                if (nSelectedJets_JESUp > 5 && nJets80_JESUp >= 2) box_JESUp = LooseLeptonSixJet;
                else if (nSelectedJets_JESUp > 3 && nJets80_JESUp >= 2) box_JESUp = LooseLeptonFourJet;
                else if (nJets80_JESUp >= 2) box_JESUp = LooseLeptonDiJet;

                if (nSelectedJets_JESDown > 5 && nJets80_JESDown >= 2) box_JESDown = LooseLeptonSixJet;
                else if (nSelectedJets_JESDown > 3 && nJets80_JESDown >= 2) box_JESDown = LooseLeptonFourJet;
                else if (nJets80_JESDown >= 2) box_JESDown = LooseLeptonDiJet;

                if (nSelectedJets_JERUp > 5 && nJets80_JERUp >= 2) box_JERUp = LooseLeptonSixJet;
                else if (nSelectedJets_JERUp > 3 && nJets80_JERUp >= 2) box_JERUp = LooseLeptonFourJet;
                else if (nJets80_JERUp >= 2) box_JERUp = LooseLeptonDiJet;

                if (nSelectedJets_JERDown > 5 && nJets80_JERDown >= 2) box_JERDown = LooseLeptonSixJet;
                else if (nSelectedJets_JERDown > 3 && nJets80_JERDown >= 2) box_JERDown = LooseLeptonFourJet;
                else if (nJets80_JERDown >= 2) box_JERDown = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger){
                if (nSelectedJets_JESUp > 5 && nJets80_JESUp >= 2) box_JESUp = SixJet;
                else if (nSelectedJets_JESUp > 3 && nJets80_JESUp >= 2) box_JESUp = FourJet;
                else if (nJets80_JESUp >= 2) box_JESUp = DiJet;

                if (nSelectedJets_JESDown > 5 && nJets80_JESDown >= 2) box_JESDown = SixJet;
                else if (nSelectedJets_JESDown > 3 && nJets80_JESDown >= 2) box_JESDown = FourJet;
                else if (nJets80_JESDown >= 2) box_JESDown = DiJet;

                if (nSelectedJets_JERUp > 5 && nJets80_JERUp >= 2) box_JERUp = SixJet;
                else if (nSelectedJets_JERUp > 3 && nJets80_JERUp >= 2) box_JERUp = FourJet;
                else if (nJets80_JERUp >= 2) box_JERUp = DiJet;

                if (nSelectedJets_JERDown > 5 && nJets80_JERDown >= 2) box_JERDown = SixJet;
                else if (nSelectedJets_JERDown > 3 && nJets80_JERDown >= 2) box_JERDown = FourJet;
                else if (nJets80_JERDown >= 2) box_JERDown = DiJet;
            }
        }

        //Nominal box
        if (passedDileptonTrigger && nTightElectrons > 0 && nTightMuons > 0){
            box = MuEle;
        }
        else if (passedDileptonTrigger && nTightMuons > 1){
            box = MuMu;
        }
        else if (passedDileptonTrigger && nTightElectrons > 1){
            box = EleEle;
        }
        else if (passedSingleLeptonTrigger && nTightMuons > 0){
            if (nSelectedJets > 5) box = MuSixJet;
            else if (nSelectedJets > 3) box = MuFourJet;
            else box = MuJet;
        }
        else if (passedSingleLeptonTrigger && nTightElectrons > 0){
            if (nSelectedJets > 5) box = EleSixJet;
            else if (nSelectedJets > 3) box = EleFourJet;
            else box = EleJet;
        }
        else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0 && nJets80 >= 2){
            if (nSelectedJets > 5) box = LooseLeptonSixJet;
            else if (nSelectedJets > 3) box = LooseLeptonFourJet;
            else box = LooseLeptonDiJet;
        }
        else if (passedHadronicTrigger && nJets80 >= 2){
            if (nSelectedJets > 5) box = SixJet;
            else if (nSelectedJets > 3) box = FourJet;
            else box = DiJet;
        }

        /////////////////////////////////
        //Scale and PDF variations
        /////////////////////////////////

        sf_facScaleUp = (*scaleWeights)[1]/genWeight;
        sf_facScaleDown = (*scaleWeights)[2]/genWeight;
        sf_renScaleUp = (*scaleWeights)[3]/genWeight;
        sf_renScaleDown = (*scaleWeights)[6]/genWeight;
        sf_facRenScaleUp = (*scaleWeights)[4]/genWeight;
        sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;

        SumScaleWeights->Fill(0.0, sf_facScaleUp);
        SumScaleWeights->Fill(1.0, sf_facScaleDown);
        SumScaleWeights->Fill(2.0, sf_renScaleUp);
        SumScaleWeights->Fill(3.0, sf_renScaleDown);
        SumScaleWeights->Fill(4.0, sf_facRenScaleUp);
        SumScaleWeights->Fill(5.0, sf_facRenScaleDown);

        //for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) 
        //{
            //SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
        //}

        /////////////////////////////////
        //Apply scale factors
        /////////////////////////////////

        //Nominal event weight
        if(!isData){
            weight *= pileupWeight; 
            if (passedSingleLeptonTrigger) {
                weight *= muonEffCorrFactor;
                weight *= muonTrigCorrFactor;
                weight *= eleEffCorrFactor;
                weight *= eleTrigCorrFactor;
            }
            else if(passedHadronicTrigger) {
                weight *= vetoMuonEffCorrFactor;
                weight *= vetoEleEffCorrFactor;
                weight *= tauEffCorrFactor;
                weight *= hadronicTrigCorrFactor;
            }
            weight *= btagCorrFactor;    
        }

        //Fill normalization histogram
        NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
        SumWeights->Fill(1.0, weight);

        /////////////////////////////////
        //SMS information
        /////////////////////////////////

        bool parsedLHE = false;
        if(isFastsimSMS && lheComments->size() > 0){
            //parse lhe comment string to get gluino and LSP masses
            stringstream parser((*lheComments)[lheComments->size()-1]);
            string item;
            getline(parser, item, '_'); //prefix
            if(getline(parser, item, '_')){ //gluino mass 
                mGluino = atoi(item.c_str());
                if(getline(parser, item, '_')){ //LSP mass 
                    mLSP = atoi(item.c_str());
                    pair<int,int> smsPair = make_pair(mGluino, mLSP);
                    parsedLHE = true;
                    if (smsFiles.count(smsPair) == 0){ //create file and tree
                        //format file name
                        string thisFileName = outFileName;
                        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
                        thisFileName += "_" + to_string(mGluino) + "_" + to_string(mLSP) + ".root";
                        
                        smsFiles[smsPair] = new TFile(thisFileName.c_str(), "recreate");
                        smsTrees[smsPair] = razorTree->CloneTree(0);
                        smsNEvents[smsPair] = new TH1F(Form("NEvents%d%d", mGluino, mLSP), "NEvents", 1,0.5,1.5);
                        smsSumWeights[smsPair] = new TH1F(Form("SumWeights%d%d", mGluino, mLSP), "SumWeights", 1,0.5,1.5);
                        smsSumScaleWeights[smsPair] = new TH1F(Form("SumScaleWeights%d%d", mGluino, mLSP), "SumScaleWeights", 6,-0.5,5.5);
                        //smsSumPdfWeights[smsPair] = new TH1F(Form("SumPdfWeights%d%d", mGluino, mLSP), "SumPdfWeights", NUM_PDF_WEIGHTS,-0.5,NUM_PDF_WEIGHTS-0.5);
                        cout << "Created new output file " << thisFileName << endl;
                    }
                    //Fill NEvents hist 
                    smsNEvents[smsPair]->Fill(1.0, genWeight);
                    smsSumWeights[smsPair]->Fill(1.0, weight);

                    smsSumScaleWeights[smsPair]->Fill(0.0, sf_facScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(1.0, sf_facScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(2.0, sf_renScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(3.0, sf_renScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(4.0, sf_facRenScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(5.0, sf_facRenScaleDown);

                    //for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) 
                    //{
                        //smsSumPdfWeights[smsPair]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
                    //}
                }
            }
        }

        /////////////////////////////////
        //Baseline cuts
        /////////////////////////////////

        //Razor
        if (MR < 300 && MR_JESUp < 300 && MR_JESDown < 300 && MR_JERUp < 300 && MR_JERDown < 300 && MR_MESUp < 300 && MR_MESDown < 300 && MR_EESUp < 300 && MR_EESDown < 300) {
            continue;
        }
        if (Rsq < 0.15 && Rsq_JESUp < 0.15 && Rsq_JESDown < 0.15 && Rsq_JERUp < 0.15 && Rsq_JERDown < 0.15 && Rsq_MESUp < 0.15 && Rsq_MESDown < 0.15 && Rsq_EESUp < 0.15 && Rsq_EESDown < 0.15) {
            continue;
        }

        //Continue if this event is not in any box
        if(box == NONE && box_JESUp == NONE && box_JESDown == NONE && box_JERUp == NONE && box_JERDown == NONE && box_MESUp == NONE && box_MESDown == NONE && box_EESUp == NONE && box_EESDown == NONE) {
            continue; 
        }

        //Trigger
        if(!passedDileptonTrigger && !passedSingleLeptonTrigger && !passedHadronicTrigger) {
            continue;
        }

        /////////////////////////////////
        //Noise filters
        /////////////////////////////////

        if(isData){
            if(!Flag_HBHENoiseFilter || !Flag_HBHEIsoNoiseFilter || !Flag_goodVertices || !Flag_eeBadScFilter) {
                continue;
            }
        }

        //Fill tree
        if(!isFastsimSMS){
            razorTree->Fill();
        }
        else if(parsedLHE){
            pair<int,int> smsPair = make_pair(mGluino, mLSP);
            smsTrees[smsPair]->Fill();
        }

    }//end of event loop

    if(!isFastsimSMS){
        cout << "Writing output tree..." << endl;
        outFile->cd();
        razorTree->Write();
        NEvents->Write();
        SumWeights->Write();
        SumScaleWeights->Write();
        //SumPdfWeights->Write();
    }
    else{
        for(auto &filePtr : smsFiles){
            cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
            filePtr.second->cd();
            smsTrees[filePtr.first]->Write();
            smsNEvents[filePtr.first]->Write("NEvents");
            smsSumWeights[filePtr.first]->Write("SumWeights");
            smsSumScaleWeights[filePtr.first]->Write("SumScaleWeights");
            //smsSumPdfWeights[filePtr.first]->Write("SumPdfWeights");
        }
    }

    outFile->Close();
}
