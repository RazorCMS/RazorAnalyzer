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
    if(!isData){
        TFile *eleEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.root");
        eleTightEfficiencyHist = (TH2D*)eleEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
        eleVetoEfficiencyHist = (TH2D*)eleEfficiencyFile->Get("ElectronEff_Veto_Fullsim");
        assert(eleTightEfficiencyHist);
        TFile *muEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.root");
        muTightEfficiencyHist = (TH2D*)muEfficiencyFile->Get("MuonEff_Tight_Fullsim");
        muVetoEfficiencyHist = (TH2D*)muEfficiencyFile->Get("MuonEff_Veto_Fullsim");
        assert(muTightEfficiencyHist);
        TFile *tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
        tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");
        assert(tauLooseEfficiencyHist);
        TFile *btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/BTagEffFastsimToFullsimCorrectionFactors.root");
        btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("BTagEff_Medium_Fullsim");
        assert(btagMediumEfficiencyHist);
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
    //For pdf uncertainties
    //For jet uncertainties
    float MR_JESUp, Rsq_JESUp, dPhiRazor_JESUp, leadingJetPt_JESUp, subleadingJetPt_JESUp; 
    float MR_JESDown, Rsq_JESDown, dPhiRazor_JESDown, leadingJetPt_JESDown, subleadingJetPt_JESDown;
    float MR_JERUp, Rsq_JERUp, dPhiRazor_JERUp, leadingJetPt_JERUp, subleadingJetPt_JERUp;
    float MR_JERDown, Rsq_JERDown, dPhiRazor_JERDown, leadingJetPt_JERDown, subleadingJetPt_JERDown;
    int nSelectedJets_JESUp, nSelectedJets_JESDown, nSelectedJets_JERUp, nSelectedJets_JERDown;
    int nBTaggedJets_JESUp, nBTaggedJets_JESDown, nBTaggedJets_JERUp, nBTaggedJets_JERDown;
    int nJets80_JESUp, nJets80_JESDown, nJets80_JERUp, nJets80_JERDown;
    RazorBox box_JESUp, box_JESDown, box_JERUp, box_JERDown;
    float mT_JESUp, mT_JESDown, mT_JERUp, mT_JERDown;
    float mTLoose_JESUp, mTLoose_JESDown, mTLoose_JERUp, mTLoose_JERDown;
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

    if (!isData) {    
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
        razorTree->Branch("sf_facScaleUp", &sf_facScaleUp, "sf_facScaleUp/F");
        razorTree->Branch("sf_facScaleDown", &sf_facScaleDown, "sf_facScaleDown/F");
        razorTree->Branch("sf_renScaleUp", &sf_renScaleUp, "sf_renScaleUp/F");
        razorTree->Branch("sf_renScaleDown", &sf_renScaleDown, "sf_renScaleDown/F");
        razorTree->Branch("sf_facRenScaleUp", &sf_facRenScaleUp, "sf_facRenScaleUp/F");
        razorTree->Branch("sf_facRenScaleDown", &sf_facRenScaleDown, "sf_facRenScaleDown/F");
        razorTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights); //get PDF weights directly from RazorEvents
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

        int nVetoMuons = 0;
        int nLooseMuons = 0;
        int nTightMuons = 0;
        TLorentzVector leadingTightMu;
        //Scale factor
        float muonEffCorrFactor = 1.0;
        float muonTrigCorrFactor = 1.0;
        float vetoMuonEffCorrFactor = 1.0;
        //Cut parameters
        const float MUON_VETO_CUT = 5;
        const float MUON_LOOSE_CUT = 20;
        //Loop muons
        for (int i = 0; i < nMuons; i++){

            //Baseline cuts
            if (muonPt[i] < MUON_VETO_CUT) continue;
            if (abs(muonEta[i]) > 2.4) continue;

            //tight lepton efficiency scale factor
            if (!isData && RazorAnalyzer::matchesGenMuon(muonEta[i], muonPhi[i]) && passedSingleLeptonTrigger && muonPt[i] > MUON_LOOSE_CUT) {
                double effTight = muTightEfficiencyHist->GetBinContent(
                        muTightEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)),
                        muTightEfficiencyHist->GetYaxis()->FindFixBin(fabs(muonEta[i]))); 
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

                if (isTightMuon(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                } 
                else {
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                }
                muonEffCorrFactor *= tmpTightSF;
                sf_muonEffUp *= tmpTightSFUp/tmpTightSF;
                sf_muonEffDown *= tmpTightSFDown/tmpTightSF;
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
            //NOTE: implemented for single lepton trigger only!  
            if(!isData && muonPt[i] >= MUON_LOOSE_CUT){
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
		  if (passedSingleLeptonTrigger && isTightMuon(i) && muonPt[i] >= MUON_LOOSE_CUT) {
		    double singleMuonTriggerEfficiencyFromFullsim = 
		      muTrigEffFromFullsimHist->GetBinContent( muTrigEffFromFullsimHist->GetXaxis()->FindFixBin(fabs(muonEta[i])),
							       muTrigEffFromFullsimHist->GetYaxis()->FindFixBin(fmax(fmin(muonPt[i],999.9),15.01))); 
		    muonTrigCorrFactor *= singleMuonTriggerEfficiencyFromFullsim;
		  }		  
		}
            }

            //TLorentzVector for this muon
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

            //Veto selection
            if (isVetoMuon(i)){
                nVetoMuons++;
                GoodLeptons.push_back(thisMuon); 
            }
            //Loose selection
            if (isLooseMuon(i) && muonPt[i] >= MUON_LOOSE_CUT) nLooseMuons++;
            //Tight selection
            if (isTightMuon(i) && muonPt[i] >= MUON_LOOSE_CUT){
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
        int nVetoElectrons = 0;
        int nLooseElectrons = 0;
        int nTightElectrons = 0;
        TLorentzVector leadingTightEle; //used for mT calculation
        float eleEffCorrFactor = 1.0;
        float vetoEleEffCorrFactor = 1.0;
        float eleTrigCorrFactor = 1.0;
        //Cut parameters
        const float ELE_VETO_CUT = 5;
        const float ELE_LOOSE_CUT = 25;
        //Loop electrons
        for (int i = 0; i < nElectrons; i++){

            //Baseline cuts
            if (elePt[i] < ELE_VETO_CUT) continue;
            if (fabs(eleEta[i]) > 2.5) continue;

            //Calculate MC->Data scale factors
            if (!isData && RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i]) && passedSingleLeptonTrigger && elePt[i] > ELE_LOOSE_CUT) {
                //Tight scale factor
                double effTight = eleTightEfficiencyHist->GetBinContent(
                        eleTightEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)),
                        eleTightEfficiencyHist->GetYaxis()->FindFixBin(fabs(eleEta[i]))); 
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

                if (isTightElectron(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                } 
                else { 
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                }
                eleEffCorrFactor *= tmpTightSF;
                sf_eleEffUp *= tmpTightSFUp/tmpTightSF;
                sf_eleEffDown *= tmpTightSFDown/tmpTightSF;
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
            //NOTE: implemented for single lepton trigger only!
            if(!isData && elePt[i] > ELE_LOOSE_CUT){
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
		  if (passedSingleLeptonTrigger && isTightElectron(i) && elePt[i] > ELE_LOOSE_CUT){
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

            //TLorentzVector for this electron
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

            //Veto selection
            if (isVetoElectron(i)){
                nVetoElectrons++;
                GoodLeptons.push_back(thisElectron);            
            }
            //Loose selection
            if (isLooseElectron(i) && elePt[i] > ELE_LOOSE_CUT) nLooseElectrons++;
            if (isTightElectron(i) && elePt[i] > ELE_LOOSE_CUT){
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

        int nLooseTaus = 0;
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
        float hadronicTrigCorrFactor = 0.972; //flat trigger scale factor
        //Loop jets
        for (int i = 0; i < nJets; i++){

            //Remove overlaps
            double deltaR = -1;
            for (auto& lep : GoodLeptons){
                double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if (deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

            //Apply Jet ID
            if (!jetPassIDTight[i]) continue;

            //Get jet energy correction
            double tmpRho = fixedGridRhoFastjetAll;
            double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], JetCorrector);   

            //Get L1-only jet energy correction
            double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], 
                    JetCorrector, 0);   

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

            //TLorentzVector for this jet
            double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
            double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;
            TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);
            TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);

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
                    sf_btagUp *= jet_scalefactorUp;
                    sf_btagDown *= jet_scalefactorDown;
                    // cout << "b-tag scale factor: " << jet_scalefactor << " (" << jetCorrPt << ", " << jetEta[i] << ")" << endl;
                    // cout << "up: " << jet_scalefactorUp << endl;
                    // cout << "down: " << jet_scalefactorDown << endl;
                }
                else {
                    btagCorrFactor *= (1/effMedium - jet_scalefactor) / (1/effMedium - 1);
                    sf_btagUp *= (1/effMedium - jet_scalefactorUp) / (1/effMedium - 1);
                    sf_btagDown *= (1/effMedium - jet_scalefactorDown) / (1/effMedium - 1);
                    // cout << "b-tag scale factor: " << (1/effMedium - jet_scalefactor)/(1/effMedium-1) << " (" << jetCorrPt << ", " << jetEta[i] << ")" << endl;
                    // cout << "up: " << (1/effMedium - jet_scalefactorUp)/(1/effMedium-1) << endl;
                    // cout << "down: " << (1/effMedium - jet_scalefactorDown)/(1/effMedium-1) << endl;
                }
            } 

            //Apply pileup jet ID 
            //UNDER CONSTRUCTION (No working point yet for Run2)
            //int level = 2; //loose jet ID
            //if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;

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
        //Get leading and subleading jet pt for JES/JER up/down
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

        //Get MET
        double PFMetCustomType1CorrectedX = metPt*cos(metPhi) + MetX_Type1Corr;
        double PFMetCustomType1CorrectedY = metPt*sin(metPhi) + MetY_Type1Corr;
        TLorentzVector PFMETCustomType1Corrected; 
        PFMETCustomType1Corrected.SetPxPyPzE(PFMetCustomType1CorrectedX, PFMetCustomType1CorrectedY, 0, 
                sqrt( pow(PFMetCustomType1CorrectedX,2) + pow(PFMetCustomType1CorrectedY,2)));  
        TLorentzVector MyMET = PFMETCustomType1Corrected; //This is the MET that will be used below.
        //TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
        //TLorentzVector MyMET = PFMETType1; //This is the MET that will be used below.

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

        //Get correct box under up/down JES/JER systematic
        if(!isData){
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

        sf_facScaleUp = (*scaleWeights)[1];
        sf_facScaleDown = (*scaleWeights)[2];
        sf_renScaleUp = (*scaleWeights)[3];
        sf_renScaleDown = (*scaleWeights)[6];
        sf_facRenScaleUp = (*scaleWeights)[4];
        sf_facRenScaleDown = (*scaleWeights)[8];

        SumScaleWeights->Fill(0.0, sf_facScaleUp);
        SumScaleWeights->Fill(1.0, sf_facScaleDown);
        SumScaleWeights->Fill(2.0, sf_renScaleUp);
        SumScaleWeights->Fill(3.0, sf_renScaleDown);
        SumScaleWeights->Fill(4.0, sf_facRenScaleUp);
        SumScaleWeights->Fill(5.0, sf_facRenScaleDown);

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
            //weight *= btagCorrFactor;    
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
                        cout << "Created new output file " << thisFileName << endl;
                    }
                    //Fill NEvents hist 
                    smsNEvents[smsPair]->Fill(1.0);
                    smsSumWeights[smsPair]->Fill(1.0, weight);

                    smsSumScaleWeights[smsPair]->Fill(0.0, sf_facScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(1.0, sf_facScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(2.0, sf_renScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(3.0, sf_renScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(4.0, sf_facRenScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(5.0, sf_facRenScaleDown);
                }
            }
        }

        /////////////////////////////////
        //Baseline cuts
        /////////////////////////////////

        //Razor
        if (MR < 300 && MR_JESUp < 300 && MR_JESDown < 300 && MR_JERUp < 300 && MR_JERDown < 300) continue;
        if (Rsq < 0.15 && Rsq_JESUp < 0.15 && Rsq_JESDown < 0.15 && Rsq_JERUp < 0.15 && Rsq_JERDown < 0.15) continue;

        //Continue if this event is not in any box
        if(box == NONE && box_JESUp == NONE && box_JESDown == NONE && box_JERUp == NONE && box_JERDown == NONE) continue; 

        //Trigger
        if(!passedDileptonTrigger && !passedSingleLeptonTrigger && !passedHadronicTrigger) continue;

        /////////////////////////////////
        //Noise filters
        /////////////////////////////////

        if(!isFastsimSMS){
	  if(!Flag_HBHENoiseFilter) continue;
	  if(!Flag_HBHETightNoiseFilter) continue;
	  if(!Flag_CSCTightHaloFilter) continue;
	  if(!Flag_goodVertices) continue;
	  if(!Flag_eeBadScFilter) continue;
	  if(!Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
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
    }
    else{
        for(auto &filePtr : smsFiles){
            cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
            filePtr.second->cd();
            smsTrees[filePtr.first]->Write();
            smsNEvents[filePtr.first]->Write("NEvents");
            smsSumWeights[filePtr.first]->Write("SumWeights");
            smsSumScaleWeights[filePtr.first]->Write("SumScaleWeights");
        }
    }

    outFile->Close();
}
