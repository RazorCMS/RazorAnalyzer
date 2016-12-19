#include "RazorHelper.h"

// Constructor
RazorHelper::RazorHelper(std::string tag_, bool isData_, bool isFastsim_): 
        tag(tag_), isData(isData_), isFastsim(isFastsim_) {
    std::cout << "RazorHelper initializing with tag " << tag << std::endl;

    eleVetoEffSFMinPt = -1;
    muVetoEffSFMinPt = -1;

    // check that CMSSW is set up
    loadCMSSWPath();
    if (cmsswPath == "") {
        loadTag_Null();
        return;
    }

    // tag for 2015 data
    else if (tag == "Razor2015") {
        loadTag_Razor2015();
    }

    // tag for 2015 76X ReReco data
    else if (tag == "Razor2015_76X") {
        loadTag_Razor2015_76X();
    }

    // tag for 2016 80X PromptReco data
    else if (tag == "Razor2016_MoriondRereco") {
        loadTag_Razor2016_MoriondRereco();
    }

    // tag for 2016G 80X data
    else if (tag == "Razor2016G_80X") {
        loadTag_Razor2016G_80X();
    }

    // tag for 2016G unblinded 80X data
    else if (tag == "Razor2016G_SUSYUnblind_80X") {
      eleVetoEffSFMinPt = 15.01;
      muVetoEffSFMinPt = 15.01;
      loadTag_Razor2016G_SUSYUnblind_80X();
    }

    // tag for 2016 ICHEP 80X data
    else if (tag == "Razor2016_ICHEP_80X") {
        loadTag_Razor2016_ICHEP_80X();
    }

    // tag not found
    else {
        std::cout << "Error in RazorHelper::RazorHelper : specified tag " << tag << " is not supported!" << std::endl;
        loadTag_Null();
        return;
    }

}

// Destructor -- close all of the open TFile objects
RazorHelper::~RazorHelper() {
    if (pileupWeightFile) {
        pileupWeightFile->Close();
        delete pileupWeightFile;
    }
    if (eleTightEfficiencyFile) {
        eleTightEfficiencyFile->Close();
        delete eleTightEfficiencyFile;
    }
    if (eleVetoEfficiencyFile) {
        eleVetoEfficiencyFile->Close();
        delete eleVetoEfficiencyFile;
    }
    if (eleEffSFFile) {
        eleEffSFFile->Close();
        delete eleEffSFFile;
    }
    if (vetoEleEffSFFile) {
        vetoEleEffSFFile->Close();
        delete vetoEleEffSFFile;
    }
    if (muTightEfficiencyFile) {
        muTightEfficiencyFile->Close();
        delete muTightEfficiencyFile;
    }
    if (muVetoEfficiencyFile) {
        muVetoEfficiencyFile->Close();
        delete muVetoEfficiencyFile;
    }
    if (muEffSFFile) {
        muEffSFFile->Close();
        delete muEffSFFile;
    }
    if (vetoMuEffSFFile) {
        vetoMuEffSFFile->Close();
        delete vetoMuEffSFFile;
    }
    if (tauEfficiencyFile) {
        tauEfficiencyFile->Close();
        delete tauEfficiencyFile;
    }
    if (btagEfficiencyFile) {
        btagEfficiencyFile->Close();
        delete btagEfficiencyFile;
    }
    if (btagCharmEfficiencyFile) {
        btagCharmEfficiencyFile->Close();
        delete btagCharmEfficiencyFile;
    }
    if (btagLightJetsEfficiencyFile) {
        btagLightJetsEfficiencyFile->Close();
        delete btagLightJetsEfficiencyFile;
    }
    if (eleTrigSFFile) {
        eleTrigSFFile->Close();
        delete eleTrigSFFile;
    }
    if (muTrigSFFile) {
        muTrigSFFile->Close();
        delete muTrigSFFile;
    }
    if (eleTrigEffFile) {
        eleTrigEffFile->Close();
        delete eleTrigEffFile;
    }
    if (muTrigEffFile) {
        muTrigEffFile->Close();
        delete muTrigEffFile;
    }
    for (uint i=0;i < JetCorrector.size() ; i++) {
      if (JetCorrector[i]) delete JetCorrector[i];
      if (JetResolutionParameters[i]) delete JetResolutionParameters[i];
      if (JetResolutionCalculator[i]) delete JetResolutionCalculator[i];
    }  
    if (btagcalib) delete btagcalib;
    if (btagreader) delete btagreader;
    if (btagreader_up) delete btagreader_up;
    if (btagreader_do) delete btagreader_do;
    if (btagreaderMistag) delete btagreaderMistag;
    if (btagreaderMistag_up) delete btagreaderMistag_up;
    if (btagreaderMistag_do) delete btagreaderMistag_do;
    if (btagreaderfastsim) delete btagreaderfastsim;
    if (btagreaderfastsim_up) delete btagreaderfastsim_up;
    if (btagreaderfastsim_do) delete btagreaderfastsim_do;
}

// Retrieves CMSSW_BASE and stores in variable cmsswPath
void RazorHelper::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in RazorHelper::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}

void RazorHelper::loadTag_Null() {
    std::cout << "Warning: initializing all RazorHelper files and histograms to 0" << std::endl;

    // pileup weights
    pileupWeightFile = 0;
    pileupWeightHist = 0;
    pileupWeightSysUpHist = 0;
    pileupWeightSysDownHist = 0;

    // electron efficiencies and scale factors
    eleTightEfficiencyFile = 0;
    eleVetoEfficiencyFile = 0;
    eleEffSFFile = 0;
    vetoEleEffSFFile = 0;
    eleTightEfficiencyHist = 0;
    eleVetoEfficiencyHist = 0;
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = 0;
    eleVetoEffSFHist = 0;

    // muon efficiencies and scale factors
    muTightEfficiencyFile = 0;
    muVetoEfficiencyFile = 0;
    muEffSFFile = 0;
    vetoMuEffSFFile = 0;
    muTightEfficiencyHist = 0;
    muVetoEfficiencyHist = 0;
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = 0;
    muVetoEffSFHist = 0;

    // tau efficiencies and scale factors
    tauEfficiencyFile = 0;
    tauLooseEfficiencyHist = 0;

    // b-tag efficiencies and scale factors
    btagEfficiencyFile = 0;
    btagCharmEfficiencyFile = 0;
    btagLightJetsEfficiencyFile = 0;
    btagMediumEfficiencyHist = 0;
    btagMediumCharmEfficiencyHist = 0;
    btagMediumLightJetsEfficiencyHist = 0;

    // single lepton trigger scale factors
    eleTrigSFFile = 0;
    eleTrigSFHist = 0;
    muTrigSFFile = 0; 
    muTrigSFHist = 0;

    // single lepton trigger efficiencies
    eleTrigEffFile = 0;
    eleTrigEffHist = 0;
    muTrigEffFile = 0; 
    muTrigEffHist = 0;

}

////////////////////////////////////////////////
//  2015 
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2015() {
    loadPileup_Razor2015();
    loadLepton_Razor2015();
    loadJECs_Razor2015();
    loadBTag_Razor2015();
    loadTrigger_Razor2015();
}

void RazorHelper::loadPileup_Razor2015() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open(
            Form("%s/src/RazorAnalyzer/data/PileupReweight_Spring15MCTo2015Data.root", cmsswPath.c_str()));
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2015(){

    // electron efficiencies and scale factors
    std::cout << "RazorHelper: loading electron efficiency histograms" << std::endl;
    // For 2015 the lepton efficiency histograms for tight and veto IDs in fullsim are both
    // stored in the fastsim-to-fullsim scale factor files. 
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.root");
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.root");
    eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightElectronSelectionEffDenominatorReco_2015Final_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_VetoElectronSelectionEffDenominatorReco_2015Final_Golden.root");
    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronEff_Veto_Fullsim");
    eleTightEffFastsimSFHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffFastsimSFHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronVeto_FastsimScaleFactor");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
    // In 2015 it was not necessary to apply any reco tracking scale factors
    // (this was introduced in 2016 to adjust for the impact of the 'HIP' effect)
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    std::cout << "RazorHelper: loading muon efficiency histograms" << std::endl;
    // For 2015 the lepton efficiency histograms for tight and veto IDs in fullsim are both
    // stored in the fastsim-to-fullsim scale factor files. 
    muTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.root");
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.root");
    muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightMuonSelectionEffDenominatorReco_2015Final_Golden.root"); 
    vetoMuEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_VetoMuonSelectionEffDenominatorReco_2015Final_Golden.root"); 
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("MuonEff_Tight_Fullsim");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("MuonEff_Veto_Fullsim");
    muTightEffFastsimSFHist = (TH2D*)muTightEfficiencyFile->Get("MuonTight_FastsimScaleFactor");
    muVetoEffFastsimSFHist = (TH2D*)muVetoEfficiencyFile->Get("MuonVeto_FastsimScaleFactor");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
    // In 2015 it was not necessary to apply any reco tracking scale factors
    // (this was introduced in 2016 to adjust for the impact of the 'HIP' effect)
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadJECs_Razor2015() {
    std::cout << "RazorHelper: loading jet energy correction constants" << std::endl;
    // initialize
    std::string jecPathname = cmsswPath + "/src/RazorAnalyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    
    if (isData) {
      
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else if (isFastsim) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fastsim_MCRUN2_74_V9_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}


void RazorHelper::loadBTag_Razor2015() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/BTagEffFastsimToFullsimCorrectionFactors.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/CharmJetBTagEffFastsimToFullsimCorrectionFactors.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/LightJetBTagEffFastsimToFullsimCorrectionFactors.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("BTagEff_Medium_Fullsim");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("BTagEff_Medium_Fullsim");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("BTagEff_Medium_Fullsim");

    std::string bTagPathname = cmsswPath + "/src/RazorAnalyzer/data/ScaleFactors/";
    // Fullsim
    btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2.csv",bTagPathname.c_str()));
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "mujets",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", Form("%s/CSV_13TEV_Combined_20_11_2015.csv",bTagPathname.c_str()));
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central"); 
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");  
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");  
}

void RazorHelper::loadTrigger_Razor2015() {
    // single lepton trigger scale factors
    std::cout << "RazorHelper: loading trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015Final_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015Final_Golden.root"); 
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleElectronTriggerEfficiencyFromFullsim.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");
    muTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleMuonTriggerEfficiencyFromFullsim.root"); 
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = std::vector<int>(8);
    hadronicTriggerNums = std::vector<int>(11);
    if( isData ) {
        dileptonTriggerNums = { 41,43,30,31,47,48,49,50 };
        singleLeptonTriggerNums = std::vector<int>(13);
        singleLeptonTriggerNums = { 2,7,12,11,15,22,23,24,25,26,27,28,29 };
        hadronicTriggerNums = { 134,135,136,137,138,139,140,141,142,143,144 };
    }
    else {
        dileptonTriggerNums = { 41,43,30,31,47,48,49,50 };
        singleLeptonTriggerNums = std::vector<int>(11);
        singleLeptonTriggerNums = { 2,7,12,11,15,18,19,20,21,28,29 };
        hadronicTriggerNums = { 134,135,136,137,138,139,140,141,142,143,144 };
    }
}

////////////////////////////////////////////////
//  2015 76X ReReco
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2015_76X() {
    loadPileup_Razor2015_76X();
    loadLepton_Razor2015(); //we'll use the same here for now
    loadPhoton_Razor2015_76X();
    loadJECs_Razor2015_76X();   //
    loadBTag_Razor2015();   //we'll use the same here for now, but this needs to be updated
    loadTrigger_Razor2015();//use the same for now
}

void RazorHelper::loadPileup_Razor2015_76X() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open(
            Form("%s/src/RazorAnalyzer/data/PileupReweight2015_7_6.root", cmsswPath.c_str()));
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadPhoton_Razor2015_76X(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2015/efficiency_results_PhoLooseEffDenominatorReco_2015.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");    
}

void RazorHelper::loadJECs_Razor2015_76X() {
    std::cout << "RazorHelper: loading jet energy correction constants" << std::endl;
    // initialize
    // load JEC parameters
    std::string jecPathname = cmsswPath + "/src/RazorAnalyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();

    if (isData) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else if (isFastsim) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fastsim_MCRUN2_74_V9_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));  
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}


////////////////////////////////////////////////
//  2016 Moriond Rereco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2016_MoriondRereco() {
    loadPileup_Razor2016_MoriondRereco();
    loadLepton_Razor2016_MoriondRereco();
    loadPhoton_Razor2016_MoriondRereco();
    loadBTag_Razor2016_MoriondRereco();
    loadTrigger_Razor2016_MoriondRereco();
    loadJECs_Razor2016_MoriondRereco();
}


void RazorHelper::loadPileup_Razor2016_MoriondRereco() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open(
            Form("%s/src/RazorAnalyzer/data/PileupWeights/PileupReweight_2016_36p2ifb.root", cmsswPath.c_str()));
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016_MoriondRereco(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root"); 
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root"); 
    eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    eleGSFTrackEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiencySF_muEleTracking_2016_average.root");
    eleGSFTrackEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have ID scale factors for Fastsim yet.
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorGen");
    eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("h2_scaleFactorsEGamma");

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root"); 
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root"); 
    muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root"); 
    vetoMuEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root"); 
    muTrackEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiencySF_muEleTracking_2016_average.root");
    muTrackEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("Efficiency_PtEta"); 
    // We don't have ID scale factors for Fastsim yet.
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorGen");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}


void RazorHelper::loadPhoton_Razor2016_MoriondRereco(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoLooseEffDenominatorReco_2016_Rereco.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");    
}


void RazorHelper::loadBTag_Razor2016_MoriondRereco() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    std::string bTagPathname = cmsswPath + "/src/RazorAnalyzer/data/ScaleFactors/";
    // Fullsim
    btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2_ichep.csv",bTagPathname.c_str()));
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "mujets",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", Form("%s/CSV_13TEV_Combined_20_11_2015.csv",bTagPathname.c_str()));
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central"); 
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");  
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");  
}

void RazorHelper::loadTrigger_Razor2016_MoriondRereco() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_Rereco_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");

    muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_Rereco_Golden.root"); 
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");
    
    eleTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016_Rereco_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016_Rereco_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root"); 
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root"); 
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    diphotonTrigLeadingLegEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root"); 
    diphotonTrigLeadingLegEffSFHist = (TH2D*)diphotonTrigLeadingLegEffSFFile->Get("ScaleFactor_PhoHLTLeadingLegEffDenominatorLoose");
    diphotonTrigTrailingLegEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root"); 
    diphotonTrigTrailingLegEffSFHist = (TH2D*)diphotonTrigTrailingLegEffSFFile->Get("ScaleFactor_PhoHLTTrailingLegEffDenominatorLoose");

    //get trigger numbers
    // (we are using the same list of trigger numbers for data and MC for 2016)
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}

void RazorHelper::loadJECs_Razor2016_MoriondRereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Spring16_23Sep2016_V1." << std::endl;
    // initialize
    std::string jecPathname = cmsswPath + "/src/RazorAnalyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
 
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016BCDV1_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016BCDV1_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016BCDV1_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016BCDV1_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Spring16_23Sep2016BCDV1_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersE = std::vector<JetCorrectorParameters> ();
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016EV1_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016EV1_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016EV1_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016EV1_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersE = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorE = new FactorizedJetCorrector(correctionParametersE);
      std::string jecUncPathE = jecPathname+"/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016EV1_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncE = new JetCorrectionUncertainty(jecUncPathE);
      SimpleJetResolution* JetResolutionCalculatorE = new SimpleJetResolution(*JetResolutionParametersE);

      correctionParameters.push_back(correctionParametersE);
      JetResolutionParameters.push_back(JetResolutionParametersE);
      JetCorrector.push_back( JetCorrectorE );
      JetResolutionCalculator.push_back(JetResolutionCalculatorE);
      jecUnc.push_back(jecUncE);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 277420 ));

      //IOV: 2016F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016FV1_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016FV1_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016FV1_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016FV1_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      std::string jecUncPathF = jecPathname+"/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016FV1_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(jecUncPathF);
      SimpleJetResolution* JetResolutionCalculatorF = new SimpleJetResolution(*JetResolutionParametersF);

      correctionParameters.push_back(correctionParametersF);
      JetResolutionParameters.push_back(JetResolutionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorF);
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 277772, 276801 ));

      //IOV: 2016GH
      std::vector<JetCorrectorParameters> correctionParametersGH = std::vector<JetCorrectorParameters> ();
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016GHV1_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016GHV1_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016GHV1_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016GHV1_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersGH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorGH = new FactorizedJetCorrector(correctionParametersGH);
      std::string jecUncPathGH = jecPathname+"/Spring16_23Sep2016_V1_DATA/Spring16_23Sep2016GHV1_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncGH = new JetCorrectionUncertainty(jecUncPathGH);
      SimpleJetResolution* JetResolutionCalculatorGH = new SimpleJetResolution(*JetResolutionParametersGH);

      correctionParameters.push_back(correctionParametersGH);
      JetResolutionParameters.push_back(JetResolutionParametersGH);
      JetCorrector.push_back( JetCorrectorGH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorGH);
      jecUnc.push_back(jecUncGH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 99999999 ));

    }
    else if (isFastsim) {
      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
   }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_MC/Spring16_23Sep2016V1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_MC/Spring16_23Sep2016V1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_23Sep2016_V1_MC/Spring16_23Sep2016V1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Spring16_23Sep2016_V1_MC/Spring16_23Sep2016V1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
  
}



////////////////////////////////////////////////
//  2016 Prompt Reco
////////////////////////////////////////////////

void RazorHelper::loadPhoton_Razor2016(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoLooseEffDenominatorReco_2016_ICHEP.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");    
}


void RazorHelper::loadBTag_Razor2016() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    std::string bTagPathname = cmsswPath + "/src/RazorAnalyzer/data/ScaleFactors/";
    // Fullsim
    btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2_ichep.csv",bTagPathname.c_str()));
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "mujets",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", Form("%s/CSV_13TEV_Combined_20_11_2015.csv",bTagPathname.c_str()));
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central"); 
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");  
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");  
}

void RazorHelper::loadJECs_Razor2016() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Spring16_25nsV6." << std::endl;
    // initialize
    std::string jecPathname = cmsswPath + "/src/RazorAnalyzer/data/JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();

    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersE = std::vector<JetCorrectorParameters> ();
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersE = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorE = new FactorizedJetCorrector(correctionParametersE);
      std::string jecUncPathE = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncE = new JetCorrectionUncertainty(jecUncPathE);
      SimpleJetResolution* JetResolutionCalculatorE = new SimpleJetResolution(*JetResolutionParametersE);

      correctionParameters.push_back(correctionParametersE);
      JetResolutionParameters.push_back(JetResolutionParametersE);
      JetCorrector.push_back( JetCorrectorE );
      JetResolutionCalculator.push_back(JetResolutionCalculatorE);
      jecUnc.push_back(jecUncE);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 277420 ));

      //IOV: 2016F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      std::string jecUncPathF = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(jecUncPathF);
      SimpleJetResolution* JetResolutionCalculatorF = new SimpleJetResolution(*JetResolutionParametersF);

      correctionParameters.push_back(correctionParametersF);
      JetResolutionParameters.push_back(JetResolutionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorF);
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 277772, 276801 ));

      //IOV: 2016GH
      std::vector<JetCorrectorParameters> correctionParametersGH = std::vector<JetCorrectorParameters> ();
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersGH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorGH = new FactorizedJetCorrector(correctionParametersGH);
      std::string jecUncPathGH = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncGH = new JetCorrectionUncertainty(jecUncPathGH);
      SimpleJetResolution* JetResolutionCalculatorGH = new SimpleJetResolution(*JetResolutionParametersGH);

      correctionParameters.push_back(correctionParametersGH);
      JetResolutionParameters.push_back(JetResolutionParametersGH);
      JetCorrector.push_back( JetCorrectorGH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorGH);
      jecUnc.push_back(jecUncGH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 99999999 ));

    }
    else if (isFastsim) {
       std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
 
}

////////////////////////////////////////////////
//  2016 ICHEP
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016_ICHEP_80X() {
    loadPileup_Razor2016_ICHEP();
    loadLepton_Razor2016_ICHEP();
    loadPhoton_Razor2016(); // same as 2016 inclusive
    loadBTag_Razor2016(); // same as 2016 inclusive
    loadTrigger_Razor2016_ICHEP();
    loadJECs_Razor2016(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016_ICHEP() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open(
            Form("%s/src/RazorAnalyzer/data/PileupWeights/PileupReweight2016_ICHEP.root", cmsswPath.c_str()));
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016_ICHEP(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root"); 
    eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016_ICHEP_Golden.root");
    eleGSFTrackEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden_12p9_ICHEP/efficiencySF_GsfTracking_2016ICHEP_Golden.root");
    eleGSFTrackEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    // We did not produce veto electron efficiencies or scale factors for this run period
    eleVetoEfficiencyFile = 0;
    eleVetoEfficiencyHist = 0;
    vetoEleEffSFFile = 0;
    eleVetoEffSFHist = 0;
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("EGamma_SF2D");
    eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root"); 
    muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016_ICHEP_Golden.root"); 
    muTrackEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden_12p9_ICHEP/efficiencySF_TrackReconstruction_2016ICHEP.root");
    muTrackEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");

    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    // We did not produce veto electron efficiencies or scale factors for this run period
    muVetoEfficiencyFile = 0;
    muVetoEfficiencyHist = 0; 
    vetoMuEffSFFile = 0;
    muVetoEffSFHist = 0;
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadTrigger_Razor2016_ICHEP() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_ICHEP_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_ICHEP_Golden.root"); 
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016_ICHEP_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016_ICHEP_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}



////////////////////////////////////////////////
//  2016 G
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016G_80X() {
    loadPileup_Razor2016G();
    loadLepton_Razor2016G();
    loadPhoton_Razor2016(); // same as 2016 inclusive
    loadBTag_Razor2016(); // same as 2016 inclusive
    loadTrigger_Razor2016G();
    loadJECs_Razor2016(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016G() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupWeights/PileupReweight2016G.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016G(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root"); 
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root"); 
    eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016G_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorReco_2016G_Golden.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root"); 
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root"); 
    muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016G_Golden.root"); 
    vetoMuEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorReco_2016G_Golden.root"); 
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("Efficiency_PtEta"); 
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadTrigger_Razor2016G() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016G_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016G_Golden.root"); 
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016G_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016G_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}

////////////////////////////////////////////////
//  2016 G Unblinded Data
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016G_SUSYUnblind_80X() {
    loadPileup_Razor2016G_SUSYUnblind();
    loadLepton_Razor2016G_SUSYUnblind();
    loadPhoton_Razor2016(); // same as 2016 inclusive
    loadBTag_Razor2016(); // same as 2016 inclusive
    loadTrigger_Razor2016G_SUSYUnblind();
    loadJECs_Razor2016(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016G_SUSYUnblind() {
    // pileup weights
    // LAST UPDATED: 3 November 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupWeights/PileupReweight2016G_SUSYUnblind.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016G_SUSYUnblind(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 3 November 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root"); 
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root"); 
    eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root"); 
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root"); 
    muEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root"); 
    vetoMuEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorGen_2016G_SUSYUnblind_Golden.root"); 
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("Efficiency_PtEta"); 
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadTrigger_Razor2016G_SUSYUnblind() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016G_Rereco_SUSYUnblind_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016G_SUSYUnblind_Golden.root"); 
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016G_Rereco_SUSYUnblind_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016G_SUSYUnblind_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root"); 
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}


////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////

double RazorHelper::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

// Get scale factor from a histogram with pt on the x-axis and eta on the y-axis
double RazorHelper::lookupPtEtaScaleFactor(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetXaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetXaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinContent(
                hist->GetXaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) ),
                hist->GetYaxis()->FindFixBin( etaToUse )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}

// Get scale factor uncertainty for a histogram with pt on the x-axis and eta on the y-axis
double RazorHelper::lookupPtEtaScaleFactorError(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetXaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetXaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinError(
                hist->GetXaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) ),
                hist->GetYaxis()->FindFixBin( etaToUse )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}

// Get scale factor from histogram with eta on the x-axis and pt on the y-axis
double RazorHelper::lookupEtaPtScaleFactor(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetYaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetYaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinContent(
                hist->GetXaxis()->FindFixBin( etaToUse ),
                hist->GetYaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}

// Gets the correct event-level scale factor depending on whether the object passes or fails selection
double RazorHelper::getPassOrFailScaleFactor(double eff, double sf, bool passes) {
  if (passes) {
    return sf;
  }
  else if (eff*sf >= 1.0 || eff >= 1.0) { //provide safety against infinite or negative weights
    return 1.0;
  }
  else if (eff == 0 || sf == 0) { //provide safety against infinite or negative weights
    return 0.0;
  }
  return (1/eff - sf) / (1/eff - 1);
}

// Helper function for computing scale factors from histograms
// The smear parameter is used to add a fractional uncertainty in quadrature with the one from the histogram
// Returns a vector of five scale factors: 
// 1) central value of scale factor
// 2,3) scale factor plus/minus uncertainty from fullsim scale factor
// 4,5) scale factor plus/minus uncertainty from fastsim scale factor
std::vector<double> RazorHelper::getLeptonScaleFactors(TH2D *effHist, TH2D *sfHist, 
        TH2D *fastsimHist, double pt, double eta, bool passes, double smear) {
    double eff = lookupPtEtaScaleFactor( effHist, pt, eta );
    double eff_fastsimUp = eff;
    double eff_fastsimDown = eff;
    if (isFastsim) { //correct efficiency for Fastsim
        double sf = lookupPtEtaScaleFactor( fastsimHist, pt, eta );
        double sfErr = lookupPtEtaScaleFactorError( fastsimHist, pt, eta );
        eff *= sf; 
        eff_fastsimUp *= (sf + sfErr);
        eff_fastsimDown *= (sf - sfErr);
    }
    double effSF = lookupPtEtaScaleFactor( sfHist, pt, eta );
    double effSFErr = lookupPtEtaScaleFactorError( sfHist, pt, eta );
    // add smearing factor to error
    effSFErr = sqrt(effSFErr*effSFErr + (smear*effSF)*(smear*effSF));
    double effSFUp = effSF + effSFErr;
    double effSFDown = effSF - effSFErr;

    double tmpSF = getPassOrFailScaleFactor( eff, effSF, passes );
    double tmpSFUp = getPassOrFailScaleFactor( eff, effSFUp, passes );
    double tmpSFDown = getPassOrFailScaleFactor( eff, effSFDown, passes );
    double tmpSF_fastsimUp = getPassOrFailScaleFactor( eff_fastsimUp, effSF, passes );
    double tmpSF_fastsimDown = getPassOrFailScaleFactor( eff_fastsimDown, effSF, passes );

    std::vector<double> out { tmpSF, tmpSFUp, tmpSFDown, tmpSF_fastsimUp, tmpSF_fastsimDown };
    return out;
}

// This function updates the provided scale factor values in-place by multiplying them by the 
// SFs corresponding to the provided values of pt and eta.  
void RazorHelper::updateScaleFactors(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist, 
        float pt, float eta, bool passes, float &sf, float &sfUp, float &sfDown, 
        float &sfFastsimUp, float &sfFastsimDown, float smear) {

    //retrieve scale factors and deal correctly with passing/failing selection
    std::vector<double> scaleFactors = getLeptonScaleFactors( effHist, sfHist, fastsimHist, pt, eta, passes, smear );

    //propagate to input values
    sf *= scaleFactors[0];
    sfUp *= scaleFactors[1]/scaleFactors[0];
    sfDown *= scaleFactors[2]/scaleFactors[0];
    sfFastsimUp *= scaleFactors[3]/scaleFactors[0];
    sfFastsimDown *= scaleFactors[4]/scaleFactors[0];
}

void RazorHelper::updateTightMuonScaleFactors(float pt, float eta, bool isTight, 
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( muTightEfficiencyHist, muTightEffSFHist, muTightEffFastsimSFHist, 
            pt, eta, isTight, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown );

}

void RazorHelper::updateVetoMuonScaleFactors(float pt, float eta, bool isVeto, 
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {
    updateScaleFactors( muVetoEfficiencyHist, muVetoEffSFHist, muVetoEffFastsimSFHist, 
            fmax( pt, muVetoEffSFMinPt), eta, isVeto, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown );

}

void RazorHelper::updateTightElectronScaleFactors(float pt, float eta, bool isTight, 
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( eleTightEfficiencyHist, eleTightEffSFHist, eleTightEffFastsimSFHist, 
            pt, eta, isTight, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown, 0.02 );

}

void RazorHelper::updateVetoElectronScaleFactors(float pt, float eta, bool isVeto, 
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {
        
    updateScaleFactors( eleVetoEfficiencyHist, eleVetoEffSFHist, eleVetoEffFastsimSFHist, 
            fmax( pt, eleVetoEffSFMinPt), eta, isVeto, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown, 0.02 );

}

// Computes a single lepton scale factor
double RazorHelper::getLeptonScaleFactor(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist, 
        double pt, double eta, bool passes) {
    double eff = lookupPtEtaScaleFactor( effHist, pt, eta );
    if (isFastsim) { //correct efficiency for Fastsim
        double sf = lookupPtEtaScaleFactor( fastsimHist, pt, eta );
        eff *= sf; 
    }
    double effSF = lookupPtEtaScaleFactor( sfHist, pt, eta );

    return getPassOrFailScaleFactor( eff, effSF, passes );
}

double RazorHelper::getTightMuonScaleFactor(float pt, float eta, bool isTight) {
    return getLeptonScaleFactor( muTightEfficiencyHist, muTightEffSFHist, muTightEffFastsimSFHist,
            pt, eta, isTight );
}

double RazorHelper::getVetoMuonScaleFactor(float pt, float eta, bool isVeto) {
    return getLeptonScaleFactor( muVetoEfficiencyHist, muVetoEffSFHist, muVetoEffFastsimSFHist,
            fmax( pt, muVetoEffSFMinPt), eta, isVeto );
}

double RazorHelper::getMuonTrackScaleFactor(float pt, float eta, bool isReconstructed) {
    double eff = lookupPtEtaScaleFactor( muTrackEffHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    double sf = lookupPtEtaScaleFactor( muTrackEffSFHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    return getPassOrFailScaleFactor( eff, sf, isReconstructed );
}

double RazorHelper::getTightElectronScaleFactor(float pt, float eta, bool isTight) {
    return getLeptonScaleFactor( eleTightEfficiencyHist, eleTightEffSFHist, eleTightEffFastsimSFHist,
            pt, eta, isTight );
}

double RazorHelper::getVetoElectronScaleFactor(float pt, float eta, bool isVeto) {
    return getLeptonScaleFactor( eleVetoEfficiencyHist, eleVetoEffSFHist, eleVetoEffFastsimSFHist,
            fmax( pt, eleVetoEffSFMinPt), eta, isVeto );
}

double RazorHelper::getEleGSFTrackScaleFactor(float pt, float eta, bool isReconstructed) {
    double eff = lookupPtEtaScaleFactor( eleGSFTrackEffHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    // note that the electron GSF tracking efficiency histogram has eta on the x-axis
    double sf = lookupEtaPtScaleFactor( eleGSFTrackEffSFHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    return getPassOrFailScaleFactor( eff, sf, isReconstructed );
}


double RazorHelper::getPhotonScaleFactor(float pt, float eta) {
  double sf = 1.0;
  if (phoLooseEffSFHist) sf = lookupPtEtaScaleFactor( phoLooseEffSFHist, pt, eta, 20.01, 99.9 ); 
  else { std::cout << "[WARNING] Could not load phoLooseEffSFHist.\n"; }
  return sf;
}

// Returns the trigger scale factor corresponding to the given values of pt and eta.  
// For fastsim events, the scale factor will be multiplied by the trigger efficiency 
// given in fastsimHist.
double RazorHelper::getTriggerScaleFactor(TH2D *sfHist, TH2D *fastsimHist, float pt, float eta, 
        bool isTight, bool passedTrigger, float fastsimPtCut, float ptCut) {
    double trigSF = lookupPtEtaScaleFactor( sfHist, pt, eta, 199.9, ptCut );
    if (isFastsim) {
        // note that the fastsim trigger histograms have pt on the y-axis
        trigSF *= lookupEtaPtScaleFactor( fastsimHist, pt, eta, 999.9, fastsimPtCut );
    }
    if (passedTrigger && isTight){
        return trigSF;
    }
    return 1.0;
}

double RazorHelper::getTriggerEfficiency(TH2D *effHist, float pt, float eta, 
        bool isTight, bool passedTrigger, float ptCut) {
    double trigEff = lookupEtaPtScaleFactor( effHist, pt, eta, 199.9, ptCut );
    if (passedTrigger && isTight){
        return trigEff;
    }
    return 1.0;
}

double RazorHelper::getSingleMuTriggerScaleFactor(float pt, float eta, bool isTight, bool passedTrigger) {
    return getTriggerScaleFactor( muTrigSFHist, muTrigEffHist, pt, eta, 
            isTight, passedTrigger, 15.01, 20.01 );
}

double RazorHelper::getSingleMuTriggerEfficiency(float pt, float eta, bool isTight, bool passedTrigger) {
    if (muTrigEffHist) {
        return getTriggerEfficiency(muTrigEffHist, pt, eta, isTight, passedTrigger, 20.01);
    }
    std::cout << "[WARNING] Single muon trigger efficiency requested, but no histogram is available.  Returning 0" << std::endl;
    return 0;
}

double RazorHelper::getSingleEleTriggerEfficiency(float pt, float eta, bool isTight, bool passedTrigger) {
    if (eleTrigEffHist) {
        return getTriggerEfficiency(eleTrigEffHist, pt, eta, isTight, passedTrigger, 25.01);
    }
    std::cout << "[WARNING] Single electron trigger efficiency requested, but no histogram is available.  Returning 0" << std::endl;
    return 0;
}

double RazorHelper::getSingleEleTriggerScaleFactor(float pt, float eta, bool isTight, bool passedTrigger) {
    return getTriggerScaleFactor( eleTrigSFHist, eleTrigEffHist, pt, eta, 
            isTight, passedTrigger, 25.01, 25.01 );
}

// Update the provided trigger scale factors in-place by multiplying them by the SFs 
// corresponding to the given values of pt and eta
void RazorHelper::updateTriggerScaleFactors(TH2D *sfHist, TH2D *fastsimHist, 
        float pt, float eta, bool isTight, bool passedTrigger, float &sf, float &sfUp, float &sfDown, 
        float fastsimPtCut, float extraSyst) {
    double trigSF = lookupPtEtaScaleFactor( sfHist, pt, eta );
    double trigSFErr = lookupPtEtaScaleFactorError( sfHist, pt, eta );
    double trigSFUp = trigSF + trigSFErr + extraSyst;
    double trigSFDown = trigSF - trigSFErr - extraSyst;
    if (passedTrigger && isTight){
        sf *= trigSF;
        sfUp *= trigSFUp/trigSF;
        sfDown *= trigSFDown/trigSF;
    }
    if (isFastsim) {
        if (passedTrigger && isTight) {
            // note that the fastsim trigger histograms have pt on the y-axis
            double trigSFFastsim = lookupEtaPtScaleFactor( fastsimHist, pt, eta, 999.9, fastsimPtCut );
            sf *= trigSFFastsim;
        }		  
    }
}

void RazorHelper::updateSingleMuTriggerScaleFactors(float pt, float eta, bool isTight, 
        bool passedTrigger, float &sf, float &sfUp, float &sfDown) {
    updateTriggerScaleFactors( muTrigSFHist, muTrigEffHist, pt, eta, isTight, 
            passedTrigger, sf, sfUp, sfDown, 15.01 );
}

void RazorHelper::updateSingleEleTriggerScaleFactors(float pt, float eta, bool isTight, 
        bool passedTrigger, float &sf, float &sfUp, float &sfDown) {
    updateTriggerScaleFactors( eleTrigSFHist, eleTrigEffHist, pt, eta, isTight, 
            passedTrigger, sf, sfUp, sfDown, 25.01, 0.02 );
}

double RazorHelper::getDiphotonTrigLeadingLegEff(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigLeadingLegEffHist) sf = lookupEtaPtScaleFactor( diphotonTrigLeadingLegEffHist, pt, eta, 20.01, 99.9 ); 
  return sf; 
}

double RazorHelper::getDiphotonTrigLeadingLegEffSF(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigLeadingLegEffSFHist) sf = lookupPtEtaScaleFactor( diphotonTrigLeadingLegEffSFHist, pt, eta, 20.01, 99.9 ); 
  return sf; 
}

double RazorHelper::getDiphotonTrigTrailingLegEff(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigTrailingLegEffHist) sf = lookupEtaPtScaleFactor( diphotonTrigTrailingLegEffHist, pt, eta, 20.01, 99.9 ); 
  return sf; 
}

double RazorHelper::getDiphotonTrigTrailingLegEffSF(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigTrailingLegEffSFHist) sf = lookupPtEtaScaleFactor( diphotonTrigTrailingLegEffSFHist, pt, eta, 20.01, 99.9 ); 
  return sf; 
}

// Retrieve jet energy uncertainty as a function of pt and eta
double RazorHelper::getJecUnc( float pt, float eta , int run) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}

// Get one b-tag scale factor
double RazorHelper::getBTagScaleFactor(float pt, float eta, int flavor, bool isCSVM) {
    // Get efficiency and jet flavor
    double effMedium = 0;
    BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
    if ( abs(flavor) == 5) {
        effMedium = lookupPtEtaScaleFactor( btagMediumEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_B;
    } 
    else if ( abs(flavor) == 4) {
        effMedium = lookupPtEtaScaleFactor( btagMediumCharmEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_C;
    } 
    else {
        effMedium = lookupPtEtaScaleFactor( btagMediumLightJetsEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_UDSG;
    }

    // Get scale factor
    if (pt >= 670) pt = 669; //670 is the largest pt range listed in the CSV text file
    double jet_scalefactor = -1;
    if ( abs(flavor) == 5 || abs(flavor) == 4 ) {
        jet_scalefactor = btagreader->eval(jetType, eta, pt); 
    }
    else {
        jet_scalefactor = 0.907317;
    }

    //correct efficiency for Fastsim
    //Do this only for b-jets for now
    if (isFastsim && abs(flavor) == 5) { 
        jet_scalefactor *= btagreaderfastsim->eval(jetType, eta, pt);
    }

    return getPassOrFailScaleFactor( effMedium, jet_scalefactor, isCSVM );
}

// Update the provided b-tag scale factors in-place by multiplying by the SFs corresponding
// to the given values of pt, eta, and flavor.  
void RazorHelper::updateBTagScaleFactors(float pt, float eta, int flavor, bool isCSVM,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown, 
        float &sfMistagUp, float &sfMistagDown) {
    // Get efficiency and jet flavor
    double effMedium = 0;
    BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
    if ( abs(flavor) == 5) {
        effMedium = lookupPtEtaScaleFactor( btagMediumEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_B;
    } 
    else if ( abs(flavor) == 4) {
        effMedium = lookupPtEtaScaleFactor( btagMediumCharmEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_C;
    } 
    else {
        effMedium = lookupPtEtaScaleFactor( btagMediumLightJetsEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_UDSG;
    }

    // Get scale factor
    if (pt >= 670) pt = 669; //670 is the largest pt range listed in the CSV text file
    double jet_scalefactor = -1;
    double jet_scalefactorUp = -1;
    double jet_scalefactorDown = -1;  
    if ( abs(flavor) == 5 || abs(flavor) == 4 ) {
        jet_scalefactor = btagreader->eval(jetType, eta, pt); 
        jet_scalefactorUp = btagreader_up->eval(jetType, eta, pt);
        jet_scalefactorDown = btagreader_do->eval(jetType, eta, pt);
    }
    else {
        jet_scalefactor = 0.907317;
        jet_scalefactorUp = 1.257317;		  
        jet_scalefactorDown = 0.557317;		
    }
    double jet_scalefactorFastsimUp = jet_scalefactor;
    double jet_scalefactorFastsimDown = jet_scalefactor;

    //correct efficiency for Fastsim
    //Do this only for b-jets for now
    if (isFastsim && abs(flavor) == 5) { 
        double jet_scalefactorFastsim = btagreaderfastsim->eval(jetType, eta, pt);
        jet_scalefactor *= jet_scalefactorFastsim;
        jet_scalefactorUp *= jet_scalefactorFastsim;
        jet_scalefactorDown *= jet_scalefactorFastsim;
        jet_scalefactorFastsimUp *= btagreaderfastsim_up->eval(jetType, eta, pt);
        jet_scalefactorFastsimDown *= btagreaderfastsim_do->eval(jetType, eta, pt);
    }

    //apply and propagate scale factor
    double tmpSF = getPassOrFailScaleFactor( effMedium, jet_scalefactor, isCSVM );
    sf *= tmpSF;
    if (abs(flavor) == 5 || abs(flavor) == 4) {
        sfUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorUp, isCSVM ) / tmpSF;
        sfDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorDown, isCSVM ) / tmpSF;
        sfFastsimUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorFastsimUp, isCSVM ) / tmpSF;
        sfFastsimDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorFastsimDown, isCSVM ) / tmpSF;
    } 
    else {
        sfMistagUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorUp, isCSVM ) / tmpSF;
        sfMistagDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorDown, isCSVM ) / tmpSF;
    }
}

// top pt reweighting from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
float RazorHelper::getTopPtWeight( float ptT, float ptTbar ) {
    // do not extrapolate correction beyond pt = 400 GeV
    if( ptT > 400 ) ptT = 400;
    if( ptTbar > 400 ) ptTbar = 400;

    float a = 0.156; // from 8 TeV, still recommended for 13 TeV use
    float b = -0.00137; //from 8 TeV, still recommended for 13 TeV use

    // weight from top
    float weightT = exp( a + b*ptT );
    // weight from antitop
    float weightTbar = exp( a + b*ptTbar );

    // combine into total weight
    return sqrt( weightT * weightTbar );
}

// electron scale corrections, derived privately on 2016 data
float RazorHelper::getElectronScaleCorrection( float eta ) {
    if ( fabs(eta) < 0.4 ) {
        return 1/0.993387;
    }
    else if ( fabs(eta) < 0.8 ) {
        return 1/0.993516;
    }
    else if ( fabs(eta) < 1.4442 ) {
        return 1/0.990877;
    }
    else {
        return 1/0.998137;
    }
}

float RazorHelper::getElectronResCorrection( float eta ) {
    if ( fabs(eta) < 0.4 ) {
        return 0.0137622;
    }
    else if ( fabs(eta) < 0.8 ) {
        return 0.961012;
    }
    else if ( fabs(eta) < 1.4442 ) {
        return 1.01832;
    }
    else {
        return 1.29274;
    }
}

float RazorHelper::getCorrectedElectronPt( float pt, float eta ) {
    return gRandom->Gaus( pt*getElectronScaleCorrection(eta), getElectronResCorrection(eta) );
}
