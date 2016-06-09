// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef RazorHelper_H
#define RazorHelper_H

#include <iostream>
#include <string>
#include <sys/stat.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"

class RazorHelper {

    public:
        // Constructor takes a string specifying which set of files to load
        RazorHelper(std::string tag_, bool isData_, bool isFastsim_); 
        virtual ~RazorHelper();

        // for pileup
        double getPileupWeight(int NPU);
        double getPileupWeightUp(int NPU);
        double getPileupWeightDown(int NPU);

        // for lepton
        void processTightMuonScaleFactors(float pt, float eta, bool isTight, 
            float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown);
        void processVetoMuonScaleFactors(float pt, float eta, bool isVeto, 
            float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown);
        void processTightElectronScaleFactors(float pt, float eta, bool isTight, 
            float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown);
        void processVetoElectronScaleFactors(float pt, float eta, bool isVeto, 
            float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown);

        // for trigger
        std::vector<int> getDileptonTriggerNums() { return dileptonTriggerNums; }
        std::vector<int> getSingleLeptonTriggerNums() { return singleLeptonTriggerNums; }
        std::vector<int> getHadronicTriggerNums() { return hadronicTriggerNums; }
        void processSingleMuTriggerScaleFactors(float pt, float eta, bool isTight, bool passedTrigger,
            float &sf, float &sfUp, float &sfDown);
        void processSingleEleTriggerScaleFactors(float pt, float eta, bool isTight, bool passedTrigger,
            float &sf, float &sfUp, float &sfDown);

        // for JECs
        FactorizedJetCorrector *getJetCorrector() { return JetCorrector; }
        double getJecUnc( float pt, float eta );

        // for btag 
        void processBTagScaleFactors(float pt, float eta, int flavor, bool isCSVM,
                float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown, 
                float &sfMistagUp, float &sfMistagDown);

    private:
        // member functions
        void loadTag_Razor2015(); // Final set of files used in 2015
        void loadTag_Razor2015_76X(); // Configuration for 2015 ReReco 
        void loadTag_Razor2016_80X(); // Evolving configuration for 2016 PromptReco
        void loadTag_Null(); // Default when tag is not provided
        void loadCMSSWPath();
        double lookupPtEtaScaleFactor(TH2D *hist, double pt, double eta, double ptmin=10.01, double ptmax=199.9);
        double lookupPtEtaScaleFactorError(TH2D *hist, double pt, double eta, double ptmin=10.01, double ptmax=199.9);
        double lookupEtaPtScaleFactor(TH2D *hist, double pt, double eta, double ptmin=10.01, double ptmax=199.9);
        double getPassOrFailScaleFactor(double eff, double sf, bool passes);
        std::vector<double> getLeptonScaleFactors(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist, 
                double pt, double eta, bool passes, double smear=0.0);
        void processScaleFactors(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist, float pt, 
                float eta, bool passes, float &sf, float &sfUp, float &sfDown, 
                float &sfFastsimUp, float &sfFastsimDown, float smear=0.0);
        void processTriggerScaleFactors(TH2D *sfHist, TH2D *fastsimHist, 
            float pt, float eta, bool isTight, bool passedTrigger, float &sf, float &sfUp, float &sfDown, 
            float fastsimPtCut = 10.01, float extraSyst = 0.);

        // for Razor2015 74X tag
        void loadPileup_Razor2015();
        void loadLepton_Razor2015();
        void loadTrigger_Razor2015();
        void loadJECs_Razor2015();
        void loadBTag_Razor2015();

        // for Razor2015 76X ReReco tag
        void loadPileup_Razor2015_76X();
        void loadLepton_Razor2015_76X();
        void loadTrigger_Razor2015_76X();
        void loadJECs_Razor2015_76X();
        void loadBTag_Razor2015_76X();

        // for Razor2016 80X tag
        void loadTrigger_Razor2016();
         
        // member data
        std::string tag;
        bool isData;
        bool isFastsim;
        std::string cmsswPath;

        // for pileup reweighting
        TFile *pileupWeightFile;
        TH1F *pileupWeightHist;
        TH1F *pileupWeightSysUpHist;
        TH1F *pileupWeightSysDownHist;

        // for electrons
        TFile *eleEfficiencyFile;
        TFile *eleEffSFFile;
        TFile *vetoEleEffSFFile;
        TH2D *eleTightEfficiencyHist;
        TH2D *eleVetoEfficiencyHist;
        TH2D *eleTightEffFastsimSFHist;
        TH2D *eleVetoEffFastsimSFHist;
        TH2D *eleTightEffSFHist;
        TH2D *eleVetoEffSFHist;

        // for muons
        TFile *muEfficiencyFile;
        TFile *muEffSFFile;
        TFile *vetoMuEffSFFile;
        TH2D *muTightEfficiencyHist;
        TH2D *muVetoEfficiencyHist;
        TH2D *muTightEffFastsimSFHist;
        TH2D *muVetoEffFastsimSFHist;
        TH2D *muTightEffSFHist;
        TH2D *muVetoEffSFHist;

        // for taus
        TFile *tauEfficiencyFile;
        TH2D *tauLooseEfficiencyHist;

        // for triggers
        TFile *eleTrigSFFile;
        TFile *muTrigSFFile;
        TFile *eleTrigEffFromFullsimFile;
        TFile *muTrigEffFromFullsimFile;
        TH2D *eleTrigSFHist;
        TH2D *muTrigSFHist;
        TH2D *eleTrigEffFromFullsimHist;
        TH2D *muTrigEffFromFullsimHist;
        std::vector<int> dileptonTriggerNums;
        std::vector<int> singleLeptonTriggerNums;
        std::vector<int> hadronicTriggerNums;

        // for jet energy corrections
        FactorizedJetCorrector *JetCorrector;
        JetCorrectionUncertainty *jecUnc;

        // for b-tag
        TFile *btagEfficiencyFile;
        TFile *btagCharmEfficiencyFile;
        TFile *btagLightJetsEfficiencyFile;
        TH2D *btagMediumEfficiencyHist;
        TH2D *btagMediumCharmEfficiencyHist;
        TH2D *btagMediumLightJetsEfficiencyHist;
        BTagCalibration *btagcalib;
        BTagCalibration *btagcalibfastsim;
        BTagCalibrationReader *btagreader;
        BTagCalibrationReader *btagreader_up;
        BTagCalibrationReader *btagreader_do;
        BTagCalibrationReader *btagreaderMistag;
        BTagCalibrationReader *btagreaderMistag_up;
        BTagCalibrationReader *btagreaderMistag_do;
        BTagCalibrationReader *btagreaderfastsim;
        BTagCalibrationReader *btagreaderfastsim_up;
        BTagCalibrationReader *btagreaderfastsim_do;
};

#endif
