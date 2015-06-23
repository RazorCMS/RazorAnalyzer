//Macro to plot yields in data and MC in the razor control regions and compare with the prediction from the fit method

#include <iostream>
#include <map>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "assert.h"
#include "math.h"

using namespace std;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

//compute b-tagging scale factor
double getBTagMediumScaleFactor(double b1Pt, bool b1Pass, double b2Pt, bool b2Pass){
    double btagScaleFactor = 1.0;
    double bjet1EventScaleFactor = 1.0;
    double bjet2EventScaleFactor = 1.0;
    if (b1Pt > 20) {
        double bjet1ScaleFactor = 0.938887 + 0.00017124 * b1Pt + (-2.76366e-07) * b1Pt * b1Pt ;
        double MCEff = 1.0;
        if (b1Pt < 50) MCEff = 0.65;
        else if (b1Pt < 80) MCEff = 0.70;
        else if (b1Pt < 120) MCEff = 0.73;
        else if (b1Pt < 210) MCEff = 0.73;
        else MCEff = 0.66;				 
        if (b1Pass) bjet1EventScaleFactor = bjet1ScaleFactor;
        else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
    }
    if (b2Pt > 20) {
        double bjet2ScaleFactor = 0.938887 + 0.00017124 * b2Pt + (-2.76366e-07) * b2Pt * b2Pt ;
        double MCEff = 1.0;
        if (b2Pt < 50) MCEff = 0.65;
        else if (b2Pt < 80) MCEff = 0.70;
        else if (b2Pt < 120) MCEff = 0.73;
        else if (b2Pt < 210) MCEff = 0.73;
        else MCEff = 0.66;		 
        if (b2Pass) bjet2EventScaleFactor = bjet2ScaleFactor;
        else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
    }
    btagScaleFactor = bjet1EventScaleFactor * bjet1EventScaleFactor;
    return btagScaleFactor;
}


void CompareDataFitCRs(){
    gROOT->SetBatch();

    int lumiInData = 19700; //in pb
    int lumiInMC = 1; //luminosity used to normalize input MC ntuples

    //set color palette 
    const Int_t NCont = 101;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //for plots
    float nMRBins = 8;
    float nRsqBins = 7;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5};

    //selection cuts for each control region
    //TODO: fill in cuts
    map<string, string> cuts;
    cuts["TTBarSingleLepton"] = "";
    cuts["TTBarDilepton"] = "";
    cuts["WSingleLepton"] = "";
    cuts["ZLLDilepton"] = "";
    cuts["ZNuNuFromDY"] = "hlt_dimuon && nLooseMuons == 2 && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    cuts["ZNuNuFromW"] = "hlt_singlemu && nBTaggedJets == 0 && nTightMuons == 1 && nLooseMuons == 1 && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100";
    cuts["ZNuNuFromGamma"] = "hlt_photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1 && leadingPhotonPt > 80";

    //get input files
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    string baseDir = "root://eoscms://store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/";
    string fileSuffix = ".root";

    mcfiles["DYJets"] = new TFile(Form("%s/DYJetsRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["TTJets"] = new TFile(Form("%s/TTJetsRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["WJets"] = new TFile(Form("%s/WJetsRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["ZJetsNuNu"] = new TFile(Form("%s/ZJetsNuNuRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["SingleTop"] = new TFile(Form("%s/SingleTopRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["TTV"] = new TFile(Form("%s/TTVRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["VV"] = new TFile(Form("%s/VVRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["QCD"] = new TFile(Form("%s/QCDRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["GJets"] = new TFile(Form("%s/GJetsRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["VG"] = new TFile(Form("%s/VGRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    mcfiles["TTG"] = new TFile(Form("%s/TTGRun1%s", baseDir.c_str(), fileSuffix.c_str()));

    //data
    datafiles["HTMHT"] = new TFile(Form("%s/HTMHTRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    datafiles["SingleMuon"] = new TFile(Form("%s/SingleMuonRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    datafiles["SingleElectron"] = new TFile(Form("%s/SingleElectronRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    datafiles["DoubleMuon"] = new TFile(Form("%s/DoubleMuonRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    datafiles["DoubleElectron"] = new TFile(Form("%s/DoubleElectronRun1%s", baseDir.c_str(), fileSuffix.c_str()));
    datafiles["Photon"] = new TFile(Form("%s/PhotonRun1%s", baseDir.c_str(), fileSuffix.c_str()));

    //assign datasets to control regions
    map<string, vector<string> > controlRegionMC;
    controlRegionMC["TTBarSingleLepton"] = vector<string> {"TTJets", "WJets", "DYJets", "SingleTop", "QCD", "VV", "TTV"};
    controlRegionMC["TTBarDilepton"] = vector<string> {"TTJets", "DYJets", "WJets", "SingleTop", "QCD", "VV", "TTV"};
    controlRegionMC["WSingleLepton"] = vector<string> {"WJets", "TTJets", "DYJets", "SingleTop", "QCD", "VV", "TTV"};
    controlRegionMC["ZNuNuFromDY"] = vector<string> {"DYJets", "TTJets", "WJets", "SingleTop", "QCD", "VV", "TTV"};
    controlRegionMC["ZNuNuFromW"] = vector<string> {"WJets", "TTJets", "DYJets", "SingleTop", "QCD", "VV", "TTV"};
    controlRegionMC["ZNuNuFromGamma"] = vector<string> {"GJets", "QCD", "VG", "TTG"};

    map<string, vector<string> > controlRegionData;
    controlRegionData["TTBarSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["TTBarDilepton"] = vector<string> {"DoubleMuon", "DoubleElectron"};
    controlRegionData["WSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["ZNuNuFromDY"] = vector<string> {"DoubleMuon"};
    controlRegionData["ZNuNuFromW"] = vector<string> {"SingleMuon"};
    controlRegionData["ZNuNuFromGamma"] = vector<string> {"Photon"};

    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    float weight;
    float leadingMuonPt, leadingPhotonPt, recoZmass, subleadingMuonPt, mTLepMet;
    float leadingTightMuonPt;
    float bjet1Pt, bjet2Pt;
    bool bjet1PassMedium, bjet2PassMedium;
    bool hlt_singlemu, hlt_dimuon, hlt_photon, hlt_razor;
    bool passedHLTPhoton50, passedHLTPhoton75, passedHLTPhoton90, passedHLTPhoton135, passedHLTPhoton150;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
    float genZPt, genWPt;
    int nPU_mean;
    float MR, Rsq, MR_noZ, Rsq_noZ, MR_noW, Rsq_noW, MR_noPho, Rsq_noPho;
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("MR_noW", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        mctrees[file.first]->SetBranchStatus("MR_noZ", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noZ", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noZ", 1);
        mctrees[file.first]->SetBranchStatus("MR_noPho", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("numJets80", 1);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("leadingMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("subleadingMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        mctrees[file.first]->SetBranchStatus("recoZmass", 1);
        mctrees[file.first]->SetBranchStatus("mTLepMet", 1);
        mctrees[file.first]->SetBranchStatus("nPU_mean", 1);
        mctrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        mctrees[file.first]->SetBranchStatus("nLooseElectrons", 1);
        mctrees[file.first]->SetBranchStatus("nTightMuons", 1);
        mctrees[file.first]->SetBranchStatus("hlt_photon", 1);
        mctrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        mctrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        mctrees[file.first]->SetBranchStatus("hlt_razor", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("genZpt", 1);
        mctrees[file.first]->SetBranchStatus("genWpt", 1);
        mctrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        mctrees[file.first]->SetBranchStatus("bjet2Pt", 1);

        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("MR_noZ", &MR_noZ);
        mctrees[file.first]->SetBranchAddress("Rsq_noZ", &Rsq_noZ);
        mctrees[file.first]->SetBranchAddress("MR_noW", &MR_noW);
        mctrees[file.first]->SetBranchAddress("Rsq_noW", &Rsq_noW);
        mctrees[file.first]->SetBranchAddress("MR_noPho", &MR_noPho);
        mctrees[file.first]->SetBranchAddress("Rsq_noPho", &Rsq_noPho);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        mctrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        mctrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        mctrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
        mctrees[file.first]->SetBranchAddress("nPU_mean", &nPU_mean);
        mctrees[file.first]->SetBranchAddress("genZpt", &genZPt);
        mctrees[file.first]->SetBranchAddress("genWpt", &genWPt);
        mctrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        mctrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
        mctrees[file.first]->SetBranchAddress("hlt_singlemu", &hlt_singlemu);
        mctrees[file.first]->SetBranchAddress("hlt_dimuon", &hlt_dimuon);
        mctrees[file.first]->SetBranchAddress("hlt_razor", &hlt_razor);
        mctrees[file.first]->SetBranchAddress("hlt_photon", &hlt_photon);
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        datatrees[file.first]->SetBranchStatus("*", 0);

        datatrees[file.first]->SetBranchStatus("MR_noW", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        datatrees[file.first]->SetBranchStatus("MR_noZ", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noZ", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noZ", 1);
        datatrees[file.first]->SetBranchStatus("MR_noPho", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        datatrees[file.first]->SetBranchStatus("MR", 1);
        datatrees[file.first]->SetBranchStatus("Rsq", 1);
        datatrees[file.first]->SetBranchStatus("numJets80", 1);
        datatrees[file.first]->SetBranchStatus("leadingMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("subleadingMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        datatrees[file.first]->SetBranchStatus("recoZmass", 1);
        datatrees[file.first]->SetBranchStatus("mTLepMet", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton50", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton75", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton90", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton135", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton150", 1);
        datatrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        datatrees[file.first]->SetBranchStatus("nLooseElectrons", 1);
        datatrees[file.first]->SetBranchStatus("nTightMuons", 1);
        datatrees[file.first]->SetBranchStatus("hlt_photon", 1);
        datatrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        datatrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        datatrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        datatrees[file.first]->SetBranchStatus("hlt_razor", 1);
        datatrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        datatrees[file.first]->SetBranchStatus("bjet2Pt", 1);
        datatrees[file.first]->SetBranchStatus("Flag_HBHENoiseFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_CSCTightHaloFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_eeBadScFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); // enable

        datatrees[file.first]->SetBranchAddress("MR", &MR);
        datatrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        datatrees[file.first]->SetBranchAddress("MR_noZ", &MR_noZ);
        datatrees[file.first]->SetBranchAddress("Rsq_noZ", &Rsq_noZ);
        datatrees[file.first]->SetBranchAddress("MR_noW", &MR_noW);
        datatrees[file.first]->SetBranchAddress("Rsq_noW", &Rsq_noW);
        datatrees[file.first]->SetBranchAddress("MR_noPho", &MR_noPho);
        datatrees[file.first]->SetBranchAddress("Rsq_noPho", &Rsq_noPho);
        datatrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        datatrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        datatrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        datatrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton50", &passedHLTPhoton50);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton75", &passedHLTPhoton75);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton90", &passedHLTPhoton90);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton135", &passedHLTPhoton135);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton150", &passedHLTPhoton150);
        datatrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        datatrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
        datatrees[file.first]->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter); // enable 
        datatrees[file.first]->SetBranchAddress("hlt_singlemu", &hlt_singlemu);
        datatrees[file.first]->SetBranchAddress("hlt_dimuon", &hlt_dimuon);
        datatrees[file.first]->SetBranchAddress("hlt_razor", &hlt_razor);
        datatrees[file.first]->SetBranchAddress("hlt_photon", &hlt_photon);
    }
    //luminosities collected by the various photon triggers
    float lumi_HLTPhoton50  = 1.353e0 + 4.921e0 + 7.947e0 + 8.131e0;
    float lumi_HLTPhoton75  = 8.111e0 + 2.953e1 + 4.768e1 + 4.879e1;
    float lumi_HLTPhoton90  = 1.622e1 + 6.408e1 + 1.010e2 + 9.948e1;
    float lumi_HLTPhoton135 = 8.893e2 + 1.476e2 + 5.429e3 + 7.318e3;
    float lumi_HLTPhoton150 = 8.893e2 + 4.429e3 + 7.152e3 + 7.318e3;

    //load TTbar scale factor histograms
    TFile *SFFileTTJets = new TFile("data/ScaleFactors/Run1/TTBarSingleLeptonScaleFactors.root");
    TH2F *SFHistTTJets = (TH2F*)SFFileTTJets->Get("TTBarSingleLeptonScaleFactor");
    float SFmaxMRTTJets = SFHistTTJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqTTJets = SFHistTTJets->GetYaxis()->GetXmax() - 0.01;

    //load WJets scale factor histogram
    TFile *SFFileWJets = new TFile("data/ScaleFactors/Run1/WJetsSingleLeptonScaleFactors.root");
    TH2F *SFHistWJets = (TH2F*)SFFileWJets->Get("WJetsSingleLeptonScaleFactor");
    float SFmaxMRWJets = SFHistWJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqWJets = SFHistWJets->GetYaxis()->GetXmax() - 0.01;

    //load DYJets scale factor histogram
    TFile *SFFileDYJets = new TFile("data/ScaleFactors/Run1/ZToLLScaleFactors.root");
    TH2F *SFHistDYJets = (TH2F*)SFFileDYJets->Get("ZToLLDileptonScaleFactor");
    float SFmaxMRDYJets = SFHistDYJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqDYJets = SFHistDYJets->GetYaxis()->GetXmax() - 0.01;

    //load ZNuNu scale factor histograms
    TFile *SFFileZJetsNuNu = new TFile("data/ScaleFactors/Run1/ZInvisibleScaleFactorsRun1.root");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("DYJetsScaleFactors");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("WJetsScaleFactors");
    TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("GJetsScaleFactors");
    float SFmaxMRZJetsNuNu = SFHistZJetsNuNu->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqZJetsNuNu = SFHistZJetsNuNu->GetYaxis()->GetXmax() - 0.01;

    //load ZNuNu-->DYJets weighting factors
    TFile *ZNuNuToDYWeightFile = new TFile("data/ScaleFactors/Run1/ZNuNuToDYScaleFactorsRun1.root");
    TH2F *ZNuNuToDYWeightHist = (TH2F*)ZNuNuToDYWeightFile->Get("razormcDYJets");
    float maxMRZNuNuToDY = ZNuNuToDYWeightHist->GetXaxis()->GetXmax() - 1;
    float maxRsqZNuNuToDY = ZNuNuToDYWeightHist->GetYaxis()->GetXmax() - 0.01;

    //load muon efficiency scale factor histogram
    TFile muIdSFFile("data/ScaleFactors/MuonEfficiencies_ID_Run2012ReReco_53X.root");
    TFile muIsoSFFile("data/ScaleFactors/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
    //TODO: use muon efficiency scale factors
    float singleMuTriggerSF = 0.97;
    float doubleMuTriggerSF = 0.97;
    float doubleMuNormalizationSF = 0.97;

    //load pileup reweighting histogram
    TFile *pileupWeightFile = new TFile("data/Run1PileupWeights.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
    assert(pileupWeightHist);

    ///////////////////////////////////////////////////////////
    // Get MC distributions
    ///////////////////////////////////////////////////////////
    map<string, map<string, TH2F*> > razorHistosMC;
    map<string, map<string, TH2F*> > razorErrorHistosMC; //store sum(w^2*error(SF)^2)
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;

        //make histograms, and make TTreeFormulas for selection cuts
        map<string, TTreeFormula*> cutsFormulas;
        for(auto &cut : cuts){
            razorHistosMC[cut.first][tree.first] = new TH2F(Form("razormc%s%s", cut.first.c_str(), tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorErrorHistosMC[cut.first][tree.first] = new TH2F(Form("razorErrormc%s%s", tree.first.c_str(), cut.second.c_str()), "sum(w^2*error(SF)^2); MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosMC[cut.first][tree.first]->Sumw2();
            razorErrorHistosMC[cut.first][tree.first]->Sumw2();
            cutsFormulas[cut.first] = new TTreeFormula(Form("%s%sCutsFormula", cut.first.c_str(), tree.first.c_str()), cuts[cut.first].c_str(), tree.second);
            cutsFormulas[cut.first]->GetNdata();
        }

        //loop over entries
        uint nEntries = tree.second->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            float eventWeight = weight*lumiInData*1.0/lumiInMC;

            //PU reweighting
            eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPU_mean));
            float sysErrorSquared = 0.0;

            for(auto &cutf : cutsFormulas){
                //selection cuts
                bool passesSelection = cutf.second->EvalInstance();
                if(!passesSelection) continue;

                //check if this sample is used for this control region
                bool found = false;
                for(auto &name : controlRegionMC[cutf.first]){
                    if(tree.first == name) found = true;
                }
                if(!found) continue;

                ///////////////////////////////////////////////////////////
                // Apply scale factors
                ///////////////////////////////////////////////////////////

                //btagging scale factor
                if(cutf.first == "TTBarSingleLepton" || cutf.first == "TTBarDilepton" || cutf.first == "WSingleLepton"){
                    eventWeight *= getBTagMediumScaleFactor(bjet1Pt, bjet1PassMedium, bjet2Pt, bjet2PassMedium);
                }

                //trigger scale factors
                if(hlt_singlemu && (cutf.first == "TTBarSingleLepton" || cutf.first == "WSingleLepton" || cutf.first == "ZNuNuFromW")){
                    eventWeight *= singleMuTriggerSF;
                }
                else if(hlt_dimuon && (cutf.first == "TTBarDilepton" || cutf.first == "ZLLDilepton" || cutf.first == "ZNuNuFromDY")){
                    eventWeight *= doubleMuTriggerSF;
                    eventWeight *= doubleMuNormalizationSF;
                }

                //TODO: muon and electron ID scale factors

                //Data/MC scale factors
                //TTJets SF
                if(tree.first == "TTJets"){
                    double SFTTJets = SFHistTTJets->GetBinContent(SFHistTTJets->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                    double SFErrorTTJets = SFHistTTJets->GetBinError(SFHistTTJets->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                    if(SFTTJets < 1e5){
                        eventWeight *= SFTTJets;
                        sysErrorSquared += weight*weight*SFErrorTTJets*SFErrorTTJets;
                    }
                    else{
                        //cout << "Warning: TTJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //WJets SF
                else if(tree.first == "WJets"){
                    double SFWJets = SFHistWJets->GetBinContent(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                    double SFErrorWJets = SFHistWJets->GetBinError(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                    if(SFWJets < 1e5){
                        eventWeight *= SFWJets;
                        sysErrorSquared += weight*weight*SFErrorWJets*SFErrorWJets;
                    }
                    else{
                        //cout << "Warning: WJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //DYJets SF
                else if(tree.first == "DYJets"){
                    double SFDYJets = SFHistDYJets->GetBinContent(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                    double SFErrorDYJets = SFHistDYJets->GetBinError(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                    if(SFDYJets < 1e5){
                        eventWeight *= SFDYJets;
                        sysErrorSquared += weight*weight*SFErrorDYJets*SFErrorDYJets;
                    }
                    else{
                        //cout << "Warning: DYJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //ZNuNu SF
                else if(tree.first == "ZJetsNuNu"){
                    double SFZJetsNuNu = SFHistZJetsNuNu->GetBinContent(SFHistZJetsNuNu->FindFixBin(min(MR, SFmaxMRZJetsNuNu), min(Rsq, SFmaxRsqZJetsNuNu)));
                    double SFErrorZJetsNuNu = SFHistZJetsNuNu->GetBinError(SFHistZJetsNuNu->FindFixBin(min(MR, SFmaxMRZJetsNuNu), min(Rsq, SFmaxRsqZJetsNuNu)));
                    if(SFZJetsNuNu < 1e5){
                        eventWeight *= SFZJetsNuNu;
                        sysErrorSquared += weight*weight*SFErrorZJetsNuNu*SFErrorZJetsNuNu;
                    }
                    else{
                        //cout << "Warning: ZJetsNuNu scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }

                ///////////////////////////////////////////////////////////
                // Fill histograms
                ///////////////////////////////////////////////////////////

                //ZNuNuFromDY CR
                if(cutf.first == "ZNuNuFromDY"){
                    razorHistosMC[cutf.first][tree.first]->Fill(MR_noZ, Rsq_noZ, eventWeight);
                    razorErrorHistosMC[cutf.first][tree.first]->Fill(MR_noZ, Rsq_noZ, sysErrorSquared);
                }
                //ZNuNuFromW CR
                else if(cutf.first == "ZNuNuFromW"){
                    razorHistosMC[cutf.first][tree.first]->Fill(MR_noW, Rsq_noW, eventWeight);
                    razorErrorHistosMC[cutf.first][tree.first]->Fill(MR_noW, Rsq_noW, sysErrorSquared);
                }
                //ZNuNuFromGamma CR
                else if(cutf.first == "ZNuNuFromGamma"){
                    razorHistosMC[cutf.first][tree.first]->Fill(MR_noPho, Rsq_noPho, eventWeight);
                    razorErrorHistosMC[cutf.first][tree.first]->Fill(MR_noPho, Rsq_noPho, sysErrorSquared);
                }
                //other CRs
                else{
                    razorHistosMC[cutf.first][tree.first]->Fill(MR, Rsq, eventWeight);
                    razorErrorHistosMC[cutf.first][tree.first]->Fill(MR, Rsq, sysErrorSquared);
                }
            }
        }
    }

    //update errors to take into account systematic uncertainties
    for(auto &tree : mctrees){
        for(auto &cut : cuts){
            for(int i = 0; i < razorHistosMC[cut.first][tree.first]->GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosMC[cut.first][tree.first]->GetNbinsY()+1; j++){
                    double squaredError = razorErrorHistosMC[cut.first][tree.first]->GetBinContent(i, j);
                    razorHistosMC[cut.first][tree.first]->SetBinError(i, j, sqrt(pow(razorHistosMC[cut.first][tree.first]->GetBinError(i, j), 2) + squaredError));
                }
            }
        }
    }

    //for rare background processes, include a 20% uncertainty on the total yield in each bin, summed in quadrature with the statistical uncertainty
    double sysErrorFrac = 0.2;
    //for QCD, assign a 100% uncertainty
    double qcdErrorFrac = 1.0;
    for(auto &tree : mctrees){
        //only do this for rare processes 
        if(tree.first == "DYJets" || tree.first == "WJets" || tree.first == "ZJetsNuNu" || tree.first == "TTJets") continue; 
        for(auto &cut : cuts){
            for(int i = 0; i < razorHistosMC[cut.first][tree.first]->GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosMC[cut.first][tree.first]->GetNbinsY()+1; j++){
                    double error = 0.0;
                    if(tree.first == "QCD"){
                        error = qcdErrorFrac*razorHistosMC[cut.first][tree.first]->GetBinContent(i, j);
                    }
                    else{
                        error = sysErrorFrac*razorHistosMC[cut.first][tree.first]->GetBinContent(i, j);
                    }
                    razorHistosMC[cut.first][tree.first]->SetBinError(i, j, sqrt(pow(razorHistosMC[cut.first][tree.first]->GetBinError(i, j), 2) + error*error));
                }
            }
        }
    }


    ///////////////////////////////////////////////////////////
    // Get data distributions
    ///////////////////////////////////////////////////////////
    map<string, TH2F*> razorHistosData; //apply only efficiency and acceptance corrections
    for(auto &cut : cuts){
        razorHistosData[cut.first] = new TH2F(Form("razordata%s", cut.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        razorHistosData[cut.first]->Sumw2();
    }
    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;

        //make TTreeFormulas for selection cuts
        map<string, TTreeFormula*> cutsFormulas;
        for(auto &cut : cuts){
            cutsFormulas[cut.first] = new TTreeFormula(Form("%s%sCutsFormula", cut.first.c_str(), tree.first.c_str()), cuts[cut.first].c_str(), tree.second);
            cutsFormulas[cut.first]->GetNdata();
        }

        //loop over entries
        uint nEntries = tree.second->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i);

            //noise filters
            if(!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) continue;

            for(auto &cutf : cutsFormulas){
                //apply selection cuts
                bool passesSelection = cutf.second->EvalInstance();
                if(!passesSelection) continue;

                //check if this sample is used for this control region
                bool found = false;
                for(auto &name : controlRegionData[cutf.first]){
                    if(tree.first == name) found = true;
                }
                if(!found) continue;

                float eventWeight = 1.0;

                //reweight for photon triggers
                if(cutf.first == "ZNuNuFromGamma"){
                    double triggerWeightRestricted = 0.0;
                    //get weight if associate each photon trigger with a particular pt range
                    if(passedHLTPhoton150 && leadingPhotonPt > 165){ 
                        triggerWeightRestricted = 1.0;
                    }
                    else if(passedHLTPhoton135 && leadingPhotonPt > 150 && leadingPhotonPt < 165){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton135;
                    }
                    else if(passedHLTPhoton90 && leadingPhotonPt > 100 && leadingPhotonPt < 150){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton90;
                    }
                    else if(passedHLTPhoton75 && leadingPhotonPt > 90 && leadingPhotonPt < 100){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton75; 
                    }
                    else if(passedHLTPhoton50 && leadingPhotonPt < 90){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton50;
                    }
                    eventWeight *= triggerWeightRestricted;

                    if(leadingPhotonPt>5000) continue; // reject noise
                }

                ///////////////////////////////////////////////////////////
                // Fill histograms
                ///////////////////////////////////////////////////////////

                //ZNuNuFromDY CR
                if(cutf.first == "ZNuNuFromDY"){
                    razorHistosData[cutf.first]->Fill(MR_noZ, Rsq_noZ, eventWeight);
                }
                //ZNuNuFromW CR
                else if(cutf.first == "ZNuNuFromW"){
                    razorHistosData[cutf.first]->Fill(MR_noW, Rsq_noW, eventWeight);
                }
                //ZNuNuFromGamma CR
                else if(cutf.first == "ZNuNuFromGamma"){
                    razorHistosData[cutf.first]->Fill(MR_noPho, Rsq_noPho, eventWeight);
                }
                //other CRs
                else{
                    razorHistosData[cutf.first]->Fill(MR, Rsq, eventWeight);
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////
    // Make plots
    ///////////////////////////////////////////////////////////
    TCanvas c("c", "c", 800, 600);
    for(auto &cut : cuts){
        gStyle->SetPaintTextFormat("1.0f");
        c.SetLogx();

        razorHistosMC[cut.first]["QCD"]->SetFillColor(33);
        razorHistosMC[cut.first]["ZJetsNuNu"]->SetFillColor(kCyan+1);
        razorHistosMC[cut.first]["WJets"]->SetFillColor(kRed+1);
        razorHistosMC[cut.first]["TTJets"]->SetFillColor(kGreen+3);
        razorHistosMC[cut.first]["DYJets"]->SetFillColor(kAzure);
        razorHistosMC[cut.first]["SingleTop"]->SetFillColor(kBlue+3);
        razorHistosMC[cut.first]["TTV"]->SetFillColor(kSpring);
        razorHistosMC[cut.first]["VV"]->SetFillColor(kViolet+2);
        razorHistosMC[cut.first]["GJets"]->SetFillColor(8);
        razorHistosMC[cut.first]["VG"]->SetFillColor(38);
        razorHistosMC[cut.first]["TTG"]->SetFillColor(7);

        //create legend
        TLegend *RazorLegend = new TLegend(0.6, 0.6, 0.9, 0.9);
        for(auto &name : controlRegionMC[cut.first]){
            RazorLegend->AddEntry(razorHistosMC[cut.first][name], name.c_str());
        }
        RazorLegend->AddEntry(razorHistosData[cut.first], "2012 Data");

        delete RazorLegend;

        //plot slices of MR and Rsq
        for(int i = 0; i < nRsqBins; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorHistosData[cut.first]->ProjectionX(Form("ThisRsqSliceData%d%s", i, cut.first.c_str()), i+1, i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqSliceMC", Form("MR (%.2f < Rsq < %.2f), %s Box", RsqBinLowEdges[i], RsqBinLowEdges[i+1], cut.first.c_str()));
            for(auto &hist : razorHistosMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second->ProjectionX(Form("hist%s%d%s", hist.first.c_str(), i, cut.first.c_str()), i+1, i+1);
                thisHist->SetFillColor(razorHistosMC[cut.first][hist.first]->GetFillColor());
                ThisRsqSliceMCMap[hist.first] = thisHist;
            }
            for(auto &name : controlRegionMC[cut.first]){
                ThisRsqSliceMC->Add(ThisRsqSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", "MRExclusiveSlice"+to_string(i)+cut.first, true);
        }
        for(int i = 0; i < nMRBins; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorHistosData[cut.first]->ProjectionY(Form("ThisMRSliceData%d%s", i, cut.first.c_str()), i+1, i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (%0.f < MR < %.0f), %s Box", MRBinLowEdges[i], MRBinLowEdges[i+1], cut.first.c_str()));
            for(auto &hist : razorHistosMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second->ProjectionY(Form("hist%s%d%s", hist.first.c_str(), i, cut.first.c_str()), i+1, i+1);
                thisHist->SetFillColor(razorHistosMC[cut.first][hist.first]->GetFillColor());
                ThisMRSliceMCMap[hist.first] = thisHist;
            }
            for(auto &name : controlRegionMC[cut.first]){
                ThisMRSliceMC->Add(ThisMRSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", "RsqExclusiveSlice"+to_string(i)+cut.first, true);
        }
        //inclusive slices
        for(int i = 0; i < nRsqBins; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorHistosData[cut.first]->ProjectionX(Form("ThisRsqIncSliceData%d%s", i, cut.first.c_str()), i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqIncSliceMC", Form("MR (Rsq > %.2f), %s Box", RsqBinLowEdges[i], cut.first.c_str()));
            for(auto &hist : razorHistosMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second->ProjectionX(Form("histinc%s%d%s", hist.first.c_str(), i, cut.first.c_str()), i+1);
                thisHist->SetFillColor(razorHistosMC[cut.first][hist.first]->GetFillColor());
                ThisRsqSliceMCMap[hist.first] = thisHist;
            }
            for(auto &name : controlRegionMC[cut.first]){
                ThisRsqSliceMC->Add(ThisRsqSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", "MRInclusiveSlice"+to_string(i)+cut.first, true);
        }
        for(int i = 0; i < nMRBins; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorHistosData[cut.first]->ProjectionY(Form("ThisMRSliceIncData%d%s", i, cut.first.c_str()), i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (MR > %.0f), %s Box", MRBinLowEdges[i], cut.first.c_str()));
            for(auto &hist : razorHistosMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second->ProjectionY(Form("histinc%s%d%s", hist.first.c_str(), i, cut.first.c_str()), i+1);
                thisHist->SetFillColor(razorHistosMC[cut.first][hist.first]->GetFillColor());
                ThisMRSliceMCMap[hist.first] = thisHist;
            }
            for(auto &name : controlRegionMC[cut.first]){
                ThisMRSliceMC->Add(ThisMRSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", "RsqInclusiveSlice"+to_string(i)+cut.first, true);
        }
    }
}

int main(){
    CompareDataFitCRs();
    return 0;
}

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX){
    TCanvas c("c", "c", 800, 600);
    c.Clear();
    c.cd();
    TPad pad1("pad1","pad1",0,0.4,1,1);
    pad1.SetBottomMargin(0);
    pad1.SetLogy();
    if(logX) pad1.SetLogx();
    pad1.Draw();
    pad1.cd();
    mcStack->Draw("hist");
    mcStack->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    mcStack->GetYaxis()->SetLabelSize(0.03);
    mcStack->GetYaxis()->SetTitleOffset(0.45);
    mcStack->GetYaxis()->SetTitleSize(0.05);
    mcStack->SetMinimum(0.1);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(1);
    dataHist->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    dataHist->Draw("pesame");
    pad1.Modified();
    gPad->Update();
    //make ratio histogram
    TList * histList = (TList*)mcStack->GetHists();
    TIter next(histList);
    TH1 *mcTotal = (TH1*) histList->First()->Clone();
    //mcTotal->Sumw2();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();

    string histoName = dataHist->GetName() ;
    if(histoName.find("mr") != std::string::npos  )
    {
        cout<<"Number of events in data: "<<dataHist->Integral()<<" "<<printString<<endl;
        cout<<"Number of events in MC: "<<mcTotal->Integral()<<" "<<endl;
    }

    //dataOverMC->Sumw2();
    dataOverMC->Divide(mcTotal);
    dataOverMC->GetXaxis()->SetTitle(xaxisTitle.c_str());
    dataOverMC->GetYaxis()->SetTitle("Data / MC");
    dataOverMC->SetMinimum(0.5);
    dataOverMC->SetMaximum(1.5);
    dataOverMC->GetXaxis()->SetLabelSize(0.1);
    dataOverMC->GetYaxis()->SetLabelSize(0.08);
    dataOverMC->GetYaxis()->SetTitleOffset(0.35);
    dataOverMC->GetXaxis()->SetTitleOffset(1.00);
    dataOverMC->GetYaxis()->SetTitleSize(0.08);
    dataOverMC->GetXaxis()->SetTitleSize(0.08);
    dataOverMC->SetStats(0);
    leg->Draw();
    c.cd();
    TPad pad2("pad2","pad2",0,0.0,1,0.4);
    pad2.SetTopMargin(0);
    pad2.SetTopMargin(0.008);
    pad2.SetBottomMargin(0.25);
    pad2.SetGridy();
    if(logX) pad2.SetLogx();
    pad2.Draw();
    pad2.cd();
    dataOverMC->Draw("pe");
    pad2.Modified();
    gPad->Update();
    c.Print(Form("%s.pdf", printString.c_str()));
    // c.Print(Form("%s.root", printString.c_str()));
}
