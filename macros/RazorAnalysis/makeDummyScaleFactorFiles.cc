#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

void makeDummyScaleFactorFiles(){

    //make dummy pileup reweighting histogram
    TFile *pileupFile = new TFile("DummyRun2PileupWeights.root", "recreate");
    TH1F *pileupHist = new TH1F("PUWeight_Run2", "PUWeight_Run2", 1, 0, 100);
    pileupHist->SetBinContent(1, 1.0);
    pileupHist->Write();
    pileupFile->Close();

    //make dummy muon scale factor histogram
    TFile *muonFile = new TFile("DummyRun2MuonWeights.root", "recreate");
    TH2F *muonLooseHist = new TH2F("MuonWeight_Run2_Loose", "MuonWeight_Run2_Loose", 1, -3, 3, 1, 0, 200);
    TH2F *muonTightHist = new TH2F("MuonWeight_Run2_Tight", "MuonWeight_Run2_Tight", 1, -3, 3, 1, 0, 200);
    muonLooseHist->SetBinContent(1, 1, 1.0);
    muonTightHist->SetBinContent(1, 1, 1.0);
    muonLooseHist->SetBinError(1, 1, .01);
    muonTightHist->SetBinError(1, 1, .01);
    muonLooseHist->Write();
    muonTightHist->Write();
    muonFile->Close();

    //make dummy electron scale factor histogram
    TFile *electronFile = new TFile("DummyRun2EleWeights.root", "recreate");
    TH2F *electronLooseHist = new TH2F("EleWeight_Run2_Loose", "EleWeight_Run2_Loose", 1, -3, 3, 1, 0, 200);
    TH2F *electronTightHist = new TH2F("EleWeight_Run2_Tight", "EleWeight_Run2_Tight", 1, -3, 3, 1, 0, 200);
    electronLooseHist->SetBinContent(1, 1, 1.0);
    electronTightHist->SetBinContent(1, 1, 1.0);
    electronLooseHist->SetBinError(1, 1, .01);
    electronTightHist->SetBinError(1, 1, .01);
    electronLooseHist->Write();
    electronTightHist->Write();
    electronFile->Close();
}
