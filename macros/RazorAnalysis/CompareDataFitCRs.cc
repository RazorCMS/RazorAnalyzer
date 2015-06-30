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

#include "RazorAnalyzer/include/ControlSampleEvents.h"

using namespace std;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void SetHistogramColor(TH1 *hist, string name){
    if(name == "QCD") hist->SetFillColor(33);
    if(name == "ZJetsNuNu") hist->SetFillColor(kCyan+1);
    if(name == "WJets") hist->SetFillColor(kRed+1);
    if(name == "TTJets") hist->SetFillColor(kGreen+3);
    if(name == "DYJets") hist->SetFillColor(kAzure);
    if(name == "SingleTop") hist->SetFillColor(kBlue+3);
    if(name == "TTV") hist->SetFillColor(kSpring);
    if(name == "VV") hist->SetFillColor(kViolet+2);
    if(name == "GJets") hist->SetFillColor(8);
    if(name == "VG") hist->SetFillColor(38);
    if(name == "TTG") hist->SetFillColor(7);
}

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
    //(note: invariant mass cuts are handled within the event loop)
    map<string, string> cuts;
    cuts["TTBarSingleLepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && MET > 30 && lep1MT > 30 && lep1MT < 100 && lep1.Pt() > 25 && NBJetsMedium >= 1 && (HLTDecision[0] || HLTDecision[1] || HLTDecision[8] || HLTDecision[9]) && MR > 300 && Rsq > 0.15";
    cuts["TTBarDilepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1.Pt() > 25 && lep2.Pt() > 25 && lep1PassLoose && lep2PassLoose && (HLTDecision[3] || HLTDecision[4] || HLTDecision[12] || HLTDecision[6] || HLTDecision[7]) && MET > 40 && NBJetsMedium >= 1 && MR > 300 && Rsq > 0.15";
    cuts["WSingleLepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1.Pt() > 25 && MET > 30 && NBJetsMedium == 0 && lep1MT > 30 && lep1MT < 100 && (HLTDecision[0] || HLTDecision[1] || HLTDecision[8] || HLTDecision[9]) && MR > 300 && Rsq > 0.15";
    cuts["ZLLDilepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && NBJetsMedium == 0 && (lep1.Pt() > 25 && lep2.Pt() > 25 && lep1PassLoose && lep2PassLoose) && (HLTDecision[3] || HLTDecision[4] || HLTDecision[12]) && MR > 300 && Rsq > 0.15";
    //cuts["ZNuNuFromDY"] = "MR_noZ > 300";
    //cuts["ZNuNuFromW"] = "MR_noW > 300";
    //cuts["ZNuNuFromGamma"] = "MR_noPho > 300";
    //cuts["ZNuNuFromDY"] = "HLT_Dimuon && nLooseMuons == 2 && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    //cuts["ZNuNuFromW"] = "HLT_SingleMu && nBTaggedJets == 0 && nTightMuons == 1 && nLooseMuons == 1 && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100";
    //cuts["ZNuNuFromGamma"] = "HLT_Photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1 && pho1.Pt() > 80";

    //get input files
    map<string, map<string, string> > mcfilenames;
    map<string, map<string, string> > datafilenames;
    string baseDir = "root://eoscms://store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/";
    string mcSuffix = "_1pb_weighted.root";
    string dataSuffix = "_GoodLumi.root";

    mcfilenames["SingleLeptonRazorSkim"] = map<string, string>();
    mcfilenames["DilptonRazorSkim"] = map<string, string>();
    datafilenames["SingleLeptonRazorSkim"] = map<string, string>();
    datafilenames["DileptonRazorSkim"] = map<string, string>();

    mcfilenames["SingleLeptonRazorSkim"]["DYJets"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_DYJetsToLL_HTBinned"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["TTJets"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_TTJets"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["WJets"] = baseDir+"SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_WJetsToLNu_HTBinned"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["TTV"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_TTV"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["VV"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_VV"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["QCD"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_QCD"+mcSuffix;
    mcfilenames["SingleLeptonRazorSkim"]["SingleTop"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_SingleTop"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["DYJets"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DYJetsToLL_HTBinned"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["TTJets"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTJets"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["WJets"] = baseDir+"DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_WJetsToLNu_HTBinned"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["TTV"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTV"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["VV"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_VV"+mcSuffix;
    mcfilenames["DileptonRazorSkim"]["SingleTop"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_SingleTop"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["GJets"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_GJets_HTBinned"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["VG"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_VG"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["TTG"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_TTG"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["QCD"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_QCD"+mcSuffix;

    //data
    datafilenames["SingleLeptonRazorSkim"]["SingleMuon"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_Data_SingleMu"+dataSuffix;
    datafilenames["SingleLeptonRazorSkim"]["SingleElectron"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_Data_SingleElectron"+dataSuffix;
    datafilenames["DileptonRazorSkim"]["DoubleMuon"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleMuParked"+dataSuffix;
    datafilenames["DileptonRazorSkim"]["DoubleElectron"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleElectron"+dataSuffix;
    //datafilenames["PhotonRazorSkim"]["Photon"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_Data_Photon"+dataSuffix;

    //assign datasets to control regions
    map<string, vector<string> > controlRegionMC;
    controlRegionMC["TTBarSingleLepton"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "WJets", "TTJets"};
    controlRegionMC["TTBarDilepton"] = vector<string> {"TTV", "VV", "SingleTop", "WJets", "DYJets", "TTJets"};
    controlRegionMC["WSingleLepton"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "TTJets", "WJets"};
    controlRegionMC["ZLLDilepton"] = vector<string> {"TTV", "VV", "SingleTop", "WJets", "TTJets", "DYJets"};
    //controlRegionMC["ZNuNuFromDY"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "WJets", "TTJets", "DYJets"};
    //controlRegionMC["ZNuNuFromW"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "TTJets", "WJets"};
    //controlRegionMC["ZNuNuFromGamma"] = vector<string> {"TTG", "VG", "QCD", "GJets"};

    map<string, vector<string> > controlRegionData;
    controlRegionData["TTBarSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["TTBarDilepton"] = vector<string> {"DoubleMuon", "DoubleElectron"};
    controlRegionData["WSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["ZLLDilepton"] = vector<string> {"DoubleMuon", "DoubleElectron"};
    //controlRegionData["ZNuNuFromDY"] = vector<string> {"DoubleMuon"};
    //controlRegionData["ZNuNuFromW"] = vector<string> {"SingleMuon"};
    //controlRegionData["ZNuNuFromGamma"] = vector<string> {"Photon"};

    map<string, string> controlRegionSkim;
    controlRegionSkim["TTBarSingleLepton"] = "SingleLeptonRazorSkim";
    controlRegionSkim["TTBarDilepton"] = "DileptonRazorSkim";
    controlRegionSkim["WSingleLepton"] = "SingleLeptonRazorSkim";
    controlRegionSkim["ZLLDilepton"] = "DileptonRazorSkim";
    //controlRegionSkim["ZNuNuFromDY"] = "DileptonRazorSkim";
    //controlRegionSkim["ZNuNuFromW"] = "SingleLeptonRazorSkim";
    //controlRegionSkim["ZNuNuFromGamma"] = "PhotonRazorSkim";

    //get trees
    map<string, map<string, ControlSampleEvents*> > mcevents;
    map<string, map<string, ControlSampleEvents*> > dataevents;
    mcevents["SingleLeptonRazorSkim"] = map<string, ControlSampleEvents*>();
    mcevents["DileptonRazorSkim"] = map<string, ControlSampleEvents*>();
    //mcevents["PhotonRazorSkim"] = map<string, ControlSampleEvents*>();
    for(auto &skim : mcfilenames){
        for(auto &file : mcfilenames[skim.first]){
            mcevents[skim.first][file.first] = new ControlSampleEvents;
            mcevents[skim.first][file.first]->LoadTree(file.second.c_str());
        }
        for(auto &file : datafilenames[skim.first]){
            dataevents[skim.first][file.first] = new ControlSampleEvents;
            dataevents[skim.first][file.first]->LoadTree(file.second.c_str());
        }
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
    TTreeFormula* cutsFormula;
    ControlSampleEvents *curTree;
    for(auto &region : controlRegionMC){ //region.first is the CR name, region.second is the list of MC samples associated with it
        cout << "Filling MC histograms for control region " << region.first << endl;
        for(auto &sample : region.second){
            cout << "   Filling MC histograms: " << sample << endl;
            string theSkim = controlRegionSkim[region.first];
            curTree = mcevents[theSkim][sample];

            //make histograms, and make TTreeFormulas for selection cuts
            razorHistosMC[region.first][sample] = new TH2F(Form("razormc%s%s", region.first.c_str(), sample.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorErrorHistosMC[region.first][sample] = new TH2F(Form("razorErrormc%s%s", sample.c_str(), region.first.c_str()), "sum(w^2*error(SF)^2); MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosMC[region.first][sample]->Sumw2();
            razorErrorHistosMC[region.first][sample]->Sumw2();
            cutsFormula = new TTreeFormula(Form("%s%sCutsFormula", region.first.c_str(), sample.c_str()), cuts[region.first].c_str(), curTree->tree_);
            cutsFormula->GetNdata();

            //loop over entries
            uint nEntries = curTree->tree_->GetEntries();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                curTree->tree_->GetEntry(i); 

                float eventWeight = curTree->weight*lumiInData*1.0/lumiInMC;

                //PU reweighting
                eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(curTree->NPU_0));
                float sysErrorSquared = 0.0;

                //selection cuts
                bool passesSelection = cutsFormula->EvalInstance();
                if(!passesSelection) continue;

                //extra selection on dilepton mass
                if(region.first == "TTBarDilepton"){
                    float mLL = (curTree->lep1 + curTree->lep2).M();
                    if(mLL > 76 && mLL < 106) continue;
                }
                if(region.first == "ZLLDilepton"){
                    float mLL = (curTree->lep1 + curTree->lep2).M();
                    if(mLL < 60 || mLL > 120) continue;
                }

                ///////////////////////////////////////////////////////////
                // Apply scale factors
                ///////////////////////////////////////////////////////////

                //btagging scale factor
                if(region.first == "TTBarSingleLepton" || region.first == "TTBarDilepton" || region.first == "WSingleLepton"){
                    eventWeight *= getBTagMediumScaleFactor(curTree->bjet1.Pt(), curTree->bjet1PassMedium, curTree->bjet2.Pt(), curTree->bjet2PassMedium);
                }

                //trigger scale factors: single muon
                if((curTree->HLTDecision[0] || curTree->HLTDecision[1]) && (region.first == "TTBarSingleLepton" || region.first == "WSingleLepton" || region.first == "ZNuNuFromW")){
                    eventWeight *= singleMuTriggerSF;
                }
                //trigger scale factors: double muon
                else if((curTree->HLTDecision[3] || curTree->HLTDecision[4]) && (region.first == "TTBarDilepton" || region.first == "ZLLDilepton" || region.first == "ZNuNuFromDY" || region.first == "ZLLDilepton")){
                    eventWeight *= doubleMuTriggerSF;
                    eventWeight *= doubleMuNormalizationSF;
                }

                //TODO: muon and electron ID scale factors

                //Data/MC scale factors
                if(region.first != "ZNuNuFromDY" && region.first != "ZNuNuFromW" && region.first != "ZNuNuFromGamma"){
                    //TTJets SF
                    if(sample == "TTJets"){
                        double SFTTJets = SFHistTTJets->GetBinContent(SFHistTTJets->FindFixBin(min(curTree->MR, SFmaxMRTTJets), min(curTree->Rsq, SFmaxRsqTTJets)));
                        double SFErrorTTJets = SFHistTTJets->GetBinError(SFHistTTJets->FindFixBin(min(curTree->MR, SFmaxMRTTJets), min(curTree->Rsq, SFmaxRsqTTJets)));
                        if(SFTTJets < 1e5){
                            eventWeight *= SFTTJets;
                            sysErrorSquared += curTree->weight*curTree->weight*SFErrorTTJets*SFErrorTTJets;
                        }
                        else{
                            //cout << "Warning: TTJets scale factor is Inf!" << endl;
                            eventWeight = 0;
                            sysErrorSquared = 0;
                        }
                    }
                    //WJets SF
                    else if(sample == "WJets"){
                        double SFWJets = SFHistWJets->GetBinContent(SFHistWJets->FindFixBin(min(curTree->MR, SFmaxMRWJets), min(curTree->Rsq, SFmaxRsqWJets)));
                        double SFErrorWJets = SFHistWJets->GetBinError(SFHistWJets->FindFixBin(min(curTree->MR, SFmaxMRWJets), min(curTree->Rsq, SFmaxRsqWJets)));
                        if(SFWJets < 1e5){
                            eventWeight *= SFWJets;
                            sysErrorSquared += curTree->weight*curTree->weight*SFErrorWJets*SFErrorWJets;
                        }
                        else{
                            //cout << "Warning: WJets scale factor is Inf!" << endl;
                            eventWeight = 0;
                            sysErrorSquared = 0;
                        }
                    }
                    //DYJets SF
                    else if(sample == "DYJets"){
                        double SFDYJets = SFHistDYJets->GetBinContent(SFHistDYJets->FindFixBin(min(curTree->MR, SFmaxMRDYJets), min(curTree->Rsq, SFmaxRsqDYJets)));
                        double SFErrorDYJets = SFHistDYJets->GetBinError(SFHistDYJets->FindFixBin(min(curTree->MR, SFmaxMRDYJets), min(curTree->Rsq, SFmaxRsqDYJets)));
                        if(SFDYJets < 1e5){
                            eventWeight *= SFDYJets;
                            sysErrorSquared += curTree->weight*curTree->weight*SFErrorDYJets*SFErrorDYJets;
                        }
                        else{
                            //cout << "Warning: DYJets scale factor is Inf!" << endl;
                            eventWeight = 0;
                            sysErrorSquared = 0;
                        }
                    }
                    //ZNuNu SF
                    else if(sample == "ZJetsNuNu"){
                        double SFZJetsNuNu = SFHistZJetsNuNu->GetBinContent(SFHistZJetsNuNu->FindFixBin(min(curTree->MR, SFmaxMRZJetsNuNu), min(curTree->Rsq, SFmaxRsqZJetsNuNu)));
                        double SFErrorZJetsNuNu = SFHistZJetsNuNu->GetBinError(SFHistZJetsNuNu->FindFixBin(min(curTree->MR, SFmaxMRZJetsNuNu), min(curTree->Rsq, SFmaxRsqZJetsNuNu)));
                        if(SFZJetsNuNu < 1e5){
                            eventWeight *= SFZJetsNuNu;
                            sysErrorSquared += curTree->weight*curTree->weight*SFErrorZJetsNuNu*SFErrorZJetsNuNu;
                        }
                        else{
                            //cout << "Warning: ZJetsNuNu scale factor is Inf!" << endl;
                            eventWeight = 0;
                            sysErrorSquared = 0;
                        }
                    }
                }

                ///////////////////////////////////////////////////////////
                // Fill histograms
                ///////////////////////////////////////////////////////////

                //ZNuNuFromDY CR
                if(region.first == "ZNuNuFromDY"){
                    //razorHistosMC[region.first][sample]->Fill(curTree->MR_noZ, curTree->Rsq_noZ, eventWeight);
                    //razorErrorHistosMC[region.first][sample]->Fill(curTree->MR_noZ, curTree->Rsq_noZ, sysErrorSquared);
                }
                //ZNuNuFromW CR
                else if(region.first == "ZNuNuFromW"){
                    //razorHistosMC[region.first][sample]->Fill(curTree->MR_noW, curTree->Rsq_noW, eventWeight);
                    //razorErrorHistosMC[region.first][sample]->Fill(curTree->MR_noW, curTree->Rsq_noW, sysErrorSquared);
                }
                //ZNuNuFromGamma CR
                else if(region.first == "ZNuNuFromGamma"){
                    //razorHistosMC[region.first][sample]->Fill(curTree->MR_noPho, curTree->Rsq_noPho, eventWeight);
                    //razorErrorHistosMC[region.first][sample]->Fill(curTree->MR_noPho, curTree->Rsq_noPho, sysErrorSquared);
                }
                //other CRs
                else{
                    razorHistosMC[region.first][sample]->Fill(curTree->MR, curTree->Rsq, eventWeight);
                    razorErrorHistosMC[region.first][sample]->Fill(curTree->MR, curTree->Rsq, sysErrorSquared);
                }
            } //end of event loop

            //update errors to take into account systematic uncertainties
            for(int i = 0; i < razorHistosMC[region.first][sample]->GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosMC[region.first][sample]->GetNbinsY()+1; j++){
                    double squaredError = razorErrorHistosMC[region.first][sample]->GetBinContent(i, j);
                    razorHistosMC[region.first][sample]->SetBinError(i, j, sqrt(pow(razorHistosMC[region.first][sample]->GetBinError(i, j), 2) + squaredError));
                }
            }

            //for rare background processes, include a 20% uncertainty on the total yield in each bin, summed in quadrature with the statistical uncertainty
            double sysErrorFrac = 0.2;
            //for QCD, assign a 100% uncertainty
            double qcdErrorFrac = 1.0;
            //only do this for rare processes 
            if(sample != "DYJets" && sample != "WJets" && sample != "ZJetsNuNu" && sample != "TTJets"){
                for(int i = 0; i < razorHistosMC[region.first][sample]->GetNbinsX()+1; i++){
                    for(int j = 0; j < razorHistosMC[region.first][sample]->GetNbinsY()+1; j++){
                        double error = 0.0;
                        if(sample == "QCD"){
                            error = qcdErrorFrac*razorHistosMC[region.first][sample]->GetBinContent(i, j);
                        }
                        else{
                            error = sysErrorFrac*razorHistosMC[region.first][sample]->GetBinContent(i, j);
                        }
                        razorHistosMC[region.first][sample]->SetBinError(i, j, sqrt(pow(razorHistosMC[region.first][sample]->GetBinError(i, j), 2) + error*error));
                    }
                }
            } //end if
        } //end of loop over samples
    } //end of loop over control regions

    ///////////////////////////////////////////////////////////
    // Get data distributions
    ///////////////////////////////////////////////////////////
    map<string, TH2F*> razorHistosData; //apply only efficiency and acceptance corrections
    for(auto &region : controlRegionData){
        cout << "Filling data histograms for control region " << region.first << endl;
        for(auto &sample : region.second){
            cout << "   Filling data histograms: " << sample << endl;
            string theSkim = controlRegionSkim[region.first];
            curTree = dataevents[theSkim][sample];

            razorHistosData[region.first] = new TH2F(Form("razordata%s%s", region.first.c_str(), sample.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosData[region.first]->Sumw2();

            //make TTreeFormulas for selection cuts
            TTreeFormula* cutsFormula;
            cutsFormula = new TTreeFormula(Form("%s%sCutsFormula", region.first.c_str(), sample.c_str()), cuts[region.first].c_str(), curTree->tree_);
            cutsFormula->GetNdata();

            //loop over entries
            uint nEntries = curTree->tree_->GetEntries();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                curTree->tree_->GetEntry(i);

                //noise filters
                if(!(curTree->Flag_HBHENoiseFilter) || !(curTree->Flag_CSCTightHaloFilter) || !(curTree->Flag_eeBadScFilter) ) continue;

                //apply selection cuts
                bool passesSelection = cutsFormula->EvalInstance();
                if(!passesSelection) continue;

                //extra selection on dilepton mass
                if(region.first == "TTBarDilepton"){
                    float mLL = (curTree->lep1 + curTree->lep2).M();
                    if(mLL > 76 && mLL < 106) continue;
                }
                if(region.first == "ZLLDilepton"){
                    float mLL = (curTree->lep1 + curTree->lep2).M();
                    if(mLL < 60 || mLL > 120) continue;
                }

                float eventWeight = 1.0;

                //reweight for photon triggers
                if(region.first == "ZNuNuFromGamma"){
                    double triggerWeightRestricted = 0.0;
                    //get weight if associate each photon trigger with a particular pt range
                    if(curTree->HLT_Photon150 && curTree->pho1.Pt() > 165){ 
                        triggerWeightRestricted = 1.0;
                    }
                    else if(curTree->HLT_Photon135 && curTree->pho1.Pt() > 150 && curTree->pho1.Pt() < 165){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton135;
                    }
                    else if(curTree->HLT_Photon90 && curTree->pho1.Pt() > 100 && curTree->pho1.Pt() < 150){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton90;
                    }
                    else if(curTree->HLT_Photon75 && curTree->pho1.Pt() > 90 && curTree->pho1.Pt() < 100){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton75; 
                    }
                    else if(curTree->HLT_Photon50 && curTree->pho1.Pt() < 90){
                        triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton50;
                    }
                    eventWeight *= triggerWeightRestricted;

                    if(curTree->pho1.Pt()>5000) continue; // reject noise
                }

                ///////////////////////////////////////////////////////////
                // Fill histograms
                ///////////////////////////////////////////////////////////

                //ZNuNuFromDY CR
                if(region.first == "ZNuNuFromDY"){
                    //razorHistosData[region.first]->Fill(curTree->MR_noZ, curTree->Rsq_noZ, eventWeight);
                }
                //ZNuNuFromW CR
                else if(region.first == "ZNuNuFromW"){
                    //razorHistosData[region.first]->Fill(curTree->MR_noW, curTree->Rsq_noW, eventWeight);
                }
                //ZNuNuFromGamma CR
                else if(region.first == "ZNuNuFromGamma"){
                    //razorHistosData[region.first]->Fill(curTree->MR_noPho, curTree->Rsq_noPho, eventWeight);
                }
                //other CRs
                else{
                    razorHistosData[region.first]->Fill(curTree->MR, curTree->Rsq, eventWeight);
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////
    // Make plots
    ///////////////////////////////////////////////////////////
    cout << "Saving output plots..." << endl;
    TCanvas c("c", "c", 800, 600);
    for(auto &cut : cuts){
        cout << "Control region: " << cut.first << endl;
        gStyle->SetPaintTextFormat("1.0f");
        c.SetLogx();

        //create legend
        TLegend *RazorLegend = new TLegend(0.6, 0.6, 0.9, 0.9);
        for(int i = controlRegionMC[cut.first].size() - 1; i >= 0; i--){
            string name = controlRegionMC[cut.first][i];
            SetHistogramColor(razorHistosMC[cut.first][name], name);
            RazorLegend->AddEntry(razorHistosMC[cut.first][name], name.c_str());
        }
        razorHistosData[cut.first]->SetMarkerStyle(20);
        razorHistosData[cut.first]->SetMarkerSize(1);
        RazorLegend->AddEntry(razorHistosData[cut.first], "2012 Data");

        //plot slices of MR and Rsq
        for(int i = 0; i < nRsqBins; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorHistosData[cut.first]->ProjectionX(Form("ThisRsqSliceData%d%s", i, cut.first.c_str()), i+1, i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqSliceMC", Form("MR (%.2f < Rsq < %.2f), %s Box", RsqBinLowEdges[i], RsqBinLowEdges[i+1], cut.first.c_str()));
            for(auto &sample : controlRegionMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cut.first][sample]->ProjectionX(Form("hist%s%d%s", sample.c_str(), i, cut.first.c_str()), i+1, i+1);
                SetHistogramColor(thisHist, sample);
                ThisRsqSliceMCMap[sample] = thisHist;
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
            for(auto &sample : controlRegionMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cut.first][sample]->ProjectionY(Form("hist%s%d%s", sample.c_str(), i, cut.first.c_str()), i+1, i+1);
                SetHistogramColor(thisHist, sample);
                ThisMRSliceMCMap[sample] = thisHist;
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
            for(auto &sample : controlRegionMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cut.first][sample]->ProjectionX(Form("histinc%s%d%s", sample.c_str(), i, cut.first.c_str()), i+1);
                SetHistogramColor(thisHist, sample);
                ThisRsqSliceMCMap[sample] = thisHist;
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
            for(auto &sample : controlRegionMC[cut.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cut.first][sample]->ProjectionY(Form("histinc%s%d%s", sample.c_str(), i, cut.first.c_str()), i+1);
                SetHistogramColor(thisHist, sample);
                ThisMRSliceMCMap[sample] = thisHist;
            }
            for(auto &name : controlRegionMC[cut.first]){
                ThisMRSliceMC->Add(ThisMRSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", "RsqInclusiveSlice"+to_string(i)+cut.first, true);
        }

        delete RazorLegend;
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
