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
#include "include/MacroHelper.h"

using namespace std;

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

void CompareDataFitCRs(){
    //bool makeScaleFactors = true; //compute the scale factors for each control region
    bool makeScaleFactors = false; 

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

    //get input files
    map<string, map<string, string> > mcfilenames;
    map<string, map<string, string> > datafilenames;
    string baseDir = "root://eoscms://store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/";
    string mcSuffix = "_1pb_weighted.root";
    string dataSuffix = "_GoodLumi.root";

    mcfilenames["SingleLeptonRazorSkim"] = map<string, string>();
    mcfilenames["DileptonRazorSkim"] = map<string, string>();
    mcfilenames["ZNuNuDileptonRazorSkim"] = map<string, string>();
    datafilenames["SingleLeptonRazorSkim"] = map<string, string>();
    datafilenames["DileptonRazorSkim"] = map<string, string>();
    datafilenames["ZNuNuDileptonRazorSkim"] = map<string, string>();

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
    mcfilenames["ZNuNuDileptonRazorSkim"]["DYJets"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_DYJetsToLL_HTBinned"+mcSuffix;
    mcfilenames["ZNuNuDileptonRazorSkim"]["TTJets"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_TTJets"+mcSuffix;
    mcfilenames["ZNuNuDileptonRazorSkim"]["WJets"] = baseDir+"ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_WJetsToLNu_HTBinned"+mcSuffix;
    mcfilenames["ZNuNuDileptonRazorSkim"]["TTV"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_TTV"+mcSuffix;
    mcfilenames["ZNuNuDileptonRazorSkim"]["VV"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_VV"+mcSuffix;
    mcfilenames["ZNuNuDileptonRazorSkim"]["SingleTop"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_SingleTop"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["GJets"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_GJets_HTBinned"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["VG"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_VG"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["TTG"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_TTG"+mcSuffix;
    //mcfilenames["PhotonRazorSkim"]["QCD"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_QCD"+mcSuffix;

    //data
    datafilenames["SingleLeptonRazorSkim"]["SingleMuon"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_Data_SingleMu"+dataSuffix;
    datafilenames["SingleLeptonRazorSkim"]["SingleElectron"] = baseDir+"/SingleLeptonRazorSkim/RunOneRazorControlRegions_SingleLeptonRazorSkim_Data_SingleElectron"+dataSuffix;
    datafilenames["DileptonRazorSkim"]["DoubleMuon"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleMuParked"+dataSuffix;
    datafilenames["DileptonRazorSkim"]["DoubleElectron"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleElectron"+dataSuffix;
    datafilenames["DileptonRazorSkim"]["MuE"] = baseDir+"/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_MuEG"+dataSuffix;
    datafilenames["ZNuNuDileptonRazorSkim"]["DoubleMuon"] = baseDir+"/ZNuNuDileptonSkim/RunOneRazorControlRegions_ZNuNuDileptonSkim_Data_DoubleMuParked"+dataSuffix;
    //datafilenames["PhotonRazorSkim"]["Photon"] = baseDir+"/PhotonRazorSkim/RunOneRazorControlRegions_PhotonRazorSkim_Data_Photon"+dataSuffix;

    //assign datasets to control regions
    map<string, vector<string> > controlRegionMC;
    controlRegionMC["TTBarSingleLepton"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "WJets", "TTJets"};
    controlRegionMC["TTBarDilepton"] = vector<string> {"TTV", "VV", "SingleTop", "WJets", "DYJets", "TTJets"};
    controlRegionMC["WSingleLepton"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "TTJets", "WJets"};
    controlRegionMC["ZLLDilepton"] = vector<string> {"TTV", "VV", "SingleTop", "WJets", "TTJets", "DYJets"};
    controlRegionMC["ZNuNuDilepton"] = vector<string> {"TTV", "VV", "SingleTop", "WJets", "TTJets", "DYJets"};
    //controlRegionMC["ZNuNuSingleLepton"] = vector<string> {"TTV", "VV", "QCD", "SingleTop", "DYJets", "TTJets", "WJets"};
    //controlRegionMC["ZNuNuPhoton"] = vector<string> {"TTG", "VG", "QCD", "GJets"};

    map<string, vector<string> > controlRegionData;
    controlRegionData["TTBarSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["TTBarDilepton"] = vector<string> {"DoubleMuon", "DoubleElectron", "MuE"};
    controlRegionData["WSingleLepton"] = vector<string> {"SingleMuon", "SingleElectron"};
    controlRegionData["ZLLDilepton"] = vector<string> {"DoubleMuon", "DoubleElectron"};
    controlRegionData["ZNuNuDilepton"] = vector<string> {"DoubleMuon"};
    //controlRegionData["ZNuNuSingleLepton"] = vector<string> {"SingleMuon"};
    //controlRegionData["ZNuNuPhoton"] = vector<string> {"Photon"};

    map<string, string> controlRegionSkim;
    controlRegionSkim["TTBarSingleLepton"] = "SingleLeptonRazorSkim";
    controlRegionSkim["TTBarDilepton"] = "DileptonRazorSkim";
    controlRegionSkim["WSingleLepton"] = "SingleLeptonRazorSkim";
    controlRegionSkim["ZLLDilepton"] = "DileptonRazorSkim";
    controlRegionSkim["ZNuNuDilepton"] = "ZNuNuDileptonRazorSkim";
    //controlRegionSkim["ZNuNuSingleLepton"] = "SingleLeptonRazorSkim";
    //controlRegionSkim["ZNuNuPhoton"] = "PhotonRazorSkim";

    //get trees
    map<string, map<string, ControlSampleEvents*> > mcevents;
    map<string, map<string, ControlSampleEvents*> > dataevents;
    mcevents["SingleLeptonRazorSkim"] = map<string, ControlSampleEvents*>();
    mcevents["DileptonRazorSkim"] = map<string, ControlSampleEvents*>();
    mcevents["ZNuNuDileptonRazorSkim"] = map<string, ControlSampleEvents*>();
    //mcevents["PhotonRazorSkim"] = map<string, ControlSampleEvents*>();
    map<string, int> kTreeTypes;
    kTreeTypes["SingleLeptonRazorSkim"] = 1;
    kTreeTypes["DileptonRazorSkim"] = 5;
    kTreeTypes["ZNuNuDileptonRazorSkim"] = 7;
    for(auto &skim : mcfilenames){
        for(auto &file : mcfilenames[skim.first]){
            mcevents[skim.first][file.first] = new ControlSampleEvents;
            mcevents[skim.first][file.first]->LoadTree(file.second.c_str(), kTreeTypes[skim.first]);
        }
        for(auto &file : datafilenames[skim.first]){
            dataevents[skim.first][file.first] = new ControlSampleEvents;
            dataevents[skim.first][file.first]->LoadTree(file.second.c_str(), kTreeTypes[skim.first]);
        }
    }

    //luminosities collected by the various photon triggers
    float lumi_HLTPhoton50  = 1.353e0 + 4.921e0 + 7.947e0 + 8.131e0;
    float lumi_HLTPhoton75  = 8.111e0 + 2.953e1 + 4.768e1 + 4.879e1;
    float lumi_HLTPhoton90  = 1.622e1 + 6.408e1 + 1.010e2 + 9.948e1;
    float lumi_HLTPhoton135 = 8.893e2 + 1.476e2 + 5.429e3 + 7.318e3;
    float lumi_HLTPhoton150 = 8.893e2 + 4.429e3 + 7.152e3 + 7.318e3;

    map<string, TH2F*> SFHists;
    //load TTbar scale factor histograms
    TFile *SFFileTTJets = new TFile("data/ScaleFactors/Run1/TTBarSingleLeptonScaleFactors.root");
    SFHists["TTJets"] = (TH2F*)SFFileTTJets->Get("TTBarSingleLeptonScaleFactor");
    //TFile *SFFileTTJets = new TFile("data/ScaleFactors/Run1/TTBarDileptonScaleFactors.root");
    //SFHists["TTJets"] = (TH2F*)SFFileTTJets->Get("TTBarDileptonScaleFactor");

    //load WJets scale factor histogram
    TFile *SFFileWJets = new TFile("data/ScaleFactors/Run1/WJetsSingleLeptonScaleFactors.root");
    SFHists["WJets"] = (TH2F*)SFFileWJets->Get("WJetsSingleLeptonScaleFactor");

    //load DYJets scale factor histogram
    TFile *SFFileDYJets = new TFile("data/ScaleFactors/Run1/ZToLLScaleFactors.root");
    SFHists["DYJets"] = (TH2F*)SFFileDYJets->Get("ZToLLDileptonScaleFactor");

    //load ZNuNu scale factor histograms
    TFile *SFFileZJetsNuNu = new TFile("mcScaleFactorsRunOne.root"); //use the scale factors derived using this macro
    SFHists["ZJetsNuNu"] = (TH2F*)SFFileZJetsNuNu->Get("ZNuNuDileptonScaleFactors");

    //TFile *SFFileZJetsNuNu = new TFile("data/ScaleFactors/Run1/ZInvisibleScaleFactorsRun1.root");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("DYJetsScaleFactors");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("WJetsScaleFactors");
    //SFHists["ZJetsNuNu"] = (TH2F*)SFFileZJetsNuNu->Get("GJetsScaleFactors");

    //load ZNuNu-->DYJets weighting factors
    TFile *ZNuNuToDYWeightFile = new TFile("data/ScaleFactors/Run1/ZNuNuToDYScaleFactorsRun1.root");
    TH2F *ZNuNuToDYWeightHist = (TH2F*)ZNuNuToDYWeightFile->Get("razormcDYJets");
    float maxMRZNuNuToDY = ZNuNuToDYWeightHist->GetXaxis()->GetXmax() - 1;
    float maxRsqZNuNuToDY = ZNuNuToDYWeightHist->GetYaxis()->GetXmax() - 0.01;

    //load muon efficiency scale factor histogram
    TFile muIdSFFile("data/ScaleFactors/MuonEfficiencies_ID_Run2012ReReco_53X.root");
    TFile muIsoSFFile("data/ScaleFactors/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
    //TODO: use muon efficiency scale factors

    //load pileup reweighting histogram
    TFile *pileupWeightFile = new TFile("data/Run1PileupWeights.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
    assert(pileupWeightHist);

    ///////////////////////////////////////////////////////////
    // Get MC distributions
    ///////////////////////////////////////////////////////////
    map<string, map<string, TH2F*> > razorHistosMC;
    map<string, map<string, TH2F*> > razorErrorHistosMC; //store sum(w^2*error(SF)^2)
    ControlSampleEvents *curTree;
    for(auto &region : controlRegionMC){ //region.first is the CR name, region.second is the list of MC samples associated with it
        cout << "Filling MC histograms for control region " << region.first << endl;
        for(auto &sample : region.second){
            cout << "   Filling MC histograms: " << sample << endl;
            string theSkim = controlRegionSkim[region.first];
            curTree = mcevents[theSkim][sample];

            //make histograms
            razorHistosMC[region.first][sample] = new TH2F(Form("razormc%s%s", region.first.c_str(), sample.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorErrorHistosMC[region.first][sample] = new TH2F(Form("razorErrormc%s%s", sample.c_str(), region.first.c_str()), "sum(w^2*error(SF)^2); MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosMC[region.first][sample]->Sumw2();
            razorErrorHistosMC[region.first][sample]->Sumw2();

            //loop over entries
            uint nEntries = curTree->tree_->GetEntries();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                curTree->tree_->GetEntry(i); 

                float eventWeight = curTree->weight*lumiInData*1.0/lumiInMC;
                float sysErrorSquared = 0.0;

                //selection cuts
                if(!curTree->inControlSample(region.first)) continue;

                // Apply scale factors
                eventWeight *= curTree->getMCCorrection(pileupWeightHist, region.first);

                //Data/MC scale factors
                if(!makeScaleFactors && region.first != "ZNuNuDilepton" && region.first != "ZNuNuSingleLepton" && region.first != "ZNuNuPhoton"){
                    if(sample == "TTJets" || sample == "WJets" || sample == "DYJets" || sample == "ZJetsNuNu"){
                        pair<double, double> sfAndErr = getDataMCSFAndError(SFHists[sample], curTree->MR, curTree->Rsq);
                        eventWeight *= sfAndErr.first; //multiply event weight by scale factor
                        sysErrorSquared += curTree->weight*curTree->weight*sfAndErr.second*sfAndErr.second; //add (w*sigma)^2 to the systematic uncertainty
                    }
                }

                ///////////////////////////////////////////////////////////
                // Fill histograms
                ///////////////////////////////////////////////////////////

                //ZNuNuDilepton CR
                if(region.first == "ZNuNuDilepton"){
                    razorHistosMC[region.first][sample]->Fill(curTree->MR_NoZ, curTree->Rsq_NoZ, eventWeight);
                    razorErrorHistosMC[region.first][sample]->Fill(curTree->MR_NoZ, curTree->Rsq_NoZ, sysErrorSquared);
                }
                //ZNuNuSingleLepton CR
                else if(region.first == "ZNuNuSingleLepton"){
                    //razorHistosMC[region.first][sample]->Fill(curTree->MR_noW, curTree->Rsq_noW, eventWeight);
                    //razorErrorHistosMC[region.first][sample]->Fill(curTree->MR_noW, curTree->Rsq_noW, sysErrorSquared);
                }
                //ZNuNuPhoton CR
                else if(region.first == "ZNuNuPhoton"){
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
    map<string, string> signalNames; //desired sample to isolate in each control region
    signalNames["TTBarSingleLepton"] = "TTJets";
    signalNames["TTBarDilepton"] = "TTJets";
    signalNames["WSingleLepton"] = "WJets";
    signalNames["ZLLDilepton"] = "DYJets";
    signalNames["ZNuNuDilepton"] = "DYJets";
    //signalNames["ZNuNuSingleLepton"] = "WJets";
    //signalNames["ZNuNuPhoton"] = "GJets";
    TFile *outSFFile;
    if(makeScaleFactors) outSFFile = new TFile("mcScaleFactorsRunOne.root", "RECREATE");
    for(auto &region : controlRegionData){
        cout << "Filling data histograms for control region " << region.first << endl;
        razorHistosData[region.first] = new TH2F(Form("razordata%s", region.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        razorHistosData[region.first]->Sumw2();

        for(auto &sample : region.second){
            cout << "   Filling data histograms: " << sample << endl;
            string theSkim = controlRegionSkim[region.first];
            curTree = dataevents[theSkim][sample];

            //loop over entries
            uint nEntries = curTree->tree_->GetEntries();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                curTree->tree_->GetEntry(i);

                //noise filters
                if(!(curTree->Flag_HBHENoiseFilter) || !(curTree->Flag_CSCTightHaloFilter) || !(curTree->Flag_eeBadScFilter) ) continue;

                //apply selection cuts
                if(!curTree->inControlSample(region.first)) continue;

                float eventWeight = 1.0;

                //reweight for photon triggers
                if(region.first == "ZNuNuPhoton"){
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

                //ZNuNuDilepton CR
                if(region.first == "ZNuNuDilepton"){
                    razorHistosData[region.first]->Fill(curTree->MR_NoZ, curTree->Rsq_NoZ, eventWeight);
                }
                //ZNuNuSingleLepton CR
                else if(region.first == "ZNuNuSingleLepton"){
                    //razorHistosData[region.first]->Fill(curTree->MR_noW, curTree->Rsq_noW, eventWeight);
                }
                //ZNuNuPhoton CR
                else if(region.first == "ZNuNuPhoton"){
                    //razorHistosData[region.first]->Fill(curTree->MR_noPho, curTree->Rsq_noPho, eventWeight);
                }
                //other CRs
                else{
                    razorHistosData[region.first]->Fill(curTree->MR, curTree->Rsq, eventWeight);
                }
            } //end event loop
        } //end loop over datasets 
        //create background subtracted histogram in each control region
        if(makeScaleFactors){
            outSFFile->cd();
            for(auto &sample : controlRegionMC[region.first]){ //loop over samples
                if(sample != signalNames[region.first]){ //subtract contribution from this process
                    cout << "Removing " << sample << " from " << region.first << " data" << endl;
                    razorHistosData[region.first]->Add(razorHistosMC[region.first][sample], -1);
                } //end if
            } //end loop over samples
            //create data/MC scale factor histogram and save it
            TH2F *dataOverMC = (TH2F*)razorHistosData[region.first]->Clone(Form("%sScaleFactors", region.first.c_str()));
            dataOverMC->Divide(razorHistosMC[region.first][signalNames[region.first]]);
            dataOverMC->Write();
        } //end if
    } //end loop over control regions

    ///////////////////////////////////////////////////////////
    // Make plots
    ///////////////////////////////////////////////////////////
    cout << "Saving output plots..." << endl;
    TCanvas c("c", "c", 800, 600);
    for(auto &cr : controlRegionMC){
        cout << "Control region: " << cr.first << endl;
        gStyle->SetPaintTextFormat("1.0f");
        c.SetLogx();

        //create legend
        TLegend *RazorLegend = new TLegend(0.6, 0.6, 0.9, 0.9);
        cout << "Building legend: ";
        for(int i = cr.second.size() - 1; i >= 0; i--){
            string name = cr.second[i];
            SetHistogramColor(razorHistosMC[cr.first][name], name);
            RazorLegend->AddEntry(razorHistosMC[cr.first][name], name.c_str());
            cout << name.c_str() << " ";
        }
        cout << endl;
        razorHistosData[cr.first]->SetMarkerStyle(20);
        razorHistosData[cr.first]->SetMarkerSize(1);
        RazorLegend->AddEntry(razorHistosData[cr.first], "2012 Data");

        //plot slices of MR and Rsq
        for(int i = 0; i < nRsqBins; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorHistosData[cr.first]->ProjectionX(Form("ThisRsqSliceData%d%s", i, cr.first.c_str()), i+1, i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqSliceMC", Form("MR (%.2f < Rsq < %.2f), %s Box", RsqBinLowEdges[i], RsqBinLowEdges[i+1], cr.first.c_str()));
            for(auto &sample : cr.second){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cr.first][sample]->ProjectionX(Form("hist%s%d%s", sample.c_str(), i, cr.first.c_str()), i+1, i+1);
                SetHistogramColor(thisHist, sample);
                ThisRsqSliceMCMap[sample] = thisHist;
            }
            for(auto &name : cr.second){
                ThisRsqSliceMC->Add(ThisRsqSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", "MRExclusiveSlice"+to_string(i)+cr.first, true);
        }
        for(int i = 0; i < nMRBins; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorHistosData[cr.first]->ProjectionY(Form("ThisMRSliceData%d%s", i, cr.first.c_str()), i+1, i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (%0.f < MR < %.0f), %s Box", MRBinLowEdges[i], MRBinLowEdges[i+1], cr.first.c_str()));
            for(auto &sample : cr.second){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cr.first][sample]->ProjectionY(Form("hist%s%d%s", sample.c_str(), i, cr.first.c_str()), i+1, i+1);
                SetHistogramColor(thisHist, sample);
                ThisMRSliceMCMap[sample] = thisHist;
            }
            for(auto &name : cr.second){
                ThisMRSliceMC->Add(ThisMRSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", "RsqExclusiveSlice"+to_string(i)+cr.first, true);
        }
        //inclusive slices
        for(int i = 0; i < nRsqBins; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorHistosData[cr.first]->ProjectionX(Form("ThisRsqIncSliceData%d%s", i, cr.first.c_str()), i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqIncSliceMC", Form("MR (Rsq > %.2f), %s Box", RsqBinLowEdges[i], cr.first.c_str()));
            for(auto &sample : cr.second){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cr.first][sample]->ProjectionX(Form("histinc%s%d%s", sample.c_str(), i, cr.first.c_str()), i+1);
                SetHistogramColor(thisHist, sample);
                ThisRsqSliceMCMap[sample] = thisHist;
            }
            for(auto &name : cr.second){
                ThisRsqSliceMC->Add(ThisRsqSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", "MRInclusiveSlice"+to_string(i)+cr.first, true);
        }
        for(int i = 0; i < nMRBins; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorHistosData[cr.first]->ProjectionY(Form("ThisMRSliceIncData%d%s", i, cr.first.c_str()), i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (MR > %.0f), %s Box", MRBinLowEdges[i], cr.first.c_str()));
            for(auto &sample : cr.second){
                TH1F *thisHist;
                thisHist = (TH1F*)razorHistosMC[cr.first][sample]->ProjectionY(Form("histinc%s%d%s", sample.c_str(), i, cr.first.c_str()), i+1);
                SetHistogramColor(thisHist, sample);
                ThisMRSliceMCMap[sample] = thisHist;
            }
            for(auto &name : cr.second){
                ThisMRSliceMC->Add(ThisMRSliceMCMap[name]);
            }
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", "RsqInclusiveSlice"+to_string(i)+cr.first, true);
        }

        delete RazorLegend;
    }
}

int main(){
    CompareDataFitCRs();
    return 0;
}
