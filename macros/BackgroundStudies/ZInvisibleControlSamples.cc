//Macro to predict the MR and Rsq distributions of the Z->nu nu background in the razor search using Z->mu mu, W->mu nu, and Gamma+Jets events

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

using namespace std;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void ZInvisibleControlSamples(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //choose which sample to normalize to
    string normalizeTo = "ZJets";

    //for plots
    float MetMin = 0.;
    float MetMax = 1000;
    float nMetBins = 20;
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.80, 1.5};

    //upper bounds of reweighing histograms
    float maxPhotonPt = 999; 
    float maxMuonPt = 999;
    float maxZPt = 2999;

    //decide to reweigh by MET or by MR and Rsq
    //(reweighByRazor = false to reweigh by MET, true to reweigh by MR and Rsq)
    bool reweighByRazor = true; 

    map<string, string> suffixes;
    suffixes["DYJets"] = "_noZ";
    suffixes["WJets"] = "_noW";
    suffixes["GJets"] = "_noPho";
    suffixes["ZJets"] = "";
    suffixes["Top"] = "_noW";
    suffixes["TopForDY"] = "_noZ";
    suffixes["EMQCD"] = "_noPho";

    map<string, string> cuts;
    cuts["DYJets"] = "nLooseMuons == 2 && hlt_dimuon && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    cuts["WJets"] = "nBTaggedJets == 0 && nTightMuons == 1 && nLooseMuons == 1 && hlt_singlemu && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100";
    cuts["GJets"] = "hlt_photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1";
    cuts["ZJets"] = "nLooseMuons == 0 && nLooseElectrons == 0 && hlt_razor && MR > 300 && Rsq > 0.15 && numJets80 > 1";
    cuts["Top"] = cuts["WJets"];
    cuts["TopForDY"] = cuts["DYJets"];
    cuts["EMQCD"] = cuts["GJets"];

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["DYJets"] = new TFile("DYJetsRun1_19700pb.root");
    mcfiles["WJets"] = new TFile("WJetsRun1_19700pb.root");
    mcfiles["GJets"] = new TFile("GJetsRun1_19700pb.root");
    mcfiles["ZJets"] = new TFile("ZJetsRun1_19700pb.root");
    mcfiles["Top"] = new TFile("TopBackgroundsRun1_19700pb.root");
    mcfiles["TopForDY"] = mcfiles["Top"];
    mcfiles["EMQCD"] = new TFile("PhotonBackgroundsRun1_19700pb.root");
    datafiles["DYJets"] = new TFile("DoubleMuRun1.root");
    datafiles["WJets"] = new TFile("SingleMuRun1.root");
    datafiles["GJets"] = new TFile("PhotonRun1.root");
    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, float> mets;
    map<string, float> mrs;
    map<string, float> rsqs;
    float weight;
    float leadingMuonPt, leadingMuonEta, leadingPhotonPt, leadingPhotonEta, recoZpt, recoZeta, recoZmass, subleadingMuonPt, subleadingMuonEta, mTLepMet;
    float leadingTightMuonPt, leadingTightMuonEta;
    float hlt_photon_weight;
    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        mctrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingTightMuonEta", &leadingTightMuonEta);
        mctrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        mctrees[file.first]->SetBranchAddress("subleadingMuonEta", &subleadingMuonEta);
        mctrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        mctrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        mctrees[file.first]->SetBranchAddress("recoZpt", &recoZpt);
        mctrees[file.first]->SetBranchAddress("recoZeta", &recoZeta);
        mctrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        mctrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
    }
    //mctrees["GJets"]->SetBranchAddress("hlt_photon_weight", &hlt_photon_weight);
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        datatrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        datatrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonEta", &leadingTightMuonEta);
        datatrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        datatrees[file.first]->SetBranchAddress("subleadingMuonEta", &subleadingMuonEta);
        datatrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        datatrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        datatrees[file.first]->SetBranchAddress("recoZpt", &recoZpt);
        datatrees[file.first]->SetBranchAddress("recoZeta", &recoZeta);
        datatrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        datatrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
    }
        datatrees["GJets"]->SetBranchAddress("hlt_photon_weight", &hlt_photon_weight);

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiency.root");
    TH2F muonLooseEffHisto = *(TH2F *)effFile.Get("MuonEfficiency");
    TH2F muonTightEffHisto = *(TH2F *)effFile.Get("MuonEfficiencyTight");
    TH2F photonEffHisto = *(TH2F *)effFile.Get("PhotonEfficiency");
    TH2F zAccHisto = *(TH2F *)effFile.Get("MuonAcceptance");
    
    //Step 1: Get the distributions to reweigh by: MET, MR, Rsq
    map<string, TH1F> metHistosForReweighing;
    map<string, TH2F> razorHistosForReweighing;
    map<string, TH1F> MRHistosBeforeReweighing;
    map<string, TH1F> RsqHistosBeforeReweighing;
    TH1F mcPhotonPt("mcPhotonPt", "Photon Pt; photon pt", 100, 0, 1000);
    mcPhotonPt.Sumw2();
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        metHistosForReweighing[tree.first] = TH1F(Form("metmc%s", tree.first.c_str()), "MET (GeV); MET(GeV)", nMetBins, MetMin, MetMax);
        razorHistosForReweighing[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        MRHistosBeforeReweighing[tree.first] = TH1F(Form("mrmc%s", tree.first.c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        RsqHistosBeforeReweighing[tree.first] = TH1F(Form("rsqmc%s", tree.first.c_str()), "; Rsq", nRsqBins, RsqBinLowEdges);
        MRHistosBeforeReweighing[tree.first].Sumw2();
        RsqHistosBeforeReweighing[tree.first].Sumw2();
        razorHistosForReweighing[tree.first].Sumw2();
        uint nEntries = tree.second->GetEntries();
        //make TTreeFormula for selection cuts
        TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
        cutsFormula.GetNdata();
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            //apply selection cuts
            bool passesSelection = cutsFormula.EvalInstance();
            if(!passesSelection) continue;

            float eventWeight = weight;
            //reweigh according to selection efficiency and acceptance
            if(tree.first == "GJets"){
                double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                if(effFactor > 1e-5) eventWeight /= effFactor;
                else{ 
                    eventWeight = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                }
                mcPhotonPt.Fill(leadingPhotonPt, eventWeight);
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) eventWeight /= effFactor;
                else{ 
                    eventWeight = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) eventWeight /= effFactor;
                else{
                    //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                    eventWeight = 0;
                }
            }
            //fill each quantity
            metHistosForReweighing[tree.first].Fill(mets[tree.first], eventWeight);
            razorHistosForReweighing[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
            MRHistosBeforeReweighing[tree.first].Fill(mrs[tree.first], eventWeight);
            RsqHistosBeforeReweighing[tree.first].Fill(rsqs[tree.first], eventWeight);
        }
    }

    //Step 2: Sanity check: apply the reweighing factors to MC
    map<string, TH2F> razorHistosMC;
    for(auto &tree : mctrees){
        cout << "Filling reweighed MC histograms: " << tree.first << endl;
        razorHistosMC[tree.first] = TH2F(Form("razorMCReweighed%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        razorHistosMC[tree.first].Sumw2();
        uint nEntries = tree.second->GetEntries();
        TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
        cutsFormula.GetNdata();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i);

            //apply selection cuts
            bool passesSelection = cutsFormula.EvalInstance();
            if(!passesSelection) continue;

            float reweighFactor = weight;
            //reweigh by efficiency and acceptance
            if(tree.first == "GJets"){
                double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{
                    reweighFactor = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{ 
                    reweighFactor = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{
                    //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                    reweighFactor = 0;
                }
            }

            if(reweighByRazor){ //reweigh by MR and Rsq
                //get the factor to reweigh by
                float denominator = razorHistosForReweighing[tree.first].GetBinContent(razorHistosForReweighing[tree.first].FindBin(mrs[tree.first], rsqs[tree.first]));
                float numerator = razorHistosForReweighing[normalizeTo].GetBinContent(razorHistosForReweighing[normalizeTo].FindBin(mrs[tree.first], rsqs[tree.first]));
                if(denominator > 0){
                    reweighFactor *= numerator / denominator;
                }
            } 
            else{ //reweigh by MET
                //get the factor to reweigh by
                float denominator = metHistosForReweighing[tree.first].GetBinContent(metHistosForReweighing[tree.first].FindBin(mets[tree.first]));
                float numerator = metHistosForReweighing[normalizeTo].GetBinContent(metHistosForReweighing[normalizeTo].FindBin(mets[tree.first]));
                if(denominator > 0){
                    reweighFactor *= numerator / denominator;    
                }
            }
            razorHistosMC[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);
        }
    }

    //Step 3: make data distributions
    map<string, TH2F> razorHistosDataBeforeReweighing; //apply only efficiency and acceptance corrections
    map<string, TH1F> MRHistosDataBeforeReweighing;
    map<string, TH1F> RsqHistosDataBeforeReweighing;
    TH1F dataPhotonPt("dataPhotonPt", "Photon Pt; photon pt", 100, 0, 1000);
    dataPhotonPt.Sumw2();
    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        razorHistosDataBeforeReweighing[tree.first] = TH2F(Form("razordatabeforereweighing%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        MRHistosDataBeforeReweighing[tree.first] = TH1F(Form("mrdata%s", tree.first.c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        RsqHistosDataBeforeReweighing[tree.first] = TH1F(Form("rsqdata%s", tree.first.c_str()), "; Rsq (GeV)", nRsqBins, RsqBinLowEdges);
        razorHistosDataBeforeReweighing[tree.first].Sumw2();
        MRHistosDataBeforeReweighing[tree.first].Sumw2();
        RsqHistosDataBeforeReweighing[tree.first].Sumw2();
        uint nEntries = tree.second->GetEntries();
        TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
        cutsFormula.GetNdata();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i);

            //apply selection cuts
            bool passesSelection = cutsFormula.EvalInstance();
            if(!passesSelection) continue;

            float reweighFactor = 1.0;
            //reweigh by efficiency and acceptance
            if(tree.first == "GJets"){
                double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{ 
                    reweighFactor = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                }
                //multiply by trigger weight
                reweighFactor *= hlt_photon_weight;
                dataPhotonPt.Fill(leadingPhotonPt, reweighFactor);
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{ 
                    reweighFactor = 0;
                    //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{
                    //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                    reweighFactor = 0;
                }
            }
            razorHistosDataBeforeReweighing[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);
            MRHistosDataBeforeReweighing[tree.first].Fill(mrs[tree.first], reweighFactor);
            RsqHistosDataBeforeReweighing[tree.first].Fill(rsqs[tree.first], reweighFactor);
        }
    }
    //subtract the top background from the WJets histogram
    razorHistosDataBeforeReweighing["WJets"] = razorHistosDataBeforeReweighing["WJets"] - razorHistosForReweighing["Top"];
    //subtract the top background from the DYJets histogram
    razorHistosDataBeforeReweighing["DYJets"] = razorHistosDataBeforeReweighing["DYJets"] - razorHistosForReweighing["TopForDY"];
    //subtract the QCD background from the photon+jets histogram
    razorHistosDataBeforeReweighing["GJets"] = razorHistosDataBeforeReweighing["GJets"] - razorHistosForReweighing["EMQCD"];
    map<string, TH2F> razorHistosData;
    //rescale the WJets MC histogram to the WJets data histogram
    razorHistosForReweighing["WJets"].Scale(razorHistosDataBeforeReweighing["WJets"].Integral()*1.0/razorHistosForReweighing["WJets"].Integral());

    //Step 4: apply reweighing factors to data
    if(reweighByRazor){
        for(auto &tree : datatrees){
            cout << "Making weighted data histograms: " << tree.first << endl;
            razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosData[tree.first].Sumw2();
            for(int i = 0; i < razorHistosData[tree.first].GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosData[tree.first].GetNbinsY()+1; j++){
                    float numerator = razorHistosForReweighing[normalizeTo].GetBinContent(i, j);
                    float denominator = razorHistosForReweighing[tree.first].GetBinContent(i, j);
                    float numeratorError = razorHistosForReweighing[normalizeTo].GetBinError(i, j);
                    float denominatorError = razorHistosForReweighing[tree.first].GetBinError(i, j);
                    if(denominator > 0){
                        razorHistosData[tree.first].SetBinContent(i, j, razorHistosDataBeforeReweighing[tree.first].GetBinContent(i, j)*numerator/denominator);
                        //compute the uncertainty on the bin, folding in uncertainties on the scale factor numerator/denominator
                        razorHistosData[tree.first].SetBinError(i, j, sqrt(pow(razorHistosDataBeforeReweighing[tree.first].GetBinError(i, j)*numerator/denominator, 2) + pow(razorHistosDataBeforeReweighing[tree.first].GetBinContent(i, j)*numeratorError/denominator, 2) + pow(razorHistosDataBeforeReweighing[tree.first].GetBinContent(i, j)*numerator*denominatorError/(denominator*denominator), 2)));
                    }
                    else{
                        razorHistosData[tree.first].SetBinContent(i, j, 0.);
                        razorHistosData[tree.first].SetBinError(i, j, 0.);
                    }
                }
            }   
        }
    }
    else{ 
        for(auto &tree : datatrees){
            cout << "Making weighted data histograms: " << tree.first << endl;
            razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosData[tree.first].Sumw2();

            uint nEntries = tree.second->GetEntries();
            TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
            cutsFormula.GetNdata();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                tree.second->GetEntry(i);

                //apply selection cuts
                bool passesSelection = cutsFormula.EvalInstance();
                if(!passesSelection) continue;

                float reweighFactor = 1.0;
                //reweigh by efficiency and acceptance
                if(tree.first == "GJets"){
                    double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                    if(effFactor > 1e-5) reweighFactor /= effFactor;
                    else{ 
                        reweighFactor = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                    }
                    //multiply by trigger weight
                    reweighFactor *= hlt_photon_weight;
                }
                else if(tree.first == "WJets"){
                    double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                    if(effFactor > 1e-5) reweighFactor /= effFactor;
                    else{ 
                        reweighFactor = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                    }
                }
                else if(tree.first == "DYJets"){
                    double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                    effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                    effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                    if(effFactor > 1e-5) reweighFactor /= effFactor;
                    else{
                        //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                        reweighFactor = 0;
                    }
                }
                //reweigh by MET
                //get the factor to reweigh by
                float denominator = metHistosForReweighing[tree.first].GetBinContent(metHistosForReweighing[tree.first].FindBin(mets[tree.first]));
                float numerator = metHistosForReweighing[normalizeTo].GetBinContent(metHistosForReweighing[normalizeTo].FindBin(mets[tree.first]));
                if(denominator > 0){
                    reweighFactor *= numerator / denominator;    
                }

                razorHistosData[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);
            }
        }
    }

    TFile outfile("controlSampleHistograms.root", "recreate");
    TCanvas c("c", "c", 800, 600);
    c.SetLogx();
    //print "step 1" histograms used for reweighing
    c.SetLogz();
    for(auto &hist : razorHistosForReweighing){
        cout << "Before reweighing, MC histo for " << hist.first << " entries: " << hist.second.Integral() << endl;
    }
    for(auto &hist : razorHistosDataBeforeReweighing){
        cout << "Before reweighing, Data histo for " << hist.first << " entries: " << hist.second.Integral() << endl;
    }
    for(auto &hist : razorHistosMC){
        cout << "After reweighing, MC histo for " << hist.first << " entries: " << hist.second.Integral() << endl;
    }
    for(auto &hist : razorHistosData){
        cout << "After reweighing, Data histo for " << hist.first << " entries: " << hist.second.Integral() << endl;
    }
    for(auto &hist : razorHistosForReweighing){
        hist.second.SetTitle(Form("MC reweighed by efficiency and acceptance, for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleMCHistogram%s.pdf", hist.first.c_str()));
        c.Print(Form("controlSampleMCHistogram%s.root", hist.first.c_str()));
        hist.second.Write();
    }
    //print "step 2" histograms
    for(auto &hist : razorHistosMC){
        hist.second.SetTitle(Form("MC with all corrections applied, for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleReweighedMCHistogram%s.pdf", hist.first.c_str()));
        c.Print(Form("controlSampleReweighedMCHistogram%s.root", hist.first.c_str()));
        hist.second.Write();
    }
    //print "step 3" razor histograms from data
    for(auto &hist : razorHistosData){
        hist.second.SetTitle(Form("Prediction for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleHistogram%s.pdf", hist.first.c_str()));
        c.Print(Form("controlSampleHistogram%s.root", hist.first.c_str()));
        hist.second.Write();
    }
    //print razor histograms from data before reweighing by MR/Rsq/MET
    for(auto &hist : razorHistosDataBeforeReweighing){
        hist.second.SetTitle(Form("Data before reweighing by MR and Rsq, for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleHistogramBeforeReweighing%s.pdf", hist.first.c_str()));
        c.Print(Form("controlSampleHistogramBeforeReweighing%s.root", hist.first.c_str()));
        hist.second.Write();
    }
    //print MR histograms, comparing data to MC
    c.SetLogy();
    //WJets
    THStack SingleMuonMC("SingleMuonMC", "MR in 1-muon control sample");
    MRHistosBeforeReweighing["WJets"].SetLineColor(kOrange+10);
    MRHistosBeforeReweighing["WJets"].SetFillColor(kOrange+10);
    MRHistosBeforeReweighing["Top"].SetLineColor(kViolet-5);
    MRHistosBeforeReweighing["Top"].SetFillColor(kViolet-5);
    SingleMuonMC.Add(&MRHistosBeforeReweighing["Top"]);
    SingleMuonMC.Add(&MRHistosBeforeReweighing["WJets"]);
    MRHistosDataBeforeReweighing["WJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighing["WJets"].SetMarkerSize(1);
    TLegend *SingleMuonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighing["WJets"], "WJets MC");
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighing["Top"], "Top MC (ttbar and single top)");
    SingleMuonLegend->AddEntry(&MRHistosDataBeforeReweighing["WJets"], "2012 Data, Single Muon CS");
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighing["WJets"], &SingleMuonMC, SingleMuonLegend, "MR (GeV)", "controlSampleMRBackgroundSingleMuon", true);
    //DYJets
    THStack DoubleMuonMC("DoubleMuonMC", "MR in 2-muon control sample");
    MRHistosBeforeReweighing["DYJets"].SetLineColor(kAzure);
    MRHistosBeforeReweighing["DYJets"].SetFillColor(kAzure);
    MRHistosBeforeReweighing["TopForDY"].SetLineColor(kViolet-5);
    MRHistosBeforeReweighing["TopForDY"].SetFillColor(kViolet-5);
    DoubleMuonMC.Add(&MRHistosBeforeReweighing["TopForDY"]);
    DoubleMuonMC.Add(&MRHistosBeforeReweighing["DYJets"]);
    MRHistosDataBeforeReweighing["DYJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighing["DYJets"].SetMarkerSize(1);
    TLegend *DoubleMuonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighing["DYJets"], "DYJets MC");
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighing["TopForDY"], "Top MC (ttbar and single top)");
    DoubleMuonLegend->AddEntry(&MRHistosDataBeforeReweighing["DYJets"], "2012 Data, Double Muon CS");
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighing["DYJets"], &DoubleMuonMC, DoubleMuonLegend, "MR (GeV)", "controlSampleMRBackgroundDoubleMuon", true);
    //Gamma+Jets
    THStack PhotonMC("PhotonMC", "MR in photon control sample");
    MRHistosBeforeReweighing["GJets"].SetLineColor(kTeal+10);
    MRHistosBeforeReweighing["GJets"].SetFillColor(kTeal+10);
    MRHistosBeforeReweighing["EMQCD"].SetLineColor(kGreen+3);
    MRHistosBeforeReweighing["EMQCD"].SetFillColor(kGreen+3);
    PhotonMC.Add(&MRHistosBeforeReweighing["EMQCD"]);
    PhotonMC.Add(&MRHistosBeforeReweighing["GJets"]);
    MRHistosDataBeforeReweighing["GJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighing["GJets"].SetMarkerSize(1);
    TLegend *PhotonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    PhotonLegend->AddEntry(&MRHistosBeforeReweighing["GJets"], "GJets MC");
    PhotonLegend->AddEntry(&MRHistosBeforeReweighing["EMQCD"], "QCD MC (EM enriched)");
    PhotonLegend->AddEntry(&MRHistosDataBeforeReweighing["GJets"], "2012 Data, Photon CS");
    PhotonLegend->Draw();
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighing["GJets"], &PhotonMC, PhotonLegend, "MR (GeV)", "controlSampleMRBackgroundPhoton", true);

    //print Rsq histograms, comparing data to MC
    c.SetLogy();
    c.SetLogx(kFALSE);
    //WJets
    THStack SingleMuonRsqMC("SingleMuonMC", "Rsq in 1-muon control sample");
    RsqHistosBeforeReweighing["WJets"].SetLineColor(kOrange+10);
    RsqHistosBeforeReweighing["WJets"].SetFillColor(kOrange+10);
    RsqHistosBeforeReweighing["Top"].SetLineColor(kViolet-5);
    RsqHistosBeforeReweighing["Top"].SetFillColor(kViolet-5);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighing["Top"]);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighing["WJets"]);
    RsqHistosDataBeforeReweighing["WJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighing["WJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighing["WJets"], &SingleMuonRsqMC, SingleMuonLegend, "Rsq", "controlSampleRsqBackgroundSingleMuon", false);
    //DYJets
    THStack DoubleMuonRsqMC("DoubleMuonMC", "Rsq in 2-muon control sample");
    RsqHistosBeforeReweighing["DYJets"].SetLineColor(kAzure);
    RsqHistosBeforeReweighing["DYJets"].SetFillColor(kAzure);
    RsqHistosBeforeReweighing["TopForDY"].SetLineColor(kViolet-5);
    RsqHistosBeforeReweighing["TopForDY"].SetFillColor(kViolet-5);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighing["TopForDY"]);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighing["DYJets"]);
    RsqHistosDataBeforeReweighing["DYJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighing["DYJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighing["DYJets"], &DoubleMuonRsqMC, DoubleMuonLegend, "Rsq", "controlSampleRsqBackgroundDoubleMuon", false);
    //Gamma+Jets
    THStack PhotonRsqMC("PhotonMC", "Rsq in photon control sample");
    RsqHistosBeforeReweighing["GJets"].SetLineColor(kTeal+10);
    RsqHistosBeforeReweighing["GJets"].SetFillColor(kTeal+10);
    RsqHistosBeforeReweighing["EMQCD"].SetLineColor(kGreen+3);
    RsqHistosBeforeReweighing["EMQCD"].SetFillColor(kGreen+3);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighing["EMQCD"]);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighing["GJets"]);
    RsqHistosDataBeforeReweighing["GJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighing["GJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighing["GJets"], &PhotonRsqMC, PhotonLegend, "Rsq", "controlSampleRsqBackgroundPhoton", false);

    //quantify agreement between DYJets and WJets predictions
    c.SetLogy(false);
    c.SetLogz(false);
    c.SetLogx();
    gStyle->SetPaintTextFormat("1.2f");
    TH2F *DYWComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYWComparisonHist");
    TH2F *DYWSigmaComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYWSigmaComparisonHist");
    for(int i = 0; i < DYWComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYWComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (WJets - DYJets)/DYJets
            DYWComparisonHist->SetBinContent(i, j, (razorHistosData["WJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
            //set bin content to (WJets - DYJets)/(error on difference)
            float sigma1 = razorHistosData["WJets"].GetBinError(i, j);
            float sigma2 = razorHistosData["DYJets"].GetBinError(i, j);
            DYWSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["WJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
        }
    }
    DYWComparisonHist->SetTitle("(WJets Prediction - DYJets Prediction)/DYJets Prediction");
    DYWComparisonHist->GetXaxis()->SetTitle("MR");
    DYWComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYWComparisonHist->SetStats(0);
    DYWComparisonHist->SetMinimum(-1.0);
    DYWComparisonHist->SetMaximum(1.0);
    DYWComparisonHist->Draw("colz");
    DYWComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonDYW.pdf");
    c.Print("controlSampleHistogramComparisonDYW.root");
    DYWComparisonHist->Write();
    DYWSigmaComparisonHist->SetTitle("(WJets Prediction - DYJets Prediction)/#sigma_{W - DY}");
    DYWSigmaComparisonHist->GetXaxis()->SetTitle("MR");
    DYWSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYWSigmaComparisonHist->SetStats(0);
    DYWSigmaComparisonHist->SetMinimum(-3);
    DYWSigmaComparisonHist->SetMaximum(3);
    DYWSigmaComparisonHist->Draw("colz");
    DYWSigmaComparisonHist->Draw("same,text");
    c.Print("controlSampleSigmaHistogramComparisonDYW.pdf");
    c.Print("controlSampleSigmaHistogramComparisonDYW.root");
    DYWSigmaComparisonHist->Write();

    //do the same for DYJets vs GJets
    TH2F *DYGComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYGComparisonHist");
    TH2F *DYGSigmaComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYGSigmaComparisonHist");
    for(int i = 0; i < DYGComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYGComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - DYJets)/DYJets
            DYGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
            //set bin content to (GJets - DYJets)/(error on difference)
            float sigma1 = razorHistosData["GJets"].GetBinError(i, j);
            float sigma2 = razorHistosData["DYJets"].GetBinError(i, j);
            DYGSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
        }
    }
    DYGComparisonHist->SetTitle("(GJets Prediction - DYJets Prediction)/DYJets Prediction");
    DYGComparisonHist->GetXaxis()->SetTitle("MR");
    DYGComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYGComparisonHist->SetStats(0);
    DYGComparisonHist->SetMinimum(-1.0);
    DYGComparisonHist->SetMaximum(1.0);
    DYGComparisonHist->Draw("colz");
    DYGComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonDYG.pdf");
    c.Print("controlSampleHistogramComparisonDYG.root");
    DYGComparisonHist->Write();
    DYGSigmaComparisonHist->SetTitle("(GJets Prediction - DYJets Prediction)/#sigma_{G - DY}");
    DYGSigmaComparisonHist->GetXaxis()->SetTitle("MR");
    DYGSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYGSigmaComparisonHist->SetStats(0);
    DYGSigmaComparisonHist->SetMinimum(-3);
    DYGSigmaComparisonHist->SetMaximum(3);
    DYGSigmaComparisonHist->Draw("colz");
    DYGSigmaComparisonHist->Draw("same,text");
    c.Print("controlSampleSigmaHistogramComparisonDYG.pdf");
    c.Print("controlSampleSigmaHistogramComparisonDYG.root");
    DYGSigmaComparisonHist->Write();

    //and for WJets vs GJets
    TH2F *WGComparisonHist = (TH2F*)razorHistosData["WJets"].Clone("WGComparisonHist");
    TH2F *WGSigmaComparisonHist = (TH2F*)razorHistosData["WJets"].Clone("WGSigmaComparisonHist");
    for(int i = 0; i < WGComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < WGComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - WJets)/WJets
            WGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["WJets"].GetBinContent(i, j))/razorHistosData["WJets"].GetBinContent(i, j));
            //set bin content to (GJets - WJets)/(error on difference)
            float sigma1 = razorHistosData["GJets"].GetBinError(i, j);
            float sigma2 = razorHistosData["WJets"].GetBinError(i, j);
            WGSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["WJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
        }
    }
    WGComparisonHist->SetTitle("(GJets Prediction - WJets Prediction)/WJets Prediction");
    WGComparisonHist->GetXaxis()->SetTitle("MR");
    WGComparisonHist->GetYaxis()->SetTitle("Rsq");
    WGComparisonHist->SetStats(0);
    WGComparisonHist->SetMinimum(-1.0);
    WGComparisonHist->SetMaximum(1.0);
    WGComparisonHist->Draw("colz");
    WGComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonWG.pdf");
    c.Print("controlSampleHistogramComparisonWG.root");
    WGComparisonHist->Write();
    WGSigmaComparisonHist->SetTitle("(GJets Prediction - WJets Prediction)/#sigma_{G - W}");
    WGSigmaComparisonHist->GetXaxis()->SetTitle("MR");
    WGSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
    WGSigmaComparisonHist->SetStats(0);
    WGSigmaComparisonHist->SetMinimum(-3);
    WGSigmaComparisonHist->SetMaximum(3);
    WGSigmaComparisonHist->Draw("colz");
    WGSigmaComparisonHist->Draw("same,text");
    c.Print("controlSampleSigmaHistogramComparisonWG.pdf");
    c.Print("controlSampleSigmaHistogramComparisonWG.root");
    WGSigmaComparisonHist->Write();

    //plot the photon pt distribution in data and MC
    mcPhotonPt.SetStats(0);
    mcPhotonPt.SetLineColor(kViolet);
    dataPhotonPt.SetStats(0);
    mcPhotonPt.Draw();
    dataPhotonPt.Draw("pesame");
    c.Print("controlSamplePhotonPt.pdf");
    c.Print("controlSamplePhotonPt.root");
}

int main(){
    ZInvisibleControlSamples();
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
    mcTotal->Sumw2();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();
    dataOverMC->Sumw2();
    dataOverMC->Divide(mcTotal);
    dataOverMC->GetXaxis()->SetTitle(xaxisTitle.c_str());
    dataOverMC->GetYaxis()->SetTitle("Data / MC");
    //dataOverMC->SetMinimum(0.7);
    //dataOverMC->SetMaximum(1.3);
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
    for(int i = 0; i < dataOverMC->GetNbinsX()+1; i++) cout << dataOverMC->GetBinContent(i) << endl;
    pad2.Modified();
    gPad->Update();
    c.Print(Form("%s.pdf", printString.c_str()));
    c.Print(Form("%s.root", printString.c_str()));
}
