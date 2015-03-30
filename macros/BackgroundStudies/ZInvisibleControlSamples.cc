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

using namespace std;

void ZInvisibleControlSamples(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.2f");

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

    map<string, string> cuts;
    cuts["DYJets"] = "hlt_dimuon && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    cuts["WJets"] = "hlt_singlemu && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100";
    cuts["GJets"] = "hlt_photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1";
    cuts["ZJets"] = "hlt_razor && MR > 300 && Rsq > 0.15 && numJets80 > 1";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["DYJets"] = new TFile("DYJetsRun1_19700pb.root");
    mcfiles["WJets"] = new TFile("WJetsRun1_19700pb.root");
    mcfiles["GJets"] = new TFile("GJetsRun1_19700pb.root");
    mcfiles["ZJets"] = new TFile("ZJetsRun1_19700pb.root");
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

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiency.root");
    TH2F muonLooseEffHisto = *(TH2F *)effFile.Get("MuonEfficiency");
    TH2F muonTightEffHisto = *(TH2F *)effFile.Get("MuonEfficiencyTight");
    TH2F photonEffHisto = *(TH2F *)effFile.Get("PhotonEfficiency");
    TH2F zAccHisto = *(TH2F *)effFile.Get("MuonAcceptance");
    
    //Step 1: Get the distributions to reweigh by: MET, MR, Rsq
    map<string, TH1F> metHistosForReweighing;
    map<string, TH2F> razorHistosForReweighing;
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        metHistosForReweighing[tree.first] = TH1F(Form("metmc%s", tree.first.c_str()), "MET (GeV); MET(GeV)", nMetBins, MetMin, MetMax);
        razorHistosForReweighing[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
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
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) eventWeight /= effFactor;
                else{ 
                    eventWeight = 0;
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) eventWeight /= effFactor;
                else{
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                    eventWeight = 0;
                }
            }
            //fill each quantity
            metHistosForReweighing[tree.first].Fill(mets[tree.first], eventWeight);
            razorHistosForReweighing[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
        }
    }

    //Step 2: Sanity check: apply the reweighing factors to MC
    map<string, TH2F> razorHistosMC;
    for(auto &tree : mctrees){
        cout << "Filling reweighed MC histograms: " << tree.first << endl;
        razorHistosMC[tree.first] = TH2F(Form("razorMCReweighed%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
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
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{ 
                    reweighFactor = 0;
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
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

    //Step 3: Apply the reweighing factors to data
    map<string, TH2F> razorHistosData;
    map<string, TH2F> razorHistosDataBeforeReweighing; //apply only efficiency and acceptance corrections
    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        razorHistosDataBeforeReweighing[tree.first] = TH2F(Form("razordatabeforereweighing%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
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
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "WJets"){
                double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{ 
                    reweighFactor = 0;
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                }
            }
            else if(tree.first == "DYJets"){
                double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                effFactor *= zAccHisto.GetBinContent(zAccHisto.FindBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                if(effFactor > 1e-5) reweighFactor /= effFactor;
                else{
                    cout << "Warning: efficiency histogram gives 0; setting event weight to 0" << endl;
                    reweighFactor = 0;
                }
            }
            razorHistosDataBeforeReweighing[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);

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
            razorHistosData[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);
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
        c.Print(Form("controlSampleHistogramBeforeReweighing%s.pdf", hist.first.c_str()));
        c.Print(Form("controlSampleHistogramBeforeReweighing%s.root", hist.first.c_str()));
        hist.second.Write();
    }

    //quantify agreement between DYJets and WJets predictions
    c.SetLogz(false);
    TH2F *DYWComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYWComparisonHist");
    for(int i = 0; i < DYWComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYWComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (WJets - DYJets)/DYJets
            DYWComparisonHist->SetBinContent(i, j, (razorHistosData["WJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
        }
    }
    DYWComparisonHist->SetTitle("(WJets Prediction - DYJets Prediction)/DYJets Prediction");
    DYWComparisonHist->GetXaxis()->SetTitle("MR");
    DYWComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYWComparisonHist->SetStats(0);
    DYWComparisonHist->Draw("colz");
    DYWComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonDYW.pdf");
    c.Print("controlSampleHistogramComparisonDYW.root");
    DYWComparisonHist->Write();

    //do the same for DYJets vs GJets
    TH2F *DYGComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYGComparisonHist");
    for(int i = 0; i < DYGComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYGComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - DYJets)/DYJets
            DYGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
        }
    }
    DYGComparisonHist->SetTitle("(GJets Prediction - DYJets Prediction)/DYJets Prediction");
    DYGComparisonHist->GetXaxis()->SetTitle("MR");
    DYGComparisonHist->GetYaxis()->SetTitle("Rsq");
    DYGComparisonHist->SetStats(0);
    DYGComparisonHist->Draw("colz");
    DYGComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonDYG.pdf");
    c.Print("controlSampleHistogramComparisonDYG.root");
    DYGComparisonHist->Write();

    //and for WJets vs GJets
    TH2F *WGComparisonHist = (TH2F*)razorHistosData["WJets"].Clone("WGComparisonHist");
    for(int i = 0; i < WGComparisonHist->GetNbinsX()+1; i++){
        for(int j = 0; j < WGComparisonHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - WJets)/WJets
            WGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["WJets"].GetBinContent(i, j))/razorHistosData["WJets"].GetBinContent(i, j));
        }
    }
    WGComparisonHist->SetTitle("(GJets Prediction - WJets Prediction)/WJets Prediction");
    WGComparisonHist->GetXaxis()->SetTitle("MR");
    WGComparisonHist->GetYaxis()->SetTitle("Rsq");
    WGComparisonHist->SetStats(0);
    WGComparisonHist->Draw("colz");
    WGComparisonHist->Draw("same,text");
    c.Print("controlSampleHistogramComparisonWG.pdf");
    c.Print("controlSampleHistogramComparisonWG.root");
    WGComparisonHist->Write();
}

int main(){
    ZInvisibleControlSamples();
    return 0;
}
