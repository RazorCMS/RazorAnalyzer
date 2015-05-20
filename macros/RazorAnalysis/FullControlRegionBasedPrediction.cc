//Runs on the output of the RazorInclusive analyzer and gives the MC-based background prediction in each bin of the MR-Rsq plane

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

void FullControlRegionBasedPrediction(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 101;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //define MR and Rsq binning
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.8, 1.5};

    //get input files -- output of RazorInclusive analyzer
    //NOTE: all data-MC correction factors should already be applied EXCEPT for the hadronic recoil scale factors obtained from the control regions
    map<string, TFile*> mcfiles;
    TFile *datafile;
    //main backgrounds
    mcfiles["DYJets"] = new TFile("RazorDYJetsRun1_19700pb_weighted.root");
    mcfiles["WJets"] = new TFile("RazorWJetsRun1_19700pb_weighted.root");
    mcfiles["ZJetsNuNu"] = new TFile("RazorZJetsNuNuRun1_19700pb_weighted.root");
    mcfiles["TTJets"] = new TFile("RazorTTJetsRun1_19700pb_weighted.root");
    mcfiles["SingleTop"] = new TFile("RazorSingleTopRun1_19700pb_weighted.root");
    mcfiles["QCD"] = new TFile("RazorQCDRun1_19700pb_weighted.root");
    //rare backgrounds
    mcfiles["TTW"] = new TFile("RazorTTWJetsRun1_19700pb_weighted.root");
    mcfiles["TTZ"] = new TFile("RazorTTZJetsRun1_19700pb_weighted.root");
    //TODO: add all other background processes!

    //data
    datafile = new TFile("RazorInclusiveRun1Data.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    TTree *datatree;
    float weight;
    float MR, Rsq, dPhiRazor;
    int nBTaggedJets, box;
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("box", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("dPhiRazor", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("box", &box);
        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("dPhiRazor", &dPhiRazor);
        mctrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
    }
    datatree = (TTree*)datafile->Get("RazorInclusive");
    datatree->SetBranchStatus("*", 0);
    datatree->SetBranchStatus("box", 1);
    datatree->SetBranchStatus("MR", 1);
    datatree->SetBranchStatus("Rsq", 1);
    datatree->SetBranchStatus("dPhiRazor", 1);
    datatree->SetBranchStatus("nBTaggedJets", 1);

    datatree->SetBranchAddress("box", &box);
    datatree->SetBranchAddress("MR", &MR);
    datatree->SetBranchAddress("Rsq", &Rsq);
    datatree->SetBranchAddress("dPhiRazor", &dPhiRazor);
    datatree->SetBranchAddress("nBTaggedJets", &nBTaggedJets);

    //load TTbar scale factor histograms
    TFile *SFFileTTBar = new TFile("data/ScaleFactors/Run1/TTBarDileptonScaleFactors.root");
    TH2F *SFHistTTBar = (TH2F*)SFFileTTBar->Get("TTBarDileptonScaleFactor");
    float SFmaxMRTTJets = SFHistTTBar->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqTTJets = SFHistTTBar->GetYaxis()->GetXmax() - 0.01;

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
    TH2F *SFHistZJetsNuNuFromDY = (TH2F*)SFFileZJetsNuNu->Get("DYJetsScaleFactors");
    TH2F *SFHistZJetsNuNuFromW = (TH2F*)SFFileZJetsNuNu->Get("WJetsScaleFactors");
    TH2F *SFHistZJetsNuNuFromGamma = (TH2F*)SFFileZJetsNuNu->Get("GJetsScaleFactors");
    float SFmaxMRZJetsNuNuFromDY = SFHistZJetsNuNuFromDY->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqZJetsNuNuFromDY = SFHistZJetsNuNuFromDY->GetYaxis()->GetXmax() - 0.01;
    float SFmaxMRZJetsNuNuFromW = SFHistZJetsNuNuFromW->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqZJetsNuNuFromW = SFHistZJetsNuNuFromW->GetYaxis()->GetXmax() - 0.01;
    float SFmaxMRZJetsNuNuFromGamma = SFHistZJetsNuNuFromGamma->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqZJetsNuNuFromGamma = SFHistZJetsNuNuFromGamma->GetYaxis()->GetXmax() - 0.01;

    //Step 1: Get the predictions from each MC process
    map<string, TH2F> razorHistosMC;
    map<string, TH1F> MRHistosMC;
    map<string, TH1F> RsqHistosMC;
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;

        //set up histograms
        razorHistosMC[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        MRHistosMC[tree.first] = TH1F(Form("mrmc%s", tree.first.c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        RsqHistosMC[tree.first] = TH1F(Form("rsqmc%s", tree.first.c_str()), "; Rsq", nRsqBins, RsqBinLowEdges);
        MRHistosMC[tree.first].Sumw2();
        RsqHistosMC[tree.first].Sumw2();
        razorHistosMC[tree.first].Sumw2();

        uint nEntries = tree.second->GetEntries();
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            float eventWeight = weight;

            //TTJets SF
            if(tree.first == "TTJets"){
                double SFTTJets = SFHistTTBar->GetBinContent(SFHistTTBar->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                if(SFTTJets > 1e-5){
                    eventWeight *= SFTTJets;
                }
                else{
                    cout << "Warning: TTJets scale factor is zero!" << endl;
                }
            }
            //WJets SF
            else if(tree.first == "WJets"){
                double SFWJets = SFHistWJets->GetBinContent(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                if(SFWJets > 1e-5){
                    eventWeight *= SFWJets;
                }
                else{
                    cout << "Warning: WJets scale factor is zero!" << endl;
                }
            }
            //DYJets SF
            else if(tree.first == "DYJets"){
                double SFDYJets = SFHistDYJets->GetBinContent(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                if(SFDYJets > 1e-5){
                    eventWeight *= SFDYJets;
                }
                else{
                    cout << "Warning: DYJets scale factor is zero!" << endl;
                }
            }

            //ZNuNu SF
            //TODO: combine the three predictions for ZNuNu?  currently use Gamma+Jets prediction
            else if(tree.first == "ZJetsNuNu"){
                double SFZJetsNuNu = SFHistZJetsNuNuFromGamma->GetBinContent(SFHistZJetsNuNuFromGamma->FindFixBin(min(MR, SFmaxMRZJetsNuNuFromGamma), min(Rsq, SFmaxRsqZJetsNuNuFromGamma)));
                if(SFZJetsNuNu > 1e-5){
                    eventWeight *= SFZJetsNuNu;
                }
                else{
                    cout << "Warning: ZJetsNuNu scale factor is zero!" << endl;
                }
            }

            //TODO: keep track of uncertainties correctly

            //fill each quantity
            razorHistosMC[tree.first].Fill(MR, Rsq, eventWeight);
            MRHistosMC[tree.first].Fill(MR, eventWeight);
            RsqHistosMC[tree.first].Fill(Rsq, eventWeight);
        }
    }

    //Step 2: make data distributions
    cout << "Filling data histograms" << endl;

    //create histograms
    TH2F razorData("razordata", "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
    TH1F MRData("mrdata", "; MR (GeV)", nMRBins, MRBinLowEdges);
    TH1F RsqData("rsqdata", "; Rsq (GeV)", nRsqBins, RsqBinLowEdges);
    TH2F razorDataUncorrected("razordataUncorrected", "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
    TH1F MRDataUncorrected("mrdataUncorrected", "; MR (GeV)", nMRBins, MRBinLowEdges);
    TH1F RsqDataUncorrected("rsqdataUncorrected", "; Rsq (GeV)", nRsqBins, RsqBinLowEdges);
    razorData.Sumw2();
    MRData.Sumw2();
    RsqData.Sumw2();
    razorDataUncorrected.Sumw2();
    MRDataUncorrected.Sumw2();
    RsqDataUncorrected.Sumw2();

    uint nEntries = datatree->GetEntries();
    for(uint i = 0; i < nEntries; i++){
        //get entry
        datatree->GetEntry(i);

        float eventWeight = 1.0;

        razorData.Fill(MR, Rsq, eventWeight);
        MRData.Fill(MR, eventWeight);
        RsqData.Fill(Rsq, eventWeight);
    }
    //for rare background processes, include a 20% uncertainty on the total yield in each bin, summed in quadrature with the statistical uncertainty
    double sysErrorFrac = 0.2;
    //for QCD, assign a 100% uncertainty
    double qcdErrorFrac = 1.0;
    for(auto &tree : mctrees){
        //only do this for rare processes 
        if(tree.first == "DYJets" || tree.first == "WJets" || tree.first == "ZJetsNuNu" || tree.first == "TTJets") continue; 
        for(int i = 0; i < razorHistosMC[tree.first].GetNbinsX()+1; i++){
            for(int j = 0; j < razorHistosMC[tree.first].GetNbinsY()+1; j++){
                double error = 0.0;
                if(tree.first == "QCD"){
                    error = qcdErrorFrac*razorHistosMC[tree.first].GetBinContent(i, j);
                }
                else{
                    error = sysErrorFrac*razorHistosMC[tree.first].GetBinContent(i, j);
                }
                razorHistosMC[tree.first].SetBinError(i, j, sqrt(pow(razorHistosMC[tree.first].GetBinError(i, j), 2) + error*error));
            }
        }
    }

    //make plots
    TCanvas c("c", "c", 800, 600);
    c.SetLogx();
    //print MC histograms
    c.SetLogz();
    for(auto &hist : razorHistosMC){
        hist.second.SetTitle(Form("MC for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("razorInclusiveMCHistogram%s.pdf", hist.first.c_str()));
        c.Print(Form("razorInclusiveMCHistogram%s.root", hist.first.c_str()));
    }
    //print data histogram
    razorData.SetTitle("Data");
    razorData.GetXaxis()->SetTitle("MR");
    razorData.GetYaxis()->SetTitle("Rsq");
    razorData.SetStats(0);
    razorData.Draw("colz");
    razorData.Draw("same,text");
    c.Print("razorInclusiveDataHistogram.pdf");
    c.Print("razorInclusiveDataHistogram.root");
    
    //print MR and Rsq 1D histograms, comparing data to MC
    c.SetLogy();
    THStack MRTotalRazorMC("TotalRazorMC", "MR");
    THStack RsqTotalRazorMC("TotalRazorMC", "Rsq");

    //format MC histograms
    MRHistosMC["DYJets"].SetFillColor(kAzure);
    MRHistosMC["WJets"].SetFillColor(kOrange+10);
    MRHistosMC["ZJetsNuNu"].SetFillColor(38);
    MRHistosMC["TTJets"].SetFillColor(kViolet-6);
    MRHistosMC["SingleTop"].SetFillColor(kViolet-5);
    MRHistosMC["QCD"].SetFillColor(kBlack);
    MRHistosMC["TTW"].SetFillColor(kRed+2);
    MRHistosMC["TTZ"].SetFillColor(kOrange-3);
    MRTotalRazorMC.Add(&MRHistosMC["TTZ"]);
    MRTotalRazorMC.Add(&MRHistosMC["TTW"]);
    MRTotalRazorMC.Add(&MRHistosMC["SingleTop"]);
    MRTotalRazorMC.Add(&MRHistosMC["QCD"]);
    MRTotalRazorMC.Add(&MRHistosMC["DYJets"]);
    MRTotalRazorMC.Add(&MRHistosMC["TTJets"]);
    MRTotalRazorMC.Add(&MRHistosMC["WJets"]);
    MRTotalRazorMC.Add(&MRHistosMC["ZJetsNuNu"]);
    //TODO: include all backgrounds
    MRData.SetMarkerStyle(20);
    MRData.SetMarkerSize(1);
    RsqHistosMC["DYJets"].SetFillColor(kAzure);
    RsqHistosMC["WJets"].SetFillColor(kOrange+10);
    RsqHistosMC["ZJetsNuNu"].SetFillColor(38);
    RsqHistosMC["TTJets"].SetFillColor(kViolet-6);
    RsqHistosMC["SingleTop"].SetFillColor(kViolet-5);
    RsqHistosMC["QCD"].SetFillColor(kBlack);
    RsqHistosMC["TTW"].SetFillColor(kRed+2);
    RsqHistosMC["TTZ"].SetFillColor(kOrange-3);
    RsqTotalRazorMC.Add(&RsqHistosMC["TTZ"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["TTW"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["SingleTop"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["QCD"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["DYJets"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["TTJets"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["WJets"]);
    RsqTotalRazorMC.Add(&RsqHistosMC["ZJetsNuNu"]);
    //TODO: include all backgrounds
    RsqData.SetMarkerStyle(20);
    RsqData.SetMarkerSize(1);

    //create legend
    TLegend *RazorLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    RazorLegend->AddEntry(&MRHistosMC["WJets"], "WJets MC");
    RazorLegend->AddEntry(&MRHistosMC["DYJets"], "DYJets MC");
    RazorLegend->AddEntry(&MRHistosMC["ZJetsNuNu"], "ZJetsNuNu MC");
    RazorLegend->AddEntry(&MRHistosMC["TTJets"], "TTJets MC");
    RazorLegend->AddEntry(&MRHistosMC["Top"], "Single Top MC");
    RazorLegend->AddEntry(&MRHistosMC["TTW"], "TTW MC");
    RazorLegend->AddEntry(&MRHistosMC["TTZ"], "TTZ MC");
    RazorLegend->AddEntry(&MRHistosMC["QCD"], "QCD MC");
    //TODO: include all backgrounds
    RazorLegend->AddEntry(&MRData, "2012 Data");
    DrawDataVsMCRatioPlot(&MRData, &MRTotalRazorMC, RazorLegend, "MR (GeV)", "razorInclusiveMRBackground", true);
    c.SetLogx(kFALSE);
    DrawDataVsMCRatioPlot(&RsqData, &RsqTotalRazorMC, RazorLegend, "Rsq (GeV)", "razorInclusiveRsqBackground", true);

    gStyle->SetPaintTextFormat("1.2f");
    c.SetLogy(false);
    c.SetLogz(false);
    c.SetLogx();
    //TODO: quantify data-MC agreement
}

int main(){
    FullControlRegionBasedPrediction();
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
    //mcTotal->Sumw2();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();
    //dataOverMC->Sumw2();
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
    pad2.Modified();
    gPad->Update();
    c.Print(Form("%s.pdf", printString.c_str()));
    c.Print(Form("%s.root", printString.c_str()));
}
