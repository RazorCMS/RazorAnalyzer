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
    bool doSFCorrections = true; //apply TT, W, Z, DY scale factors
    bool doMiscCorrections = true; //apply lepton efficiency, b-tagging, ... scale factors
    gROOT->SetBatch();

    float dPhiRazorCut = 4; //cut on the angle between the two razor hemispheres

    //set color palette 
    const Int_t NCont = 101;
    gStyle->SetNumberContours(NCont);

    //define MR and Rsq binning
    int nMRBins = 10;
    int nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.8, 1.5};

    //define cuts for 1D MR and Rsq plots
    float MRCutFor1DPlots = 400;
    float RsqCutFor1DPlots = 0.25;

    //get input files -- output of RazorInclusive analyzer
    int lumiInData = 19700; //in /pb
    int lumiInMC = 1; //luminosity used to normalize MC ntuples
    string mcPrefix = "";
    if(doMiscCorrections){
        //NOTE: all data-MC correction factors should already be applied EXCEPT for the hadronic recoil scale factors obtained from the control regions 
        mcPrefix = "eos/cms/store/group/phys_susy/razor/run2/RunOneRazorInclusive/done/MC_WithCorrectionFactors/";//location of MC ntuples
    }
    else{
        mcPrefix = "eos/cms/store/group/phys_susy/razor/run2/RunOneRazorInclusive/done/MC_NoCorrectionFactors/";//location of MC ntuples
    }
    string dataPrefix = "eos/cms/store/group/phys_susy/razor/run2/RunOneRazorInclusive/done/Data/"; //location of data ntuples

    map<string, TFile*> mcfiles;
    mcfiles["DYJets"] = new TFile(Form("%s/RazorInclusive_DYJetsToLL_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["WJets"] = new TFile(Form("%s/RazorInclusive_WJetsToLNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["ZJetsNuNu"] = new TFile(Form("%s/RazorInclusive_ZJetsToNuNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTJets"] = new TFile(Form("%s/RazorInclusive_TTJets_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["SingleTop"] = new TFile(Form("%s/RazorInclusive_SingleTop_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["QCD"] = new TFile(Form("%s/RazorInclusive_QCD_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTV"] = new TFile(Form("%s/RazorInclusive_TTV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["VV"] = new TFile(Form("%s/RazorInclusive_VV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTTT"] = new TFile(Form("%s/RazorInclusive_TTTT_TuneZ2star_8TeV-madgraph-tauola_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));

    //data
    TFile *datafile;
    datafile = new TFile(Form("%s/RazorInclusive_Data_HTMHTParked_Run2012_GoodLumi.root", dataPrefix.c_str()));

    //get trees and set branches
    map<string, TTree*> mctrees;
    TTree *datatree;
    float weight;
    float MR, Rsq, dPhiRazor;
    int nBTaggedJets, nSelectedJets, box;
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("box", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("dPhiRazor", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("nSelectedJets", 1);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("box", &box);
        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("dPhiRazor", &dPhiRazor);
        mctrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
        mctrees[file.first]->SetBranchAddress("nSelectedJets", &nSelectedJets);
    }
    datatree = (TTree*)datafile->Get("RazorInclusive");
    datatree->SetBranchStatus("*", 0);
    datatree->SetBranchStatus("box", 1);
    datatree->SetBranchStatus("MR", 1);
    datatree->SetBranchStatus("Rsq", 1);
    datatree->SetBranchStatus("dPhiRazor", 1);
    datatree->SetBranchStatus("nBTaggedJets", 1);
    datatree->SetBranchStatus("nSelectedJets", 1);

    datatree->SetBranchAddress("box", &box);
    datatree->SetBranchAddress("MR", &MR);
    datatree->SetBranchAddress("Rsq", &Rsq);
    datatree->SetBranchAddress("dPhiRazor", &dPhiRazor);
    datatree->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
    datatree->SetBranchAddress("nSelectedJets", &nSelectedJets);

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

    //load ZNuNu-->DYJets weighting factors
    TFile *ZNuNuToDYWeightFile = new TFile("data/ScaleFactors/Run1/ZNuNuToDYScaleFactorsRun1.root");
    TH2F *ZNuNuToDYWeightHist = (TH2F*)ZNuNuToDYWeightFile->Get("razormcDYJets");
    float maxMRZNuNuToDY = ZNuNuToDYWeightHist->GetXaxis()->GetXmax() - 1;
    float maxRsqZNuNuToDY = ZNuNuToDYWeightHist->GetYaxis()->GetXmax() - 0.01;

    vector<int> boxes;
    vector<string> boxNames;
    boxes.push_back(8);
    boxNames.push_back("MultiJet");
    boxes.push_back(10);
    boxNames.push_back("DiJet");
    int minNBTags = 1; //TODO: bin in nBTags instead of cutting
    //loop over boxes
    for(uint iBox = 0; iBox < boxes.size(); iBox++){
        gStyle->SetPaintTextFormat("1.0f");
        cout << "Analyzing " << boxNames[iBox] << " Box " << endl;
        //Step 1: Get the predictions from each MC process
        map<string, TH2F> razorHistosMC;
        map<string, TH1F> MRHistosMC;
        map<string, TH1F> RsqHistosMC;
        for(auto &tree : mctrees){
            cout << "Filling MC histograms: " << tree.first << endl;

            //set up histograms
            razorHistosMC[tree.first] = TH2F(Form("razormc%s%s", tree.first.c_str(), boxNames[iBox].c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            MRHistosMC[tree.first] = TH1F(Form("mrmc%s%s", tree.first.c_str(), boxNames[iBox].c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
            RsqHistosMC[tree.first] = TH1F(Form("rsqmc%s%s", tree.first.c_str(), boxNames[iBox].c_str()), "; Rsq", nRsqBins, RsqBinLowEdges);
            MRHistosMC[tree.first].Sumw2();
            RsqHistosMC[tree.first].Sumw2();
            razorHistosMC[tree.first].Sumw2();

            uint nEntries = tree.second->GetEntries();
            //loop over entries
            for(uint i = 0; i < nEntries; i++){
                //get entry
                tree.second->GetEntry(i); 

                //enforce correct box and number of B-tags
                if(box != boxes[iBox]) continue;
                if(nBTaggedJets < minNBTags) continue;

                //cut on MR and Rsq
                if(MR < 300 || Rsq < 0.15) continue;
                //cut on dPhiRazor
                if(fabs(dPhiRazor) > dPhiRazorCut) continue;

                float eventWeight = weight*lumiInData*1.0/lumiInMC;

                //Data/MC scale factors
                if(doSFCorrections){
                    //TTJets SF
                    if(tree.first == "TTJets"){
                        double SFTTJets = SFHistTTBar->GetBinContent(SFHistTTBar->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                        if(SFTTJets > 1e-5){
                            eventWeight *= SFTTJets;
                        }
                        else{
                            //cout << "Warning: TTJets scale factor is zero!" << endl;
                        }
                    }
                    //WJets SF
                    else if(tree.first == "WJets"){
                        double SFWJets = SFHistWJets->GetBinContent(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                        if(SFWJets > 1e-5){
                            eventWeight *= SFWJets;
                        }
                        else{
                            //cout << "Warning: WJets scale factor is zero!" << endl;
                        }
                    }
                    //DYJets SF
                    else if(tree.first == "DYJets"){
                        double SFDYJets = SFHistDYJets->GetBinContent(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                        if(SFDYJets > 1e-5){
                            eventWeight *= SFDYJets;
                        }
                        else{
                            //cout << "Warning: DYJets scale factor is zero!" << endl;
                        }
                    }
                    //ZNuNu SF
                    //TODO: combine the three predictions for ZNuNu?  currently use Gamma+Jets prediction
                    else if(tree.first == "ZJetsNuNu"){
                        double SFZJetsNuNu = SFHistZJetsNuNuFromDY->GetBinContent(SFHistZJetsNuNuFromDY->FindFixBin(min(MR, SFmaxMRZJetsNuNuFromDY), min(Rsq, SFmaxRsqZJetsNuNuFromDY)));
                        if(SFZJetsNuNu > 1e-5){
                            eventWeight *= SFZJetsNuNu;
                        }
                        else{
                            //cout << "Warning: ZJetsNuNu scale factor is zero!" << endl;
                        }

                        //scale ZNuNu so it looks like DYJets
                        double SFZNuNuToDY = ZNuNuToDYWeightHist->GetBinContent(ZNuNuToDYWeightHist->FindFixBin(min(MR, maxMRZNuNuToDY), min(Rsq, maxRsqZNuNuToDY)));
                        if(SFZNuNuToDY > 1e-5){
                            eventWeight *= SFZNuNuToDY;
                        }
                        else{
                            //cout << "Warning: ZNuNuToDY scale factor is zero!" << endl;
                        }
                    }
                }

                //TODO: keep track of uncertainties on the scale factors

                //fill each quantity
                razorHistosMC[tree.first].Fill(MR, Rsq, eventWeight);
                if(Rsq > RsqCutFor1DPlots) MRHistosMC[tree.first].Fill(MR, eventWeight);
                if(MR > MRCutFor1DPlots) RsqHistosMC[tree.first].Fill(Rsq, eventWeight);
            }
        }

        //Step 2: make data distributions
        cout << "Filling data histograms" << endl;

        //create histograms
        TH2F razorData(Form("razordata%s", boxNames[iBox].c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        TH1F MRData(Form("mrdata%s", boxNames[iBox].c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        TH1F RsqData(Form("rsqdata%s", boxNames[iBox].c_str()), "; Rsq (GeV)", nRsqBins, RsqBinLowEdges);
        razorData.Sumw2();
        MRData.Sumw2();
        RsqData.Sumw2();

        uint nEntries = datatree->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            datatree->GetEntry(i);

            //enforce correct box and number of B-tags
            if(box != boxes[iBox]) continue;
            if(nBTaggedJets < minNBTags) continue;

            //cut on MR and Rsq
            if(MR < 300 || Rsq < 0.15) continue;
            //cut on dPhiRazor
            if(fabs(dPhiRazor) > dPhiRazorCut) continue;

            float eventWeight = 1.0;

            razorData.Fill(MR, Rsq, eventWeight);
            if(Rsq > RsqCutFor1DPlots) MRData.Fill(MR, eventWeight);
            if(MR > MRCutFor1DPlots) RsqData.Fill(Rsq, eventWeight);
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
        //TODO: also draw MC before scale factors, for comparison
        TCanvas c("c", "c", 800, 600);
        c.SetLogx();

        //total MC histogram
        TH2F TotalRazorMC("TotalRazorMC", "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        //print MC histograms
        c.SetLogz();
        for(auto &hist : razorHistosMC){
            hist.second.SetTitle(Form("MC for %s, %s Box", hist.first.c_str(), boxNames[iBox].c_str()));
            hist.second.GetXaxis()->SetTitle("MR");
            hist.second.GetYaxis()->SetTitle("Rsq");
            hist.second.SetStats(0);
            hist.second.Draw("colz");
            hist.second.Draw("same,text");
            if(doSFCorrections){
                c.Print(Form("razorInclusiveMCHistogram%s%s.pdf", hist.first.c_str(), boxNames[iBox].c_str()));
                c.Print(Form("razorInclusiveMCHistogram%s%s.root", hist.first.c_str(), boxNames[iBox].c_str()));
            }
            else{
                c.Print(Form("razorInclusiveMCHistogram%s%sNoSFCorr.pdf", hist.first.c_str(), boxNames[iBox].c_str()));
                c.Print(Form("razorInclusiveMCHistogram%s%sNoSFCorr.root", hist.first.c_str(), boxNames[iBox].c_str()));
            }

            //add to total histogram
            TotalRazorMC = TotalRazorMC + hist.second;
        }
        TotalRazorMC.SetTitle(Form("Total MC, %s Box", boxNames[iBox].c_str()));
        TotalRazorMC.SetStats(0);
        TotalRazorMC.Draw("colz");
        TotalRazorMC.Draw("same,text");
        if(doSFCorrections){
            c.Print(Form("razorInclusiveMCHistogramTotal%s.pdf", boxNames[iBox].c_str()));
            c.Print(Form("razorInclusiveMCHistogramTotal%s.root", boxNames[iBox].c_str()));
        }
        else{
            c.Print(Form("razorInclusiveMCHistogramTotal%sNoSFCorr.pdf", boxNames[iBox].c_str()));
            c.Print(Form("razorInclusiveMCHistogramTotal%sNoSFCorr.root", boxNames[iBox].c_str()));
        }

        //print data histogram
        razorData.SetTitle(Form("Data, %s Box", boxNames[iBox].c_str()));
        razorData.GetXaxis()->SetTitle("MR");
        razorData.GetYaxis()->SetTitle("Rsq");
        razorData.SetStats(0);
        razorData.Draw("colz");
        razorData.Draw("same,text");
        c.Print(Form("razorInclusiveDataHistogram%s.pdf", boxNames[iBox].c_str()));
        c.Print(Form("razorInclusiveDataHistogram%s.root", boxNames[iBox].c_str()));

        //print MR and Rsq 1D histograms, comparing data to MC
        c.SetLogy();
        THStack MRTotalRazorMC("MRTotalRazorMC", Form("MR (Rsq > %.2f), %s Box", RsqCutFor1DPlots, boxNames[iBox].c_str()));
        THStack RsqTotalRazorMC("RsqTotalRazorMC", Form("Rsq (MR > %.0f), %s Box", MRCutFor1DPlots, boxNames[iBox].c_str()));

        //format MC histograms
        MRHistosMC["QCD"].SetFillColor(33);
        MRHistosMC["ZJetsNuNu"].SetFillColor(kCyan+1);
        MRHistosMC["WJets"].SetFillColor(kRed+1);
        MRHistosMC["TTJets"].SetFillColor(kGreen+3);
        MRHistosMC["DYJets"].SetFillColor(kAzure);
        MRHistosMC["SingleTop"].SetFillColor(kBlue+3);
        MRHistosMC["TTV"].SetFillColor(kSpring);
        MRHistosMC["VV"].SetFillColor(kViolet+2);
        MRHistosMC["TTTT"].SetFillColor(kRed+4);
        MRTotalRazorMC.Add(&MRHistosMC["TTTT"]);
        MRTotalRazorMC.Add(&MRHistosMC["VV"]);
        MRTotalRazorMC.Add(&MRHistosMC["TTV"]);
        MRTotalRazorMC.Add(&MRHistosMC["SingleTop"]);
        MRTotalRazorMC.Add(&MRHistosMC["DYJets"]);
        MRTotalRazorMC.Add(&MRHistosMC["TTJets"]);
        MRTotalRazorMC.Add(&MRHistosMC["WJets"]);
        MRTotalRazorMC.Add(&MRHistosMC["ZJetsNuNu"]);
        MRTotalRazorMC.Add(&MRHistosMC["QCD"]);
        MRData.SetMarkerStyle(20);
        MRData.SetMarkerSize(1);
        RsqHistosMC["QCD"].SetFillColor(33);
        RsqHistosMC["ZJetsNuNu"].SetFillColor(kCyan+1);
        RsqHistosMC["WJets"].SetFillColor(kRed+1);
        RsqHistosMC["TTJets"].SetFillColor(kGreen+3);
        RsqHistosMC["DYJets"].SetFillColor(kAzure);
        RsqHistosMC["SingleTop"].SetFillColor(kBlue+3);
        RsqHistosMC["TTV"].SetFillColor(kSpring);
        RsqHistosMC["VV"].SetFillColor(kViolet+2);
        RsqHistosMC["TTTT"].SetFillColor(kRed+4);
        RsqTotalRazorMC.Add(&RsqHistosMC["TTTT"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["VV"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["TTV"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["SingleTop"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["DYJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["TTJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["WJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["ZJetsNuNu"]);
        RsqTotalRazorMC.Add(&RsqHistosMC["QCD"]);
        RsqData.SetMarkerStyle(20);
        RsqData.SetMarkerSize(1);

        //create legend
        TLegend *RazorLegend = new TLegend(0.6, 0.6, 0.9, 0.9);
        RazorLegend->AddEntry(&MRHistosMC["WJets"], "WJets MC");
        RazorLegend->AddEntry(&MRHistosMC["DYJets"], "DYJets MC");
        RazorLegend->AddEntry(&MRHistosMC["ZJetsNuNu"], "ZJetsNuNu MC");
        RazorLegend->AddEntry(&MRHistosMC["TTJets"], "TTJets MC");
        RazorLegend->AddEntry(&MRHistosMC["SingleTop"], "Single Top MC");
        RazorLegend->AddEntry(&MRHistosMC["VV"], "VV MC");
        RazorLegend->AddEntry(&MRHistosMC["TTV"], "TTV MC");
        RazorLegend->AddEntry(&MRHistosMC["TTTT"], "TTTT MC");
        RazorLegend->AddEntry(&MRHistosMC["QCD"], "QCD MC");
        RazorLegend->AddEntry(&MRData, "2012 Data");
        if(doSFCorrections){
            DrawDataVsMCRatioPlot(&MRData, &MRTotalRazorMC, RazorLegend, "MR (GeV)", "razorInclusiveMRBackground"+boxNames[iBox], true);
            c.SetLogx(kFALSE);
            DrawDataVsMCRatioPlot(&RsqData, &RsqTotalRazorMC, RazorLegend, "Rsq (GeV)", "razorInclusiveRsqBackground"+boxNames[iBox], true);
        }
        else{
            DrawDataVsMCRatioPlot(&MRData, &MRTotalRazorMC, RazorLegend, "MR (GeV)", "razorInclusiveMRBackground"+boxNames[iBox]+"NoSFCorr", true);
            c.SetLogx(kFALSE);
            DrawDataVsMCRatioPlot(&RsqData, &RsqTotalRazorMC, RazorLegend, "Rsq (GeV)", "razorInclusiveRsqBackground"+boxNames[iBox]+"NoSFCorr", true);
        }

        gStyle->SetPaintTextFormat("1.2f");
        c.SetLogy(false);
        c.SetLogz(false);
        c.SetLogx();

        //data/MC
        TH2F DataOverMCHist = *((TH2F*)razorData.Clone(Form("DataOverMCHist%s", boxNames[iBox].c_str())));
        DataOverMCHist.Divide(&TotalRazorMC);
        DataOverMCHist.SetStats(0);
        DataOverMCHist.SetMinimum(0.1);
        DataOverMCHist.SetMaximum(3.0);
        DataOverMCHist.SetTitle(Form("Data/MC, %s Box", boxNames[iBox].c_str()));
        DataOverMCHist.Draw("colz");
        DataOverMCHist.Draw("same,text");
        if(doSFCorrections){
            c.Print(Form("razorInclusiveDataOverMC%s.pdf", boxNames[iBox].c_str()));
        }
        else{
            c.Print(Form("razorInclusiveDataOverMC%sNoSFCorr.pdf", boxNames[iBox].c_str()));
        }

        //TODO: quantify data-MC agreement in nSigmas
        delete RazorLegend;
    } //end of loop over boxes
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
    c.Print(Form("%s.root", printString.c_str()));
}
