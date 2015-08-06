//Misc functions used for razor analysis macros

#ifndef MacroHelper_H
#define MacroHelper_H

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"

#include "include/RazorAnalyzer.h"

//returns the Data/MC scale factor, and its error, at the desired value of MR and Rsq
pair<double, double> getDataMCSFAndError(TH2* sfHist, float MR, float Rsq){
    float maxMR = sfHist->GetXaxis()->GetXmax() - 1;
    float maxRsq = sfHist->GetYaxis()->GetXmax() - 0.01;
    double sf = sfHist->GetBinContent(sfHist->FindFixBin(min(MR, maxMR), min(Rsq, maxRsq)));
    double sfErr = sfHist->GetBinError(sfHist->FindFixBin(min(MR, maxMR), min(Rsq, maxRsq)));
    if(sf < 1e5){
        return std::make_pair(sf, sfErr);
    }
    else{
        std::cout << "Warning: Data/MC scale factor is Inf!  Returning 0" << endl;
        return std::make_pair(0.0, 0.0);
    }
}

//plotting macro for data and MC, including ratio plot
//parameters: 
//dataHist: TH1F with data points (to draw as black dots)
//mcStack: a THStack containing all MC histograms (to draw as filled bar histograms)
//leg: TLegend containing all of the entries in the desired order
//printString: name of output file
//fitHist: optional fit result histogram (to draw as a solid line)
void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX, TH1F *fitHist=0){
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
    mcStack->SetMaximum(max(mcStack->GetMaximum(), dataHist->GetMaximum()));
    mcStack->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    mcStack->GetYaxis()->SetLabelSize(0.03);
    mcStack->GetYaxis()->SetTitleOffset(0.45);
    mcStack->GetYaxis()->SetTitleSize(0.05);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(1);
    dataHist->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    dataHist->Draw("pesame");
    if(fitHist){
        fitHist->SetLineWidth(2);
        fitHist->Draw("lsame");
    }
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
    c.Print(Form("%sLog.gif", printString.c_str()));
    //c.Print(Form("%s.root", printString.c_str()));
    pad1.SetLogy(kFALSE);
    pad1.Modified();
    gPad->Update();
    c.Print(Form("%sLinear.gif", printString.c_str()));
}

//check if the given box is a muon box
bool isSingleMuonBox(RazorAnalyzer::RazorBox box){
    if(box == RazorAnalyzer::MuSixJet || box == RazorAnalyzer::MuFourJet || box == RazorAnalyzer::MuJet || box == RazorAnalyzer::MuMultiJet) return true;
    return false;
}

//check if the given box is a electron box
bool isSingleElectronBox(RazorAnalyzer::RazorBox box){
    if(box == RazorAnalyzer::EleSixJet || box == RazorAnalyzer::EleFourJet || box == RazorAnalyzer::EleJet || box == RazorAnalyzer::EleMultiJet) return true;
    return false;
}

//check if the given box is a lepton box (excluding LooseLepton boxes)
bool isLeptonicBox(RazorAnalyzer::RazorBox box){
    if(box == RazorAnalyzer::MuJet
            || box == RazorAnalyzer::MuFourJet
            || box == RazorAnalyzer::MuSixJet
            || box == RazorAnalyzer::MuMultiJet
            || box == RazorAnalyzer::EleJet
            || box == RazorAnalyzer::EleFourJet
            || box == RazorAnalyzer::EleSixJet
            || box == RazorAnalyzer::EleMultiJet
            || box == RazorAnalyzer::MuMu
            || box == RazorAnalyzer::MuEle
            || box == RazorAnalyzer::EleEle) return true;
    return false;
}

#endif
