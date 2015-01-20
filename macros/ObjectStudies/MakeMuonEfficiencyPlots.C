//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeMuonEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/MuonNtuple/MuonNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Muon")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TLegend.h>                
#include <TGraphAsymmErrors.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/MuonTree.h"

#endif


bool PassSelection( MuonTree* MuTree ) {

  bool pass = false;

  //if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && MuTree->fMuPFIso04*MuTree->fMuPt < 11) {
  //if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose) {
  // if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04 < 0.4) {
  //   pass = true;
  // }
  
  //improved isolation

  // //**********************************
  // //Tight Selection
  // //**********************************
  // if (MuTree->fMuPt > 0 && MuTree->fMuIsTight && fabs(MuTree->fMuIP3dSig)<4 && fabs(MuTree->fMuD0) < 0.2 && MuTree->fMuPFIso04 < 0.12) {
  //   pass = true;
  // }
 
  // //**********************************
  // //Loose Selection
  // //**********************************
  // if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04 < 0.2) {
  //   pass = true;
  // }
  
  // //**********************************
  // //Veto Selection
  // //**********************************
  // if (MuTree->fMuPt > 20) {
  //   if (MuTree->fMuPt > 0 
  // 	&& MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 &&  MuTree->fMuPFIso04 < 0.4) {
  //     pass = true;
  //   }
  // } else {
  //   if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04*MuTree->fMuPt < 10) {
  //     pass = true;
  //   }
  // }

  //reco only
  //if (MuTree->fMuPt > 0) pass = true;

  return pass;

}


void plotMuonEfficiency() {

  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileRelIso = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_NumeratorRelIso0p4_DenominatorLooseIDAndIPCut.root","READ");
  TFile *fileImprovedIso = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_NumeratorImprovedIso_DenominatorLooseIDAndIPCut.root","READ");
  TFile *fileVeto = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_VetoSelection.root","READ");
  TFile *fileLoose = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_LooseSelection.root","READ");
  TFile *fileTight = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_TightSelection.root","READ");
  TFile *fileFakesVeto = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_FakeMuon_VetoSelection.root","READ");
  TFile *fileFakesLoose = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_FakeMuon_LooseSelection.root","READ");
  TFile *fileFakesTight = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_FakeMuon_TightSelection.root","READ");

  TGraphAsymmErrors* effPtRelIso = (TGraphAsymmErrors*)fileRelIso->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtImprovedIso = (TGraphAsymmErrors*)fileImprovedIso->Get("Efficiency_Pt");

  TGraphAsymmErrors* effPtVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effEtaVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effNpvVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_NPV");

  TGraphAsymmErrors* effFakePtVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakeEtaVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeNpvVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_NPV");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtRelIso, "Relative Iso < 0.4", "LP");
  legend->AddEntry(effPtImprovedIso, "Improved Veto Isolation", "LP");

  effPtRelIso->SetLineWidth(3);
  effPtRelIso->SetLineColor(kBlack);
  effPtRelIso->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effPtRelIso->GetYaxis()->SetTitle("Isolation Efficiency");
  effPtRelIso->GetYaxis()->SetTitleOffset(1.2);

  effPtRelIso->SetLineWidth(3);
  effPtRelIso->SetLineColor(kRed);
  effPtImprovedIso->SetLineWidth(3);
  effPtImprovedIso->SetLineColor(kBlue);

  effPtRelIso->Draw("AP");
  effPtImprovedIso->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonIsolationEfficiencyVsPt.gif");
  cv->SaveAs("MuonIsolationEfficiencyVsPt.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtVeto, "Veto", "LP");
  legend->AddEntry(effPtLoose, "Loose", "LP");
  legend->AddEntry(effPtTight, "Tight", "LP");

  effPtVeto->SetLineWidth(3);
  effPtVeto->SetLineColor(kBlack);
  effPtVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effPtVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effPtVeto->GetYaxis()->SetTitleOffset(1.2);

  effPtLoose->SetLineWidth(3);
  effPtLoose->SetLineColor(kBlue);
  effPtTight->SetLineWidth(3);
  effPtTight->SetLineColor(kRed);

  effPtVeto->Draw("AP");
  effPtLoose->Draw("Psame");
  effPtTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsPt.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsPt.pdf");


 
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaVeto, "Veto", "LP");
  legend->AddEntry(effEtaLoose, "Loose", "LP");
  legend->AddEntry(effEtaTight, "Tight", "LP");

  effEtaVeto->SetLineWidth(3);
  effEtaVeto->SetLineColor(kBlack);
  effEtaVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaVeto->GetYaxis()->SetTitleOffset(1.2);

  effEtaLoose->SetLineWidth(3);
  effEtaLoose->SetLineColor(kBlue);
  effEtaTight->SetLineWidth(3);
  effEtaTight->SetLineColor(kRed);

  effEtaVeto->Draw("AP");
  effEtaLoose->Draw("Psame");
  effEtaTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsEta.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsEta.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNpvVeto, "Veto", "LP");
  legend->AddEntry(effNpvLoose, "Loose", "LP");
  legend->AddEntry(effNpvTight, "Tight", "LP");

  effNpvVeto->SetLineWidth(3);
  effNpvVeto->SetLineColor(kBlack);
  effNpvVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effNpvVeto->GetYaxis()->SetTitleOffset(1.2);
  effNpvVeto->GetXaxis()->SetRangeUser(5,35);
  effNpvVeto->GetYaxis()->SetRangeUser(0.5,1.0);

  effNpvLoose->SetLineWidth(3);
  effNpvLoose->SetLineColor(kBlue);
  effNpvTight->SetLineWidth(3);
  effNpvTight->SetLineColor(kRed);

  effNpvVeto->Draw("AP");
  effNpvLoose->Draw("Psame");
  effNpvTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsNpv.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsNpv.pdf");






  //***************************************************************
  //Fake Muons : Efficiency Vs Pt
  //***************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakePtVeto, "Veto", "LP");
  legend->AddEntry(effFakePtLoose, "Loose", "LP");
  legend->AddEntry(effFakePtTight, "Tight", "LP");

  effFakePtVeto->SetLineWidth(3);
  effFakePtVeto->SetLineColor(kBlack);
  effFakePtVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effFakePtVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakePtVeto->GetYaxis()->SetTitleOffset(1.35);

  effFakePtLoose->SetLineWidth(3);
  effFakePtLoose->SetLineColor(kBlue);
  effFakePtTight->SetLineWidth(3);
  effFakePtTight->SetLineColor(kRed);

  effFakePtVeto->Draw("AP");
  effFakePtLoose->Draw("Psame");
  effFakePtTight->Draw("Psame");
  
  effFakePtVeto->GetYaxis()->SetRangeUser(0,0.05);

  legend->Draw();  
  //cv->SetLogy();
  cv->SaveAs("MuonSelectionFakeRateVsPt.gif");
  cv->SaveAs("MuonSelectionFakeRateVsPt.pdf");

  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.7,0.90,0.85);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeEtaVeto, "Veto", "LP");
  legend->AddEntry(effFakeEtaLoose, "Loose", "LP");
  legend->AddEntry(effFakeEtaTight, "Tight", "LP");

  effFakeEtaVeto->SetLineWidth(3);
  effFakeEtaVeto->SetLineColor(kBlack);
  effFakeEtaVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effFakeEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeEtaVeto->GetYaxis()->SetTitleOffset(1.35);

  effFakeEtaLoose->SetLineWidth(3);
  effFakeEtaLoose->SetLineColor(kBlue);
  effFakeEtaTight->SetLineWidth(3);
  effFakeEtaTight->SetLineColor(kRed);

  effFakeEtaVeto->GetYaxis()->SetRangeUser(0,0.025);

  effFakeEtaVeto->Draw("AP");
  effFakeEtaLoose->Draw("Psame");
  effFakeEtaTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionFakeRateVsEta.gif");
  cv->SaveAs("MuonSelectionFakeRateVsEta.pdf");




 cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.65,0.90,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeNpvVeto, "Veto", "LP");
  legend->AddEntry(effFakeNpvLoose, "Loose", "LP");
  legend->AddEntry(effFakeNpvTight, "Tight", "LP");

  effFakeNpvVeto->SetLineWidth(3);
  effFakeNpvVeto->SetLineColor(kBlack);
  effFakeNpvVeto->GetXaxis()->SetTitle("Muon p_{T} [Gev/C]");
  effFakeNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeNpvVeto->GetYaxis()->SetTitleOffset(1.35);
  effFakeNpvVeto->GetXaxis()->SetRangeUser(5,35);
  effFakeNpvVeto->GetYaxis()->SetRangeUser(0.0,0.04);

  effFakeNpvLoose->SetLineWidth(3);
  effFakeNpvLoose->SetLineColor(kBlue);
  effFakeNpvTight->SetLineWidth(3);
  effFakeNpvTight->SetLineColor(kRed);

  effFakeNpvVeto->Draw("AP");
  effFakeNpvLoose->Draw("Psame");
  effFakeNpvTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionFakeRateVsNpv.gif");
  cv->SaveAs("MuonSelectionFakeRateVsNpv.pdf");






}




//=== MAIN MACRO ================================================================================================= 

void MakeMuonEfficiencyPlots(const string inputfile, int option = -1, string label = "") {

  // plotMuonEfficiency();
  // return;

  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Muon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Muon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Muon Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Muon Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Muon Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Muon Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Muon Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Muon Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Muon Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Muon Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Muon Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Muon Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  MuonTree *MuTree = new MuonTree;
  MuTree->LoadTree(inputfile.c_str());
  MuTree->InitTree(MuonTree::kMuTreeLight);

  cout << "Total Entries: " << MuTree->tree_->GetEntries() << "\n";
  int nentries = MuTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < MuTree->tree_->GetEntries(); ientry++) {       	
    MuTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (MuTree->fMuGenPt < 5) continue;
    if (abs(MuTree->fMuGenEta) > 2.4) continue;

    if (!(MuTree->fMuPt > 5)) continue;

    //Pass ID cuts   
    // if (!(MuTree->fMuPt > 0 && MuTree->fMuIsLoose  && fabs(MuTree->fMuIP3dSig)<4)) continue;


    //**** PT - ETA ****
    histDenominatorPtEta->Fill(MuTree->fMuGenPt,MuTree->fMuGenEta);
    if(PassSelection(MuTree)) {
      histNumeratorPtEta->Fill(MuTree->fMuGenPt,MuTree->fMuGenEta);
    }


    //**** PT ****
      histDenominatorPt->Fill(MuTree->fMuGenPt);

      //Numerator
      if(PassSelection(MuTree)) {
        histNumeratorPt->Fill(MuTree->fMuGenPt);        
      }


    //**** Eta ****
    if (fabs(MuTree->fMuGenPt) > 30) {
      histDenominatorEta->Fill(MuTree->fMuGenEta);

      //Numerator
      if(PassSelection(MuTree)) {
        histNumeratorEta->Fill(MuTree->fMuGenEta);        
      }

    }

    //**** Phi ****
    if (fabs(MuTree->fMuGenEta) < 2.4) {
      histDenominatorPhi->Fill(MuTree->fMuGenPhi);

      //Numerator
      if(PassSelection(MuTree)) {
        histNumeratorPhi->Fill(MuTree->fMuGenPhi);        
      }

    }

    //**** Rho ****
    if (fabs(MuTree->fMuGenEta) < 2.4) {
      histDenominatorRho->Fill(MuTree->fRho);

      //Numerator
      if(PassSelection(MuTree)) {
        histNumeratorRho->Fill(MuTree->fRho);        
      }

    }
    //**** Npv ****
    if (fabs(MuTree->fMuGenEta) < 2.4) {
      histDenominatorNpv->Fill(MuTree->fNVertices);

      //Numerator
      if(PassSelection(MuTree)) {
        histNumeratorNpv->Fill(MuTree->fNVertices);        
      }

    }

    // //**** Npu ****
    // if (fabs(MuTree->fMuGenEta) < 2.4) {
    //   histDenominatorNpu->Fill(MuTree->);

    //   //Numerator
    //   if(PassSelection(MuTree)) {
    //     histNumeratorNpu->Fill(MuTree->);        
    //   }

    // }


  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "Efficiency_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "Efficiency_Npu" , vector<double>() ,  -99, -99, 0, 1);  
  TH2F *efficiency_pteta = createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Pt.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Eta.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Phi.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Rho.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npv.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npu.gif").c_str());


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");
  file->WriteTObject(efficiency_eta, "Efficiency_Eta", "WriteDelete");
  file->WriteTObject(efficiency_phi, "Efficiency_Phi", "WriteDelete");
  file->WriteTObject(efficiency_rho, "Efficiency_Rho", "WriteDelete");
  file->WriteTObject(efficiency_npv, "Efficiency_NPV", "WriteDelete");
  file->WriteTObject(efficiency_npu, "Efficiency_NPU", "WriteDelete");
  file->WriteTObject(efficiency_pteta, "Efficiency_PtEta", "WriteDelete");

  file->Close();
  delete file;       

}
