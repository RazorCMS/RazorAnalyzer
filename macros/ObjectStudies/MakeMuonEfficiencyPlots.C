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
#include <TGraphAsymmErrors.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/MuonTree.h"

#endif


bool PassSelection( MuonTree* MuTree ) {

  bool pass = false;
  //if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && MuTree->fMuPFIso04*MuTree->fMuPt < 11) {
  //if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose) {
  //if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && MuTree->fMuPFIso04 < 0.4) {
  
  //improved isolation
  if (MuTree->fMuPt > 20) {
    if (MuTree->fMuPt > 0 
	&& MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 &&  MuTree->fMuPFIso04 < 0.4) {
      pass = true;
    }
  } else {
    if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose  && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04*MuTree->fMuPt < 10) {
      pass = true;
    }
  }

  //reco only
  //if (MuTree->fMuPt > 0) pass = true;

  return pass;

}

//=== MAIN MACRO ================================================================================================= 

void MakeMuonEfficiencyPlots(const string inputfile, int option = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Electron Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Electron Npu; Number of Events", 50, 0 , 100);

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
