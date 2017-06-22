//================================================================================================
//
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
#include "RazorAnalyzer/include/JetTree.h"

#endif


bool PassSelection( JetTree* JetTree ) {

  bool pass = false;

  //Medium WP
  if (JetTree->fJetCISV > 0.890) {
    pass = true;
  }

  return pass;
}


void plotBTaggingEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileCSVMediumBJets = new TFile("Efficiency_BJets_25ns.root","READ");
  TFile *fileCSVMediumCJets = new TFile("Efficiency_CharmJets_25ns.root","READ");
  TFile *fileCSVMediumLightJets = new TFile("Efficiency_LightJets_25ns.root","READ");
 
  TGraphAsymmErrors* effPtCSVMediumBJets = (TGraphAsymmErrors*)fileCSVMediumBJets->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtCSVMediumCJets = (TGraphAsymmErrors*)fileCSVMediumCJets->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtCSVMediumLightJets = (TGraphAsymmErrors*)fileCSVMediumLightJets->Get("Efficiency_Pt");

  //*********************************************************************
  //B-tag efficiency Vs Pt
  //*********************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.30,0.64,0.60,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtCSVMediumBJets, "B-Jets", "LP");

  effPtCSVMediumBJets->SetLineWidth(3);
  effPtCSVMediumBJets->SetLineColor(kBlack);
  effPtCSVMediumBJets->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  effPtCSVMediumBJets->GetYaxis()->SetTitle("Selection Efficiency");
  effPtCSVMediumBJets->GetYaxis()->SetTitleOffset(1.2);

  effPtCSVMediumBJets->Draw("AP");
  
  legend->Draw();

  cv->SaveAs("BTaggingEfficiencyVsPt.png");

  //*********************************************************************
  //Mis Tag efficiency Vs Pt
  //*********************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.30,0.64,0.60,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtCSVMediumCJets, "Charm Jets", "LP");
  legend->AddEntry(effPtCSVMediumLightJets, "Light Jets", "LP");

  effPtCSVMediumCJets->SetLineWidth(3);
  effPtCSVMediumCJets->SetLineColor(kBlue);
  effPtCSVMediumCJets->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  effPtCSVMediumCJets->GetYaxis()->SetTitle("Selection Efficiency");
  effPtCSVMediumCJets->GetYaxis()->SetTitleOffset(1.2);
  effPtCSVMediumCJets->GetYaxis()->SetRangeUser(0,0.3);

  effPtCSVMediumLightJets->SetLineWidth(3);
  effPtCSVMediumLightJets->SetLineColor(kRed);

  effPtCSVMediumCJets->Draw("AP");
  effPtCSVMediumLightJets->Draw("Psame");
  
  legend->Draw();

  cv->SaveAs("BTaggingMistagEfficiencyVsPt.png");

}

//=== MAIN MACRO ================================================================================================= 

void ProduceBTaggingEfficiencyPlots(const string inputfile, int option = -1, string label = "") {


  string Label = "";
  if (label != "") Label = "_" + label;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  JetTree *jetTree = new JetTree;
  jetTree->LoadTree(inputfile.c_str());
  jetTree->InitTree();
  int NCounts = 0;

  cout << "Total Entries: " << jetTree->tree_->GetEntries() << "\n";
  int nentries = jetTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < jetTree->tree_->GetEntries(); ientry++) {       	

    jetTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Selection options
    if (option == 5) {
      if (abs(jetTree->fJetPartonFlavor) != 5) continue;
    }
    if (option == 4) {
      if (abs(jetTree->fJetPartonFlavor) != 4) continue;
    }
    if (option == 0) {
      if (!(abs(jetTree->fJetPartonFlavor) == 21 || (abs(jetTree->fJetPartonFlavor) >= 1 && abs(jetTree->fJetPartonFlavor) <= 3))) continue;
    }

    //Cuts
    if (jetTree->fJetGenPt < 30) continue;
    if (abs(jetTree->fJetGenEta) > 2.4) continue;

    NCounts++;

    //**** PT ****
    histDenominatorPt->Fill(jetTree->fJetGenPt);
    
    //Numerator
    if(PassSelection(jetTree)) {
      histNumeratorPt->Fill(jetTree->fJetGenPt);        
    }
    
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");

  file->Close();
  delete file;       

}

void MakeBTaggingEff_tutorial( int Option = 0) {
 
    ProduceBTaggingEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/temp/JetNtuple_Prompt_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.Job0Of637.root", 5 , "BJets_25ns");
    ProduceBTaggingEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/temp/JetNtuple_Prompt_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.Job0Of637.root", 4 , "CharmJets_25ns");
    ProduceBTaggingEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/temp/JetNtuple_Prompt_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.Job0Of637.root", 0 , "LightJets_25ns");

    plotBTaggingEfficiency();
  
}
