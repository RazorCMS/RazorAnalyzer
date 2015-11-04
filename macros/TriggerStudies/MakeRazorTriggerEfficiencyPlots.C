//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
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
#include <TLegend.h>                
#include <TLatex.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"

#endif




bool PassSelection( bool *HLTDecision, int wp ) {

  bool pass = false;

 //trigger
  if (wp == 1) {
    pass = HLTDecision[94];
  }
  if (wp == 2) {
    pass = HLTDecision[95];
  }
  if (wp == 3) {
    pass = HLTDecision[96];
  }
  if (wp == 4) {
    pass = HLTDecision[97];
  }
  if (wp == 10) {
    pass = HLTDecision[94] || HLTDecision[95];
  }
  if (wp == 11) {
    pass = HLTDecision[96] || HLTDecision[97];
  }

  // if (wp == 5) {
  //   pass = HLTDecision[90];
  // }
  // if (wp == 6) {
  //   pass = HLTDecision[91];
  // }
  // if (wp == 7) {
  //   pass = HLTDecision[92];
  // }
  // if (wp == 8) {
  //   pass = HLTDecision[93];
  // }
  // if (wp == 15) {
  //   pass = HLTDecision[90] || HLTDecision[91];
  // }
  // if (wp == 16) {
  //   pass = HLTDecision[92] || HLTDecision[93];
  // }


  // //trigger
  // if (wp == 1) {
  //   pass = HLTDecision[138];
  // }
  // if (wp == 2) {
  //   pass = HLTDecision[139];
  // }
  // if (wp == 3) {
  //   pass = HLTDecision[140];
  // }
  // if (wp == 4) {
  //   pass = HLTDecision[141];
  // }
  // if (wp == 10) {
  //   pass = HLTDecision[138] || HLTDecision[139];
  // }
  // if (wp == 11) {
  //   pass = HLTDecision[140] || HLTDecision[141];
  // }

  if (wp == 5) {
    pass = HLTDecision[134];
  }
  if (wp == 6) {
    pass = HLTDecision[135];
  }
  if (wp == 7) {
    pass = HLTDecision[136];
  }
  if (wp == 8) {
    pass = HLTDecision[137];
  }
  if (wp == 15) {
    pass = HLTDecision[134] || HLTDecision[135];
  }
  if (wp == 16) {
    pass = HLTDecision[136] || HLTDecision[137];
  }

  return pass;  
}





void plotRazorTriggerEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;
  TLatex *boxLabel = 0;

  TFile *fileTrigger_RsqMR260_Rsq0p09_MR200_TTJets = new TFile("Efficiency_RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns.root","READ");
  TFile *fileTrigger_RsqMR300_Rsq0p09_MR200_TTJets = new TFile("Efficiency_RazorTrigger_RsqMR300_Rsq0p09_MR200_All_TTJets_25ns.root","READ");
  TFile *fileTrigger_RsqMR260_Rsq0p09_MR200_WJets = new TFile("Efficiency_RazorTrigger_RsqMR260_Rsq0p09_MR200_All_WJets_25ns.root","READ");
  TFile *fileTrigger_RsqMR300_Rsq0p09_MR200_WJets = new TFile("Efficiency_RazorTrigger_RsqMR300_Rsq0p09_MR200_All_WJets_25ns.root","READ");
  TFile *fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = new TFile("Efficiency_RazorTrigger_RsqMR240_Rsq0p09_MR200_All_AllLeptonData_2015D.root","READ");
  TFile *fileTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData = new TFile("Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_AllLeptonData_2015D.root","READ");

  TH2F* effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets = (TH2F*)fileTrigger_RsqMR260_Rsq0p09_MR200_TTJets->Get("Efficiency_MRRsq");
  TH2F* effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets = (TH2F*)fileTrigger_RsqMR300_Rsq0p09_MR200_TTJets->Get("Efficiency_MRRsq");
  TH2F* effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets = (TH2F*)fileTrigger_RsqMR260_Rsq0p09_MR200_WJets->Get("Efficiency_MRRsq");
  TH2F* effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets = (TH2F*)fileTrigger_RsqMR300_Rsq0p09_MR200_WJets->Get("Efficiency_MRRsq");
  TH2F* effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = (TH2F*)fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_MRRsq");
  TH2F* effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData = (TH2F*)fileTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_MRRsq");


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->SetStats(0);
  effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->GetYaxis()->SetRangeUser(0,0.8);
  
  effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"Hadronic & Loose Lepton Categories ( t#bar{t} Monte Carlo )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR260_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.pdf");


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->SetStats(0);
  effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"Hadronic & Loose Lepton Categories ( t#bar{t} Monte Carlo )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR300_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.pdf");


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->SetStats(0);
  effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"1-Lepton Categories ( W+Jets Monte Carlo )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR260_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.pdf");


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->SetStats(0);
  effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"1-Lepton Lepton Categories ( W+Jets Monte Carlo )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR300_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.pdf");


  //*******************************************************************************  
  // Data Efficiency Plots
  //*******************************************************************************


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetStats(0);
  effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetRangeUser(0,0.8);
  effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"1-Lepton Categories ( Lepton Trigger Data )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR240_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR240_Rsq0p09_MR200_AllLeptonData_2015B.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR240_Rsq0p09_MR200_AllLeptonData_2015B.pdf");


  cv = new TCanvas("cv","cv", 800,600);
  cv->SetRightMargin(0.15);
  effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->SetStats(0);
  effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
  effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->Draw("colz");

  boxLabel = new TLatex();
  boxLabel->SetNDC();
  boxLabel->SetTextSize(0.040);
  boxLabel->SetTextFont(42);
  boxLabel->SetTextColor(kBlack);
  boxLabel->DrawLatex(0.2,0.92,"1-Lepton Lepton Categories ( Lepton Trigger Data )");
  boxLabel->DrawLatex(0.4,0.85,"RsqMR270_Rsq0p09_MR200 Trigger");
  boxLabel->Draw();
 
  cv->SaveAs("RazorTriggerEfficiency_RsqMR270_Rsq0p09_MR200_AllLeptonData_2015B.gif");
  cv->SaveAs("RazorTriggerEfficiency_RsqMR270_Rsq0p09_MR200_AllLeptonData_2015B.pdf");

 
}



//=== MAIN MACRO ================================================================================================= 

void ProduceRazorTriggerEfficiencyPlots(const string inputfile, int wp, int option = -1, string label = "") {
 
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorMR = new TH1F ("histDenominatorMR",";M_{R} [GeV/c^{2}]; Number of Events", 40, 0, 1000);
  TH1F *histNumeratorMR = new TH1F ("histNumeratorMR",";M_{R} [GeV/c^{2}]; Number of Events", 40, 0 , 1000);
  TH1F *histDenominatorRsq = new TH1F ("histDenominatorRsq",";R^{2}; Number of Events", 30, 0 , 1.5);
  TH1F *histNumeratorRsq = new TH1F ("histNumeratorRsq",";R^{2}; Number of Events", 30, 0 , 1.5);

  TH2F *histDenominatorMRRsq = 0;
  TH2F *histNumeratorMRRsq = 0;
  // if (option >= 10) {
  //   histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 20, 0 , 1000, 15, 0, 1.5);
  //   histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 20, 0 , 1000, 15, 0, 1.5);
  // } else {
  histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);
  histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);
    // histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 40, 0 , 1000, 30, 0, 1.5);
    // histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 40, 0 , 1000, 30, 0, 1.5);
    
    vector<double> MRBins;
    vector<double> RsqBins;
    // MRBins.push_back(200);
    // MRBins.push_back(225); 
    // MRBins.push_back(250); 
    // MRBins.push_back(275); 
    // MRBins.push_back(300); 
    // MRBins.push_back(325); 
    // MRBins.push_back(350); 
    // MRBins.push_back(375); 
    MRBins.push_back(400); 
    // MRBins.push_back(425); 
    // MRBins.push_back(450); 
    // MRBins.push_back(500); 
    // MRBins.push_back(550); 
    // MRBins.push_back(600); 
    // MRBins.push_back(800); 
    MRBins.push_back(1000); 
    // RsqBins.push_back(0.10);
    // RsqBins.push_back(0.15);
    // RsqBins.push_back(0.20);
    RsqBins.push_back(0.25);
    // RsqBins.push_back(0.30);
    // RsqBins.push_back(0.40);
    // RsqBins.push_back(0.60);
    RsqBins.push_back(10);

    histNumeratorMRRsq = rebin(histNumeratorMRRsq, MRBins, RsqBins);
    histDenominatorMRRsq = rebin(histDenominatorMRRsq, MRBins, RsqBins);


  // }

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  TFile* inputFile = new TFile(inputfile.c_str(),"READ");
  assert(inputFile);
  TTree* tree = 0;
  tree = (TTree*)inputFile->Get("RazorInclusive");
  
  float weight = 0;
  int nvtx = 0;
  int box = -1;
  int nBTaggedJets = 0;
  int nSelectedJets = 0;
  int nJets80 = 0;
  float dPhiRazor = 0;
  float MR = 0;
  float Rsq = 0;
  float met = 0;
  bool  HLTDecision[150];
  UInt_t run = 0;
  UInt_t event = 0;

  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("box",&box);
  tree->SetBranchAddress("nVtx",&nvtx);
  tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
  tree->SetBranchAddress("nSelectedJets",&nSelectedJets);
  tree->SetBranchAddress("nJets80",&nJets80);
  tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
  tree->SetBranchAddress("MR",&MR);
  tree->SetBranchAddress("Rsq",&Rsq);
  tree->SetBranchAddress("met",&met);
  tree->SetBranchAddress("HLTDecision",&HLTDecision);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);

  cout << "Total Entries: " << tree->GetEntries() << "\n";

  //for duplicate event checking
  map<pair<uint,uint>, bool > processedRunEvents;

  for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {       	
    tree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;


    //Cuts
    if (!(nJets80 >= 2)) continue;

    //Select Hadronic boxes
    if (option == 0) {
      //hadronic boxes
      //if (!(box == 9 || box == 10 || box == 11 || box == 12 || box == 13 || box == 14)) continue;

      //leptonic boxes
      // if (!(box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box ==  8)) continue;      
      if (!(box == 6 || box == 7 || box ==  8)) continue;      
    }

    //Select 1L and 2L boxes
    if (option == 1 || option == 11) {
      //1L and 2L boxes
      bool passedDileptonTrigger = bool( HLTDecision[41] || HLTDecision[43] 
				    || HLTDecision[30] || HLTDecision[31] 
				    || HLTDecision[47] || HLTDecision[48] || HLTDecision[49] || HLTDecision[50] );
      bool passedSingleLeptonTrigger = bool(HLTDecision[2] || HLTDecision[7] || HLTDecision[12] || HLTDecision[11] || HLTDecision[15]
				       || HLTDecision[22] || HLTDecision[23] || HLTDecision[24] || HLTDecision[25] || 
				       HLTDecision[26] || HLTDecision[27] ||
				       HLTDecision[28] || HLTDecision[29]);      
      if (!(passedDileptonTrigger || passedSingleLeptonTrigger)) continue;
      // if (!(box == 0 || box == 1 || box == 2 || box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box ==  8)) continue;
      
      //only 1L boxes
      //if (!(box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box ==  8)) continue;

      if (!(box == 6 || box == 7 || box ==  8)) continue;


    }


    //Use events triggered by HT
    if (option == 2) {
      if (!(HLTDecision[88] || HLTDecision[89])) continue;      
    }

    //Use events triggered by HT
    if (option == 3) {
      if (!(HLTDecision[98] || HLTDecision[99] || HLTDecision[100])) continue;      
    }

    //Use events triggered by MET
    if (option == 4) {
      if (!(HLTDecision[85] )) continue;      
    }


    //**** MR - Rsq ****
    histDenominatorMRRsq->Fill(MR,Rsq);
    if(PassSelection(HLTDecision,wp)) {
      //cout << MR << " " << Rsq << " " << weight << "\n";
      histNumeratorMRRsq->Fill(MR,Rsq);
    }

    //Cuts
    if (!(MR > 400 && Rsq > 0.25)) continue;

    //Remove double counted events        
    if(!(processedRunEvents.find(make_pair(run, event)) == processedRunEvents.end())) {
      continue;
    } else {
      processedRunEvents[make_pair(run, event)] = true;
    }


    //**** MR ****
    if (Rsq>0.25) { 
      histDenominatorMR->Fill(MR);
      //Numerator
      if(PassSelection(HLTDecision,wp)) {
	histNumeratorMR->Fill(MR);
      }
    }

    //**** Rsq ****
    if (MR>400) {
      histDenominatorRsq->Fill(Rsq);      
      //Numerator
      if(PassSelection(HLTDecision,wp)) {
	histNumeratorRsq->Fill(Rsq);
      }
    }
 
     
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_MR = createEfficiencyGraph(histNumeratorMR, histDenominatorMR, "Efficiency_MR" , vector<double>() ,  -99, -99, 0.9, 1);
  TGraphAsymmErrors *efficiency_Rsq = createEfficiencyGraph(histNumeratorRsq, histDenominatorRsq, "Efficiency_Rsq" , vector<double>() ,  -99, -99, 0.9, 1);
  TH2F *efficiency_MRRsq = createEfficiencyHist2D(histNumeratorMRRsq, histDenominatorMRRsq, "Efficiency_MRRsq" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_MR->Draw("AP");
  efficiency_MR->SetTitle("");
  efficiency_MR->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_MR->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  efficiency_MR->GetYaxis()->SetTitle("Efficiency");
  efficiency_MR->GetYaxis()->SetTitleOffset(1.2);
  efficiency_MR->SetLineWidth(3);  
  cv->SaveAs(("Efficiency"+Label+"_MR.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_Rsq->Draw("AP");
  efficiency_Rsq->SetTitle("");
  efficiency_Rsq->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_Rsq->GetXaxis()->SetTitle("R^{2}");
  efficiency_Rsq->GetYaxis()->SetTitle("Efficiency");
  efficiency_Rsq->GetYaxis()->SetTitleOffset(1.2);
  efficiency_Rsq->SetLineWidth(3);  
  cv->SaveAs(("Efficiency"+Label+"_Rsq.gif").c_str());



  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_MR, "Efficiency_MR", "WriteDelete");
  file->WriteTObject(efficiency_Rsq, "Efficiency_Rsq", "WriteDelete");
  file->WriteTObject(efficiency_MRRsq, "Efficiency_MRRsq", "WriteDelete");

  file->Close();
  delete file;       

}

void MakeRazorTriggerEfficiencyPlots( int option = 0) {

  if (option == 1) {  
    //***************************************
    // TTbar MC : Use Hadronic Boxes
    //***************************************
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 1, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_TTJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_4jet_TTJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 3, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_TTJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 4, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_4jet_TTJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 11, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_TTJets_25ns");   
 
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns");   

    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns"); 
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns"); 

    //***************************************
    // W+Jets MC: Use Single Lepton Boxes
    //***************************************
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 1, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_WJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 2, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_4jet_WJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 3, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_WJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 4, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_4jet_WJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 10, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_WJets_25ns");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 11, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_WJets_25ns");   
  
  }

  if (option == 2) {
 
    //***************************************
    // Try all different lepton datasets
    //***************************************
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_MuonEG_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_MuonEG_2015B");   
     // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_SingleElectron_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_SingleElectron_2015B");   
     // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_SingleMuon_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_SingleMuon_2015B");   
     // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_DoubleMuon_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_DoubleMuon_2015B");   
     // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_DoubleEG_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_DoubleEG_2015B");   
     // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_Tau_Run2015B_GoodLumi.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_Tau_2015B");   

  

    //***************************************
    // Use combined data
    //***************************************
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_AllLeptonData_Run2015B_GoodLumi.root", 5, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData_2015B");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_AllLeptonData_Run2015B_GoodLumi.root", 6, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_4jet_AllLeptonData_2015B");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_AllLeptonData_Run2015B_GoodLumi.root", 7, 11, "RazorTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData_2015B");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_AllLeptonData_Run2015B_GoodLumi.root", 8, 11, "RazorTrigger_RsqMR270_Rsq0p09_MR200_4jet_AllLeptonData_2015B");   
    //ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/oldData/RazorInclusive_CombinedLeptonic_Run2015D_GoodLumiGolden.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_AllLeptonData_2015D");   
    // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/oldData/RazorInclusive_CombinedLeptonic_Run2015D_GoodLumiGolden.root", 16, 11, "RazorTrigger_RsqMR270_Rsq0p09_MR200_All_AllLeptonData_2015D");   

    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data/RazorInclusive_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_AllLeptonData_2015D");   


  }

  if (option == 3) {
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_JetHT_Run2015B_GoodLumi.root", 15, 2, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_HT800Data_2015B");   
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_JetHT_Run2015B_GoodLumi.root", 16, 2, "RazorTrigger_RsqMR270_Rsq0p09_MR200_All_HT800Data_2015B");   
  }

  if (option == 4) {
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_HTMHT_Run2015B_GoodLumi.root", 15, 3, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_RsqData_2015B");   
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_HTMHT_Run2015B_GoodLumi.root", 16, 3, "RazorTrigger_RsqMR270_Rsq0p09_MR200_All_RsqData_2015B");   
  }

  if (option == 5) {
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_MET_Run2015B_GoodLumi.root", 15, 4, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_METData_2015B");   
    ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Run2015B/RazorInclusive_MET_Run2015B_GoodLumi.root", 16, 4, "RazorTrigger_RsqMR270_Rsq0p09_MR200_All_METData_2015B");   
  }

  plotRazorTriggerEfficiency();

}
