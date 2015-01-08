//================================================================================================
//
// Simple Example
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
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"

#endif

int color[6] = {kRed, kGreen+2, kBlue, kMagenta, kBlack, kCyan};
//=== MAIN MACRO ================================================================================================= 


void RunSelectWJetsControlSample( vector<string> inputfiles, vector<string> processLabels, int option = -1, int finalState = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;
  double lumi = 5000;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************

  double MRBins[8] = {300, 400, 500, 750, 1000, 1500, 2000, 4000};
  double RsqBins[10] = {0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5};

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histMT;
  vector<TH1F*> histNBJetsLoose;
  vector<TH1F*> histNBJetsMedium;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 0, 3000));
    histRsq.push_back(new TH1F(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0, 2.0));
    histMT.push_back(new TH1F(Form("histMT_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 1000));
    histNBJetsLoose.push_back(new TH1F(Form("histNBJetsLoose_%s",processLabels[i].c_str()), "; Number of B-tagged Jets; Number of Events", 6, -0.5, 5.5));
    histNBJetsMedium.push_back(new TH1F(Form("histNBJetsMedium_%s",processLabels[i].c_str()), "; Number of B-tagged Jets; Number of Events", 6, -0.5, 5.5));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", 7, MRBins, 9, RsqBins));
  }
 
  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {
    ControlSampleEvents *events = new ControlSampleEvents;
    events->LoadTree(inputfiles[i].c_str());

    cout << "process: " << processLabels[i] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
    for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
      events->tree_->GetEntry(ientry);
      
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      
      //if (ientry > 1000000) break;

      //******************************
      //Selection Cuts 
      //******************************
      if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13))	     
	  ) continue;

      //lepton selection
      if (! (events->lep1.Pt() > 30 && events->lep1PassTight)
	  ) continue;

      //2nd lepton veto
      if (! (events->lep2.Pt() <= 0 
	     ||	   
	     ( (events->lep1+events->lep2).M() < 66 || (events->lep1+events->lep2).M() > 116) )
	  )
	continue;

      //MT cuts
      if (!(events->lep1MT > 30 && events->lep1MT < 100)) continue;

      //Razor signal region cuts
      if (!(events->MR > 300 && events->Rsq > 0.1)) continue;
    
      //******************************
      //Final States
      //******************************
      //e only
      if (finalState == 0 &&
	  !(abs(events->lep1Type) == 11)
	  ) continue;
      
      if (finalState == 1 &&
	  !(abs(events->lep1Type) == 13)
	  ) continue;
            
      //******************************
      //Options
      //******************************
      //Default Control Region: No CSV Medium B-tags
      if (option == 0 &&
	  events->NBJetsMedium > 0
	  ) continue;

      //1 or less CSV Loose B-tags
      if (option == 1 &&
	  events->NBJetsLoose > 1)
	continue;

      //inclusive in b-tags
      if (option == 2) {
      }


      //******************************
      //Fill histograms
      //******************************
      histMR[i]->Fill(events->MR,events->weight*lumi);
      histRsq[i]->Fill(events->Rsq,events->weight*lumi);
      histMT[i]->Fill(events->lep1MT,events->weight*lumi);      
      histNBJetsLoose[i]->Fill(events->NBJetsLoose, events->weight*lumi);    
      histNBJetsMedium[i]->Fill(events->NBJetsMedium, events->weight*lumi);    
      histMRVsRsq[i]->Fill(events->MR,events->Rsq, events->weight*lumi);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Compute Expected Statistical Uncertainty
  //==============================================================================================================
  TH2F *statUnc_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("statUnc_MRVsRsq"));
  for (int i=0; i<statUnc_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<statUnc_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
      if (statUnc_MRVsRsq->GetBinContent(i,j) > 1) {
  	statUnc_MRVsRsq->SetBinContent(i,j,1.0 / sqrt(statUnc_MRVsRsq->GetBinContent(i,j)));  	
      } else {
  	statUnc_MRVsRsq->SetBinContent(i,j,0);
      }
    }
  }

 
  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR = new THStack();
  for (int i = histMR.size()-1; i >= 0; --i) {
   histMR[i]->SetFillColor(color[i]);
    histMR[i]->SetFillStyle(1001);

    if ( histMR[i]->Integral() > 0) {
      stackMR->Add(histMR[i]);
    }
  }
  for (uint i = 0 ; i < histMR.size(); ++i) {
    legend->AddEntry(histMR[i],(processLabels[i]).c_str(), "F");
  }

  if (stackMR->GetHists()->GetEntries() > 0) {
    stackMR->Draw();
    stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_WJetsControlRegion_MR%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //Rsq
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq = new THStack();
  for (int i = histRsq.size()-1; i >= 0; --i) {
    histRsq[i]->SetFillColor(color[i]);
    histRsq[i]->SetFillStyle(1001);

    if ( histRsq[i]->Integral() > 0) {
      stackRsq->Add(histRsq[i]);
    }
  }
  for (uint i = 0 ; i < histRsq.size(); ++i) {
    legend->AddEntry(histRsq[i],(processLabels[i]).c_str(), "F");
  }

  if (stackRsq->GetHists()->GetEntries() > 0) {
    stackRsq->Draw();
    stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_WJetsControlRegion_Rsq%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //MT
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMT = new THStack();
  for (int i = histMT.size()-1; i >= 0; --i) {
    histMT[i]->SetFillColor(color[i]);
    histMT[i]->SetFillStyle(1001);

    if ( histMT[i]->Integral() > 0) {
      stackMT->Add(histMT[i]);
    }
  }
  for (uint i = 0 ; i < histMT.size(); ++i) {
    legend->AddEntry(histMT[i],(processLabels[i]).c_str(), "F");
  }

  if (stackMT->GetHists()->GetEntries() > 0) {
    stackMT->Draw();
    stackMT->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMT->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMT->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMT->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_WJetsControlRegion_MT%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //NBJetsLoose
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNBJetsLoose = new THStack();
  for (int i = histNBJetsLoose.size()-1; i >= 0; --i) {
    histNBJetsLoose[i]->SetFillColor(color[i]);
    histNBJetsLoose[i]->SetFillStyle(1001);

    if ( histNBJetsLoose[i]->Integral() > 0) {
      stackNBJetsLoose->Add(histNBJetsLoose[i]);
    }
  }
  for (uint i = 0 ; i < histNBJetsLoose.size(); ++i) {
    legend->AddEntry(histNBJetsLoose[i],(processLabels[i]).c_str(), "F");
  }

  if (stackNBJetsLoose->GetHists()->GetEntries() > 0) {
    stackNBJetsLoose->Draw();
    stackNBJetsLoose->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNBJetsLoose->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackNBJetsLoose->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNBJetsLoose->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_WJetsControlRegion_NBJetsLoose%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //NBJetsMedium
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNBJetsMedium = new THStack();
  for (int i = histNBJetsMedium.size()-1; i >= 0; --i) {
    histNBJetsMedium[i]->SetFillColor(color[i]);
    histNBJetsMedium[i]->SetFillStyle(1001);

    if ( histNBJetsMedium[i]->Integral() > 0) {
      stackNBJetsMedium->Add(histNBJetsMedium[i]);
    }
  }
  for (uint i = 0 ; i < histNBJetsMedium.size(); ++i) {
    legend->AddEntry(histNBJetsMedium[i],(processLabels[i]).c_str(), "F");
  }

  if (stackNBJetsMedium->GetHists()->GetEntries() > 0) {
    stackNBJetsMedium->Draw();
    stackNBJetsMedium->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNBJetsMedium->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackNBJetsMedium->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNBJetsMedium->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_WJetsControlRegion_NBJetsMedium%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //MR Vs Rsq
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  histMRVsRsq[0]->SetStats(false);
  histMRVsRsq[0]->Draw("colz");
  cv->SaveAs(Form("Razor_WJetsControlRegion_EventCountsWJets_MRVsRsq%s.gif",Label.c_str()));

  cv = new TCanvas("cv","cv", 800,600);
  statUnc_MRVsRsq->SetStats(false);
  statUnc_MRVsRsq->Draw("colz");
  cv->SaveAs(Form("Razor_WJetsControlRegion_StatUncWJets_MRVsRsq%s.gif",Label.c_str()));

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("WJetsControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMT[i], Form("histMT_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBJetsLoose[i], Form("histNBJetsLoose_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBJetsMedium[i], Form("histNBJetsMedium_%s",processLabels[i].c_str()), "WriteDelete");
  }
  file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");
  file->WriteTObject(stackMT, "stackMT", "WriteDelete");
  file->WriteTObject(stackNBJetsLoose, "stackNBJetsLoose", "WriteDelete");
  file->WriteTObject(stackNBJetsMedium, "stackNBJetsMedium", "WriteDelete");

  for(int i=0; i<int(histMRVsRsq.size()); i++) {
    file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  }
  file->WriteTObject(statUnc_MRVsRsq,"statUnc_MRVsRsq_WJets","WriteDelete");

  file->Close();
  delete file;       

}


void SelectWJetsControlSample( int option = -1, string label = "") {

  vector<string> inputfiles;
  vector<string> processLabels;

  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_LeptonPlusMTCutSkim_WJetsToLNu_HT100ToInf_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_LeptonPlusMTCutSkim_TTJets_25ns_weighted.root");  
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_LeptonPlusMTCutSkim_DYJetsToLL_HT100ToInf_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_LeptonPlusMTCutSkim_QCDHT100ToInf_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_LeptonPlusMTCutSkim_SingleTop_25ns_weighted.root"); 
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_Multiboson_25ns_weighted.root");

  processLabels.push_back("WJetsToLNu");
  processLabels.push_back("TTJets");
  processLabels.push_back("DYJetsToLL");
  processLabels.push_back("QCD");
  processLabels.push_back("SingleTop");
  processLabels.push_back("Multiboson");

  RunSelectWJetsControlSample(inputfiles,processLabels,0,-1,"ZeroMediumBTags_all");
  // RunSelectWJetsControlSample(inputfiles,processLabels,0,0,"ZeroMediumBTags_e");
  // RunSelectWJetsControlSample(inputfiles,processLabels,0,1,"ZeroMediumBTags_mu");

  RunSelectWJetsControlSample(inputfiles,processLabels,1,-1,"OneOrZeroLooseBTags_all");

  RunSelectWJetsControlSample(inputfiles,processLabels,2,-1,"InclusiveBTags_all");

}
