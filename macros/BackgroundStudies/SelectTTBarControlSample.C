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
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"

#endif

int color[5] = {kRed, kGreen+2, kBlue, kMagenta, kBlack};
//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarControlSample( vector<string> inputfiles, vector<string> processLabels, int option = -1, string label = "") {
  
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
  vector<TH1F*> histDileptonMass;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 0, 3000));
    histRsq.push_back(new TH1F(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0, 2.0));
    histDileptonMass.push_back(new TH1F(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 1000));
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
      if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     &&
	     (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     )
	  ) continue;

      //lepton selection
      if (! (events->lep1.Pt() > 20 && events->lep2.Pt() > 20
	     && events->lep1PassLoose && events->lep2PassLoose)
	  ) continue;

      //b-tagging
      if ( !( (events->bjet1PassLoose || events->bjet2PassLoose)
	      && events->bjet1.Pt() > 30 && events->bjet2.Pt() > 30)
	   ) continue;

      //dilepton mass cut
      if ( (events->lep1+events->lep2).M() < 20) continue;

      //Z-mass window cut
      if ( abs(events->lep1Type) == abs(events->lep2Type) 
	   && 
	   (events->lep1+events->lep2).M() > 76 && (events->lep1+events->lep2).M() < 106
	   ) continue;
  
      //Razor signal region cuts
      if (!(events->MR > 300 && events->Rsq > 0.1)) continue;
    
      //******************************
      //Options
      //******************************
      //e-mu only
      if (option == 0 &&
	  !((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13) ||
	    (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11))
	  ) continue;
      
      //ee or mumu only
      if (option == 1 && 
	  !(abs(events->lep1Type) == abs(events->lep2Type))
	  ) continue;
      
      //ee only
      if (option == 2 &&
	  !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 11)
	  ) continue;

      //mumu only
      if (option == 3 &&
	  !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 13)
	  ) continue;
      

      //******************************
      //Fill histograms
      //******************************
      histMR[i]->Fill(events->MR,events->weight*lumi);
      histRsq[i]->Fill(events->Rsq,events->weight*lumi);
      histDileptonMass[i]->Fill((events->lep1+events->lep2).M(),events->weight*lumi);      
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
    cv->SaveAs(Form("Razor_TTBarControlRegion_MR%s.gif",Label.c_str()));
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
    cv->SaveAs(Form("Razor_TTBarControlRegion_Rsq%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //DileptonMass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonMass = new THStack();
  for (int i = histDileptonMass.size()-1; i >= 0; --i) {
    histDileptonMass[i]->SetFillColor(color[i]);
    histDileptonMass[i]->SetFillStyle(1001);

    if ( histDileptonMass[i]->Integral() > 0) {
      stackDileptonMass->Add(histDileptonMass[i]);
    }
  }
  for (uint i = 0 ; i < histDileptonMass.size(); ++i) {
    legend->AddEntry(histDileptonMass[i],(processLabels[i]).c_str(), "F");
  }

  if (stackDileptonMass->GetHists()->GetEntries() > 0) {
    stackDileptonMass->Draw();
    stackDileptonMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonMass->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarControlRegion_DileptonMass%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //MR Vs Rsq
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  histMRVsRsq[0]->SetStats(false);
  histMRVsRsq[0]->Draw("colz");
  cv->SaveAs(Form("Razor_TTBarControlRegion_EventCountsTTJets_MRVsRsq%s.gif",Label.c_str()));

  cv = new TCanvas("cv","cv", 800,600);
  statUnc_MRVsRsq->SetStats(false);
  statUnc_MRVsRsq->Draw("colz");
  cv->SaveAs(Form("Razor_TTBarControlRegion_StatUncTTJets_MRVsRsq%s.gif",Label.c_str()));

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TTBarControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
  }

  for(int i=0; i<int(histMRVsRsq.size()); i++) {
    file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  }
  file->WriteTObject(statUnc_MRVsRsq,"statUnc_MRVsRsq_TTJets","WriteDelete");

  file->Close();
  delete file;       

}


void SelectTTBarControlSample( int option = -1, string label = "") {

  vector<string> inputfiles;
  vector<string> processLabels;

  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_DileptonSkim_TTJets_25ns_weighted.root");  
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_DileptonSkim_DYJetsToLL_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorControlRegions/RazorControlRegions_DileptonSkim_WJetsToLNu_25ns_weighted.root");

  processLabels.push_back("TTJets");
  processLabels.push_back("DYJetsToLL");
  processLabels.push_back("WJetsToLNu");

  RunSelectTTBarControlSample(inputfiles,processLabels,-1,"all");
  //RunSelectTTBarControlSample(inputfiles,processLabels,0,"emu");
  //RunSelectTTBarControlSample(inputfiles,processLabels,1,"eemumu");
  //RunSelectTTBarControlSample(inputfiles,processLabels,2,"ee");
  //RunSelectTTBarControlSample(inputfiles,processLabels,3,"mumu");

}
