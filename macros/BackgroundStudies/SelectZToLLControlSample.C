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
#include <TH1D.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 
#include <TRandom3.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"

#endif

void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv", 800,700);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.50,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1D *histDataOverMC = (TH1D*)hist[0]->Clone("histDataOverMC");

  if (hasData) {
    for (int i = hist.size()-1; i >= 1; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  } else {
    for (int i = hist.size()-1; i >= 0; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  }

  for (uint i = 0 ; i < hist.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stack->GetHists()->GetEntries() > 0) {
    stack->Draw("hist");
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetLineWidth(2);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("e1same");
    }
    legend->Draw();
  }
  cv->cd();
  cv->Update();


  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->Draw();
  pad2->cd();
    
  for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
    double data = 0;
    if (hasData) {
      data = hist[0]->GetBinContent(b);
    }
    double MC = 0;
    double MCErrSqr = 0;
    if (hasData) {
      for (uint i = 1 ; i < hist.size(); ++i) {
	MC += hist[i]->GetBinContent(b);
	MCErrSqr += pow(hist[i]->GetBinError(b),2);
      }
    } else {
      MC = 1;
    }
      
    if (MC > 0) {
      histDataOverMC->SetBinContent(b, data / MC);
      histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
    } else {
      histDataOverMC->SetBinContent(b, 0);
      histDataOverMC->SetBinError(b, 0);
    }
    //cout << "bin " << b << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
  }

  histDataOverMC->GetYaxis()->SetTitle("Data/MC");
  histDataOverMC->GetYaxis()->SetNdivisions(306);
  histDataOverMC->GetYaxis()->SetTitleSize(0.10);
  histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
  histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
  histDataOverMC->GetYaxis()->SetLabelSize(0.10);
  histDataOverMC->GetXaxis()->SetLabelSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleSize(0.15);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("e1");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s.gif",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s_Logy.gif",varName.c_str(),label.c_str()));

}


// void CompareMetResolution () {

//   //Compare MET resolution for different MET's
//   TFile *fileNoMetCorr = new TFile("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim_PFMetNoCorrections/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root","READ");
//   TFile *fileType1Corr40GeVJetsNoPUJets = new TFile("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim_Type1Met40GeVNoPileupJets/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root","READ");
//   TFile *fileType1Corr20GeVJetsNoPUJets = new TFile("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim_Type1Met20GeVNoPileupJets/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root","READ");
//   TFile *fileType1Corr20GeVJetsWithPUJets = new TFile("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim_Type1Met20GeVInclPileupJets/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root","READ");

//   TTree *treeNoMetCorr = (TTree*)fileNoMetCorr->Get("ControlSampleEvent");
//   TTree *treeType1Corr40GeVJetsNoPUJets = (TTree*)fileType1Corr40GeVJetsNoPUJets->Get("ControlSampleEvent");
//   TTree *treeType1Corr20GeVJetsNoPUJets = (TTree*)fileType1Corr20GeVJetsNoPUJets->Get("ControlSampleEvent");
//   TTree *treeType1Corr20GeVJetsWithPUJets = (TTree*)fileType1Corr20GeVJetsWithPUJets->Get("ControlSampleEvent");

//   treeNoMetCorr->Draw("MET>>hNoMetCorr(100,0,300)","lep1.Pt() > 35 && lep2.Pt() > 35");
//   treeType1Corr40GeVJetsNoPUJets->Draw("MET>>hType1Corr40GeVJetsNoPUJets(100,0,300)","lep1.Pt() > 35 && lep2.Pt() > 35");
//   treeType1Corr20GeVJetsNoPUJets->Draw("MET>>hType1Corr20GeVJetsNoPUJets(100,0,300)","lep1.Pt() > 35 && lep2.Pt() > 35");
//   treeType1Corr20GeVJetsWithPUJets->Draw("MET>>hType1Corr20GeVJetsWithPUJets(100,0,300)","lep1.Pt() > 35 && lep2.Pt() > 35");

//   TLegend* l = new TLegend(0.5,0.7,0.7,0.85);
//   l->AddEntry(hNoMetCorr,"No Corrections","L");
//   l->AddEntry(hType1Corr40GeVJetsNoPUJets,"Type 1 Corr ( 40 GeV Jets)","L");
//   l->AddEntry(hType1Corr20GeVJetsNoPUJets,"Type 1 Corr ( 20 GeV Jets)","L");
//   l->AddEntry(hType1Corr20GeVJetsWithPUJets,"Type 1 Corr ( 20 GeV Jets incl Pileup)","L");

//   TCanvas *cv = new TCanvas("cv","cv", 800,600);
//   hNoMetCorr->SetLineColor(kBlue);
//   hType1Corr40GeVJetsNoPUJets->SetLineColor(kRed);
//   hType1Corr20GeVJetsNoPUJets->SetLineColor(kGreen+2);
//   hType1Corr20GeVJetsWithPUJets->SetLineColor(kViolet);
//   hNoMetCorr->Draw();
//   hType1Corr40GeVJetsNoPUJets->Draw("same");
//   hType1Corr20GeVJetsNoPUJets->Draw("same");
//   hType1Corr20GeVJetsWithPUJets->Draw("same");
  
// }

//=== MAIN MACRO ================================================================================================= 


void RunSelectZToLLControlSample( string datafile, vector<string> bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;
  TRandom3 *random = new TRandom3;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/Run1/ElectronSelection_Run2012ReReco_53X.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  assert(eleEffSFHist);

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 5;
  const int NRsqBins = 4;
  double MRBins[NMRBins] = {300, 350, 400, 450, 550};
  double RsqBins[NRsqBins] = {0.05, 0.075, 0.10, 1.5};

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;
  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back("Data");
    color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }

  vector<TH1D*> histNVtx;
  vector<TH1D*> histMR;
  vector<TH1D*> histRsq;
  vector<TH1D*> histDileptonMass;
  vector<TH1D*> histMET;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histNVtx.push_back(new TH1D(Form("histNVtx_%s",processLabels[i].c_str()), "; Number of Primary Vertices; R^{2}; Number of Events", 100, -0.5, 99.5));
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 400, 2500));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0.05, 1.0));
    histDileptonMass.push_back(new TH1D(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 200));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV/c] ; Number of Events", 100, 0, 1000));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 15, -0.5, 14.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 10, -0.5,9.5));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histNVtx[i]->Sumw2();
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histDileptonMass[i]->Sumw2();
    histMET[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNJets80[i]->Sumw2();
    histNBtags[i]->Sumw2();   
    histMRVsRsq[i]->Sumw2();
  }
 
  double dataYield = 0;
  double MCYield = 0;

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

      double puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(events->NPU_0));
      double weight = lumi * events->weight * puWeight;

      //******************************
      //Trigger Selection
      //******************************
      bool passTrigger = false;

      //DiMuon Triggers: Mu17Mu8 , Mu17TkMu8
      if (events->HLTDecision[3] ==true || events->HLTDecision[4] ==true) passTrigger = true;

      //DiElectron Triggers:
      if (events->HLTDecision[12] ==true) passTrigger = true;

      if (!passTrigger) continue;

      //******************************
      //Selection Cuts 
      //******************************
      if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     &&
	     (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     )
	  ) continue;

      //lepton selection
      if (! (events->lep1.Pt() > 25 && events->lep2.Pt() > 25
	     && events->lep1PassLoose && events->lep2PassLoose)
	  ) continue;

      //dilepton mass cut
      if ( (events->lep1+events->lep2).M() < 20) continue;

      // BTag Veto
      if ( !( events->NBJetsMedium == 0)) continue;

      //Razor signal region cuts
      if (option == "TwoJet80") {
	if (!(events->NJets80 >= 2 )) continue;
      }
      
      if (option == "MR300Rsq0p05") {
	if (!(events->MR > 300 && events->Rsq > 0.05 )) continue;
      }
      
    
      //******************************
      //Options
      //******************************
      //e-mu only
      if (channelOption == 0 &&
	  !((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13) ||
	    (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11))
	  ) continue;
      
      //ee or mumu only
      if (channelOption == 1 && 
	  !(abs(events->lep1Type) == abs(events->lep2Type))
	  ) continue;
      
      //ee only
      if (channelOption == 2 &&
	  !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 11)
	  ) continue;

      //mumu only
      if (channelOption == 3 &&
	  !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 13)
	  ) continue;
      

      //******************************
      //Apply Scale Factors
      //******************************
      if (!(hasData && i==0)) {
      	double triggerEffScaleFactor = 1.0;
      	if ( (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 13)) triggerEffScaleFactor = 0.970; 
	
	double leptonEffScaleFactor = 1.0;
	if (abs(events->lep1Type) == 11) {
	  leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fabs(events->lep1.Eta())) , 
							       eleEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),10.01)));	 
	}
	if (abs(events->lep2Type) == 11) {
	  leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fabs(events->lep2.Eta())) , 
							       eleEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),10.01)));
	}
	
	//extra scale factor for MC to agree with data yield
	double normalizationScaleFactor = 1.0;
	if ( (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 13)) normalizationScaleFactor = 0.970;
	if ( (abs(events->lep1Type) == 11 && abs(events->lep2Type) == 11)) normalizationScaleFactor = 0.970;
	
      	weight *= leptonEffScaleFactor;
      	weight *= triggerEffScaleFactor;
	weight *= normalizationScaleFactor;
      }


      //******************************
      //Fill histograms
      //******************************
      if (hasData && i==0) {
	if ((events->lep1+events->lep2).M() > 80 && (events->lep1+events->lep2).M() < 100) {
	  histMR[i]->Fill(events->MR);
	  histRsq[i]->Fill(events->Rsq);
	  histDileptonMass[i]->Fill((events->lep1+events->lep2).M());      
	  histMET[i]->Fill(events->MET);  
	  histNJets40[i]->Fill(events->NJets40);  
	  histNJets80[i]->Fill(events->NJets80);  
	  histNBtags[i]->Fill(events->NBJetsMedium);  
	  histMRVsRsq[i]->Fill(events->MR,events->Rsq);
	}
	if ((events->lep1+events->lep2).M() > 60 && (events->lep1+events->lep2).M() < 120) {
	  dataYield += 1.0;
	}
      } else {
	if ((events->lep1+events->lep2).M() > 80 && (events->lep1+events->lep2).M() < 100) {
	  histMR[i]->Fill(events->MR, weight );
	  histRsq[i]->Fill(events->Rsq, weight );
	  histDileptonMass[i]->Fill((events->lep1+events->lep2).M(), weight );      
	  
	  // histMET[i]->Fill(events->MET + fmax(0,random->Gaus(0,1.4)), weight);  
	  histMET[i]->Fill(events->MET, weight);  
	  histNJets40[i]->Fill(events->NJets40, weight);  
	  histNJets80[i]->Fill(events->NJets80, weight);  
	  histNBtags[i]->Fill(events->NBJetsMedium, weight);  
	  histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);
	}
	if ((events->lep1+events->lep2).M() > 60 && (events->lep1+events->lep2).M() < 120) {
	  MCYield += weight;
	}
      }
    }
  }

  cout << "here1\n";

  //--------------------------------------------------------------------------------------------------------------
  // Compute Expected Statistical Uncertainty
  //==============================================================================================================
  // TH2F *statUnc_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("statUnc_MRVsRsq"));
  // for (int i=0; i<statUnc_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
  //   for (int j=0; j<statUnc_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
  //     if (statUnc_MRVsRsq->GetBinContent(i,j) > 1) {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,1.0 / sqrt(statUnc_MRVsRsq->GetBinContent(i,j)));  	
  //     } else {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,0);
  //     }
  //   }
  // }

   //--------------------------------------------------------------------------------------------------------------
  // Subtract Non WJets Bkg
  //==============================================================================================================
  TH2F *DataMinusBkg_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("DataMinusBkg_MRVsRsq"));
  TH2F *MCToDataScaleFactor_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("MCToDataScaleFactor_MRVsRsq"));
  if (hasData) {

    for (int i=0; i<DataMinusBkg_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
      for (int j=0; j<DataMinusBkg_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
      
	double data = histMRVsRsq[0]->GetBinContent(i,j);
	double mc = 0; 
	double mc_StatErr = 0; 
	double bkg = 0;
	double bkg_StatErrSqr = 0;
	double bkg_SysErrSqr = 0;

	for (uint k=1; k < inputfiles.size(); ++k) {

	  if (processLabels[k] == "DY") {
	    mc = histMRVsRsq[k]->GetBinContent(i,j);
	    mc_StatErr = sqrt(histMRVsRsq[k]->GetBinError(i,j));
	    continue;
	  }

	  double systematicUncertainty = 0;
	  if (processLabels[k] == "VV") systematicUncertainty = 0.2;
	  if (processLabels[k] == "SingleTop") systematicUncertainty = 0.2;
	  if (processLabels[k] == "TT+V") systematicUncertainty = 0.2;
	  if (processLabels[k] == "WJets") systematicUncertainty = 0.2;
	  if (processLabels[k] == "TTJets") systematicUncertainty = 0.2;
 
	  bkg += histMRVsRsq[k]->GetBinContent(i,j);
	  bkg_StatErrSqr += pow(histMRVsRsq[k]->GetBinError(i,j),2);
	  bkg_SysErrSqr += pow( histMRVsRsq[k]->GetBinContent(i,j) * systematicUncertainty, 2);
	}

	DataMinusBkg_MRVsRsq->SetBinContent(i,j, data - bkg );
	double dataMinusBkgTotalErr = sqrt(data + bkg_StatErrSqr + bkg_SysErrSqr);
	DataMinusBkg_MRVsRsq->SetBinError(i,j, dataMinusBkgTotalErr );


	cout << "Bin " << DataMinusBkg_MRVsRsq->GetXaxis()->GetBinCenter(i) << " " << DataMinusBkg_MRVsRsq->GetYaxis()->GetBinCenter(j) << " : "
	     << data << " " << mc << " " << bkg << " " << mc_StatErr << " " << bkg_StatErrSqr << " " << bkg_SysErrSqr << "\n";

	MCToDataScaleFactor_MRVsRsq->SetBinContent(i,j, (data - bkg)/mc );
	MCToDataScaleFactor_MRVsRsq->SetBinError(i,j, ((data - bkg)/mc)*sqrt( pow(mc_StatErr/mc,2) + pow(dataMinusBkgTotalErr/(data-bkg),2)) );

      }
    }

   for (int i=0; i<DataMinusBkg_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
      for (int j=0; j<DataMinusBkg_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
	cout << "Bin " << DataMinusBkg_MRVsRsq->GetXaxis()->GetBinCenter(i) << " " << DataMinusBkg_MRVsRsq->GetYaxis()->GetBinCenter(j) << " : " << MCToDataScaleFactor_MRVsRsq->GetBinContent(i,j) << " +/- " << MCToDataScaleFactor_MRVsRsq->GetBinError(i,j) << "\n";	
      }
   }
  }







  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Yield inside Z Mass window 60-120\n";
  cout << "Data: " << dataYield << "\n";
  cout << "MC: " << MCYield << "\n";

  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================
  //PlotDataAndStackedBkg( histNVtx, processLabels, color, hasData, "NVtx", Label);
  PlotDataAndStackedBkg( histMR, processLabels, color, hasData, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, hasData, "Rsq", Label);
  PlotDataAndStackedBkg( histDileptonMass, processLabels, color, hasData, "DileptonMass", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, hasData, "MET", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, hasData, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, hasData, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, hasData, "NBtags", Label);
  


  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("ZToLLControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
 
 
  for(int i=0; i<int(histMRVsRsq.size()); i++) {
    file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  }
  //file->WriteTObject(statUnc_MRVsRsq,"statUnc_MRVsRsq_TTJets","WriteDelete");

  file->WriteTObject(DataMinusBkg_MRVsRsq, "DataMinusBkg_MRVsRsq", "WriteDelete");
  file->WriteTObject(MCToDataScaleFactor_MRVsRsq, "MCToDataScaleFactor_MRVsRsq", "WriteDelete");
  file->Close();
  delete file;       

}







void SelectZToLLControlSample( int option = 0) {

  string datafile = "";
  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> colors;

  //Inclusive sample
  if (option == 11) {
    datafile = "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root";
  } else if (option == 10) {
    datafile = "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DoubleElectron_GoodLumi.root";
  }
  if (option == 10 || option == 11) {
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_TTJets_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_VV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_SingleTop_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_TTV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted.root");
  }

  // //Razor Skim sample
  if (option == 1) {
    datafile = "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleMuParked_GoodLumi.root";
  } else if (option == 0) {
    datafile = "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_Data_DoubleElectron_GoodLumi.root";
  }
  if (option == 0 || option == 1) {
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DYJetsToLL_HTBinned_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTJets_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_VV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_SingleTop_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_WJetsToLNu_HTBinned_1pb_weighted.root");
  }

  processLabels.push_back("DY");
  processLabels.push_back("TTJets");  
  processLabels.push_back("VV");
  processLabels.push_back("SingleTop");
  processLabels.push_back("TT+V");
  processLabels.push_back("WJets");

  colors.push_back(kGreen+2);
  colors.push_back(kRed);
  colors.push_back(kOrange+1);
  colors.push_back(kBlack);
  colors.push_back(kGray);
  colors.push_back(kBlue);
 

  //*********************************************************************
  //EE Control Region
  //*********************************************************************
  if (option == 0) {
    RunSelectZToLLControlSample(datafile, inputfiles,processLabels,colors, 19789,"MR300Rsq0p05", 2,  "MR300Rsq0p05_ee");
  }
  if (option == 10) {
    RunSelectZToLLControlSample(datafile, inputfiles,processLabels,colors, 19789,"Inclusive", 2,  "Inclusive_ee");
  }
  //*********************************************************************
  //MM Control Region
  //*********************************************************************
  if ( option == 1) {
    RunSelectZToLLControlSample(datafile, inputfiles,processLabels, colors, 19751, "MR300Rsq0p05", 3, "MR300Rsq0p05_mumu");
  }
  if ( option == 11) {
    RunSelectZToLLControlSample(datafile, inputfiles,processLabels, colors, 19751, "Inclusive", 3, "Inclusive_mumu");    
  }
  
 
}


//**********************
//Mu-Mu Yields
//**********************
//Inclusive
//Data: 7.67962e+06
//MC: 7.87451e+06

//**********************
//E-E Yields
//**********************
// Inclusive
//Data: 5.42906e+06
//MC: 5.56361e+06

