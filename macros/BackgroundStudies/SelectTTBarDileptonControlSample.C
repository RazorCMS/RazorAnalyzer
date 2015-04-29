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


void PlotDataAndStackedBkg( vector<TH1F*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

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

  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1F *histDataOverMC = (TH1F*)hist[0]->Clone("histDataOverMC");

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
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
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
  histDataOverMC->GetXaxis()->SetTitleSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.2);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("e1");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s.gif",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s_Logy.gif",varName.c_str(),label.c_str()));

}



//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarDileptonControlSample( string datafile, vector<string> bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  // TFile *file = TFile::Open("data.root", "UPDATE");
  // file->cd();
  
  // TTree *tree = new TTree("tree", "tree");
  // tree->Branch("MR",&MR,"MR/F");
  // tree->Branch("Rsq",&Rsq,"Rsq/F");
  // tree->Branch("mll",&mll,"mll/F");



  bool printdebug = false;

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/ElectronSelection_Run2012ReReco_53X.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  assert(eleEffSFHist);

  TFile *DYScaleFactorsFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/ZToLLScaleFactors.root", "READ");
  TH2D *DYScaleFactorsHist = (TH2D*)DYScaleFactorsFile->Get("ZToLLDileptonScaleFactor");
  assert(DYScaleFactorsHist);


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 7;
  const int NRsqBins = 7;
  double MRBins[NMRBins] = {300, 350, 400, 450, 500, 550, 700};
  double RsqBins[NRsqBins] = {0.15,0.175,0.20,0.225,0.25,0.30,1.5};

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

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histMuonPt;
  vector<TH1F*> histElectronPt;
  vector<TH1F*> histDileptonDeltaPhi;
  vector<TH1F*> histDileptonPt;
  vector<TH1F*> histDileptonEta;
  vector<TH1F*> histMET;
  vector<TH1F*> histDileptonMass;
  vector<TH1F*> histDileptonCharge;
  vector<TH1F*> histNJets40;
  vector<TH1F*> histNJets80;
  vector<TH1F*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 50, 0, 1500));
    histRsq.push_back(new TH1F(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 25, 0, 0.5));
    histMuonPt.push_back(new TH1F(Form("histMuonPt_%s",processLabels[i].c_str()), "; MuonPt [GeV/c] ; Number of Events", 20, 0, 400));
    histElectronPt.push_back(new TH1F(Form("histElectronPt_%s",processLabels[i].c_str()), "; ElectronPt [GeV/c] ; Number of Events", 20, 0, 400));
    histDileptonDeltaPhi.push_back(new TH1F(Form("histDileptonDeltaPhi_%s",processLabels[i].c_str()), "; DileptonDeltaPhi [GeV/c] ; Number of Events", 25, 0, 3.1416));
    histDileptonPt.push_back(new TH1F(Form("histDileptonPt_%s",processLabels[i].c_str()), "; DileptonPt [GeV/c] ; Number of Events", 50, 0, 500));
    histDileptonEta.push_back(new TH1F(Form("histDileptonEta_%s",processLabels[i].c_str()), "; DileptonEta [GeV/c] ; Number of Events", 25, -10, 10));
    histMET.push_back(new TH1F(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histDileptonMass.push_back(new TH1F(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 1000));
    histDileptonCharge.push_back(new TH1F(Form("histDileptonCharge_%s",processLabels[i].c_str()), "; Charge; Number of Events", 5, -2.5, 2.5));
    histNJets40.push_back(new TH1F(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 10, -0.5, 9.5));
    histNJets80.push_back(new TH1F(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1F(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 5, -0.5, 4.5));
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histMuonPt[i]->Sumw2();
    histElectronPt[i]->Sumw2();
    histDileptonDeltaPhi[i]->Sumw2();
    histDileptonPt[i]->Sumw2();
    histDileptonEta[i]->Sumw2();
    histDileptonMass[i]->Sumw2();
    histDileptonCharge[i]->Sumw2();
    histMET[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNBtags[i]->Sumw2();  

    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histMRVsRsq[i]->Sumw2();
  }
 
  double dataYield = 0;
  double MCYield = 0;

  vector<pair<UInt_t,UInt_t> > RunAndEvent;

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
      //remove duplicate events in data
      //******************************
      // if ((hasData && i==0)) {
      // 	bool isDuplicateEvent = false;
      // 	for(int k=0; k<RunAndEvent.size(); ++k) {	
      // 	  if (RunAndEvent[k].first == events->run && RunAndEvent[k].second == events->event) {
      // 	    //cout << "Found duplicate: " << events->run << " " << events->event << "\n";
      // 	    isDuplicateEvent = true;
      // 	    break;
      // 	  } 
      // 	}
      // 	if (isDuplicateEvent) {
      // 	  continue;
      // 	} else {
      // 	  //cout << "Add " << events->run << " " << events->event << "\n";
      // 	  RunAndEvent.push_back( pair<UInt_t,UInt_t>(events->run, events->event));
      // 	}
      // }


      //******************************
      //Trigger Selection
      //******************************
      bool passTrigger = false;

      //DiMuon Triggers: Mu17Mu8 , Mu17TkMu8      
      if (events->HLTDecision[3] ==true || events->HLTDecision[4] ==true) passTrigger = true;
      //DiElectron Triggers:
      if (events->HLTDecision[12] ==true) passTrigger = true;
      //MuEG Triggers: Mu17Ele8 , Mu8Ele17
      if (events->HLTDecision[6] ==true || events->HLTDecision[7] ==true) passTrigger = true;

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
      if ( (events->lep1+events->lep2).M() < 50) continue;

      //MET cut
      if (option == "Met>40" ||option == "topEnhanced" || option == "TwoJet80" || option == "MR300Rsq0p15"  || option == "BumpRegion" || option == "LooserBumpRegion"
	  ) {
	if ( events->MET < 40 ) continue;
      }

      //b-tagging      
      if (option == "topEnhanced" || option == "TwoJet80" || option == "MR300Rsq0p15"  || option == "BumpRegion" || option == "LooserBumpRegion"
	  ) {
	// if ( !( (events->jet1PassCSVLoose || events->jet2PassCSVLoose)
	// 	&& events->jet1.Pt() > 30 && events->jet2.Pt() > 30)
	//      ) continue;
	if (!(events->NBJetsMedium >= 1)) continue;	
      }
      
      //Z-mass window cut
      if ( abs(events->lep1Type) == abs(events->lep2Type) 
      	   && 
      	   (events->lep1+events->lep2).M() > 76 && (events->lep1+events->lep2).M() < 106
      	   ) continue;
      

      // //Razor signal region cuts
      if (option == "TwoJet80") {
	if (!(events->MR > 0 )) continue;
      }
      
      if (option == "MR300Rsq0p15" || option == "BumpRegion") {
	if (!(events->MR > 300 && events->Rsq > 0.15 )) continue;
      }
      
// Bin 500 -0.01875 : 0 0 0 0 0 0
// Bin 500 0.175 : 32 38.9447 2.43988 1.17867 0.732361 0.111041
// Bin 500 0.225 : 13 15.8852 2.23035 0.936241 0.701604 0.0887993
// Bin 500 0.275 : 7 8.10963 0.648978 0.794302 0.194775 0.00819985
// Bin 500 0.9 : 6 5.76329 2.09984 0.728553 0.933627 0.132999
// Bin 500 0.175 : 0.759027 +/- 0.148937
// Bin 500 0.225 : 0.677968 +/- 0.237164
// Bin 500 0.275 : 0.783146 +/- 0.339717
// Bin 500 0.9 : 0.676724 +/- 0.469115


      //******************************
      //ChannelOptions
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
	if ( (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 13)) triggerEffScaleFactor = 0.9700; 

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
	if ( ((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13)) || 
	     ((abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11))
	     ) normalizationScaleFactor = 0.970;


	//b-tagging scale factors
	double btagScaleFactor = 1.0;
	double bjet1EventScaleFactor = 1.0;
	double bjet2EventScaleFactor = 1.0;
	if (events->bjet1.Pt() > 20) {
	  double bjet1ScaleFactor = 0.938887 + 0.00017124 * events->bjet1.Pt() + (-2.76366e-07) * events->bjet1.Pt() * events->bjet1.Pt() ;
	  double MCEff = 1.0;
	  if (events->bjet1.Pt() < 50) MCEff = 0.65;
	  else if (events->bjet1.Pt() < 80) MCEff = 0.70;
	  else if (events->bjet1.Pt() < 120) MCEff = 0.73;
	  else if (events->bjet1.Pt() < 210) MCEff = 0.73;
	  else MCEff = 0.66;				 
	  if (events->bjet1PassMedium) bjet1EventScaleFactor = bjet1ScaleFactor;
	  else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
	}
	if (events->bjet2.Pt() > 20) {
	  double bjet2ScaleFactor = 0.938887 + 0.00017124 * events->bjet2.Pt() + (-2.76366e-07) * events->bjet2.Pt() * events->bjet1.Pt() ;
	  double MCEff = 1.0;
	  if (events->bjet2.Pt() < 50) MCEff = 0.65;
	  else if (events->bjet2.Pt() < 80) MCEff = 0.70;
	  else if (events->bjet2.Pt() < 120) MCEff = 0.73;
	  else if (events->bjet2.Pt() < 210) MCEff = 0.73;
	  else MCEff = 0.66;				 
	  if (events->bjet2PassMedium) bjet2EventScaleFactor = bjet2ScaleFactor;
	  else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
	}
	btagScaleFactor = bjet1EventScaleFactor * bjet2EventScaleFactor;

	// cout << events->NBJetsMedium << " : " << events->bjet1.Pt() << " " << events->bjet1PassMedium << " " << events->bjet2.Pt() << " " << events->bjet2PassMedium  
	//      << " : " << bjet1EventScaleFactor << " " << bjet2EventScaleFactor << " " << btagScaleFactor << " "
	//      <<  " \n";


	weight = weight * leptonEffScaleFactor;
	weight = weight * triggerEffScaleFactor;
	weight = weight * normalizationScaleFactor;
	weight = weight * btagScaleFactor;

	//cout << weight << " " << leptonEffScaleFactor << " " << triggerEffScaleFactor << " " <<  events->weight << " " << puWeight << "\n";

	// if ((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13) ||
	//     (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11)) weight = weight * (733.0 / 902.0);
	// if (abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 13) weight = weight * (302.0 / 409.0);
	// if (abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 11) weight = weight * (263.0 / 318.0);

	if (processLabels[i] == "DY") {
	  weight *= DYScaleFactorsHist->GetBinContent( DYScaleFactorsHist->GetXaxis()->FindFixBin(fmin(fmax(events->MR,300.1),549.9)) ,  
						       DYScaleFactorsHist->GetYaxis()->FindFixBin(fmin(fmax(events->Rsq,0.0501),1.499)) );	 
	}
       
      }


      //******************************
      //Fill histograms
      //******************************
      if (hasData && i==0) {
	dataYield += 1.0;
	histDileptonMass[i]->Fill((events->lep1+events->lep2).M());      
	histMET[i]->Fill(fmin(events->MET,399.9));

	if (abs(events->lep1Type) == 13) {
	  histMuonPt[i]->Fill(events->lep1.Pt());
	  histElectronPt[i]->Fill(events->lep2.Pt());
	} else {
	  histElectronPt[i]->Fill(events->lep1.Pt());
	  histMuonPt[i]->Fill(events->lep2.Pt());
	}
	histDileptonDeltaPhi[i]->Fill( acos(cos(events->lep1.Phi() - events->lep2.Phi())) );
	histDileptonPt[i]->Fill( (events->lep1+events->lep2).Pt() );
	histDileptonEta[i]->Fill( (events->lep1+events->lep2).Eta() );
	histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type));
	histNJets40[i]->Fill( events->NJets40 );
	histNJets80[i]->Fill( events->NJets80 );
	histNBtags[i]->Fill( events->NBJetsLoose );

	if (events->MR > 0) {
	  histMR[i]->Fill(events->MR);
	  histRsq[i]->Fill(events->Rsq);
	}

	histMRVsRsq[i]->Fill(events->MR,events->Rsq);

      } else {
	MCYield += weight;
	histDileptonMass[i]->Fill((events->lep1+events->lep2).M(), weight );      
	histMET[i]->Fill(fmin(events->MET,399.9), weight );
	if (abs(events->lep1Type) == 13) {
	  histMuonPt[i]->Fill(events->lep1.Pt(), weight);
	  histElectronPt[i]->Fill(events->lep2.Pt(), weight);
	} else {
	  histElectronPt[i]->Fill(events->lep1.Pt(), weight);
	  histMuonPt[i]->Fill(events->lep2.Pt(), weight);
	}
	histDileptonDeltaPhi[i]->Fill( acos(cos(events->lep1.Phi() - events->lep2.Phi())) , weight);
	histDileptonPt[i]->Fill( (events->lep1+events->lep2).Pt() , weight);
	histDileptonEta[i]->Fill( (events->lep1+events->lep2).Eta() , weight);

	histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type) , weight);
	histNJets40[i]->Fill( events->NJets40 , weight);
	histNJets80[i]->Fill( events->NJets80 , weight);
	histNBtags[i]->Fill( events->NBJetsLoose , weight);
	
	if (events->MR > 0) {
	  histMR[i]->Fill(events->MR, weight );
	  histRsq[i]->Fill(events->Rsq, weight );
	}

	histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);

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
  // Subtract Non TTBar Bkg
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

	  if (processLabels[k] == "TTJets") {
	    mc = histMRVsRsq[k]->GetBinContent(i,j);
	    mc_StatErr = sqrt(histMRVsRsq[k]->GetBinError(i,j));
	    continue;
	  }

	  double systematicUncertainty = 0;
	  if (processLabels[k] == "VV") systematicUncertainty = 0.2;
	  if (processLabels[k] == "SingleTop") systematicUncertainty = 0.2;
	  if (processLabels[k] == "TT+V") systematicUncertainty = 0.2;
	  if (processLabels[k] == "DY") systematicUncertainty = 0.2;
 
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
  PlotDataAndStackedBkg( histMR, processLabels, color, hasData, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, hasData, "Rsq", Label);
  PlotDataAndStackedBkg( histDileptonMass, processLabels, color, hasData, "DileptonMass", Label);
  PlotDataAndStackedBkg( histMuonPt, processLabels, color, hasData, "MuonPt", Label);
  PlotDataAndStackedBkg( histElectronPt, processLabels, color, hasData, "ElectronPt", Label);
  PlotDataAndStackedBkg( histDileptonDeltaPhi, processLabels, color, hasData, "DileptonDeltaPhi", Label);
  PlotDataAndStackedBkg( histDileptonPt, processLabels, color, hasData, "DileptonPt", Label);
  PlotDataAndStackedBkg( histDileptonEta, processLabels, color, hasData, "DileptonEta", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, hasData, "MET", Label);
  PlotDataAndStackedBkg( histDileptonCharge, processLabels, color, hasData, "DileptonCharge", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, hasData, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, hasData, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, hasData, "NBtags", Label);
 


  // //*******************************************************************************************
  // //MR Vs Rsq
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // histMRVsRsq[0]->SetStats(false);
  // histMRVsRsq[0]->Draw("colz");
  // cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_EventCountsTTJets_MRVsRsq%s.gif",Label.c_str()));

  // cv = new TCanvas("cv","cv", 800,600);
  // statUnc_MRVsRsq->SetStats(false);
  // statUnc_MRVsRsq->Draw("colz");
  // cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_StatUncTTJets_MRVsRsq%s.gif",Label.c_str()));

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Yield inside Z Mass window 60-120\n";
  cout << "Data: " << dataYield << "\n";
  cout << "MC: " << MCYield << "\n";

  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("TTBarDileptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonCharge[i], Form("histDileptonCharge_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  }


  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonCharge[i], Form("histDileptonCharge_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
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







void SelectTTBarDileptonControlSample( int option = 0) {

  string datafile = "";
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_DoubleElectron_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_MuEG_GoodLumi.root";

  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_MuEG_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_HTAndHTMHTParked_GoodLumi.root";
  //string datafile = "";
  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> colors;

  //No Skims  
  if (option == 10 || option == 20) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_MuEG_GoodLumi.root";
  if (option == 11) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DoubleElectron_GoodLumi.root";
  if (option == 12) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root";

  if (option == 10 || option == 11 || option == 12 || option == 20) {
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_TTJets_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_VV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_SingleTop_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_TTV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkim/RunOneRazorControlRegions_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted.root");
  }

  //MR300 Skim
  // if (option == 0) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_MuEG_GoodLumi.root";
  if (option == 0) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_MuEG_GoodLumi.root";
  if (option == 1) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DoubleElectron_GoodLumi.root";
  if (option == 2) datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DoubleMuParked_GoodLumi.root";

  if (option == 0 || option == 1 || option == 2) {
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTJets_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_VV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_SingleTop_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_TTV_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DYJetsToLL_HTBinned_1pb_weighted.root");
    inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_WJetsToLNu_HTBinned_1pb_weighted.root");  
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_1pb_weighted.root");
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_VV_1pb_weighted.root");
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_SingleTop_1pb_weighted.root");
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_TTV_1pb_weighted.root");
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonRazorSkim/RunOneRazorControlRegions_DileptonRazorSkim_DYJetsToLL_HTBinned_1pb_weighted.root");
    // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/DileptonSkimNoPUJetID/RunOneRazorControlRegions_DileptonSkim_WJetsToLNu_HTBinned_1pb_weighted.root");  
  }

  processLabels.push_back("TTJets");  
  processLabels.push_back("VV");
  processLabels.push_back("SingleTop");
  processLabels.push_back("TT+V");
  processLabels.push_back("DY");
  processLabels.push_back("WJets");

  colors.push_back(kRed);
  colors.push_back(kOrange+1);
  colors.push_back(kBlack);
  colors.push_back(kGray);
  colors.push_back(kGreen+2);
  colors.push_back(kBlue);
 

  //*********************************************************************
  //E-Mu Control Region
  //*********************************************************************
  // RunSelectTTBarDileptonControlSample(datafile, inputfiles, processLabels, colors, 19780,"Inclusive",0,"Inclusive_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"Met>40",0,"MetGreaterThan40_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"topEnhanced",0,"TopEnhanced_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"TwoJet80",0,"TwoJet80_emu");
  
  if (option == 0) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"MR300Rsq0p15",0,"MR300Rsq0p15_emu");
  if (option == 10) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"topEnhanced",0,"TopEnhanced_emu");
  if (option == 20) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"Inclusive",0,"Inclusive_emu");

  //*********************************************************************
  //E-E Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"Inclusive",2,"Inclusive_ee");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"topEnhanced",2,"TopEnhanced_ee");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"TwoJet80",2,"TwoJet80_ee");
  if (option == 1) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"MR300Rsq0p15",2,"MR300Rsq0p15_ee");
  if (option == 11) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"topEnhanced",2,"TopEnhanced_ee");

  //*********************************************************************
  //MuMu Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"Inclusive",3,"Inclusive_mumu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"topEnhanced",3,"TopEnhanced_mumu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"TwoJet80",3,"TwoJet80_mumu");
  if (option == 2) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"MR300Rsq0p15",3,"MR300Rsq0p15_mumu");
  if (option == 12) RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"topEnhanced",3,"TopEnhanced_mumu");





}


//**********************
//E-Mu Yields
//**********************
//Inclusive
// Data: 80236
// MC: 79823.6

//Met>40
//Data: 44398
// MC: 44765.7

//TwoJet80
// Data: 4579
// MC: 5488.65

//MR300Rsq0p15
// Data: 579
// MC: 678.504

//**********************
//Mu-Mu Yields
//**********************
//Inclusive
//Data: 7.67962e+06
//MC: 7.87451e+06

//TwoJet80
//Data: 5731
//MC: 6126.44

//**********************
//E-E Yields
//**********************
// Inclusive
//Data: 5.42906e+06
//MC: 5.56361e+06

//TwoJet80
//Data: 263
//MC: 318



// Z->MM
// Data: 7.47318e+06
// MC: 7.62973e+06 Powheg
// MC: 7.62596e+06 Madgraph

// Z->EE
//Data: 5.27829e+06
//MC: 5.37693e+06 Madgraph


