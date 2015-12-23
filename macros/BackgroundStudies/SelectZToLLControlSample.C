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
#include <TLatex.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

#endif

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1D* NormalizeHist(TH1D *originalHist) {
  TH1D* hist = (TH1D*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}



void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

  double DataMean = 0;
  double DataRMS = 0;
  double MCMean = 0;
  double MCRMS = 0;


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
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.0);
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetMarkerStyle(20);      
      hist[0]->SetMarkerSize(1);
      hist[0]->SetLineWidth(1);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("pesame");
      DataMean = hist[0]->GetMean();
      DataRMS = hist[0]->GetRMS();
    }
    legend->Draw();

    for (uint i = 0 ; i < hist.size(); ++i) {
      if (processLabels[i] == "Data") {
	DataMean = hist[i]->GetMean();
	DataRMS = hist[i]->GetRMS();
      } 
      if (processLabels[i] == "DY") {
	MCMean = hist[i]->GetMean();
	MCRMS = hist[i]->GetRMS();
      } 
    }
  }

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  lumi_13TeV = "2.1 fb^{-1}";
  //lumi_13TeV = "Run257396-257400";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

  // TLatex *StatLabels  = new TLatex;
  // StatLabels->SetTextSize(0.03);
  // StatLabels->DrawLatexNDC(0.8,0.8, Form("Data Mean: %.1f",DataMean));
  // StatLabels->DrawLatexNDC(0.8,0.76, Form("Data RMS: %.1f",DataRMS));
  // StatLabels->DrawLatexNDC(0.8,0.72, Form("MC Mean: %.1f",MCMean));
  // StatLabels->DrawLatexNDC(0.8,0.68, Form("MC RMS: %.1f",MCRMS));


  cv->cd();
  cv->Update();


  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
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
  histDataOverMC->SetLineColor(kBlack);
  histDataOverMC->SetMarkerStyle(20);      
  histDataOverMC->SetMarkerSize(1);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("pe");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("Razor_ZToLLCR_%s%s_Logy.pdf",varName.c_str(),label.c_str()));


 

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


void RunSelectZToLLControlSample( vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {

  //****************************************
  //Set Plot Style
  //****************************************
  setTDRStyle();

  string Label = "";
  if (label != "") Label = "_" + label;
  TRandom3 *random = new TRandom3;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  TFile *NVtxWeightFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/NVtxReweight_ZToMuMu_2015D_2093ipb.root", "READ");
  TH1D *NVtxWeightHist = (TH1D*)NVtxWeightFile->Get("NVtxReweight");
  assert(NVtxWeightHist);

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupReweight_Spring15MCTo2015Data.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
  assert(pileupWeightHist);

  // TFile *eleEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20151013_PR_2015D_Golden_1264/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root","READ");
  TFile *eleEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightElectronSelectionEffDenominatorReco_2015Final_Golden.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
  assert(eleEffSFHist);

  // TFile *muonEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightMuonSelectionEffDenominatorReco_2015Final_Golden.root","READ");
  TFile *muonEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightMuonSelectionEffDenominatorReco_2015Final_Golden.root","READ");
  TH2D *muonEffSFHist = (TH2D*)muonEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
  assert(muonEffSFHist);

  TFile *eleTriggerEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015Final_Golden.root","READ");
  TH2D *eleTriggerEffSFHist = (TH2D*)eleTriggerEffSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
  assert(eleTriggerEffSFHist);

  TFile *muonTriggerEffSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015Final_Golden.root","READ");
  TH2D *muonTriggerEffSFHist = (TH2D*)muonTriggerEffSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");
  assert(muonTriggerEffSFHist);



  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 5;
  const int NRsqBins = 4;
  double MRBins[NMRBins] = {300, 350, 400, 450, 550};
  double RsqBins[NRsqBins] = {0.05, 0.075, 0.10, 1.5};

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;
  
  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);
  
  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
    inputfiles.push_back(bkgfiles[i]);
    processLabels.push_back(bkgLabels[i]);
    color.push_back(bkgColors[i]);
  }

  TH1D* histNVtxData = new TH1D( "histNVtxData", ";Number of Reconstructed Primary Vertices;Number of Events", 40, -0.5,39.5);
  TH1D* histNVtxAllBkg = new TH1D( "histNVtxAllBkg", ";Number of Reconstructed Primary Vertices;Number of Events", 40, -0.5,39.5);

  vector<TH1D*> histLepLeadPt;
  vector<TH1D*> histLepSubleadPt;
  vector<TH1D*> histLepPt;
  vector<TH1D*> histLepEta;
  vector<TH1D*> histNVtx;
  vector<TH1D*> histMR;
  vector<TH1D*> histRsq;
  vector<TH1D*> histDileptonMass;
  vector<TH1D*> histMET;
  vector<TH1D*> histMETPhi;
  vector<TH1D*> histMETNoHF;
  vector<TH1D*> histMETRaw;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histLepLeadPt.push_back(new TH1D(Form("histLepLeadPt_%s",processLabels[i].c_str()), "; Leading Lepton p_{T} [GeV/c]; Number of Events ", 50, 0, 200));
    histLepSubleadPt.push_back(new TH1D(Form("histLepSubleadPt_%s",processLabels[i].c_str()), "; Sub-leading Lepton p_{T} [GeV/c]; Number of Events ", 50, 0, 200));
    histLepPt.push_back(new TH1D(Form("histLepPt_%s",processLabels[i].c_str()), "; Lepton p_{T} [GeV/c]; Number of Events ", 50, 0, 200));
    histLepEta.push_back(new TH1D(Form("histLepEta_%s",processLabels[i].c_str()), "; Lepton #eta; Number of Events ", 50, -2.5, 2.5));

    histNVtx.push_back(new TH1D(Form("histNVtx_%s",processLabels[i].c_str()), "; Number of Primary Vertices; R^{2}; Number of Events", 100, -0.5, 99.5));
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Events / 50 GeV", 50, 200, 2700));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 20, 0.00, 1.0));
    histDileptonMass.push_back(new TH1D(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 60, 60, 120));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV/c] ; Number of Events", 50, 0, 200));
    histMETPhi.push_back(new TH1D(Form("histMETPhi_%s",processLabels[i].c_str()), "; METPhi [GeV/c] ; Number of Events", 50, -3.2, 3.1));
    histMETNoHF.push_back(new TH1D(Form("histMETNoHF_%s",processLabels[i].c_str()), "; MET (NoHF) [GeV/c] ; Number of Events", 50, 0, 200));
    histMETRaw.push_back(new TH1D(Form("histMETRaw_%s",processLabels[i].c_str()), "; MET (Raw) [GeV/c] ; Number of Events", 50, 0, 200));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 15, -0.5, 14.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 10, -0.5,9.5));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histLepLeadPt[i]->Sumw2();
    histLepSubleadPt[i]->Sumw2();
    histLepPt[i]->Sumw2();
    histLepEta[i]->Sumw2();
    histNVtx[i]->Sumw2();
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histDileptonMass[i]->Sumw2();
    histMET[i]->Sumw2();
    histMETPhi[i]->Sumw2();
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
    for (uint j=0; j < inputfiles[i].size(); ++j) {

      //for duplicate event checking
      map<pair<uint,uint>, bool > processedRunEvents;

      ControlSampleEvents *events = new ControlSampleEvents;
      events->LoadTree(inputfiles[i][j].c_str(), ControlSampleEvents::kTreeType_Dilepton_Full);
  
      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;
      
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
	events->tree_->GetEntry(ientry);
      
	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      
	//if (ientry > 1000000) break;

	double puWeight = 1;      
	double weight = 1;
	if (!isData) {
	  //puWeight = NVtxWeightHist->GetBinContent(NVtxWeightHist->GetXaxis()->FindFixBin(events->NPV));
	  puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(events->NPU_0));
	  weight = lumi * events->weight * puWeight;
	}


	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	//Use Single Lepton Triggers
	// if ( events->HLTDecision[2] || 
	//      events->HLTDecision[7] ||  
	//      events->HLTDecision[11] ||
	//      events->HLTDecision[12] ||
	//      events->HLTDecision[15]
	//      )  
	//   passTrigger = true;

	if (isData) {
	  if ( events->HLTDecision[22] || events->HLTDecision[23] || events->HLTDecision[24] || 
	       events->HLTDecision[25] || events->HLTDecision[26] ||
	       events->HLTDecision[27] || events->HLTDecision[28] || events->HLTDecision[29]	  
	       ) passTrigger = true;
	} else {
	  if ( events->HLTDecision[18] || events->HLTDecision[19] || events->HLTDecision[20] || 
	       events->HLTDecision[21] ||
	       events->HLTDecision[28] || events->HLTDecision[29]	  
	       ) passTrigger = true;
	}

	//DiMuon Triggers: Mu17Mu8 , Mu17TkMu8
	if (events->HLTDecision[39] || events->HLTDecision[41] || events->HLTDecision[43] ) passTrigger = true;

	// //DiElectron Triggers:
	// if ( events->HLTDecision[28] 
	// 	  || events->HLTDecision[29]	
	// 	  ) passTrigger = true;

	// //Razor Hadronic Triggers
	// if (isData) {
	// 	if ( events->HLTDecision[132] || events->HLTDecision[133] ) passTrigger = true;
	// } else {
	// 	if ( events->HLTDecision[136] || events->HLTDecision[137] ) passTrigger = true;
	// }

	if (!passTrigger) continue;

	//******************************
	//Selection Cuts 
	//******************************
	if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       &&
	       (abs(events->lep2Type) == 11 || abs(events->lep2Type) == 13)
	       )
	    ) continue;

	//******************************
	//apply scale corrections
	//******************************
	TLorentzVector l1;
	TLorentzVector l2;
	l1.SetPtEtaPhiM( events->lep1.Pt() , events->lep1.Eta(), events->lep1.Phi(), events->lep1.M());
	l2.SetPtEtaPhiM( events->lep2.Pt() , events->lep2.Eta(), events->lep2.Phi(), events->lep2.M());
     
	//Single lepton triggered data requires leading lepton pt cut
	if ( !( l1.Pt() > 30 || l2.Pt() > 30)) continue;
	if (! (l1.Pt() > 20 && l2.Pt() > 20
	       && events->lep1PassTight && events->lep2PassTight)
	    ) continue;


	// if(! (fabs(l1.Eta()) < 1.42 && fabs(l2.Eta()) < 1.42 )) continue;
	//if(! (fabs(l1.Eta()) > 1.5 && fabs(l2.Eta()) > 1.5)) continue;

	//dilepton mass cut
	if ( (events->lep1+events->lep2).M() < 50) continue;

	// BTag Veto
	//if ( !( events->NBJetsMedium == 0)) continue;


	//Check for duplicate data events
	if (isData) {
	  //cout << "event: " << events->run << " " << events->event << "\n";
	  if(!(processedRunEvents.find(make_pair(events->run, events->event)) == processedRunEvents.end())) {
	    cout << "Duplicate event: " << events->run << " " << events->event << "\n";
	    continue;
	  } else {
	    processedRunEvents[make_pair(events->run, events->event)] = true;
	    //cout << processedRunEvents.size() << "\n";
	  }
	}


	// //Razor signal region cuts
	// if (option == "TwoJet80" || option == "MR300Rsq0p05" || option == "MR300") {
	//   if (!(events->NJets80 >= 2 )) continue;
	// }
      
	// if (option == "MR300Rsq0p05") {
	//   if (!(events->MR > 300 && events->Rsq > 0.05 )) continue;
	// }
      
	// if (option == "MR300") {
	//   if (!(events->MR > 300 )) continue;
	// }
      
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
	if (!isData) {
	  double triggerEffScaleFactor = 1.0;
	
	  double leptonEffScaleFactor = 1.0;
	  if (abs(events->lep1Type) == 11  ) {
	    leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),20.01)),
								 eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep1.Eta()))
								 );	 
	  }
	  if (abs(events->lep2Type) == 11 ) {
	    leptonEffScaleFactor *= eleEffSFHist->GetBinContent(  eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep2.Pt(),199.9),20.01)),
								  eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep2.Eta()))
								  );
	  }
	  if (abs(events->lep1Type) == 13) {
	    leptonEffScaleFactor *= muonEffSFHist->GetBinContent( muonEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),20.01)),
								  muonEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep1.Eta()))
								  );	 
	  }
	  if (abs(events->lep2Type) == 13) {
	    leptonEffScaleFactor *= muonEffSFHist->GetBinContent(  muonEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep2.Pt(),199.9),20.01)),
								   muonEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep2.Eta()))
								   );
	  }
	
	  // cout << events->lep1.Pt() << " " << events->lep1.Eta() << " : " << eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),15.01)) << " "
	  //      << eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep1.Eta())) << " : "
	  //      << eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep1.Pt(),199.9),15.01)),
	  // 				     eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep1.Eta()))
	  // 				     )
	  //      << " : " 
	  //      << abs(events->lep1Type) << " "
	  //      << "\n";
	  // cout << events->lep2.Pt() << " " << events->lep2.Eta() << " : " << eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep2.Pt(),199.9),15.01)) << " "
	  //      << eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep2.Eta())) << " : "
	  //      << eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(events->lep2.Pt(),199.9),15.01)),
	  // 				     eleEffSFHist->GetYaxis()->FindFixBin(fabs(events->lep2.Eta()))
	  // 				     )
	  //      << " : " 
	  //      << abs(events->lep2Type) << " "
	  //      << "\n";
	  // cout << "lepton SF: " << leptonEffScaleFactor << " " << weight << " \n";
	    
 

	  weight *= leptonEffScaleFactor;
	  weight *= triggerEffScaleFactor;

	}


	//******************************
	//Fill histograms
	//******************************
	if (isData) {
	  histNVtxData->Fill(events->NPV);
	  if ((l1+l2).M() > 60 && (l1+l2).M() < 120) {
	    histLepLeadPt[i]->Fill(fmax(l1.Pt(),l2.Pt()));
	    histLepSubleadPt[i]->Fill(fmin(l1.Pt(),l2.Pt()));
	    histLepPt[i]->Fill(l1.Pt());
	    histLepPt[i]->Fill(l2.Pt());
	    histLepEta[i]->Fill(l1.Eta());
	    histLepEta[i]->Fill(l2.Eta());

	    histNVtx[i]->Fill(events->NPV);	  
	    histMR[i]->Fill(events->MR);
	    histRsq[i]->Fill(events->Rsq);
	    histDileptonMass[i]->Fill((l1+l2).M());      
	    histMET[i]->Fill(events->MET);  
	    histMETPhi[i]->Fill(events->METPhi);  
	    histMETNoHF[i]->Fill(events->METnoHF);  
	    histMETRaw[i]->Fill(events->METRaw);  
	    histNJets40[i]->Fill(events->NJets40);  
	    histNJets80[i]->Fill(events->NJets80);  
	    histNBtags[i]->Fill(events->NBJetsMedium);  
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq);
	  }
	  if ((l1+l2).M() > 60 && (l1+l2).M() < 120) {
	    dataYield += 1.0;
	  }
	} else {
	  histNVtxAllBkg->Fill(events->NPV, weight);
	  if ((l1+l2).M() > 60 && (l1+l2).M() < 120) {
	    histLepLeadPt[i]->Fill(fmax(l1.Pt(),l2.Pt()),weight);
	    histLepSubleadPt[i]->Fill(fmin(l1.Pt(),l2.Pt()),weight);
	    histLepPt[i]->Fill(l1.Pt(),weight);
	    histLepPt[i]->Fill(l2.Pt(),weight);
	    histLepEta[i]->Fill(l1.Eta(),weight);
	    histLepEta[i]->Fill(l2.Eta(),weight);

	    histNVtx[i]->Fill(events->NPV,weight);
	    histMR[i]->Fill(events->MR, weight );
	    histRsq[i]->Fill(events->Rsq, weight );
	    histDileptonMass[i]->Fill((l1+l2).M(), weight );      
	  
	    histMET[i]->Fill(events->MET, weight);  
	    histMETPhi[i]->Fill(events->METPhi, weight);  
	    histMETNoHF[i]->Fill(events->METnoHF, weight);  
	    histMETRaw[i]->Fill(events->METRaw, weight);  
	    histNJets40[i]->Fill(events->NJets40, weight);  
	    histNJets80[i]->Fill(events->NJets80, weight);  
	    histNBtags[i]->Fill(events->NBJetsMedium, weight);  
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);
	  }
	  if ((l1+l2).M() > 60 && (l1+l2).M() < 120) {
	    MCYield += weight;
	  }
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


  //*******************************************************************************************
  //Make NVtx Reweighting Function
  //*******************************************************************************************
  TH1D *NVtxDataNormalized = NormalizeHist( histNVtxData );
  TH1D *NVtxBkgNormalized = NormalizeHist( histNVtxAllBkg );

  TH1D *NVtxReweight = (TH1D*)NVtxDataNormalized->Clone("NVtxReweight");  
  for (int i=1; i<NVtxReweight->GetXaxis()->GetNbins()+1; i++) {

    double data = 0;
    double bkg = 0;
    if (NVtxBkgNormalized->GetBinContent(i) > 0) {
      NVtxReweight->SetBinContent(i,NVtxDataNormalized->GetBinContent(i)/NVtxBkgNormalized->GetBinContent(i));
    } else if (NVtxDataNormalized->GetBinContent(i) == 0){
      NVtxReweight->SetBinContent(i,0.0);
    } else {
      NVtxReweight->SetBinContent(i,1.0);
    }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Subtract Non WJets Bkg
  //==============================================================================================================
  TH2F *DataMinusBkg_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("DataMinusBkg_MRVsRsq"));
  TH2F *MCToDataScaleFactor_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("MCToDataScaleFactor_MRVsRsq"));

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
  PlotDataAndStackedBkg( histLepLeadPt, processLabels, color, true, "LepLeadPt", Label);
  PlotDataAndStackedBkg( histLepSubleadPt, processLabels, color, true, "LepSubleadPt", Label);
  PlotDataAndStackedBkg( histLepPt, processLabels, color, true, "LepPt", Label);
  PlotDataAndStackedBkg( histLepEta, processLabels, color, true, "LepEta", Label);

  PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
  PlotDataAndStackedBkg( histDileptonMass, processLabels, color, true, "DileptonMass", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET", Label);
  PlotDataAndStackedBkg( histMETPhi, processLabels, color, true, "METPhi", Label);
  PlotDataAndStackedBkg( histMETNoHF, processLabels, color, true, "METNoHF", Label);
  PlotDataAndStackedBkg( histMETRaw, processLabels, color, true, "METRaw", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
  PlotDataAndStackedBkg( histNVtx, processLabels, color, true, "NVtx", Label);
  


  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open("NVtxReweight_ZToMuMu.root", "UPDATE");
  file->cd();
  file->WriteTObject(NVtxReweight, "NVtxReweight", "WriteDelete");
  file->Close();
  delete file;


  file = TFile::Open(("ZToLLControlRegionPlots"+Label+".root").c_str(), "UPDATE");
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

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  //Inclusive sample
  if (option == 1 || option == 11) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleMuon_Run2015D_GoodLumiGolden.root");
    // datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_New/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DoubleMuon_Run2015D_Oct05ReMiniAOD_GoodLumiGolden.root");
    // datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_New/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DoubleMuon_Run2015D_PRv4_GoodLumiGolden.root");

  } else if (option == 0 || option == 10) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleElectron_Run2015D_GoodLumiGolden.root");
  } else if (option == 2 || option == 12) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleMuon_Run2015D_GoodLumiGolden.root");
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleElectron_Run2015D_GoodLumiGolden.root"); 
  }

  vector<string> bkgfiles_dy;
  bkgfiles_dy.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_New/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");

  vector<string> bkgfiles_ttbar;
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");
  vector<string> bkgfiles_vv;
  bkgfiles_vv.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted.root");
  vector<string> bkgfiles_singletop;
  bkgfiles_singletop.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted.root");
  vector<string> bkgfiles_wjets;
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");

  bkgfiles.push_back(bkgfiles_dy);
  bkgfiles.push_back(bkgfiles_ttbar);
  bkgfiles.push_back(bkgfiles_vv);
  bkgfiles.push_back(bkgfiles_singletop);
  // bkgfiles.push_back(bkgfiles_wjets);
  
  processLabels.push_back("DY");
  processLabels.push_back("TTJets");  
  processLabels.push_back("VV");
  processLabels.push_back("SingleTop");
  // processLabels.push_back("WJets");
  // processLabels.push_back("TT+V");

  colors.push_back(kGreen+2);
  colors.push_back(kAzure+10);
  colors.push_back(kGray);
  colors.push_back(kBlue);
  // colors.push_back(kRed);
  // colors.push_back(kOrange+1);


   double lumi = 0;
   lumi = 2137;

  //*********************************************************************
  //EE Control Region
  //*********************************************************************
  if (option == 0) {
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels,colors, lumi,"MR300Rsq0p05", 2,  "MR300Rsq0p05_ee");
  }
  if (option == 10) {
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels,colors, lumi,"Inclusive", 2,  "Inclusive_ee");
  }
  //*********************************************************************
  //MM Control Region
  //*********************************************************************
  if ( option == 1) {
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels, colors, lumi, "MR300Rsq0p05", 3, "MR300Rsq0p05_mumu");
  }
  if ( option == 11) {
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels, colors, lumi, "Inclusive", 3, "Inclusive_mumu");    
  }
  
  //*********************************************************************
  //EE+MM Control Region
  //*********************************************************************
  if ( option == 2) {    
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels, colors, lumi, "MR300", 1, "MR300_all");
  }
  if ( option == 12) {
    RunSelectZToLLControlSample(datafiles, bkgfiles,processLabels, colors, lumi, "Inclusive", 1, "Inclusive_all");    
  }


}


//**********************
//Mu-Mu Yields
//**********************
// Single Mu Triggers ( 30 - 20 )
// Data: 1.52257e+06
// MC: 1.43049e+06 (new PU weights, only DY)

// Double Mu Triggers (30-20)
// Data: 1.36959e+06
// MC: 1.37669e+06

// DY: 99.4%
// TTBAR: 0.3%
// ST: 0.03%
// W+jets: nothing
// VV: 0.2%

//Only IsoMu
//Data: 1.51388e+06
//MC: 1.42669e+06 (DY only)

//Only IsoTkMu
//Data: 1.49898e+06
//MC: 1.42098e+06

//30-30
//Data: 1.19584e+06
//MC: 1.1206e+06

// Re-miniAOD (575 ipb)
// Data: 396431
// MC: 385646

// 1016
// Data: 726978
// MC: 681420


// 1562 
// Data: 1.12677e+06
// MC: 1.04762e+06
//1562 Double Mu
//Data: 1.01021e+06
//MC: 1.00009e+06



//**********************
//E-E Yields
//**********************
// Single Ele Triggers ( 30 - 20 )
// Data: 931286
// MC: 920790 ( Nov result )
//MC: 927400 (only DY, new PU weights)

