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
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

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

  legend = new TLegend(0.75,0.50,0.90,0.84);
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
    }
    legend->Draw();
  }

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  // lumi_13TeV = "42 pb^{-1}";
  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

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
  histDataOverMC->GetYaxis()->SetRangeUser(0.0,3.0);
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
  cv->SaveAs(Form("Razor_PhotonControlRegion_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("Razor_PhotonControlRegion_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_PhotonControlRegion_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("Razor_PhotonControlRegion_%s%s_Logy.pdf",varName.c_str(),label.c_str()));


 

}



//=== MAIN MACRO ================================================================================================= 


void RunSelectPhotonControlSample(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 5;
  const int NRsqBins = 5;
  double MRBins[NMRBins] = {300, 400, 500, 700, 4000};
  double RsqBins[NRsqBins] = {0.15, 0.25, 0.30, 0.41, 1.5};

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

  vector<TH1D*> histMR;
  vector<TH1D*> histRsq;
  vector<TH1D*> histPhotonPt;
  vector<TH1D*> histPhotonEta;
  vector<TH1D*> histMET;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH1D*> histMRRsqUnrolled;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 40, 300, 2300));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 50, 0.0, 1.50));
    histMRRsqUnrolled.push_back(new TH1D(Form("histMRRsqUnrolled_%s",processLabels[i].c_str()), "; Bin Number ; Number of Events", NMRBins*NRsqBins, 0, NMRBins*NRsqBins));
    histPhotonPt.push_back(new TH1D(Form("histPhotonPt_%s",processLabels[i].c_str()), "; Photon p_{T} [GeV/c] ; Number of Events", 80, 0, 400));    
    histPhotonEta.push_back(new TH1D(Form("histPhotonEta_%s",processLabels[i].c_str()), "; Photon #eta ; Number of Events", 50, -3, 3));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 10, -0.5, 9.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 5, -0.5, 4.5));
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histMRRsqUnrolled[i]->Sumw2();
    histPhotonPt[i]->Sumw2();
    histPhotonEta[i]->Sumw2();
    histMET[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNJets80[i]->Sumw2();
    histNBtags[i]->Sumw2();  

    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histMRVsRsq[i]->Sumw2();
  }
 

  TH1D* PhotonPt_HLTPho22 = new TH1D( "PhotonPt_HLTPho22", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho30 = new TH1D( "PhotonPt_HLTPho30", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho36 = new TH1D( "PhotonPt_HLTPho36", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho50 = new TH1D( "PhotonPt_HLTPho50", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho75 = new TH1D( "PhotonPt_HLTPho75", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho90 = new TH1D( "PhotonPt_HLTPho90", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho120 = new TH1D( "PhotonPt_HLTPho120", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);
  TH1D* PhotonPt_HLTPho165 = new TH1D( "PhotonPt_HLTPho165", "; Photon p_{T} [GeV/c]; Number of Events", 150, 0, 300);

  double YieldPho22_35To40 = 0;
  double YieldPho30_35To40 = 0;
  double YieldPho30_42To50 = 0;
  double YieldPho36_42To50 = 0;
  double YieldPho36_58To70 = 0;
  double YieldPho50_58To70 = 0;
  double YieldPho50_85To95 = 0;
  double YieldPho75_85To95 = 0;
  double YieldPho75_105To115 = 0;
  double YieldPho90_105To115 = 0;
  double YieldPho90_135To145 = 0;
  double YieldPho120_135To145 = 0;
  double YieldPho120_185To200 = 0;
  double YieldPho165_185To200 = 0;


  double dataYield = 0;
  double MCYield = 0;


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {

    //for duplicate event checking
    map<pair<uint,uint>, bool > processedRunEvents;

    for (uint j=0; j < inputfiles[i].size(); ++j) {
      ControlSampleEvents *events = new ControlSampleEvents;
      events->LoadTree(inputfiles[i][j].c_str(), ControlSampleEvents::kTreeType_Photon_Full);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;
    
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
	events->tree_->GetEntry(ientry);
      

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	double puWeight = 1;      
	double weight = 1;
	if (!isData) {
	  weight = lumi * events->weight;
	}

	if (isnan(events->weight) || isinf(events->weight)) {
	  continue;
	  cout << "...bad event: " << weight << " " << "\n";
	}

	//******************************
	// remove double-counting fakes
	//******************************
	if (!isData) {
	  bool isFake = false;	   
	  if(abs(events->pho1_motherID) > 5 && events->pho1_motherID != 21 && events->pho1_motherID != 2212) isFake = true;	  
	  if(processLabels[i]  == "GJets" && isFake) continue;
	  if(processLabels[i] == "QCD" && !isFake) continue;
	}

	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	if (isData) {	  
	  double dataWeight = 1;
	  if (events->pho1.Pt() > 185) {
	    dataWeight = 1;
	    if (events->HLTDecision[93]) passTrigger = true;
	  } 
	  else if (events->pho1.Pt() > 135) {
	    dataWeight = events->HLTPrescale[92];
	    if (events->HLTDecision[92]) passTrigger = true;
	  } else if (events->pho1.Pt() > 105) {
	    dataWeight = events->HLTPrescale[91];
	    if (events->HLTDecision[91]) passTrigger = true;
	  } else if (events->pho1.Pt() > 85) {
	    dataWeight = events->HLTPrescale[90];
	    if (events->HLTDecision[90]) passTrigger = true;
	  } else {
	    dataWeight = events->HLTPrescale[89];
	    if (events->HLTDecision[89]) passTrigger = true;
	  } 

	  weight = dataWeight;

	} else {
	  if (
	       events->HLTDecision[88] || events->HLTDecision[89] || events->HLTDecision[90] 
	       || events->HLTDecision[91] || events->HLTDecision[92] || 
	      events->HLTDecision[93]
	      ) passTrigger = true;
	}

	if (!passTrigger) continue;

	//******************************
	//Selection Cuts 
	//******************************
	//Photon selection
	if (! (events->pho1.Pt() > 50)) continue;
  
	//By default it's EG Loose selection + chargedIso < 2.5 GeV.
	//We can tighten the cuts if we wish...

	if (option == "MR300Rsq0p15" ) {
	  if (!(events->NJets80 >= 2 && events->MR_NoPho > 300 && events->Rsq_NoPho > 0.15 )) continue;	
	}
      
	//MET Filters
	if (!(events->Flag_HBHENoiseFilter && events->Flag_goodVertices && events->Flag_eeBadScFilter)) continue;


	// //******************************
	// //Apply Scale Factors
	// //******************************
	// if (!isData) {
	//   double razorSF = 1.0;
	//   if (processLabels[i] == "TTJets") {
	//     razorSF = TTBarSFHist->GetBinContent( TTBarSFHist->GetXaxis()->FindFixBin( fmin(fmax(events->MR,400.1),3999)),
	// 					  TTBarSFHist->GetYaxis()->FindFixBin( fmin(fmax(events->Rsq,0.151),1.49))
	// 					  );
	//   }
	//   if (processLabels[i] == "WJets") {
	//     razorSF = WJetsSFHist->GetBinContent( WJetsSFHist->GetXaxis()->FindFixBin( fmin(fmax(events->MR,400.1),3999)),
	// 					  WJetsSFHist->GetYaxis()->FindFixBin( fmin(fmax(events->Rsq,0.151),1.49))
	// 					  );
	//   }
	//   // cout << processLabels[i] << " " << events->MR << " " << events->Rsq << " " << razorSF << "\n";
	//   weight *= razorSF;
	// }


	//******************************
	//Fill histograms
	//******************************		
	int MRBin = histMRVsRsq[i]->GetXaxis()->FindFixBin( events->MR );
	int RsqBin = histMRVsRsq[i]->GetYaxis()->FindFixBin( events->Rsq );
	int MRRsqBin = NRsqBins*(MRBin-1) + RsqBin;

	if (isData) {

	  if (passTrigger) {
	    dataYield += weight;
	    histPhotonPt[i]->Fill(events->pho1.Pt(), weight);
	    histPhotonEta[i]->Fill(events->pho1.Eta(), weight);
	    histMET[i]->Fill(events->MET, weight);
	    histNJets80[i]->Fill(events->NJets80_NoPho, weight);
	    histNJets40[i]->Fill(events->NJets_NoPho, weight);
	    histNBtags[i]->Fill(events->NBJetsMedium, weight);
	    histMR[i]->Fill(events->MR_NoPho, weight);
	    histRsq[i]->Fill(events->Rsq_NoPho, weight);
	    histMRRsqUnrolled[i]->Fill(MRRsqBin+0.5, weight);
	  }

	 
	  //************************************************************
	  //Histograms to compute scale factors for each trigger path
	  //************************************************************		
	  if (events->HLTDecision[88]) {
	    PhotonPt_HLTPho36->Fill( events->pho1.Pt(), events->HLTPrescale[88]);
	    if (events->pho1.Pt() > 42 && events->pho1.Pt() < 50) YieldPho36_42To50 += events->HLTPrescale[88];
	    if (events->pho1.Pt() > 58 && events->pho1.Pt() < 70) YieldPho36_58To70 += events->HLTPrescale[88];	    
	  }
	  if (events->HLTDecision[89]) {
	    PhotonPt_HLTPho50->Fill( events->pho1.Pt(), events->HLTPrescale[89]);    
	    if (events->pho1.Pt() > 58 && events->pho1.Pt() < 70) YieldPho50_58To70 += events->HLTPrescale[89];	    
	    if (events->pho1.Pt() > 85 && events->pho1.Pt() < 95) YieldPho50_85To95 += events->HLTPrescale[89];	    
	  }
	  if (events->HLTDecision[90]) {
	    PhotonPt_HLTPho75->Fill( events->pho1.Pt(), events->HLTPrescale[90]);
	    if (events->pho1.Pt() > 85 && events->pho1.Pt() < 95) YieldPho75_85To95 += events->HLTPrescale[90];	    
	    if (events->pho1.Pt() > 105 && events->pho1.Pt() < 115) YieldPho75_105To115 += events->HLTPrescale[90];	    
	  }
	  if (events->HLTDecision[91]) {
	    PhotonPt_HLTPho90->Fill( events->pho1.Pt(), events->HLTPrescale[91]);
	    if (events->pho1.Pt() > 105 && events->pho1.Pt() < 115) YieldPho90_105To115 += events->HLTPrescale[91];	    
	    if (events->pho1.Pt() > 135 && events->pho1.Pt() < 145) YieldPho90_135To145 += events->HLTPrescale[91];	    
	  }
	  if (events->HLTDecision[92]) {
	    PhotonPt_HLTPho120->Fill( events->pho1.Pt(), events->HLTPrescale[92]);
	    if (events->pho1.Pt() > 135 && events->pho1.Pt() < 145) YieldPho120_135To145 += events->HLTPrescale[92];	    
	    if (events->pho1.Pt() > 185 && events->pho1.Pt() < 200) YieldPho120_185To200 += events->HLTPrescale[92];	    
	  }
	  if (events->HLTDecision[93]) {
	    PhotonPt_HLTPho165->Fill( events->pho1.Pt(), events->HLTPrescale[93]);
	    if (events->pho1.Pt() > 185 && events->pho1.Pt() < 200) YieldPho165_185To200 += events->HLTPrescale[93];	    
	  }


	} else {
	  MCYield += weight;
	  histPhotonPt[i]->Fill(events->pho1.Pt(), weight);
	  histPhotonEta[i]->Fill(events->pho1.Eta(), weight);
	  histMET[i]->Fill(events->MET, weight);
	  histNJets80[i]->Fill(events->NJets80_NoPho, weight);
	  histNJets40[i]->Fill(events->NJets_NoPho, weight);
	  histNBtags[i]->Fill(events->NBJetsMedium, weight);
	  histMR[i]->Fill(events->MR_NoPho, weight);
	  histRsq[i]->Fill(events->Rsq_NoPho, weight);
	  histMRRsqUnrolled[i]->Fill(MRRsqBin+0.5, weight);
	}


      } //loop over events
    } //loop over input files
  } //loop over input file groups


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
  PlotDataAndStackedBkg( histPhotonPt, processLabels, color, true, "PhotonPt", Label);
  PlotDataAndStackedBkg( histPhotonEta, processLabels, color, true, "PhotonEta", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET_NoPho", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
  PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
  PlotDataAndStackedBkg( histMRRsqUnrolled, processLabels, color, true, "MRRsqUnrolled", Label);
  

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  // cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Selected Event Yield \n";
  cout << "Data: " << dataYield << "\n";
  cout << "MC: " << MCYield << "\n";

  cout << "YieldPho22_35To40 : " << YieldPho22_35To40 << "\n";
  cout << "YieldPho30_35To40 : " << YieldPho30_35To40 << "\n";
  cout << "Ratio : " << YieldPho30_35To40 / YieldPho22_35To40 << "\n\n";
  cout << "YieldPho30_42To50 : " << YieldPho30_42To50 << "\n";
  cout << "YieldPho36_42To50 : " << YieldPho36_42To50 << "\n";
  cout << "Ratio : " << YieldPho36_42To50 / YieldPho30_42To50 << "\n\n";
  cout << "YieldPho36_58To70 : " << YieldPho36_58To70 << "\n";
  cout << "YieldPho50_58To70 : " << YieldPho50_58To70 << "\n";
  cout << "Ratio : " << YieldPho50_58To70 / YieldPho36_58To70 << "\n\n";
  cout << "YieldPho50_85To95 : " << YieldPho50_85To95 << "\n";
  cout << "YieldPho75_85To95 : " << YieldPho75_85To95 << "\n";
  cout << "Ratio : " << YieldPho75_85To95 / YieldPho50_85To95 << "\n\n";
  cout << "YieldPho75_105To115 : " << YieldPho75_105To115 << "\n";
  cout << "YieldPho90_105To115 : " << YieldPho90_105To115 << "\n";
  cout << "Ratio : " << YieldPho90_105To115 / YieldPho75_105To115 << "\n\n";
  cout << "YieldPho90_135To145 : " << YieldPho90_135To145 << "\n";
  cout << "YieldPho120_135To145 : " << YieldPho120_135To145 << "\n";
  cout << "Ratio : " << YieldPho120_135To145 / YieldPho90_135To145 << "\n\n";
  cout << "YieldPho120_185To200 : " << YieldPho120_185To200 << "\n";
  cout << "YieldPho165_185To200 : " << YieldPho165_185To200 << "\n";
  cout << "Ratio : " << YieldPho165_185To200 / YieldPho120_185To200 << "\n\n";

  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("PhotonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  file->WriteTObject( PhotonPt_HLTPho22, "PhotonPt_HLTPho22" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho30, "PhotonPt_HLTPho30" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho36, "PhotonPt_HLTPho36" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho50, "PhotonPt_HLTPho50" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho75, "PhotonPt_HLTPho75" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho90, "PhotonPt_HLTPho90" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho120, "PhotonPt_HLTPho120" , "WriteDelete");
  file->WriteTObject( PhotonPt_HLTPho165, "PhotonPt_HLTPho165" , "WriteDelete");

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histPhotonPt[i], Form("histPhotonPt_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histPhotonEta[i], Form("histPhotonEta_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMET[i], Form("histMET_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets80[i], Form("histNJets80_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMRRsqUnrolled[i], Form("histMRRsqUnrolled_%s",processLabels[i].c_str()), "WriteDelete");
  }


  // for(int i=0; i<int(histMRVsRsq.size()); i++) {
  //   file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  // }

  file->Close();
  delete file;       

}







void SelectPhotonControlSample( int option = 0) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";


  //No Skims  
  if (option >= 10) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden.root");     
  } else {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden_RazorSkim.root");     
  }

  vector<string> bkgfiles_gjets;
  vector<string> bkgfiles_qcd;
  vector<string> bkgfiles_other;

  if (option >= 10) {
    bkgfiles_gjets.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8_1pb_weighted.root");    
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root");     
  } else {
    bkgfiles_gjets.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_RazorSkim.root");    
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_QCD_HTBinned_1pb_weighted_RazorSkim.root"); 
    bkgfiles_other.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_Other_1pb_weighted_RazorSkim.root"); 
  }
   

  bkgfiles.push_back(bkgfiles_gjets);
  bkgfiles.push_back(bkgfiles_qcd);
  bkgfiles.push_back(bkgfiles_other);

  processLabels.push_back("GJets");  
  processLabels.push_back("QCD");
  processLabels.push_back("Other");

  colors.push_back(kOrange);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  
  double lumi = 2185;

  //*********************************************************************
  //GJets Control Region
  //*********************************************************************
   if (option == 0) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",0,"MR300Rsq0p15");
   }
   if (option == 10) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",0,"Inclusive");
   }


}



//**********************
//With Photon ID + Iso Cuts
//**********************
// YieldPho36_58To70 : 1.30574e+07
// YieldPho50_58To70 : 1.30185e+07
// Ratio : 0.997014

// YieldPho50_85To95 : 2.09166e+06
// YieldPho75_85To95 : 1.9931e+06
// Ratio : 0.95288

// YieldPho75_105To115 : 767340
// YieldPho90_105To115 : 768620
// Ratio : 1.00167

// YieldPho90_135To145 : 239080
// YieldPho120_135To145 : 230440
// Ratio : 0.963861

// YieldPho120_185To200 : 71070
// YieldPho165_185To200 : 70911
// Ratio : 0.997763

//**********************
//After Razor Cuts
//**********************
// YieldPho36_58To70 : 18000
// YieldPho50_58To70 : 10620
// Ratio : 0.59

// YieldPho50_85To95 : 10520
// YieldPho75_85To95 : 9120
// Ratio : 0.86692

// YieldPho75_105To115 : 11520
// YieldPho90_105To115 : 12480
// Ratio : 1.08333

// YieldPho90_135To145 : 7990
// YieldPho120_135To145 : 7545
// Ratio : 0.944305

// YieldPho120_185To200 : 7230
// YieldPho165_185To200 : 6711
// Ratio : 0.928216
