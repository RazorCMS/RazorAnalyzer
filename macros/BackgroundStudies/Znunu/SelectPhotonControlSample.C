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
    double dataErrSqr = 0;
    if (hasData) {
      data = hist[0]->GetBinContent(b);
      dataErrSqr = pow(hist[0]->GetBinError(b),2);
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
      //histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
      histDataOverMC->SetBinError(b, (data / MC)*sqrt(dataErrSqr/pow(data,2) + MCErrSqr/pow(MC,2) ));
    } else {
      histDataOverMC->SetBinContent(b, 0);
      histDataOverMC->SetBinError(b, 0);
    }
    cout << "bin " << b << " : " << data << " +/- " << sqrt(dataErrSqr) << " / " << MC << " +/- " << sqrt(MCErrSqr)
	 << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
  }

  histDataOverMC->GetYaxis()->SetTitle("Data/MC");
  histDataOverMC->GetYaxis()->SetNdivisions(306);
  histDataOverMC->GetYaxis()->SetTitleSize(0.10);
  histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
  //  histDataOverMC->GetYaxis()->SetRangeUser(0.0,3.0);
  histDataOverMC->GetYaxis()->SetRangeUser(0.0,2.0);
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


void RunSelectPhotonControlSample(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, 
				    int SFOption, string option, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;

  TFile *SFInputFile = 0;
  if (SFOption == 0) SFInputFile = 0;
  else if (SFOption == 1) SFInputFile = TFile::Open("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_GJetsInv.root", "READ");
  else if (SFOption == 2) SFInputFile = TFile::Open("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", "READ");
  TH2F *InputSFHist = 0;
  if (SFInputFile) {
    InputSFHist = (TH2F*)SFInputFile->Get("GJetsInvScaleFactors");
  }
  

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 6;
  const int NRsqBins = 6;
  const int NNJetBins = 4;
  const int NNBTagBins = 5;
  double MRBins[NMRBins] = {400, 500, 700, 900, 1200, 4000};
  double RsqBins[NRsqBins] = {0.25, 0.30, 0.41, 0.52, 0.64, 1.5};
  double NJetBins[NNJetBins] = {0,4,7,20};
  double NBTagBins[NNBTagBins] = {0,1,2,3,4};

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
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", NMRBins-1, MRBins));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", NRsqBins-1, RsqBins));
    histMRRsqUnrolled.push_back(new TH1D(Form("histMRRsqUnrolled_%s",processLabels[i].c_str()), "; Bin Number ; Number of Events", (NMRBins-1)*(NRsqBins-1), 0, (NMRBins-1)*(NRsqBins-1)));
    histPhotonPt.push_back(new TH1D(Form("histPhotonPt_%s",processLabels[i].c_str()), "; Photon p_{T} [GeV/c] ; Number of Events", 80, 0, 400));    
    histPhotonEta.push_back(new TH1D(Form("histPhotonEta_%s",processLabels[i].c_str()), "; Photon #eta ; Number of Events", 50, -3, 3));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", NNJetBins-1, NJetBins));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", NNBTagBins-1,NBTagBins));
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
	if (!isData && !(option != "Inclusive" && processLabels[i] == "QCD")) {
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
	  if(processLabels[i]  == "GJets" && events->minDRGenPhotonToParton < 0.4) continue;
	  if(processLabels[i]  == "GJetsFrag" && isFake) continue;	  
	  if(processLabels[i]  == "GJetsFrag" && events->minDRGenPhotonToParton >= 0.4) continue;
	  if(processLabels[i]  == "GJets" || processLabels[i]  == "GJetsFrag") weight *= 1.55; //Overall K-Factor & photon efficiency correction

	  if (option == "Inclusive") {
	    if(processLabels[i] == "QCD" && !isFake) continue;
	    if(processLabels[i]  == "QCD") weight *= 1.37; //Overall K-Factor & photon fake rate correction
	  } else {
	    if(processLabels[i] == "QCD") {
	      weight = 0.05;
	    }
	  }
	}

	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	if (isData || (option != "Inclusive" && processLabels[i] == "QCD") ) {	  
	  double dataWeight = 1;
	  if (events->pho1.Pt() > 185) {
	    dataWeight = 1;
	    if (events->HLTDecision[93] && events->pho1HLTFilter[22]) passTrigger = true;
	  }
	  else if (events->pho1.Pt() > 135) {
	    dataWeight = events->HLTPrescale[92];
	    if (events->HLTDecision[92] && events->pho1HLTFilter[23]) passTrigger = true;
	  } else if (events->pho1.Pt() > 105) {
	    dataWeight = events->HLTPrescale[91];
	    if (events->HLTDecision[91] && events->pho1HLTFilter[24]) passTrigger = true;
	  } else if (events->pho1.Pt() > 85) {
	    dataWeight = events->HLTPrescale[90];
	    if (events->HLTDecision[90] && events->pho1HLTFilter[25]) passTrigger = true;
	  } else if (events->pho1.Pt() > 58){
	    dataWeight = events->HLTPrescale[89];
	    if (events->HLTDecision[89] && events->pho1HLTFilter[26]) passTrigger = true;
	  } else {
	    dataWeight = events->HLTPrescale[88];
	    if (events->HLTDecision[88] && events->pho1HLTFilter[27]) passTrigger = true;
	  } 
	  
	  if (isData) {
	    weight = dataWeight;
	  }

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
	if (fabs(events->pho1.Eta()) < 1.5) {
	  if (! (events->pho1_sigmaietaieta < 0.0103)) continue;
	} else {
	  if (! (events->pho1_sigmaietaieta < 0.0271)) continue;
	}
	if (! (events->pho1_chargediso < 2.5)) continue;
	if (! (events->pho1.Pt() > 50)) continue;
  
	//By default it's EG Loose selection + chargedIso < 2.5 GeV.
	//We can tighten the cuts if we wish...

	if (option == "MR300Rsq0p15" ) {
	  if (!(events->NJets80_NoPho >= 2 && events->MR_NoPho > 400 && events->Rsq_NoPho > 0.25 )) continue;	
	}
    	if (option == "MR300Rsq0p15_4Jet" ) {
	  if (!(events->NJets80_NoPho >= 2 && events->NJets_NoPho >= 4 && events->MR_NoPho > 400 && events->Rsq_NoPho > 0.25 )) continue;	
	}
      
	//MET Filters
	if (!(events->Flag_HBHENoiseFilter && events->Flag_goodVertices && events->Flag_eeBadScFilter)) continue;


	// //******************************
	// //Apply Scale Factors
	// //******************************
	if (!isData && !(option != "Inclusive" && processLabels[i] == "QCD") ) {
	   double razorSF = 1.0;
	   if (processLabels[i] == "GJets" && InputSFHist) {
	     razorSF = InputSFHist->GetBinContent( InputSFHist->GetXaxis()->FindFixBin( fmin(fmax(events->MR_NoPho,400.1),3999)),
						   InputSFHist->GetYaxis()->FindFixBin( fmin(fmax(events->Rsq_NoPho,0.251),1.49))
						   );
	   }
	   weight *= razorSF;
	 }



	//******************************
	//Fill histograms
	//******************************		
	int MRBin = histMRVsRsq[i]->GetXaxis()->FindFixBin( events->MR_NoPho );
	int RsqBin = histMRVsRsq[i]->GetYaxis()->FindFixBin( events->Rsq_NoPho );
	int MRRsqBin = (NRsqBins-1)*(MRBin-1) + RsqBin-1 ;

	if (isData) {

	  if (passTrigger) {
	    dataYield += weight;
	    histPhotonPt[i]->Fill(events->pho1.Pt(), weight);
	    histPhotonEta[i]->Fill(events->pho1.Eta(), weight);
	    histMET[i]->Fill(events->MET, weight);
	    histNJets80[i]->Fill(events->NJets80_NoPho, weight);
	    histNJets40[i]->Fill(events->NJets_NoPho, weight);
	    histNBtags[i]->Fill(fmin(events->NBJetsMedium,3), weight);
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
  // Calculate Scale Factor For GJets
  //==============================================================================================================
  if ( SFOption == 0 ) {

    TH2F* HistSF = (TH2F*)histMRVsRsq[0]->Clone("GJetsInvScaleFactors");
    TH2F* HistSFUp = (TH2F*)histMRVsRsq[0]->Clone("GJetsInvScaleFactorsUp");
    TH2F* HistSFDown = (TH2F*)histMRVsRsq[0]->Clone("GJetsInvScaleFactorsDown");

    for (int i=1; i <= (NMRBins-1)*(NRsqBins-1); ++i) {
    
      int bin_i = int ((i-1) / (NMRBins-1)) + 1;
      int bin_j = bin_j = ((i-1) % (NMRBins-1)) + 1;
  
    
      cout << "Bin: " << i << " -> " << bin_i << " , " << bin_j << "\n";

      double N_data = 0;
      double N_QCD = 0;
      double NErr_QCD = 0;
      double N_GJetsFrag = 0;
      double NErr_GJetsFrag = 0;
      double N_GJetsDirect = 0;
      double NErr_GJetsDirect = 0;
      double N_Other = 0;
      double NErr_Other = 0;

      for (uint j=0; j < inputfiles.size(); ++j) {
	if ( processLabels[j] == "Data") {
	  N_data = histMRRsqUnrolled[j]->GetBinContent(i);
	  N_QCD = 0.05 * N_data; //Fake Fraction is 5% with 5% systematic
	  NErr_QCD = 0.05 * N_data;
	}
	if ( processLabels[j] == "GJets") {
	  N_GJetsDirect = histMRRsqUnrolled[j]->GetBinContent(i);
	  NErr_GJetsDirect = histMRRsqUnrolled[j]->GetBinError(i);
	}
	if ( processLabels[j] == "GJetFrag") {
	  N_GJetsFrag = histMRRsqUnrolled[j]->GetBinContent(i);
	  NErr_GJetsFrag = histMRRsqUnrolled[j]->GetBinError(i);
	}      
	if ( processLabels[j] == "Other") {
	  N_Other = histMRRsqUnrolled[j]->GetBinContent(i);
	  NErr_Other = sqrt( pow(histMRRsqUnrolled[j]->GetBinError(i),2) + pow( 0.50 * N_Other , 2)); 
	}      
      }
     
      double N_dataMinusBkg = N_data - N_QCD - N_GJetsFrag - NErr_Other;
      double NErr_dataMinusBkg = sqrt( N_data + pow(NErr_QCD,2) + pow(NErr_GJetsFrag,2) + pow( NErr_Other, 2));
      double SF = N_dataMinusBkg / N_GJetsDirect;
      double SFErr = sqrt( ( pow(NErr_dataMinusBkg,2) * pow(N_GJetsDirect,2) + pow(NErr_GJetsDirect,2)*pow(N_dataMinusBkg,2)) / (pow(N_GJetsDirect,4)));

      cout << "Data-Bkg : " << N_dataMinusBkg << " +/- " << NErr_dataMinusBkg << "\n";
      cout << "MC (GJets) : " << N_GJetsDirect << " +/- " << NErr_GJetsDirect << "\n";
      cout << "SF : " << SF << " +/- " << SFErr << "\n";
      cout << "\n\n";

      HistSF->SetBinContent( bin_i, bin_j, SF);
      HistSF->SetBinError( bin_i, bin_j, SFErr);

    }
  
    TFile *SFFile = TFile::Open("RazorScaleFactors_GJets.root", "UPDATE");
    SFFile->WriteTObject(HistSF, "GJetsInvScaleFactors", "WriteDelete");
    SFFile->Close();
  }

  if ( SFOption == 2 ) {

    TH1D* HistSysUnc_NBTag = (TH1D*)histNBtags[0]->Clone("GJetsInvScaleFactors_Sys_NBTag");
    for (int i=1; i < HistSysUnc_NBTag->GetXaxis()->GetNbins() + 1; ++i) {
      double N_data = 0;
      double N_QCD = 0;
      double N_GJetsFrag = 0;
      double N_GJetsDirect = 0;
      double N_Other = 0;
      for (uint j=0; j < inputfiles.size(); ++j) {
	if ( processLabels[j] == "Data") {
	  N_data = histNBtags[j]->GetBinContent(i);
	  N_QCD = 0.05 * N_data; //Fake Fraction is 5% with 5% systematic
	}
	if ( processLabels[j] == "GJets") {
	  N_GJetsDirect = histNBtags[j]->GetBinContent(i);
	}
	if ( processLabels[j] == "GJetFrag") {
	  N_GJetsFrag = histNBtags[j]->GetBinContent(i);
	}      
	if ( processLabels[j] == "Other") {
	  N_Other = histNBtags[j]->GetBinContent(i);
	}
      }
      double N_dataMinusQCD = N_data * 0.95;
      double N_MCNoQCD = N_GJetsDirect + N_GJetsFrag + N_Other;
      double SFSystematics = fabs ( 1 - N_dataMinusQCD / N_MCNoQCD);
      HistSysUnc_NBTag->SetBinContent( i , SFSystematics );
      HistSysUnc_NBTag->SetBinError( i , 0 );

      cout << "Data-QCD : " << N_dataMinusQCD << "\n";
      cout << "MC (excluding QCD) : " << N_GJetsDirect << " + " << N_GJetsFrag << " + " << N_Other << " = " << N_MCNoQCD << "\n";
      cout << "Systematic Uncertainty : " << SFSystematics << "\n";
      cout << "\n\n";

    }
    TFile *SFFile = TFile::Open("RazorZNuNuBTagClosureTests.root", "UPDATE");
    SFFile->WriteTObject(HistSysUnc_NBTag, "ZNuNuBTagClosureSysUnc", "WriteDelete");
    SFFile->Close();
  }


  //--------------------------------------------------------------------------------------------------------------
  // Plot Event Density for MR-Rsq spectrum 
  //==============================================================================================================
  for (int i=1; i <= (NMRBins-1)*(NRsqBins-1); ++i) {
    
    int bin_i = int ((i-1) / (NMRBins-1)) + 1;
    int bin_j = bin_j = ((i-1) % (NMRBins-1)) + 1;        
    cout << "Bin: " << i << " -> " << bin_i << " , " << bin_j << " : " 
	 << "MR: " << histMRVsRsq[0]->GetXaxis()->GetBinLowEdge(bin_i) << " - " << histMRVsRsq[0]->GetXaxis()->GetBinUpEdge(bin_i) << " , "
	 << "Rsq: " << histMRVsRsq[0]->GetYaxis()->GetBinLowEdge(bin_j) << " - " << histMRVsRsq[0]->GetYaxis()->GetBinUpEdge(bin_j) << " : "
	 << (histMRVsRsq[0]->GetXaxis()->GetBinUpEdge(bin_i) - histMRVsRsq[0]->GetXaxis()->GetBinLowEdge(bin_i))  << " * "
	 << (histMRVsRsq[0]->GetYaxis()->GetBinUpEdge(bin_j) - histMRVsRsq[0]->GetYaxis()->GetBinLowEdge(bin_j))
	 << "\n";
    double binVolume = (histMRVsRsq[0]->GetXaxis()->GetBinUpEdge(bin_i) - histMRVsRsq[0]->GetXaxis()->GetBinLowEdge(bin_i)) *
      (histMRVsRsq[0]->GetYaxis()->GetBinUpEdge(bin_j) - histMRVsRsq[0]->GetYaxis()->GetBinLowEdge(bin_j));
         
    for (uint j=0; j < inputfiles.size(); ++j) {
      histMRRsqUnrolled[j]->SetBinContent( i , histMRRsqUnrolled[j]->GetBinContent(i) / binVolume );
      histMRRsqUnrolled[j]->SetBinError( i , histMRRsqUnrolled[j]->GetBinError(i) / binVolume );
    }        

  }
    

  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  // PlotDataAndStackedBkg( histPhotonPt, processLabels, color, true, "PhotonPt", Label);
  // PlotDataAndStackedBkg( histPhotonEta, processLabels, color, true, "PhotonEta", Label);
  // PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET_NoPho", Label);
  // PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  // PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  // PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
  // PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  // PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
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
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden.root");     
  } else {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RazorSkim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden_RazorSkim.root");     
  }

  vector<string> bkgfiles_gjets;
  vector<string> bkgfiles_gjetsFrag;
  vector<string> bkgfiles_qcd;
  vector<string> bkgfiles_other;

  if (option >= 10) {
    bkgfiles_gjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8_1pb_weighted.root");    
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root"); 
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_1pb_weighted.root");     
  } else {
    bkgfiles_gjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RazorSkim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_RazorSkim.root");    
    bkgfiles_gjetsFrag.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RazorSkim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_RazorSkim.root");    
    bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RazorSkim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden_RazorSkim.root"); 
    bkgfiles_other.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RazorSkim/RunTwoRazorControlRegions_PhotonFull_Other_1pb_weighted_RazorSkim.root"); 
  }
   

  bkgfiles.push_back(bkgfiles_gjets);
  bkgfiles.push_back(bkgfiles_gjetsFrag);
  bkgfiles.push_back(bkgfiles_qcd);
  bkgfiles.push_back(bkgfiles_other);

  processLabels.push_back("GJets");  
  processLabels.push_back("GJetsFrag");  
  processLabels.push_back("QCD");
  processLabels.push_back("Other");

  colors.push_back(kOrange);
  colors.push_back(kOrange+4);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  
  double lumi = 2185;

  //*********************************************************************
  //GJets Control Region
  //*********************************************************************
   if (option == 0) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,0,"MR300Rsq0p15","PhotonCR");
   }
   if (option == 1) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,1,"MR300Rsq0p15","PhotonCRTestSF");
   }
   if (option == 2) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,2,"MR300Rsq0p15_4Jet","PhotonCRTestSFMultiJet");
   }
   if (option == 10) {
     RunSelectPhotonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,0,"Inclusive","Inclusive");
   }


}



//**********************
//With Photon ID + Iso Cuts
//**********************
// YieldPho36_58To70 : 6.79226e+06
// YieldPho50_58To70 : 6.90576e+06
// Ratio : 1.01671

// YieldPho50_85To95 : 1.21875e+06
// YieldPho75_85To95 : 1.22005e+06
// Ratio : 1.00107

// YieldPho75_105To115 : 488278
// YieldPho90_105To115 : 489032
// Ratio : 1.00154

// YieldPho90_135To145 : 161048
// YieldPho120_135To145 : 160721
// Ratio : 0.99797

// YieldPho120_185To200 : 53750
// YieldPho165_185To200 : 53762
// Ratio : 1.00022


//**********************
//After Razor Cuts
//**********************
// YieldPho36_58To70 : 4500
// YieldPho50_58To70 : 4510
// Ratio : 1.00222

// YieldPho50_85To95 : 6936
// YieldPho75_85To95 : 6320
// Ratio : 0.911188

// YieldPho75_105To115 : 6750
// YieldPho90_105To115 : 7023
// Ratio : 1.04044

// YieldPho90_135To145 : 5926
// YieldPho120_135To145 : 6004
// Ratio : 1.01316

// YieldPho120_185To200 : 5481
// YieldPho165_185To200 : 5475
// Ratio : 0.998905

