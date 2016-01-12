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
#include <TFractionFitter.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

#endif


void RunPhotonControlSample_Fit(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;

  TFile *FakeFile = new TFile("FakePhotonTemplate-new.root", "READ");
  TH1F *fakeHist = (TH1F*)FakeFile->Get("h1");

  TFile *PromptFile = new TFile("PromptPhotonTemplate-new3.root", "READ");
  TH1F *promptHist = (TH1F*)PromptFile->Get("h1");

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  // const int NMRBins = 5;
  // const int NRsqBins = 5;
  double MRBins[] = {300, 400, 500, 700, 900, 1200, 4000};
  double RsqBins[] = {0.15, 0.2, 0.25, 0.30, 0.41, 0.52, 1.5};

  const int NMRBins = sizeof(MRBins)/sizeof(double)-1;
  const int NRsqBins = sizeof(RsqBins)/sizeof(double)-1;

  cout<<"Number of MR bins: "<<NMRBins<<" , Number of Rsq bins: "<<NRsqBins<<endl;

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

  vector<TH1D*> histSigmaIetaIeta_EB;
  vector<TH1D*> histSigmaIetaIetaTemplate_EB;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histSigmaIetaIeta_EB.push_back(new TH1D(Form("histSigmaIetaIeta_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", 100, 0, 0.05));
    histSigmaIetaIetaTemplate_EB.push_back(new TH1D(Form("histSigmaIetaIetaTemplate_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", 100, 0, 0.05));

    histSigmaIetaIeta_EB[i]->Sumw2();
    histSigmaIetaIetaTemplate_EB[i]->Sumw2();
  }
 
  TH1D *histSigmaIetaIeta_EB_binned[NMRBins][NRsqBins];
  TH1D *histSigmaIetaIeta_FakesTemplate_EB_binned[NMRBins];
  TH1D *histSigmaIetaIeta_PromptTemplate_EB_binned[NMRBins];

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histSigmaIetaIeta_EB_binned[a][b] = new TH1D (Form("histSigmaIetaIeta_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", 100, 0, 0.05);
      histSigmaIetaIeta_EB_binned[a][b] -> Sumw2();
    }
  }

  for (int a=0; a<NMRBins; a++) {
    histSigmaIetaIeta_FakesTemplate_EB_binned[a] = new TH1D (Form("histSigmaIetaIeta_FakesTemplate_EB_binned_%d", a), "; Sigma_ietaieta; Number of Events", 100, 0, 0.05);
    histSigmaIetaIeta_FakesTemplate_EB_binned[a] -> Sumw2();
 
    histSigmaIetaIeta_PromptTemplate_EB_binned[a] = new TH1D (Form("histSigmaIetaIeta_PromptTemplate_EB_binned_%d", a), "; Sigma_ietaieta; Number of Events", 100, 0, 0.05);
    histSigmaIetaIeta_PromptTemplate_EB_binned[a] -> Sumw2();
 }

  TH2D *histSigmaIetaIeta_EB_MRRsq= new TH2D("histSigmaIetaIeta_EB_MRRsq", "; MR; Rsq", NMRBins, MRBins, NRsqBins, RsqBins);
  TH1D *histweights= new TH1D("histweights", "; weights", 1000, -2., 10000.);

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
      // for(UInt_t ientry=0; ientry < 1000000; ientry++) {       	
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
  
	if (isData) {
	  if (passTrigger) {
	    
	    if( fabs(events->pho1.Eta()) < 1.44 && events->pho1_chargediso > 2.5 && events->pho1_sigmaietaieta < 0.015 ){
	      histSigmaIetaIetaTemplate_EB[i]->Fill(events->pho1_sigmaietaieta);
	      
	      histweights->Fill(weight);
	      
	      // create binned fakes template
	      for(int ii = 0; ii<NMRBins; ii++)
		{
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] )
		    histSigmaIetaIeta_FakesTemplate_EB_binned[ii]->Fill(events->pho1_sigmaietaieta);
		}
	    }

	    // select photons in EB
	    // if(fabs(events->pho1.Eta()) < 1.44 && events->pho1.Pt()>190 && events->pho1_sigmaietaieta < 0.015 && events->pho1_chargediso < 2.5) {
	    if(fabs(events->pho1.Eta()) < 1.44 && events->pho1_sigmaietaieta < 0.015 && events->pho1_chargediso < 2.5) {
	      histSigmaIetaIeta_EB[i]->Fill(events->pho1_sigmaietaieta);
	      
	      for(int ii = 0; ii<NMRBins; ii++)
		for(int jj = 0; jj<NRsqBins; jj++)
		  {
		    if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] )
		      if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
			histSigmaIetaIeta_EB_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
			// cout<<MRBins[ii]<<" "<<MRBins[ii+1]<<" "<<RsqBins[jj]<< " "<<RsqBins[jj+1]<<endl;
		      }
		  }
	    }
	  }
	}
	else
	  {
	    if( fabs(events->pho1.Eta()) < 1.44 && events->pho1_chargediso < 2.5 ) {
	      histSigmaIetaIetaTemplate_EB[i]->Fill(events->pho1_sigmaietaieta, weight);
	      
	      // create binned prompt template
	      for(int ii = 0; ii<NMRBins; ii++)
		{
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] )
		    histSigmaIetaIeta_PromptTemplate_EB_binned[ii]->Fill(events->pho1_sigmaietaieta, weight);
		}
	    }
	  }
	//By default it's EG Loose selection + chargedIso < 2.5 GeV.
	//We can tighten the cuts if we wish...

	if (option == "MR300Rsq0p15" ) {
	  if (!(events->NJets80 >= 2 && events->MR_NoPho > 300 && events->Rsq_NoPho > 0.15 )) continue;	
	}
      
	//MET Filters
	if (!(events->Flag_HBHENoiseFilter && events->Flag_goodVertices && events->Flag_eeBadScFilter)) continue;

      } //loop over events
    } //loop over input files
  } //loop over input file groups

  // TH1F *data;                              //data histogram
  // TH1F *mc0;                               // first MC histogram
  // TH1F *mc1;                               // second MC histogram
  // TH1F *mc2;                               // third MC histogram
  
  // cout<<"before scaling: "<<histSigmaIetaIetaTemplate_EB[0]->Integral()<<endl;
  // histSigmaIetaIetaTemplate_EB[0]->Scale(histSigmaIetaIeta_EB[0]->Integral()/histSigmaIetaIetaTemplate_EB[0]->Integral());
  // cout<<"after scaling: "<<histSigmaIetaIetaTemplate_EB[0]->Integral()<<endl;

  TObjArray *mc = new TObjArray(2);        // template histograms are put in this array
  mc->Add(histSigmaIetaIetaTemplate_EB[0]); //fakes
  // mc->Add(fakeHist); //fakes
  mc->Add(promptHist); //prompt
  
  TFractionFitter* fit = new TFractionFitter(histSigmaIetaIeta_EB[0], mc); // initialise
  fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  fit->SetRangeX(1,50);                 // use only the first 15 bins in the fit

   Double_t p0, p1, errP0, errP1;
   Int_t status = fit->Fit();               // perform the fit
   std::cout << "fit status: " << status << std::endl;
   if (status == 0) {                       // check on fit status
     TCanvas * c1 = new TCanvas("c", "c", 600, 600) ;
     c1->cd();
     TPad pad1("pad1","pad1",0,0.4,1,1);
     pad1.SetBottomMargin(0);
     pad1.SetLogy();
     pad1.Draw();
     pad1.cd();
     
     TH1F* result = (TH1F*) fit->GetPlot();
     result->SetLineColor(2);
     
     fit->GetResult( 0, p0, errP0);
     cout<<"Fakes Fraction: "<< p0 <<"+/-"<<errP0<<endl;
     fit->GetResult( 1, p1, errP1);
     cout<<"Prompt Fraction: "<< p1 <<"+/-"<<errP1<<endl;

     histSigmaIetaIeta_EB[0]->Draw("Ep");
     result->Draw("samehist");

     TH1F* resultfake = (TH1F*) fit->GetMCPrediction(0);

     resultfake->Scale((p0*histSigmaIetaIeta_EB[0]->Integral())/resultfake->Integral()); // normalize the Template of the fakes
     resultfake->SetLineColor(3);

     TH1F* resultprompt = (TH1F*) fit->GetMCPrediction(1);
     resultprompt->Scale((p1*histSigmaIetaIeta_EB[0]->Integral())/resultprompt->Integral()); // normalize the Template of the prompts
     resultprompt->SetLineColor(4);

     Double_t TotalFakes  = resultfake->Integral(0, resultfake->GetXaxis()->FindBin(0.0105));
     Double_t TotalPrompt = resultprompt->Integral(0, resultprompt->GetXaxis()->FindBin(0.0105));
     Double_t Total       = result->Integral(0, result->GetXaxis()->FindBin(0.0105));

     cout<<"Fake amount passing sigma_ietaieta cut: "<<TotalFakes<<endl;
     cout<<"Prompt amount passing sigma_ietaieta cut: "<<TotalPrompt<<endl;
     cout<<"Total amount passing sigma_ietaieta cut: "<<Total<<endl;
     
     resultfake->Draw("samehist");
     resultprompt->Draw("samehist");
 
     pad1.Modified();
     gPad->Update();
     
     c1->cd();

     TH1F *resultcopy = (TH1F*)result->Clone();
     resultcopy -> Divide(histSigmaIetaIeta_EB[0]);
     resultcopy->GetYaxis()->SetTitle("Total Fit/Data");
     resultcopy->SetMinimum(0.2);
     resultcopy->SetMaximum(2.0);
     
     TPad pad2("pad2","pad2",0,0.0,1,0.4);
     pad2.SetTopMargin(0);
     pad2.SetTopMargin(0.008);
     pad2.SetBottomMargin(0.25);
     pad2.SetGridy();
     pad2.Draw();
     pad2.cd();
     resultcopy->Draw("pe");
     pad2.Modified();
     gPad->Update();

     c1->SetLogy();
     c1->SaveAs("Fit.png");
  }


   // now do the fits in the MR/Rsq bins   
   for(int ii=0; ii<NMRBins; ii++)
     for(int jj=0; jj<NRsqBins; jj++) {
       
       cout<<"BINNED FIT: "<<ii<<" "<<jj<<endl;
       
       // TObjArray *mc_binned = new TObjArray(2);        // template histograms are put in this array
       // mc_binned->Add(histSigmaIetaIeta_FakesTemplate_EB_binned[ii]); //fakes
       // mc_binned->Add(histSigmaIetaIeta_PromptTemplate_EB_binned[ii]); //prompt
       // mc_binned->Add(promptHist); //prompt
       
       TFractionFitter* fit_binned = new TFractionFitter(histSigmaIetaIeta_EB_binned[ii][jj], mc); // initialise
       
       Double_t p0_binned, p1_binned, errP0_binned, errP1_binned;

       Int_t status_binned = fit_binned->Fit();               // perform the fit

       std::cout << "fit status binned: "<<status_binned<<" "<< ii<<" "<< jj << std::endl;
       if (status_binned == 0) {                       // check on fit status
   	 TCanvas c("c", "c", 800, 600);
   	 c.Clear();
   	 c.cd();
   	 TPad pad1("pad1","pad1",0,0.4,1,1);
   	 pad1.SetBottomMargin(0);
   	 pad1.SetLogy();
   	 pad1.Draw();
   	 pad1.cd();
	 
   	 TH1F* result = (TH1F*) fit_binned->GetPlot();
	 result->Scale((histSigmaIetaIeta_EB_binned[ii][jj]->Integral())/result->Integral());
   	 result->SetLineColor(2);
	 
   	 fit_binned->GetResult( 0, p0_binned, errP0_binned);
   	 cout<<"Fakes Fraction: "<< ii<<" "<< jj << " "<<p0_binned <<"+/-"<<errP0_binned<<endl;
   	 fit_binned->GetResult( 1, p1_binned, errP1_binned);
   	 cout<<"Prompt Fraction: "<< ii<<" "<< jj << " "<< p1_binned <<"+/-"<<errP1_binned<<endl;
	 
   	 histSigmaIetaIeta_EB_binned[ii][jj]->Draw("Ep");
   	 result->Draw("samehist");
	 
   	 TH1F* resultfake = (TH1F*) fit_binned->GetMCPrediction(0);
	 
   	 resultfake->Scale((p0_binned*histSigmaIetaIeta_EB_binned[ii][jj]->Integral())/resultfake->Integral()); // normalize the Template of the fakes
   	 resultfake->SetLineColor(3);
	 
   	 TH1F* resultprompt = (TH1F*) fit_binned->GetMCPrediction(1);
   	 resultprompt->Scale((p1_binned*histSigmaIetaIeta_EB_binned[ii][jj]->Integral())/resultprompt->Integral()); // normalize the Template of the prompts
   	 resultprompt->SetLineColor(4);
	 
   	 Double_t TotalFakes  = resultfake->Integral(0, resultfake->GetXaxis()->FindBin(0.0105));
   	 Double_t TotalPrompt = resultprompt->Integral(0, resultprompt->GetXaxis()->FindBin(0.0105));
   	 Double_t Total       = result->Integral(0, result->GetXaxis()->FindBin(0.0105));
	 
   	 cout<<"Fake amount passing sigma_ietaieta cut: "<<TotalFakes<<endl;
   	 cout<<"Prompt amount passing sigma_ietaieta cut: "<<TotalPrompt<<endl;
   	 cout<<"Total amount passing sigma_ietaieta cut: "<<Total<<endl;
	 
   	 histSigmaIetaIeta_EB_MRRsq->SetBinContent(ii+1, jj+1, TotalFakes/Total);
	 
   	 resultfake->Draw("samehist");
   	 resultprompt->Draw("samehist");
   	 pad1.Modified();
   	 gPad->Update();

   	 c.cd();

   	 TH1F *resultcopy = (TH1F*)result->Clone();
   	 resultcopy -> Divide(histSigmaIetaIeta_EB_binned[ii][jj]);
   	 resultcopy->GetYaxis()->SetTitle("Total Fit/Data");
   	 resultcopy->SetMinimum(0.2);
   	 resultcopy->SetMaximum(2.0);

   	 TPad pad2("pad2","pad2",0,0.0,1,0.4);
   	 pad2.SetTopMargin(0);
   	 pad2.SetTopMargin(0.008);
   	 pad2.SetBottomMargin(0.25);
   	 pad2.SetGridy();
   	 pad2.Draw();
   	 pad2.cd();
   	 resultcopy->Draw("pe");
   	 pad2.Modified();
   	 gPad->Update();

   	 c.SaveAs(Form("Fit_%d_%d.png", ii, jj));
       }
     }   

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("PhotonControlRegionPlots"+Label+".root").c_str(), "RECREATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histSigmaIetaIeta_EB[i], Form("histSigmaIetaIeta_EB_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIetaTemplate_EB[i], Form("histSigmaIetaIetaTemplate_EB_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EB_MRRsq, "histSigmaIetaIeta_EB_MRRsq", "WriteDelete");
    file->WriteTObject(histweights, "histweights", "WriteDelete");
    
    for(int ii=0; ii<NMRBins; ii++) {
      for(int jj=0; jj<NRsqBins; jj++) {
	file->WriteTObject(histSigmaIetaIeta_EB_binned[ii][jj], Form("histSigmaIetaIeta_EB_binned_%s_%d_%d",processLabels[i].c_str(), ii, jj), "WriteDelete");
      }
      file->WriteTObject(histSigmaIetaIeta_FakesTemplate_EB_binned[ii], Form("histSigmaIetaIeta_FakesTemplate_EB_binned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
      file->WriteTObject(histSigmaIetaIeta_PromptTemplate_EB_binned[ii], Form("histSigmaIetaIeta_PromptTemplate_EB_binned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
    }
  }

  file->Close();
  delete file;       
}

void PhotonControlSample_Fit( int option = 10) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";


  //No Skims  
  if (option >= 10) {
    datafiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_PRv4_GoodLumiGolden.root");
    //  
  } else {
    datafiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_Run2015D_GoodLumiGolden.root");     
  }

  vector<string> bkgfiles_gjets;
  // vector<string> bkgfiles_qcd;
  // vector<string> bkgfiles_other;

  if (option >= 10) {
    // bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");    
    // bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");    
    // bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");    
    // bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");    
    // bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");    
  // } else {
  //   bkgfiles_gjets.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_RazorSkim.root");    
  //   bkgfiles_qcd.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_QCD_HTBinned_1pb_weighted_RazorSkim.root"); 
  //   bkgfiles_other.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final/RazorSkim/RunTwoRazorControlRegions_PhotonFull_Other_1pb_weighted_RazorSkim.root"); 
  }
   

  bkgfiles.push_back(bkgfiles_gjets);
  // bkgfiles.push_back(bkgfiles_qcd);
  // bkgfiles.push_back(bkgfiles_other);

  processLabels.push_back("GJets");  
  // processLabels.push_back("QCD");
  // processLabels.push_back("Other");

  colors.push_back(kOrange);
  // colors.push_back(kMagenta);
  // colors.push_back(kCyan);
  
  double lumi = 2185;

  //*********************************************************************
  //GJets Control Region
  //*********************************************************************
   if (option == 0) {
     RunPhotonControlSample_Fit(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",0,"MR300Rsq0p15");
   }
   if (option == 10) {
     RunPhotonControlSample_Fit(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",0,"Inclusive");
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


/* The template for prompts was made with this snippet on command line
TChain * chain = new TChain("ControlSampleEvent");
 chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
h1 = new TH1F("h1","h1",100,0.,0.05);
h1->Sumw2();

chain->Draw("pho1_sigmaietaieta>>h1","weight*((HLTDecision[88]==1 || HLTDecision[89]==1 || HLTDecision[90]==1 || HLTDecision[91] ==1 ||HLTDecision[92]==1 || HLTDecision[93]==1)&&fabs(pho1.Eta())<1.44&&pho1_chargediso<2.5&&pho1_sigmaietaieta<0.015)");
h1->SaveAs("a.root")
*/
