//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakePhotonEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/PhotonNtuple/PhotonNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Photon")'
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
#include "RazorAnalyzer/include/PhotonTree.h"

#endif


bool PassSelection( PhotonTree* PhoTree , int wp = 0 ) {

  bool pass = false;

  //**********************************
  //Tight Selection
  //**********************************
  if (wp == 2 && PhoTree->fPhoIsTight ) {
    pass = true;
  }

  //**********************************
  //Medium Selection
  //**********************************
  if (wp == 1 && PhoTree->fPhoIsMedium ) {
    pass = true;
  }

  //**********************************
  //Loose Selection
  //**********************************
  if (wp == 0 && PhoTree->fPhoIsLoose ) {
    pass = true;
  }
 

  return pass;

}


TGraphAsymmErrors* getEffGraph( string filename, string graphname) {
  
  TFile *f = new TFile(filename.c_str(),"READ");
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)f->Get(graphname.c_str());
  f->Close();
  delete f;
  if (!graph) {
    cout << "Graph " << graphname << " from " << filename << " could not be retrieved\n";
    assert(false);
  }
  return graph;
}


void plotPhotonEfficiency() {

  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Loose.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Medium.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Tight.root","Efficiency_Pt");

  TGraphAsymmErrors* effVsPt_ttH_25ns_Loose = getEffGraph("Efficiency_ttH_25ns_Spring15_Loose.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_25ns_Medium = getEffGraph("Efficiency_ttH_25ns_Spring15_Medium.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_25ns_Tight = getEffGraph("Efficiency_ttH_25ns_Spring15_Tight.root","Efficiency_Pt");

  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Loose.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Medium.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Tight.root","Efficiency_Eta");

  TGraphAsymmErrors* effVsEta_ttH_25ns_Loose = getEffGraph("Efficiency_ttH_25ns_Spring15_Loose.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_ttH_25ns_Medium = getEffGraph("Efficiency_ttH_25ns_Spring15_Medium.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_ttH_25ns_Tight = getEffGraph("Efficiency_ttH_25ns_Spring15_Tight.root","Efficiency_Eta");

  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Loose.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Medium.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Tight.root","Efficiency_NPV");

  TGraphAsymmErrors* effVsNPV_ttH_25ns_Loose = getEffGraph("Efficiency_ttH_25ns_Spring15_Loose.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_ttH_25ns_Medium = getEffGraph("Efficiency_ttH_25ns_Spring15_Medium.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_ttH_25ns_Tight = getEffGraph("Efficiency_ttH_25ns_Spring15_Tight.root","Efficiency_NPV");


  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_NotCloseToParton_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Spring15_Loose_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_25ns_NotCloseToParton_Loose = getEffGraph("Efficiency_ttH_25ns_Spring15_Loose_NotCloseToParton.root","Efficiency_Pt");




   //****************************************************************************
  //Make Plots
  //****************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //****************************************************************************
  //POG WP Efficiencies
  //****************************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsPt_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsPt_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsPt_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsPt_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);

  effVsPt_GJetFlat_50ns_Loose->Draw("AP");
  effVsPt_GJetFlat_50ns_Medium->Draw("Psame");
  effVsPt_GJetFlat_50ns_Tight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt.gif");
  cv->SaveAs("PhotonEfficiencyVsPt.pdf");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsEta_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsEta_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsEta_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsEta_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsEta_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon #eta");
  effVsEta_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsEta_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsEta_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsEta_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsEta_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);

  effVsEta_GJetFlat_50ns_Loose->Draw("AP");
  effVsEta_GJetFlat_50ns_Medium->Draw("Psame");
  effVsEta_GJetFlat_50ns_Tight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsEta.gif");
  cv->SaveAs("PhotonEfficiencyVsEta.pdf");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsNPV_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsNPV_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effVsNPV_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsNPV_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsNPV_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsNPV_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsNPV_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);


  effVsNPV_GJetFlat_50ns_Loose->Draw("AP");
  effVsNPV_GJetFlat_50ns_Medium->Draw("Psame");
  effVsNPV_GJetFlat_50ns_Tight->Draw("Psame");
  effVsNPV_GJetFlat_50ns_Loose->GetXaxis()->SetRangeUser(0,40);
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsNPV.gif");
  cv->SaveAs("PhotonEfficiencyVsNPV.pdf");


  //****************************************************************************
  //Compare GJet with ttH
  //****************************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_Loose, "#gamma+Jet Flat Loose WP", "LP");
  legend->AddEntry(effVsPt_ttH_25ns_Loose, "ttH#rightarrow#gamma#gamma Loose WP", "LP");

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsPt_ttH_25ns_Loose->SetLineWidth(3);
  effVsPt_ttH_25ns_Loose->SetLineColor(kBlue);

  effVsPt_GJetFlat_50ns_Loose->Draw("AP");
  effVsPt_ttH_25ns_Loose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_Loose.gif");
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_Loose.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_NotCloseToParton_Loose, "#gamma+Jet Flat Loose WP", "LP");
  legend->AddEntry(effVsPt_ttH_25ns_NotCloseToParton_Loose, "ttH#rightarrow#gamma#gamma Loose WP", "LP");

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineColor(kRed);
  effVsPt_ttH_25ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_ttH_25ns_NotCloseToParton_Loose->SetLineColor(kBlue);

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->Draw("AP");
  effVsPt_ttH_25ns_NotCloseToParton_Loose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_NotCloseToParton_Loose.gif");
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_NotCloseToParton_Loose.pdf");




 return;







}




//=== MAIN MACRO ================================================================================================= 

void ProducePhotonEfficiencyPlots(const string inputfile, int wp = 0, int option = -1, bool usePhotonNotNearParton = false, string label = "") {

   // plotPhotonEfficiency();
   // return;

  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Photon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Photon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Photon #eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Photon #eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Photon #phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Photon #phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Number of Reconstructed Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Number of Reconstructed Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Number of Pileup Interactions; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Number of Pileup Interactions; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  PhotonTree *PhoTree = new PhotonTree;
  PhoTree->LoadTree(inputfile.c_str());
  PhoTree->InitTree(PhotonTree::kPhotonTreeLight);

  cout << "Total Entries: " << PhoTree->tree_->GetEntries() << "\n";
  int nentries = PhoTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < PhoTree->tree_->GetEntries(); ientry++) {       	
    PhoTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (PhoTree->fPhoGenPt < 5) continue;
    if (abs(PhoTree->fPhoGenEta) > 2.5) continue;

    if (!(PhoTree->fPhoPt > 20)) continue;

    if (usePhotonNotNearParton) {
      if (PhoTree->fDRToClosestParton < 1.0) continue;
    }


    if (option==0) {

      if (PhoTree->fPhoGenPt < 20) continue;      
      

      //**** PT - ETA ****
      histDenominatorPtEta->Fill(PhoTree->fPhoGenPt,PhoTree->fPhoGenEta);
      if(PassSelection(PhoTree, wp)) {
	histNumeratorPtEta->Fill(PhoTree->fPhoGenPt,PhoTree->fPhoGenEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(PhoTree->fPhoGenPt);

      //Numerator
      if(PassSelection(PhoTree, wp)) {
        histNumeratorPt->Fill(PhoTree->fPhoGenPt);        
      }


      //**** Eta ****
      if (fabs(PhoTree->fPhoGenPt) > 30) {
	histDenominatorEta->Fill(PhoTree->fPhoGenEta);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorEta->Fill(PhoTree->fPhoGenEta);        
	}

      }

      //**** Phi ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorPhi->Fill(PhoTree->fPhoGenPhi);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorPhi->Fill(PhoTree->fPhoGenPhi);        
	}

      }

      //**** Rho ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorRho->Fill(PhoTree->fRho);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorRho->Fill(PhoTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorNpv->Fill(PhoTree->fNVertices);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorNpv->Fill(PhoTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(PhoTree->fPhoGenEta) < 2.5) {
      //   histDenominatorNpu->Fill(PhoTree->);

      //   //Numerator
      //   if(PassSelection(PhoTree, wp)) {
      //     histNumeratorNpu->Fill(PhoTree->);        
      //   }

      // }

    } //end if option  == 0



    if (option==1) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(PhoTree->fPhoPt,PhoTree->fPhoEta);
      if(PassSelection(PhoTree, wp)) {
	histNumeratorPtEta->Fill(PhoTree->fPhoPt,PhoTree->fPhoEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(PhoTree->fPhoPt);

      //Numerator
      if(PassSelection(PhoTree, wp)) {
        histNumeratorPt->Fill(PhoTree->fPhoPt);        
      }


      //**** Eta ****
      if (fabs(PhoTree->fPhoPt) > 30) {
	histDenominatorEta->Fill(PhoTree->fPhoEta);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorEta->Fill(PhoTree->fPhoEta);        
	}

      }

      //**** Phi ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorPhi->Fill(PhoTree->fPhoPhi);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorPhi->Fill(PhoTree->fPhoPhi);        
	}

      }

      //**** Rho ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorRho->Fill(PhoTree->fRho);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorRho->Fill(PhoTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorNpv->Fill(PhoTree->fNVertices);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorNpv->Fill(PhoTree->fNVertices);        
	}

      }

    
    }



  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "Efficiency_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "Efficiency_Npu" , vector<double>() ,  -99, -99, 0, 1);  
  TH2F *efficiency_pteta = createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Pt.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Eta.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Phi.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Rho.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npv.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npu.gif").c_str());


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");
  file->WriteTObject(efficiency_eta, "Efficiency_Eta", "WriteDelete");
  file->WriteTObject(efficiency_phi, "Efficiency_Phi", "WriteDelete");
  file->WriteTObject(efficiency_rho, "Efficiency_Rho", "WriteDelete");
  file->WriteTObject(efficiency_npv, "Efficiency_NPV", "WriteDelete");
  file->WriteTObject(efficiency_npu, "Efficiency_NPU", "WriteDelete");
  file->WriteTObject(efficiency_pteta, "Efficiency_PtEta", "WriteDelete");

  file->Close();
  delete file;       

}




void MakePhotonEfficiencyPlots( int option = 0) {

  if (option == 1) {
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 0, 0, false, "ttH_25ns_Spring15_Loose");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 1, 0, false, "ttH_25ns_Spring15_Medium");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 2, 0, false, "ttH_25ns_Spring15_Tight");

    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 0, 0, true, "ttH_25ns_Spring15_Loose_NotCloseToParton");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 1, 0, true, "ttH_25ns_Spring15_Medium_NotCloseToParton");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_Spring15_25ns.root", 2, 0, true, "ttH_25ns_Spring15_Tight_NotCloseToParton");

    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 0, 0, false, "GJetFlat_50ns_Spring15_Loose");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 1, 0, false, "GJetFlat_50ns_Spring15_Medium");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 2, 0, false, "GJetFlat_50ns_Spring15_Tight");

    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 0, 0, true, "GJetFlat_50ns_Spring15_Loose_NotCloseToParton");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 1, 0, true, "GJetFlat_50ns_Spring15_Medium_NotCloseToParton");
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_PHYS14_50ns.root", 2, 0, true, "GJetFlat_50ns_Spring15_Tight_NotCloseToParton");
  }	



  plotPhotonEfficiency();

}
