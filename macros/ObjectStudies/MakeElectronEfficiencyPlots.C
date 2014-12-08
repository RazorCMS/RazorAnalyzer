//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
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

#include "RazorAnalyzer/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/ElectronTree.h"

#endif

Bool_t passCSA14VetoID( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.02
	   && fabs(eleTree->fEleDPhiIn) < 0.2579
	   && eleTree->fEleSigmaIEtaIEta < 0.0125
	   && eleTree->fEleHoverE < 0.2564
	   && fabs(eleTree->fEleD0) < 0.025
	   && fabs(eleTree->fEleDZ) < 0.5863
	   && eleTree->fEleOneOverEMinusOneOverP < 0.1508	   
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.0141
	   && fabs(eleTree->fEleDPhiIn) < 0.2591
	   && eleTree->fEleSigmaIEtaIEta < 0.0371 
	   && eleTree->fEleHoverE < 0.1335 
	  && fabs(eleTree->fEleD0) < 0.2232
	   && fabs(eleTree->fEleDZ) <  0.9513
	   && eleTree->fEleOneOverEMinusOneOverP <  0.1542
	  && eleTree->fElePassConversion
	  && eleTree->fEleNMissHits < 3
	  ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passCSA14VetoIso( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( eleTree->fElePFIso04 < 0.3313   
	   ) {
	pass = true;
      }
    } else {
      if (  eleTree->fElePFIso04 < 0.3816
	  ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passImprovedIso( ElectronTree *eleTree) {

    bool pass = false;
    if(eleTree->fElePt > 20) {
      pass = passCSA14VetoIso(eleTree);
    } else {
      pass = bool(eleTree->fElePFIso04*eleTree->fElePt < 8);
    }
    return pass;
}

Bool_t passIDMVANonTrigVeto( ElectronTree *eleTree) {

  Int_t subdet = 0;
  if (fabs(eleTree->fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(eleTree->fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (eleTree->fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.0;
  if (MVABin == 1) MVACut = 0.0;
  if (MVABin == 2) MVACut = 0.3; 
  if (MVABin == 3) MVACut = -0.3;
  if (MVABin == 4) MVACut = -0.3;
  if (MVABin == 5) MVACut = -0.3;  

    bool pass = false;
    if (eleTree->fIDMVANonTrig > MVACut 
	) {
      pass = true;
    }   
    return pass;
}




bool PassSelection( ElectronTree* eleTree ) {

  bool pass = false;

  if ( eleTree->fElePt > 0 && 
       //passCSA14VetoID(eleTree) &&
       passIDMVANonTrigVeto(eleTree) && 

       //(0 == 0)
       //(eleTree->fElePFIso04*eleTree->fElePt < 8)
       //passCSA14VetoIso(eleTree)
       passImprovedIso(eleTree)
       ) {		  	   
    pass = true;
  }
  return pass;
  
}



//=== MAIN MACRO ================================================================================================= 

void MakeElectronEfficiencyPlots(const string inputfile, int option = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Electron Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Electron Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 50, 0 , 200, 50, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  ElectronTree *EleTree = new ElectronTree;
  EleTree->LoadTree(inputfile.c_str());
  EleTree->InitTree(ElectronTree::kEleTreeLight);

  cout << "Total Entries: " << EleTree->tree_->GetEntries() << "\n";
  int nentries = EleTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < EleTree->tree_->GetEntries(); ientry++) {       	
    EleTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (EleTree->fEleGenPt < 5) continue;
    if (abs(EleTree->fEleGenEta) > 2.4) continue;



    //**** PT - ETA ****
    histDenominatorPtEta->Fill(EleTree->fEleGenPt,EleTree->fEleGenEta);
    if(PassSelection(EleTree)) {
      histNumeratorPtEta->Fill(EleTree->fEleGenPt,EleTree->fEleGenEta);
    }


    //**** PT ****
      histDenominatorPt->Fill(EleTree->fEleGenPt);

      //Numerator
      if(PassSelection(EleTree)) {
        histNumeratorPt->Fill(EleTree->fEleGenPt);        
      }


    //**** Eta ****
    if (fabs(EleTree->fEleGenPt) > 30) {
      histDenominatorEta->Fill(EleTree->fEleGenEta);

      //Numerator
      if(PassSelection(EleTree)) {
        histNumeratorEta->Fill(EleTree->fEleGenEta);        
      }

    }

    //**** Phi ****
    if (fabs(EleTree->fEleGenEta) < 2.4) {
      histDenominatorPhi->Fill(EleTree->fEleGenPhi);

      //Numerator
      if(PassSelection(EleTree)) {
        histNumeratorPhi->Fill(EleTree->fEleGenPhi);        
      }

    }

    //**** Rho ****
    if (fabs(EleTree->fEleGenEta) < 2.4) {
      histDenominatorRho->Fill(EleTree->fRho);

      //Numerator
      if(PassSelection(EleTree)) {
        histNumeratorRho->Fill(EleTree->fRho);        
      }

    }
    //**** Npv ****
    if (fabs(EleTree->fEleGenEta) < 2.4) {
      histDenominatorNpv->Fill(EleTree->fNVertices);

      //Numerator
      if(PassSelection(EleTree)) {
        histNumeratorNpv->Fill(EleTree->fNVertices);        
      }

    }

    // //**** Npu ****
    // if (fabs(EleTree->fEleGenEta) < 2.4) {
    //   histDenominatorNpu->Fill(EleTree->);

    //   //Numerator
    //   if(PassSelection(EleTree)) {
    //     histNumeratorNpu->Fill(EleTree->);        
    //   }

    // }


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
