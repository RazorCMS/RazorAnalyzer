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
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
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



void MakePileupReweight() {



  TFile *pileupTargetFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupTarget2015.root", "READ");
  TH1F *pileupTargetHist = (TH1F*)pileupTargetFile->Get("pileup");
  assert(pileupTargetHist);

  TFile *pileupSourceFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupSource_Spring15_25ns.root", "READ");
  TH1F *pileupSourceHist = (TH1F*)pileupSourceFile->Get("PileupSourceHist");
  assert(pileupSourceHist);


  //*******************************************************************************************
  //Make NVtx Reweighting Function
  //*******************************************************************************************
  TH1F *PileupTargetNormalized = NormalizeHist( pileupTargetHist );
  TH1F *PileupSourceNormalized = NormalizeHist( pileupSourceHist );

  TH1F *PileupReweight = new TH1F ("PileupReweight",";NPU;Weight",50,-0.5,49.5);

  for (int i=1; i<PileupReweight->GetXaxis()->GetNbins()+1; i++) {

    double data = 0;
    double bkg = 0;
    if (PileupSourceNormalized->GetBinContent(i) > 0) {
      PileupReweight->SetBinContent(i,PileupTargetNormalized->GetBinContent(i)/PileupSourceNormalized->GetBinContent(i));
    } else if (PileupTargetNormalized->GetBinContent(i) == 0){
      PileupReweight->SetBinContent(i,0.0);
    } else {
      if (i == 1) {
	PileupReweight->SetBinContent(i,1);
      } else {
	PileupReweight->SetBinContent(i,PileupReweight->GetBinContent(i-1));
      }
    }

    cout << "Bin " << i << " : " << PileupReweight->GetBinCenter(i) << " : " << PileupTargetNormalized->GetBinContent(i) << " / " << PileupSourceNormalized->GetBinContent(i) << " = " << PileupReweight->GetBinContent(i) << "\n";
  }



  TFile *file = TFile::Open("PileupReweight2015.root", "UPDATE");
  file->cd();
  file->WriteTObject(PileupReweight, "PileupReweight", "WriteDelete");
  file->Close();
  delete file;


 

}





