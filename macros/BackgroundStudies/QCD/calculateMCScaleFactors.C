//--------------------------------------------------------------
//
// make ratio plots for Razor sideband and dijet control region
//
//--------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TEventList.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2Poly.h>
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include <CalStyleRemix.hh>

//enum {WJet=0, ZJet, TTJet};
void processFile(TString inputFileName, TString outputFileName, int physProc);
void calculateMCScaleFactors() {

  //processFile("/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/clean/CMSSW_7_1_5/src/RazorAnalyzer/Backgrounds/Signal/FullRazorInclusive_Razor2016_MoriondRereco_WJets_1pb_weighted.root", "WJets.root", 0);

  processFile("/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/clean/CMSSW_7_1_5/src/RazorAnalyzer/Backgrounds/Signal/FullRazorInclusive_Razor2016_MoriondRereco_ZInv_1pb_weighted.root", "ZInv.root", 1);

  //processFile("/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/clean/CMSSW_7_1_5/src/RazorAnalyzer/Backgrounds/Signal/FullRazorInclusive_Razor2016_MoriondRereco_TTJetsInclusive_1pb_weighted.root", "TTJets.root", 2);


	      
}

void processFile(TString inputFileName, TString outputFileName, int physProc) {

  TFile *inputFile = TFile::Open(inputFileName,"read");
  TTree *inputTree = (TTree*) inputFile->Get("RazorInclusive");

  float MR,Rsq,dPhiRazor,leadingJetPt;
  int nSelectedJets,nBTaggedJets, NISRJets;
  float met, metOverCaloMet;
  int box;
  float weight;

  inputTree->SetBranchAddress("MR",             &MR);
  inputTree->SetBranchAddress("Rsq",            &Rsq);
  inputTree->SetBranchAddress("dPhiRazor",      &dPhiRazor);
  inputTree->SetBranchAddress("met",            &met);
  inputTree->SetBranchAddress("weight",         &weight);
  inputTree->SetBranchAddress("box",            &box);

  inputTree->SetBranchAddress("leadingJetPt",   &leadingJetPt);
  inputTree->SetBranchAddress("nSelectedJets",  &nSelectedJets);
  inputTree->SetBranchAddress("nBTaggedJets",   &nBTaggedJets);
  inputTree->SetBranchAddress("NISRJets",       &NISRJets);
  inputTree->SetBranchAddress("metOverCaloMet", &metOverCaloMet);

  inputTree->GetEntry(0);

  TString cut="(box==11||box==12||box==14)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)";

  TFile *fScaleFactors = TFile::Open("RazorScaleFactors_Razor2016_MoriondRereco.root");
  TFile *fNJetScaleFactors = TFile::Open("RazorNJetsScaleFactors_Razor2016_MoriondRereco.root");

  TH2Poly *hScaleFactors=0; 
  TH1F *hNJetScaleFactors=0; 

  if (physProc==0) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("WJetsInvScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("WJetsScaleFactors");
  }
  else if (physProc==1) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("WJetsInvScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("WJetsInvScaleFactors");
  }
  else if (physProc==2) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("TTJetsScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("TTJetsScaleFactors");
  }

  TFile *outputFile = new TFile(outputFileName,"recreate");
  TTree *outputTree = (TTree*) inputTree->GetTree()->CloneTree(0);

  float mcScaleFactor;
  outputTree->Branch("mcScaleFactor",&mcScaleFactor, "mcScaleFactor/F");

  inputTree->Draw(">>elist1",cut);
  
  TEventList *list = (TEventList*)gDirectory->Get("elist1");

  for (Int_t i=0; i<list->GetN(); i++) {
  //for (UInt_t i=0; i<10; i++) {
    inputTree->GetEntry(list->GetEntry(i));

    double tNJets=min((double)nSelectedJets, hNJetScaleFactors->GetXaxis()->GetXmax()*0.999);
    tNJets=max(tNJets, hNJetScaleFactors->GetXaxis()->GetXmin()*1.001);

    double tMR=min((double)MR, hScaleFactors->GetXaxis()->GetXmax()*0.999);
    tMR=max(tMR, hScaleFactors->GetXaxis()->GetXmin()*1.001);

    double tRsq=min((double)Rsq, hScaleFactors->GetYaxis()->GetXmax()*0.999);
    tRsq=max(tRsq, hScaleFactors->GetYaxis()->GetXmin()*1.001);

    double scaleFactor = hScaleFactors->GetBinContent(hScaleFactors->FindBin(tMR, tRsq));
    double njetScaleFactor = hNJetScaleFactors->GetBinContent(hNJetScaleFactors->FindFixBin(tNJets));

    double nisrCorr=1;
    if (physProc==2) {
      if (NISRJets==0) nisrCorr=1;
      else if (NISRJets==1) nisrCorr=0.920;
      else if (NISRJets==2) nisrCorr=0.821;
      else if (NISRJets==3) nisrCorr=0.715;
      else if (NISRJets==4) nisrCorr=0.662;
      else if (NISRJets==5) nisrCorr=0.561;
      else nisrCorr=0.511;
    }

    if (scaleFactor==0 || scaleFactor>2) {
      cout << scaleFactor << endl;
      scaleFactor=1;
    }
    if (njetScaleFactor==0 || njetScaleFactor>2) {
      cout << njetScaleFactor << endl;
      njetScaleFactor=1;
    }

    mcScaleFactor=scaleFactor*njetScaleFactor*nisrCorr;

    outputTree->Fill();
  }

  outputFile->Write();
  outputFile->Close();

}
