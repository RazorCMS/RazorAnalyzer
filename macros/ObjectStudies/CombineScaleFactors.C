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
#include <TStyle.h>                
#include <TGraphAsymmErrors.h>                

#endif





//=== MAIN MACRO ================================================================================================= 

void CombineScaleFactors() {


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  TFile *recoSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiencySF_muEleTracking_2016_average.root","READ");
  TH2F *histRecoMuonSF = (TH2F*)recoSFFile->Get("muon");
  TH2F *histRecoElectronSF = (TH2F*)recoSFFile->Get("h2_scaleFactorsEGamma");
 
  TFile *tightMuonSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016_26p4_Golden.root","READ");
  TH2F *histTightMuonSF = (TH2F*)tightMuonSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");

  TFile *vetoMuonSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorReco_2016_26p4_Golden.root","READ");
  TH2F *histVetoMuonSF = (TH2F*)vetoMuonSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");

  TFile *tightElectronSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016_26p4_Golden.root","READ");
  TH2F *histTightElectronSF = (TH2F*)tightElectronSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");

  TFile *vetoElectronSFFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorReco_2016_26p4_Golden.root","READ");
  TH2F *histVetoElectronSF = (TH2F*)vetoElectronSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");

  //tight electrons 
  const int NPtBinsTightElectron = 5;
  const int NEtaBinsTightElectron = 15;
  double ptBinsTightElectron[NPtBinsTightElectron+1] = {20, 30, 40, 50, 70, 7000};
  double etaBinsTightElectron[NEtaBinsTightElectron+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.444, 1.566, 1.63, 1.8, 2.0, 2.2, 2.3, 2.4, 2.5};

  //veto electrons
  const int NPtBinsVetoElectron = 7;
  const int NEtaBinsVetoElectron = 15;
  double ptBinsVetoElectron[NPtBinsVetoElectron+1] = {10, 15, 20, 30, 40, 50, 70, 7000};
  double etaBinsVetoElectron[NEtaBinsVetoElectron+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.444, 1.566, 1.63, 1.8, 2.0, 2.2, 2.3, 2.4, 2.5};


  //tight muons
  const int NPtBinsTightMuon = 5;
  const int NEtaBinsTightMuon = 9;
  double ptBinsTightMuon[NPtBinsTightMuon+1] = {20, 30, 40, 50, 70, 7000};
  double etaBinsTightMuon[NEtaBinsTightMuon+1] = {0,0.5,0.6,1.0,1.1,1.5,1.6,2.0,2.1,2.4};
 
  //veto muons
  const int NPtBinsVetoMuon = 8;
  const int NEtaBinsVetoMuon = 9;
  double ptBinsVetoMuon[NPtBinsVetoMuon+1] = {10, 15, 20, 30, 40, 50, 70, 100, 7000};
  double etaBinsVetoMuon[NEtaBinsVetoMuon+1] = {0,0.5,0.6,1.0,1.1,1.5,1.6,2.0,2.1,2.4};



  TH2F *outputVetoElectronSF = new TH2F( "ScaleFactor_VetoElectronSelectionEffDenominatorGen", ";p_{T} [GeV/c] ; #eta; Scale Factor", NPtBinsVetoElectron, ptBinsVetoElectron, NEtaBinsVetoElectron, etaBinsVetoElectron);
  for (int i=1; i< NPtBinsVetoElectron+1; i++) {
    for (int j=1; j< NEtaBinsVetoElectron+1; j++) {
      double tmpPt = outputVetoElectronSF->GetXaxis()->GetBinCenter(i);
      double tmpEta = outputVetoElectronSF->GetYaxis()->GetBinCenter(j);

      double recoEffSF = 0.5 * histRecoElectronSF->GetBinContent( histRecoElectronSF->GetXaxis()->FindFixBin(tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30))
	+ 0.5 * histRecoElectronSF->GetBinContent( histRecoElectronSF->GetXaxis()->FindFixBin(-1*tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30));
      double recoEffSFErr = histRecoElectronSF->GetBinError( histRecoElectronSF->GetXaxis()->FindFixBin(tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30));

      double sf = histVetoElectronSF->GetBinContent( histVetoElectronSF->GetXaxis()->FindFixBin(tmpPt),  histVetoElectronSF->GetYaxis()->FindFixBin(tmpEta)) * recoEffSF;
      double sfErr = sf* sqrt ( pow( histVetoElectronSF->GetBinError( histVetoElectronSF->GetXaxis()->FindFixBin(tmpPt),  histVetoElectronSF->GetYaxis()->FindFixBin(tmpEta)) / histVetoElectronSF->GetBinContent( histVetoElectronSF->GetXaxis()->FindFixBin(tmpPt),  histVetoElectronSF->GetYaxis()->FindFixBin(tmpEta)) , 2) + pow(recoEffSFErr / recoEffSF,2) );
      outputVetoElectronSF->SetBinContent( i,j, sf);
      outputVetoElectronSF->SetBinError( i,j, sfErr);
      cout << i << " " << j << " : " << tmpPt << " " << tmpEta << " : " 
	   << histVetoElectronSF->GetBinContent( histVetoElectronSF->GetXaxis()->FindFixBin(tmpPt),  histVetoElectronSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- "
	   << histVetoElectronSF->GetBinError( histVetoElectronSF->GetXaxis()->FindFixBin(tmpPt),  histVetoElectronSF->GetYaxis()->FindFixBin(tmpEta)) << " * "
	   << recoEffSF << " +/- " << recoEffSFErr
	   << " = "
	   << sf << " +/- " << sfErr << "\n";
    }
  }

  TH2F *outputTightElectronSF = new TH2F( "ScaleFactor_TightElectronSelectionEffDenominatorGen", ";p_{T} [GeV/c] ; #eta; Scale Factor", NPtBinsTightElectron, ptBinsTightElectron, NEtaBinsTightElectron, etaBinsTightElectron);
  for (int i=1; i< NPtBinsTightElectron+1; i++) {
    for (int j=1; j< NEtaBinsTightElectron+1; j++) {
      double tmpPt = outputTightElectronSF->GetXaxis()->GetBinCenter(i);
      double tmpEta = outputTightElectronSF->GetYaxis()->GetBinCenter(j);

      double recoEffSF = 0.5 * histRecoElectronSF->GetBinContent( histRecoElectronSF->GetXaxis()->FindFixBin(tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30))
	+ 0.5 * histRecoElectronSF->GetBinContent( histRecoElectronSF->GetXaxis()->FindFixBin(-1*tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30));
      double recoEffSFErr = histRecoElectronSF->GetBinError( histRecoElectronSF->GetXaxis()->FindFixBin(tmpEta), histRecoElectronSF->GetYaxis()->FindFixBin(30));

      double sf = histTightElectronSF->GetBinContent( histTightElectronSF->GetXaxis()->FindFixBin(tmpPt),  histTightElectronSF->GetYaxis()->FindFixBin(tmpEta)) * recoEffSF;
      double sfErr = sf* sqrt ( pow( histTightElectronSF->GetBinError( histTightElectronSF->GetXaxis()->FindFixBin(tmpPt),  histTightElectronSF->GetYaxis()->FindFixBin(tmpEta)) / histTightElectronSF->GetBinContent( histTightElectronSF->GetXaxis()->FindFixBin(tmpPt),  histTightElectronSF->GetYaxis()->FindFixBin(tmpEta)) , 2) + pow(recoEffSFErr / recoEffSF,2) );
      outputTightElectronSF->SetBinContent( i,j, sf);
      outputTightElectronSF->SetBinError( i,j, sfErr);
      cout << i << " " << j << " : " << tmpPt << " " << tmpEta << " : " 
	   << histTightElectronSF->GetBinContent( histTightElectronSF->GetXaxis()->FindFixBin(tmpPt),  histTightElectronSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- "
	   << histTightElectronSF->GetBinError( histTightElectronSF->GetXaxis()->FindFixBin(tmpPt),  histTightElectronSF->GetYaxis()->FindFixBin(tmpEta)) << " * "
	   << recoEffSF << " +/- " << recoEffSFErr
	   << " = "
	   << sf << " +/- " << sfErr << "\n";
    }
  }

  TH2F *outputTightMuonSF = new TH2F( "ScaleFactor_TightMuonSelectionEffDenominatorGen", ";p_{T} [GeV/c] ; #eta; Scale Factor", NPtBinsTightMuon, ptBinsTightMuon, NEtaBinsTightMuon, etaBinsTightMuon);
  for (int i=1; i< NPtBinsTightMuon+1; i++) {
    for (int j=1; j< NEtaBinsTightMuon+1; j++) {
      double tmpPt = outputTightMuonSF->GetXaxis()->GetBinCenter(i);
      double tmpEta = outputTightMuonSF->GetYaxis()->GetBinCenter(j);

      double recoEffSF = 0.5 * histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(30), histRecoMuonSF->GetYaxis()->FindFixBin(tmpEta))
	+ 0.5 * histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(30), histRecoMuonSF->GetYaxis()->FindFixBin(-1*tmpEta));
      
      double sf = histTightMuonSF->GetBinContent( histTightMuonSF->GetXaxis()->FindFixBin(tmpPt),  histTightMuonSF->GetYaxis()->FindFixBin(tmpEta)) * recoEffSF;
      double sfErr = sf* sqrt ( pow( histTightMuonSF->GetBinError( histTightMuonSF->GetXaxis()->FindFixBin(tmpPt),  histTightMuonSF->GetYaxis()->FindFixBin(tmpEta)) / histTightMuonSF->GetBinContent( histTightMuonSF->GetXaxis()->FindFixBin(tmpPt),  histTightMuonSF->GetYaxis()->FindFixBin(tmpEta)) , 2) + pow(0.003,2) );
      outputTightMuonSF->SetBinContent( i,j, sf);
      outputTightMuonSF->SetBinError( i,j, sfErr);
      cout << i << " " << j << " : " << tmpPt << " " << tmpEta << " : " 
	   << histTightMuonSF->GetBinContent( histTightMuonSF->GetXaxis()->FindFixBin(tmpPt),  histTightMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- "
	   << histTightMuonSF->GetBinError( histTightMuonSF->GetXaxis()->FindFixBin(tmpPt),  histTightMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " * "
	   << histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(tmpPt), histRecoMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- 0.003 "
	   << " = "
	   << sf << " +/- " << sfErr << "\n";
    }
  }

  TH2F *outputVetoMuonSF = new TH2F( "ScaleFactor_VetoMuonSelectionEffDenominatorGen", ";p_{T} [GeV/c] ; #eta; Scale Factor", NPtBinsVetoMuon, ptBinsVetoMuon, NEtaBinsVetoMuon, etaBinsVetoMuon);
  for (int i=1; i< NPtBinsVetoMuon+1; i++) {
    for (int j=1; j< NEtaBinsVetoMuon+1; j++) {
      double tmpPt = outputVetoMuonSF->GetXaxis()->GetBinCenter(i);
      double tmpEta = outputVetoMuonSF->GetYaxis()->GetBinCenter(j);

      double recoEffSF = 0.5 * histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(30), histRecoMuonSF->GetYaxis()->FindFixBin(tmpEta))
	+ 0.5 * histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(30), histRecoMuonSF->GetYaxis()->FindFixBin(-1*tmpEta));
      
      double sf = histVetoMuonSF->GetBinContent( histVetoMuonSF->GetXaxis()->FindFixBin(tmpPt),  histVetoMuonSF->GetYaxis()->FindFixBin(tmpEta)) * recoEffSF;
      double sfErr = sf* sqrt ( pow( histVetoMuonSF->GetBinError( histVetoMuonSF->GetXaxis()->FindFixBin(tmpPt),  histVetoMuonSF->GetYaxis()->FindFixBin(tmpEta)) / histVetoMuonSF->GetBinContent( histVetoMuonSF->GetXaxis()->FindFixBin(tmpPt),  histVetoMuonSF->GetYaxis()->FindFixBin(tmpEta)) , 2) + pow(0.003,2) );
      outputVetoMuonSF->SetBinContent( i,j, sf);
      outputVetoMuonSF->SetBinError( i,j, sfErr);
      cout << i << " " << j << " : " << tmpPt << " " << tmpEta << " : " 
	   << histVetoMuonSF->GetBinContent( histVetoMuonSF->GetXaxis()->FindFixBin(tmpPt),  histVetoMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- "
	   << histVetoMuonSF->GetBinError( histVetoMuonSF->GetXaxis()->FindFixBin(tmpPt),  histVetoMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " * "
	   << histRecoMuonSF->GetBinContent( histRecoMuonSF->GetXaxis()->FindFixBin(tmpPt), histRecoMuonSF->GetYaxis()->FindFixBin(tmpEta)) << " +/- 0.003 "
	   << " = "
	   << sf << " +/- " << sfErr << "\n";
    }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("efficiency_results_TightElectronSelectionEffDenominatorGen_2016_26p4_Golden.root", "RECREATE");
  file->cd();
  file->WriteTObject(outputTightElectronSF, "ScaleFactor_TightElectronSelectionEffDenominatorGen", "WriteDelete");
  file->Close();
  delete file;       

  file = TFile::Open("efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_26p4_Golden.root", "RECREATE");
  file->cd();
  file->WriteTObject(outputVetoElectronSF, "ScaleFactor_VetoElectronSelectionEffDenominatorGen", "WriteDelete");
  file->Close();
  delete file;       

  file = TFile::Open("efficiency_results_TightMuonSelectionEffDenominatorGen_2016_26p4_Golden.root", "RECREATE");
  file->cd();
  file->WriteTObject(outputTightMuonSF, "ScaleFactor_TightMuonSelectionEffDenominatorGen", "WriteDelete");
  file->Close();
  delete file; 

  file = TFile::Open("efficiency_results_VetoMuonSelectionEffDenominatorGen_2016_26p4_Golden.root", "RECREATE");
  file->cd();
  file->WriteTObject(outputVetoMuonSF, "ScaleFactor_VetoMuonSelectionEffDenominatorGen", "WriteDelete");
  file->Close();
  delete file;       


}



