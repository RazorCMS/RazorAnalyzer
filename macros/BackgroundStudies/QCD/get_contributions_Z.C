#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include <CalStyleRemix.hh>

void get_contributions_Z() {

  //--------------------------------------------------------------
  //
  // setup
  //
  //--------------------------------------------------------------

  TCanvas *c = MakeCanvas("c","c",800,600);

  TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/DoubleMuon_Run2015D_Golden.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/DYJetsToLL_amcatnlo_2137pb_weighted.root","read");
  TFile *fVV = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/VV.root","read");
  TFile *fTT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_leptonic.root","read");

  TString cut_str="(box==1)*(zMass>60 && zMass<120)*passedDileptonTrigger*weight*puWeight*(NJets40>1)*(MR>300)";

  //TString cut_str_mr1="*(MR>500 && MR<600)";
  //TString cut_str_mr2="*(MR>600 && MR<700)";
  //TString cut_str_mr3="*(MR>700 && MR<800)";
  //TString cut_str_mr4="*(MR>800 && MR<900)";
  //TString cut_str_mr5="*(MR>900 && MR<1000)";

  TString cut_str_mr1="*(Rsq > 0 && Rsq < 0.05)";
  TString cut_str_mr2="*(Rsq > 0.05 && Rsq < 0.1)";
  TString cut_str_mr3="*(Rsq > 0.1 && Rsq < 0.15)";
  TString cut_str_mr4="*(Rsq > 0.15 && Rsq < 0.2)";
  TString cut_str_mr5="*(Rsq > 0.2 && Rsq < 1.0)";

  // draw Z mass
  Float_t nbin=30, xmin=60, xmax=120;
  TH1F *hMZ_D = new TH1F("hMZ_D", "hMZ_D", nbin, xmin, xmax); hMZ_D->Sumw2();
  TH1F *hMZ_Z = new TH1F("hMZ_Z", "hMZ_Z", nbin, xmin, xmax); hMZ_Z->Sumw2();
  TH1F *hMZ_VV = new TH1F("hMZ_VV", "hMZ_VV", nbin, xmin, xmax); hMZ_VV->Sumw2();
  TH1F *hMZ_TT = new TH1F("hMZ_TT", "hMZ_TT", nbin, xmin, xmax); hMZ_TT->Sumw2();

  // draw MR
  nbin=15; xmin=500; xmax=2000;

  TH1F *hMR_D = new TH1F("hMR_D", "hMR_D", nbin, xmin, xmax); hMR_D->Sumw2();
  TH1F *hMR_Z = new TH1F("hMR_Z", "hMR_Z", nbin, xmin, xmax); hMR_Z->Sumw2();
  TH1F *hMR_VV = new TH1F("hMR_VV", "hMR_VV", nbin, xmin, xmax); hMR_VV->Sumw2();
  TH1F *hMR_TT = new TH1F("hMR_TT", "hMR_TT", nbin, xmin, xmax); hMR_TT->Sumw2();

  TH1F *hMR_pass_D = new TH1F("hMR_pass_D", "hMR_pass_D", nbin, xmin, xmax); hMR_pass_D->Sumw2();
  TH1F *hMR_pass_Z = new TH1F("hMR_pass_Z", "hMR_pass_Z", nbin, xmin, xmax); hMR_pass_Z->Sumw2();
  TH1F *hMR_pass_VV = new TH1F("hMR_pass_VV", "hMR_pass_VV", nbin, xmin, xmax); hMR_pass_VV->Sumw2();
  TH1F *hMR_pass_TT = new TH1F("hMR_pass_TT", "hMR_pass_TT", nbin, xmin, xmax); hMR_pass_TT->Sumw2();

  TH1F *hMR_fail_D = new TH1F("hMR_fail_D", "hMR_fail_D", nbin, xmin, xmax); hMR_fail_D->Sumw2();
  TH1F *hMR_fail_Z = new TH1F("hMR_fail_Z", "hMR_fail_Z", nbin, xmin, xmax); hMR_fail_Z->Sumw2();
  TH1F *hMR_fail_VV = new TH1F("hMR_fail_VV", "hMR_fail_VV", nbin, xmin, xmax); hMR_fail_VV->Sumw2();
  TH1F *hMR_fail_TT = new TH1F("hMR_fail_TT", "hMR_fail_TT", nbin, xmin, xmax); hMR_fail_TT->Sumw2();

  // draw Rsq
  nbin=10; xmin=0.0; xmax=1.0;

  TH1F *hRsq_D = new TH1F("hRsq_D", "hRsq_D", nbin, xmin, xmax); hRsq_D->Sumw2();
  TH1F *hRsq_Z = new TH1F("hRsq_Z", "hRsq_Z", nbin, xmin, xmax); hRsq_Z->Sumw2();
  TH1F *hRsq_VV = new TH1F("hRsq_VV", "hRsq_VV", nbin, xmin, xmax); hRsq_VV->Sumw2();
  TH1F *hRsq_TT = new TH1F("hRsq_TT", "hRsq_TT", nbin, xmin, xmax); hRsq_TT->Sumw2();

  TH1F *hRsq_pass_D = new TH1F("hRsq_pass_D", "hRsq_pass_D", nbin, xmin, xmax); hRsq_pass_D->Sumw2();
  TH1F *hRsq_pass_Z = new TH1F("hRsq_pass_Z", "hRsq_pass_Z", nbin, xmin, xmax); hRsq_pass_Z->Sumw2();
  TH1F *hRsq_pass_VV = new TH1F("hRsq_pass_VV", "hRsq_pass_VV", nbin, xmin, xmax); hRsq_pass_VV->Sumw2();
  TH1F *hRsq_pass_TT = new TH1F("hRsq_pass_TT", "hRsq_pass_TT", nbin, xmin, xmax); hRsq_pass_TT->Sumw2();

  TH1F *hRsq_fail_D = new TH1F("hRsq_fail_D", "hRsq_fail_D", nbin, xmin, xmax); hRsq_fail_D->Sumw2();
  TH1F *hRsq_fail_Z = new TH1F("hRsq_fail_Z", "hRsq_fail_Z", nbin, xmin, xmax); hRsq_fail_Z->Sumw2();
  TH1F *hRsq_fail_VV = new TH1F("hRsq_fail_VV", "hRsq_fail_VV", nbin, xmin, xmax); hRsq_fail_VV->Sumw2();
  TH1F *hRsq_fail_TT = new TH1F("hRsq_fail_TT", "hRsq_fail_TT", nbin, xmin, xmax); hRsq_fail_TT->Sumw2();

  // draw dPhiRazor
  nbin=8; 
  Float_t xbins[9] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.0, 3.2};

  TH1F *hDPhiR_D = new TH1F("hDPhiR_D", "hDPhiR_D", nbin, &xbins[0]); hDPhiR_D->Sumw2();
  TH1F *hDPhiR_Z = new TH1F("hDPhiR_Z", "hDPhiR_Z", nbin, &xbins[0]); hDPhiR_Z->Sumw2();
  TH1F *hDPhiR_VV = new TH1F("hDPhiR_VV", "hDPhiR_VV", nbin, &xbins[0]); hDPhiR_VV->Sumw2();
  TH1F *hDPhiR_TT = new TH1F("hDPhiR_TT", "hDPhiR_TT", nbin, &xbins[0]); hDPhiR_TT->Sumw2();

  TH1F *hDPhiR_mr1_D = new TH1F("hDPhiR_mr1_D", "hDPhiR_mr1_D", nbin, &xbins[0]); hDPhiR_mr1_D->Sumw2();
  TH1F *hDPhiR_mr1_Z = new TH1F("hDPhiR_mr1_Z", "hDPhiR_mr1_Z", nbin, &xbins[0]); hDPhiR_mr1_Z->Sumw2();
  TH1F *hDPhiR_mr1_TT = new TH1F("hDPhiR_mr1_TT", "hDPhiR_mr1_TT", nbin, &xbins[0]); hDPhiR_mr1_TT->Sumw2();
  TH1F *hDPhiR_mr1_VV = new TH1F("hDPhiR_mr1_VV", "hDPhiR_mr1_VV", nbin, &xbins[0]); hDPhiR_mr1_VV->Sumw2();

  TH1F *hDPhiR_mr2_D = new TH1F("hDPhiR_mr2_D", "hDPhiR_mr2_D", nbin, &xbins[0]); hDPhiR_mr2_D->Sumw2();
  TH1F *hDPhiR_mr2_Z = new TH1F("hDPhiR_mr2_Z", "hDPhiR_mr2_Z", nbin, &xbins[0]); hDPhiR_mr2_Z->Sumw2();
  TH1F *hDPhiR_mr2_TT = new TH1F("hDPhiR_mr2_TT", "hDPhiR_mr2_TT", nbin, &xbins[0]); hDPhiR_mr2_TT->Sumw2();
  TH1F *hDPhiR_mr2_VV = new TH1F("hDPhiR_mr2_VV", "hDPhiR_mr2_VV", nbin, &xbins[0]); hDPhiR_mr2_VV->Sumw2();

  TH1F *hDPhiR_mr3_D = new TH1F("hDPhiR_mr3_D", "hDPhiR_mr3_D", nbin, &xbins[0]); hDPhiR_mr3_D->Sumw2();
  TH1F *hDPhiR_mr3_Z = new TH1F("hDPhiR_mr3_Z", "hDPhiR_mr3_Z", nbin, &xbins[0]); hDPhiR_mr3_Z->Sumw2();
  TH1F *hDPhiR_mr3_TT = new TH1F("hDPhiR_mr3_TT", "hDPhiR_mr3_TT", nbin, &xbins[0]); hDPhiR_mr3_TT->Sumw2();
  TH1F *hDPhiR_mr3_VV = new TH1F("hDPhiR_mr3_VV", "hDPhiR_mr3_VV", nbin, &xbins[0]); hDPhiR_mr3_VV->Sumw2();

  TH1F *hDPhiR_mr4_D = new TH1F("hDPhiR_mr4_D", "hDPhiR_mr4_D", nbin, &xbins[0]); hDPhiR_mr4_D->Sumw2();
  TH1F *hDPhiR_mr4_Z = new TH1F("hDPhiR_mr4_Z", "hDPhiR_mr4_Z", nbin, &xbins[0]); hDPhiR_mr4_Z->Sumw2();
  TH1F *hDPhiR_mr4_TT = new TH1F("hDPhiR_mr4_TT", "hDPhiR_mr4_TT", nbin, &xbins[0]); hDPhiR_mr4_TT->Sumw2();
  TH1F *hDPhiR_mr4_VV = new TH1F("hDPhiR_mr4_VV", "hDPhiR_mr4_VV", nbin, &xbins[0]); hDPhiR_mr4_VV->Sumw2();

  TH1F *hDPhiR_mr5_D = new TH1F("hDPhiR_mr5_D", "hDPhiR_mr5_D", nbin, &xbins[0]); hDPhiR_mr5_D->Sumw2();
  TH1F *hDPhiR_mr5_Z = new TH1F("hDPhiR_mr5_Z", "hDPhiR_mr5_Z", nbin, &xbins[0]); hDPhiR_mr5_Z->Sumw2();
  TH1F *hDPhiR_mr5_TT = new TH1F("hDPhiR_mr5_TT", "hDPhiR_mr5_TT", nbin, &xbins[0]); hDPhiR_mr5_TT->Sumw2();
  TH1F *hDPhiR_mr5_VV = new TH1F("hDPhiR_mr5_VV", "hDPhiR_mr5_VV", nbin, &xbins[0]); hDPhiR_mr5_VV->Sumw2();

  // open trees
  TTree *tD = (TTree*) fD->Get("QCDTree");
  TTree *tZ = (TTree*) fZ->Get("QCDTree");
  TTree *tVV = (TTree*) fVV->Get("QCDTree");
  TTree *tTT = (TTree*) fTT->Get("QCDTree");

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  // draw MZ 
  tD->Draw("zMass>>hMZ_D", cut_str);
  tZ->Draw("zMass>>hMZ_Z", cut_str);
  tVV->Draw("zMass>>hMZ_VV", cut_str);
  tTT->Draw("zMass>>hMZ_TT", cut_str);

  // draw MR (again)
  tD->Draw("MR>>hMR_D", cut_str);
  tZ->Draw("MR>>hMR_Z", cut_str);
  tVV->Draw("MR>>hMR_VV", cut_str);
  tTT->Draw("MR>>hMR_TT", cut_str);
  
  tD->Draw("MR>>hMR_pass_D", cut_str+"*(abs(dPhiRazor)>2.8)");
  tZ->Draw("MR>>hMR_pass_Z", cut_str+"*(abs(dPhiRazor)>2.8)");
  tVV->Draw("MR>>hMR_pass_VV", cut_str+"*(abs(dPhiRazor)>2.8)");
  tTT->Draw("MR>>hMR_pass_TT", cut_str+"*(abs(dPhiRazor)>2.8)");
  
  tD->Draw("MR>>hMR_fail_D", cut_str+"*(abs(dPhiRazor)<2.8)");
  tZ->Draw("MR>>hMR_fail_Z", cut_str+"*(abs(dPhiRazor)<2.8)");
  tVV->Draw("MR>>hMR_fail_VV", cut_str+"*(abs(dPhiRazor)<2.8)");
  tTT->Draw("MR>>hMR_fail_TT", cut_str+"*(abs(dPhiRazor)<2.8)");

  // draw Rsq (again)
  tD->Draw("Rsq>>hRsq_D", cut_str);
  tZ->Draw("Rsq>>hRsq_Z", cut_str);
  tVV->Draw("Rsq>>hRsq_VV", cut_str);
  tTT->Draw("Rsq>>hRsq_TT", cut_str);

  tD->Draw("Rsq>>hRsq_pass_D", cut_str+"*(abs(dPhiRazor)>2.8)");
  tZ->Draw("Rsq>>hRsq_pass_Z", cut_str+"*(abs(dPhiRazor)>2.8)");
  tVV->Draw("Rsq>>hRsq_pass_VV", cut_str+"*(abs(dPhiRazor)>2.8)");
  tTT->Draw("Rsq>>hRsq_pass_TT", cut_str+"*(abs(dPhiRazor)>2.8)");
  
  tD->Draw("Rsq>>hRsq_fail_D", cut_str+"*(abs(dPhiRazor)<2.8)");
  tZ->Draw("Rsq>>hRsq_fail_Z", cut_str+"*(abs(dPhiRazor)<2.8)");
  tVV->Draw("Rsq>>hRsq_fail_VV", cut_str+"*(abs(dPhiRazor)<2.8)");
  tTT->Draw("Rsq>>hRsq_fail_TT", cut_str+"*(abs(dPhiRazor)<2.8)");

  // draw dPhiR (again)
  tD->Draw("dPhiRazor>>hDPhiR_D", cut_str);
  tZ->Draw("dPhiRazor>>hDPhiR_Z", cut_str);
  tVV->Draw("dPhiRazor>>hDPhiR_VV", cut_str);
  tTT->Draw("dPhiRazor>>hDPhiR_TT", cut_str);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr1_D", cut_str+cut_str_mr1);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr1_Z", cut_str+cut_str_mr1);
  tVV->Draw("abs(dPhiRazor)>>hDPhiR_mr1_VV", cut_str+cut_str_mr1);
  tTT->Draw("abs(dPhiRazor)>>hDPhiR_mr1_TT", cut_str+cut_str_mr1);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr2_D", cut_str+cut_str_mr2);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr2_Z", cut_str+cut_str_mr2);
  tVV->Draw("abs(dPhiRazor)>>hDPhiR_mr2_VV", cut_str+cut_str_mr2);
  tTT->Draw("abs(dPhiRazor)>>hDPhiR_mr2_TT", cut_str+cut_str_mr2);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr3_D", cut_str+cut_str_mr3);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr3_Z", cut_str+cut_str_mr3);
  tVV->Draw("abs(dPhiRazor)>>hDPhiR_mr3_VV", cut_str+cut_str_mr3);
  tTT->Draw("abs(dPhiRazor)>>hDPhiR_mr3_TT", cut_str+cut_str_mr3);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr4_D", cut_str+cut_str_mr4);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr4_Z", cut_str+cut_str_mr4);
  tVV->Draw("abs(dPhiRazor)>>hDPhiR_mr4_VV", cut_str+cut_str_mr4);
  tTT->Draw("abs(dPhiRazor)>>hDPhiR_mr4_TT", cut_str+cut_str_mr4);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr5_D", cut_str+cut_str_mr5);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr5_Z", cut_str+cut_str_mr5);
  tVV->Draw("abs(dPhiRazor)>>hDPhiR_mr5_VV", cut_str+cut_str_mr5);
  tTT->Draw("abs(dPhiRazor)>>hDPhiR_mr5_TT", cut_str+cut_str_mr5);

  //--------------------------------------------------------------
  //
  // configure draw options
  //
  //--------------------------------------------------------------

  //data
  InitData(hMZ_D,"","",kBlack);
  InitData(hMR_D,"","",kBlack);
  InitData(hMR_pass_D,"","",kBlack);
  InitData(hMR_fail_D,"","",kBlack);
  InitData(hRsq_D,"","",kBlack);
  InitData(hRsq_pass_D,"","",kBlack);
  InitData(hRsq_fail_D,"","",kBlack);
  InitData(hDPhiR_D,"","",kBlack);

  InitData(hDPhiR_mr1_D,"","",kBlack);
  InitData(hDPhiR_mr2_D,"","",kBlack);
  InitData(hDPhiR_mr3_D,"","",kBlack);
  InitData(hDPhiR_mr4_D,"","",kBlack);
  InitData(hDPhiR_mr5_D,"","",kBlack);

  // z
  InitHist(hMZ_Z,"","",kBlue+1);
  InitHist(hMR_Z,"","",kBlue+1);
  InitHist(hMR_pass_Z,"","",kBlue+1);
  InitHist(hMR_fail_Z,"","",kBlue+1);
  InitHist(hRsq_Z,"","",kBlue+1);
  InitHist(hRsq_pass_Z,"","",kBlue+1);
  InitHist(hRsq_fail_Z,"","",kBlue+1);
  InitHist(hDPhiR_Z,"","",kBlue+1);

  InitHist(hDPhiR_mr1_Z,"","",kBlue+1);
  InitHist(hDPhiR_mr2_Z,"","",kBlue+1);
  InitHist(hDPhiR_mr3_Z,"","",kBlue+1);
  InitHist(hDPhiR_mr4_Z,"","",kBlue+1);
  InitHist(hDPhiR_mr5_Z,"","",kBlue+1);

  // vv
  InitHist(hMZ_VV,"","",kViolet+3);
  InitHist(hMR_VV,"","",kViolet+3);
  InitHist(hMR_pass_VV,"","",kViolet+3);
  InitHist(hMR_fail_VV,"","",kViolet+3);
  InitHist(hRsq_VV,"","",kViolet+3);
  InitHist(hRsq_pass_VV,"","",kViolet+3);
  InitHist(hRsq_fail_VV,"","",kViolet+3);
  InitHist(hDPhiR_VV,"","",kViolet+3);

  InitHist(hDPhiR_mr1_VV,"","",kViolet+3);
  InitHist(hDPhiR_mr2_VV,"","",kViolet+3);
  InitHist(hDPhiR_mr3_VV,"","",kViolet+3);
  InitHist(hDPhiR_mr4_VV,"","",kViolet+3);
  InitHist(hDPhiR_mr5_VV,"","",kViolet+3);

  // tt
  InitHist(hMZ_TT,"","",kGreen+2);
  InitHist(hMR_TT,"","",kGreen+2);
  InitHist(hMR_pass_TT,"","",kGreen+2);
  InitHist(hMR_fail_TT,"","",kGreen+2);
  InitHist(hRsq_TT,"","",kGreen+2);
  InitHist(hRsq_pass_TT,"","",kGreen+2);
  InitHist(hRsq_fail_TT,"","",kGreen+2);
  InitHist(hDPhiR_TT,"","",kGreen+2);

  InitHist(hDPhiR_mr1_TT,"","",kGreen+2);
  InitHist(hDPhiR_mr2_TT,"","",kGreen+2);
  InitHist(hDPhiR_mr3_TT,"","",kGreen+2);
  InitHist(hDPhiR_mr4_TT,"","",kGreen+2);
  InitHist(hDPhiR_mr5_TT,"","",kGreen+2);
  
  TLegend *leg = new TLegend(0.67,0.57,0.89,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  leg->AddEntry(hMR_D, "Data", "pel");
  leg->AddEntry(hMR_Z, "Zll+jets", "f");  
  leg->AddEntry(hMZ_VV, "VV", "f");  
  leg->AddEntry(hMZ_TT, "TT", "f");  

  hMZ_TT->Add(hMZ_VV); 
  hMR_TT->Add(hMR_VV); 
  hMR_pass_TT->Add(hMR_pass_VV); 
  hMR_fail_TT->Add(hMR_fail_VV); 
  hRsq_TT->Add(hRsq_VV); 
  hRsq_pass_TT->Add(hRsq_pass_VV); 
  hRsq_fail_TT->Add(hRsq_fail_VV); 
  hDPhiR_TT->Add(hDPhiR_VV); 
  
  hDPhiR_mr1_TT->Add(hDPhiR_mr1_VV); 
  hDPhiR_mr2_TT->Add(hDPhiR_mr2_VV); 
  hDPhiR_mr3_TT->Add(hDPhiR_mr3_VV); 
  hDPhiR_mr4_TT->Add(hDPhiR_mr4_VV); 
  hDPhiR_mr5_TT->Add(hDPhiR_mr5_VV); 
  
  hMZ_Z->Add(hMZ_TT); 
  hMR_Z->Add(hMR_TT); 
  hMR_pass_Z->Add(hMR_pass_TT); 
  hMR_fail_Z->Add(hMR_fail_TT); 
  hRsq_Z->Add(hRsq_TT); 
  hRsq_pass_Z->Add(hRsq_pass_TT);
  hRsq_fail_Z->Add(hRsq_fail_TT);
  hDPhiR_Z->Add(hDPhiR_TT); 

  hDPhiR_mr1_Z->Add(hDPhiR_mr1_TT); 
  hDPhiR_mr2_Z->Add(hDPhiR_mr2_TT); 
  hDPhiR_mr3_Z->Add(hDPhiR_mr3_TT); 
  hDPhiR_mr4_Z->Add(hDPhiR_mr4_TT); 
  hDPhiR_mr5_Z->Add(hDPhiR_mr5_TT); 

  Float_t scale=hMZ_D->Integral()/hMZ_Z->Integral();
  hMZ_Z->Scale(scale);
  hMZ_TT->Scale(scale);
  hMZ_VV->Scale(scale);

  scale=hMR_D->Integral()/hMR_Z->Integral();
  hMR_Z->Scale(scale);
  hMR_TT->Scale(scale);
  hMR_VV->Scale(scale);

  scale=hMR_pass_D->Integral()/hMR_pass_Z->Integral();
  hMR_pass_Z->Scale(scale);
  hMR_pass_TT->Scale(scale);
  hMR_pass_VV->Scale(scale);

  scale=hMR_fail_D->Integral()/hMR_fail_Z->Integral();
  hMR_fail_Z->Scale(scale);
  hMR_fail_TT->Scale(scale);
  hMR_fail_VV->Scale(scale);

  scale=hRsq_D->Integral()/hRsq_Z->Integral();
  hRsq_Z->Scale(scale);
  hRsq_TT->Scale(scale);
  hRsq_VV->Scale(scale);

  scale=hRsq_pass_D->Integral()/hRsq_pass_Z->Integral();
  hRsq_pass_Z->Scale(scale);
  hRsq_pass_TT->Scale(scale);
  hRsq_pass_VV->Scale(scale);

  scale=hRsq_fail_D->Integral()/hRsq_fail_Z->Integral();
  hRsq_fail_Z->Scale(scale);
  hRsq_fail_TT->Scale(scale);
  hRsq_fail_VV->Scale(scale);

  scale=hDPhiR_D->Integral()/hDPhiR_Z->Integral();
  hDPhiR_Z->Scale(scale);
  hDPhiR_TT->Scale(scale);
  hDPhiR_VV->Scale(scale);

  scale=hDPhiR_mr1_D->Integral()/hDPhiR_mr1_Z->Integral();
  hDPhiR_mr1_Z->Scale(scale);
  hDPhiR_mr1_TT->Scale(scale);
  hDPhiR_mr1_VV->Scale(scale);

  scale=hDPhiR_mr2_D->Integral()/hDPhiR_mr2_Z->Integral();
  hDPhiR_mr2_Z->Scale(scale);
  hDPhiR_mr2_TT->Scale(scale);
  hDPhiR_mr2_VV->Scale(scale);

  scale=hDPhiR_mr3_D->Integral()/hDPhiR_mr3_Z->Integral();
  hDPhiR_mr3_Z->Scale(scale);
  hDPhiR_mr3_TT->Scale(scale);
  hDPhiR_mr3_VV->Scale(scale);

  scale=hDPhiR_mr4_D->Integral()/hDPhiR_mr4_Z->Integral();
  hDPhiR_mr4_Z->Scale(scale);
  hDPhiR_mr4_TT->Scale(scale);
  hDPhiR_mr4_VV->Scale(scale);

  scale=hDPhiR_mr5_D->Integral()/hDPhiR_mr5_Z->Integral();
  hDPhiR_mr5_Z->Scale(scale);
  hDPhiR_mr5_TT->Scale(scale);
  hDPhiR_mr5_VV->Scale(scale);


  //--------------------------------------------------------------
  //
  // draw
  //
  //--------------------------------------------------------------

  hMZ_Z->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMZ_Z->GetMaximum(), hMZ_D->GetMaximum()));
  hMZ_Z->GetXaxis()->SetTitle("M_{Z}");
  hMZ_Z->GetYaxis()->SetTitle("Events");
  hMZ_Z->SetTitle("");
  hMZ_Z->Draw("hist");
  hMZ_TT->Draw("hist same");
  hMZ_VV->Draw("hist same");
  hMZ_D->Draw("pe same");
  leg->Draw();
  c->SaveAs("MZ_zjets.png");

  hMR_Z->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_Z->GetMaximum(), hMR_D->GetMaximum()));
  hMR_Z->GetXaxis()->SetTitle("M_{R}");
  hMR_Z->GetYaxis()->SetTitle("Events");
  hMR_Z->SetTitle("");
  hMR_Z->Draw("hist");
  hMR_TT->Draw("hist same");
  hMR_VV->Draw("hist same");
  hMR_D->Draw("pe same");
  leg->Draw();
  c->SaveAs("MR_zjets.png");

  hMR_pass_Z->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_pass_Z->GetMaximum(), hMR_pass_D->GetMaximum()));
  hMR_pass_Z->GetXaxis()->SetTitle("M_{R}");
  hMR_pass_Z->GetYaxis()->SetTitle("Events");
  hMR_pass_Z->SetTitle("");
  hMR_pass_Z->Draw("hist");
  hMR_pass_TT->Draw("hist same");
  hMR_pass_VV->Draw("hist same");
  hMR_pass_D->Draw("pe same");
  leg->Draw();
  c->SaveAs("MR_pass_zjets.png");

  hMR_fail_Z->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_fail_Z->GetMaximum(), hMR_fail_D->GetMaximum()));
  hMR_fail_Z->GetXaxis()->SetTitle("M_{R}");
  hMR_fail_Z->GetYaxis()->SetTitle("Events");
  hMR_fail_Z->SetTitle("");
  hMR_fail_Z->Draw("hist");
  hMR_fail_TT->Draw("hist same");
  hMR_fail_VV->Draw("hist same");
  hMR_fail_D->Draw("pe same");
  leg->Draw();
  c->SaveAs("MR_fail_zjets.png");

  c->SetLogy(1);
  hRsq_Z->GetYaxis()->SetRangeUser(0.01,100*TMath::Max(hRsq_Z->GetMaximum(), hRsq_D->GetMaximum()));
  hRsq_Z->GetXaxis()->SetTitle("R^{2}");
  hRsq_Z->GetYaxis()->SetTitle("Events");
  hRsq_Z->SetTitle("");
  hRsq_Z->Draw("hist");
  hRsq_TT->Draw("hist same");
  hRsq_VV->Draw("hist same");
  hRsq_D->Draw("same ep");
  leg->Draw();
  c->SaveAs("Rsq_zjets.png");

  hRsq_pass_Z->GetYaxis()->SetRangeUser(0.01,100*TMath::Max(hRsq_pass_Z->GetMaximum(), hRsq_pass_D->GetMaximum()));
  hRsq_pass_Z->GetXaxis()->SetTitle("R^{2}");
  hRsq_pass_Z->GetYaxis()->SetTitle("Events");
  hRsq_pass_Z->SetTitle("");
  hRsq_pass_Z->Draw("hist");
  hRsq_pass_TT->Draw("hist same");
  hRsq_pass_VV->Draw("hist same");
  hRsq_pass_D->Draw("same ep");
  leg->Draw();
  c->SaveAs("Rsq_pass_zjets.png");

  hRsq_fail_Z->GetYaxis()->SetRangeUser(0.01,100*TMath::Max(hRsq_fail_Z->GetMaximum(), hRsq_fail_D->GetMaximum()));
  hRsq_fail_Z->GetXaxis()->SetTitle("R^{2}");
  hRsq_fail_Z->GetYaxis()->SetTitle("Events");
  hRsq_fail_Z->SetTitle("");
  hRsq_fail_Z->Draw("hist");
  hRsq_fail_TT->Draw("hist same");
  hRsq_fail_VV->Draw("hist same");
  hRsq_fail_D->Draw("same ep");
  leg->Draw();
  c->SaveAs("Rsq_fail_zjets.png");

  c->SetLogy(0);
  leg->SetX1NDC(0.39); leg->SetX2NDC(0.61);

  hDPhiR_Z->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_Z->GetMaximum(), hDPhiR_D->GetMaximum()));
  hDPhiR_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_Z->SetTitle("");
  hDPhiR_Z->Draw("hist");
  hDPhiR_TT->Draw("hist same");
  hDPhiR_VV->Draw("hist same");
  hDPhiR_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_zjets.png");

  hDPhiR_mr1_Z->GetYaxis()->SetRangeUser(0.0,2*TMath::Max(hDPhiR_mr1_Z->GetMaximum(), hDPhiR_mr1_D->GetMaximum()));
  hDPhiR_mr1_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_mr1_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_mr1_Z->SetTitle("0 < R^2 < 0.05");
  //hDPhiR_mr1_Z->SetTitle("500 < MR < 600");
  hDPhiR_mr1_Z->Draw("hist");
  hDPhiR_mr1_TT->Draw("hist same");
  hDPhiR_mr1_VV->Draw("hist same");
  hDPhiR_mr1_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_rsq1_zjets.png");
  //c->SaveAs("DPhiR_mr1_zjets.png");

  hDPhiR_mr2_Z->GetYaxis()->SetRangeUser(0.0,2*TMath::Max(hDPhiR_mr2_Z->GetMaximum(), hDPhiR_mr2_D->GetMaximum()));
  hDPhiR_mr2_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_mr2_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_mr2_Z->SetTitle("0.05 < R^2 < 0.1");
  //hDPhiR_mr2_Z->SetTitle("600 < MR < 700");
  hDPhiR_mr2_Z->Draw("hist");
  hDPhiR_mr2_TT->Draw("hist same");
  hDPhiR_mr2_VV->Draw("hist same");
  hDPhiR_mr2_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_rsq2_zjets.png");
  //c->SaveAs("DPhiR_mr2_zjets.png");

  hDPhiR_mr3_Z->GetYaxis()->SetRangeUser(0.0,2*TMath::Max(hDPhiR_mr3_Z->GetMaximum(), hDPhiR_mr3_D->GetMaximum()));
  hDPhiR_mr3_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_mr3_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_mr3_Z->SetTitle("0.1 < R^2 < 0.15");
  hDPhiR_mr3_Z->SetTitle("700 < MR < 800");
  hDPhiR_mr3_Z->Draw("hist");
  hDPhiR_mr3_TT->Draw("hist same");
  hDPhiR_mr3_VV->Draw("hist same");
  hDPhiR_mr3_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_rsq3_zjets.png");
  //c->SaveAs("DPhiR_mr3_zjets.png");

  hDPhiR_mr4_Z->GetYaxis()->SetRangeUser(0.0,2*TMath::Max(hDPhiR_mr4_Z->GetMaximum(), hDPhiR_mr4_D->GetMaximum()));
  hDPhiR_mr4_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_mr4_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_mr4_Z->SetTitle("0.15 < R^2 < 0.2");
  //hDPhiR_mr4_Z->SetTitle("800 < MR < 900");
  hDPhiR_mr4_Z->Draw("hist");
  hDPhiR_mr4_TT->Draw("hist same");
  hDPhiR_mr4_VV->Draw("hist same");
  hDPhiR_mr4_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_rsq4_zjets.png");
  //c->SaveAs("DPhiR_mr4_zjets.png");

  hDPhiR_mr5_Z->GetYaxis()->SetRangeUser(0.0,2*TMath::Max(hDPhiR_mr5_Z->GetMaximum(), hDPhiR_mr5_D->GetMaximum()));
  hDPhiR_mr5_Z->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_mr5_Z->GetYaxis()->SetTitle("Events");
  hDPhiR_mr5_Z->SetTitle("0.2 < R^2 < 1.0");
  //hDPhiR_mr5_Z->SetTitle("900 < MR < 1000");
  hDPhiR_mr5_Z->Draw("hist");
  hDPhiR_mr5_TT->Draw("hist same");
  hDPhiR_mr5_VV->Draw("hist same");
  hDPhiR_mr5_D->Draw("same pe");
  leg->Draw();
  c->SaveAs("DPhiR_rsq5_zjets.png");
  //c->SaveAs("DPhiR_mr5_zjets.png");

}
