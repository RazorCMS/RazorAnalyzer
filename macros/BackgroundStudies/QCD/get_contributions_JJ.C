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

void get_contributions_JJ() {

  //--------------------------------------------------------------
  //
  // setup
  //
  //--------------------------------------------------------------

  TCanvas *c = MakeCanvas("c","c",800,600);

  TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/JetHT_Run2015D_PRv4_Golden_skim.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/QCD_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_skim.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/ZJetsToNuNu_13TeV-madgraph_2137pb_skim.root","read");

  TString cut_str="weight*(box==100)*puWeight*(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)";
 
  TString cut_str_mr1="*(MR>500 && MR<600)";
  TString cut_str_mr2="*(MR>600 && MR<700)";
  TString cut_str_mr3="*(MR>700 && MR<800)";
  TString cut_str_mr4="*(MR>800 && MR<900)";
  TString cut_str_mr5="*(MR>900 && MR<1000)";

  TString cut_str_dat=cut_str+"*(HLTDecision[105]*HLTPrescale[105])";

  // draw MR
  Float_t nbin=5, xmin=500, xmax=1000;

  TH1F *hMR_D = new TH1F("hMR_D", "hMR_D", nbin, xmin, xmax); hMR_D->Sumw2();
  TH1F *hMR_Q = new TH1F("hMR_Q", "hMR_Q", nbin, xmin, xmax); hMR_Q->Sumw2();
  TH1F *hMR_T = new TH1F("hMR_T", "hMR_T", nbin, xmin, xmax); hMR_T->Sumw2();
  TH1F *hMR_W = new TH1F("hMR_W", "hMR_W", nbin, xmin, xmax); hMR_W->Sumw2();
  TH1F *hMR_Z = new TH1F("hMR_Z", "hMR_Z", nbin, xmin, xmax); hMR_Z->Sumw2();

  TH1F *hMR_pass_D = new TH1F("hMR_pass_D", "hMR_pass_D", nbin, xmin, xmax); hMR_pass_D->Sumw2();
  TH1F *hMR_pass_Q = new TH1F("hMR_pass_Q", "hMR_pass_Q", nbin, xmin, xmax); hMR_pass_Q->Sumw2();
  TH1F *hMR_pass_T = new TH1F("hMR_pass_T", "hMR_pass_T", nbin, xmin, xmax); hMR_pass_T->Sumw2();
  TH1F *hMR_pass_W = new TH1F("hMR_pass_W", "hMR_pass_W", nbin, xmin, xmax); hMR_pass_W->Sumw2();
  TH1F *hMR_pass_Z = new TH1F("hMR_pass_Z", "hMR_pass_Z", nbin, xmin, xmax); hMR_pass_Z->Sumw2();

  TH1F *hMR_fail_D = new TH1F("hMR_fail_D", "hMR_fail_D", nbin, xmin, xmax); hMR_fail_D->Sumw2();
  TH1F *hMR_fail_Q = new TH1F("hMR_fail_Q", "hMR_fail_Q", nbin, xmin, xmax); hMR_fail_Q->Sumw2();
  TH1F *hMR_fail_T = new TH1F("hMR_fail_T", "hMR_fail_T", nbin, xmin, xmax); hMR_fail_T->Sumw2();
  TH1F *hMR_fail_W = new TH1F("hMR_fail_W", "hMR_fail_W", nbin, xmin, xmax); hMR_fail_W->Sumw2();
  TH1F *hMR_fail_Z = new TH1F("hMR_fail_Z", "hMR_fail_Z", nbin, xmin, xmax); hMR_fail_Z->Sumw2();

  // draw Rsq
  nbin=17; xmin=0.15; xmax=1.0;

  TH1F *hRsq_D = new TH1F("hRsq_D", "hRsq_D", nbin, xmin, xmax); hRsq_D->Sumw2();
  TH1F *hRsq_Q = new TH1F("hRsq_Q", "hRsq_Q", nbin, xmin, xmax); hRsq_Q->Sumw2();
  TH1F *hRsq_T = new TH1F("hRsq_T", "hRsq_T", nbin, xmin, xmax); hRsq_T->Sumw2();
  TH1F *hRsq_W = new TH1F("hRsq_W", "hRsq_W", nbin, xmin, xmax); hRsq_W->Sumw2();
  TH1F *hRsq_Z = new TH1F("hRsq_Z", "hRsq_Z", nbin, xmin, xmax); hRsq_Z->Sumw2();

  TH1F *hRsq_pass_D = new TH1F("hRsq_pass_D", "hRsq_pass_D", nbin, xmin, xmax); hRsq_pass_D->Sumw2();
  TH1F *hRsq_pass_Q = new TH1F("hRsq_pass_Q", "hRsq_pass_Q", nbin, xmin, xmax); hRsq_pass_Q->Sumw2();
  TH1F *hRsq_pass_T = new TH1F("hRsq_pass_T", "hRsq_pass_T", nbin, xmin, xmax); hRsq_pass_T->Sumw2();
  TH1F *hRsq_pass_W = new TH1F("hRsq_pass_W", "hRsq_pass_W", nbin, xmin, xmax); hRsq_pass_W->Sumw2();
  TH1F *hRsq_pass_Z = new TH1F("hRsq_pass_Z", "hRsq_pass_Z", nbin, xmin, xmax); hRsq_pass_Z->Sumw2();

  TH1F *hRsq_fail_D = new TH1F("hRsq_fail_D", "hRsq_fail_D", nbin, xmin, xmax); hRsq_fail_D->Sumw2();
  TH1F *hRsq_fail_Q = new TH1F("hRsq_fail_Q", "hRsq_fail_Q", nbin, xmin, xmax); hRsq_fail_Q->Sumw2();
  TH1F *hRsq_fail_T = new TH1F("hRsq_fail_T", "hRsq_fail_T", nbin, xmin, xmax); hRsq_fail_T->Sumw2();
  TH1F *hRsq_fail_W = new TH1F("hRsq_fail_W", "hRsq_fail_W", nbin, xmin, xmax); hRsq_fail_W->Sumw2();
  TH1F *hRsq_fail_Z = new TH1F("hRsq_fail_Z", "hRsq_fail_Z", nbin, xmin, xmax); hRsq_fail_Z->Sumw2();

  // draw dPhiRazor
  nbin=8; 
  Float_t xbins[9] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.0, 3.2};

  TH1F *hDPhiR_D = new TH1F("hDPhiR_D", "hDPhiR_D", nbin, &xbins[0]); hDPhiR_D->Sumw2();
  TH1F *hDPhiR_Q = new TH1F("hDPhiR_Q", "hDPhiR_Q", nbin, &xbins[0]); hDPhiR_Q->Sumw2();
  TH1F *hDPhiR_T = new TH1F("hDPhiR_T", "hDPhiR_T", nbin, &xbins[0]); hDPhiR_T->Sumw2();
  TH1F *hDPhiR_W = new TH1F("hDPhiR_W", "hDPhiR_W", nbin, &xbins[0]); hDPhiR_W->Sumw2();
  TH1F *hDPhiR_Z = new TH1F("hDPhiR_Z", "hDPhiR_Z", nbin, &xbins[0]); hDPhiR_Z->Sumw2();

  TH1F *hDPhiR_mr1_D = new TH1F("hDPhiR_mr1_D", "hDPhiR_mr1_D", nbin, &xbins[0]); hDPhiR_mr1_D->Sumw2();
  TH1F *hDPhiR_mr1_Q = new TH1F("hDPhiR_mr1_Q", "hDPhiR_mr1_Q", nbin, &xbins[0]); hDPhiR_mr1_Q->Sumw2();
  TH1F *hDPhiR_mr1_T = new TH1F("hDPhiR_mr1_T", "hDPhiR_mr1_T", nbin, &xbins[0]); hDPhiR_mr1_T->Sumw2();
  TH1F *hDPhiR_mr1_W = new TH1F("hDPhiR_mr1_W", "hDPhiR_mr1_W", nbin, &xbins[0]); hDPhiR_mr1_W->Sumw2();
  TH1F *hDPhiR_mr1_Z = new TH1F("hDPhiR_mr1_Z", "hDPhiR_mr1_Z", nbin, &xbins[0]); hDPhiR_mr1_Z->Sumw2();
  
  TH1F *hDPhiR_mr2_D = new TH1F("hDPhiR_mr2_D", "hDPhiR_mr2_D", nbin, &xbins[0]); hDPhiR_mr2_D->Sumw2();
  TH1F *hDPhiR_mr2_Q = new TH1F("hDPhiR_mr2_Q", "hDPhiR_mr2_Q", nbin, &xbins[0]); hDPhiR_mr2_Q->Sumw2();
  TH1F *hDPhiR_mr2_T = new TH1F("hDPhiR_mr2_T", "hDPhiR_mr2_T", nbin, &xbins[0]); hDPhiR_mr2_T->Sumw2();
  TH1F *hDPhiR_mr2_W = new TH1F("hDPhiR_mr2_W", "hDPhiR_mr2_W", nbin, &xbins[0]); hDPhiR_mr2_W->Sumw2();
  TH1F *hDPhiR_mr2_Z = new TH1F("hDPhiR_mr2_Z", "hDPhiR_mr2_Z", nbin, &xbins[0]); hDPhiR_mr2_Z->Sumw2();
  
  TH1F *hDPhiR_mr3_D = new TH1F("hDPhiR_mr3_D", "hDPhiR_mr3_D", nbin, &xbins[0]); hDPhiR_mr3_D->Sumw2();
  TH1F *hDPhiR_mr3_Q = new TH1F("hDPhiR_mr3_Q", "hDPhiR_mr3_Q", nbin, &xbins[0]); hDPhiR_mr3_Q->Sumw2();
  TH1F *hDPhiR_mr3_T = new TH1F("hDPhiR_mr3_T", "hDPhiR_mr3_T", nbin, &xbins[0]); hDPhiR_mr3_T->Sumw2();
  TH1F *hDPhiR_mr3_W = new TH1F("hDPhiR_mr3_W", "hDPhiR_mr3_W", nbin, &xbins[0]); hDPhiR_mr3_W->Sumw2();
  TH1F *hDPhiR_mr3_Z = new TH1F("hDPhiR_mr3_Z", "hDPhiR_mr3_Z", nbin, &xbins[0]); hDPhiR_mr3_Z->Sumw2();
  
  TH1F *hDPhiR_mr4_D = new TH1F("hDPhiR_mr4_D", "hDPhiR_mr4_D", nbin, &xbins[0]); hDPhiR_mr4_D->Sumw2();
  TH1F *hDPhiR_mr4_Q = new TH1F("hDPhiR_mr4_Q", "hDPhiR_mr4_Q", nbin, &xbins[0]); hDPhiR_mr4_Q->Sumw2();
  TH1F *hDPhiR_mr4_T = new TH1F("hDPhiR_mr4_T", "hDPhiR_mr4_T", nbin, &xbins[0]); hDPhiR_mr4_T->Sumw2();
  TH1F *hDPhiR_mr4_W = new TH1F("hDPhiR_mr4_W", "hDPhiR_mr4_W", nbin, &xbins[0]); hDPhiR_mr4_W->Sumw2();
  TH1F *hDPhiR_mr4_Z = new TH1F("hDPhiR_mr4_Z", "hDPhiR_mr4_Z", nbin, &xbins[0]); hDPhiR_mr4_Z->Sumw2();
  
  TH1F *hDPhiR_mr5_D = new TH1F("hDPhiR_mr5_D", "hDPhiR_mr5_D", nbin, &xbins[0]); hDPhiR_mr5_D->Sumw2();
  TH1F *hDPhiR_mr5_Q = new TH1F("hDPhiR_mr5_Q", "hDPhiR_mr5_Q", nbin, &xbins[0]); hDPhiR_mr5_Q->Sumw2();
  TH1F *hDPhiR_mr5_T = new TH1F("hDPhiR_mr5_T", "hDPhiR_mr5_T", nbin, &xbins[0]); hDPhiR_mr5_T->Sumw2();
  TH1F *hDPhiR_mr5_W = new TH1F("hDPhiR_mr5_W", "hDPhiR_mr5_W", nbin, &xbins[0]); hDPhiR_mr5_W->Sumw2();
  TH1F *hDPhiR_mr5_Z = new TH1F("hDPhiR_mr5_Z", "hDPhiR_mr5_Z", nbin, &xbins[0]); hDPhiR_mr5_Z->Sumw2();

  // open trees
  TTree *tD = (TTree*) fD->Get("QCDTree");
  TTree *tQ = (TTree*) fQ->Get("QCDTree");
  TTree *tT = (TTree*) fT->Get("QCDTree");
  TTree *tW = (TTree*) fW->Get("QCDTree");
  TTree *tZ = (TTree*) fZ->Get("QCDTree");

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  // draw MR (again)
  tD->Draw("MR>>hMR_D", cut_str_dat);
  tQ->Draw("MR>>hMR_Q", cut_str);
  tT->Draw("MR>>hMR_T", cut_str);
  tW->Draw("MR>>hMR_W", cut_str);
  tZ->Draw("MR>>hMR_Z", cut_str);

  tD->Draw("MR>>hMR_pass_D", cut_str_dat+"*(abs(dPhiRazor)>2.8)");
  tQ->Draw("MR>>hMR_pass_Q", cut_str+"*(abs(dPhiRazor)>2.8)");
  tT->Draw("MR>>hMR_pass_T", cut_str+"*(abs(dPhiRazor)>2.8)");
  tW->Draw("MR>>hMR_pass_W", cut_str+"*(abs(dPhiRazor)>2.8)");
  tZ->Draw("MR>>hMR_pass_Z", cut_str+"*(abs(dPhiRazor)>2.8)");

  tD->Draw("MR>>hMR_fail_D", cut_str_dat+"*(abs(dPhiRazor)<2.8)");
  tQ->Draw("MR>>hMR_fail_Q", cut_str+"*(abs(dPhiRazor)<2.8)");
  tT->Draw("MR>>hMR_fail_T", cut_str+"*(abs(dPhiRazor)<2.8)");
  tW->Draw("MR>>hMR_fail_W", cut_str+"*(abs(dPhiRazor)<2.8)");
  tZ->Draw("MR>>hMR_fail_Z", cut_str+"*(abs(dPhiRazor)<2.8)");

  // draw Rsq (again)
  tD->Draw("Rsq>>hRsq_D", cut_str_dat);
  tQ->Draw("Rsq>>hRsq_Q", cut_str);
  tT->Draw("Rsq>>hRsq_T", cut_str);
  tW->Draw("Rsq>>hRsq_W", cut_str);
  tZ->Draw("Rsq>>hRsq_Z", cut_str);

  tD->Draw("Rsq>>hRsq_pass_D", cut_str_dat+"*(abs(dPhiRazor)>2.8)");
  tQ->Draw("Rsq>>hRsq_pass_Q", cut_str+"*(abs(dPhiRazor)>2.8)");
  tT->Draw("Rsq>>hRsq_pass_T", cut_str+"*(abs(dPhiRazor)>2.8)");
  tW->Draw("Rsq>>hRsq_pass_W", cut_str+"*(abs(dPhiRazor)>2.8)");
  tZ->Draw("Rsq>>hRsq_pass_Z", cut_str+"*(abs(dPhiRazor)>2.8)");

  tD->Draw("Rsq>>hRsq_fail_D", cut_str_dat+"*(abs(dPhiRazor)<2.8)");
  tQ->Draw("Rsq>>hRsq_fail_Q", cut_str+"*(abs(dPhiRazor)<2.8)");
  tT->Draw("Rsq>>hRsq_fail_T", cut_str+"*(abs(dPhiRazor)<2.8)");
  tW->Draw("Rsq>>hRsq_fail_W", cut_str+"*(abs(dPhiRazor)<2.8)");
  tZ->Draw("Rsq>>hRsq_fail_Z", cut_str+"*(abs(dPhiRazor)<2.8)");

  // draw dPhiR (again)
  tD->Draw("dPhiRazor>>hDPhiR_D", cut_str);
  tQ->Draw("dPhiRazor>>hDPhiR_Q", cut_str);
  tT->Draw("dPhiRazor>>hDPhiR_T", cut_str);
  tW->Draw("dPhiRazor>>hDPhiR_W", cut_str);
  tZ->Draw("dPhiRazor>>hDPhiR_Z", cut_str);

  tD->Draw("dPhiRazor>>hDPhiR_mr1_D", cut_str_dat+cut_str_mr1);
  tQ->Draw("dPhiRazor>>hDPhiR_mr1_Q", cut_str+cut_str_mr1);
  tT->Draw("dPhiRazor>>hDPhiR_mr1_T", cut_str+cut_str_mr1);
  tW->Draw("dPhiRazor>>hDPhiR_mr1_W", cut_str+cut_str_mr1);
  tZ->Draw("dPhiRazor>>hDPhiR_mr1_Z", cut_str+cut_str_mr1);

  tD->Draw("dPhiRazor>>hDPhiR_mr2_D", cut_str_dat+cut_str_mr2);
  tQ->Draw("dPhiRazor>>hDPhiR_mr2_Q", cut_str+cut_str_mr2);
  tT->Draw("dPhiRazor>>hDPhiR_mr2_T", cut_str+cut_str_mr2);
  tW->Draw("dPhiRazor>>hDPhiR_mr2_W", cut_str+cut_str_mr2);
  tZ->Draw("dPhiRazor>>hDPhiR_mr2_Z", cut_str+cut_str_mr2);

  tD->Draw("dPhiRazor>>hDPhiR_mr3_D", cut_str_dat+cut_str_mr3);
  tQ->Draw("dPhiRazor>>hDPhiR_mr3_Q", cut_str+cut_str_mr3);
  tT->Draw("dPhiRazor>>hDPhiR_mr3_T", cut_str+cut_str_mr3);
  tW->Draw("dPhiRazor>>hDPhiR_mr3_W", cut_str+cut_str_mr3);
  tZ->Draw("dPhiRazor>>hDPhiR_mr3_Z", cut_str+cut_str_mr3);

  tD->Draw("dPhiRazor>>hDPhiR_mr4_D", cut_str_dat+cut_str_mr4);
  tQ->Draw("dPhiRazor>>hDPhiR_mr4_Q", cut_str+cut_str_mr4);
  tT->Draw("dPhiRazor>>hDPhiR_mr4_T", cut_str+cut_str_mr4);
  tW->Draw("dPhiRazor>>hDPhiR_mr4_W", cut_str+cut_str_mr4);
  tZ->Draw("dPhiRazor>>hDPhiR_mr4_Z", cut_str+cut_str_mr4);

  tD->Draw("dPhiRazor>>hDPhiR_mr5_D", cut_str_dat+cut_str_mr5);
  tQ->Draw("dPhiRazor>>hDPhiR_mr5_Q", cut_str+cut_str_mr5);
  tT->Draw("dPhiRazor>>hDPhiR_mr5_T", cut_str+cut_str_mr5);
  tW->Draw("dPhiRazor>>hDPhiR_mr5_W", cut_str+cut_str_mr5);
  tZ->Draw("dPhiRazor>>hDPhiR_mr5_Z", cut_str+cut_str_mr5);

  //--------------------------------------------------------------
  //
  // configure draw options
  //
  //--------------------------------------------------------------

  //data
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

  // qcd
  InitHist(hMR_Q,"","",kMagenta);
  InitHist(hMR_pass_Q,"","",kMagenta);
  InitHist(hMR_fail_Q,"","",kMagenta);
  InitHist(hRsq_Q,"","",kMagenta);
  InitHist(hRsq_pass_Q,"","",kMagenta);
  InitHist(hRsq_fail_Q,"","",kMagenta);
  InitHist(hDPhiR_Q,"","",kMagenta);

  InitHist(hDPhiR_mr1_Q,"","",kMagenta);
  InitHist(hDPhiR_mr2_Q,"","",kMagenta);
  InitHist(hDPhiR_mr3_Q,"","",kMagenta);
  InitHist(hDPhiR_mr4_Q,"","",kMagenta);
  InitHist(hDPhiR_mr5_Q,"","",kMagenta);

  // tt
  InitHist(hMR_T,"","",kGreen+2);
  InitHist(hMR_pass_T,"","",kGreen+2);
  InitHist(hMR_fail_T,"","",kGreen+2);
  InitHist(hRsq_T,"","",kGreen+2);
  InitHist(hRsq_pass_T,"","",kGreen+2);
  InitHist(hRsq_fail_T,"","",kGreen+2);
  InitHist(hDPhiR_T,"","",kGreen+2);

  InitHist(hDPhiR_mr1_T,"","",kGreen+2);
  InitHist(hDPhiR_mr2_T,"","",kGreen+2);
  InitHist(hDPhiR_mr3_T,"","",kGreen+2);
  InitHist(hDPhiR_mr4_T,"","",kGreen+2);
  InitHist(hDPhiR_mr5_T,"","",kGreen+2);

  // z
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

  // w
  InitHist(hMR_W,"","",kRed+1);
  InitHist(hMR_pass_W,"","",kRed+1);
  InitHist(hMR_fail_W,"","",kRed+1);
  InitHist(hRsq_W,"","",kRed+1);
  InitHist(hRsq_pass_W,"","",kRed+1);
  InitHist(hRsq_fail_W,"","",kRed+1);
  InitHist(hDPhiR_W,"","",kRed+1);

  InitHist(hDPhiR_mr1_W,"","",kRed+1);
  InitHist(hDPhiR_mr2_W,"","",kRed+1);
  InitHist(hDPhiR_mr3_W,"","",kRed+1);
  InitHist(hDPhiR_mr4_W,"","",kRed+1);
  InitHist(hDPhiR_mr5_W,"","",kRed+1);

  hMR_Q->Add(hMR_T);
  hMR_Q->Add(hMR_W);
  hMR_Q->Add(hMR_Z);
  hMR_T->Add(hMR_W);
  hMR_T->Add(hMR_Z);
  hMR_W->Add(hMR_Z);

  Float_t scale=hMR_D->Integral()/hMR_Q->Integral();
  hMR_Q->Scale(scale);
  hMR_T->Scale(scale);
  hMR_W->Scale(scale);
  hMR_Z->Scale(scale);
  
  hMR_pass_Q->Add(hMR_pass_T);
  hMR_pass_Q->Add(hMR_pass_W);
  hMR_pass_Q->Add(hMR_pass_Z);
  hMR_pass_T->Add(hMR_pass_W);
  hMR_pass_T->Add(hMR_pass_Z);
  hMR_pass_W->Add(hMR_pass_Z);
  
  //scale=hMR_pass_D->Integral()/hMR_pass_Q->Integral();
  hMR_pass_Q->Scale(scale);
  hMR_pass_T->Scale(scale);
  hMR_pass_W->Scale(scale);
  hMR_pass_Z->Scale(scale);
  
  hMR_fail_Q->Add(hMR_fail_T);
  hMR_fail_Q->Add(hMR_fail_W);
  hMR_fail_Q->Add(hMR_fail_Z);
  hMR_fail_T->Add(hMR_fail_W);
  hMR_fail_T->Add(hMR_fail_Z);
  hMR_fail_W->Add(hMR_fail_Z);
  
  //scale=hMR_fail_D->Integral()/hMR_fail_Q->Integral();
  hMR_fail_Q->Scale(scale);
  hMR_fail_T->Scale(scale);
  hMR_fail_W->Scale(scale);
  hMR_fail_Z->Scale(scale);
  
  hRsq_Q->Add(hRsq_T);
  hRsq_Q->Add(hRsq_W);
  hRsq_Q->Add(hRsq_Z);
  hRsq_T->Add(hRsq_W);
  hRsq_T->Add(hRsq_Z);
  hRsq_W->Add(hRsq_Z);
  
  scale=hRsq_D->Integral()/hRsq_Q->Integral();
  hRsq_Q->Scale(scale);
  hRsq_T->Scale(scale);
  hRsq_W->Scale(scale);
  hRsq_Z->Scale(scale);

  hRsq_pass_Q->Add(hRsq_pass_T);
  hRsq_pass_Q->Add(hRsq_pass_W);
  hRsq_pass_Q->Add(hRsq_pass_Z);
  hRsq_pass_T->Add(hRsq_pass_W);
  hRsq_pass_T->Add(hRsq_pass_Z);
  hRsq_pass_W->Add(hRsq_pass_Z);

  //scale=hRsq_pass_D->Integral()/hRsq_pass_Q->Integral();
  hRsq_pass_Q->Scale(scale);
  hRsq_pass_T->Scale(scale);
  hRsq_pass_W->Scale(scale);
  hRsq_pass_Z->Scale(scale);

  hRsq_fail_Q->Add(hRsq_fail_T);
  hRsq_fail_Q->Add(hRsq_fail_W);
  hRsq_fail_Q->Add(hRsq_fail_Z);
  hRsq_fail_T->Add(hRsq_fail_W);
  hRsq_fail_T->Add(hRsq_fail_Z);
  hRsq_fail_W->Add(hRsq_fail_Z);

  //scale=hRsq_fail_D->Integral()/hRsq_fail_Q->Integral();
  hRsq_fail_Q->Scale(scale);
  hRsq_fail_T->Scale(scale);
  hRsq_fail_W->Scale(scale);
  hRsq_fail_Z->Scale(scale);

  hDPhiR_Q->Add(hDPhiR_T);
  hDPhiR_Q->Add(hDPhiR_W);
  hDPhiR_Q->Add(hDPhiR_Z);
  hDPhiR_T->Add(hDPhiR_W);
  hDPhiR_T->Add(hDPhiR_Z);
  hDPhiR_W->Add(hDPhiR_Z);

  scale=hDPhiR_D->Integral()/hDPhiR_Q->Integral();
  hDPhiR_Q->Scale(scale);
  hDPhiR_T->Scale(scale);
  hDPhiR_W->Scale(scale);
  hDPhiR_Z->Scale(scale);

  hDPhiR_mr1_Q->Add(hDPhiR_mr1_T);
  hDPhiR_mr1_Q->Add(hDPhiR_mr1_W);
  hDPhiR_mr1_Q->Add(hDPhiR_mr1_Z);
  hDPhiR_mr1_T->Add(hDPhiR_mr1_W);
  hDPhiR_mr1_T->Add(hDPhiR_mr1_Z);
  hDPhiR_mr1_W->Add(hDPhiR_mr1_Z);

  //scale=hDPhiR_mr1_D->Integral()/hDPhiR_mr1_Q->Integral();
  hDPhiR_mr1_Q->Scale(scale);
  hDPhiR_mr1_T->Scale(scale);
  hDPhiR_mr1_W->Scale(scale);
  hDPhiR_mr1_Z->Scale(scale);

  hDPhiR_mr2_Q->Add(hDPhiR_mr2_T);
  hDPhiR_mr2_Q->Add(hDPhiR_mr2_W);
  hDPhiR_mr2_Q->Add(hDPhiR_mr2_Z);
  hDPhiR_mr2_T->Add(hDPhiR_mr2_W);
  hDPhiR_mr2_T->Add(hDPhiR_mr2_Z);
  hDPhiR_mr2_W->Add(hDPhiR_mr2_Z);

  //scale=hDPhiR_mr2_D->Integral()/hDPhiR_mr2_Q->Integral();
  hDPhiR_mr2_Q->Scale(scale);
  hDPhiR_mr2_T->Scale(scale);
  hDPhiR_mr2_W->Scale(scale);
  hDPhiR_mr2_Z->Scale(scale);
  
  hDPhiR_mr3_Q->Add(hDPhiR_mr3_T);
  hDPhiR_mr3_Q->Add(hDPhiR_mr3_W);
  hDPhiR_mr3_Q->Add(hDPhiR_mr3_Z);
  hDPhiR_mr3_T->Add(hDPhiR_mr3_W);
  hDPhiR_mr3_T->Add(hDPhiR_mr3_Z);
  hDPhiR_mr3_W->Add(hDPhiR_mr3_Z);
  
  //scale=hDPhiR_mr3_D->Integral()/hDPhiR_mr3_Q->Integral();
  hDPhiR_mr3_Q->Scale(scale);
  hDPhiR_mr3_T->Scale(scale);
  hDPhiR_mr3_W->Scale(scale);
  hDPhiR_mr3_Z->Scale(scale);

  hDPhiR_mr4_Q->Add(hDPhiR_mr4_T);
  hDPhiR_mr4_Q->Add(hDPhiR_mr4_W);
  hDPhiR_mr4_Q->Add(hDPhiR_mr4_Z);
  hDPhiR_mr4_T->Add(hDPhiR_mr4_W);
  hDPhiR_mr4_T->Add(hDPhiR_mr4_Z);
  hDPhiR_mr4_W->Add(hDPhiR_mr4_Z);

  //scale=hDPhiR_mr4_D->Integral()/hDPhiR_mr4_Q->Integral();
  hDPhiR_mr4_Q->Scale(scale);
  hDPhiR_mr4_T->Scale(scale);
  hDPhiR_mr4_W->Scale(scale);
  hDPhiR_mr4_Z->Scale(scale);
  
  hDPhiR_mr5_Q->Add(hDPhiR_mr5_T);
  hDPhiR_mr5_Q->Add(hDPhiR_mr5_W);
  hDPhiR_mr5_Q->Add(hDPhiR_mr5_Z);
  hDPhiR_mr5_T->Add(hDPhiR_mr5_W);
  hDPhiR_mr5_T->Add(hDPhiR_mr5_Z);
  hDPhiR_mr5_W->Add(hDPhiR_mr5_Z);

  //scale=hDPhiR_mr5_D->Integral()/hDPhiR_mr5_Q->Integral();
  hDPhiR_mr5_Q->Scale(scale);
  hDPhiR_mr5_T->Scale(scale);
  hDPhiR_mr5_W->Scale(scale);
  hDPhiR_mr5_Z->Scale(scale);

  TLegend *leg = new TLegend(0.67,0.57,0.89,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  leg->AddEntry(hMR_D, "Data", "pel");
  leg->AddEntry(hMR_Q, "QCD", "f");
  leg->AddEntry(hMR_T, "TT+jets", "f");
  leg->AddEntry(hMR_W, "W+jets", "f");
  leg->AddEntry(hMR_Z, "Zvv+jets", "f");  

  //--------------------------------------------------------------
  //
  // draw
  //
  //--------------------------------------------------------------

  hMR_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_Q->GetMaximum(), hMR_D->GetMaximum()));
  hMR_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_Q->GetYaxis()->SetTitle("Events");
  hMR_Q->SetTitle("");
  hMR_Q->Draw("hist");
  hMR_T->Draw("histsame");
  hMR_W->Draw("histsame");
  hMR_Z->Draw("histsame");
  hMR_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_dijet.png");

  hMR_pass_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_pass_Q->GetMaximum(), hMR_pass_D->GetMaximum()));
  hMR_pass_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_pass_Q->GetYaxis()->SetTitle("Events");
  hMR_pass_Q->SetTitle("");
  hMR_pass_Q->Draw("hist");
  hMR_pass_T->Draw("histsame");
  hMR_pass_W->Draw("histsame");
  hMR_pass_Z->Draw("histsame");
  hMR_pass_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_pass_dijet.png");

  hMR_fail_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_fail_Q->GetMaximum(), hMR_fail_D->GetMaximum()));
  hMR_fail_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_fail_Q->GetYaxis()->SetTitle("Events");
  hMR_fail_Q->SetTitle("");
  hMR_fail_Q->Draw("hist");
  hMR_fail_T->Draw("histsame");
  hMR_fail_W->Draw("histsame");
  hMR_fail_Z->Draw("histsame");
  hMR_fail_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_fail_dijet.png");

  c->SetLogy(1);
  hRsq_Q->GetYaxis()->SetRangeUser(0.001,1.2*TMath::Max(hRsq_Q->GetMaximum(), hRsq_D->GetMaximum()));
  hRsq_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_Q->GetYaxis()->SetTitle("Events");
  hRsq_Q->SetTitle("");
  hRsq_Q->Draw("hist");
  hRsq_T->Draw("histsame");
  hRsq_W->Draw("histsame");
  hRsq_Z->Draw("histsame");
  hRsq_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_dijet.png");

  hRsq_pass_Q->GetYaxis()->SetRangeUser(0.001,1.2*TMath::Max(hRsq_pass_Q->GetMaximum(), hRsq_pass_D->GetMaximum()));
  hRsq_pass_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_pass_Q->GetYaxis()->SetTitle("Events");
  hRsq_pass_Q->SetTitle("");
  hRsq_pass_Q->Draw("hist");
  hRsq_pass_T->Draw("histsame");
  hRsq_pass_W->Draw("histsame");
  hRsq_pass_Z->Draw("histsame");
  hRsq_pass_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_pass_dijet.png");

  hRsq_fail_Q->GetYaxis()->SetRangeUser(0.001,1.2*TMath::Max(hRsq_fail_Q->GetMaximum(), hRsq_fail_D->GetMaximum()));
  hRsq_fail_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_fail_Q->GetYaxis()->SetTitle("Events");
  hRsq_fail_Q->SetTitle("");
  hRsq_fail_Q->Draw("hist");
  hRsq_fail_T->Draw("histsame");
  hRsq_fail_W->Draw("histsame");
  hRsq_fail_Z->Draw("histsame");
  hRsq_fail_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_fail_dijet.png");
  c->SetLogy(0);
  leg->SetX1NDC(0.39); leg->SetX2NDC(0.61);

  hDPhiR_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_Q->GetMaximum(), hDPhiR_D->GetMaximum()));
  hDPhiR_Q->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_Q->SetTitle("");
  hDPhiR_Q->Draw("hist");
  hDPhiR_T->Draw("histsame");
  hDPhiR_W->Draw("histsame");
  hDPhiR_Z->Draw("histsame");
  hDPhiR_D->Draw("same e");

  leg->Draw();

  c->SaveAs("DPhiR_dijet.png");

  hDPhiR_mr1_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr1_Q->GetMaximum(), hDPhiR_mr1_D->GetMaximum()));
  hDPhiR_mr1_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr1_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr1_Q->SetTitle("500 < MR < 600");
  hDPhiR_mr1_Q->Draw("hist");
  hDPhiR_mr1_T->Draw("histsame");
  hDPhiR_mr1_W->Draw("histsame");
  hDPhiR_mr1_Z->Draw("histsame");
  hDPhiR_mr1_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr1_dijet.png");

  hDPhiR_mr2_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr2_Q->GetMaximum(), hDPhiR_mr2_D->GetMaximum()));
  hDPhiR_mr2_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr2_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr2_Q->SetTitle("600 < MR < 700");
  hDPhiR_mr2_Q->Draw("hist");
  hDPhiR_mr2_T->Draw("histsame");
  hDPhiR_mr2_W->Draw("histsame");
  hDPhiR_mr2_Z->Draw("histsame");
  hDPhiR_mr2_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr2_dijet.png");
  
  hDPhiR_mr3_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr3_Q->GetMaximum(), hDPhiR_mr3_D->GetMaximum()));
  hDPhiR_mr3_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr3_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr3_Q->SetTitle("700 < MR < 800");
  hDPhiR_mr3_Q->Draw("hist");
  hDPhiR_mr3_T->Draw("histsame");
  hDPhiR_mr3_W->Draw("histsame");
  hDPhiR_mr3_Z->Draw("histsame");
  hDPhiR_mr3_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr3_dijet.png");
  
  hDPhiR_mr4_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr4_Q->GetMaximum(), hDPhiR_mr4_D->GetMaximum()));
  hDPhiR_mr4_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr4_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr4_Q->SetTitle("800 < MR < 900");
  hDPhiR_mr4_Q->Draw("hist");
  hDPhiR_mr4_T->Draw("histsame");
  hDPhiR_mr4_W->Draw("histsame");
  hDPhiR_mr4_Z->Draw("histsame");
  hDPhiR_mr4_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr4_dijet.png");
  
  hDPhiR_mr5_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr5_Q->GetMaximum(), hDPhiR_mr5_D->GetMaximum()));
  hDPhiR_mr5_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr5_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr5_Q->SetTitle("900 < MR < 1000");
  hDPhiR_mr5_Q->Draw("hist");
  hDPhiR_mr5_T->Draw("histsame");
  hDPhiR_mr5_W->Draw("histsame");
  hDPhiR_mr5_Z->Draw("histsame");
  hDPhiR_mr5_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr5_dijet.png");

}
