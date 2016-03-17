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
#include <TAxis.h>
#include <TMath.h>
#include <TH3F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
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

void makepaperplots() {

  enum {mr=0, ljpt, rsq};
  enum {razor=0, dijet};
  //--------------------------------------------------------------
  //
  // RAZOR-TRIGGERED region setup
  //
  //--------------------------------------------------------------

  // for mr
  Int_t binIn=mr;
  const Int_t nbinx=1, nbiny=7, nbinz=2;
  Float_t xmin=0.15, xmax=0.25; 
  Float_t ymin=400,  ymax=3000; 
  Float_t zmin=0,    zmax=2;
  Float_t xbins[nbinx+1] = {xmin, xmax};
  Float_t ybins[nbiny+1] = {ymin, 500, 600, 700, 800, 900, 1000, ymax};
  Float_t zbins[nbinz+1] = {zmin, 1, zmax};

  TString pname = "npf_vs_mr_razor_fit.pdf";

  // for rsq
  //Int_t binIn=rsq;
  ////const Int_t nbinx=1, nbiny=10, nbinz=2;
  //const Int_t nbinx=1, nbiny=8, nbinz=2;
  //Float_t xmin=400, xmax=3000; 
  //Float_t ymin=0.15,  ymax=0.35; 
  //Float_t zmin=0,    zmax=2;
  //Float_t xbins[nbinx+1] = {xmin, xmax};
  ////Float_t ybins[nbiny+1] = { ymin, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, ymax};
  //Float_t ybins[nbiny+1] = { ymin, 0.16, 0.17, 0.18, 0.19, 0.20, 0.225, 0.25, ymax};
  //Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //TString pname = "npf_vs_rsq_razor.png";

  TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/HTMHT_Run2015D_Golden.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/QCD_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_skim.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/ZJetsToNuNu_13TeV-madgraph_2137pb_skim.root","read");
  
  TString cut_str="weight*(box==11||box==12)*puWeight*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)*(Rsq<0.25)";
  TString cut_str_dat=cut_str+"*(Rsq<0.25)";

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  TH3F *dPhiPF_Q = new TH3F("dPhiPF_Q", "dPhiPF_Q", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Q->Sumw2();
  TH3F *dPhiPF_D = new TH3F("dPhiPF_D", "dPhiPF_D", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_D->Sumw2();
  TH3F *dPhiPF_T = new TH3F("dPhiPF_T", "dPhiPF_T", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_T->Sumw2();
  TH3F *dPhiPF_W = new TH3F("dPhiPF_W", "dPhiPF_W", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_W->Sumw2();
  TH3F *dPhiPF_Z = new TH3F("dPhiPF_Z", "dPhiPF_Z", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Z->Sumw2();

  Float_t ybins2[nbiny+3] = {ymin-50, ymin, 500, 600, 700, 800, 900, 1000, ymax, ymax+50};
  TH1F *fxn_plus_err = new TH1F("fxn_plus_err","fxn_plus_err", nbiny+2, &ybins2[0]); fxn_plus_err->Sumw2();

  // open trees
  TTree *tQ = (TTree*) fQ->Get("QCDTree");
  TTree *tD = (TTree*) fD->Get("QCDTree");
  TTree *tT = (TTree*) fT->Get("QCDTree");
  TTree *tW = (TTree*) fW->Get("QCDTree");
  TTree *tZ = (TTree*) fZ->Get("QCDTree");

  TCanvas *c = MakeCanvas("c","c",800,600);

  if (binIn==mr) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_Q", cut_str);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_T", cut_str);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_W", cut_str);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_Z", cut_str);
  } else if (binIn==rsq) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Q", cut_str);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_T", cut_str);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_W", cut_str);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Z", cut_str);
  }

  TF1 *qcd_fxn = new TF1("qcd_fxn","[0]*x^[1]+[2]", 400, 3000);
  qcd_fxn->SetParameter(0,3.1e7);
  qcd_fxn->SetParameter(1,-3.1);
  qcd_fxn->SetParameter(2,0.062);

  fxn_plus_err->SetBinContent(1,qcd_fxn->Eval(ymin-25));
  fxn_plus_err->SetBinError(1,qcd_fxn->Eval(ymin-25)*0.87);

  fxn_plus_err->SetBinContent(nbiny+2,qcd_fxn->Eval(ymax+25));
  fxn_plus_err->SetBinError(nbiny+2,qcd_fxn->Eval(ymax+25)*0.87);

  Double_t wtf=0, lesswtf=0;

  vector<TGraphAsymmErrors *> qcd_with_mr; 
  for (Int_t i=0; i<dPhiPF_Q->GetNbinsX(); i++) {
    qcd_with_mr.push_back(new TGraphAsymmErrors());

    Int_t k=0;
  
    for (Int_t j=0; j<dPhiPF_Q->GetNbinsY(); j++) {
  
      Float_t rsq=dPhiPF_Q->GetYaxis()->GetBinCenter(j+1);
  
      Float_t dnP=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1));
      Float_t dnF=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,2))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
  
      Float_t nPF=0;
      if (dnP>0 && dnF>0) {
	nPF=dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
	wtf+=qcd_fxn->Eval(rsq)*(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));
	lesswtf+=(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));

	Float_t drsq=0.5*dPhiPF_Q->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);

	fxn_plus_err->SetBinContent(k+2,qcd_fxn->Eval(rsq));
	fxn_plus_err->SetBinError(k+2,qcd_fxn->Eval(rsq)*0.87);
  
	qcd_with_mr[i]->SetPoint(k, rsq, nPF);
	qcd_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
  }

  vector<TGraphAsymmErrors *> dat_with_mr; 
  Float_t xrsq=0, drsq=0, dnP=0, dnF=0, dnP_T=0, dnF_T=0, dnP_W=0, dnF_W=0, dnP_Z=0, dnF_Z=0, nPass=0, nFail=0, nPF=0;
  for (Int_t i=0; i<dPhiPF_D->GetNbinsX(); i++) {
    dat_with_mr.push_back(new TGraphAsymmErrors());

    Int_t k=0;

    for (Int_t j=0; j<dPhiPF_D->GetNbinsY(); j++) {
      xrsq=0; drsq=0;
      dnP=0; dnF=0; 
      dnP_T=0; dnF_T=0; 
      dnP_W=0; dnF_W=0; 
      dnP_Z=0; dnF_Z=0; 
      nPass=0; nFail=0; nPF=0;

      xrsq=dPhiPF_D->GetYaxis()->GetBinCenter(j+1);

      dnP=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(i+1,j+1,1))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1));
      dnF=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(i+1,j+1,2))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2));

      dnP=dnP*dnP;
      dnF=dnF*dnF;

      dnP_T=dPhiPF_T->GetBinError(dPhiPF_T->GetBin(i+1,j+1,1))/dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1));
      dnF_T=dPhiPF_T->GetBinError(dPhiPF_T->GetBin(i+1,j+1,2))/dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2));

      dnP_W=dPhiPF_W->GetBinError(dPhiPF_W->GetBin(i+1,j+1,1))/dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1));
      dnF_W=dPhiPF_W->GetBinError(dPhiPF_W->GetBin(i+1,j+1,2))/dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2));

      dnP_Z=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,1))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      dnF_Z=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,2))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));

      nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      
      nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));
      
      dnP+=dnP_T*dnP_T+dnP_W*dnP_W+dnP_Z*dnP_Z;
      dnF+=dnF_T*dnF_T+dnF_W*dnF_W+dnP_Z*dnP_Z;

      if (nPass>0 && nFail>0) {
	nPF=nPass/nFail;

	Float_t drsq=0.5*dPhiPF_D->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP + dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP + dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
	
	dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
  }


  Int_t i=0;

  qcd_fxn->SetLineColor(kBlue+1);

  qcd_with_mr[i]->SetMarkerStyle(24);
  //qcd_with_mr[i]->SetMarkerStyle(20);
  //qcd_with_mr[i]->SetLineColor(kRed);
  //qcd_with_mr[i]->SetMarkerColor(kRed);
  dat_with_mr[i]->SetMarkerStyle(21);
  
  if (binIn==mr) qcd_with_mr[i]->GetXaxis()->SetTitle("M_{R} [ GeV ]");  
  else if (binIn==rsq) qcd_with_mr[i]->GetXaxis()->SetTitle("R^{2}");
  qcd_with_mr[i]->GetYaxis()->SetTitle("Translation Factor #zeta");
  qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.0);
  if (binIn==mr) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 0.5);
  qcd_with_mr[i]->GetXaxis()->SetNdivisions(508);
  qcd_with_mr[i]->GetYaxis()->SetNdivisions(508);
  qcd_with_mr[i]->SetTitle("");
  qcd_with_mr[i]->Draw("ap e1");
  dat_with_mr[i]->Draw("p e1 same");

  fxn_plus_err->SetFillColor(kAzure+7);
  fxn_plus_err->SetFillStyle(3254);
  fxn_plus_err->SetMarkerStyle(0);
  fxn_plus_err->SetLineColor(kBlue+1);

  fxn_plus_err->Draw("same f e3");
  qcd_with_mr[i]->Draw("same p e1");
  dat_with_mr[i]->Draw("same p e1");
  qcd_fxn->Draw("same l");


  TLegend *leg = new TLegend(0.50,0.70,0.75,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);

  leg->AddEntry(dat_with_mr[0], "Data-(t#bar{t}+EWK)", "pel");
  leg->AddEntry(qcd_with_mr[0], "QCD MC", "pel");
  leg->AddEntry(fxn_plus_err, "Functional Form", "lf");


  leg->Draw();
  c->SaveAs(pname);

}
