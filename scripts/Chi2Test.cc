#include <iostream>
//ROOT INCLUDES
#include <TTree.h>
#include <TLatex.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBox.h>
#include <TSystem.h>
//ROOFIT INCLUDES
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooExtendPdf.h>
#include <RooStats/SPlot.h>
#include <RooStats/ModelConfig.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooBernstein.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooRandom.h>

void Chi2Test()
{
  TFile* f = new TFile("../BinnedFitResults_MultiJet.root", "READ");
  
  //gSystem->Load("../python/lib/libRazorRun2.so");
  RooWorkspace* ws = (RooWorkspace*)f->Get("wMultiJet");
  ws->Print();
  

  //RooRealVar mr("MR","M_{R}", 100, 4000,"GeV");
  RooRealVar* th1x = ws->var("th1x");
  
  RooRealVar* mr = ws->var("MR");
  mr->setBins(20);
  mr->setRange("low", 103, 120);
  mr->setRange("high", 135, 160);
  mr->setRange("Full", 100, 4000);

  RooRealVar* rsq = ws->var("Rsq");
  rsq->setBins(20);
  rsq->setRange("low", 0, 2);
  rsq->setRange("high", 0, 2);
  rsq->setRange("Full", 0, 2);
  
  RooRealVar* nBtag = ws->var("nBtag");
  nBtag->setBins(20);
  //nBtag->setRange("low", 0, 10);
  //nBtag->setRange("high", 0, 10);
  //nBtag->setRange("Full", 0, 10);
  
  RooDataSet* dataset = (RooDataSet*)ws->data("RMRTree");
  
  RooDataHist* dhist  = (RooDataHist*)ws->data("data_obs");
  //RooDataSet*  dset2  = (RooDataSet*)dhist->reduce(" nBtag == 0 ");
  //RooDataHist* dhist2 = (RooDataHist*)dhist->binnedClone();

  
  RooAbsPdf* pdf_0b = ws->pdf("MultiJet_TTj0b");
  RooAbsPdf* pdf_0b_copy = ws->pdf("MultiJet_TTj0b");
  RooRealVar* Ntot_TTj0b_MultiJet = ws->var("Ntot_TTj0b_MultiJet");
  RooAddPdf* my_0b_pdf = new RooAddPdf( "my_0b_pdf","", RooArgList(*pdf_0b), RooArgList(*Ntot_TTj0b_MultiJet) );
  RooDataHist* gdhist;

  TH1F* hdata;
  TH1F* hpdf;
  TH1F* h_chi2 = new TH1F("h_chi2", "chi2", 200, -100, 100);
  
  TTree* outTree = new TTree("outTree", "outTree of chi2");
  double chi2;
  outTree->Branch("chi2", &chi2, "chi2/F");
  
  for ( int i = 0; i < 10000; i++ )
    {
      gdhist = pdf_0b->generateBinned(*th1x, 1000, kFALSE);
      RooFitResult* bres = my_0b_pdf->fitTo( *gdhist, RooFit::Strategy(0), RooFit::Extended(kTRUE), RooFit::Save(kTRUE) );
      hdata = (TH1F*)gdhist->createHistogram("gdth1", *th1x);
      hpdf = (TH1F*)my_0b_pdf->createHistogram("pdf_th1", *th1x);
      chi2 = 0;
      for ( int j = 1; j <= hdata->GetNbinsX(); j++ )
	{
	  if ( hpdf->GetBinContent(j) > 0 ) chi2 += pow(hdata->GetBinContent(j)-hpdf->GetBinContent(j),2)/hpdf->GetBinContent(j);
	  //if ( hpdf->GetBinContent(j) > 0 ) chi2 += (hdata->GetBinContent(j)-hpdf->GetBinContent(j))/sqrt(hpdf->GetBinContent(j));
	  //std::cout << "-----------" << chi2 << std::endl;
	}
      h_chi2->Fill( chi2 );
      outTree->Fill( );
    }

  //RooAbsPdf* pdf_0b = ws->pdf("MultiJet_TTj0b");
  RooPlot *fmr = th1x->frame();
  //dataset->plotOn(fmr);
  gdhist->plotOn(fmr);
  //my_0b_pdf->plotOn(fmr);
  //pdf_0b->plotOn(fmr, RooFit::LineColor(kRed));
  //fmr->Draw();

  //hdata->Draw();
  hpdf->SetLineColor(kRed);
  //hpdf->Draw("same");
  h_chi2->Draw();

  TFile* fout = new TFile("outTree.root", "RECREATE");
  outTree->Write();
  fout->Close();
}
