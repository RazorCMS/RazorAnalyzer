#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

void plotScaleFactor() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","READ");

  TH2F *ttbarNominal = (TH2F*)inf->Get("TTJetsScaleFactors");
  TH2F *ttbarUp = (TH2F*)inf->Get("TTJetsScaleFactorsUp");
  TH2F *ttbarDown = (TH2F*)inf->Get("TTJetsScaleFactorsDown");
  TH2F *wNominal = (TH2F*)inf->Get("WJetsScaleFactors");
  TH2F *wUp = (TH2F*)inf->Get("WJetsScaleFactorsUp");
  TH2F *wDown = (TH2F*)inf->Get("WJetsScaleFactorsDown");
  TH2F *wInvNominal = (TH2F*)inf->Get("WJetsInvScaleFactors");
  TH2F *wInvUp = (TH2F*)inf->Get("WJetsInvScaleFactorsUp");
  TH2F *wInvDown = (TH2F*)inf->Get("WJetsInvScaleFactorsDown");
  TH2F *GJetInvNominal = (TH2F*)inf->Get("GJetsInvScaleFactors");


  TCanvas *cv = 0;
  gStyle->SetPaintTextFormat("4.2f");

  //****************************************************
  //Plot GJetsInv Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  GJetInvNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  GJetInvNominal->GetXaxis()->SetRangeUser(400,4000);
  GJetInvNominal->GetYaxis()->SetRangeUser(0.25,1.5);
  GJetInvNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  GJetInvNominal->GetZaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetZaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetXaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetXaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetXaxis()->SetTitleOffset(0.8);
  GJetInvNominal->GetYaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetYaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetYaxis()->SetTitleOffset(0.8);
  GJetInvNominal->SetStats(false);
  GJetInvNominal->SetMaximum(1.8);
  GJetInvNominal->SetMinimum(0.35);


  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("GJetsInvScaleFactor_CorrectedToMultiJet.png");
  cv->SaveAs("GJetsInvScaleFactor_CorrectedToMultiJet.pdf");


  //****************************************************
  //Plot WJetsInv Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  wInvNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  wInvNominal->GetXaxis()->SetRangeUser(400,4000);
  wInvNominal->GetYaxis()->SetRangeUser(0.25,1.5);
  wInvNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  wInvNominal->GetZaxis()->SetLabelSize(0.05);
  wInvNominal->GetZaxis()->SetTitleSize(0.05);
  wInvNominal->GetXaxis()->SetLabelSize(0.05);
  wInvNominal->GetXaxis()->SetTitleSize(0.05);
  wInvNominal->GetXaxis()->SetTitleOffset(0.8);
  wInvNominal->GetYaxis()->SetLabelSize(0.05);
  wInvNominal->GetYaxis()->SetTitleSize(0.05);
  wInvNominal->GetYaxis()->SetTitleOffset(0.8);
  wInvNominal->SetStats(false);
  wInvNominal->SetMaximum(1.8);
  wInvNominal->SetMinimum(0.35);


  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("WJetsInvScaleFactor_CorrectedToMultiJet.png");
  cv->SaveAs("WJetsInvScaleFactor_CorrectedToMultiJet.pdf");




 //****************************************************
  //Plot WJets Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  wNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  wNominal->GetXaxis()->SetRangeUser(400,4000);
  wNominal->GetYaxis()->SetRangeUser(0.25,1.5);
  wNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  wNominal->GetZaxis()->SetLabelSize(0.05);
  wNominal->GetZaxis()->SetTitleSize(0.05);
  wNominal->GetXaxis()->SetLabelSize(0.05);
  wNominal->GetXaxis()->SetTitleSize(0.05);
  wNominal->GetXaxis()->SetTitleOffset(0.8);
  wNominal->GetYaxis()->SetLabelSize(0.05);
  wNominal->GetYaxis()->SetTitleSize(0.05);
  wNominal->GetYaxis()->SetTitleOffset(0.8);
  wNominal->SetStats(false);
  wNominal->SetMaximum(1.8);
  wNominal->SetMinimum(0.35);

  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("WJetsScaleFactor_CorrectedToMultiJet.png");
  cv->SaveAs("WJetsScaleFactor_CorrectedToMultiJet.pdf");



 //****************************************************
  //Plot TTBar Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  ttbarNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  ttbarNominal->GetXaxis()->SetRangeUser(400,4000);
  ttbarNominal->GetYaxis()->SetRangeUser(0.25,1.5);
  ttbarNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  ttbarNominal->GetZaxis()->SetLabelSize(0.05);
  ttbarNominal->GetZaxis()->SetTitleSize(0.05);
  ttbarNominal->GetXaxis()->SetLabelSize(0.05);
  ttbarNominal->GetXaxis()->SetTitleSize(0.05);
  ttbarNominal->GetXaxis()->SetTitleOffset(0.8);
  ttbarNominal->GetYaxis()->SetLabelSize(0.05);
  ttbarNominal->GetYaxis()->SetTitleSize(0.05);
  ttbarNominal->GetYaxis()->SetTitleOffset(0.8);
  ttbarNominal->SetStats(false);
  // ttbarNominal->SetMaximum(10000);
  // ttbarNominal->SetMinimum(0.0);
  ttbarNominal->SetMaximum(1.8);
  ttbarNominal->SetMinimum(0.35);

  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("TTBarScaleFactor_CorrectedToMultiJet.png");
  cv->SaveAs("TTBarScaleFactor_CorrectedToMultiJet.pdf");






}




void plotGJetsScaleFactorSystematics() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","READ");

  TH2F *ttbarNominal = (TH2F*)inf->Get("TTJetsScaleFactors");
  TH2F *ttbarUp = (TH2F*)inf->Get("TTJetsScaleFactorsUp");
  TH2F *ttbarDown = (TH2F*)inf->Get("TTJetsScaleFactorsDown");
  TH2F *wNominal = (TH2F*)inf->Get("WJetsScaleFactors");
  TH2F *wUp = (TH2F*)inf->Get("WJetsScaleFactorsUp");
  TH2F *wDown = (TH2F*)inf->Get("WJetsScaleFactorsDown");
  TH2F *wInvNominal = (TH2F*)inf->Get("WJetsInvScaleFactors");
  TH2F *wInvUp = (TH2F*)inf->Get("WJetsInvScaleFactorsUp");
  TH2F *wInvDown = (TH2F*)inf->Get("WJetsInvScaleFactorsDown");
  TH2F *GJetInvNominal = (TH2F*)inf->Get("GJetsInvScaleFactors");


  TCanvas *cv = 0;
  gStyle->SetPaintTextFormat("4.2f");

  //****************************************************
  //Systematic Uncertainty
  //****************************************************
  TH2F *GJetsSystematicUnc = (TH2F*)GJetInvNominal->Clone("GJetsSystematicUnc");
  for ( int i = 1; i<GJetsSystematicUnc->GetXaxis()->GetNbins()+1; ++i) {
    for ( int j = 1; j<GJetsSystematicUnc->GetYaxis()->GetNbins()+1; ++j) {
      double gjet = GJetInvNominal->GetBinContent(i,j);
      double wjet = wInvNominal->GetBinContent(wInvNominal->GetXaxis()->FindFixBin(GJetInvNominal->GetXaxis()->GetBinCenter(i)),
					       wInvNominal->GetYaxis()->FindFixBin(GJetInvNominal->GetYaxis()->GetBinCenter(j)));            
      GJetsSystematicUnc->SetBinContent(i,j, fabs(gjet - wjet)/gjet );
    }
  }

  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(1);
  GJetsSystematicUnc->Draw("colztext");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  GJetsSystematicUnc->GetXaxis()->SetRangeUser(400,4000);
  GJetsSystematicUnc->GetYaxis()->SetRangeUser(0.25,1.5);
  GJetsSystematicUnc->GetZaxis()->SetTitle("Systematic Uncertainty");
  GJetsSystematicUnc->GetZaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetZaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetTitleOffset(0.8);
  GJetsSystematicUnc->GetYaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetYaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetYaxis()->SetTitleOffset(0.8);
  GJetsSystematicUnc->SetStats(false);
  GJetsSystematicUnc->SetMaximum(1.0);
  GJetsSystematicUnc->SetMinimum(0.0);

  lumi_13TeV = "2.2 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("GJetsVsWJetsSystematic.png");
  cv->SaveAs("GJetsVsWJetsSystematic.pdf");


  //****************************************************
  //Make Up and Down Scale Factor Histogram
  //****************************************************
  TH2F *GJetsScaleFactor_Down = (TH2F*)GJetInvNominal->Clone("GJetsInvScaleFactors_Down");
  for ( int i = 1; i<GJetsScaleFactor_Down->GetXaxis()->GetNbins()+1; ++i) {
    for ( int j = 1; j<GJetsScaleFactor_Down->GetYaxis()->GetNbins()+1; ++j) {
      double gjet = GJetInvNominal->GetBinContent(i,j);
      double wjet = wInvNominal->GetBinContent(wInvNominal->GetXaxis()->FindFixBin(GJetInvNominal->GetXaxis()->GetBinCenter(i)),
					       wInvNominal->GetYaxis()->FindFixBin(GJetInvNominal->GetYaxis()->GetBinCenter(j)));            
      GJetsScaleFactor_Down->SetBinContent(i,j, gjet - (wjet - gjet) );
      cout << "Bin " << i << " " << j << " : " << gjet << " , " << wjet << " , " <<  gjet - (wjet - gjet) << "\n";
    }
  }

  TFile *outf = new TFile("RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","UPDATE");
  outf->WriteTObject(ttbarNominal);
  outf->WriteTObject(ttbarUp);
  outf->WriteTObject(ttbarDown);
  outf->WriteTObject(wNominal);
  outf->WriteTObject(wUp);
  outf->WriteTObject(wDown);
  outf->WriteTObject(wInvNominal);
  outf->WriteTObject(wInvUp);
  outf->WriteTObject(wInvDown);
  outf->WriteTObject(GJetInvNominal);
  outf->WriteTObject(GJetsScaleFactor_Down);
  outf->Close();
  

}



void plotScaleFactorHistograms() {

  //plotScaleFactor();
  plotGJetsScaleFactorSystematics();

}
