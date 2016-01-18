void doCopy() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_TTJets.root","READ");
  TFile *inf2 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_WJets.root","READ");
  TFile *inf3 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_WJetsInv.root","READ");
  TFile *inf4 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_GJetsInv.root","READ");
  TFile *outf = new TFile("RazorScaleFactors_Inclusive_NEW.root","RECREATE");

  TH2F *ttbarNominal = (TH2F*)inf->Get("TTJetsScaleFactors");
  TH2F *ttbarUp = (TH2F*)inf->Get("TTJetsScaleFactorsUp");
  TH2F *ttbarDown = (TH2F*)inf->Get("TTJetsScaleFactorsDown");
  TH2F *wNominal = (TH2F*)inf2->Get("WJetsScaleFactors");
  TH2F *wUp = (TH2F*)inf2->Get("WJetsScaleFactorsUp");
  TH2F *wDown = (TH2F*)inf2->Get("WJetsScaleFactorsDown");
  TH2F *wInvNominal = (TH2F*)inf3->Get("WJetsInvScaleFactors");
  TH2F *wInvUp = (TH2F*)inf3->Get("WJetsInvScaleFactorsUp");
  TH2F *wInvDown = (TH2F*)inf3->Get("WJetsInvScaleFactorsDown");
  TH2F *GJetInvNominal = (TH2F*)inf4->Get("GJetsInvScaleFactors");

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
  outf->Close();
  inf->Close();

}

void doModify( TH2F *f, double correction, double error) {
  
  for (int i=0; i<f->GetXaxis()->GetNbins()+2; ++i) {
    for (int j=0; j<f->GetXaxis()->GetNbins()+2; ++j) {
      double val = correction*f->GetBinContent(i,j);
      double err = val * sqrt( pow(f->GetBinError(i,j)/f->GetBinContent(i,j),2) 
			       + pow( error/correction , 2));
      f->SetBinContent(i,j, val);
      f->SetBinError(i,j, err);
    }
  }
  
}

void addMultiJetCorrection() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root","READ");
  assert(inf);
  TFile *outf = new TFile("RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","RECREATE");

  TH2F *ttbarNominal = (TH2F*)inf->Get("TTJetsScaleFactors");
  TH2F *ttbarUp = (TH2F*)inf->Get("TTJetsScaleFactorsUp");
  TH2F *ttbarDown = (TH2F*)inf->Get("TTJetsScaleFactorsDown");
  TH2F *wNominal = (TH2F*)inf->Get("WJetsScaleFactors");
  TH2F *wUp = (TH2F*)inf->Get("WJetsScaleFactorsUp");
  TH2F *wDown = (TH2F*)inf->Get("WJetsScaleFactorsDown");
  TH2F *wInvNominal = (TH2F*)inf->Get("WJetsInvScaleFactors");
  TH2F *wInvUp = (TH2F*)inf->Get("WJetsInvScaleFactorsUp");
  TH2F *wInvDown = (TH2F*)inf->Get("WJetsInvScaleFactorsDown");
  TH2F *GJetsInvNominal = (TH2F*)inf->Get("GJetsInvScaleFactors");

  doModify(ttbarNominal, 0.9, 0.1);
  doModify(ttbarUp, 0.9, 0.1);
  doModify(ttbarDown, 0.9, 0.1);
  doModify(wNominal, 0.9, 0.1);
  doModify(wUp, 0.9, 0.1);
  doModify(wDown, 0.9, 0.1);
  doModify(wInvNominal, 0.9, 0.1);
  doModify(wInvUp, 0.9, 0.1);
  doModify(wInvDown, 0.9, 0.1);
  doModify(GJetsInvNominal, 0.87, 0.05);

  outf->WriteTObject(ttbarNominal);
  outf->WriteTObject(ttbarUp);
  outf->WriteTObject(ttbarDown);
  outf->WriteTObject(wNominal);
  outf->WriteTObject(wUp);
  outf->WriteTObject(wDown);
  outf->WriteTObject(wInvNominal);
  outf->WriteTObject(wInvUp);
  outf->WriteTObject(wInvDown);
  outf->WriteTObject(GJetsInvNominal);
  outf->Close();
  inf->Close();

}


void copyScaleFactorHistograms() {

  //doCopy();
  addMultiJetCorrection();

}
