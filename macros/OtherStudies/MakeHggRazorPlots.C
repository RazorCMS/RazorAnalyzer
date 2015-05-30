
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>
#include <map>
#include <iostream>

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kGray, kOrange+1, kGray+3, kBlack, kBlack};


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

double deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TRandom3 *random = new TRandom3(0);
  double intLumi = 4000; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processLabels;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back(dataLabel);
  } else {
    hasData = true;
    inputfiles.push_back("");
    processLabels.push_back("Hypothetical Data");    
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1F* histMRAllBkg =  new TH1F( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 100, 0, 3000);
  TH1F* histRsqAllBkg =  new TH1F( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", 100, 0, 1.5);
  TH1F* histMRAllBkg_AfterDPhiCut =  new TH1F("MRAllBkg_AfterDPhiCut", ";M_{R} [GeV/c^{2}];Number of Events", 100, 0, 3000);
  TH1F* histRsqAllBkg_AfterDPhiCut =  new TH1F( "RsqAllBkg_AfterDPhiCut", ";M_{R} [GeV/c^{2}];Number of Events", 100, 0, 1.5);
  histMRAllBkg->SetStats(false);
  histMRAllBkg_AfterDPhiCut->SetStats(false);
  histRsqAllBkg->SetStats(false);
  histRsqAllBkg_AfterDPhiCut->SetStats(false);
  

  vector<TH1F*> histNJets80;
  vector<TH1F*> histNJets60;
  vector<TH1F*> histNJets40;
  vector<TH1F*> histNJets30;
  vector<TH1F*> histPhoton1Pt;
  vector<TH1F*> histPhoton2Pt;
  vector<TH1F*> histJet1Pt;
  vector<TH1F*> histJet2Pt;
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histDPhiRazor;
  vector<TH1F*> histMgg;
  vector<TH1F*> histPtggPeakRegion;
  vector<TH1F*> histPtggSidebandRegion;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
    histMgg.push_back( new TH1F( Form("Mgg_%s",processLabels[i].c_str()), ";M_{#gamma#gamma} [GeV/c^{2}];Number of Events", 20, 100, 160));
    if (!hasData || i != 0) histMgg[i]->SetFillColor(color[i]);
    if (hasData && i==0) histMgg[i]->SetLineWidth(3);
    histMgg[i]->SetLineColor(color[i]);    
    histMgg[i]->SetStats(false);    
    histMgg[i]->Sumw2();

    histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 3000));
    if (!hasData || i != 0) histMR[i]->SetFillColor(color[i]);
    if (hasData && i==0) histMR[i]->SetLineWidth(3);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 1.5));
    if (!hasData || i != 0) histRsq[i]->SetFillColor(color[i]);
    if (hasData && i==0) histRsq[i]->SetLineWidth(3);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histNJets80.push_back( new TH1F( Form("NJets80_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5));
    if (!hasData || i != 0) histNJets80[i]->SetFillColor(color[i]);
    if (hasData && i==0) histNJets80[i]->SetLineWidth(3);
    histNJets80[i]->SetLineColor(color[i]);    
    histNJets80[i]->SetStats(false);    
    histNJets80[i]->Sumw2();

    histNJets60.push_back( new TH1F( Form("NJets60_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5));
    if (!hasData || i != 0) histNJets60[i]->SetFillColor(color[i]);
    if (hasData && i==0) histNJets60[i]->SetLineWidth(3);
    histNJets60[i]->SetLineColor(color[i]);    
    histNJets60[i]->SetStats(false);    
    histNJets60[i]->Sumw2();

    histNJets40.push_back( new TH1F( Form("NJets40_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5));
    if (!hasData || i != 0) histNJets40[i]->SetFillColor(color[i]);
    if (hasData && i==0) histNJets40[i]->SetLineWidth(3);
    histNJets40[i]->SetLineColor(color[i]);    
    histNJets40[i]->SetStats(false);    
    histNJets40[i]->Sumw2();

    histNJets30.push_back( new TH1F( Form("NJets30_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5));
    if (!hasData || i != 0) histNJets30[i]->SetFillColor(color[i]);
    if (hasData && i==0) histNJets30[i]->SetLineWidth(3);
    histNJets30[i]->SetLineColor(color[i]);    
    histNJets30[i]->SetStats(false);    
    histNJets30[i]->Sumw2();

    histPhoton1Pt.push_back( new TH1F( Form("Photon1Pt_%s",processLabels[i].c_str()), ";Leading Photon p_{T} [GeV/c];Number of Events", 50, 0, 200));
    if (!hasData || i != 0) histPhoton1Pt[i]->SetFillColor(color[i]);
    if (hasData && i==0) histPhoton1Pt[i]->SetLineWidth(3);
    histPhoton1Pt[i]->SetLineColor(color[i]);    
    histPhoton1Pt[i]->SetStats(false);    
    histPhoton1Pt[i]->Sumw2();

    histPhoton2Pt.push_back( new TH1F( Form("Photon2Pt_%s",processLabels[i].c_str()), ";Second Photon p_{T} [GeV/c];Number of Events", 50, 0, 200));
    if (!hasData || i != 0) histPhoton2Pt[i]->SetFillColor(color[i]);
    if (hasData && i==0) histPhoton2Pt[i]->SetLineWidth(3);
    histPhoton2Pt[i]->SetLineColor(color[i]);    
    histPhoton2Pt[i]->SetStats(false);    
    histPhoton2Pt[i]->Sumw2();

    histJet1Pt.push_back( new TH1F( Form("Jet1Pt_%s",processLabels[i].c_str()), ";Leading Jet p_{T} [GeV/c];Number of Events", 125, 0, 500));
    if (!hasData || i != 0) histJet1Pt[i]->SetFillColor(color[i]);
    if (hasData && i==0) histJet1Pt[i]->SetLineWidth(3);
    histJet1Pt[i]->SetLineColor(color[i]);    
    histJet1Pt[i]->SetStats(false);    
    histJet1Pt[i]->Sumw2();

    histJet2Pt.push_back( new TH1F( Form("Jet2Pt_%s",processLabels[i].c_str()), ";Second Jet p_{T} [GeV/c];Number of Events", 125, 0, 500));
    if (!hasData || i != 0) histJet2Pt[i]->SetFillColor(color[i]);
    if (hasData && i==0) histJet2Pt[i]->SetLineWidth(3);
    histJet2Pt[i]->SetLineColor(color[i]);    
    histJet2Pt[i]->SetStats(false);    
    histJet2Pt[i]->Sumw2();

    histPtggPeakRegion.push_back( new TH1F( Form("PtggPeakRegion_%s",processLabels[i].c_str()), ";p_{T #gamma#gamma} [GeV/c];Number of Events", 20, 0, 500));
    if (!hasData || i != 0) histPtggPeakRegion[i]->SetFillColor(color[i]);
    if (hasData && i==0) histPtggPeakRegion[i]->SetLineWidth(3);
    histPtggPeakRegion[i]->SetLineColor(color[i]);    
    histPtggPeakRegion[i]->SetStats(false);    
    histPtggPeakRegion[i]->Sumw2();

    histPtggSidebandRegion.push_back( new TH1F( Form("PtggSidebandRegion_%s",processLabels[i].c_str()), ";p_{T #gamma#gamma} [GeV/c];Number of Events", 20, 0, 500));
    if (!hasData || i != 0) histPtggSidebandRegion[i]->SetFillColor(color[i]);
    if (hasData && i==0) histPtggSidebandRegion[i]->SetLineWidth(3);
    histPtggSidebandRegion[i]->SetLineColor(color[i]);    
    histPtggSidebandRegion[i]->SetStats(false);    
    histPtggSidebandRegion[i]->Sumw2();

    histDPhiRazor.push_back( new TH1F( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    if (!hasData || i != 0) histDPhiRazor[i]->SetFillColor(color[i]);
    if (hasData && i==0) histDPhiRazor[i]->SetLineWidth(3);
    histDPhiRazor[i]->SetLineColor(color[i]);
    histDPhiRazor[i]->SetStats(false); 
 }

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************



  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HighRes");
    assert(tree);

    float w = 0;
    int njets = 0;
    float jetpt[40]; //[njets]
    float jeteta[40]; //[njets]
    float jetphi[40]; //[njets]
    float pho1pt = 0;
    float pho1eta = 0;
    float pho1phi = 0;
    float pho2pt = 0;
    float pho2eta = 0;
    float pho2phi = 0;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float mgg = 0;
    float ptgg = 0;
    float Rsq = 0;
    uint run = 0;
    uint event = 0;

    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("weight",&w);
    tree->SetBranchAddress("n_Jets",&njets);
    tree->SetBranchAddress("jet_Pt",&jetpt);
    tree->SetBranchAddress("jet_Eta",&jeteta);
    tree->SetBranchAddress("jet_Phi",&jetphi);
    tree->SetBranchAddress("pho1Pt",&pho1pt);
    tree->SetBranchAddress("Pho1Eta",&pho1eta);
    tree->SetBranchAddress("pho1Phi",&pho1phi);
    tree->SetBranchAddress("pho2Pt",&pho2pt);
    tree->SetBranchAddress("Pho2Eta",&pho2eta);
    tree->SetBranchAddress("pho2Phi",&pho2phi);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("mGammaGamma",&mgg);
    tree->SetBranchAddress("pTGammaGamma",&ptgg);
    tree->SetBranchAddress("t1Rsq",&Rsq);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      double weight = 1;

      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       

      //Options
      if (option == 0 ) {
      }

      //apply selection cuts
      if (!(pho1pt > 25 && pho2pt > 25)) continue;
      if (!(pho1pt > 40 || pho2pt > 40)) continue;
      if (!( fabs(pho1eta) < 1.44 && fabs(pho2eta) < 1.44 ) ) continue;
      if (!(MR > 350 && Rsq > 0.035)) continue;

      //count jets
      int njets80 = 0;
      int njets60 = 0;
      int njets40 = 0;
      double leadjetpt = 0;
      double secondjetpt = 0;
      for(int j=0; j < njets; ++j) {

	if (fabs(jeteta[j]) >=  3.0) continue;

	//cout << "jet " << j << " : " << jetpt[j] << "\n";
	if ( deltaR( jeteta[j], pho1eta, jetphi[j], pho1phi) < 0.4) {
	  //cout << "overlap : " << pho1pt << " " << jetpt[j] << " " << deltaR( jeteta[j], pho1eta, jetphi[j], pho1phi) << "\n";
	  continue; 
	}
	if ( deltaR( jeteta[j], pho2eta, jetphi[j], pho2phi) < 0.4) {
	  //cout << "overlap2 : " << pho2pt << " " << jetpt[j] << " " << deltaR( jeteta[j], pho2eta, jetphi[j], pho2phi) << "\n";
	  continue; 
	}

	if (jetpt[j] > leadjetpt) {
	  secondjetpt = leadjetpt;
	  leadjetpt = jetpt[j];
	} else if (jetpt[j] > secondjetpt) {
	  secondjetpt = jetpt[j];
	}

      	if (jetpt[j] > 80 && fabs(jeteta[j]) < 3.0) {
      	  njets80++;
      	}
      	if (jetpt[j] > 60 && fabs(jeteta[j]) < 3.0) {
      	  njets60++;
      	}
      	if (jetpt[j] > 40 && fabs(jeteta[j]) < 3.0) {
      	  njets40++;
      	}
      }
      //cout << "lead jets: " << leadjetpt << " " << secondjetpt << "\n";
      cout << run << " " << event << " : " << mgg << " " << ptgg << " " << MR << " " << Rsq << " : " << pho1eta << " " << pho2eta << " \n";

      histMgg[i]->Fill(mgg, weight);
      if (mgg > 125 && mgg < 135) {
	histRsq[i]->Fill(Rsq, weight);	  
	histMR[i]->Fill(MR, weight);
	histPtggPeakRegion[i]->Fill(ptgg, weight);
	histNJets80[i]->Fill(njets80, weight);
	histNJets60[i]->Fill(njets60, weight);
	histNJets40[i]->Fill(njets40, weight);
	histNJets30[i]->Fill(njets, weight);
	histJet1Pt[i]->Fill(fmin(fmax(0.1,leadjetpt),499), weight);
	histJet2Pt[i]->Fill(fmin(fmax(0.1,secondjetpt),499), weight);
	if (pho1pt > pho2pt) {
	  histPhoton1Pt[i]->Fill(pho1pt, weight);
	  histPhoton2Pt[i]->Fill(pho2pt, weight);
	} else {
	  histPhoton1Pt[i]->Fill(pho2pt, weight);
	  histPhoton2Pt[i]->Fill(pho1pt, weight);
	}
      } else {
	histPtggSidebandRegion[i]->Fill(ptgg, weight);
      }
  
    }

    inputFile->Close();
    delete inputFile;
  
  }
  

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //Mgg
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMgg = new THStack();

  if (hasData) {
    for (Int_t i = histMgg.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMgg[i]->GetNbinsX()+1; j++) {
  	intError += histMgg[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMgg[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMgg[i]->Integral() > 0) {
  	stackMgg->Add(histMgg[i]);
      }
    }    
  } else {
    for (Int_t i = histMgg.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMgg[i]->GetNbinsX()+1; j++) {
  	intError += histMgg[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMgg[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMgg[i]->Integral() > 0) {
  	stackMgg->Add(histMgg[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histMgg.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histMgg[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMgg[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackMgg->Draw("hist");
  // stackMgg->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMgg->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackMgg->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMgg->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackMgg->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histMgg[0]->SetLineColor(kBlack);
    histMgg[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histMgg[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s", "19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Mgg%s.gif",Label.c_str()));

 

 //*******************************************************************************************
  //NJets80
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNJets80 = new THStack();

  if (hasData) {
    for (Int_t i = histNJets80.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets80[i]->GetNbinsX()+1; j++) {
  	intError += histNJets80[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets80[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets80[i]->Integral() > 0) {
  	stackNJets80->Add(histNJets80[i]);
      }
    }    
  } else {
    for (Int_t i = histNJets80.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets80[i]->GetNbinsX()+1; j++) {
  	intError += histNJets80[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets80[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets80[i]->Integral() > 0) {
  	stackNJets80->Add(histNJets80[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histNJets80.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNJets80[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histNJets80[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackNJets80->Draw("hist");
  // stackNJets80->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNJets80->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackNJets80->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNJets80->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackNJets80->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histNJets80[0]->SetLineColor(kBlack);
    histNJets80[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histNJets80[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("NJets80%s.gif",Label.c_str()));

  //*******************************************************************************************
  //NJets60
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNJets60 = new THStack();

  if (hasData) {
    for (Int_t i = histNJets60.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets60[i]->GetNbinsX()+1; j++) {
  	intError += histNJets60[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets60[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets60[i]->Integral() > 0) {
  	stackNJets60->Add(histNJets60[i]);
      }
    }    
  } else {
    for (Int_t i = histNJets60.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets60[i]->GetNbinsX()+1; j++) {
  	intError += histNJets60[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets60[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets60[i]->Integral() > 0) {
  	stackNJets60->Add(histNJets60[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histNJets60.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNJets60[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histNJets60[i],processLabels[i].c_str(), "F");
    }
  }
  


  if (hasData) {
    histNJets60[0]->SetLineColor(kBlack);
    histNJets60[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histNJets60[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("NJets60%s.gif",Label.c_str()));

  //*******************************************************************************************
  //NJets40
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNJets40 = new THStack();

  if (hasData) {
    for (Int_t i = histNJets40.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets40[i]->GetNbinsX()+1; j++) {
  	intError += histNJets40[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets40[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets40[i]->Integral() > 0) {
  	stackNJets40->Add(histNJets40[i]);
      }
    }    
  } else {
    for (Int_t i = histNJets40.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets40[i]->GetNbinsX()+1; j++) {
  	intError += histNJets40[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets40[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets40[i]->Integral() > 0) {
  	stackNJets40->Add(histNJets40[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histNJets40.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNJets40[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histNJets40[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackNJets40->Draw("hist");
  // stackNJets40->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNJets40->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackNJets40->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNJets40->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackNJets40->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histNJets40[0]->SetLineColor(kBlack);
    histNJets40[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histNJets40[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("NJets40%s.gif",Label.c_str()));

  //*******************************************************************************************
  //NJets30
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNJets30 = new THStack();

  if (hasData) {
    for (Int_t i = histNJets30.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets30[i]->GetNbinsX()+1; j++) {
  	intError += histNJets30[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets30[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets30[i]->Integral() > 0) {
  	stackNJets30->Add(histNJets30[i]);
      }
    }    
  } else {
    for (Int_t i = histNJets30.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histNJets30[i]->GetNbinsX()+1; j++) {
  	intError += histNJets30[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histNJets30[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histNJets30[i]->Integral() > 0) {
  	stackNJets30->Add(histNJets30[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histNJets30.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNJets30[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histNJets30[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackNJets30->Draw("hist");
  // stackNJets30->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNJets30->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackNJets30->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNJets30->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackNJets30->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histNJets30[0]->SetLineColor(kBlack);
    histNJets30[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histNJets30[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("NJets30%s.gif",Label.c_str()));

 

  //*******************************************************************************************
  //Photon1Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackPhoton1Pt = new THStack();

  if (hasData) {
    for (Int_t i = histPhoton1Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPhoton1Pt[i]->GetNbinsX()+1; j++) {
  	intError += histPhoton1Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPhoton1Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPhoton1Pt[i]->Integral() > 0) {
  	stackPhoton1Pt->Add(histPhoton1Pt[i]);
      }
    }    
  } else {
    for (Int_t i = histPhoton1Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPhoton1Pt[i]->GetNbinsX()+1; j++) {
  	intError += histPhoton1Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPhoton1Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPhoton1Pt[i]->Integral() > 0) {
  	stackPhoton1Pt->Add(histPhoton1Pt[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histPhoton1Pt.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histPhoton1Pt[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histPhoton1Pt[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackPhoton1Pt->Draw("hist");
  // stackPhoton1Pt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackPhoton1Pt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackPhoton1Pt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackPhoton1Pt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackPhoton1Pt->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histPhoton1Pt[0]->SetLineColor(kBlack);
    histPhoton1Pt[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histPhoton1Pt[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Photon1Pt%s.gif",Label.c_str()));

 
  //*******************************************************************************************
  //Photon2Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackPhoton2Pt = new THStack();

  if (hasData) {
    for (Int_t i = histPhoton2Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPhoton2Pt[i]->GetNbinsX()+1; j++) {
  	intError += histPhoton2Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPhoton2Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPhoton2Pt[i]->Integral() > 0) {
  	stackPhoton2Pt->Add(histPhoton2Pt[i]);
      }
    }    
  } else {
    for (Int_t i = histPhoton2Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPhoton2Pt[i]->GetNbinsX()+1; j++) {
  	intError += histPhoton2Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPhoton2Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPhoton2Pt[i]->Integral() > 0) {
  	stackPhoton2Pt->Add(histPhoton2Pt[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histPhoton2Pt.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histPhoton2Pt[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histPhoton2Pt[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackPhoton2Pt->Draw("hist");
  // stackPhoton2Pt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackPhoton2Pt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackPhoton2Pt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackPhoton2Pt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackPhoton2Pt->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histPhoton2Pt[0]->SetLineColor(kBlack);
    histPhoton2Pt[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histPhoton2Pt[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Photon2Pt%s.gif",Label.c_str()));

 

  //*******************************************************************************************
  //Jet1Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackJet1Pt = new THStack();

  if (hasData) {
    for (Int_t i = histJet1Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histJet1Pt[i]->GetNbinsX()+1; j++) {
  	intError += histJet1Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histJet1Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histJet1Pt[i]->Integral() > 0) {
  	stackJet1Pt->Add(histJet1Pt[i]);
      }
    }    
  } else {
    for (Int_t i = histJet1Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histJet1Pt[i]->GetNbinsX()+1; j++) {
  	intError += histJet1Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histJet1Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histJet1Pt[i]->Integral() > 0) {
  	stackJet1Pt->Add(histJet1Pt[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histJet1Pt.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histJet1Pt[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histJet1Pt[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackJet1Pt->Draw("hist");
  // stackJet1Pt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackJet1Pt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackJet1Pt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackJet1Pt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackJet1Pt->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histJet1Pt[0]->SetLineColor(kBlack);
    histJet1Pt[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histJet1Pt[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Jet1Pt%s.gif",Label.c_str()));

 
  //*******************************************************************************************
  //Jet2Pt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackJet2Pt = new THStack();

  if (hasData) {
    for (Int_t i = histJet2Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histJet2Pt[i]->GetNbinsX()+1; j++) {
  	intError += histJet2Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histJet2Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histJet2Pt[i]->Integral() > 0) {
  	stackJet2Pt->Add(histJet2Pt[i]);
      }
    }    
  } else {
    for (Int_t i = histJet2Pt.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histJet2Pt[i]->GetNbinsX()+1; j++) {
  	intError += histJet2Pt[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histJet2Pt[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histJet2Pt[i]->Integral() > 0) {
  	stackJet2Pt->Add(histJet2Pt[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histJet2Pt.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histJet2Pt[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histJet2Pt[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackJet2Pt->Draw("hist");
  // stackJet2Pt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackJet2Pt->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackJet2Pt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackJet2Pt->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackJet2Pt->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histJet2Pt[0]->SetLineColor(kBlack);
    histJet2Pt[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histJet2Pt[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Jet2Pt%s.gif",Label.c_str()));

 


  //*******************************************************************************************
  //PtggPeakRegion
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackPtggPeakRegion = new THStack();

  if (hasData) {
    for (Int_t i = histPtggPeakRegion.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPtggPeakRegion[i]->GetNbinsX()+1; j++) {
  	intError += histPtggPeakRegion[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPtggPeakRegion[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPtggPeakRegion[i]->Integral() > 0) {
  	stackPtggPeakRegion->Add(histPtggPeakRegion[i]);
      }
    }    
  } else {
    for (Int_t i = histPtggPeakRegion.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPtggPeakRegion[i]->GetNbinsX()+1; j++) {
  	intError += histPtggPeakRegion[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPtggPeakRegion[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPtggPeakRegion[i]->Integral() > 0) {
  	stackPtggPeakRegion->Add(histPtggPeakRegion[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histPtggPeakRegion.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histPtggPeakRegion[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histPtggPeakRegion[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackPtggPeakRegion->Draw("hist");
  // stackPtggPeakRegion->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackPtggPeakRegion->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackPtggPeakRegion->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackPtggPeakRegion->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackPtggPeakRegion->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histPtggPeakRegion[0]->SetLineColor(kBlack);
    histPtggPeakRegion[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histPtggPeakRegion[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("PtggPeakRegion%s.gif",Label.c_str()));

   //*******************************************************************************************
  //PtggSidebandRegion
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.70,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackPtggSidebandRegion = new THStack();

  if (hasData) {
    for (Int_t i = histPtggSidebandRegion.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPtggSidebandRegion[i]->GetNbinsX()+1; j++) {
  	intError += histPtggSidebandRegion[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPtggSidebandRegion[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPtggSidebandRegion[i]->Integral() > 0) {
  	stackPtggSidebandRegion->Add(histPtggSidebandRegion[i]);
      }
    }    
  } else {
    for (Int_t i = histPtggSidebandRegion.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histPtggSidebandRegion[i]->GetNbinsX()+1; j++) {
  	intError += histPtggSidebandRegion[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histPtggSidebandRegion[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histPtggSidebandRegion[i]->Integral() > 0) {
  	stackPtggSidebandRegion->Add(histPtggSidebandRegion[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histPtggSidebandRegion.size()); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histPtggSidebandRegion[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histPtggSidebandRegion[i],processLabels[i].c_str(), "F");
    }
  }
  

  //stackPtggSidebandRegion->Draw("hist");
  // stackPtggSidebandRegion->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackPtggSidebandRegion->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackPtggSidebandRegion->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackPtggSidebandRegion->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackPtggSidebandRegion->GetHistogram()->SetLineColor(kBlack);

  if (hasData) {
    histPtggSidebandRegion[0]->SetLineColor(kBlack);
    histPtggSidebandRegion[0]->Draw("e1same");
    cout << processLabels[0] << " : " << histPtggSidebandRegion[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("PtggSidebandRegion%s.gif",Label.c_str()));

 

  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
 
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMgg[i], Form("histMgg_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets80[i], Form("histNJets80_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets60[i], Form("histNJets60_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets30[i], Form("histNJets30_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histPhoton1Pt[i], Form("histPhoton1Pt_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histPhoton2Pt[i], Form("histPhoton2Pt_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histJet1Pt[i], Form("histJet1Pt_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histJet2Pt[i], Form("histJet2Pt_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histPtggPeakRegion[i], Form("histPtggPeakRegion_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histPtggSidebandRegion[i], Form("histPtggSidebandRegion_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
  // file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


 void MakeHggRazorPlots() {

   string datafile = "";
   string dataLabel = "";

   // string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1qqqq_2J_mGl-1400_mLSP-100_1pb_weighted.root";  
   // string dataLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   // string dataLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
   //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   //string dataLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   

   dataLabel = "data";
   //Main Analysis Categories
   //datafile = "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run1/HggRazor/HggRazorNtuple/DoublePhotonDataRunABCD_Trigger_GoodLumi.root";

   //HighRes - LowRes categories only
   datafile = "/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorHgg/HggRazorSimpleCat_Filtered.root";


   vector<string> bkgfiles;
   vector<string> bkgLabels;
   
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_TTJets_1pb_weighted.root");  
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_DYJetsToLL_1pb_weighted.root");
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_WJetsToLNu_1pb_weighted.root");
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_ZJetsToNuNu_1pb_weighted.root");
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_QCD_1pb_weighted.root"); 
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_SingleTop_1pb_weighted.root"); 
   // bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVM/HbbRazor_Multiboson_1pb_weighted.root"); 
   
   // bkgLabels.push_back("TTJets");
   // bkgLabels.push_back("DYJetsToLL");
   // bkgLabels.push_back("WJetsToLNu");
   // bkgLabels.push_back("ZJetsToNuNu");
   // bkgLabels.push_back("QCD");
   // bkgLabels.push_back("SingleTop");
   // bkgLabels.push_back("Other");

   RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,0,1,"MR350Rsq0p035","H#rightarrow#gamma#gamma Razor");
   
 
 }
 
