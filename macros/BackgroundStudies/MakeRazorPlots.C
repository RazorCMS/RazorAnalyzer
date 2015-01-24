
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <map>
#include <iostream>

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kBlack, kOrange+1, kGray, kBlack, kBlack};


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


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( string signalfile, string signalLabel,  vector<string> bkgfiles,vector<string> bkgLabels, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  double intLumi = 4000; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processLabels;

  bool hasSignal = false;
  if (signalfile != "") {
    hasSignal = true;
    inputfiles.push_back(signalfile);
    processLabels.push_back(signalLabel);
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
  

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq; 
  vector<TH1F*> histDPhiRazor;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
    histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 3000));
    if (!hasSignal || i != 0) histMR[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histMR[i]->SetLineWidth(3);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 1.5));
    if (!hasSignal || i != 0) histRsq[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histRsq[i]->SetLineWidth(3);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histDPhiRazor.push_back( new TH1F( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    if (!hasSignal || i != 0) histDPhiRazor[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histDPhiRazor[i]->SetLineWidth(3);
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
    tree = (TTree*)inputFile->Get("RazorInclusive");
  // if (box == 0) {
    //   tree = (TTree*)inputFile->Get("MultiJet");
    // } else if (box == 1) {
    //   tree = (TTree*)inputFile->Get("LooseLeptonMultiJet");
    // } else if (box == 2) {
    //   tree = (TTree*)inputFile->Get("MuMultiJet");
    // } else if (box == 3) {
    //   tree = (TTree*)inputFile->Get("EleMultiJet");
    // }
 
    float weight = 0;
    int box = -1;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);


    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       


      if (intLumi*weight > 100) continue;

      //Box Options
      if (option == 0 ) {
	if (!(nBTaggedJets == 0)) continue;
      }
      if (option == 1 ) {
	if (!(nBTaggedJets >= 1)) continue;
      }

      if (boxOption == 0) {
	if (box != 8) continue;
      } else if (boxOption == 1) {
	if (box != 7) continue;
      } else if (boxOption == 2) {
	if (box != 3) continue;
      } else if (boxOption == 3) {
	if (box != 5) continue;
      } else if (boxOption == 11) {
	if (!(box == 7 || box == 3 || box == 5)) continue;
      }

      //apply baseline cuts
      if (!(MR > 300 && Rsq > 0.1)) continue;

      if (!hasSignal || i>0) {
	histMRAllBkg->Fill(MR, intLumi*weight);
	histRsqAllBkg->Fill(Rsq, intLumi*weight);
	if (fabs(dPhiRazor) < 2.7) {
	  histMRAllBkg_AfterDPhiCut->Fill(MR, intLumi*weight);
	  histRsqAllBkg_AfterDPhiCut->Fill(Rsq, intLumi*weight);
	}
      }

      if (
	  //Rsq > 0.1 && MR > 500
	  Rsq > 0.25 && MR > 1000
	  ) {
	histDPhiRazor[i]->Fill(dPhiRazor, intLumi*weight);	
	
	if (fabs(dPhiRazor) < 2.7) {
	  histMR[i]->Fill(MR, intLumi*weight);
	  histRsq[i]->Fill(Rsq, intLumi*weight);	  
	}
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
  //MR
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR = new THStack();

  if (hasSignal) {
    for (Int_t i = histMR.size()-1; i >= 1; i--) {
      double intError = 0;
      for(int j=1; j < histMR[i]->GetNbinsX()+1; j++) {
	intError += histMR[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMR[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMR[i]->Integral() > 0) {
	stackMR->Add(histMR[i]);
      }
    }    
  } else {
    for (Int_t i = histMR.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMR[i]->GetNbinsX()+1; j++) {
	intError += histMR[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMR[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMR[i]->Integral() > 0) {
	stackMR->Add(histMR[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histMR.size()); ++i) {
    if (hasSignal && i==0) {
      legend->AddEntry(histMR[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMR[i],processLabels[i].c_str(), "F");
    }
  }
  

  stackMR->Draw("hist");
  stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());

  if (hasSignal) {
    histMR[0]->Draw("histsame");
    cout << processLabels[0] << " : " << histMR[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 13 TeV, #int L = %d fb^{-1}, %s",int(intLumi/1000), latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("MR%s.gif",Label.c_str()));

 

  // //*******************************************************************************************
  // //Rsq
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // legend = new TLegend(0.60,0.54,0.90,0.84);
  // legend->SetTextSize(0.03);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);

  // THStack *stackRsq = new THStack();

  // if (hasSignal) {
  //   for (Int_t i = histRsq.size()-1; i >= 1; i--) {
  //     cout << processLabels[i] << " : " << histRsq[i]->GetSumOfWeights() << "\n";
  //     if ( histRsq[i]->Integral() > 0) {
  // 	stackRsq->Add(histRsq[i]);
  //     }
  //   }    
  // } else {
  //   for (Int_t i = histRsq.size()-1; i >= 0; i--) {
  //     cout << processLabels[i] << " : " << histRsq[i]->GetSumOfWeights() << "\n";
  //     if ( histRsq[i]->Integral() > 0) {
  // 	stackRsq->Add(histRsq[i]);
  //     }
  //   }
  // }
  // for (Int_t i = 0 ; i < int(histRsq.size()); ++i) {
  //   if (hasSignal && i==0) {
  //     legend->AddEntry(histRsq[i],processLabels[i].c_str(), "L");
  //   } else {
  //     legend->AddEntry(histRsq[i],processLabels[i].c_str(), "F");
  //   }
  // }
  

  // stackRsq->Draw();
  // stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());

  // if (hasSignal) {
  //   histRsq[0]->Draw("sameL");
  // }

  // legend->Draw();
  // cv->SaveAs(Form("Rsq%s.gif",Label.c_str()));

 

  // //*******************************************************************************************
  // //DPhiRazor
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // legend = new TLegend(0.15,0.54,0.65,0.84);
  // legend->SetTextSize(0.03);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);

  // THStack *stackDPhiRazor = new THStack();

  // if (hasSignal) {
  //   for (Int_t i = histDPhiRazor.size()-1; i >= 1; i--) {
  //     cout << processLabels[i] << " : " << histDPhiRazor[i]->GetSumOfWeights() << "\n";
  //     if ( histDPhiRazor[i]->Integral() > 0) {
  // 	stackDPhiRazor->Add(histDPhiRazor[i]);
  //     }
  //   }    
  // } else {
  //   for (Int_t i = histDPhiRazor.size()-1; i >= 0; i--) {
  //     cout << processLabels[i] << " : " << histDPhiRazor[i]->GetSumOfWeights() << "\n";
  //     if ( histDPhiRazor[i]->Integral() > 0) {
  // 	stackDPhiRazor->Add(histDPhiRazor[i]);
  //     }
  //   }
  // }
  // for (Int_t i = 0 ; i < int(histDPhiRazor.size()); ++i) {
  //   if (hasSignal && i==0) {
  //     legend->AddEntry(histDPhiRazor[i],processLabels[i].c_str(), "L");
  //   } else {
  //     legend->AddEntry(histDPhiRazor[i],processLabels[i].c_str(), "F");
  //   }
  // }
  

  // stackDPhiRazor->Draw();
  // stackDPhiRazor->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDPhiRazor->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackDPhiRazor->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDPhiRazor->GetHists()->At(0)))->GetYaxis()->GetTitle());
  // stackDPhiRazor->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);

  // if (hasSignal) {
  //   histDPhiRazor[0]->Draw("sameL");
  // }

  // legend->Draw();
  // cv->SaveAs(Form("DPhiRazor%s.gif",Label.c_str()));

 
  // //*******************************************************************************************
  // //MR Before and After DPhi Cut
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // legend = new TLegend(0.50,0.54,0.90,0.84);
  // legend->SetTextSize(0.03);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);

  // histMRAllBkg = NormalizeHist(histMRAllBkg);
  // histMRAllBkg_AfterDPhiCut = NormalizeHist(histMRAllBkg_AfterDPhiCut);

  // histMRAllBkg->SetLineColor(kBlue);
  // histMRAllBkg_AfterDPhiCut->SetLineColor(kRed);
  // histMRAllBkg->GetYaxis()->SetTitle("Fraction of Events");
  // histMRAllBkg->GetYaxis()->SetTitleOffset(1.2);
  // histMRAllBkg_AfterDPhiCut->GetYaxis()->SetTitle("Fraction of Events");

  // legend->AddEntry(histMRAllBkg, "No #Delta#phi_{Razor} cut", "L");
  // legend->AddEntry(histMRAllBkg_AfterDPhiCut, "#Delta#phi_{Razor} < 2.7 cut", "L");

  // histMRAllBkg->Draw("hist");
  // histMRAllBkg_AfterDPhiCut->Draw("histsame");

  // legend->Draw();
  // cv->SetLogy();
  // cv->SaveAs(Form("MRBeforeAfterDPhiCut.gif",Label.c_str()));

  // //*******************************************************************************************
  // //Rsq Before and After DPhi Cut
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // legend = new TLegend(0.50,0.54,0.90,0.84);
  // legend->SetTextSize(0.03);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);

  // histRsqAllBkg = NormalizeHist(histRsqAllBkg);
  // histRsqAllBkg_AfterDPhiCut = NormalizeHist(histRsqAllBkg_AfterDPhiCut);

  // histRsqAllBkg->SetLineColor(kBlue);
  // histRsqAllBkg_AfterDPhiCut->SetLineColor(kRed);
  // histRsqAllBkg->GetYaxis()->SetTitle("Fraction of Events");
  // histRsqAllBkg->GetYaxis()->SetTitleOffset(1.2);
  // histRsqAllBkg_AfterDPhiCut->GetYaxis()->SetTitle("Fraction of Events");

  // legend->AddEntry(histRsqAllBkg, "No #Delta#phi_{Razor} cut", "L");
  // legend->AddEntry(histRsqAllBkg_AfterDPhiCut, "#Delta#phi_{Razor} < 2.7 cut", "L");

  // histRsqAllBkg->Draw("hist");
  // histRsqAllBkg_AfterDPhiCut->Draw("histsame");

  // legend->Draw();
  // cv->SetLogy();
  // cv->SaveAs(Form("RsqBeforeAfterDPhiCut.gif",Label.c_str()));


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
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
  file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


 void MakeRazorPlots() {

    // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/newWithElectronD0Cut/Inclusive/RazorAnalysis_SMS-T1qqqq_2J_mGl-1400_mLSP-100_25ns_weighted.root";  
    // string signalLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/newWithElectronD0Cut/Inclusive/RazorAnalysis_SMS-T1bbbb_2J_mGl-1500_mLSP-100_25ns_weighted.root";  
   // string signalLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
     string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/newWithElectronD0Cut/Inclusive/RazorAnalysis_SMS-T1tttt_2J_mGl-1500_mLSP-100_25ns_weighted.root";  
      string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   
   vector<string> bkgfiles;
   vector<string> bkgLabels;
   
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_TTJets_25ns_weighted.root");  
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_DYJetsToLL_HT100ToInf_25ns_weighted.root");
    bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_WJetsToLNu_HT100ToInf_25ns_weighted.root");
    bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_ZJetsToNuNu_HT100ToInf_25ns_weighted.root");
    bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_QCDHT100ToInf_25ns_weighted.root"); 
    bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_SingleTop_25ns_weighted.root"); 
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_Multiboson_25ns_weighted.root"); 
    
   bkgLabels.push_back("TTJets");
   bkgLabels.push_back("DYJetsToLL");
   bkgLabels.push_back("WJetsToLNu");
   bkgLabels.push_back("ZJetsToNuNu");
   bkgLabels.push_back("QCD");
   bkgLabels.push_back("SingleTop");
   bkgLabels.push_back("Other");
   
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,0,"T1qqqq_MultiJet_ZeroBTags", "MultiJet Box 0 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"T1bbbb_MultiJet_OneOrMoreBTags","MultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"T1tttt_MultiJet_OneOrMoreBTags","MultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,1,1,"T1tttt_LooseLeptonMultiJet_OneOrMoreBTags","LooseLeptonMultiJet Box #geq 1 b-tag" );
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,2,1,"T1tttt_MuMultiJet_OneOrMoreBTags","MuMultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,3,1,"T1tttt_EleMultiJet_OneOrMoreBTags","EleMultiJet Box #geq 1 b-tag");
   RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,11,1,"T1tttt_LeptonMultiJet_OneOrMoreBTags","Lepton+MultiJet Box #geq 1 b-tag");
 
   //RunMakeRazorPlots("","",bkgfiles,bkgLabels,8,0,"QCD_MultiJet_ZeroBTag");
   //RunMakeRazorPlots("","",bkgfiles,bkgLabels,8,1,"QCD_MultiJet_OneOrMoreBTag");
 
 }
 
