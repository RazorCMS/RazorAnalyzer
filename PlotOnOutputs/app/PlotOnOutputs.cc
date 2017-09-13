#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>

//#include <tree.hh>

const bool _debug = true;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./PlotOnOytputs inputList outputFile\n" << std::endl;
  srand(time(NULL));
  gROOT->Reset();
  gStyle->SetOptStat(0);
  if ( argc == 3 )
  { 
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        } else {
                while(getline(file,line)) ++n;
                std::cout << "n = " << n << "\n" << std::endl;
        }
        //std::string outputFile = argv[2];
  }

  std::ifstream ifs( argv[1], std::ifstream::in );
  assert(ifs);
  //while(std::getline(ifs,line)) n++;

  TFile* fin[n]; 
  TTree* tree[n]; 

  int i = 0;
  std::string process, rootFileName;
  std::string labelProcess[n];
  while ( ifs.good() ){
          ifs >> process >> rootFileName;
          if ( ifs.eof() ) continue;
          if ( process.find("#") != std::string::npos ) continue;
          if ( _debug ) std::cout << process << " " << rootFileName << std::endl;
          labelProcess[i] = process;
          fin[i] = new TFile( rootFileName.c_str(), "READ");
          //assert( fin[i] );
          if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;

          //------------------------
          //Getting TTree and Histos
          //------------------------
          tree[i] = (TTree*)fin[i]->Get("HggRazorLeptons");

          i++;
          //std::cout << "i = " << i << "\n" << std::endl;
  }

  
  //********************************************************
  //Print output
  //********************************************************
  std::string outputFile = argv[2];
  TFile* fout = new TFile( outputFile.c_str(), "RECREATE");

  int m = 2*n;
  std::cout << "m = " << m << "\n" << std::endl;

  TH1F* hist_MR[m];
  TH1F* hist_Rsq[m]; 
  TH1F* hist_HPt[m]; //Higgs Pt
  TH1F* hist_pTGammaGamma[m]; //reconstructed Higgs Pt
  TH1F* hist_HT[m];  //HT = sum of photon pT  + jet pT
  TH1F* hist_MET[m];

  //MR
  TCanvas* c_MR = new TCanvas( "c_MR", "c_MR", 800, 700 );
  TLegend* leg_MR = new TLegend(0.7,0.7,0.9,0.9);
  c_MR->SetHighLightColor(2);
  c_MR->SetFillColor(0);
  c_MR->SetBorderMode(0);
  c_MR->SetBorderSize(2);
  c_MR->SetLeftMargin( leftMargin );
  c_MR->SetRightMargin( 1.6*rightMargin );
  c_MR->SetTopMargin( topMargin );
  c_MR->SetBottomMargin( bottomMargin );
  c_MR->SetFrameBorderMode(0);
  c_MR->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",100,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i]->Draw(draw_MR_str.c_str(),"box > -1 && box <5","HSAME");
          hist_MR[i]->SetLineColor(i+1);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_MR->AddEntry(hist_MR[i],add_entry_str.c_str(),"l");
          std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",100,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i-n]->Draw(draw_MR_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_MR[i]->SetLineColor(i+1);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_MR->AddEntry(hist_MR[i],add_entry_str.c_str(),"l");
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MR[0]->SetTitle(";MR;#Events");
  hist_MR[0]->GetYaxis()->SetRangeUser(0,0.09);
  leg_MR->SetBorderSize(0);
  leg_MR->Draw("SAME");
  leg_MR->Write();
  c_MR->Write();
  c_MR->Update();
  c_MR->SaveAs("c_MR.pdf");

  //Rsq
  TCanvas* c_Rsq = new TCanvas( "c_Rsq", "c_Rsq", 800, 700 );
  TLegend* leg_Rsq = new TLegend(0.7,0.7,0.9,0.9);
  c_Rsq->SetHighLightColor(2);
  c_Rsq->SetFillColor(0);
  c_Rsq->SetBorderMode(0);
  c_Rsq->SetBorderSize(2);
  c_Rsq->SetLeftMargin( leftMargin );
  c_Rsq->SetRightMargin( 1.6*rightMargin );
  c_Rsq->SetTopMargin( topMargin );
  c_Rsq->SetBottomMargin( bottomMargin );
  c_Rsq->SetFrameBorderMode(0);
  c_Rsq->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",100,0,1.3);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i]->Draw(draw_Rsq_str.c_str(),"box > -1 && box < 5","HSAME");
          hist_Rsq[i]->SetLineColor(i+1);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_Rsq->AddEntry(hist_Rsq[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",100,0,1.3);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i-n]->Draw(draw_Rsq_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_Rsq[i]->SetLineColor(i+1);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_Rsq->AddEntry(hist_Rsq[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_Rsq[0]->SetTitle(";Rsq;#Events");
  hist_Rsq[0]->GetYaxis()->SetRangeUser(0,0.08);
  leg_Rsq->SetBorderSize(0);
  leg_Rsq->Draw("SAME");
  leg_Rsq->Write();
  c_Rsq->Write();
  c_Rsq->Update();
  c_Rsq->SaveAs("c_Rsq.pdf");

  //HPt
  TCanvas* c_HPt = new TCanvas( "c_HPt", "c_HPt", 800, 700 );
  TLegend* leg_HPt = new TLegend(0.7,0.7,0.9,0.9);
  c_HPt->SetHighLightColor(2);
  c_HPt->SetFillColor(0);
  c_HPt->SetBorderMode(0);
  c_HPt->SetBorderSize(2);
  c_HPt->SetLeftMargin( leftMargin );
  c_HPt->SetRightMargin( 1.6*rightMargin );
  c_HPt->SetTopMargin( topMargin );
  c_HPt->SetBottomMargin( bottomMargin );
  c_HPt->SetFrameBorderMode(0);
  c_HPt->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_HPt_str = "hist_HPt_"+Convert.str();
          hist_HPt[i] = new TH1F(hist_HPt_str.c_str(),"",100,0,1000);
          std::string draw_HPt_str = "HPt>>"+hist_HPt_str;
          tree[i]->Draw(draw_HPt_str.c_str(),"box > -1 && box < 5","HSAME");
          hist_HPt[i]->SetLineColor(i+1);
          double norm = hist_HPt[i]->Integral();
          hist_HPt[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_HPt->AddEntry(hist_HPt[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_HPt_str = "hist_HPt_"+Convert.str();
          hist_HPt[i] = new TH1F(hist_HPt_str.c_str(),"",100,0,1000);
          std::string draw_HPt_str = "HPt>>"+hist_HPt_str;
          tree[i-n]->Draw(draw_HPt_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_HPt[i]->SetLineColor(i+1);
          double norm = hist_HPt[i]->Integral();
          hist_HPt[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_HPt->AddEntry(hist_HPt[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_HPt[0]->SetTitle(";HPt;#Events");
  hist_HPt[0]->GetYaxis()->SetRangeUser(0,0.08);
  leg_HPt->SetBorderSize(0);
  leg_HPt->Draw("SAME");
  leg_HPt->Write();
  c_HPt->Write();
  c_HPt->Update();
  c_HPt->SaveAs("c_HPt.pdf");

  //pTGammaGamma
  TCanvas* c_pTGammaGamma = new TCanvas( "c_pTGammaGamma", "c_pTGammaGamma", 800, 700 );
  TLegend* leg_pTGammaGamma = new TLegend(0.7,0.7,0.9,0.9);
  c_pTGammaGamma->SetHighLightColor(2);
  c_pTGammaGamma->SetFillColor(0);
  c_pTGammaGamma->SetBorderMode(0);
  c_pTGammaGamma->SetBorderSize(2);
  c_pTGammaGamma->SetLeftMargin( leftMargin );
  c_pTGammaGamma->SetRightMargin( 1.6*rightMargin );
  c_pTGammaGamma->SetTopMargin( topMargin );
  c_pTGammaGamma->SetBottomMargin( bottomMargin );
  c_pTGammaGamma->SetFrameBorderMode(0);
  c_pTGammaGamma->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",100,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"box > -1 && box < 5","HSAME");
          hist_pTGammaGamma[i]->SetLineColor(i+1);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_pTGammaGamma->AddEntry(hist_pTGammaGamma[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",100,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i-n]->Draw(draw_pTGammaGamma_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_pTGammaGamma[i]->SetLineColor(i+1);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_pTGammaGamma->AddEntry(hist_pTGammaGamma[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[0]->SetTitle(";pTGammaGamma;#Events");
  hist_pTGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.12);
  leg_pTGammaGamma->SetBorderSize(0);
  leg_pTGammaGamma->Draw("SAME");
  leg_pTGammaGamma->Write();
  c_pTGammaGamma->Write();
  c_pTGammaGamma->Update();
  c_pTGammaGamma->SaveAs("c_pTGammaGamma.pdf");

  //HT
  TCanvas* c_HT = new TCanvas( "c_HT", "c_HT", 800, 700 );
  TLegend* leg_HT = new TLegend(0.7,0.7,0.9,0.9);
  c_HT->SetHighLightColor(2);
  c_HT->SetFillColor(0);
  c_HT->SetBorderMode(0);
  c_HT->SetBorderSize(2);
  c_HT->SetLeftMargin( leftMargin );
  c_HT->SetRightMargin( 1.6*rightMargin );
  c_HT->SetTopMargin( topMargin );
  c_HT->SetBottomMargin( bottomMargin );
  c_HT->SetFrameBorderMode(0);
  c_HT->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_HT_str = "hist_HT_"+Convert.str();
          hist_HT[i] = new TH1F(hist_HT_str.c_str(),"",100,0,2000);
          std::string draw_HT_str = "HT>>"+hist_HT_str;
          tree[i]->Draw(draw_HT_str.c_str(),"box > -1 && box < 5","HSAME");
          hist_HT[i]->SetLineColor(i+1);
          double norm = hist_HT[i]->Integral();
          hist_HT[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_HT->AddEntry(hist_HT[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_HT_str = "hist_HT_"+Convert.str();
          hist_HT[i] = new TH1F(hist_HT_str.c_str(),"",100,0,2000);
          std::string draw_HT_str = "HT>>"+hist_HT_str;
          tree[i-n]->Draw(draw_HT_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_HT[i]->SetLineColor(i+1);
          double norm = hist_HT[i]->Integral();
          hist_HT[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_HT->AddEntry(hist_HT[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_HT[0]->SetTitle(";HT;#Events");
  hist_HT[0]->GetYaxis()->SetRangeUser(0,0.12);
  leg_HT->SetBorderSize(0);
  leg_HT->Draw("SAME");
  leg_HT->Write();
  c_HT->Write();
  c_HT->Update();
  c_HT->SaveAs("c_HT.pdf");

  //MET
  TCanvas* c_MET = new TCanvas( "c_MET", "c_MET", 800, 700 );
  TLegend* leg_MET = new TLegend(0.7,0.7,0.9,0.9);
  c_MET->SetHighLightColor(2);
  c_MET->SetFillColor(0);
  c_MET->SetBorderMode(0);
  c_MET->SetBorderSize(2);
  c_MET->SetLeftMargin( leftMargin );
  c_MET->SetRightMargin( 1.6*rightMargin );
  c_MET->SetTopMargin( topMargin );
  c_MET->SetBottomMargin( bottomMargin );
  c_MET->SetFrameBorderMode(0);
  c_MET->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MET_str = "hist_MET_"+Convert.str();
          hist_MET[i] = new TH1F(hist_MET_str.c_str(),"",100,0,1000);
          std::string draw_MET_str = "MET>>"+hist_MET_str;
          tree[i]->Draw(draw_MET_str.c_str(),"box > -1 && box < 5","HSAME");
          hist_MET[i]->SetLineColor(i+1);
          double norm = hist_MET[i]->Integral();
          hist_MET[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i]+" leptons";
          leg_MET->AddEntry(hist_MET[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MET_str = "hist_MET_"+Convert.str();
          hist_MET[i] = new TH1F(hist_MET_str.c_str(),"",100,0,1000);
          std::string draw_MET_str = "MET>>"+hist_MET_str;
          tree[i-n]->Draw(draw_MET_str.c_str(),"box > 4 && box < 10","HSAME");
          hist_MET[i]->SetLineColor(i+1);
          double norm = hist_MET[i]->Integral();
          hist_MET[i]->Scale(1/norm);
          std::string add_entry_str = labelProcess[i-n]+" inclusive";
          leg_MET->AddEntry(hist_MET[i],add_entry_str.c_str(),"l");
          //std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MET[0]->SetTitle(";MET;#Events");
  hist_MET[0]->GetYaxis()->SetRangeUser(0,0.15);
  leg_MET->SetBorderSize(0);
  leg_MET->Draw("SAME");
  leg_MET->Write();
  c_MET->Write();
  c_MET->Update();
  c_MET->SaveAs("c_MET.pdf");

  fout->Close();


}
