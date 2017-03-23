//================================================================================================
//
// Skim
//
//________________________________________________________________________________________________
//

//#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include <THStack.h> 
#include <TKey.h> 
#include <TApplication.h>
//#endif


//=== MAIN MACRO ================================================================================================= 

int DoMixSamples( string inputfile1, string inputfile2, double weightFactor1, double weightFactor2, string outputfile) {

  //create output file
  TFile *outputFile = new TFile(outputfile.c_str(), "RECREATE");

  //loop over all TTrees in the file
  TFile *inputFile1 = TFile::Open(inputfile1.c_str(), "READ");
  assert(inputFile1);
  inputFile1->cd();
  inputFile1->Purge(); //purge unwanted TTree cycles in file
  TFile *inputFile2 = TFile::Open(inputfile2.c_str(), "READ");
  assert(inputFile2);
  inputFile2->cd();
  inputFile2->Purge(); //purge unwanted TTree cycles in file

  //  TIter nextkey(inputFile1->GetListOfKeys());
  //TKey *key;


  TTree *inputTree1 = (TTree*)inputFile1->Get("HggRazor");
  TTree *inputTree2 = (TTree*)inputFile2->Get("HggRazor");
  cout << "Processing tree " << inputTree1->GetName() << endl;
  
  //create new tree
  outputFile->cd();
  TTree *outputTree1 = inputTree1->CloneTree(0);  
  TTree *outputTree2 = inputTree2->CloneTree(0);  
  cout << "Events in the ntuple: " << inputTree1->GetEntries() << endl;
  
  uint run;
  uint event;
  float weight;
 
  //First File
  inputTree1->SetBranchAddress("weight", &weight);
    
  for (int n=0;n<inputTree1->GetEntries();n++) { 
    if (n%1000000==0) cout << "Processed Event " << n << "\n";
    inputTree1->GetEntry(n);
    weight = weight * weightFactor1;
    outputTree1->Fill(); 
  }
    
  //Second File
  inputTree2->SetBranchAddress("weight", &weight);	    
  for (int n=0;n<inputTree2->GetEntries();n++) { 
    if (n%1000000==0) cout << "Processed Event " << n << "\n";
    inputTree2->GetEntry(n);
    weight = weight * weightFactor2;
    outputTree2->Fill(); 
  }
    
  TList *list = new TList;
  list->Add(outputTree1);
  list->Add(outputTree2);
  TTree *outputTree = TTree::MergeTrees(list);
  outputTree->SetName("HggRazor");
  outputTree->Write();

  cout << "Number of Input Events From File 1: " << inputTree1->GetEntries() << "\n";
  cout << "Number of Input Events From File 2: " << inputTree2->GetEntries() << "\n";
  cout << "Number of Output Events In File: " << outputTree->GetEntries() << "\n";



  //Merge the Histograms
  TH1F *NEvents1 = (TH1F*)inputFile1->Get("NEvents");
  TH1F *SumWeights1 = (TH1F*)inputFile1->Get("SumWeights");
  TH1F *SumScaleWeights1 = (TH1F*)inputFile1->Get("SumScaleWeights");
  TH1F *SumPdfWeights1 = (TH1F*)inputFile1->Get("SumPdfWeights");
  TH1F *NISRJets1 = (TH1F*)inputFile1->Get("NISRJets");
  TH1F *PtISR1 = (TH1F*)inputFile1->Get("PtISR");
  TH1F *NPV1 = (TH1F*)inputFile1->Get("NPV");

  TH1F *NEvents2 = (TH1F*)inputFile1->Get("NEvents");
  TH1F *SumWeights2 = (TH1F*)inputFile1->Get("SumWeights");
  TH1F *SumScaleWeights2 = (TH1F*)inputFile1->Get("SumScaleWeights");
  TH1F *SumPdfWeights2 = (TH1F*)inputFile1->Get("SumPdfWeights");
  TH1F *NISRJets2 = (TH1F*)inputFile1->Get("NISRJets");
  TH1F *PtISR2 = (TH1F*)inputFile1->Get("PtISR");
  TH1F *NPV2 = (TH1F*)inputFile1->Get("NPV");

  TH1F *NEvents = (TH1F*)NEvents1->Clone("NEvents");
  TH1F *SumWeights = (TH1F*)SumWeights1->Clone("SumWeights");
  TH1F *SumScaleWeights = (TH1F*)SumScaleWeights1->Clone("SumScaleWeights");
  TH1F *SumPdfWeights = (TH1F*)SumPdfWeights1->Clone("SumPdfWeights");
  TH1F *NISRJets = (TH1F*)NISRJets1->Clone("NISRJets");
  TH1F *PtISR = (TH1F*)PtISR1->Clone("PtISR");
  TH1F *NPV = (TH1F*)NPV1->Clone("NPV");

  for (int i=0; i< NEvents->GetXaxis()->GetNbins()+2; i++) {
    NEvents->SetBinContent( i, weightFactor1 * NEvents1->GetBinContent(i) + 
			    weightFactor2 * NEvents2->GetBinContent(i));
  }

  for (int i=0; i< SumWeights->GetXaxis()->GetNbins()+2; i++) {
    SumWeights->SetBinContent( i, weightFactor1 * SumWeights1->GetBinContent(i) + 
			    weightFactor2 * SumWeights2->GetBinContent(i));
  }

  for (int i=0; i< SumScaleWeights->GetXaxis()->GetNbins()+2; i++) {
    SumScaleWeights->SetBinContent( i, weightFactor1 * SumScaleWeights1->GetBinContent(i) + 
			    weightFactor2 * SumScaleWeights2->GetBinContent(i));
  }

  for (int i=0; i< SumPdfWeights->GetXaxis()->GetNbins()+2; i++) {
    SumPdfWeights->SetBinContent( i, weightFactor1 * SumPdfWeights1->GetBinContent(i) + 
			    weightFactor2 * SumPdfWeights2->GetBinContent(i));
  }

  for (int i=0; i< NISRJets->GetXaxis()->GetNbins()+2; i++) {
    NISRJets->SetBinContent( i, weightFactor1 * NISRJets1->GetBinContent(i) + 
			    weightFactor2 * NISRJets2->GetBinContent(i));
  }

  for (int i=0; i< PtISR->GetXaxis()->GetNbins()+2; i++) {
    PtISR->SetBinContent( i, weightFactor1 * PtISR1->GetBinContent(i) + 
			    weightFactor2 * PtISR2->GetBinContent(i));
  }

  for (int i=0; i< NPV->GetXaxis()->GetNbins()+2; i++) {
    NPV->SetBinContent( i, weightFactor1 * NPV1->GetBinContent(i) + 
			    weightFactor2 * NPV2->GetBinContent(i));
  }

  NEvents->Write(); 
  SumWeights->Write(); 
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  NISRJets->Write();
  PtISR->Write();
  NPV->Write();
	
  inputFile1->Close();
  inputFile2->Close();
  cout << "Closing output file." << endl;
  outputFile->Close();
  delete outputFile;
  delete list;
  //gApplication->Terminate();
  return 0;
}


void MixSignalSamples() {

  int samples[36] = { 127,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000};

  // for (int i=29; i<36; i++) {
  //   cout << "Sample : " << i << " : " << samples[i] << "\n";
  //   DoMixSamples(Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHH_%i_1pb_weighted.root",samples[i]),Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHZ_%i_1pb_weighted.root",samples[i]), 0.25, 0.5, Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHHHZ_BR5050_%i_1pb_weighted.root",samples[i]));    
  // }


 // for (int i=29; i<36; i++) {
 //    cout << "Sample : " << i << " : " << samples[i] << "\n";
 //    DoMixSamples(Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHH_%i_1pb_weighted.root",samples[i]),Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHZ_%i_1pb_weighted.root",samples[i]), 0.5625, 0.375, Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHHHZ_BR7525_%i_1pb_weighted.root",samples[i]));    
 //  }

 for (int i=29; i<36; i++) {
    cout << "Sample : " << i << " : " << samples[i] << "\n";
    DoMixSamples(Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHH_%i_1pb_weighted.root",samples[i]),Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHZ_%i_1pb_weighted.root",samples[i]), 0.0625, 0.375, Form("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/2016/V3p12_PhotonCorrDec06_JECSep23V3_20170219/FastsimSignal/combined/SMS-TChiHHHZ_BR2575_%i_1pb_weighted.root",samples[i]));    
  }



}


