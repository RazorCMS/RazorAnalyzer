#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include <string>
#include <iostream>
#include <math.h>

void TMVA_vtx_prepareTree(string inputFilename = "../HggRazorUpgradeTiming_PU0_NoTiming.root", string outputFilename = "HggRazorUpgradeTiming_PU0_NoTiming_vtx.root")
{
	float ptasym = 0.;
  	float ptbal = 0.;
  	float logsumpt2 = 0.;
  	float pull_conv = 0.;
  	float nConv = 0.;
  	float weight = 1.0;

	Int_t nPVAll = 0;
	Int_t isMatchPv[200];
  	Int_t isMaxbdtPv[200];
   	float ptasym_all[200];// = {0.};
  	float ptbal_all[200];// = {0.};
  	float logsumpt2_all[200];// = {0.};
  	float pull_conv_all[200];// = {0.};
  	float nConv_all[200];// = {0.};

	for(int i=0;i<200;i++)
  	{
        isMatchPv[i]=0;
        isMaxbdtPv[i]=0;
        ptasym_all[i]=0.0;
        ptbal_all[i]=0.0;
        logsumpt2_all[i]=0.0;
        pull_conv_all[i]=0.0;
        nConv_all[i]=0.0;
  	}


   TFile *file_in = new TFile(inputFilename.c_str(),"READ");
   TTree *tree_in = (TTree*)file_in->Get("HggRazor");
    
   tree_in->SetBranchAddress( "nPVAll", 	&nPVAll);
   tree_in->SetBranchAddress( "isMatchPv",   	isMatchPv);
   tree_in->SetBranchAddress( "ptasym",   	ptasym_all);
   tree_in->SetBranchAddress( "ptbal",   	ptbal_all);
   tree_in->SetBranchAddress( "logsumpt2",   	logsumpt2_all);
   tree_in->SetBranchAddress( "limPullToConv",  pull_conv_all);
   tree_in->SetBranchAddress( "nConv",   	nConv_all);

   TFile *file_out = new TFile(outputFilename.c_str(),"RECREATE");
   TTree *tree_out_S = new TTree("TreeS", "Tree_Signal");
   TTree *tree_out_B = new TTree("TreeB", "Tree_Background");

   tree_out_S->Branch( "ptasym",	&ptasym, "ptasym/F");
   tree_out_S->Branch( "ptbal", 	&ptbal, "ptbal/F");
   tree_out_S->Branch( "logsumpt2", 	&logsumpt2, "logsumpt2/F");
   tree_out_S->Branch( "limPullToConv", &pull_conv, "limPullToConv/F");
   tree_out_S->Branch( "nConv", 	&nConv, "nConv/F");
   //tree_out_S->Branch( "weight", 	&weight, "weight/F");

   tree_out_B->Branch( "ptasym",	&ptasym, "ptasym/F");
   tree_out_B->Branch( "ptbal", 	&ptbal, "ptbal/F");
   tree_out_B->Branch( "logsumpt2", 	&logsumpt2, "logsumpt2/F");
   tree_out_B->Branch( "limPullToConv", &pull_conv, "limPullToConv/F");
   tree_out_B->Branch( "nConv", 	&nConv, "nConv/F");
   tree_out_B->Branch( "weight", 	&weight, "weight/F");


   int NEntries =  tree_in->GetEntries();
   for(int i=0;i<NEntries;i++)
   {
	tree_in->GetEntry(i);
	for(int j=0;j<nPVAll;j++)
	{	
		ptasym = ptasym_all[j];
		ptbal = ptbal_all[j];
		logsumpt2 = logsumpt2_all[j];
		pull_conv = pull_conv_all[j];
		nConv = nConv_all[j];
		if(!isinf(logsumpt2))
		{
		if(isMatchPv[j]==0)
		{
			tree_out_B->Fill();
		}		
		else
		{
			tree_out_S->Fill();
		}
		}
	}
   }

   tree_out_B->Write();
   tree_out_S->Write();	

}
