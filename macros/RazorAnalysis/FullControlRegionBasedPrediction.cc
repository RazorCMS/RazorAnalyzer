//Runs on the output of the RazorInclusive analyzer and gives the MC-based background prediction in each bin of the MR-Rsq plane

#include <iostream>
#include <map>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"

#include "include/RazorAnalyzer.h"

using namespace std;

//check if the given box is a muon box
bool isSingleMuonBox(RazorAnalyzer::RazorBox box){
    if(box == RazorAnalyzer::MuSixJet || box == RazorAnalyzer::MuFourJet || box == RazorAnalyzer::MuJet || box == RazorAnalyzer::MuMultiJet) return true;
    return false;
}
//check if the given box is a electron box
bool isSingleElectronBox(RazorAnalyzer::RazorBox box){
    if(box == RazorAnalyzer::EleSixJet || box == RazorAnalyzer::EleFourJet || box == RazorAnalyzer::EleJet || box == RazorAnalyzer::EleMultiJet) return true;
    return false;
}

//define MR and Rsq binning

//int NMRBINS = 10;
//float MRBINLOWEDGES[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
int NRSQBINS = 8;
float RSQBINLOWEDGES[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.8, 1.5};

//even binning
int NMRBINS = 20;
float MRBINLOWEDGES[] = {300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600, 630, 660, 690, 720, 750, 780, 810, 840, 870, 900};
//float MRBINLOWEDGES[] = {200, 230, 260, 290, 320, 350, 380, 410, 440, 470, 500, 530, 560, 590, 620, 650, 680, 710, 740, 770, 800};
//float MRBINLOWEDGES[] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200};

//ttbar single lepton SF bins
//int NMRBINS = 8;
//int NRSQBINS = 7;
//float MRBINLOWEDGES[] = {300, 350, 400, 450, 500, 550, 700, 900, 4000};
//float RSQBINLOWEDGES[] = {0.15,0.175,0.20,0.225,0.25,0.30,0.41, 1.5};

//wjets single lepton SF bins
//const int NMRBINS = 9;
//const int NRSQBINS = 8;
//double MRBINLOWEDGES[] = {300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000};
//double RSQBINLOWEDGES[] = {0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5};  

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void FullControlRegionBasedPrediction(){
    bool doSFCorrections = true; //apply TT, W, Z, DY scale factors
    bool doMiscCorrections = true; //apply lepton efficiency, b-tagging, ... scale factors
    //bool scaleZNuNuToDY = true; //scale Z->nunu MC razor variable distribution to match that of DY+Jets MC
    bool scaleZNuNuToDY = false; //scale Z->nunu MC razor variable distribution to match that of DY+Jets MC
    gROOT->SetBatch();

    //////////////////////////////////////////////////
    //Define baseline cuts
    //////////////////////////////////////////////////

    bool doDPhiRazorCut = true;
    //bool doDPhiRazorCut = false;
    float dPhiRazorCut = 2.7; //cut on the angle between the two razor hemispheres

    bool doMetCut = false;
    //bool doMetCut = true;
    float metCut = 30;

    //bool doMTCut = true;
    bool doMTCut = false;
    float mTLowerCut = 30;
    float mTUpperCut = 100;

    //bool doLeptonPtCut = true;
    bool doLeptonPtCut = false;
    float leptonPtCut = 30;

    //bool bTagsInclusive = true; //true = require >= minNBTags, false = require = minNBTags
    bool bTagsInclusive = false; //true = require >= minNBTags, false = require = minNBTags
    int minNBTags = 0; //TODO: bin in nBTags instead of cutting

    //set color palette 
    const Int_t NCont = 101;
    gStyle->SetNumberContours(NCont);

    //define cuts for 1D MR and Rsq plots
    float MRCutFor1DPlots = 400;
    float RsqCutFor1DPlots = 0.25;

    //get input files -- output of RazorInclusive analyzer
    int lumiInData = 19700; //in /pb
    int lumiInMC = 1; //luminosity used to normalize MC ntuples
    string mcPrefix;
    if(doMiscCorrections){
        //NOTE: all data-MC correction factors should already be applied EXCEPT for the hadronic recoil scale factors obtained from the control regions 
        mcPrefix = "/afs/cern.ch/work/d/duanders/run2Studies/CMSSW_7_3_0_pre1/src/RazorAnalyzer/normtest/June10/";
        //mcPrefix = "eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorInclusive/done/MC_WithCorrectionFactors"; //location of MC ntuples
    }
    else{
        mcPrefix = "eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorInclusive/done/MC_NoCorrectionFactors/";//location of MC ntuples
    }
    //string dataPrefix = "/afs/cern.ch/work/d/duanders/run2Studies/CMSSW_7_3_0_pre1/src/RazorAnalyzer/normtest/June10/"; //location of data ntuples
    string dataPrefix = "eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorInclusive/done/June10/"; //location of data ntuples

    map<string, TFile*> mcfiles;
    //mcfiles["DYJets"] = new TFile(Form("%s/RazorInclusive_DYJetsToLL_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["WJets"] = new TFile(Form("%s/RazorInclusive_WJetsToLNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["ZJetsNuNu"] = new TFile(Form("%s/RazorInclusive_ZJetsToNuNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["TTJets"] = new TFile(Form("%s/RazorInclusive_TTJets_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["SingleTop"] = new TFile(Form("%s/RazorInclusive_SingleTop_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["QCD"] = new TFile(Form("%s/RazorInclusive_QCD_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["TTV"] = new TFile(Form("%s/RazorInclusive_TTV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["VV"] = new TFile(Form("%s/RazorInclusive_VV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    //mcfiles["TTTT"] = new TFile(Form("%s/RazorInclusive_TTTT_TuneZ2star_8TeV-madgraph-tauola_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["DYJets"] = new TFile(Form("%s/DYJetsToLL_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["WJets"] = new TFile(Form("%s/WJetsToLNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["ZJetsNuNu"] = new TFile(Form("%s/ZJetsToNuNu_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTJets"] = new TFile(Form("%s/TTJets_All_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["SingleTop"] = new TFile(Form("%s/SingleTop_All_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["QCD"] = new TFile(Form("%s/QCD_HTBinned_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTV"] = new TFile(Form("%s/TTV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["VV"] = new TFile(Form("%s/VV_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));
    mcfiles["TTTT"] = new TFile(Form("%s/TTTT_TuneZ2star_8TeV-madgraph-tauola_%dpb_weighted.root", mcPrefix.c_str(), lumiInMC));

    //data
    //TFile *datafile;
    //datafile = new TFile(Form("%s/RazorInclusive_Data_HTMHTParked_Run2012_GoodLumi.root", dataPrefix.c_str()));
    map<string, TFile*> datafiles;
    vector<string> datanames{"HTMHT", "SingleMu", "SingleElectron", "DoubleMuParked", "DoubleElectron", "MuEG"};
    for(auto &name : datanames) datafiles[name] = new TFile(Form("%s/Data_%s_GoodLumi.root", dataPrefix.c_str(), name.c_str()));

    //get trees and set branches
    map<string, TTree*> mctrees;
    //TTree *datatree;
    map<string, TTree*> datatrees;
    float weight;
    float MR, Rsq, dPhiRazor, met, mT, leadingTightMuPt, leadingTightElePt;
    int nBTaggedJets, nSelectedJets, box;
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("box", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("dPhiRazor", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("nSelectedJets", 1);
        mctrees[file.first]->SetBranchStatus("met", 1);
        mctrees[file.first]->SetBranchStatus("mT", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightMuPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightElePt", 1);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("box", &box);
        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("dPhiRazor", &dPhiRazor);
        mctrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
        mctrees[file.first]->SetBranchAddress("nSelectedJets", &nSelectedJets);
        mctrees[file.first]->SetBranchAddress("met", &met);
        mctrees[file.first]->SetBranchAddress("mT", &mT);
        mctrees[file.first]->SetBranchAddress("leadingTightMuPt", &leadingTightMuPt);
        mctrees[file.first]->SetBranchAddress("leadingTightElePt", &leadingTightElePt);
    }
    //datatree = (TTree*)datafile->Get("RazorInclusive");
    //datatree->SetBranchStatus("*", 0);
    //datatree->SetBranchStatus("box", 1);
    //datatree->SetBranchStatus("MR", 1);
    //datatree->SetBranchStatus("Rsq", 1);
    //datatree->SetBranchStatus("dPhiRazor", 1);
    //datatree->SetBranchStatus("nBTaggedJets", 1);
    //datatree->SetBranchStatus("nSelectedJets", 1);
    //datatree->SetBranchStatus("met", 1);

    //datatree->SetBranchAddress("box", &box);
    //datatree->SetBranchAddress("MR", &MR);
    //datatree->SetBranchAddress("Rsq", &Rsq);
    //datatree->SetBranchAddress("dPhiRazor", &dPhiRazor);
    //datatree->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
    //datatree->SetBranchAddress("met", &met);
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        datatrees[file.first]->SetBranchStatus("*", 0);
        datatrees[file.first]->SetBranchStatus("box", 1);
        datatrees[file.first]->SetBranchStatus("MR", 1);
        datatrees[file.first]->SetBranchStatus("Rsq", 1);
        datatrees[file.first]->SetBranchStatus("dPhiRazor", 1);
        datatrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        datatrees[file.first]->SetBranchStatus("nSelectedJets", 1);
        datatrees[file.first]->SetBranchStatus("met", 1);
        datatrees[file.first]->SetBranchStatus("mT", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightElePt", 1);

        datatrees[file.first]->SetBranchAddress("box", &box);
        datatrees[file.first]->SetBranchAddress("MR", &MR);
        datatrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        datatrees[file.first]->SetBranchAddress("dPhiRazor", &dPhiRazor);
        datatrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
        datatrees[file.first]->SetBranchAddress("met", &met);
        datatrees[file.first]->SetBranchAddress("mT", &mT);
        datatrees[file.first]->SetBranchAddress("leadingTightMuPt", &leadingTightMuPt);
        datatrees[file.first]->SetBranchAddress("leadingTightElePt", &leadingTightElePt);
    }

    //load TTbar scale factor histograms
    TFile *SFFileTTJets = new TFile("data/ScaleFactors/Run1/TTBarSingleLeptonScaleFactors.root");
    TH2F *SFHistTTJets = (TH2F*)SFFileTTJets->Get("TTBarSingleLeptonScaleFactor");
    float SFmaxMRTTJets = SFHistTTJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqTTJets = SFHistTTJets->GetYaxis()->GetXmax() - 0.01;
    //cout << "TTJets " << SFmaxMRTTJets << " " << SFmaxRsqTTJets << endl;

    //load WJets scale factor histogram
    TFile *SFFileWJets = new TFile("data/ScaleFactors/Run1/WJetsSingleLeptonScaleFactors.root");
    TH2F *SFHistWJets = (TH2F*)SFFileWJets->Get("WJetsSingleLeptonScaleFactor");
    float SFmaxMRWJets = SFHistWJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqWJets = SFHistWJets->GetYaxis()->GetXmax() - 0.01;
    //cout << "WJets " << SFmaxMRWJets << " " << SFmaxRsqWJets << endl;

    //load DYJets scale factor histogram
    TFile *SFFileDYJets = new TFile("data/ScaleFactors/Run1/ZToLLScaleFactors.root");
    TH2F *SFHistDYJets = (TH2F*)SFFileDYJets->Get("ZToLLDileptonScaleFactor");
    float SFmaxMRDYJets = SFHistDYJets->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqDYJets = SFHistDYJets->GetYaxis()->GetXmax() - 0.01;
    //cout << "DYJets " << SFmaxMRDYJets << " " << SFmaxRsqDYJets << endl;

    //load ZNuNu scale factor histograms
    TFile *SFFileZJetsNuNu = new TFile("data/ScaleFactors/Run1/ZInvisibleScaleFactorsRun1.root");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("DYJetsScaleFactors");
    //TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("WJetsScaleFactors");
    TH2F *SFHistZJetsNuNu = (TH2F*)SFFileZJetsNuNu->Get("GJetsScaleFactors");
    float SFmaxMRZJetsNuNu = SFHistZJetsNuNu->GetXaxis()->GetXmax() - 1;
    float SFmaxRsqZJetsNuNu = SFHistZJetsNuNu->GetYaxis()->GetXmax() - 0.01;
    //cout << "ZJetsNuNu " << SFmaxMRZJetsNuNu << " " << SFmaxRsqZJetsNuNu << endl;

    //load ZNuNu-->DYJets weighting factors
    TFile *ZNuNuToDYWeightFile = new TFile("data/ScaleFactors/Run1/ZNuNuToDYScaleFactorsRun1.root");
    TH2F *ZNuNuToDYWeightHist = (TH2F*)ZNuNuToDYWeightFile->Get("razormcDYJets");
    float maxMRZNuNuToDY = ZNuNuToDYWeightHist->GetXaxis()->GetXmax() - 1;
    float maxRsqZNuNuToDY = ZNuNuToDYWeightHist->GetYaxis()->GetXmax() - 0.01;
    //cout << "ZNuNuToDY " << maxMRZNuNuToDY << " " << maxRsqZNuNuToDY << endl;

    //declare which boxes to check
    map<RazorAnalyzer::RazorBox, string> boxes;
    boxes[RazorAnalyzer::MuEle] = "MuEle";
    boxes[RazorAnalyzer::MuMu] = "MuMu";
    boxes[RazorAnalyzer::EleEle] = "EleEle";
    boxes[RazorAnalyzer::MuSixJet] = "MuSixJet";
    boxes[RazorAnalyzer::MuFourJet] = "MuFourJet";
    boxes[RazorAnalyzer::MuJet] = "MuJet";
    boxes[RazorAnalyzer::EleSixJet] = "EleSixJet";
    boxes[RazorAnalyzer::EleFourJet] = "EleFourJet";
    boxes[RazorAnalyzer::EleJet] = "EleJet";
    boxes[RazorAnalyzer::LooseLeptonSixJet] = "LooseLeptonSixJet";
    boxes[RazorAnalyzer::LooseLeptonFourJet] = "LooseLeptonFourJet";
    boxes[RazorAnalyzer::LooseLeptonDiJet] = "LooseLeptonDiJet";
    boxes[RazorAnalyzer::SixJet] = "SixJet";
    boxes[RazorAnalyzer::FourJet] = "FourJet";
    boxes[RazorAnalyzer::DiJet] = "DiJet";
    boxes[RazorAnalyzer::NONE] = "MultiJetPlusLooseLeptonMultiJet"; //(a temporary hack to combine three boxes)
    //boxes[RazorAnalyzer::NONE] = "MultiJet"; //(a temporary hack to combine three boxes)
    //boxes[RazorAnalyzer::NONE] = "WJetsSingleLepton"; //(a temporary hack to combine three boxes)
    //boxes[RazorAnalyzer::NONE] = "TTJetsSingleLepton"; //(a temporary hack to combine three boxes)
    //associate each box with a dataset
    map<RazorAnalyzer::RazorBox, string> boxDatasets;
    boxDatasets[RazorAnalyzer::MuEle] = "MuEG";
    boxDatasets[RazorAnalyzer::MuMu] = "DoubleMuParked";
    boxDatasets[RazorAnalyzer::EleEle] = "DoubleElectron";
    boxDatasets[RazorAnalyzer::MuSixJet] = "SingleMu";
    boxDatasets[RazorAnalyzer::MuFourJet] = "SingleMu";
    boxDatasets[RazorAnalyzer::MuJet] = "SingleMu";
    boxDatasets[RazorAnalyzer::EleSixJet] = "SingleElectron";
    boxDatasets[RazorAnalyzer::EleFourJet] = "SingleElectron";
    boxDatasets[RazorAnalyzer::EleJet] = "SingleElectron";
    boxDatasets[RazorAnalyzer::LooseLeptonSixJet] = "HTMHT";
    boxDatasets[RazorAnalyzer::LooseLeptonFourJet] = "HTMHT";
    boxDatasets[RazorAnalyzer::LooseLeptonDiJet] = "HTMHT";
    boxDatasets[RazorAnalyzer::SixJet] = "HTMHT";
    boxDatasets[RazorAnalyzer::FourJet] = "HTMHT";
    boxDatasets[RazorAnalyzer::DiJet] = "HTMHT";

    //make directories for plots
    string plotDir = "/afs/cern.ch/work/d/duanders/public/plotsNoZNuNuSF0BtagDPhiCut";
    struct stat st;
    if (stat(plotDir.c_str(), &st) == -1) {
        mkdir(plotDir.c_str(), 0777);
    }
    if (stat(Form("%s/RazorInclusive", plotDir.c_str()), &st) == -1) {
        mkdir(Form("%s/RazorInclusive", plotDir.c_str()), 0777);
    }
    for(auto &ibox : boxes){
        //check if directory exists
        if (stat(Form("%s/RazorInclusive/%s", plotDir.c_str(), ibox.second.c_str()), &st) == -1) {
            mkdir(Form("%s/RazorInclusive/%s", plotDir.c_str(), ibox.second.c_str()), 0777);
        }
    }       
    string plotPath = plotDir+"/RazorInclusive";

    //////////////////////////////////////////////////
    //Make histograms for each box, and fill them
    //////////////////////////////////////////////////

    //////////////////////////////////////////////////
    //Step 1: Get the predictions from each MC process
    //////////////////////////////////////////////////
    map<RazorAnalyzer::RazorBox, map<string, TH2F> > razorHistosMC;
    map<RazorAnalyzer::RazorBox, map<string, TH1F> > MRHistosMC;
    map<RazorAnalyzer::RazorBox, map<string, TH1F> > RsqHistosMC;
    map<RazorAnalyzer::RazorBox, map<string, TH2F> > razorErrorHistosMC;
    map<RazorAnalyzer::RazorBox, map<string, TH1F> > MRErrorHistosMC;
    map<RazorAnalyzer::RazorBox, map<string, TH1F> > RsqErrorHistosMC;
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;

        //set up histograms
        for(auto &ibox : boxes){
            razorHistosMC[ibox.first][tree.first] = TH2F(Form("razormc%s%s", tree.first.c_str(), ibox.second.c_str()), "; MR (GeV); Rsq", NMRBINS, MRBINLOWEDGES, NRSQBINS, RSQBINLOWEDGES);
            MRHistosMC[ibox.first][tree.first] = TH1F(Form("mrmc%s%s", tree.first.c_str(), ibox.second.c_str()), "; MR (GeV)", NMRBINS, MRBINLOWEDGES);
            RsqHistosMC[ibox.first][tree.first] = TH1F(Form("rsqmc%s%s", tree.first.c_str(), ibox.second.c_str()), "; Rsq", NRSQBINS, RSQBINLOWEDGES);
            MRHistosMC[ibox.first][tree.first].Sumw2();
            RsqHistosMC[ibox.first][tree.first].Sumw2();
            razorHistosMC[ibox.first][tree.first].Sumw2();
            //histograms to hold sum(w^2*error(SF)^2) for each bin
            razorErrorHistosMC[ibox.first][tree.first] = TH2F(Form("razorErrormc%s%s", tree.first.c_str(), ibox.second.c_str()), "sum(w^2*error(SF)^2); MR (GeV); Rsq", NMRBINS, MRBINLOWEDGES, NRSQBINS, RSQBINLOWEDGES);
            MRErrorHistosMC[ibox.first][tree.first] = TH1F(Form("mrErrormc%s%s", tree.first.c_str(), ibox.second.c_str()), "sum(w^2*error(SF)^2); MR (GeV)", NMRBINS, MRBINLOWEDGES);
            RsqErrorHistosMC[ibox.first][tree.first] = TH1F(Form("rsqErrormc%s%s", tree.first.c_str(), ibox.second.c_str()), "sum(w^2*error(SF)^2); Rsq", NRSQBINS, RSQBINLOWEDGES);
        }

        uint nEntries = tree.second->GetEntries();
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 
            RazorAnalyzer::RazorBox razorbox = static_cast<RazorAnalyzer::RazorBox>(box);

            //enforce correct number of B-tags
            if(bTagsInclusive){
                if(nBTaggedJets < minNBTags) continue;
            }
            else {
                if(nBTaggedJets != minNBTags) continue;
            }
            //cut on MR and Rsq
            if(MR < MRBINLOWEDGES[0] || Rsq < RSQBINLOWEDGES[0]) continue;
            //cut on dPhiRazor
            if(doDPhiRazorCut && fabs(dPhiRazor) > dPhiRazorCut) continue;
            //cut on met
            if(doMetCut && met < metCut) continue;
            //cut on mT
            if(doMTCut && (isSingleMuonBox(razorbox) || isSingleElectronBox(razorbox))){
                if(mT < mTLowerCut || mT > mTUpperCut) continue;
            }
            //cut on lepton pt
            if(doLeptonPtCut){
                if(isSingleMuonBox(razorbox) && leadingTightMuPt < leptonPtCut) continue;
                if(isSingleElectronBox(razorbox) && leadingTightElePt < leptonPtCut) continue;
            }

            float eventWeight = weight*lumiInData*1.0/lumiInMC;
            float sysErrorSquared = 0.0;

            //Data/MC scale factors
            if(doSFCorrections){
                //TTJets SF
                if(tree.first == "TTJets"){
                    double SFTTJets = SFHistTTJets->GetBinContent(SFHistTTJets->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                    double SFErrorTTJets = SFHistTTJets->GetBinError(SFHistTTJets->FindFixBin(min(MR, SFmaxMRTTJets), min(Rsq, SFmaxRsqTTJets)));
                    if(SFTTJets < 1e5){
                        eventWeight *= SFTTJets;
                        sysErrorSquared += weight*weight*SFErrorTTJets*SFErrorTTJets;
                        //cout << "TTJets SF: " << SFTTJets << " (MR, Rsq) = (" << MR << ", " << Rsq << ") " << endl;
                    }
                    else{
                        //cout << "Warning: TTJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //WJets SF
                else if(tree.first == "WJets"){
                    double SFWJets = SFHistWJets->GetBinContent(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                    double SFErrorWJets = SFHistWJets->GetBinError(SFHistWJets->FindFixBin(min(MR, SFmaxMRWJets), min(Rsq, SFmaxRsqWJets)));
                    if(SFWJets < 1e5){
                        eventWeight *= SFWJets;
                        sysErrorSquared += weight*weight*SFErrorWJets*SFErrorWJets;
                        //cout << "WJets SF: " << SFWJets << " (MR, Rsq) = (" << MR << ", " << Rsq << ") " << endl;
                    }
                    else{
                        //cout << "Warning: WJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //DYJets SF
                else if(tree.first == "DYJets"){
                    double SFDYJets = SFHistDYJets->GetBinContent(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                    double SFErrorDYJets = SFHistDYJets->GetBinError(SFHistDYJets->FindFixBin(min(MR, SFmaxMRDYJets), min(Rsq, SFmaxRsqDYJets)));
                    if(SFDYJets < 1e5){
                        eventWeight *= SFDYJets;
                        sysErrorSquared += weight*weight*SFErrorDYJets*SFErrorDYJets;
                        //cout << "DYJets SF: " << SFDYJets << " (MR, Rsq) = (" << MR << ", " << Rsq << ") " << endl;
                    }
                    else{
                        //cout << "Warning: DYJets scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }
                }
                //ZNuNu SF
                //TODO: combine the three predictions for ZNuNu
                /*else if(tree.first == "ZJetsNuNu"){
                    double SFZJetsNuNu = SFHistZJetsNuNu->GetBinContent(SFHistZJetsNuNu->FindFixBin(min(MR, SFmaxMRZJetsNuNu), min(Rsq, SFmaxRsqZJetsNuNu)));
                    double SFErrorZJetsNuNu = SFHistZJetsNuNu->GetBinError(SFHistZJetsNuNu->FindFixBin(min(MR, SFmaxMRZJetsNuNu), min(Rsq, SFmaxRsqZJetsNuNu)));
                    if(SFZJetsNuNu < 1e5){
                        eventWeight *= SFZJetsNuNu;
                        sysErrorSquared += weight*weight*SFErrorZJetsNuNu*SFErrorZJetsNuNu;
                        //cout << "ZJetsNuNu SF: " << SFZJetsNuNu << " (MR, Rsq) = (" << MR << ", " << Rsq << ") " << endl;
                    }
                    else{
                        //cout << "Warning: ZJetsNuNu scale factor is Inf!" << endl;
                        eventWeight = 0;
                        sysErrorSquared = 0;
                    }

                    //scale ZNuNu so it looks like DYJets
                    if(scaleZNuNuToDY){
                        double SFZNuNuToDY = ZNuNuToDYWeightHist->GetBinContent(ZNuNuToDYWeightHist->FindFixBin(min(MR, maxMRZNuNuToDY), min(Rsq, maxRsqZNuNuToDY)));
                        double SFErrorZNuNuToDY = ZNuNuToDYWeightHist->GetBinError(ZNuNuToDYWeightHist->FindFixBin(min(MR, maxMRZNuNuToDY), min(Rsq, maxRsqZNuNuToDY)));
                        if(SFZNuNuToDY < 1e5){
                            //cout << "ZNuNuToDY SF: " << SFZNuNuToDY << " (MR, Rsq) = (" << MR << ", " << Rsq << ") " << endl;
                            eventWeight *= SFZNuNuToDY;
                            sysErrorSquared += weight*weight*SFErrorZNuNuToDY*SFErrorZNuNuToDY;
                        }
                        else{
                            //cout << "Warning: ZNuNuToDY scale factor is Inf!" << endl;
                            eventWeight = 0;
                            sysErrorSquared = 0;
                        }
                    }
                }*/
            }

            //single lepton trigger scale factor
            if(isSingleMuonBox(razorbox) || isSingleElectronBox(razorbox)){
                eventWeight = eventWeight*0.97;
            }

            //fill each quantity
            if(boxes.find(razorbox) == boxes.end()){
                cout << "Box " << razorbox << " not in the list!" << endl;
                continue;
            }
            razorHistosMC[razorbox][tree.first].Fill(MR, Rsq, eventWeight);
            razorErrorHistosMC[razorbox][tree.first].Fill(MR, Rsq, sysErrorSquared);
            if(Rsq > RsqCutFor1DPlots){
                MRHistosMC[razorbox][tree.first].Fill(MR, eventWeight);
                MRErrorHistosMC[razorbox][tree.first].Fill(MR, sysErrorSquared);
            }
            if(MR > MRCutFor1DPlots){
                RsqHistosMC[razorbox][tree.first].Fill(Rsq, eventWeight);
                RsqErrorHistosMC[razorbox][tree.first].Fill(Rsq, sysErrorSquared);
            }

            //combined Single Lepton boxes
            //if(isSingleMuonBox(razorbox) || isSingleElectronBox(razorbox)){
            //MultiJet box
            //if(razorbox == RazorAnalyzer::FourJet || razorbox == RazorAnalyzer::SixJet){
            //LooseLeptonMultiJet+MultiJet box
            if(razorbox == RazorAnalyzer::FourJet || razorbox == RazorAnalyzer::SixJet || razorbox == RazorAnalyzer::LooseLeptonFourJet || razorbox == RazorAnalyzer::LooseLeptonSixJet){
                razorHistosMC[RazorAnalyzer::NONE][tree.first].Fill(MR, Rsq, eventWeight);
                razorErrorHistosMC[RazorAnalyzer::NONE][tree.first].Fill(MR, Rsq, sysErrorSquared);
                if(Rsq > RsqCutFor1DPlots){
                    MRHistosMC[RazorAnalyzer::NONE][tree.first].Fill(MR, eventWeight);
                    MRErrorHistosMC[RazorAnalyzer::NONE][tree.first].Fill(MR, sysErrorSquared);
                }
                if(MR > MRCutFor1DPlots){
                    RsqHistosMC[RazorAnalyzer::NONE][tree.first].Fill(Rsq, eventWeight);
                    RsqErrorHistosMC[RazorAnalyzer::NONE][tree.first].Fill(Rsq, sysErrorSquared);
                }
            }
        }
    }
    //update errors to take into account systematic uncertainties
    for(auto &tree : mctrees){
        for(auto &ibox : boxes){
            for(int i = 0; i < razorHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosMC[ibox.first][tree.first].GetNbinsY()+1; j++){
                    double squaredError = razorErrorHistosMC[ibox.first][tree.first].GetBinContent(i, j);
                    razorHistosMC[ibox.first][tree.first].SetBinError(i, j, sqrt(pow(razorHistosMC[ibox.first][tree.first].GetBinError(i, j), 2) + squaredError));
                }
            }
            for(int i = 0; i < MRHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                double squaredError = MRErrorHistosMC[ibox.first][tree.first].GetBinContent(i);
                MRHistosMC[ibox.first][tree.first].SetBinError(i, sqrt(pow(MRHistosMC[ibox.first][tree.first].GetBinError(i), 2) + squaredError));
            }
            for(int i = 0; i < RsqHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                double squaredError = RsqErrorHistosMC[ibox.first][tree.first].GetBinContent(i);
                RsqHistosMC[ibox.first][tree.first].SetBinError(i, sqrt(pow(RsqHistosMC[ibox.first][tree.first].GetBinError(i), 2) + squaredError));
            }
        }
    }

    //for rare background processes, include a 20% uncertainty on the total yield in each bin, summed in quadrature with the statistical uncertainty
    double sysErrorFrac = 0.2;
    //for QCD, assign a 100% uncertainty
    double qcdErrorFrac = 1.0;
    for(auto &tree : mctrees){
        //only do this for rare processes 
        if(tree.first == "DYJets" || tree.first == "WJets" || tree.first == "ZJetsNuNu" || tree.first == "TTJets") continue; 
        for(auto &ibox : boxes){
            for(int i = 0; i < razorHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                for(int j = 0; j < razorHistosMC[ibox.first][tree.first].GetNbinsY()+1; j++){
                    double error = 0.0;
                    if(tree.first == "QCD"){
                        error = qcdErrorFrac*razorHistosMC[ibox.first][tree.first].GetBinContent(i, j);
                    }
                    else{
                        error = sysErrorFrac*razorHistosMC[ibox.first][tree.first].GetBinContent(i, j);
                    }
                    razorHistosMC[ibox.first][tree.first].SetBinError(i, j, sqrt(pow(razorHistosMC[ibox.first][tree.first].GetBinError(i, j), 2) + error*error));
                }
            }
            for(int i = 0; i < MRHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                double error = 0.0;
                if(tree.first == "QCD"){
                    error = qcdErrorFrac*MRHistosMC[ibox.first][tree.first].GetBinContent(i);
                }
                else{
                    error = sysErrorFrac*MRHistosMC[ibox.first][tree.first].GetBinContent(i);
                }
                MRHistosMC[ibox.first][tree.first].SetBinError(i, sqrt(pow(MRHistosMC[ibox.first][tree.first].GetBinError(i), 2) + error*error));
            }
            for(int i = 0; i < RsqHistosMC[ibox.first][tree.first].GetNbinsX()+1; i++){
                double error = 0.0;
                if(tree.first == "QCD"){
                    error = qcdErrorFrac*RsqHistosMC[ibox.first][tree.first].GetBinContent(i);
                }
                else{
                    error = sysErrorFrac*RsqHistosMC[ibox.first][tree.first].GetBinContent(i);
                }
                RsqHistosMC[ibox.first][tree.first].SetBinError(i, sqrt(pow(RsqHistosMC[ibox.first][tree.first].GetBinError(i), 2) + error*error));
            }
        }
    }

    //////////////////////////////////////////////////
    //Step 2: make data distributions
    //////////////////////////////////////////////////

    //create histograms
    map<RazorAnalyzer::RazorBox, TH2F> razorData;
    map<RazorAnalyzer::RazorBox, TH1F> MRData;
    map<RazorAnalyzer::RazorBox, TH1F> RsqData;
    for(auto &ibox : boxes){
        razorData[ibox.first] = TH2F(Form("razordata%s", ibox.second.c_str()), "; MR (GeV); Rsq", NMRBINS, MRBINLOWEDGES, NRSQBINS, RSQBINLOWEDGES);
        MRData[ibox.first] = TH1F(Form("mrdata%s", ibox.second.c_str()), "; MR (GeV)", NMRBINS, MRBINLOWEDGES);
        RsqData[ibox.first] = TH1F(Form("rsqdata%s", ibox.second.c_str()), "; Rsq (GeV)", NRSQBINS, RSQBINLOWEDGES);
        razorData[ibox.first].Sumw2();
        MRData[ibox.first].Sumw2();
        RsqData[ibox.first].Sumw2();
    }

    for(auto &tree : datatrees){ 
        cout << "Filling data histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i);
            RazorAnalyzer::RazorBox razorbox = static_cast<RazorAnalyzer::RazorBox>(box);

            //enforce correct box
            if(tree.first != boxDatasets[razorbox]) continue;

            //enforce correct number of B-tags
            if(bTagsInclusive){
                if(nBTaggedJets < minNBTags) continue;
            }
            else {
                if(nBTaggedJets != minNBTags) continue;
            }
            //cut on MR and Rsq
            if(MR < MRBINLOWEDGES[0] || Rsq < RSQBINLOWEDGES[0]) continue;
            //cut on dPhiRazor
            if(doDPhiRazorCut && fabs(dPhiRazor) > dPhiRazorCut) continue;
            //cut on met
            if(doMetCut && met < metCut) continue;
            //cut on mT
            if(doMTCut && (isSingleMuonBox(razorbox) || isSingleElectronBox(razorbox))){
                if(mT < mTLowerCut || mT > mTUpperCut) continue;
            }
            //cut on lepton pt
            if(doLeptonPtCut){
                if(isSingleMuonBox(razorbox) && leadingTightMuPt < leptonPtCut) continue;
                if(isSingleElectronBox(razorbox) && leadingTightElePt < leptonPtCut) continue;
            }

            float eventWeight = 1.0;

            if(boxes.find(razorbox) == boxes.end()){
                cout << "Box " << razorbox << " not in the list!" << endl;
                continue;
            }
            razorData[razorbox].Fill(MR, Rsq, eventWeight);
            if(Rsq > RsqCutFor1DPlots) MRData[razorbox].Fill(MR, eventWeight);
            if(MR > MRCutFor1DPlots) RsqData[razorbox].Fill(Rsq, eventWeight);

            //combined Single Lepton boxes
            //if(isSingleMuonBox(razorbox) || isSingleElectronBox(razorbox)){
            //MultiJet box
            //if(razorbox == RazorAnalyzer::FourJet || razorbox == RazorAnalyzer::SixJet){
            //LooseLeptonMultiJet+MultiJet box
            if(razorbox == RazorAnalyzer::FourJet || razorbox == RazorAnalyzer::SixJet || razorbox == RazorAnalyzer::LooseLeptonFourJet || razorbox == RazorAnalyzer::LooseLeptonSixJet){
                razorData[RazorAnalyzer::NONE].Fill(MR, Rsq, eventWeight);
                if(Rsq > RsqCutFor1DPlots){
                    MRData[RazorAnalyzer::NONE].Fill(MR, eventWeight);
                }
                if(MR > MRCutFor1DPlots){
                    RsqData[RazorAnalyzer::NONE].Fill(Rsq, eventWeight);
                }
            }
        }
    }

    //////////////////////////////////////////////////
    //make plots
    //////////////////////////////////////////////////
    TCanvas c("c", "c", 800, 600);
    for(auto &ibox : boxes){
        gStyle->SetPaintTextFormat("1.0f");
        c.SetLogx();

        //total MC histogram
        TH2F TotalRazorMC("TotalRazorMC", "; MR (GeV); Rsq", NMRBINS, MRBINLOWEDGES, NRSQBINS, RSQBINLOWEDGES);
        //print MC histograms
        c.SetLogz();
        for(auto &hist : razorHistosMC[ibox.first]){
            hist.second.SetTitle(Form("MC for %s, %s Box", hist.first.c_str(), ibox.second.c_str()));
            hist.second.GetXaxis()->SetTitle("MR");
            hist.second.GetYaxis()->SetTitle("Rsq");
            hist.second.SetStats(0);
            hist.second.Draw("colz");
            hist.second.Draw("same,text");
            //c.Print(Form("%s/%s/MCHistogram%s%s.gif", plotPath.c_str(), ibox.second.c_str(), hist.first.c_str(), ibox.second.c_str()));
            //c.Print(Form("%s/%s/MCHistogram%s%s.root", plotPath.c_str(), ibox.second.c_str(), hist.first.c_str(), ibox.second.c_str()));

            //add to total histogram
            TotalRazorMC = TotalRazorMC + hist.second;
        }
        TotalRazorMC.SetTitle(Form("Total MC, %s Box", ibox.second.c_str()));
        TotalRazorMC.SetStats(0);
        TotalRazorMC.Draw("colz");
        TotalRazorMC.Draw("same,text");
        c.Print(Form("%s/%s/MCHistogramTotal%s.gif", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));
        //c.Print(Form("%s/%s/MCHistogramTotal%s.root", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));

        //print data histogram
        razorData[ibox.first].SetTitle(Form("Data, %s Box", ibox.second.c_str()));
        razorData[ibox.first].GetXaxis()->SetTitle("MR");
        razorData[ibox.first].GetYaxis()->SetTitle("Rsq");
        razorData[ibox.first].SetStats(0);
        razorData[ibox.first].Draw("colz");
        razorData[ibox.first].Draw("same,text");
        c.Print(Form("%s/%s/DataHistogram%s.gif", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));
        //c.Print(Form("%s/%s/DataHistogram%s.root", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));

        //print MR and Rsq 1D histograms, comparing data to MC
        c.SetLogy();
        THStack MRTotalRazorMC("MRTotalRazorMC", Form("MR (Rsq > %.2f), %s Box", RsqCutFor1DPlots, ibox.second.c_str()));
        THStack RsqTotalRazorMC("RsqTotalRazorMC", Form("Rsq (MR > %.0f), %s Box", MRCutFor1DPlots, ibox.second.c_str()));

        //format MC histograms
        MRHistosMC[ibox.first]["QCD"].SetFillColor(33);
        MRHistosMC[ibox.first]["ZJetsNuNu"].SetFillColor(kCyan+1);
        MRHistosMC[ibox.first]["WJets"].SetFillColor(kRed+1);
        MRHistosMC[ibox.first]["TTJets"].SetFillColor(kGreen+3);
        MRHistosMC[ibox.first]["DYJets"].SetFillColor(kAzure);
        MRHistosMC[ibox.first]["SingleTop"].SetFillColor(kBlue+3);
        MRHistosMC[ibox.first]["TTV"].SetFillColor(kSpring);
        MRHistosMC[ibox.first]["VV"].SetFillColor(kViolet+2);
        MRHistosMC[ibox.first]["TTTT"].SetFillColor(kRed+4);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["TTTT"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["VV"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["TTV"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["SingleTop"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["DYJets"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["TTJets"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["WJets"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["ZJetsNuNu"]);
        MRTotalRazorMC.Add(&MRHistosMC[ibox.first]["QCD"]);
        MRData[ibox.first].SetMarkerStyle(20);
        MRData[ibox.first].SetMarkerSize(1);
        RsqHistosMC[ibox.first]["QCD"].SetFillColor(33);
        RsqHistosMC[ibox.first]["ZJetsNuNu"].SetFillColor(kCyan+1);
        RsqHistosMC[ibox.first]["WJets"].SetFillColor(kRed+1);
        RsqHistosMC[ibox.first]["TTJets"].SetFillColor(kGreen+3);
        RsqHistosMC[ibox.first]["DYJets"].SetFillColor(kAzure);
        RsqHistosMC[ibox.first]["SingleTop"].SetFillColor(kBlue+3);
        RsqHistosMC[ibox.first]["TTV"].SetFillColor(kSpring);
        RsqHistosMC[ibox.first]["VV"].SetFillColor(kViolet+2);
        RsqHistosMC[ibox.first]["TTTT"].SetFillColor(kRed+4);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["TTTT"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["VV"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["TTV"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["SingleTop"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["DYJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["TTJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["WJets"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["ZJetsNuNu"]);
        RsqTotalRazorMC.Add(&RsqHistosMC[ibox.first]["QCD"]);
        RsqData[ibox.first].SetMarkerStyle(20);
        RsqData[ibox.first].SetMarkerSize(1);

        //create legend
        TLegend *RazorLegend = new TLegend(0.6, 0.6, 0.9, 0.9);
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["WJets"], "WJets MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["DYJets"], "DYJets MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["ZJetsNuNu"], "ZJetsNuNu MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["TTJets"], "TTJets MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["SingleTop"], "Single Top MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["VV"], "VV MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["TTV"], "TTV MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["TTTT"], "TTTT MC");
        RazorLegend->AddEntry(&MRHistosMC[ibox.first]["QCD"], "QCD MC");
        RazorLegend->AddEntry(&MRData[ibox.first], "2012 Data");
        //DrawDataVsMCRatioPlot(&MRData[ibox.first], &MRTotalRazorMC, RazorLegend, "MR (GeV)", plotPath+"/"+ibox.second+"/MR"+ibox.second, true);
        //DrawDataVsMCRatioPlot(&RsqData[ibox.first], &RsqTotalRazorMC, RazorLegend, "Rsq (GeV)", plotPath+"/"+ibox.second+"/Rsq"+ibox.second, true);

        gStyle->SetPaintTextFormat("1.2f");
        c.SetLogy(false);
        c.SetLogz(false);
        c.SetLogx();

        //plot slices of MR and Rsq
        for(int i = 0; i < NRSQBINS; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorData[ibox.first].ProjectionX(Form("ThisRsqSliceData%d%s", i, ibox.second.c_str()), i+1, i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqSliceMC", Form("MR (%.2f < Rsq < %.2f), %s Box", RSQBINLOWEDGES[i], RSQBINLOWEDGES[i+1], ibox.second.c_str()));
            for(auto &hist : razorHistosMC[ibox.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second.ProjectionX(Form("hist%s%d%s", hist.first.c_str(), i, ibox.second.c_str()), i+1, i+1);
                thisHist->SetFillColor(MRHistosMC[ibox.first][hist.first].GetFillColor());
                ThisRsqSliceMCMap[hist.first] = thisHist;
            }
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTTT"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["VV"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTV"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["SingleTop"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["DYJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["WJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["ZJetsNuNu"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["QCD"]);
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", plotPath+"/"+ibox.second+"/MRExclusiveSlice"+to_string(i)+ibox.second, true);
        }
        for(int i = 0; i < NMRBINS; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorData[ibox.first].ProjectionY(Form("ThisMRSliceData%d%s", i, ibox.second.c_str()), i+1, i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (%0.f < MR < %.0f), %s Box", MRBINLOWEDGES[i], MRBINLOWEDGES[i+1], ibox.second.c_str()));
            for(auto &hist : razorHistosMC[ibox.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second.ProjectionY(Form("hist%s%d%s", hist.first.c_str(), i, ibox.second.c_str()), i+1, i+1);
                thisHist->SetFillColor(RsqHistosMC[ibox.first][hist.first].GetFillColor());
                ThisMRSliceMCMap[hist.first] = thisHist;
            }
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTTT"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["VV"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTV"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["SingleTop"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["DYJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["WJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["ZJetsNuNu"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["QCD"]);
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", plotPath+"/"+ibox.second+"/RsqExclusiveSlice"+to_string(i)+ibox.second, true);
        }
        //inclusive slices
        for(int i = 0; i < NRSQBINS; i++){
            map<string, TH1F*> ThisRsqSliceMCMap;    
            TH1F *ThisRsqSliceData = (TH1F*)razorData[ibox.first].ProjectionX(Form("ThisRsqIncSliceData%d%s", i, ibox.second.c_str()), i+1);
            THStack *ThisRsqSliceMC = new THStack("ThisRsqIncSliceMC", Form("MR (Rsq > %.2f), %s Box", RSQBINLOWEDGES[i], ibox.second.c_str()));
            for(auto &hist : razorHistosMC[ibox.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second.ProjectionX(Form("histinc%s%d%s", hist.first.c_str(), i, ibox.second.c_str()), i+1);
                thisHist->SetFillColor(MRHistosMC[ibox.first][hist.first].GetFillColor());
                ThisRsqSliceMCMap[hist.first] = thisHist;
            }
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTTT"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["VV"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTV"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["SingleTop"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["DYJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["TTJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["WJets"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["ZJetsNuNu"]);
            ThisRsqSliceMC->Add(ThisRsqSliceMCMap["QCD"]);
            DrawDataVsMCRatioPlot(ThisRsqSliceData, ThisRsqSliceMC, RazorLegend, "MR (GeV)", plotPath+"/"+ibox.second+"/MRInclusiveSlice"+to_string(i)+ibox.second, true);
        }
        for(int i = 0; i < NMRBINS; i++){
            map<string, TH1F*> ThisMRSliceMCMap;    
            TH1F *ThisMRSliceData = (TH1F*)razorData[ibox.first].ProjectionY(Form("ThisMRSliceIncData%d%s", i, ibox.second.c_str()), i+1);
            THStack *ThisMRSliceMC = new THStack("ThisMRSliceMC", Form("Rsq (MR > %.0f), %s Box", MRBINLOWEDGES[i], ibox.second.c_str()));
            for(auto &hist : razorHistosMC[ibox.first]){
                TH1F *thisHist;
                thisHist = (TH1F*)hist.second.ProjectionY(Form("histinc%s%d%s", hist.first.c_str(), i, ibox.second.c_str()), i+1);
                thisHist->SetFillColor(RsqHistosMC[ibox.first][hist.first].GetFillColor());
                ThisMRSliceMCMap[hist.first] = thisHist;
            }
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTTT"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["VV"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTV"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["SingleTop"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["DYJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["TTJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["WJets"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["ZJetsNuNu"]);
            ThisMRSliceMC->Add(ThisMRSliceMCMap["QCD"]);
            DrawDataVsMCRatioPlot(ThisMRSliceData, ThisMRSliceMC, RazorLegend, "Rsq", plotPath+"/"+ibox.second+"/RsqInclusiveSlice"+to_string(i)+ibox.second, true);
        }

        //data/MC
        TH2F DataOverMCHist = *((TH2F*)razorData[ibox.first].Clone(Form("DataOverMCHist%s", ibox.second.c_str())));
        DataOverMCHist.Divide(&TotalRazorMC);
        DataOverMCHist.SetStats(0);
        DataOverMCHist.SetMinimum(0.1);
        DataOverMCHist.SetMaximum(3.0);
        DataOverMCHist.SetTitle(Form("Data/MC, %s Box", ibox.second.c_str()));
        DataOverMCHist.Draw("colz");
        DataOverMCHist.Draw("same,text");
        c.Print(Form("%s/%s/DataOverMC%s.gif", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));
        //(data-MC)/sigma
        TH2F DataMCSigmaComparison = *((TH2F*)razorData[ibox.first].Clone(Form("DataMCSigmaComparison%s", ibox.second.c_str())));
        for(int i = 0; i < razorData[ibox.first].GetNbinsX()+1; i++){
            for(int j = 0; j < razorData[ibox.first].GetNbinsY()+1; j++){
                double diffError = sqrt(pow(razorData[ibox.first].GetBinError(i, j), 2) + pow(TotalRazorMC.GetBinError(i, j), 2));
                DataMCSigmaComparison.SetBinContent(i, j, (razorData[ibox.first].GetBinContent(i, j) - TotalRazorMC.GetBinContent(i, j))/diffError);
            }
        }
        DataMCSigmaComparison.SetStats(0);
        DataMCSigmaComparison.SetMinimum(-3.0);
        DataMCSigmaComparison.SetMaximum(3.0);
        DataMCSigmaComparison.SetTitle(Form("(Data-MC)/#sigma_{Data-MC}, %s Box", ibox.second.c_str()));
        DataMCSigmaComparison.Draw("colz");
        DataMCSigmaComparison.Draw("same,text");
            c.Print(Form("%s/%s/DataMCSigma%s.gif", plotPath.c_str(), ibox.second.c_str(), ibox.second.c_str()));

        delete RazorLegend;
    }
}

int main(){
    FullControlRegionBasedPrediction();
    return 0;
}

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX){
    TCanvas c("c", "c", 800, 600);
    c.Clear();
    c.cd();
    TPad pad1("pad1","pad1",0,0.4,1,1);
    pad1.SetBottomMargin(0);
    pad1.SetLogy();
    if(logX) pad1.SetLogx();
    pad1.Draw();
    pad1.cd();
    mcStack->Draw("hist");
    mcStack->SetMaximum(max(mcStack->GetMaximum(), dataHist->GetMaximum()));
    mcStack->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    mcStack->GetYaxis()->SetLabelSize(0.03);
    mcStack->GetYaxis()->SetTitleOffset(0.45);
    mcStack->GetYaxis()->SetTitleSize(0.05);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(1);
    dataHist->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    dataHist->Draw("pesame");
    pad1.Modified();
    gPad->Update();
    //make ratio histogram
    TList * histList = (TList*)mcStack->GetHists();
    TIter next(histList);
    TH1 *mcTotal = (TH1*) histList->First()->Clone();
    //mcTotal->Sumw2();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();
    //dataOverMC->Sumw2();
    dataOverMC->Divide(mcTotal);
    dataOverMC->GetXaxis()->SetTitle(xaxisTitle.c_str());
    dataOverMC->GetYaxis()->SetTitle("Data / MC");
    dataOverMC->SetMinimum(0.5);
    dataOverMC->SetMaximum(1.5);
    dataOverMC->GetXaxis()->SetLabelSize(0.1);
    dataOverMC->GetYaxis()->SetLabelSize(0.08);
    dataOverMC->GetYaxis()->SetTitleOffset(0.35);
    dataOverMC->GetXaxis()->SetTitleOffset(1.00);
    dataOverMC->GetYaxis()->SetTitleSize(0.08);
    dataOverMC->GetXaxis()->SetTitleSize(0.08);
    dataOverMC->SetStats(0);
    leg->Draw();
    c.cd();
    TPad pad2("pad2","pad2",0,0.0,1,0.4);
    pad2.SetTopMargin(0);
    pad2.SetTopMargin(0.008);
    pad2.SetBottomMargin(0.25);
    pad2.SetGridy();
    if(logX) pad2.SetLogx();
    pad2.Draw();
    pad2.cd();
    dataOverMC->Draw("pe");
    pad2.Modified();
    gPad->Update();
    c.Print(Form("%sLog.gif", printString.c_str()));
    //c.Print(Form("%s.root", printString.c_str()));
    pad1.SetLogy(kFALSE);
    pad1.Modified();
    gPad->Update();
    c.Print(Form("%sLinear.gif", printString.c_str()));
}
