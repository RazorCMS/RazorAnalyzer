//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 13 12:47:37 2014 by ROOT version 5.32/00
// from TTree RazorEvents/selected miniAOD information
// found on file: razorNtupleAna.root
//////////////////////////////////////////////////////////

#ifndef RazorEvents_h
#define RazorEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class RazorEvents {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nPV;
   Int_t           runNum;
   Int_t           lumiNum;
   Int_t           eventNum;
   Float_t         pvX;
   Float_t         pvY;
   Float_t         pvZ;
   Int_t           nBunchXing;
   Int_t           BunchXing[16];   //[nBunchXing]
   Int_t           nPU[16];   //[nBunchXing]
   Float_t         nPUmean[16];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[7];   //[nMuons]
   Float_t         muonPt[7];   //[nMuons]
   Float_t         muonEta[7];   //[nMuons]
   Float_t         muonPhi[7];   //[nMuons]
   Int_t           muonCharge[7];   //[nMuons]
   Bool_t          muonIsLoose[7];   //[nMuons]
   Bool_t          muonIsTight[7];   //[nMuons]
   Float_t         muon_d0[7];   //[nMuons]
   Float_t         muon_dZ[7];   //[nMuons]
   Float_t         muon_ip3d[7];   //[nMuons]
   Float_t         muon_ip3dSignificance[7];   //[nMuons]
   UShort_t        muonType[7];   //[nMuons]
   Float_t         muon_sumChargedHadronPt[7];   //[nMuons]
   Float_t         muon_sumChargedParticlePt[7];   //[nMuons]
   Float_t         muon_sumNeutralHadronEt[7];   //[nMuons]
   Float_t         muon_sumPhotonEt[7];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[6];   //[nElectrons]
   Float_t         elePt[6];   //[nElectrons]
   Float_t         eleEta[6];   //[nElectrons]
   Float_t         elePhi[6];   //[nElectrons]
   Float_t         eleCharge[6];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[6];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[6];   //[nElectrons]
   Float_t         eleR9[6];   //[nElectrons]
   Float_t         ele_dEta[6];   //[nElectrons]
   Float_t         ele_dPhi[6];   //[nElectrons]
   Float_t         ele_HoverE[6];   //[nElectrons]
   Float_t         ele_d0[6];   //[nElectrons]
   Float_t         ele_dZ[6];   //[nElectrons]
   Float_t         ele_sumChargedHadronPt[6];   //[nElectrons]
   Float_t         ele_sumNeutralHadronEt[6];   //[nElectrons]
   Float_t         ele_sumPhotonEt[6];   //[nElectrons]
   Int_t           ele_MissHits[6];   //[nElectrons]
   Int_t           ele_ConvRejec[6];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[6];   //[nElectrons]
   Float_t         ele_RegressionE[6];   //[nElectrons]
   Float_t         ele_CombineP4[6];   //[nElectrons]
   Int_t           nTaus;
   Float_t         tauE[6];   //[nTaus]
   Float_t         tauPt[6];   //[nTaus]
   Float_t         tauEta[6];   //[nTaus]
   Float_t         tauPhi[6];   //[nTaus]
   Int_t           nPhotons;
   Float_t         phoE[3];   //[nPhotons]
   Float_t         phoPt[3];   //[nPhotons]
   Float_t         phoEta[3];   //[nPhotons]
   Float_t         phoPhi[3];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[3];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[3];   //[nPhotons]
   Float_t         phoR9[3];   //[nPhotons]
   Float_t         pho_HoverE[3];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[3];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[3];   //[nPhotons]
   Float_t         pho_sumPhotonEt[3];   //[nPhotons]
   Int_t           pho_isConversion[3];   //[nPhotons]
   Float_t         pho_RegressionE[3];   //[nPhotons]
   Float_t         pho_IDMVA[3];   //[nPhotons]
   Int_t           nJets;
   Float_t         jetE[17];   //[nJets]
   Float_t         jetPt[17];   //[nJets]
   Float_t         jetEta[17];   //[nJets]
   Float_t         jetPhi[17];   //[nJets]
   Float_t         jetCSV[17];   //[nJets]
   Float_t         jetCISV[17];   //[nJets]
   Float_t         jetMass[17];   //[nJets]
   Float_t         jetJetArea[17];   //[nJets]
   Float_t         jetPileupE[17];   //[nJets]
   Float_t         jetPileupId[17];   //[nJets]
   Int_t           nFatJets;
   Float_t         fatJetE[5];   //[nFatJets]
   Float_t         fatJetPt[5];   //[nFatJets]
   Float_t         fatJetEta[5];   //[nFatJets]
   Float_t         fatJetPhi[5];   //[nFatJets]
   Float_t         metPt;
   Float_t         metPhi;
   Float_t         sumMET;
   Float_t         MR;
   Float_t         RSQ;
   Int_t           nGenJets;
   Float_t         genJetE[21];   //[nGenJets]
   Float_t         genJetPt[21];   //[nGenJets]
   Float_t         genJetEta[21];   //[nGenJets]
   Float_t         genJetPhi[21];   //[nGenJets]
   Float_t         genMetPt;
   Float_t         genMetPhi;
   UShort_t        nGenParticle;
   Int_t           motherIndex[158];   //[nGenParticle]
   Int_t           gParticleId[158];   //[nGenParticle]
   Int_t           gParticleStatus[158];   //[nGenParticle]
   Float_t         gParticleE[158];   //[nGenParticle]
   Float_t         gParticlePt[158];   //[nGenParticle]
   Float_t         gParticleEta[158];   //[nGenParticle]
   Float_t         gParticlePhi[158];   //[nGenParticle]

   // List of branches
   TBranch        *b_nPV;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_lumiNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_pvX;   //!
   TBranch        *b_pvY;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_nBunchXing;   //!
   TBranch        *b_BunchXing;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUmean;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_muonE;   //!
   TBranch        *b_muonPt;   //!
   TBranch        *b_muonEta;   //!
   TBranch        *b_muonPhi;   //!
   TBranch        *b_muonCharge;   //!
   TBranch        *b_muonIsLoose;   //!
   TBranch        *b_muonIsTight;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_dZ;   //!
   TBranch        *b_muon_ip3d;   //!
   TBranch        *b_muon_ip3dSignificance;   //!
   TBranch        *b_muonType;   //!
   TBranch        *b_muon_sumChargedHadronPt;   //!
   TBranch        *b_muon_sumChargedParticlePt;   //!
   TBranch        *b_muon_sumNeutralHadronEt;   //!
   TBranch        *b_muon_sumPhotonEt;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_eleE;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleSigmaIetaIeta;   //!
   TBranch        *b_eleFull5x5SigmaIetaIeta;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_ele_dEta;   //!
   TBranch        *b_ele_dPhi;   //!
   TBranch        *b_ele_HoverE;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_dZ;   //!
   TBranch        *b_ele_sumChargedHadronPt;   //!
   TBranch        *b_ele_sumNeutralHadronEt;   //!
   TBranch        *b_ele_sumPhotonEt;   //!
   TBranch        *b_ele_MissHits;   //!
   TBranch        *b_ele_ConvRejec;   //!
   TBranch        *b_ele_OneOverEminusOneOverP;   //!
   TBranch        *b_ele_RegressionE;   //!
   TBranch        *b_ele_CombineP4;   //!
   TBranch        *b_nTaus;   //!
   TBranch        *b_tauE;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSigmaIetaIeta;   //!
   TBranch        *b_phoFull5x5SigmaIetaIeta;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_pho_HoverE;   //!
   TBranch        *b_pho_sumChargedHadronPt;   //!
   TBranch        *b_pho_sumNeutralHadronEt;   //!
   TBranch        *b_pho_sumPhotonEt;   //!
   TBranch        *b_pho_isConversion;   //!
   TBranch        *b_pho_RegressionE;   //!
   TBranch        *b_pho_IDMVA;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCSV;   //!
   TBranch        *b_jetCISV;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetJetArea;   //!
   TBranch        *b_jetPileupE;   //!
   TBranch        *b_jetPileupId;   //!
   TBranch        *b_nFatJets;   //!
   TBranch        *b_fatJetE;   //!
   TBranch        *b_fatJetPt;   //!
   TBranch        *b_fatJetEta;   //!
   TBranch        *b_fatJetPhi;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_sumMET;   //!
   TBranch        *b_MR;   //!
   TBranch        *b_RSQ;   //!
   TBranch        *b_nGenJets;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genMetPt;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_nGenParticle;   //!
   TBranch        *b_motherIndex;   //!
   TBranch        *b_gParticleId;   //!
   TBranch        *b_gParticleStatus;   //!
   TBranch        *b_gParticleE;   //!
   TBranch        *b_gParticlePt;   //!
   TBranch        *b_gParticleEta;   //!
   TBranch        *b_gParticlePhi;   //!

   RazorEvents(TTree *tree=0);
   virtual ~RazorEvents();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RazorEvents_cxx
RazorEvents::RazorEvents(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("razorNtupleAna.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("razorNtupleAna.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("razorNtupleAna.root:/ntuples");
      dir->GetObject("RazorEvents",tree);

   }
   Init(tree);
}

RazorEvents::~RazorEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RazorEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RazorEvents::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RazorEvents::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("lumiNum", &lumiNum, &b_lumiNum);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("pvX", &pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", &pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("nBunchXing", &nBunchXing, &b_nBunchXing);
   fChain->SetBranchAddress("BunchXing", BunchXing, &b_BunchXing);
   fChain->SetBranchAddress("nPU", nPU, &b_nPU);
   fChain->SetBranchAddress("nPUmean", nPUmean, &b_nPUmean);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("muonE", muonE, &b_muonE);
   fChain->SetBranchAddress("muonPt", muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonEta", muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonCharge", muonCharge, &b_muonCharge);
   fChain->SetBranchAddress("muonIsLoose", muonIsLoose, &b_muonIsLoose);
   fChain->SetBranchAddress("muonIsTight", muonIsTight, &b_muonIsTight);
   fChain->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_dZ", muon_dZ, &b_muon_dZ);
   fChain->SetBranchAddress("muon_ip3d", muon_ip3d, &b_muon_ip3d);
   fChain->SetBranchAddress("muon_ip3dSignificance", muon_ip3dSignificance, &b_muon_ip3dSignificance);
   fChain->SetBranchAddress("muonType", muonType, &b_muonType);
   fChain->SetBranchAddress("muon_sumChargedHadronPt", muon_sumChargedHadronPt, &b_muon_sumChargedHadronPt);
   fChain->SetBranchAddress("muon_sumChargedParticlePt", muon_sumChargedParticlePt, &b_muon_sumChargedParticlePt);
   fChain->SetBranchAddress("muon_sumNeutralHadronEt", muon_sumNeutralHadronEt, &b_muon_sumNeutralHadronEt);
   fChain->SetBranchAddress("muon_sumPhotonEt", muon_sumPhotonEt, &b_muon_sumPhotonEt);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("eleE", eleE, &b_eleE);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleSigmaIetaIeta", eleSigmaIetaIeta, &b_eleSigmaIetaIeta);
   fChain->SetBranchAddress("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, &b_eleFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("eleR9", eleR9, &b_eleR9);
   fChain->SetBranchAddress("ele_dEta", ele_dEta, &b_ele_dEta);
   fChain->SetBranchAddress("ele_dPhi", ele_dPhi, &b_ele_dPhi);
   fChain->SetBranchAddress("ele_HoverE", ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_d0", ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dZ", ele_dZ, &b_ele_dZ);
   fChain->SetBranchAddress("ele_sumChargedHadronPt", ele_sumChargedHadronPt, &b_ele_sumChargedHadronPt);
   fChain->SetBranchAddress("ele_sumNeutralHadronEt", ele_sumNeutralHadronEt, &b_ele_sumNeutralHadronEt);
   fChain->SetBranchAddress("ele_sumPhotonEt", ele_sumPhotonEt, &b_ele_sumPhotonEt);
   fChain->SetBranchAddress("ele_MissHits", ele_MissHits, &b_ele_MissHits);
   fChain->SetBranchAddress("ele_ConvRejec", ele_ConvRejec, &b_ele_ConvRejec);
   fChain->SetBranchAddress("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, &b_ele_OneOverEminusOneOverP);
   fChain->SetBranchAddress("ele_RegressionE", ele_RegressionE, &b_ele_RegressionE);
   fChain->SetBranchAddress("ele_CombineP4", ele_CombineP4, &b_ele_CombineP4);
   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("tauE", tauE, &b_tauE);
   fChain->SetBranchAddress("tauPt", tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEta", tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoPt", phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSigmaIetaIeta", phoSigmaIetaIeta, &b_phoSigmaIetaIeta);
   fChain->SetBranchAddress("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, &b_phoFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("pho_HoverE", pho_HoverE, &b_pho_HoverE);
   fChain->SetBranchAddress("pho_sumChargedHadronPt", pho_sumChargedHadronPt, &b_pho_sumChargedHadronPt);
   fChain->SetBranchAddress("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, &b_pho_sumNeutralHadronEt);
   fChain->SetBranchAddress("pho_sumPhotonEt", pho_sumPhotonEt, &b_pho_sumPhotonEt);
   fChain->SetBranchAddress("pho_isConversion", pho_isConversion, &b_pho_isConversion);
   fChain->SetBranchAddress("pho_RegressionE", pho_RegressionE, &b_pho_RegressionE);
   fChain->SetBranchAddress("pho_IDMVA", pho_IDMVA, &b_pho_IDMVA);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCSV", jetCSV, &b_jetCSV);
   fChain->SetBranchAddress("jetCISV", jetCISV, &b_jetCISV);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetJetArea", jetJetArea, &b_jetJetArea);
   fChain->SetBranchAddress("jetPileupE", jetPileupE, &b_jetPileupE);
   fChain->SetBranchAddress("jetPileupId", jetPileupId, &b_jetPileupId);
   fChain->SetBranchAddress("nFatJets", &nFatJets, &b_nFatJets);
   fChain->SetBranchAddress("fatJetE", fatJetE, &b_fatJetE);
   fChain->SetBranchAddress("fatJetPt", fatJetPt, &b_fatJetPt);
   fChain->SetBranchAddress("fatJetEta", fatJetEta, &b_fatJetEta);
   fChain->SetBranchAddress("fatJetPhi", fatJetPhi, &b_fatJetPhi);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("sumMET", &sumMET, &b_sumMET);
   fChain->SetBranchAddress("MR", &MR, &b_MR);
   fChain->SetBranchAddress("RSQ", &RSQ, &b_RSQ);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("genJetE", genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genMetPt", &genMetPt, &b_genMetPt);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_nGenParticle);
   fChain->SetBranchAddress("motherIndex", motherIndex, &b_motherIndex);
   fChain->SetBranchAddress("gParticleId", gParticleId, &b_gParticleId);
   fChain->SetBranchAddress("gParticleStatus", gParticleStatus, &b_gParticleStatus);
   fChain->SetBranchAddress("gParticleE", gParticleE, &b_gParticleE);
   fChain->SetBranchAddress("gParticlePt", gParticlePt, &b_gParticlePt);
   fChain->SetBranchAddress("gParticleEta", gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", gParticlePhi, &b_gParticlePhi);
   Notify();
}

Bool_t RazorEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RazorEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RazorEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RazorEvents_cxx
