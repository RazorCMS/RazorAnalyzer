//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 11 17:22:56 2014 by ROOT version 5.34/10
// from TTree RazorEvents/selected miniAOD information
// found on file: /afs/cern.ch/user/s/sixie/eos/cms/store/user/sixie/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/crab_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Spring14miniaod-PU20bx25_POSTLS170_V5-v1_RazorNtupleV1p2_10Nov2014_V1/141111_015750/0000/razorNtupleAna_1.root
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
   Float_t         muonE[22];   //[nMuons]
   Float_t         muonPt[22];   //[nMuons]
   Float_t         muonEta[22];   //[nMuons]
   Float_t         muonPhi[22];   //[nMuons]
   Int_t           muonCharge[22];   //[nMuons]
   Bool_t          muonIsLoose[22];   //[nMuons]
   Bool_t          muonIsTight[22];   //[nMuons]
   Float_t         muon_d0[22];   //[nMuons]
   Float_t         muon_dZ[22];   //[nMuons]
   Float_t         muon_ip3d[22];   //[nMuons]
   Float_t         muon_ip3dSignificance[22];   //[nMuons]
   UShort_t        muonType[22];   //[nMuons]
   UInt_t          muonQuality[22];   //[nMuons]
   Float_t         muon_relIso04DBetaCorr[22];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[8];   //[nElectrons]
   Float_t         elePt[8];   //[nElectrons]
   Float_t         eleEta[8];   //[nElectrons]
   Float_t         elePhi[8];   //[nElectrons]
   Float_t         eleCharge[8];   //[nElectrons]
   Float_t         eleEta_SC[8];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[8];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[8];   //[nElectrons]
   Float_t         eleR9[8];   //[nElectrons]
   Float_t         ele_dEta[8];   //[nElectrons]
   Float_t         ele_dPhi[8];   //[nElectrons]
   Float_t         ele_HoverE[8];   //[nElectrons]
   Float_t         ele_d0[8];   //[nElectrons]
   Float_t         ele_dZ[8];   //[nElectrons]
   Float_t         ele_relIsoDBetaCorr[8];   //[nElectrons]
   Int_t           ele_MissHits[8];   //[nElectrons]
   Bool_t          ele_PassConvVeto[8];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[8];   //[nElectrons]
   Float_t         ele_IDMVATrig[8];   //[nElectrons]
   Float_t         ele_IDMVANonTrig[8];   //[nElectrons]
   Float_t         ele_RegressionE[8];   //[nElectrons]
   Float_t         ele_CombineP4[8];   //[nElectrons]
   Int_t           nTaus;
   Float_t         tauE[9];   //[nTaus]
   Float_t         tauPt[9];   //[nTaus]
   Float_t         tauEta[9];   //[nTaus]
   Float_t         tauPhi[9];   //[nTaus]
   Bool_t          tau_IsLoose[9];   //[nTaus]
   Bool_t          tau_isMedium[9];   //[nTaus]
   Bool_t          tau_isTight[9];   //[nTaus]
   Bool_t          tau_passEleVetoLoose[9];   //[nTaus]
   Bool_t          tau_passEleVetoMedium[9];   //[nTaus]
   Bool_t          tau_passEleVetoTight[9];   //[nTaus]
   Bool_t          tau_passMuVetoLoose[9];   //[nTaus]
   Bool_t          tau_passMuVetoMedium[9];   //[nTaus]
   Bool_t          tau_passMuVetoTight[9];   //[nTaus]
   UInt_t          tau_ID[9];   //[nTaus]
   Float_t         tau_combinedIsoDeltaBetaCorr3Hits[9];   //[nTaus]
   Float_t         tau_eleVetoMVA[9];   //[nTaus]
   Int_t           tau_eleVetoCategory[9];   //[nTaus]
   Float_t         tau_muonVetoMVA[9];   //[nTaus]
   Float_t         tau_isoMVAnewDMwLT[9];   //[nTaus]
   Float_t         tau_isoMVAnewDMwoLT[9];   //[nTaus]
   Float_t         tau_leadCandPt[9];   //[nTaus]
   Int_t           tau_leadCandID[9];   //[nTaus]
   Float_t         tau_leadChargedHadrCandPt[9];   //[nTaus]
   Int_t           tau_leadChargedHadrCandID[9];   //[nTaus]
   Int_t           nPhotons;
   Float_t         phoE[5];   //[nPhotons]
   Float_t         phoPt[5];   //[nPhotons]
   Float_t         phoEta[5];   //[nPhotons]
   Float_t         phoPhi[5];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[5];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[5];   //[nPhotons]
   Float_t         phoR9[5];   //[nPhotons]
   Float_t         pho_HoverE[5];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[5];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[5];   //[nPhotons]
   Float_t         pho_sumPhotonEt[5];   //[nPhotons]
   Int_t           pho_isConversion[5];   //[nPhotons]
   Float_t         pho_RegressionE[5];   //[nPhotons]
   Float_t         pho_IDMVA[5];   //[nPhotons]
   Int_t           nJets;
   Float_t         jetE[20];   //[nJets]
   Float_t         jetPt[20];   //[nJets]
   Float_t         jetEta[20];   //[nJets]
   Float_t         jetPhi[20];   //[nJets]
   Float_t         jetCSV[20];   //[nJets]
   Float_t         jetCISV[20];   //[nJets]
   Float_t         jetMass[20];   //[nJets]
   Float_t         jetJetArea[20];   //[nJets]
   Float_t         jetPileupE[20];   //[nJets]
   Float_t         jetPileupId[20];   //[nJets]
   Int_t           nFatJets;
   Float_t         fatJetE[7];   //[nFatJets]
   Float_t         fatJetPt[7];   //[nFatJets]
   Float_t         fatJetEta[7];   //[nFatJets]
   Float_t         fatJetPhi[7];   //[nFatJets]
   Float_t         metPt;
   Float_t         metPhi;
   Float_t         sumMET;
   Float_t         MR;
   Float_t         RSQ;
   Int_t           nGenJets;
   Float_t         genJetE[29];   //[nGenJets]
   Float_t         genJetPt[29];   //[nGenJets]
   Float_t         genJetEta[29];   //[nGenJets]
   Float_t         genJetPhi[29];   //[nGenJets]
   Float_t         genMetPt;
   Float_t         genMetPhi;
   UShort_t        nGenParticle;
   Int_t           gParticleMotherId[52];   //[nGenParticle]
   Int_t           gParticleMotherIndex[52];   //[nGenParticle]
   Int_t           gParticleId[52];   //[nGenParticle]
   Int_t           gParticleStatus[52];   //[nGenParticle]
   Float_t         gParticleE[52];   //[nGenParticle]
   Float_t         gParticlePt[52];   //[nGenParticle]
   Float_t         gParticleEta[52];   //[nGenParticle]
   Float_t         gParticlePhi[52];   //[nGenParticle]

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
   TBranch        *b_muonQuality;   //!
   TBranch        *b_muon_relIso04DBetaCorr;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_eleE;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEta_SC;   //!
   TBranch        *b_eleSigmaIetaIeta;   //!
   TBranch        *b_eleFull5x5SigmaIetaIeta;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_ele_dEta;   //!
   TBranch        *b_ele_dPhi;   //!
   TBranch        *b_ele_HoverE;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_dZ;   //!
   TBranch        *b_ele_relIsoDBetaCorr;   //!
   TBranch        *b_ele_MissHits;   //!
   TBranch        *b_ele_PassConvVeto;   //!
   TBranch        *b_ele_OneOverEminusOneOverP;   //!
   TBranch        *b_ele_IDMVATrig;   //!
   TBranch        *b_ele_IDMVANonTrig;   //!
   TBranch        *b_ele_RegressionE;   //!
   TBranch        *b_ele_CombineP4;   //!
   TBranch        *b_nTaus;   //!
   TBranch        *b_tauE;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tau_IsLoose;   //!
   TBranch        *b_tau_isMedium;   //!
   TBranch        *b_tau_isTight;   //!
   TBranch        *b_tau_passEleVetoLoose;   //!
   TBranch        *b_tau_passEleVetoMedium;   //!
   TBranch        *b_tau_passEleVetoTight;   //!
   TBranch        *b_tau_passMuVetoLoose;   //!
   TBranch        *b_tau_passMuVetoMedium;   //!
   TBranch        *b_tau_passMuVetoTight;   //!
   TBranch        *b_tau_ID;   //!
   TBranch        *b_tau_combinedIsoDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_eleVetoMVA;   //!
   TBranch        *b_tau_eleVetoCategory;   //!
   TBranch        *b_tau_muonVetoMVA;   //!
   TBranch        *b_tau_isoMVAnewDMwLT;   //!
   TBranch        *b_tau_isoMVAnewDMwoLT;   //!
   TBranch        *b_tau_leadCandPt;   //!
   TBranch        *b_tau_leadCandID;   //!
   TBranch        *b_tau_leadChargedHadrCandPt;   //!
   TBranch        *b_tau_leadChargedHadrCandID;   //!
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
   TBranch        *b_gParticleMotherId;   //!
   TBranch        *b_gParticleMotherIndex;   //!
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
   fChain->SetBranchAddress("muonQuality", muonQuality, &b_muonQuality);
   fChain->SetBranchAddress("muon_relIso04DBetaCorr", muon_relIso04DBetaCorr, &b_muon_relIso04DBetaCorr);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("eleE", eleE, &b_eleE);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEta_SC", eleEta_SC, &b_eleEta_SC);
   fChain->SetBranchAddress("eleSigmaIetaIeta", eleSigmaIetaIeta, &b_eleSigmaIetaIeta);
   fChain->SetBranchAddress("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, &b_eleFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("eleR9", eleR9, &b_eleR9);
   fChain->SetBranchAddress("ele_dEta", ele_dEta, &b_ele_dEta);
   fChain->SetBranchAddress("ele_dPhi", ele_dPhi, &b_ele_dPhi);
   fChain->SetBranchAddress("ele_HoverE", ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_d0", ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dZ", ele_dZ, &b_ele_dZ);
   fChain->SetBranchAddress("ele_relIsoDBetaCorr", ele_relIsoDBetaCorr, &b_ele_relIsoDBetaCorr);
   fChain->SetBranchAddress("ele_MissHits", ele_MissHits, &b_ele_MissHits);
   fChain->SetBranchAddress("ele_PassConvVeto", ele_PassConvVeto, &b_ele_PassConvVeto);
   fChain->SetBranchAddress("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, &b_ele_OneOverEminusOneOverP);
   fChain->SetBranchAddress("ele_IDMVATrig", ele_IDMVATrig, &b_ele_IDMVATrig);
   fChain->SetBranchAddress("ele_IDMVANonTrig", ele_IDMVANonTrig, &b_ele_IDMVANonTrig);
   fChain->SetBranchAddress("ele_RegressionE", ele_RegressionE, &b_ele_RegressionE);
   fChain->SetBranchAddress("ele_CombineP4", ele_CombineP4, &b_ele_CombineP4);
   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("tauE", tauE, &b_tauE);
   fChain->SetBranchAddress("tauPt", tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEta", tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tau_IsLoose", tau_IsLoose, &b_tau_IsLoose);
   fChain->SetBranchAddress("tau_isMedium", tau_isMedium, &b_tau_isMedium);
   fChain->SetBranchAddress("tau_isTight", tau_isTight, &b_tau_isTight);
   fChain->SetBranchAddress("tau_passEleVetoLoose", tau_passEleVetoLoose, &b_tau_passEleVetoLoose);
   fChain->SetBranchAddress("tau_passEleVetoMedium", tau_passEleVetoMedium, &b_tau_passEleVetoMedium);
   fChain->SetBranchAddress("tau_passEleVetoTight", tau_passEleVetoTight, &b_tau_passEleVetoTight);
   fChain->SetBranchAddress("tau_passMuVetoLoose", tau_passMuVetoLoose, &b_tau_passMuVetoLoose);
   fChain->SetBranchAddress("tau_passMuVetoMedium", tau_passMuVetoMedium, &b_tau_passMuVetoMedium);
   fChain->SetBranchAddress("tau_passMuVetoTight", tau_passMuVetoTight, &b_tau_passMuVetoTight);
   fChain->SetBranchAddress("tau_ID", tau_ID, &b_tau_ID);
   fChain->SetBranchAddress("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, &b_tau_combinedIsoDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_eleVetoMVA", tau_eleVetoMVA, &b_tau_eleVetoMVA);
   fChain->SetBranchAddress("tau_eleVetoCategory", tau_eleVetoCategory, &b_tau_eleVetoCategory);
   fChain->SetBranchAddress("tau_muonVetoMVA", tau_muonVetoMVA, &b_tau_muonVetoMVA);
   fChain->SetBranchAddress("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, &b_tau_isoMVAnewDMwLT);
   fChain->SetBranchAddress("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, &b_tau_isoMVAnewDMwoLT);
   fChain->SetBranchAddress("tau_leadCandPt", tau_leadCandPt, &b_tau_leadCandPt);
   fChain->SetBranchAddress("tau_leadCandID", tau_leadCandID, &b_tau_leadCandID);
   fChain->SetBranchAddress("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, &b_tau_leadChargedHadrCandPt);
   fChain->SetBranchAddress("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, &b_tau_leadChargedHadrCandID);
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
   fChain->SetBranchAddress("gParticleMotherId", gParticleMotherId, &b_gParticleMotherId);
   fChain->SetBranchAddress("gParticleMotherIndex", gParticleMotherIndex, &b_gParticleMotherIndex);
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
