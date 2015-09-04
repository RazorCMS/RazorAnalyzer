#include "RazorAnalyzer.h"

using namespace std;

RazorAnalyzer::RazorAnalyzer(TTree *tree) : RazorEvents(tree)
{
    //turn off all branches
    fChain->SetBranchStatus("*", 0);
}

RazorAnalyzer::~RazorAnalyzer()
{

}

//NOTE: the functions below need to be maintained by hand.  If variables are added or removed from the ntuple, these functions need to be updated to reflect the changes.

void RazorAnalyzer::EnableAll(){
    EnableEventInfo();
    EnableMuons();
    EnableElectrons();
    EnableTaus();
    EnablePhotons();
    EnableJets();
    EnableFatJets();
    EnableMet();
    EnablePileup();
    EnableMC();
    EnableGenParticles();
    EnableRazor();
}

void RazorAnalyzer::EnableEventInfo(){
    fChain->SetBranchStatus("nPV", 1);
    fChain->SetBranchStatus("pvX", 1);
    fChain->SetBranchStatus("pvY", 1);
    fChain->SetBranchStatus("pvZ", 1);
    fChain->SetBranchStatus("isData", 1);
    fChain->SetBranchStatus("runNum", 1);
    fChain->SetBranchStatus("lumiNum", 1);
    fChain->SetBranchStatus("eventNum", 1);
    fChain->SetBranchStatus("fixedGridRhoAll", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetAll", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetAllCalo", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralCalo", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralChargedPileUp", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralNeutral", 1);
    fChain->SetBranchStatus("HLTDecision", 1);
    fChain->SetBranchStatus("HLTPrescale", 1);
}

void RazorAnalyzer::EnablePileup(){
    fChain->SetBranchStatus("nBunchXing", 1);
    fChain->SetBranchStatus("BunchXing", 1);
    fChain->SetBranchStatus("nPU", 1);
    fChain->SetBranchStatus("nPUmean", 1);
}

void RazorAnalyzer::EnableMuons(){
    fChain->SetBranchStatus("nMuons", 1);
    fChain->SetBranchStatus("muonE", 1);
    fChain->SetBranchStatus("muonPt", 1);
    fChain->SetBranchStatus("muonEta", 1);
    fChain->SetBranchStatus("muonPhi", 1);
    fChain->SetBranchStatus("muonCharge", 1);
    fChain->SetBranchStatus("muonIsLoose", 1);
    fChain->SetBranchStatus("muonIsMedium", 1);
    fChain->SetBranchStatus("muonIsTight", 1);
    fChain->SetBranchStatus("muon_d0", 1);
    fChain->SetBranchStatus("muon_dZ", 1);
    fChain->SetBranchStatus("muon_ip3d", 1);
    fChain->SetBranchStatus("muon_ip3dSignificance", 1);
    fChain->SetBranchStatus("muonType", 1);
    fChain->SetBranchStatus("muonQuality", 1);
    fChain->SetBranchStatus("muon_pileupIso", 1);
    fChain->SetBranchStatus("muon_chargedIso", 1);
    fChain->SetBranchStatus("muon_photonIso", 1);
    fChain->SetBranchStatus("muon_neutralHadIso", 1);
    fChain->SetBranchStatus("muon_ptrel", 1);
    fChain->SetBranchStatus("muon_chargedMiniIso", 1);
    fChain->SetBranchStatus("muon_photonAndNeutralHadronMiniIso", 1);
    fChain->SetBranchStatus("muon_chargedPileupMiniIso", 1);
    fChain->SetBranchStatus("muon_activityMiniIsoAnnulus", 1);
    fChain->SetBranchStatus("muon_passSingleMuTagFilter", 1);
    fChain->SetBranchStatus("muon_passHLTFilter", 1);
}

void RazorAnalyzer::EnableElectrons(){
    fChain->SetBranchStatus("nElectrons", 1);
    fChain->SetBranchStatus("eleE", 1);
    fChain->SetBranchStatus("elePt", 1);
    fChain->SetBranchStatus("eleEta", 1);
    fChain->SetBranchStatus("elePhi", 1);
    fChain->SetBranchStatus("eleCharge", 1);
    fChain->SetBranchStatus("eleEta_SC", 1);
    fChain->SetBranchStatus("eleSigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleFull5x5SigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleR9", 1);
    fChain->SetBranchStatus("ele_dEta", 1);
    fChain->SetBranchStatus("ele_dPhi", 1);
    fChain->SetBranchStatus("ele_HoverE", 1);
    fChain->SetBranchStatus("ele_d0", 1);
    fChain->SetBranchStatus("ele_dZ", 1);
    fChain->SetBranchStatus("ele_ip3d", 1);
    fChain->SetBranchStatus("ele_ip3dSignificance", 1);
    fChain->SetBranchStatus("ele_pileupIso", 1);
    fChain->SetBranchStatus("ele_chargedIso", 1);
    fChain->SetBranchStatus("ele_photonIso", 1);
    fChain->SetBranchStatus("ele_neutralHadIso", 1);
    fChain->SetBranchStatus("ele_MissHits", 1);
    fChain->SetBranchStatus("ele_PassConvVeto", 1);
    fChain->SetBranchStatus("ele_OneOverEminusOneOverP", 1);
    fChain->SetBranchStatus("ele_IDMVATrig", 1);
    fChain->SetBranchStatus("ele_IDMVANonTrig", 1);
    fChain->SetBranchStatus("ele_RegressionE", 1);
    fChain->SetBranchStatus("ele_CombineP4", 1);
    fChain->SetBranchStatus("ele_ptrel", 1);
    fChain->SetBranchStatus("ele_chargedMiniIso", 1);
    fChain->SetBranchStatus("ele_photonAndNeutralHadronMiniIso", 1);
    fChain->SetBranchStatus("ele_chargedPileupMiniIso", 1);
    fChain->SetBranchStatus("ele_activityMiniIsoAnnulus", 1);
    fChain->SetBranchStatus("ele_passSingleEleTagFilter", 1); 
    fChain->SetBranchStatus("ele_passTPOneTagFilter", 1); 
    fChain->SetBranchStatus("ele_passTPTwoTagFilter", 1); 
    fChain->SetBranchStatus("ele_passTPOneProbeFilter", 1); 
    fChain->SetBranchStatus("ele_passTPTwoProbeFilter", 1); 
    fChain->SetBranchStatus("ele_passHLTFilter", 1); 
}

void RazorAnalyzer::EnableTaus(){
    fChain->SetBranchStatus("nTaus", 1);
    fChain->SetBranchStatus("tauE", 1);
    fChain->SetBranchStatus("tauPt", 1);
    fChain->SetBranchStatus("tauEta", 1);
    fChain->SetBranchStatus("tauPhi", 1);
    fChain->SetBranchStatus("tau_IsLoose", 1);
    fChain->SetBranchStatus("tau_IsMedium", 1);
    fChain->SetBranchStatus("tau_IsTight", 1);
    fChain->SetBranchStatus("tau_passEleVetoLoose", 1);
    fChain->SetBranchStatus("tau_passEleVetoMedium", 1);
    fChain->SetBranchStatus("tau_passEleVetoTight", 1);
    fChain->SetBranchStatus("tau_passMuVetoLoose", 1);
    fChain->SetBranchStatus("tau_passMuVetoMedium", 1);
    fChain->SetBranchStatus("tau_passMuVetoTight", 1);
    fChain->SetBranchStatus("tau_ID", 1);
    fChain->SetBranchStatus("tau_combinedIsoDeltaBetaCorr3Hits", 1);
    fChain->SetBranchStatus("tau_eleVetoMVA", 1);
    fChain->SetBranchStatus("tau_eleVetoCategory", 1);
    fChain->SetBranchStatus("tau_muonVetoMVA", 1);
    fChain->SetBranchStatus("tau_isoMVAnewDMwLT", 1);
    fChain->SetBranchStatus("tau_isoMVAnewDMwoLT", 1);
    fChain->SetBranchStatus("tau_leadCandPt", 1);
    fChain->SetBranchStatus("tau_leadCandID", 1);
    fChain->SetBranchStatus("tau_leadChargedHadrCandPt", 1);
    fChain->SetBranchStatus("tau_leadChargedHadrCandID", 1);
}

void RazorAnalyzer::EnableIsoPFCandidates(){
    fChain->SetBranchStatus("nIsoPFCandidates", 1);
    fChain->SetBranchStatus("isoPFCandidatePt", 1);
    fChain->SetBranchStatus("isoPFCandidateEta", 1);
    fChain->SetBranchStatus("isoPFCandidatePhi", 1);
    fChain->SetBranchStatus("isoPFCandidateIso04", 1);
    fChain->SetBranchStatus("isoPFCandidateD0", 1);
    fChain->SetBranchStatus("isoPFCandidatePdgId", 1);  
}

void RazorAnalyzer::EnablePhotons(){
    fChain->SetBranchStatus("nPhotons", 1);
    fChain->SetBranchStatus("phoE", 1);
    fChain->SetBranchStatus("phoPt", 1);
    fChain->SetBranchStatus("phoEta", 1);
    fChain->SetBranchStatus("phoPhi", 1);
    fChain->SetBranchStatus("phoSigmaIetaIeta", 1);
    fChain->SetBranchStatus("phoFull5x5SigmaIetaIeta", 1);
    fChain->SetBranchStatus("phoR9", 1);
    fChain->SetBranchStatus("pho_HoverE", 1);
    fChain->SetBranchStatus("pho_sumChargedHadronPt", 1);
    fChain->SetBranchStatus("pho_sumNeutralHadronEt", 1);
    fChain->SetBranchStatus("pho_sumPhotonEt", 1);
    fChain->SetBranchStatus("pho_sumWorstVertexChargedHadronPt", 1);
    fChain->SetBranchStatus("pho_isConversion", 1);
    fChain->SetBranchStatus("pho_passEleVeto", 1);
    fChain->SetBranchStatus("pho_RegressionE", 1);
    fChain->SetBranchStatus("pho_RegressionEUncertainty", 1);
    fChain->SetBranchStatus("pho_IDMVA", 1);
    fChain->SetBranchStatus("pho_superClusterEta", 1);
    fChain->SetBranchStatus("pho_superClusterPhi", 1);
    fChain->SetBranchStatus("pho_hasPixelSeed", 1);
    fChain->SetBranchStatus("pho_passHLTFilter", 1);
}

void RazorAnalyzer::EnableJets(){
    fChain->SetBranchStatus("nJets", 1);
    fChain->SetBranchStatus("jetE", 1);
    fChain->SetBranchStatus("jetPt", 1);
    fChain->SetBranchStatus("jetEta", 1);
    fChain->SetBranchStatus("jetPhi", 1);
    fChain->SetBranchStatus("jetCSV", 1);
    fChain->SetBranchStatus("jetCISV", 1);
    fChain->SetBranchStatus("jetMass", 1);
    fChain->SetBranchStatus("jetJetArea", 1);
    fChain->SetBranchStatus("jetPileupE", 1);
    fChain->SetBranchStatus("jetPileupId", 1);
    fChain->SetBranchStatus("jetPileupIdFlag", 1);
    fChain->SetBranchStatus("jetPassIDLoose", 1);
    fChain->SetBranchStatus("jetPassIDTight", 1);
    fChain->SetBranchStatus("jetPassMuFrac", 1);
    fChain->SetBranchStatus("jetPassEleFrac", 1);
    fChain->SetBranchStatus("jetPartonFlavor", 1);
    fChain->SetBranchStatus("jetHadronFlavor", 1);
    fChain->SetBranchStatus("jetChargedEMEnergyFraction", 1);
    fChain->SetBranchStatus("jetNeutralEMEnergyFraction", 1);
    fChain->SetBranchStatus("jetChargedHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetNeutralHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetMuonEnergyFraction", 1);
    fChain->SetBranchStatus("jetHOEnergyFraction", 1);
    fChain->SetBranchStatus("jetHFHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetHFEMEnergyFraction", 1);
}

void RazorAnalyzer::EnableFatJets(){
    fChain->SetBranchStatus("nFatJets", 1);
    fChain->SetBranchStatus("fatJetE", 1);
    fChain->SetBranchStatus("fatJetPt", 1);
    fChain->SetBranchStatus("fatJetEta", 1);
    fChain->SetBranchStatus("fatJetPhi", 1);
}

void RazorAnalyzer::EnableMet(){
    fChain->SetBranchStatus("metPt", 1);
    fChain->SetBranchStatus("metPhi", 1);
    fChain->SetBranchStatus("metNoHFPt", 1);
    fChain->SetBranchStatus("metNoHFPhi", 1);
    fChain->SetBranchStatus("metType0Pt", 1);
    fChain->SetBranchStatus("metType0Phi", 1);
    fChain->SetBranchStatus("metType1Pt", 1);
    fChain->SetBranchStatus("metType1Phi", 1);
    fChain->SetBranchStatus("metType0Plus1Pt", 1);
    fChain->SetBranchStatus("metType0Plus1Phi", 1);
    fChain->SetBranchStatus("sumMET", 1);
    fChain->SetBranchStatus("Flag_HBHENoiseFilter", 1);
    fChain->SetBranchStatus("Flag_CSCTightHaloFilter", 1);
    fChain->SetBranchStatus("Flag_hcalLaserEventFilter", 1);
    fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
    fChain->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 1);
    fChain->SetBranchStatus("Flag_goodVertices", 1);
    fChain->SetBranchStatus("Flag_trackingFailureFilter", 1);
    fChain->SetBranchStatus("Flag_eeBadScFilter", 1);
    fChain->SetBranchStatus("Flag_ecalLaserCorrFilter", 1);
    fChain->SetBranchStatus("Flag_trkPOGFilters", 1);
    fChain->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1);
    fChain->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1);
    fChain->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1);
    fChain->SetBranchStatus("Flag_METFilters", 1);    
    fChain->SetBranchStatus("Flag_EcalDeadCellEvent", 1);    
    fChain->SetBranchStatus("Flag_IsNotDeadEcalCluster", 1);    
    fChain->SetBranchStatus("Flag_EcalDeadDR", 1);    
    fChain->SetBranchStatus("Flag_EcalBoundaryDR", 1);    
}

void RazorAnalyzer::EnableRazor(){
    fChain->SetBranchStatus("MR", 1);
    fChain->SetBranchStatus("Rsq", 1);
}

void RazorAnalyzer::EnableMC(){
    fChain->SetBranchStatus("nGenJets", 1);
    fChain->SetBranchStatus("genJetE", 1);
    fChain->SetBranchStatus("genJetPt", 1);
    fChain->SetBranchStatus("genJetEta", 1);
    fChain->SetBranchStatus("genJetPhi", 1);
    fChain->SetBranchStatus("genMetPt", 1);
    fChain->SetBranchStatus("genMetPhi", 1);
    fChain->SetBranchStatus("genVertexZ", 1);
    fChain->SetBranchStatus("genWeight", 1);
    fChain->SetBranchStatus("genSignalProcessID", 1);
    fChain->SetBranchStatus("genQScale", 1);
    fChain->SetBranchStatus("genAlphaQCD", 1);
    fChain->SetBranchStatus("genAlphaQED", 1);
}

void RazorAnalyzer::EnableGenParticles(){
    fChain->SetBranchStatus("nGenParticle", 1);
    fChain->SetBranchStatus("gParticleMotherId", 1);
    fChain->SetBranchStatus("gParticleMotherIndex", 1);
    fChain->SetBranchStatus("gParticleId", 1);
    fChain->SetBranchStatus("gParticleStatus", 1);
    fChain->SetBranchStatus("gParticleE", 1);
    fChain->SetBranchStatus("gParticlePt", 1);
    fChain->SetBranchStatus("gParticleEta", 1);
    fChain->SetBranchStatus("gParticlePhi", 1);
}
