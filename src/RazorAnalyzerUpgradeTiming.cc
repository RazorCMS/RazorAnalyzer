#include "RazorAnalyzerUpgradeTiming.h"
#include "TLorentzVector.h"
#include "Hemisphere.hh"
#include "Davismt2.hh"

using namespace std;

RazorAnalyzerUpgradeTiming::RazorAnalyzerUpgradeTiming(TTree *tree) : RazorEventsUpgradeTiming(tree)
{
    //turn off all branches
    fChain->SetBranchStatus("*", 0);
}

RazorAnalyzerUpgradeTiming::~RazorAnalyzerUpgradeTiming()
{

}

void RazorAnalyzerUpgradeTiming::Analyze(bool isData, bool useTiming, bool usePhoChi2, bool useOddEvent, int option, string outputFileName, string label) {
    cout << "Analyze method called on base RazorAnalyzerUpgradeTiming instance.  Parameters were: " << isData <<"   "<< useTiming << " " <<usePhoChi2<<"   "<< useOddEvent << "  " << option << " " << outputFileName << " " << label << endl;
}

//NOTE: the functions below need to be maintained by hand.  If variables are added or removed from the ntuple, these functions need to be updated to reflect the changes.

void RazorAnalyzerUpgradeTiming::EnableAll(){
    EnableEventInfo();
    EnablePVAll();
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
}

void RazorAnalyzerUpgradeTiming::EnableEventInfo(){
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
//    fChain->SetBranchStatus("HLTDecision", 1);
//    fChain->SetBranchStatus("HLTPrescale", 1);
}

void RazorAnalyzerUpgradeTiming::EnablePVAll() {
//    fChain->SetBranchStatus("beamSpotX",1);
//    fChain->SetBranchStatus("beamSpotY",1);
//    fChain->SetBranchStatus("beamSpotZ",1);
    fChain->SetBranchStatus("nPVAll",1);
    fChain->SetBranchStatus("pvIndex",1);
    fChain->SetBranchStatus("pvAllX",1);
    fChain->SetBranchStatus("pvAllY",1);
    fChain->SetBranchStatus("pvAllZ",1);
    fChain->SetBranchStatus("pvAllT",1);
    fChain->SetBranchStatus("pvAllLogSumPtSq",1);
    fChain->SetBranchStatus("pvNtrack_dt",1);
    fChain->SetBranchStatus("pvNtrack",1);
    fChain->SetBranchStatus("pvAllSumPt",1);
    fChain->SetBranchStatus("pvAllSumPt_dt",1);
    fChain->SetBranchStatus("pvAllSumPx",1);
    fChain->SetBranchStatus("pvAllSumPy",1);
    fChain->SetBranchStatus("pvAllLogSumPtSq_dt",1);
    fChain->SetBranchStatus("pvAllSumPx_dt",1);
    fChain->SetBranchStatus("pvAllSumPy_dt",1);

}

void RazorAnalyzerUpgradeTiming::EnablePileup(){
    fChain->SetBranchStatus("nBunchXing", 1);
    fChain->SetBranchStatus("BunchXing", 1);
    fChain->SetBranchStatus("nPU", 1);
    fChain->SetBranchStatus("nPUmean", 1);
}

void RazorAnalyzerUpgradeTiming::EnableMuons(){
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
    fChain->SetBranchStatus("muon_validFractionTrackerHits", 1);
    fChain->SetBranchStatus("muon_isGlobal", 1);
    fChain->SetBranchStatus("muon_normChi2", 1);
    fChain->SetBranchStatus("muon_chi2LocalPosition", 1);
    fChain->SetBranchStatus("muon_kinkFinder", 1);
    fChain->SetBranchStatus("muon_segmentCompatability", 1);
    fChain->SetBranchStatus("muonIsICHEPMedium", 1);
    fChain->SetBranchStatus("muon_passSingleMuTagFilter", 1);
    fChain->SetBranchStatus("muon_passHLTFilter", 1);
}

void RazorAnalyzerUpgradeTiming::EnableElectrons(){
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

void RazorAnalyzerUpgradeTiming::EnableTaus(){
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
    fChain->SetBranchStatus("tau_chargedIsoPtSum", 1);
    fChain->SetBranchStatus("tau_neutralIsoPtSum", 1);
    fChain->SetBranchStatus("tau_puCorrPtSum", 1);
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

void RazorAnalyzerUpgradeTiming::EnableIsoPFCandidates(){
    fChain->SetBranchStatus("nIsoPFCandidates", 1);
    fChain->SetBranchStatus("isoPFCandidatePt", 1);
    fChain->SetBranchStatus("isoPFCandidateEta", 1);
    fChain->SetBranchStatus("isoPFCandidatePhi", 1);
    fChain->SetBranchStatus("isoPFCandidateIso04", 1);
    fChain->SetBranchStatus("isoPFCandidateD0", 1);
    fChain->SetBranchStatus("isoPFCandidatePdgId", 1);  
}

void RazorAnalyzerUpgradeTiming::EnablePhotons(){
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
    fChain->SetBranchStatus("pho_pfIsoChargedHadronIso", 1);
    fChain->SetBranchStatus("pho_pfIsoChargedHadronIsoWrongVtx", 1);
    fChain->SetBranchStatus("pho_pfIsoNeutralHadronIso", 1);
    fChain->SetBranchStatus("pho_pfIsoPhotonIso", 1);
    fChain->SetBranchStatus("pho_pfIsoModFrixione", 1);
    fChain->SetBranchStatus("pho_pfIsoSumPUPt", 1);    
    fChain->SetBranchStatus("pho_isConversion", 1);
    fChain->SetBranchStatus("pho_passEleVeto", 1);
    fChain->SetBranchStatus("pho_RegressionE", 1);
    fChain->SetBranchStatus("pho_RegressionEUncertainty", 1);
    fChain->SetBranchStatus("pho_IDMVA", 1);
    fChain->SetBranchStatus("pho_superClusterEnergy", 1);
    fChain->SetBranchStatus("pho_superClusterRawEnergy", 1);
    fChain->SetBranchStatus("pho_superClusterEta", 1);
    fChain->SetBranchStatus("pho_superClusterPhi", 1);
    fChain->SetBranchStatus("pho_superClusterSeedX", 1);
    fChain->SetBranchStatus("pho_superClusterSeedY", 1);
    fChain->SetBranchStatus("pho_superClusterSeedZ", 1);
    fChain->SetBranchStatus("pho_superClusterSeedT", 1);
    fChain->SetBranchStatus("pho_superClusterX", 1);
    fChain->SetBranchStatus("pho_superClusterY", 1);
    fChain->SetBranchStatus("pho_superClusterZ", 1);
    fChain->SetBranchStatus("pho_hasPixelSeed", 1);
    fChain->SetBranchStatus("pho_passHLTFilter", 1);
    fChain->SetBranchStatus("pho_convType", 1);
    fChain->SetBranchStatus("pho_convTrkZ", 1);
    fChain->SetBranchStatus("pho_convTrkClusZ", 1);
    fChain->SetBranchStatus("pho_vtxSumPx", 1);
    fChain->SetBranchStatus("pho_vtxSumPy", 1);
}

void RazorAnalyzerUpgradeTiming::EnableJets(){
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
    fChain->SetBranchStatus("jetAllMuonPt", 1);
    fChain->SetBranchStatus("jetAllMuonEta", 1);
    fChain->SetBranchStatus("jetAllMuonPhi", 1);
    fChain->SetBranchStatus("jetAllMuonM", 1);
}

void RazorAnalyzerUpgradeTiming::EnableFatJets(){
    fChain->SetBranchStatus("nFatJets", 1);
    fChain->SetBranchStatus("fatJetE", 1);
    fChain->SetBranchStatus("fatJetPt", 1);
    fChain->SetBranchStatus("fatJetEta", 1);
    fChain->SetBranchStatus("fatJetPhi", 1);
}

void RazorAnalyzerUpgradeTiming::EnableMet(){
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
    fChain->SetBranchStatus("metPuppiPt", 1);
    fChain->SetBranchStatus("metPuppiPhi", 1);
    fChain->SetBranchStatus("metCaloPt", 1);
    fChain->SetBranchStatus("metCaloPhi", 1);
    fChain->SetBranchStatus("sumMET", 1);
    fChain->SetBranchStatus("Flag_HBHENoiseFilter", 1);
    fChain->SetBranchStatus("Flag_HBHETightNoiseFilter", 1);
    fChain->SetBranchStatus("Flag_HBHEIsoNoiseFilter", 1);
    //fChain->SetBranchStatus("Flag_badChargedCandidateFilter", 1);
    //fChain->SetBranchStatus("Flag_badMuonFilter", 1);
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

    fChain->SetBranchStatus("metType1PtJetResUp", 1);
    fChain->SetBranchStatus("metType1PtJetResDown", 1);
    fChain->SetBranchStatus("metType1PtJetEnUp", 1);
    fChain->SetBranchStatus("metType1PtJetEnDown", 1);
    fChain->SetBranchStatus("metType1PtMuonEnUp", 1);
    fChain->SetBranchStatus("metType1PtMuonEnDown", 1);
    fChain->SetBranchStatus("metType1PtElectronEnUp", 1);
    fChain->SetBranchStatus("metType1PtElectronEnDown", 1);
    fChain->SetBranchStatus("metType1PtTauEnUp", 1);
    fChain->SetBranchStatus("metType1PtTauEnDown", 1);
    fChain->SetBranchStatus("metType1PtUnclusteredEnUp", 1);
    fChain->SetBranchStatus("metType1PtUnclusteredEnDown", 1);
    fChain->SetBranchStatus("metType1PtPhotonEnUp", 1);
    fChain->SetBranchStatus("metType1PtPhotonEnDown", 1);
    fChain->SetBranchStatus("metType1PhiJetResUp", 1);
    fChain->SetBranchStatus("metType1PhiJetResDown", 1);
    fChain->SetBranchStatus("metType1PhiJetEnUp", 1);
    fChain->SetBranchStatus("metType1PhiJetEnDown", 1);
    fChain->SetBranchStatus("metType1PhiMuonEnUp", 1);
    fChain->SetBranchStatus("metType1PhiMuonEnDown", 1);
    fChain->SetBranchStatus("metType1PhiElectronEnUp", 1);
    fChain->SetBranchStatus("metType1PhiElectronEnDown", 1);
    fChain->SetBranchStatus("metType1PhiTauEnUp", 1);
    fChain->SetBranchStatus("metType1PhiTauEnDown", 1);
    fChain->SetBranchStatus("metType1PhiUnclusteredEnUp", 1);
    fChain->SetBranchStatus("metType1PhiUnclusteredEnDown", 1);
    fChain->SetBranchStatus("metType1PhiPhotonEnUp", 1);
    fChain->SetBranchStatus("metType1PhiPhotonEnDown", 1);
}

void RazorAnalyzerUpgradeTiming::EnableRazor(){
    fChain->SetBranchStatus("MR", 1);
    fChain->SetBranchStatus("Rsq", 1);
}

void RazorAnalyzerUpgradeTiming::EnableMC(){
    fChain->SetBranchStatus("nGenJets", 1);
    fChain->SetBranchStatus("genJetE", 1);
    fChain->SetBranchStatus("genJetPt", 1);
    fChain->SetBranchStatus("genJetEta", 1);
    fChain->SetBranchStatus("genJetPhi", 1);
    fChain->SetBranchStatus("genMetPt", 1);
    fChain->SetBranchStatus("genMetPhi", 1);
    fChain->SetBranchStatus("genVertexZ", 1);
    fChain->SetBranchStatus("genVertexX", 1);
    fChain->SetBranchStatus("genVertexY", 1);
    fChain->SetBranchStatus("genVertexT", 1);
    fChain->SetBranchStatus("genWeight", 1);
    fChain->SetBranchStatus("genSignalProcessID", 1);
    fChain->SetBranchStatus("genQScale", 1);
    fChain->SetBranchStatus("genAlphaQCD", 1);
    fChain->SetBranchStatus("genAlphaQED", 1);
    //fChain->SetBranchStatus("lheComments", 1);
    fChain->SetBranchStatus("scaleWeights", 1);
    fChain->SetBranchStatus("pdfWeights", 1);
    fChain->SetBranchStatus("allTrackdT", 1);
    fChain->SetBranchStatus("allTrackPt", 1);
    fChain->SetBranchStatus("allTrackdZ", 1);
    fChain->SetBranchStatus("allTrackPvIndex", 1);
    fChain->SetBranchStatus("alphasWeights", 1);
}

void RazorAnalyzerUpgradeTiming::EnableGenParticles(){
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

//////////////////////////////
//ELECTRON
//////////////////////////////

float RazorAnalyzerUpgradeTiming::GetElectronScaleCorrection( double pt, double eta ) {  
  double scaleCorr = 1.0;
  // if ( pt > 0 && fabs(eta) < 1.5) {
  //   scaleCorr = 1.015;
  // } else {
  //   scaleCorr = 1.05;
  // }
  //some dummy code so that it compiler doesn't complain
  if (pt > 0 || eta > 0) {
    scaleCorr = 1.0;
  }
  return scaleCorr;
}

float RazorAnalyzerUpgradeTiming::GetElectronEffectiveAreaMean(int i, bool use25nsCuts ){ 

    double effArea = 0.0;
    //Effective areas below are for the sum of Neutral Hadrons + Photons
    if (use25nsCuts) {
        // These are the Spring15 25ns effective areas reported in this presentation:
        // https://indico.cern.ch/event/369239/contributions/874575/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
        if (fabs(eleEta_SC[i]) < 1.0) {
            effArea = 0.1752;
        } else if (fabs(eleEta_SC[i]) < 1.479) {
            effArea = 0.1862;	
        } else if (fabs(eleEta_SC[i]) < 2.0) {
            effArea = 0.1411;	
        } else if (fabs(eleEta_SC[i]) < 2.2) {
            effArea = 0.1534;	
        } else if (fabs(eleEta_SC[i]) < 2.3) {
            effArea = 0.1903;	
        } else if (fabs(eleEta_SC[i]) < 2.4) {
            effArea = 0.2243;	
        } else if (fabs(eleEta_SC[i]) < 2.5) {
            effArea = 0.2687;	
        }
        return effArea;
    }
    else {
        // These are the Spring15 50ns effective areas reported in this presentation:
        // https://indico.cern.ch/event/369235/contributions/874560/attachments/734635/1007867/Rami_EffAreas.pdf
        if (fabs(eleEta_SC[i]) < 0.8) {
            effArea = 0.0973;
        } else if (fabs(eleEta_SC[i]) < 1.3) {
            effArea = 0.0954;
        } else if (fabs(eleEta_SC[i]) < 2.0) {
            effArea = 0.0632;	
        } else if (fabs(eleEta_SC[i]) < 2.2) {
            effArea = 0.0727;	
        } else {
            effArea = 0.1337;	
        } 
    }
    return effArea;
}

float RazorAnalyzerUpgradeTiming::GetElectronEffectiveArea90(int i){ 

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (fabs(eleEta_SC[i]) < 1.0) {
    effArea = 0.1752;
  } else if (fabs(eleEta_SC[i]) < 1.479) {
    effArea = 0.1862;	
  } else if (fabs(eleEta_SC[i]) < 2.0) {
    effArea = 0.1411;	
  } else if (fabs(eleEta_SC[i]) < 2.2) {
    effArea = 0.1534;	
  } else if (fabs(eleEta_SC[i]) < 2.3) {
    effArea = 0.1903;	
  } else if (fabs(eleEta_SC[i]) < 2.4) {
    effArea = 0.2243;	
  } else if (fabs(eleEta_SC[i]) < 2.5) {
    effArea = 0.2687;	
  }
  return effArea;
}

bool RazorAnalyzerUpgradeTiming::isEGammaPOGVetoElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGVetoElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGVetoElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isEGammaPOGLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGLooseElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isEGammaPOGMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGMediumElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isEGammaPOGTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGTightElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isVetoElectron(int i, bool applyID, bool applyIso){
  return isMVANonTrigVetoElectron(i, applyID, applyIso);
}

bool RazorAnalyzerUpgradeTiming::isLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMVANonTrigVetoElectron(int i, bool applyID, bool applyIso){

  bool pass = true;
  if (applyID) {
    if (!passMVANonTrigVetoElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passMVANonTrigVetoElectronIso(i)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGVetoElectronID(int i, bool use25nsCuts){
    // Veto ID recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;
    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( fabs(ele_dEta[i]) < 0.00749
                && fabs(ele_dPhi[i]) < 0.228
                && eleFull5x5SigmaIetaIeta[i] < 0.0115
                && ele_HoverE[i] < 0.356
                && fabs(ele_d0[i]) < 0.05
                && fabs(ele_dZ[i]) < 0.10
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.299
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 2
           ) {
            pass = true;
        }
    } else {
        if (fabs(ele_dEta[i]) < 0.00895
                && fabs(ele_dPhi[i]) < 0.213
                && eleFull5x5SigmaIetaIeta[i] < 0.037
                && ele_HoverE[i] < 0.211
                && fabs(ele_d0[i]) < 0.1
                && fabs(ele_dZ[i]) < 0.2
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.15
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 3
           ) {
            pass = true;
        }
    } 
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGLooseElectronID(int i, bool use25nsCuts){
    // Loose ID recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;
    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( fabs(ele_dEta[i]) < 0.00477
                && fabs(ele_dPhi[i]) < 0.222
                && eleFull5x5SigmaIetaIeta[i] < 0.011
                && ele_HoverE[i] < 0.298
                && fabs(ele_d0[i]) < 0.05
                && fabs(ele_dZ[i]) < 0.10
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.241
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } else {
        if (fabs(ele_dEta[i]) < 0.00868
                && fabs(ele_dPhi[i]) < 0.213
                && eleFull5x5SigmaIetaIeta[i] < 0.0314
                && ele_HoverE[i] < 0.101
                && fabs(ele_d0[i]) < 0.1
                && fabs(ele_dZ[i]) < 0.2
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.14
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } 
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGMediumElectronID(int i, bool use25nsCuts){
    // Medium ID recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;
    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( fabs(ele_dEta[i]) < 0.00311
                && fabs(ele_dPhi[i]) < 0.103
                && eleFull5x5SigmaIetaIeta[i] < 0.00998
                && ele_HoverE[i] < 0.253
                && fabs(ele_d0[i]) < 0.05
                && fabs(ele_dZ[i]) < 0.10
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.134
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } else {
        if (fabs(ele_dEta[i]) < 0.00609
                && fabs(ele_dPhi[i]) < 0.045
                && eleFull5x5SigmaIetaIeta[i] < 0.0298
                && ele_HoverE[i] < 0.0878
                && fabs(ele_d0[i]) < 0.1
                && fabs(ele_dZ[i]) < 0.2
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.13
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } 
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGTightElectronID(int i, bool use25nsCuts){
    // Tight ID recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;
    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( fabs(ele_dEta[i]) < 0.00308
                && fabs(ele_dPhi[i]) < 0.0816
                && eleFull5x5SigmaIetaIeta[i] < 0.00998
                && ele_HoverE[i] < 0.0414
                && fabs(ele_d0[i]) < 0.05
                && fabs(ele_dZ[i]) < 0.10
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.0129
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } else {
        if (fabs(ele_dEta[i]) < 0.00605
                && fabs(ele_dPhi[i]) < 0.0394
                && eleFull5x5SigmaIetaIeta[i] < 0.0292
                && ele_HoverE[i] < 0.0641
                && fabs(ele_d0[i]) < 0.1
                && fabs(ele_dZ[i]) < 0.2
                && fabs(ele_OneOverEminusOneOverP[i]) < 0.0129
                && ele_PassConvVeto[i]
                && ele_MissHits[i] <= 1
           ) {
            pass = true;
        }
    } 
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passMVANonTrigVetoElectronID(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = -0.11;
  if (subdet == 1 && ptBin == 0) MVACut = -0.55;
  if (subdet == 2 && ptBin == 0) MVACut = -0.60;
  if (subdet == 0 && ptBin == 1) MVACut = -0.16;
  if (subdet == 1 && ptBin == 1) MVACut = -0.65;
  if (subdet == 2 && ptBin == 1) MVACut = -0.74;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut
      && fabs(ele_ip3dSignificance[i]) < 4
      ) {
    pass = true;
  }

  return pass;

}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGVetoElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.175
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.159
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGLooseElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0994
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.107
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGMediumElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0695
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0821
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passEGammaPOGTightElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0588
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0571
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzerUpgradeTiming::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (  ( (elePt[i] > 20 && (ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.2 )
	  ||
	   (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) < 5)
	  )
	) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::passHZZElectronPreselection(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;
  
  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut   
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      ) {
    pass = true;
  }

  return pass;
}


bool RazorAnalyzerUpgradeTiming::isHZZElectron(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;
  
  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut   
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      && fabs(ele_d0[i]) < 0.05
      && ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.4)      
      ) {
    pass = true;
  }

  return pass;
}




bool RazorAnalyzerUpgradeTiming::matchTagElectronHLTFilters(int i){
  bool match = false;
  if ( 
      //Data filters
      ele_passHLTFilter[i][1] 
      || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13] 
      || ele_passHLTFilter[i][49] || ele_passHLTFilter[i][53] || ele_passHLTFilter[i][57] || ele_passHLTFilter[i][60]      
      //MC filters
      || ele_passHLTFilter[i][64] 
      || ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
       ) {
    match = true;
  }

  return match;
}

bool RazorAnalyzerUpgradeTiming::matchProbeElectronHLTFilters(int i){
  bool match = false;
  if ( 
      ele_passHLTFilter[i][50] || ele_passHLTFilter[i][51]     
      || ele_passHLTFilter[i][61]     || ele_passHLTFilter[i][62]     
       ) {
    match = true;
  }
  
  return match;
}

bool RazorAnalyzerUpgradeTiming::matchProbeSCHLTFilters(int i){
  bool match = false;
  if ( 
      ele_passHLTFilter[i][54] || ele_passHLTFilter[i][55] 
      || ele_passHLTFilter[i][58]    || ele_passHLTFilter[i][59]    
       ) {
    match = true;
  }
  
  return match;
}

bool RazorAnalyzerUpgradeTiming::matchElectronHLTFilters(int i, string HLTFilter, string analysisTag) {
  if (analysisTag == "2015") return matchElectronHLTFilters2015(i, HLTFilter);
  else if (analysisTag == "2016") return matchElectronHLTFilters2016(i, HLTFilter);
  else {
    cout << "Analysis Tag " << analysisTag << " is not supported. Returning false.\n";
    return false;
  }
}


bool RazorAnalyzerUpgradeTiming::matchElectronHLTFilters2016(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleElectron") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13] 
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele23Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] 	
	 ) {
      match = true;
    }
  }
  if (HLTFilter == "Ele27Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][5] || ele_passHLTFilter[i][7]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele27Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][6]  || ele_passHLTFilter[i][8]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele32Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][12] 
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele32Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][13] 
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele105") {
    if ( ele_passHLTFilter[i][17] ) match = true;
  }
  if (HLTFilter == "Ele115") {
    if ( ele_passHLTFilter[i][19] ) match = true;
  }
   
  if (HLTFilter == "DoubleElectronLeg1") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][24] || ele_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][25] || ele_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DoubleElectronLeg2DZ") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][26] || ele_passHLTFilter[i][30]
	 ) {
      match = true;
    }
  }
   
  return match;
}

bool RazorAnalyzerUpgradeTiming::matchElectronHLTFilters2015(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleElectron") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13] 
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele23Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] 
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][64]
	 ) {
      match = true;
    }
  }
  if (HLTFilter == "Ele27Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][5] 
	//MC filters
	|| ele_passHLTFilter[i][8]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele27Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][6] 
	//MC filters
	|| ele_passHLTFilter[i][10]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele32Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][13] 
	//MC filters
	|| ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele105") {
    if ( ele_passHLTFilter[i][17] ) match = true;
  }
  if (HLTFilter == "Ele115") {
    if ( ele_passHLTFilter[i][19] ) match = true;
  }
   
  if (HLTFilter == "DoubleElectronLeg1") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][24] || ele_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][25] || ele_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DoubleElectronLeg2DZ") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][26] || ele_passHLTFilter[i][30]
	 ) {
      match = true;
    }
  }
   
  return match;
}

//////////////////////////////
// MUON
//////////////////////////////

float RazorAnalyzerUpgradeTiming::GetMuonEffectiveAreaMean(int i, string type ){ 

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (type == "neutral") {
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0735;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0619;	
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0465;	
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.0433;	
    } else {
      effArea = 0.0577;	
    } 
  } 
  if (type == "charged") {
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0106;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0096;	
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0079;	
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.0058;	
    } else {
      effArea = 0.0053;	
    } 
  } 
  return effArea;
}



bool RazorAnalyzerUpgradeTiming::isMuonPOGLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMuonPOGMediumMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsICHEPMedium[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMuonPOGTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsTight[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}   

bool RazorAnalyzerUpgradeTiming::isVetoMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!(
	  ( muonPt[i] > 20 && (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)) )/muonPt[i] < 0.2 )
	  ||
	  ( muonPt[i] <= 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral") )) < 10 )
	  )) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)))/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzerUpgradeTiming::isTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsICHEPMedium[i] && fabs(muon_ip3dSignificance[i]) < 4 && fabs(muon_d0[i]) < 0.2)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)) )/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}   




bool RazorAnalyzerUpgradeTiming::passHZZMuonPreselection(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_d0[i]) < 0.5
      && fabs(muon_dZ[i]) < 1.0
      ) {
    pass = true;
  }
  return pass;
}


bool RazorAnalyzerUpgradeTiming::isHZZMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_d0[i]) < 0.5
      && fabs(muon_dZ[i]) < 1.0
      && fabs(muon_ip3dSignificance[i]) < 4
      && ( ( (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.4 ))
      ) {
    pass = true;
  }
  return pass;
}



bool RazorAnalyzerUpgradeTiming::matchTagMuonHLTFilters(int i){
  return  matchMuonHLTFilters(i, "SingleMuon");
}


bool RazorAnalyzerUpgradeTiming::matchMuonHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleMuon") {
    if ( 
	//Data filters
	muon_passHLTFilter[i][1]      // IsoMu17_eta2p1
	|| muon_passHLTFilter[i][3]   // IsoMu20
	|| muon_passHLTFilter[i][4]   // HLT_IsoMu20_eta2p1
	|| muon_passHLTFilter[i][6]   // HLT_IsoMu24_eta2p1
	|| muon_passHLTFilter[i][8]   // IsoMu27
	|| muon_passHLTFilter[i][9]   // IsoTkMu20
	|| muon_passHLTFilter[i][10]  // IsoTkMu20_eta2p1 
	|| muon_passHLTFilter[i][11]  // IsoTkMu24_eta2p1
	|| muon_passHLTFilter[i][12]  // IsoTkMu27
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "IsoMu20") {
    if ( muon_passHLTFilter[i][3] ) match = true;    
  }
  if (HLTFilter == "IsoMu22") {
    if ( muon_passHLTFilter[i][7] || muon_passHLTFilter[i][9] ) match = true;    
  }
  if (HLTFilter == "IsoMu24") {
    if ( muon_passHLTFilter[i][11] || muon_passHLTFilter[i][12] ) match = true;    
  }
  if (HLTFilter == "IsoTkMu20") {
    if ( muon_passHLTFilter[i][16] ) match = true;    
  }
  if (HLTFilter == "IsoTkMu22") {
    if ( muon_passHLTFilter[i][17] || muon_passHLTFilter[i][18] ) match = true;    
  }
  if (HLTFilter == "IsoTkMu24") {
    if ( muon_passHLTFilter[i][19] ) match = true;    
  }
  if (HLTFilter == "IsoMu27") {
    if ( muon_passHLTFilter[i][14] ) match = true;    
  }
  if (HLTFilter == "IsoTkMu27") {
    if ( muon_passHLTFilter[i][21] ) match = true;    
  }
  if (HLTFilter == "Mu50") {
    if ( muon_passHLTFilter[i][22] ) match = true;    
  }
  if (HLTFilter == "Mu55") {
    if ( muon_passHLTFilter[i][25] ) match = true;    
  }
  if (HLTFilter == "Mu45_eta2p1") {
    if ( muon_passHLTFilter[i][26] ) match = true;    
  }
  if (HLTFilter == "TkMu50") {
    if ( muon_passHLTFilter[i][23] ) match = true;    
  }
   

  if (HLTFilter == "DoubleMuonLeg1") {
    if ( 
	//Data filters : 17_8 Triggers
	muon_passHLTFilter[i][19] || muon_passHLTFilter[i][27] 
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleMuonLeg2") {
    if ( 
	//Data filters : 17_8 triggers
	muon_passHLTFilter[i][25] || muon_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DoubleMuonLeg2DZ") {
    if ( 
	//Data filters : 17_8 triggers
	muon_passHLTFilter[i][26] || muon_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }
   
  return match;
}

//////////////////////////////
// TAU
//////////////////////////////

bool RazorAnalyzerUpgradeTiming::isLooseTau(int i){
  bool pass = false;
  if (tau_IsLoose[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMediumTau(int i){
  bool pass = false;
  if (tau_IsMedium[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isTightTau(int i){
  bool pass = false;
  if (tau_IsTight[i]) {
    pass = true;
  }

  return pass;
}

//////////////////////////////
// JET 
//////////////////////////////

//From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
bool RazorAnalyzerUpgradeTiming::isCSVL(int i, string dataset){
  if (dataset == "74X") {
    return jetCISV[i] > 0.605;
  } else {
    return jetCISV[i] > 0.460;
  }
}

bool RazorAnalyzerUpgradeTiming::isCSVM(int i, string dataset){
  if (dataset == "74X") {
    return jetCISV[i] > 0.890;
  } else {
    return jetCISV[i] > 0.800;
  }
}

bool RazorAnalyzerUpgradeTiming::isCSVT(int i, string dataset){
  if (dataset == "74X") {
    return jetCISV[i] > 0.970;
  } else {
    return jetCISV[i] > 0.935;
  }
}



//Jet Energy Corrections
double RazorAnalyzerUpgradeTiming::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 FactorizedJetCorrector *jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {
  if (!jetcorrector) {
    cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector->setJetEta(jetEta);
  jetcorrector->setJetPt(jetRawPt);
  jetcorrector->setJetPhi(jetPhi);
  jetcorrector->setJetE(jetE);
  jetcorrector->setRho(rho);
  jetcorrector->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector->getSubCorrections();

  if (printDebug) cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}

//compute the smeared jet pt (if option = "up" or "down", will change the smear factor by +/- 1 sigma )
//NOTE: these are Run 1 recommendations and should be replaced as soon as a Run 2 prescription is available.  
double RazorAnalyzerUpgradeTiming::JetEnergySmearingFactor( double jetPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, TRandom3 *random ) {

  std::vector<float> fJetEta, fJetPtNPU;
  fJetEta.push_back(jetEta);  
  fJetPtNPU.push_back(jetPt); 
  fJetPtNPU.push_back(NPU); 
  double MCJetResolution = JetResolutionCalculator->resolution(fJetEta,fJetPtNPU);
  
  double c = 1;
  if (fabs(jetEta) < 0.5) c = 1.079;
  else if(fabs(jetEta) < 1.1) c = 1.099;
  else if(fabs(jetEta) < 1.7) c = 1.121;
  else if(fabs(jetEta) < 2.3) c = 1.208;
  else if(fabs(jetEta) < 2.8) c = 1.254;
  else if(fabs(jetEta) < 3.2) c = 1.395;
  else if(fabs(jetEta) < 5.0) c = 1.056;

  double sigma = sqrt( c*c - 1) * MCJetResolution;

  return fmax( 1.0 + random->Gaus(0, sigma) , 0);

}

//return smearing factor for up/down shift in JER 
double RazorAnalyzerUpgradeTiming::UpDownJetEnergySmearingFactor(double unsmearedPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, double smearedPt, string option){
    //get jet resolution
    std::vector<float> fJetEta, fJetPtNPU;
    fJetEta.push_back(jetEta);  
    fJetPtNPU.push_back(unsmearedPt); 
    fJetPtNPU.push_back(NPU); 
    double MCJetResolution = JetResolutionCalculator->resolution(fJetEta,fJetPtNPU);

    //get sigma used to smear the jet
    double c = 1;
    if (fabs(jetEta) < 0.5) c = 1.079;
    else if(fabs(jetEta) < 1.1) c = 1.099;
    else if(fabs(jetEta) < 1.7) c = 1.121;
    else if(fabs(jetEta) < 2.3) c = 1.208;
    else if(fabs(jetEta) < 2.8) c = 1.254;
    else if(fabs(jetEta) < 3.2) c = 1.395;
    else if(fabs(jetEta) < 5.0) c = 1.056;
    double sigma = sqrt( c*c - 1) * MCJetResolution;
    //get number of sigmas the jet was smeared
    double z = (smearedPt / unsmearedPt - 1) /sigma;

    if(option == "up"){ //get c plus 1 sigma
        double cUp = 1.0;
        if (fabs(jetEta) < 0.5) cUp = 1.105;
        else if(fabs(jetEta) < 1.1) cUp = 1.127;
        else if(fabs(jetEta) < 1.7) cUp = 1.150;
        else if(fabs(jetEta) < 2.3) cUp = 1.254;
        else if(fabs(jetEta) < 2.8) cUp = 1.316;
        else if(fabs(jetEta) < 3.2) cUp = 1.458;
        else if(fabs(jetEta) < 5.0) cUp = 1.247;
        double sigmaUp = sqrt( cUp*cUp - 1) * MCJetResolution;
        return 1.0 + z*sigmaUp;
    }
    else if(option == "down"){ //get c minus 1 sigma
        double cDown = 1.0;
        if (fabs(jetEta) < 0.5) cDown = 1.053;
        else if(fabs(jetEta) < 1.1) cDown = 1.071;
        else if(fabs(jetEta) < 1.7) cDown = 1.092;
        else if(fabs(jetEta) < 2.3) cDown = 1.162;
        else if(fabs(jetEta) < 2.8) cDown = 1.192;
        else if(fabs(jetEta) < 3.2) cDown = 1.332;
        else if(fabs(jetEta) < 5.0) cDown = 0.865;
        double sigmaDown = sqrt( cDown*cDown - 1) * MCJetResolution;
        return 1.0 + z*sigmaDown;
    }
    else{ 
        std::cout << "Error in UpDownJetEnergySmear: please specify option='up' or 'down'.  Returning 1.0" << std::endl;
    }
    return 1.0;
}


//b-tagging scale factors from https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFb-pt_WITHttbar_payload_EPS13.txt
//(if option = "up" or "down", will change the scale factor by +/- 1 sigma)
//NOTE: These are the Run 1 recommended scale factors and should be replaced as soon as Run 2 factors are available.
double RazorAnalyzerUpgradeTiming::BTagScaleFactor(double jetPt, bool CSVM, string option){
    double tmpBTagCorrFactor = 1.0;
    //nominal correction factor
    double tmpCorrFactor = 0.938887 + 0.00017124 * jetPt + (-2.76366e-07) * jetPt * jetPt ;

    if(option == "up" || option == "down"){
        double uncertainty = 0.0;
        if (jetPt < 30) uncertainty = 0.0415707;
        else if (jetPt < 40) uncertainty = 0.0204209;
        else if (jetPt < 50) uncertainty = 0.0223227;
        else if (jetPt < 60) uncertainty = 0.0206655;
        else if (jetPt < 70) uncertainty = 0.0199325;
        else if (jetPt < 80) uncertainty = 0.0174121;
        else if (jetPt < 100) uncertainty = 0.0202332;
        else if (jetPt < 120) uncertainty = 0.0182446;
        else if (jetPt < 160) uncertainty = 0.0159777;
        else if (jetPt < 210) uncertainty = 0.0218531;
        else if (jetPt < 260) uncertainty = 0.0204688;
        else if (jetPt < 320) uncertainty = 0.0265191;
        else if (jetPt < 400) uncertainty = 0.0313175;
        else if (jetPt < 500) uncertainty = 0.0415417;
        else if (jetPt < 600) uncertainty = 0.0740446;
        else if (jetPt < 800) uncertainty = 0.0596716;
        else uncertainty = 2*0.0596716;

        if (option == "up") tmpCorrFactor += uncertainty;
        else if (option == "down") tmpCorrFactor -= uncertainty;
    }

    double MCEff = 1.0;
    if (jetPt < 50) MCEff = 0.65;
    else if (jetPt < 80) MCEff = 0.70;
    else if (jetPt < 120) MCEff = 0.73;
    else if (jetPt < 210) MCEff = 0.73;
    else MCEff = 0.66;

    //If pass CSV Medium
    if (CSVM) {
        tmpBTagCorrFactor = tmpCorrFactor;
    } else {
        tmpBTagCorrFactor = (1/MCEff - tmpCorrFactor) / (1/MCEff - 1);
    }
    return tmpBTagCorrFactor;
}

//////////////////////////////
// PHOTON
//////////////////////////////

bool RazorAnalyzerUpgradeTiming::photonPassesElectronVeto(int i){
    //use presence of a pixel seed as proxy for an electron veto
    return (pho_passEleVeto[i]);
}


void RazorAnalyzerUpgradeTiming::getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0157;
      effAreaNHad  = 0.0143;
      effAreaPho   = 0.0725;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0143;
      effAreaNHad  = 0.0210;
      effAreaPho   = 0.0604;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0115;
      effAreaNHad  = 0.0148;
      effAreaPho   = 0.0320;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0094;
      effAreaNHad  = 0.0082;
      effAreaPho   = 0.0512;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0095;
      effAreaNHad  = 0.0124;
      effAreaPho   = 0.0766;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0068;
      effAreaNHad  = 0.0186;
      effAreaPho   = 0.0949;
    }
  else
    {
      effAreaChHad = 0.0053;
      effAreaNHad  = 0.0320;
      effAreaPho   = 0.1160;
    }
};

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
void RazorAnalyzerUpgradeTiming::getPhotonEffArea90( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0360;
      effAreaNHad  = 0.0597;
      effAreaPho   = 0.1210;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0377;
      effAreaNHad  = 0.0807;
      effAreaPho   = 0.1107;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0306;
      effAreaNHad  = 0.0629;
      effAreaPho   = 0.0699;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0283;
      effAreaNHad  = 0.0197;
      effAreaPho   = 0.1056;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0254;
      effAreaNHad  = 0.0184;
      effAreaPho   = 0.1457;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0217;
      effAreaNHad  = 0.0284;
      effAreaPho   = 0.1719;
    }
  else
    {
      effAreaChHad = 0.0167;
      effAreaNHad  = 0.0591;
      effAreaPho   = 0.1998;
    }
};



//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzerUpgradeTiming::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    //get the effective areas. results are passed to variables by reference
    if (useEffectiveArea90) 
      {
	getPhotonEffArea90( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }
    else 
      {
	getPhotonEffAreaRun2( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = fmax(pho_pfIsoChargedHadronIso[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    
    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;
    
    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassLooseIDWithoutEleVeto(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01042) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.02683) pass = false;    
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassMediumIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01012) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.02678) pass = false;    
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassTightIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01012) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.02649) pass = false;    
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}


bool RazorAnalyzerUpgradeTiming::photonPassLooseID(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::photonPassMediumID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::photonPassTightID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassLooseIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 1.325, 4.50 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 2.554 + 0.0047*phoPt[i], true );
    } else {
      return photonPassesIsolation(i, 1.293, 4.187 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 3.86 + 0.0034*phoPt[i], true);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassMediumIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 0.789, 2.364 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 0.425 + 0.0047*phoPt[i], true);
    } else {
      return photonPassesIsolation(i, 0.447, 1.765 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 3.15 + 0.0034*phoPt[i], true);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzerUpgradeTiming::photonPassTightIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 0.227, 1.691 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 0.346 + 0.0047*phoPt[i], true);
    } else {
      return photonPassesIsolation(i, 0.146, 0.432 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 2.75 + 0.0034*phoPt[i], true);
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

bool RazorAnalyzerUpgradeTiming::isLoosePhoton(int i, bool use25nsCuts){

  bool pass = true;
  if(!isLoosePhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMediumPhoton(int i, bool use25nsCuts){
  bool pass = true;

  if(!isMediumPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isTightPhoton(int i, bool use25nsCuts){
  bool pass = true;
  if (!isTightPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isLoosePhotonWithoutEleVeto(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isMediumPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzerUpgradeTiming::isTightPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso(i,use25nsCuts)) pass = false;

  return pass;
}


bool RazorAnalyzerUpgradeTiming::matchPhotonHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg1") {
    if ( 
	//Data filters
	pho_passHLTFilter[i][8] 
	//MC filters

	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg2") {
    if ( 
	//Data filters
	pho_passHLTFilter[i][9] 
	 ) {
      match = true;
    }
  }
   
  return match;
}

TLorentzVector RazorAnalyzerUpgradeTiming::GetCorrectedMomentum( TVector3 vtx, TVector3 phoPos, double phoE )
{
  TVector3 phoDir = phoPos - vtx;
  TVector3 phoP3  = phoDir.Unit()*phoE;
  return TLorentzVector( phoP3, phoE);
};

bool RazorAnalyzerUpgradeTiming::photonPassLooseIDWithoutEleVetoExo15004( int i )
{
  bool pass = true;
  if ( fabs(pho_superClusterEta[i]) < 1.4442 )
    {
      if( pho_HoverE[i] > 0.05 ) pass = false;
      if( phoFull5x5SigmaIetaIeta[i] > 0.0105 ) pass = false;
    }
  else
    {
      if( pho_HoverE[i] > 0.05 ) pass = false;
      if( phoFull5x5SigmaIetaIeta[i] > 0.0280 ) pass = false;
    }
  return pass;
};
 
bool RazorAnalyzerUpgradeTiming::photonPassesIsolationExo15004(int i, double PFChHadIsoCut, double PFPhotIsoCut )
{
  double effAreaPhotons = 0.0;
  double eta = pho_superClusterEta[i];
  getPhotonEffAreaExo15004( eta, effAreaPhotons );
  
  //Rho corrected PF charged hadron isolation
  //double PFIsoCorrected_ChHad = pho_sumChargedHadronPt[i];//No PU correction (Caltech Original)
  double PFIsoCorrected_ChHad = pho_pfIsoChargedHadronIso[i];//(Exo15004 default pfIso)
  if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    
  //Rho corrected PF photon isolation
  //double PFIsoCorrected_Photons = pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons;//PU corr can go neg!(Caltech Original)
  double PFIsoCorrected_Photons = pho_pfIsoPhotonIso[i] - fixedGridRhoAll*effAreaPhotons;//PU corr can go neg!(Exo15004 default pfIso)
  if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;
  
  //photon passed all cuts
  return true;
};

void RazorAnalyzerUpgradeTiming::getPhotonEffAreaExo15004( float eta, double& effAreaPho )
{
  if ( fabs( eta ) < 0.9 )
    {
      effAreaPho = 0.17;
    }
  else if ( fabs( eta ) < 1.4442 )
    {
      effAreaPho = 0.14;
    }
  else if ( fabs( eta ) < 2.0 )
    {
      effAreaPho = 0.11;
    }
  else if ( fabs( eta ) < 2.2 )
    {
      effAreaPho = 0.14;
    }
  else
    {
      effAreaPho = 0.22;
    }
};

bool RazorAnalyzerUpgradeTiming::photonPassLooseIsoExo15004(int i)
{
  if( fabs(pho_superClusterEta[i]) < 1.4442 )
    {
      return photonPassesIsolationExo15004(i, 5, (2.75 - 2.5) + 0.0045*phoPt[i] );
    } 
  else if ( fabs(pho_superClusterEta[i]) < 2.0 )
    {
      return photonPassesIsolationExo15004(i, 5, (2.0 - 2.5) + 0.003*phoPt[i] );
    }
  else
    {
      return photonPassesIsolationExo15004(i, 5, (2.0 - 2.5) + 0.003*phoPt[i] );
    }
};

//////////////////////////////
// GEN
//////////////////////////////

//Finds closes gen electron and returns index to gParticle
int RazorAnalyzerUpgradeTiming::findClosestGenElectron(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 11) continue; 
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }
  
  return matchedIndex;
};


//Finds closes gen muon and returns index to gParticle
int RazorAnalyzerUpgradeTiming::findClosestGenMuon(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 13) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }
  
  return matchedIndex;
};

//Finds closest gen jet and returns index to gParticle
int RazorAnalyzerUpgradeTiming::findClosestGenJet(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenJets; j++){
    //    if (gParticleStatus[j] != 1) continue;
    //if (abs(gParticleId[j]) != 13) continue;

    if ( deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < 0.3
         && deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < minDR
         ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, genJetEta[j], genJetPhi[j]);
    }

  }

  return matchedIndex;
};



//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzerUpgradeTiming::isGenTau(int index){
  return ( abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2 &&
	   (abs(gParticleMotherId[index]) == 23 ||abs(gParticleMotherId[index]) == 24) 
	   );
};


//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzerUpgradeTiming::isGenLeptonicTau(int index){
  if (abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2 
      && (abs(gParticleMotherId[index]) == 24 || abs(gParticleMotherId[index]) == 23)
      ) {    

    for(int k = 0; k < nGenParticle; k++){
      if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == index) {
	return true;
      }
    }
  }
  return false;
};
 
//Finds closes gen tau and returns index to gParticle
int RazorAnalyzerUpgradeTiming::findClosestGenTau(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (abs(gParticleId[j]) != 15) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }
  
  return matchedIndex;
};

int RazorAnalyzerUpgradeTiming::findClosestRecoTau(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nTaus; j++){
    if ( deltaR(eta, phi, tauEta[j], tauPhi[j]) < 0.1
	 && deltaR(eta, phi, tauEta[j], tauPhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR(eta, phi, tauEta[j], tauPhi[j]);
    }
  }

  return matchedIndex;
};

//Finds closest gen tau and checks its mother, if matched return pdgID
int RazorAnalyzerUpgradeTiming::GetTauMatchedID(double eta, double phi){
  int matchedID = 0;
  int matchedIndex = findClosestGenTau(eta, phi);

  //find muon if no tau was found
  if (matchedIndex < 0) matchedIndex = findClosestGenMuon(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //find electron if no tau or muon was found
  if (matchedIndex < 0) matchedIndex = findClosestGenElectron(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //if nothing was found
  if(matchedIndex < 0) return 0;//No Match -> ID == 0
  if (gParticleMotherId[matchedIndex] > 50) {
    matchedID = gParticleMotherId[matchedIndex];
  } else if (abs(gParticleMotherId[matchedIndex]) == 23 || 
	     abs(gParticleMotherId[matchedIndex]) == 24) {
    matchedID = gParticleId[matchedIndex];
  }
  
  return matchedID;
};

//Returns index of the closest parton. If no match is found returns zero.
int RazorAnalyzerUpgradeTiming::findClosestParton(float eta, float phi){
  float minDRToParton = 9999;
  int partonIndex = -1;
  for(int j = 0; j < nGenParticle; j++){
    //only look for outgoing partons                                                                    
    if  (!( ( (abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21)
	    && gParticleStatus[j] == 23)
	 ) continue;
    double tmpDR = deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]);
    if ( tmpDR < minDRToParton ) {
      minDRToParton = tmpDR;
      partonIndex = j;
    }
  }
  
  return partonIndex;
};


//Checks if a gen muon is in the given eta and phi direction
bool RazorAnalyzerUpgradeTiming::matchesGenMuon(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 && 
	 abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 15) 
	 ) {
      result = true;
      break;
    }    
  }
  return result;
};

//Checks if a gen electron is in the given eta and phi direction
bool RazorAnalyzerUpgradeTiming::matchesGenElectron(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 && 
	 abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 15) 
	 ) {
      result = true;
      break;
    }    
  }
  return result;
};

//////////////////////////////
// MISC 
//////////////////////////////

double RazorAnalyzerUpgradeTiming::deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzerUpgradeTiming::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}

TLorentzVector RazorAnalyzerUpgradeTiming::makeTLorentzVector(double pt, double eta, double phi, double energy){
    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt, eta, phi, energy);
    return vec;
}

TLorentzVector RazorAnalyzerUpgradeTiming::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass){
    TLorentzVector vec;
    vec.SetPtEtaPhiM(pt, eta, phi, mass);
    return vec;
}

double RazorAnalyzerUpgradeTiming::GetAlphaT(vector<TLorentzVector> jets) 
{   
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;
    double alphaT = 0;
    
    if(nJets < 2) return alphaT;
   
    //int nComb = pow(2, nJets); // # possible combinations
    int nComb = 14;
    if(nJets<14)
	{ 
    	nComb = pow(2, nJets); // # possible combinations
    	}
    // steal from the getHemispheres method

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }
    
    //step 2: Select combination that mininize |ET1 - ET2|
    double eMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;

    for(size_t i=0; i < possibleHem1s.size(); i++)
    {
        double eTemp = fabs(possibleHem1s[i].Et() - possibleHem2s[i].Et());
        if (eMin < 0 || eTemp < eMin)
        {
            eMin = eTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }
    
    float MhtX = 0., MhtY = 0.;
    float HT = 0.; 
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    float MHT = MyMHT.Pt();

    // Calculate alphaT
    alphaT = 0.5 * (1-eMin/HT)/sqrt(1-pow(MHT/HT,2));

    return alphaT;  
}

double RazorAnalyzerUpgradeTiming::GetDPhiMin(vector<TLorentzVector> jets)
    // This variable is used in the alphaT analysis
{
    //int nJets = jets.size();
    double dPhiMin = -1.;
    float HT = 0.;
    float MhtX = 0.;
    float MhtY = 0.;
    // Search for min dPhi between recomputed missing HT and test jets
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }
    TLorentzVector MyMHT;
    MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    for (auto& obj : jets)
    {
    // Recompute MHT by ignoring a test jet 
        float recomputedMHTX = MhtX - obj.Px();
        float recomputedMHTY = MhtY - obj.Py();
        TLorentzVector recomputedMHT;
        recomputedMHT.SetPxPyPzE(-recomputedMHTX, -recomputedMHTY, 0, sqrt(pow(recomputedMHTX,2) + pow(recomputedMHTY,2)));
        double phiTemp = fabs(recomputedMHT.Phi() - obj.Phi());
        if (dPhiMin < 0 || phiTemp < dPhiMin)   dPhiMin = phiTemp;
    }

    return dPhiMin;
}

vector<TLorentzVector> RazorAnalyzerUpgradeTiming::getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided
        TLorentzVector emptyHem1, emptyHem2;
        vector<TLorentzVector> emptyHemsOut;
        emptyHemsOut.push_back(emptyHem1);
        emptyHemsOut.push_back(emptyHem2);
        return emptyHemsOut;
    }

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    //int nComb = pow(2, nJets);
    int nComb = 14;
    if(nJets<14)
	{ 
    	nComb = pow(2, nJets); // # possible combinations
    	}

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }

    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

std::vector< std::vector<int> > RazorAnalyzerUpgradeTiming::getHemispheresV2( std::vector<TLorentzVector> jets )
{
  //returns vector with original indices to jets
  int nJets = jets.size();
  vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations                                                              
  std::vector< std::vector<int> > index1;
  vector<TLorentzVector> possibleHem2s;
  std::vector< std::vector<int> > index2;

  if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided                                                           
    std::vector<int> emptyIndex1, emptyIndex2;
    std::vector< std::vector<int> > void_return;
    void_return.push_back( emptyIndex1 );
    void_return.push_back( emptyIndex2 );
    return void_return;
  }
  
  //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc                                              
  int nComb = pow(2, nJets);
  //std::cout << "njets: " << nJets << " ncomb: " << nComb << std::endl;
  //step 1: store all possible partitions of the input jets                                                                                   
  int j_count;
  for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
    //std::cout << "=iter: " << i << std::endl; 
    TLorentzVector j_temp1, j_temp2;
    std::vector<int> tmp_index1, tmp_index2;
    int itemp = i;
    j_count = nComb/2;
    int count = 0;
    while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2               
      //std::cout << "j_count: " << j_count << " itemp: " << itemp << " count: " << count << std::endl; 
      if(itemp/j_count == 1){
	j_temp1 += jets[count];
	tmp_index1.push_back( count );
      } else {
	j_temp2 += jets[count];
	tmp_index2.push_back( count );
      }
      itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count                                                            
      j_count /= 2;
      count++;
    }
    possibleHem1s.push_back(j_temp1);
    index1.push_back( tmp_index1 );
    possibleHem2s.push_back(j_temp2);
    index2.push_back( tmp_index2 );
  }
  
  //step 2: choose the partition that minimizes m1^2 + m2^2                                                                                   
  double mMin = -1;
  TLorentzVector myHem1;
  TLorentzVector myHem2;
  int partition_index = -1;
  for(size_t i=0; i < possibleHem1s.size(); i++){
    double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
    if(mMin < 0 || mTemp < mMin){
      mMin = mTemp;
      myHem1 = possibleHem1s[i];
      myHem2 = possibleHem2s[i];
      partition_index = i;
    }
  }

  //return the hemispheres in decreasing order of pt                                                                                          
  vector<TLorentzVector> hemsOut;
  std::vector< std::vector<int> > index_out;
  if(myHem1.Pt() > myHem2.Pt()){
    hemsOut.push_back(myHem1);
    hemsOut.push_back(myHem2);
    index_out.push_back( index1[partition_index] );
    index_out.push_back( index2[partition_index] );
  } else {
    hemsOut.push_back(myHem2);
    hemsOut.push_back(myHem1);
    index_out.push_back( index2[partition_index] );
    index_out.push_back( index1[partition_index] );
  }
  
  return index_out;
};


double RazorAnalyzerUpgradeTiming::computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

double RazorAnalyzerUpgradeTiming::computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}


double RazorAnalyzerUpgradeTiming::GetMT( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  //return sqrt( visible.M2() + 2.0*( vis.Pt()*met.Pt() - vis.Dot( met ) ) );
  return sqrt( 2.0*( vis.Pt()*met.Pt() - vis.Dot( met ) ) );
};

double RazorAnalyzerUpgradeTiming::GetMT( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetMT( visible, _met );
};


double RazorAnalyzerUpgradeTiming::GetMTEnergy( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  //return sqrt( visible.M2() + 2.0*( visible.E()*met.Pt() - vis.Dot( met ) ) );
  return sqrt( 2.0*( visible.E()*met.Pt() - vis.Dot( met ) ) );
};

double RazorAnalyzerUpgradeTiming::GetMTEnergy( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetMTEnergy( visible, _met );
};


double RazorAnalyzerUpgradeTiming::GetDphi( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  return vis.DeltaPhi( met );
};

double RazorAnalyzerUpgradeTiming::GetDphi( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetDphi( visible, _met );
};

//auxiliary functions for RazorInclusive and MatchedRazorInclusive analyses
bool RazorAnalyzerUpgradeTiming::passesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool RazorAnalyzerUpgradeTiming::passesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    // if(MR < 300 || Rsq < 0.15) passes = false;
    // if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}

//Checks if ToSubtract matches any particles in Collection, and subtracts the momentum of ToSubtract from the closest one
int RazorAnalyzerUpgradeTiming::SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch){
    //if Collection has any elements within R<deltaRMatch of the vector ToSubtract, find the closest such element
    //otherwise, return -1
    double closestDR = -1;
    int closestDRIndex = -1;
    for(uint i = 0; i < Collection.size(); i++){
        double thisDR = Collection[i].DeltaR(ToSubtract);
        if(closestDR < 0 || thisDR < closestDR){
            closestDR = thisDR;
            closestDRIndex = i;
        }
    }
    if(closestDR < 0 || closestDR > deltaRMatch){ //if we didn't look at any objects or we didn't get one within 0.4, return -1
        return -1;
    }
    
    //subtract the momentum (magnitude) of ToSubtract from that of the closest vector in Collection

    //if we're subtracting away everything, just set the vector to 0 and return the index of the changed vector
    if(ToSubtract.P() >= Collection[closestDRIndex].P()){ 
        Collection[closestDRIndex].SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
        return closestDRIndex;
    }
    //otherwise, scale the 4-momentum by the appropriate factor
    double scalePFactor = (Collection[closestDRIndex].P() - ToSubtract.P())/Collection[closestDRIndex].P(); // = new P / old P
    Collection[closestDRIndex].SetPxPyPzE(Collection[closestDRIndex].Px()*scalePFactor, Collection[closestDRIndex].Py()*scalePFactor, Collection[closestDRIndex].Pz()*scalePFactor, Collection[closestDRIndex].E()*scalePFactor);
    return closestDRIndex;
}

double RazorAnalyzerUpgradeTiming::calcMT2(float testMass, bool massive, vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association)
{
  //computes MT2 using a test mass of testMass, with hemispheres made massless if massive is set to false
  //hemispheres are clustered by finding the grouping of input jets that minimizes the Lund distance
  
  if(jets.size() < 2) return -9999; //need at least two jets for the calculation
  vector<float> px, py, pz, E;
  for(uint i = 0; i < jets.size(); i++){
    //push 4vector components onto individual lists
    px.push_back(jets[i].Px());
    py.push_back(jets[i].Py());
    pz.push_back(jets[i].Pz());
    E.push_back(jets[i].E());
  }
  
  //form the hemispheres using the provided Hemisphere class
  Hemisphere* hemis = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemis->getGrouping();
  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);
        
  //make the hemisphere vectors
  for(uint i=0; i<px.size(); ++i){
    if(grouping[i]==1){
      pseudojet1.SetPx(pseudojet1.Px() + px[i]);
      pseudojet1.SetPy(pseudojet1.Py() + py[i]);
      pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
      pseudojet1.SetE( pseudojet1.E()  + E[i]);   
    }else if(grouping[i] == 2){
      pseudojet2.SetPx(pseudojet2.Px() + px[i]);
      pseudojet2.SetPy(pseudojet2.Py() + py[i]);
      pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
      pseudojet2.SetE( pseudojet2.E()  + E[i]);
    }
  }
  delete hemis;
  
  //now compute MT2 using the Davismt2 class
  
  //these arrays contain (mass, px, py) for the pseudojets and the MET
  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());
  
  pa[0] = static_cast<double> (massive ? pseudojet1.M() : 0);
  pa[1] = static_cast<double> (pseudojet1.Px());
  pa[2] = static_cast<double> (pseudojet1.Py());
  
  pb[0] = static_cast<double> (massive ? pseudojet2.M() : 0);
  pb[1] = static_cast<double> (pseudojet2.Px());
  pb[2] = static_cast<double> (pseudojet2.Py());
  
  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testMass);
  Float_t MT2=mt2->get_mt2();
  delete mt2;
  return MT2;
};

