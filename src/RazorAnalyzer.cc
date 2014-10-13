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
    fChain->SetBranchStatus("runNum", 1);
    fChain->SetBranchStatus("lumiNum", 1);
    fChain->SetBranchStatus("eventNum", 1);
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
    fChain->SetBranchStatus("muonIsTight", 1);
    fChain->SetBranchStatus("muon_d0", 1);
    fChain->SetBranchStatus("muon_dZ", 1);
    fChain->SetBranchStatus("muon_ip3d", 1);
    fChain->SetBranchStatus("muon_ip3dSignificance", 1);
    fChain->SetBranchStatus("muonType", 1);
    fChain->SetBranchStatus("muon_sumChargedHadronPt", 1);
    fChain->SetBranchStatus("muon_sumChargedParticlePt", 1);
    fChain->SetBranchStatus("muon_sumNeutralHadronEt", 1);
    fChain->SetBranchStatus("muon_sumPhotonEt", 1);
}

void RazorAnalyzer::EnableElectrons(){
    fChain->SetBranchStatus("nElectrons", 1);
    fChain->SetBranchStatus("eleE", 1);
    fChain->SetBranchStatus("elePt", 1);
    fChain->SetBranchStatus("eleEta", 1);
    fChain->SetBranchStatus("elePhi", 1);
    fChain->SetBranchStatus("eleCharge", 1);
    fChain->SetBranchStatus("eleSigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleFull5x5SigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleR9", 1);
    fChain->SetBranchStatus("ele_dEta", 1);
    fChain->SetBranchStatus("ele_dPhi", 1);
    fChain->SetBranchStatus("ele_HoverE", 1);
    fChain->SetBranchStatus("ele_d0", 1);
    fChain->SetBranchStatus("ele_dZ", 1);
    fChain->SetBranchStatus("ele_sumChargedHadronPt", 1);
    fChain->SetBranchStatus("ele_sumNeutralHadronEt", 1);
    fChain->SetBranchStatus("ele_sumPhotonEt", 1);
    fChain->SetBranchStatus("ele_MissHits", 1);
    fChain->SetBranchStatus("ele_ConvRejec", 1);
    fChain->SetBranchStatus("ele_OneOverEminusOneOverP", 1);
    fChain->SetBranchStatus("ele_RegressionE", 1);
    fChain->SetBranchStatus("ele_CombineP4", 1);
}

void RazorAnalyzer::EnableTaus(){
    fChain->SetBranchStatus("nTaus", 1);
    fChain->SetBranchStatus("tauE", 1);
    fChain->SetBranchStatus("tauPt", 1);
    fChain->SetBranchStatus("tauEta", 1);
    fChain->SetBranchStatus("tauPhi", 1);
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
    fChain->SetBranchStatus("pho_isConversion", 1);
    fChain->SetBranchStatus("pho_RegressionE", 1);
    fChain->SetBranchStatus("pho_IDMVA", 1);
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
    fChain->SetBranchStatus("sumMET", 1);
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
}

void RazorAnalyzer::EnableGenParticles(){
    fChain->SetBranchStatus("nGenParticle", 1);
    fChain->SetBranchStatus("motherIndex", 1);
    fChain->SetBranchStatus("gParticleStatus", 1);
    fChain->SetBranchStatus("gParticleId", 1);
    fChain->SetBranchStatus("gParticleE", 1);
    fChain->SetBranchStatus("gParticlePt", 1);
    fChain->SetBranchStatus("gParticleEta", 1);
    fChain->SetBranchStatus("gParticlePhi", 1);
}
