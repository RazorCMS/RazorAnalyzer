#include "RazorAnalyzer.h"

bool RazorAnalyzer::isVetoMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      && ( ( muonPt[i] > 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.4 )
	   ||
	   ( muonPt[i] <= 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) < 10 ) 
	   )
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::isLooseMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.2
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::isTightMuon(int i){
  bool pass = false;
  if (
      muonIsTight[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      && fabs(muon_d0[i]) < 0.2
      && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12     
      ) {
    pass = true;
  }
  return pass;
}   



bool RazorAnalyzer::passVetoMuonID(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passLooseMuonID(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passTightMuonID(int i){
  bool pass = false;
  if (
      muonIsTight[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      && abs(muon_d0[i]) < 0.2
      ) {
    pass = true;
  }
  return pass;
}   

bool RazorAnalyzer::passVetoMuonIso(int i){
  bool pass = false;
  if ( ( ( muonPt[i] > 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.4 )
	 ||
	 ( muonPt[i] <= 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i]))  < 10 ) 
	 )
       ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passLooseMuonIso(int i){
  bool pass = false;
  if ( (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.2
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passTightMuonIso(int i){
  bool pass = false;
  if ( (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12     
      ) {
    pass = true;
  }
  return pass;
}   
