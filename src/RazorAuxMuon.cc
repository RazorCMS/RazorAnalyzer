#include "RazorAnalyzer.h"

bool RazorAnalyzer::isMuonPOGLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMuonPOGMediumMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsMedium[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMuonPOGTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsTight[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}   

bool RazorAnalyzer::isVetoMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!(
	  ( muonPt[i] > 20 && (muon_chargedMiniIso[i] + muon_photonAndNeutralHadronMiniIso[i] - 0.5*muon_chargedPileupMiniIso[i])/muonPt[i] < 0.2 )
	  ||
	  ( muonPt[i] <= 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) < 10 ) 	  
	  )) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + muon_photonAndNeutralHadronMiniIso[i] - 0.5*muon_chargedPileupMiniIso[i])/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsMedium[i] && fabs(muon_ip3dSignificance[i]) < 4 && fabs(muon_d0[i]) < 0.2)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + muon_photonAndNeutralHadronMiniIso[i] - 0.5*muon_chargedPileupMiniIso[i])/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}   




bool RazorAnalyzer::passHZZMuonPreselection(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_d0[i]) < 0.5
      && fabs(muon_dZ[i]) < 1.0
      ) {
    pass = true;
  }
  return pass;
}


bool RazorAnalyzer::isHZZMuon(int i){
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



bool RazorAnalyzer::matchTagMuonHLTFilters(int i){
  return  matchMuonHLTFilters(i, "SingleMuon");
}


bool RazorAnalyzer::matchMuonHLTFilters(int i, string HLTFilter){
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
  if (HLTFilter == "IsoTkMu20") {
    if ( muon_passHLTFilter[i][9] ) match = true;    
  }
  if (HLTFilter == "IsoMu27") {
    if ( muon_passHLTFilter[i][8] ) match = true;    
  }
  if (HLTFilter == "IsoTkMu27") {
    if ( muon_passHLTFilter[i][12] ) match = true;    
  }
  if (HLTFilter == "Mu50") {
    if ( muon_passHLTFilter[i][14] ) match = true;    
  }
  if (HLTFilter == "Mu55") {
    if ( muon_passHLTFilter[i][15] ) match = true;    
  }
  if (HLTFilter == "Mu45_eta2p1") {
    if ( muon_passHLTFilter[i][16] ) match = true;    
  }
  if (HLTFilter == "Mu50_eta2p1") {
    if ( muon_passHLTFilter[i][17] ) match = true;    
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
