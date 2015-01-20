#include "RazorAnalyzer.h"

//TODO: implement these!

bool RazorAnalyzer::isVetoMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && fabs(muon_ip3dSignificance[i]) < 4
      && ( ( muonPt[i] > 20 && muon_relIso04DBetaCorr[i] < 0.4 )
	   ||
	   ( muonPt[i] <= 20 && muon_relIso04DBetaCorr[i]*muonPt[i] < 10 ) 
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
      && muon_relIso04DBetaCorr[i] < 0.2
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
      && muon_relIso04DBetaCorr[i] < 0.12     
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
  if ( ( ( muonPt[i] > 20 && muon_relIso04DBetaCorr[i] < 0.4 )
	 ||
	 ( muonPt[i] <= 20 && muon_relIso04DBetaCorr[i]*muonPt[i] < 10 ) 
	 )
       ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passLooseMuonIso(int i){
  bool pass = false;
  if ( muon_relIso04DBetaCorr[i] < 0.2
      ) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passTightMuonIso(int i){
  bool pass = false;
  if ( muon_relIso04DBetaCorr[i] < 0.12     
      ) {
    pass = true;
  }
  return pass;
}   
