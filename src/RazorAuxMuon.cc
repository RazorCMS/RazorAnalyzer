#include "RazorAnalyzer.h"

//TODO: implement these!

bool RazorAnalyzer::isVetoMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
      && muon_relIso04DBetaCorr[i] < 0.4
      ) {
    pass = true;
  }
  pass = true;
  return pass;
}

bool RazorAnalyzer::isLooseMuon(int i){
  bool pass = false;
  if (muonIsLoose[i] 
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
      && abs(muon_d0[i]) < 0.2
      && muon_relIso04DBetaCorr[i] < 0.12     
      ) {
    pass = true;
  }
  return pass;
}   
