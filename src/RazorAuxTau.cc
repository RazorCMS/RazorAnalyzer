#include "RazorAnalyzer.h"

bool RazorAnalyzer::isLooseTau(int i){
  bool pass = false;
  if (tau_IsLoose[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzer::isMediumTau(int i){
  bool pass = false;
  if (tau_IsMedium[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzer::isTightTau(int i){
  bool pass = false;
  if (tau_IsTight[i]) {
    pass = true;
  }

  return pass;
}
