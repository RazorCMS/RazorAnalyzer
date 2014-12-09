#include "RazorAnalyzer.h"

//B-tag working points for 53X Summer13 studies at 8 TeV (with 22Jan2013 ReReco Data)
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
bool RazorAnalyzer::isCSVL(int i){
    return jetCSV[i] > 0.244;
}

bool RazorAnalyzer::isCSVM(int i){
    return jetCSV[i] > 0.679;
}

bool RazorAnalyzer::isCSVT(int i){
    return jetCSV[i] > 0.898;
}

//Jet Energy Corrections
double RazorAnalyzer::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 FactorizedJetCorrector *jetcorrector,  
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
    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";
  
  return cumulativeCorrection;

}
