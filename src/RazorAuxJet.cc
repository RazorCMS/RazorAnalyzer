#include "RazorAnalyzer.h"

//B-tag working points for 53X Summer13 studies at 8 TeV (with 22Jan2013 ReReco Data)
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
bool RazorAnalyzer::isOldCSVL(int i){
    return jetCSV[i] > 0.244;
}

bool RazorAnalyzer::isOldCSVM(int i){
    return jetCSV[i] > 0.679;
}

bool RazorAnalyzer::isOldCSVT(int i){
    return jetCSV[i] > 0.898;
}

//From https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging
//Preliminary operating points derived from ttbar events
bool RazorAnalyzer::isCSVL(int i){
    return jetCISV[i] > 0.423;
}

bool RazorAnalyzer::isCSVM(int i){
    return jetCISV[i] > 0.814;
}

bool RazorAnalyzer::isCSVT(int i){
    return jetCISV[i] > 0.941;
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

double RazorAnalyzer::JetEnergySmearingFactor( double jetPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, TRandom3 *random) {

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
