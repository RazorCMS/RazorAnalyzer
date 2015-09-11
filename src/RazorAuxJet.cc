#include "RazorAnalyzer.h"

//From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
//Preliminary operating points derived from ttbar events
bool RazorAnalyzer::isCSVL(int i){
    return jetCISV[i] > 0.605;
}

bool RazorAnalyzer::isCSVM(int i){
    return jetCISV[i] > 0.890;
}

bool RazorAnalyzer::isCSVT(int i){
    return jetCISV[i] > 0.970;
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

//compute the smeared jet pt (if option = "up" or "down", will change the smear factor by +/- 1 sigma )
//NOTE: these are Run 1 recommendations and should be replaced as soon as a Run 2 prescription is available.  
double RazorAnalyzer::JetEnergySmearingFactor( double jetPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, TRandom3 *random, string option ) {

  std::vector<float> fJetEta, fJetPtNPU;
  fJetEta.push_back(jetEta);  
  fJetPtNPU.push_back(jetPt); 
  fJetPtNPU.push_back(NPU); 
  double MCJetResolution = JetResolutionCalculator->resolution(fJetEta,fJetPtNPU);
  
  double c = 1;
  if(option == "up"){ //get c plus 1 sigma
      if (fabs(jetEta) < 0.5) c = 1.105;
      else if(fabs(jetEta) < 1.1) c = 1.127;
      else if(fabs(jetEta) < 1.7) c = 1.150;
      else if(fabs(jetEta) < 2.3) c = 1.254;
      else if(fabs(jetEta) < 2.8) c = 1.316;
      else if(fabs(jetEta) < 3.2) c = 1.458;
      else if(fabs(jetEta) < 5.0) c = 1.247;
  }
  else if(option == "down"){ //get c minus 1 sigma
      if (fabs(jetEta) < 0.5) c = 1.053;
      else if(fabs(jetEta) < 1.1) c = 1.071;
      else if(fabs(jetEta) < 1.7) c = 1.092;
      else if(fabs(jetEta) < 2.3) c = 1.162;
      else if(fabs(jetEta) < 2.8) c = 1.192;
      else if(fabs(jetEta) < 3.2) c = 1.332;
      else if(fabs(jetEta) < 5.0) c = 0.865;
  }
  else{ //get nominal value of c
      if (fabs(jetEta) < 0.5) c = 1.079;
      else if(fabs(jetEta) < 1.1) c = 1.099;
      else if(fabs(jetEta) < 1.7) c = 1.121;
      else if(fabs(jetEta) < 2.3) c = 1.208;
      else if(fabs(jetEta) < 2.8) c = 1.254;
      else if(fabs(jetEta) < 3.2) c = 1.395;
      else if(fabs(jetEta) < 5.0) c = 1.056;
  }

  double sigma = sqrt( c*c - 1) * MCJetResolution;

  return fmax( 1.0 + random->Gaus(0, sigma) , 0);

}

//b-tagging scale factors from https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFb-pt_WITHttbar_payload_EPS13.txt
//(if option = "up" or "down", will change the scale factor by +/- 1 sigma)
//NOTE: These are the Run 1 recommended scale factors and should be replaced as soon as Run 2 factors are available.
double RazorAnalyzer::BTagScaleFactor(double jetPt, bool CSVM, string option){
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
