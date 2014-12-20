#include "RazorAnalyzer.h"
#include "TLorentzVector.h"

//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzer::isGenTau(int index){
  return ( abs(gParticleId[index]) == 15 && 
	   (abs(gParticleMotherId[index]) == 23 ||abs(gParticleMotherId[index]) == 24) 
	   );
};
  
//Finds closes gen tau and returns index to gParticle
int RazorAnalyzer::findClosestGenTau(double eta, double phi) {
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

int RazorAnalyzer::findClosestRecoTau(double eta, double phi) {
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
int RazorAnalyzer::GetTauMatchedID(double eta, double phi){
  int matchedIndex = findClosestGenTau(eta, phi);
  if(matchedIndex < 0)return 0;//No Match -> ID == 0
  int matchedID = 0;
  if (gParticleMotherId[matchedIndex] > 50) {
    matchedID = gParticleMotherId[matchedIndex];
  } else if (abs(gParticleMotherId[matchedIndex]) == 23 || 
	     abs(gParticleMotherId[matchedIndex]) == 24) {
    matchedID = gParticleId[matchedIndex];
  }
  
  return matchedID;
};

//Returns index of the closest parton. If no match is found returns zero.
int RazorAnalyzer::findClosestParton(float eta, float phi){
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
