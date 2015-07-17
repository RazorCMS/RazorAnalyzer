#include "RazorAnalyzer.h"
#include "TLorentzVector.h"



//Finds closes gen electron and returns index to gParticle
int RazorAnalyzer::findClosestGenElectron(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 11) continue; 
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }
  
  return matchedIndex;
};


//Finds closes gen muon and returns index to gParticle
int RazorAnalyzer::findClosestGenMuon(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 13) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }
  
  return matchedIndex;
};

//Finds closest gen jet and returns index to gParticle
int RazorAnalyzer::findClosestGenJet(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenJets; j++){
    //    if (gParticleStatus[j] != 1) continue;
    //if (abs(gParticleId[j]) != 13) continue;

    if ( deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < 0.3
         && deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < minDR
         ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, genJetEta[j], genJetPhi[j]);
    }

  }

  return matchedIndex;
};



//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzer::isGenTau(int index){
  return ( abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2 &&
	   (abs(gParticleMotherId[index]) == 23 ||abs(gParticleMotherId[index]) == 24) 
	   );
};


//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzer::isGenLeptonicTau(int index){
  if (abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2 
      && (abs(gParticleMotherId[index]) == 24 || abs(gParticleMotherId[index]) == 23)
      ) {    

    for(int k = 0; k < nGenParticle; k++){
      if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == index) {
	return true;
      }
    }
  }
  return false;
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
  int matchedID = 0;
  int matchedIndex = findClosestGenTau(eta, phi);

  //find muon if no tau was found
  if (matchedIndex < 0) matchedIndex = findClosestGenMuon(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //find electron if no tau or muon was found
  if (matchedIndex < 0) matchedIndex = findClosestGenElectron(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //if nothing was found
  if(matchedIndex < 0) return 0;//No Match -> ID == 0
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


//Checks if a gen muon is in the given eta and phi direction
bool RazorAnalyzer::matchesGenMuon(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 && 
	 abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24) 
	 ) {
      result = true;
      break;
    }    
  }
  return result;
};

//Checks if a gen electron is in the given eta and phi direction
bool RazorAnalyzer::matchesGenElectron(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 && 
	 abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24) 
	 ) {
      result = true;
      break;
    }    
  }
  return result;
};
