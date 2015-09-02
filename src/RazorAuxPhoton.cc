//LOCAL INCLUDES
#include "RazorAnalyzer.h"


bool RazorAnalyzer::photonPassesElectronVeto(int i){
    //use presence of a pixel seed as proxy for an electron veto
    return (pho_passEleVeto[i]);
}


void RazorAnalyzer::getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0157;
      effAreaNHad  = 0.0143;
      effAreaPho   = 0.0725;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0143;
      effAreaNHad  = 0.0210;
      effAreaPho   = 0.0604;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0115;
      effAreaNHad  = 0.0148;
      effAreaPho   = 0.0320;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0094;
      effAreaNHad  = 0.0082;
      effAreaPho   = 0.0512;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0095;
      effAreaNHad  = 0.0124;
      effAreaPho   = 0.0766;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0068;
      effAreaNHad  = 0.0186;
      effAreaPho   = 0.0949;
    }
  else
    {
      effAreaChHad = 0.0053;
      effAreaNHad  = 0.0320;
      effAreaPho   = 0.1160;
    }
};

void RazorAnalyzer::getPhotonEffArea95( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.141;
      effAreaNHad  = 0.075;
      effAreaPho   = 0.052;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.117;
      effAreaNHad  = 0.092;
      effAreaPho   = 0.052;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.087;
      effAreaNHad  = 0.069;
      effAreaPho   = 0.041;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.135;
      effAreaNHad  = 0.032;
      effAreaPho   = 0.038;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.186;
      effAreaNHad  = 0.038;
      effAreaPho   = 0.033;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.215;
      effAreaNHad  = 0.051;
      effAreaPho   = 0.029;
    }
  else
    {
      effAreaChHad = 0.250;
      effAreaNHad  = 0.087;
      effAreaPho   = 0.025;
    }
};



//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea95){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    //get the effective areas. results are passed to variables by reference
    if (useEffectiveArea95) {
      getPhotonEffArea95( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
    } else {
      getPhotonEffAreaRun2( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
    }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    
    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;
    
    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = max(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::photonPassLooseIDWithoutEleVeto(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0106) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0285) pass = false;    
    }
  } 

  //50 ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0277) pass = false;    
    }
  }

  return pass;
}

bool RazorAnalyzer::photonPassMediumIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0264) pass = false;    
    }
  }

  //50 ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
    }
  }

  return pass;
}

bool RazorAnalyzer::photonPassTightIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0099) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0264) pass = false;    
    }
  }

  //50 ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){    
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
    } else { 
      if(pho_HoverE[i] > 0.05) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
    }
  }

  return pass;
}


bool RazorAnalyzer::photonPassLooseID(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassMediumID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassTightID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


bool RazorAnalyzer::photonPassLooseIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 2.73, 1.85+exp( 0.0050*phoPt[i] + 0.2907 ), 1.34 + 0.0038*phoPt[i], true );
    } else {
      return photonPassesIsolation(i, 2.52, 18.51+exp( 0.0041*phoPt[i] + 0.8933 ), 0.95 + 0.0047*phoPt[i], true);
    }
  } 

  //50ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 2.44, 2.57+exp( 0.0044*phoPt[i] + 0.5809 ), 1.92 + 0.0043*phoPt[i], false );
    } else {
      return photonPassesIsolation(i, 1.84, 4.00+exp( 0.0040*phoPt[i] + 0.9402 ), 2.15 + 0.0041*phoPt[i], false);
    }
  }
}

bool RazorAnalyzer::photonPassMediumIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 0.70, 1.01+exp( 0.0050*phoPt[i] + 0.2907 ), 0.68 + 0.0038*phoPt[i], true);
    } else {
      return photonPassesIsolation(i, 0.58, 0.59+exp( 0.0041*phoPt[i] + 0.8933 ), 0.72 + 0.0047*phoPt[i], true);
    }
  }

  //50ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 1.31, 0.60+exp( 0.0044*phoPt[i] + 0.5809 ), 1.33 + 0.0043*phoPt[i], false);
    } else {
      return photonPassesIsolation(i, 1.25, 1.65+exp( 0.0040*phoPt[i] + 0.9402 ), 1.02 + 0.0041*phoPt[i], false);
    }
  }
}

bool RazorAnalyzer::photonPassTightIso(int i, bool use25nsCuts){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 0.09, 0.99+exp( 0.0050*phoPt[i] + 0.2907 ), 0.26 + 0.0038*phoPt[i], true);
    } else {
      return photonPassesIsolation(i, 0.02, 0.01+exp( 0.0041*phoPt[i] + 0.8933 ), 0.18 + 0.0047*phoPt[i], true);
    }
  } 

  //50ns cuts below
  else {
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 0.91, 0.33+exp( 0.0044*phoPt[i] + 0.5809 ), 0.61 + 0.0043*phoPt[i], false);
    } else {
      return photonPassesIsolation(i, 0.65, 0.93+exp( 0.0040*phoPt[i] + 0.9402 ), 0.54 + 0.0041*phoPt[i], false);
    }
  }
}

bool RazorAnalyzer::isLoosePhoton(int i, bool use25nsCuts){

  bool pass = true;
  if(!isLoosePhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhoton(int i, bool use25nsCuts){
  bool pass = true;

  if(!isMediumPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhoton(int i, bool use25nsCuts){
  bool pass = true;
  if (!isTightPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isLoosePhotonWithoutEleVeto(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso(i,use25nsCuts)) pass = false;

  return pass;
}


bool RazorAnalyzer::matchPhotonHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg1") {
    if ( 
	//Data filters
	pho_passHLTFilter[i][8] 
	//MC filters

	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg2") {
    if ( 
	//Data filters
	pho_passHLTFilter[i][9] 
	 ) {
      match = true;
    }
  }
   
  return match;
}
