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



//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    //get the effective areas. results are passed to variables by reference
    getPhotonEffAreaRun2( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);

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

bool RazorAnalyzer::photonPassLooseID(int i){

  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0277) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassMediumID(int i){
  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassTightID(int i){
  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


bool RazorAnalyzer::photonPassLooseIso(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 2.44, 2.57+exp( 0.0044*phoPt[i] + 0.5809 ), 1.92 + 0.0043*phoPt[i] );
    }
    //endcap photons
    return photonPassesIsolation(i, 1.84, 4.00+exp( 0.0040*phoPt[i] + 0.9402 ), 2.15 + 0.0041*phoPt[i]);
}

bool RazorAnalyzer::photonPassMediumIso(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 1.31, 0.60+exp( 0.0044*phoPt[i] + 0.5809 ), 1.33 + 0.0043*phoPt[i] );
    }
    //endcap photons
    return photonPassesIsolation(i, 1.25, 1.65+exp( 0.0040*phoPt[i] + 0.9402 ), 1.02 + 0.0041*phoPt[i]);
}

bool RazorAnalyzer::photonPassTightIso(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 0.91, 0.33+exp( 0.0044*phoPt[i] + 0.5809 ), 0.61 + 0.0043*phoPt[i] );
    }
    //endcap photons
    return photonPassesIsolation(i, 0.65, 0.93+exp( 0.0040*phoPt[i] + 0.9402 ), 0.54 + 0.0041*phoPt[i]);
}

bool RazorAnalyzer::isLoosePhoton(int i){

  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0277) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;
  if(!photonPassLooseIso(i)) pass = false;

  //cout << "Loose Photon ? : " << bool(pho_HoverE[i] > 0.05) << " " << bool(phoFull5x5SigmaIetaIeta[i] > 0.0103) << " " << photonPassesElectronVeto(i) << " " << photonPassLooseIso(i) << "\n";

  return pass;
}

bool RazorAnalyzer::isMediumPhoton(int i){
  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;
  if(!photonPassMediumIso(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhoton(int i){
  bool pass = true;

  if(fabs(pho_superClusterEta[i]) < 1.479){    
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0100) pass = false;   
  } else { 
    if(pho_HoverE[i] > 0.05) pass = false;
    if(phoFull5x5SigmaIetaIeta[i] > 0.0267) pass = false;    
  }

  if(!photonPassesElectronVeto(i)) pass = false;
  if(!photonPassTightIso(i)) pass = false;

  return pass;
}

