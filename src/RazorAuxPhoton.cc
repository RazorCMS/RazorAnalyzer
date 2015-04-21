//LOCAL INCLUDES
#include "RazorAuxPhoton.hh"
#include "RazorAnalyzer.h"

bool RazorAnalyzer::photonPassesElectronVeto(int i){
    //use presence of a pixel seed as proxy for an electron veto
    return !(pho_hasPixelSeed[i]);
}

bool RazorAnalyzer::passEleVetoRun1( int i )
{
  return pho_passEleVeto[i];
};

//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;
    if(fabs(pho_superClusterEta[i]) < 1.0){
        effAreaChargedHadrons = 0.0130;
        effAreaNeutralHadrons = 0.0056;
        effAreaPhotons = 0.0896;
    }
    else if(fabs(pho_superClusterEta[i]) < 1.479){
        effAreaChargedHadrons = 0.0096;
        effAreaNeutralHadrons = 0.0107;
        effAreaPhotons = 0.0762;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.0){
        effAreaChargedHadrons = 0.0107;
        effAreaNeutralHadrons = 0.0019;
        effAreaPhotons = 0.0383;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.2){
        effAreaChargedHadrons = 0.0077;
        effAreaNeutralHadrons = 0.0011;
        effAreaPhotons = 0.0534;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.3){
        effAreaChargedHadrons = 0.0088;
        effAreaNeutralHadrons = 0.0077;
        effAreaPhotons = 0.0846;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.4){
        effAreaChargedHadrons = 0.0065;
        effAreaNeutralHadrons = 0.0178;
        effAreaPhotons = 0.1032;
    }
    else{
        effAreaChargedHadrons = 0.0030;
        effAreaNeutralHadrons = 0.1675;
        effAreaPhotons = 0.1598;
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

void RazorAnalyzer::getPhotonEffAreaRun1( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.012;
      effAreaNHad  = 0.030;
      effAreaPho   = 0.148;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.010;
      effAreaNHad  = 0.057;
      effAreaPho   = 0.130;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.014;
      effAreaNHad  = 0.039;
      effAreaPho   = 0.112;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.012;
      effAreaNHad  = 0.015;
      effAreaPho   = 0.216;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.016;
      effAreaNHad  = 0.024;
      effAreaPho   = 0.262;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.020;
      effAreaNHad  = 0.039;
      effAreaPho   = 0.260;
    }
  else
    {
      effAreaChHad = 0.012;
      effAreaNHad  = 0.072;
      effAreaPho   = 0.266;
    }
};

bool RazorAnalyzer::photonPassIsoRun1( int i )
{
  //Define variables
  float eta = pho_superClusterEta[i];
  float pt  = phoPt[i]; 
  bool _isEB = false;
  if ( fabs( eta ) < 1.44 ) _isEB = true;
  //get effective area for isolation calculations                                                
  double effAreaChargedHadrons = 0.0;
  double effAreaNeutralHadrons = 0.0;
  double effAreaPhotons = 0.0;
  getPhotonEffAreaRun1( eta, effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons );
  /*
  std::cout << "chA:" << effAreaChargedHadrons << " nA: " << effAreaNeutralHadrons << " phoA: " << effAreaPhotons << std::endl;
  std::cout << "chHad: " << pho_sumChargedHadronPt[i] << std::endl;
  std::cout << "nHad: " << pho_sumNeutralHadronEt[i] << std::endl;
  std::cout << "pho: " << pho_sumPhotonEt[i] << std::endl;
  std::cout << "pu: " << fixedGridRhoAll << std::endl;
  */
  //Compute Photon Isolation
  //Rho corrected PF charged hadron isolation
  double PFIsoCorrected_chHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoAll*effAreaChargedHadrons, 0.);
  //std::cout << "chHad Iso: " << PFIsoCorrected_chHad << std::endl;
  //Rho corrected PF neutral hadron isolation                                                                              
  double PFIsoCorrected_nHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoAll*effAreaNeutralHadrons, 0.);
  //std::cout << "nHad Iso: " <<PFIsoCorrected_nHad <<std::endl;
  //Rho corrected PF photon isolation                                                              
  double PFIsoCorrected_pho = max(pho_sumPhotonEt[i] - fixedGridRhoAll*effAreaPhotons, 0.);
  //std::cout << "pho Iso: " << PFIsoCorrected_pho <<std::endl;
  
  //Apply Isolation
  if ( _isEB )
    {
      //std::cout << "==EB==" << std::endl;
      //std::cout << "pass chHad iso?" << std::endl;
      //std::cout << "chHad Cut: " << EB_GetPFchHadIsoCut() << std::endl;
      if ( PFIsoCorrected_chHad > EB_GetPFchHadIsoCut() ) return false;
      //std::cout<< "pass nHad iso?" << std::endl;
      //std::cout << "nHad Cut: " << EB_GetPFnHadIsoCut( pt ) << std::endl;
      if ( PFIsoCorrected_nHad > EB_GetPFnHadIsoCut( pt ) ) return false;
      //std::cout << "Pho Cut: " << EB_GetPFphoIsoCut( pt ) << std::endl;
      //std::cout<< "pass pho iso?" <<std::endl;
      if ( PFIsoCorrected_pho > EB_GetPFphoIsoCut( pt ) ) return false;
      //std::cout<< "EB Yes" <<std::endl;
    }
  else
    {
      //std::cout << "==EE==" << std::endl;
      //std::cout << "pass chHad iso?" << std::endl;
      //std::cout << "chHad Cut: " << EE_GetPFchHadIsoCut() << std::endl;
      if ( PFIsoCorrected_chHad > EE_GetPFchHadIsoCut() ) return false;
      //std::cout << "nHad Cut: " << EE_GetPFnHadIsoCut( pt ) << std::endl;
      //std::cout<< "pass nHad iso?" <<std::endl;
      if ( PFIsoCorrected_nHad > EE_GetPFnHadIsoCut( pt ) ) return false;
      //std::cout<< "EE Yes" <<std::endl;
      //No Photn Isolation requirement;
    }
  
  return true;
};

//Cut Based Photon ID WP 90/85 (EB/EE) recommendation from EGamma: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
bool RazorAnalyzer::isGoodPhotonRun1( int i, bool _iso = false)
{
  bool _isEB = false;
  if ( fabs( pho_superClusterEta[i] ) < 1.44 ) _isEB = true;
  if ( _isEB )
    {
      //EB
      //std::cout<< "=====EB====" << std::endl;
      //std::cout << "pass Ele Veto?" << std::endl;
      if ( EB_EleVeto && !passEleVetoRun1( i ) ) return false;
      //std::cout << "yes!, pass HoverE?" << std::endl;
      if ( pho_HoverE[i] > EB_HoverECut ) return false;// HoverE Cut
      //std::cout << "yes!, pass ietaieta?" << std::endl;
      if ( phoSigmaIetaIeta[i] > EB_SigmaIetaIetaCut ) return false;// SigmaIetaIeta Cut
      //if ( phoFull5x5SigmaIetaIeta[i] > EB_SigmaIetaIetaCut ) return false;
      //std::cout << "yes!" << std::endl;
      if ( _iso && !photonPassIsoRun1( i ) ) return false;//Apply Isolation if flag (_iso) is true
      //if ( _iso) std::cout << "yes, passed ISO" << std::endl;
    }
  else
    {
      //EE
      //std::cout << "=====EE====" << std::endl;
      //std::cout << "pass Ele Veto?" << std::endl;
      if ( EE_EleVeto && !passEleVetoRun1( i ) ) return false;
      //std::cout<< "yes!, pass HoverE?" << std::endl;
      if ( pho_HoverE[i] > EE_HoverECut ) return false;// HoverE Cut 
      //std::cout<< "yes!, pass ietaieta?" << std::endl;
      if ( phoSigmaIetaIeta[i] > EE_SigmaIetaIetaCut ) return false;// SigmaIetaIeta Cut
      //if ( phoFull5x5SigmaIetaIeta[i] > EE_SigmaIetaIetaCut ) return false;
      //std::cout<< "yes!" << std::endl;
      if ( _iso && !photonPassIsoRun1( i ) ) return false;//Apply Isolation if flag (_iso) is true
      //if ( _iso) std::cout<< "yes, passed ISO" << std::endl;
    }

  return true;
};


bool RazorAnalyzer::photonPassesRunOneIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;
    if(fabs(pho_superClusterEta[i]) < 1.0){
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.030;
        effAreaPhotons = 0.148;
    }
    else if(fabs(pho_superClusterEta[i]) < 1.479){
        effAreaChargedHadrons = 0.010;
        effAreaNeutralHadrons = 0.057;
        effAreaPhotons = 0.130;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.0){
        effAreaChargedHadrons = 0.014;
        effAreaNeutralHadrons = 0.039;
        effAreaPhotons = 0.112;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.2){
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.015;
        effAreaPhotons = 0.216;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.3){
        effAreaChargedHadrons = 0.016;
        effAreaNeutralHadrons = 0.024;
        effAreaPhotons = 0.262;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.4){
        effAreaChargedHadrons = 0.020;
        effAreaNeutralHadrons = 0.039;
        effAreaPhotons = 0.260;
    }
    else{
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.072;
        effAreaPhotons = 0.266;
    }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    
    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;
    
    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = max(pho_sumPhotonEt[i] - fixedGridRhoAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoFull5x5SigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //Isolation
    if(!photonPassesIsolation(i, PFChHadIsoCut, PFNeuHadIsoCut, PFPhotIsoCut)) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesRunOneCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoSigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //Isolation
    if(!photonPassesRunOneIsolation(i, PFChHadIsoCut, PFNeuHadIsoCut, PFPhotIsoCut)) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesCutsBasedPhotonIDNoIsoCuts(int i, double HoverECut, double SigmaIetaIetaCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoFull5x5SigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::isLoosePhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.048, 0.0106, 2.56, 3.74+0.0025*phoPt[i], 2.68+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.069, 0.0266, 3.12, 17.11+0.0118*phoPt[i], 2.70+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isMediumPhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.032, 0.0101, 1.90, 2.96+0.0025*phoPt[i], 1.39+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.0166, 0.0264, 1.95, 4.42+0.0118*phoPt[i], 1.89+0.0059*phoPt[i]);
}

//photon ID with electron veto and HoverE and sigmaIetaIeta cuts only
bool RazorAnalyzer::isMediumPhotonNoIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonIDNoIsoCuts(i, 0.032, 0.0101);
    }
    //endcap photons
    return passesCutsBasedPhotonIDNoIsoCuts(i, 0.0166,0.0264);
}

bool RazorAnalyzer::photonPassesLooseIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 2.56, 3.74+0.0025*phoPt[i], 2.68+0.001*phoPt[i] );
    }
    //endcap photons
    return photonPassesIsolation(i, 3.12, 17.11+0.0118*phoPt[i], 2.70+0.0059*phoPt[i]);
}

bool RazorAnalyzer::photonPassesMediumIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 1.90, 2.96+0.0025*phoPt[i], 1.39+0.001*phoPt[i]);
    }
    //endcap photons
    return photonPassesIsolation(i, 1.95, 4.42+0.0118*phoPt[i], 1.89+0.0059*phoPt[i]);
}

bool RazorAnalyzer::photonPassesTightIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 1.86, 2.64+0.0025*phoPt[i], 1.20+0.001*phoPt[i]);
    }
    //endcap photons
    return photonPassesIsolation(i, 1.68, 4.42+0.0118*phoPt[i], 1.03+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isTightPhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.011, 0.0099, 1.86, 2.64+0.0025*phoPt[i], 1.20+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.015, 0.0263, 1.68, 4.42+0.0118*phoPt[i], 1.03+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isMediumRunOnePhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesRunOneCutsBasedPhotonID(i, 0.05, 0.011, 1.5, 1.0+0.04*phoPt[i], 0.7*0.005*phoPt[i]);
    }
    //endcap photons
    return passesRunOneCutsBasedPhotonID(i, 0.05, 0.033, 1.2, 1.5+0.04*phoPt[i], 1.0+0.005*phoPt[i]);
}
