//LOCAL INCLUDES
#include "HggRazorAuxPhoton.hh"
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

bool RazorAnalyzer::photonPassIsoRun1( int i , bool _debug )
{
  //Define variables
  float eta = pho_superClusterEta[i];
  //float pt  = phoE[i]/cosh(phoEta[i]); 
  float pt  = phoPt[i];//default pt
  bool _isEB = false;
  if ( fabs( eta ) < 1.48 ) _isEB = true;
  //get effective area for isolation calculations                                                
  double effAreaChargedHadrons = 0.0;
  double effAreaNeutralHadrons = 0.0;
  double effAreaPhotons = 0.0;
  getPhotonEffAreaRun1( eta, effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons );
  if ( _debug )
    {
      std::cout << "chA:" << effAreaChargedHadrons << " nA: " << effAreaNeutralHadrons << " phoA: " << effAreaPhotons << std::endl;
      std::cout << "chHad: " << pho_sumChargedHadronPt[i] << std::endl;
      std::cout << "nHad: " << pho_sumNeutralHadronEt[i] << std::endl;
      std::cout << "pho: " << pho_sumPhotonEt[i] << std::endl;
      std::cout << "pu: " << fixedGridRhoAll << std::endl;
    }
  //Compute Photon Isolation
  //Rho corrected PF charged hadron isolation
  double PFIsoCorrected_chHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoAll*effAreaChargedHadrons, 0.);
  //Rho corrected PF neutral hadron isolation                                                                              
  double PFIsoCorrected_nHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoAll*effAreaNeutralHadrons, 0.);
  //Rho corrected PF photon isolation                                                              
  double PFIsoCorrected_pho = max(pho_sumPhotonEt[i] - fixedGridRhoAll*effAreaPhotons, 0.);
  if ( _debug )
    {
      std::cout << "chHad Iso: " << PFIsoCorrected_chHad << std::endl;
      std::cout << "nHad Iso: " <<PFIsoCorrected_nHad <<std::endl;
      std::cout << "pho Iso: " << PFIsoCorrected_pho <<std::endl; 
    }
  
  //Apply Isolation
  if ( _isEB )
    {
      if ( PFIsoCorrected_chHad > EB_GetPFchHadIsoCut() )
	{
	  if ( _debug ) std::cout << "EB, Iso: failed chHadIso, cut @ " << EB_GetPFchHadIsoCut() << std::endl;
	  return false;
	}
      if ( PFIsoCorrected_nHad > EB_GetPFnHadIsoCut( pt ) )
	{
	  if ( _debug ) std::cout << "EB, Iso: failed nHadIso, cut @ " << EB_GetPFnHadIsoCut( pt ) << std::endl;
	  return false;
	}
      if ( PFIsoCorrected_pho > EB_GetPFphoIsoCut( pt ) )
	{
	  if ( _debug ) std::cout << "EB, Iso: failed phoIso, cut @ " << EB_GetPFphoIsoCut( pt ) << std::endl;
	  return false;
	}
    }
  else
    {
      if ( PFIsoCorrected_chHad > EE_GetPFchHadIsoCut() )
	{
	  if ( _debug ) std::cout << "EE, Iso: failed chHadIso, cut @ " << EE_GetPFchHadIsoCut() << std::endl;
	  return false;
	}
      if ( PFIsoCorrected_nHad > EE_GetPFnHadIsoCut( pt ) )
	{
	  if ( _debug ) std::cout << "EE, Iso: failed nHadIso, cut @ " << EE_GetPFnHadIsoCut( pt ) << std::endl;
	  return false;
	}
    }
  
  return true;
};

//Cut Based Photon ID WP 90/85 (EB/EE) recommendation from EGamma: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
bool RazorAnalyzer::isGoodPhotonRun1( int i, bool _iso = false, bool _debug = false)
{
  bool _isEB = false;
  if ( fabs( pho_superClusterEta[i] ) < 1.48 ) _isEB = true;
  if ( _isEB )
    {
      //EB
      if ( EB_EleVeto && !passEleVetoRun1( i ) ) 
	{
	  if ( _debug ) std::cout << "EB, ID: failed CSEV" << std::endl;
	  return false;
	}
      if ( pho_HoverE[i] > EB_HoverECut )
	{
	  if ( _debug ) std::cout << "EB, ID: failed HoveE: " << pho_HoverE[i] << std::endl;
	  return false;// HoverE Cut
	}
      if ( phoSigmaIetaIeta[i] > EB_SigmaIetaIetaCut )
	{
	  if ( _debug ) std::cout << "EB, ID: failed SigmaIetaIeta: " << phoSigmaIetaIeta[i] <<  std::endl;
	  return false;// SigmaIetaIeta Cut
	}
      if ( _iso && !photonPassIsoRun1( i, _debug ) ) return false;//Apply Isolation if flag (_iso) is true
      //if ( _iso) std::cout << "yes, passed ISO" << std::endl;
    }
  else
    {
      //EE
      if ( EE_EleVeto && !passEleVetoRun1( i ) )
	{
	  if ( _debug ) std::cout << "EE, ID: failed CSEV" << std::endl;
	  return false;
	}
      if ( pho_HoverE[i] > EE_HoverECut )
	{
	  if ( _debug ) std::cout << "EE, ID: failed HoveE: " << pho_HoverE[i] << std::endl;
	  return false;// HoverE Cut 
	}
      if ( phoSigmaIetaIeta[i] > EE_SigmaIetaIetaCut )
	{
	  if ( _debug ) std::cout << "EE, ID: failed SigmaIetaIeta: " << phoSigmaIetaIeta[i] <<  std::endl;
	  return false;// SigmaIetaIeta Cut
	}
      if ( _iso && !photonPassIsoRun1( i, _debug ) ) return false;//Apply Isolation if flag (_iso) is true
    }

  return true;
};

// R u n 2   H g g   P h o t o n   ID
//-----------------------------------
bool RazorAnalyzer::isGoodPhotonRun2( int i, bool _iso, WP wp, bool _debug )
{
  bool _isEB = false;
  if ( fabs( pho_superClusterEta[i] ) < 1.48 ) _isEB = true;
  if ( _isEB )
    {
      // L o o s e  W o r k i n g   P o i n t
      //-------------------------------------
      if ( wp == WP::Loose )
	{
	  if ( pho_HoverE[i] > HoverE_EB[0] )
	    {
	      if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
	      return false;// HoverE Cut
	    }
	  if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EB[0] )
	    {
	      if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
	    }
	}
      // M e d i u m  W o r k i n g   P o i n t
      //-------------------------------------
      if ( wp == WP::Medium )
        {
          if ( pho_HoverE[i] > HoverE_EB[1] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
              return false;// HoverE Cut
            }
          if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EB[1] )
            {
	      if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
	}
      // T i g h t  W o r k i n g   P o i n t
      //-------------------------------------
      if ( wp == WP::Tight )
        {
          if ( pho_HoverE[i] > HoverE_EB[2] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
              return false;// HoverE Cut
            }
          if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EB[2] )
            {
	      if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
	}
    }
  else
    {
      // L o o s e  W o r k i n g   P o i n t                    
      //------------------------------------- 
      if ( wp == WP::Loose )
        {
          if ( pho_HoverE[i] > HoverE_EE[0] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
              return false;// HoverE Cut
            }
          if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EB[0] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
        }
      // M e d i u m  W o r k i n g   P o i n t 
      //-------------------------------------         
      if ( wp == WP::Medium )
        {
          if ( pho_HoverE[i] > HoverE_EE[1] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
              return false;// HoverE Cut
            }
          if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EE[1] )
            {
              if ( _debug ) std::cout << "EB, ID: failed SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
        }
      // T i g h t  W o r k i n g   P o i n t
      //-------------------------------------
      if ( wp == WP::Tight )
        {
          if ( pho_HoverE[i] > HoverE_EE[2] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 HoverE: " << pho_HoverE[i] << std::endl;
              return false;// HoverE Cut
            }
          if ( phoSigmaIetaIeta[i] > SigmaIetaIeta_EE[2] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
        }
    }
  
  if ( _iso && photonPassIsoRun2( i, wp, _debug ) )return false;
  return true;
};

//R u n 2  H g g R a z o r  I s o l a t  i o n
//--------------------------------------------
bool RazorAnalyzer::photonPassIsoRun2( int i, WP wp ,bool _debug )
{
  float eta = pho_superClusterEta[i];
  float pt  = phoPt[i];//default pt
  
  double effAreaChargedHadrons = -1.;
  double effAreaNeutralHadrons = -1.;
  double effAreaPhotons        = -1.;
  GetEA( eta, effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons );
  
  //Compute Photon Isolation
  //Rho corrected PF charged hadron isolation
  double PFIsoCorrected_chHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoAll*effAreaChargedHadrons, 0.);
  //Rho corrected PF neutral hadron isolation
  double PFIsoCorrected_nHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoAll*effAreaNeutralHadrons, 0.);
  //Rho corrected PF photon isolation
  double PFIsoCorrected_pho = max(pho_sumPhotonEt[i] - fixedGridRhoAll*effAreaPhotons, 0.);
  
  if ( _debug )
    {
      std::cout << "run2 chHad Iso: " << PFIsoCorrected_chHad << std::endl;
      std::cout << "run2 nHad Iso: " <<PFIsoCorrected_nHad <<std::endl;
      std::cout << "run2 pho Iso: " << PFIsoCorrected_pho <<std::endl;
    }
  
  if ( PFIsoCorrected_chHad > GetchHadIsoCut( wp, eta ) ) return false;
  if ( PFIsoCorrected_nHad  > GetnHadIsoCut( wp, eta, pt ) ) return false;
  if ( PFIsoCorrected_pho   > GetPhotonIsoCut( wp, eta, pt ) ) return false;
  
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

bool RazorAnalyzer::isTightRunOnePhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesRunOneCutsBasedPhotonID(i, 0.05, 0.011, 0.7, 0.4+0.04*phoPt[i], 0.5*0.005*phoPt[i]);
    }
    //endcap photons
    return passesRunOneCutsBasedPhotonID(i, 0.05, 0.031, 0.5, 1.5+0.04*phoPt[i], 1.0+0.005*phoPt[i]);
}
