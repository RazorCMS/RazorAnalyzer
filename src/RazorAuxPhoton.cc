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



void RazorAnalyzer::getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0234;
      effAreaNHad  = 0.0053;
      effAreaPho   = 0.0896;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0189;
      effAreaNHad  = 0.0103;
      effAreaPho   = 0.0762;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0171;
      effAreaNHad  = 0.0057;
      effAreaPho   = 0.0383;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0129;
      effAreaNHad  = 0.0070;
      effAreaPho   = 0.0534;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0110;
      effAreaNHad  = 0.0152;
      effAreaPho   = 0.0846;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0074;
      effAreaNHad  = 0.0230;
      effAreaPho   = 0.1032;
    }
  else
    {
      effAreaChHad = 0.0035;
      effAreaNHad  = 0.1709;
      effAreaPho   = 0.1598;
    }
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
  getPhotonEffAreaRun2( eta, effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons );
  
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
      std::cout << "run2 chHad Iso: " << PFIsoCorrected_chHad << std::endl;
      std::cout << "run2 nHad Iso: " <<PFIsoCorrected_nHad <<std::endl;
      std::cout << "run2 pho Iso: " << PFIsoCorrected_pho <<std::endl;
    }
  
  if ( PFIsoCorrected_chHad > GetchHadIsoCut( wp, eta ) ) return false;
  if ( PFIsoCorrected_nHad  > GetnHadIsoCut( wp, eta, pt ) ) return false;
  if ( PFIsoCorrected_pho   > GetPhotonIsoCut( wp, eta, pt ) ) return false;
  
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
	  if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EB[0] )
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
          if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EB[1] )
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
          if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EB[2] )
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
          if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EE[0] )
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
          if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EE[1] )
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
          if ( phoFull5x5SigmaIetaIeta[i] > SigmaIetaIeta_EE[2] )
            {
              if ( _debug ) std::cout << "EB, ID: failed run2 SigmaIetaIeta: " << phoSigmaIetaIeta[i] << std::endl;
              return false;// SigmaIetaIeta Cut
            }
        }
    }
  
  if ( _iso && !photonPassIsoRun2( i, wp, _debug ) )return false;
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

double RazorAnalyzer::getPhotonScaleCorrectionRunOne7TeV(int run, double eta, double r9, double et) {

  double corr1 = 1.0;
  if ( fabs(eta) < 1 && r9 < 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9961;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9971;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9974;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9976;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9964;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9971;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9962;
  else if ( fabs(eta) < 1 && r9 < 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9961;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9929;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9939;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9942;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9944;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9932;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9939;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9929;
  else if ( fabs(eta) < 1 && r9 >= 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9929;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9981;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9973;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9981;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9980;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9952;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9965;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9964;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9965;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9881;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9874;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9881;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9881;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9852;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9866;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9865;
  else if ( fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9866;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9959;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9971;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9971;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9971;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9983;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9973;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9967;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9960;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9942;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9954;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9954;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9954;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9966;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9956;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9950;
  else if ( fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9943;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9985;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9980;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9991;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9991;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9981;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9992;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 177776 && run <= 178723) corr1 = 1.0004;
  else if ( fabs(eta) >= 2 && r9 < 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9993;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 160431 && run <= 165547) corr1 = 0.9940;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 165548 && run <= 167042) corr1 = 0.9936;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 167043 && run <= 172400) corr1 = 0.9947;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 172401 && run <= 173663) corr1 = 0.9947;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 173664 && run <= 176840) corr1 = 0.9937;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 176841 && run <= 177775) corr1 = 0.9948;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 177776 && run <= 178723) corr1 = 0.9960;
  else if ( fabs(eta) >= 2 && r9 >= 0.94 && run >= 178724 && run <= 180252) corr1 = 0.9949;

  double corr2 = 1.0;
  if ( fabs(eta) < 1  && r9 < 0.94 && 20 < et && et <=  33) corr2 =  1.0015;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 33 < et && et <=  39) corr2 =  1.0008;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 39 < et && et <=  45) corr2 =  1.0010;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 45 < et && et <=  50) corr2 =  1.0006;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 50 < et && et <=  58) corr2 =  0.9999;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 58 < et && et <=  100) corr2 =  0.9982;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 20 < et && et <=  35) corr2 =  1.0017;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 35 < et && et <=  43) corr2 =  1.0011;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 43 < et && et <=  50) corr2 =  1.0011;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 50 < et && et <=  55) corr2 =  1.0008;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 55 < et && et <=  100) corr2 =  0.9991;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 20 < et && et <=  33) corr2 =  1.0041;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 33 < et && et <=  39) corr2 =  1.0018;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 39 < et && et <=  45) corr2 =  0.9994;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 45 < et && et <=  50) corr2 =  0.9970;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 50 < et && et <=  58) corr2 =  0.9960;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 58 < et && et <=  100) corr2 =  0.9951;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 20 < et && et <= 40) corr2 =  1.0013;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 40 < et && et <= 50) corr2 =  1.0016;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 50 < et && et <= 100) corr2 = 1.0007;
  
  double corr = corr1 * corr2;
  return corr;
}

double RazorAnalyzer::getPhotonScaleCorrectionRunOne8TeV(int run, double eta, double r9, double et) {

  double corr1 = 1.0;

  if (fabs(eta) < 1 && r9 < 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9922;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9989;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9931;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9922;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9929;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9935;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9929;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9932;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9938;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9940;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9931;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9936;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9942;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9936;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9943;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9938;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9934;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9936;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9933;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9935;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9938;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9940;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9942;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9940;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9947;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9937;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9943;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9945;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9944;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9947;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9951;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9949;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9944;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9958;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9935;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9939;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9938;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9939;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9938;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9942;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9949;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9946;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9945;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9945;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9941;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9952;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9949;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9950;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9951;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9947;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9952;
  } else if (fabs(eta) < 1 && r9 >= 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9894;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9961;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9902;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9894;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9900;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9907;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9901;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9904;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9910;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9912;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9903;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9908;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9914;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9908;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9915;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9910;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9906;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9907;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9904;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9907;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9910;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9912;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9913;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9912;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9919;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9909;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9915;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9917;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9916;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9919;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9923;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9921;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9916;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9930;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9906;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9910;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9910;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9911;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9909;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9914;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9921;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9918;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9917;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9917;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9913;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9923;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9921;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9922;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9922;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9919;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9924;
  } else if (fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 < 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9982;
    else if ( run >= 190782 && run <= 191042) corr1 = 1.0014;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9963;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9982;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9970;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9975;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9973;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9979;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9992;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9976;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9968;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9993;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9983;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9980;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9983;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9977;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9970;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9970;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9973;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9976;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9985;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9983;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9970;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9984;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9971;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9969;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9992;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9983;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9974;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9977;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9989;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9994;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9962;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9936;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9991;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9999;
    else if ( run >= 204563 && run <= 205085) corr1 = 1.0006;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9991;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9995;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9982;
    else if ( run >= 205826 && run <= 206207) corr1 = 1.0012;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9995;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9995;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9997;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9998;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9997;
    else if ( run >= 207221 && run <= 207315) corr1 = 1.0006;
    else if ( run >= 207316 && run <= 207489) corr1 = 1.0001;
    else if ( run >= 207490 && run <= 207919) corr1 = 1.0006;
    else if ( run >= 207920 && run <= 208351) corr1 = 1.0000;
    else if ( run >= 208352 && run <= 208686) corr1 = 1.0011;
  } else if (fabs(eta) >= 1 && fabs(eta) <= 1.479 && r9 >= 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9876;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9909;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9858;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9876;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9865;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9870;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9867;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9874;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9887;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9871;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9863;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9888;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9878;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9874;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9877;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9871;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9864;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9865;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9868;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9870;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9879;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9877;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9865;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9879;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9866;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9864;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9886;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9878;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9869;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9871;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9883;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9888;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9857;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9831;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9886;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9894;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9901;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9886;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9890;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9877;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9908;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9890;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9890;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9892;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9893;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9892;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9901;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9896;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9901;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9895;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9906;
  } else if (fabs(eta) > 1.479 && fabs(eta) < 2 && r9 < 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9919;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9917;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9944;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9917;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9906;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9918;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9905;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9940;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9963;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9948;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9930;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9929;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9962;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9942;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9944;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9923;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9940;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9920;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9921;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9916;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9909;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9938;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9924;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9942;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9924;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9914;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9920;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9919;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9929;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9914;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9950;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9920;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9924;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9769;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9914;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9931;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9934;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9932;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9916;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9930;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9945;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9906;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9957;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9949;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9910;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9935;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9934;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9944;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9930;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9930;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9934;
  } else if (fabs(eta) > 1.479 && fabs(eta) < 2 && r9 >= 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9856;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9854;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9880;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9853;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9842;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9855;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9841;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9876;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9900;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9885;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9867;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9866;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9898;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9878;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9880;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9859;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9877;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9856;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9857;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9853;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9845;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9875;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9860;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9878;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9861;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9850;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9857;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9856;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9866;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9850;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9887;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9857;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9860;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9705;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9851;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9868;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9871;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9868;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9852;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9866;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9882;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9842;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9894;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9886;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9846;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9872;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9870;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9881;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9866;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9867;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9871;
  } else if (fabs(eta) > 2 && r9 < 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9862;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9810;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9845;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9841;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9843;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9855;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9836;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9838;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9853;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9866;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9851;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9856;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9846;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9846;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9831;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9841;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9868;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9855;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9858;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9853;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9853;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9873;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9867;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9858;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9864;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9883;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9876;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9870;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9883;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9877;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9874;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9861;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9865;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9778;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9827;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9864;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9839;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9832;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9816;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9842;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9866;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9844;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9873;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9841;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9839;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9863;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9848;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9841;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9819;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9863;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9859;
  } else if (fabs(eta) > 2 && r9 >= 0.94) {
    if ( run >= 190645 && run <= 190781) corr1 = 0.9812;
    else if ( run >= 190782 && run <= 191042) corr1 = 0.9760;
    else if ( run >= 191043 && run <= 191720) corr1 = 0.9795;
    else if ( run >= 191721 && run <= 193833) corr1 = 0.9791;
    else if ( run >= 193834 && run <= 194116) corr1 = 0.9793;
    else if ( run >= 194117 && run <= 194427) corr1 = 0.9805;
    else if ( run >= 194428 && run <= 194618) corr1 = 0.9786;
    else if ( run >= 194619 && run <= 194789) corr1 = 0.9788;
    else if ( run >= 194790 && run <= 195111) corr1 = 0.9803;
    else if ( run >= 195112 && run <= 195377) corr1 = 0.9817;
    else if ( run >= 195378 && run <= 195398) corr1 = 0.9801;
    else if ( run >= 195399 && run <= 195657) corr1 = 0.9806;
    else if ( run >= 195658 && run <= 195918) corr1 = 0.9797;
    else if ( run >= 195919 && run <= 196198) corr1 = 0.9796;
    else if ( run >= 196199 && run <= 196356) corr1 = 0.9781;
    else if ( run >= 196357 && run <= 198115) corr1 = 0.9791;
    else if ( run >= 198116 && run <= 198940) corr1 = 0.9818;
    else if ( run >= 198941 && run <= 199317) corr1 = 0.9805;
    else if ( run >= 199318 && run <= 199428) corr1 = 0.9808;
    else if ( run >= 199429 && run <= 199697) corr1 = 0.9803;
    else if ( run >= 199698 && run <= 199832) corr1 = 0.9803;
    else if ( run >= 199833 && run <= 199960) corr1 = 0.9823;
    else if ( run >= 199961 && run <= 200151) corr1 = 0.9817;
    else if ( run >= 200152 && run <= 200490) corr1 = 0.9809;
    else if ( run >= 200491 && run <= 200991) corr1 = 0.9814;
    else if ( run >= 200992 && run <= 201201) corr1 = 0.9833;
    else if ( run >= 201202 && run <= 201624) corr1 = 0.9827;
    else if ( run >= 201625 && run <= 201707) corr1 = 0.9820;
    else if ( run >= 201708 && run <= 202059) corr1 = 0.9833;
    else if ( run >= 202060 && run <= 202204) corr1 = 0.9827;
    else if ( run >= 202205 && run <= 202332) corr1 = 0.9825;
    else if ( run >= 202333 && run <= 202972) corr1 = 0.9812;
    else if ( run >= 202973 && run <= 203002) corr1 = 0.9816;
    else if ( run >= 203003 && run <= 203852) corr1 = 0.9728;
    else if ( run >= 203853 && run <= 204099) corr1 = 0.9777;
    else if ( run >= 204100 && run <= 204562) corr1 = 0.9814;
    else if ( run >= 204563 && run <= 205085) corr1 = 0.9789;
    else if ( run >= 205086 && run <= 205310) corr1 = 0.9782;
    else if ( run >= 205311 && run <= 205617) corr1 = 0.9766;
    else if ( run >= 205618 && run <= 205825) corr1 = 0.9792;
    else if ( run >= 205826 && run <= 206207) corr1 = 0.9816;
    else if ( run >= 206208 && run <= 206389) corr1 = 0.9794;
    else if ( run >= 206390 && run <= 206483) corr1 = 0.9824;
    else if ( run >= 206484 && run <= 206597) corr1 = 0.9791;
    else if ( run >= 206598 && run <= 206896) corr1 = 0.9789;
    else if ( run >= 206897 && run <= 207220) corr1 = 0.9813;
    else if ( run >= 207221 && run <= 207315) corr1 = 0.9798;
    else if ( run >= 207316 && run <= 207489) corr1 = 0.9791;
    else if ( run >= 207490 && run <= 207919) corr1 = 0.9769;
    else if ( run >= 207920 && run <= 208351) corr1 = 0.9813;
    else if ( run >= 208352 && run <= 208686) corr1 = 0.9809;
  }

  double corr2 = 1.0;
  if ( fabs(eta) < 1  && r9 < 0.94 && 20 < et && et <=  33) corr2 =  1.0015;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 33 < et && et <=  39) corr2 =  1.0008;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 39 < et && et <=  45) corr2 =  1.0010;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 45 < et && et <=  50) corr2 =  1.0006;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 50 < et && et <=  58) corr2 =  0.9999;
  else if ( fabs(eta) < 1  && r9 < 0.94 && 58 < et && et <=  100) corr2 =  0.9982;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 20 < et && et <=  35) corr2 =  1.0017;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 35 < et && et <=  43) corr2 =  1.0011;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 43 < et && et <=  50) corr2 =  1.0011;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 50 < et && et <=  55) corr2 =  1.0008;
  else if ( fabs(eta) < 1  && r9 >= 0.94 && 55 < et && et <=  100) corr2 =  0.9991;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 20 < et && et <=  33) corr2 =  1.0041;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 33 < et && et <=  39) corr2 =  1.0018;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 39 < et && et <=  45) corr2 =  0.9994;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 45 < et && et <=  50) corr2 =  0.9970;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 50 < et && et <=  58) corr2 =  0.9960;
  else if ( fabs(eta) > 1  && r9 < 0.94 && 58 < et && et <=  100) corr2 =  0.9951;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 20 < et && et <= 40) corr2 =  1.0013;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 40 < et && et <= 50) corr2 =  1.0016;
  else if ( fabs(eta) > 1  && r9 >= 0.94 && 50 < et && et <= 100) corr2 = 1.0007;

  double corr = corr1 * corr2;
  return corr;

}

