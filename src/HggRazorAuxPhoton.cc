//C++ INCLUDES
#include <iostream>
#include <math.h>
//ROOT INCLUDES
//LOCAL INCLUDES
#include "HggRazorAuxPhoton.hh"

//std::map< EA, std::vector<float> > effAreas;
float effArea[3][7];

float EB_GetPFchHadIsoCut( ){ return EB_PFchHadIsoCut; };
float EB_GetPFnHadIsoCut( float pt ){ return (EB_PFnHadIsoConst + EB_PFnHadIsoSlope*pt); };
float EB_GetPFphoIsoCut( float pt ){ return (EB_PFphoIsoConst + EB_PFphoIsoSlope*pt); };

float EE_GetPFchHadIsoCut( ){ return EE_PFchHadIsoCut; };
float EE_GetPFnHadIsoCut( float pt ){ return (EE_PFnHadIsoConst + EE_PFnHadIsoSlope*pt); };


//R u n 2  F u n c t i o n s
//--------------------------
void InitEffArea( )
{
  //F i l l i n g   c h H a d
  //-------------------------
  effArea[0][0] = 0.0234;
  effArea[0][1] = 0.0189;
  effArea[0][2] = 0.0171;
  effArea[0][3] = 0.0129;
  effArea[0][4] = 0.0110;
  effArea[0][5] = 0.0074;
  effArea[0][6] = 0.0035;

  //F i l l i n g   n H a d
  //-----------------------                            
  effArea[1][0] = 0.0053;
  effArea[1][1] = 0.0103;
  effArea[1][2] = 0.0057;
  effArea[1][3] = 0.0070;
  effArea[1][4] = 0.0152;
  effArea[1][5] = 0.0232;
  effArea[1][6] = 0.1709;
  
  // F i l l i n g   p h o t o n
  //----------------------------
  effArea[2][0] = 0.0780;
  effArea[2][1] = 0.0629;
  effArea[2][2] = 0.0264;
  effArea[2][3] = 0.0462;
  effArea[2][4] = 0.0740;
  effArea[2][5] = 0.0924;
  effArea[2][6] = 0.1484;
};

// G e t   E f f e c t i v e   A r e a   B i n
//--------------------------------------------
bool GetEA( float eta, double& EAcHad, double& EAnHad, double& EAphoton  )
{
  int eabin = n_EA_Threshold;
  for( int i = 0; i < n_EA_Threshold; i++ )
    {
      if ( fabs( eta ) < EA_Threshold[i] )
	{
	  eabin = i;
	  break;
	}
    }
  
  EAcHad = effArea[0][eabin];
  EAnHad = effArea[1][eabin];
  EAphoton = effArea[2][eabin];
  if ( EAcHad <= .0 || EAnHad <= .0 || EAphoton <= .0 ) return false;
  return true;
};

// I s o l a t i o n s
//--------------------
float GetchHadIsoCut( WP wp, float eta )
{
  bool _isEB = true;
  if ( fabs( eta ) > 1.479 ) _isEB = false;
  if ( _isEB )
    {
      if ( wp == WP::Loose ) return 2.67;
      if ( wp == WP::Medium )return 1.79;
      if ( wp == WP::Tight ) return 1.66;
    }
  else
    {
      if ( wp == WP::Loose ) return 1.79;
      if ( wp == WP::Medium )return 1.09;
      if ( wp == WP::Tight ) return 1.04;
    }
  
  return -1.;
};

float GetnHadIsoCut( WP wp, float eta, float pt )
{
  bool _isEB = true;
  if ( fabs( eta ) > 1.479 ) _isEB = false;
  if ( _isEB )
    {
      if ( wp == WP::Loose ) return 7.23 + exp( 0.0028*pt + 0.5408 );
      if ( wp == WP::Medium )return 0.16 + exp( 0.0028*pt + 0.5408 );
      if ( wp == WP::Tight ) return 0.14 + exp( 0.0028*pt + 0.5408 );
    }
  else
    {
      if ( wp == WP::Loose ) return 8.89 + 0.01725*pt;
      if ( wp == WP::Medium )return 4.31 + 0.01720*pt;
      if ( wp == WP::Tight ) return 3.89 + 0.01720*pt;
    }

  return -1.;
};

float GetPhotonIsoCut( WP wp, float eta, float pt )
{
  bool _isEB = true;
  if ( fabs( eta ) > 1.479 ) _isEB = false;
  if ( _isEB )
    {
      if ( wp == WP::Loose ) return 2.11 + 0.0014*pt;
      if ( wp == WP::Medium )return 1.90 + 0.0014*pt;
      if ( wp == WP::Tight ) return 1.40 + 0.0014*pt;
    }
  else
    {
      if ( wp == WP::Loose ) return 3.09 + 0.0091*pt;
      if ( wp == WP::Medium )return 1.90 + 0.0091*pt;
      if ( wp == WP::Tight ) return 1.40 + 0.0091*pt;
    }

  return -1.;
};
