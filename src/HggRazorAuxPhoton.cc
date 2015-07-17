//C++ INCLUDES
#include <iostream>
#include <math.h>
//ROOT INCLUDES
//LOCAL INCLUDES
#include "HggRazorAuxPhoton.hh"

//std::map< EA, std::vector<float> > effAreas;
// float effArea[3][7];

float EB_GetPFchHadIsoCut( ){ return EB_PFchHadIsoCut; };
float EB_GetPFchHadIsoCut_VL( ){ return 2.0*EB_PFchHadIsoCut; };
float EB_GetPFchHadIsoCut_M( ){ return 1.5; };

float EB_GetPFnHadIsoCut( float pt ){ return (EB_PFnHadIsoConst + EB_PFnHadIsoSlope*pt); };
float EB_GetPFnHadIsoCut_VL( float pt ){ return (2.*EB_PFphoIsoConst + EB_PFnHadIsoSlope*pt); };
float EB_GetPFnHadIsoCut_M( float pt ){ return (1.0 + EB_PFnHadIsoSlope*pt); };

float EB_GetPFphoIsoCut( float pt ){ return (EB_PFphoIsoConst + EB_PFphoIsoSlope*pt); };
float EB_GetPFphoIsoCut_VL( float pt ){ return (2.*EB_PFphoIsoConst + EB_PFphoIsoSlope*pt); };
float EB_GetPFphoIsoCut_M( float pt ){ return (0.7 + EB_PFphoIsoSlope*pt); };

float EE_GetPFchHadIsoCut( ){ return EE_PFchHadIsoCut; };
float EE_GetPFchHadIsoCut_VL( ){ return 2*EE_PFchHadIsoCut; };
float EE_GetPFchHadIsoCut_M( ){ return 1.2; };

float EE_GetPFnHadIsoCut( float pt ){ return (EE_PFnHadIsoConst + EE_PFnHadIsoSlope*pt); };
float EE_GetPFnHadIsoCut_VL( float pt ){ return (2.*EE_PFnHadIsoConst + EE_PFnHadIsoSlope*pt); };
float EE_GetPFnHadIsoCut_M( float pt ){ return (1.5 + EE_PFnHadIsoSlope*pt); };

float EE_GetPFphoIsoCut_M( float pt ){ return (1.0 + 0.005*pt); };

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
