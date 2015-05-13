/*
Defines 2012 Cut Based Photon ID
Defines PhotonCandidate struct
*/

#ifndef RAZOR_AUX_PHOTON_HH
#define RAZOR_AUX_PHOTON_HH 1
//C++ INCLUDES
#include <iostream>
#include <math.h>
//ROOT INCLUDES
//LOCAL INCLUDES

constexpr bool  EB_EleVeto           = true;
constexpr float EB_HoverECut         = 0.05;
constexpr float EB_SigmaIetaIetaCut  = 0.012;
constexpr float EB_PFchHadIsoCut     = 2.6;
constexpr float EB_PFnHadIsoConst    = 3.5;
constexpr float EB_PFnHadIsoSlope    = 0.04;
constexpr float EB_PFphoIsoConst     = 1.3;
constexpr float EB_PFphoIsoSlope     = 0.005;

constexpr bool  EE_EleVeto           = true;
constexpr float EE_HoverECut         = 0.05;
constexpr float EE_SigmaIetaIetaCut  = 0.034;
constexpr float EE_PFchHadIsoCut     = 2.3;
constexpr float EE_PFnHadIsoConst    = 2.9;
constexpr float EE_PFnHadIsoSlope    = 0.04;

constexpr float EB_GetPFchHadIsoCut( ){ return EB_PFchHadIsoCut; };
float EB_GetPFnHadIsoCut( float pt ){ return (EB_PFnHadIsoConst + EB_PFnHadIsoSlope*pt); };
float EB_GetPFphoIsoCut( float pt ){ return (EB_PFphoIsoConst + EB_PFphoIsoSlope*pt); };

constexpr float EE_GetPFchHadIsoCut( ){ return EE_PFchHadIsoCut; };
float EE_GetPFnHadIsoCut( float pt ){ return (EE_PFnHadIsoConst + EE_PFnHadIsoSlope*pt); };

//------------------------------------------------------
// R u n 2   P h o t o n   I D   f o r   H g g R a z o r
//------------------------------------------------------

// E f f e c t i v e   A r e a   T h r e s h o l d s 
// -------------------------------------------------
constexpr int n_EA_Threshold = 6;
constexpr float EA_Threshold[] = { 1.0, 1.479, 2.0 , 2.2, 2.3, 2.4 };
enum class EA { chHad, nHad, photon};

// E f f e c t i v e   A r e a   V a l u e s 
// -----------------------------------------
float effArea[3][7];
void InitEffArea()
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
  //-------------------------
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
int GetEAbin( float eta )
{
  for( int i = 0; i < n_EA_Threshold; i++ )
    {
      if ( fabs( eta ) < EA_Threshold[i] ) return i;
    }
  return n_EA_Threshold;
};


// P h o t o n   I D   W o r k i n g   P o in t s
// ----------------------------------------------
enum class WP { Loose, Medium, Tight };
//EB
float HoverE_EB[] = { 0.0280, 0.0120, 0.0100 };
float SigmaIetaIeta_EB[] = { 0.01070, 0.0100, 0.0100 };
//EE
float HoverE_EE[] = { 0.093, 0.023, 0.015 };
float SigmaIetaIeta_EE[] = { 0.0272, 0.0267, 0.0265 };

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
};

#endif
