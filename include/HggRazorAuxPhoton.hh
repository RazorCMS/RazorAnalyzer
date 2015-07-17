/*
Defines 2012 Cut Based Photon ID
Defines PhotonCandidate struct
*/

#ifndef HGGRAZOR_AUX_PHOTON_HH
#define HGGRAZOR_AUX_PHOTON_HH 1
//C++ INCLUDES
#include <iostream>
#include <math.h>
#include <map>
#include <vector>
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

float EB_GetPFchHadIsoCut( );
float EB_GetPFchHadIsoCut_VL( );
float EB_GetPFchHadIsoCut_M( );

float EB_GetPFnHadIsoCut( float pt );
float EB_GetPFnHadIsoCut_VL( float pt );
float EB_GetPFnHadIsoCut_M( float pt );

float EB_GetPFphoIsoCut( float pt );
float EB_GetPFphoIsoCut_VL( float pt );
float EB_GetPFphoIsoCut_M( float pt );

float EE_GetPFchHadIsoCut( );
float EE_GetPFchHadIsoCut_VL( );
float EE_GetPFchHadIsoCut_M( );

float EE_GetPFnHadIsoCut( float pt );
float EE_GetPFnHadIsoCut_VL( float pt );
float EE_GetPFnHadIsoCut_M( float pt );

float EE_GetPFphoIsoCut_M( float pt );

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
//std::map< EA, std::vector<float> > effAreas;
//float effArea[3][7] = {{}};
void InitEffArea();

// P h o t o n   I D   W o r k i n g   P o in t s
// ----------------------------------------------
enum class WP { VeryLoose, Loose, Medium, Tight };
//EB
constexpr float HoverE_EB[] = { 0.0280, 0.0120, 0.0100 };
constexpr float SigmaIetaIeta_EB[] = { 0.01070, 0.0100, 0.0100 };
//EE
constexpr float HoverE_EE[] = { 0.093, 0.023, 0.015 };
constexpr float SigmaIetaIeta_EE[] = { 0.0272, 0.0267, 0.0265 };

// I s o l a t i o n s
//--------------------
float GetchHadIsoCut( WP wp, float eta );
float GetnHadIsoCut( WP wp, float eta, float pt );
float GetPhotonIsoCut( WP wp, float eta, float pt );

#endif
