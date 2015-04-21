/*
Defines 2012 Cut Based Photon ID
Defines PhotonCandidate struct
*/

#ifndef RAZOR_AUX_PHOTON_HH
#define RAZOR_AUX_PHOTON_HH 1
//C++ INCLUDES
#include <iostream>
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

/*
struct PhotonCandidate
{
  float phoPt;
  float phoEta;
  float phoPho;
  float phoSigmaIetaIeta;
  float phoR9;
  float phoHoverE;
  float phosumChargedHadronPt;
  float phosumNeutralHadronEt;
  float phosumPhotonEt;
  float phopassEleVeto;
};
*/
#endif
