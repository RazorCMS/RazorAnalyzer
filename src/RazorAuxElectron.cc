#include "RazorAnalyzer.h"

float RazorAnalyzer::GetElectronScaleCorrection( double pt, double eta ) {  
  double scaleCorr = 1.0;
  // if ( pt > 0 && fabs(eta) < 1.5) {
  //   scaleCorr = 1.015;
  // } else {
  //   scaleCorr = 1.05;
  // }
  //some dummy code so that it compiler doesn't complain
  if (pt > 0 || eta > 0) {
    scaleCorr = 1.0;
  }
  return scaleCorr;
}

float RazorAnalyzer::GetElectronEffectiveAreaMean(int i, bool use25nsCuts ){ 

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (use25nsCuts) {
    if (fabs(eleEta_SC[i]) < 1.0) {
      effArea = 0.0960;
    } else if (fabs(eleEta_SC[i]) < 1.479) {
      effArea = 0.0947;	
    } else if (fabs(eleEta_SC[i]) < 2.0) {
      effArea = 0.0580;	
    } else if (fabs(eleEta_SC[i]) < 2.2) {
      effArea = 0.0688;	
    } else if (fabs(eleEta_SC[i]) < 2.3) {
      effArea = 0.0967;	
    } else if (fabs(eleEta_SC[i]) < 2.4) {
      effArea = 0.1195;	
    } else if (fabs(eleEta_SC[i]) < 2.5) {
      effArea = 0.1475;	
    } 
  } else {
    if (fabs(eleEta_SC[i]) < 0.8) {
      effArea = 0.0973;
    } else if (fabs(eleEta_SC[i]) < 1.3) {
      effArea = 0.0954;
    } else if (fabs(eleEta_SC[i]) < 2.0) {
      effArea = 0.0632;	
    } else if (fabs(eleEta_SC[i]) < 2.2) {
      effArea = 0.0727;	
    } else {
      effArea = 0.1337;	
    } 
  }
  return effArea;
}

float RazorAnalyzer::GetElectronEffectiveArea90(int i){ 

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (fabs(eleEta_SC[i]) < 1.0) {
    effArea = 0.1752;
  } else if (fabs(eleEta_SC[i]) < 1.479) {
    effArea = 0.1862;	
  } else if (fabs(eleEta_SC[i]) < 2.0) {
    effArea = 0.1411;	
  } else if (fabs(eleEta_SC[i]) < 2.2) {
    effArea = 0.1534;	
  } else if (fabs(eleEta_SC[i]) < 2.3) {
    effArea = 0.1903;	
  } else if (fabs(eleEta_SC[i]) < 2.4) {
    effArea = 0.2243;	
  } else if (fabs(eleEta_SC[i]) < 2.5) {
    effArea = 0.2687;	
  }
  return effArea;
}

bool RazorAnalyzer::isEGammaPOGVetoElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGVetoElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGVetoElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGLooseElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGMediumElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGTightElectronIso(i,use25nsCuts)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isVetoElectron(int i, bool applyID, bool applyIso){
  return isMVANonTrigVetoElectron(i, applyID, applyIso);
}

bool RazorAnalyzer::isLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMVANonTrigVetoElectron(int i, bool applyID, bool applyIso){

  bool pass = true;
  if (applyID) {
    if (!passMVANonTrigVetoElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passMVANonTrigVetoElectronIso(i)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGVetoElectronID(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0152
	   && fabs(ele_dPhi[i]) < 0.216
	   && eleFull5x5SigmaIetaIeta[i] < 0.0114
	   && ele_HoverE[i] < 0.181
	   && fabs(ele_d0[i]) < 0.0564
	   && fabs(ele_dZ[i]) < 0.472
	   && ele_OneOverEminusOneOverP[i] < 0.207
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.0113
	  && fabs(ele_dPhi[i]) < 0.237
	  && eleFull5x5SigmaIetaIeta[i] < 0.0352
	  && ele_HoverE[i] < 0.116
	  && fabs(ele_d0[i]) < 0.222
	  && fabs(ele_dZ[i]) < 0.921
	  && ele_OneOverEminusOneOverP[i] < 0.174
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 3
	  ) {
	pass = true;
      }
    } 
  } 

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0126
	   && fabs(ele_dPhi[i]) < 0.107
	   && eleFull5x5SigmaIetaIeta[i] < 0.012
	   && ele_HoverE[i] < 0.186
	   && fabs(ele_d0[i]) < 0.0621
	   && fabs(ele_dZ[i]) < 0.613
	   && ele_OneOverEminusOneOverP[i] < 0.239
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.0109
	  && fabs(ele_dPhi[i]) < 0.219
	  && eleFull5x5SigmaIetaIeta[i] < 0.0339
	  && ele_HoverE[i] < 0.0962
	  && fabs(ele_d0[i]) < 0.279
	  && fabs(ele_dZ[i]) < 0.947
	  && ele_OneOverEminusOneOverP[i] < 0.141
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 3
	  ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronID(int i, bool use25nsCuts){
  bool pass = false;
  
  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0105
	   && fabs(ele_dPhi[i]) < 0.115
	   && eleFull5x5SigmaIetaIeta[i] < 0.0103
	   && ele_HoverE[i] < 0.104
	   && fabs(ele_d0[i]) < 0.0261
	   && fabs(ele_dZ[i]) < 0.41
	   && ele_OneOverEminusOneOverP[i] < 0.102
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00814
	  && fabs(ele_dPhi[i]) < 0.182
	  && eleFull5x5SigmaIetaIeta[i] < 0.0301
	  && ele_HoverE[i] < 0.0897
	  && fabs(ele_d0[i]) < 0.118
	  && fabs(ele_dZ[i]) < 0.822
	  && ele_OneOverEminusOneOverP[i] < 0.126
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    } 
  }

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00976
	   && fabs(ele_dPhi[i]) < 0.0929
	   && eleFull5x5SigmaIetaIeta[i] < 0.0105
	   && ele_HoverE[i] < 0.0765
	   && fabs(ele_d0[i]) < 0.0227
	   && fabs(ele_dZ[i]) < 0.379
	   && ele_OneOverEminusOneOverP[i] < 0.184
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00952
	  && fabs(ele_dPhi[i]) < 0.181
	  && eleFull5x5SigmaIetaIeta[i] < 0.0318
	  && ele_HoverE[i] < 0.0824
	  && fabs(ele_d0[i]) < 0.242
	  && fabs(ele_dZ[i]) < 0.921
	  && ele_OneOverEminusOneOverP[i] < 0.125
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronID(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0103
	   && fabs(ele_dPhi[i]) < 0.0336
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0876
	   && fabs(ele_d0[i]) < 0.0118
	   && fabs(ele_dZ[i]) < 0.373
	   && ele_OneOverEminusOneOverP[i] < 0.0174
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00733
	  && fabs(ele_dPhi[i]) < 0.114
	  && eleFull5x5SigmaIetaIeta[i] < 0.0283
	  && ele_HoverE[i] < 0.0678
	  && fabs(ele_d0[i]) < 0.0739
	  && fabs(ele_dZ[i]) < 0.602
	  && ele_OneOverEminusOneOverP[i] < 0.0898
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  } 

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0094
	   && fabs(ele_dPhi[i]) < 0.0296
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0372
	   && fabs(ele_d0[i]) < 0.0151
	   && fabs(ele_dZ[i]) < 0.238
	   && ele_OneOverEminusOneOverP[i] < 0.118	 
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00773
	  && fabs(ele_dPhi[i]) < 0.148
	  && eleFull5x5SigmaIetaIeta[i] < 0.0287
	  && ele_HoverE[i] < 0.0546
	  && fabs(ele_d0[i]) < 0.0535
	  && fabs(ele_dZ[i]) < 0.572
	  && ele_OneOverEminusOneOverP[i] < 0.104
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronID(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00926
	   && fabs(ele_dPhi[i]) < 0.0336
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0597
	   && fabs(ele_d0[i]) < 0.0111	   
	   && fabs(ele_dZ[i]) < 0.466
	   && ele_OneOverEminusOneOverP[i] < 0.012
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00724
	  && fabs(ele_dPhi[i]) < 0.0918
	  && eleFull5x5SigmaIetaIeta[i] < 0.0279
	  && ele_HoverE[i] < 0.0615
	  && fabs(ele_d0[i]) < 0.0351
	  && fabs(ele_dZ[i]) < 0.417
	  && ele_OneOverEminusOneOverP[i] < 0.00999
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    } 
  } 

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0095
	   && fabs(ele_dPhi[i]) < 0.0291
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0372
	   && fabs(ele_d0[i]) < 0.0144
	   && fabs(ele_dZ[i]) < 0.323
	   && ele_OneOverEminusOneOverP[i] < 0.0174
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00762
	  && fabs(ele_dPhi[i]) < 0.0439
	  && eleFull5x5SigmaIetaIeta[i] < 0.0287
	  && ele_HoverE[i] < 0.0544
	  && fabs(ele_d0[i]) < 0.0377
	  && fabs(ele_dZ[i]) < 0.571
	  && ele_OneOverEminusOneOverP[i] < 0.01
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronID(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = -0.11;
  if (subdet == 1 && ptBin == 0) MVACut = -0.55;
  if (subdet == 2 && ptBin == 0) MVACut = -0.60;
  if (subdet == 0 && ptBin == 1) MVACut = -0.16;
  if (subdet == 1 && ptBin == 1) MVACut = -0.65;
  if (subdet == 2 && ptBin == 1) MVACut = -0.74;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut
      && fabs(ele_ip3dSignificance[i]) < 4
      ) {
    pass = true;
  }

  return pass;

}

bool RazorAnalyzer::passEGammaPOGVetoElectronIso(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveArea90(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.126
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveArea90(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.144
	   ) {
	pass = true;
      }
    }
  }

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.161
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.193
	   ) {
	pass = true;
      }
    }
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronIso(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
     if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveArea90(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0893
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveArea90(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.121
	   ) {
	pass = true;
      }
    } 
 } 

  //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.118
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.118
	   ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronIso(int i, bool use25nsCuts){
  bool pass = false;

  if(use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0766
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0678
	   ) {
	pass = true;
      }
    } 
  }

 //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0987
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0902
	   ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronIso(int i, bool use25nsCuts){
  bool pass = false;

  if (use25nsCuts) {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0354
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0646
	   ) {
	pass = true;
      }
    } 
  } 

 //50ns cuts below
  else {
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0468
	   ) {
	pass = true;
      }
    } else {
      if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0759
	   ) {
	pass = true;
      }
    } 
  }
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (  ( (elePt[i] > 20 && (ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.2 )
	  ||
	   (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) < 5)
	  )
	) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passHZZElectronPreselection(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;
  
  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut   
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      ) {
    pass = true;
  }

  return pass;
}


bool RazorAnalyzer::isHZZElectron(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;
  
  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut   
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      && fabs(ele_d0[i]) < 0.05
      && ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.4)      
      ) {
    pass = true;
  }

  return pass;
}




bool RazorAnalyzer::matchTagElectronHLTFilters(int i){
  bool match = false;
  if ( 
      //Data filters
      ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13] 
      || ele_passHLTFilter[i][49] || ele_passHLTFilter[i][53] || ele_passHLTFilter[i][57] || ele_passHLTFilter[i][60]
      //MC filters
      || ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
       ) {
    match = true;
  }

  return match;
}

bool RazorAnalyzer::matchProbeElectronHLTFilters(int i){
  bool match = false;
  if ( 
      ele_passHLTFilter[i][50] || ele_passHLTFilter[i][51]     
      || ele_passHLTFilter[i][61]     || ele_passHLTFilter[i][62]     
       ) {
    match = true;
  }
  
  return match;
}

bool RazorAnalyzer::matchProbeSCHLTFilters(int i){
  bool match = false;
  if ( 
      ele_passHLTFilter[i][54] || ele_passHLTFilter[i][55] 
      || ele_passHLTFilter[i][58]    || ele_passHLTFilter[i][59]    
       ) {
    match = true;
  }
  
  return match;
}

bool RazorAnalyzer::matchElectronHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleElectron") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13] 
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele23Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][1] 
	//MC filters
	|| ele_passHLTFilter[i][3]
	 ) {
      match = true;
    }
  }
  if (HLTFilter == "Ele27Loose") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][5] 
	//MC filters
	|| ele_passHLTFilter[i][8]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele27Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][6] 
	//MC filters
	|| ele_passHLTFilter[i][10]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele32Tight") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][13] 
	//MC filters
	|| ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "Ele105") {
    if ( ele_passHLTFilter[i][17] ) match = true;
  }
  if (HLTFilter == "Ele115") {
    if ( ele_passHLTFilter[i][19] ) match = true;
  }
   
  if (HLTFilter == "DoubleElectronLeg1") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][24] || ele_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][25] || ele_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }
   
  if (HLTFilter == "DoubleElectronLeg2DZ") {
    if ( 
	//Data filters
	ele_passHLTFilter[i][26] || ele_passHLTFilter[i][30]
	 ) {
      match = true;
    }
  }
   
  return match;
}
