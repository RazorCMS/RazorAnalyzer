#include "RazorAnalyzer.h"

float RazorAnalyzer::GetEffectiveAreaMean(int i){ 

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
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
  return effArea;
}

float RazorAnalyzer::GetEffectiveArea90(int i){ 

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

bool RazorAnalyzer::isEGammaPOGVetoElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGVetoElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGVetoElectronIso(i)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGLooseElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGLooseElectronIso(i)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGMediumElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGMediumElectronIso(i)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isEGammaPOGTightElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGTightElectronIso(i)) pass = false;
  } 
  return pass;
}

bool RazorAnalyzer::isVetoElectron(int i, bool applyID, bool applyIso){
  return isMVANonTrigVetoElectron(i, applyID, applyIso);
}

bool RazorAnalyzer::isLooseElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!(ele_miniiso[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMediumElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!(ele_miniiso[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isTightElectron(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i)) pass = false;
  }
  if (applyIso) {
    if (!(ele_miniiso[i] < 0.1)) pass = false;
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

bool RazorAnalyzer::passEGammaPOGVetoElectronID(int i){
  bool pass = false;
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
  return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronID(int i){
  bool pass = false;
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
  return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronID(int i){
  bool pass = false;
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
  return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronID(int i){
  bool pass = false;
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
  if (subdet == 0 && ptBin == 0) MVACut = -0.1;
  if (subdet == 1 && ptBin == 0) MVACut = -0.75;
  if (subdet == 2 && ptBin == 0) MVACut = -0.1;
  if (subdet == 0 && ptBin == 1) MVACut = -0.5;
  if (subdet == 1 && ptBin == 1) MVACut = -0.8;
  if (subdet == 2 && ptBin == 1) MVACut = -0.3;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut
      &&
      ( (fabs(eleEta_SC[i]) < 1.479 && fabs(ele_d0[i]) < 0.02)
	||
	(fabs(eleEta_SC[i]) >= 1.479 && fabs(ele_d0[i]) < 0.1)
	)
      ) {
    pass = true;
  }

  return pass;

}

bool RazorAnalyzer::passEGammaPOGVetoElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.161
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.193
	) {
      pass = true;
    }
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.118
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.118
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0987
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0902
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0468
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0759
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  if (  ( (elePt[i] > 20 && ele_miniiso[i] < 0.1 )
	  ||
	  (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) < 5)
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
      || ele_passHLTFilter[i][49] || ele_passHLTFilter[i][53]
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
       ) {
    match = true;
  }
  
  return match;
}

bool RazorAnalyzer::matchProbeSCHLTFilters(int i){
  bool match = false;
  if ( 
      ele_passHLTFilter[i][54] || ele_passHLTFilter[i][55]     
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
