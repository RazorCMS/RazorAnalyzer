#include "RazorAnalyzer.h"

bool RazorAnalyzer::isVetoElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.016315
	 && fabs(ele_dPhi[i]) < 0.252044
	 && eleFull5x5SigmaIetaIeta[i] < 0.011100
	 && ele_HoverE[i] < 0.345843
	 && fabs(ele_d0[i]) < 0.060279
	 && fabs(ele_dZ[i]) < 0.800538
	 && ele_OneOverEminusOneOverP[i] < 0.248070
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.164369
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 2
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.010671
	 && fabs(ele_dPhi[i]) < 0.245263
	 && eleFull5x5SigmaIetaIeta[i] < 0.033987
	 && ele_HoverE[i] < 0.134691
	 && fabs(ele_d0[i]) < 0.273097
	 && fabs(ele_dZ[i]) < 0.885860
	 && ele_OneOverEminusOneOverP[i] < 0.157160
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.212604
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] < 3
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::isLooseElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.012442
	 && fabs(ele_dPhi[i]) < 0.072624
	 && eleFull5x5SigmaIetaIeta[i] < 0.010557
	 && ele_HoverE[i] < 0.121476
	 && fabs(ele_d0[i]) < 0.022664
	 && fabs(ele_dZ[i]) < 0.173670
	 && ele_OneOverEminusOneOverP[i] < 0.221803
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.120026
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.010654
	 && fabs(ele_dPhi[i]) < 0.145129
	 && eleFull5x5SigmaIetaIeta[i] < 0.032602
	 && ele_HoverE[i] < 0.131862
	 && fabs(ele_d0[i]) < 0.097358
	 && fabs(ele_dZ[i]) < 0.198444
	 && ele_OneOverEminusOneOverP[i] < 0.142283
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.162914
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::isTightElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.006574
	 && fabs(ele_dPhi[i]) < 0.022868
	 && eleFull5x5SigmaIetaIeta[i] < 0.010181
	 && ele_HoverE[i] < 0.037553
	 && fabs(ele_d0[i]) < 0.009924
	 && fabs(ele_dZ[i]) < 0.015310
	 && ele_OneOverEminusOneOverP[i] < 0.131191
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.074355
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.005681
	 && fabs(ele_dPhi[i]) < 0.032046
	 && eleFull5x5SigmaIetaIeta[i] < 0.028766
	 && ele_HoverE[i] < 0.081902
	 && fabs(ele_d0[i]) < 0.027261
	 && fabs(ele_dZ[i]) < 0.147154
	 && ele_OneOverEminusOneOverP[i] < 0.106055
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.090185
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::isMVANonTrigVetoElectron(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.0;
  if (subdet == 1 && ptBin == 0) MVACut = 0.6;
  if (subdet == 0 && ptBin == 1) MVACut = -0.3;
  if (subdet == 1 && ptBin == 1) MVACut = 0.5;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut
      &&
      ( (fabs(eleEta_SC[i]) < 1.479 && fabs(ele_d0[i]) < 0.0166)
      	||
      	(fabs(eleEta_SC[i]) >= 1.479 && fabs(ele_d0[i]) < 0.098)
      	)
      && ( (elePt[i] > 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.3)
	   ||
	   (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i]))  < 5)
	   )
      ) {
    pass = true;
  }


  return pass;

}




bool RazorAnalyzer::passLooseElectronID(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.012442
	 && fabs(ele_dPhi[i]) < 0.072624
	 && eleFull5x5SigmaIetaIeta[i] < 0.010557
	 && ele_HoverE[i] < 0.121476
	 && fabs(ele_d0[i]) < 0.0221803
	 && fabs(ele_dZ[i]) < 0.173670
	 && ele_OneOverEminusOneOverP[i] < 0.221803
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.010654
	 && fabs(ele_dPhi[i]) < 0.145129
	 && eleFull5x5SigmaIetaIeta[i] < 0.032602
	 && ele_HoverE[i] < 0.131862
	 && fabs(ele_d0[i]) < 0.097358
	 && fabs(ele_dZ[i]) < 0.198444
	 && ele_OneOverEminusOneOverP[i] < 0.142283
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passTightElectronID(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.006574
	 && fabs(ele_dPhi[i]) < 0.022868
	 && eleFull5x5SigmaIetaIeta[i] < 0.010181
	 && ele_HoverE[i] < 0.037553
	 && fabs(ele_d0[i]) < 0.009924
	 && fabs(ele_dZ[i]) < 0.015310
	 && ele_OneOverEminusOneOverP[i] < 0.131191
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.005681
	 && fabs(ele_dPhi[i]) < 0.032046
	 && eleFull5x5SigmaIetaIeta[i] < 0.028766
	 && ele_HoverE[i] < 0.081902
	 && fabs(ele_d0[i]) < 0.027261
	 && fabs(ele_dZ[i]) < 0.147154
	 && ele_OneOverEminusOneOverP[i] < 0.106055
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronID(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.0;
  if (subdet == 1 && ptBin == 0) MVACut = 0.6;
  if (subdet == 0 && ptBin == 1) MVACut = -0.3;
  if (subdet == 1 && ptBin == 1) MVACut = 0.5;

  bool pass = false;
  if (ele_IDMVANonTrig[i] > MVACut
      &&
      ( (fabs(eleEta_SC[i]) < 1.479 && fabs(ele_d0[i]) < 0.0166)
	||
	(fabs(eleEta_SC[i]) >= 1.479 && fabs(ele_d0[i]) < 0.098)
	)
      ) {
    pass = true;
  }


  return pass;

}


bool RazorAnalyzer::passLooseElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.120026
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.162914
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passTightElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.074355
	) {
      pass = true;
    }
  } else {
    if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.090185
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  if (  ( (elePt[i] > 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.3)
	  ||
	  (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) < 5)
	  )
	) {
    pass = true;
  }

  return pass;

}


bool RazorAnalyzer::isRunOneLooseElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.007
	 && fabs(ele_dPhi[i]) < 0.15
	 && eleSigmaIetaIeta[i] < 0.01
	 && ele_HoverE[i] < 0.12
	 && fabs(ele_d0[i]) < 0.02
	 && fabs(ele_dZ[i]) < 0.2
	 && ele_OneOverEminusOneOverP[i] < 0.05
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.15
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] <= 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.009
	 && fabs(ele_dPhi[i]) < 0.10
	 && eleSigmaIetaIeta[i] < 0.03
	 && ele_HoverE[i] < 0.10
	 && fabs(ele_d0[i]) < 0.02
	 && fabs(ele_dZ[i]) < 0.2
	 && ele_OneOverEminusOneOverP[i] < 0.05
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.15
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] <= 1
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::isRunOneTightElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.004
	 && fabs(ele_dPhi[i]) < 0.03
	 && eleSigmaIetaIeta[i] < 0.01
	 && ele_HoverE[i] < 0.12
	 && fabs(ele_d0[i]) < 0.02
	 && fabs(ele_dZ[i]) < 0.2
	 && ele_OneOverEminusOneOverP[i] < 0.05
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.10
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] <= 0
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.005
	 && fabs(ele_dPhi[i]) < 0.02
	 && eleFull5x5SigmaIetaIeta[i] < 0.03
	 && ele_HoverE[i] < 0.10
	 && fabs(ele_d0[i]) < 0.02
	 && fabs(ele_dZ[i]) < 0.1
	 && ele_OneOverEminusOneOverP[i] < 0.05
	 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.10
	&& ele_PassConvVeto[i]
	 && ele_MissHits[i] <= 0
	) {
      pass = true;
    }
  } 
  return pass;
}
