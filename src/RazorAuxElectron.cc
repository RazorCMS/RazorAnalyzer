#include "RazorAnalyzer.h"

bool RazorAnalyzer::isVetoElectron(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.02
	 && fabs(ele_dPhi[i]) < 0.2579
	 && eleFull5x5SigmaIetaIeta[i] < 0.0125
	 && ele_HoverE[i] < 0.2564
	 && fabs(ele_d0[i]) < 0.025
	 && fabs(ele_dZ[i]) < 0.5863
	 && ele_OneOverEminusOneOverP[i] < 0.1508
	 && ele_relIsoDBetaCorr[i] < 0.3313
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 2
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.0141
	 && fabs(ele_dPhi[i]) < 0.2591
	 && eleFull5x5SigmaIetaIeta[i] < 0.0371
	 && ele_HoverE[i] < 0.1335
	 && fabs(ele_d0[i]) < 0.2232
	 && fabs(ele_dZ[i]) < 0.9513
	 && ele_OneOverEminusOneOverP[i] < 0.1542
	 && ele_relIsoDBetaCorr[i] < 0.3816
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
    if ( fabs(ele_dEta[i]) < 0.0181
	 && fabs(ele_dPhi[i]) < 0.0936
	 && eleFull5x5SigmaIetaIeta[i] < 0.0123
	 && ele_HoverE[i] < 0.141
	 && fabs(ele_d0[i]) < 0.0166
	 && fabs(ele_dZ[i]) < 0.54342
	 && ele_OneOverEminusOneOverP[i] < 0.1353
	 && ele_relIsoDBetaCorr[i] < 0.24
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.0124
	 && fabs(ele_dPhi[i]) < 0.0642
	 && eleFull5x5SigmaIetaIeta[i] < 0.035
	 && ele_HoverE[i] < 0.1115
	 && fabs(ele_d0[i]) < 0.098
	 && fabs(ele_dZ[i]) < 0.9187
	 && ele_OneOverEminusOneOverP[i] < 0.1443
	 && ele_relIsoDBetaCorr[i] < 0.3529
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
    if ( fabs(ele_dEta[i]) < 0.0091
	 && fabs(ele_dPhi[i]) < 0.031
	 && eleFull5x5SigmaIetaIeta[i] < 0.0106
	 && ele_HoverE[i] < 0.0532
	 && fabs(ele_d0[i]) < 0.0126
	 && fabs(ele_dZ[i]) < 0.0116
	 && ele_OneOverEminusOneOverP[i] < 0.0609
	 && ele_relIsoDBetaCorr[i] < 0.1649
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.0106
	 && fabs(ele_dPhi[i]) < 0.0359
	 && eleFull5x5SigmaIetaIeta[i] < 0.0305
	 && ele_HoverE[i] < 0.0835
	 && fabs(ele_d0[i]) < 0.0163
	 && fabs(ele_dZ[i]) < 0.5999
	 && ele_OneOverEminusOneOverP[i] < 0.1126
	 && ele_relIsoDBetaCorr[i] < 0.2075
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
      && ( (elePt[i] > 20 && ele_relIsoDBetaCorr[i] < 0.3)
	   ||
	   (elePt[i] <= 20 && ele_relIsoDBetaCorr[i]*elePt[i] < 5)
	   )
      ) {
    pass = true;
  }


  return pass;

}




bool RazorAnalyzer::passLooseElectronID(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( fabs(ele_dEta[i]) < 0.0181
	 && fabs(ele_dPhi[i]) < 0.0936
	 && eleFull5x5SigmaIetaIeta[i] < 0.0123
	 && ele_HoverE[i] < 0.141
	 && fabs(ele_d0[i]) < 0.0166
	 && fabs(ele_dZ[i]) < 0.54342
	 && ele_OneOverEminusOneOverP[i] < 0.1353
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.0124
	 && fabs(ele_dPhi[i]) < 0.0642
	 && eleFull5x5SigmaIetaIeta[i] < 0.035
	 && ele_HoverE[i] < 0.1115
	 && fabs(ele_d0[i]) < 0.098
	 && fabs(ele_dZ[i]) < 0.9187
	 && ele_OneOverEminusOneOverP[i] < 0.1443
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
    if ( fabs(ele_dEta[i]) < 0.0091
	 && fabs(ele_dPhi[i]) < 0.031
	 && eleFull5x5SigmaIetaIeta[i] < 0.0106
	 && ele_HoverE[i] < 0.0532
	 && fabs(ele_d0[i]) < 0.0126
	 && fabs(ele_dZ[i]) < 0.0116
	 && ele_OneOverEminusOneOverP[i] < 0.0609
	 && ele_PassConvVeto[i]
	 && ele_MissHits[i] < 1
	) {
      pass = true;
    }
  } else {
    if (fabs(ele_dEta[i]) < 0.0106
	 && fabs(ele_dPhi[i]) < 0.0359
	 && eleFull5x5SigmaIetaIeta[i] < 0.0305
	 && ele_HoverE[i] < 0.0835
	 && fabs(ele_d0[i]) < 0.0163
	 && fabs(ele_dZ[i]) < 0.5999
	 && ele_OneOverEminusOneOverP[i] < 0.1126
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
    if ( ele_relIsoDBetaCorr[i] < 0.24	
	) {
      pass = true;
    }
  } else {
    if ( ele_relIsoDBetaCorr[i] < 0.3529
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passTightElectronIso(int i){
  bool pass = false;
  if(fabs(eleEta_SC[i]) < 1.479) {
    if ( ele_relIsoDBetaCorr[i] < 0.1649
	) {
      pass = true;
    }
  } else {
    if ( ele_relIsoDBetaCorr[i] < 0.2075
	) {
      pass = true;
    }
  } 
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  if (  ( (elePt[i] > 20 && ele_relIsoDBetaCorr[i] < 0.3)
	  ||
	  (elePt[i] <= 20 && ele_relIsoDBetaCorr[i]*elePt[i] < 5)
	  )
	) {
    pass = true;
  }

  return pass;

}
