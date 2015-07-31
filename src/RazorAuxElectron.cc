#include "RazorAnalyzer.h"


bool RazorAnalyzer::isEGammaPOGVetoElectron(int i){
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

bool RazorAnalyzer::isEGammaPOGLooseElectron(int i){
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

bool RazorAnalyzer::isEGammaPOGTightElectron(int i){
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

bool RazorAnalyzer::isVetoElectron(int i){
  return isMVANonTrigVetoElectron(i);
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
	 && ele_miniiso[i] < 0.1
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
	 && ele_miniiso[i] < 0.1
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
	 && ele_miniiso[i] < 0.1
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
	 && ele_miniiso[i] < 0.1
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
      && ( (elePt[i] > 20 && ele_miniiso[i] < 0.2)
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
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = -0.143;
  if (subdet == 1 && ptBin == 0) MVACut = -0.519;
  if (subdet == 2 && ptBin == 0) MVACut = -0.180;
  if (subdet == 0 && ptBin == 1) MVACut = -0.304;
  if (subdet == 1 && ptBin == 1) MVACut = -0.526;
  if (subdet == 2 && ptBin == 1) MVACut = -0.693;

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


bool RazorAnalyzer::passEGammaPOGLooseElectronIso(int i){
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

bool RazorAnalyzer::passEGammaPOGTightElectronIso(int i){
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

bool RazorAnalyzer::passLooseElectronIso(int i){
  bool pass = false;
  if ( ele_miniiso[i] < 0.1 ){
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passTightElectronIso(int i){
  bool pass = false;
  if ( ele_miniiso[i] < 0.1 ){
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){
 
  bool pass = false;
  if (  ( (elePt[i] > 20 && ele_miniiso[i] < 0.2 )
	  ||
	  (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) < 5)
	  )
	) {
    pass = true;
  }

  return pass;

}


bool RazorAnalyzer::passRunOneHZZElectronPreselection(int i){

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


bool RazorAnalyzer::isRunOneHZZElectron(int i){

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

double RazorAnalyzer::getElectronEfficiencyRunOne( string selectionType, double pt, double eta) {

  double result = 0;
  if (selectionType == "tight") {
    if (fabs(eta) < 0.8) {
      if (pt >= 5 && pt < 15) {
	result =  0.413;
      } else if (pt >= 15 && pt < 20) {
	result =  0.561;
      } else if (pt >= 20 && pt < 30) {
	result =  0.703;
      } else if (pt >= 30 && pt < 40) {
	result =  0.797;
      } else if (pt >= 40 && pt < 50) {
	result =  0.844;
      } else if (pt >= 50) {
	result =  0.862;
      } 
    } else if ( fabs(eta) < 1.44) {
      if (pt >= 5 && pt < 15) {
	result =  0.410;
      } else if (pt >= 15 && pt < 20) {
	result =  0.523;
      } else if (pt >= 20 && pt < 30) {
	result =  0.652;
      } else if (pt >= 30 && pt < 40) {
	result =  0.773;
      } else if (pt >= 40 && pt < 50) {
	result =  0.840;
      } else if (pt >= 50) {
	result =  0.865;
      } 
    } else if ( fabs(eta) < 1.57) {
      if (pt >= 5 && pt < 15) {
	result =  0.233;
      } else if (pt >= 15 && pt < 20) {
	result =  0.279;
      } else if (pt >= 20 && pt < 30) {
	result =  0.366;
      } else if (pt >= 30 && pt < 40) {
	result =  0.521;
      } else if (pt >= 40 && pt < 50) {
	result =  0.686;
      } else if (pt >= 50) {
	result =  0.713;
      } 
    } else if ( fabs(eta) < 2.0) {
      if (pt >= 5 && pt < 15) {
	result =  0.180;
      } else if (pt >= 15 && pt < 20) {
	result =  0.290;
      } else if (pt >= 20 && pt < 30) {
	result =  0.462;
      } else if (pt >= 30 && pt < 40) {
	result =  0.589;
      } else if (pt >= 40 && pt < 50) {
	result =  0.667;
      } else if (pt >= 50) {
	result =  0.710;
      } 
    } else if ( fabs(eta) < 2.5) {
      if (pt >= 5 && pt < 15) {
	result =  0.188;
      } else if (pt >= 15 && pt < 20) {
	result =  0.305;
      } else if (pt >= 20 && pt < 30) {
	result =  0.458;
      } else if (pt >= 30 && pt < 40) {
	result =  0.564;
      } else if (pt >= 40 && pt < 50) {
	result =  0.638;
      } else if (pt >= 50) {
	result =  0.675;
      } 
    }
  }

  if (selectionType == "loose") {
    if (fabs(eta) < 0.8) {
      if (pt >= 5 && pt < 15) {
	result =  0.488;
      } else if (pt >= 15 && pt < 20) {
	result =  0.642;
      } else if (pt >= 20 && pt < 30) {
	result =  0.774;
      } else if (pt >= 30 && pt < 40) {
	result =  0.851;
      } else if (pt >= 40 && pt < 50) {
	result =  0.883;
      } else if (pt >= 50) {
	result =  0.894;
      } 
    } else if ( fabs(eta) < 1.44) {
      if (pt >= 5 && pt < 15) {
	result =  0.566;
      } else if (pt >= 15 && pt < 20) {
	result =  0.687;
      } else if (pt >= 20 && pt < 30) {
	result =  0.789;
      } else if (pt >= 30 && pt < 40) {
	result =  0.872;
      } else if (pt >= 40 && pt < 50) {
	result =  0.908;
      } else if (pt >= 50) {
	result =  0.917;
      } 
    } else if ( fabs(eta) < 1.57) {
      if (pt >= 5 && pt < 15) {
	result =  0.415;
      } else if (pt >= 15 && pt < 20) {
	result =  0.475;
      } else if (pt >= 20 && pt < 30) {
	result =  0.584;
      } else if (pt >= 30 && pt < 40) {
	result =  0.737;
      } else if (pt >= 40 && pt < 50) {
	result =  0.846;
      } else if (pt >= 50) {
	result =  0.852;
      } 
    } else if ( fabs(eta) < 2.0) {
      if (pt >= 5 && pt < 15) {
	result =  0.306;
      } else if (pt >= 15 && pt < 20) {
	result =  0.485;
      } else if (pt >= 20 && pt < 30) {
	result =  0.680;
      } else if (pt >= 30 && pt < 40) {
	result =  0.783;
      } else if (pt >= 40 && pt < 50) {
	result =  0.827;
      } else if (pt >= 50) {
	result =  0.838;
      } 
    } else if ( fabs(eta) < 2.5) {
      if (pt >= 5 && pt < 15) {
	result =  0.315;
      } else if (pt >= 15 && pt < 20) {
	result =  0.471;
      } else if (pt >= 20 && pt < 30) {
	result =  0.632;
      } else if (pt >= 30 && pt < 40) {
	result =  0.718;
      } else if (pt >= 40 && pt < 50) {
	result =  0.767;
      } else if (pt >= 50) {
	result =  0.780;
      } 
    }
  } //end if loose

  return result;

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
