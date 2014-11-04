#include "RazorAnalyzer.h"

bool RazorAnalyzer::isLooseElectron(int i){
    return true;
}

bool RazorAnalyzer::isTightElectron(int i){
  bool pass = false;
  if(eleEta[i] < 1.479) {
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


bool RazorAnalyzer::isVetoElectron(int i){
  bool pass = false;
  if(eleEta[i] < 1.479) {
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
