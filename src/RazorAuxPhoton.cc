#include "RazorAnalyzer.h"

bool RazorAnalyzer::photonPassesElectronVeto(int i){
    //use presence of a pixel seed as proxy for an electron veto
    return !(pho_hasPixelSeed[i]);
}

//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;
    if(fabs(pho_superClusterEta[i]) < 1.0){
        effAreaChargedHadrons = 0.0130;
        effAreaNeutralHadrons = 0.0056;
        effAreaPhotons = 0.0896;
    }
    else if(fabs(pho_superClusterEta[i]) < 1.479){
        effAreaChargedHadrons = 0.0096;
        effAreaNeutralHadrons = 0.0107;
        effAreaPhotons = 0.0762;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.0){
        effAreaChargedHadrons = 0.0107;
        effAreaNeutralHadrons = 0.0019;
        effAreaPhotons = 0.0383;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.2){
        effAreaChargedHadrons = 0.0077;
        effAreaNeutralHadrons = 0.0011;
        effAreaPhotons = 0.0534;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.3){
        effAreaChargedHadrons = 0.0088;
        effAreaNeutralHadrons = 0.0077;
        effAreaPhotons = 0.0846;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.4){
        effAreaChargedHadrons = 0.0065;
        effAreaNeutralHadrons = 0.0178;
        effAreaPhotons = 0.1032;
    }
    else{
        effAreaChargedHadrons = 0.0030;
        effAreaNeutralHadrons = 0.1675;
        effAreaPhotons = 0.1598;
    }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    
    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;
    
    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = max(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::photonPassesRunOneIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;
    if(fabs(pho_superClusterEta[i]) < 1.0){
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.030;
        effAreaPhotons = 0.148;
    }
    else if(fabs(pho_superClusterEta[i]) < 1.479){
        effAreaChargedHadrons = 0.010;
        effAreaNeutralHadrons = 0.057;
        effAreaPhotons = 0.130;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.0){
        effAreaChargedHadrons = 0.014;
        effAreaNeutralHadrons = 0.039;
        effAreaPhotons = 0.112;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.2){
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.015;
        effAreaPhotons = 0.216;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.3){
        effAreaChargedHadrons = 0.016;
        effAreaNeutralHadrons = 0.024;
        effAreaPhotons = 0.262;
    }
    else if(fabs(pho_superClusterEta[i]) < 2.4){
        effAreaChargedHadrons = 0.020;
        effAreaNeutralHadrons = 0.039;
        effAreaPhotons = 0.260;
    }
    else{
        effAreaChargedHadrons = 0.012;
        effAreaNeutralHadrons = 0.072;
        effAreaPhotons = 0.266;
    }
    cout << "Photon areas: " << effAreaChargedHadrons << " " << effAreaNeutralHadrons << " " << effAreaPhotons << endl;

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = max(pho_sumChargedHadronPt[i] - fixedGridRhoAll*effAreaChargedHadrons, 0.);
    cout << "charged " << PFIsoCorrected_ChHad << endl;
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;
    cout << "passed charged" << endl;
    
    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = max(pho_sumNeutralHadronEt[i] - fixedGridRhoAll*effAreaNeutralHadrons, 0.);
    cout << "neutral " << PFIsoCorrected_NeuHad << endl;
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;
    cout << "passed neutral had" << endl;
    
    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = max(pho_sumPhotonEt[i] - fixedGridRhoAll*effAreaPhotons, 0.);
    cout << "pho " << PFIsoCorrected_Photons << endl;
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;
    cout << "passed photon" << endl;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoFull5x5SigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //Isolation
    if(!photonPassesIsolation(i, PFChHadIsoCut, PFNeuHadIsoCut, PFPhotIsoCut)) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesRunOneCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoSigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //Isolation
    if(!photonPassesRunOneIsolation(i, PFChHadIsoCut, PFNeuHadIsoCut, PFPhotIsoCut)) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::passesCutsBasedPhotonIDNoIsoCuts(int i, double HoverECut, double SigmaIetaIetaCut){
    //electron veto
    if(!photonPassesElectronVeto(i)) return false;

    //HoverE
    if(pho_HoverE[i] > HoverECut) return false;
    
    //SigmaIetaIeta
    if(phoFull5x5SigmaIetaIeta[i] > SigmaIetaIetaCut) return false;

    //photon passed all cuts
    return true;
}

bool RazorAnalyzer::isLoosePhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.048, 0.0106, 2.56, 3.74+0.0025*phoPt[i], 2.68+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.069, 0.0266, 3.12, 17.11+0.0118*phoPt[i], 2.70+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isMediumPhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.032, 0.0101, 1.90, 2.96+0.0025*phoPt[i], 1.39+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.0166, 0.0264, 1.95, 4.42+0.0118*phoPt[i], 1.89+0.0059*phoPt[i]);
}

//photon ID with electron veto and HoverE and sigmaIetaIeta cuts only
bool RazorAnalyzer::isMediumPhotonNoIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonIDNoIsoCuts(i, 0.032, 0.0101);
    }
    //endcap photons
    return passesCutsBasedPhotonIDNoIsoCuts(i, 0.0166,0.0264);
}

bool RazorAnalyzer::photonPassesLooseIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 2.56, 3.74+0.0025*phoPt[i], 2.68+0.001*phoPt[i] );
    }
    //endcap photons
    return photonPassesIsolation(i, 3.12, 17.11+0.0118*phoPt[i], 2.70+0.0059*phoPt[i]);
}

bool RazorAnalyzer::photonPassesMediumIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 1.90, 2.96+0.0025*phoPt[i], 1.39+0.001*phoPt[i]);
    }
    //endcap photons
    return photonPassesIsolation(i, 1.95, 4.42+0.0118*phoPt[i], 1.89+0.0059*phoPt[i]);
}

bool RazorAnalyzer::photonPassesTightIsoCuts(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation(i, 1.86, 2.64+0.0025*phoPt[i], 1.20+0.001*phoPt[i]);
    }
    //endcap photons
    return photonPassesIsolation(i, 1.68, 4.42+0.0118*phoPt[i], 1.03+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isTightPhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesCutsBasedPhotonID(i, 0.011, 0.0099, 1.86, 2.64+0.0025*phoPt[i], 1.20+0.001*phoPt[i]);
    }
    //endcap photons
    return passesCutsBasedPhotonID(i, 0.015, 0.0263, 1.68, 4.42+0.0118*phoPt[i], 1.03+0.0059*phoPt[i]);
}

bool RazorAnalyzer::isMediumRunOnePhoton(int i){
    //barrel photons
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return passesRunOneCutsBasedPhotonID(i, 0.05, 0.011, 1.5, 1.0+0.04*phoPt[i], 0.7*0.005*phoPt[i]);
    }
    //endcap photons
    return passesRunOneCutsBasedPhotonID(i, 0.05, 0.033, 1.2, 1.5+0.04*phoPt[i], 1.0+0.005*phoPt[i]);
}
