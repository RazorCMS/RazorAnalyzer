#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

using namespace std;

void RazorAnalyzer::PhotonNtupler(string outputFilename)
{
    //create output file
    TFile outFile(outputFilename.c_str(), "RECREATE");
    TTree *PhotonNtuple = new TTree("PhotonNtuple", "Selected photon information");

    //photon branches carried over from ntuple
    PhotonNtuple->Branch("nPhotons", &nPhotons, "nPhotons/I");
    PhotonNtuple->Branch("phoPt", phoPt, "phoPt[nPhotons]/F");
    PhotonNtuple->Branch("phoEta", phoEta, "phoEta[nPhotons]/F");
    PhotonNtuple->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
    PhotonNtuple->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta,"phoFull5x5SigmaIetaIeta[nPhotons]/F");
    PhotonNtuple->Branch("phoR9", phoR9,"phoR9[nPhotons]/F");
    PhotonNtuple->Branch("pho_HoverE", pho_HoverE,"pho_HoverE[nPhotons]/F");
    PhotonNtuple->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt,"pho_sumChargedHadronPt[nPhotons]/F");
    PhotonNtuple->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt,"pho_sumNeutralHadronEt[nPhotons]/F");
    PhotonNtuple->Branch("pho_sumPhotonEt", pho_sumPhotonEt,"pho_sumPhotonEt[nPhotons]/F");
    PhotonNtuple->Branch("pho_isConversion", pho_isConversion,"pho_isConversion[nPhotons]/O");
    PhotonNtuple->Branch("pho_passEleVeto", pho_passEleVeto,"pho_passEleVeto[nPhotons]/O");
    PhotonNtuple->Branch("pho_RegressionE", pho_RegressionE,"pho_RegressionE[nPhotons]/F");
    PhotonNtuple->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty,"pho_RegressionEUncertainty[nPhotons]/F");
    PhotonNtuple->Branch("pho_IDMVA", pho_IDMVA,"pho_IDMVA[nPhotons]/F");
    PhotonNtuple->Branch("pho_superClusterEta", pho_superClusterEta,"pho_superClusterEta[nPhotons]/F");
    PhotonNtuple->Branch("pho_hasPixelSeed", pho_hasPixelSeed,"pho_hasPixelSeed[nPhotons]/O");
    PhotonNtuple->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
    PhotonNtuple->Branch("nPV", &nPV, "nPV/I");
    //extra photon branches
    int nGenPhotons;
    float genPhotonPt[200];
    float genPhotonEta[200];
    float genPhotonPhi[200];
    float genPhotonE[200];
    PhotonNtuple->Branch("nGenPhotons", &nGenPhotons,"nGenPhotons/I");
    PhotonNtuple->Branch("genPhotonPt", genPhotonPt, "genPhotonPt[nGenPhotons]/F");
    PhotonNtuple->Branch("genPhotonEta", genPhotonEta, "genPhotonEta[nGenPhotons]/F");
    PhotonNtuple->Branch("genPhotonPhi", genPhotonPhi, "genPhotonPhi[nGenPhotons]/F");
    PhotonNtuple->Branch("genPhotonE", genPhotonE, "genPhotonE[nGenPhotons]/F");

    bool phoMatchesGen[40];
    int phoMatchingGenPhotonIndex[40];
    float deltaEOverEBest[40];
    bool phoIsLoose[40];
    bool phoIsMedium[40];
    bool phoIsTight[40];
    bool phoIsBarrel[40];
    bool phoIsIsolatedLoose[40];
    bool phoIsIsolatedMedium[40];
    bool phoIsIsolatedTight[40];
    float phoChHadIsolation[40];
    float phoNeuHadIsolation[40];
    float phoPhotIsolation[40];
    PhotonNtuple->Branch("phoMatchesGen", phoMatchesGen, "phoMatchesGen[nPhotons]/O");
    PhotonNtuple->Branch("phoMatchingGenPhotonIndex", phoMatchingGenPhotonIndex, "phoMatchingGenPhotonIndex[nPhotons]/I");
    PhotonNtuple->Branch("deltaEOverEBest", deltaEOverEBest, "deltaEOverEBest[nPhotons]/F");
    PhotonNtuple->Branch("phoIsLoose", phoIsLoose, "phoIsLoose[nPhotons]/O");
    PhotonNtuple->Branch("phoIsMedium", phoIsMedium, "phoIsMedium[nPhotons]/O");
    PhotonNtuple->Branch("phoIsTight", phoIsTight, "phoIsTight[nPhotons]/O");
    PhotonNtuple->Branch("phoIsBarrel", phoIsBarrel, "phoIsBarrel[nPhotons]/O");
    PhotonNtuple->Branch("phoIsIsolatedLoose", phoIsIsolatedLoose, "phoIsIsolatedLoose[nPhotons]/O");
    PhotonNtuple->Branch("phoIsIsolatedMedium", phoIsIsolatedMedium, "phoIsIsolatedMedium[nPhotons]/O");
    PhotonNtuple->Branch("phoIsIsolatedTight", phoIsIsolatedTight, "phoIsIsolatedTight[nPhotons]/O");
    PhotonNtuple->Branch("phoChHadIsolation", phoChHadIsolation, "phoChHadIsolation[nPhotons]/F");
    PhotonNtuple->Branch("phoNeuHadIsolation", phoNeuHadIsolation, "phoNeuHadIsolation[nPhotons]/F");
    PhotonNtuple->Branch("phoPhotIsolation", phoPhotIsolation, "phoPhotIsolation[nPhotons]/F");

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //reset tree variables
        for(int i = 0; i < 40; i++){
            phoChHadIsolation[i] = -1;
            phoNeuHadIsolation[i] = -1;
            phoPhotIsolation[i] = -1;
        }
        nGenPhotons = 0;
        for(int i = 0; i < 200; i++){
            genPhotonPt[i] = -999;
            genPhotonEta[i] = -999;
            genPhotonPhi[i] = -999;
        }

        //select gen photons
        for(int g = 0; g < nGenParticle; g++){
           if(gParticleStatus[g] != 1) continue;
           if(gParticleId[g] != 22) continue;
           if(gParticleE[g] < 1.) continue;

           //fill gen photon info
           genPhotonPt[nGenPhotons] = gParticlePt[g];
           genPhotonEta[nGenPhotons] = gParticleEta[g];
           genPhotonPhi[nGenPhotons] = gParticlePhi[g];
           genPhotonE[nGenPhotons] = gParticleE[g];
           nGenPhotons++;
        }

        //fill photon information
        for(int i = 0; i < nPhotons; i++){
            //try to match photon with a gen photon
            phoMatchesGen[i] = false;
            phoMatchingGenPhotonIndex[i] = -1;
            deltaEOverEBest[i] = 999;
            for(int g = 0; g < nGenPhotons; g++){
                if(deltaR(phoEta[i], phoPhi[i], genPhotonEta[g], genPhotonPhi[g]) > 0.2) continue;
                float deltaEOverE = fabs(pho_RegressionE[i] - genPhotonE[g])/genPhotonE[g];
                if(deltaEOverE > 1.) continue;

                phoMatchesGen[i] = true;

                if(deltaEOverE < deltaEOverEBest[i]){
                    deltaEOverEBest[i] = deltaEOverE;
                    phoMatchingGenPhotonIndex[i] = g;
                }
            }

            phoIsLoose[i] = isLoosePhoton(i);
            phoIsMedium[i] = isMediumPhoton(i);
            phoIsTight[i] = isTightPhoton(i);
            if(fabs(pho_superClusterEta[i]) < 1.479) phoIsBarrel[i] = true;
            else phoIsBarrel[i] = false;

            phoIsIsolatedLoose[i] = photonPassesLooseIsoCuts(i);
            phoIsIsolatedMedium[i] = photonPassesMediumIsoCuts(i);
            phoIsIsolatedTight[i] = photonPassesTightIsoCuts(i);

            //get effective area for isolation calculations
            double effAreaChargedHadrons = 0.0;
            double effAreaNeutralHadrons = 0.0;
            double effAreaPhotons = 0.0;
            if(fabs(pho_superClusterEta[i]) < 1.0){
                effAreaChargedHadrons = 0.0080;
                effAreaNeutralHadrons = 0.0126;
                effAreaPhotons = 0.0982;
            }
            else if(fabs(pho_superClusterEta[i]) < 1.479){
                effAreaChargedHadrons = 0.0079;
                effAreaNeutralHadrons = 0.0237;
                effAreaPhotons = 0.0857;
            }
            else if(fabs(pho_superClusterEta[i]) < 2.0){
                effAreaChargedHadrons = 0.0080;
                effAreaNeutralHadrons = 0.0;
                effAreaPhotons = 0.0484;
            }
            else if(fabs(pho_superClusterEta[i]) < 2.2){
                effAreaChargedHadrons = 0.0048;
                effAreaNeutralHadrons = 0.0;
                effAreaPhotons = 0.0668;
            }
            else if(fabs(pho_superClusterEta[i]) < 2.3){
                effAreaChargedHadrons = 0.0029;
                effAreaNeutralHadrons = 0.0;
                effAreaPhotons = 0.0868;
            }
            else if(fabs(pho_superClusterEta[i]) < 2.4){
                effAreaChargedHadrons = 0.0036;
                effAreaNeutralHadrons = 0.0;
                effAreaPhotons = 0.0982;
            }
            else{
                effAreaChargedHadrons = 0.0016;
                effAreaNeutralHadrons = 0.0769;
                effAreaPhotons = 0.1337;
            }

            //compute isolation quantities
            phoChHadIsolation[i] = max(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
            phoNeuHadIsolation[i] = max(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
            phoPhotIsolation[i] = max(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
        }
        PhotonNtuple->Fill();
    }
    PhotonNtuple->Write();
    outFile.Close();
}
