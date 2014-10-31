#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes

using namespace std;

enum RazorBox {
    MUELE, 
    MUMU,
    ELEELE,
    MUMULTIJET,
    MUJET,
    ELEMULTIJET,
    ELEJET,
    MULTIJET,
    TWOBJET,
    ONEBJET,
    ZEROBJET,
    NONE
};

bool passesHadronicRazorBaseline(double MR, double Rsq);
bool passesLeptonicRazorBaseline(double MR, double Rsq);

void RazorAnalyzer::RazorInclusive()
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    TFile outFile("RazorInclusive.root", "RECREATE");
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

    //tree variables
    int nSelectedJets, nBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nSelectedTaus;
    float theMR;
    float theRsq;
    RazorBox box;
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nSelectedTaus", &nSelectedTaus, "nSelectedTaus/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("box", &box, "box/I");

    //begin loop
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        //begin event
        if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //reset tree variables
        nSelectedJets = 0;
        nBTaggedJets = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nSelectedTaus = 0;
        theMR = -1;
        theRsq = -1;
        box = NONE;

        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        
        vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
        for(int i = 0; i < nMuons; i++){
            if(!isLooseMuon(i)) continue;  
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;

            nLooseMuons++;
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
            GoodLeptons.push_back(thisMuon);

            if(isTightMuon(i)){ 
                nTightMuons++;
            }
        }
        for(int i = 0; i < nElectrons; i++){
            if(!isLooseElectron(i)) continue; 
            if(elePt[i] < 10) continue;

            nLooseElectrons++;
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            GoodLeptons.push_back(thisElectron);

            if(isTightElectron(i)){ 
                nTightElectrons++;
            }
        }
        for(int i = 0; i < nTaus; i++){
            if(!isSelectedTau(i)) continue; 

            nSelectedTaus++;
        }
        
        vector<TLorentzVector> GoodJets;
        int numJetsAbove80GeV = 0;
        for(int i = 0; i < nJets; i++){
            if(jetPt[i] < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //exclude selected muons and electrons from the jet collection
            double deltaR = -1;
            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            for(auto& lep : GoodLeptons){
                double thisDR = thisJet.DeltaR(lep);
                if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
            
            if(jetPt[i] > 80) numJetsAbove80GeV++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVM(i)){ 
                nBTaggedJets++;
            }
        }
        if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);

        vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
        theMR = computeMR(hemispheres[0], hemispheres[1]); 
        theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);

        //MuEle Box
        if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = MUELE;
                razorTree->Fill();
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = MUMU;
                razorTree->Fill();
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = ELEELE;
                razorTree->Fill();
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = MUMULTIJET;
                razorTree->Fill();
            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = MUJET;
                razorTree->Fill();
            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = ELEMULTIJET;
                razorTree->Fill();
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                box = ELEJET;
                razorTree->Fill();
            }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesHadronicRazorBaseline(theMR, theRsq)){  
                box = MULTIJET;
                razorTree->Fill();
            }
        }
        //2BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 1){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                box = TWOBJET;
                razorTree->Fill();
            }
        }
        //1BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                box = ONEBJET;
                razorTree->Fill();
            }
        }
        //0BJetBox
        else if(passedHadronicTrigger){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                box = ZEROBJET;
                razorTree->Fill();
            }
        }
    }//end of event loop

    cout << "Writing output trees..." << endl;
    razorTree->Write();

    outFile.Close();
}

bool passesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool passesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 300 || Rsq < 0.15) passes = false;
    if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}
