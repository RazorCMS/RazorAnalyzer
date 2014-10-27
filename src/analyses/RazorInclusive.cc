#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes

using namespace std;

bool passesHadronicRazorBaseline(double MR, double Rsq);
bool passesLeptonicRazorBaseline(double MR, double Rsq);

void RazorAnalyzer::RazorInclusive()
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    TFile outFile("RazorInclusive.root", "RECREATE");
    map<string, TTree*> razorBoxes;

    vector<string> boxNames;
    boxNames.push_back("MuEle");
    boxNames.push_back("MuMu");
    boxNames.push_back("EleEle");
    boxNames.push_back("MuMultiJet");
    boxNames.push_back("MuJet");
    boxNames.push_back("EleMultiJet");
    boxNames.push_back("EleJet");
    boxNames.push_back("MultiJet");
    boxNames.push_back("2BJet");
    boxNames.push_back("1BJet");
    boxNames.push_back("0BJet");
    for(int i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //tree variables
    int nSelectedJets, nBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nSelectedTaus;
    float theMR;
    float theRsq;
    for(auto& box : razorBoxes){
        box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        box.second->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
        box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        box.second->Branch("nSelectedTaus", &nSelectedTaus, "nSelectedTaus/I");
        box.second->Branch("MR", &theMR, "MR/F");
        box.second->Branch("Rsq", &theRsq, "Rsq/F");
    }
    
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
                razorBoxes["MuEle"]->Fill();
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["MuMu"]->Fill();
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["EleEle"]->Fill();
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["MuMultiJet"]->Fill();

            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["MuJet"]->Fill();

            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["EleMultiJet"]->Fill();
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["EleJet"]->Fill();
            }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesHadronicRazorBaseline(theMR, theRsq)){  
                razorBoxes["MultiJet"]->Fill();
            }
        }
        //2BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 1){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["2BJet"]->Fill();
            }
        }
        //1BJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["1BJet"]->Fill();
            }
        }
        //0BJetBox
        else if(passedHadronicTrigger){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                razorBoxes["0BJet"]->Fill();
            }
        }
    }//end of event loop

    cout << "Writing output trees..." << endl;
    for(auto& box : razorBoxes) box.second->Write();

    outFile.Close();
}

bool passesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 400 || Rsq < 0.25) passes = false;
    if(MR < 550 && Rsq < 0.3) passes = false;
    return passes;
}

bool passesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 300 || Rsq < 0.15) passes = false;
    if(MR < 450 && Rsq < 0.2) passes = false;
    return passes;
}
