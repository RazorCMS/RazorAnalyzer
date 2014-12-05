#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

enum HggRazorBox {
    HighPt,
    Hbb,
    Zbb,
    HighRes,
    LowRes
};

void RazorAnalyzer::HggRazor(string outFileName, bool combineTrees)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    if (outFileName.empty()){
        cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << endl;
        outFileName = "HggRazor.root";
    }
    TFile outFile(outFileName.c_str(), "RECREATE");

    //one tree to hold all events
    TTree *razorTree = new TTree("HggRazor", "Info on selected razor inclusive events");

    //separate trees for individual boxes
    map<string, TTree*> razorBoxes;
    vector<string> boxNames;
    boxNames.push_back("HighPt");
    boxNames.push_back("Hbb");
    boxNames.push_back("Zbb");
    boxNames.push_back("HighRes");
    boxNames.push_back("LowRes");
    for(size_t i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nLooseBTaggedJets, nMediumBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
    float theMR;
    float theRsq;
    float mGammaGamma, pTGammaGamma;
    float mbbZ, mbbH;
    float sigmaEOverE1, sigmaEOverE2;
    HggRazorBox box;

    //set branches on big tree
    if(combineTrees){
        razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
        razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
        razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
        razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
        razorTree->Branch("mbbZ", &mbbZ, "mbbZ/F");
        razorTree->Branch("mbbH", &mbbH, "mbbH/F");
        razorTree->Branch("sigmaEOverE1", &sigmaEOverE1, "sigmaEOverE1/F");
        razorTree->Branch("sigmaEOverE2", &sigmaEOverE2, "sigmaEOverE2/F");
        razorTree->Branch("box", &box, "box/I");
    }
    //set branches on all trees
    else{ 
        for(auto& box : razorBoxes){
            box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
            box.second->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
            box.second->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
            box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
            box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
            box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
            box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
            box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
            box.second->Branch("MR", &theMR, "MR/F");
            box.second->Branch("Rsq", &theRsq, "Rsq/F");
            box.second->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
            box.second->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
            box.second->Branch("mbbZ", &mbbZ, "mbbZ/F");
            box.second->Branch("mbbH", &mbbH, "mbbH/F");
            box.second->Branch("sigmaEOverE1", &sigmaEOverE1, "sigmaEOverE1/F");
            box.second->Branch("sigmaEOverE2", &sigmaEOverE2, "sigmaEOverE2/F");
        }
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

        //fill normalization histogram
        NEvents->Fill(1.0);

        //reset tree variables
        nSelectedJets = 0;
        nLooseBTaggedJets = 0;
        nMediumBTaggedJets = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nTightTaus = 0;
        theMR = -1;
        theRsq = -1;
        mGammaGamma = -1;
        pTGammaGamma = -1;
        mbbZ = 0;
        mbbH = 0;
        sigmaEOverE1 = -1;
        sigmaEOverE2 = -1;
        if(combineTrees) box = LowRes;

        //TODO: triggers!
        bool passedDiphotonTrigger = true;

        if(!passedDiphotonTrigger) continue;

        //muon selection
        for(int i = 0; i < nMuons; i++){
            if(!isLooseMuon(i)) continue;  
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;

            nLooseMuons++;

            if(isTightMuon(i)){ 
                nTightMuons++;
            }
        }
        for(int i = 0; i < nElectrons; i++){
            if(!isLooseElectron(i)) continue; 
            if(elePt[i] < 10) continue;

            nLooseElectrons++;

            if(isTightElectron(i)){ 
                nTightElectrons++;
            }
        }
        for(int i = 0; i < nTaus; i++){
            if(!isTightTau(i)) continue; 

            nTightTaus++;
        }

        vector<TLorentzVector> GoodJets;
        for(int i = 0; i < nJets; i++){
            if(jetPt[i] < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //TODO: exclude selected photons from the jet collection

            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVL(i)){
                nLooseBTaggedJets++;
            }
            if(isCSVM(i)){ 
                nMediumBTaggedJets++;
            }
        }
        //TODO: if there are two loose bjets and one medium bjet, compute mbbH, mbbZ

        //TODO: photon selection

        //Compute the razor variables using the selected jets and the diphoton system
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        //TODO: add photons
        TLorentzVector PFMET = makeTLorentzVector(metPt, 0, metPhi, 0);

        vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
        theMR = computeMR(hemispheres[0], hemispheres[1]); 
        theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);

        //TODO: sort events into boxes
        //HighPt Box
        if(true){
            if(combineTrees){
                box = HighPt;
                razorTree->Fill();
            }
            else razorBoxes["HighPt"]->Fill();
        }
        //Hbb Box
        else if(true){
            if(combineTrees){
                box = Hbb;
                razorTree->Fill();
            }
            else razorBoxes["Hbb"]->Fill();
        }

        //TODO: everything else

    }//end of event loop

    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();

    outFile.Close();
}
