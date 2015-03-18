//Produces trees for estimating lepton selection efficiency (using gen-level information) and for studying the use of DY+Jets, W+Jets, and G+Jets control regions to model the distribution of the Z->invisible background in the razor analysis.

#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};

void RazorAnalyzer::RazorPhotonStudy(string outputfilename, bool isData)
{
    //****************************************************//
    //            Initialization of the tree              //
    //****************************************************//

    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorPhotonStudy.root";
    TFile outFile(outfilename.c_str(), "RECREATE");

    //one tree to hold all events
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nBTaggedJets;
    int nVetoMuons, nLooseMuons, nTightMuons;
    int nVetoElectrons, nLooseElectrons, nTightElectrons;
    int nLooseTaus, nMediumTaus, nTightTaus;
    int nGenMuons, nGenElectrons, nGenTauMuons, nGenTauElectrons, nGenTaus, nGenPhotons;
    float theMR, MR_noPho, MR_noZ, MR_noW, MR_noGenZ;
    float theRsq, Rsq_noPho, Rsq_noZ, Rsq_noW, Rsq_noGenZ; 
    float met,genmet,met_noPho, met_noZ, met_noW, genZmass, met_noGenZ;
    float leadingGenMuonPt, leadingGenPhotonPt;
    float leadingGenMuonEta, leadingGenPhotonEta;
    float leadingGenMuonPhi, leadingGenPhotonPhi;
    float leadingGenMuonE, leadingGenPhotonE;
    float subleadingGenMuonPt, subleadingGenPhotonPt;
    float subleadingGenMuonEta, subleadingGenPhotonEta;
    float subleadingGenMuonPhi, subleadingGenPhotonPhi;
    float subleadingGenMuonE, subleadingGenPhotonE;
    float leadingMuonPt, leadingPhotonPt;
    float leadingMuonEta, leadingPhotonEta;
    float leadingMuonPhi, leadingPhotonPhi;
    float leadingMuonE, leadingPhotonE;
    float subleadingMuonPt, subleadingPhotonPt;
    float subleadingMuonEta, subleadingPhotonEta;
    float subleadingMuonPhi, subleadingPhotonPhi;
    float subleadingMuonE, subleadingPhotonE;
    float metphi, genmetphi, metphi_noZ, metphi_noW, metphi_noPho, metphi_noGenZ;
    float HT, HT_noZ, HT_noW, HT_noPho, HT_noGenZ;
    int numJets, numJets_noZ, numJets_noW, numJets_noPho, numJets_noGenZ; 
    int numJets80, numJets80_noZ, numJets80_noW, numJets80_noPho, numJets80_noGenZ; 
    float genZpt, recoZpt, genZeta, recoZeta, genZphi, recoZphi, recoZmass, genWpt, recoWpt, genWeta, genWphi, recoWphi;
    bool leadGenMuonIsFound, leadGenPhotonIsFound;
    bool leadGenMuonIsFoundTight;
    float ptMatchingLeadGenMuon, ptMatchingLeadGenPhoton; //pt of the object matching the gen particle
    int nSelectedPhotons;    
    float mTLepMet;

    //set branches on big tree
    if(!isData){
        razorTree->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
        razorTree->Branch("nGenElectrons", &nGenElectrons, "nGenElectrons/I");
        razorTree->Branch("nGenTauMuons", &nGenTauMuons, "nGenTauMuons/I");
        razorTree->Branch("nGenTauElectrons", &nGenTauElectrons, "nGenTauElectrons/I");
        razorTree->Branch("nGenTaus", &nGenTaus, "nGenTaus/I");
        razorTree->Branch("nGenPhotons", &nGenPhotons, "nGenPhotons/I");
        razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
        razorTree->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
        razorTree->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
        razorTree->Branch("leadingGenMuonPhi", &leadingGenMuonPhi, "leadingGenMuonPhi/F");
        razorTree->Branch("leadingGenMuonE", &leadingGenMuonE, "leadingGenMuonE/F");
        razorTree->Branch("leadingGenPhotonPt", &leadingGenPhotonPt, "leadingGenPhotonPt/F");
        razorTree->Branch("leadingGenPhotonEta", &leadingGenPhotonEta, "leadingGenPhotonEta/F");
        razorTree->Branch("leadingGenPhotonPhi", &leadingGenPhotonPhi, "leadingGenPhotonPhi/F");
        razorTree->Branch("leadingGenPhotonE", &leadingGenPhotonE, "leadingGenPhotonE/F");
        razorTree->Branch("subleadingGenMuonPt", &subleadingGenMuonPt, "subleadingGenMuonPt/F");
        razorTree->Branch("subleadingGenMuonEta", &subleadingGenMuonEta, "subleadingGenMuonEta/F");
        razorTree->Branch("subleadingGenMuonPhi", &subleadingGenMuonPhi, "subleadingGenMuonPhi/F");
        razorTree->Branch("subleadingGenMuonE", &subleadingGenMuonE, "subleadingGenMuonE/F");
        razorTree->Branch("subleadingGenPhotonPt", &subleadingGenPhotonPt, "subleadingGenPhotonPt/F");
        razorTree->Branch("subleadingGenPhotonEta", &subleadingGenPhotonEta, "subleadingGenPhotonEta/F");
        razorTree->Branch("subleadingGenPhotonPhi", &subleadingGenPhotonPhi, "subleadingGenPhotonPhi/F");
        razorTree->Branch("subleadingGenPhotonE", &subleadingGenPhotonE, "subleadingGenPhotonE/F");
        razorTree->Branch("genZpt", &genZpt, "genZpt/F");
        razorTree->Branch("genZeta", &genZeta, "genZeta/F");
        razorTree->Branch("genZphi", &genZphi, "genZphi/F");
        razorTree->Branch("genZmass", &genZmass, "genZmass/F");
        razorTree->Branch("genWpt", &genWpt, "genWpt/F");
        razorTree->Branch("genWeta", &genWeta, "genWeta/F");
        razorTree->Branch("genWphi", &genWphi, "genWphi/F");
        razorTree->Branch("MR_noGenZ", &MR_noGenZ, "MR_noGenZ/F");
        razorTree->Branch("Rsq_noGenZ", &Rsq_noGenZ, "Rsq_noGenZ/F");
        razorTree->Branch("met_noGenZ", &met_noGenZ, "met_noGenZ/F");
        razorTree->Branch("metphi_noGenZ", &metphi_noGenZ, "metphi_noGenZ/F");
        razorTree->Branch("genmet", &genmet, "genmet/F");
        razorTree->Branch("genmetphi", &genmetphi, "genmetphi/F");
        razorTree->Branch("HT_noGenZ", &HT_noGenZ, "HT_noGenZ/F");
        razorTree->Branch("numJets_noGenZ", &numJets_noGenZ, "numJets_noGenZ/I");
        razorTree->Branch("numJets80_noGenZ", &numJets80_noGenZ, "numJets80_noGenZ/I");
        razorTree->Branch("leadGenMuonIsFound", &leadGenMuonIsFound, "leadGenMuonIsFound/O");
        razorTree->Branch("leadGenPhotonIsFound", &leadGenPhotonIsFound, "leadGenPhotonIsFound/O");
        razorTree->Branch("leadGenMuonIsFoundTight", &leadGenMuonIsFoundTight, "leadGenMuonIsFoundTight/O");
        razorTree->Branch("ptMatchingLeadGenMuon", &ptMatchingLeadGenMuon, "ptMatchingLeadGenMuon/F");
        razorTree->Branch("ptMatchingLeadGenPhoton", &ptMatchingLeadGenPhoton, "ptMatchingLeadGenPhoton/F");
    }
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("nMediumTaus", &nMediumTaus, "nMediumTaus/I");
    razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
    razorTree->Branch("leadingMuonPt", &leadingMuonPt, "leadingMuonPt/F");
    razorTree->Branch("leadingMuonEta", &leadingMuonEta, "leadingMuonEta/F");
    razorTree->Branch("leadingMuonPhi", &leadingMuonPhi, "leadingMuonPhi/F");
    razorTree->Branch("leadingMuonE", &leadingMuonE, "leadingMuonE/F");
    razorTree->Branch("leadingPhotonPt", &leadingPhotonPt, "leadingPhotonPt/F");
    razorTree->Branch("leadingPhotonEta", &leadingPhotonEta, "leadingPhotonEta/F");
    razorTree->Branch("leadingPhotonPhi", &leadingPhotonPhi, "leadingPhotonPhi/F");
    razorTree->Branch("leadingPhotonE", &leadingPhotonE, "leadingPhotonE/F");
    razorTree->Branch("subleadingPhotonPt", &subleadingPhotonPt, "subleadingPhotonPt/F");
    razorTree->Branch("subleadingPhotonEta", &subleadingPhotonEta, "subleadingPhotonEta/F");
    razorTree->Branch("subleadingPhotonPhi", &subleadingPhotonPhi, "subleadingPhotonPhi/F");
    razorTree->Branch("subleadingPhotonE", &subleadingPhotonE, "subleadingPhotonE/F");
    razorTree->Branch("subleadingMuonPt", &subleadingMuonPt, "subleadingMuonPt/F");
    razorTree->Branch("subleadingMuonEta", &subleadingMuonEta, "subleadingMuonEta/F");
    razorTree->Branch("subleadingMuonPhi", &subleadingMuonPhi, "subleadingMuonPhi/F");
    razorTree->Branch("subleadingMuonE", &subleadingMuonE, "subleadingMuonE/F");
    razorTree->Branch("recoZpt", &recoZpt, "recoZpt/F");
    razorTree->Branch("recoZeta", &recoZeta, "recoZeta/F");
    razorTree->Branch("recoZphi", &recoZphi, "recoZphi/F");
    razorTree->Branch("recoZmass", &recoZmass, "recoZmass/F");
    razorTree->Branch("recoWpt", &recoWpt, "recoWpt/F");
    razorTree->Branch("recoWphi", &recoWphi, "recoWphi/F");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("MR_noZ", &MR_noZ, "MR_noZ/F");
    razorTree->Branch("MR_noW", &MR_noW, "MR_noW/F");
    razorTree->Branch("MR_noPho", &MR_noPho, "MR_noPho/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("Rsq_noPho", &Rsq_noPho, "Rsq_noPho/F");
    razorTree->Branch("Rsq_noZ", &Rsq_noZ, "Rsq_noZ/F");
    razorTree->Branch("Rsq_noW", &Rsq_noW, "Rsq_noW/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("metphi", &metphi, "metphi/F");
    razorTree->Branch("met_noPho", &met_noPho, "met_noPho/F");
    razorTree->Branch("met_noZ", &met_noZ, "met_noZ/F");
    razorTree->Branch("met_noW", &met_noW, "met_noW/F");
    razorTree->Branch("metphi_noZ", &metphi_noZ, "metphi_noZ/F");
    razorTree->Branch("metphi_noW", &metphi_noW, "metphi_noW/F");
    razorTree->Branch("metphi_noPho", &metphi_noPho, "metphi_noPho/F");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("HT_noZ", &HT_noZ, "HT_noZ/F");
    razorTree->Branch("HT_noW", &HT_noW, "HT_noW/F");
    razorTree->Branch("HT_noPho", &HT_noPho, "HT_noPho/F");
    razorTree->Branch("numJets", &numJets, "numJets/I");
    razorTree->Branch("numJets_noZ", &numJets_noZ, "numJets_noZ/I");
    razorTree->Branch("numJets_noW", &numJets_noW, "numJets_noW/I");
    razorTree->Branch("numJets_noPho", &numJets_noPho, "numJets_noPho/I");
    razorTree->Branch("numJets80", &numJets80, "numJets80/I");
    razorTree->Branch("numJets80_noZ", &numJets80_noZ, "numJets80_noZ/I");
    razorTree->Branch("numJets80_noW", &numJets80_noW, "numJets80_noW/I");
    razorTree->Branch("numJets80_noPho", &numJets80_noPho, "numJets80_noPho/I");
    razorTree->Branch("mTLepMet", &mTLepMet, "mTLepMet/F");

    //****************************************************//
    //            Begin the event loop                    //
    //****************************************************//

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        //****************************************************//
        //               Initialize the event                 //
        //****************************************************//
        if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //fill normalization histogram
        NEvents->Fill(1.0);

        //reset tree variables
        if(!isData){
            nGenMuons = 0;
            nGenElectrons = 0;
            nGenTauMuons = 0;
            nGenTauElectrons = 0;
            nGenTaus = 0;
            nGenPhotons = 0;
            MR_noGenZ = -1;
            genZpt = -1;
            genZeta = -99;
            genZphi = -99;
            genZmass = -1;
            genWpt = -1;
            genWeta = -999;
            genWphi = -999;
            met_noGenZ = -1.;
            metphi_noGenZ = -99.;
            genmet = genMetPt;
            genmetphi = genMetPhi;
            Rsq_noGenZ = -1;
            HT_noGenZ = 0;
            numJets_noGenZ = 0;
            numJets80_noGenZ = 0;
            leadGenMuonIsFound = false;
            leadGenPhotonIsFound = false;
            leadGenMuonIsFoundTight = false;
            ptMatchingLeadGenMuon = -1;
            ptMatchingLeadGenPhoton = -1;
            leadingGenMuonPt = 0;
            leadingGenPhotonPt = 0;
            leadingGenMuonEta = -999;
            leadingGenPhotonEta = -999;
            leadingGenMuonPhi = -999;
            leadingGenPhotonPhi = -999;
            leadingGenMuonE = 0;
            leadingGenPhotonE = 0;
            subleadingGenMuonPt = 0;
            subleadingGenPhotonPt = 0;
            subleadingGenMuonEta = -999;
            subleadingGenPhotonEta = -999;
            subleadingGenMuonPhi = -999;
            subleadingGenPhotonPhi = -999;
            subleadingGenMuonE = 0;
            subleadingGenPhotonE = 0;
        }
        nSelectedJets = 0;
        nBTaggedJets = 0;
        nVetoMuons = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nVetoElectrons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nLooseTaus = 0;
        nMediumTaus = 0;
        nTightTaus = 0;
        nSelectedPhotons = 0;
        theMR = -1;
        MR_noZ = -1;
        MR_noW = -1;
        MR_noPho = -1;
        recoZpt = -1;
        recoZeta = -999;
        recoZphi = -999;
        recoZmass = -1;
        recoWpt = -1;
        recoWphi = -999;
        met = metPt;
        met_noPho = -1.;
        met_noZ = -1.;
        met_noW = -1.;
        metphi_noZ = -99.;
        metphi_noW = -99.;
        metphi_noPho = -99.;
        metphi = -99.;
        leadingMuonPt = -1;
        leadingMuonEta = -999;
        leadingMuonPhi = -999;
        leadingMuonE = -999;
        leadingPhotonPt = -1;
        leadingPhotonEta = -999;
        leadingPhotonPhi = -999;
        leadingPhotonE = -999;
        subleadingMuonPt = -1;
        subleadingMuonEta = -999;
        subleadingMuonPhi = -999;
        subleadingMuonE = -999;
        subleadingPhotonPt = -1;
        subleadingPhotonEta = -999;
        subleadingPhotonPhi = -999;
        subleadingPhotonE = -999;
        theRsq = -1;
        Rsq_noPho = -1;
        Rsq_noZ = -1;
        Rsq_noW = -1;
        HT = 0;
        HT_noPho = 0;
        HT_noZ = 0;
        HT_noW = 0;
        numJets = 0;
        numJets_noPho = 0;
        numJets_noZ = 0;
        numJets_noW = 0;
        numJets80 = 0;
        numJets80_noPho = 0;
        numJets80_noZ = 0;
        numJets80_noW = 0;
        mTLepMet = -1;

        //****************************************************//
        //               Select gen particles                 //
        //****************************************************//
        vector<TLorentzVector> GoodGenMuons; //for removing gen muons from jet collection later
        if(!isData){
            for(int j = 0; j < nGenParticle; j++){
                //electrons
                if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1) {
                    if (  (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23) ) {
                        nGenElectrons++;
                    }
                    if ( abs(gParticleMotherId[j]) == 15) {
                        nGenTauElectrons++;
                    }
                }
                //muons
                if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1) {
                    if ( (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)) {
                        nGenMuons++;
                        TLorentzVector thisGenMuon = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
                        GoodGenMuons.push_back(thisGenMuon);
                        if (gParticlePt[j] > leadingGenMuonPt) {
                            //make leading gen muon into subleading
                            subleadingGenMuonPt = leadingGenMuonPt;
                            subleadingGenMuonEta = leadingGenMuonEta;
                            subleadingGenMuonPhi = leadingGenMuonPhi;
                            subleadingGenMuonE = leadingGenMuonE; 
                            //make this the leading gen muon
                            leadingGenMuonPt = gParticlePt[j];
                            leadingGenMuonEta = gParticleEta[j];
                            leadingGenMuonPhi = gParticlePhi[j];
                            leadingGenMuonE = gParticleE[j];
                        }
                        else if(gParticlePt[j] > subleadingGenMuonPt){
                            //make this the subleading gen muon
                            subleadingGenMuonPt = gParticlePt[j];
                            subleadingGenMuonEta = gParticleEta[j];
                            subleadingGenMuonPhi = gParticlePhi[j];
                            subleadingGenMuonE = gParticleE[j];
                        }
                    }
                    if ( abs(gParticleMotherId[j]) == 15) {
                        nGenTauMuons++;
                    }
                }
                //taus
                if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
                        && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
                   ) nGenTaus++;
                //photons
                if (abs(gParticleId[j]) == 22 && 
                        ( (abs(gParticleMotherId[j]) >= 1 && abs(gParticleMotherId[j]) <= 5) || 
                          (abs(gParticleMotherId[j]) == 21) || (abs(gParticleMotherId[j]) == 2212) ) && 
                        abs(gParticleStatus[j]) == 1){	     
                    nGenPhotons++; 
                    if(gParticlePt[j] > leadingGenPhotonPt){
                        //make leading gen photon into subleading
                        subleadingGenPhotonPt = leadingGenPhotonPt;
                        subleadingGenPhotonEta = leadingGenPhotonEta;
                        subleadingGenPhotonPhi = leadingGenPhotonPhi;
                        subleadingGenPhotonE = leadingGenPhotonE;
                        //make this the leading gen photon
                        leadingGenPhotonPt = gParticlePt[j];
                        leadingGenPhotonEta = gParticleEta[j];
                        leadingGenPhotonPhi = gParticlePhi[j];
                        leadingGenPhotonE = gParticleE[j];
                    }
                    else if(gParticlePt[j] > subleadingGenPhotonPt){
                        //make this the subleading photon
                        subleadingGenPhotonPt = gParticlePt[j];
                        subleadingGenPhotonEta = gParticleEta[j];
                        subleadingGenPhotonPhi = gParticlePhi[j];
                        subleadingGenPhotonE = gParticleE[j];
                    }
                }
            }

            // gen level Z pt
            for(int j = 0; j < nGenParticle; j++){
                if(gParticleStatus[j] != 22) continue; //gen-level Z and W have pythia8 status 22
                TLorentzVector boson = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]); 
                if(abs(gParticleId[j]) == 23){ //Z boson
                    genZpt = gParticlePt[j];
                    genZeta = gParticleEta[j];
                    genZphi = gParticlePhi[j];
                    genZmass = boson.M();
                }
                else if(abs(gParticleId[j]) == 24){ //W boson
                    genWpt = gParticlePt[j];
                    genWeta = gParticleEta[j];
                    genWphi = gParticlePhi[j];
                }
            }
        }

        //****************************************************//
        //               Select muons                         //
        //****************************************************//
        vector<TLorentzVector> GoodMuons; 
        vector<TLorentzVector> GoodMuonsTight;
        for(int i = 0; i < nMuons; i++){

            if(!isLooseMuon(i)) continue;
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

            if(isVetoMuon(i)) nVetoMuons++;
            if(isTightMuon(i)){
                nTightMuons++;
                GoodMuonsTight.push_back(thisMuon);
            }
            nLooseMuons++;
            
            GoodMuons.push_back(thisMuon);

            //check if this muon is leading or subleading
            if(muonPt[i] > leadingMuonPt){
                //make leading muon into subleading
                subleadingMuonPt = leadingMuonPt;
                subleadingMuonEta = leadingMuonEta;
                subleadingMuonPhi = leadingMuonPhi;
                subleadingMuonE = leadingMuonE;
                //make this the leading muon
                leadingMuonPt = muonPt[i];
                leadingMuonEta = muonEta[i];
                leadingMuonPhi = muonPhi[i];
                leadingMuonE = muonE[i];
            }
            else if(muonPt[i] > subleadingMuonPt){
                //make this the subleading muon
                subleadingMuonPt = muonPt[i];
                subleadingMuonEta = muonEta[i];
                subleadingMuonPhi = muonPhi[i];
                subleadingMuonE = muonE[i];
            }

        }

        //****************************************************//
        //               Select electrons                     //
        //****************************************************//
        for(int i = 0; i < nElectrons; i++){

            if(!isLooseElectron(i)) continue;
            if(elePt[i] < 10) continue;
            if(fabs(eleEta[i]) > 2.5) continue;
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

            if(isMVANonTrigVetoElectron(i)) nVetoElectrons++;
            if(isTightElectron(i)){
                nTightElectrons++;
            }
            nLooseElectrons++;
        }

        //****************************************************//
        //               Select taus                          //
        //****************************************************//
        for(int i = 0; i < nTaus; i++){

            if(isLooseTau(i)){
                nLooseTaus++;
            }
            if(isMediumTau(i)){
                nMediumTaus++;
            }
            if(isTightTau(i)){
                nTightTaus++;
            }

        }

        //****************************************************//
        //               Select jets                          //
        //****************************************************//
        vector<TLorentzVector> GoodJets; //will contain leptons above 40 GeV in addition to jets
        for(int i = 0; i < nJets; i++){
            if(jetPt[i] < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            if(jetPt[i] > 80) numJets80++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVM(i)){ 
                nBTaggedJets++;
            }
        }
        sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());

        if(numJets80 < 2) continue; //event fails to have two 80 GeV jets
        if(GoodJets.size() > 15) continue; //TODO: instead do jet ID

        //****************************************************//
        //     Compute the razor variables and HT, nJets      //
        //****************************************************//
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
        metphi = metPhi;

        //count jets and compute HT
        numJets = GoodJets.size();
        for(auto& pf : GoodJets) HT += pf.Pt();

        // compute R and MR
        if(GoodJets.size() >= 2){
            vector<TLorentzVector> hemispheres = getHemispheres(GoodJets);
            theMR = computeMR(hemispheres[0], hemispheres[1]); 
            theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
        }

        //****************************************************//
        //             Select photons                         //
        //****************************************************//
        vector<TLorentzVector> GoodPhotons;
        vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
        int nPhotonsAbove40GeV = 0;
        for(int i = 0; i < nPhotons; i++){
            if(!isMediumPhoton(i)) continue;
            if(phoPt[i] < 10) continue;
            if(fabs(phoEta[i]) > 2.5) continue;

            if(phoPt[i] > 40) nPhotonsAbove40GeV++;
            TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], pho_RegressionE[i]);
            GoodPhotons.push_back(thisPhoton);
            GoodPhotonSigmaE.push_back(pho_RegressionEUncertainty[i]);
            nSelectedPhotons++;

            //check if this photon is leading or subleading
            if(phoPt[i] > leadingPhotonPt){
                //make leading photon into subleading
                subleadingPhotonPt = leadingPhotonPt;
                subleadingPhotonEta = leadingPhotonEta;
                subleadingPhotonPhi = leadingPhotonPhi;
                subleadingPhotonE = leadingPhotonE;
                //make this the leading photon
                leadingPhotonPt = phoPt[i];
                leadingPhotonEta = phoEta[i];
                leadingPhotonPhi = phoPhi[i];
                leadingPhotonE = phoE[i];
            }
            else if(phoPt[i] > subleadingPhotonPt){
                //make this the subleading photon
                subleadingPhotonPt = phoPt[i];
                subleadingPhotonEta = phoEta[i];
                subleadingPhotonPhi = phoPhi[i];
                subleadingPhotonE = phoE[i];
            }
        }

        if(GoodMuons.size() == 0 && GoodPhotons.size() == 0) continue; //don't save event if no muons or photons

        //****************************************************//
        //        Match gen-level and reco objects            //
        //****************************************************//
        if(!isData){
            if(nGenMuons > 0){
                //see if we selected a muon matching the leading gen muon
                TLorentzVector leadGenMuon = makeTLorentzVector(leadingGenMuonPt, leadingGenMuonEta, leadingGenMuonPhi, leadingGenMuonE);
                for(auto& mu : GoodMuons){
                    float thisDeltaR = mu.DeltaR(leadGenMuon);
                    if(thisDeltaR < 0.1){ //muon matches leading gen muon
                        leadGenMuonIsFound = true;
                        ptMatchingLeadGenMuon = mu.Pt(); //fill pt of the muon matching the gen muon
                        break;
                    }
                }
                for(auto& mu : GoodMuonsTight){
                    float thisDeltaR = mu.DeltaR(leadGenMuon);
                    if(thisDeltaR < 0.1){ //muon matches leading gen muon
                        leadGenMuonIsFoundTight = true;
                        break;
                    }
                }
            }
            if(nGenPhotons > 0){
                //see if we selected a photon matching the leading gen photon
                TLorentzVector leadGenPhoton = makeTLorentzVector(leadingGenPhotonPt, leadingGenPhotonEta, leadingGenPhotonPhi, leadingGenPhotonE);
                for(auto& pho : GoodPhotons){
                    float thisDeltaR = pho.DeltaR(leadGenPhoton);
                    if(thisDeltaR < 0.1){ //photon matches leading gen photon
                        leadGenPhotonIsFound = true;
                        ptMatchingLeadGenPhoton = pho.Pt(); //fill pt of the photon matching the gen photon
                        break;
                    }
                }
            }
        }

        //****************************************************//
        //    Compute razor vars for DY, W, Gamma samples     //
        //****************************************************//
        //photons
        if(GoodPhotons.size()>0){
            sort(GoodPhotons.begin(), GoodPhotons.end(), greater_than_pt());

            //compute MET with leading photon added
            TLorentzVector m1 = GoodPhotons[0];
            TLorentzVector m2 = PFMET;
            TLorentzVector photonPlusMet_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), 0.0);

            met_noPho = photonPlusMet_perp.Pt();
            metphi_noPho = photonPlusMet_perp.Phi();

            //remove leading photon from collection of selected jets
            vector<TLorentzVector> GoodJetsNoLeadPhoton = GoodJets;
            int subtractedIndex = SubtractParticleFromCollection(GoodPhotons[0], GoodJetsNoLeadPhoton);
            if(subtractedIndex >= 0){
                if(GoodJetsNoLeadPhoton[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoLeadPhoton.erase(GoodJetsNoLeadPhoton.begin()+subtractedIndex);
                }
            }
            //count the number of jets above 80 GeV now
            for(auto& jet : GoodJetsNoLeadPhoton){
                if(jet.Pt() > 80) numJets80_noPho++;
            }
            
            //count jets and compute HT
            numJets_noPho = GoodJetsNoLeadPhoton.size();
            for(auto& pf : GoodJetsNoLeadPhoton) HT_noPho += pf.Pt();

            if(GoodJetsNoLeadPhoton.size() >= 2){
                //remake the hemispheres using the new jet collection
                vector<TLorentzVector> hemispheresNoLeadPhoton = getHemispheres(GoodJetsNoLeadPhoton);
                TLorentzVector PFMET_NOPHO = makeTLorentzVectorPtEtaPhiM(met_noPho, 0, metphi_noPho, 0);
                MR_noPho = computeMR(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1]); 
                Rsq_noPho = computeRsq(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1], PFMET_NOPHO);
            }
        } 
        else{ //save some info even if no photons are found
            numJets_noPho = numJets;
            numJets80_noPho = numJets80;
            met_noPho = met;
            metphi_noPho = metphi;
            HT_noPho = HT;
        }

        // Muons for Z

        //remove selected muons from collection of selected jets and add them to the MET
        vector<TLorentzVector> GoodJetsNoMuons = GoodJets;
        TLorentzVector TotalMuonVec;
        for(auto& mu : GoodMuons){
            TotalMuonVec = TotalMuonVec + mu; //add this muon's momentum to the sum
            int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoMuons);
            if(subtractedIndex >= 0){
                if(GoodJetsNoMuons[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoMuons.erase(GoodJetsNoMuons.begin()+subtractedIndex);
                }
            }
        }
        //remove selected TIGHT muons from collection of selected jets and add them to the MET
        vector<TLorentzVector> GoodJetsNoTightMuons = GoodJets;
        TLorentzVector TotalTightMuonVec;
        for(auto& mu : GoodMuonsTight){
            TotalTightMuonVec = TotalTightMuonVec + mu; //add this muon's momentum to the sum
            int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoTightMuons);
            if(subtractedIndex >= 0){
                if(GoodJetsNoTightMuons[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoTightMuons.erase(GoodJetsNoTightMuons.begin()+subtractedIndex);
                }
            }
        }

        //do the same for GEN muons
        vector<TLorentzVector> GoodJetsNoGenMuons = GoodJets;
        TLorentzVector TotalGenMuonVec;
        if(!isData){
            for(auto& mu : GoodGenMuons){
                TotalGenMuonVec = TotalGenMuonVec + mu;
                int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoGenMuons);
                if(subtractedIndex >= 0){
                    if(GoodJetsNoGenMuons[subtractedIndex].Pt() < 40){ //erase this jet
                        GoodJetsNoGenMuons.erase(GoodJetsNoGenMuons.begin()+subtractedIndex);
                    }
                }
            }
        }

        //make the MET vector with the muons (or gen muons) added
        TLorentzVector ZPlusMet_perp = makeTLorentzVector((TotalMuonVec + PFMET).Pt(), 0., (TotalMuonVec + PFMET).Phi(), 0.);
        met_noZ = ZPlusMet_perp.Pt();
        metphi_noZ = ZPlusMet_perp.Phi();

        TLorentzVector WPlusMet_perp = makeTLorentzVector((TotalTightMuonVec + PFMET).Pt(), 0., (TotalTightMuonVec + PFMET).Phi(), 0.);
        met_noW = WPlusMet_perp.Pt();
        metphi_noW = WPlusMet_perp.Phi(); 

        TLorentzVector ZPlusMetGen_perp;
        if(!isData){
            ZPlusMetGen_perp = makeTLorentzVectorPtEtaPhiM((TotalGenMuonVec + PFMET).Pt(), 0., (TotalGenMuonVec + PFMET).Phi(), 0.);
            met_noGenZ = ZPlusMetGen_perp.Pt();
            metphi_noGenZ = ZPlusMetGen_perp.Phi();
        }

        //count jets and compute HT
        //Z
        numJets_noZ = GoodJetsNoMuons.size();
        for(auto& jet : GoodJetsNoMuons){
            HT_noZ += jet.Pt();
            if(jet.Pt() > 80) numJets80_noZ++;
        }
        //W
        numJets_noW = GoodJetsNoTightMuons.size();
        for(auto& jet : GoodJetsNoTightMuons){
            HT_noW += jet.Pt();
            if(jet.Pt() > 80) numJets80_noW++;
        }

        //get reco Z information
        if(GoodMuons.size() >= 1){
            recoZpt = TotalMuonVec.Pt();
            recoZeta = TotalMuonVec.Eta();
            recoZphi = TotalMuonVec.Phi();
            recoZmass = TotalMuonVec.M();
        }

        //compute reco Z information and razor variables for DY
        if(numJets_noZ > 1)
        {
            vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodJetsNoMuons);
            Rsq_noZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], ZPlusMet_perp);
            MR_noZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
        }
        if(!isData){
            //Gen Z
            numJets_noGenZ = GoodJetsNoGenMuons.size();
            for(auto& jet : GoodJetsNoGenMuons){
                HT_noGenZ += jet.Pt();
                if(jet.Pt() > 80) numJets80_noGenZ++;
            }
            //razor variables using GEN muons
            if(numJets_noGenZ > 1)
            {
                vector<TLorentzVector> hemispheresNoGenZ = getHemispheres(GoodJetsNoGenMuons);
                Rsq_noGenZ = computeRsq(hemispheresNoGenZ[0], hemispheresNoGenZ[1], ZPlusMetGen_perp);
                MR_noGenZ = computeMR(hemispheresNoGenZ[0], hemispheresNoGenZ[1]); 
            }
        }
        //razor variables using tight muons (for W)
        if(numJets_noW > 1){
            vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodJetsNoTightMuons);
            Rsq_noW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], WPlusMet_perp);
            MR_noW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 
        }

        //for W, also get the transverse mass of the first tight muon and the MET
        if(GoodMuonsTight.size() > 0) 
        {
            TLorentzVector m1 = GoodMuonsTight[0];
            TLorentzVector m2 = PFMET;
            double deltaPhiLepMet = m1.DeltaPhi(m2);
            mTLepMet = sqrt(2*m2.Pt()*m1.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) ); //transverse mass calculation

            //store reco W information
            recoWpt = (m1+m2).Pt();
            recoWphi = (m1+m2).Phi();
        }

        if(theMR < 300 && MR_noZ < 300 && MR_noW < 300 && MR_noPho < 300) continue;
        if(theRsq < 0.15 && Rsq_noZ < 0.15 && Rsq_noW < 0.15 && Rsq_noPho < 0.15) continue;

        razorTree->Fill();
    }//end of event loop

    cout << "Writing output tree..." << endl;
    razorTree->Write();
    NEvents->Write();

    outFile.Close();
}


