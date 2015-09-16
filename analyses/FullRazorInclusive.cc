#define RazorAnalyzer_cxx
#include "RazorAnalyzer.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

void RazorAnalyzer::FullRazorInclusive(string outFileName, bool isData)
{
    /////////////////////////////////
    //Basic setup
    /////////////////////////////////

    cout << "Initializing..." << endl;
    TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    //Output file
    if (outFileName.empty()){
        cout << "FullRazorInclusive: Output filename not specified!" << endl << "Using default output name FullRazorInclusive.root" << endl;
        outFileName = "FullRazorInclusive.root";
    }
    TFile outFile(outFileName.c_str(), "RECREATE");

    //Output tree
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

    //Histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    /////////////////////////////////
    //Pileup Weights
    /////////////////////////////////

    TFile *pileupWeightFile = 0;
    TH1D *pileupWeightHist = 0;
    // UNDER CONSTRUCTION (no Run 2 PU weights file available yet)
    if(!isData){
        pileupWeightFile = TFile::Open("data/ScaleFactors/Placeholders/DummyRun2PileupWeights.root");
        pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run2");
        assert(pileupWeightHist);
    }

    /////////////////////////////////
    //Lepton Efficiency Correction Factors
    /////////////////////////////////

    TH2D *eleLooseEffSFHist = 0;
    TH2D *eleTightEffSFHist = 0;
    TH2D *muLooseEffSFHist = 0;
    TH2D *muTightEffSFHist = 0;
    // UNDER CONSTRUCTION (no Run 2 scale factor files available yet)
    if(!isData){
        TFile *eleEffSFFile = TFile::Open("data/ScaleFactors/Placeholders/DummyRun2EleWeights.root");
        eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("EleWeight_Run2_Loose");
        assert(eleLooseEffSFHist);
        eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("EleWeight_Run2_Tight");
        assert(eleTightEffSFHist);
        TFile *muEffSFFile = TFile::Open("data/ScaleFactors/Placeholders/DummyRun2MuonWeights.root"); 
        muLooseEffSFHist = (TH2D*)muEffSFFile->Get("MuonWeight_Run2_Loose"); 
        assert(muLooseEffSFHist);
        muTightEffSFHist = (TH2D*)muEffSFFile->Get("MuonWeight_Run2_Tight");
        assert(muTightEffSFHist);
    }

    /////////////////////////////////
    //Jet Energy Corrections
    /////////////////////////////////

    //Get directory for JEC files 
    char* cmsswPath;
    cmsswPath = getenv("CMSSW_BASE");
    string pathname;
    if (cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
    else {
        cout << "Warning: CMSSW_BASE not detected.  Looking for JEC parameters in data/JEC";
        pathname = "data/JEC/";
    }

    //Get JEC files
    cout << "Getting JEC parameters from " << pathname << endl;
    std::vector<JetCorrectorParameters> correctionParameters;
    if (isData) {
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
    } 
    else {
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
        correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    }
    //Set up JEC machinery 
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
    JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
    SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);

    //Get JEC uncertainty file and set up JetCorrectionUncertainty
    string jecUncPath;
    if (isData) {
        jecUncPath = pathname+"/Summer15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
    }
    else {
        jecUncPath = pathname+"/Summer15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(jecUncPath);

    /////////////////////////////////
    //Tree Initialization
    /////////////////////////////////

    //Basic tree variables
    int nVtx, nSelectedJets, nBTaggedJets; 
    float dPhiRazor, MR, Rsq, mT, mTLoose, 
          leadingJetPt, subleadingJetPt, leadingTightMuPt, leadingTightElePt;
    float weight = 1.0;
    RazorBox box;
    //For lepton efficiency scale factor uncertainty
    float sf_muonEffUp, sf_muonEffDown;
    float sf_eleEffUp, sf_eleEffDown;
    //For btag scale factor uncertainty
    float sf_btagUp, sf_btagDown;
    //For jet uncertainties
    float MR_JESUp, Rsq_JESUp, dPhiRazor_JESUp, leadingJetPt_JESUp, subleadingJetPt_JESUp; 
    float MR_JESDown, Rsq_JESDown, dPhiRazor_JESDown, leadingJetPt_JESDown, subleadingJetPt_JESDown;
    float MR_JERUp, Rsq_JERUp, dPhiRazor_JERUp, leadingJetPt_JERUp, subleadingJetPt_JERUp;
    float MR_JERDown, Rsq_JERDown, dPhiRazor_JERDown, leadingJetPt_JERDown, subleadingJetPt_JERDown;
    int nSelectedJets_JESUp, nSelectedJets_JESDown, nSelectedJets_JERUp, nSelectedJets_JERDown;
    int nBTaggedJets_JESUp, nBTaggedJets_JESDown, nBTaggedJets_JERUp, nBTaggedJets_JERDown;
    RazorBox box_JESUp, box_JESDown, box_JERUp, box_JERDown;

    //Set branches
    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("MR", &MR, "MR/F");
    razorTree->Branch("Rsq", &Rsq, "Rsq/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("mT", &mT, "mT/F");
    razorTree->Branch("mTLoose", &mTLoose, "mTLoose/F");//for LooseLepton boxes
    razorTree->Branch("leadingTightMuPt", &leadingTightMuPt, "leadingTightMuPt/F");
    razorTree->Branch("leadingTightElePt", &leadingTightElePt, "leadingTightElePt/F");
    razorTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
    razorTree->Branch("subleadingJetPt", &subleadingJetPt, "subleadingJetPt/F");
    razorTree->Branch("box", &box, "box/I");
    if (!isData) {    
        razorTree->Branch("weight", &weight, "weight/F");
        razorTree->Branch("sf_muonEffUp", &sf_muonEffUp, "sf_muonEffUp/F");
        razorTree->Branch("sf_muonEffDown", &sf_muonEffDown, "sf_muonEffDown/F");
        razorTree->Branch("sf_eleEffUp", &sf_eleEffUp, "sf_eleEffUp/F");
        razorTree->Branch("sf_eleEffDown", &sf_eleEffDown, "sf_eleEffDown/F");
        razorTree->Branch("sf_btagUp", &sf_btagUp, "sf_btagUp/F");
        razorTree->Branch("sf_btagDown", &sf_btagDown, "sf_btagDown/F");
        razorTree->Branch("MR_JESUp", &MR_JESUp, "MR_JESUp/F");
        razorTree->Branch("Rsq_JESUp", &Rsq_JESUp, "Rsq_JESUp/F");
        razorTree->Branch("dPhiRazor_JESUp", &dPhiRazor_JESUp, "dPhiRazor_JESUp/F");
        razorTree->Branch("leadingJetPt_JESUp", &leadingJetPt_JESUp, "leadingJetPt_JESUp/F");
        razorTree->Branch("subleadingJetPt_JESUp", &subleadingJetPt_JESUp, "subleadingJetPt_JESUp/F");
        razorTree->Branch("nSelectedJets_JESUp", &nSelectedJets_JESUp, "nSelectedJets_JESUp/I");
        razorTree->Branch("nBTaggedJets_JESUp", &nBTaggedJets_JESUp, "nBTaggedJets_JESUp/I");
        razorTree->Branch("box_JESUp", &box_JESUp, "box_JESUp/I");
        razorTree->Branch("MR_JESDown", &MR_JESDown, "MR_JESDown/F");
        razorTree->Branch("Rsq_JESDown", &Rsq_JESDown, "Rsq_JESDown/F");
        razorTree->Branch("dPhiRazor_JESDown", &dPhiRazor_JESDown, "dPhiRazor_JESDown/F");
        razorTree->Branch("leadingJetPt_JESDown", &leadingJetPt_JESDown, "leadingJetPt_JESDown/F");
        razorTree->Branch("subleadingJetPt_JESDown", &subleadingJetPt_JESDown, "subleadingJetPt_JESDown/F");
        razorTree->Branch("nSelectedJets_JESDown", &nSelectedJets_JESDown, "nSelectedJets_JESDown/I");
        razorTree->Branch("nBTaggedJets_JESDown", &nBTaggedJets_JESDown, "nBTaggedJets_JESDown/I");
        razorTree->Branch("box_JESDown", &box_JESDown, "box_JESDown/I");
        razorTree->Branch("MR_JERUp", &MR_JERUp, "MR_JERUp/F");
        razorTree->Branch("Rsq_JERUp", &Rsq_JERUp, "Rsq_JERUp/F");
        razorTree->Branch("dPhiRazor_JERUp", &dPhiRazor_JERUp, "dPhiRazor_JERUp/F");
        razorTree->Branch("leadingJetPt_JERUp", &leadingJetPt_JERUp, "leadingJetPt_JERUp/F");
        razorTree->Branch("subleadingJetPt_JERUp", &subleadingJetPt_JERUp, "subleadingJetPt_JERUp/F");
        razorTree->Branch("nSelectedJets_JERUp", &nSelectedJets_JERUp, "nSelectedJets_JERUp/I");
        razorTree->Branch("nBTaggedJets_JERUp", &nBTaggedJets_JERUp, "nBTaggedJets_JERUp/I");
        razorTree->Branch("box_JERUp", &box_JERUp, "box_JERUp/I");
        razorTree->Branch("MR_JERDown", &MR_JERDown, "MR_JERDown/F");
        razorTree->Branch("Rsq_JERDown", &Rsq_JERDown, "Rsq_JERDown/F");
        razorTree->Branch("dPhiRazor_JERDown", &dPhiRazor_JERDown, "dPhiRazor_JERDown/F");
        razorTree->Branch("leadingJetPt_JERDown", &leadingJetPt_JERDown, "leadingJetPt_JERDown/F");
        razorTree->Branch("subleadingJetPt_JERDown", &subleadingJetPt_JERDown, "subleadingJetPt_JERDown/F");
        razorTree->Branch("nSelectedJets_JERDown", &nSelectedJets_JERDown, "nSelectedJets_JERDown/I");
        razorTree->Branch("nBTaggedJets_JERDown", &nBTaggedJets_JERDown, "nBTaggedJets_JERDown/I");
        razorTree->Branch("box_JERDown", &box_JERDown, "box_JERDown/I");
    } 
    else {
        razorTree->Branch("run", &runNum, "run/i");
        razorTree->Branch("lumi", &lumiNum, "lumi/i");
        razorTree->Branch("event", &eventNum, "event/i");
    }

    /////////////////////////////////
    //Event loop
    /////////////////////////////////

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        /////////////////////////////////
        //Begin event
        /////////////////////////////////

        //Initialize
        if (jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);

        //Fill normalization histogram
        NEvents->Fill(1.0);

        //Reset tree variables
        nVtx = nPV;
        nSelectedJets = 0;
        nBTaggedJets = 0;
        MR = -1;
        Rsq = -1;
        dPhiRazor = -9;
        mT = -1;
        mTLoose = -1;
        leadingJetPt = -1;
        subleadingJetPt = -1;
        leadingTightMuPt = -1;
        leadingTightElePt = -1;
        box = NONE;
        weight = 1.0;
        if(!isData){
            sf_muonEffUp = 1.0;
            sf_muonEffDown = 1.0;
            sf_eleEffUp = 1.0;
            sf_eleEffDown = 1.0;
            sf_btagUp = 1.0;
            sf_btagDown = 1.0;
            MR_JESUp = -1;
            Rsq_JESUp = -1;
            dPhiRazor_JESUp = -9;
            leadingJetPt_JESUp = -1;
            subleadingJetPt_JESUp = -1;
            nSelectedJets_JESUp = 0;
            nBTaggedJets_JESUp = 0;
            box_JESUp = NONE;
            MR_JESDown = -1;
            Rsq_JESDown = -1;
            dPhiRazor_JESDown = -9;
            leadingJetPt_JESDown = -1;
            subleadingJetPt_JESDown = -1;
            nSelectedJets_JESDown = 0;
            nBTaggedJets_JESDown = 0;
            box_JESDown = NONE;
            MR_JERUp = -1;
            Rsq_JERUp = -1;
            dPhiRazor_JERUp = -9;
            leadingJetPt_JERUp = -1;
            subleadingJetPt_JERUp = -1;
            nSelectedJets_JERUp = 0;
            nBTaggedJets_JERUp = 0;
            box_JERUp = NONE;
            MR_JERDown = -1;
            Rsq_JERDown = -1;
            dPhiRazor_JERDown = -9;
            leadingJetPt_JERDown = -1;
            subleadingJetPt_JERDown = -1;
            nSelectedJets_JERDown = 0;
            nBTaggedJets_JERDown = 0;
            box_JERDown = NONE;
        }

        //Reset TLorentzVector collections
        vector<TLorentzVector> GoodJets; //jets used to compute hemispheres
        vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
        //For systematics
        vector<TLorentzVector> GoodJetsJESUp;
        vector<TLorentzVector> GoodJetsJESDown;
        vector<TLorentzVector> GoodJetsJERUp;
        vector<TLorentzVector> GoodJetsJERDown;

        /////////////////////////////////
        //Trigger
        /////////////////////////////////

        bool passedDileptonTrigger = false;
        bool passedSingleLeptonTrigger = false;
        bool passedHadronicTrigger= false;

        if (isData) {
            passedDileptonTrigger = true;    
            passedSingleLeptonTrigger = HLTDecision[1] || HLTDecision[2] || HLTDecision[8] ||
                HLTDecision[20] || HLTDecision[22] || HLTDecision[24] || HLTDecision[25]  ;
            passedHadronicTrigger = true;
        } else {
            passedDileptonTrigger = true;  
            passedSingleLeptonTrigger = HLTDecision[1] || HLTDecision[2] || HLTDecision[8] ||
                HLTDecision[17] || HLTDecision[18] || HLTDecision[19] || HLTDecision[24]|| HLTDecision[25] ;
            passedHadronicTrigger = true;
        }

        if(!passedDileptonTrigger && !passedSingleLeptonTrigger && !passedHadronicTrigger) continue;

        /////////////////////////////////
        //Noise filters
        /////////////////////////////////
        //UNDER CONSTRUCTION

        /////////////////////////////////
        //Pileup reweighting
        /////////////////////////////////

        double NPU = 0;
        double pileupWeight = 1.0;
        if(!isData){
            //Get number of PU interactions
            for (int i = 0; i < nBunchXing; i++) {
                if (BunchXing[i] == 0) {
                    NPU = nPUmean[i];
                }
            }
            pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
        }

        /////////////////////////////////
        //Muon selection
        /////////////////////////////////

        int nVetoMuons = 0;
        int nLooseMuons = 0;
        int nTightMuons = 0;
        TLorentzVector leadingTightMu;
        //Scale factor
        float muonEffCorrFactor = 1.0;
        //Cut parameters
        const float MUON_VETO_CUT = 5;
        const float MUON_LOOSE_CUT = 30;
        //Loop muons
        for (int i = 0; i < nMuons; i++){

            //Baseline cuts
            if (muonPt[i] < MUON_VETO_CUT) continue;
            if (abs(muonEta[i]) > 2.4) continue;

            //Calculate MC->Data scale factors
            if (!isData && RazorAnalyzer::matchesGenMuon(muonEta[i], muonPhi[i])) {	
                //UNDER CONSTRUCTION (no efficiencies or scale factors available yet
                double effTight = 0.9; //NOTE: placeholder value
                double effTightSF = muTightEffSFHist->GetBinContent( 
                        muTightEffSFHist->GetXaxis()->FindFixBin(fabs(muonEta[i])) , 
                        muTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)));
                double effTightSFErr = muTightEffSFHist->GetBinError( 
                        muTightEffSFHist->GetXaxis()->FindFixBin(fabs(muonEta[i])) , 
                        muTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(muonPt[i],199.9),10.01)));
                double effTightSFUp = effTightSF + effTightSFErr;
                double effTightSFDown = effTightSF - effTightSFErr;
                double tmpTightSF = 1.0;
                double tmpTightSFUp = 1.0;
                double tmpTightSFDown = 1.0;
                if (isTightMuon(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                } 
                else {
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                }

                muonEffCorrFactor *= tmpTightSF;
                sf_muonEffUp *= tmpTightSFUp;
                sf_muonEffDown *= tmpTightSFDown;
            }

            //TLorentzVector for this muon
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

            //Veto selection
            if (isVetoMuon(i)){
                nVetoMuons++;
                GoodLeptons.push_back(thisMuon); 
            }
            //Loose selection
            if (isLooseMuon(i) && muonPt[i] >= MUON_LOOSE_CUT) nLooseMuons++;
            //Tight selection
            if (isTightMuon(i) && muonPt[i] >= MUON_LOOSE_CUT){
                nTightMuons++;
                if (muonPt[i] > leadingTightMuPt){
                    leadingTightMu = thisMuon;
                    leadingTightMuPt = muonPt[i];
                }
            }
        }

        /////////////////////////////////
        //Electron selection
        /////////////////////////////////
        int nVetoElectrons = 0;
        int nLooseElectrons = 0;
        int nTightElectrons = 0;
        TLorentzVector leadingTightEle; //used for mT calculation
        float eleEffCorrFactor = 1.0;
        //Cut parameters
        const float ELE_VETO_CUT = 5;
        const float ELE_LOOSE_CUT = 30;
        //Loop electrons
        for (int i = 0; i < nElectrons; i++){

            //Baseline cuts
            if (elePt[i] < ELE_VETO_CUT) continue;
            if (fabs(eleEta[i]) > 2.5) continue;

            //Calculate MC->Data scale factors
            if (!isData && RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i])) {
                //UNDER CONSTRUCTION (no scale factors or efficiencies available for run 2)
                //Tight scale factor
                double effTight = 0.9; //NOTE: placeholder value
                double effTightSF = eleTightEffSFHist->GetBinContent( 
                        eleTightEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
                        eleTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
                double effTightSFErr = eleTightEffSFHist->GetBinError( 
                        eleTightEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
                        eleTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
                double effTightSFUp = effTightSF + effTightSFErr;
                double effTightSFDown = effTightSF - effTightSFErr;
                double tmpTightSF = 1.0;
                double tmpTightSFUp = 1.0;
                double tmpTightSFDown = 1.0;
                if (isTightElectron(i)) {
                    tmpTightSF = effTightSF;
                    tmpTightSFUp = effTightSFUp;
                    tmpTightSFDown = effTightSFDown;
                } 
                else {
                    tmpTightSF = (1/effTight - effTightSF) / (1/effTight - 1);
                    tmpTightSFUp = (1/effTight - effTightSFUp) / (1/effTight - 1);
                    tmpTightSFDown = (1/effTight - effTightSFDown) / (1/effTight - 1);
                }

                eleEffCorrFactor *= tmpTightSF;
                sf_eleEffUp *= tmpTightSFUp;
                sf_eleEffDown *= tmpTightSFDown;
            }

            //Remove overlaps
            bool overlap = false;
            for (auto& lep : GoodLeptons){
                if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;

            //TLorentzVector for this electron
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

            //Veto selection
            if (isVetoElectron(i)){
                nVetoElectrons++;
                GoodLeptons.push_back(thisElectron);            
            }
            //Loose selection
            if (isLooseElectron(i) && elePt[i] > ELE_LOOSE_CUT) nLooseElectrons++;
            if (isTightElectron(i) && elePt[i] > ELE_LOOSE_CUT){
                nTightElectrons++;
                if (elePt[i] > leadingTightElePt){
                    leadingTightEle = thisElectron;
                    leadingTightElePt = elePt[i];
                }
            }
        }

        /////////////////////////////////
        //Tau selection
        /////////////////////////////////

        int nLooseTaus = 0;
        const float TAU_LOOSE_CUT = 20;
        //Loop taus
        for (int i = 0; i < nTaus; i++){	 
            //Baseline cuts
            if (tauPt[i] < TAU_LOOSE_CUT) continue;
            if (fabs(tauEta[i]) > 2.4) continue;

            //remove overlaps
            bool overlap = false;
            for (auto& lep : GoodLeptons){
                if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;

            //Loose selection
            if (isLooseTau(i)){
                nLooseTaus++;
                TLorentzVector thisTau = makeTLorentzVectorPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
                GoodLeptons.push_back(thisTau);  
            }
        }

        /////////////////////////////////
        //Jet selection
        /////////////////////////////////

        //BTag scale factor
        float btagCorrFactor = 1.0;
        //Propagate jet uncertainties to MET
        float MetXCorr_JESUp = 0;
        float MetYCorr_JESUp = 0;
        float MetXCorr_JESDown = 0;
        float MetYCorr_JESDown = 0;
        float MetXCorr_JERUp = 0;
        float MetYCorr_JERUp = 0;
        float MetXCorr_JERDown = 0;
        float MetYCorr_JERDown = 0;
        //Loop jets
        for (int i = 0; i < nJets; i++){

            //Remove overlaps
            double deltaR = -1;
            for (auto& lep : GoodLeptons){
                double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if (deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

            //Apply Jet ID
            if (!jetPassIDTight[i]) continue;

            //Get jet energy correction
            double tmpRho = fixedGridRhoFastjetAll;
            double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], JetCorrector);   

            //Get jet energy resolution correction, with up/down variants
            double jetEnergySmearFactor = 1.0;
            double jetEnergySmearFactorUp = 1.0;
            double jetEnergySmearFactorDown = 1.0;
            if (!isData) {
                jetEnergySmearFactor = JetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random);
                jetEnergySmearFactorUp = JetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random, "up");
                jetEnergySmearFactorDown = JetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random, "down");
            }

            //TLorentzVector for this jet
            double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
            double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;
            TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);

            //Apply b-tagging correction factor 
            //UNDER CONSTRUCTION (no b-tagging corrections for run 2 yet)
            if (!isData && abs(jetPartonFlavor[i]) == 5 && jetCorrPt > 20) {
                btagCorrFactor *= BTagScaleFactor(jetCorrPt, isCSVM(i));
                sf_btagUp *= BTagScaleFactor(jetCorrPt, isCSVM(i), "up");
                sf_btagDown *= BTagScaleFactor(jetCorrPt, isCSVM(i), "down");
            }

            //Apply pileup jet ID 
            //UNDER CONSTRUCTION (No working point yet for Run2)
            //int level = 2; //loose jet ID
            //if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;

            //Cut on jet eta
            if (fabs(jetEta[i]) > 3.0) continue;

            //Get uncertainty on JEC and JER
            if(!isData){
                jecUnc->setJetEta(jetEta[i]);
                jecUnc->setJetPt(jetCorrPt);
                double unc = jecUnc->getUncertainty(true);
                double jetPtJESUp = jetCorrPt*(1+unc);
                double jetPtJESDown = jetCorrPt/(1+unc);
                double jetPtJERUp = jetPt[i]*JEC*jetEnergySmearFactorUp;
                double jetPtJERDown = jetPt[i]*JEC*jetEnergySmearFactorDown;
                double jetEJESUp = jetCorrE*(1+unc);
                double jetEJESDown = jetCorrE/(1+unc);
                double jetEJERUp = jetE[i]*JEC*jetEnergySmearFactorUp;
                double jetEJERDown = jetE[i]*JEC*jetEnergySmearFactorDown;
                TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
                TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
                TLorentzVector thisJetJERUp = makeTLorentzVector(jetPtJERUp, jetEta[i], jetPhi[i], jetEJERUp);
                TLorentzVector thisJetJERDown = makeTLorentzVector(jetPtJERDown, jetEta[i], jetPhi[i], jetEJERDown);

                //Propagate uncertainties to the MET
                if (jetPtJESUp > 20) {
                    MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
                    MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
                }
                if (jetPtJESDown > 20) {
                    MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
                    MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
                }
                if (jetPtJERUp > 20) {
                    MetXCorr_JERUp += -1 * (thisJetJERUp.Px() - thisJet.Px());
                    MetYCorr_JERUp += -1 * (thisJetJERUp.Py() - thisJet.Py());
                }
                if (jetPtJERDown > 20) {
                    MetXCorr_JERDown += -1 * (thisJetJERDown.Px() - thisJet.Px());
                    MetYCorr_JERDown += -1 * (thisJetJERDown.Py() - thisJet.Py());
                }
                //Record jets that pass the cut
                if(jetPtJESUp > 40){
                    GoodJetsJESUp.push_back(thisJetJESUp);
                    nSelectedJets_JESUp++;
                    if (isCSVM(i)) nBTaggedJets_JESUp++;
                }
                if(jetPtJESDown > 40){
                    GoodJetsJESDown.push_back(thisJetJESDown);
                    nSelectedJets_JESDown++;
                    if (isCSVM(i)) nBTaggedJets_JESDown++;
                }
                if(jetPtJERUp > 40){
                    GoodJetsJERUp.push_back(thisJetJERUp);
                    nSelectedJets_JERUp++;
                    if (isCSVM(i)) nBTaggedJets_JERUp++;
                }
                if(jetPtJERDown > 40){
                    GoodJetsJERDown.push_back(thisJetJERDown);
                    nSelectedJets_JERDown++;
                    if (isCSVM(i)) nBTaggedJets_JERDown++;
                }
            }

            //Cut on jet pt
            if (jetCorrPt < 40) continue;

            //Record this jet
            GoodJets.push_back(thisJet);
            nSelectedJets++;
            if (isCSVM(i)){ 
                nBTaggedJets++;
            }
        }

        //Get leading and subleading jet pt
        for (auto &jet : GoodJets){
            if (jet.Pt() > leadingJetPt){
                subleadingJetPt = leadingJetPt;
                leadingJetPt = jet.Pt();
            }
            else if (jet.Pt() > subleadingJetPt){
                subleadingJetPt = jet.Pt();
            }
        }
        //Get leading and subleading jet pt for JES/JER up/down
        if (!isData){
            for (auto &jet : GoodJetsJESUp){
                if (jet.Pt() > leadingJetPt_JESUp){
                    subleadingJetPt_JESUp = leadingJetPt_JESUp;
                    leadingJetPt_JESUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JESUp){
                    subleadingJetPt_JESUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJESDown){
                if (jet.Pt() > leadingJetPt_JESDown){
                    subleadingJetPt_JESDown = leadingJetPt_JESDown;
                    leadingJetPt_JESDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JESDown){
                    subleadingJetPt_JESDown = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJERUp){
                if (jet.Pt() > leadingJetPt_JERUp){
                    subleadingJetPt_JERUp = leadingJetPt_JERUp;
                    leadingJetPt_JERUp = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JERUp){
                    subleadingJetPt_JERUp = jet.Pt();
                }
            }
            for (auto &jet : GoodJetsJERDown){
                if (jet.Pt() > leadingJetPt_JERDown){
                    subleadingJetPt_JERDown = leadingJetPt_JERDown;
                    leadingJetPt_JERDown = jet.Pt();
                }
                else if (jet.Pt() > subleadingJetPt_JERDown){
                    subleadingJetPt_JERDown = jet.Pt();
                }
            }
        }

        /////////////////////////////////
        //Compute razor variables and mT
        /////////////////////////////////

        //Combine jet and lepton collections
        for (auto& lep : GoodLeptons) {
            GoodJets.push_back(lep);
            if(!isData){
                GoodJetsJESUp.push_back(lep);
                GoodJetsJESDown.push_back(lep);
                GoodJetsJERUp.push_back(lep);
                GoodJetsJERDown.push_back(lep);
            }
        }

        //Get MET
        TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
        TLorentzVector MyMET = PFMETType1; //This is the MET that will be used below.

        //Compute razor variables and dPhiRazor
        vector<TLorentzVector> hemispheres = getHemispheres(GoodJets);
        MR = computeMR(hemispheres[0], hemispheres[1]); 
        Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
        dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());

        //Propagate up/down jet uncertainties to MET and recompute razor variables
        if(!isData){
            //JES up
            float PFMetXJESUp = MyMET.Px() + MetXCorr_JESUp;
            float PFMetYJESUp = MyMET.Py() + MetYCorr_JESUp;
            TLorentzVector PFMET_JESUp(PFMetXJESUp, PFMetYJESUp, 0, sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) )); 
            vector<TLorentzVector> hemispheres_JESUp = getHemispheres(GoodJetsJESUp);
            MR_JESUp = computeMR(hemispheres_JESUp[0], hemispheres_JESUp[1]); 
            Rsq_JESUp = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], PFMET_JESUp);
            dPhiRazor_JESUp = deltaPhi(hemispheres_JESUp[0].Phi(),hemispheres_JESUp[1].Phi());

            //JES down
            float PFMetXJESDown = MyMET.Px() + MetXCorr_JESDown;
            float PFMetYJESDown = MyMET.Py() + MetYCorr_JESDown;
            TLorentzVector PFMET_JESDown(PFMetXJESDown, PFMetYJESDown, 0, sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) )); 
            vector<TLorentzVector> hemispheres_JESDown = getHemispheres(GoodJetsJESDown);
            MR_JESDown = computeMR(hemispheres_JESDown[0], hemispheres_JESDown[1]); 
            Rsq_JESDown = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], PFMET_JESDown);
            dPhiRazor_JESDown = deltaPhi(hemispheres_JESDown[0].Phi(),hemispheres_JESDown[1].Phi());

            //JER up
            float PFMetXJERUp = MyMET.Px() + MetXCorr_JERUp;
            float PFMetYJERUp = MyMET.Py() + MetYCorr_JERUp;
            TLorentzVector PFMET_JERUp(PFMetXJERUp, PFMetYJERUp, 0, sqrt( pow(PFMetXJERUp,2) + pow(PFMetYJERUp,2) )); 
            vector<TLorentzVector> hemispheres_JERUp = getHemispheres(GoodJetsJERUp);
            MR_JERUp = computeMR(hemispheres_JERUp[0], hemispheres_JERUp[1]); 
            Rsq_JERUp = computeRsq(hemispheres_JERUp[0], hemispheres_JERUp[1], PFMET_JERUp);
            dPhiRazor_JERUp = deltaPhi(hemispheres_JERUp[0].Phi(),hemispheres_JERUp[1].Phi());

            //JER down
            float PFMetXJERDown = MyMET.Px() + MetXCorr_JERDown;
            float PFMetYJERDown = MyMET.Py() + MetYCorr_JERDown;
            TLorentzVector PFMET_JERDown(PFMetXJERDown, PFMetYJERDown, 0, sqrt( pow(PFMetXJERDown,2) + pow(PFMetYJERDown,2) )); 
            vector<TLorentzVector> hemispheres_JERDown = getHemispheres(GoodJetsJERDown);
            MR_JERDown = computeMR(hemispheres_JERDown[0], hemispheres_JERDown[1]); 
            Rsq_JERDown = computeRsq(hemispheres_JERDown[0], hemispheres_JERDown[1], PFMET_JERDown);
            dPhiRazor_JERDown = deltaPhi(hemispheres_JERDown[0].Phi(),hemispheres_JERDown[1].Phi());
        }

        //Compute transverse mass 
        if (nTightMuons + nTightElectrons > 0){
            TLorentzVector leadingLepton;
            if (leadingTightMuPt > leadingTightElePt) leadingLepton = leadingTightMu;
            else leadingLepton = leadingTightEle;

            float deltaPhiLepMet = leadingLepton.DeltaPhi(MyMET);
            mT = sqrt(2*leadingLepton.Pt()*MyMET.Pt()*(1.0 - cos(deltaPhiLepMet))); 
        }
        //Transverse mass with leading lepton, regardless of quality
        if (GoodLeptons.size() > 0){
            //get the highest pt lepton
            float maxLepPt = -1;
            int maxLepIndex = -1;
            for (uint i = 0; i < GoodLeptons.size(); i++){
                if (GoodLeptons[i].Pt() > maxLepPt){
                    maxLepPt = GoodLeptons[i].Pt();
                    maxLepIndex = i;
                }
            }
            if (maxLepIndex >= 0){
                float deltaPhiLepMet = GoodLeptons[maxLepIndex].DeltaPhi(MyMET);
                mTLoose = sqrt(2*GoodLeptons[maxLepIndex].Pt()*MyMET.Pt()*(1.0 - cos(deltaPhiLepMet)));
            }
        }

        /////////////////////////////////
        //Make baseline razor cuts
        /////////////////////////////////
        if (MR < 300 && MR_JESUp < 300 && MR_JESDown < 300 && MR_JERUp < 300 && MR_JERDown < 300) continue;
        if (Rsq < 0.15 && Rsq_JESUp < 0.15 && Rsq_JESDown < 0.15 && Rsq_JERUp < 0.15 && Rsq_JERDown < 0.15) continue;

        /////////////////////////////////
        //Apply scale factors
        /////////////////////////////////

        //Nominal event weight
        weight *= pileupWeight;
        weight *= muonEffCorrFactor;
        weight *= eleEffCorrFactor;
        weight *= btagCorrFactor;    

        /////////////////////////////////
        //Categorize into boxes
        /////////////////////////////////

        //Get correct box under up/down JES/JER systematic
        if(!isData){
            if(passedDileptonTrigger && nTightElectrons > 0 && nLooseMuons > 0){
                box_JESUp = MuEle;
                box_JESDown = MuEle;
                box_JERUp = MuEle;
                box_JERDown = MuEle;
            }
            else if(passedDileptonTrigger && nTightMuons > 0 && nLooseMuons > 1){
                box_JESUp = MuMu;
                box_JESDown = MuMu;
                box_JERUp = MuMu;
                box_JERDown = MuMu;
            }
            else if(passedDileptonTrigger && nTightElectrons > 0 && nLooseElectrons > 1){
                box_JESUp = EleEle;
                box_JESDown = EleEle;
                box_JERUp = EleEle;
                box_JERDown = EleEle;
            }
            else if (passedSingleLeptonTrigger && nTightMuons > 0){
                if (nSelectedJets_JESUp > 5) box_JESUp = MuSixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = MuFourJet;
                else box_JESUp = MuJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = MuSixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = MuFourJet;
                else box_JESDown = MuJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = MuSixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = MuFourJet;
                else box_JERUp = MuJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = MuSixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = MuFourJet;
                else box_JERDown = MuJet;
            }
            else if (passedSingleLeptonTrigger && nTightElectrons > 0){
                if (nSelectedJets_JESUp > 5) box_JESUp = EleSixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = EleFourJet;
                else box_JESUp = EleJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = EleSixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = EleFourJet;
                else box_JESDown = EleJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = EleSixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = EleFourJet;
                else box_JERUp = EleJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = EleSixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = EleFourJet;
                else box_JERDown = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0){
                if (nSelectedJets_JESUp > 5) box_JESUp = LooseLeptonSixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = LooseLeptonFourJet;
                else box_JESUp = LooseLeptonDiJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = LooseLeptonSixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = LooseLeptonFourJet;
                else box_JESDown = LooseLeptonDiJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = LooseLeptonSixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = LooseLeptonFourJet;
                else box_JERUp = LooseLeptonDiJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = LooseLeptonSixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = LooseLeptonFourJet;
                else box_JERDown = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger){
                if (nSelectedJets_JESUp > 5) box_JESUp = SixJet;
                else if (nSelectedJets_JESUp > 3) box_JESUp = FourJet;
                else box_JESUp = DiJet;

                if (nSelectedJets_JESDown > 5) box_JESDown = SixJet;
                else if (nSelectedJets_JESDown > 3) box_JESDown = FourJet;
                else box_JESDown = DiJet;

                if (nSelectedJets_JERUp > 5) box_JERUp = SixJet;
                else if (nSelectedJets_JERUp > 3) box_JERUp = FourJet;
                else box_JERUp = DiJet;

                if (nSelectedJets_JERDown > 5) box_JERDown = SixJet;
                else if (nSelectedJets_JERDown > 3) box_JERDown = FourJet;
                else box_JERDown = DiJet;
            }
        }

        //Nominal box
        if (passedDileptonTrigger && nTightElectrons > 0 && nLooseMuons > 0){
            box = MuEle;
        }
        else if (passedDileptonTrigger && nTightMuons > 0 && nLooseMuons > 1){
            box = MuMu;
        }
        else if (passedDileptonTrigger && nTightElectrons > 0 && nLooseElectrons > 1){
            box = EleEle;
        }
        else if (passedSingleLeptonTrigger && nTightMuons > 0){
            if (nSelectedJets > 5) box = MuSixJet;
            else if (nSelectedJets > 3) box = MuFourJet;
            else box = MuJet;
        }
        else if (passedSingleLeptonTrigger && nTightElectrons > 0){
            if (nSelectedJets > 5) box = EleSixJet;
            else if (nSelectedJets > 3) box = EleFourJet;
            else box = EleJet;
        }
        else if (passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0){
            if (nSelectedJets > 5) box = LooseLeptonSixJet;
            else if (nSelectedJets > 3) box = LooseLeptonFourJet;
            else box = LooseLeptonDiJet;
        }
        else if (passedHadronicTrigger){
            if (nSelectedJets > 5) box = SixJet;
            else if (nSelectedJets > 3) box = FourJet;
            else box = DiJet;
        }

        //Continue if this event is not in any box
        if(box == NONE && box_JESUp == NONE && box_JESDown == NONE && box_JERUp == NONE && box_JERDown == NONE) continue; 

        //Fill tree
        razorTree->Fill();

    }//end of event loop

    cout << "Writing output tree..." << endl;
    outFile.cd();
    razorTree->Write();
    NEvents->Write();

    outFile.Close();
}
