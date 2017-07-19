#include "DelayedPhotonAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns


float DelayedPhotonAnalyzer::getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID) {
  float timeCalib = 1.0;
  
  int N_entries = tree->GetEntries(); 
  int i_entry=0;
  for(uint i=0;i<start_run.size();i++) {
    if(run>= start_run[i] && run<= end_run[i]) {
      i_entry = i;
      break;
    }
  }

  
  if(i_entry> N_entries) return timeCalib;
  tree->GetEntry(i_entry);
  std::vector<int>::iterator p_id;
  p_id = std::find(detID_all->begin(), detID_all->end(), detID);
  if (p_id == detID_all->end()) return timeCalib;
  uint idx = std::distance(detID_all->begin(), p_id);
  
  if(idx<=IC_time_all->size()) timeCalib = IC_time_all->at(idx);  
  
  return timeCalib;
};

float DelayedPhotonAnalyzer::getPedestalNoise(TTree *tree, vector <uint> & start_time, vector <uint> & end_time, uint time, uint detID) {
  float pedestalNoise = 1.0;
  
  int N_entries = tree->GetEntries();
  int i_entry=0;
  for(uint i=0;i<start_time.size();i++) {
    if(time>= start_time[i] && time<= end_time[i])
    {
      i_entry = i;
      break;
    }
  }
  
  if(i_entry> N_entries) return pedestalNoise;
  tree->GetEntry(i_entry);
  std::vector<int>::iterator p_id;
  p_id = std::find(detID_all->begin(), detID_all->end(), detID);
  if (p_id == detID_all->end()) return pedestalNoise;
  uint idx = std::distance(detID_all->begin(), p_id);
  
  if(idx<=rms_G12_all->size()) pedestalNoise = rms_G12_all->at(idx);  
  
  return pedestalNoise;
};


float DelayedPhotonAnalyzer::getADCToGeV( uint run, int isEBOrEE) {
  double ADCToGeV = 0;
  //EB
  if (isEBOrEE == 0) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.039680;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.039798;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.039436;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.039298;
  }   
  //EE
  else if (isEBOrEE == 1) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.067230;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.067370;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.066764;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.065957;
  }
  return ADCToGeV;
}


void DelayedPhotonAnalyzer::Analyze(bool isData, int option, string outFileName, string label) {

  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  //bool doPhotonScaleCorrection = true;

  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;
  

  //*****************************************************************************
  //Load Intercalibration constants
  //*****************************************************************************
  vector <uint> start_run;//start run of all IOV 
  vector <uint> end_run;//end run of all IOV
  vector <uint> start_run_rereco;// for SepRereco tags
  vector <uint> end_run_rereco;// for SepRereco tags
  start_run_tmp=0; 
  end_run_tmp=0;
  IC_time_all=0;
  detID_all=0;

  TFile f_timeCalib("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalTimeCalibConstants_Legacy2016_v1/EcalTimeCalibConstants_Legacy2016_v1.root","READ");
  TTree *tree_timeCalib = (TTree*)f_timeCalib.Get("timeCalib");
  
  tree_timeCalib->SetBranchAddress("start_run", &start_run_tmp);
  tree_timeCalib->SetBranchAddress("end_run", &end_run_tmp);
  tree_timeCalib->SetBranchAddress("IC_time", &IC_time_all);
  tree_timeCalib->SetBranchAddress("detID", &detID_all);
  
  int N_entries_timeCalib = tree_timeCalib->GetEntries();
  
  for(int i=0;i<N_entries_timeCalib;i++) {
    tree_timeCalib->GetEntry(i);
    start_run.push_back(start_run_tmp);
    end_run.push_back(end_run_tmp);
  }


  TFile f_timeCalib_rereco("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalTimeCalibConstants_v08_offline/tree_EcalTimeCalibConstants_v08_offline.root","READ");
  TTree *tree_timeCalib_rereco = (TTree*)f_timeCalib_rereco.Get("timeCalib");
  
  tree_timeCalib_rereco->SetBranchAddress("start_run", &start_run_tmp);
  tree_timeCalib_rereco->SetBranchAddress("end_run", &end_run_tmp);
  tree_timeCalib_rereco->SetBranchAddress("IC_time", &IC_time_all);
  tree_timeCalib_rereco->SetBranchAddress("detID", &detID_all);
  
  int N_entries_timeCalib_rereco = tree_timeCalib_rereco->GetEntries();
  
  for(int i=0;i<N_entries_timeCalib_rereco;i++) {
    tree_timeCalib_rereco->GetEntry(i);
    start_run_rereco.push_back(start_run_tmp);
    end_run_rereco.push_back(end_run_tmp);
  }


  //*****************************************************************************
  //Load Pedestals
  //*****************************************************************************
  vector <uint> start_time;//start run of all IOV 
  vector <uint> end_time;//end run of all IOV
  start_time_tmp=0; 
  end_time_tmp=0;
  rms_G12_all=0;
  detID_all=0 ;

  // TFile f_pedestal("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1.root","READ");
  // TTree *tree_pedestal = (TTree*)f_pedestal.Get("pedestal");
  
  // tree_pedestal->SetBranchAddress("start_time_second", &start_time_tmp);
  // tree_pedestal->SetBranchAddress("end_time_second", &end_time_tmp);
  // tree_pedestal->SetBranchAddress("rms_G12", &rms_G12_all);
  // tree_pedestal->SetBranchAddress("detID", &detID_all);
  
  // int N_entries_pedestal = tree_pedestal->GetEntries();
  
  // cout << "Total Pedestal IOVs: " << N_entries_pedestal << "\n";
  // for(int i=0;i<N_entries_pedestal;i++) {
  //   cout << "Loading Pedestal IOV " << i << "\n";
  //   tree_pedestal->GetEntry(i);
  //   start_time.push_back(start_time_tmp);
  //   end_time.push_back(end_time_tmp);
  // }

  // //test 
  // uint test_time = 1464000000;
  // //cout<<"EB test..."<<endl;
  // for(int ieta=-85;ieta<=85 && ieta!=0;ieta++) {
  //   for(int iphi=1;iphi<=360;iphi++) {
  //     int detID = detID_from_iEtaiPhi(ieta, iphi, true, false);
  //     cout<<test_time<<"  "<<ieta<<"  "<<iphi<<"  "<<detID;      
  //     float pedestalRMS = getPedestalNoise(tree_pedestal, start_time,end_time,test_time, detID);
  //     cout << "   " << pedestalRMS << endl;      
  //   }
  // }
  
  
  //*****************************************************************************
  //Open Output File
  //*****************************************************************************
  if ( outFileName.empty() ) {
    std::cout << "DelayedPhotonRazor: Output filename not specified!" << endl << "Using default output name DelayedPhotonRazor.root" << std::endl;
    outFileName = "DelayedPhotonRazor.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );




  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("DelayedPhotonAnalyzer", "Info on selected razor inclusive events");

  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float mass;
  float t1, t2;
  float t1_seed, t2_seed;
  float t1calib_seed, t2calib_seed;
  float t1raw_seed, t2raw_seed;
  float phoNumber;
  float pho1Energy, pho1_Pt, pho1_Eta, pho1_Phi;
  float pho2Energy, pho2_Pt, pho2_Eta, pho2_Phi;
  float jetEnergy, jet_Pt, jet_Eta, jet_Phi;
  float met_Pt, met_Phi, sum_MET;

  float pho2SeedIY;
  float pho2SeedIX;
  float pho1SeedIY;
  float pho1SeedIX;
  float pho2SeedIPhi;
  float pho2SeedIEta;
  float pho1SeedIPhi;
  float pho1SeedIEta;
  float phoEta_SC[1000];
  bool pho2IsEB;
  bool pho1IsEB;
  float calibratedSeedHitTime;


  //pho_hasPixelSeed, pho_passHLTFilter, pho_convType, pho_convTrkZ, pho_convTrkClusZ, pho_vtxSumPx, pho_vtxSumPy, pho_seedRecHitSwitchToGain6, pho_seedRecHitSwitchToGain1, pho_anyRecHitSwitchToGain6, 
  //pho_anyRecHitSwitchToGain1, pho_SeedRechitIndex, phoSigmaIetaIeta, phoFill5x5SigmaIetaIeta, phoR9, pho_HoverE, pho_sumChargedHadronPtAllVertices, pho_sumChargedHadronPt, pho_sumNeutralHadronEt, 
  //pho_sumPhotonEt, pho_sumWorstVertexChargedHadronPt, pho_pfIsoChargedHadronIso, pho_pfIsoChargedHadronIsoWrongVtx, pho_pfIsoNeutralHadronIso, pho_pfIsoPhotonIso, pho_pfIsoModFrixione, 
  //pho_isConversion, pho_passEleVeto, pho_RegressionE, pho_RegressionEUncertainty, pho_IDMVA, 

  int NPU;
  //int nPV;
  unsigned int run, lumi, event;


  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("nPV", &nPV, "nPV/i");
  outputTree->Branch("mass", &mass, "mass/F");
  outputTree->Branch("t1", &t1, "t1/F");
  outputTree->Branch("t2", &t2, "t2/F");
  outputTree->Branch("t1_seed", &t1_seed, "t1_seed/F");
  outputTree->Branch("t2_seed", &t2_seed, "t2_seed/F");
  outputTree->Branch("t1calib_seed", &t1calib_seed, "t1calib_seed/F");
  outputTree->Branch("t2calib_seed", &t2calib_seed, "t2calib_seed/F");
  outputTree->Branch("t1raw_seed", &t1raw_seed, "t1raw_seed/F");
  outputTree->Branch("t2raw_seed", &t2raw_seed, "t2raw_seed/F");

  outputTree->Branch("phoNumber", &phoNumber, "phoNumber/F");
  outputTree->Branch("pho1Energy", &pho1Energy, "pho1Energy/F");
  outputTree->Branch("pho1_Pt", &pho1_Pt, "pho1_Pt/F");
  outputTree->Branch("pho1_Eta", &pho1_Eta, "pho1_Eta/F");
  outputTree->Branch("pho1_Phi", &pho1_Phi, "pho1_Phi/F");  

  outputTree->Branch("pho2Energy", &pho2Energy, "pho2Energy/F");
  outputTree->Branch("pho2_Pt", &pho2_Pt, "pho2_Pt/F");
  outputTree->Branch("pho2_Eta", &pho2_Eta, "pho2_Eta/F");
  outputTree->Branch("pho2_Phi", &pho2_Phi, "pho2_Phi/F");  

  outputTree->Branch("jetEnergy", &jetEnergy, "jetEnergy/F");
  outputTree->Branch("jet_Pt", &jet_Pt, "jet_Pt/F");
  outputTree->Branch("jet_Eta", &jet_Eta, "jet_Eta/F");
  outputTree->Branch("jet_Phi", &jet_Phi, "jet_Phi/F");

  outputTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
  outputTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
  outputTree->Branch("sum_MET", &sum_MET, "sum_MET/F");

  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);


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


    //initialize branches
    weight = 0;
    pileupWeight = 0;
    pileupWeightUp = 0;
    pileupWeightDown = 0;
    mass = 0;
    t1 = -999;
    t2 = -999;
    t1_seed = -999;
    t2_seed = -999;
    t1calib_seed = -999;
    t2calib_seed = -999;
    t1raw_seed = -999;
    t2raw_seed = -999;

    phoNumber = -999;
    pho1Energy = -999;
    pho1_Pt = -999;
    pho1_Eta = -999;
    pho1_Phi = -999;
    
    pho2Energy = -999;
    pho2_Pt = -999;
    pho2_Eta = -999;
    pho2_Phi = -999;

    jetEnergy = -999;
    jet_Pt = -999;
    jet_Eta = -999;
    jet_Phi = -999;

    met_Pt = -999;
    met_Phi = -999;
    sum_MET = -999;

    NPU = 0;   

    //fill normalization histogram
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;

    //get NPU
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
        NPU = nPUmean[i];
      }
    }

    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    
    double pvX = 0;

    int nPho = 0;
    TLorentzVector pho1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector pho2 = makeTLorentzVector(0,0,0,0);
    double pho1_time = 0;
    double pho2_time= 0;
    double pho1_seedtime = 0;
    double pho2_seedtime = 0;
    double pho1_seedtimeCalib = 0;
    double pho2_seedtimeCalib = 0;
    double pho1_seedtimeraw = 0;
    double pho2_seedtimeraw = 0;
    double MET_Pt_event = 0;
    double MET_Phi_event = 0;
    double MET_event = 0;

    for(int i = 0; i < nPhotons; i++) {
      // apply cuts
      if(phoPt[i] < 25) continue; 
      if(fabs(phoEta[i]) > 2.5) continue;
      if(fabs(phoEta[i]) > 1.4442 && fabs(phoEta[i]) < 1.566) continue; //the eta range for photon, this takes care of the gap between barrel and endcap
      //if(!(isEGammaPOGTightElectron(i))) continue;
      
      nPho++;
      TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], phoE[i]);
      
      //rough definition
      uint seedhitIndex =  (*pho_SeedRechitIndex)[i];

      cout<<"seedhitIndex: "<<seedhitIndex<<endl;
      cout<<"ecalRechit_ID size: "<<ecalRechit_ID->size()<<endl;
      cout<<"ecalRechit_ID: "<<(*ecalRechit_ID)[seedhitIndex]<<endl;

      bool isEBOrEE = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 ); //barrel vs. endcap
      double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      //apply intercalibration
      double IC_time_SeptRereco = getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double IC_time_LagacyRereco = getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double calibratedSeedHitTime = rawSeedHitTime + IC_time_LagacyRereco - IC_time_SeptRereco;

      //apply TOF correction
      double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
      

      double TOFrawSeedHitTime = rawSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;


      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << " : \n" 
      //  << "  runNum: " << runNum << "  detID: " << (*ecalRechit_ID)[seedhitIndex] << "  IC_time_SeptRereco: " << IC_time_SeptRereco << "  IC_time_LagacyRereco: " << IC_time_LagacyRereco << " \n"
      //    << " " << rawSeedHitTime << " -> " << calibratedSeedHitTime << " -> " << TOFCorrectedSeedHitTime << " "
      //  << "\n";
 

      double tmpSumWeightedTime = 0;
      double tmpSumWeight = 0;

      for (uint k=0; k<(*pho_EcalRechitIndex)[i].size(); ++k) {

	cout << metPt << endl;
        
        uint rechitIndex = (*pho_EcalRechitIndex)[i][k];

	MET_event = sumMET;

        MET_Pt_event = metPt;
        MET_Phi_event = metPhi;
      
        double rawT = (*ecalRechit_T)[rechitIndex];
        //apply intercalibration
        //apply TOF correction
        double corrT = rawT + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

        // double pedNoise = getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex]);
        double pedNoise = 1;
        double ADCToGeV = getADCToGeV(runNum, isEBOrEE);
        double sigmaE = pedNoise * ADCToGeV;
  
        float C_EB = 1;
        float N_EB = 1;
  
        double sigmaT = N_EB / ((*ecalRechit_E)[rechitIndex] / sigmaE) + sqrt(2) * C_EB;
        tmpSumWeightedTime += corrT * ( 1.0 / (sigmaT*sigmaT) );
        tmpSumWeight += ( 1.0 / (sigmaT*sigmaT) );
        // cout << "\n";
      }
      double weightedTime = tmpSumWeightedTime / tmpSumWeight;

            
      if (thisPhoton.Pt() > pho1.Pt()) {
        pho1 = thisPhoton;
        pho1_time = weightedTime;
        //pho1_seedtime = TOFCorrectedSeedHitTime;
        pho1_seedtime = TOFrawSeedHitTime;
        pho1_seedtimeCalib = calibratedSeedHitTime;
        pho1_seedtimeraw = rawSeedHitTime;
        pho1IsEB = bool( phoEta_SC[i] < 1.5 );
        if (pho1IsEB) { //if the photon is in the barrel
        pho1SeedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        pho1SeedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        pho1SeedIX = -999;
        pho1SeedIY = -999;
      } 
      else {
        pho1SeedIEta = -999;
        pho1SeedIPhi = -999;
        pho1SeedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
        pho1SeedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
      }
    } 
    else if (thisPhoton.Pt() > pho2.Pt()) {
      pho2 = thisPhoton;
      pho2_time = weightedTime;
      //pho2_seedtime = TOFCorrectedSeedHitTime; 
      pho2_seedtime = TOFrawSeedHitTime; 
      pho2_seedtimeCalib = calibratedSeedHitTime;
      pho2_seedtimeraw = rawSeedHitTime;
      pho2IsEB = bool( phoEta_SC[i] < 1.5 );
      if (pho2IsEB) {
        pho2SeedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        pho2SeedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        pho2SeedIX = -999;
        pho2SeedIY = -999;
      } 
      else {
        pho2SeedIEta = -999;
        pho2SeedIPhi = -999;
        pho2SeedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
        pho2SeedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
      }
    } 
  }
    
  if (nPho >= 2) {
    pho1Energy = pho1.E();
    pho1_Pt = pho1.Pt();
    pho1_Eta = pho1.Eta();
    pho1_Phi = pho1.Phi();
    pho2Energy = pho2.E();
    pho2_Pt = pho2.Pt();
    pho2_Eta = pho2.Eta();
    pho2_Phi = pho2.Phi();

    mass = (pho1+pho2).M();
    t1 = pho1_time;
    t2 = pho2_time;
    t1_seed = pho1_seedtime;
    t2_seed = pho2_seedtime;
    t1calib_seed = pho1_seedtimeCalib;
    t2calib_seed = pho2_seedtimeCalib;
    t1raw_seed = pho1_seedtimeraw;
    t2raw_seed = pho2_seedtimeraw; 
    met_Pt = MET_Pt_event;
    met_Phi =  MET_Phi_event;
    sum_MET = MET_event;

    //cout << "ele2: " << ele2.Pt() << " " << ele2_seedtime << "\n";
  }

  //Fill Event
  //if (mass > 60 && mass < 120) {
  if ( t1_seed > 0 && t2_seed > 0) {
    outputTree->Fill();
  }

}//end of event loop
  
cout << "Writing output trees..." << endl;
outputTree->Write();
NEvents->Write();

outFile->Close();
}
