#include "DelayedPhotonAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"
#include "TLorentzVector.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns


TLorentzVector DelayedPhotonAnalyzer::intersectPoint(float x0,float y0,float z0,float px,float py,float pz,float R)
{
  TLorentzVector sol;


  float x1,y1,z1;
  float x2,y2,z2;

  if((px*px*py*py*R*R + py*py*py*py*R*R - py*py*py*py*x0*x0 + 2.0*px*py*py*py*x0*y0 - px*px*py*py*y0*y0) > 0.0 )
  {
  y1 = (-px*py*x0 + px*px*y0)/(px*px+py*py) - pow((px*px*py*py*R*R + py*py*py*py*R*R - py*py*py*py*x0*x0 + 2.0*px*py*py*py*x0*y0 - px*px*py*py*y0*y0),0.5) / (px*px+py*py);
  x1 =  x0 + (px/py) * (y1-y0);
  z1 =  z0 + (pz/py) * (y1-y0);

  y2 = (-px*py*x0 + px*px*y0)/(px*px+py*py) + pow((px*px*py*py*R*R + py*py*py*py*R*R - py*py*py*py*x0*x0 + 2.0*px*py*py*py*x0*y0 - px*px*py*py*y0*y0),0.5) / (px*px+py*py); 
  x2 =  x0 + (px/py) * (y2-y0);
  z2 =  z0 + (pz/py) * (y2-y0);
  }
 
  else
  {
   x1=0,y1=0,z1=0,x2=0,y2=0,z2=0; 
  }


  if( (z1-z0)*pz > 0.0 ) 
  {
	sol.SetXYZM(x1,y1,z1,0.0);
  }

  else if( (z2-z0)*pz > 0.0 )
  {
	sol.SetXYZM(x2,y2,z2,0.0);
  }

  else
  {
	sol.SetXYZM(0,0,0,0.0);
  }	
//  cout<<"solution: "<<sol.X()<<"  "<<sol.Y()<<"  "<<sol.Z()<<endl;
  return sol;
};

float DelayedPhotonAnalyzer::getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID) {
  float timeCalib = 1.0;

  // accessing variables for the seed hit information
  std::vector<float> *ecalRechit_X;
  std::vector<float> *ecalRechit_Y;
  std::vector<float> *ecalRechit_Z;
  
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
  float t1_seed_genV, t2_seed_genV;
  float TOF_total1;
  float TOF_total2;
  float TOF_total1_genV;
  float TOF_total2_genV;
  float TOF_initial;
  float TOF_neu1, TOF_neu2;
  float TOF_neu1_RF, TOF_neu2_RF;
  float TOF_pho1, TOF_pho2;
  //float TOF_pho1_RF, TOF_pho2_RF;
  float t1calib_seed, t2calib_seed;
  float t1raw_seed, t2raw_seed;
  float phoNumber;
  float pho1Energy, pho1_Pt, pho1_Eta, pho1_Phi;
  float pho2Energy, pho2_Pt, pho2_Eta, pho2_Phi;
  float jetEnergy, jet_Pt, jet_Eta, jet_Phi;
  float met_Pt, met_Phi, sum_MET;
  float diff_gpt1, diff_gpt2;
  float deltaR_gpt1, deltaR_gpt2;

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
  outputTree->Branch("pvX", &pvX, "pvX/F");
  outputTree->Branch("pvY", &pvY, "pvY/F");
  outputTree->Branch("pvZ", &pvZ, "pvZ/F");
  outputTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  outputTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  outputTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  outputTree->Branch("t1", &t1, "t1/F");
  outputTree->Branch("t2", &t2, "t2/F");
  outputTree->Branch("t1_seed", &t1_seed, "t1_seed/F");
  outputTree->Branch("t2_seed", &t2_seed, "t2_seed/F");
  outputTree->Branch("t1_seed_genV", &t1_seed_genV, "t1_seed_genV/F");
  outputTree->Branch("t2_seed_genV", &t2_seed_genV, "t2_seed_genV/F");
  outputTree->Branch("t1calib_seed", &t1calib_seed, "t1calib_seed/F");
  outputTree->Branch("t2calib_seed", &t2calib_seed, "t2calib_seed/F");
  outputTree->Branch("t1raw_seed", &t1raw_seed, "t1raw_seed/F");
  outputTree->Branch("t2raw_seed", &t2raw_seed, "t2raw_seed/F");
  outputTree->Branch("TOF_total1", &TOF_total1, "TOF_total1/F");
  outputTree->Branch("TOF_total2", &TOF_total2, "TOF_total2/F");
  outputTree->Branch("TOF_total1_genV", &TOF_total1_genV, "TOF_total1_genV/F");
  outputTree->Branch("TOF_total2_genV", &TOF_total2_genV, "TOF_total2_genV/F");
  outputTree->Branch("TOF_initial", &TOF_initial, "TOF_initial/F");
  outputTree->Branch("TOF_neu1", &TOF_neu1, "TOF_neu1/F");
  outputTree->Branch("TOF_neu1_RF", &TOF_neu1_RF, "TOF_neu1_RF/F");
  outputTree->Branch("TOF_neu2", &TOF_neu2, "TOF_neu2/F");
  outputTree->Branch("TOF_neu2_RF", &TOF_neu2_RF, "TOF_neu2_RF/F");
  outputTree->Branch("TOF_pho1", &TOF_pho1, "TOF_pho1/F");
  //outputTree->Branch("TOF_pho1_RF", &TOF_pho1_RF, "TOF_pho1_RF/F");
  outputTree->Branch("TOF_pho2", &TOF_pho2, "TOF_pho2/F");
  //outputTree->Branch("TOF_pho2_RF", &TOF_pho2_RF, "TOF_pho2_RF/F");

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

  outputTree->Branch("diff_gpt1", &diff_gpt1, "diff_gpt1/F");
  outputTree->Branch("diff_gpt2", &diff_gpt2, "diff_gpt2/F");
  outputTree->Branch("deltaR_gpt1", &deltaR_gpt1, "deltaR_gpt1/F");
  outputTree->Branch("deltaR_gpt2", &deltaR_gpt2, "deltaR_gpt2/F");

  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  float matched = 0;
  float unmatched = 0;


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
    t1_seed_genV = -999;
    t2_seed_genV = -999;
    t1calib_seed = -999;
    t2calib_seed = -999;
    t1raw_seed = -999;
    t2raw_seed = -999;

    TOF_total1 = -999;
    TOF_total2 = -999;
    TOF_total1_genV = -999;
    TOF_total2_genV = -999;
    TOF_initial = 0; // currently set to 0 to assume gen vertex is at 0
    TOF_neu1 = -999;
    TOF_neu2 = -999;
    TOF_pho1 = -999;
    TOF_pho2 = -999;

    TOF_neu1_RF = -999;
    TOF_neu2_RF = -999;
    //TOF_pho1_RF = -999;
    //TOF_pho2_RF = -999;


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

    diff_gpt1 = -999;
    diff_gpt2 = -999;
    deltaR_gpt1 = -999;
    deltaR_gpt2 = -999;

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
    
    //double pvX = 0;

    int nPho = 0;
    TLorentzVector pho1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector pho2 = makeTLorentzVector(0,0,0,0);
    double pho1_time = 0;
    double pho2_time= 0;
    double pho1_seedtime = 0;
    double pho2_seedtime = 0;
    double pho1_seedtime_genV = 0;
    double pho2_seedtime_genV = 0;
    double pho1_seedtime_recoV = 0;
    double pho2_seedtime_recoV = 0;
    double pho1_seedtimeCalib = 0;
    double pho2_seedtimeCalib = 0;
    double pho1_seedtimeraw = 0;
    double pho2_seedtimeraw = 0;
    double MET_Pt_event = 0;
    double MET_Phi_event = 0;
    double MET_event = 0;

    // XYZ rechit where photon is detected
    float pho1SeedX = 0;
    float pho1SeedY = 0;
    float pho1SeedZ = 0;
    float pho2SeedX = 0;
    float pho2SeedY = 0;
    float pho2SeedZ = 0;


    for(int i = 0; i < nPhotons; i++) { //photon loop
      // apply cuts
//      if(phoPt[i] < 25) continue; 
//      if(fabs(phoEta[i]) > 2.5) continue;
//      if(fabs(phoEta[i]) > 1.4442 && fabs(phoEta[i]) < 1.566) continue; //the eta range for photon, this takes care of the gap between barrel and endcap
      //if(!(isEGammaPOGTightElectron(i))) continue;
      
      nPho++;
      TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], phoE[i]);
      
      //rough definition
      uint seedhitIndex =  (*pho_SeedRechitIndex)[i];
    
      cout<<"reco Photon - "<<i<<" : seedX = "<<(*ecalRechit_X)[seedhitIndex]<<" : seedY = "<<(*ecalRechit_Y)[seedhitIndex]<<" : seedZ = "<<(*ecalRechit_Z)[seedhitIndex]<<"  pT = "<<phoPt[i]<<"  Energy = "<<phoE[i]<<endl;

      //cout<<"seedhitIndex: "<<seedhitIndex<<endl;
      //cout<<"ecalRechit_ID size: "<<ecalRechit_ID->size()<<endl;
      //cout<<"ecalRechit_ID: "<<(*ecalRechit_ID)[seedhitIndex]<<endl;

      bool isEBOrEE = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 ); //barrel vs. endcap
      double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      //apply intercalibration
      double IC_time_SeptRereco = getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double IC_time_LagacyRereco = getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double calibratedSeedHitTime = rawSeedHitTime + IC_time_LagacyRereco - IC_time_SeptRereco;

      //apply TOF correction
      double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
      

      double TOFrawSeedHitTime = rawSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT; //this is TOF correction with reco V

      //generator xyz information
      float primary_x1;
      float primary_y1;
      float primary_z1;

      primary_x1 = genVertexX;
      primary_y1 = genVertexY; 
      primary_z1 = genVertexZ;

      double TOFrawSeedHitTime_genV = rawSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-primary_x1,2)+pow((*ecalRechit_Y)[seedhitIndex]-primary_y1,2)+pow((*ecalRechit_Z)[seedhitIndex]-primary_z1,2)))/SPEED_OF_LIGHT;

      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << " : \n" 
      //  << "  runNum: " << runNum << "  detID: " << (*ecalRechit_ID)[seedhitIndex] << "  IC_time_SeptRereco: " << IC_time_SeptRereco << "  IC_time_LagacyRereco: " << IC_time_LagacyRereco << " \n"
      //    << " " << rawSeedHitTime << " -> " << calibratedSeedHitTime << " -> " << TOFCorrectedSeedHitTime << " "
      //  << "\n";

      double tmpSumWeightedTime = 0;
      double tmpSumWeight = 0;

      for (uint k=0; k<(*pho_EcalRechitIndex)[i].size(); ++k) {

	//cout << metPt << endl;
        
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

        // finds the rechit position of the first photon
  	pho1SeedX = (*ecalRechit_X)[seedhitIndex];
  	pho1SeedY = (*ecalRechit_Y)[seedhitIndex];
  	pho1SeedZ = (*ecalRechit_Z)[seedhitIndex];

        pho1_time = weightedTime;
        //pho1_seedtime = TOFCorrectedSeedHitTime;
        pho1_seedtime = TOFrawSeedHitTime;
        pho1_seedtime_genV = TOFrawSeedHitTime_genV;
        pho1_seedtimeCalib = calibratedSeedHitTime;
        pho1_seedtimeraw = rawSeedHitTime;
        pho1IsEB = bool( phoEta_SC[i] < 1.5 );
	if (pho1IsEB) { //if the photon is in the barrel
        	pho1SeedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        	pho1SeedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
        	pho1SeedIX = -999;
        	pho1SeedIY = -999;

		//cout << "rawSeedHitTime:  " << rawSeedHitTime << "  TOFrawSeedHitTime:  " << TOFrawSeedHitTime << endl;
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
      pho2SeedX = (*ecalRechit_X)[seedhitIndex];
      pho2SeedY = (*ecalRechit_Y)[seedhitIndex];
      pho2SeedZ = (*ecalRechit_Z)[seedhitIndex];

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
  } //end photon loop

      
  if (nPho >= 2) { 
    cout << "THIS IS THE 2 PHOTON LOOP" << endl;

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
    t1_seed_genV = pho1_seedtime_genV;
    t2_seed_genV = pho2_seedtime_genV;
    t1calib_seed = pho1_seedtimeCalib;
    t2calib_seed = pho2_seedtimeCalib;
    t1raw_seed = pho1_seedtimeraw;
    t2raw_seed = pho2_seedtimeraw; 
    met_Pt = MET_Pt_event;
    met_Phi =  MET_Phi_event;
    sum_MET = MET_event;

    //cout << "ele2: " << ele2.Pt() << " " << ele2_seedtime << "\n";


    //TOF of neutralino and photon
    bool foundN1 = 0;
    bool foundN2 = 0; 
    int pho1_index = 0;
    int pho2_index = 0;
    int neu1_index = 0;
    int neu2_index = 0;

    // finding the neutralino and photon index
    for(unsigned int i = 0; i < nGenParticle; i++){ //this is gen particle loop within event and photon loop
	//cout << "GEN PARTICLE LOOP NOW " << endl;
	//if ( foundN1 == 0 && gParticleId[i] == 1000022 &&  gParticleMotherIndex[i]==0){ // neutralino is first particle
	if ( foundN1 == 0 && gParticleId[i] == 1000022 &&  gParticleMotherId[i]==2000011 ){ // neutralino comes from slepton
		cout << "FOUND NEUTRALINO" << endl;
		cout << "particle ID   "<<gParticleId[i]<<"  mother ID  "<<gParticleMotherId[i] <<"  mother index  "<<gParticleMotherIndex[i] << endl;
        	neu1_index = i;
        	foundN1 = 1;
        	for(unsigned int j = 0; j < nGenParticle; j++){
		  	if ( gParticleId[j] == 22 && gParticleMotherIndex[j] == neu1_index ){
        	        pho1_index = j;
        	        //cout << "neutralino 1 index  "<<neu1_index<< "  photon 1 index  "<<pho1_index<<endl;
        	        }
        	}
        }
 	//else if ( foundN1 == 1 && foundN2 == 0 && gParticleId[i] == 1000022 &&  gParticleMotherIndex[i]==0 ){
 	else if ( foundN1 == 1 && foundN2 == 0 && gParticleId[i] == 1000022 &&  gParticleMotherId[i]==2000011 ){
        	neu2_index = i;
        	foundN2 = 1;
        	for(unsigned int j = 0; j < nGenParticle; j++){
			if ( gParticleId[j] == 22 && gParticleMotherIndex[j] == neu2_index ){
                	pho2_index = j;
                	//cout << "neutralino 2 index  "<<neu2_index<< "  photon 2 index  "<<pho2_index<<endl;
                	}
          	}
       	}
  }



  if(foundN1==1 && foundN2==1){

  float decay_x1;
  float decay_y1;
  float decay_z1;
  float decay_x2;
  float decay_y2;
  float decay_z2;

  // this is where the neutralino decays (photon is created)
  decay_x1 = gParticleDecayVertexX[neu1_index];
  decay_y1 = gParticleDecayVertexY[neu1_index];
  decay_z1 = gParticleDecayVertexZ[neu1_index];

  decay_x2 = gParticleDecayVertexX[neu2_index];
  decay_y2 = gParticleDecayVertexY[neu2_index];
  decay_z2 = gParticleDecayVertexZ[neu2_index];


  // need to match up the photon index and the reco photon - this is done based on momentum
  // pho1_Pt is reco level, gpho1_Pt is gen level information

  float gpho1_Pt = gParticlePt[pho1_index];
  float gpho2_Pt = gParticlePt[pho2_index];
 
  float diffpt11 = fabs(pho1_Pt - gpho1_Pt);
  float diffpt21 = fabs(pho1_Pt - gpho2_Pt);
  float diffpt12 = fabs(pho2_Pt - gpho1_Pt);
  float diffpt22 = fabs(pho2_Pt - gpho2_Pt);
  

  TLorentzVector genSeed1;
  TLorentzVector genSeed2;
  genSeed1 = intersectPoint(decay_x1, decay_y1, decay_z1, gParticlePx[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5), gParticlePy[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5), gParticlePz[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5),129.0); // using intersection function written above, radius as 129 cm


  cout<<"genSeed1: input : "<<decay_x1<<"   "<<decay_y1<<"   "<<decay_z1<<"   "<<
	gParticlePx[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5)<<"  "<<
	gParticlePy[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5)<<"   "<<
	gParticlePz[pho1_index]/pow(gParticlePx[pho1_index]*gParticlePx[pho1_index] + gParticlePy[pho1_index]+gParticlePy[pho1_index]+gParticlePz[pho1_index]*gParticlePz[pho1_index],0.5)<<"   "<<
	"  output "<<genSeed1.X()<<"  "<<genSeed1.Y()<<"   "<<genSeed1.Z()<<endl;

  genSeed2 = intersectPoint(decay_x2, decay_y2, decay_z2, gParticlePx[pho2_index]/pow(gParticlePx[pho2_index]*gParticlePx[pho2_index] + gParticlePy[pho2_index]+gParticlePy[pho2_index]+gParticlePz[pho2_index]*gParticlePz[pho2_index],0.5), gParticlePy[pho2_index]/pow(gParticlePx[pho2_index]*gParticlePx[pho2_index] + gParticlePy[pho2_index]+gParticlePy[pho2_index]+gParticlePz[pho2_index]*gParticlePz[pho2_index],0.5), gParticlePz[pho2_index]/pow(gParticlePx[pho2_index]*gParticlePx[pho2_index] + gParticlePy[pho2_index]+gParticlePy[pho2_index]+gParticlePz[pho2_index]*gParticlePz[pho2_index],0.5),129.0);

 
  
  TLorentzVector recoPho1LV, recoPho2LV, genPho1LV, genPho2LV;
  //recoPho1LV.SetPtEtaPhiE(pho1_Pt, pho1_Eta,pho1_Phi,pho1Energy);
  //recoPho2LV.SetPtEtaPhiE(pho2_Pt, pho2_Eta,pho2_Phi,pho2Energy);
  recoPho1LV.SetXYZM(pho1SeedX-decay_x1,pho1SeedY-decay_y1, pho1SeedZ-decay_z1 ,0.0);
  recoPho2LV.SetXYZM(pho2SeedX-decay_x2,pho2SeedY-decay_y2, pho2SeedZ-decay_z2 ,0.0);
  genPho1LV.SetPxPyPzE(gParticlePx[pho1_index], gParticlePy[pho1_index], gParticlePz[pho1_index], gParticleE[pho1_index]);
  genPho2LV.SetPxPyPzE(gParticlePx[pho2_index], gParticlePy[pho2_index], gParticlePz[pho2_index], gParticleE[pho2_index]);

  //float gpho1_Eta = gParticleEta[pho1_index];
  //float gpho2_Eta = gParticleEta[pho2_index];
  //float gpho1_Phi = gParticlePhi[pho1_index];
  //float gpho2_Phi = gParticlePhi[pho2_index];

  //float diff11_Eta = pho1_Eta - gpho1_Eta;
  //float diff12_Eta = pho2_Eta - gpho1_Eta;
  //float diff11_Phi = pho1_Phi - gpho1_Phi;
  //float diff12_Phi = pho2_Phi - gpho1_Phi;

  cout <<recoPho1LV.Eta() << "  reco1 eta  "<<genSeed1.Eta()<<"  gen1 eta  "<<recoPho2LV.Eta() << "  reco2 eta  "<<genSeed2.Eta()<<"  gen2 eta  "<< endl;
  cout <<recoPho1LV.Phi() << "  reco1 phi  "<<genSeed1.Phi()<<"  gen1 phi  "<<recoPho2LV.Phi() << "  reco2 eta  "<<genSeed2.Phi()<<"  gen2 phi  "<< endl;
  cout <<genSeed1.X()<<"  gen1 X  "<<genSeed1.Y()<<"  gen1 Y  "<<genSeed1.Z()<<"  gen1 Z  "<< endl;
  cout <<genSeed2.X()<<"  gen2 X  "<<genSeed2.Y()<<"  gen2 Y  "<<genSeed2.Z()<<"  gen2 Z  "<< endl;

  float diffR11 = pow(pow(recoPho1LV.Eta()-genSeed1.Eta(),2.0)+pow(recoPho1LV.Phi()-genSeed1.Phi(),2.0),0.5);
  float diffR21 = pow(pow(recoPho1LV.Eta()-genSeed2.Eta(),2.0)+pow(recoPho1LV.Phi()-genSeed2.Phi(),2.0),0.5);
  float diffR12 = pow(pow(recoPho2LV.Eta()-genSeed1.Eta(),2.0)+pow(recoPho2LV.Phi()-genSeed1.Phi(),2.0),0.5); 
  float diffR22 = pow(pow(recoPho2LV.Eta()-genSeed2.Eta(),2.0)+pow(recoPho2LV.Phi()-genSeed2.Phi(),2.0),0.5); 
 
  float distanceX_1p;
  float distanceX_2p;
  float distanceY_1p;
  float distanceY_2p;
  float distanceZ_1p;
  float distanceZ_2p;

  float distance1_p;
  float distance2_p;

  cout << "diffR11  "<<diffR11<<"  diffR12  "<<diffR12<<"  diffR21  "<<diffR21<<"  diffR22  "<<diffR22<<endl;

  if ( diffR11 < diffR12 && 0 < diffR11 && diffR11 <0.400 && 0 < diffR22 && diffR22 < 0.400 ){ // cut made on diffRerence between momenta
	cout << "passed diff cuts"<< endl;
  //if ( diff11_Eta < diff12_Eta && diff11_Phi < diff12_Phi ) {
  	// then gpho1 and pho1 are matched up - find distance photon travels
        distanceX_1p = pho1SeedX - decay_x1;
  	distanceY_1p = pho1SeedY - decay_y1;
  	distanceZ_1p = pho1SeedZ - decay_z1;
        distanceX_2p = pho2SeedX - decay_x2;
  	distanceY_2p = pho2SeedY - decay_y2;
  	distanceZ_2p = pho2SeedZ - decay_z2;
	diff_gpt1 = diffpt11;
	diff_gpt2 = diffpt22;

	deltaR_gpt1 = diffR11;//pow(pow(recoPho1LV.Eta()-genSeed1.Eta(),2.0)+pow(recoPho1LV.Phi()-genSeed1.Phi(),2.0),0.5);
	deltaR_gpt2 = diffR22;//pow(pow(recoPho2LV.Eta()-genSeed2.Eta(),2.0)+pow(recoPho2LV.Phi()-genSeed2.Phi(),2.0),0.5);

	distance1_p = pow( (pow(distanceX_1p,2) + pow(distanceY_1p,2) + pow(distanceZ_1p,2)),0.5) ;
	distance2_p = pow( (pow(distanceX_2p,2) + pow(distanceY_2p,2) + pow(distanceZ_2p,2)),0.5) ;
	matched += 1;
  }
  else if ( diffR12 < diffR11 && 0 < diffR12 && diffR12 < 0.100 && 0 < diffR21 && diffR21 <0.100 ) {
  //if ( diff11_Eta > diff12_Eta && diff11_Phi > diff12_Phi ) {
  	// then gpho1 and pho2 are matched up - find distance photon travels
        distanceX_1p = pho2SeedX - decay_x1;
  	distanceY_1p = pho2SeedY - decay_y1;
  	distanceZ_1p = pho2SeedZ - decay_z1;
        distanceX_2p = pho1SeedX - decay_x2;
  	distanceY_2p = pho1SeedY - decay_y2;
  	distanceZ_2p = pho1SeedZ - decay_z2;
	diff_gpt1 = diffpt12;
	diff_gpt2 = diffpt21;

        deltaR_gpt1 = diffR12;//pow(pow(recoPho2LV.Eta()-genSeed1.Eta(),2.0)+pow(recoPho2LV.Phi()-genSeed1.Phi(),2.0),0.5);
	deltaR_gpt2 = diffR21;//pow(pow(recoPho1LV.Eta()-genSeed2.Eta(),2.0)+pow(recoPho1LV.Phi()-genSeed2.Phi(),2.0),0.5);


	distance1_p = pow( pow(distanceX_1p,2) + pow(distanceY_1p,2) + pow(distanceZ_1p,2),0.5) ;
	distance2_p = pow( pow(distanceX_2p,2) + pow(distanceY_2p,2) + pow(distanceZ_2p,2),0.5) ;
	matched += 1;
  }
  else {
	unmatched += 1;  // find efficiency of matching gen photons
  }


  if (diff_gpt1 == -999 && diff_gpt2 == -999 ) {
	unmatched += 1;
  }
  //else {
	// += 1;
  //}
  //if (diff_gpt2 == -999) {
	//unmatched += 1;
  //}
  //else {
	//matched += 1;
  //}
  //cout << "matched:  " << matched << "  unmatched:  " << unmatched << endl;

  //const float SPEED_OF_LIGHT = 29.9792; //cm/ns

  // time of flight calculations for the photon
  // find the total distance the photon traveled
  //float distance1_p = pow( pow(distanceX_1p,2) + pow(distanceY_1p,2) + pow(distanceZ_1p,2),0.5) ;
  //float distance2_p = pow( pow(distanceX_2p,2) + pow(distanceY_2p,2) + pow(distanceZ_2p,2),0.5) ;

  float time1_p = distance1_p / (SPEED_OF_LIGHT) ;
  float time2_p = distance2_p / (SPEED_OF_LIGHT) ;

  if ( 0 <= time1_p ) {
  	TOF_pho1 = time1_p; // fill the ROOT branch histogram with the TOF of first photon
  }
  if ( 0 <= time2_p ) {
  	TOF_pho2 = time2_p; // fill the ROOT branch histogram with the TOF of first photon
  }


  float primary_x1;
  float primary_y1;
  float primary_z1;
  float NeutraPt1;
  float NeutraPt2;
  float NeutraEta1;
  float NeutraEta2;

  float creation_time; // this will be genVertexT, time the neutralino is created

  primary_x1 = genVertexX;
  primary_y1 = genVertexY; 
  primary_z1 = genVertexZ;
  NeutraPt1 = gParticlePt[neu1_index]; // transverse momentum of neutralino particles when created
  NeutraPt2 = gParticlePt[neu2_index]; // transverse momentum of neutralino particles when created
  NeutraEta1 = gParticleEta[neu1_index]; // eta value used to find momentum
  NeutraEta2 = gParticleEta[neu2_index]; // eta value used to find momentum

  //cout << genVertexT << "   generated time vertex" << endl;
  //cout << genVertexZ << "   generated vertex Z" << endl;

  //creation_time = genVertexT;

  // creation time is added to neutralino TOF and photon TOF to find total hit time
  //TOF_initial = creation_time;
   
  //cout << TOF_initial << "   TOF_initial value" << endl;
  
  float distanceX1;
  float distanceX2;
  float distanceY1;
  float distanceY2;
  float distanceZ1;
  float distanceZ2;

  // distance calculations for the neutralino
  distanceX1 = decay_x1 - primary_x1;
  distanceX2 = decay_x2 - primary_x1;

  distanceY1 = decay_y1 - primary_y1;
  distanceY2 = decay_y2 - primary_y1;

  distanceZ1 = decay_z1 - primary_z1;
  distanceZ2 = decay_z2 - primary_z1;

  // now do the calculation for the time of flight for the neutralino
  // find the total distance the particle traveled
  float distance1 = pow( (pow(distanceX1,2) + pow(distanceY1,2) + pow(distanceZ1,2)),0.5) ;
  float distance2 = pow( (pow(distanceX2,2) + pow(distanceY2,2) + pow(distanceZ2,2)),0.5) ;
  // find the momentum from the transverse momentum and the eta value
  float momentum1 = NeutraPt1 * cosh(NeutraEta1);
  float momentum2 = NeutraPt2 * cosh(NeutraEta2);

  float mass = 1000; // in GeV

  float time1 = distance1 / (SPEED_OF_LIGHT*momentum1) * pow((pow(mass,2) + pow(momentum1,2)),0.5);
  float time2 = distance2 / (SPEED_OF_LIGHT*momentum2) * pow((pow(mass,2) + pow(momentum2,2)),0.5);

  TOF_neu1 = time1; // fill the ROOT branch histogram with the TOF of first neutralino
  TOF_neu1_RF = time1*mass*pow((pow(mass,2) + pow(momentum1,2)),-0.5);
  TOF_total1 = TOF_initial + TOF_neu1 + TOF_pho1 - pow(pho1SeedX*pho1SeedX + pho1SeedY*pho1SeedY + pho1SeedZ*pho1SeedZ, 0.5) / SPEED_OF_LIGHT; // this is with 000 vertex

  float deltax = pho1SeedX - primary_x1;
  float deltay = pho1SeedY - primary_y1;
  float deltaz = pho1SeedZ - primary_z1;

  TOF_total1_genV = TOF_initial + TOF_neu1 + TOF_pho1 - pow((pho1SeedX-primary_x1)*(pho1SeedX-primary_x1) + (pho1SeedY-primary_y1)*(pho1SeedY-primary_y1) + (pho1SeedZ-primary_z1)*(pho1SeedZ-primary_z1), 0.5) / SPEED_OF_LIGHT; 
  float L_gen_crys1 = pow( (pow(deltax,2) + pow(deltay,2) + pow(deltaz,2) ), 0.5);


  TOF_neu2 = time2; // fill the ROOT branch histogram with the TOF of first neutralino
  TOF_neu2_RF = time2*mass*pow((pow(mass,2) + pow(momentum2,2)),-0.5);
  TOF_total2 = TOF_initial + TOF_neu2 + TOF_pho2 - pow(pho2SeedX*pho2SeedX + pho2SeedY*pho2SeedY + pho2SeedZ*pho2SeedZ, 0.5) / SPEED_OF_LIGHT; 
	
  float deltax2 = pho2SeedX - primary_x1;
  float deltay2 = pho2SeedY - primary_y1;
  float deltaz2 = pho2SeedZ - primary_z1;

  float L_gen_crys2 = pow( (pow(deltax2,2) + pow(deltay2,2) + pow(deltaz2,2) ), 0.5);
  TOF_total2_genV = TOF_initial + TOF_neu2 + TOF_pho2 - pow((pho2SeedX-primary_x1)*(pho2SeedX-primary_x1) + (pho2SeedY-primary_y1)*(pho2SeedY-primary_y1) + (pho2SeedZ-primary_z1)*(pho2SeedZ-primary_z1), 0.5) / SPEED_OF_LIGHT; 
      


  if( diffR11 < diffR12 && 0 < diffR11 && diffR11 < 0.100 && 0 < diffR22 && diffR22 < 0.100 ) { // pt matching, photon1 and photon1
	cout <<"EVENT -GEN  genSeedXYZ: "<<genSeed1.X() <<"  "<< genSeed1.Y() <<"  "<<genSeed1.Z() <<"   deltaR(gen,reco): "<<deltaR_gpt1<<" genPhoPt: "<<gpho1_Pt<<" pT diff  "<< diffpt11 <<"  genPhoPx: "<<gParticlePx[pho1_index]<< "  genPhoPy: "<<gParticlePy[pho1_index]<< "  genPhoPz: "<<gParticlePz[pho1_index]<<"  genPhoE: "<<gParticleE[pho1_index]<<"    pho1 - gen vertex:  " << creation_time << "  neutralino dist:  "<<distance1<< "  neutralino TOF:  " << time1 << "  photon dist:  "<<distance1_p<< "  photon TOF:  " << time1_p << "  initialTOF:  " << TOF_initial << "  TOF_neu1:  " << TOF_neu1 << "  TOF_pho1:  " << TOF_pho1 << "  TOF_total1:  "<< TOF_initial + TOF_neu1 + TOF_pho1 <<"  TOF_total1-L/c:  " << TOF_total1 <<"  TOF_total1-L/c (gen):  " << TOF_total1_genV <<"  neutralino decay x (gen):  "<<decay_x1<< "  neutralino decay y (gen):  "<<decay_y1<< "  neutralino decay z (gen):  "<<decay_z1<< "  Length gen to crystal:  " <<L_gen_crys1 <<"  RECO  seedHitRawTime:   "<<t1raw_seed <<"  seedHitTimeTOF: "<<t1_seed<<"  seedHitTimeTOF genV: "<<t1_seed_genV<< "  pho1SeedX:  "<<pho1SeedX<<"  pho1SeedY:  "<<pho1SeedY<<"  pho1SeedZ:  "<<pho1SeedZ<< endl; 
	cout <<"EVENT -GEN  genSeedXYZ: "<<genSeed2.X() <<"  "<< genSeed2.Y() <<"  "<<genSeed2.Z() <<"   deltaR(gen,reco): "<<deltaR_gpt2<<" genPhoPt: "<<gpho2_Pt<<" pT diff  "<< diffpt22 <<"  genPhoPx: "<<gParticlePx[pho2_index]<< "  genPhoPy: "<<gParticlePy[pho2_index]<< "  genPhoPz: "<<gParticlePz[pho2_index]<< "  genPhoE: "<<gParticleE[pho2_index]<<"   pho2 - gen vertex:  " << creation_time << "  neutralino dist:  "<<distance2<< "  neutralino TOF:  " << time2 << "  photon dist:  "<<distance2_p<< "  photon TOF:  " << time2_p << "  initialTOF:  " << TOF_initial << "  TOF_neu2:  " << TOF_neu2 << "  TOF_pho2:  " << TOF_pho2 << "TOF_total2:  "<< TOF_initial + TOF_neu2 + TOF_pho2 <<"  TOF_total2-L/c:  " << TOF_total2 <<"  TOF_total2-L/c (gen):  " << TOF_total2_genV <<"  neutralino decay x (gen):  "<<decay_x2<< "  neutralino decay y (gen):  "<<decay_y2<< "  neutralino decay z (gen):  "<<decay_z2 << "  Length gen to crystal:  " <<L_gen_crys2<<"  RECO  seedHitRawTime:   "<<t2raw_seed <<"  seedHitTimeTOF: "<<t2_seed<<"  seedHitTimeTOF genV: "<<t2_seed_genV<< "  pho2SeedX:  "<<pho2SeedX<<"  pho2SeedY:  "<<pho2SeedY<<"  pho2SeedZ:  "<<pho2SeedZ<< endl; 
  }

  else if ( diffR12 < diffR11 && 0 < diffR12 && diffR12 < 0.100 && 0 < diffR21 && diffR21 < 0.100 ) { // matching photon1 (gen) and photon2 (reco)
	cout <<"EVENT -GEN  genSeedXYZ: "<<genSeed1.X() <<"  "<< genSeed1.Y() <<"  "<<genSeed1.Z() <<"   deltaR(gen,reco): "<<deltaR_gpt1<<" genPhoPt: "<<gpho1_Pt<<" pT diff  "<< diffpt12 <<"  genPhoPx: "<<gParticlePx[pho1_index]<< "  genPhoPy: "<<gParticlePy[pho1_index]<< "  genPhoPz: "<<gParticlePz[pho1_index]<< "  genPhoE: "<<gParticleE[pho1_index]<<" pho1 - gen vertex:  " << creation_time << "  neutralino dist:  "<<distance1<< "  neutralino TOF:  " << time1 << "  photon dist:  "<<distance1_p<< "  photon TOF:  " << time1_p << "  initialTOF:  " << TOF_initial << "  TOF_neu1:  " << TOF_neu1 << "  TOF_pho1:  " << TOF_pho1 <<  "TOF_total1:  "<< TOF_initial + TOF_neu1 + TOF_pho1 <<"  TOF_total1-L/c:  " << TOF_total2  <<"  TOF_total1-L/c (gen):  " << TOF_total2_genV <<"  neutralino decay x (gen):  "<<decay_x1<< "  neutralino decay y (gen):  "<<decay_y1<< "  neutralino decay z (gen):  "<<decay_z1<< "  Length gen to crystal:  " <<L_gen_crys2<<"  RECO  seedHitRawTime:   "<<t2raw_seed <<"  seedHitTimeTOF: "<<t2_seed<<"  seedHitTimeTOF genV: "<<t2_seed_genV<< "  pho2SeedX:  "<<pho2SeedX<<"  pho2SeedY:  "<<pho2SeedY<<"  pho2SeedZ:  "<<pho2SeedZ<< endl; 
	cout <<"EVENT -GEN  genSeedXYZ: "<<genSeed2.X() <<"  "<< genSeed2.Y() <<"  "<<genSeed2.Z() <<"   deltaR(gen,reco): "<<deltaR_gpt2<<" genPhoPt: "<<gpho2_Pt<<" pT diff  "<< diffpt21 <<"  genPhoPx: "<<gParticlePx[pho2_index]<< "  genPhoPy: "<<gParticlePy[pho2_index]<< "  genPhoPz: "<<gParticlePz[pho2_index]<< "  genPhoE: "<<gParticleE[pho2_index]<<" pho2 - gen vertex:  " << creation_time << "  neutralino dist:  "<<distance2<< "  neutralino TOF:  " << time2 << "  photon dist:  "<<distance2_p<< "  photon TOF:  " << time2_p << "  initialTOF:  " << TOF_initial << "  TOF_neu2:  " << TOF_neu2 << "  TOF_pho2:  " << TOF_pho2 <<  "TOF_total2:  "<< TOF_initial + TOF_neu2 + TOF_pho2 <<"  TOF_total2-L/c:  " << TOF_total1  <<"  TOF_total2-L/c (gen):  " << TOF_total1_genV <<"  neutralino decay x (gen):  "<<decay_x2<< "  neutralino decay y (gen):  "<<decay_y2<< "  neutralino decay z (gen):  "<<decay_z2<< "  Length gen to crystal:  " <<L_gen_crys1<<"  RECO  seedHitRawTime:   "<<t1raw_seed <<"  seedHitTimeTOF: "<<t1_seed<<"  seedHitTimeTOF genV: "<<t1_seed_genV<< "  pho1SeedX:  "<<pho1SeedX<<"  pho1SeedY:  "<<pho1SeedY<<"  pho1SeedZ:  "<<pho1SeedZ<< endl; 
  }

      //<<"time1: "<<time1<<" time2: "<<time2<<endl; 
      //cout<<"momentum1: "<<momentum1<<" momentum2: "<<momentum2<<endl; 
      //cout<<"distance1: "<<distance1<<" distance2: "<<distance2<<endl; 
      //cout<<"distanceX1: "<<distanceX1<<" distanceX2: "<<distanceX2<<endl; 
      //cout<<"distanceY1: "<<distanceY1<<" distanceY2: "<<distanceY2<<endl; 
      //cout<<"distanceZ1: "<<distanceZ1<<" distanceZ2: "<<distanceZ2<<endl; 
      //cout<<"px1: "<<primary_x1<<"  dx1: "<<decay_x1<<"py1: "<<primary_y1<<"  dy1: "<<decay_y1<<"pz1: "<<primary_z1<<"  dz1: "<<decay_z1<<endl;
      //cout<<"NeutraPt: "<<NeutraPt<<"  NeutraEta: "<<NeutraEta<<endl;
      //cout<<"distances: "<<distanceY1<<"  x direction:"<<distanceX1<<endl;
      //cout << "distance: "<<distance1<<"  momentum: "<<momentum<<"time of flight: " << time1 << endl;
}


  //Fill Event
  //if (mass > 60 && mass < 120) {
//  if ( t1_seed > 0 && t2_seed > 0) {
  {
  outputTree->Fill();
  }

  }
}//end of event loop
  
cout << "Writing output trees..." << endl;
outputTree->Write();
NEvents->Write();

outFile->Close();
}
