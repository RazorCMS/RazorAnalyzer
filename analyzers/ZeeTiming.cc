#include "ZeeTiming.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns


float ZeeTiming::getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID) {
  float timeCalib = 1.0;
  
  int N_entries = tree->GetEntries();
  int i_entry=0;
  for(uint i=0;i<start_run.size();i++)
    {
      if(run>= start_run[i] && run<= end_run[i])
	{
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

float ZeeTiming::getPedestalNoise(TTree *tree, vector <uint> & start_time, vector <uint> & end_time, uint time, uint detID) {
  float pedestalNoise = 1.0;
  
  int N_entries = tree->GetEntries();
  int i_entry=0;
  for(uint i=0;i<start_time.size();i++)
    {
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


float ZeeTiming::getADCToGeV( uint run, int isEBOrEE) {
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


void ZeeTiming::Analyze(bool isData, int option, string outFileName, string label)
{

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
  detID_all=0 ;

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
  TFile f_pedestal("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1_G12rmsonly.root","READ");
  TTree *tree_pedestal = (TTree*)f_pedestal.Get("pedestal");
  
  tree_pedestal->SetBranchAddress("start_time_second", &start_time_tmp);
  tree_pedestal->SetBranchAddress("end_time_second", &end_time_tmp);
  tree_pedestal->SetBranchAddress("rms_G12", &rms_G12_all);
  tree_pedestal->SetBranchAddress("detID", &detID_all);
  
  int N_entries_pedestal = tree_pedestal->GetEntries();
  
  cout << "Total Pedestal IOVs: " << N_entries_pedestal << "\n";
  for(int i=0;i<N_entries_pedestal;i++) {
    cout << "Loading Pedestal IOV " << i << "\n";
    tree_pedestal->GetEntry(i);
    start_time.push_back(start_time_tmp);
    end_time.push_back(end_time_tmp);
  }

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
    std::cout << "ZeeTiming: Output filename not specified!" << endl << "Using default output name ZeeTiming.root" << std::endl;
    outFileName = "ZeeTiming.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );




  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("ZeeTiming", "Info on selected razor inclusive events");

  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float mass;
  float t1_TOF2, t2_TOF2;
  float t1, t2;
  float t1_seed, t2_seed;
  float seed1_TC_lagacy, seed2_TC_lagacy;
  float seed1_TC_sept, seed2_TC_sept;
  float seed1_pedestal, seed2_pedestal;
  float seed1_transpCorr, seed2_transpCorr;
  float t1calib_seed, t2calib_seed;
  float t1calib_seed_sept, t2calib_seed_sept;
  float t1raw_seed, t2raw_seed;
  float ele1E, ele1Pt, ele1Eta, ele1Phi;
  float ele2E, ele2Pt, ele2Eta, ele2Phi;
  int  ele1SeedIEta, ele1SeedIPhi, ele1SeedIX, ele1SeedIY;
  bool ele1IsEB;
  int  ele2SeedIEta, ele2SeedIPhi, ele2SeedIX, ele2SeedIY; 
  bool ele2IsEB;
  int NPU;
  vector<float> *ecalElectronRechit_E;
  vector<float> *ecalElectronRechit_rawT;
  vector<float> *ecalElectronRechit_calibT_sept;
  vector<float> *ecalElectronRechit_calibT_lagacy;
  vector<float> *ecalElectronRechit_calibT_lagacy_TOF;
  vector<float> *ecalElectronRechit_Eta;
  vector<float> *ecalElectronRechit_Phi;
  vector<float> *ecalElectronRechit_transpCorr;
	
  ecalElectronRechit_E = new std::vector<float>; ecalElectronRechit_E->clear();
  ecalElectronRechit_rawT = new std::vector<float>; ecalElectronRechit_rawT->clear();
  ecalElectronRechit_calibT_sept = new std::vector<float>; ecalElectronRechit_calibT_sept->clear();
  ecalElectronRechit_calibT_lagacy = new std::vector<float>; ecalElectronRechit_calibT_lagacy->clear();
  ecalElectronRechit_calibT_lagacy_TOF = new std::vector<float>; ecalElectronRechit_calibT_lagacy_TOF->clear();
  ecalElectronRechit_Eta = new std::vector<float>; ecalElectronRechit_Eta->clear();
  ecalElectronRechit_Phi = new std::vector<float>; ecalElectronRechit_Phi->clear();
  ecalElectronRechit_transpCorr = new std::vector<float>; ecalElectronRechit_transpCorr->clear();
  //int nPV;
  unsigned int run, lumi, event;


  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("eventTime", &eventTime, "eventTime/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("nPV", &nPV, "nPV/i");
  outputTree->Branch("mass", &mass, "mass/F");
  outputTree->Branch("t1", &t1, "t1/F");
  outputTree->Branch("t2", &t2, "t2/F");
  outputTree->Branch("t1_TOF2", &t1_TOF2, "t1_TOF2/F");
  outputTree->Branch("t2_TOF2", &t2_TOF2, "t2_TOF2/F");
  outputTree->Branch("t1_seed", &t1_seed, "t1_seed/F");
  outputTree->Branch("t2_seed", &t2_seed, "t2_seed/F");
  outputTree->Branch("seed1_TC_lagacy", &seed1_TC_lagacy, "seed1_TC_lagacy/F");
  outputTree->Branch("seed2_TC_lagacy", &seed2_TC_lagacy, "seed2_TC_lagacy/F");
  outputTree->Branch("seed2_TC_sept", &seed2_TC_sept, "seed2_TC_sept/F");
  outputTree->Branch("seed1_TC_sept", &seed1_TC_sept, "seed1_TC_sept/F");
  outputTree->Branch("seed1_pedestal", &seed1_pedestal, "seed1_pedestal/F");
  outputTree->Branch("seed2_pedestal", &seed2_pedestal, "seed2_pedestal/F");
  outputTree->Branch("seed1_transpCorr", &seed1_transpCorr, "seed1_transpCorr/F");
  outputTree->Branch("seed2_transpCorr", &seed2_transpCorr, "seed2_transpCorr/F");
  outputTree->Branch("t1calib_seed", &t1calib_seed, "t1calib_seed/F");
  outputTree->Branch("t2calib_seed", &t2calib_seed, "t2calib_seed/F");
  outputTree->Branch("t1calib_seed_sept", &t1calib_seed_sept, "t1calib_seed_sept/F");
  outputTree->Branch("t2calib_seed_sept", &t2calib_seed_sept, "t2calib_seed_sept/F");
  outputTree->Branch("t1raw_seed", &t1raw_seed, "t1raw_seed/F");
  outputTree->Branch("t2raw_seed", &t2raw_seed, "t2raw_seed/F");
  outputTree->Branch("ele1E", &ele1E, "ele1E/F");
  outputTree->Branch("ele1Pt", &ele1Pt, "ele1Pt/F");
  outputTree->Branch("ele1Eta", &ele1Eta, "ele1Eta/F");
  outputTree->Branch("ele1Phi", &ele1Phi, "ele1Phi/F");
  outputTree->Branch("ele1IsEB", &ele1IsEB, "ele1IsEB/O");
  outputTree->Branch("ele1SeedIEta", &ele1SeedIEta, "ele1SeedIEta/I");
  outputTree->Branch("ele1SeedIPhi", &ele1SeedIPhi, "ele1SeedIPhi/I");
  outputTree->Branch("ele1SeedIX", &ele1SeedIX, "ele1SeedIX/I");
  outputTree->Branch("ele1SeedIY", &ele1SeedIY, "ele1SeedIY/I");
  outputTree->Branch("ele2E", &ele2E, "ele2E/F");
  outputTree->Branch("ele2Pt", &ele2Pt, "ele2Pt/F");
  outputTree->Branch("ele2Eta", &ele2Eta, "ele2Eta/F");
  outputTree->Branch("ele2Phi", &ele2Phi, "ele2Phi/F");
  outputTree->Branch("ele2IsEB", &ele2IsEB, "ele2IsEB/O");
  outputTree->Branch("ele2SeedIEta", &ele2SeedIEta, "ele2SeedIEta/I");
  outputTree->Branch("ele2SeedIPhi", &ele2SeedIPhi, "ele2SeedIPhi/I");
  outputTree->Branch("ele2SeedIX", &ele2SeedIX, "ele2SeedIX/I");
  outputTree->Branch("ele2SeedIY", &ele2SeedIY, "ele2SeedIY/I");
  outputTree->Branch("ecalElectronRechit_E", "std::vector<float>",&ecalElectronRechit_E);
  outputTree->Branch("ecalElectronRechit_rawT", "std::vector<float>",&ecalElectronRechit_rawT);
  outputTree->Branch("ecalElectronRechit_calibT_sept", "std::vector<float>",&ecalElectronRechit_calibT_sept);
  outputTree->Branch("ecalElectronRechit_calibT_lagacy", "std::vector<float>",&ecalElectronRechit_calibT_lagacy);
  outputTree->Branch("ecalElectronRechit_calibT_lagacy_TOF", "std::vector<float>",&ecalElectronRechit_calibT_lagacy_TOF);
  outputTree->Branch("ecalElectronRechit_Eta", "std::vector<float>",&ecalElectronRechit_Eta);
  outputTree->Branch("ecalElectronRechit_Phi", "std::vector<float>",&ecalElectronRechit_Phi);
  outputTree->Branch("ecalElectronRechit_transpCorr", "std::vector<float>",&ecalElectronRechit_transpCorr);

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
    t1_TOF2 = -999;
    t2_TOF2 = -999;
    t1_seed = -999;
    t2_seed = -999;
    seed1_TC_lagacy = -999;
    seed2_TC_lagacy = -999;
    seed1_TC_sept = -999;
    seed2_TC_sept = -999;
    seed1_pedestal = -999;
    seed2_pedestal = -999;
    seed1_transpCorr = -999;
    seed2_transpCorr = -999;
    t1calib_seed = -999;
    t2calib_seed = -999;
    t1calib_seed_sept = -999;
    t2calib_seed_sept = -999;
    t1raw_seed = -999;
    t2raw_seed = -999;
    ele1E = -999; 
    ele1Pt = -999;
    ele1Eta = -999;
    ele1Phi = -999;
    ele2E = -999;
    ele2Pt = -999;
    ele2Eta = -999;
    ele2Phi = -999;
    ele1SeedIEta = -999;
    ele1SeedIPhi = -999;
    ele1SeedIX = -999;
    ele1SeedIY = -999;
    ele1IsEB = 0;
    ele2SeedIEta = -999;
    ele2SeedIPhi = -999;
    ele2SeedIX = -999;
    ele2SeedIY = -999; 
    ele2IsEB = 0;
    NPU = 0;   
    ecalElectronRechit_E->clear();
    ecalElectronRechit_rawT->clear();
    ecalElectronRechit_calibT_sept->clear();
    ecalElectronRechit_calibT_lagacy->clear();
    ecalElectronRechit_calibT_lagacy_TOF->clear();
    ecalElectronRechit_Eta->clear();
    ecalElectronRechit_Phi->clear();
    ecalElectronRechit_transpCorr->clear();


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


    int nEle = 0;
    TLorentzVector ele1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector ele2 = makeTLorentzVector(0,0,0,0);
    double ele1_time = 0;
    double ele2_time= 0;
    double ele1_time_TOF2 = 0;
    double ele2_time_TOF2 = 0;
    double ele1_seedtime = 0;
    double ele2_seedtime = 0;
    double ele1_seedtimeCalib = 0;
    double ele2_seedtimeCalib = 0;
    double ele1_seedtimeCalib_sept = 0;
    double ele2_seedtimeCalib_sept = 0;
    double ele1_seedtimeraw = 0;
    double ele2_seedtimeraw = 0;
    for(int i = 0; i < nElectrons; i++){
      if(elePt[i] < 35) continue;
      if(fabs(eleEta[i]) > 2.5) continue;
      if(fabs(eleEta[i]) > 1.4442 && fabs(eleEta[i]) < 1.566) continue;
      if(!(isEGammaPOGTightElectron(i))) continue;
      
      nEle++;
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
          
      //rough definition
      uint seedhitIndex =  (*ele_SeedRechitIndex)[i];
      bool isEBOrEE = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 );
      double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      //apply intercalibration
      double IC_time_SeptRereco = getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double IC_time_LagacyRereco = getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]);
      double calibratedSeedHitTime = rawSeedHitTime + IC_time_LagacyRereco - IC_time_SeptRereco;
      double calibratedSeedHitTime_sept = rawSeedHitTime;
      double rawSeedHitTime_no_sept = rawSeedHitTime - IC_time_SeptRereco;

      //apply TOF correction
      double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;


      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << " : \n" 
      //	<< "  runNum: " << runNum << "  detID: " << (*ecalRechit_ID)[seedhitIndex] << "  IC_time_SeptRereco: " << IC_time_SeptRereco << "  IC_time_LagacyRereco: " << IC_time_LagacyRereco << " \n"
      //   	<< " " << rawSeedHitTime << " -> " << calibratedSeedHitTime << " -> " << TOFCorrectedSeedHitTime << " "
      //	<< "\n";
 


      double seedPedNoise = getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex]);
   
      double tmpSumWeightedTime = 0;
      double tmpSumWeightedTime_TOF2 = 0;
      double tmpSumWeight = 0;
      double tmpSumWeight_TOF2 = 0;

      for (uint k=0; k<(*ele_EcalRechitIndex)[i].size(); ++k) {	
	uint rechitIndex = (*ele_EcalRechitIndex)[i][k];
    
      	  
	double rawT_this = (*ecalRechit_T)[rechitIndex];
	//apply intercalibration
	double IC_time_SeptRereco_this = getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[rechitIndex]);
        double IC_time_LagacyRereco_this = getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[rechitIndex]);
        double calibratedSeedHitTime_this = rawT_this + IC_time_LagacyRereco_this - IC_time_SeptRereco_this;

	//apply TOF correction
	double corrT = calibratedSeedHitTime_this + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

  	ecalElectronRechit_E->push_back((*ecalRechit_E)[rechitIndex]); 
        ecalElectronRechit_Eta->push_back((*ecalRechit_Eta)[rechitIndex]); 
        ecalElectronRechit_Phi->push_back((*ecalRechit_Phi)[rechitIndex]); 
        ecalElectronRechit_transpCorr->push_back((*ecalRechit_transpCorr)[rechitIndex]); 
        ecalElectronRechit_rawT->push_back(rawT_this - IC_time_SeptRereco_this); 
        ecalElectronRechit_calibT_sept->push_back(rawT_this); 
        ecalElectronRechit_calibT_lagacy->push_back(calibratedSeedHitTime_this); 
        ecalElectronRechit_calibT_lagacy_TOF->push_back(corrT); 

	double pedNoise = getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[rechitIndex]);
	//double pedNoise = 1;
	double ADCToGeV = getADCToGeV(runNum, isEBOrEE);
	double sigmaE = pedNoise * ADCToGeV;
	double sigmaT = N_EB / ((*ecalRechit_E)[rechitIndex] / sigmaE) + sqrt(2) * C_EB;
	tmpSumWeightedTime += corrT * ( 1.0 / (sigmaT*sigmaT) );
	tmpSumWeightedTime_TOF2 += calibratedSeedHitTime_this * ( 1.0 / (sigmaT*sigmaT) );
	tmpSumWeight += ( 1.0 / (sigmaT*sigmaT) );
	tmpSumWeight_TOF2 += ( 1.0 / (sigmaT*sigmaT) );
	// cout << "\n";
      }
      double weightedTime = tmpSumWeightedTime / tmpSumWeight;
      double weightedTime_TOF2 = tmpSumWeightedTime_TOF2 / tmpSumWeight_TOF2 + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

            
      if (thisElectron.Pt() > ele1.Pt()) {
	ele1 = thisElectron;
	seed1_TC_lagacy = IC_time_LagacyRereco;
	seed1_TC_sept = IC_time_SeptRereco;
	seed1_pedestal = seedPedNoise; 
	seed1_transpCorr = (*ecalRechit_transpCorr)[seedhitIndex]; 
	ele1_time = weightedTime;
	ele1_time_TOF2 = weightedTime_TOF2;
	ele1_seedtime = TOFCorrectedSeedHitTime;
	ele1_seedtimeCalib = calibratedSeedHitTime;
	ele1_seedtimeCalib_sept = calibratedSeedHitTime_sept;
	ele1_seedtimeraw = rawSeedHitTime_no_sept;
	ele1IsEB = bool( eleEta_SC[i] < 1.5 );
	if (ele1IsEB) {
	  ele1SeedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele1SeedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele1SeedIX = -999;
	  ele1SeedIY = -999;
	} else {
	  ele1SeedIEta = -999;
	  ele1SeedIPhi = -999;
	  ele1SeedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele1SeedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	}

      } else if (thisElectron.Pt() > ele2.Pt()) {
	ele2 = thisElectron;
	seed2_TC_lagacy = IC_time_LagacyRereco;
	seed2_TC_sept = IC_time_SeptRereco;
	seed2_pedestal = seedPedNoise; 
	seed2_transpCorr = (*ecalRechit_transpCorr)[seedhitIndex]; 
	ele2_time = weightedTime;
	ele2_time_TOF2 = weightedTime_TOF2;
	ele2_seedtime = TOFCorrectedSeedHitTime; 
	ele2_seedtimeCalib = calibratedSeedHitTime;
	ele2_seedtimeCalib_sept = calibratedSeedHitTime_sept;
 	ele2_seedtimeraw = rawSeedHitTime_no_sept;
 	ele2IsEB = bool( eleEta_SC[i] < 1.5 );
	if (ele2IsEB) {
	  ele2SeedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele2SeedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele2SeedIX = -999;
	  ele2SeedIY = -999;
	} else {
	  ele2SeedIEta = -999;
	  ele2SeedIPhi = -999;
	  ele2SeedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele2SeedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	}

      }	
    }
    
    if (nEle >= 2) {
      ele1E = ele1.E();
      ele1Pt = ele1.Pt();
      ele1Eta = ele1.Eta();
      ele1Phi = ele1.Phi();
      ele2E = ele2.E();
      ele2Pt = ele2.Pt();
      ele2Eta = ele2.Eta();
      ele2Phi = ele2.Phi();

      mass = (ele1+ele2).M();
      t1 = ele1_time;
      t1_TOF2 = ele1_time_TOF2;
      t2 = ele2_time;
      t2_TOF2 = ele2_time_TOF2;
      t1_seed = ele1_seedtime;
      t2_seed = ele2_seedtime;
      t1calib_seed = ele1_seedtimeCalib;
      t2calib_seed = ele2_seedtimeCalib;
      t1calib_seed_sept = ele1_seedtimeCalib_sept;
      t2calib_seed_sept = ele2_seedtimeCalib_sept;
      t1raw_seed = ele1_seedtimeraw;
      t2raw_seed = ele2_seedtimeraw;
       //cout << "ele2: " << ele2.Pt() << " " << ele2_seedtime << "\n";
    }

    //Fill Event
    if (mass > 60 && mass < 120) {
      outputTree->Fill();
    }

  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  outputTree->Write();
  NEvents->Write();

  outFile->Close();
}
