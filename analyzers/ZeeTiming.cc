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

  //*****************************************************************************
  //Load Pedestals
  //*****************************************************************************
  vector <uint> start_time;//start run of all IOV 
  vector <uint> end_time;//end run of all IOV
  start_time_tmp=0; 
  end_time_tmp=0;
  rms_G12_all=0;
  detID_all=0 ;

  TFile f_pedestal("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1.root","READ");
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
    std::cout << "HggRazor: Output filename not specified!" << endl << "Using default output name HggRazor.root" << std::endl;
    outFileName = "HggRazor.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );




  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("ZeeTiming", "Info on selected razor inclusive events");

  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float mass;
  float t1, t2;
  float t1_seed, t2_seed;
  float t1calib_seed, t2calib_seed;
  float t1raw_seed, t2raw_seed;
  float ele1E, ele1Pt, ele1Eta, ele1Phi;
  float ele2E, ele2Pt, ele2Eta, ele2Phi;
  int  ele1SeedIEta, ele1SeedIPhi, ele1SeedIX, ele1SeedIY;
  bool ele1IsEB;
  int  ele2SeedIEta, ele2SeedIPhi, ele2SeedIX, ele2SeedIY; 
  bool ele2IsEB;
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
  // outputTree->Branch("ele1E", &ele1E, "ele1E/F");
  outputTree->Branch("ele1Pt", &ele1Pt, "ele1Pt/F");
  // outputTree->Branch("ele1Eta", &ele1Eta, "ele1Eta/F");
  // outputTree->Branch("ele1Phi", &ele1Phi, "ele1Phi/F");
  outputTree->Branch("ele1IsEB", &ele1IsEB, "ele1IsEB/O");
  outputTree->Branch("ele1SeedIEta", &ele1SeedIEta, "ele1SeedIEta/I");
  outputTree->Branch("ele1SeedIPhi", &ele1SeedIPhi, "ele1SeedIPhi/I");
  outputTree->Branch("ele1SeedIX", &ele1SeedIX, "ele1SeedIX/I");
  outputTree->Branch("ele1SeedIY", &ele1SeedIY, "ele1SeedIY/I");
  // outputTree->Branch("ele2E", &ele2E, "ele2E/F");
  outputTree->Branch("ele2Pt", &ele2Pt, "ele2Pt/F");
  // outputTree->Branch("ele2Eta", &ele2Eta, "ele2Eta/F");
  // outputTree->Branch("ele2Phi", &ele2Phi, "ele2Phi/F");
  outputTree->Branch("ele2IsEB", &ele2IsEB, "ele2IsEB/O");
  outputTree->Branch("ele2SeedIEta", &ele2SeedIEta, "ele2SeedIEta/I");
  outputTree->Branch("ele2SeedIPhi", &ele2SeedIPhi, "ele2SeedIPhi/I");
  outputTree->Branch("ele2SeedIX", &ele2SeedIX, "ele2SeedIX/I");
  outputTree->Branch("ele2SeedIY", &ele2SeedIY, "ele2SeedIY/I");

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
    double ele1_seedtime = 0;
    double ele2_seedtime = 0;
    double ele1_seedtimeCalib = 0;
    double ele2_seedtimeCalib = 0;
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
      bool isEBOrEE = bool( eleEta_SC[i] < 1.5 );

      uint seedhitIndex =  (*ele_SeedRechitIndex)[i];
      double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      //apply intercalibration
      double calibratedSeedHitTime = rawSeedHitTime + getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]);

      //apply TOF correction
      double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;


      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << " : " 
      // 	   << (*ecalRechit_ID)[seedhitIndex] << " " << rawSeedHitTime << " -> " << calibratedSeedHitTime << " -> " << TOFCorrectedSeedHitTime << " "
      // 	   << "\n";
 

   
      double tmpSumWeightedTime = 0;
      double tmpSumWeight = 0;

      for (uint k=0; k<(*ele_EcalRechitIndex)[i].size(); ++k) {
      	
	uint rechitIndex = (*ele_EcalRechitIndex)[i][k];
		  
	double rawT = (*ecalRechit_T)[rechitIndex];
	//apply intercalibration
	//apply TOF correction
	double corrT = rawT + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

	double pedNoise = getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex]);
	double ADCToGeV = getADCToGeV(runNum, isEBOrEE);
	double sigmaE = pedNoise * ADCToGeV;
	double sigmaT = N_EB / ((*ecalRechit_E)[rechitIndex] / sigmaE) + sqrt(2) * C_EB;
	tmpSumWeightedTime += corrT * ( 1.0 / (sigmaT*sigmaT) );
	tmpSumWeight += ( 1.0 / (sigmaT*sigmaT) );
	// cout << "\n";
      }
      double weightedTime = tmpSumWeightedTime / tmpSumWeight;

            
      if (thisElectron.Pt() > ele1.Pt()) {
	ele1 = thisElectron;
	ele1_time = weightedTime;
	ele1_seedtime = TOFCorrectedSeedHitTime;
	ele1_seedtimeCalib = calibratedSeedHitTime;
	ele1_seedtimeraw = rawSeedHitTime;
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
	ele2_time = weightedTime;
	ele2_seedtime = TOFCorrectedSeedHitTime; 
	ele2_seedtimeCalib = calibratedSeedHitTime;
 	ele2_seedtimeraw = rawSeedHitTime;
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
      t2 = ele2_time;
      t1_seed = ele1_seedtime;
      t2_seed = ele2_seedtime;
      t1calib_seed = ele1_seedtimeCalib;
      t2calib_seed = ele2_seedtimeCalib;
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
