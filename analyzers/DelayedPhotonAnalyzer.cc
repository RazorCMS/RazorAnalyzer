#include "DelayedPhotonAnalyzer.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns


TVector3 DelayedPhotonAnalyzer::intersectPoint(float x0,float y0,float z0,float px,float py,float pz,float R)
{
  TVector3 sol;


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
	sol.SetXYZ(x1,y1,z1);
  }

  else if( (z2-z0)*pz > 0.0 )
  {
	sol.SetXYZ(x2,y2,z2);
  }

  else
  {
	sol.SetXYZ(0,0,0);
  }	
  return sol;
};

float DelayedPhotonAnalyzer::getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID) {
  float timeCalib = 1.0;

  // accessing variables for the seed hit information
  
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

  isData = false;//////remember to delete this when we are dealing with data.....

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
  float TOF_neu1, TOF_neu2;
  float TOF_neu1_RF, TOF_neu2_RF;
  float TOF_pho1, TOF_pho2;
  float t1calib_seed, t2calib_seed;
  float t1raw_seed, t2raw_seed;
  float phoNumber;
  float pho1Energy, pho1_Pt, pho1_Eta, pho1_Phi, pho1_angle_xtal;
  float pho2Energy, pho2_Pt, pho2_Eta, pho2_Phi, pho2_angle_xtal;
  bool pho1_isStandardPhoton, pho2_isStandardPhoton;
  float jetEnergy, jet_Pt, jet_Eta, jet_Phi;
  float met_Pt, met_Phi, sum_MET;
  float deltaPt_pho1, deltaPt_pho2;
  float deltaR_pho1, deltaR_pho2;
  float deltaEta_pho1, deltaEta_pho2;
  float deltaPhi_pho1, deltaPhi_pho2;
  float reco_eta1, reco_eta2;
  float gen_eta1, gen_eta2;
  float R1, R2;
  float ZD1, ZD2;

  float pho1SeedE;
  float pho1SeedEta;
  float pho1SeedPhi;
  float pho2SeedE;
  float pho2SeedEta;
  float pho2SeedPhi;

  int NPU;
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
  outputTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  outputTree->Branch("TOF_neu1", &TOF_neu1, "TOF_neu1/F");
  outputTree->Branch("TOF_neu1_RF", &TOF_neu1_RF, "TOF_neu1_RF/F");
  outputTree->Branch("TOF_neu2", &TOF_neu2, "TOF_neu2/F");
  outputTree->Branch("TOF_neu2_RF", &TOF_neu2_RF, "TOF_neu2_RF/F");
  outputTree->Branch("TOF_pho1", &TOF_pho1, "TOF_pho1/F");
  outputTree->Branch("TOF_pho2", &TOF_pho2, "TOF_pho2/F");

  outputTree->Branch("phoNumber", &phoNumber, "phoNumber/F");
  outputTree->Branch("pho1Energy", &pho1Energy, "pho1Energy/F");
  outputTree->Branch("pho1SeedE", &pho1SeedE, "pho1SeedE/F");
  outputTree->Branch("pho1SeedEta", &pho1SeedEta, "pho1SeedEta/F");
  outputTree->Branch("pho1SeedPhi", &pho1SeedPhi, "pho1SeedPhi/F");
  outputTree->Branch("pho1_isStandardPhoton", &pho1_isStandardPhoton, "pho1_isStandardPhoton/O");
  outputTree->Branch("pho1_Pt", &pho1_Pt, "pho1_Pt/F");
  outputTree->Branch("pho1_angle_xtal", &pho1_angle_xtal, "pho1_angle_xtal/F");
  outputTree->Branch("pho1_Eta", &pho1_Eta, "pho1_Eta/F");
  outputTree->Branch("pho1_Phi", &pho1_Phi, "pho1_Phi/F");  

  outputTree->Branch("pho2Energy", &pho2Energy, "pho2Energy/F");
  outputTree->Branch("pho2SeedE", &pho2SeedE, "pho2SeedE/F");
  outputTree->Branch("pho2SeedEta", &pho2SeedEta, "pho2SeedEta/F");
  outputTree->Branch("pho2SeedPhi", &pho2SeedPhi, "pho2SeedPhi/F");
  outputTree->Branch("pho2_isStandardPhoton", &pho2_isStandardPhoton, "pho2_isStandardPhoton/O");
  outputTree->Branch("pho2_Pt", &pho2_Pt, "pho2_Pt/F");
  outputTree->Branch("pho2_angle_xtal", &pho2_angle_xtal, "pho2_angle_xtal/F");
  outputTree->Branch("pho2_Eta", &pho2_Eta, "pho2_Eta/F");
  outputTree->Branch("pho2_Phi", &pho2_Phi, "pho2_Phi/F");  

  outputTree->Branch("jetEnergy", &jetEnergy, "jetEnergy/F");
  outputTree->Branch("jet_Pt", &jet_Pt, "jet_Pt/F");
  outputTree->Branch("jet_Eta", &jet_Eta, "jet_Eta/F");
  outputTree->Branch("jet_Phi", &jet_Phi, "jet_Phi/F");

  outputTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
  outputTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
  outputTree->Branch("sum_MET", &sum_MET, "sum_MET/F");

  outputTree->Branch("deltaPt_pho1", &deltaPt_pho1, "deltaPt_pho1/F");
  outputTree->Branch("deltaPt_pho2", &deltaPt_pho2, "deltaPt_pho2/F");
  outputTree->Branch("deltaR_pho1", &deltaR_pho1, "deltaR_pho1/F");
  outputTree->Branch("deltaR_pho2", &deltaR_pho2, "deltaR_pho2/F");
  outputTree->Branch("deltaEta_pho1", &deltaEta_pho1, "deltaEta_pho1/F");
  outputTree->Branch("deltaEta_pho2", &deltaEta_pho2, "deltaEta_pho2/F");
  outputTree->Branch("deltaPhi_pho1", &deltaPhi_pho1, "deltaPhi_pho1/F");
  outputTree->Branch("deltaPhi_pho2", &deltaPhi_pho2, "deltaPhi_pho2/F");
  outputTree->Branch("reco_eta1", &reco_eta1, "reco_eta1/F");
  outputTree->Branch("reco_eta2", &reco_eta2, "reco_eta2/F");
  outputTree->Branch("gen_eta1", &gen_eta1, "gen_eta1/F");
  outputTree->Branch("gen_eta2", &gen_eta2, "gen_eta2/F");

  outputTree->Branch("R1", &R1, "R1/F");
  outputTree->Branch("R2", &R2, "R2/F");
  outputTree->Branch("ZD1", &ZD1, "ZD1/F");
  outputTree->Branch("ZD2", &ZD2, "ZD2/F");

  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);


  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<10000;jentry++) {
    //begin event
    if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
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
    TOF_neu1 = -999;
    TOF_neu2 = -999;
    TOF_pho1 = -999;
    TOF_pho2 = -999;

    TOF_neu1_RF = -999;
    TOF_neu2_RF = -999;

    phoNumber = -999;
    pho1Energy = -999;
    pho1SeedE = -999;
    pho1SeedEta = -999;
    pho1SeedPhi = -999;
    pho1_Pt = -999;
    pho1_angle_xtal = -999;
    pho1_isStandardPhoton = true;
    pho1_Eta = -999;
    pho1_Phi = -999;
    
    pho2Energy = -999;
    pho2SeedE = -999;
    pho2SeedEta = -999;
    pho2SeedPhi = -999;
    pho2_Pt = -999;
    pho2_angle_xtal = -999;
    pho2_isStandardPhoton = true;
    pho2_Eta = -999;
    pho2_Phi = -999;

    jetEnergy = -999;
    jet_Pt = -999;
    jet_Eta = -999;
    jet_Phi = -999;

    met_Pt = -999;
    met_Phi = -999;
    sum_MET = -999;

    deltaPt_pho1 = -999;
    deltaPt_pho2 = -999;
    deltaR_pho1 = -999;
    deltaR_pho2 = -999;
    deltaEta_pho1 = -999;
    deltaEta_pho2 = -999;
    deltaPhi_pho1 = -999;
    deltaPhi_pho2 = -999;
    reco_eta1 = -999;
    reco_eta2 = -999;
    gen_eta1 = -999;
    gen_eta2 = -999;

    R1 = -999;
    R2 = -999;
    ZD1 = -999;
    ZD2 = -999;

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

  for(int ind_pho = 0; ind_pho < nPhotons; ind_pho++) 
  { //photon loop
      	// apply cuts
      	if(phoPt[ind_pho] < 25) continue; 
      	if(fabs(phoEta[ind_pho]) > 2.5) continue;
      	if(fabs(phoEta[ind_pho]) > 1.4442 && fabs(phoEta[ind_pho]) < 1.566) continue; //the eta range for photon, this takes care of the gap between barrel and endcap
      	//if(!(isEGammaPOGTightElectron(i))) continue;
      
      	nPho++;
      	TLorentzVector thisPhoton = makeTLorentzVector(phoPt[ind_pho], phoEta[ind_pho], phoPhi[ind_pho], phoE[ind_pho]);
      
      	//rough definition
      	uint seedhitIndex =  (*pho_SeedRechitIndex)[ind_pho];
    
      	//cout<<"reco Photon - "<<i<<" : seedX = "<<(*ecalRechit_X)[seedhitIndex]<<" : seedY = "<<(*ecalRechit_Y)[seedhitIndex]<<" : seedZ = "<<(*ecalRechit_Z)[seedhitIndex]<<"  pT = "<<phoPt[ind_pho]<<"  Energy = "<<phoE[ind_pho]<<endl;

      	//cout<<"seedhitIndex: "<<seedhitIndex<<endl;
      	//cout<<"ecalRechit_ID size: "<<ecalRechit_ID->size()<<endl;
      	//cout<<"ecalRechit_ID: "<<(*ecalRechit_ID)[seedhitIndex]<<endl;

      	bool isEBOrEE = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 ); //barrel vs. endcap
      	double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      	//apply intercalibration2      
      	double IC_time_SeptRereco = isData ? getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[seedhitIndex]) : 0;
      	double IC_time_LagacyRereco = isData ? getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]) : 0;
      	double calibratedSeedHitTime = rawSeedHitTime + IC_time_LagacyRereco - IC_time_SeptRereco;

      	//apply TOF correction
      	double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
      
      	//generator xyz information

      	double TOFCorrectedSeedHitTime_genV = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-genVertexX,2)+pow((*ecalRechit_Y)[seedhitIndex]-genVertexY,2)+pow((*ecalRechit_Z)[seedhitIndex]-genVertexZ,2)))/SPEED_OF_LIGHT;

      	double tmpSumWeightedTime = 0;
      	double tmpSumWeight = 0;

      	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k) 
	{
		//cout << metPt << endl;
        	uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
		MET_event = sumMET;
        	MET_Pt_event = metPt;
        	MET_Phi_event = metPhi;
      
        	double rawT = (*ecalRechit_T)[rechitIndex];
        	//apply intercalibration
        	double corrT = rawT + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

        	// double pedNoise = getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex]);
        	double pedNoise = 1;
        	double ADCToGeV = isData ? getADCToGeV(runNum, isEBOrEE) : 1;
        	double sigmaE = pedNoise * ADCToGeV;
  
        	float C_EB = isData ? 1 : 0;
        	float N_EB = 1;
  
        	double sigmaT = N_EB / ((*ecalRechit_E)[rechitIndex] / sigmaE) + sqrt(2) * C_EB;
        	tmpSumWeightedTime += corrT * ( 1.0 / (sigmaT*sigmaT) );
        	tmpSumWeight += ( 1.0 / (sigmaT*sigmaT) );
        	// cout << "\n";
      	}
      	double weightedTime = tmpSumWeightedTime / tmpSumWeight;
         
      	if (thisPhoton.Pt() > pho1.Pt()) 
	{ // find two highest momentum photons
        	pho1 = thisPhoton;
        	// finds the rechit position of the first photon
  		pho1SeedX = (*ecalRechit_X)[seedhitIndex];
  		pho1SeedY = (*ecalRechit_Y)[seedhitIndex];
  		pho1SeedZ = (*ecalRechit_Z)[seedhitIndex];
  	
		pho1SeedE = (*ecalRechit_E)[seedhitIndex];
		pho1SeedEta = (*ecalRechit_Eta)[seedhitIndex];
		pho1SeedPhi = (*ecalRechit_Phi)[seedhitIndex];
        
        	pho1_isStandardPhoton = pho_isStandardPhoton[ind_pho];
        	pho1_time = weightedTime;
        	pho1_seedtime = TOFCorrectedSeedHitTime;
        	pho1_seedtime_genV = TOFCorrectedSeedHitTime_genV;
        	pho1_seedtimeCalib = calibratedSeedHitTime;
        	pho1_seedtimeraw = rawSeedHitTime;
    	} 
    	else if (thisPhoton.Pt() > pho2.Pt()) 
	{
      		pho2 = thisPhoton;
      		pho2SeedX = (*ecalRechit_X)[seedhitIndex];
      		pho2SeedY = (*ecalRechit_Y)[seedhitIndex];
      		pho2SeedZ = (*ecalRechit_Z)[seedhitIndex];

      		pho2SeedE = (*ecalRechit_E)[seedhitIndex]; 
      		pho2SeedEta = (*ecalRechit_Eta)[seedhitIndex];
      		pho2SeedPhi = (*ecalRechit_Phi)[seedhitIndex];

      		pho2_isStandardPhoton = pho_isStandardPhoton[ind_pho];

      		pho2_time = weightedTime;
      		pho2_seedtime = TOFCorrectedSeedHitTime; 
      		pho2_seedtimeCalib = calibratedSeedHitTime;
      		pho2_seedtimeraw = rawSeedHitTime;
	}    
 } //end photon loop

 if (nPho >= 2) 
   { 
    	//cout << "THIS IS THE 2 PHOTON LOOP" << endl;
    
	bool isMatched = false;

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

    	bool foundN1 = false;
    	bool foundN2 = false; 
    	int pho1_index = 0;
    	int pho2_index = 0;
    	int neu1_index = 0;
    	int neu2_index = 0;

    	// finding the neutralino and photon index
    	for(unsigned int i = 0; i < nGenParticle; i++)
	{ //this is gen particle loop within event and photon loop
		if ( !foundN1 && gParticleId[i] == 22 && gParticleMotherId[i] == 1000022 )
		{ //finds a photon from a neutralino
			pho1_index = i;
			neu1_index = gParticleMotherIndex[i];
			foundN1 = true;
		}
		else if ( foundN1 && !foundN2 && gParticleId[i] == 22 && gParticleMotherId[i] == 1000022 ) 
		{
			pho2_index = i;
			neu2_index = gParticleMotherIndex[i];
			foundN2 = true;
		}
    	}

  	if(foundN1==1 && foundN2==1)
	{

  		R1 = pow(gParticleDecayVertexX[neu1_index]*gParticleDecayVertexX[neu1_index] + gParticleDecayVertexY[neu1_index]*gParticleDecayVertexY[neu1_index] + gParticleDecayVertexZ[neu1_index]*gParticleDecayVertexZ[neu1_index],0.5);
  		R2 = pow(gParticleDecayVertexX[neu2_index]*gParticleDecayVertexX[neu2_index] + gParticleDecayVertexY[neu2_index]*gParticleDecayVertexY[neu2_index] + gParticleDecayVertexZ[neu2_index]*gParticleDecayVertexZ[neu2_index],0.5);
 
  		float decay_x1 = gParticleDecayVertexX[neu1_index];
  		float decay_y1 = gParticleDecayVertexY[neu1_index];
  		float decay_z1 = gParticleDecayVertexZ[neu1_index];
  		float decay_x2 = gParticleDecayVertexX[neu2_index];
  		float decay_y2 = gParticleDecayVertexY[neu2_index];
  		float decay_z2 = gParticleDecayVertexZ[neu2_index];

  		ZD1 = decay_z1;
  		ZD2 = decay_z2;

  		// need to match up the photon index and the reco photon - this is done based on momentum
  		// pho1_Pt is reco level, gpho1_Pt is gen level information
  		float gpho1_Pt = gParticlePt[pho1_index];
  		float gpho2_Pt = gParticlePt[pho2_index];
  		float deltaPt11 = fabs(pho1_Pt - gpho1_Pt);
  		float deltaPt21 = fabs(pho1_Pt - gpho2_Pt);
  		float deltaPt12 = fabs(pho2_Pt - gpho1_Pt);
  		float deltaPt22 = fabs(pho2_Pt - gpho2_Pt);

  		TVector3 genSeed1;
  		TVector3 genSeed2;

  		float norm1 = pow((pow(gParticlePx[pho1_index],2)+pow(gParticlePy[pho1_index],2)+pow(gParticlePz[pho1_index],2)),0.5);
  		float px1 = (gParticlePx[pho1_index]) / norm1;
  		float py1 = (gParticlePy[pho1_index]) / norm1;
  		float pz1 = (gParticlePz[pho1_index]) / norm1;
  		genSeed1 = intersectPoint(decay_x1, decay_y1, decay_z1, px1, py1, pz1, 129.7); // using intersection function written above, radius as 129.7 cm
  		float norm2 = pow((pow(gParticlePx[pho2_index],2)+pow(gParticlePy[pho2_index],2)+pow(gParticlePz[pho2_index],2)),0.5);
  		float px2 = (gParticlePx[pho2_index]) / norm2;
  		float py2 = (gParticlePy[pho2_index]) / norm2;
  		float pz2 = (gParticlePz[pho2_index]) / norm2;
  		genSeed2 = intersectPoint(decay_x2, decay_y2, decay_z2, px2, py2, pz2, 129.0); // using intersection function written above, radius as 129 cm

  		TVector3 recoSeed1(pho1SeedX,pho1SeedY,pho1SeedZ);
  		TVector3 recoSeed2(pho2SeedX,pho2SeedY,pho2SeedZ);

  		float deltaR11 = genSeed1.DeltaR(recoSeed1); 
  		float deltaR12 = genSeed1.DeltaR(recoSeed2);
  		float deltaR21 = genSeed2.DeltaR(recoSeed1);
  		float deltaR22 = genSeed2.DeltaR(recoSeed2);

  		float deltaEta11 = recoSeed1.Eta()-genSeed1.Eta();
  		float deltaPhi11 = recoSeed1.Phi()-genSeed1.Phi();
  		float deltaEta12 = recoSeed2.Eta()-genSeed1.Eta();
  		float deltaPhi12 = recoSeed2.Phi()-genSeed1.Phi();

  		float deltaEta21 = recoSeed1.Eta()-genSeed2.Eta();
  		float deltaPhi21 = recoSeed1.Phi()-genSeed2.Phi();
  		float deltaEta22 = recoSeed2.Eta()-genSeed2.Eta();
  		float deltaPhi22 = recoSeed2.Phi()-genSeed2.Phi();

  		reco_eta1 = recoSeed1.Eta();
  		reco_eta2 = recoSeed2.Eta();
  		
		bool is1To1 = false;
  //cout << "deltaR11  "<<deltaR11<<"  deltaR12  "<<deltaR12<<"  deltaR21  "<<deltaR21<<"  deltaR22  "<<deltaR22<<endl;

  		if ( deltaR11 < deltaR12 && deltaPt11 < deltaPt12) 
		{
			is1To1 = true;
			isMatched = true;
		}
  		else if(deltaR12 < deltaR11 && deltaPt12 < deltaPt11)
		{
			is1To1 = false;
			isMatched = true;
		}
		
		if(isMatched)
		{
			gen_eta1 = is1To1 ? genSeed1.Eta() : genSeed2.Eta();
			gen_eta2 = is1To1 ? genSeed2.Eta() : genSeed1.Eta();
			deltaR_pho1 = is1To1 ? deltaR11 : deltaR21;
			deltaEta_pho1 = is1To1 ? deltaEta11 : deltaEta21;
			deltaPhi_pho1 = is1To1 ? deltaPhi11 : deltaPhi21;
			deltaPt_pho1 = is1To1 ? deltaPt11 : deltaPt21;
			deltaR_pho2 = is1To1 ? deltaR22 : deltaR12;
			deltaEta_pho2 = is1To1 ? deltaEta22 : deltaEta12;
			deltaPhi_pho2 = is1To1 ? deltaPhi22 : deltaPhi12;
			deltaPt_pho2 = is1To1 ? deltaPt22 : deltaPt12;

			float mass = 1000.0;
			float p_neu1 = is1To1 ? gParticlePt[neu1_index]/cosh(gParticleEta[neu1_index]) : gParticlePt[neu2_index]/cosh(gParticleEta[neu2_index]) ;
			float p_neu2 = is1To1 ? gParticlePt[neu2_index]/cosh(gParticleEta[neu2_index]) : gParticlePt[neu1_index]/cosh(gParticleEta[neu1_index]) ;
		
			TVector3 point_genPV(genVertexX,genVertexY,genVertexZ);	
			TVector3 point_decayV1(is1To1 ? decay_x1: decay_x2,is1To1 ? decay_y1: decay_y2, is1To1 ? decay_z1: decay_z2);
			TVector3 point_decayV2(is1To1 ? decay_x2: decay_x1,is1To1 ? decay_y2: decay_y1, is1To1 ? decay_z2: decay_z1);

			TOF_neu1 = (point_decayV1-point_genPV).Mag() / (SPEED_OF_LIGHT*p_neu1) * pow((pow(mass,2) + pow(p_neu1,2)),0.5);
			TOF_neu2 = (point_decayV2-point_genPV).Mag() / (SPEED_OF_LIGHT*p_neu2) * pow((pow(mass,2) + pow(p_neu2,2)),0.5);
			TOF_neu1_RF = TOF_neu1*mass*pow((pow(mass,2) + pow(p_neu1,2)),-0.5);
			TOF_neu2_RF = TOF_neu2*mass*pow((pow(mass,2) + pow(p_neu2,2)),-0.5);

			TOF_pho1 = (recoSeed1 - point_decayV1).Mag() / SPEED_OF_LIGHT ;
			TOF_pho2 = (recoSeed2 - point_decayV2).Mag() / SPEED_OF_LIGHT ;

			TOF_total1 = genVertexT + TOF_neu1 + TOF_pho1 - recoSeed1.Mag() / SPEED_OF_LIGHT;
			TOF_total1_genV = genVertexT + TOF_neu1 + TOF_pho1 - (recoSeed1 - point_genPV).Mag() / SPEED_OF_LIGHT;
			TOF_total2 = genVertexT + TOF_neu2 + TOF_pho2 - recoSeed2.Mag() / SPEED_OF_LIGHT;
			TOF_total2_genV = genVertexT + TOF_neu2 + TOF_pho2 - (recoSeed2 - point_genPV).Mag() / SPEED_OF_LIGHT;
			
			pho1_angle_xtal = recoSeed1.Angle(recoSeed1 - point_decayV1); 
			pho2_angle_xtal = recoSeed2.Angle(recoSeed2 - point_decayV2); 
			
			outputTree->Fill();		
		}//if isMatched
	}//if gen found
  }//if nPho>=2
}//event loop

cout << "Writing output trees..." << endl;
outputTree->Write();
NEvents->Write();
outFile->Close();

}//analyzer function


