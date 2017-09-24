//LOCLAL INCLUDES
#include "DelayedPhotonAnalyzer.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "RazorHelper.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns
const float EB_R = 129.7;
const float EE_Z = 317.0;
const double JET_CUT = 30.;

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


float DelayedPhotonAnalyzer::getADCToGeV( uint run, int isFromEB) {
  double ADCToGeV = 0;
  //EB
  if (isFromEB == 1) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.039680;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.039798;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.039436;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.039298;
  }   
  //EE
  else if (isFromEB == 0) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.067230;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.067370;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.066764;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.065957;
  }
  return ADCToGeV;
}


void DelayedPhotonAnalyzer::Analyze(bool isData, int option, string outFileName, string label) {

  //isData = false;//////remember to delete this when we are dealing with data.....

  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  //bool doPhotonScaleCorrection = true;

  string analysisTag = "Razor2016_MoriondRereco";
  if ( label != "") analysisTag = label;

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
  if(isData)
  { 
	  for(int i=0;i<N_entries_timeCalib;i++) {
	    tree_timeCalib->GetEntry(i);
	    start_run.push_back(start_run_tmp);
	    end_run.push_back(end_run_tmp);
	  }
  }


  TFile f_timeCalib_rereco("/eos/cms/store/group/phys_susy/razor/EcalTiming/EcalTimeCalibConstants_v08_offline/tree_EcalTimeCalibConstants_v08_offline.root","READ");
  TTree *tree_timeCalib_rereco = (TTree*)f_timeCalib_rereco.Get("timeCalib");
  
  tree_timeCalib_rereco->SetBranchAddress("start_run", &start_run_tmp);
  tree_timeCalib_rereco->SetBranchAddress("end_run", &end_run_tmp);
  tree_timeCalib_rereco->SetBranchAddress("IC_time", &IC_time_all);
  tree_timeCalib_rereco->SetBranchAddress("detID", &detID_all);
  
  int N_entries_timeCalib_rereco = tree_timeCalib_rereco->GetEntries();
  if(isData)
  { 
	  for(int i=0;i<N_entries_timeCalib_rereco;i++) {
	    tree_timeCalib_rereco->GetEntry(i);
	    start_run_rereco.push_back(start_run_tmp);
	    end_run_rereco.push_back(end_run_tmp);
	  }
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
 
  if(isData)
  { 
	  cout << "Total Pedestal IOVs: " << N_entries_pedestal << "\n";
	  for(int i=0;i<N_entries_pedestal;i++) {
	    cout << "Loading Pedestal IOV " << i << "\n";
	    tree_pedestal->GetEntry(i);
	    start_time.push_back(start_time_tmp);
	    end_time.push_back(end_time_tmp);
	  }
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
    std::cout << "DelayedPhotonAnalyzer: Output filename not specified!" << endl << "Using default output name DelayedPhoton.root" << std::endl;
    outFileName = "DelayedPhoton.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );


  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("DelayedPhoton", "Delayed photon events");


  int NPU;
  unsigned int run, lumi, event;
  float genVertexTime = 0.0;//genVertexT;
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float TOF_total1;
  float TOF_total2;
  float TOF_total1_genV;
  float TOF_total2_genV;
  float TOF_neu1, TOF_neu2;
  float TOF_neu1_RF, TOF_neu2_RF;
  float TOF_pho1, TOF_pho2;
  
  int n_Photons;
  float pho1E, pho1Pt, pho1Eta, pho1Phi, pho1SeedE, pho1SeedPt, pho1SeedEta, pho1SeedPhi, pho1SC_E, pho1SC_Pt, pho1SC_Eta, pho1SC_Phi, pho1angle_xtal, pho1SigmaIetaIeta, pho1R9, pho1HoverE, pho1sumChargedHadronPt, pho1sumNeutralHadronEt, pho1sumPhotonEt, pho1sigmaEOverE, pho1SeedTimeRaw, pho1SeedTimeCalib, pho1SeedTimeCalibTOF, pho1SeedTimeGenV, pho1ClusterTime, pho1Sminor, pho1Smajor, pho1Setaeta, pho1Sphiphi, pho1Setaphi, pho1GenE, pho1GenPt, pho1GenEta, pho1GenPhi;
  float pho2E, pho2Pt, pho2Eta, pho2Phi, pho2SeedE, pho2SeedPt, pho2SeedEta, pho2SeedPhi, pho2SC_E, pho2SC_Pt, pho2SC_Eta, pho2SC_Phi, pho2angle_xtal, pho2SigmaIetaIeta, pho2R9, pho2HoverE, pho2sumChargedHadronPt, pho2sumNeutralHadronEt, pho2sumPhotonEt, pho2sigmaEOverE, pho2SeedTimeRaw, pho2SeedTimeCalib, pho2SeedTimeCalibTOF, pho2SeedTimeGenV, pho2ClusterTime, pho2Sminor, pho2Smajor, pho2Setaeta, pho2Sphiphi, pho2Setaphi, pho2GenE, pho2GenPt, pho2GenEta, pho2GenPhi;
  bool pho1passEleVeto, pho1passIsoLoose, pho1passIsoMedium, pho1passIsoTight, pho1isStandardPhoton;
  bool pho2passEleVeto, pho2passIsoLoose, pho2passIsoMedium, pho2passIsoTight, pho2isStandardPhoton;

  int n_Jets;
  float MET, t1MET;//, MET_JESUp, MET_JESDown, t1MET_JESUp, t1MET_JESDown;
  float HT;

  float deltaPt_pho1, deltaPt_pho2;
  float deltaR_pho1, deltaR_pho2;
  float deltaEta_pho1, deltaEta_pho2;
  float deltaPhi_pho1, deltaPhi_pho2;
  float reco_eta1, reco_eta2;
  float gen_eta1, gen_eta2;
  float R1, R2;
  float ZD1, ZD2;

  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("nPV", &nPV, "nPV/i");
  outputTree->Branch("pvX", &pvX, "pvX/F");
  outputTree->Branch("pvY", &pvY, "pvY/F");
  outputTree->Branch("pvZ", &pvZ, "pvZ/F");

  outputTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  outputTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  outputTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  outputTree->Branch("genVertexTime", &genVertexTime, "genVertexTime/F");
  
  outputTree->Branch("TOF_total1", &TOF_total1, "TOF_total1/F");
  outputTree->Branch("TOF_total2", &TOF_total2, "TOF_total2/F");
  outputTree->Branch("TOF_total1_genV", &TOF_total1_genV, "TOF_total1_genV/F");
  outputTree->Branch("TOF_total2_genV", &TOF_total2_genV, "TOF_total2_genV/F");
  outputTree->Branch("TOF_neu1", &TOF_neu1, "TOF_neu1/F");
  outputTree->Branch("TOF_neu1_RF", &TOF_neu1_RF, "TOF_neu1_RF/F");
  outputTree->Branch("TOF_neu2", &TOF_neu2, "TOF_neu2/F");
  outputTree->Branch("TOF_neu2_RF", &TOF_neu2_RF, "TOF_neu2_RF/F");
  outputTree->Branch("TOF_pho1", &TOF_pho1, "TOF_pho1/F");
  outputTree->Branch("TOF_pho2", &TOF_pho2, "TOF_pho2/F");

  outputTree->Branch("n_Photons", &n_Photons, "n_Photons/I"); // 1 or 2
  outputTree->Branch("pho1E", &pho1E, "pho1E/F");
  outputTree->Branch("pho1Pt", &pho1Pt, "pho1Pt/F");
  outputTree->Branch("pho1Eta", &pho1Eta, "pho1Eta/F");
  outputTree->Branch("pho1Phi", &pho1Phi, "pho1Phi/F");  
  outputTree->Branch("pho1SeedE", &pho1SeedE, "pho1SeedE/F");
  outputTree->Branch("pho1SeedPt", &pho1SeedPt, "pho1SeedPt/F");
  outputTree->Branch("pho1SeedEta", &pho1SeedEta, "pho1SeedEta/F");
  outputTree->Branch("pho1SeedPhi", &pho1SeedPhi, "pho1SeedPhi/F");
  outputTree->Branch("pho1SC_E", &pho1SC_E, "pho1SC_E/F");
  outputTree->Branch("pho1SC_Pt", &pho1SC_Pt, "pho1SC_Pt/F");
  outputTree->Branch("pho1SC_Eta", &pho1SC_Eta, "pho1SC_Eta/F");
  outputTree->Branch("pho1SC_Phi", &pho1SC_Phi, "pho1SC_Phi/F");
  outputTree->Branch("pho1isStandardPhoton", &pho1isStandardPhoton, "pho1isStandardPhoton/O");
  outputTree->Branch("pho1angle_xtal", &pho1angle_xtal, "pho1angle_xtal/F");
  outputTree->Branch("pho1SigmaIetaIeta", &pho1SigmaIetaIeta, "pho1SigmaIetaIeta/F");
  outputTree->Branch("pho1R9", &pho1R9, "pho1R9/F");
  outputTree->Branch("pho1HoverE", &pho1HoverE, "pho1HoverE/F");
  outputTree->Branch("pho1sumChargedHadronPt", &pho1sumChargedHadronPt, "pho1sumChargedHadronPt/F");
  outputTree->Branch("pho1sumNeutralHadronEt", &pho1sumNeutralHadronEt, "pho1sumNeutralHadronEt/F");
  outputTree->Branch("pho1sumPhotonEt", &pho1sumPhotonEt, "pho1sumPhotonEt/F");
  outputTree->Branch("pho1sigmaEOverE", &pho1sigmaEOverE, "pho1sigmaEOverE/F");
  outputTree->Branch("pho1passEleVeto", &pho1passEleVeto, "pho1passEleVeto/O");
  outputTree->Branch("pho1passIsoLoose", &pho1passIsoLoose, "pho1passIsoLoose/O");
  outputTree->Branch("pho1passIsoMedium", &pho1passIsoMedium, "pho1passIsoMedium/O");
  outputTree->Branch("pho1passIsoTight", &pho1passIsoTight, "pho1passIsoTight/O");
  outputTree->Branch("pho1SeedTimeRaw", &pho1SeedTimeRaw, "pho1SeedTimeRaw/F");
  outputTree->Branch("pho1SeedTimeCalib", &pho1SeedTimeCalib, "pho1SeedTimeCalib/F");
  outputTree->Branch("pho1SeedTimeCalibTOF", &pho1SeedTimeCalibTOF, "pho1SeedTimeCalibTOF/F");
  outputTree->Branch("pho1SeedTimeGenV", &pho1SeedTimeGenV, "pho1SeedTimeGenV/F");
  outputTree->Branch("pho1ClusterTime", &pho1ClusterTime, "pho1ClusterTime/F");
  outputTree->Branch("pho1Sminor", &pho1Sminor, "pho1Sminor/F");
  outputTree->Branch("pho1Smajor", &pho1Smajor, "pho1Smajor/F");
  outputTree->Branch("pho1Setaeta", &pho1Setaeta, "pho1Setaeta/F");
  outputTree->Branch("pho1Sphiphi", &pho1Sphiphi, "pho1Sphiphi/F");
  outputTree->Branch("pho1Setaphi", &pho1Setaphi, "pho1Setaphi/F");
  outputTree->Branch("pho1GenE", &pho1GenE, "pho1GenE/F");
  outputTree->Branch("pho1GenPt", &pho1GenPt, "pho1GenPt/F");
  outputTree->Branch("pho1GenEta", &pho1GenEta, "pho1GenEta/F");
  outputTree->Branch("pho1GenPhi", &pho1GenPhi, "pho1GenPhi/F");

  outputTree->Branch("pho2E", &pho2E, "pho2E/F");
  outputTree->Branch("pho2Pt", &pho2Pt, "pho2Pt/F");
  outputTree->Branch("pho2Eta", &pho2Eta, "pho2Eta/F");
  outputTree->Branch("pho2Phi", &pho2Phi, "pho2Phi/F");  
  outputTree->Branch("pho2SeedE", &pho2SeedE, "pho2SeedE/F");
  outputTree->Branch("pho2SeedPt", &pho2SeedPt, "pho2SeedPt/F");
  outputTree->Branch("pho2SeedEta", &pho2SeedEta, "pho2SeedEta/F");
  outputTree->Branch("pho2SeedPhi", &pho2SeedPhi, "pho2SeedPhi/F");
  outputTree->Branch("pho2SC_E", &pho2SC_E, "pho2SC_E/F");
  outputTree->Branch("pho2SC_Pt", &pho2SC_Pt, "pho2SC_Pt/F");
  outputTree->Branch("pho2SC_Eta", &pho2SC_Eta, "pho2SC_Eta/F");
  outputTree->Branch("pho2SC_Phi", &pho2SC_Phi, "pho2SC_Phi/F");
  outputTree->Branch("pho2isStandardPhoton", &pho2isStandardPhoton, "pho2isStandardPhoton/O");
  outputTree->Branch("pho2angle_xtal", &pho2angle_xtal, "pho2angle_xtal/F");
  outputTree->Branch("pho2SigmaIetaIeta", &pho2SigmaIetaIeta, "pho2SigmaIetaIeta/F");
  outputTree->Branch("pho2R9", &pho2R9, "pho2R9/F");
  outputTree->Branch("pho2HoverE", &pho2HoverE, "pho2HoverE/F");
  outputTree->Branch("pho2sumChargedHadronPt", &pho2sumChargedHadronPt, "pho2sumChargedHadronPt/F");
  outputTree->Branch("pho2sumNeutralHadronEt", &pho2sumNeutralHadronEt, "pho2sumNeutralHadronEt/F");
  outputTree->Branch("pho2sumPhotonEt", &pho2sumPhotonEt, "pho2sumPhotonEt/F");
  outputTree->Branch("pho2sigmaEOverE", &pho2sigmaEOverE, "pho2sigmaEOverE/F");
  outputTree->Branch("pho2passEleVeto", &pho2passEleVeto, "pho2passEleVeto/O");
  outputTree->Branch("pho2passIsoLoose", &pho2passIsoLoose, "pho2passIsoLoose/O");
  outputTree->Branch("pho2passIsoMedium", &pho2passIsoMedium, "pho2passIsoMedium/O");
  outputTree->Branch("pho2passIsoTight", &pho2passIsoTight, "pho2passIsoTight/O");
  outputTree->Branch("pho2SeedTimeRaw", &pho2SeedTimeRaw, "pho2SeedTimeRaw/F");
  outputTree->Branch("pho2SeedTimeCalib", &pho2SeedTimeCalib, "pho2SeedTimeCalib/F");
  outputTree->Branch("pho2SeedTimeCalibTOF", &pho2SeedTimeCalibTOF, "pho2SeedTimeCalibTOF/F");
  outputTree->Branch("pho2SeedTimeGenV", &pho2SeedTimeGenV, "pho2SeedTimeGenV/F");
  outputTree->Branch("pho2ClusterTime", &pho2ClusterTime, "pho2ClusterTime/F");
  outputTree->Branch("pho2Sminor", &pho2Sminor, "pho2Sminor/F");
  outputTree->Branch("pho2Smajor", &pho2Smajor, "pho2Smajor/F");
  outputTree->Branch("pho2Setaeta", &pho2Setaeta, "pho2Setaeta/F");
  outputTree->Branch("pho2Sphiphi", &pho2Sphiphi, "pho2Sphiphi/F");
  outputTree->Branch("pho2Setaphi", &pho2Setaphi, "pho2Setaphi/F");
  outputTree->Branch("pho2GenE", &pho2GenE, "pho2GenE/F");
  outputTree->Branch("pho2GenPt", &pho2GenPt, "pho2GenPt/F");
  outputTree->Branch("pho2GenEta", &pho2GenEta, "pho2GenEta/F");
  outputTree->Branch("pho2GenPhi", &pho2GenPhi, "pho2GenPhi/F");

  outputTree->Branch("n_Jets", &n_Jets, "n_Jets/I");

  outputTree->Branch("MET", &MET, "MET/F");
//  outputTree->Branch("MET_JESUp", &MET_JESUp, "MET_JESUp/F");
//  outputTree->Branch("MET_JESDown", &MET_JESDown, "MET_JESDown/F");
  outputTree->Branch("t1MET", &t1MET, "t1MET/F");
//  outputTree->Branch("t1MET_JESUp", &t1MET_JESUp, "t1MET_JESUp/F");
//  outputTree->Branch("t1MET_JESDown", &t1MET_JESDown, "t1MET_JESDown/F");
  outputTree->Branch("HT", &HT, "HT/F");

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
    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    genVertexTime = 0.0;//genVertexT
    NPU = 0;   
    weight = 0;
    pileupWeight = 0;
    pileupWeightUp = 0;
    pileupWeightDown = 0;
    
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

    n_Photons = 0;
    pho1E = -999, pho1Pt = -999, pho1Eta = -999, pho1Phi = -999, pho1SeedE = -999, pho1SeedPt = -999, pho1SeedEta = -999, pho1SeedPhi = -999, pho1SC_E = -999, pho1SC_Pt = -999, pho1SC_Eta = -999, pho1SC_Phi = -999, pho1angle_xtal = -999, pho1SigmaIetaIeta = -999, pho1R9 = -999, pho1HoverE = -999, pho1sumChargedHadronPt = -999, pho1sumNeutralHadronEt = -999, pho1sumPhotonEt = -999, pho1sigmaEOverE = -999, pho1SeedTimeRaw = -999, pho1SeedTimeCalib = -999, pho1SeedTimeCalibTOF = -999, pho1SeedTimeGenV = -999, pho1ClusterTime = -999, pho1Sminor = -999, pho1Smajor = -999, pho1Setaeta = -999, pho1Sphiphi = -999, pho1Setaphi = -999, pho1GenE = -999, pho1GenPt = -999, pho1GenEta = -999, pho1GenPhi = -999;
    pho2E = -999, pho2Pt = -999, pho2Eta = -999, pho2Phi = -999, pho2SeedE = -999, pho2SeedPt = -999, pho2SeedEta = -999, pho2SeedPhi = -999, pho2SC_E = -999, pho2SC_Pt = -999, pho2SC_Eta = -999, pho2SC_Phi = -999, pho2angle_xtal = -999, pho2SigmaIetaIeta = -999, pho2R9 = -999, pho2HoverE = -999, pho2sumChargedHadronPt = -999, pho2sumNeutralHadronEt = -999, pho2sumPhotonEt = -999, pho2sigmaEOverE = -999, pho2SeedTimeRaw = -999, pho2SeedTimeCalib = -999, pho2SeedTimeCalibTOF = -999, pho2SeedTimeGenV = -999, pho2ClusterTime = -999, pho2Sminor = -999, pho2Smajor = -999, pho2Setaeta = -999, pho2Sphiphi = -999, pho2Setaphi = -999, pho2GenE = -999, pho2GenPt = -999, pho2GenEta = -999, pho2GenPhi = -999;
    pho1passEleVeto = false, pho1passIsoLoose = false, pho1passIsoMedium = false, pho1passIsoTight = false, pho1isStandardPhoton = false;
    pho2passEleVeto = false, pho2passIsoLoose = false, pho2passIsoMedium = false, pho2passIsoTight = false, pho2isStandardPhoton = false;

    n_Jets = 0;
    MET = -999, t1MET = -999;//MET_JESUp = -999, MET_JESDown = -999, t1MET_JESUp = -999, t1MET_JESDown = -999;
    HT = -999;

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


    //fill normalization histogram
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;

    //get NPU
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
        NPU = nPUmean[i];
      }
    }
   
    int nPho = 0;
    TLorentzVector pho1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector pho2 = makeTLorentzVector(0,0,0,0);

    // XYZ rechit where photon is detected
    float pho1SeedX = 0;
    float pho1SeedY = 0;
    float pho1SeedZ = 0;
    float pho2SeedX = 0;
    float pho2SeedY = 0;
    float pho2SeedZ = 0;

    TVector3 vtx( pvX, pvY, pvZ );


  for(int ind_pho = 0; ind_pho < nPhotons; ind_pho++) 
  { //photon loop
      	// apply cuts
      	if(phoPt[ind_pho] < 25) continue; // basic Pt cut
      	if(fabs(phoEta[ind_pho]) > 2.5) continue; // tracker region
      	if(fabs(phoEta[ind_pho]) > 1.4442 && fabs(phoEta[ind_pho]) < 1.566) continue; //the eta range for photon, this takes care of the gap between barrel and endcap
      	//if(!(isEGammaPOGTightElectron(i))) continue;
      	
	nPho++;
   	//photon cluster
      	TLorentzVector thisPhoton = makeTLorentzVector(phoPt[ind_pho], phoEta[ind_pho], phoPhi[ind_pho], phoE[ind_pho]);
      
	//photon super cluster
	TVector3 phoPos;
       	if ( fabs( pho_superClusterEta[ind_pho] ) < 1.479 )
       	{
        	phoPos.SetXYZ( EB_R*cos( pho_superClusterPhi[ind_pho]), EB_R*sin( pho_superClusterPhi[ind_pho] ), EB_R*sinh( pho_superClusterEta[ind_pho] ) );
      	}
        else
        {
              	double R = fabs( EE_Z/sinh( pho_superClusterEta[ind_pho] ) );
              	if ( pho_superClusterEta[ind_pho] > .0 )
                {
                  	phoPos.SetXYZ( R*cos( pho_superClusterPhi[ind_pho] ), R*sin( pho_superClusterPhi[ind_pho] ), EE_Z);
                }
              	else
                {
                  	phoPos.SetXYZ( R*cos( pho_superClusterPhi[ind_pho] ), R*sin( pho_superClusterPhi[ind_pho] ), -EE_Z);
                }
     	}
   	TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_RegressionE[ind_pho] );
 
      	//rough definition
      	uint seedhitIndex =  (*pho_SeedRechitIndex)[ind_pho];
    
      	//cout<<"reco Photon - "<<i<<" : seedX = "<<(*ecalRechit_X)[seedhitIndex]<<" : seedY = "<<(*ecalRechit_Y)[seedhitIndex]<<" : seedZ = "<<(*ecalRechit_Z)[seedhitIndex]<<"  pT = "<<phoPt[ind_pho]<<"  Energy = "<<phoE[ind_pho]<<endl;
      	//cout<<"seedhitIndex: "<<seedhitIndex<<endl;
      	//cout<<"ecalRechit_ID size: "<<ecalRechit_ID->size()<<endl;
      	//cout<<"ecalRechit_ID: "<<(*ecalRechit_ID)[seedhitIndex]<<endl;

      	bool isFromEB = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 ); //barrel vs. endcap
      	double rawSeedHitTime =  (*ecalRechit_T)[seedhitIndex];

      	//apply intercalibration2      
      	double IC_time_SeptRereco = isData ? getTimeCalibConstant(tree_timeCalib_rereco, start_run_rereco,end_run_rereco,runNum, (*ecalRechit_ID)[seedhitIndex]) : 0;
      	double IC_time_LagacyRereco = isData ? getTimeCalibConstant(tree_timeCalib, start_run,end_run,runNum, (*ecalRechit_ID)[seedhitIndex]) : 0;
      	double calibratedSeedHitTime = rawSeedHitTime + IC_time_LagacyRereco - IC_time_SeptRereco;

      	//apply TOF correction
      	double TOFCorrectedSeedHitTime = calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
      
      	//generator xyz information

      	double TOFCorrectedSeedHitTime_genV = isData ? TOFCorrectedSeedHitTime : ( calibratedSeedHitTime + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-genVertexX,2)+pow((*ecalRechit_Y)[seedhitIndex]-genVertexY,2)+pow((*ecalRechit_Z)[seedhitIndex]-genVertexZ,2)))/SPEED_OF_LIGHT );

      	double tmpSumWeightedTime = 0;
      	double tmpSumWeight = 0;

	double etaAverage = 0.0;
	double phiAverage = 0.0;
	double mTotalWeight = 0.0;
	double tmpSumE = 0.0;
	double phoSetaeta = 0.0;//second moments of eta eta
	double phoSphiphi = 0.0;
	double phoSetaphi = 0.0;
	double phoSminor = 0.0;
	double phoSmajor = 0.0;
	
      	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k) 
	{
		//cout << metPt << endl;
        	uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
      
        	double rawT = (*ecalRechit_T)[rechitIndex];
        	//apply intercalibration
        	double corrT = rawT + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

        	double pedNoise = isData ? (getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex])) : 1.0;
        	//double pedNoise = 1;
        	double ADCToGeV = isData ? getADCToGeV(runNum, isFromEB) : 1;
        	double sigmaE = pedNoise * ADCToGeV;
  
        	float C_EB = isData ? 1 : 0;
        	float N_EB = 1;
  
        	double sigmaT = N_EB / ((*ecalRechit_E)[rechitIndex] / sigmaE) + sqrt(2) * C_EB;
        	tmpSumWeightedTime += corrT * ( 1.0 / (sigmaT*sigmaT) );
        	tmpSumWeight += ( 1.0 / (sigmaT*sigmaT) );
        	// cout << "\n";
        	tmpSumE += (*ecalRechit_E)[rechitIndex];	
		
      	}
	
	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k)
        {
		uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
		double thisWeight = max(4.2+log(((*ecalRechit_E)[rechitIndex])/tmpSumE),0.0);
		mTotalWeight += thisWeight;
		//float thisIEtaIX = (*ecalRechit_Eta)[rechitIndex];
		//float thisIPhiIY = (*ecalRechit_Phi)[rechitIndex];
		float thisIPhiIY =  1.0 * iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);
		float thisIEtaIX =  1.0 * iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);

		etaAverage += thisWeight * (thisIEtaIX) ;	
		phiAverage += thisWeight * (thisIPhiIY) ;	
	}
	
	etaAverage = etaAverage / mTotalWeight;
	phiAverage = phiAverage / mTotalWeight;

	//cout<<"DEBUG rechit: ";	
	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k)
        {
                uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
                double thisWeight = max(4.2+log(((*ecalRechit_E)[rechitIndex])/tmpSumE),0.0);
		//float thisIEtaIX = (*ecalRechit_Eta)[rechitIndex];
		//float thisIPhiIY = (*ecalRechit_Phi)[rechitIndex];
		float thisIPhiIY =  1.0 * iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);
		float thisIEtaIX =  1.0 * iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);


		//computing the moments of eta and phi
		//cout<<" [ i "<<k<<" w "<<thisWeight<<" dd "<<(thisIEtaIX - etaAverage) * (thisIEtaIX - etaAverage)<<"  s  "<<thisWeight * (thisIEtaIX - etaAverage) * (thisIEtaIX - etaAverage)<<" ] ";
		phoSetaeta += thisWeight * (thisIEtaIX - etaAverage) * (thisIEtaIX - etaAverage);
		phoSphiphi += thisWeight * (thisIPhiIY - phiAverage) * (thisIPhiIY - phiAverage);
		phoSetaphi += thisWeight * (thisIEtaIX - etaAverage) * (thisIPhiIY - phiAverage);
        }
	//cout<<"  ===>   setaeta "<<phoSetaeta<<endl;
	
      	double weightedTime = tmpSumWeightedTime / tmpSumWeight;

	if(tmpSumE>0.0)
	{
		phoSetaeta = phoSetaeta / mTotalWeight;	
		phoSphiphi = phoSphiphi / mTotalWeight;	
		phoSetaphi = phoSetaphi / mTotalWeight;	
		phoSminor = 0.5 * (phoSetaeta + phoSphiphi - pow(pow(phoSetaeta - phoSphiphi,2.0) + 4.0*pow(phoSetaphi,2.0),0.5));
		phoSmajor = 0.5 * (phoSetaeta + phoSphiphi + pow(pow(phoSetaeta - phoSphiphi,2.0) + 4.0*pow(phoSetaphi,2.0),0.5));
	}
         
      	if (thisPhoton.Pt() > pho1.Pt()) 
	{ // find two highest momentum photons
        	pho1 = thisPhoton;
	
		pho1E = thisPhoton.E();
		pho1Pt = thisPhoton.Pt();
		pho1Eta = thisPhoton.Eta();
		pho1Phi = thisPhoton.Phi();
		pho1SeedE = (*ecalRechit_E)[seedhitIndex];
		pho1SeedEta = (*ecalRechit_Eta)[seedhitIndex];
		pho1SeedPhi = (*ecalRechit_Phi)[seedhitIndex];
		pho1SeedPt = pho1SeedE/cosh(pho1SeedEta);
		pho1SC_E = phoSC.E();
		pho1SC_Pt = phoSC.Pt();
		pho1SC_Eta = phoSC.Eta();
		pho1SC_Phi = phoSC.Phi();
		pho1SigmaIetaIeta = phoFull5x5SigmaIetaIeta[ind_pho];
		pho1R9 = phoR9[ind_pho];
		pho1HoverE = pho_HoverE[ind_pho];
		pho1sumChargedHadronPt = pho_pfIsoChargedHadronIso[ind_pho];
		pho1sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[ind_pho];
		pho1sumPhotonEt = pho_pfIsoPhotonIso[ind_pho];
		pho1sigmaEOverE = pho_RegressionEUncertainty[ind_pho]/pho_RegressionE[ind_pho];
        	pho1SeedTimeRaw = rawSeedHitTime;
        	pho1SeedTimeCalib = calibratedSeedHitTime;
        	pho1SeedTimeCalibTOF = TOFCorrectedSeedHitTime;
        	pho1SeedTimeGenV = TOFCorrectedSeedHitTime_genV;
        	pho1ClusterTime = weightedTime;
		pho1Sminor = phoSminor;
		pho1Smajor = phoSmajor;
		pho1Setaeta = phoSetaeta;
		pho1Sphiphi = phoSphiphi;
		pho1Setaphi = phoSetaphi;
		pho1passEleVeto = pho_passEleVeto[ind_pho];
		pho1passIsoLoose = photonPassLooseIso(ind_pho);
		pho1passIsoMedium = photonPassMediumIso(ind_pho);
		pho1passIsoTight = photonPassTightIso(ind_pho);
        	pho1isStandardPhoton = pho_isStandardPhoton[ind_pho];

  		pho1SeedX = (*ecalRechit_X)[seedhitIndex];
  		pho1SeedY = (*ecalRechit_Y)[seedhitIndex];
  		pho1SeedZ = (*ecalRechit_Z)[seedhitIndex];
    	} 
    	else if (thisPhoton.Pt() > pho2.Pt()) 
	{
      		pho2 = thisPhoton;
      	
		pho2E = thisPhoton.E();
		pho2Pt = thisPhoton.Pt();
		pho2Eta = thisPhoton.Eta();
		pho2Phi = thisPhoton.Phi();
		pho2SeedE = (*ecalRechit_E)[seedhitIndex];
		pho2SeedEta = (*ecalRechit_Eta)[seedhitIndex];
		pho2SeedPhi = (*ecalRechit_Phi)[seedhitIndex];
		pho2SeedPt = pho2SeedE/cosh(pho2SeedEta);
		pho2SC_E = phoSC.E();
		pho2SC_Pt = phoSC.Pt();
		pho2SC_Eta = phoSC.Eta();
		pho2SC_Phi = phoSC.Phi();
		pho2SigmaIetaIeta = phoFull5x5SigmaIetaIeta[ind_pho];
		pho2R9 = phoR9[ind_pho];
		pho2HoverE = pho_HoverE[ind_pho];
		pho2sumChargedHadronPt = pho_pfIsoChargedHadronIso[ind_pho];
		pho2sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[ind_pho];
		pho2sumPhotonEt = pho_pfIsoPhotonIso[ind_pho];
		pho2sigmaEOverE = pho_RegressionEUncertainty[ind_pho]/pho_RegressionE[ind_pho];
        	pho2SeedTimeRaw = rawSeedHitTime;
        	pho2SeedTimeCalib = calibratedSeedHitTime;
        	pho2SeedTimeCalibTOF = TOFCorrectedSeedHitTime;
        	pho2SeedTimeGenV = TOFCorrectedSeedHitTime_genV;
        	pho2ClusterTime = weightedTime;
		pho2Sminor = phoSminor;
		pho2Smajor = phoSmajor;
		pho2Setaeta = phoSetaeta;
		pho2Sphiphi = phoSphiphi;
		pho2Setaphi = phoSetaphi;
		pho2passEleVeto = pho_passEleVeto[ind_pho];
		pho2passIsoLoose = photonPassLooseIso(ind_pho);
		pho2passIsoMedium = photonPassMediumIso(ind_pho);
		pho2passIsoTight = photonPassTightIso(ind_pho);
        	pho2isStandardPhoton = pho_isStandardPhoton[ind_pho];

  		pho2SeedX = (*ecalRechit_X)[seedhitIndex];
  		pho2SeedY = (*ecalRechit_Y)[seedhitIndex];
  		pho2SeedZ = (*ecalRechit_Z)[seedhitIndex];
	}    
 } //end photon loop

 if(nPho == 0) continue; 

HT = 0.0;
HT = pho1Pt;
if(nPho>=2) HT += pho2Pt;

//jet loop
/*
RazorHelper *helper = 0;
helper = new RazorHelper(analysisTag, isData, false); 
std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
*/

for(int i = 0; i < nJets; i++)
{

	double JEC = 1.0;//JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
                         //                        fixedGridRhoAll, jetJetArea[i], runNum,
                         //                        JetCorrectorIOV, JetCorrector );
	
      	TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

       	if( thisJet.Pt() < JET_CUT ) continue;//According to the April 1st 2015 AN
        if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
        if ( !jetPassIDLoose[i] ) continue;

	double deltaRJetPhoton = 0.0;
	if(nPho==1) deltaRJetPhoton = thisJet.DeltaR( pho1 );
	if(nPho>=2) deltaRJetPhoton = min( thisJet.DeltaR( pho1 ), thisJet.DeltaR( pho2 ) );
          
	if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN
	n_Jets++;
	HT += thisJet.Pt();
}

 MET = metPt;
 t1MET = metType1Pt;
	
 
 //fill the output tree
 if (nPho >= 1) // require at least one photon
 { 
    	//cout << "THIS IS THE 2 PHOTON LOOP" << endl;
 	if(nPho == 1) n_Photons = 1;
 	if(nPho > 1) n_Photons = 2;
	bool isMatched = false;

	if(!isData && nPho>=2) //for two neutralino -> photon + gravitino signal MC studies only
	{
		bool foundN1 = false;
		bool foundN2 = false; 
		int pho1index = 0;
		int pho2index = 0;
		int neu1_index = 0;
		int neu2_index = 0;
		// finding the neutralino and photon index
		for(int ind_gen = 0; ind_gen < nGenParticle; ind_gen++)
		{ //this is gen particle loop within event and photon loop
			if ( !foundN1 && gParticleId[ind_gen] == 22 && gParticleMotherId[ind_gen] == 1000022 )
			{ //finds a photon from a neutralino
				pho1index = ind_gen;
				neu1_index = gParticleMotherIndex[ind_gen];
				foundN1 = true;
			}
			else if ( foundN1 && !foundN2 && gParticleId[ind_gen] == 22 && gParticleMotherId[ind_gen] == 1000022 ) 
			{
				pho2index = ind_gen;
				neu2_index = gParticleMotherIndex[ind_gen];
				foundN2 = true;
			}
		}
	
		//bool insideECAL = false;
		//if((gParticleDecayVertexX[neu1_index]*gParticleDecayVertexX[neu1_index]+gParticleDecayVertexY[neu1_index]*gParticleDecayVertexY[neu1_index] < EB_R*EB_R) && abs(gParticleDecayVertexZ[neu1_index])<300.0 && (gParticleDecayVertexX[neu2_index]*gParticleDecayVertexX[neu2_index]+gParticleDecayVertexY[neu2_index]*gParticleDecayVertexY[neu2_index] < EB_R*EB_R) && abs(gParticleDecayVertexZ[neu2_index])<300.0 ) insideECAL = true;
		//if(foundN1==1 && foundN2==1 && insideECAL)
		if(foundN1==1 && foundN2==1)
		{

			float decay_x1 = gParticleDecayVertexX[neu1_index];
			float decay_y1 = gParticleDecayVertexY[neu1_index];
			float decay_z1 = gParticleDecayVertexZ[neu1_index];
			float decay_x2 = gParticleDecayVertexX[neu2_index];
			float decay_y2 = gParticleDecayVertexY[neu2_index];
			float decay_z2 = gParticleDecayVertexZ[neu2_index];

			// need to match up the photon index and the reco photon - this is done based on momentum
			// pho1Pt is reco level, gpho1Pt is gen level information
			float gpho1Pt = gParticlePt[pho1index];
			float gpho2Pt = gParticlePt[pho2index];
			float deltaPt11 = fabs(pho1Pt - gpho1Pt);
			float deltaPt21 = fabs(pho1Pt - gpho2Pt);
			float deltaPt12 = fabs(pho2Pt - gpho1Pt);
			float deltaPt22 = fabs(pho2Pt - gpho2Pt);

			TVector3 genSeed1;
			TVector3 genSeed2;

			float norm1 = pow((pow(gParticlePx[pho1index],2)+pow(gParticlePy[pho1index],2)+pow(gParticlePz[pho1index],2)),0.5);
			float px1 = (gParticlePx[pho1index]) / norm1;
			float py1 = (gParticlePy[pho1index]) / norm1;
			float pz1 = (gParticlePz[pho1index]) / norm1;
			genSeed1 = intersectPoint(decay_x1, decay_y1, decay_z1, px1, py1, pz1, EB_R); // using intersection function written above, radius as 129.7 cm
			float norm2 = pow((pow(gParticlePx[pho2index],2)+pow(gParticlePy[pho2index],2)+pow(gParticlePz[pho2index],2)),0.5);
			float px2 = (gParticlePx[pho2index]) / norm2;
			float py2 = (gParticlePy[pho2index]) / norm2;
			float pz2 = (gParticlePz[pho2index]) / norm2;
			genSeed2 = intersectPoint(decay_x2, decay_y2, decay_z2, px2, py2, pz2, EB_R); // using intersection function written above, radius as 129 cm

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
				pho1GenE = is1To1 ? gParticleE[pho1index] : gParticleE[pho2index];
				pho1GenPt = is1To1 ? gParticlePt[pho1index] : gParticlePt[pho2index];
				pho1GenEta = is1To1 ? gParticleEta[pho1index] : gParticleEta[pho2index];
				pho1GenPhi = is1To1 ? gParticlePhi[pho1index] : gParticlePhi[pho2index];
					
				pho2GenE = is1To1 ? gParticleE[pho2index] : gParticleE[pho1index];
				pho2GenPt = is1To1 ? gParticlePt[pho2index] : gParticlePt[pho1index];
				pho2GenEta = is1To1 ? gParticleEta[pho2index] : gParticleEta[pho1index];
				pho2GenPhi = is1To1 ? gParticlePhi[pho2index] : gParticlePhi[pho1index];
					
				R1 = is1To1 ? pow(decay_x1*decay_x1 + decay_y1*decay_y1, 0.5) : pow(decay_x2*decay_x2 + decay_y2*decay_y2, 0.5) ; 
				R2 = is1To1 ? pow(decay_x2*decay_x2 + decay_y2*decay_y2, 0.5) : pow(decay_x1*decay_x1 + decay_y1*decay_y1, 0.5) ; 
				ZD1 = is1To1 ? decay_z1 : decay_z2 ;
				ZD2 = is1To1 ? decay_z2 : decay_z1 ;

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

				float massNeu = 1000.0;
				//float p_neu1 = is1To1 ? (gParticlePt[neu1_index]*cosh(gParticleEta[neu1_index])) : (gParticlePt[neu2_index]*cosh(gParticleEta[neu2_index]) );
				float p_neu1 = is1To1 ? pow(gParticlePx[neu1_index]*gParticlePx[neu1_index]+gParticlePy[neu1_index]*gParticlePy[neu1_index]+gParticlePz[neu1_index]*gParticlePz[neu1_index],0.5) : pow(gParticlePx[neu2_index]*gParticlePx[neu2_index]+gParticlePy[neu2_index]*gParticlePy[neu2_index]+gParticlePz[neu2_index]*gParticlePz[neu2_index],0.5);
				//float p_neu2 = is1To1 ? (gParticlePt[neu2_index]*cosh(gParticleEta[neu2_index])) : (gParticlePt[neu1_index]*cosh(gParticleEta[neu1_index]) );
				float p_neu2 = is1To1 ? pow(gParticlePx[neu2_index]*gParticlePx[neu2_index]+gParticlePy[neu2_index]*gParticlePy[neu2_index]+gParticlePz[neu2_index]*gParticlePz[neu2_index],0.5) : pow(gParticlePx[neu1_index]*gParticlePx[neu1_index]+gParticlePy[neu1_index]*gParticlePy[neu1_index]+gParticlePz[neu1_index]*gParticlePz[neu1_index],0.5);
			
				TVector3 point_genPV(genVertexX,genVertexY,genVertexZ);	
				TVector3 point_decayV1(is1To1 ? decay_x1: decay_x2,is1To1 ? decay_y1: decay_y2, is1To1 ? decay_z1: decay_z2);
				TVector3 point_decayV2(is1To1 ? decay_x2: decay_x1,is1To1 ? decay_y2: decay_y1, is1To1 ? decay_z2: decay_z1);

				TOF_neu1 = (point_decayV1-point_genPV).Mag() / (SPEED_OF_LIGHT*p_neu1) * pow((pow(massNeu,2) + pow(p_neu1,2)),0.5);
				TOF_neu2 = (point_decayV2-point_genPV).Mag() / (SPEED_OF_LIGHT*p_neu2) * pow((pow(massNeu,2) + pow(p_neu2,2)),0.5);
				TOF_neu1_RF = TOF_neu1*massNeu*pow((pow(massNeu,2) + pow(p_neu1,2)),-0.5);
				TOF_neu2_RF = TOF_neu2*massNeu*pow((pow(massNeu,2) + pow(p_neu2,2)),-0.5);

				TOF_pho1 = (recoSeed1 - point_decayV1).Mag() / SPEED_OF_LIGHT ;
				TOF_pho2 = (recoSeed2 - point_decayV2).Mag() / SPEED_OF_LIGHT ;

				if(abs(genVertexTime) < 100.)
				{
					TOF_total1 = genVertexTime + TOF_neu1 + TOF_pho1 - recoSeed1.Mag() / SPEED_OF_LIGHT;
					TOF_total1_genV = genVertexTime + TOF_neu1 + TOF_pho1 - (recoSeed1 - point_genPV).Mag() / SPEED_OF_LIGHT;
					TOF_total2 = genVertexTime + TOF_neu2 + TOF_pho2 - recoSeed2.Mag() / SPEED_OF_LIGHT;
				TOF_total2_genV = genVertexTime + TOF_neu2 + TOF_pho2 - (recoSeed2 - point_genPV).Mag() / SPEED_OF_LIGHT;
				}
				else
				{
					TOF_total1 = TOF_neu1 + TOF_pho1 - recoSeed1.Mag() / SPEED_OF_LIGHT;
					TOF_total1_genV = TOF_neu1 + TOF_pho1 - (recoSeed1 - point_genPV).Mag() / SPEED_OF_LIGHT;
					TOF_total2 = TOF_neu2 + TOF_pho2 - recoSeed2.Mag() / SPEED_OF_LIGHT;
					TOF_total2_genV = TOF_neu2 + TOF_pho2 - (recoSeed2 - point_genPV).Mag() / SPEED_OF_LIGHT;
				}
				
				pho1angle_xtal = recoSeed1.Angle(recoSeed1 - point_decayV1); 
				pho2angle_xtal = recoSeed2.Angle(recoSeed2 - point_decayV2); 
				
			}//if isMatched
		}//if gen found
	}//if !isData

	outputTree->Fill();		

   }//if nPho>=1
}//event loop

cout << "Writing output trees..." << endl;
outputTree->Write();
NEvents->Write();
outFile->Close();

}//analyzer function


