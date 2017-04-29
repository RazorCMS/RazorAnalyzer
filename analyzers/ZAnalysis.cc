//LOCAL INCLUDES
#include "ZAnalysis.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"
#include "AngleConversion.h"

using namespace std;


const Double_t MASS_LOW  = 40;
const Double_t MASS_HIGH = 200;
const Double_t PT_CUT    = 22;
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS = 0.105658369;

const Int_t BOSON_ID  = 23;
const Int_t LEPTON_ID = 13;
const int NUM_PDF_WEIGHTS = 60;

enum { eNone=0, eMuMu=1, eEleEle=2 };  // event category enum

//Testing branching and merging
void ZAnalysis::Analyze(bool isData, int option, string outFileName, string label)
{
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);

  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;

  std::cout << "[INFO]: option = " << option << std::endl;
  std::cout << "[INFO]: analysisTag --> " << analysisTag << std::endl;

  
  if ( outFileName.empty() )
    {
      std::cout << "ZAnalysis: Output filename not specified!" << endl << "Using default output name ZAnalysis.root" << std::endl;
      outFileName = "ZAnalysis.root";
    }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outTree = new TTree("ZAnalysis", "Info on selected razor inclusive events");
  
  //Get CMSSW Directory
  // char* cmsswPath;
  // cmsswPath = getenv("CMSSW_BASE");

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, false);
  else if (analysisTag == "Razor2016_80X") helper = new RazorHelper("Razor2016_80X", isData, false);
  else helper = new RazorHelper(analysisTag, isData, false);
  

  //--------------------------------
  //Including Jet Energy Corrections
  //--------------------------------
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();


  //----------
  //pu histo
  //----------
  TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  
  //histogram containing total number of processed events (for normalization)
  TH1F *histNPV = new TH1F("NPV", "NPV", 2, -0.5, 1.5);
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

  //--------------
  //tree variables
  //--------------
  int NPU;
  bool matchGen;
  int category;

  // UInt_t  id_1, id_2;
  // Double_t x_1, x_2, xPDF_1, xPDF_2;
  // Double_t scalePDF, weightPDF;

  TLorentzVector *genV=0;
  float genVPt, genVPhi, genVy, genVMass;
  float genWeight;
  float PUWeight, PUWeightUp, PUWeightDown;
  float weight;

  float met, metPhi, sumEt, u1, u2;
  float tkMet, tkMetPhi, tkSumEt, tkU1, tkU2;
  float mvaMet, mvaMetPhi, mvaSumEt, mvaU1, mvaU2;
  float puppiMet, puppiMetPhi, puppiSumEt, puppiU1, puppiU2;

  int q1, q2;
  TLorentzVector *dilep;
  TLorentzVector *lep1;
  TLorentzVector *lep2;

  float trkIso1, emIso1, hadIso1, pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1; 
  float d01, dz1;
  float muNchi21, nPixHits1, nTkLayers1, nMatch1, nValidHits1, typeBits1;
  TLorentzVector *sta1;

  float trkIso2, emIso2, hadIso2, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;  
  float d02, dz2;
  float muNchi22, nPixHits2, nTkLayers2, nMatch2, nValidHits2, typeBits2;
  TLorentzVector *sta2;

  //PDF SF
  std::vector<float> sf_pdf;
  

  //------------------------
  //set branches on big tree
  //------------------------
 
    outTree->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",     &lumiNum,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",      &eventNum,     "evtNum/i");      // event number
    outTree->Branch("matchGen",    &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    outTree->Branch("category",    &category,   "category/i");    // dilepton category
    // outTree->Branch("id_1",        &id_1,       "id_1/i");        // PDF info -- parton ID for parton 1
    // outTree->Branch("id_2",        &id_2,       "id_2/i");        // PDF info -- parton ID for parton 2
    // outTree->Branch("x_1",         &x_1,        "x_1/d");         // PDF info -- x for parton 1
    // outTree->Branch("x_2",         &x_2,        "x_2/d");         // PDF info -- x for parton 2
    // outTree->Branch("xPDF_1",      &xPDF_1,     "xPDF_1/d");      // PDF info -- x*F for parton 1
    // outTree->Branch("xPDF_2",      &xPDF_2,     "xPDF_2/d");      // PDF info -- x*F for parton 2
    // outTree->Branch("scalePDF",    &scalePDF,   "scalePDF/d");    // PDF info -- energy scale of parton interaction
    // outTree->Branch("weightPDF",   &weightPDF,  "weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("npv",         &nPV,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",         &NPU,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genV",        "TLorentzVector",  &genV);     // GEN boson 4-vector (signal MC)
    outTree->Branch("genVPt",      &genVPt,     "genVPt/F");      // GEN boson pT (signal MC)
    outTree->Branch("genVPhi",     &genVPhi,    "genVPhi/F");     // GEN boson phi (signal MC)
    outTree->Branch("genVy",       &genVy,      "genVy/F");       // GEN boson rapidity (signal MC)
    outTree->Branch("genVMass",    &genVMass,   "genVMass/F");    // GEN boson mass (signal MC)
    outTree->Branch("genWeight",   &genWeight,  "genWeight/F");
    outTree->Branch("PUWeight",    &PUWeight,   "PUWeight/F");
    outTree->Branch("PUWeightUp",  &PUWeightUp, "PUWeightUp/F");
    outTree->Branch("PUWeightDown",  &PUWeightDown, "PUWeightDown/F");
    outTree->Branch("weight",      &weight,     "weight/F");    // event weight per 1/fb (MC)
    outTree->Branch("met",         &met,        "met/F");         // MET
    outTree->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
    outTree->Branch("sumEt",       &sumEt,      "sumEt/F");       // Sum ET
    outTree->Branch("u1",          &u1,         "u1/F");          // parallel component of recoil
    outTree->Branch("u2",          &u2,         "u2/F");          // perpendicular component of recoil
    outTree->Branch("tkMet",       &tkMet,      "tkMet/F");       // MET (track MET)
    outTree->Branch("tkMetPhi",    &tkMetPhi,   "tkMetPhi/F");    // phi(MET) (track MET)
    outTree->Branch("tkSumEt",     &tkSumEt,    "tkSumEt/F");     // Sum ET (track MET)
    outTree->Branch("tkU1",        &tkU1,       "tkU1/F");        // parallel component of recoil (track MET)
    outTree->Branch("tkU2",        &tkU2,       "tkU2/F");        // perpendicular component of recoil (track MET)
    outTree->Branch("mvaMet",      &mvaMet,     "mvaMet/F");      // MVA MET
    outTree->Branch("mvaMetPhi",   &mvaMetPhi,  "mvaMetPhi/F");   // phi(MVA MET)
    outTree->Branch("mvaSumEt",    &mvaSumEt,   "mvaSumEt/F");    // Sum ET (mva MET)
    outTree->Branch("mvaU1",       &mvaU1,      "mvaU1/F");       // parallel component of recoil (mva MET)
    outTree->Branch("mvaU2",       &mvaU2,      "mvaU2/F");       // perpendicular component of recoil (mva MET) 
    outTree->Branch("puppiMet",    &puppiMet,   "puppiMet/F");      // Puppi MET
    outTree->Branch("puppiMetPhi", &puppiMetPhi,"puppiMetPhi/F");   // phi(Puppi MET)
    outTree->Branch("puppiSumEt",  &puppiSumEt, "puppiSumEt/F");    // Sum ET (Puppi MET)
    outTree->Branch("puppiU1",     &puppiU1,    "puppiU1/F");       // parallel component of recoil (Puppi MET)
    outTree->Branch("puppiU2",     &puppiU2,    "puppiU2/F");       // perpendicular component of recoil (Puppi MET)
    outTree->Branch("q1",          &q1,         "q1/I");          // charge of tag lepton
    outTree->Branch("q2",          &q2,         "q2/I");          // charge of probe lepton
    outTree->Branch("dilep",       "TLorentzVector", &dilep);     // di-lepton 4-vector
    outTree->Branch("lep1",        "TLorentzVector", &lep1);      // tag lepton 4-vector
    outTree->Branch("lep2",        "TLorentzVector", &lep2);      // probe lepton 4-vector
    ///// muon specific /////
    outTree->Branch("trkIso1",     &trkIso1,     "trkIso1/F");       // track isolation of tag lepton
    outTree->Branch("trkIso2",     &trkIso2,     "trkIso2/F");       // track isolation of probe lepton
    outTree->Branch("emIso1",      &emIso1,      "emIso1/F");        // ECAL isolation of tag lepton
    outTree->Branch("emIso2",      &emIso2,      "emIso2/F");        // ECAL isolation of probe lepton
    outTree->Branch("hadIso1",     &hadIso1,     "hadIso1/F");       // HCAL isolation of tag lepton
    outTree->Branch("hadIso2",     &hadIso2,     "hadIso2/F");       // HCAL isolation of probe lepton
    outTree->Branch("pfChIso1",    &pfChIso1,    "pfChIso1/F");      // PF charged hadron isolation of tag lepton
    outTree->Branch("pfChIso2",    &pfChIso2,    "pfChIso2/F");      // PF charged hadron isolation of probe lepton
    outTree->Branch("pfGamIso1",   &pfGamIso1,   "pfGamIso1/F");     // PF photon isolation of tag lepton
    outTree->Branch("pfGamIso2",   &pfGamIso2,   "pfGamIso2/F");     // PF photon isolation of probe lepton
    outTree->Branch("pfNeuIso1",   &pfNeuIso1,   "pfNeuIso1/F");     // PF neutral hadron isolation of tag lepton
    outTree->Branch("pfNeuIso2",   &pfNeuIso2,   "pfNeuIso2/F");     // PF neutral hadron isolation of probe lepton
    outTree->Branch("pfCombIso1",  &pfCombIso1,  "pfCombIso1/F");    // PF combined isolation of tag lepton
    outTree->Branch("pfCombIso2",  &pfCombIso2,  "pfCombIso2/F");    // PF combined isolation of probe lepton    
    outTree->Branch("d01",         &d01,         "d01/F");           // transverse impact parameter of tag lepton
    outTree->Branch("d02",         &d02,         "d02/F");           // transverse impact parameter of probe lepton	 
    outTree->Branch("dz1",         &dz1,         "dz1/F");           // longitudinal impact parameter of tag lepton
    outTree->Branch("dz2",         &dz2,         "dz2/F");           // longitudinal impact parameter of probe lepton	 
    outTree->Branch("muNchi21",    &muNchi21,    "muNchi21/F");      // muon fit normalized chi^2 of tag lepton
    outTree->Branch("muNchi22",    &muNchi22,    "muNchi22/F");      // muon fit normalized chi^2 of probe lepton
    outTree->Branch("nPixHits1",   &nPixHits1,	 "nPixHits1/i");     // number of pixel hits of tag muon
    outTree->Branch("nPixHits2",   &nPixHits2,	 "nPixHits2/i");     // number of pixel hits of probe muon
    outTree->Branch("nTkLayers1",  &nTkLayers1,  "nTkLayers1/i");    // number of tracker layers of tag muon
    outTree->Branch("nTkLayers2",  &nTkLayers2,  "nTkLayers2/i");    // number of tracker layers of probe muon
    outTree->Branch("nMatch1",     &nMatch1,	 "nMatch1/i");       // number of matched segments of tag muon
    outTree->Branch("nMatch2",     &nMatch2,	 "nMatch2/i");       // number of matched segments of probe muon 
    outTree->Branch("nValidHits1", &nValidHits1, "nValidHits1/i");   // number of valid muon hits of tag muon
    outTree->Branch("nValidHits2", &nValidHits2, "nValidHits2/i");   // number of valid muon hits of probe muon
    outTree->Branch("typeBits1",   &typeBits1,   "typeBits1/i");     // muon type of tag muon
    outTree->Branch("typeBits2",   &typeBits2,   "typeBits2/i");     // muon type of probe muon
    outTree->Branch("sta1",        "TLorentzVector", &sta1);         // tag standalone muon 4-vector
    outTree->Branch("sta2",        "TLorentzVector", &sta2);         // probe standalone muon 4-vector


  //begin loop
  if ( fChain == 0 ) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::cout << "[INFO]: Total Entries = " << fChain->GetEntries() << "\n";
  for ( Long64_t jentry=0; jentry < nentries; jentry++ ) {
      //begin event
      if( jentry % 10000 == 0 ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      Long64_t ientry = LoadTree( jentry );
      if ( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
    
      //fill normalization histogram    
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      weight = genWeight;
      SumWeights->Fill(1.0, weight);
      
      //reset tree variables
      runNum = -999;
      lumiNum = -999;
      eventNum = -999;
      matchGen = -999;
      category = -999;
      // id_1 = -999;
      // id_2 = -999;
      // x_1 = -999;
      // x_2 = -999;
      // xPDF_1 = -999;
      // xPDF_2 = -999;
      // scalePDF = -999;
      // weightPDF = -999;
      nPV = -999;
      NPU = -999;
      genV = 0;
      genVPt = -999;
      genVPhi = -999;
      genVy = -999;
      genVMass = -999;
      genWeight = -999;
      PUWeight = -999;
      weight = 0;
      met = -999;
      metPhi = -999;
      sumEt = -999;
      u1 = -999;
      u2 = -999;
      tkMet = -999;
      tkMetPhi = -999;
      tkSumEt = -999;
      tkU1 = -999;
      tkU2 = -999;
      mvaMet = -999;
      mvaMetPhi = -999;
      mvaSumEt = -999;
      mvaU1 = -999;
      mvaU2 = -999;
      puppiMet = -999;
      puppiMetPhi = -999;
      puppiSumEt = -999;
      puppiU1 = -999;
      puppiU2 = -999;
      q1 = -999;
      q2 = -999;
      dilep = 0;
      lep1 = 0;
      lep2 = 0;
      trkIso1 = -999;
      trkIso2 = -999;
      emIso1 = -999;
      emIso2 = -999;
      hadIso1 = -999;
      hadIso2 = -999;
      pfChIso1 = -999;
      pfChIso2 = -999;
      pfGamIso1 = -999;
      pfGamIso2 = -999;
      pfNeuIso1 = -999;
      pfNeuIso2 = -999;
      pfCombIso1 = -999;
      pfCombIso2 = -999;
      d01 = -999;
      d02 = -999;
      dz1 = -999;
      dz2 = -999;
      muNchi21 = -999;
      muNchi22 = -999;
      nPixHits1 = -999;
      nPixHits2 = -999;
      nTkLayers1 = -999;
      nTkLayers2 = -999;
      nMatch1 = -999;
      nMatch2 = -999;
      nValidHits1 = -999;
      nValidHits2 = -999;
      typeBits1 = -999;
      typeBits2 = -999;
      sta1 = 0;
      sta2 = 0;

 

      //------------------
      //Pileup reweighting
      //------------------
      PUWeight = 1.0;
      if( !isData ) {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
	puhisto->Fill(NPU);
	PUWeight = helper->getPileupWeight(NPU);
	PUWeightUp = helper->getPileupWeightUp(NPU) / PUWeight;
	PUWeightDown = helper->getPileupWeightDown(NPU) / PUWeight;	
      }
      
      /////////////////////////////////
      //Scale and PDF variations
      /////////////////////////////////

      if ( (*scaleWeights).size() >= 9 ) 
	{
	  // sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
	  // sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
	  // sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
	  // sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
	  // sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
	  // sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;

	  
	  SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
	  SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
	  SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
	  SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
	  SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
	  SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
	}
      
      sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
      for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) {
	sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
      }
      

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************

      //-------------------------------
      //1) Look for Zmm Candidate
      //-------------------------------
      int lep1Type = 0;
      double lep1Pt = 0;
      double lep1Eta = -999;
      double lep1Phi = -999;
      int lep2Type = 0;
      double lep2Pt = 0;
      double lep2Eta = -999;
      double lep2Phi = -999;

      double bestDiMuonMass = 0;
      TLorentzVector ZCandidate;
      for( int i = 0; i < nMuons; i++ )	{
	  if(!isMuonPOGLooseMuon(i)) continue;  
	  if(muonPt[i] < 15) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
	  for( int j = i+1; j < nMuons; j++ )	{
	    if(!isVetoMuon(j)) continue;  
	    if(muonPt[j] < 15) continue;
	    if(abs(muonEta[j]) > 2.4) continue;
	    
	    TLorentzVector tmpMuon1;
	    tmpMuon1.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i],0.1057);
	    TLorentzVector tmpMuon2;
	    tmpMuon2.SetPtEtaPhiM(muonPt[j],muonEta[j], muonPhi[j],0.1057);
	    double tmpMass = (tmpMuon1+tmpMuon2).M();	    


	    if ( fabs(tmpMass - 91.2) < fabs(bestDiMuonMass - 91.2))  {
	      bestDiMuonMass = tmpMass;
	      category = eMuMu;
	      lep1Type = 13 * -1 * muonCharge[i];
	      lep1Pt = muonPt[i];
	      lep1Eta = muonEta[i];
	      lep1Phi = muonPhi[i];
	      lep2Type = 13 * -1 * muonCharge[j];
	      lep2Pt = muonPt[j];
	      lep2Eta = muonEta[j];
	      lep2Phi = muonPhi[j];
	      ZCandidate = tmpMuon1 + tmpMuon2;
	      matchGen = matchesGenMuon(lep1Eta,lep1Phi) && matchesGenMuon(lep2Eta,lep2Phi);	      
	    } //if better Z mass match
	  }// loop 2nd muon
      } //loop 1st muon


      //-------------------------------
      //2) Look for Zee Candidate
      //-------------------------------
      if (category == eNone) {
	double bestDielectronMass = 0;
	for( int i = 0; i < nElectrons; i++ )	{
	  if(!isVetoElectron(i)) continue;  
	  if(elePt[i] < 20) continue;
	  if(abs(eleEta[i]) > 2.4) continue;
	  for( int j = i+1; j < nElectrons; j++ )	{
	    if(!isVetoElectron(j)) continue;  
	    if(elePt[j] < 20) continue;
	    if(abs(eleEta[j]) > 2.4) continue;
	    
	    TLorentzVector tmpElectron1;
	    tmpElectron1.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i],0.000511);
	    TLorentzVector tmpElectron2;
	    tmpElectron2.SetPtEtaPhiM(elePt[j],eleEta[j], elePhi[j],0.000511);
	    double tmpMass = (tmpElectron1+tmpElectron2).M();	    


	    if ( fabs(tmpMass - 91.2) < fabs(bestDielectronMass - 91.2)){
	      bestDielectronMass = tmpMass;
	      category = eEleEle;
	      lep1Type = 11 * -1 * eleCharge[i];
	      lep1Pt = elePt[i];
	      lep1Eta = eleEta[i];
	      lep1Phi = elePhi[i];
	      lep2Type = 11 * -1 * eleCharge[j];
	      lep2Pt = elePt[j];
	      lep2Eta = eleEta[j];
	      lep2Phi = elePhi[j];
	      ZCandidate = tmpElectron1 + tmpElectron2;
	      matchGen = matchesGenElectron(lep1Eta,lep1Phi) && matchesGenElectron(lep2Eta,lep2Phi);	      
	    } //if better Z mass match
	  } //loop 2nd electron
	} //loop 1st electron
      }

      lep1 = new TLorentzVector;
      lep2 = new TLorentzVector;
      dilep = new TLorentzVector;
      if ( abs(lep1Type) == 11) {
	lep1->SetPtEtaPhiM(lep1Pt, lep1Eta, lep1Phi,0.000511);
	lep2->SetPtEtaPhiM(lep2Pt, lep2Eta, lep2Phi,0.000511);
      } else {
	lep1->SetPtEtaPhiM(lep1Pt, lep1Eta, lep1Phi,0.1057);
	lep2->SetPtEtaPhiM(lep2Pt, lep2Eta, lep2Phi,0.1057);
      }
      dilep->SetPtEtaPhiM( (*lep1+*lep2).Pt(), (*lep1+*lep2).Eta(), (*lep1+*lep2).Phi(), (*lep1+*lep2).M());
      q1 = -1 * lep1Type / abs(lep1Type);
      q2 = -1 * lep2Type / abs(lep2Type);

      //*************************************************************************
      //MET Variables
      //*************************************************************************
      met = metType1Pt;
      metPhi = metType1Phi;

 
      //******************************************************
      //Filters
      //******************************************************
      //2-L filter
      if (category == eNone) continue;

      //Fill Event
      outTree->Fill();

      //end of event loop
  }
  
  std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;
  
  std::cout << "[INFO]: Writing output trees..." << std::endl;    
  outFile->cd();
  outTree->Write();
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  histNPV->Write();
  puhisto->Write();
  
  outFile->Close();
  delete helper;

}
