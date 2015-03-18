#ifndef ControlSampleEvents_H
#define ControlSampleEvents_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"
#include <Rtypes.h>
#include "TLorentzVector.h"

  class ControlSampleEvents {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== Process IDs  ====
      //*******************************************
      enum BkgProcessId { kData = 0,
			  kQCD = 1,
			  kWJets = 2,
			  kZJets = 3,
			  kTTJets = 4,
			  kSingleT = 5,
			  kVV = 6,
			  kHiggs = 7,
			  kSUSY = 99,
			  kUnknown = 999
      };

      /// variables
      Float_t                 weight;
      UInt_t                  run;
      UInt_t                  lumi;
      UInt_t                  event;
      UInt_t                  processID;
      UInt_t                  NPU_0;
      UInt_t                  NPU_Minus1;
      UInt_t                  NPU_Plus1;
      TLorentzVector          genlep1;
      TLorentzVector          genlep2;
      Int_t                   genlep1Type;
      Int_t                   genlep2Type;
      Bool_t                  HLTDecision[100];
      TLorentzVector          lep1;
      TLorentzVector          lep2;
      Int_t                   lep1Type;
      Int_t                   lep2Type;
      Int_t                   lep1MatchedGenLepIndex;
      Int_t                   lep2MatchedGenLepIndex;      
      Bool_t                  lep1PassVeto;
      Bool_t                  lep1PassLoose;
      Bool_t                  lep1PassTight;
      Bool_t                  lep1PassVetoID;
      Bool_t                  lep1PassLooseID;
      Bool_t                  lep1PassTightID;
      Bool_t                  lep1PassVetoIso;
      Bool_t                  lep1PassLooseIso;
      Bool_t                  lep1PassTightIso;
      Float_t                 lep1MinDRToBJet;
      Bool_t                  lep2PassVeto;
      Bool_t                  lep2PassLoose;
      Bool_t                  lep2PassTight;
      Bool_t                  lep2PassVetoID;
      Bool_t                  lep2PassLooseID;
      Bool_t                  lep2PassTightID;
      Bool_t                  lep2PassVetoIso;
      Bool_t                  lep2PassLooseIso;
      Bool_t                  lep2PassTightIso;
      Float_t                 lep2MinDRToBJet;
      TLorentzVector          bjet1;
      TLorentzVector          bjet2;
      Bool_t                  bjet1PassLoose;
      Bool_t                  bjet1PassMedium;
      Bool_t                  bjet1PassTight;
      Bool_t                  bjet2PassLoose;
      Bool_t                  bjet2PassMedium;
      Bool_t                  bjet2PassTight;
      TLorentzVector          jet1;
      TLorentzVector          jet2;      
      Bool_t                  jet1PassCSVLoose;
      Bool_t                  jet1PassCSVMedium;
      Bool_t                  jet1PassCSVTight;
      Bool_t                  jet2PassCSVLoose;
      Bool_t                  jet2PassCSVMedium;
      Bool_t                  jet2PassCSVTight;
      Float_t                 MR;
      Float_t                 Rsq;
      Float_t                 MR_NoDilepton;
      Float_t                 Rsq_NoDilepton;
      Float_t                 MR_NoLeadJet;
      Float_t                 Rsq_NoLeadJet;
      Float_t                 MET;
      Float_t                 MET_NoDilepton;
      Float_t                 MET_NoLeadJet;
      Float_t                 minDPhi;
      Float_t                 minDPhiN;
      Float_t                 dPhiRazor;
      UInt_t                  NJets40;
      UInt_t                  NBJetsLoose;
      UInt_t                  NBJetsMedium;
      UInt_t                  NBJetsTight;
      Float_t                 HT;
      Float_t                 lep1MT;
      Bool_t                  Flag_HBHENoiseFilter;//
      Bool_t                  Flag_CSCTightHaloFilter;
      Bool_t                  Flag_hcalLaserEventFilter; //
      Bool_t                  Flag_EcalDeadCellTriggerPrimitiveFilter; //
      Bool_t                  Flag_goodVertices;
      Bool_t                  Flag_trackingFailureFilter;//
      Bool_t                  Flag_eeBadScFilter; //
      Bool_t                  Flag_ecalLaserCorrFilter;
      Bool_t                  Flag_trkPOGFilters; //
      Bool_t                  Flag_trkPOG_manystripclus53X;
      Bool_t                  Flag_trkPOG_toomanystripclus53X;
      Bool_t                  Flag_trkPOG_logErrorTooManyClusters;
      Bool_t                  Flag_METFilters;
     

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      ControlSampleEvents()  {
	genlep1Ptr  = &genlep1;
	genlep2Ptr  = &genlep2;
	lep1Ptr     = &lep1;
	lep2Ptr     = &lep2;
	bjet1Ptr    = &bjet1;
	bjet2Ptr    = &bjet2;       
	jet1Ptr     = &jet1;
	jet2Ptr     = &jet2;       
      };

      /// default destructor
      ~ControlSampleEvents(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
	weight               = 0.0;
	run                  = 0.0;
	lumi                 = 0.0;
	event                = 0.0;
	processID            = ControlSampleEvents::kUnknown;
	NPU_0                = 0.0;
	NPU_Minus1           = 0.0;
	NPU_Plus1            = 0.0;
	genlep1              = TLorentzVector();
	genlep2              = TLorentzVector();
	genlep1Type          = 0.0;
	genlep2Type          = 0.0;
	for(int i=0;i<100;++i) HLTDecision[i] = false;
	lep1                 = TLorentzVector();
	lep2                 = TLorentzVector();
	lep1Type             = 0.0;
	lep2Type             = 0.0;
	lep1MatchedGenLepIndex = -1;
	lep2MatchedGenLepIndex = -1;
	lep1PassVeto         = 0.0;
	lep1PassLoose        = 0.0;
	lep1PassTight        = 0.0;
	lep1PassVetoID       = 0.0;
	lep1PassLooseID      = 0.0;
	lep1PassTightID      = 0.0;
	lep1PassVetoIso      = 0.0;
	lep1PassLooseIso     = 0.0;
	lep1PassTightIso     = 0.0;
	lep1MinDRToBJet      = 0.0;
	lep2PassVeto         = 0.0;
	lep2PassLoose        = 0.0;
	lep2PassTight        = 0.0;
	lep2PassVetoID       = 0.0;
	lep2PassLooseID      = 0.0;
	lep2PassTightID      = 0.0;
	lep2PassVetoIso      = 0.0;
	lep2PassLooseIso     = 0.0;
	lep2PassTightIso     = 0.0;
	lep2MinDRToBJet      = 0.0;
	bjet1                = TLorentzVector();
	bjet2                = TLorentzVector();
	bjet1PassLoose       = 0.0;
	bjet1PassMedium      = 0.0;
	bjet1PassTight       = 0.0;
	bjet2PassLoose       = 0.0;
	bjet2PassMedium      = 0.0;
	bjet2PassTight       = 0.0;
	jet1                 = TLorentzVector();
	jet2                 = TLorentzVector();
	jet1PassCSVLoose     = 0.0;
	jet1PassCSVMedium    = 0.0;
	jet1PassCSVTight     = 0.0;
	jet2PassCSVLoose     = 0.0;
	jet2PassCSVMedium    = 0.0;
	jet2PassCSVTight     = 0.0;
	MR                   = 0.0;
	Rsq                  = 0.0;
	MR_NoDilepton        = 0.0;
	Rsq_NoDilepton       = 0.0;
	MR_NoLeadJet         = 0.0;
	Rsq_NoLeadJet        = 0.0;
	MET                  = 0.0;
	MET_NoDilepton       = 0.0;
	MET_NoLeadJet        = 0.0;
	minDPhi              = 0.0;
	minDPhiN             = 0.0;
	dPhiRazor            = 0.0;
	NJets40              = 0.0;
	NBJetsLoose          = 0.0;
	NBJetsMedium         = 0.0;
	NBJetsTight          = 0.0;
	HT                   = 0.0;      
	lep1MT               = 0.0;  
	Flag_HBHENoiseFilter = 0.0;//
	Flag_CSCTightHaloFilter = 0.0;
	Flag_hcalLaserEventFilter = 0.0; //
	Flag_EcalDeadCellTriggerPrimitiveFilter = 0.0; //
	Flag_goodVertices = 0.0;
	Flag_trackingFailureFilter = 0.0;//
	Flag_eeBadScFilter = 0.0; //
	Flag_ecalLaserCorrFilter = 0.0;
	Flag_trkPOGFilters = 0.0; //
	Flag_trkPOG_manystripclus53X = 0.0;
	Flag_trkPOG_toomanystripclus53X = 0.0;
	Flag_trkPOG_logErrorTooManyClusters = 0.0;
	Flag_METFilters = 0.0;
      }
    
      /// load a ControlSampleEvents
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("ControlSampleEvent"));
	InitTree();
        assert(tree_);
      }
    
      /// create a ControlSampleEvents
      void CreateTree(){
        tree_ = new TTree("ControlSampleEvent","ControlSampleEvent");
        f_ = 0;

        //book the branches
	tree_->Branch("weight",&weight,"weight/F");
	tree_->Branch("run",&run,"run/i");
	tree_->Branch("lumi",&lumi,"lumi/i");
	tree_->Branch("event",&event,"event/i");
	tree_->Branch("processID",&processID,"processID/i");
	tree_->Branch("NPU_0",&NPU_0,"NPU_0/i");
	tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
	tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
	tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[100]/O");
	tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
	tree_->Branch("genlep2Type",&genlep2Type,"genlep2Type/I");
	tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
	tree_->Branch("lep2Type",&lep2Type,"lep2Type/I");
	tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
	tree_->Branch("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex,"lep2MatchedGenLepIndex/I");
	tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
	tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
	tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
	tree_->Branch("lep1PassVetoID",&lep1PassVetoID,"lep1PassVetoID/O");
	tree_->Branch("lep1PassLooseID",&lep1PassLooseID,"lep1PassLooseID/O");
	tree_->Branch("lep1PassTightID",&lep1PassTightID,"lep1PassTightID/O");
	tree_->Branch("lep1PassVetoIso",&lep1PassVetoIso,"lep1PassVetoIso/O");
	tree_->Branch("lep1PassLooseIso",&lep1PassLooseIso,"lep1PassLooseIso/O");
	tree_->Branch("lep1PassTightIso",&lep1PassTightIso,"lep1PassTightIso/O");
	tree_->Branch("lep1MinDRToBJet",&lep1MinDRToBJet,"lep1MinDRToBJet/F");
	tree_->Branch("lep2PassVeto",&lep2PassVeto,"lep2PassVeto/O");
	tree_->Branch("lep2PassLoose",&lep2PassLoose,"lep2PassLoose/O");
	tree_->Branch("lep2PassTight",&lep2PassTight,"lep2PassTight/O");
	tree_->Branch("lep2PassVetoID",&lep2PassVetoID,"lep2PassVetoID/O");
	tree_->Branch("lep2PassLooseID",&lep2PassLooseID,"lep2PassLooseID/O");
	tree_->Branch("lep2PassTightID",&lep2PassTightID,"lep2PassTightID/O");
	tree_->Branch("lep2PassVetoIso",&lep2PassVetoIso,"lep2PassVetoIso/O");
	tree_->Branch("lep2PassLooseIso",&lep2PassLooseIso,"lep2PassLooseIso/O");
	tree_->Branch("lep2PassTightIso",&lep2PassTightIso,"lep2PassTightIso/O");
	tree_->Branch("lep2MinDRToBJet",&lep2MinDRToBJet,"lep2MinDRToBJet/F");
	tree_->Branch("bjet1PassLoose",&bjet1PassLoose,"bjet1PassLoose/O");
	tree_->Branch("bjet1PassMedium",&bjet1PassMedium,"bjet1PassMedium/O");
	tree_->Branch("bjet1PassTight",&bjet1PassTight,"bjet1PassTight/O");
	tree_->Branch("bjet2PassLoose",&bjet2PassLoose,"bjet2PassLoose/O");
	tree_->Branch("bjet2PassMedium",&bjet2PassMedium,"bjet2PassMedium/O");
	tree_->Branch("bjet2PassTight",&bjet2PassTight,"bjet2PassTight/O");
	tree_->Branch("jet1PassCSVLoose",&jet1PassCSVLoose,"jet1PassCSVLoose/O");
	tree_->Branch("jet1PassCSVMedium",&jet1PassCSVMedium,"jet1PassCSVMedium/O");
	tree_->Branch("jet1PassCSVTight",&jet1PassCSVTight,"jet1PassCSVTight/O");
	tree_->Branch("jet2PassCSVLoose",&jet2PassCSVLoose,"jet2PassCSVLoose/O");
	tree_->Branch("jet2PassCSVMedium",&jet2PassCSVMedium,"jet2PassCSVMedium/O");
	tree_->Branch("jet2PassCSVTight",&jet2PassCSVTight,"jet2PassCSVTight/O");
	tree_->Branch("MR",&MR,"MR/F");
	tree_->Branch("Rsq",&Rsq,"Rsq/F");
	tree_->Branch("MR_NoDilepton",&MR_NoDilepton,"MR_NoDilepton/F");
	tree_->Branch("Rsq_NoDilepton",&Rsq_NoDilepton,"Rsq_NoDilepton/F");
	tree_->Branch("MR_NoLeadJet",&MR_NoLeadJet,"MR_NoLeadJet/F");
	tree_->Branch("Rsq_NoLeadJet",&Rsq_NoLeadJet,"Rsq_NoLeadJet/F");
	tree_->Branch("MET",&MET,"MET/F");
	tree_->Branch("MET_NoDilepton",&MET_NoDilepton,"MET_NoDilepton/F");
	tree_->Branch("MET_NoLeadJet",&MET_NoLeadJet,"MET_NoLeadJet/F");
	tree_->Branch("minDPhi",&minDPhi,"minDPhi/F");
	tree_->Branch("minDPhiN",&minDPhiN,"minDPhiN/F");
	tree_->Branch("dPhiRazor",&dPhiRazor,"dPhiRazor/F");
	tree_->Branch("NJets40",&NJets40,"NJets40/i");
	tree_->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
	tree_->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
	tree_->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
	tree_->Branch("HT",&HT,"HT/F");
	tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
	tree_->Branch("genlep1", "TLorentzVector", &genlep1Ptr);
	tree_->Branch("genlep2", "TLorentzVector", &genlep2Ptr);
	tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
	tree_->Branch("lep2",    "TLorentzVector", &lep2Ptr);
	tree_->Branch("bjet1",   "TLorentzVector", &bjet1Ptr);
	tree_->Branch("bjet2",   "TLorentzVector", &bjet2Ptr);
	tree_->Branch("jet1",    "TLorentzVector", &jet1Ptr);
	tree_->Branch("jet2",    "TLorentzVector", &jet2Ptr);
	tree_->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
	tree_->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter,"Flag_CSCTightHaloFilter/O");
	tree_->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter,"Flag_hcalLaserEventFilter/O");
	tree_->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
	tree_->Branch("Flag_goodVertices", &Flag_goodVertices,"Flag_goodVertices/O");
	tree_->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter,"Flag_trackingFailureFilter/O");
	tree_->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
	tree_->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
	tree_->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters,"Flag_trkPOGFilters/O");
	tree_->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X,"Flag_trkPOG_manystripclus53X/O");
	tree_->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X,"Flag_trkPOG_toomanystripclus53X/O");
	tree_->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters,"Flag_trkPOG_logErrorTooManyClusters/O");
	tree_->Branch("Flag_METFilters", &Flag_METFilters,"Flag_METFilters/O");	
      } 

      // initialze a ControlSampleEvents
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

	tree_->SetBranchAddress("weight",&weight);
	tree_->SetBranchAddress("run",&run);
	tree_->SetBranchAddress("lumi",&lumi);
	tree_->SetBranchAddress("event",&event);
	tree_->SetBranchAddress("processID",&processID);
	tree_->SetBranchAddress("NPU_0",&NPU_0);
	tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
	tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
	tree_->SetBranchAddress("genlep1Type",&genlep1Type);
	tree_->SetBranchAddress("genlep2Type",&genlep2Type);
	tree_->SetBranchAddress("HLTDecision",&HLTDecision);
	tree_->SetBranchAddress("lep1Type",&lep1Type);
	tree_->SetBranchAddress("lep2Type",&lep2Type);
	tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
	tree_->SetBranchAddress("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex);
	tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
	tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
	tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
	tree_->SetBranchAddress("lep1PassVetoID",&lep1PassVetoID);
	tree_->SetBranchAddress("lep1PassLooseID",&lep1PassLooseID);
	tree_->SetBranchAddress("lep1PassTightID",&lep1PassTightID);
	tree_->SetBranchAddress("lep1PassVetoIso",&lep1PassVetoIso);
	tree_->SetBranchAddress("lep1PassLooseIso",&lep1PassLooseIso);
	tree_->SetBranchAddress("lep1PassTightIso",&lep1PassTightIso);
	tree_->SetBranchAddress("lep1MinDRToBJet",&lep1MinDRToBJet);
	tree_->SetBranchAddress("lep2PassVeto",&lep2PassVeto);
	tree_->SetBranchAddress("lep2PassLoose",&lep2PassLoose);
	tree_->SetBranchAddress("lep2PassTight",&lep2PassTight);
	tree_->SetBranchAddress("lep2PassVetoID",&lep2PassVetoID);
	tree_->SetBranchAddress("lep2PassLooseID",&lep2PassLooseID);
	tree_->SetBranchAddress("lep2PassTightID",&lep2PassTightID);
	tree_->SetBranchAddress("lep2PassVetoIso",&lep2PassVetoIso);
	tree_->SetBranchAddress("lep2PassLooseIso",&lep2PassLooseIso);
	tree_->SetBranchAddress("lep2PassTightIso",&lep2PassTightIso);
	tree_->SetBranchAddress("lep2MinDRToBJet",&lep2MinDRToBJet);
	tree_->SetBranchAddress("bjet1PassLoose",&bjet1PassLoose);
	tree_->SetBranchAddress("bjet1PassMedium",&bjet1PassMedium);
	tree_->SetBranchAddress("bjet1PassTight",&bjet1PassTight);
	tree_->SetBranchAddress("bjet2PassLoose",&bjet2PassLoose);
	tree_->SetBranchAddress("bjet2PassMedium",&bjet2PassMedium);
	tree_->SetBranchAddress("bjet2PassTight",&bjet2PassTight);
	tree_->SetBranchAddress("jet1PassCSVLoose",&jet1PassCSVLoose);
	tree_->SetBranchAddress("jet1PassCSVMedium",&jet1PassCSVMedium);
	tree_->SetBranchAddress("jet1PassCSVTight",&jet1PassCSVTight);
	tree_->SetBranchAddress("jet2PassCSVLoose",&jet2PassCSVLoose);
	tree_->SetBranchAddress("jet2PassCSVMedium",&jet2PassCSVMedium);
	tree_->SetBranchAddress("jet2PassCSVTight",&jet2PassCSVTight);
	tree_->SetBranchAddress("MR",&MR);
	tree_->SetBranchAddress("Rsq",&Rsq);
	tree_->SetBranchAddress("MR_NoDilepton",&MR_NoDilepton);
	tree_->SetBranchAddress("Rsq_NoDilepton",&Rsq_NoDilepton);
	tree_->SetBranchAddress("MR_NoLeadJet",&MR_NoLeadJet);
	tree_->SetBranchAddress("Rsq_NoLeadJet",&Rsq_NoLeadJet);
	tree_->SetBranchAddress("MET",&MET);
	tree_->SetBranchAddress("MET_NoDilepton",&MET_NoDilepton);
	tree_->SetBranchAddress("MET_NoLeadJet",&MET_NoLeadJet);
	tree_->SetBranchAddress("minDPhi",&minDPhi);
	tree_->SetBranchAddress("minDPhiN",&minDPhiN);
	tree_->SetBranchAddress("dPhiRazor",&dPhiRazor);
	tree_->SetBranchAddress("NJets40",&NJets40);
	tree_->SetBranchAddress("NBJetsLoose",&NBJetsLoose);
	tree_->SetBranchAddress("NBJetsMedium",&NBJetsMedium);
	tree_->SetBranchAddress("NBJetsTight",&NBJetsTight);
	tree_->SetBranchAddress("HT",&HT);
	tree_->SetBranchAddress("lep1MT",&lep1MT);	
	tree_->SetBranchAddress("genlep1",&genlep1Ptr);
	tree_->SetBranchAddress("genlep2",&genlep2Ptr);
	tree_->SetBranchAddress("lep1",&lep1Ptr);
	tree_->SetBranchAddress("lep2",&lep2Ptr);
	tree_->SetBranchAddress("bjet1",&bjet1Ptr);
	tree_->SetBranchAddress("bjet2",&bjet2Ptr);
	tree_->SetBranchAddress("jet1",&jet1Ptr);
	tree_->SetBranchAddress("jet2",&jet2Ptr);
	tree_->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
	tree_->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter);
	tree_->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter);
	tree_->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
	tree_->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
	tree_->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter);
	tree_->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
	tree_->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter);
	tree_->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters);
	tree_->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X);
	tree_->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X);
	tree_->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters);
	tree_->SetBranchAddress("Flag_METFilters", &Flag_METFilters);
	
        gErrorIgnoreLevel = currentState;
      }

    private:
      TLorentzVector* genlep1Ptr;
      TLorentzVector* genlep2Ptr;
      TLorentzVector* lep1Ptr;
      TLorentzVector* lep2Ptr;
      TLorentzVector* bjet1Ptr;
      TLorentzVector* bjet2Ptr;
      TLorentzVector* jet1Ptr;
      TLorentzVector* jet2Ptr;
      
  }; 


#endif

